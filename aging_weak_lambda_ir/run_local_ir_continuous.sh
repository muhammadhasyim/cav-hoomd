#!/usr/bin/env bash
# Continuous dual-GPU equilibrium IR production (400 ps, one MD job per GPU).
# Each finished run is registered for incremental DACF averaging by ir_continuous_daemon.py.
set -euo pipefail

REPO="${REPO:-/scratch/mh7373/projects/cav-hoomd}"
PY="${PY:-python}"
SCRIPT="$REPO/examples/05_advanced_run.py"
BASE="$REPO/aging_weak_lambda_ir"
IC="$REPO/aging_weak_lambda_ic_library/init-500-from-ic-library.gsd"
LOGDIR="$BASE/logs/ir_continuous"
QUEUE="$LOGDIR/ir_queue.tsv"
LOCKDIR="$LOGDIR/queue_locks"
REGISTRY="$LOGDIR/active_runs.jsonl"

RUNTIME_PS="${RUNTIME_PS:-400.0}"
N_REPLICAS="${N_REPLICAS:-2}"
GPUS="${GPUS:-0,1}"
PER_GPU="${PER_GPU:-1}"
DIPOLE_PERIOD_PS="${DIPOLE_PERIOD_PS:-0.002}"
CYCLES="${CYCLES:-0}"
DRY_RUN="${DRY_RUN:-0}"

mkdir -p "$LOGDIR" "$LOCKDIR" "$BASE/derived/continuous"

if [[ ! -f "$IC" ]]; then
  echo "ERROR: missing cavity IC: $IC" >&2
  exit 1
fi

write_queue() {
  local cycles="$1"
  "$PY" - "$BASE" "$N_REPLICAS" "$QUEUE" "$cycles" <<'PY'
import sys
from pathlib import Path

base = Path(sys.argv[1])
n_replicas = int(sys.argv[2])
queue_path = Path(sys.argv[3])
cycles = int(sys.argv[4])

lambdas = [
    (0.0, "0"),
    (0.01, "0p01"),
    (0.016667, "0p016667"),
    (0.023333, "0p023333"),
    (0.03, "0p03"),
]

rows = []
cycle_count = max(cycles, 1)
for _cycle in range(cycle_count):
    for lam, tag in lambdas:
        lam_dir = base / f"lambda{tag}"
        lam_dir.mkdir(parents=True, exist_ok=True)
        (lam_dir / "runs").mkdir(parents=True, exist_ok=True)
        for replica in range(n_replicas):
            rows.append(f"{lam}\t{tag}\t{replica}\t{_cycle}\n")

queue_path.write_text("".join(rows), encoding="utf-8")
print(f"queued {len(rows)} jobs ({cycle_count} cycle(s), replicas 0..{n_replicas - 1})")
PY
}

if [[ "$CYCLES" == "0" ]]; then
  write_queue 1
else
  write_queue "$CYCLES"
fi

njobs=$(wc -l <"$QUEUE" | tr -d ' ')
echo "$(date -u +%H:%M:%SZ) IR continuous queue size=$njobs GPUS=$GPUS runtime=${RUNTIME_PS}ps dipole=${DIPOLE_PERIOD_PS}ps"

if [[ "$DRY_RUN" == "1" ]]; then
  echo "DRY_RUN=1: queue written to $QUEUE"
  head -n 5 "$QUEUE"
  exit 0
fi

register_run() {
  local run_id="$1" run_dir="$2" lam="$3" tag="$4" replica="$5" gpu="$6"
  "$PY" - "$REGISTRY" "$run_id" "$run_dir" "$lam" "$tag" "$replica" "$gpu" <<'PY'
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

registry = Path(sys.argv[1])
payload = {
    "run_id": sys.argv[2],
    "run_dir": sys.argv[3],
    "lam": float(sys.argv[4]),
    "tag": sys.argv[5],
    "replica": int(sys.argv[6]),
    "gpu": int(sys.argv[7]),
    "started_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    "completed": False,
}
registry.parent.mkdir(parents=True, exist_ok=True)
with registry.open("a", encoding="utf-8") as handle:
    handle.write(json.dumps(payload) + "\n")
PY
}

mark_completed() {
  local run_id="$1"
  "$PY" - "$REGISTRY" "$run_id" <<'PY'
import json
import sys
from pathlib import Path

registry = Path(sys.argv[1])
run_id = sys.argv[2]
rows = []
if registry.is_file():
    for line in registry.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        row = json.loads(line)
        if row["run_id"] == run_id:
            row["completed"] = True
        rows.append(row)
registry.write_text(
    "".join(json.dumps(row) + "\n" for row in rows),
    encoding="utf-8",
)
PY
}

claim_next() {
  (
    flock 9
    if [[ ! -s "$QUEUE" ]]; then
      write_queue 1 >/dev/null
    fi
    while IFS=$'\t' read -r lam tag replica cycle; do
      claim="$LOCKDIR/claim_${tag}_${replica}_${cycle}.lock"
      [[ -e "$claim" ]] && continue
      grep -v $'^'"${lam}"$'\t'"${tag}"$'\t'"${replica}"$'\t'"${cycle}"'$' "$QUEUE" >"${QUEUE}.tmp" || true
      mv "${QUEUE}.tmp" "$QUEUE"
      printf '%s\t%s\t%s\t%s\n' "$lam" "$tag" "$replica" "$cycle" >"$claim"
      echo "$lam $tag $replica $cycle"
      exit 0
    done <"$QUEUE"
    exit 1
  ) 9>"$LOCKDIR/queue.flock"
}

run_one() {
  local gpu="$1" lam="$2" tag="$3" replica="$4" cycle="$5"
  local run_id run_dir
  run_id="$(date -u +%Y%m%dT%H%M%SZ)_r${replica}_c${cycle}_gpu${gpu}"
  run_dir="$BASE/lambda${tag}/runs/${run_id}"
  mkdir -p "$run_dir"
  register_run "$run_id" "$run_dir" "$lam" "$tag" "$replica" "$gpu"
  echo "$(date -u +%H:%M:%SZ) start gpu=$gpu lam=$lam replica=$replica cycle=$cycle run_id=$run_id"
  (
    cd "$run_dir"
    CUDA_VISIBLE_DEVICES="$gpu" "$PY" "$SCRIPT" \
      --molecular-bath bussi --cavity-bath langevin \
      --finite-q \
      --lambda-coupling "$lam" --coupling-type constant \
      --temperature 100.0 --frequency 1560.0 \
      --runtime "$RUNTIME_PS" \
      --input-gsd "$IC" --frame "$replica" \
      --device GPU --gpu-id 0 \
      --fixed-timestep --timestep 1.0 \
      --enable-dipole-fdr \
      --dipole-fdr-output-period-ps "$DIPOLE_PERIOD_PS" \
      --dipole-fdr-max-correlation-time-ps "$RUNTIME_PS" \
      --hdf5-output-period-ps "$DIPOLE_PERIOD_PS" \
      --enable-energy-tracker --energy-output-period-ps 2.0 \
      --disable-gsd --disable-temp-tracker \
      --console-output-period-ps 20.0 \
      --replicas 0 --seed $((replica + 1 + cycle * 1000)) \
      --electrostatics reaction_field \
      --eps-rf 80.0 --coulomb-rcut 16.0
  ) >"$LOGDIR/${run_id}.out" 2>"$LOGDIR/${run_id}.err"
  local rc=$?
  rm -f "$LOCKDIR/claim_${tag}_${replica}_${cycle}.lock"
  if [[ $rc -eq 0 ]]; then
    mark_completed "$run_id"
  fi
  echo "$(date -u +%H:%M:%SZ) done gpu=$gpu lam=$lam replica=$replica cycle=$cycle run_id=$run_id rc=$rc"
  return $rc
}

worker() {
  local gpu="$1"
  while true; do
    local job
    if ! job=$(claim_next); then
      echo "$(date -u +%H:%M:%SZ) gpu=$gpu queue empty after refill attempt; sleeping 30s"
      sleep 30
      continue
    fi
    read -r lam tag replica cycle <<<"$job"
    run_one "$gpu" "$lam" "$tag" "$replica" "$cycle" || true
  done
}

IFS=',' read -ra GPU_ARR <<<"$GPUS"
pids=()
for gpu in "${GPU_ARR[@]}"; do
  for ((slot=0; slot<PER_GPU; slot++)); do
    worker "$gpu" >"$LOGDIR/worker_gpu${gpu}_slot${slot}.log" 2>&1 &
    pids+=($!)
  done
done

printf '%s\n' "${pids[@]}" >"$LOGDIR/worker.pids"
echo "started ${#pids[@]} continuous workers on GPUs ${GPUS}: ${pids[*]}"
wait "${pids[@]}" || true
