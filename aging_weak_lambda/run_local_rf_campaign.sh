#!/usr/bin/env bash
# Local 4-GPU aging_weak_lambda production with reaction-field electrostatics.
# Queue order: replica-major, lambda ascending (r0 all λ, then r1 all λ, ...).
set -euo pipefail

REPO=/workspace/refs/cav-hoomd
PY=python3
SCRIPT=$REPO/examples/05_advanced_run.py
BASE=$REPO/aging_weak_lambda
IC=$REPO/aging_weak_lambda_ic_library/init-500-from-ic-library.gsd
LOGDIR=$BASE/logs/rf_campaign
QUEUE=$LOGDIR/replica_queue.tsv
LOCKDIR=$LOGDIR/queue_locks
RUNTIME_PS=${RUNTIME_PS:-1500.0}
SWITCH_PS=${SWITCH_PS:-200.0}
N_REPLICAS=${N_REPLICAS:-500}
PER_GPU=${PER_GPU:-8}
GPUS=${GPUS:-0,1,2,3}
MIN_H5_BYTES=${MIN_H5_BYTES:-500000000}
FKT_MAX_REFS=${FKT_MAX_REFS:-8}
SAMPLE_PS=${SAMPLE_PS:-2.0}

export PYTHONPATH=/opt/hoomd-md:${PYTHONPATH:-}

mkdir -p "$LOGDIR" "$LOCKDIR"

if [[ ! -f "$IC" ]]; then
  echo "ERROR: missing cavity IC: $IC" >&2
  exit 1
fi

"$PY" - "$BASE" "$N_REPLICAS" "$MIN_H5_BYTES" "$QUEUE" <<'PY'
import sys
from pathlib import Path

base = Path(sys.argv[1])
n_replicas = int(sys.argv[2])
min_bytes = int(sys.argv[3])
queue_path = Path(sys.argv[4])

lambdas = [
    (0.0, "0"),
    (0.01, "0p01"),
    (0.016667, "0p016667"),
    (0.023333, "0p023333"),
    (0.03, "0p03"),
]

def is_complete(lam_dir: Path, replica: int) -> bool:
    if not lam_dir.is_dir():
        return False
    patterns = ("cavity_coupling_*", "coupling_*")
    for pattern in patterns:
        for run_dir in lam_dir.glob(pattern):
            h5 = run_dir / f"observables_replica_{replica}.h5"
            fkt_last = run_dir / f"prod-{replica}_fkt_ref_007.txt"
            if (
                h5.is_file()
                and h5.stat().st_size >= min_bytes
                and fkt_last.is_file()
            ):
                return True
    return False

rows = []
# Replica-major, lambda ascending
for replica in range(n_replicas):
    for lam, tag in lambdas:
        lam_dir = base / f"lambda{tag}"
        lam_dir.mkdir(parents=True, exist_ok=True)
        if not is_complete(lam_dir, replica):
            rows.append(f"{lam}\t{tag}\t{replica}\n")

queue_path.write_text("".join(rows), encoding="utf-8")
print(f"queued {len(rows)} jobs (replica-major, lambda ascending)")
PY

njobs=$(wc -l <"$QUEUE" | tr -d ' ')
echo "$(date -u +%H:%M:%SZ) RF campaign queue size=$njobs PER_GPU=$PER_GPU GPUS=$GPUS runtime=${RUNTIME_PS}ps sample=${SAMPLE_PS}ps"

random_frame() {
  "$PY" - "$1" "$2" <<'PY'
import sys
replica = int(sys.argv[1])
lam = float(sys.argv[2])
print((replica * 10007 + int(lam * 1e6)) % 500)
PY
}

claim_next() {
  (
    flock 9
    while IFS=$'\t' read -r lam tag replica; do
      claim="$LOCKDIR/claim_${tag}_${replica}.lock"
      [[ -e "$claim" ]] && continue
      grep -v $'^'"${lam}"$'\t'"${tag}"$'\t'"${replica}"'$' "$QUEUE" >"${QUEUE}.tmp" || true
      mv "${QUEUE}.tmp" "$QUEUE"
      printf '%s\t%s\t%s\n' "$lam" "$tag" "$replica" >"$claim"
      echo "$lam $tag $replica"
      exit 0
    done <"$QUEUE"
    exit 1
  ) 9>"$LOCKDIR/queue.flock"
}

run_one() {
  local gpu="$1" lam="$2" tag="$3" replica="$4"
  local cwd frame
  cwd="$BASE/lambda${tag}"
  frame=$(random_frame "$replica" "$lam")
  echo "$(date -u +%H:%M:%SZ) start gpu=$gpu lam=$lam replica=$replica frame=$frame"
  (
    cd "$cwd"
    CUDA_VISIBLE_DEVICES="$gpu" "$PY" "$SCRIPT" \
      --molecular-bath bussi --cavity-bath langevin \
      --lambda-coupling "$lam" --coupling-type step \
      --temperature 100.0 --frequency 1560.0 \
      --runtime "$RUNTIME_PS" --switch-time "$SWITCH_PS" \
      --input-gsd "$IC" --frame "$frame" \
      --device GPU --gpu-id 0 \
      --fixed-timestep --timestep 1.0 \
      --enable-energy-tracker --energy-output-period-ps "$SAMPLE_PS" \
      --hdf5-output-period-ps "$SAMPLE_PS" \
      --enable-fkt --fkt-kmag 6.02 --fkt-wavevectors 50 \
      --fkt-ref-interval 200.0 --fkt-max-refs "$FKT_MAX_REFS" \
      --fkt-output-period-ps "$SAMPLE_PS" \
      --disable-gsd --disable-temp-tracker \
      --console-output-period-ps 100.0 \
      --replicas "$replica" --seed $((replica + 1)) \
      --electrostatics reaction_field \
      --eps-rf 80.0 --coulomb-rcut 16.0
  ) >"$LOGDIR/rf_lam${tag}_r${replica}_gpu${gpu}.out" \
    2>"$LOGDIR/rf_lam${tag}_r${replica}_gpu${gpu}.err"
  local rc=$?
  rm -f "$LOCKDIR/claim_${tag}_${replica}.lock"
  echo "$(date -u +%H:%M:%SZ) done gpu=$gpu lam=$lam replica=$replica rc=$rc"
  return $rc
}

worker() {
  local gpu="$1"
  while true; do
    local job
    if ! job=$(claim_next); then
      echo "$(date -u +%H:%M:%SZ) gpu=$gpu queue empty; exit"
      break
    fi
    read -r lam tag replica <<<"$job"
    run_one "$gpu" "$lam" "$tag" "$replica" || true
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
echo "started ${#pids[@]} workers on GPUs ${GPUS}: ${pids[*]}"
wait "${pids[@]}" || true
echo "$(date -u +%H:%M:%SZ) RF campaign launcher finished"
