#!/usr/bin/env bash
# Local dual-GPU equilibrium IR smoke at aging_weak_lambda couplings.
# Constant coupling + finite-q + velocity thermalization from shared IC library.
set -euo pipefail

REPO="${REPO:-/scratch/mh7373/projects/cav-hoomd}"
PY="${PY:-python}"
SCRIPT="$REPO/examples/05_advanced_run.py"
BASE="$REPO/aging_weak_lambda_ir"
IC="$REPO/aging_weak_lambda_ic_library/init-500-from-ic-library.gsd"
LOGDIR="$BASE/logs/ir_smoke"
QUEUE="$LOGDIR/ir_queue.tsv"
LOCKDIR="$LOGDIR/queue_locks"

RUNTIME_PS="${RUNTIME_PS:-100.0}"
N_REPLICAS="${N_REPLICAS:-2}"
GPUS="${GPUS:-0,1}"
PER_GPU="${PER_GPU:-1}"
DIPOLE_PERIOD_PS="${DIPOLE_PERIOD_PS:-0.002}"
DRY_RUN="${DRY_RUN:-0}"

mkdir -p "$LOGDIR" "$LOCKDIR"

if [[ ! -f "$IC" ]]; then
  echo "ERROR: missing cavity IC: $IC" >&2
  exit 1
fi

"$PY" - "$BASE" "$N_REPLICAS" "$QUEUE" <<'PY'
import sys
from pathlib import Path

base = Path(sys.argv[1])
n_replicas = int(sys.argv[2])
queue_path = Path(sys.argv[3])

lambdas = [
    (0.0, "0"),
    (0.01, "0p01"),
    (0.016667, "0p016667"),
    (0.023333, "0p023333"),
    (0.03, "0p03"),
]

rows = []
for lam, tag in lambdas:
    lam_dir = base / f"lambda{tag}"
    lam_dir.mkdir(parents=True, exist_ok=True)
    for replica in range(n_replicas):
        rows.append(f"{lam}\t{tag}\t{replica}\n")

queue_path.write_text("".join(rows), encoding="utf-8")
print(f"queued {len(rows)} jobs (lambda-major, replicas 0..{n_replicas - 1})")
PY

njobs=$(wc -l <"$QUEUE" | tr -d ' ')
echo "$(date -u +%H:%M:%SZ) IR smoke queue size=$njobs GPUS=$GPUS runtime=${RUNTIME_PS}ps dipole=${DIPOLE_PERIOD_PS}ps"

if [[ "$DRY_RUN" == "1" ]]; then
  echo "DRY_RUN=1: queue written to $QUEUE"
  cat "$QUEUE"
  exit 0
fi

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
  local cwd
  cwd="$BASE/lambda${tag}"
  echo "$(date -u +%H:%M:%SZ) start gpu=$gpu lam=$lam replica=$replica frame=$replica"
  (
    cd "$cwd"
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
      --replicas "$replica" --seed $((replica + 1)) \
      --electrostatics reaction_field \
      --eps-rf 80.0 --coulomb-rcut 16.0
  ) >"$LOGDIR/ir_lam${tag}_r${replica}_gpu${gpu}.out" \
    2>"$LOGDIR/ir_lam${tag}_r${replica}_gpu${gpu}.err"
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
echo "$(date -u +%H:%M:%SZ) IR smoke launcher finished"
