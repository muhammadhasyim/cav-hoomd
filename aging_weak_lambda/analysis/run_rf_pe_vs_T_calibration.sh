#!/usr/bin/env bash
# Reaction-field PE-vs-T calibration on local GPUs (one MD job per GPU).
set -euo pipefail

REPO="${REPO:-/scratch/mh7373/projects/cav-hoomd}"
PY="${PY:-/scratch/mh7373/miniforge3/envs/hoomd/bin/python}"
SCRIPT="$REPO/examples/05_advanced_run.py"
ANALYZE="$REPO/aging_weak_lambda/analysis/analyze_rf_pe_vs_T_calibration.py"
OUT="${OUT:-$REPO/aging_weak_lambda/analysis/pe_vs_T_calib_rf}"
IC="$REPO/aging_weak_lambda_ic_library/init-500-from-ic-library.gsd"
LOGDIR="$OUT/logs"
QUEUE="$LOGDIR/temperature_queue.txt"
LOCKDIR="$LOGDIR/locks"
RUNTIME_PS="${RUNTIME_PS:-2000.0}"
DISCARD_PS="${DISCARD_PS:-1000.0}"
GPUS="${GPUS:-0,1}"
PER_GPU="${PER_GPU:-1}"
MIN_END_PS="${MIN_END_PS:-1999.0}"

TEMPERATURES=(300 280 260 240 220 200 180 160 140 120 100 80)

mkdir -p "$OUT" "$LOGDIR" "$LOCKDIR"

if [[ ! -f "$IC" ]]; then
  echo "ERROR: missing init GSD: $IC" >&2
  exit 1
fi

is_complete() {
  local temp_k="$1"
  local h5="$OUT/temp_${temp_k}K/no_cavity/observables_replica_0.h5"
  [[ -f "$h5" ]] || return 1
  "$PY" - "$h5" "$MIN_END_PS" <<'PY'
import sys
from pathlib import Path
import h5py
import numpy as np

path = Path(sys.argv[1])
min_end = float(sys.argv[2])
with h5py.File(path, "r") as handle:
    time = handle["time"][:]
valid = time[time > 0.0]
if valid.size == 0:
    raise SystemExit(1)
if float(valid[-1]) < min_end:
    raise SystemExit(1)
PY
}

build_queue() {
  : >"$QUEUE"
  local temp_k
  for temp_k in "${TEMPERATURES[@]}"; do
    if is_complete "$temp_k"; then
      echo "$(date -u +%H:%M:%SZ) skip complete temp_${temp_k}K"
      continue
    fi
    printf '%s\n' "$temp_k" >>"$QUEUE"
  done
  echo "$(date -u +%H:%M:%SZ) queue size=$(wc -l <"$QUEUE" | tr -d ' ')"
}

claim_next() {
  (
    flock 9
    local temp_k
    while IFS= read -r temp_k; do
      [[ -n "$temp_k" ]] || continue
      local claim="$LOCKDIR/claim_${temp_k}K.lock"
      [[ -e "$claim" ]] && continue
      grep -vx "$temp_k" "$QUEUE" >"${QUEUE}.tmp" || true
      mv "${QUEUE}.tmp" "$QUEUE"
      touch "$claim"
      echo "$temp_k"
      exit 0
    done <"$QUEUE"
    exit 1
  ) 9>"$LOCKDIR/queue.flock"
}

run_one() {
  local gpu="$1"
  local temp_k="$2"
  local run_dir="$OUT/temp_${temp_k}K"
  local log="$LOGDIR/temp_${temp_k}K_gpu${gpu}.log"
  mkdir -p "$run_dir"
  {
    echo "START T=${temp_k}K gpu=${gpu} host=$(hostname) date=$(date -Iseconds)"
    nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv,noheader || true
    cd "$run_dir"
    CUDA_VISIBLE_DEVICES="$gpu" "$PY" "$SCRIPT" \
      --molecular-bath bussi --no-cavity \
      --temperature "$temp_k" --runtime "$RUNTIME_PS" \
      --input-gsd "$IC" --frame -1 \
      --device GPU --gpu-id 0 \
      --fixed-timestep --timestep 1.0 \
      --enable-energy-tracker --energy-output-period-ps 1.0 \
      --hdf5-output-period-ps 0.01 \
      --disable-gsd --disable-temp-tracker \
      --console-output-period-ps 100 \
      --replicas 0 --seed "$((1000 + temp_k))" \
      --electrostatics reaction_field --eps-rf 80.0 --coulomb-rcut 16.0
    echo "END T=${temp_k}K gpu=${gpu} rc=$? date=$(date -Iseconds)"
  } >"$log" 2>&1
  local rc=$?
  rm -f "$LOCKDIR/claim_${temp_k}K.lock"
  return "$rc"
}

worker() {
  local gpu="$1"
  while true; do
    local temp_k
    if ! temp_k=$(claim_next); then
      echo "$(date -u +%H:%M:%SZ) gpu=$gpu queue empty; exit"
      break
    fi
    echo "$(date -u +%H:%M:%SZ) gpu=$gpu starting temp_${temp_k}K"
    run_one "$gpu" "$temp_k" || echo "$(date -u +%H:%M:%SZ) gpu=$gpu temp_${temp_k}K failed rc=$?"
  done
}

build_queue

IFS=',' read -ra GPU_ARR <<<"$GPUS"
pids=()
for gpu in "${GPU_ARR[@]}"; do
  for ((slot=0; slot<PER_GPU; slot++)); do
    worker "$gpu" >"$LOGDIR/worker_gpu${gpu}_slot${slot}.log" 2>&1 &
    pids+=($!)
  done
done

for pid in "${pids[@]}"; do
  wait "$pid"
done

echo "$(date -u +%H:%M:%SZ) all workers finished; running analyzer"
"$PY" "$ANALYZE" \
  --root "$OUT" \
  --discard-ps "$DISCARD_PS" \
  --prod-end-ps "$RUNTIME_PS"

echo "$(date -u +%H:%M:%SZ) calibration complete"
