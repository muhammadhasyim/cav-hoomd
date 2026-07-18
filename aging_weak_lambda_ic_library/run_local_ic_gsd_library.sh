#!/usr/bin/env bash
# No-cavity IC library at 100 K: 8 replicas, 100 ns, GSD every 200 ps (500 frames each).
# GSD only — no Fkt / energy tracker / dipole.
# Restarts from per-replica checkpoints in ic_checkpoints/ (saved from truncated prod-*.gsd).
set -euo pipefail

PY=/scratch/mh7373/miniforge3/envs/hoomd/bin/python
REPO=/scratch/mh7373/projects/cav-hoomd
BASE=$REPO/aging_weak_lambda_ic_library
SCRIPT=$REPO/examples/05_advanced_run.py
ICDIR=$BASE/ic_checkpoints
LOGDIR=$BASE/logs
N_REPLICAS=${N_REPLICAS:-8}
PER_GPU=${PER_GPU:-4}
RUNTIME_PS=${RUNTIME_PS:-100000.0}   # 500 frames * 200 ps
GSD_PERIOD_PS=${GSD_PERIOD_PS:-200.0}
TEMP_K=${TEMP_K:-100.0}

mkdir -p "$LOGDIR" "$BASE/replicas" "$ICDIR"

run_one() {
  local gpu="$1" replica="$2"
  local outdir="$BASE/replicas/replica_${replica}"
  local ic_gsd="$ICDIR/replica_${replica}.gsd"
  mkdir -p "$outdir"
  if [[ ! -f "$ic_gsd" ]]; then
    echo "missing checkpoint: $ic_gsd" >&2
    return 1
  fi
  echo "$(date -u +%H:%M:%SZ) start gpu=$gpu replica=$replica ic=$ic_gsd T=${TEMP_K}K runtime=${RUNTIME_PS}ps gsd_period=${GSD_PERIOD_PS}ps append=1"
  (
    cd "$outdir"
    CUDA_VISIBLE_DEVICES="$gpu" "$PY" "$SCRIPT" \
      --no-cavity \
      --molecular-bath bussi \
      --temperature "$TEMP_K" \
      --runtime "$RUNTIME_PS" \
      --input-gsd "$ic_gsd" \
      --frame -1 \
      --device GPU \
      --gpu-id 0 \
      --fixed-timestep \
      --timestep 1.0 \
      --gsd-output-period-ps "$GSD_PERIOD_PS" \
      --console-output-period-ps 1000.0 \
      --replicas "$replica" \
      --seed $((1000 + replica))
  ) >"$LOGDIR/replica_${replica}_gpu${gpu}.out" \
    2>"$LOGDIR/replica_${replica}_gpu${gpu}.err"
  local rc=$?
  echo "$(date -u +%H:%M:%SZ) done gpu=$gpu replica=$replica rc=$rc"
  return $rc
}

# Round-robin replicas across GPUs with PER_GPU slots each.
pids=()
slot_gpu0=0
slot_gpu1=0
for ((r=0; r<N_REPLICAS; r++)); do
  if (( slot_gpu0 < PER_GPU )); then
    gpu=0
    slot_gpu0=$((slot_gpu0 + 1))
  elif (( slot_gpu1 < PER_GPU )); then
    gpu=1
    slot_gpu1=$((slot_gpu1 + 1))
  else
    wait -n || true
    gpu=0
  fi
  run_one "$gpu" "$r" &
  pids+=($!)
done

printf '%s\n' "${pids[@]}" >"$LOGDIR/ic_library.pids"
echo "started ${#pids[@]} IC-library jobs: ${pids[*]}"
echo "checkpoints: $ICDIR/replica_*.gsd"
echo "outputs under $BASE/replicas/replica_*/no_cavity/prod-*.gsd (append mode, no truncate)"
wait "${pids[@]}" || true
echo "$(date -u +%H:%M:%SZ) IC library launcher finished"
