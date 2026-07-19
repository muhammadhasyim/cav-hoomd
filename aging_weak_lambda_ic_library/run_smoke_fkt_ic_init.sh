#!/usr/bin/env bash
# 50-replica lambda=0 smoke test using init-500-from-ic-library.gsd (--frame replica).
# F(k,t) at |k|=1.00 a.u., zero coupling, production-like bath settings.
set -euo pipefail

PY=/scratch/mh7373/miniforge3/envs/hoomd/bin/python
REPO=/scratch/mh7373/projects/cav-hoomd
BASE=$REPO/aging_weak_lambda_ic_library/smoke_fkt_ic_init
SCRIPT=$REPO/examples/05_advanced_run.py
INIT_GSD=$REPO/aging_weak_lambda_ic_library/init-500-from-ic-library.gsd
LOGDIR=$BASE/logs
RUNDIR=$BASE/lambda0
COUPLING_DIR=coupling_000000e+00_switch_200.0ps
N_REPLICAS=${N_REPLICAS:-50}
PER_GPU=${PER_GPU:-2}
RUNTIME_PS=${RUNTIME_PS:-300.0}
FKT_KMAG=${FKT_KMAG:-1.0}

mkdir -p "$LOGDIR" "$RUNDIR"

run_one() {
  local gpu="$1" replica="$2"
  local outdir="$RUNDIR/$COUPLING_DIR"
  mkdir -p "$outdir"
  local fkt_ref="$outdir/prod-${replica}_fkt_ref_000.txt"
  if [[ -f "$fkt_ref" ]] && [[ $(wc -l <"$fkt_ref") -gt 100 ]]; then
    echo "$(date -u +%H:%M:%SZ) skip replica=$replica (fkt present)"
    return 0
  fi
  echo "$(date -u +%H:%M:%SZ) start gpu=$gpu replica=$replica frame=$replica k=$FKT_KMAG"
  (
    cd "$RUNDIR"
    CUDA_VISIBLE_DEVICES="$gpu" "$PY" "$SCRIPT" \
      --molecular-bath bussi \
      --cavity-bath langevin \
      --lambda-coupling 0.0 \
      --coupling-type step \
      --temperature 100.0 \
      --frequency 1560.0 \
      --runtime "$RUNTIME_PS" \
      --switch-time 200.0 \
      --input-gsd "$INIT_GSD" \
      --frame "$replica" \
      --device GPU \
      --gpu-id 0 \
      --fixed-timestep \
      --timestep 1.0 \
      --enable-fkt \
      --fkt-kmag "$FKT_KMAG" \
      --fkt-wavevectors 50 \
      --fkt-ref-interval 200.0 \
      --fkt-max-refs 5 \
      --fkt-output-period-ps 1.0 \
      --disable-gsd \
      --disable-hdf5-output \
      --disable-temp-tracker \
      --console-output-period-ps 200.0 \
      --replicas "$replica" \
      --seed $((replica + 1))
  ) >"$LOGDIR/replica_${replica}_gpu${gpu}.out" \
    2>"$LOGDIR/replica_${replica}_gpu${gpu}.err"
  local rc=$?
  echo "$(date -u +%H:%M:%SZ) done gpu=$gpu replica=$replica rc=$rc"
  return $rc
}

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
    if (( slot_gpu0 < PER_GPU )); then
      gpu=0
      slot_gpu0=$((slot_gpu0 + 1))
    else
      gpu=1
      slot_gpu1=$((slot_gpu1 + 1))
    fi
  fi
  run_one "$gpu" "$r" &
  pids+=($!)
done

printf '%s\n' "${pids[@]}" >"$LOGDIR/smoke.pids"
echo "started ${#pids[@]} smoke jobs: ${pids[*]}"
wait "${pids[@]}" || true
echo "$(date -u +%H:%M:%SZ) smoke launcher finished"
