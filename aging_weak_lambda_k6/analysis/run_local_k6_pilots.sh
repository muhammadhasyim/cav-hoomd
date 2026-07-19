#!/usr/bin/env bash
# Run k=6.02 pilot replicas locally on available GPUs (gl061).
# Two workers (GPU 0 and GPU 1), 50 replicas each for lambda=0 and 0.03.
set -euo pipefail

PY=/scratch/mh7373/miniforge3/envs/hoomd/bin/python
REPO=/scratch/mh7373/projects/cav-hoomd
BASE=$REPO/aging_weak_lambda_k6
SCRIPT=$REPO/examples/05_advanced_run.py
GSD=$REPO/examples/init-0.gsd
LOGDIR=$BASE/analysis/local_pilot_logs
N_REPLICAS=${N_REPLICAS:-50}
mkdir -p "$LOGDIR" \
  "$BASE/lambda0" \
  "$BASE/lambda0p03"

run_one() {
  local gpu="$1" lam="$2" replica="$3" tag="$4"
  local outdir cwd
  if [[ "$tag" == "0" ]]; then
    cwd=$BASE/lambda0
  else
    cwd=$BASE/lambda0p03
  fi
  local run_dir
  run_dir=$(ls -d "$cwd"/cavity_coupling_*_switch_200.0ps 2>/dev/null | head -1 || true)
  if [[ -n "${run_dir}" ]]; then
    local h5="$run_dir/observables_replica_${replica}.h5"
    if [[ -f "$h5" ]] && [[ $(stat -c%s "$h5") -ge 1380000000 ]]; then
      echo "$(date -u +%H:%M:%S) skip complete lam=$lam replica=$replica"
      return 0
    fi
    rm -f "$run_dir/.replica_${replica}.lock"
  fi
  echo "$(date -u +%H:%M:%S) start gpu=$gpu lam=$lam replica=$replica"
  (
    cd "$cwd"
    CUDA_VISIBLE_DEVICES="$gpu" "$PY" "$SCRIPT" \
      --molecular-bath bussi --cavity-bath langevin \
      --lambda-coupling "$lam" --coupling-type step \
      --temperature 100.0 --frequency 1560.0 \
      --runtime 2500.0 --switch-time 200.0 \
      --input-gsd "$GSD" --frame -1 \
      --device GPU --gpu-id 0 \
      --fixed-timestep --timestep 1.0 \
      --enable-energy-tracker --energy-output-period-ps 1.0 \
      --enable-fkt --fkt-kmag 6.02 --fkt-wavevectors 50 \
      --fkt-ref-interval 200.0 --fkt-max-refs 13 \
      --fkt-output-period-ps 1.0 --disable-gsd \
      --console-output-period-ps 100.0 \
      --replicas "$replica" --seed $((replica + 1))
  ) >"$LOGDIR/lam${tag}_r${replica}_gpu${gpu}.out" 2>"$LOGDIR/lam${tag}_r${replica}_gpu${gpu}.err"
  local rc=$?
  echo "$(date -u +%H:%M:%S) done gpu=$gpu lam=$lam replica=$replica rc=$rc"
  return $rc
}

worker() {
  local gpu="$1"
  shift
  local jobs=("$@")
  for job in "${jobs[@]}"; do
    IFS=: read -r lam tag replica <<<"$job"
    run_one "$gpu" "$lam" "$replica" "$tag" || true
  done
}

# Interleave lambda=0 and 0.03 across replicas 0..N-1
jobs_gpu0=()
jobs_gpu1=()
for ((r=0; r<N_REPLICAS; r++)); do
  if (( r % 2 == 0 )); then
    jobs_gpu0+=("0.0:0:$r")
    jobs_gpu1+=("0.03:0p03:$r")
  else
    jobs_gpu0+=("0.03:0p03:$r")
    jobs_gpu1+=("0.0:0:$r")
  fi
done

echo "Launching local k6 pilots: N_REPLICAS=$N_REPLICAS on GPUs 0 and 1"
worker 0 "${jobs_gpu0[@]}" >"$LOGDIR/worker_gpu0.log" 2>&1 &
pid0=$!
worker 1 "${jobs_gpu1[@]}" >"$LOGDIR/worker_gpu1.log" 2>&1 &
pid1=$!
echo "workers pid0=$pid0 pid1=$pid1"
echo "$pid0" >"$LOGDIR/worker_gpu0.pid"
echo "$pid1" >"$LOGDIR/worker_gpu1.pid"
wait "$pid0" "$pid1"
echo "All local pilot workers finished"
