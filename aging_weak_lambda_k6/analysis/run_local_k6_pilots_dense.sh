#!/usr/bin/env bash
# Dense local k=6.02 pilot: multiple HOOMD replicas per GPU via a shared queue.
set -euo pipefail

PY=/scratch/mh7373/miniforge3/envs/hoomd/bin/python
REPO=/scratch/mh7373/projects/cav-hoomd
BASE=$REPO/aging_weak_lambda_k6
SCRIPT=$REPO/examples/05_advanced_run.py
GSD=$REPO/examples/init-0.gsd
LOGDIR=$BASE/analysis/local_pilot_logs
N_REPLICAS=${N_REPLICAS:-50}
PER_GPU=${PER_GPU:-4}
MIN_H5_BYTES=${MIN_H5_BYTES:-1380000000}
QUEUE=$LOGDIR/replica_queue.tsv
LOCKDIR=$LOGDIR/queue_locks

mkdir -p "$LOGDIR" "$LOCKDIR" "$BASE/lambda0" "$BASE/lambda0p03"

is_complete() {
  local tag="$1" replica="$2"
  local cwd run_dir h5
  if [[ "$tag" == "0" ]]; then cwd=$BASE/lambda0; else cwd=$BASE/lambda0p03; fi
  run_dir=$(ls -d "$cwd"/cavity_coupling_*_switch_200.0ps 2>/dev/null | head -1 || true)
  [[ -n "$run_dir" ]] || return 1
  h5="$run_dir/observables_replica_${replica}.h5"
  [[ -f "$h5" ]] && [[ $(stat -c%s "$h5") -ge $MIN_H5_BYTES ]]
}

# Build queue: all incomplete (lam, tag, replica) jobs.
: >"$QUEUE"
for ((r=0; r<N_REPLICAS; r++)); do
  for pair in "0.0:0" "0.03:0p03"; do
    IFS=: read -r lam tag <<<"$pair"
    if ! is_complete "$tag" "$r"; then
      # Skip if another process already holds a run lock / h5 is actively growing
      # with a live lock from the packed convention.
      cwd=$BASE/lambda0
      [[ "$tag" == "0p03" ]] && cwd=$BASE/lambda0p03
      run_dir=$(ls -d "$cwd"/cavity_coupling_*_switch_200.0ps 2>/dev/null | head -1 || true)
      if [[ -n "$run_dir" && -f "$run_dir/.replica_${r}.lock" ]]; then
        continue
      fi
      # Also skip if a local claim lock exists
      if [[ -f "$LOCKDIR/claim_${tag}_${r}.lock" ]]; then
        continue
      fi
      printf '%s\t%s\t%s\n' "$lam" "$tag" "$r" >>"$QUEUE"
    fi
  done
done

njobs=$(wc -l <"$QUEUE" | tr -d ' ')
echo "$(date -u +%H:%M:%SZ) dense queue size=$njobs PER_GPU=$PER_GPU"

claim_next() {
  local line lam tag replica claim
  (
    flock 9
    while IFS=$'\t' read -r lam tag replica; do
      claim="$LOCKDIR/claim_${tag}_${replica}.lock"
      if [[ -e "$claim" ]]; then
        continue
      fi
      if is_complete "$tag" "$replica"; then
        continue
      fi
      # Remove this line from queue
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
  if [[ "$tag" == "0" ]]; then cwd=$BASE/lambda0; else cwd=$BASE/lambda0p03; fi
  echo "$(date -u +%H:%M:%SZ) start gpu=$gpu lam=$lam replica=$replica"
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
  ) >"$LOGDIR/dense_lam${tag}_r${replica}_gpu${gpu}.out" \
    2>"$LOGDIR/dense_lam${tag}_r${replica}_gpu${gpu}.err"
  local rc=$?
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

# Mark currently running replicas (from pgrep) so we do not double-claim them.
while read -r line; do
  lam=$(echo "$line" | sed -n 's/.*--lambda-coupling \([^ ]*\).*/\1/p')
  rep=$(echo "$line" | sed -n 's/.*--replicas \([^ ]*\).*/\1/p')
  if [[ "$lam" == "0.0" || "$lam" == "0" ]]; then tag=0; else tag=0p03; fi
  if [[ -n "$rep" ]]; then
    touch "$LOCKDIR/claim_${tag}_${rep}.lock"
    echo "preserve running claim tag=$tag replica=$rep"
  fi
done < <(pgrep -af '05_advanced_run.py' | rg 'aging_weak_lambda_k6|lambda-coupling' || true)

# Rebuild queue excluding preserved claims
: >"$QUEUE"
for ((r=0; r<N_REPLICAS; r++)); do
  for pair in "0.0:0" "0.03:0p03"; do
    IFS=: read -r lam tag <<<"$pair"
    if is_complete "$tag" "$r"; then continue; fi
    if [[ -f "$LOCKDIR/claim_${tag}_${r}.lock" ]]; then continue; fi
    printf '%s\t%s\t%s\n' "$lam" "$tag" "$r" >>"$QUEUE"
  done
done
echo "$(date -u +%H:%M:%SZ) queue after preserve: $(wc -l <"$QUEUE") jobs"

pids=()
for gpu in 0 1; do
  for ((slot=0; slot<PER_GPU; slot++)); do
    # Leave one slot free on each GPU if a legacy sequential worker still owns a sim
    if [[ $slot -eq 0 ]] && pgrep -af "CUDA_VISIBLE_DEVICES=$gpu" >/dev/null 2>&1; then
      :
    fi
    worker "$gpu" >"$LOGDIR/dense_worker_gpu${gpu}_slot${slot}.log" 2>&1 &
    pids+=($!)
  done
done

printf '%s\n' "${pids[@]}" >"$LOGDIR/dense_worker.pids"
echo "started ${#pids[@]} dense workers: ${pids[*]}"
wait "${pids[@]}" || true
echo "$(date -u +%H:%M:%SZ) dense launcher finished"
