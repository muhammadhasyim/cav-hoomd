#!/usr/bin/env bash
# Poll the k=6.02 pilot until enough complete replicas exist, then analyze.
set -euo pipefail

PY=/scratch/mh7373/miniforge3/envs/hoomd/bin/python
BASE=/scratch/mh7373/projects/cav-hoomd/aging_weak_lambda_k6
ANALYSIS=$BASE/analysis
MIN_REPLICAS=${MIN_REPLICAS:-20}
POLL_SEC=${POLL_SEC:-300}
JOB_ID=${JOB_ID:-14156210}

mkdir -p "$ANALYSIS"
LOG=$ANALYSIS/watch_pilot.log

count_complete() {
  local dir="$1"
  "$PY" - "$dir" <<'PY'
from pathlib import Path
import sys
run = Path(sys.argv[1])
min_bytes = 1_380_000_000
n = sum(1 for p in run.glob("observables_replica_*.h5") if p.is_file() and p.stat().st_size >= min_bytes)
print(n)
PY
}

echo "==== watch start $(date -u +%Y%m%dT%H%M%SZ) job=$JOB_ID min=$MIN_REPLICAS ====" | tee -a "$LOG"

while true; do
  d0=$(ls -d "$BASE"/lambda0/cavity_coupling_*_switch_200.0ps 2>/dev/null | head -1 || true)
  d3=$(ls -d "$BASE"/lambda0p03/cavity_coupling_*_switch_200.0ps 2>/dev/null | head -1 || true)
  n0=0
  n3=0
  if [[ -n "${d0}" ]]; then n0=$(count_complete "$d0"); fi
  if [[ -n "${d3}" ]]; then n3=$(count_complete "$d3"); fi
  echo "$(date -u +%Y%m%dT%H%M%SZ) complete_h5 lambda0=$n0 lambda0p03=$n3" | tee -a "$LOG"
  squeue -j "$JOB_ID" -h -o '%T %r' 2>/dev/null | head -5 | tee -a "$LOG" || true

  if [[ "$n0" -ge "$MIN_REPLICAS" && "$n3" -ge "$MIN_REPLICAS" ]]; then
    echo "Threshold reached; running analysis" | tee -a "$LOG"
    "$PY" "$ANALYSIS/run_k6_pilot_analysis.py" --min-replicas "$MIN_REPLICAS" 2>&1 | tee -a "$LOG"
    "$PY" /scratch/mh7373/projects/cav-hoomd/aging_weak_lambda/analysis/ic_equilibration_audit.py \
      --fkt-kmag 6.02 \
      --output "$ANALYSIS/ic_equilibration_audit_report.txt" 2>&1 | tee -a "$LOG"
    echo "==== watch done $(date -u +%Y%m%dT%H%M%SZ) ====" | tee -a "$LOG"
    exit 0
  fi

  # Exit if the array vanished and we still lack data.
  if ! squeue -j "$JOB_ID" -h >/dev/null 2>&1; then
    if [[ "$n0" -lt "$MIN_REPLICAS" || "$n3" -lt "$MIN_REPLICAS" ]]; then
      echo "Job $JOB_ID no longer queued and threshold not met" | tee -a "$LOG"
      exit 4
    fi
  fi
  sleep "$POLL_SEC"
done
