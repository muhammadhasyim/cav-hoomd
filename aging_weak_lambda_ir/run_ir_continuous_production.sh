#!/usr/bin/env bash
# Launch continuous 400 ps IR production on local GPUs + incremental DACF daemon.
# Usage:
#   nohup bash aging_weak_lambda_ir/run_ir_continuous_production.sh \
#     > aging_weak_lambda_ir/logs/ir_continuous/nohup_production.log 2>&1 &
set -euo pipefail

REPO="${REPO:-/scratch/mh7373/projects/cav-hoomd}"
LOGDIR="$REPO/aging_weak_lambda_ir/logs/ir_continuous"
mkdir -p "$LOGDIR"

source /scratch/mh7373/miniforge3/etc/profile.d/conda.sh
conda activate hoomd
cd "$REPO"

GPUS="${GPUS:-0,1}"
RUNTIME_PS="${RUNTIME_PS:-400.0}"
N_REPLICAS="${N_REPLICAS:-2}"
POLL_INTERVAL_S="${POLL_INTERVAL_S:-120.0}"
IR_METHOD="${IR_METHOD:-dct}"

echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) IR continuous production start"
echo "GPUS=$GPUS RUNTIME_PS=$RUNTIME_PS N_REPLICAS=$N_REPLICAS POLL_INTERVAL_S=$POLL_INTERVAL_S"

python aging_weak_lambda_ir/ir_continuous_daemon.py \
  --poll-interval-s "$POLL_INTERVAL_S" \
  --method "$IR_METHOD" \
  >"$LOGDIR/daemon.stdout.log" 2>"$LOGDIR/daemon.stderr.log" &
daemon_pid=$!
echo "$daemon_pid" >"$LOGDIR/daemon.pid"
echo "started DACF daemon pid=$daemon_pid"

GPUS="$GPUS" RUNTIME_PS="$RUNTIME_PS" N_REPLICAS="$N_REPLICAS" \
  bash aging_weak_lambda_ir/run_local_ir_continuous.sh

echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) IR continuous workers exited; stopping daemon"
kill "$daemon_pid" 2>/dev/null || true
