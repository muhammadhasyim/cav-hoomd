#!/usr/bin/env bash
# Full equilibrium IR smoke: 10 GPU jobs + IR spectrum post-processing.
set -euo pipefail

REPO="${REPO:-/scratch/mh7373/projects/cav-hoomd}"
LOGDIR="$REPO/aging_weak_lambda_ir/logs/ir_smoke"
mkdir -p "$LOGDIR"

source /scratch/mh7373/miniforge3/etc/profile.d/conda.sh
conda activate hoomd
cd "$REPO"

echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) IR smoke production start"
bash aging_weak_lambda_ir/run_local_ir_smoke.sh

echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) IR analysis start"
python aging_weak_lambda_ir/analyze_ir_from_dipole.py

echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) IR smoke production complete"
