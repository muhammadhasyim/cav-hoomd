#!/bin/bash
# ============================================================
# Weak-coupling aging campaign — cav-hoomd, ADAPTIVE timestep
# Mirrors submit_aging_campaign.sh but uses adaptive dt
# (default 05_advanced_run.py parameters: error_tol=5.0,
#  initial_fraction=1e-5, time_constant=50 ps)
#
# With these params dt≈0.025 fs steady-state → ~36 h/replica;
# time limit set to 47:59:00 so all 2500 ps complete.
#
# Usage:
#   ./submit_aging_adaptive.sh [--dry-run] [--test] [--lam <v>]
#   ./submit_aging_adaptive.sh --start <n> --end <m>   # partial range
# ============================================================

set -euo pipefail

# ── Physical parameters (identical to fixed-dt campaign) ──────────────
LAMBDAS=(0.0 0.01 0.016667 0.023333 0.03)
N_REPLICAS=1000
TEMPERATURE=100.0
FREQUENCY=1560.0
RUNTIME_PS=2500.0
SWITCH_TIME_PS=200.0
FKT_KMAG=1.0085
FKT_REF_INTERVAL=200.0
FKT_MAX_REFS=13
FKT_OUTPUT_PERIOD=1.0
ENERGY_OUTPUT_PERIOD=1.0

# ── Adaptive timestep parameters (05_advanced_run.py defaults) ─────────
ERROR_TOLERANCE=5.0      # target error tolerance (a.u.²)
INITIAL_FRACTION=1e-5    # initial fraction of error_tol (shock dampening start)
TIME_CONSTANT_PS=50.0    # exponential recovery time constant (ps)

# ── Paths ──────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CAV_HOOMD_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
EXAMPLES_DIR="$CAV_HOOMD_DIR/examples"
INPUT_GSD="$CAV_HOOMD_DIR/examples/init-0.gsd"
OUTPUT_BASE="/scratch/mh7373/projects/cav-hoomd/aging_weak_lambda_adaptive"
CONDA_ENV="hoomd"
MINIFORGE="/scratch/mh7373/miniforge3"

# ── SLURM settings ─────────────────────────────────────────────────────
PARTITION="a100_chemistry"
ACCOUNT="torch_pr_283_chemistry"
GRES="gpu:a100:1"
CPUS=4
MEM="16G"
TIME_LIMIT="47:59:00"   # just under 48 h (gpu48 QOS) — enough for ~36 h adaptive run
CONCURRENT=4

# ── Argument parsing ───────────────────────────────────────────────────
DRY_RUN=false
TEST_MODE=false
ONLY_LAM=""
RANGE_START=0
RANGE_END=$((N_REPLICAS - 1))

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)  DRY_RUN=true;  shift ;;
        --test)     TEST_MODE=true; shift ;;
        --lam)      ONLY_LAM="$2"; shift 2 ;;
        --start)    RANGE_START="$2"; shift 2 ;;
        --end)      RANGE_END="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ "$TEST_MODE" == "true" ]]; then
    RANGE_START=0; RANGE_END=3
    echo "*** TEST MODE — replicas 0-3 per λ ***"
fi

if [[ "$DRY_RUN" == "true" ]]; then
    echo "*** DRY RUN — no jobs submitted ***"
fi

# ── Setup ──────────────────────────────────────────────────────────────
mkdir -p "$OUTPUT_BASE"
mkdir -p "$SCRIPT_DIR/logs"
mkdir -p "$SCRIPT_DIR/generated"

SUBMITTED_LOG="$SCRIPT_DIR/submitted_aging_adaptive.txt"

echo "============================================================"
echo "Weak-coupling aging campaign — ADAPTIVE timestep"
echo "============================================================"
echo "Date:          $(date)"
echo "λ values:      ${LAMBDAS[*]}"
echo "Replicas:      ${RANGE_START}–${RANGE_END}"
echo "Runtime:       $RUNTIME_PS ps  (switch at $SWITCH_TIME_PS ps)"
echo "T, freq:       $TEMPERATURE K, $FREQUENCY cm⁻¹"
echo "Adaptive:      tol=${ERROR_TOLERANCE}, frac=${INITIAL_FRACTION}, τ=${TIME_CONSTANT_PS} ps"
echo "Output:        $OUTPUT_BASE"
echo "Time limit:    $TIME_LIMIT"
echo "============================================================"

TOTAL_JOBS=0
ARRAY_END=$RANGE_END

for LAM in "${LAMBDAS[@]}"; do
    if [[ -n "$ONLY_LAM" && "$ONLY_LAM" != "$LAM" ]]; then
        continue
    fi

    if [[ "$LAM" == "0.0" || "$LAM" == "0" ]]; then
        LAM_TAG="0"
    else
        LAM_TAG="${LAM//./p}"
    fi
    JOB_NAME="cav-adp-lam${LAM_TAG}"
    LAM_DIR="$OUTPUT_BASE/lambda${LAM_TAG}"
    mkdir -p "$LAM_DIR"

    SBATCH_FILE="$SCRIPT_DIR/generated/adaptive_lam${LAM_TAG}_r${RANGE_START}-${RANGE_END}.sbatch"

    cat > "$SBATCH_FILE" << EOF
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --partition=${PARTITION}
#SBATCH --account=${ACCOUNT}
#SBATCH --gres=${GRES}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --time=${TIME_LIMIT}
#SBATCH --array=${RANGE_START}-${ARRAY_END}%${CONCURRENT}
#SBATCH --output=${SCRIPT_DIR}/logs/${JOB_NAME}_%A_%a.out
#SBATCH --error=${SCRIPT_DIR}/logs/${JOB_NAME}_%A_%a.err
#SBATCH --mail-type=NONE

source ${MINIFORGE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

REPLICA=\${SLURM_ARRAY_TASK_ID}
cd "${LAM_DIR}"

echo "===================================================================="
echo "cav-hoomd aging campaign — ADAPTIVE timestep"
echo "  λ          = ${LAM}"
echo "  Replica    = \${REPLICA}"
echo "  Node       = \$(hostname)"
echo "  SLURM job  = \${SLURM_JOB_ID}_\${SLURM_ARRAY_TASK_ID}"
echo "  Date       = \$(date)"
echo "  error_tol  = ${ERROR_TOLERANCE}  initial_frac = ${INITIAL_FRACTION}  τ = ${TIME_CONSTANT_PS} ps"
echo "===================================================================="
nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader
echo "===================================================================="

COUPLING_FLAG=""
if [[ "${LAM}" == "0.0" ]]; then
    COUPLING_FLAG="--coupling 0.0"
else
    COUPLING_FLAG="--coupling ${LAM}"
fi

python "${EXAMPLES_DIR}/05_advanced_run.py" \\
    --molecular-bath bussi \\
    --cavity-bath langevin \\
    \${COUPLING_FLAG} \\
    --temperature ${TEMPERATURE} \\
    --frequency ${FREQUENCY} \\
    --runtime ${RUNTIME_PS} \\
    --switch-time ${SWITCH_TIME_PS} \\
    --input-gsd "${INPUT_GSD}" \\
    --frame -1 \\
    --device GPU \\
    --gpu-id 0 \\
    --error-tolerance ${ERROR_TOLERANCE} \\
    --initial-fraction ${INITIAL_FRACTION} \\
    --time-constant-ps ${TIME_CONSTANT_PS} \\
    --enable-energy-tracker \\
    --energy-output-period-ps ${ENERGY_OUTPUT_PERIOD} \\
    --enable-fkt \\
    --fkt-kmag ${FKT_KMAG} \\
    --fkt-wavevectors 50 \\
    --fkt-ref-interval ${FKT_REF_INTERVAL} \\
    --fkt-max-refs ${FKT_MAX_REFS} \\
    --fkt-output-period-ps ${FKT_OUTPUT_PERIOD} \\
    --gsd-output-period-ps 999999 \\
    --console-output-period-ps 100.0

RC=\$?
if [[ \$RC -eq 0 ]]; then
    echo "SUCCESS: λ=${LAM} replica \${REPLICA} completed"
else
    echo "FAILED:  λ=${LAM} replica \${REPLICA} exit code \${RC}" >&2
    exit \$RC
fi
EOF

    chmod +x "$SBATCH_FILE"
    echo ""
    echo "λ=${LAM}  →  ${JOB_NAME}  array ${RANGE_START}-${ARRAY_END}%${CONCURRENT}"

    if [[ "$DRY_RUN" == "false" ]]; then
        JOB_ID=$(sbatch --parsable "$SBATCH_FILE")
        echo "  Submitted: job $JOB_ID"
        echo "${JOB_ID}:lam=${LAM}:replicas ${RANGE_START}-${RANGE_END}" >> "$SUBMITTED_LOG"
        TOTAL_JOBS=$((TOTAL_JOBS + 1))
        sleep 0.5
    else
        echo "  [DRY RUN] Would submit: $SBATCH_FILE"
        TOTAL_JOBS=$((TOTAL_JOBS + 1))
    fi
done

echo ""
echo "============================================================"
echo "Submitted: $TOTAL_JOBS job arrays"
echo "Logged to: $SUBMITTED_LOG"
echo ""
echo "Monitor:   squeue -u \$USER | grep cav-adp"
echo "Output:    $OUTPUT_BASE"
echo ""
echo "To submit replicas 4-999 once queue slots open:"
echo "  ./submit_aging_adaptive.sh --start 4 --end 999"
echo "============================================================"
