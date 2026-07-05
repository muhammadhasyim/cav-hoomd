#!/bin/bash
# ============================================================
# Weak-coupling aging campaign — cav-hoomd
# Matches the OpenMM aging_weak_lambda N=1000 campaign
# (same λ grid, T, frequency, runtime, switch time)
# ============================================================
#
# Usage:
#   ./submit_aging_campaign.sh [--dry-run] [--test] [--lam <val>]
#
#   --dry-run   Print sbatch commands without submitting
#   --test      Submit 4 replicas per λ instead of 1000 (quick smoke test)
#   --lam <v>   Submit only for a single λ value (e.g. --lam 0.01)
#
# ============================================================

set -euo pipefail

# ── Physical parameters (matching OpenMM config.py) ──────────────────
LAMBDAS=(0.0 0.01 0.016667 0.023333 0.03)
N_REPLICAS=1000          # replicas per λ (SLURM array 0–999)
TEMPERATURE=100.0        # K
FREQUENCY=1560.0         # cm⁻¹ (resonant with A–A stretch)
RUNTIME_PS=2500.0        # total simulation time
SWITCH_TIME_PS=200.0     # cavity coupling turns on at this time
FKT_KMAG=1.0085          # k-magnitude in Bohr⁻¹ (matches OpenMM FKT_KMAG_AU)
FKT_REF_INTERVAL=200.0   # ps between F(k,t) reference updates
FKT_MAX_REFS=13          # maximum number of F(k,t) references kept
FKT_OUTPUT_PERIOD=1.0    # ps between F(k,t) data points
ENERGY_OUTPUT_PERIOD=1.0 # ps between energy data points

# ── Paths ─────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CAV_HOOMD_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
EXAMPLES_DIR="$CAV_HOOMD_DIR/examples"
INPUT_GSD="$CAV_HOOMD_DIR/cooling-0.gsd"    # 250 dimers, T=100K, box=40 Bohr
OUTPUT_BASE="/scratch/mh7373/projects/cav-hoomd/aging_weak_lambda"
CONDA_ENV="hoomd"                            # micromamba env with HOOMD 5.4.0 + plugin
MINIFORGE="/scratch/mh7373/miniforge3"

# ── SLURM settings ────────────────────────────────────────────────────
PARTITION="a100_chem"
GRES="gpu:a100:1"
CPUS=4
MEM="16G"
TIME_LIMIT="24:00:00"
CONCURRENT=4             # simultaneous tasks per array (%N throttle)

# ── Argument parsing ──────────────────────────────────────────────────
DRY_RUN=false
TEST_MODE=false
ONLY_LAM=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)  DRY_RUN=true;  shift ;;
        --test)     TEST_MODE=true; shift ;;
        --lam)      ONLY_LAM="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ "$TEST_MODE" == "true" ]]; then
    N_REPLICAS=4
    echo "*** TEST MODE — 4 replicas per λ ***"
fi

if [[ "$DRY_RUN" == "true" ]]; then
    echo "*** DRY RUN — no jobs will be submitted ***"
fi

# ── Setup ─────────────────────────────────────────────────────────────
mkdir -p "$OUTPUT_BASE"
mkdir -p "$SCRIPT_DIR/logs"
mkdir -p "$SCRIPT_DIR/generated"

SUBMITTED_LOG="$SCRIPT_DIR/submitted_aging_campaign.txt"
> "$SUBMITTED_LOG"

echo "============================================================"
echo "Weak-coupling aging campaign — cav-hoomd"
echo "============================================================"
echo "Date: $(date)"
echo "λ values:      ${LAMBDAS[*]}"
echo "Replicas/λ:    $N_REPLICAS"
echo "Runtime:       $RUNTIME_PS ps  (switch at $SWITCH_TIME_PS ps)"
echo "T, freq:       $TEMPERATURE K, $FREQUENCY cm⁻¹"
echo "Output:        $OUTPUT_BASE"
echo "Input GSD:     $INPUT_GSD"
echo "HOOMD env:     $CONDA_ENV"
echo "Partition:     $PARTITION"
echo "Concurrent:    $CONCURRENT per array"
echo "============================================================"

ARRAY_END=$((N_REPLICAS - 1))
TOTAL_JOBS=0

for LAM in "${LAMBDAS[@]}"; do
    # Optional single-λ filter
    if [[ -n "$ONLY_LAM" && "$ONLY_LAM" != "$LAM" ]]; then
        continue
    fi

    # Human-readable λ tag for filenames
    if [[ "$LAM" == "0.0" || "$LAM" == "0" ]]; then
        LAM_TAG="0"
    else
        LAM_TAG="${LAM//./p}"
    fi
    JOB_NAME="cav-aging-lam${LAM_TAG}"
    LAM_DIR="$OUTPUT_BASE/lambda${LAM_TAG}"
    mkdir -p "$LAM_DIR"

    SBATCH_FILE="$SCRIPT_DIR/generated/aging_lam${LAM_TAG}.sbatch"

    cat > "$SBATCH_FILE" << EOF
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --partition=${PARTITION}
#SBATCH --gres=${GRES}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --time=${TIME_LIMIT}
#SBATCH --array=0-${ARRAY_END}%${CONCURRENT}
#SBATCH --output=${SCRIPT_DIR}/logs/${JOB_NAME}_%A_%a.out
#SBATCH --error=${SCRIPT_DIR}/logs/${JOB_NAME}_%A_%a.err
#SBATCH --mail-type=NONE

# ── Environment ────────────────────────────────────────────────
source ${MINIFORGE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

REPLICA=\${SLURM_ARRAY_TASK_ID}
LAM=${LAM}
LAM_DIR="${LAM_DIR}"
EXAMPLES_DIR="${EXAMPLES_DIR}"
INPUT_GSD="${INPUT_GSD}"

echo "================================================================"
echo "cav-hoomd aging campaign"
echo "  λ          = \$LAM"
echo "  Replica    = \$REPLICA"
echo "  Node       = \$(hostname)"
echo "  SLURM job  = \${SLURM_JOB_ID}_\${SLURM_ARRAY_TASK_ID}"
echo "  Date       = \$(date)"
echo "================================================================"
nvidia-smi --query-gpu=index,name,memory.total,driver_version --format=csv,noheader
echo "================================================================"

# ── Run in the λ output directory ──────────────────────────────
cd "\${LAM_DIR}"

COUPLING_FLAG=""
if [[ "\${LAM}" == "0.0" ]]; then
    # λ=0: cavity present but uncoupled (keeps output structure consistent)
    COUPLING_FLAG="--coupling 0.0"
else
    COUPLING_FLAG="--coupling \${LAM}"
fi

python "\${EXAMPLES_DIR}/05_advanced_run.py" \\
    --molecular-bath bussi \\
    --cavity-bath langevin \\
    \${COUPLING_FLAG} \\
    --temperature ${TEMPERATURE} \\
    --frequency ${FREQUENCY} \\
    --runtime ${RUNTIME_PS} \\
    --switch-time ${SWITCH_TIME_PS} \\
    --input-gsd "\${INPUT_GSD}" \\
    --frame -1 \\
    --device GPU \\
    --gpu-id 0 \\
    --enable-energy-tracker \\
    --energy-output-period-ps ${ENERGY_OUTPUT_PERIOD} \\
    --enable-fkt \\
    --fkt-kmag ${FKT_KMAG} \\
    --fkt-wavevectors 50 \\
    --fkt-ref-interval ${FKT_REF_INTERVAL} \\
    --fkt-max-refs ${FKT_MAX_REFS} \\
    --fkt-output-period-ps ${FKT_OUTPUT_PERIOD} \\
    --gsd-output-period-ps 999999 \\
    --console-output-period-ps 75.0 \\
    --error-tolerance 5.0 \\
    --initial-fraction 1e-5 \\
    --time-constant-ps 50.0

RC=\$?
if [[ \$RC -eq 0 ]]; then
    echo "SUCCESS: λ=\${LAM} replica \${REPLICA} completed"
else
    echo "FAILED:  λ=\${LAM} replica \${REPLICA} exit code \${RC}" >&2
    exit \$RC
fi
EOF

    chmod +x "$SBATCH_FILE"
    echo ""
    echo "λ=${LAM}  →  ${JOB_NAME}  array 0-${ARRAY_END}%${CONCURRENT}"

    if [[ "$DRY_RUN" == "false" ]]; then
        JOB_ID=$(sbatch --parsable "$SBATCH_FILE")
        echo "  Submitted: job $JOB_ID"
        echo "${JOB_ID}:lam=${LAM}:${N_REPLICAS} replicas" >> "$SUBMITTED_LOG"
        TOTAL_JOBS=$((TOTAL_JOBS + 1))
        sleep 0.5
    else
        echo "  [DRY RUN] Would submit: $SBATCH_FILE"
        TOTAL_JOBS=$((TOTAL_JOBS + 1))
    fi
done

echo ""
echo "============================================================"
echo "Total λ values submitted: $TOTAL_JOBS"
echo "Total simulations:        $((TOTAL_JOBS * N_REPLICAS))"
echo "Job IDs logged to:        $SUBMITTED_LOG"
echo ""
echo "Monitor with:"
echo "  squeue -u \$USER | grep cav-aging"
echo "  tail -f $SCRIPT_DIR/logs/cav-aging-lam0_*.out"
echo "============================================================"
