#!/bin/bash
# ============================================================
# Weak-coupling aging campaign — cav-hoomd
# Matches the OpenMM aging_weak_lambda N=1000 campaign
# (same λ grid, T, frequency, runtime, switch time)
# ============================================================
#
# Usage:
#   ./submit_aging_campaign.sh [--dry-run] [--test] [--lam <val>]
#   ./submit_aging_campaign.sh --resume [--cleanup] [--dry-run] [--lam <val>]
#
#   --dry-run   Print sbatch commands without submitting
#   --test      Submit 4 replicas per λ instead of 1000 (quick smoke test)
#   --lam <v>   Submit only for a single λ value (e.g. --lam 0.01)
#   --resume    Submit only missing replicas (uses completed HDF5 size check)
#   --cleanup   Delete partial/failed replica outputs before resume submit
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
COMPLETE_MIN_BYTES=1380000000

# ── Paths ─────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CAV_HOOMD_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
EXAMPLES_DIR="$CAV_HOOMD_DIR/examples"
INPUT_GSD="$CAV_HOOMD_DIR/examples/init-0.gsd"  # 250 dimers + cavity particle, T=100K, box=40 Bohr
OUTPUT_BASE="/scratch/mh7373/projects/cav-hoomd/aging_weak_lambda"
STATUS_SCRIPT="$SCRIPT_DIR/aging_campaign_status.py"
CONDA_ENV="hoomd"                            # micromamba env with HOOMD 5.4.0 + plugin
MINIFORGE="/scratch/mh7373/miniforge3"

# ── SLURM settings ────────────────────────────────────────────────────
PARTITION="a100_chemistry"
ACCOUNT="torch_pr_283_chemistry"
GRES="gpu:a100:1"
CPUS=4
MEM="16G"
TIME_LIMIT="24:00:00"
CONCURRENT=4             # simultaneous tasks per array (%N throttle)

# ── Argument parsing ──────────────────────────────────────────────────
DRY_RUN=false
TEST_MODE=false
RESUME_MODE=false
CLEANUP=false
ONLY_LAM=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)  DRY_RUN=true;  shift ;;
        --test)     TEST_MODE=true; shift ;;
        --resume)   RESUME_MODE=true; shift ;;
        --cleanup)  CLEANUP=true; shift ;;
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

if [[ "$RESUME_MODE" == "false" && "$CLEANUP" == "true" ]]; then
    echo "ERROR: --cleanup requires --resume"
    exit 1
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
echo "Mode:          $([[ "$RESUME_MODE" == "true" ]] && echo resume || echo full)"
echo "Cleanup first: $([[ "$CLEANUP" == "true" ]] && echo yes || echo no)"
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

TOTAL_JOBS=0
TOTAL_SIMULATIONS=0

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

    ARRAY_SPEC=""
    N_SUBMITTED=$N_REPLICAS
    RUN_DIR_NAME=""

    if [[ "$RESUME_MODE" == "true" ]]; then
        STATUS_ARGS=(
            python3 "$STATUS_SCRIPT"
            --output-base "$OUTPUT_BASE"
            --lambdas "$LAM"
            --n-replicas "$N_REPLICAS"
            --switch-time-ps "$SWITCH_TIME_PS"
            --json
        )
        if [[ "$CLEANUP" == "true" ]]; then
            STATUS_ARGS+=(--cleanup)
            if [[ "$DRY_RUN" == "true" ]]; then
                STATUS_ARGS+=(--dry-run)
            fi
        fi

        STATUS_JSON="$("${STATUS_ARGS[@]}")"
        ARRAY_SPEC=$(STATUS_JSON="$STATUS_JSON" python3 -c 'import json, os; item=json.loads(os.environ["STATUS_JSON"])[0]; print(item["slurm_array"])')
        N_COMPLETE=$(STATUS_JSON="$STATUS_JSON" python3 -c 'import json, os; print(json.loads(os.environ["STATUS_JSON"])[0]["complete"])')
        N_PARTIAL=$(STATUS_JSON="$STATUS_JSON" python3 -c 'import json, os; print(json.loads(os.environ["STATUS_JSON"])[0]["partial"])')
        N_SUBMITTED=$(STATUS_JSON="$STATUS_JSON" python3 -c 'import json, os; print(json.loads(os.environ["STATUS_JSON"])[0]["missing"])')
        RUN_DIR_NAME=$(STATUS_JSON="$STATUS_JSON" python3 -c 'import json, os; print(json.loads(os.environ["STATUS_JSON"])[0]["run_dir"].split("/")[-1])')

        echo ""
        echo "λ=${LAM}: complete=${N_COMPLETE}, partial_removed=${N_PARTIAL}, to_submit=${N_SUBMITTED}"

        if [[ -z "$ARRAY_SPEC" ]]; then
            echo "  Nothing to submit for λ=${LAM}"
            continue
        fi
    else
        ARRAY_END=$((N_REPLICAS - 1))
        ARRAY_SPEC="0-${ARRAY_END}"
        RUN_DIR_NAME=$(PYTHONPATH="${CAV_HOOMD_DIR}" python3 -c "from examples.slurm.aging_campaign_status import coupling_dir_name; print(coupling_dir_name(${LAM}, ${SWITCH_TIME_PS}))")
    fi

    if [[ "$RESUME_MODE" == "true" ]]; then
        SBATCH_SUFFIX="_resume"
    else
        SBATCH_SUFFIX=""
    fi
    SBATCH_FILE="$SCRIPT_DIR/generated/aging_lam${LAM_TAG}${SBATCH_SUFFIX}.sbatch"

    cat > "$SBATCH_FILE" << EOF
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --partition=${PARTITION}
#SBATCH --account=${ACCOUNT}
#SBATCH --gres=${GRES}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --time=${TIME_LIMIT}
#SBATCH --array=${ARRAY_SPEC}%${CONCURRENT}
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
RUN_DIR_NAME="${RUN_DIR_NAME}"
COMPLETE_MIN_BYTES=${COMPLETE_MIN_BYTES}

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

RUN_DIR="${RUN_DIR_NAME}"
H5_FILE="\${RUN_DIR}/observables_replica_\${REPLICA}.h5"
GSD_FILE="\${RUN_DIR}/prod-\${REPLICA}.gsd"

if [[ -f "\${H5_FILE}" ]]; then
    FILE_SIZE=\$(stat -c%s "\${H5_FILE}")
    if [[ "\${FILE_SIZE}" -ge "\${COMPLETE_MIN_BYTES}" ]]; then
        echo "Replica \${REPLICA} already complete (\${FILE_SIZE} bytes), skipping"
        exit 0
    fi
    echo "Removing partial replica artifacts for \${REPLICA}"
    rm -f "\${H5_FILE}" "\${GSD_FILE}"
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
    --fixed-timestep --timestep 1.0 \\
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
    echo "SUCCESS: λ=\${LAM} replica \${REPLICA} completed"
else
    echo "FAILED:  λ=\${LAM} replica \${REPLICA} exit code \${RC}" >&2
    exit \$RC
fi
EOF

    chmod +x "$SBATCH_FILE"
    echo ""
    echo "λ=${LAM}  →  ${JOB_NAME}  array ${ARRAY_SPEC}%${CONCURRENT}"

    if [[ "$DRY_RUN" == "false" ]]; then
        JOB_ID=$(sbatch --parsable "$SBATCH_FILE")
        echo "  Submitted: job $JOB_ID"
        echo "${JOB_ID}:lam=${LAM}:${N_SUBMITTED} replicas mode=${RESUME_MODE}" >> "$SUBMITTED_LOG"
        TOTAL_JOBS=$((TOTAL_JOBS + 1))
        TOTAL_SIMULATIONS=$((TOTAL_SIMULATIONS + N_SUBMITTED))
        sleep 0.5
    else
        echo "  [DRY RUN] Would submit: $SBATCH_FILE"
        TOTAL_JOBS=$((TOTAL_JOBS + 1))
        TOTAL_SIMULATIONS=$((TOTAL_SIMULATIONS + N_SUBMITTED))
    fi
done

echo ""
echo "============================================================"
echo "Total λ values submitted: $TOTAL_JOBS"
echo "Total simulations queued: $TOTAL_SIMULATIONS"
echo "Job IDs logged to:        $SUBMITTED_LOG"
echo ""
echo "Monitor with:"
echo "  squeue -u \$USER | grep cav-aging"
echo "  tail -f $SCRIPT_DIR/logs/cav-aging-lam0_*.out"
echo "============================================================"
