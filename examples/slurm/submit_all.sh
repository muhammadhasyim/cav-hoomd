#!/bin/bash
# ============================================================
# Finite-Size Scaling Study - Master Submission Script
# NYU Greene HPC
# ============================================================
#
# This script generates and submits all SLURM job arrays for the
# finite-size scaling study.
#
# Usage:
#   ./submit_all.sh [options]
#
# Options:
#   --dry-run       Print commands without submitting
#   --test          Submit only 2 replicas per system size for testing
#   --output-dir    Override default output directory
#
# ============================================================

set -e

# ============================================================
# Configuration
# ============================================================

# Physical parameters
COLLECTIVE_COUPLING=0.07      # Target collective coupling g (a.u.)
RUNTIME_PS=3000.0             # Total simulation time (ps)
SWITCH_TIME_PS=1000.0         # Time to turn on coupling (ps)
TEMPERATURE=100.0             # Temperature (K)
FREQUENCY=2000.0              # Cavity frequency (cm^-1)

# Adaptive timestepping parameters
# NOTE: Lower error_tolerance (1.0) needed for stronger coupling (g=0.07)
# This prevents cavity particle from escaping due to large timesteps
ERROR_TOLERANCE=1.0
INITIAL_FRACTION=1e-5
TIME_CONSTANT_PS=50.0

# Output parameters
DISABLE_GSD=true              # Disable GSD trajectory output
ENERGY_OUTPUT_PERIOD=0.1      # Energy output period (ps)
FKT_OUTPUT_PERIOD=0.1         # F(k,t) output period (ps)
FKT_REF_INTERVAL=200.0        # F(k,t) reference interval (ps)

# System sizes and replica counts
# Format: "N_MOLECULES:REPLICAS"
declare -a SYSTEM_CONFIGS=(
    "100:1250"
    "250:500"
    "500:250"
    "1000:125"
    "2500:50"
    "5000:25"
    "10000:13"
)

# SLURM limits
MAX_ARRAY_SIZE=250            # Max tasks per job array (conservative)

# Environment settings
CONDA_ENV="hoomd"
CAV_HOOMD_DIR="/home/mh7373/GitRepos/cav-hoomd"
OUTPUT_DIR="${SCRATCH}/finite_size_scaling_g${COLLECTIVE_COUPLING}_${RUNTIME_PS}ps"

# ============================================================
# Parse command line arguments
# ============================================================

DRY_RUN=false
TEST_MODE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --test)
            TEST_MODE=true
            shift
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# ============================================================
# Setup
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEMPLATE_FILE="${SCRIPT_DIR}/finite_size_scaling_template.sbatch"
JOB_IDS_FILE="${SCRIPT_DIR}/submitted_jobs.txt"

echo "============================================================"
echo "Finite-Size Scaling Study - SLURM Submission"
echo "============================================================"
echo "Date: $(date)"
echo "User: ${USER}"
echo ""
echo "Configuration:"
echo "  Collective coupling: ${COLLECTIVE_COUPLING} a.u."
echo "  Runtime: ${RUNTIME_PS} ps"
echo "  Switch time: ${SWITCH_TIME_PS} ps"
echo "  Output directory: ${OUTPUT_DIR}"
echo "  Conda environment: ${CONDA_ENV}"
echo "  cav-hoomd directory: ${CAV_HOOMD_DIR}"
echo ""

if [ "${DRY_RUN}" = true ]; then
    echo "*** DRY RUN MODE - No jobs will be submitted ***"
    echo ""
fi

if [ "${TEST_MODE}" = true ]; then
    echo "*** TEST MODE - Only 2 replicas per system size ***"
    echo ""
fi

# Create output directory
if [ "${DRY_RUN}" = false ]; then
    mkdir -p "${OUTPUT_DIR}"
    mkdir -p "${SCRIPT_DIR}/logs"
    mkdir -p "${SCRIPT_DIR}/generated"
fi

# Clear previous job IDs file
> "${JOB_IDS_FILE}"

# ============================================================
# Submit job arrays
# ============================================================

TOTAL_JOBS=0
TOTAL_TASKS=0

for config in "${SYSTEM_CONFIGS[@]}"; do
    # Parse configuration
    N_MOLECULES="${config%%:*}"
    REPLICAS="${config##*:}"
    
    # Override replicas in test mode
    if [ "${TEST_MODE}" = true ]; then
        REPLICAS=2
    fi
    
    echo "============================================================"
    echo "System size N=${N_MOLECULES}: ${REPLICAS} replicas"
    echo "============================================================"
    
    # Calculate number of batches needed
    NUM_BATCHES=$(( (REPLICAS + MAX_ARRAY_SIZE - 1) / MAX_ARRAY_SIZE ))
    
    for ((BATCH_ID=1; BATCH_ID<=NUM_BATCHES; BATCH_ID++)); do
        # Calculate array range for this batch
        BATCH_OFFSET=$(( (BATCH_ID - 1) * MAX_ARRAY_SIZE ))
        REMAINING=$(( REPLICAS - BATCH_OFFSET ))
        BATCH_SIZE=$(( REMAINING < MAX_ARRAY_SIZE ? REMAINING : MAX_ARRAY_SIZE ))
        ARRAY_END=$(( BATCH_SIZE - 1 ))
        ARRAY_RANGE="0-${ARRAY_END}"
        
        echo "  Batch ${BATCH_ID}/${NUM_BATCHES}: replicas ${BATCH_OFFSET}-$((BATCH_OFFSET + BATCH_SIZE - 1))"
        
        # Generate sbatch file for this batch
        SBATCH_FILE="${SCRIPT_DIR}/generated/N${N_MOLECULES}_batch${BATCH_ID}.sbatch"
        
        if [ "${DRY_RUN}" = false ]; then
            # Export environment variables and generate sbatch file
            cat > "${SBATCH_FILE}" << EOF
#!/bin/bash
#SBATCH --job-name=cav_N${N_MOLECULES}_b${BATCH_ID}
#SBATCH --output=${SCRIPT_DIR}/logs/cav_N${N_MOLECULES}_b${BATCH_ID}_%A_%a.out
#SBATCH --error=${SCRIPT_DIR}/logs/cav_N${N_MOLECULES}_b${BATCH_ID}_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:v100:1
#SBATCH --mem=16GB
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --array=${ARRAY_RANGE}

# Environment variables
export N_MOLECULES=${N_MOLECULES}
export BATCH_ID=${BATCH_ID}
export BATCH_OFFSET=${BATCH_OFFSET}
export COLLECTIVE_COUPLING=${COLLECTIVE_COUPLING}
export RUNTIME_PS=${RUNTIME_PS}
export SWITCH_TIME_PS=${SWITCH_TIME_PS}
export ERROR_TOLERANCE=${ERROR_TOLERANCE}
export INITIAL_FRACTION=${INITIAL_FRACTION}
export TIME_CONSTANT_PS=${TIME_CONSTANT_PS}
export ENERGY_OUTPUT_PERIOD=${ENERGY_OUTPUT_PERIOD}
export DISABLE_GSD=${DISABLE_GSD}
export ENABLE_FKT=true
export FKT_OUTPUT_PERIOD=${FKT_OUTPUT_PERIOD}
export FKT_REF_INTERVAL=${FKT_REF_INTERVAL}
export OUTPUT_DIR="${OUTPUT_DIR}"
export CONDA_ENV="${CONDA_ENV}"
export CAV_HOOMD_DIR="${CAV_HOOMD_DIR}"

# Source the template script
source "${TEMPLATE_FILE}"
EOF
            chmod +x "${SBATCH_FILE}"
            
            # Submit the job
            JOB_ID=$(sbatch --parsable "${SBATCH_FILE}")
            echo "    Submitted: Job ID ${JOB_ID}"
            echo "${JOB_ID}:N${N_MOLECULES}:batch${BATCH_ID}:${BATCH_OFFSET}-$((BATCH_OFFSET + BATCH_SIZE - 1))" >> "${JOB_IDS_FILE}"
            
            TOTAL_JOBS=$((TOTAL_JOBS + 1))
            TOTAL_TASKS=$((TOTAL_TASKS + BATCH_SIZE))
        else
            echo "    [DRY RUN] Would submit: ${SBATCH_FILE}"
            echo "    [DRY RUN] Array range: ${ARRAY_RANGE}"
            TOTAL_JOBS=$((TOTAL_JOBS + 1))
            TOTAL_TASKS=$((TOTAL_TASKS + BATCH_SIZE))
        fi
    done
done

echo ""
echo "============================================================"
echo "Submission Summary"
echo "============================================================"
echo "Total job arrays submitted: ${TOTAL_JOBS}"
echo "Total tasks (replicas): ${TOTAL_TASKS}"
echo "Job IDs saved to: ${JOB_IDS_FILE}"
echo ""
echo "Monitor with:"
echo "  squeue -u ${USER} | grep cav_"
echo "  ./check_status.sh"
echo ""
echo "Output will be written to:"
echo "  ${OUTPUT_DIR}"
echo "============================================================"
