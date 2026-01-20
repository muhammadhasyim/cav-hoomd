#!/bin/bash
#SBATCH --job-name=cav_N${N_MOLECULES}_b${BATCH_ID}
#SBATCH --output=logs/cav_N${N_MOLECULES}_b${BATCH_ID}_%A_%a.out
#SBATCH --error=logs/cav_N${N_MOLECULES}_b${BATCH_ID}_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:v100:1
#SBATCH --mem=16GB
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --array=${ARRAY_RANGE}

# ============================================================
# Finite-Size Scaling Study - SLURM Job Array Template
# NYU Greene HPC
# ============================================================
#
# This script runs a single replica of the finite-size scaling study.
# It is designed to be used as a job array task.
#
# Required environment variables (set by submit_all.sh):
#   N_MOLECULES     - Number of molecules in the system
#   BATCH_ID        - Batch number for this job array
#   BATCH_OFFSET    - Starting replica ID for this batch
#   COLLECTIVE_COUPLING - Target collective coupling g (a.u.)
#   OUTPUT_DIR      - Base output directory
#   CONDA_ENV       - Conda environment name
#   CAV_HOOMD_DIR   - Path to cav-hoomd repository
#
# ============================================================

set -e

# Print job information
echo "============================================================"
echo "SLURM Job Array Task"
echo "============================================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "============================================================"

# Calculate actual replica ID
REPLICA_ID=$((BATCH_OFFSET + SLURM_ARRAY_TASK_ID))
echo "System size N: ${N_MOLECULES}"
echo "Batch ID: ${BATCH_ID}"
echo "Batch offset: ${BATCH_OFFSET}"
echo "Replica ID: ${REPLICA_ID}"
echo "============================================================"

# Load conda environment
echo "Loading conda environment: ${CONDA_ENV}"
source ~/.bashrc
conda activate ${CONDA_ENV}

# Verify CUDA is available
echo "CUDA devices:"
nvidia-smi --query-gpu=index,name,memory.total --format=csv
echo "============================================================"

# Set working directory
cd ${CAV_HOOMD_DIR}/examples

# Check if this replica has already completed successfully
REPLICA_OUTPUT_DIR="${OUTPUT_DIR}/N_${N_MOLECULES}/replica_${REPLICA_ID}"
REPLICA_STATUS_FILE="${REPLICA_OUTPUT_DIR}/status.txt"

if [ -f "${REPLICA_STATUS_FILE}" ]; then
    PREV_STATUS=$(cat "${REPLICA_STATUS_FILE}")
    if [ "${PREV_STATUS}" == "success" ]; then
        echo "Replica ${REPLICA_ID} already completed successfully. Skipping."
        exit 0
    else
        echo "Previous run status: ${PREV_STATUS}. Re-running..."
    fi
fi

# Create output directory for this replica
mkdir -p "${REPLICA_OUTPUT_DIR}"

# Run the simulation
echo "Starting simulation..."
echo "Command: python finite_size_scaling_study.py \\"
echo "    --single-replica-mode \\"
echo "    --replica-id ${REPLICA_ID} \\"
echo "    --n-molecules-single ${N_MOLECULES} \\"
echo "    --collective-coupling ${COLLECTIVE_COUPLING} \\"
echo "    --runtime-ps ${RUNTIME_PS} \\"
echo "    --switch-time-ps ${SWITCH_TIME_PS} \\"
echo "    --device GPU \\"
echo "    --output-dir ${OUTPUT_DIR}"
echo "============================================================"

# Trap errors and write status
trap 'echo "failure" > "${REPLICA_STATUS_FILE}"; echo "Job failed with exit code $?"' ERR

# Build command with optional flags
CMD="python finite_size_scaling_study.py \
    --single-replica-mode \
    --replica-id ${REPLICA_ID} \
    --n-molecules-single ${N_MOLECULES} \
    --collective-coupling ${COLLECTIVE_COUPLING} \
    --runtime-ps ${RUNTIME_PS} \
    --switch-time-ps ${SWITCH_TIME_PS} \
    --device GPU \
    --gpu-id 0 \
    --output-dir ${OUTPUT_DIR} \
    --molecular-thermostat bussi \
    --cavity-thermostat langevin \
    --molecular-thermostat-tau 5.0 \
    --cavity-thermostat-tau 5.0 \
    --error-tolerance ${ERROR_TOLERANCE:-5.0} \
    --initial-fraction ${INITIAL_FRACTION:-1e-5} \
    --time-constant-ps ${TIME_CONSTANT_PS:-50.0} \
    --energy-output-period-ps ${ENERGY_OUTPUT_PERIOD:-0.1} \
    --enable-hdf5-output"

# Add F(k,t) if enabled
if [ "${ENABLE_FKT:-true}" = "true" ]; then
    CMD="${CMD} --enable-fkt"
    CMD="${CMD} --fkt-output-period-ps ${FKT_OUTPUT_PERIOD:-0.1}"
    CMD="${CMD} --fkt-ref-interval-ps ${FKT_REF_INTERVAL:-200.0}"
fi

# Disable GSD if requested
if [ "${DISABLE_GSD:-true}" = "true" ]; then
    CMD="${CMD} --disable-gsd"
fi

# Run the Python script
eval ${CMD}

# Mark as successful
echo "success" > "${REPLICA_STATUS_FILE}"

echo "============================================================"
echo "Simulation completed successfully!"
echo "Output directory: ${REPLICA_OUTPUT_DIR}"
echo "End time: $(date)"
echo "============================================================"
