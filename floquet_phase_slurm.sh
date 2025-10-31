#!/bin/bash
#SBATCH --job-name=floquet_phase
#SBATCH --output=slurm_logs/floquet_%A_%a.out
#SBATCH --error=slurm_logs/floquet_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=cpu
#SBATCH --array=0-19

# Exact-Cancellation + PI Controller SLURM Job Array Script
# ====================================================================
# Runs multiple replicas of the same lambda coupling value
# Each replica is run as a separate SLURM array job on CPU
# 
# Usage:
#   sbatch floquet_phase_slurm.sh <lambda_coupling> <initial_temp> <runtime>
#
# Example:
#   sbatch floquet_phase_slurm.sh 0.025 300.0 20000.0
#
# This will launch 20 replicas (0-19) each with lambda=0.025
# ====================================================================

# Get parameters from command line or use defaults
LAMBDA_COUPLING=${1:-0.025}
INITIAL_TEMP=${2:-300.0}
RUNTIME=${3:-20000.0}

# The replica ID is the SLURM array task ID
REPLICA=${SLURM_ARRAY_TASK_ID}

# Create log directory if it doesn't exist
mkdir -p slurm_logs

echo "========================================================================"
echo "SLURM Job Information"
echo "========================================================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "Started: $(date)"
echo ""
echo "Simulation Parameters:"
echo "  Lambda coupling: ${LAMBDA_COUPLING}"
echo "  Initial temperature: ${INITIAL_TEMP} K"
echo "  Runtime: ${RUNTIME} ps"
echo "  Replica: ${REPLICA}"
echo "  CPUs: ${SLURM_CPUS_PER_TASK}"
echo "========================================================================"
echo ""

# Load necessary modules (adjust based on your cluster)
# module load python/3.12
# module load cuda/12.1  # If you want GPU support
# source activate your_environment

# Set OpenMP threads
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Change to submission directory
cd ${SLURM_SUBMIT_DIR}

# Run the simulation
python3 18_unified_cavity_dynamics.py \
  --lambda-coupling ${LAMBDA_COUPLING} \
  --coupling-type square \
  --square-period 10.0 \
  --square-duty-cycle 0.1 \
  --square-start-time 10.0 \
  --square-stop-time ${RUNTIME} \
  --temperature ${INITIAL_TEMP} \
  --frequency 1560.0 \
  --molecular-tau 1.0 \
  --cavity-tau 1.0 \
  --runtime ${RUNTIME} \
  --device CPU \
  --seed 42 \
  --error-tolerance 5.0 \
  --initial-fraction 1e-5 \
  --time-constant-ps 50.0 \
  --enable-energy-tracker \
  --enable-temp-tracker \
  --enable-hdf5-output \
  --energy-output-period 0.1 \
  --temp-tracker-output-period 0.1 \
  --hdf5-output-period 0.1 \
  --console-output-period 1.0 \
  --input-gsd ../cooling-0.gsd \
  --replicas ${REPLICA} \
  --temp-tracker-empirical-data-file ../potential_energy_vs_T.txt \
  --empirical-data-file ../potential_energy_vs_T.txt \
  --enable-diffeq-controller \
  --diffeq-temperature-method lj_coulombic \
  --diffeq-time-constant 1.0 \
  --diffeq-turn-on-time 10.0 \
  --diffeq-update-interval 0.0 \
  --diffeq-apply-to both \
  --diffeq-T-min 0.01

EXIT_CODE=$?

echo ""
echo "========================================================================"
echo "Job Completed"
echo "========================================================================"
echo "Exit code: ${EXIT_CODE}"
echo "Finished: $(date)"
echo "Output directory: constant_lambda${LAMBDA_COUPLING}_diffeq/"
echo "========================================================================"

exit ${EXIT_CODE}

