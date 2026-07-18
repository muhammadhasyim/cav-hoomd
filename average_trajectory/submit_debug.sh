#!/bin/bash
#SBATCH --job-name=cavitymd_debug      # Job name with debug suffix
#SBATCH --array=0-9                    # Small array for debugging (10 tasks)
#SBATCH --time=00:30:00                # 30 minutes max
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1              # CPU cores per task
#SBATCH --mem=2G                       # Reduced memory for debug
#SBATCH --output=logs/debug_output_%A_%a.out  # Debug-specific output
#SBATCH --error=logs/debug_error_%A_%a.err    # Debug-specific error

# Exit on any error
set -e

# Get coupling strength from command line argument or environment variable
COUPLING=${1:-${COUPLING:-"1.0e-4"}}

# Create logs directory if it doesn't exist
mkdir -p logs

# Debug parameters - much shorter runtime
runtime=20                    # 10 ps for quick debugging
molecular_tau=1.0
cavity_tau=1.0
switchtime=10

echo "=== DEBUG MODE ==="
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Using coupling strength: $COUPLING"
echo "Frame will be automatically set to: $SLURM_ARRAY_TASK_ID"
echo "DEBUG: Short runtime of $runtime ps"
echo "Molecular tau: $molecular_tau"
echo "Cavity tau: $cavity_tau"
echo "Working directory: $(pwd)"
echo "Python executable: $(which python)"
echo "=================="

# Create a unique seed for each array task to ensure different random states
SEED=0

python 05_advanced_run.py \
  --molecular-bath bussi \
  --cavity-bath langevin \
  --coupling $COUPLING \
  --temperature 100.0 \
  --frequency 1560.0 \
  --molecular-tau $molecular_tau \
  --cavity-tau $cavity_tau \
  --enable-energy-tracker \
  --energy-output-period-ps 0.01 \
  --gsd-output-period-ps 1.0 \
  --console-output-period-ps 1.0 \
  --max-energy-output-time 250.0 \
  --enable-fkt \
  --fkt-ref-interval 100.0 \
  --fkt-max-refs 10 \
  --fkt-output-period-ps 1.0 \
  --enable-dipole-autocorr \
  --dipole-ref-interval 100.0 \
  --dipole-max-refs 10 \
  --dipole-output-period-ps 0.01 \
  --runtime $runtime \
  --device CPU \
  --seed=0 \
  --switch-time $switchtime \
  --error-tolerance 4.0 \
  --initial-fraction 1e-20 \
  --time-constant-ps 100.0


# Check if the simulation completed successfully
if [[ $? -eq 0 ]]; then
    echo "✓ Successfully completed DEBUG task $SLURM_ARRAY_TASK_ID"
    echo "  - Coupling: $COUPLING"
    echo "  - Runtime: ${runtime}ps"
    echo "  - Seed: $SEED"
else
    echo "✗ FAILED: DEBUG task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Optional: Quick output validation
if [[ -f "energy_output.txt" ]]; then
    echo "Energy output file size: $(wc -l < energy_output.txt) lines"
fi

echo "Debug task completed at $(date)"
