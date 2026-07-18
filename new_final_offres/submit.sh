#!/bin/bash
#
# Cavity MD Simulation with Fixed Coupling and Variable Frequency
# 
# Usage: sbatch submit.sh [frequency_in_cm1]
# Example: sbatch submit.sh 1200.0
# 
# Uses a fixed dimensionless coupling constant (lambda_coupling = 0.042)
# and varies the cavity frequency. The lambda_coupling parameter automatically
# handles the frequency dependence internally in the simulation.
#
#SBATCH --job-name=cavitymd             # Job name
#SBATCH --array=0-499                    # Task range (frames 1-50 per coupling)
#SBATCH --time=50:00:00
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1               # CPU cores per task
#SBATCH --mem=4G                        # Memory per job (adjust as needed)
#SBATCH --output=logs/output_%A_%a.out  # Standard output
#SBATCH --error=logs/error_%A_%a.err    # Standard error

# Get frequency from command line argument or environment variable (in cm-1)
FREQUENCY_CM1=${1:-${FREQUENCY_CM1:-"1560.0"}}

# Fixed coupling constant (dimensionless lambda_coupling parameter)
# This value is frequency-independent when using lambda_coupling
COUPLING=0.042

# Create logs directory if it doesn't exist
mkdir -p logs

# Parameters
runtime=3000
molecular_tau=1.0 #ps
cavity_tau=1.0 #ps
switchtime=200.0 #ps

echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Using frequency: $FREQUENCY_CM1 cm-1"
echo "Fixed coupling strength (lambda): $COUPLING (dimensionless)"
echo "Frame will be automatically set to: $SLURM_ARRAY_TASK_ID"

# Run the cavity MD experiment with debug settings
SEED=0

# Run the cavity MD experiment with debug settings
python 05_advanced_run.py \
  --molecular-bath bussi \
  --cavity-bath langevin \
  --coupling $COUPLING \
  --temperature 100.0 \
  --frequency $FREQUENCY_CM1 \
  --molecular-tau $molecular_tau \
  --cavity-tau $cavity_tau \
  --enable-energy-tracker \
  --energy-output-period-ps 0.01 \
  --gsd-output-period-ps 1.0 \
  --console-output-period-ps 1.0 \
  --max-energy-output-time 250.0 \
  --enable-fkt \
  --fkt-ref-interval 200.0 \
  --fkt-max-refs 20 \
  --fkt-output-period-ps 1.0 \
  --enable-dipole-autocorr \
  --dipole-ref-interval 200.0 \
  --dipole-max-refs 20 \
  --dipole-output-period-ps 0.01 \
  --runtime $runtime \
  --device CPU \
  --seed=0 \
  --switch-time $switchtime \
  --error-tolerance 2.0 \
  --initial-fraction 1e-10 \
  --time-constant-ps 1000.0

# Check if the simulation completed successfully
if [[ $? -eq 0 ]]; then
    echo "✓ Successfully completed DEBUG task $SLURM_ARRAY_TASK_ID"
    echo "  - Frequency: $FREQUENCY_CM1 cm-1"
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
