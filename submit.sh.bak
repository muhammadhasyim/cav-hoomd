#!/bin/bash
#SBATCH --job-name=cavitymd             # Job name
#SBATCH --array=0-499                    # Task range (frames 1-50 per coupling)
#SBATCH --time=10:00:00                 # Maximum run time (10 hours)
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1               # CPU cores per task
#SBATCH --mem=4G                        # Memory per job (adjust as needed)
#SBATCH --output=logs/output_%A_%a.out  # Standard output
#SBATCH --error=logs/error_%A_%a.err    # Standard error

# Get coupling strength from command line argument or environment variable
COUPLING=${1:-${COUPLING:-"1.0e-4"}}

# Create logs directory if it doesn't exist
mkdir -p logs

# Load any necessary modules
conda init
conda activate hoomd-clean

# Parameters
runtime=1000
molecular_tau=1.0 #ps
cavity_tau=1.0 #ps

echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Using coupling strength: $COUPLING"
echo "Frame will be automatically set to: $SLURM_ARRAY_TASK_ID"

# Run the cavity MD experiment with debug settings
python 05_advanced_run.py --molecular-bath bussi --cavity-bath langevin --coupling $COUPLING --temperature 100.0 --frequency 1560 --runtime $runtime --molecular-tau $molecular_tau --cavity-tau $cavity_tau --enable-energy-tracker --energy-output-period-ps 0.01 --enable-fkt --fkt-output-period-ps 1.0 --gsd-output-period-ps 1.0 --truncate-gsd --fkt-ref-interval 100.0 --max-energy-output-time 1.0 --console-output-period-ps 0.05 --device CPU

echo "Completed task $SLURM_ARRAY_TASK_ID with coupling $COUPLING"
