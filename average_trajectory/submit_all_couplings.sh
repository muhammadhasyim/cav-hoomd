#!/bin/bash

# Script to submit multiple SLURM jobs, each for a different coupling strength

# Array of coupling strengths to test
coupling_values=(
    "1.0e-4"
    "5.0e-4"
    "1.0e-3"
)

echo "Submitting ${#coupling_values[@]} jobs for different coupling strengths..."

# Submit a job for each coupling strength
for coupling in "${coupling_values[@]}"
do
    # Submit job with coupling as argument
    job_id=$(sbatch --parsable submit.sh "$coupling")
    echo "Submitted job $job_id for coupling $coupling"
    
    # Optional: Add a small delay to avoid overwhelming the scheduler
    sleep 1
done

echo "All jobs submitted!"
echo "Monitor with: squeue -u \$USER"
echo "Cancel all with: scancel -u \$USER" 