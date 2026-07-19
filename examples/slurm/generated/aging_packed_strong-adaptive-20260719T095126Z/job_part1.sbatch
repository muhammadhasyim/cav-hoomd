#!/bin/bash
#SBATCH --job-name=cav-strong-adp-part1
#SBATCH --partition=l40s_public,h200_public,a100_chemistry
#SBATCH --account=torch_pr_283_chemistry
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=47:59:00
#SBATCH --array=0-999%24
#SBATCH --output=/scratch/mh7373/projects/cav-hoomd/examples/slurm/logs/cav-strong-adp_%A_%a.out
#SBATCH --error=/scratch/mh7373/projects/cav-hoomd/examples/slurm/logs/cav-strong-adp_%A_%a.err
#SBATCH --mail-type=NONE

set -euo pipefail
source /scratch/mh7373/miniforge3/etc/profile.d/conda.sh
conda activate hoomd
export OMP_NUM_THREADS=4

python /scratch/mh7373/projects/cav-hoomd/examples/slurm/run_packed_aging_task.py \
    --manifest /scratch/mh7373/projects/cav-hoomd/examples/slurm/generated/aging_packed_strong-adaptive-20260719T095126Z/manifest_part1.tsv \
    --task-index "${SLURM_ARRAY_TASK_ID}" \
    --output-base /scratch/mh7373/projects/cav-hoomd/aging_strong_lambda_adaptive \
    --examples-dir /scratch/mh7373/projects/cav-hoomd/examples \
    --input-gsd /scratch/mh7373/projects/cav-hoomd/aging_weak_lambda_ic_library/init-500-from-ic-library.gsd \
    --log-dir /scratch/mh7373/projects/cav-hoomd/examples/slurm/logs/packed_strong_adaptive \
    --adaptive-timestep \
    --error-tolerance 5.0 \
    --initial-fraction 1e-5 \
    --time-constant-ps 50.0 \
    --allowed-lambda 0.03 \
    --allowed-lambda 0.06 \
    --allowed-lambda 0.09 \
    --allowed-lambda 0.12
