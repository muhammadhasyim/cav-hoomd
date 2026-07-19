#!/bin/bash
#SBATCH --job-name=cav-aging-packed
#SBATCH --partition=l40s_public,h200_public,a100_chemistry
#SBATCH --account=torch_pr_283_chemistry
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=0-2499%24
#SBATCH --output=/scratch/mh7373/projects/cav-hoomd/examples/slurm/logs/cav-aging-packed_%A_%a.out
#SBATCH --error=/scratch/mh7373/projects/cav-hoomd/examples/slurm/logs/cav-aging-packed_%A_%a.err
#SBATCH --mail-type=NONE

set -euo pipefail
source /scratch/mh7373/miniforge3/etc/profile.d/conda.sh
conda activate hoomd
export OMP_NUM_THREADS=4

python /scratch/mh7373/projects/cav-hoomd/examples/slurm/run_packed_aging_task.py \
    --manifest /scratch/mh7373/projects/cav-hoomd/examples/slurm/generated/aging_packed_2000ps-r1-20260719T072243Z/manifest.tsv \
    --task-index "${SLURM_ARRAY_TASK_ID}" \
    --output-base /scratch/mh7373/projects/cav-hoomd/aging_weak_lambda \
    --examples-dir /scratch/mh7373/projects/cav-hoomd/examples \
    --input-gsd /scratch/mh7373/projects/cav-hoomd/aging_weak_lambda_ic_library/init-500-from-ic-library.gsd \
    --log-dir /scratch/mh7373/projects/cav-hoomd/examples/slurm/logs/packed
