#!/bin/bash
#SBATCH --job-name=cav-aging-r1-p1r
#SBATCH --partition=l40s_public,h200_public
#SBATCH --account=torch_pr_283_chemistry
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --array=0-1178%24
#SBATCH --output=/scratch/mh7373/projects/cav-hoomd/examples/slurm/logs/cav-aging-r1-p1r_%A_%a.out
#SBATCH --error=/scratch/mh7373/projects/cav-hoomd/examples/slurm/logs/cav-aging-r1-p1r_%A_%a.err
#SBATCH --mail-type=NONE

set -euo pipefail
if [[ "$(hostname -s)" == "ga015" ]]; then
  echo "Refusing ga015: CUDA init failures observed on this node" >&2
  exit 1
fi
source /scratch/mh7373/miniforge3/etc/profile.d/conda.sh
conda activate hoomd
export OMP_NUM_THREADS=4

python /scratch/mh7373/projects/cav-hoomd/examples/slurm/run_packed_aging_task.py \
    --manifest /scratch/mh7373/projects/cav-hoomd/examples/slurm/generated/aging_packed_20260719T051155Z-r1-defer24/manifest_part1_relaunch_exclude_ga015.tsv \
    --task-index "${SLURM_ARRAY_TASK_ID}" \
    --output-base /scratch/mh7373/projects/cav-hoomd/aging_weak_lambda \
    --examples-dir /scratch/mh7373/projects/cav-hoomd/examples \
    --input-gsd /scratch/mh7373/projects/cav-hoomd/aging_weak_lambda_ic_library/init-500-from-ic-library.gsd \
    --log-dir /scratch/mh7373/projects/cav-hoomd/examples/slurm/logs/packed
