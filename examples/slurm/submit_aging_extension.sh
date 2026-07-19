#!/bin/bash
# Submit checkpointed 1600→2000 ps aging extensions.

set -euo pipefail

LAMBDAS=(0.0 0.01 0.016667 0.023333 0.03)
MAX_REPLICA_ID=999
REPLICAS_PER_TASK=1
CONCURRENT_TASKS=24
WALLTIME="00:45:00"
CPUS=8
MEMORY="32G"
PARTITION="l40s_public,h200_public,a100_chemistry"
ACCOUNT="torch_pr_283_chemistry"
GPU_RESOURCE="gpu:1"
DRY_RUN=false
PLAN_ID=""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CAV_HOOMD_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTPUT_BASE="/scratch/mh7373/projects/cav-hoomd/aging_weak_lambda"
GENERATED_DIR="${SCRIPT_DIR}/generated"
LOG_DIR="${SCRIPT_DIR}/logs"
MINIFORGE="/scratch/mh7373/miniforge3"
CONDA_ENV="hoomd"
EXTENSION_RUNNER="${SCRIPT_DIR}/run_extension_aging_task.py"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --plan-id) PLAN_ID="$2"; shift 2 ;;
        --output-base) OUTPUT_BASE="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

if [[ -z "${PLAN_ID}" ]]; then
    PLAN_ID="$(date -u +%Y%m%dT%H%M%SZ)-ext-$$"
fi

PLAN_DIR="${GENERATED_DIR}/aging_extension_${PLAN_ID}"
mkdir -p "${PLAN_DIR}" "${LOG_DIR}" "${OUTPUT_BASE}"
MANIFEST="${PLAN_DIR}/manifest.tsv"
PLAN_JSON="${PLAN_DIR}/plan.json"
SBATCH_FILE="${PLAN_DIR}/job.sbatch"

(
    cd "${CAV_HOOMD_DIR}"
    PYTHONPATH="${CAV_HOOMD_DIR}${PYTHONPATH:+:${PYTHONPATH}}" \
        python3 -m examples.slurm.aging_extension_planner \
        --output-base "${OUTPUT_BASE}" \
        --lambdas "${LAMBDAS[@]}" \
        --max-replica-id "${MAX_REPLICA_ID}" \
        --replicas-per-task "${REPLICAS_PER_TASK}" \
        --manifest "${MANIFEST}" \
        --json
) > "${PLAN_JSON}"

TOTAL_GROUPS="$(
    python3 - "${PLAN_JSON}" <<'PY'
import json, sys
print(json.load(open(sys.argv[1], encoding="utf-8"))["total_groups"])
PY
)"

if [[ "${TOTAL_GROUPS}" -eq 0 ]]; then
    echo "No extension replicas require submission"
    exit 0
fi

ARRAY_END=$((TOTAL_GROUPS - 1))
cat > "${SBATCH_FILE}" <<EOF
#!/bin/bash
#SBATCH --job-name=cav-aging-ext
#SBATCH --partition=${PARTITION}
#SBATCH --account=${ACCOUNT}
#SBATCH --gres=${GPU_RESOURCE}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEMORY}
#SBATCH --time=${WALLTIME}
#SBATCH --array=0-${ARRAY_END}%${CONCURRENT_TASKS}
#SBATCH --output=${LOG_DIR}/cav-aging-ext_%A_%a.out
#SBATCH --error=${LOG_DIR}/cav-aging-ext_%A_%a.err
#SBATCH --mail-type=NONE

set -euo pipefail
source ${MINIFORGE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}
export OMP_NUM_THREADS=4

python ${EXTENSION_RUNNER} \\
    --manifest ${MANIFEST} \\
    --task-index "\${SLURM_ARRAY_TASK_ID}" \\
    --output-base ${OUTPUT_BASE} \\
    --examples-dir ${CAV_HOOMD_DIR}/examples \\
    --log-dir ${LOG_DIR}/extension
EOF

chmod +x "${SBATCH_FILE}"
if [[ "${DRY_RUN}" == "true" ]]; then
    echo "Dry run: ${SBATCH_FILE}"
    exit 0
fi

JOB_ID="$(sbatch --parsable "${SBATCH_FILE}")"
echo "Submitted aging extension array: ${JOB_ID}"
