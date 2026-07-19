#!/usr/bin/env bash
# Submit strong-coupling adaptive aging campaign (part 1) with part-2 watcher.
# SLURM resources mirror submit_aging_campaign.sh; only lambdas, output dir,
# and adaptive-dt runner flags differ from the weak-lambda packed campaign.
set -euo pipefail

LAMBDAS=(0.03 0.06 0.09 0.12)
TARGET_VALID=500
MAX_REPLICA_ID=999
REPLICAS_PER_TASK=1
CONCURRENT_TASKS=24
WALLTIME="02:00:00"
CPUS=8
MEMORY="32G"
PARTITION="l40s_public,h200_public,a100_chemistry"
ACCOUNT="torch_pr_283_chemistry"
GPU_RESOURCE="gpu:1"
DRY_RUN=false
PLAN_ONLY=false
PLAN_ID=""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CAV_HOOMD_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTPUT_BASE="/scratch/mh7373/projects/cav-hoomd/aging_strong_lambda_adaptive"
GENERATED_DIR="${SCRIPT_DIR}/generated"
LOG_DIR="${SCRIPT_DIR}/logs"
MINIFORGE="/scratch/mh7373/miniforge3"
CONDA_ENV="hoomd"
INPUT_GSD="${CAV_HOOMD_DIR}/aging_weak_lambda_ic_library/init-500-from-ic-library.gsd"
PACKED_RUNNER="${SCRIPT_DIR}/run_packed_aging_task.py"
ERROR_TOLERANCE=5.0
INITIAL_FRACTION=1e-5
TIME_CONSTANT_PS=50.0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --plan-only) PLAN_ONLY=true; shift ;;
        --plan-id)
            PLAN_ID="$2"
            shift 2
            ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

if [[ "${DRY_RUN}" == "true" ]]; then
    PLAN_ONLY=true
fi

if [[ -z "${PLAN_ID}" ]]; then
    PLAN_ID="strong-adaptive-$(date -u +%Y%m%dT%H%M%SZ)"
fi

mkdir -p "${GENERATED_DIR}" "${LOG_DIR}" "${OUTPUT_BASE}"
PLAN_DIR="${GENERATED_DIR}/aging_packed_${PLAN_ID}"
if [[ -e "${PLAN_DIR}" ]]; then
    echo "ERROR: plan directory exists: ${PLAN_DIR}" >&2
    exit 2
fi
mkdir "${PLAN_DIR}"

MANIFEST="${PLAN_DIR}/manifest.tsv"
PLAN_JSON="${PLAN_DIR}/plan.json"
SBATCH_FILE="${PLAN_DIR}/job.sbatch"

(
    cd "${CAV_HOOMD_DIR}"
    PYTHONPATH="${CAV_HOOMD_DIR}${PYTHONPATH:+:${PYTHONPATH}}" \
        python3 -m examples.slurm.aging_campaign_planner \
        --output-base "${OUTPUT_BASE}" \
        --lambdas "${LAMBDAS[@]}" \
        --target-valid "${TARGET_VALID}" \
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
    echo "Refusing to submit an empty plan" >&2
    exit 3
fi

ARRAY_END=$((TOTAL_GROUPS - 1))
ALLOWED_LAMBDA_CMD=""
for idx in "${!LAMBDAS[@]}"; do
    lam="${LAMBDAS[$idx]}"
    if (( idx < ${#LAMBDAS[@]} - 1 )); then
        ALLOWED_LAMBDA_CMD+="    --allowed-lambda ${lam} \\
"
    else
        ALLOWED_LAMBDA_CMD+="    --allowed-lambda ${lam}"
    fi
done

cat > "${SBATCH_FILE}" <<EOF
#!/bin/bash
#SBATCH --job-name=cav-strong-adp
#SBATCH --partition=${PARTITION}
#SBATCH --account=${ACCOUNT}
#SBATCH --gres=${GPU_RESOURCE}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEMORY}
#SBATCH --time=${WALLTIME}
#SBATCH --array=0-${ARRAY_END}%${CONCURRENT_TASKS}
#SBATCH --output=${LOG_DIR}/cav-strong-adp_%A_%a.out
#SBATCH --error=${LOG_DIR}/cav-strong-adp_%A_%a.err
#SBATCH --mail-type=NONE

set -euo pipefail
source ${MINIFORGE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}
export OMP_NUM_THREADS=4

python ${PACKED_RUNNER} \\
    --manifest ${MANIFEST} \\
    --task-index "\${SLURM_ARRAY_TASK_ID}" \\
    --output-base ${OUTPUT_BASE} \\
    --examples-dir ${CAV_HOOMD_DIR}/examples \\
    --input-gsd ${INPUT_GSD} \\
    --log-dir ${LOG_DIR}/packed_strong_adaptive \\
    --adaptive-timestep \\
    --error-tolerance ${ERROR_TOLERANCE} \\
    --initial-fraction ${INITIAL_FRACTION} \\
    --time-constant-ps ${TIME_CONSTANT_PS} \\
${ALLOWED_LAMBDA_CMD}
EOF

chmod +x "${SBATCH_FILE}"

python3 "${SCRIPT_DIR}/split_packed_manifest.py" \
    --plan-dir "${PLAN_DIR}" \
    --part1-job-name "cav-strong-adp-part1" \
    --part2-job-name "cav-strong-adp-part2"

PART1_TASKS=$(wc -l < "${PLAN_DIR}/manifest_part1.tsv")
PART2_TASKS=$(wc -l < "${PLAN_DIR}/manifest_part2.tsv")

echo "plan_dir=${PLAN_DIR}"
echo "part1_tasks=${PART1_TASKS} part2_tasks=${PART2_TASKS}"
echo "walltime=${WALLTIME} concurrent_tasks=${CONCURRENT_TASKS} replicas_per_task=${REPLICAS_PER_TASK}"

if [[ "${PLAN_ONLY}" == "true" ]]; then
    echo "Plan-only: generated ${PLAN_DIR}/job_part1.sbatch (not submitted)"
    exit 0
fi

JOB_ID="$(sbatch --parsable "${PLAN_DIR}/job_part1.sbatch" 2>&1)" || {
    echo "${JOB_ID}"
    echo "Part 1 blocked by quota; starting part2 watcher."
    nohup env PLAN_DIR="${PLAN_DIR}" PART2_TASKS="${PART2_TASKS}" \
        "${SCRIPT_DIR}/watch_submit_aging_strong_adaptive_part2.sh" \
        >> "${LOG_DIR}/watch_submit_strong_adaptive_part2.nohup.out" 2>&1 &
    echo "Started strong-adaptive part2 watcher PID $!"
    exit 0
}

echo "${JOB_ID}:strong-adaptive-part1:${PART1_TASKS}:walltime=${WALLTIME}" \
    >> "${SCRIPT_DIR}/submitted_aging_campaign.txt"
echo "Submitted part1: ${JOB_ID}"

nohup env PLAN_DIR="${PLAN_DIR}" PART2_TASKS="${PART2_TASKS}" \
    "${SCRIPT_DIR}/watch_submit_aging_strong_adaptive_part2.sh" \
    >> "${LOG_DIR}/watch_submit_strong_adaptive_part2.nohup.out" 2>&1 &
echo "Started strong-adaptive part2 watcher PID $!"
