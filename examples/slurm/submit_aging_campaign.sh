#!/bin/bash
# Submit a robust, manifest-backed packed aging campaign.

set -euo pipefail

LAMBDAS=(0.0 0.01 0.016667 0.023333 0.03)
TARGET_VALID=500
MAX_REPLICA_ID=999
REPLICAS_PER_TASK=2
CONCURRENT_TASKS=4
WALLTIME="03:00:00"
CPUS=8
MEMORY="32G"
PARTITION="a100_chemistry"
ACCOUNT="torch_pr_283_chemistry"
GPU_RESOURCE="gpu:a100:1"
DRY_RUN=false
TEST_MODE=false
PLAN_ID=""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CAV_HOOMD_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
EXAMPLES_DIR="${CAV_HOOMD_DIR}/examples"
OUTPUT_BASE="/scratch/mh7373/projects/cav-hoomd/aging_weak_lambda"
GENERATED_DIR="${SCRIPT_DIR}/generated"
LOG_DIR="${SCRIPT_DIR}/logs"
MINIFORGE="/scratch/mh7373/miniforge3"
CONDA_ENV="hoomd"
INPUT_GSD="${EXAMPLES_DIR}/init-0.gsd"
PACKED_RUNNER="${SCRIPT_DIR}/run_packed_aging_task.py"

usage() {
    printf '%s\n' \
        "Usage: $0 [options]" \
        "" \
        "Options:" \
        "  --dry-run                 Plan and generate files without sbatch" \
        "  --plan-id ID              Immutable filename identifier" \
        "  --test                    Use a four-replica local planning domain" \
        "  --lam VALUE               Plan one coupling value only" \
        "  --output-base PATH        Campaign output directory" \
        "  --generated-dir PATH      Manifest and SBATCH output directory" \
        "  --target-valid N          Valid replicas requested per coupling" \
        "  --max-replica-id N        Highest eligible replica ID" \
        "  --replicas-per-task N     Concurrent replicas per task (must be 2)" \
        "  --concurrent N            Maximum active array tasks" \
        "  --walltime HH:MM:SS       Packed task walltime" \
        "  --resume                  Accepted for backward compatibility" \
        "  --cleanup                 Cleanup occurs safely inside each task" \
        "  --help                    Show this message"
}

require_option_value() {
    if [[ $# -lt 2 || "$2" == --* ]]; then
        echo "ERROR: $1 requires a value" >&2
        exit 2
    fi
}

canonical_path() {
    python3 - "$1" <<'PY'
from pathlib import Path
import sys

print(Path(sys.argv[1]).expanduser().resolve())
PY
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --test)
            TEST_MODE=true
            shift
            ;;
        --plan-id)
            require_option_value "$@"
            PLAN_ID="$2"
            shift 2
            ;;
        --lam)
            require_option_value "$@"
            LAMBDAS=("$2")
            shift 2
            ;;
        --output-base)
            require_option_value "$@"
            OUTPUT_BASE="$2"
            shift 2
            ;;
        --generated-dir)
            require_option_value "$@"
            GENERATED_DIR="$2"
            shift 2
            ;;
        --target-valid)
            require_option_value "$@"
            TARGET_VALID="$2"
            shift 2
            ;;
        --max-replica-id)
            require_option_value "$@"
            MAX_REPLICA_ID="$2"
            shift 2
            ;;
        --replicas-per-task)
            require_option_value "$@"
            REPLICAS_PER_TASK="$2"
            shift 2
            ;;
        --concurrent)
            require_option_value "$@"
            CONCURRENT_TASKS="$2"
            shift 2
            ;;
        --walltime)
            require_option_value "$@"
            WALLTIME="$2"
            shift 2
            ;;
        --resume|--cleanup)
            shift
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ "${TEST_MODE}" == "true" ]]; then
    TARGET_VALID=4
    MAX_REPLICA_ID=3
fi

for lam in "${LAMBDAS[@]}"; do
    case "${lam}" in
        0|0.0|0.01|0.016667|0.023333|0.03)
            ;;
        *)
            echo "ERROR: ${lam} is not an aging campaign coupling" >&2
            exit 2
            ;;
    esac
done

if [[ "${REPLICAS_PER_TASK}" -ne 2 ]]; then
    echo "ERROR: --replicas-per-task must be 2 for this campaign" >&2
    exit 2
fi
if [[
    "${TARGET_VALID}" -lt 0
    || "${MAX_REPLICA_ID}" -lt 0
    || "${MAX_REPLICA_ID}" -gt 999
]]; then
    echo "ERROR: target must be nonnegative and replica IDs must be 0..999" >&2
    exit 2
fi
if [[ "${CONCURRENT_TASKS}" -le 0 ]]; then
    echo "ERROR: --concurrent must be positive" >&2
    exit 2
fi

if [[ -z "${PLAN_ID}" ]]; then
    PLAN_ID="$(date -u +%Y%m%dT%H%M%SZ)-$$"
fi
if [[ ! "${PLAN_ID}" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "ERROR: --plan-id may contain only letters, numbers, dot, _, and -" >&2
    exit 2
fi

OUTPUT_BASE="$(canonical_path "${OUTPUT_BASE}")"
GENERATED_DIR="$(canonical_path "${GENERATED_DIR}")"
LOG_DIR="$(canonical_path "${LOG_DIR}")"
mkdir -p "${GENERATED_DIR}" "${LOG_DIR}" "${OUTPUT_BASE}"

PLAN_DIR="${GENERATED_DIR}/aging_packed_${PLAN_ID}"
if ! mkdir "${PLAN_DIR}"; then
    echo "ERROR: immutable plan directory already exists: ${PLAN_DIR}" >&2
    exit 2
fi
MANIFEST="${PLAN_DIR}/manifest.tsv"
PLAN_JSON="${PLAN_DIR}/plan.json"
SBATCH_FILE="${PLAN_DIR}/job.sbatch"

PLANNER_ARGS=(
    python3 -m examples.slurm.aging_campaign_planner
    --output-base "${OUTPUT_BASE}"
    --lambdas "${LAMBDAS[@]}"
    --target-valid "${TARGET_VALID}"
    --max-replica-id "${MAX_REPLICA_ID}"
    --replicas-per-task "${REPLICAS_PER_TASK}"
    --manifest "${MANIFEST}"
    --json
)

(
    cd "${CAV_HOOMD_DIR}"
    PYTHONPATH="${CAV_HOOMD_DIR}${PYTHONPATH:+:${PYTHONPATH}}" \
        "${PLANNER_ARGS[@]}"
) > "${PLAN_JSON}"

python3 - "${PLAN_JSON}" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    plan = json.load(handle)

print("Packed aging campaign plan")
print("=" * 72)
for item in plan["couplings"]:
    selected = ",".join(str(value) for value in item["selected_replicas"])
    print(
        f"lambda={item['lam']:g}: target={item['target_valid']} "
        f"valid={item['valid']} invalid={item['invalid']} "
        f"needed={item['needed']} selected={item['selected']} "
        f"groups={item['groups']} pairs={item['pairs']} "
        f"singletons={item['singletons']}"
    )
    print(f"  selected_replicas={selected or '(none)'}")
print(f"total_array_tasks={plan['total_groups']}")
print(f"manifest={plan['manifest']}")
PY

TOTAL_GROUPS="$(
    python3 - "${PLAN_JSON}" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    print(json.load(handle)["total_groups"])
PY
)"

echo "resources=gres:${GPU_RESOURCE} cpus:${CPUS} memory:${MEMORY}"
echo "walltime=${WALLTIME} concurrent_tasks=${CONCURRENT_TASKS}"

if [[ "${TOTAL_GROUPS}" -eq 0 ]]; then
    if [[ "${DRY_RUN}" == "true" ]]; then
        echo "Dry run: no work is required and no job was submitted"
        exit 0
    fi
    echo "Refusing to submit an empty plan" >&2
    exit 3
fi

ARRAY_END=$((TOTAL_GROUPS - 1))
printf -v PROFILE_Q '%q' "${MINIFORGE}/etc/profile.d/conda.sh"
printf -v RUNNER_Q '%q' "${PACKED_RUNNER}"
printf -v MANIFEST_Q '%q' "${MANIFEST}"
printf -v OUTPUT_BASE_Q '%q' "${OUTPUT_BASE}"
printf -v EXAMPLES_DIR_Q '%q' "${EXAMPLES_DIR}"
printf -v INPUT_GSD_Q '%q' "${INPUT_GSD}"
printf -v PACKED_LOG_DIR_Q '%q' "${LOG_DIR}/packed"
cat > "${SBATCH_FILE}" <<EOF
#!/bin/bash
#SBATCH --job-name=cav-aging-packed
#SBATCH --partition=${PARTITION}
#SBATCH --account=${ACCOUNT}
#SBATCH --gres=${GPU_RESOURCE}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEMORY}
#SBATCH --time=${WALLTIME}
#SBATCH --array=0-${ARRAY_END}%${CONCURRENT_TASKS}
#SBATCH --output=${LOG_DIR}/cav-aging-packed_%A_%a.out
#SBATCH --error=${LOG_DIR}/cav-aging-packed_%A_%a.err
#SBATCH --mail-type=NONE

set -euo pipefail
source ${PROFILE_Q}
conda activate ${CONDA_ENV}
export OMP_NUM_THREADS=4

python ${RUNNER_Q} \\
    --manifest ${MANIFEST_Q} \\
    --task-index "\${SLURM_ARRAY_TASK_ID}" \\
    --output-base ${OUTPUT_BASE_Q} \\
    --examples-dir ${EXAMPLES_DIR_Q} \\
    --input-gsd ${INPUT_GSD_Q} \\
    --log-dir ${PACKED_LOG_DIR_Q}
EOF

chmod +x "${SBATCH_FILE}"
bash -n "${SBATCH_FILE}"

echo "generated_sbatch=${SBATCH_FILE}"
echo "array=0-${ARRAY_END}%${CONCURRENT_TASKS}"
if [[ "${DRY_RUN}" == "true" ]]; then
    printf 'sbatch --parsable %q\n' "${SBATCH_FILE}"
    echo "Dry run: no job submitted"
    exit 0
fi

JOB_ID="$(sbatch --parsable "${SBATCH_FILE}")"
echo "${JOB_ID}:packed:${TOTAL_GROUPS}:target=${TARGET_VALID}" \
    >> "${SCRIPT_DIR}/submitted_aging_campaign.txt"
echo "Submitted packed aging array: ${JOB_ID}"
