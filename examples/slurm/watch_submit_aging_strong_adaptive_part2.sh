#!/bin/bash
# Submit strong-adaptive part2 when MaxSubmitPU allows.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PLAN_DIR="${PLAN_DIR:?PLAN_DIR required}"
PART2_TASKS="${PART2_TASKS:?PART2_TASKS required}"
SBATCH="${PLAN_DIR}/job_part2.sbatch"
MARKER="${PLAN_DIR}/part2_submitted.marker"
LOG="${SCRIPT_DIR}/logs/watch_submit_strong_adaptive_part2.log"
SUBMITTED="${SCRIPT_DIR}/submitted_aging_campaign.txt"
POLL_S="${POLL_INTERVAL_S:-300}"
MAX_SUBMIT=2000
SAFETY=8
USER_NAME="${USER:-mh7373}"

mkdir -p "$(dirname "${LOG}")"
if [[ -f "${MARKER}" ]]; then
    echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) marker exists; exiting" >> "${LOG}"
    exit 0
fi

while true; do
    ts="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    active="$(squeue -u "${USER_NAME}" -h -t PENDING,RUNNING -r | wc -l)"
    projected=$((active + PART2_TASKS))
    limit=$((MAX_SUBMIT - SAFETY))
    echo "${ts} active_tasks=${active} projected=${projected} limit=${limit}" >> "${LOG}"
    if (( projected <= limit )); then
        out="$(sbatch --parsable "${SBATCH}" 2>&1)" || true
        echo "${ts} sbatch: ${out}" >> "${LOG}"
        if [[ "${out}" =~ ^[0-9]+$ ]]; then
            echo "${out}:strong-adaptive-part2:${PART2_TASKS}:walltime=02:00:00" >> "${SUBMITTED}"
            echo "job_id=${out}" > "${MARKER}"
            echo "${ts} submitted part2 as ${out}; exiting" >> "${LOG}"
            exit 0
        fi
    fi
    sleep "${POLL_S}"
done
