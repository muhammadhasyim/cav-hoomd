#!/bin/bash
# Submit 2000ps-r1 part2 when MaxSubmitPU allows.
# Count array tasks with ``squeue -r`` (compact squeue undercounts).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PLAN_DIR="${SCRIPT_DIR}/generated/aging_packed_2000ps-r1-20260719T072243Z"
SBATCH="${PLAN_DIR}/job_part2.sbatch"
MARKER="${PLAN_DIR}/part2_r1_submitted.marker"
LOG="${SCRIPT_DIR}/logs/watch_submit_part2_2000_r1.log"
SUBMITTED="${SCRIPT_DIR}/submitted_aging_campaign.txt"
POLL_S="${POLL_INTERVAL_S:-300}"
PART2_TASKS=1250
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
            echo "${out}:packed-2000ps-r1-part2:${PART2_TASKS}:walltime=02:00:00:replicas_per_task=1" \
                >> "${SUBMITTED}"
            echo "job_id=${out}" > "${MARKER}"
            echo "${ts} submitted part2 as ${out}; exiting" >> "${LOG}"
            nohup "${SCRIPT_DIR}/start_watch_launch_strong_adaptive.sh" \
                >> "${SCRIPT_DIR}/logs/watch_launch_strong_adaptive.nohup.out" 2>&1 &
            echo "${ts} started strong-adaptive launch watcher PID $!" >> "${LOG}"
            exit 0
        fi
    fi
    sleep "${POLL_S}"
done
