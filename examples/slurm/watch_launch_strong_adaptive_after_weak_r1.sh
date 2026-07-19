#!/bin/bash
# Wait for weak-lambda 2000ps r1 to finish, then launch strong-adaptive campaign.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WEAK_R1_PLAN_DIR="${WEAK_R1_PLAN_DIR:-${SCRIPT_DIR}/generated/aging_packed_2000ps-r1-20260719T072243Z}"
WEAK_R1_PART1_JOB="${WEAK_R1_PART1_JOB:-14250221}"
PART2_MARKER="${WEAK_R1_PLAN_DIR}/part2_r1_submitted.marker"
LAUNCH_MARKER="${SCRIPT_DIR}/generated/strong_adaptive_launched.marker"
LOG="${SCRIPT_DIR}/logs/watch_launch_strong_adaptive.log"
POLL_S="${POLL_INTERVAL_S:-300}"
USER_NAME="${USER:-mh7373}"
SUBMIT_SCRIPT="${SCRIPT_DIR}/submit_aging_strong_adaptive.sh"

mkdir -p "$(dirname "${LOG}")" "${SCRIPT_DIR}/generated"

if [[ -f "${LAUNCH_MARKER}" ]]; then
    echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) launch marker exists; exiting" >> "${LOG}"
    exit 0
fi

echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) waiting for ${PART2_MARKER}" >> "${LOG}"
while [[ ! -f "${PART2_MARKER}" ]]; do
    sleep "${POLL_S}"
done

WEAK_R1_PART2_JOB="$(
    cd "${SCRIPT_DIR}/../.."
    PYTHONPATH="$(pwd)${PYTHONPATH:+:${PYTHONPATH}}" \
        python3 - "${PART2_MARKER}" <<'PY'
from pathlib import Path
import sys
from examples.slurm.slurm_job_watch import parse_marker_job_id
print(parse_marker_job_id(Path(sys.argv[1]).read_text(encoding="utf-8")))
PY
)"

echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) part2 submitted as ${WEAK_R1_PART2_JOB}; waiting for drain" >> "${LOG}"

while true; do
    ts="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    part1_active="$(squeue -u "${USER_NAME}" -j "${WEAK_R1_PART1_JOB}" -h -t PENDING,RUNNING -r | wc -l)"
    part2_active="$(squeue -u "${USER_NAME}" -j "${WEAK_R1_PART2_JOB}" -h -t PENDING,RUNNING -r | wc -l)"
    echo "${ts} part1_job=${WEAK_R1_PART1_JOB} active=${part1_active} part2_job=${WEAK_R1_PART2_JOB} active=${part2_active}" >> "${LOG}"
    if (( part1_active == 0 && part2_active == 0 )); then
        break
    fi
    sleep "${POLL_S}"
done

launch_ts="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "${launch_ts} weak r1 drained; launching strong-adaptive" >> "${LOG}"
launch_out="$("${SUBMIT_SCRIPT}" 2>&1)" || true
echo "${launch_out}" >> "${LOG}"
echo "launched_at=${launch_ts}" > "${LAUNCH_MARKER}"
echo "${launch_out}"
