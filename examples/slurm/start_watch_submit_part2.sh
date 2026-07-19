#!/bin/bash
# Launch the part 2 submit watcher in the background.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PY:-/scratch/mh7373/miniforge3/envs/hoomd/bin/python}"
PID_FILE="${SCRIPT_DIR}/logs/watch_submit_part2.pid"
LOG_FILE="${SCRIPT_DIR}/logs/watch_submit_part2.log"

mkdir -p "${SCRIPT_DIR}/logs"

if [[ -f "${PID_FILE}" ]]; then
    OLD_PID="$(cat "${PID_FILE}")"
    if kill -0 "${OLD_PID}" 2>/dev/null; then
        echo "Watcher already running with PID ${OLD_PID}"
        echo "Log: ${LOG_FILE}"
        exit 0
    fi
fi

nohup "${PY}" "${SCRIPT_DIR}/watch_submit_aging_part2.py" > /dev/null 2>&1 &
WATCHER_PID=$!
echo "${WATCHER_PID}" > "${PID_FILE}"
echo "Started part 2 watcher PID ${WATCHER_PID}"
echo "Log: ${LOG_FILE}"
