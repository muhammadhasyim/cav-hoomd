#!/bin/bash
# Restart the 2000ps-r1 part2 submit watcher in the background.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG="${SCRIPT_DIR}/logs/watch_submit_part2_2000_r1.nohup.out"
PID_FILE="${SCRIPT_DIR}/logs/watch_submit_part2_2000_r1.pid"
MARKER="${SCRIPT_DIR}/generated/aging_packed_2000ps-r1-20260719T072243Z/part2_r1_submitted.marker"

if [[ -f "${MARKER}" ]]; then
    echo "Part2 already submitted ($(cat "${MARKER}"))"
    exit 0
fi

if [[ -f "${PID_FILE}" ]]; then
    OLD_PID="$(cat "${PID_FILE}")"
    if kill -0 "${OLD_PID}" 2>/dev/null; then
        echo "Watcher already running PID ${OLD_PID}"
        exit 0
    fi
fi

nohup "${SCRIPT_DIR}/watch_submit_aging_2000_r1_part2.sh" >> "${LOG}" 2>&1 &
echo $! > "${PID_FILE}"
echo "Started 2000ps-r1 part2 watcher PID $(cat "${PID_FILE}")"
echo "Log: ${SCRIPT_DIR}/logs/watch_submit_part2_2000_r1.log"
