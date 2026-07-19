#!/bin/bash
# Start the completion-gated strong-adaptive launcher in the background.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG="${SCRIPT_DIR}/logs/watch_launch_strong_adaptive.nohup.out"
PID_FILE="${SCRIPT_DIR}/logs/watch_launch_strong_adaptive.pid"
MARKER="${SCRIPT_DIR}/generated/strong_adaptive_launched.marker"

if [[ -f "${MARKER}" ]]; then
    echo "Strong-adaptive already launched ($(cat "${MARKER}"))"
    exit 0
fi

if [[ -f "${PID_FILE}" ]]; then
    OLD_PID="$(cat "${PID_FILE}")"
    if kill -0 "${OLD_PID}" 2>/dev/null; then
        echo "Launch watcher already running PID ${OLD_PID}"
        exit 0
    fi
fi

nohup "${SCRIPT_DIR}/watch_launch_strong_adaptive_after_weak_r1.sh" >> "${LOG}" 2>&1 &
echo $! > "${PID_FILE}"
echo "Started strong-adaptive launch watcher PID $(cat "${PID_FILE}")"
echo "Log: ${SCRIPT_DIR}/logs/watch_launch_strong_adaptive.log"
