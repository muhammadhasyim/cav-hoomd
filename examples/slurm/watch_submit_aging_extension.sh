#!/bin/bash
# Poll for checkpoint GSD files and submit extension array when tasks exist.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG="${SCRIPT_DIR}/logs/watch_extension_submit.log"
POLL_S="${POLL_INTERVAL_S:-1800}"
MARKER="${SCRIPT_DIR}/generated/.extension_submitted_marker"

mkdir -p "${SCRIPT_DIR}/logs" "$(dirname "${MARKER}")"
while true; do
  ts="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  if [[ -f "${MARKER}" ]]; then
    echo "${ts} extension already submitted ($(cat "${MARKER}"))" >> "${LOG}"
    exit 0
  fi
  out="$("${SCRIPT_DIR}/submit_aging_extension.sh" --plan-id "ext-${ts}" 2>&1)" || true
  echo "${ts} ${out}" >> "${LOG}"
  if echo "${out}" | grep -q "Submitted aging extension array:"; then
    echo "${ts} $(echo "${out}" | awk '/Submitted aging extension array:/ {print $NF}')" > "${MARKER}"
    exit 0
  fi
  sleep "${POLL_S}"
done
