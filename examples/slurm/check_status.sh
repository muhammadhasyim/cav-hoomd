#!/bin/bash
# ============================================================
# Finite-Size Scaling Study - Job Status Monitor
# NYU Greene HPC
# ============================================================
#
# This script checks the status of submitted job arrays and
# counts completed replicas.
#
# Usage:
#   ./check_status.sh [options]
#
# Options:
#   --output-dir DIR   Specify output directory to scan
#   --detailed         Show per-system-size breakdown
#
# ============================================================

set -e

# ============================================================
# Configuration
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOB_IDS_FILE="${SCRIPT_DIR}/submitted_jobs.txt"

# Default output directory (matches submit_all.sh)
COLLECTIVE_COUPLING=0.02
RUNTIME_PS=3000.0
OUTPUT_DIR="${SCRATCH}/finite_size_scaling_g${COLLECTIVE_COUPLING}_${RUNTIME_PS}ps"

# System configurations (must match submit_all.sh)
declare -a SYSTEM_SIZES=(100 250 500 1000 2500 5000 10000)
declare -A EXPECTED_REPLICAS=(
    [100]=1250
    [250]=500
    [500]=250
    [1000]=125
    [2500]=50
    [5000]=25
    [10000]=13
)

# ============================================================
# Parse command line arguments
# ============================================================

DETAILED=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --detailed)
            DETAILED=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# ============================================================
# Functions
# ============================================================

count_status() {
    local n_molecules=$1
    local status=$2
    local count=0
    
    if [ -d "${OUTPUT_DIR}/N_${n_molecules}" ]; then
        count=$(find "${OUTPUT_DIR}/N_${n_molecules}" -name "status.txt" -exec grep -l "^${status}$" {} \; 2>/dev/null | wc -l)
    fi
    echo $count
}

count_running() {
    local n_molecules=$1
    local count=0
    
    if [ -d "${OUTPUT_DIR}/N_${n_molecules}" ]; then
        # Count directories that exist but don't have a status file yet
        local total_dirs=$(find "${OUTPUT_DIR}/N_${n_molecules}" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l)
        local with_status=$(find "${OUTPUT_DIR}/N_${n_molecules}" -name "status.txt" 2>/dev/null | wc -l)
        count=$((total_dirs - with_status))
    fi
    echo $count
}

# ============================================================
# Main
# ============================================================

echo "============================================================"
echo "Finite-Size Scaling Study - Status Report"
echo "============================================================"
echo "Date: $(date)"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Check SLURM queue
echo "SLURM Queue Status:"
echo "-------------------"
QUEUE_INFO=$(squeue -u ${USER} --name="cav_*" -h -o "%j %t" 2>/dev/null || echo "")
if [ -z "${QUEUE_INFO}" ]; then
    echo "  No jobs in queue"
else
    PENDING=$(echo "${QUEUE_INFO}" | grep -c " PD$" || echo "0")
    RUNNING=$(echo "${QUEUE_INFO}" | grep -c " R$" || echo "0")
    echo "  Pending: ${PENDING}"
    echo "  Running: ${RUNNING}"
fi
echo ""

# Count completed replicas
echo "Replica Completion Status:"
echo "--------------------------"

TOTAL_EXPECTED=0
TOTAL_SUCCESS=0
TOTAL_FAILED=0
TOTAL_RUNNING=0

printf "%-12s %10s %10s %10s %10s %10s\n" "System" "Expected" "Success" "Failed" "Running" "Pending"
printf "%-12s %10s %10s %10s %10s %10s\n" "------" "--------" "-------" "------" "-------" "-------"

for n in "${SYSTEM_SIZES[@]}"; do
    expected=${EXPECTED_REPLICAS[$n]}
    success=$(count_status $n "success")
    failed=$(count_status $n "failure")
    running=$(count_running $n)
    pending=$((expected - success - failed - running))
    if [ $pending -lt 0 ]; then pending=0; fi
    
    printf "%-12s %10d %10d %10d %10d %10d\n" "N=${n}" "$expected" "$success" "$failed" "$running" "$pending"
    
    TOTAL_EXPECTED=$((TOTAL_EXPECTED + expected))
    TOTAL_SUCCESS=$((TOTAL_SUCCESS + success))
    TOTAL_FAILED=$((TOTAL_FAILED + failed))
    TOTAL_RUNNING=$((TOTAL_RUNNING + running))
done

TOTAL_PENDING=$((TOTAL_EXPECTED - TOTAL_SUCCESS - TOTAL_FAILED - TOTAL_RUNNING))
if [ $TOTAL_PENDING -lt 0 ]; then TOTAL_PENDING=0; fi

printf "%-12s %10s %10s %10s %10s %10s\n" "------" "--------" "-------" "------" "-------" "-------"
printf "%-12s %10d %10d %10d %10d %10d\n" "TOTAL" "$TOTAL_EXPECTED" "$TOTAL_SUCCESS" "$TOTAL_FAILED" "$TOTAL_RUNNING" "$TOTAL_PENDING"
echo ""

# Calculate progress
if [ $TOTAL_EXPECTED -gt 0 ]; then
    PROGRESS=$(awk "BEGIN {printf \"%.1f\", 100.0 * ${TOTAL_SUCCESS} / ${TOTAL_EXPECTED}}")
    echo "Overall Progress: ${PROGRESS}% complete (${TOTAL_SUCCESS}/${TOTAL_EXPECTED} replicas)"
else
    echo "Overall Progress: No replicas expected"
fi

# Show failed replicas if any
if [ $TOTAL_FAILED -gt 0 ]; then
    echo ""
    echo "WARNING: ${TOTAL_FAILED} replica(s) failed!"
    echo "Failed replicas:"
    for n in "${SYSTEM_SIZES[@]}"; do
        if [ -d "${OUTPUT_DIR}/N_${n}" ]; then
            find "${OUTPUT_DIR}/N_${n}" -name "status.txt" -exec grep -l "^failure$" {} \; 2>/dev/null | while read f; do
                replica_dir=$(dirname "$f")
                replica_id=$(basename "$replica_dir" | sed 's/replica_//')
                echo "  N=${n}, replica=${replica_id}"
            done
        fi
    done
fi

# Estimate time remaining
if [ $TOTAL_RUNNING -gt 0 ] && [ $TOTAL_PENDING -gt 0 ]; then
    echo ""
    echo "Estimated time remaining: $((TOTAL_PENDING * 35 / TOTAL_RUNNING)) minutes"
    echo "(assuming ~35 min/replica with ${TOTAL_RUNNING} concurrent jobs)"
fi

echo ""
echo "============================================================"

# Detailed view
if [ "${DETAILED}" = true ]; then
    echo ""
    echo "Detailed Logs (last 5 per system size):"
    echo "----------------------------------------"
    for n in "${SYSTEM_SIZES[@]}"; do
        echo ""
        echo "N=${n}:"
        ls -lt ${SCRIPT_DIR}/logs/cav_N${n}_*.out 2>/dev/null | head -5 | awk '{print "  " $NF}' || echo "  (no logs yet)"
    done
fi
