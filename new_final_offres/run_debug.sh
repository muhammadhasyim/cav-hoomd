#!/bin/bash
# Debug job submission script for cavity MD experiments
# Usage: ./run_debug.sh [coupling_value]
# Example: ./run_debug.sh 1e-3

# Default coupling if not provided
COUPLING=${1:-"1.0e-4"}

echo "Submitting DEBUG cavity MD experiment..."
echo "Coupling: $COUPLING"
echo "Array tasks: 1-5"
echo "Runtime: ~5 minutes per task"
echo "Total wall time limit: 15 minutes"
echo ""

# Submit the debug job
sbatch submit_debug.sh $COUPLING

echo "Debug job submitted! Check status with:"
echo "  squeue -u \$USER"
echo ""
echo "Monitor logs in real-time with:"
echo "  tail -f logs/debug_output_*"
echo ""
echo "Quick coupling value examples for testing:"
echo "  ./run_debug.sh 1e-3    # Weak coupling"
echo "  ./run_debug.sh 1e-2    # Medium coupling" 
echo "  ./run_debug.sh 1e-1    # Strong coupling" 