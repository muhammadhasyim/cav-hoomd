#!/bin/bash

# Unified Exponential Wave vs Zero Coupling Benchmark Script (Fair Comparison)
# Aligned with exponential_wave_command.sh configuration
# 
# This script runs two simulations sequentially to test controller effectiveness:
# 1. Exponential wave coupling (3e-4) with DiffEq + Adaptive Bath controllers
# 2. Zero coupling (0.0) with GD controller targeting the same final temperature
#
# FAIR COMPARISON DESIGN:
#   - Both use same thermostats (Bussi + Langevin)
#   - Both use same adaptive timestep (error_tolerance=2.0)
#   - Both start from same initial configuration (cooling-0.gsd)
#   - Both track molecular temperatures (T_trans, T_rot, T_vib)
#   - Exp wave uses: DiffEq + Adaptive Bath + cavity coupling
#   - Zero coupling uses: GD controller only, NO cavity coupling
#   - Question: Can GD controller reach same T without cavity help?
#
# Usage:
#   ./exponential_vs_zero_benchmark.sh [num_replicas] [device] [runtime_ps] [diffeq_turn_on_ps] [diffeq_turn_off_ps] [adaptive_turn_on_ps] [gd_turn_on_ps]
#
# Examples:
#   ./exponential_vs_zero_benchmark.sh                                        # Default: 1 replica, GPU, 3000ps
#   ./exponential_vs_zero_benchmark.sh 3 CPU                                 # 3 replicas on CPU
#   ./exponential_vs_zero_benchmark.sh 1 GPU 5000 500.0 1000.0 1000.0 500.0  # 5000ps runtime
#   ./exponential_vs_zero_benchmark.sh 5 GPU 10000 1000.0 2000.0 2000.0 1000.0 # Long run

set -e  # Exit on any error

# Help function
show_help() {
    echo "Unified Exponential Wave vs Zero Coupling Benchmark Script"
    echo "Aligned with exponential_wave_command.sh configuration"
    echo ""
    echo "Usage:"
    echo "  $0 [num_replicas] [device] [runtime_ps] [diffeq_turn_on_ps] [diffeq_turn_off_ps] [adaptive_turn_on_ps] [gd_turn_on_ps]"
    echo ""
    echo "Parameters:"
    echo "  num_replicas        Number of simulation replicas (default: 1)"
    echo "  device              Compute device: CPU or GPU (default: GPU)"
    echo "  runtime_ps          Total simulation time in picoseconds (default: 3000.0)"
    echo "  diffeq_turn_on_ps   DiffEq controller turn-on time in ps (default: 500.0)"
    echo "  diffeq_turn_off_ps  DiffEq controller turn-off time in ps (default: 1000.0)"
    echo "  adaptive_turn_on_ps Adaptive Bath controller turn-on time in ps (default: 1000.0)"
    echo "  gd_turn_on_ps       GD controller turn-on time for zero coupling (default: 500.0)"
    echo ""
    echo "Examples:"
    echo "  $0                                           # Default: 1 replica, GPU, 3000ps"
    echo "  $0 3 CPU                                    # 3 replicas on CPU"
    echo "  $0 1 GPU 5000 500.0 1000.0 1000.0 500.0     # 5000ps runtime"
    echo "  $0 5 GPU 10000 1000.0 2000.0 2000.0 1000.0  # Long run: 10000ps"
    echo ""
    echo "Controller Timeline (Exponential Wave):"
    echo "  Phase 1: DiffEq Controller + Exp Wave Coupling (diffeq_turn_on → diffeq_turn_off)"
    echo "  Phase 2: Adaptive Bath Controller (adaptive_turn_on → end)"
    echo "  Coupling: 3e-4, period=50ps, tau=1.0ps, amplitude_scale=5.0"
    echo "  Time constant: 0.1 ps (same for both DiffEq and Adaptive)"
    echo ""
    echo "Zero Coupling (Fair Comparison):"
    echo "  Identical settings to exponential wave, but:"
    echo "  - coupling = 0.0 (no cavity)"
    echo "  - GD controller with time_constant=0.1 ps (same as exp wave controllers)"
    echo "  - GD target = final temperature from exp wave simulation"
    echo "  - Question: Can GD reach same T without cavity help?"
    echo ""
    exit 0
}

# Check for help request
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_help
fi

# Configuration with defaults (aligned with exponential_wave_command.sh)
NUM_REPLICAS=${1:-1}        # Default to 1 replica if not specified
DEVICE=${2:-GPU}            # Default to GPU if not specified
RUNTIME_PS=${3:-3000.0}     # Default to 3000.0 ps runtime (matches wave command)
DIFFEQ_TURN_ON_PS=${4:-500.0}    # Default DiffEq controller turn-on time (matches exp-start-time)
DIFFEQ_TURN_OFF_PS=${5:-1000.0}  # Default DiffEq controller turn-off time (matches exp-stop-time)
ADAPTIVE_TURN_ON_PS=${6:-1000.0} # Default Adaptive Bath controller turn-on time (matches exp-stop-time)
GD_TURN_ON_PS=${7:-500.0}        # Default GD controller turn-on time for zero coupling (match wave timing)

echo "=========================================="
echo "EXPONENTIAL WAVE vs ZERO COUPLING BENCHMARK"
echo "=========================================="
echo "Configuration:"
echo "  Number of replicas: $NUM_REPLICAS"
echo "  Device: $DEVICE"
echo "  Runtime: $RUNTIME_PS ps each simulation"
echo ""
echo "Exponential Wave Simulation:"
echo "  Coupling: 3e-4 (exp wave: period=50ps, tau=1.0ps)"
echo "  Controllers: DiffEq (${DIFFEQ_TURN_ON_PS}-${DIFFEQ_TURN_OFF_PS}ps) → Adaptive Bath (${ADAPTIVE_TURN_ON_PS}ps+)"
echo "  Time constant: 0.1 ps"
echo ""
echo "Zero Coupling Simulation (Fair Comparison):"
echo "  Coupling: 0.0 (no cavity)"
echo "  Controller: GD (kinetic temp tracking, ON at ${GD_TURN_ON_PS}ps)"
echo "  Time constant: 0.1 ps (same as exp wave)"
echo "  Target: Final temperature from exp wave simulation"
echo ""

# Generate replica range (0-N)
generate_replica_range() {
    local num_replicas=$1
    if [[ $num_replicas -eq 1 ]]; then
        echo "0"
    else
        local last_replica=$((num_replicas - 1))
        echo "0-$last_replica"
    fi
}

REPLICA_RANGE=$(generate_replica_range $NUM_REPLICAS)
echo "Replica range: $REPLICA_RANGE"

# Validate timing parameters
echo ""
echo "Validating timing parameters..."

# Check that turn-off time is after turn-on time for DiffEq controller
if (( $(echo "$DIFFEQ_TURN_OFF_PS <= $DIFFEQ_TURN_ON_PS" | bc -l) )); then
    echo "ERROR: DiffEq turn-off time ($DIFFEQ_TURN_OFF_PS ps) must be after turn-on time ($DIFFEQ_TURN_ON_PS ps)"
    exit 1
fi

# Check that adaptive controller starts when DiffEq stops (or later)
if (( $(echo "$ADAPTIVE_TURN_ON_PS < $DIFFEQ_TURN_OFF_PS" | bc -l) )); then
    echo "WARNING: Adaptive Bath controller starts ($ADAPTIVE_TURN_ON_PS ps) before DiffEq stops ($DIFFEQ_TURN_OFF_PS ps)"
    echo "         This may cause controller conflicts. Consider setting adaptive turn-on >= DiffEq turn-off."
fi

# Check that all turn-on times are within the runtime
for time_var in DIFFEQ_TURN_ON_PS DIFFEQ_TURN_OFF_PS ADAPTIVE_TURN_ON_PS GD_TURN_ON_PS; do
    time_val=${!time_var}
    if (( $(echo "$time_val >= $RUNTIME_PS" | bc -l) )); then
        echo "ERROR: $time_var ($time_val ps) must be less than runtime ($RUNTIME_PS ps)"
        exit 1
    fi
done

echo "✓ Timing parameters validated successfully"

# Function to extract final temperature from CSV
extract_final_temperature() {
    local csv_file="$1"
    local temp_column="$2"
    
    if [[ ! -f "$csv_file" ]]; then
        echo "Error: CSV file $csv_file not found!" >&2
        return 1
    fi
    
    # Get the last line and extract the specified temperature column
    local final_temp=$(tail -1 "$csv_file" | cut -d',' -f"$temp_column")
    
    # Validate it's a number
    if [[ ! "$final_temp" =~ ^[0-9]+\.?[0-9]*$ ]]; then
        echo "Error: Could not extract valid temperature from $csv_file column $temp_column" >&2
        echo "Got: '$final_temp'" >&2
        return 1
    fi
    
    echo "$final_temp"
}

# Function to extract and average final temperatures across multiple replicas
extract_average_temperature_from_replicas() {
    local base_csv_pattern="$1"  # e.g., "temperature_tracker_replica_"
    local temp_column="$2"
    local num_replicas="$3"
    
    local sum=0
    local count=0
    local temps=()
    
    for ((i=0; i<num_replicas; i++)); do
        local csv_file="${base_csv_pattern}${i}.csv"
        if [[ -f "$csv_file" ]]; then
            local temp=$(extract_final_temperature "$csv_file" "$temp_column")
            if [[ -n "$temp" ]]; then
                temps+=("$temp")
                sum=$(python3 -c "print($sum + $temp)")
                count=$((count + 1))
            fi
        else
            echo "Warning: Replica $i CSV file not found: $csv_file" >&2
        fi
    done
    
    if [[ $count -eq 0 ]]; then
        echo "Error: No valid temperature data found across replicas!" >&2
        return 1
    fi
    
    local avg_temp=$(python3 -c "print(f'{$sum / $count:.2f}')")
    echo "Replica temperatures: ${temps[*]}"
    echo "Average temperature: $avg_temp K (from $count replicas)"
    echo "$avg_temp"
}

# Function to wait for simulation completion
wait_for_completion() {
    local experiment_name="$1"
    
    echo "Waiting for $experiment_name to complete..."
    
    # Monitor the simulation by checking if the process is still running
    local start_time=$(date +%s)
    while true; do
        # Check if there are any python processes running the simulation
        if ! pgrep -f "18_unified_cavity_dynamics.py" > /dev/null; then
            echo "$experiment_name completed!"
            break
        fi
        
        # Show progress every 60 seconds
        local current_time=$(date +%s)
        local elapsed=$((current_time - start_time))
        if [[ $((elapsed % 60)) -eq 0 ]] && [[ $elapsed -gt 0 ]]; then
            echo "  [$experiment_name] Running for ${elapsed}s..."
            # Count running processes
            local running_procs=$(pgrep -cf "18_unified_cavity_dynamics.py" || echo "0")
            echo "  Active python processes: $running_procs"
        fi
        
        sleep 1
    done
    
    local end_time=$(date +%s)
    local total_time=$((end_time - start_time))
    echo "$experiment_name finished in ${total_time} seconds"
}

# Step 1: Run exponential wave experiment
echo "=========================================="
echo "STEP 1: Running Exponential Wave Experiment"
echo "=========================================="
echo "Configuration:"
echo "  - Coupling: 3e-4 (exponential wave, period=50ps, tau=1.0ps)"
echo "  - Phase 1 (${DIFFEQ_TURN_ON_PS}-${DIFFEQ_TURN_OFF_PS} ps): DiffEq Controller + Exponential Wave"
echo "  - Phase 2 (${ADAPTIVE_TURN_ON_PS}+ ps): Adaptive Bath Controller (amplitude_scale=5.0)"
echo "  - Runtime: $RUNTIME_PS ps"
echo "  - Replicas: $REPLICA_RANGE"
echo "  - Molecular temps: Enabled"
echo ""

# Run exponential wave simulation directly (aligned with exponential_wave_command.sh)
python3 18_unified_cavity_dynamics.py \
  --molecular-bath bussi \
  --cavity-bath langevin \
  --coupling-type exponentialwave \
  --coupling 3e-4 \
  --exp-period 50.0 \
  --exp-tau 1.0 \
  --exp-start-time $DIFFEQ_TURN_ON_PS \
  --exp-stop-time $DIFFEQ_TURN_OFF_PS \
  --temperature 300.0 \
  --frequency 1560.0 \
  --molecular-tau 1.0 \
  --cavity-tau 1.0 \
  --runtime $RUNTIME_PS \
  --enable-energy-tracker \
  --enable-temp-tracker \
  --temp-tracker-empirical-data-file /home/mh7373/GitRepos/cav-hoomd/potential_energy_vs_T.txt \
  --enable-diffeq-controller \
  --diffeq-temperature-method lj_coulombic \
  --diffeq-time-constant 0.1 \
  --diffeq-turn-on-time $DIFFEQ_TURN_ON_PS \
  --diffeq-turn-off-time $DIFFEQ_TURN_OFF_PS \
  --diffeq-update-interval 0.0 \
  --diffeq-apply-to both \
  --diffeq-T-min 0.1 \
  --enable-adaptive-bath \
  --adaptive-bath-amplitude-scale 10.0 \
  --adaptive-bath-time-constant 0.1 \
  --adaptive-bath-turn-on-time $ADAPTIVE_TURN_ON_PS \
  --adaptive-bath-update-interval 0.0 \
  --adaptive-bath-T-min 0.1 \
  --adaptive-bath-dynamic-target \
  --adaptive-bath-apply-to both \
  --adaptive-bath-empirical-data-file /home/mh7373/GitRepos/cav-hoomd/potential_energy_vs_T.txt \
  --adaptive-bath-signal-temperature-method lj_coulombic \
  --device $DEVICE \
  --seed 42 \
  --replicas $REPLICA_RANGE \
  --energy-output-period-ps 0.01 \
  --console-output-period-ps 1.0 \
  --error-tolerance 2.0 \
  --initial-fraction 1e-6 \
  --input-gsd cooling-0.gsd \
  --time-constant-ps 10.0 \
  --enable-dynamic-coupling-detection \
  --enable-molecular-temps \
  --molecular-temps-output-period-ps 1.0 \
  --coupling-change-threshold 1e-5 &

# Wait for completion
wait_for_completion "Exponential Wave Experiment"

# Step 2: Extract final equilibrium temperature
echo ""
echo "=========================================="
echo "STEP 2: Extracting Final Equilibrium Temperature"
echo "=========================================="

# Find the exponential wave output directory (coupling 3e-4)
EXP_OUTPUT_DIR=$(find . -maxdepth 1 -name "*cavity_coupling_3eneg04*" -type d | head -1)

if [[ -z "$EXP_OUTPUT_DIR" ]]; then
    echo "ERROR: Could not find exponential wave output directory matching pattern: *cavity_coupling_3eneg04*"
    echo "Available cavity directories:"
    find . -maxdepth 1 -name "*cavity*" -type d 2>/dev/null || echo "  No cavity directories found"
    exit 1
fi

echo "Found exponential wave output directory: $EXP_OUTPUT_DIR"

# Extract final temperature from temperature tracker files
FINAL_TEMP_K=""
if [[ -f "${EXP_OUTPUT_DIR}/temperature_tracker_replica_0.csv" ]]; then
    echo "Found temperature tracker CSV files, extracting average from replicas..."
    
    # Extract cavity bath temperatures (column 5) across all replicas
    CAVITY_TEMP_BASE="${EXP_OUTPUT_DIR}/temperature_tracker_replica_"
    CAVITY_TEMP=$(extract_average_temperature_from_replicas "$CAVITY_TEMP_BASE" 5 $NUM_REPLICAS | tail -1)
    
    # Extract molecular bath temperatures (column 6) across all replicas  
    MOLECULAR_TEMP=$(extract_average_temperature_from_replicas "$CAVITY_TEMP_BASE" 6 $NUM_REPLICAS | tail -1)
    
    if [[ -n "$MOLECULAR_TEMP" && -n "$CAVITY_TEMP" ]]; then
        # Use average of molecular and cavity bath temperatures
        FINAL_TEMP_K=$(python3 -c "print(f'{(float('$MOLECULAR_TEMP') + float('$CAVITY_TEMP')) / 2.0:.2f}')")
        echo ""
        echo "Final average cavity bath temperature: ${CAVITY_TEMP} K"
        echo "Final average molecular bath temperature: ${MOLECULAR_TEMP} K" 
        echo "Overall final temperature: ${FINAL_TEMP_K} K"
    fi
else
    echo "Error: Could not find temperature tracker CSV files!"
    echo "Expected files pattern: ${EXP_OUTPUT_DIR}/temperature_tracker_replica_*.csv"
    echo "Current directory contents:"
    ls -la "$EXP_OUTPUT_DIR/" 2>/dev/null || echo "Directory not found"
    exit 1
fi

if [[ -z "$FINAL_TEMP_K" ]]; then
    echo "Error: Could not extract final temperature!"
    exit 1
fi

echo ""
echo "Extracted equilibrium temperature: $FINAL_TEMP_K K"

# Step 3: Run zero coupling experiment with extracted temperature
echo ""
echo "=========================================="
echo "STEP 3: Running Zero Coupling Experiment"
echo "=========================================="
echo "Configuration:"
echo "  - Coupling: 0.0 (no cavity coupling)"
echo "  - GD Controller: kinetic temperature tracking (ON at ${GD_TURN_ON_PS}ps)"
echo "  - Target temperature: $FINAL_TEMP_K K (from exponential wave)"
echo "  - Time constant: 0.1 ps (same as exp wave controllers)"
echo "  - Runtime: $RUNTIME_PS ps"
echo "  - Replicas: $REPLICA_RANGE"
echo "  - Molecular temps: Enabled"
echo ""

# Run zero coupling simulation with GD controller targeting extracted temperature
python3 18_unified_cavity_dynamics.py \
  --molecular-bath bussi \
  --cavity-bath langevin \
  --coupling-type constant \
  --coupling 0.0 \
  --temperature 300.0 \
  --frequency 1560.0 \
  --molecular-tau 1.0 \
  --cavity-tau 1.0 \
  --runtime $RUNTIME_PS \
  --enable-energy-tracker \
  --enable-temp-tracker \
  --temp-tracker-empirical-data-file /home/mh7373/GitRepos/cav-hoomd/potential_energy_vs_T.txt \
  --enable-gd-feedback \
  --gd-method kinetic \
  --gd-target $FINAL_TEMP_K \
  --gd-time-constant 0.1 \
  --gd-turn-on-time $GD_TURN_ON_PS \
  --gd-update-interval 0.0 \
  --gd-apply-to both \
  --gd-T-min 0.1 \
  --device $DEVICE \
  --seed 42 \
  --replicas $REPLICA_RANGE \
  --energy-output-period-ps 0.01 \
  --console-output-period-ps 1.0 \
  --error-tolerance 2.0 \
  --initial-fraction 1e-6 \
  --input-gsd cooling-0.gsd \
  --time-constant-ps 10.0 \
  --enable-dynamic-coupling-detection \
  --enable-molecular-temps \
  --molecular-temps-output-period-ps 1.0 \
  --coupling-change-threshold 1e-5 &

# Wait for completion
wait_for_completion "Zero Coupling Experiment"

# Step 4: Generate benchmark summary
echo ""
echo "=========================================="
echo "STEP 4: Benchmark Summary"
echo "=========================================="

echo "Exponential Wave vs Zero Coupling benchmark complete!"
echo ""
echo "FAIR COMPARISON SUMMARY:"
echo "========================================"
echo "Both simulations use IDENTICAL settings:"
echo "  - Same thermostats: Bussi (molecular) + Langevin (cavity)"
echo "  - Same initial temperature: 300K"
echo "  - Same runtime: $RUNTIME_PS ps"
echo "  - Same timestep strategy: adaptive (error_tolerance=2.0)"
echo "  - Same initial conditions: cooling-0.gsd"
echo "  - Same output frequencies"
echo "  - Same molecular temperature tracking"
echo "  - Same controller time constant: 0.1 ps"
echo ""
echo "KEY DIFFERENCES:"
echo "  Exponential Wave:"
echo "    - Coupling: 3e-4 (exponential wave, period=50ps, tau=1.0ps)"
echo "    - Controllers: DiffEq (${DIFFEQ_TURN_ON_PS}-${DIFFEQ_TURN_OFF_PS}ps) → Adaptive Bath (${ADAPTIVE_TURN_ON_PS}ps+)"
echo "    - Final temperature reached: $FINAL_TEMP_K K"
echo ""
echo "  Zero Coupling:"
echo "    - Coupling: 0.0 (NO cavity)"
echo "    - Controller: GD feedback (kinetic temp, ON at ${GD_TURN_ON_PS}ps)"
echo "    - Target temperature: $FINAL_TEMP_K K (same as exp wave final)"
echo ""
echo "QUESTION: Can GD controller reach the same temperature without cavity coupling?"
echo "          Or does cavity provide unique heating/cooling effects?"
echo ""
echo "Output directories:"
echo "  - Exponential Wave: $EXP_OUTPUT_DIR"
echo "  - Zero Coupling: cavity_coupling_0eneg00/ (or similar)"
echo ""
echo "Key files to analyze:"
echo "  - Temperature trackers: temperature_tracker_replica_*.csv"
echo "  - Molecular temps: molecular_temperatures_replica_*.csv"
echo "  - Energy trackers: prod-*_energy_comprehensive.txt"
echo "  - Controller outputs:"
echo "      Exp wave: diffeq_controller_*.csv, adaptive_bath_controller_*.csv"
echo "      Zero coupling: gd_feedback_replica_*.csv"
echo ""
echo "Analysis questions:"
echo "  1. Does GD controller reach $FINAL_TEMP_K K as effectively as cavity+controllers?"
echo "  2. How do convergence rates compare?"
echo "  3. How does molecular temp decomposition differ (T_trans, T_rot, T_vib)?"
echo "  4. Are energy distributions statistically different?"
echo "  5. Does cavity coupling cause unique molecular-level effects?"

echo ""
echo "Benchmark complete! 🎯"
