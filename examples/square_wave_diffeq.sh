#!/bin/bash

# Square Wave Coupling + Vanilla DiffEq Controller Test Script
# ====================================================================
# CONFIGURATION:
# - Square wave cavity coupling with period 10 ps, duty cycle 0.1 (on for 1 ps)
# - Vanilla differential equation controller (no bias correction, no PI control)
# - Simple first-order response: dT_bath/dt = -(T_bath - T_signal) / τ
# - Time constant: 1 ps
#
# CONTROL STRATEGY:
# - Direct tracking of LJ+Coulombic fictive temperature
# - No bias estimation, no PI control, no adaptive features
# - Pure differential equation feedback with fixed 1 ps time constant
# ====================================================================

# Square wave coupling parameters
LAMBDA_COUPLING=${1:-0.1}           # Maximum lambda value (dimensionless)
PERIOD_PS=${2:-10.0}                # Square wave period in ps (default: 10.0)
DUTY_CYCLE=${3:-0.1}                # Duty cycle (default: 0.1 = on for 1 ps out of 10 ps)

# DiffEq controller parameters
TIME_CONSTANT=${4:-1.0}             # Fixed time constant in ps (default: 1.0)
TARGET_TEMP=${5:-300.0}             # Initial/target temperature in K

# Simulation parameters
RUNTIME=${6:-200.0}                 # Total runtime in ps (default: 200.0 = 20 periods)
REPLICA=${7:-10}                    # Replica ID

# Calculate pulse width from duty cycle
PULSE_WIDTH=$(echo "$PERIOD_PS * $DUTY_CYCLE" | bc -l)

echo ""
echo "========================================================================"
echo "Starting Square Wave + Vanilla DiffEq Controller simulation"
echo "========================================================================"
echo "COUPLING CONFIGURATION:"
echo "  Type: Square wave (periodic)"
echo "  Maximum lambda: ${LAMBDA_COUPLING} (dimensionless)"
echo "  Period: ${PERIOD_PS} ps"
echo "  Duty cycle: ${DUTY_CYCLE} (${PULSE_WIDTH} ps on, $(echo "$PERIOD_PS - $PULSE_WIDTH" | bc -l) ps off)"
echo ""
echo "CONTROLLER CONFIGURATION:"
echo "  Type: Vanilla DiffEq (first-order)"
echo "  Signal: LJ+Coulombic fictive temperature"
echo "  Time constant: ${TIME_CONSTANT} ps"
echo "  Bias correction: DISABLED"
echo "  PI control: DISABLED"
echo "  Adaptive features: DISABLED"
echo ""
echo "SIMULATION SETTINGS:"
echo "  Initial temperature: ${TARGET_TEMP} K"
echo "  Runtime: ${RUNTIME} ps ($(echo "$RUNTIME / $PERIOD_PS" | bc -l | xargs printf "%.1f") periods)"
echo "  Replica: ${REPLICA}"
echo "========================================================================"
echo ""

# Build the command
CMD="python3 18_unified_cavity_dynamics.py \
  --lambda-coupling ${LAMBDA_COUPLING} \
  --coupling-type square \
  --square-period ${PERIOD_PS} \
  --square-duty-cycle ${DUTY_CYCLE} \
  --temperature ${TARGET_TEMP} \
  --frequency 1560.0 \
  --molecular-bath bussi \
  --cavity-bath bussi \
  --molecular-tau 1.0 \
  --cavity-tau 1.0 \
  --runtime ${RUNTIME} \
  --device GPU \
  --seed 42 \
  --error-tolerance 5.0 \
  --initial-fraction 1e-3 \
  --time-constant-ps 50.0 \
  --enable-energy-tracker \
  --enable-temp-tracker \
  --enable-molecular-temps \
  --enable-hdf5-output \
  --energy-output-period-ps 0.01 \
  --temp-tracker-output-period-ps 0.01 \
  --hdf5-output-period-ps 0.01 \
  --console-output-period-ps 1.0 \
  --input-gsd /home/mh7373/GitRepos/cav-hoomd/examples/cooling-0.gsd \
  --replicas ${REPLICA} \
  --temp-tracker-empirical-data-file /home/mh7373/GitRepos/cav-hoomd/potential_energy_vs_T.txt \
  --empirical-data-file /home/mh7373/GitRepos/cav-hoomd/potential_energy_vs_T.txt \
  --enable-diffeq-controller \
  --diffeq-temperature-method lj_coulombic \
  --diffeq-time-constant ${TIME_CONSTANT} \
  --diffeq-turn-on-time 0.0 \
  --diffeq-update-interval 0.01 \
  --diffeq-apply-to both \
  --diffeq-T-min 0.1 \
  --diffeq-disable-bias-estimation"

echo "Executing command:"
echo "$CMD"
echo ""

# Execute
eval $CMD

EXIT_CODE=$?

echo ""
echo "========================================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "Square Wave + DiffEq Controller simulation COMPLETED SUCCESSFULLY"
    echo ""
    echo "📁 Output files (in replica_${REPLICA}/):"
    echo "  - GSD trajectory: prod-${REPLICA}.gsd"
    echo "  - DiffEq controller CSV: diffeq_control_replica_${REPLICA}.csv"
    echo "  - HDF5 observables: observables_replica_${REPLICA}.h5"
    echo "  - Temperature tracker CSV: temperature_tracker_replica_${REPLICA}.csv"
    echo "  - Energy tracker CSV: prod-${REPLICA}_energy_comprehensive.txt"
    echo ""
    echo "📊 Analysis suggestions:"
    echo "  - Plot T_signal vs T_bath to see controller response"
    echo "  - Overlay square wave coupling on temperature plot"
    echo "  - Check if controller maintains stability during pulses"
    echo "  - Measure response time to coupling changes"
else
    echo "Square Wave + DiffEq Controller simulation FAILED with exit code ${EXIT_CODE}"
fi
echo "========================================================================"
echo ""

exit $EXIT_CODE

