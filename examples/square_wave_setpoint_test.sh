#!/bin/bash

# Square Wave Coupling + Simple Setpoint Controller Test
# ====================================================================
# SEQUENTIAL CONTROL STRATEGY:
# - Phase 1 (Square wave active): Run square wave coupling with fixed bath temperature
# - Phase 2 (Square wave stops): SimpleSetpointController activates and captures lj_coulombic temperature as setpoint
# - Phase 3 (Setpoint control): Bath converges kinetic temperature to captured setpoint
#
# CONFIGURATION:
# - Square Wave Coupling: Periodic on/off pulses with fixed bath temperature
# - SimpleSetpointController: Activated when square wave stops, captures lj_coulombic as setpoint
# - Control Law: dT_bath/dt = -(T_kinetic - T_setpoint) / τ
#
# PHYSICAL PRINCIPLES:  
# - During square wave: System heats/cools under coupling, bath temperature held constant
# - At square wave stop: Capture final lj_coulombic temperature as equilibrium estimate
# - During setpoint control: Drive kinetic temperature to match that equilibrium estimate
# - Final state: System should stabilize at captured equilibrium temperature
# ====================================================================

# Lambda coupling parameter (dimensionless)
LAMBDA_COUPLING=${1:-0.025}  # Default lambda coupling

# Target temperatures
TARGET_TEMP=${2:-300.0}  # Default to 300.0 K (initial bath temperature)

# Square wave parameters
SQUARE_PERIOD=${3:-10.0}       # Period in ps (default: 10.0 ps)
SQUARE_DUTY_CYCLE=${4:-0.5}   # Duty cycle 0.0-1.0 (default: 0.5 = 50% on)
SQUARE_START_TIME=${5:-10.0}   # When square wave starts (default: 10.0 ps)
SQUARE_STOP_TIME=${6:-100.0}   # When square wave stops (default: 100.0 ps)

# SimpleSetpointController parameters
SETPOINT_TIME_CONSTANT=${7:-5.0}  # Control time constant in ps (default: 5.0 ps)
SETPOINT_APPLY_TO=${8:-molecular}  # Which thermostats to control (default: molecular)

# Simulation parameters
RUNTIME=${9:-500.0}  # Total simulation time in ps (default: 500.0 ps)

# Randomly choose a replica ID
REPLICA=22

echo ""
echo "========================================================================"
echo "Starting Square Wave + Simple Setpoint Controller simulation"
echo "Lambda coupling: ${LAMBDA_COUPLING} (dimensionless)"
echo "Initial temperature: ${TARGET_TEMP} K"
echo "Square wave parameters:"
echo "  Period: ${SQUARE_PERIOD} ps"
echo "  Duty cycle: ${SQUARE_DUTY_CYCLE}"
echo "  Start time: ${SQUARE_START_TIME} ps" 
echo "  Stop time: ${SQUARE_STOP_TIME} ps"
echo "SimpleSetpointController parameters:"
echo "  Signal method: lj_coulombic (for setpoint capture)"
echo "  Time constant: ${SETPOINT_TIME_CONSTANT} ps"
echo "  Apply to: ${SETPOINT_APPLY_TO}"
echo "  Turn on time: ${SQUARE_STOP_TIME} ps (when square wave stops)"
echo "Total runtime: ${RUNTIME} ps"
echo "Assigned replica: ${REPLICA}"
echo ""
echo "CONTROL SEQUENCE:"
echo "  t=0 to ${SQUARE_START_TIME} ps: Equilibration (no coupling, no control)"
echo "  t=${SQUARE_START_TIME} to ${SQUARE_STOP_TIME} ps: Square wave coupling (fixed bath temperature)" 
echo "  t=${SQUARE_STOP_TIME}+ ps: SimpleSetpointController (captures lj_coulombic → drives kinetic)"
echo "========================================================================"

python3 18_unified_cavity_dynamics.py \
  --lambda-coupling ${LAMBDA_COUPLING} \
  --coupling-type square \
  --square-period ${SQUARE_PERIOD} \
  --square-duty-cycle ${SQUARE_DUTY_CYCLE} \
  --square-phase-offset 0.0 \
  --square-start-time ${SQUARE_START_TIME} \
  --square-stop-time ${SQUARE_STOP_TIME} \
  --temperature ${TARGET_TEMP} \
  --frequency 1560.0 \
  --molecular-tau 1.0 \
  --runtime ${RUNTIME} \
  --device GPU \
  --seed 42 \
  --error-tolerance 1.0 \
  --initial-fraction 1e-3 \
  --time-constant-ps 50.0 \
  --enable-energy-tracker \
  --enable-temp-tracker \
  --enable-hdf5-output \
  --enable-molecular-temps \
  --energy-output-period 0.01 \
  --temp-tracker-output-period 0.005 \
  --molecular-temps-output-period-ps 0.01 \
  --hdf5-output-period 0.005 \
  --console-output-period 1.0 \
  --input-gsd /home/mh7373/GitRepos/cav-hoomd/cooling-0.gsd \
  --replicas ${REPLICA} \
  --temp-tracker-empirical-data-file /home/mh7373/GitRepos/cav-hoomd/potential_energy_vs_T.txt \
  --enable-simple-setpoint-controller \
  --simple-setpoint-signal-method lj_coulombic \
  --simple-setpoint-time-constant-ps ${SETPOINT_TIME_CONSTANT} \
  --simple-setpoint-apply-to ${SETPOINT_APPLY_TO} \
  --simple-setpoint-turn-on-time-ps ${SQUARE_STOP_TIME} \
  --simple-setpoint-update-interval-ps 0.1 \
  --simple-setpoint-T-min 0.0 \
  --simple-setpoint-output-file simple_setpoint_control.csv \
  --simple-setpoint-console-output-period-ps 1.0 \
  --enable-dynamic-coupling-detection \
  --coupling-change-threshold 1e-1

echo ""
echo "========================================================================"
echo "Simulation complete!"
echo "Output files:"
echo "  - replica_${REPLICA}/simple_setpoint_control_replica_${REPLICA}.csv"
echo "  - replica_${REPLICA}/comprehensive_temperatures_replica_${REPLICA}.csv"
echo "  - replica_${REPLICA}/energy_tracker_replica_${REPLICA}.csv"
echo "  - replica_${REPLICA}/molecular_temperatures_replica_${REPLICA}.csv"
echo ""
echo "EXPECTED BEHAVIOR:"
echo "  1. t=0-${SQUARE_START_TIME} ps: System equilibrates at ${TARGET_TEMP} K" 
echo "  2. t=${SQUARE_START_TIME}-${SQUARE_STOP_TIME} ps: Square wave coupling perturbs system"
echo "  3. t=${SQUARE_STOP_TIME}+ ps: SimpleSetpointController captures lj_coulombic temperature"
echo "  4. t=${SQUARE_STOP_TIME}+ ps: Bath drives kinetic temperature → captured setpoint"
echo "========================================================================"
