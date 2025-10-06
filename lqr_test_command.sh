#!/bin/bash

# LQR Controller Test Script - MULTI-STEP SYSTEM ID VERSION
# Sequential control experiment: DiffEqController → LQR Controller → Control
# Phase 1 (10-110 ps):    DiffEqController tracks LJ+Coulombic, exponential wave coupling
# Phase 2 (110-210 ps):   LQR Multi-Step System ID
#                         - 4 sequential excitations (cool/heat/cool/return)
#                         - Each step is 25ps → total 100ps
#                         - Tests both directions + multiple magnitudes
#                         - Much better frequency coverage than single quench!
# Phase 3 (210+ ps):      LQR optimal control (signal-only, hot has inverse response)

python3 18_unified_cavity_dynamics.py \
  --molecular-bath bussi \
  --cavity-bath langevin \
  --coupling-type exponentialwave \
  --coupling 5e-4 \
  --exp-period 50.0 \
  --exp-tau 1.0 \
  --exp-start-time 10.0 \
  --exp-stop-time 110.0 \
  --exp-adaptive \
  --temperature 300.0 \
  --frequency 1560.0 \
  --molecular-tau 1.0 \
  --cavity-tau 1.0 \
  --runtime 10000.0 \
  --enable-energy-tracker \
  --enable-temp-tracker \
  --temp-tracker-empirical-data-file /home/mh7373/GitRepos/cav-hoomd/potential_energy_vs_T.txt \
  --enable-diffeq-controller \
  --diffeq-temperature-method lj_coulombic \
  --diffeq-time-constant 0.1 \
  --diffeq-turn-on-time 10.0 \
  --diffeq-turn-off-time 110.0 \
  --diffeq-update-interval 0.0 \
  --diffeq-apply-to both \
  --diffeq-T-min 2.0 \
  --enable-lqr-controller \
  --lqr-signal-method lj_coulombic \
  --lqr-hot-method harmonic \
  --lqr-target-temperature 300.0 \
  --lqr-dynamic-target \
  --lqr-dynamic-target-method lj_coulombic \
  --lqr-weight-signal 100.0 \
  --lqr-weight-hot 0.01 \
  --lqr-weight-bath 0.1 \
  --lqr-weight-integral 5.0 \
  --lqr-control-effort 10.0 \
  --lqr-process-noise-signal 2.0 \
  --lqr-process-noise-hot 10.0 \
  --lqr-measurement-noise-signal 2.0 \
  --lqr-measurement-noise-hot 30.0 \
  --lqr-system-id-mode multi_step \
  --lqr-system-id-temp 10.0 \
  --lqr-system-id-duration 100.0 \
  --lqr-periodic-system-id \
  --lqr-periodic-system-id-interval 100.0 \
  --lqr-turn-on-time 110.0 \
  --lqr-update-interval 1.0000 \
  --lqr-T-min 2.0 \
  --lqr-T-max 600.0 \
  --lqr-apply-to both \
  --lqr-empirical-data-file /home/mh7373/GitRepos/cav-hoomd/potential_energy_vs_T.txt \
  --device GPU \
  --seed 42 \
  --energy-output-period-ps 0.01 \
  --console-output-period-ps 1.0 \
  --error-tolerance 5.0 \
  --initial-fraction 1e-6 \
  --input-gsd cooling-0.gsd \
  --time-constant-ps 10.0 \
  --replicas 0 \
  --enable-dynamic-coupling-detection \
  --enable-molecular-temps \
  --molecular-temps-output-period-ps 1.0 \
  --coupling-change-threshold 1e-5

echo ""
echo "============================================"
echo "LQR Controller Test Complete!"
echo "============================================"
echo "Check output files:"
echo "  - lqr_control_replica_0.csv"
echo "  - lqr_system_params_replica_0.json"
echo "  - temperature_tracker_replica_0.csv"
echo "  - energy_tracker_replica_0.csv"
echo "============================================"

