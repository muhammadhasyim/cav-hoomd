#!/bin/bash

# LQG Controller Test - TIGHT SIGNAL+KINETIC TRACKING 🎯
# ====================================================================
# ✨ STRATEGY: Ignore T_hot, focus on tight T_signal and T_kinetic!
#
# INNOVATION: Strong cross-coupling + high individual weights force
# T_signal and T_kinetic to stay LOCKED TOGETHER at target
#
# Cost Function:
#   J = q_s·x_s² + q_kin·x_kin² + q_int·ξ² + w_sk·(x_s - x_kin)² + R·u²
#   (T_hot is IGNORED: q_h ≈ 0, no hot cross-coupling)
#
# CONFIGURATION (GENTLER CONTROL - ANTI-OSCILLATION):
# - State: [x_s, x_h, x_kin, ξ] (4D augmented, but T_h has minimal weight)
# - Weights: q_s=10.0, q_h=0.01 (ignored!), q_kin=10.0, q_int=1.0 (reduced 10x)
# - Cross-coupling: w_sk=5.0 (moderate signal↔kinetic coupling)
# - Control effort: R=10.0 (gentler than before)
# - Bath limits: [50K, 500K] (wider authority to prevent saturation)
# - System ID: 200ps duration, 50K excitation (safer than 0K!)
# - Gain scheduling: DISABLED
# - No periodic re-ID (prevents 0K crashes)
#
# EXPECTED: Smooth convergence T_signal ≈ T_kinetic ≈ T_target, T_hot free-floating!
# ====================================================================

python3 18_unified_cavity_dynamics.py \
  --molecular-bath bussi \
  --cavity-bath langevin \
  --coupling-type exponentialwave \
  --coupling 1e-4 \
  --exp-period 1.0 \
  --exp-tau 1.0 \
  --exp-start-time 10.0 \
  --exp-stop-time 210.0 \
  --temperature 300.0 \
  --frequency 1560.0 \
  --molecular-tau 1.0 \
  --cavity-tau 1.0 \
  --runtime 50000.0 \
  --enable-energy-tracker \
  --enable-temp-tracker \
  --temp-tracker-empirical-data-file /home/mh7373/GitRepos/cav-hoomd/potential_energy_vs_T.txt \
  --enable-diffeq-controller \
  --diffeq-temperature-method lj_coulombic \
  --diffeq-time-constant 0.1 \
  --diffeq-turn-on-time 10.0 \
  --diffeq-turn-off-time 210.0 \
  --diffeq-update-interval 0.01 \
  --diffeq-apply-to both \
  --diffeq-T-min 2.0 \
  --enable-lqr-controller \
  --lqr-signal-method lj_coulombic \
  --lqr-hot-method harmonic \
  --lqr-target-temperature 300.0 \
  --lqr-dynamic-target \
  --lqr-dynamic-target-method lj_coulombic \
  --lqr-weight-signal 10.0 \
  --lqr-weight-hot 0.01 \
  --lqr-weight-kinetic 10.0 \
  --lqr-weight-integral 1.0 \
  --lqr-control-effort 10.0 \
  --lqr-cross-coupling-signal-kinetic 5.0 \
  --lqr-cross-coupling-signal-hot 0.01 \
  --lqr-cross-coupling-hot-kinetic 0.01 \
  --lqr-process-noise-signal 2.0 \
  --lqr-process-noise-hot 20.0 \
  --lqr-process-noise-kinetic 2.0 \
  --lqr-measurement-noise-signal 2.0 \
  --lqr-measurement-noise-hot 30.0 \
  --lqr-measurement-noise-kinetic 2.0 \
  --lqr-system-id-mode multi_step \
  --lqr-system-id-temp 10.0 \
  --lqr-system-id-duration 200.0 \
  --lqr-track-kinetic-temp \
  --lqr-turn-on-time 210.0 \
  --lqr-update-interval 1.0 \
  --lqr-T-min 50.0 \
  --lqr-T-max 500.0 \
  --lqr-apply-to both \
  --lqr-th-filter-enabled \
  --lqr-th-filter-time-constant 2.0 \
  --lqr-disable-gain-scheduling \
  --lqr-adaptive-lqr-threshold 1000.0 \
  --lqr-empirical-data-file /home/mh7373/GitRepos/cav-hoomd/potential_energy_vs_T.txt \
  --device GPU \
  --seed 42 \
  --energy-output-period-ps 0.01 \
  --console-output-period-ps 1.0 \
  --error-tolerance 2.0 \
  --initial-fraction 1e-6 \
  --input-gsd cooling-0.gsd \
  --time-constant-ps 2.0 \
  --replicas 297 \
  --enable-dynamic-coupling-detection \
  --enable-molecular-temps \
  --molecular-temps-output-period-ps 1.0 \
  --coupling-change-threshold 1e-4

echo ""
echo "============================================"
echo "TIGHT SIGNAL+KINETIC TRACKING! 🎯🔥"
echo "============================================"
echo "Key configuration (ANTI-OSCILLATION TUNING):"
echo "  1. 3D base system: [x_s, x_h, x_kin] (T_h ignored!)"
echo "  2. Individual weights: q=[signal:10.0, hot:0.01, kin:10.0, int:1.0] (10x GENTLER)"
echo "  3. Cross-coupling: w_sk=5.0 (moderate, not aggressive)"
echo "  4. Control effort: R=10.0 (gentler control)"
echo "  5. Bath limits: [50K, 500K] (wide authority, prevents saturation)"
echo "  6. System ID: 200ps @ 50K (longer, safer than 0K)"
echo "  7. T_h low-pass filter: τ=2.0ps"
echo "  8. NO periodic re-identification (prevents crashes)"
echo "  9. Gain scheduling: DISABLED"
echo ""
echo "Expected behavior:"
echo "  - SMOOTH convergence to T_target=300K ✅"
echo "  - T_signal ≈ T_kinetic ≈ 300K (within ±5K)"
echo "  - No wild oscillations or integral windup"
echo "  - T_hot FREE-FLOATING (not controlled)"
echo "  - No particle escapes!"
echo ""
echo "Focus: STABLE, gentle temperature control! 🎯"
echo "============================================"

