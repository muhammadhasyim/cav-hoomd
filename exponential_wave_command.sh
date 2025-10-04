#!/bin/bash

# Exponential wave coupling experiment with controller switching at 550 ps
# Phase 1 (50-550 ps): DiffEqController tracks LJ+Coulombic temperature, exponential wave coupling
# Phase 2 (550+ ps): Adaptive bath controller with direct temperature adjustment
# Bath schedule: Direct proportional control based on signal temperature error
# If Ts(t) > Ttarget: Tbath = Ttarget - amplitude_scale * |Ts(t) - Ttarget| (cool down)
# If Ts(t) < Ttarget: Tbath = Ttarget + amplitude_scale * |Ts(t) - Ttarget| (heat up)
# Exponential wave: g(t) = A * exp(-t_period/tau) with period=25ps, tau=0.1ps (fast decay)
#python3 18_unified_cavity_dynamics.py \

python3 18_unified_cavity_dynamics.py \
  --molecular-bath bussi \
  --cavity-bath langevin \
  --coupling-type exponentialwave \
  --coupling 3e-4 \
  --exp-period 50.0 \
  --exp-tau 1.0 \
  --exp-start-time  500.0 \
  --exp-stop-time 1000.0 \
  --temperature 300.0 \
  --frequency 1560.0 \
  --molecular-tau 1.0 \
  --cavity-tau 1.0 \
  --runtime 3000.0 \
  --enable-energy-tracker \
  --enable-temp-tracker \
  --temp-tracker-empirical-data-file /home/mh7373/GitRepos/cav-hoomd/potential_energy_vs_T.txt \
  --enable-diffeq-controller \
  --diffeq-temperature-method lj_coulombic \
  --diffeq-time-constant 0.1 \
  --diffeq-turn-on-time 500.0 \
  --diffeq-turn-off-time 1000.0 \
  --diffeq-update-interval 0.0 \
  --diffeq-apply-to both \
  --diffeq-T-min 0.1 \
  --enable-adaptive-bath \
  --adaptive-bath-amplitude-scale 50.0 \
  --adaptive-bath-time-constant 0.1 \
  --adaptive-bath-turn-on-time 1000.0 \
  --adaptive-bath-update-interval 0.001 \
  --adaptive-bath-T-min 0.1 \
  --adaptive-bath-dynamic-target \
  --adaptive-bath-apply-to both \
  --adaptive-bath-empirical-data-file /home/mh7373/GitRepos/cav-hoomd/potential_energy_vs_T.txt \
  --adaptive-bath-signal-temperature-method lj_coulombic \
  --device GPU \
  --seed 42 \
  --energy-output-period-ps 0.01 \
  --console-output-period-ps 1.0 \
  --error-tolerance 2.0 \
  --initial-fraction 1e-6 \
  --input-gsd cooling-0.gsd \
  --time-constant-ps 10.0 \
  --replicas 325 \
  --enable-dynamic-coupling-detection \
  --enable-molecular-temps \
  --molecular-temps-output-period-ps 1.0 \
  --coupling-change-threshold 1e-5
