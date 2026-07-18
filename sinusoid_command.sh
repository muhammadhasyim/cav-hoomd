#!/bin/bash

# Sinusoidal coupling experiment with dual independent temperature control
# Controls harmonic and LJ/Coulombic modes separately to prevent energy transfer bottleneck

python3 18_unified_cavity_dynamics.py \
  --molecular-bath bussi \
  --cavity-bath langevin \
  --coupling-type periodic \
  --coupling 1e-3 \
  --periodic-period 1.0 \
  --periodic-phase-offset 1.57 \
  --periodic-start-time 10.0 \
  --periodic-stop-time 1000.0 \
  --temperature 300.0 \
  --frequency 1560.0 \
  --molecular-tau 1.0 \
  --cavity-tau 1.0 \
  --runtime 10000.0 \
  --enable-energy-tracker \
  --enable-temp-tracker \
  --temp-tracker-empirical-data-file /home/mh7373/GitRepos/cav-hoomd/final_nodiss_cavitymd/potential_energy_components_vs_temperature.txt \
  --device GPU \
  --seed 42 \
  --replicas 0 \
  --energy-output-period-ps 0.01 \
  --console-output-period-ps 1.0 \
  --error-tolerance 5.0 \
  --initial-fraction 1e-6 \
  --input-gsd /home/mh7373/GitRepos/cav-hoomd/cooling-0.gsd \
  --time-constant-ps 100.0 \
  --enable-dynamic-coupling-detection \
  --coupling-change-threshold 1e-5 \
