#!/usr/bin/env python3
"""
Universal Cavity MD Experiment Runner - All Phase 3 Features

This script provides a comprehensive framework for running cavity MD simulations
with ALL enhanced coupling modes and control features:

 ENHANCED COUPLING VARIANTS:
- CONSTANT: g(t) = constant (default)
- STEP: g(t) = 0 → g_target at t_switch (with optional decay and turn-off)
- PERIODIC: g(t) = A * sin(2π*t/T + φ) (with start/stop times)
- EXPONENTIAL: g(t) = A * exp(-t/τ) (with turn-on/off times)
- SQUARE WAVE: g(t) = periodic on/off pulses (with duty cycle)

 ENHANCED LASER DRIVE:
- F_laser = F_L * cos(ω_L * t) with precise timing control
- Start/stop times for laser pulses
- Arbitrary k-vector directions

 PI FEEDBACK CONTROLLER:
- Kinetic temperature control (from particle velocities)
- LJ+Coulombic fictive temperature (Rosenfeld scaling)
- Harmonic fictive temperature (T^3/2 scaling)
- IMC auto-tuning or manual PI parameters

 HARMONIC BOND RESET:
- One-time thermal reinitialization of bond stretch DOF for ALL bond types
- Auto-detects bond parameters (K, r0) from simulation force field
- Preserves COM motion and molecular orientation
- Option B: Exact canonical sampling at target temperature
- Triggered at specified time during simulation
- Dynamic temperature: automatically uses bath temperature at reset time

 PRODUCTION FEATURES:
- All Phase 3 enhanced features integrated
- Full GPU/CPU compatibility and consistency 
- Backward compatibility with existing workflows
- Smart folder naming with all parameters
- Comprehensive energy tracking and analysis
- Multi-replica execution support
- Advanced timing control for all features
- Real-time PI feedback with multiple temperature methods

 ENHANCED USAGE EXAMPLES:

  # Exponential decay coupling with PI feedback
  # Exponential wave coupling with adaptive amplitude
  python 18_unified_cavity_dynamics.py --coupling-type exponentialwave --coupling 2e-4 --exp-period 50.0 --exp-tau 0.1 --exp-adaptive --enable-temp-tracker
  
  # Square wave coupling with timed laser
  python 18_unified_cavity_dynamics.py --coupling-type square --coupling 1e-3 --square-period 8.0 --square-duty-cycle 0.3 --laser --laser-start-time 10.0 --laser-stop-time 30.0
  
  # Decaying square wave coupling (10% amplitude decay per period)
  python 18_unified_cavity_dynamics.py --coupling-type decaying_square --coupling 1e-3 --decaying-square-period 5.0 --decaying-square-decay-rate 0.10 --decaying-square-duty-cycle 0.75 --runtime 100
  
  # Adaptive square wave coupling with quench controller (amplitude adjusts based on harmonic temperature)
  python 18_unified_cavity_dynamics.py --coupling-type adaptive_square --coupling 1e-3 --adaptive-square-period 8.0 --adaptive-square-duty-cycle 0.75 --enable-quench-controller --quench-target-temperature 10.0 --enable-temp-tracker
  
  # Auto-stop coupling when temperatures converge (stops coupling when |T_fictive - T_kinetic| < 1K for 10ps)
  python 18_unified_cavity_dynamics.py --coupling-type decaying_square --coupling 1e-3 --enable-auto-stop --auto-stop-tol 1.0 --auto-stop-window 10.0 --enable-temp-tracker
  
  # Enhanced periodic coupling with full parameter control
  # Step coupling with LJ+Coulombic feedback
  # Quench experiment: rapid cooling from 300K to 100K at 50ps
  python 18_unified_cavity_dynamics.py --coupling 1e-3 --enable-quench-controller --quench-initial-temperature 300.0 --quench-target-temperature 100.0 --quench-time 50.0 --runtime 200
  
  # Multi-feature advanced experiment
  # Legacy examples (still supported):
  # Periodic coupling with mechanical modulation
  python 18_unified_cavity_dynamics.py --periodic --coupling 1e-3 --period 1.0 --mech-periodic --mech-frequency 100.0 --mech-magnitude 1e-4
  
  # Step coupling (legacy switch-time)
  python 18_unified_cavity_dynamics.py --coupling 1e-3 --switch-time 10.0 --enable-empirical-feedback --empirical-data-file data.txt

Created: 2025-09-25 (Phase 1)
Enhanced: 2025-09-26 (Phase 3 - Full Integration)
Author: AI Assistant + User Collaboration
Status: Production Ready - All Phase 3 Features Integrated
"""

import sys
import os
import argparse
import logging
from pathlib import Path
from typing import Optional, List, Tuple, Dict, Any
import numpy as np

# Import HOOMD and plugin modules
import hoomd
from hoomd.cavitymd.simulation import CavityMDSimulation
from hoomd.cavitymd.utils import PhysicalConstants
from hoomd.cavitymd.updaters import HarmonicBondReset

# =============================================================================
# SIMULATION FUNCTIONS (Unified: Periodic, Laser, Mechanical, Switch-time)
# =============================================================================

def run_single_experiment(molecular_thermo, cavity_thermo, finite_q,
                         coupling=None, lambda_coupling=1e-3, temperature=100.0, frequency=1560.0,
                         runtime_ps=500.0, input_gsd='molecular-0.gsd', frame=-1,
                         # Enhanced coupling variant parameters
                         coupling_variant_type='constant',
                         # Step coupling parameters (enhanced)
                         switch_time_ps=None, decay_time_constant_ps=None, step_turn_off_time_ps=None,
                         # Periodic coupling parameters (enhanced)
                         periodic_coupling=False, periodic_period_ps=1.0, periodic_phase_offset=0.0, 
                         periodic_start_time_ps=0.0, periodic_stop_time_ps=None,
                         # Exponential decay parameters
                         exponential_amplitude=None, exponential_decay_time_ps=10.0,
                         exponential_turn_on_time_ps=0.0, exponential_turn_off_time_ps=None,
                         # Square wave parameters
                         square_period_ps=2.0, square_duty_cycle=0.5, square_phase_offset=0.0,
                         square_start_time_ps=0.0, square_stop_time_ps=None,
                         # Decaying square wave parameters
                         decaying_square_period_ps=2.0, decaying_square_duty_cycle=0.5, decaying_square_phase_offset=0.0,
                         decaying_square_start_time_ps=0.0, decaying_square_stop_time_ps=None,
                         decaying_square_decay_rate=0.1, decaying_square_minimum_amplitude=1e-6,
                         # Adaptive square wave parameters
                         adaptive_square_period_ps=5.0, adaptive_square_duty_cycle=0.5, adaptive_square_phase_offset=0.0,
                         adaptive_square_start_time_ps=0.0, adaptive_square_stop_time_ps=None,
                         adaptive_square_min_amplitude=1e-8, adaptive_square_max_amplitude=1e-1,
                         # Exponential wave parameters
                         exp_period_ps=2.0, exp_tau_ps=0.5, exp_start_time_ps=0.0, exp_stop_time_ps=None, exp_adaptive=False,
                         # Composite coupling parameters
                         composite_sinusoid_amplitude=1e-4, composite_sinusoid_period=1.0, composite_sinusoid_phase=0.0,
                         composite_sinusoid_start_time=0.0, composite_sinusoid_stop_time=None,
                         composite_square_amplitude=2e-4, composite_square_period=50.0, composite_square_duty_cycle=0.02,
                         composite_square_start_time=10.0, composite_square_stop_time=1000.0, 
                         composite_square_adaptive=False, composite_max_amplitude=None,
                         # Auto-stop coupling parameters
                         enable_auto_stop=False, auto_stop_tol=1.0, auto_stop_window=10.0,
                         # Enhanced laser drive parameters  
                         laser_enabled=False, laser_frequency_cm1=1560.0, laser_amplitude=1e-5, 
                         laser_start_time_ps=0.0, laser_stop_time_ps=None, laser_kvector=None,
                         # Gradient descent feedback controller parameters
                         enable_gd_feedback=False, gd_target_temperature=100.0,
                         gd_turn_on_time_ps=0.0, gd_turn_off_time_ps=None,
                         gd_temperature_method='kinetic', gd_update_interval_ps=0.1,
                         gd_time_constant_ps=10.0, gd_apply_to='both', gd_T_min=0.0, gd_T_max=None,
                         gd_disable_effective_temp=False,
                         # Multi-signal error function parameters
                         gd_enable_multi_signal=False, gd_weight_system_target=1.0,
                         gd_weight_bath_target=0.0, gd_weight_system_bath=0.0,
                         # Dual independent feedback controller parameters
                         enable_dual_feedback=False,
                         dual_cavity_method='harmonic_equipartition',
                         dual_molecular_method='lj_coulombic_kinetic',
                         dual_cavity_target_temperature=100.0,
                         dual_molecular_target_temperature=100.0,
                         dual_cavity_time_constant_ps=5.0,
                         dual_molecular_time_constant_ps=10.0,
                         dual_turn_on_time_ps=0.0,
                         dual_turn_off_time_ps=None,
                         dual_update_interval_ps=0.1,
                         dual_cavity_T_min=0.0,
                         dual_cavity_T_max=None,
                         dual_molecular_T_min=0.0,
                         dual_molecular_T_max=None,
                         dual_cavity_dynamic_target=False,
                         dual_molecular_dynamic_target=False,
                         dual_cavity_integral_time_constant_ps=None,
                         dual_molecular_integral_time_constant_ps=None,
                         # Sinusoidal bath temperature controller parameters
                         enable_sinusoidal_bath=False,
                         sinusoidal_bath_period_ps=1.0,
                         sinusoidal_bath_amplitude_scale=0.1,
                         sinusoidal_bath_phase_offset=0.0,
                         sinusoidal_bath_target_temperature=100.0,
                         sinusoidal_bath_dynamic_target=False,
                         sinusoidal_bath_turn_on_time_ps=0.0,
                         sinusoidal_bath_turn_off_time_ps=None,
                         sinusoidal_bath_update_interval_ps=0.1,
                         sinusoidal_bath_apply_to='both',
                         sinusoidal_bath_T_min=0.1,
                         sinusoidal_bath_T_max=None,
                         sinusoidal_bath_empirical_data_file=None,
                         sinusoidal_bath_amplitude_update_interval_ps=1.0,
                         sinusoidal_bath_amplitude_temperature_method='harmonic_equipartition',
                         sinusoidal_bath_adaptive_range_mode=False,
                         # Adaptive bath temperature controller parameters
                         enable_adaptive_bath=False,
                         adaptive_bath_amplitude_scale=1.0,
                         adaptive_bath_time_constant_ps=1.0,
                         adaptive_bath_target_temperature=100.0,
                         adaptive_bath_dynamic_target=True,
                         adaptive_bath_turn_on_time_ps=0.0,
                         adaptive_bath_turn_off_time_ps=None,
                         adaptive_bath_update_interval_ps=0.1,
                         adaptive_bath_apply_to='both',
                         adaptive_bath_T_min=0.1,
                        adaptive_bath_T_max=None,
                        adaptive_bath_empirical_data_file=None,
                        adaptive_bath_signal_temperature_method='harmonic_equipartition',
                        adaptive_bath_dynamic_target_temperature_method=None,
                        # Quench controller parameters
                         enable_quench_controller=False, quench_initial_temperature=100.0,
                         quench_target_temperature=50.0, quench_time_ps=50.0, quench_apply_to='both',
                         # Offset temperature controller parameters
                         enable_offset_controller=False, offset_temperature_method='kinetic',
                         offset_temperature_offset_K=-50.0, offset_turn_on_time_ps=0.0,
                         offset_turn_off_time_ps=None, offset_update_interval_ps=0.1,
                         offset_apply_to='both', offset_T_min=0.0, offset_T_max=None,
                         # Harmonic bond reset parameters
                         enable_harmonic_reset=False, harmonic_reset_turn_on_time_ps=0.0,
                         harmonic_reset_temperature=None, harmonic_reset_seed=42,
                         # Differential equation controller parameters
                         enable_diffeq_controller=False, diffeq_temperature_method='kinetic',
                         diffeq_time_constant_ps=5.0, diffeq_time_constant_auto=False,
                         diffeq_turn_on_time_ps=0.0, diffeq_turn_off_time_ps=None, diffeq_update_interval_ps=0.1,
                        diffeq_apply_to='both', diffeq_T_min=0.0, diffeq_T_max=None,
                        diffeq_rate_limit_K_per_ps=None, diffeq_disable_bias_estimation=False,
                        # Exact-cancellation + PI control parameters for DiffEq controller
                        diffeq_enable_pi_control=False, diffeq_pi_rho=1.0, diffeq_pi_epsilon=0.5,
                        diffeq_pi_zeta=0.8, diffeq_relaxation_data_file=None,
                        diffeq_filter_window=0.0,
                        # Adaptive bias cancellation parameters for DiffEq controller
                        diffeq_enable_bias_cancellation=False, diffeq_bias_tau_b_ps=50.0,
                        diffeq_bias_tau_b_auto=False, diffeq_bias_kappa=0.01,
                       diffeq_bias_kappa_auto=False, diffeq_bias_tau_b_prefactor=5.0,
                       diffeq_bias_kappa_prefactor=50.0, diffeq_bias_calibration_time_ps=10.0,
                       # Simple Setpoint Controller parameters
                       enable_simple_setpoint_controller=False, simple_setpoint_signal_method='kinetic',
                       simple_setpoint_time_constant_ps=5.0, simple_setpoint_apply_to='both',
                       simple_setpoint_turn_on_time_ps=0.0, simple_setpoint_turn_off_time_ps=None,
                       simple_setpoint_update_interval_ps=0.1, simple_setpoint_T_min=0.0,
                       simple_setpoint_T_max=None, simple_setpoint_output_file='simple_setpoint_control.csv',
                       simple_setpoint_console_output_period_ps=1.0,
                       # Adaptive MPC controller parameters
                       enable_adaptive_mpc_controller=False, adaptive_mpc_target_temperature=100.0,
                       adaptive_mpc_turn_on_time_ps=0.0, adaptive_mpc_turn_off_time_ps=None,
                       adaptive_mpc_system_id_duration_ps=50.0, adaptive_mpc_system_id_step_duration_ps=5.0,
                       adaptive_mpc_system_id_seed=42, adaptive_mpc_update_interval_ps=0.1,
                       adaptive_mpc_prediction_horizon=10, adaptive_mpc_control_horizon=5,
                       adaptive_mpc_output_weight=100.0, adaptive_mpc_control_effort_weight=None,
                       adaptive_mpc_rate_penalty_weight=None, adaptive_mpc_lambda_min=0.0,
                       adaptive_mpc_lambda_max=1e-2, adaptive_mpc_T_bath_min=0.1,
                       adaptive_mpc_T_bath_max=500.0, adaptive_mpc_delta_lambda_max=1e-4,
                       adaptive_mpc_delta_T_bath_max=10.0, adaptive_mpc_apply_to='both',
                       adaptive_mpc_rls_forgetting_factor=0.995, adaptive_mpc_rls_initial_covariance=100.0,
                       adaptive_mpc_model_update_interval=10, adaptive_mpc_output_file='adaptive_mpc_control.csv',
                       adaptive_mpc_console_output_period_ps=1.0,
                       adaptive_mpc_regularization_param=1e-3,
                       adaptive_mpc_use_scaling=True,
                       adaptive_mpc_debug_mode=False,
                       # PID controller parameters
                       enable_pid_controller=False, pid_signal_choice='lj_coulombic',
                       pid_target_temperature=100.0, pid_self_loop=False,
                       pid_Kp=None, pid_Ti=None, pid_Td=None,
                       pid_auto_tune=True, pid_auto_tune_step_size=20.0,
                       pid_auto_tune_duration_ps=50.0, pid_turn_on_time_ps=0.0,
                       pid_turn_off_time_ps=None, pid_update_interval_ps=0.1,
                       pid_apply_to='both', pid_T_min=0.1, pid_T_max=None,
                       pid_rate_limit_K_per_ps=None, pid_derivative_filter_N=10.0,
                       pid_enable_anti_windup=True, pid_output_file='pid_control.csv',
                       pid_console_output_period_ps=1.0,
                       # BathPI controller parameters
                       enable_bath_pi_controller=False, bath_pi_apply_to='both',
                       bath_pi_K_p_molecular=0.1, bath_pi_K_i_molecular=0.01, bath_pi_K_T_molecular=0.0,
                       bath_pi_K_p_cavity=0.1, bath_pi_K_i_cavity=0.01, bath_pi_K_T_cavity=0.0,
                       bath_pi_filter_window_ps=5.0, bath_pi_flux_source='reservoir',
                       bath_pi_anti_windup_alpha=0.01,
                       bath_pi_enable_feedforward=False, bath_pi_T_nominal=None,
                       bath_pi_feedforward_tau_ps=1000.0, bath_pi_turn_on_time_ps=0.0,
                       bath_pi_turn_off_time_ps=None, bath_pi_update_interval_ps=0.1,
                       bath_pi_T_min=0.1, bath_pi_T_max=None, bath_pi_rate_limit_K_per_ps=None,
                       bath_pi_output_file='bath_pi_control.csv', bath_pi_relaxation_data_file=None,
                       # LQR controller parameters
                       enable_lqr_controller=False, lqr_signal_method='lj_coulombic',
                        lqr_hot_method='harmonic_equipartition', lqr_target_temperature=300.0,
                        lqr_dynamic_target=False, lqr_dynamic_target_method=None,
                        lqr_weight_signal=100.0, lqr_weight_hot=1.0, lqr_weight_bath=0.1,
                        lqr_weight_integral=10.0, lqr_control_effort=1.0,
                        lqr_process_noise_signal=0.1, lqr_process_noise_hot=0.1,
                        lqr_measurement_noise_signal=0.5, lqr_measurement_noise_hot=0.5,
                        lqr_system_id_mode='step', lqr_system_id_temp_K=5.0,
                        lqr_system_id_duration_ps=50.0, lqr_system_id_file='lqr_system_params.json',
                        lqr_periodic_system_id=False, lqr_periodic_system_id_interval_ps=1000.0,
                        # EKF-based adaptation
                        lqr_use_ekf_adaptation=True, lqr_ekf_update_interval=50,
                        lqr_ekf_process_noise_param=0.001, lqr_ekf_initial_covariance_param=0.1,
                        lqr_adaptive_lqr_threshold=0.05,
                        # Gain scheduling
                        lqr_enable_gain_scheduling=True, lqr_gain_schedule_far_threshold=20.0,
                        lqr_gain_schedule_near_threshold=10.0,
                        # T_h low-pass filter
                        lqr_th_filter_enabled=True, lqr_th_filter_time_constant=20.0,
                        # Gentle startup
                        lqr_gentle_startup_steps=10, lqr_gentle_startup_min_authority=0.1,
                        # Kinetic temperature tracking (3D state)
                        lqr_track_kinetic_temp=False, lqr_weight_kinetic=100.0,
                        lqr_process_noise_kinetic=2.0, lqr_measurement_noise_kinetic=2.0,
                        # Cross-coupling weights for thermal equilibration
                        lqr_cross_coupling_signal_kinetic=0.0,
                        lqr_cross_coupling_signal_hot=0.0,
                        lqr_cross_coupling_hot_kinetic=0.0,
                        # Timing and limits
                        lqr_turn_on_time_ps=0.0, lqr_turn_off_time_ps=None,
                        lqr_update_interval_ps=0.1, lqr_T_min=0.1, lqr_T_max=None,
                        lqr_apply_to='both', lqr_output_file='lqr_controller.csv',
                        lqr_empirical_data_file=None,
                        # Adaptive LQI controller parameters
                        lqr_controller_type='standard', lqr_tau_L_initial=200.0, lqr_tau_H_initial=30.0,
                        lqr_k_initial=0.01, lqr_tau_b=1.0, lqr_q_common=1000.0, lqr_q_diff=100.0,
                        lqr_q_eta_common=1e6, lqr_q_eta_diff=10.0, lqr_process_noise_drift=0.01,
                        lqr_rls_forgetting=0.998, lqr_rls_regularization=1e-6, lqr_rls_update_interval=50,
                        lqr_max_control_rate=10.0, lqr_integral_max_common=1000.0, lqr_integral_max_diff=100.0,
                        lqr_theta_change_threshold=0.05,
                        # LQG pulse controller parameters
                        enable_lqg_controller=False,
                        lqg_A_matrix=None, lqg_B_matrix=None, lqg_C_matrix=None,
                        lqg_S_matrix=None, lqg_r_target=None,
                        lqg_Qz_weights=None, lqg_Ru_weights=None, lqg_Qe_weights=None,
                        lqg_W_noise=None, lqg_V_noise=None,
                        lqg_g0_baseline=1e-4, lqg_dt=5.0,
                        lqg_delta_g_min=-0.3, lqg_delta_g_max=0.3,
                        lqg_T_bath_min=10.0, lqg_T_bath_max=500.0,
                        lqg_rate_limit_K_per_ps=1.0,
                        lqg_turn_on_time_ps=0.0, lqg_turn_off_time_ps=None,
                        lqg_update_interval_ps=5.0, lqg_console_output_period_ps=10.0,
                        lqg_c_affine_term=None, lqg_temperature_methods=None,
                        # LQG square wave parameters
                        lqg_square_period_ps=5.0, lqg_square_duty_cycle=0.5,
                        lqg_square_phase_offset=0.0, lqg_square_start_time_ps=0.0,
                        lqg_square_stop_time_ps=None, lqg_square_g_low_level=1e-6,
                        # LQG system identification parameters
                        lqg_enable_system_id=True, lqg_equilibrium_duration_ps=500.0,
                        lqg_prbs_n_unique_steps=10, lqg_prbs_amplitude=5e-5,
                        lqg_bath_n_unique_steps=8, lqg_bath_excitation_amplitude=10.0,
                        lqg_system_id_output_file='lqg_system_id.json',
                        # LQG coupling controller parameters (new coupling-type: lqg_coupling)
                        lqg_coupling_target_temperature=100.0, lqg_coupling_update_interval_ps=0.1,
                        lqg_coupling_equilibrium_duration_ps=5.0, lqg_coupling_step_duration_ps=5.0,
                        lqg_coupling_n_steps=2, lqg_coupling_lambda_min=0.0, lqg_coupling_lambda_max=1e-2,
                        lqg_coupling_process_noise_std=0.1, lqg_coupling_measurement_noise_std=0.5,
                        lqg_coupling_weight_signal=100.0, lqg_coupling_weight_harmonic=1.0,
                        lqg_coupling_weight_kinetic=1.0, lqg_coupling_weight_bath=0.1,
                        lqg_coupling_control_effort=1.0, lqg_coupling_integral_gain=2.0,
                        lqg_coupling_system_id_file='lqg_coupling_system_id.json',
                        lqg_coupling_control_file='lqg_coupling_control.csv',
                        lqg_coupling_temperature_methods=['lj_coulombic', 'harmonic_equipartition', 'kinetic'],
                        lqg_coupling_bath_temperature_method='kinetic',
                        # Temperature tracker parameters
                         enable_temp_tracker=False, temp_tracker_output_period_ps=0.1,
                         temp_tracker_empirical_data_file=None,
                         # Molecular temperature decomposition parameters
                         enable_molecular_temps=False, molecular_temps_output_period_ps=1.0,
                         # Dipole moment FDR parameters
                         enable_dipole_fdr=False, dipole_fdr_output_period_ps=0.1,
                         dipole_fdr_max_correlation_time_ps=100.0, dipole_fdr_field_direction=[0, 0, 1],
                         dipole_fdr_exclude_cavity=True, enable_dipole_response=False,
                         dipole_response_field_strength=1e-5, dipole_response_sign=1.0,
                         # Legacy parameters (backward compatibility)
                         damping_ratio=0.0,
                         # Mechanical cavity modulation parameters
                         mech_periodic=False, mech_frequency_cm1=100.0, mech_magnitude=1e-4,
                         # Thermostat parameters
                         molecular_tau=5.0, cavity_tau=5.0,
                         # Timestep parameters
                         fixed_timestep=False, timestep_fs=0.1, error_tolerance=0.01,
                         initial_fraction=1e-5, time_constant_ps=50.0,
                         # Dynamic coupling detection parameters
                         enable_dynamic_coupling_detection=True, coupling_change_threshold=1e-5,
                         # Output and tracking parameters
                         enable_energy_tracker=False, energy_output_period_ps=0.1,
                         fkt_output_period_ps=1.0, gsd_output_period_ps=50.0,
                         console_output_period_ps=1.0, enable_fkt=False,
                         fkt_kmag=1.0, fkt_wavevectors=50, fkt_ref_interval=1.0, fkt_max_refs=10,
                         enable_dipole_autocorr=False, dipole_ref_interval=1.0, dipole_max_refs=10,
                         dipole_output_period_ps=1.0, max_energy_output_time=None,
                         # Hardware parameters
                         device='CPU', gpu_id=0, truncate_gsd=False, seed=None,
                         restart_velocities=True, zero_momentum_enabled=False, zero_momentum_period_ps=1.0,
                         # Empirical feedback parameters
                         enable_empirical_feedback=False, empirical_data_file=None,
                         feedback_output_period_ps=0.1, feedback_energy_component='lj_coulombic',
                         feedback_apply_to='both', feedback_averaging_window_ps=5.0,
                         feedback_update_interval_ps=5.0, feedback_T_min=0.0, feedback_T_max=None,
                         feedback_turn_off_time_ps=None, feedback_use_direct_harmonic=True,
                         feedback_auto_adjust_molecular_tau=False, 
                         enable_auto_stopping=False, # Renamed from feedback_auto_stopping
                         auto_stop_min_time_ps=10.0, 
                         auto_stop_smoothing_window=5,
                         # Experiment metadata
                         replica=0, 
                         incavity=True) -> bool:
    """
    Run a unified cavity MD experiment with support for ALL Phase 3 enhanced features.

    Args:
        molecular_thermo (str): Molecular thermostat type ('bussi', 'langevin', 'none')
        cavity_thermo (str): Cavity thermostat type ('bussi', 'langevin', 'none')
        finite_q (bool): Use finite-q cavity mode (default: q=0 mode)
        coupling (float): Cavity coupling strength / amplitude (default: 1e-3)
        temperature (float): System temperature in Kelvin (default: 100.0)
        frequency (float): Cavity frequency in cm⁻¹ (default: 1560.0)
        runtime_ps (float): Simulation runtime in ps (default: 500.0)
        input_gsd (str): Input GSD file path (default: molecular-0.gsd)
        frame (int): Frame number to read from GSD file (default: -1, last frame)
        # Enhanced coupling variant parameters
        coupling_variant_type (str): Coupling variant type ('constant', 'step', 'periodic', 'exponential', 'square', 'decaying_square', 'adaptive_square', 'exponentialwave', 'composite', 'lqg_square')
        # Enhanced step coupling parameters
        switch_time_ps (float): Time in ps to turn off step coupling (default: None, no turn-off)
        decay_time_constant_ps (float): Exponential decay time constant in ps for time-varying coupling (default: None)
        step_turn_off_time_ps (float): Time in ps to turn off step coupling (default: None, no turn-off)
        # Enhanced periodic coupling parameters
        periodic_coupling (bool): Enable periodic coupling (default: False)
        periodic_period_ps (float): Period in ps (default: 1.0)
        periodic_phase_offset (float): Phase offset in radians (default: 0.0)
        periodic_start_time_ps (float): Start time in ps (default: 0.0)
        periodic_stop_time_ps (float): Stop time in ps (default: None)
        # Exponential decay parameters
        exponential_amplitude (float): Exponential decay amplitude (default: None, uses --coupling value)
        exponential_decay_time_ps (float): Exponential decay time constant in ps (default: 10.0)
        exponential_turn_on_time_ps (float): Exponential turn-on time in ps (default: 0.0)
        exponential_turn_off_time_ps (float): Exponential turn-off time in ps (default: None, no turn-off)
        # Square wave parameters
        square_period_ps (float): Square wave period in ps (default: 2.0)
        square_duty_cycle (float): Square wave duty cycle (0.0-1.0) (default: 0.5)
        square_phase_offset (float): Square wave phase offset in radians (default: 0.0)
        square_start_time_ps (float): Square wave start time in ps (default: 0.0)
        square_stop_time_ps (float): Square wave stop time in ps (default: None, no stop)
        # Decaying square wave parameters
        decaying_square_period_ps (float): Decaying square wave period in ps (default: 2.0)
        decaying_square_duty_cycle (float): Decaying square wave duty cycle (0.0-1.0) (default: 0.5)
        decaying_square_phase_offset (float): Decaying square wave phase offset in radians (default: 0.0)
        decaying_square_start_time_ps (float): Decaying square wave start time in ps (default: 0.0)
        decaying_square_stop_time_ps (float): Decaying square wave stop time in ps (default: None, no stop)
        decaying_square_decay_rate (float): Decay rate per period 0.0-1.0 (default: 0.1 = 10%)
        decaying_square_minimum_amplitude (float): Minimum amplitude threshold (default: 1e-6)
        # Adaptive square wave parameters
        adaptive_square_period_ps (float): Adaptive square wave period in ps (default: 5.0)
        adaptive_square_duty_cycle (float): Adaptive square wave duty cycle (0.0-1.0) (default: 0.5)
        adaptive_square_phase_offset (float): Adaptive square wave phase offset in radians (default: 0.0)
        adaptive_square_start_time_ps (float): Adaptive square wave start time in ps (default: 0.0)
        adaptive_square_stop_time_ps (float): Adaptive square wave stop time in ps (default: None, no stop)
        adaptive_square_min_amplitude (float): Minimum allowed amplitude (safety limit) (default: 1e-8)
        adaptive_square_max_amplitude (float): Maximum allowed amplitude (safety limit) (default: 1e-1)
        # Exponential wave parameters
        exp_period_ps (float): Exponential wave period in ps (default: 2.0)
        exp_tau_ps (float): Exponential decay time constant in ps (default: 0.5)
        exp_start_time_ps (float): Exponential wave start time in ps (default: 0.0)
        exp_stop_time_ps (float): Exponential wave stop time in ps (default: None, no stop)
        exp_adaptive (bool): Enable adaptive amplitude scaling based on T_bath/T_harmonic ratio
        # Composite coupling parameters
        composite_sinusoid_amplitude (float): Amplitude for sinusoidal component (default: 1e-4)
        composite_sinusoid_period (float): Period for sinusoidal component in ps (default: 1.0)
        composite_sinusoid_phase (float): Phase offset for sinusoidal component in radians (default: 0.0)
        composite_sinusoid_start_time (float): Start time for sinusoidal component in ps (default: 0.0)
        composite_sinusoid_stop_time (float): Stop time for sinusoidal component in ps (default: None, runs indefinitely)
        composite_square_amplitude (float): Amplitude for adaptive square component (default: 2e-4)
        composite_square_period (float): Period for adaptive square component in ps (default: 50.0)
        composite_square_duty_cycle (float): Duty cycle for adaptive square component (default: 0.02)
        composite_square_start_time (float): Start time for adaptive square component in ps (default: 10.0)
        composite_square_stop_time (float): Stop time for square wave component in ps (default: 1000.0)
        composite_square_adaptive (bool): Use adaptive square wave (amplitude adjusts to temperature) instead of fixed amplitude
        composite_max_amplitude (float): Maximum total amplitude for composite (default: None, no limit)
        # Auto-stop coupling parameters
        enable_auto_stop (bool): Enable automatic coupling stop when temperatures converge (default: False)
        auto_stop_tol (float): Temperature tolerance for auto-stop in K (default: 1.0)
        auto_stop_window (float): Averaging window for auto-stop in ps (default: 10.0)
        # Enhanced laser drive parameters  
        laser_enabled (bool): Enable laser with timing control F_laser = F_L*cos(ω_L*t)
        laser_frequency_cm1 (float): Laser frequency in cm⁻¹ (default: 1560.0)
        laser_amplitude (float): Laser amplitude in a.u. (default: 1e-5)
        laser_start_time_ps (float): Time to start laser drive in ps (default: 0.0)
        laser_stop_time_ps (float): Time to stop laser drive in ps (default: None, runs indefinitely)
        laser_kvector (list): Laser k-vector components [kx, ky, kz] (default: [1.0, 0.0, 0.0])
        # Gradient descent feedback controller parameters
        enable_gd_feedback (bool): Enable gradient descent feedback temperature controller
        gd_target_temperature (float): GD controller target temperature in K (default: 100.0)
        gd_turn_on_time_ps (float): GD controller turn-on time in ps (default: 0.0)
        gd_turn_off_time_ps (float): GD controller turn-off time in ps (default: None, never turn off)
        gd_temperature_method (str): GD controller temperature method (default: kinetic)
        gd_update_interval_ps (float): GD controller update interval in ps (default: 0.1)
        gd_time_constant_ps (float): GD controller time constant in ps (default: 10.0)
        gd_apply_to (str): Which thermostats to control: molecular, cavity, or both (default: both)
        gd_T_min (float): GD controller minimum temperature in K (default: 0.0)
        gd_T_max (float): GD controller maximum temperature in K (default: None)
        gd_disable_effective_temp (bool): Use raw measured temperature instead of effective temperature mixing (default: False)
        # Multi-signal error function parameters
        gd_enable_multi_signal (bool): Enable multi-signal error function combining three error signals (default: False)
        gd_weight_system_target (float): Weight for (Ts - Ttarget) error signal (default: 1.0)
        gd_weight_bath_target (float): Weight for (Tbath - Ttarget) error signal (default: 0.0)  
        gd_weight_system_bath (float): Weight for (Ts - Tbath) error signal (default: 0.0)
        # Dual independent feedback controller parameters
        enable_dual_feedback (bool): Enable dual independent feedback controller (default: False)
        dual_cavity_method (str): Temperature method for cavity bath control (default: harmonic_equipartition)
        dual_molecular_method (str): Temperature method for molecular bath control (default: lj_coulombic_kinetic)
        dual_cavity_target_temperature (float): Target temperature for cavity bath in K (default: 100.0)
        dual_molecular_target_temperature (float): Target temperature for molecular bath in K (default: 100.0)
        dual_cavity_time_constant_ps (float): Time constant for cavity bath control in ps (default: 5.0)
        dual_molecular_time_constant_ps (float): Time constant for molecular bath control in ps (default: 10.0)
        dual_turn_on_time_ps (float): Turn-on time for dual controller in ps (default: 0.0)
        dual_turn_off_time_ps (float): Turn-off time for dual controller in ps (default: None, never turn off)
        dual_update_interval_ps (float): Update interval for dual controller in ps (default: 0.1)
        dual_cavity_T_min (float): Minimum cavity bath temperature in K (default: 0.0)
        dual_cavity_T_max (float): Maximum cavity bath temperature in K (default: None)
        dual_molecular_T_min (float): Minimum molecular bath temperature in K (default: 0.0)
        dual_molecular_T_max (float): Maximum molecular bath temperature in K (default: None)
        dual_cavity_dynamic_target (bool): Use dynamic target temperature for cavity (measured at turn-on time) instead of fixed target (default: False)
        dual_molecular_dynamic_target (bool): Use dynamic target temperature for molecular (measured at turn-on time) instead of fixed target (default: False)
        dual_cavity_integral_time_constant_ps (float): Integral time constant for cavity PI controller in ps (default: None, P-only). Smaller values = more aggressive integral action
        dual_molecular_integral_time_constant_ps (float): Integral time constant for molecular PI controller in ps (default: None, P-only). Smaller values = more aggressive integral action
        # Asymmetric gain factors
        dual_cavity_heating_gain_factor (float): Proportional gain factor for cavity heating (error < 0) (default: 1.0)
        dual_cavity_cooling_gain_factor (float): Proportional gain factor for cavity cooling (error > 0) (default: 1.0)
        dual_molecular_heating_gain_factor (float): Proportional gain factor for molecular heating (error < 0) (default: 1.0)
        dual_molecular_cooling_gain_factor (float): Proportional gain factor for molecular cooling (error > 0) (default: 1.0)
        # Asymmetric integral time constants
        dual_cavity_integral_heating_time_constant (float): Integral time constant for cavity heating in ps (default: None, use base integral constant)
        dual_cavity_integral_cooling_time_constant (float): Integral time constant for cavity cooling in ps (default: None, use base integral constant)
        dual_molecular_integral_heating_time_constant (float): Integral time constant for molecular heating in ps (default: None, use base integral constant)
        dual_molecular_integral_cooling_time_constant (float): Integral time constant for molecular cooling in ps (default: None, use base integral constant)
        # Sinusoidal bath temperature controller parameters
        enable_sinusoidal_bath (bool): Enable sinusoidal bath temperature controller (default: False)
        sinusoidal_bath_period_ps (float): Period of sinusoidal temperature oscillation in ps (default: 1.0)
        sinusoidal_bath_amplitude_scale (float): Scaling factor for amplitude based on harmonic equipartition deviation (default: 0.1)
        sinusoidal_bath_phase_offset (float): Phase offset for sinusoidal function in radians (default: 0.0)
        sinusoidal_bath_target_temperature (float): Static target temperature in K (ignored if dynamic target enabled) (default: 100.0)
        sinusoidal_bath_dynamic_target (bool): Set mean temperature dynamically at turn-on time (default: False)
        sinusoidal_bath_turn_on_time_ps (float): Time to start control in ps (default: 0.0)
        sinusoidal_bath_turn_off_time_ps (float): Time to stop control in ps (default: None = never turn off)
        sinusoidal_bath_update_interval_ps (float): Interval between control updates in ps (default: 0.1)
        sinusoidal_bath_apply_to (str): Which thermostats to control (default: both)
        sinusoidal_bath_T_min (float): Minimum allowed bath temperature in K (default: 0.1)
        sinusoidal_bath_T_max (float): Maximum allowed bath temperature in K (default: None)
        sinusoidal_bath_empirical_data_file (str): Path to empirical data file for harmonic equipartition calculation
        sinusoidal_bath_amplitude_update_interval_ps (float): Interval for updating amplitude based on harmonic temperature in ps (default: 1.0)
        sinusoidal_bath_amplitude_temperature_method (str): Temperature method for amplitude calculation (default: harmonic_equipartition)
        sinusoidal_bath_adaptive_range_mode (bool): Enable adaptive range mode: sinusoid bounds fixed at target, range adapts to |Ts-Ttarget| (default: False)
        # Adaptive bath temperature controller parameters
        enable_adaptive_bath (bool): Enable adaptive bath temperature controller (default: False)
        adaptive_bath_amplitude_scale (float): Amplitude scale for bath temperature adjustment (default: 1.0)
        adaptive_bath_time_constant_ps (float): Time constant tau for GD evolution in ps (default: 1.0)
        adaptive_bath_target_temperature (float): Static target temperature in K (ignored if --adaptive-bath-dynamic-target) (default: 100.0)
        adaptive_bath_dynamic_target (bool): Set target temperature dynamically at turn-on time (default: False)
        adaptive_bath_turn_on_time_ps (float): Time to start adaptive bath control in ps (default: 0.0)
        adaptive_bath_turn_off_time_ps (float): Time to stop adaptive bath control in ps (default: None)
        adaptive_bath_update_interval_ps (float): Update interval for adaptive bath control in ps (default: 0.1)
        adaptive_bath_apply_to (str): Which thermostats to apply adaptive bath control to (default: both)
        adaptive_bath_T_min (float): Minimum allowed bath temperature in K (default: 0.1)
        adaptive_bath_T_max (float): Maximum allowed bath temperature in K (default: None)
        adaptive_bath_empirical_data_file (str): Path to empirical data file for signal temperature calculation
        adaptive_bath_signal_temperature_method (str): Temperature method for control signal (default: harmonic_equipartition)
        adaptive_bath_dynamic_target_temperature_method (str): Temperature method for setting dynamic target (default: same as signal method)
        # Quench controller parameters
        enable_quench_controller (bool): Enable quench controller for instantaneous temperature changes (default: False)
        quench_initial_temperature (float): Initial temperature before quench in K (default: 100.0)
        quench_target_temperature (float): Target temperature after quench in K (default: 50.0)
        quench_time_ps (float): Time when quench occurs in ps (default: 50.0)
        quench_apply_to (str): Which thermostats to apply quench to (default: both)
        # Offset temperature controller parameters
        enable_offset_controller (bool): Enable offset temperature controller for fixed temperature offsets (default: False)
        offset_temperature_method (str): Temperature calculation method for offset controller (default: kinetic)
        offset_temperature_offset_K (float): Temperature offset in K (can be positive or negative) (default: -50.0)
        offset_turn_on_time_ps (float): Turn-on time for offset controller in ps (default: 0.0)
        offset_turn_off_time_ps (float): Turn-off time for offset controller in ps (default: None, never turn off)
        offset_update_interval_ps (float): Update interval for offset controller in ps (default: 0.1)
        offset_apply_to (str): Apply offset to molecular, cavity, or both thermostats (default: both)
        offset_T_min (float): Minimum allowed bath temperature in K (default: 0.0)
        offset_T_max (float): Maximum allowed bath temperature in K (default: None)
        # Harmonic bond reset parameters
        enable_harmonic_reset (bool): Enable one-time harmonic bond reset (thermal sampling of bond stretch DOF for ALL bond types)
        harmonic_reset_turn_on_time_ps (float): Time in ps to trigger harmonic bond reset (default: 0.0)
        harmonic_reset_temperature (float): Target temperature for internal stretch DOF (default: None, use dynamic bath temperature at reset time)
        harmonic_reset_seed (int): Random seed for harmonic reset (default: 42)
        # Differential equation controller parameters
        enable_diffeq_controller (bool): Enable differential equation temperature controller (default: False)
        diffeq_temperature_method (str): Temperature calculation method for differential equation controller (default: kinetic)
        diffeq_time_constant_ps (float): Time constant τ in ps for differential equation dT/dt = -(T_bath - T_signal)/τ (default: 5.0)
        diffeq_time_constant_auto (bool): Use adaptive time constant τ = T[T_bath] from relaxation model instead of fixed value (default: False)
        diffeq_turn_on_time_ps (float): Turn-on time for differential equation controller in ps (default: 0.0)
        diffeq_turn_off_time_ps (float): Turn-off time for differential equation controller in ps (default: None, never turn off)
        diffeq_update_interval_ps (float): Update interval for differential equation controller in ps (default: 0.1)
        diffeq_apply_to (str): Apply differential equation control to molecular, cavity, or both thermostats (default: both)
        diffeq_T_min (float): Minimum allowed bath temperature in K (default: 0.0)
        diffeq_T_max (float): Maximum allowed bath temperature in K (default: None)
        diffeq_rate_limit_K_per_ps (float): Maximum rate of temperature change in K/ps (default: None, no limit)
        # Exact-cancellation + PI control parameters for DiffEq controller
        diffeq_enable_pi_control (bool): Enable exact-cancellation + PI control in DiffEq controller (default: False)
        diffeq_pi_rho (float): Expected disturbance rate for PI control in K/ps (default: 1.0)
        diffeq_pi_epsilon (float): Target tolerance for PI control in K (default: 0.5)
        diffeq_pi_zeta (float): Damping coefficient for PI control (default: 0.8)
        diffeq_relaxation_data_file (str): Path to relaxation time data file for exact cancellation (default: None)
        diffeq_filter_window (float): Moving window size for filtering T_signal measurements in ps (default: 0.0)
        # LQR controller parameters
        enable_lqr_controller (bool): Enable LQR optimal temperature controller with Kalman filter (default: False)
        lqr_signal_method (str): Temperature method for regulated signal (default: lj_coulombic)
        lqr_hot_method (str): Temperature method for disturbance signal (default: harmonic_equipartition)
        lqr_target_temperature (float): Fixed target temperature in K (used if not dynamic) (default: 300.0)
        lqr_dynamic_target (bool): Set target temperature dynamically from signal at turn-on (default: False)
        lqr_dynamic_target_method (str): Temperature method for dynamic target (default: same as signal method)
        lqr_weight_signal (float): LQR weight for signal regulation (higher = tighter) (default: 100.0)
        lqr_weight_hot (float): LQR weight for hot temperature (lower = allow more variation) (default: 1.0)
        lqr_weight_bath (float): [DEPRECATED - LQG uses 2D state, bath is control input] Old 3D formulation parameter (default: 0.1)
        lqr_weight_integral (float): LQR weight for integral action (default: 10.0)
        lqr_control_effort (float): LQR penalty on control effort (higher = gentler) (default: 1.0)
        lqr_process_noise_signal (float): Process noise std dev for signal (K) (default: 0.1)
        lqr_process_noise_hot (float): Process noise std dev for hot (K) (default: 0.1)
        lqr_measurement_noise_signal (float): Measurement noise std dev for signal (K) (default: 0.5)
        lqr_measurement_noise_hot (float): Measurement noise std dev for hot (K) (default: 0.5)
        lqr_system_id_mode (str): System ID mode: step (single quench), multi_step (4 sequential excitations), load (from file), skip (manual), none (adaptive no-quench) (default: step)
        lqr_system_id_temp_K (float): Quench temperature for system ID (K) (default: 5.0)
        lqr_system_id_duration_ps (float): Duration of system ID cold quench (ps) (default: 50.0)
        lqr_system_id_file (str): File to save/load system parameters (default: lqr_system_params.json)
        lqr_periodic_system_id (bool): Enable periodic re-identification of system parameters
        lqr_periodic_system_id_interval_ps (float): Time interval between periodic system IDs (ps) (default: 1000.0)
        # EKF adaptation
        lqr_use_ekf_adaptation (bool): Use EKF for online parameter estimation (default: True)
        lqr_ekf_update_interval (int): EKF update interval (control steps) (default: 50)
        lqr_ekf_process_noise_param (float): EKF parameter drift noise std (default: 0.001)
        lqr_ekf_initial_covariance_param (float): EKF initial parameter uncertainty (default: 0.1)
        lqr_adaptive_lqr_threshold (float): Relative parameter change threshold for LQR redesign (default: 0.05 = 5%)
        # Gain scheduling
        lqr_enable_gain_scheduling (bool): Enable gain scheduling based on temperature spread (default: True)
        lqr_gain_schedule_far_threshold (float): Temperature spread for full gain (K) (default: 20.0)
        lqr_gain_schedule_near_threshold (float): Temperature spread for reduced gain (K) (default: 10.0)
        # T_h low-pass filter
        lqr_th_filter_enabled (bool): Enable low-pass filtering on T_h measurements (default: True)
        lqr_th_filter_time_constant (float): T_h filter time constant in ps (default: 20.0, cutoff ~8Hz)
        # Gentle startup
        lqr_gentle_startup_steps (int): Number of steps to ramp up control authority (default: 10)
        lqr_gentle_startup_min_authority (float): Starting control authority fraction (default: 0.1 = 10%)
        # Kinetic temperature tracking (3D state)
        lqr_track_kinetic_temp (bool): Enable kinetic temperature as 3rd state variable (default: False, uses 2D system)
        lqr_weight_kinetic (float): LQR weight for kinetic temperature error (default: 100.0, same as signal temp)
        lqr_process_noise_kinetic (float): Kalman filter process noise for kinetic temp (K) (default: 2.0)
        lqr_measurement_noise_kinetic (float): Kalman filter measurement noise for kinetic temp (K) (default: 2.0)
        # Cross-coupling weights for thermal equilibration
        lqr_cross_coupling_signal_kinetic (float): Cross-coupling weight for (T_signal - T_kinetic)² penalty (default: 0.0, disabled)
        lqr_cross_coupling_signal_hot (float): Cross-coupling weight for (T_signal - T_hot)² penalty (default: 0.0, disabled)
        lqr_cross_coupling_hot_kinetic (float): Cross-coupling weight for (T_hot - T_kinetic)² penalty (default: 0.0, disabled)
        # Adaptive LQI Controller parameters
        lqr_controller_type (str): Controller type: standard (basic LQR) or adaptive_lqi (full mode-based) (default: standard)
        lqr_tau_L_initial (float): Initial guess for signal time constant (ps) for adaptive_lqi (default: 200.0)
        lqr_tau_H_initial (float): Initial guess for hot time constant (ps) for adaptive_lqi (default: 30.0)
        lqr_k_initial (float): Initial guess for coupling coefficient for adaptive_lqi (default: 0.01)
        lqr_tau_b (float): Bath actuator time constant (ps) for adaptive_lqi (default: 1.0)
        lqr_q_common (float): Common mode state weight for adaptive_lqi (default: 1000.0)
        lqr_q_diff (float): Difference mode state weight for adaptive_lqi (default: 100.0)
        lqr_q_eta_common (float): Common mode integral weight for adaptive_lqi (default: 1e6)
        lqr_q_eta_diff (float): Difference mode integral weight for adaptive_lqi (default: 10.0)
        lqr_process_noise_drift (float): Drift process noise for adaptive_lqi (default: 0.01)
        lqr_rls_forgetting (float): RLS forgetting factor for adaptive_lqi (default: 0.998)
        lqr_rls_regularization (float): RLS regularization for adaptive_lqi (default: 1e-6)
        lqr_rls_update_interval (int): RLS update interval (samples) for adaptive_lqi (default: 50)
        lqr_max_control_rate (float): Maximum control rate (K/ps) for adaptive_lqi (default: 10.0)
        lqr_integral_max_common (float): Maximum common mode integral for adaptive_lqi (default: 1000.0)
        lqr_integral_max_diff (float): Maximum difference mode integral for adaptive_lqi (default: 100.0)
        lqr_theta_change_threshold (float): Parameter change threshold for relinearization for adaptive_lqi (default: 0.05)
        # LQG pulse controller parameters
        enable_lqg_controller (bool): Enable LQG pulse-to-pulse optimal controller (default: False)
        lqg_A_matrix (str): JSON string or file path for A matrix (state transition)
        lqg_B_matrix (str): JSON string or file path for B matrix (control input)
        lqg_C_matrix (str): JSON string or file path for C matrix (measurement)
        lqg_S_matrix (str): JSON string or file path for S matrix (tracked outputs)
        lqg_r_target (str): JSON string or file path for target reference vector
        lqg_Qz_weights (str): JSON string or file path for output tracking penalty matrix
        lqg_Ru_weights (str): JSON string or file path for control effort penalty matrix
        lqg_Qe_weights (str): JSON string or file path for integral error penalty matrix
        lqg_W_noise (str): JSON string or file path for process noise covariance matrix
        lqg_V_noise (str): JSON string or file path for measurement noise covariance matrix
        lqg_g0_baseline (float): Baseline coupling strength g0 in a.u. (default: 1e-4)
        lqg_dt (float): LQG pulse period in ps (default: 5.0)
        lqg_delta_g_min (float): Minimum coupling deviation Δg in a.u. (default: -0.3)
        lqg_delta_g_max (float): Maximum coupling deviation Δg in a.u. (default: 0.3)
        lqg_T_bath_min (float): Minimum bath temperature in K (default: 10.0)
        lqg_T_bath_max (float): Maximum bath temperature in K (default: 500.0)
        lqg_turn_on_time_ps (float): Turn-on time for LQG controller in ps (default: 0.0)
        lqg_turn_off_time_ps (float): Turn-off time for LQG controller in ps (default: None)
        lqg_update_interval_ps (float): LQG update interval in ps (default: 5.0)
        lqg_console_output_period_ps (float): Console output period for LQG in ps (default: 10.0)
        lqg_c_affine_term (str): JSON string or file path for affine term c vector
        lqg_temperature_methods (str): JSON string mapping state indices to temperature methods
        lqg_square_period_ps (float): LQG square wave period in ps (default: 5.0)
        lqg_square_duty_cycle (float): LQG square wave duty cycle 0.0-1.0 (default: 0.5)
        lqg_square_phase_offset (float): LQG square wave phase offset in radians (default: 0.0)
        lqg_square_start_time_ps (float): LQG square wave start time in ps (default: 0.0)
        lqg_square_stop_time_ps (float): LQG square wave stop time in ps (default: None)
        lqg_square_g_low_level (float): LQG square wave low coupling level in a.u. (default: 1e-6)
        lqg_enable_system_id (bool): Enable automatic system identification (default behavior, can be overridden by --lqg-disable-system-id)
        lqg_disable_system_id (bool): Disable system identification (use provided matrices or seed model)
        lqg_equilibrium_duration_ps (float): Duration for equilibrium data collection in ps (default: 500.0)
        lqg_prbs_n_unique_steps (int): Number of unique coupling levels for PRBS excitation (default: 10)
        lqg_prbs_amplitude (float): PRBS excitation amplitude in a.u. (default: 5e-5)
        lqg_bath_n_unique_steps (int): Number of unique bath temperature levels for excitation (default: 8)
        lqg_bath_excitation_amplitude (float): Bath temperature excitation amplitude in K (default: 10.0)
        lqg_system_id_output_file (str): Output file for identified model (default: lqg_system_id.json)
        
        # LQG coupling controller parameters (new coupling-type: lqg_coupling)
        lqg_coupling_target_temperature (float): Target temperature for T_s regulation in K (default: 100.0)
        lqg_coupling_update_interval_ps (float): LQG coupling controller update period in ps (default: 0.1)
        lqg_coupling_equilibrium_duration_ps (float): Initial equilibrium phase duration for system ID in ps (default: 5.0)
        lqg_coupling_step_duration_ps (float): Duration for each step in system ID excitation in ps (default: 5.0)
        lqg_coupling_n_steps (int): Number of step up/down cycles for system ID (default: 2)
        lqg_coupling_lambda_min (float): Minimum allowed λ value (default: 0.0)
        lqg_coupling_lambda_max (float): Maximum allowed λ value (default: 1e-2)
        lqg_coupling_process_noise_std (float): Process noise standard deviation for Kalman filter in K (default: 0.1)
        lqg_coupling_measurement_noise_std (float): Measurement noise standard deviation for Kalman filter in K (default: 0.5)
        lqg_coupling_weight_signal (float): LQR weight for T_s regulation (higher = tighter control) (default: 100.0)
        lqg_coupling_weight_harmonic (float): LQR weight for T_harmonic observation (default: 1.0)
        lqg_coupling_weight_kinetic (float): LQR weight for T_kinetic observation (default: 1.0)
        lqg_coupling_weight_bath (float): LQR weight for T_bath observation (default: 0.1)
        lqg_coupling_control_effort (float): LQR penalty on control effort (higher = gentler) (default: 1.0)
        lqg_coupling_integral_gain (float): Integral action gain for zero steady-state error (default: 2.0)
        lqg_coupling_system_id_file (str): JSON file to save identified system parameters (default: lqg_coupling_system_id.json)
        lqg_coupling_control_file (str): CSV file to save control trajectory (default: lqg_coupling_control.csv)
        lqg_coupling_temperature_methods (list): Temperature methods for [T_s, T_harmonic, T_kinetic] (default: ['lj_coulombic', 'harmonic_equipartition', 'kinetic'])
        lqg_coupling_bath_temperature_method (str): Temperature method for T_bath observation (default: kinetic)
        
        enable_temp_tracker (bool): Enable comprehensive temperature tracker (kinetic, fictive, bath temperatures)
        temp_tracker_output_period_ps (float): Temperature tracker output period in ps (default: 0.1)
        temp_tracker_empirical_data_file (str): Empirical data file for LJ+Coul fictive temperature (default: None)
        enable_molecular_temps (bool): Enable molecular temperature decomposition (translational/rotational/vibrational)
        molecular_temps_output_period_ps (float): Molecular temperature output period in ps (default: 1.0)
        enable_dipole_fdr (bool): Enable dipole moment FDR analysis (autocorrelation tracking)
        dipole_fdr_output_period_ps (float): Dipole FDR output period in ps (default: 0.1)
        dipole_fdr_max_correlation_time_ps (float): Maximum correlation time for dipole FDR in ps (default: 100.0)
        dipole_fdr_field_direction (list): Electric field direction for FDR [x, y, z] (default: [0, 0, 1])
        dipole_fdr_exclude_cavity (bool): Exclude cavity particles from dipole calculation (default: True)
        enable_dipole_response (bool): Enable dipole response force for fork-and-clone FDR measurements
        dipole_response_field_strength (float): Electric field strength for response measurements (default: 1e-5)
        dipole_response_sign (float): Sign of electric field (+1 for plus clone, -1 for minus clone) (default: +1.0)
        mech_periodic (bool): Enable mechanical cavity length modulation: L(t) = L₀ + A*sin(ω*t)
        mech_frequency_cm1 (float): Mechanical modulation frequency in cm⁻¹ (default: 100.0)
        mech_magnitude (float): Mechanical modulation magnitude in a.u. (default: 1e-4)
        switch_time_ps (float): Time in ps when coupling and dissipation turn on (default: on from start)
        decay_time_constant_ps (float): Exponential decay time constant in ps for time-varying coupling (default: None)
        damping_ratio (float): Damping ratio zeta for dissipation calculation (default: 0.0)
        molecular_tau (float): Molecular thermostat time constant in ps (default: 5.0)
        cavity_tau (float): Cavity thermostat time constant in ps (default: 5.0)
        fixed_timestep (bool): Use fixed timestep instead of adaptive timestep
        timestep_fs (float): Fixed timestep in fs (only used with --fixed-timestep) (default: 0.1)
        error_tolerance (float): Error tolerance for adaptive timestep (ignored with --fixed-timestep) (default: 0.01)
        initial_fraction (float): Initial fraction for shock dampening (default: 1e-5)
        time_constant_ps (float): Time constant for error tolerance ramping in ps (default: 50.0)
        enable_dynamic_coupling_detection (bool): Enable dynamic coupling change detection for shock dampening (default: True)
        coupling_change_threshold (float): Threshold for detecting large coupling changes in a.u. (default: 1e-5)
        enable_energy_tracker (bool): Enable comprehensive energy tracking
        energy_output_period_ps (float): Energy output period in ps (default: 0.1)
        fkt_output_period_ps (float): F(k,t) output period in ps (default: 1.0)
        gsd_output_period_ps (float): GSD trajectory output period in ps (default: 50.0)
        console_output_period_ps (float): Console output period in ps (default: 1.0)
        enable_fkt (bool): Enable F(k,t) density correlation tracking
        fkt_kmag (float): F(k,t) k-magnitude (default: 1.0)
        fkt_wavevectors (int): Number of F(k,t) wavevectors (default: 50)
        fkt_ref_interval (float): F(k,t) reference interval in ps (default: 1.0)
        fkt_max_refs (int): Maximum F(k,t) reference states (default: 10)
        enable_dipole_autocorr (bool): Enable dipole autocorrelation tracking
        dipole_ref_interval (float): Dipole reference interval in ps (default: 1.0)
        dipole_max_refs (int): Maximum dipole reference states (default: 10)
        dipole_output_period_ps (float): Dipole output period in ps (default: 1.0)
        max_energy_output_time (float): Maximum time for energy output in ps (default: None, no limit)
        device (str): Compute device (default: CPU)
        gpu_id (int): GPU device ID (default: 0)
        truncate_gsd (bool): Truncate existing GSD file (default: append)
        seed (int): Random seed for simulation (default: replica-based deterministic seed)
        no_restart_velocities (bool): Do not restart velocities - use existing velocities from GSD file (default: restart velocities)
        zero_momentum_enabled (bool): Enable periodic momentum zeroing to prevent center-of-mass drift (default: disabled)
        zero_momentum_period_ps (float): Period for momentum zeroing in ps (default: 1.0)
        enable_empirical_feedback (bool): Enable empirical temperature feedback control based on systemic temperature
        empirical_data_file (str): Path to empirical energy-temperature data file (required if --enable-empirical-feedback)
        feedback_output_period_ps (float): Empirical feedback CSV output period in ps (default: 0.1)
        feedback_energy_component (str): Energy component for empirical feedback (default: lj_coulombic)
        feedback_apply_to (str): Which thermostats to apply empirical feedback to (default: both)
        feedback_averaging_window_ps (float): Time window for averaging fictive temperature in ps (default: 5.0)
        feedback_update_interval_ps (float): Interval between temperature updates in ps (default: 5.0)
        feedback_T_min (float): Minimum allowed temperature in K (default: 50.0)
        feedback_T_max (float): Maximum allowed temperature in K (default: 300.0)
        feedback_turn_off_time_ps (float): Time in ps to turn off empirical feedback (default: None, feedback runs until end)
        feedback_use_direct_harmonic (bool): Use direct harmonic calculation for fictive temperature (default: empirical fitting)
        feedback_auto_adjust_molecular_tau (bool): Auto-adjust molecular thermostat tau when feedback turns off (default: disabled)
        enable_auto_stopping (bool): Enable automatic stopping when fictive and harmonic temperatures converge (default: disabled)
        auto_stop_min_time_ps (float): Minimum time before auto-stop can trigger in ps (default: 10.0)
        auto_stop_smoothing_window (int): Window size for temperature smoothing in auto-stop detection (default: 5)
        replica (int): Replica number
        incavity (bool): Run without cavity coupling (molecular-only simulation)

    Returns:
        bool: True if the experiment was successful, False otherwise
    """
    try:
        # Handle backward compatibility between coupling and lambda_coupling
        if coupling is not None and lambda_coupling != 1e-3:
            raise ValueError("Cannot specify both --coupling (deprecated) and --lambda-coupling simultaneously")
        if coupling is not None:
            print("WARNING: Using deprecated --coupling parameter. Please migrate to --lambda-coupling.")
            # Convert epsilon to lambda for backward compatibility
            omegac = frequency / 219474.6  # Convert cm^-1 to a.u.
            effective_lambda_coupling = coupling / omegac
        else:
            effective_lambda_coupling = lambda_coupling
        
        # Determine controller type for directory naming
        controller_type = "nocontrol"
        if enable_pid_controller:
            if pid_self_loop:
                controller_type = "pid_selfloop"
            else:
                controller_type = "pid"
        elif enable_diffeq_controller:
            if diffeq_enable_pi_control:
                controller_type = "diffeq_pi"
            elif diffeq_enable_bias_cancellation:
                controller_type = "diffeq_biascancel"
            else:
                controller_type = "diffeq"
        elif enable_lqr_controller:
            controller_type = "lqr"
        elif enable_lqg_controller:
            controller_type = "lqg"
        elif enable_adaptive_mpc_controller:
            controller_type = "mpc"
        elif enable_gd_feedback:
            controller_type = "gd"
        elif enable_dual_feedback:
            controller_type = "dual"
        elif enable_offset_controller:
            controller_type = "offset"
        elif enable_simple_setpoint_controller:
            controller_type = "setpoint"
        elif enable_quench_controller:
            controller_type = "quench"
        
        # Create descriptive job directory name
        # Format: <coupling_type>_lambda<value>_<controller>
        # Example: constant_lambda0.025_diffeq
        # Note: If multiple replicas need the same parameters, use different output directories
        lambda_str = f"{effective_lambda_coupling:.4f}".rstrip('0').rstrip('.')
        job_dir = f"{coupling_variant_type}_lambda{lambda_str}_{controller_type}"
        
        # Create the job directory if it doesn't exist
        import os
        os.makedirs(job_dir, exist_ok=True)
        
        print(f"Output directory: {job_dir}")
        print(f"  Coupling type: {coupling_variant_type}")
        print(f"  Lambda coupling: {effective_lambda_coupling}")
        print(f"  Controller: {controller_type}")
        print(f"  Replica: {replica}")
        print()
        
        sim = CavityMDSimulation(
            # Required parameters (must be first)
            job_dir=job_dir,
            replica=replica,
            freq=frequency,
            lambda_coupling=effective_lambda_coupling,
            incavity=incavity,
            
            # Basic simulation parameters
            name='prod',
            runtime_ps=runtime_ps,
            input_gsd=input_gsd,
            frame=frame,
            temperature=temperature,
            molecular_thermostat=molecular_thermo,
            cavity_thermostat=cavity_thermo,
            finite_q=finite_q,
            
            # Enhanced coupling variant parameters
            coupling_variant_type=coupling_variant_type,
            # Enhanced periodic coupling parameters
            periodic_coupling=periodic_coupling,
            periodic_period_ps=periodic_period_ps,
            periodic_phase_offset=periodic_phase_offset,
            periodic_start_time_ps=periodic_start_time_ps,
            periodic_stop_time_ps=periodic_stop_time_ps,
            # Exponential decay parameters
            exponential_amplitude=exponential_amplitude,
            exponential_decay_time_ps=exponential_decay_time_ps,
            exponential_turn_on_time_ps=exponential_turn_on_time_ps,
            exponential_turn_off_time_ps=exponential_turn_off_time_ps,
            # Square wave parameters
            square_period_ps=square_period_ps,
            square_duty_cycle=square_duty_cycle,
            square_phase_offset=square_phase_offset,
            square_start_time_ps=square_start_time_ps,
            square_stop_time_ps=square_stop_time_ps,
            # Decaying square wave parameters
            decaying_square_period_ps=decaying_square_period_ps,
            decaying_square_duty_cycle=decaying_square_duty_cycle,
            decaying_square_phase_offset=decaying_square_phase_offset,
            decaying_square_start_time_ps=decaying_square_start_time_ps,
            decaying_square_stop_time_ps=decaying_square_stop_time_ps,
            decaying_square_decay_rate=decaying_square_decay_rate,
            decaying_square_minimum_amplitude=decaying_square_minimum_amplitude,
            # Adaptive square wave parameters
            adaptive_square_period_ps=adaptive_square_period_ps,
            adaptive_square_duty_cycle=adaptive_square_duty_cycle,
            adaptive_square_phase_offset=adaptive_square_phase_offset,
            adaptive_square_start_time_ps=adaptive_square_start_time_ps,
            adaptive_square_stop_time_ps=adaptive_square_stop_time_ps,
            adaptive_square_min_amplitude=adaptive_square_min_amplitude,
            adaptive_square_max_amplitude=adaptive_square_max_amplitude,
            # Exponential wave parameters
            exp_period_ps=exp_period_ps,
            exp_tau_ps=exp_tau_ps,
            exp_start_time_ps=exp_start_time_ps,
            exp_stop_time_ps=exp_stop_time_ps,
            exp_adaptive=exp_adaptive,
            # Composite coupling parameters
            composite_sinusoid_amplitude=composite_sinusoid_amplitude,
            composite_sinusoid_period=composite_sinusoid_period,
            composite_sinusoid_phase=composite_sinusoid_phase,
            composite_sinusoid_start_time=composite_sinusoid_start_time,
            composite_sinusoid_stop_time=composite_sinusoid_stop_time,
            composite_square_amplitude=composite_square_amplitude,
            composite_square_period=composite_square_period,
            composite_square_duty_cycle=composite_square_duty_cycle,
            composite_square_start_time=composite_square_start_time,
            composite_square_stop_time=composite_square_stop_time,
            composite_square_adaptive=composite_square_adaptive,
            composite_max_amplitude=composite_max_amplitude,
            # Auto-stop coupling parameters
            enable_auto_stop=enable_auto_stop,
            auto_stop_tol=auto_stop_tol,
            auto_stop_window=auto_stop_window,
            # Enhanced laser drive parameters
            laser_enabled=laser_enabled,
            laser_frequency_cm1=laser_frequency_cm1,
            laser_amplitude=laser_amplitude,
            laser_start_time_ps=laser_start_time_ps,
            laser_stop_time_ps=laser_stop_time_ps,
            laser_kvector=laser_kvector,
            # Gradient descent feedback controller parameters
            enable_gd_feedback=enable_gd_feedback,
            gd_target_temperature=gd_target_temperature,
            gd_turn_on_time_ps=gd_turn_on_time_ps,
            gd_turn_off_time_ps=gd_turn_off_time_ps,
            gd_temperature_method=gd_temperature_method,
            gd_update_interval_ps=gd_update_interval_ps,
            gd_time_constant_ps=gd_time_constant_ps,
            gd_apply_to=gd_apply_to,
            gd_T_min=gd_T_min,
            gd_T_max=gd_T_max,
            gd_disable_effective_temp=gd_disable_effective_temp,
            # Multi-signal error function parameters
            gd_enable_multi_signal=gd_enable_multi_signal,
            gd_weight_system_target=gd_weight_system_target,
            gd_weight_bath_target=gd_weight_bath_target,
            gd_weight_system_bath=gd_weight_system_bath,
            # Dual independent feedback controller parameters
            enable_dual_feedback=enable_dual_feedback,
            dual_cavity_method=dual_cavity_method,
            dual_molecular_method=dual_molecular_method,
            dual_cavity_target_temperature=dual_cavity_target_temperature,
            dual_molecular_target_temperature=dual_molecular_target_temperature,
            dual_cavity_time_constant_ps=dual_cavity_time_constant_ps,
            dual_molecular_time_constant_ps=dual_molecular_time_constant_ps,
            dual_turn_on_time_ps=dual_turn_on_time_ps,
            dual_turn_off_time_ps=dual_turn_off_time_ps,
            dual_update_interval_ps=dual_update_interval_ps,
            dual_cavity_T_min=dual_cavity_T_min,
            dual_cavity_T_max=dual_cavity_T_max,
            dual_molecular_T_min=dual_molecular_T_min,
            dual_molecular_T_max=dual_molecular_T_max,
            dual_cavity_dynamic_target=dual_cavity_dynamic_target,
            dual_molecular_dynamic_target=dual_molecular_dynamic_target,
            dual_cavity_integral_time_constant_ps=dual_cavity_integral_time_constant_ps,
            dual_molecular_integral_time_constant_ps=dual_molecular_integral_time_constant_ps,
            # Sinusoidal bath temperature controller parameters
            enable_sinusoidal_bath=enable_sinusoidal_bath,
            sinusoidal_bath_period_ps=sinusoidal_bath_period_ps,
            sinusoidal_bath_amplitude_scale=sinusoidal_bath_amplitude_scale,
            sinusoidal_bath_phase_offset=sinusoidal_bath_phase_offset,
            sinusoidal_bath_target_temperature=sinusoidal_bath_target_temperature,
            sinusoidal_bath_dynamic_target=sinusoidal_bath_dynamic_target,
            sinusoidal_bath_turn_on_time_ps=sinusoidal_bath_turn_on_time_ps,
            sinusoidal_bath_turn_off_time_ps=sinusoidal_bath_turn_off_time_ps,
            sinusoidal_bath_update_interval_ps=sinusoidal_bath_update_interval_ps,
            sinusoidal_bath_apply_to=sinusoidal_bath_apply_to,
            sinusoidal_bath_T_min=sinusoidal_bath_T_min,
            sinusoidal_bath_T_max=sinusoidal_bath_T_max,
            sinusoidal_bath_empirical_data_file=sinusoidal_bath_empirical_data_file,
            sinusoidal_bath_amplitude_update_interval_ps=sinusoidal_bath_amplitude_update_interval_ps,
            sinusoidal_bath_amplitude_temperature_method=sinusoidal_bath_amplitude_temperature_method,
            sinusoidal_bath_adaptive_range_mode=sinusoidal_bath_adaptive_range_mode,
            # Adaptive bath temperature controller parameters
            enable_adaptive_bath=enable_adaptive_bath,
            adaptive_bath_amplitude_scale=adaptive_bath_amplitude_scale,
            adaptive_bath_time_constant_ps=adaptive_bath_time_constant_ps,
            adaptive_bath_target_temperature=adaptive_bath_target_temperature,
            adaptive_bath_dynamic_target=adaptive_bath_dynamic_target,
            adaptive_bath_turn_on_time_ps=adaptive_bath_turn_on_time_ps,
            adaptive_bath_turn_off_time_ps=adaptive_bath_turn_off_time_ps,
            adaptive_bath_update_interval_ps=adaptive_bath_update_interval_ps,
            adaptive_bath_apply_to=adaptive_bath_apply_to,
            adaptive_bath_T_min=adaptive_bath_T_min,
            adaptive_bath_T_max=adaptive_bath_T_max,
            adaptive_bath_empirical_data_file=adaptive_bath_empirical_data_file,
            adaptive_bath_signal_temperature_method=adaptive_bath_signal_temperature_method,
            adaptive_bath_dynamic_target_temperature_method=adaptive_bath_dynamic_target_temperature_method,
            # Quench controller parameters
            enable_quench_controller=enable_quench_controller,
            quench_initial_temperature=quench_initial_temperature,
            quench_target_temperature=quench_target_temperature,
            quench_time_ps=quench_time_ps,
            quench_apply_to=quench_apply_to,
            # Offset temperature controller parameters
            enable_offset_controller=enable_offset_controller,
            offset_temperature_method=offset_temperature_method,
            offset_temperature_offset_K=offset_temperature_offset_K,
            offset_turn_on_time_ps=offset_turn_on_time_ps,
            offset_turn_off_time_ps=offset_turn_off_time_ps,
            offset_update_interval_ps=offset_update_interval_ps,
            offset_apply_to=offset_apply_to,
            offset_T_min=offset_T_min,
            offset_T_max=offset_T_max,
            # Harmonic bond reset parameters
            enable_harmonic_reset=enable_harmonic_reset,
            harmonic_reset_turn_on_time_ps=harmonic_reset_turn_on_time_ps,
            harmonic_reset_temperature=harmonic_reset_temperature,
            harmonic_reset_seed=harmonic_reset_seed,
            # Differential equation controller parameters
            enable_diffeq_controller=enable_diffeq_controller,
            diffeq_temperature_method=diffeq_temperature_method,
            diffeq_time_constant_ps=diffeq_time_constant_ps,
            diffeq_time_constant_auto=diffeq_time_constant_auto,
            diffeq_turn_on_time_ps=diffeq_turn_on_time_ps,
            diffeq_turn_off_time_ps=diffeq_turn_off_time_ps,
            diffeq_update_interval_ps=diffeq_update_interval_ps,
            diffeq_apply_to=diffeq_apply_to,
            diffeq_T_min=diffeq_T_min,
            diffeq_T_max=diffeq_T_max,
            diffeq_rate_limit_K_per_ps=diffeq_rate_limit_K_per_ps,
            # Exact-cancellation + PI control parameters for DiffEq controller
            diffeq_enable_pi_control=diffeq_enable_pi_control,
            diffeq_pi_rho=diffeq_pi_rho,
            diffeq_pi_epsilon=diffeq_pi_epsilon,
            diffeq_pi_zeta=diffeq_pi_zeta,
            diffeq_relaxation_data_file=diffeq_relaxation_data_file,
            diffeq_filter_window=diffeq_filter_window,
            # Adaptive bias cancellation parameters for DiffEq controller
            diffeq_enable_bias_cancellation=diffeq_enable_bias_cancellation,
            diffeq_bias_tau_b_ps=diffeq_bias_tau_b_ps,
            diffeq_bias_tau_b_auto=diffeq_bias_tau_b_auto,
            diffeq_bias_kappa=diffeq_bias_kappa,
            diffeq_bias_kappa_auto=diffeq_bias_kappa_auto,
            diffeq_bias_tau_b_prefactor=diffeq_bias_tau_b_prefactor,
            diffeq_bias_kappa_prefactor=diffeq_bias_kappa_prefactor,
            diffeq_bias_calibration_time_ps=diffeq_bias_calibration_time_ps,
            # Simple Setpoint Controller parameters
            enable_simple_setpoint_controller=enable_simple_setpoint_controller,
            simple_setpoint_signal_method=simple_setpoint_signal_method,
            simple_setpoint_time_constant_ps=simple_setpoint_time_constant_ps,
            simple_setpoint_apply_to=simple_setpoint_apply_to,
            simple_setpoint_turn_on_time_ps=simple_setpoint_turn_on_time_ps,
            simple_setpoint_turn_off_time_ps=simple_setpoint_turn_off_time_ps,
            simple_setpoint_update_interval_ps=simple_setpoint_update_interval_ps,
            simple_setpoint_T_min=simple_setpoint_T_min,
            simple_setpoint_T_max=simple_setpoint_T_max,
            simple_setpoint_output_file=simple_setpoint_output_file,
            simple_setpoint_console_output_period_ps=simple_setpoint_console_output_period_ps,
            # Adaptive MPC controller parameters
            enable_adaptive_mpc_controller=enable_adaptive_mpc_controller,
            adaptive_mpc_target_temperature=adaptive_mpc_target_temperature,
            adaptive_mpc_turn_on_time_ps=adaptive_mpc_turn_on_time_ps,
            adaptive_mpc_turn_off_time_ps=adaptive_mpc_turn_off_time_ps,
            adaptive_mpc_system_id_duration_ps=adaptive_mpc_system_id_duration_ps,
            adaptive_mpc_system_id_step_duration_ps=adaptive_mpc_system_id_step_duration_ps,
            adaptive_mpc_system_id_seed=adaptive_mpc_system_id_seed,
            adaptive_mpc_update_interval_ps=adaptive_mpc_update_interval_ps,
            adaptive_mpc_prediction_horizon=adaptive_mpc_prediction_horizon,
            adaptive_mpc_control_horizon=adaptive_mpc_control_horizon,
            adaptive_mpc_output_weight=adaptive_mpc_output_weight,
            adaptive_mpc_control_effort_weight=adaptive_mpc_control_effort_weight if adaptive_mpc_control_effort_weight else [1.0, 0.1],
            adaptive_mpc_rate_penalty_weight=adaptive_mpc_rate_penalty_weight if adaptive_mpc_rate_penalty_weight else [10.0, 1.0],
            adaptive_mpc_lambda_min=adaptive_mpc_lambda_min,
            adaptive_mpc_lambda_max=adaptive_mpc_lambda_max,
            adaptive_mpc_T_bath_min=adaptive_mpc_T_bath_min,
            adaptive_mpc_T_bath_max=adaptive_mpc_T_bath_max,
            adaptive_mpc_delta_lambda_max=adaptive_mpc_delta_lambda_max,
            adaptive_mpc_delta_T_bath_max=adaptive_mpc_delta_T_bath_max,
            adaptive_mpc_apply_to=adaptive_mpc_apply_to,
            adaptive_mpc_rls_forgetting_factor=adaptive_mpc_rls_forgetting_factor,
            adaptive_mpc_rls_initial_covariance=adaptive_mpc_rls_initial_covariance,
            adaptive_mpc_model_update_interval=adaptive_mpc_model_update_interval,
            adaptive_mpc_output_file=adaptive_mpc_output_file,
            adaptive_mpc_console_output_period_ps=adaptive_mpc_console_output_period_ps,
            # BathPI controller parameters
            enable_bath_pi_controller=enable_bath_pi_controller,
            bath_pi_apply_to=bath_pi_apply_to,
            bath_pi_K_p_molecular=bath_pi_K_p_molecular,
            bath_pi_K_i_molecular=bath_pi_K_i_molecular,
            bath_pi_K_T_molecular=bath_pi_K_T_molecular,
            bath_pi_K_p_cavity=bath_pi_K_p_cavity,
            bath_pi_K_i_cavity=bath_pi_K_i_cavity,
            bath_pi_K_T_cavity=bath_pi_K_T_cavity,
            bath_pi_filter_window_ps=bath_pi_filter_window_ps,
            bath_pi_flux_source=bath_pi_flux_source,
            bath_pi_anti_windup_alpha=bath_pi_anti_windup_alpha,
            bath_pi_enable_feedforward=bath_pi_enable_feedforward,
            bath_pi_T_nominal=bath_pi_T_nominal,
            bath_pi_feedforward_tau_ps=bath_pi_feedforward_tau_ps,
            bath_pi_turn_on_time_ps=bath_pi_turn_on_time_ps,
            bath_pi_turn_off_time_ps=bath_pi_turn_off_time_ps,
            bath_pi_update_interval_ps=bath_pi_update_interval_ps,
            bath_pi_T_min=bath_pi_T_min,
            bath_pi_T_max=bath_pi_T_max,
            bath_pi_rate_limit_K_per_ps=bath_pi_rate_limit_K_per_ps,
            bath_pi_output_file=bath_pi_output_file,
            bath_pi_relaxation_data_file=bath_pi_relaxation_data_file,
            # LQR controller parameters
            enable_lqr_controller=enable_lqr_controller,
            lqr_signal_method=lqr_signal_method,
            lqr_hot_method=lqr_hot_method,
            lqr_target_temperature=lqr_target_temperature,
            lqr_dynamic_target=lqr_dynamic_target,
            lqr_dynamic_target_method=lqr_dynamic_target_method,
            lqr_weight_signal=lqr_weight_signal,
            lqr_weight_hot=lqr_weight_hot,
            lqr_weight_bath=lqr_weight_bath,  # Deprecated: not used in LQG 2D formulation
            lqr_weight_integral=lqr_weight_integral,
            lqr_control_effort=lqr_control_effort,
            lqr_process_noise_signal=lqr_process_noise_signal,
            lqr_process_noise_hot=lqr_process_noise_hot,
            lqr_measurement_noise_signal=lqr_measurement_noise_signal,
            lqr_measurement_noise_hot=lqr_measurement_noise_hot,
            lqr_system_id_mode=lqr_system_id_mode,
            lqr_system_id_temp_K=lqr_system_id_temp_K,
            lqr_system_id_duration_ps=lqr_system_id_duration_ps,
            lqr_system_id_file=lqr_system_id_file,
            lqr_periodic_system_id=lqr_periodic_system_id,
            lqr_periodic_system_id_interval_ps=lqr_periodic_system_id_interval_ps,
            # EKF adaptation
            lqr_use_ekf_adaptation=lqr_use_ekf_adaptation,
            lqr_ekf_update_interval=lqr_ekf_update_interval,
            lqr_ekf_process_noise_param=lqr_ekf_process_noise_param,
            lqr_ekf_initial_covariance_param=lqr_ekf_initial_covariance_param,
            lqr_adaptive_lqr_threshold=lqr_adaptive_lqr_threshold,
            # Gain scheduling
            lqr_enable_gain_scheduling=lqr_enable_gain_scheduling,
            lqr_gain_schedule_far_threshold=lqr_gain_schedule_far_threshold,
            lqr_gain_schedule_near_threshold=lqr_gain_schedule_near_threshold,
            # T_h low-pass filter
            lqr_th_filter_enabled=lqr_th_filter_enabled,
            lqr_th_filter_time_constant=lqr_th_filter_time_constant,
            # Gentle startup
            lqr_gentle_startup_steps=lqr_gentle_startup_steps,
            lqr_gentle_startup_min_authority=lqr_gentle_startup_min_authority,
            # Kinetic temperature tracking (3D state)
            lqr_track_kinetic_temp=lqr_track_kinetic_temp,
            lqr_weight_kinetic=lqr_weight_kinetic,
            lqr_process_noise_kinetic=lqr_process_noise_kinetic,
            lqr_measurement_noise_kinetic=lqr_measurement_noise_kinetic,
            # Cross-coupling weights for thermal equilibration
            lqr_cross_coupling_signal_kinetic=lqr_cross_coupling_signal_kinetic,
            lqr_cross_coupling_signal_hot=lqr_cross_coupling_signal_hot,
            lqr_cross_coupling_hot_kinetic=lqr_cross_coupling_hot_kinetic,
            # Timing
            lqr_turn_on_time_ps=lqr_turn_on_time_ps,
            lqr_turn_off_time_ps=lqr_turn_off_time_ps,
            lqr_update_interval_ps=lqr_update_interval_ps,
            lqr_T_min=lqr_T_min,
            lqr_T_max=lqr_T_max,
            lqr_apply_to=lqr_apply_to,
            lqr_output_file=lqr_output_file,
            lqr_empirical_data_file=lqr_empirical_data_file,
            # Adaptive LQI controller parameters
            lqr_controller_type=lqr_controller_type,
            lqr_tau_L_initial=lqr_tau_L_initial,
            lqr_tau_H_initial=lqr_tau_H_initial,
            lqr_k_initial=lqr_k_initial,
            lqr_tau_b=lqr_tau_b,
            lqr_q_common=lqr_q_common,
            lqr_q_diff=lqr_q_diff,
            lqr_q_eta_common=lqr_q_eta_common,
            lqr_q_eta_diff=lqr_q_eta_diff,
            lqr_process_noise_drift=lqr_process_noise_drift,
            lqr_rls_forgetting=lqr_rls_forgetting,
            lqr_rls_regularization=lqr_rls_regularization,
            lqr_rls_update_interval=lqr_rls_update_interval,
            lqr_max_control_rate=lqr_max_control_rate,
            lqr_integral_max_common=lqr_integral_max_common,
            lqr_integral_max_diff=lqr_integral_max_diff,
            lqr_theta_change_threshold=lqr_theta_change_threshold,
            # LQG pulse controller parameters
            enable_lqg_controller=enable_lqg_controller,
            lqg_A_matrix=lqg_A_matrix,
            lqg_B_matrix=lqg_B_matrix,
            lqg_C_matrix=lqg_C_matrix,
            lqg_S_matrix=lqg_S_matrix,
            lqg_r_target=lqg_r_target,
            lqg_Qz_weights=lqg_Qz_weights,
            lqg_Ru_weights=lqg_Ru_weights,
            lqg_Qe_weights=lqg_Qe_weights,
            lqg_W_noise=lqg_W_noise,
            lqg_V_noise=lqg_V_noise,
            lqg_g0_baseline=lqg_g0_baseline,
            lqg_dt=lqg_dt,
            lqg_delta_g_min=lqg_delta_g_min,
            lqg_delta_g_max=lqg_delta_g_max,
            lqg_T_bath_min=lqg_T_bath_min,
            lqg_T_bath_max=lqg_T_bath_max,
            lqg_rate_limit_K_per_ps=lqg_rate_limit_K_per_ps,
            lqg_turn_on_time_ps=lqg_turn_on_time_ps,
            lqg_turn_off_time_ps=lqg_turn_off_time_ps,
            lqg_update_interval_ps=lqg_update_interval_ps,
            lqg_console_output_period_ps=lqg_console_output_period_ps,
            lqg_c_affine_term=lqg_c_affine_term,
            lqg_temperature_methods=lqg_temperature_methods,
            lqg_square_period_ps=lqg_square_period_ps,
            lqg_square_duty_cycle=lqg_square_duty_cycle,
            lqg_square_phase_offset=lqg_square_phase_offset,
            lqg_square_start_time_ps=lqg_square_start_time_ps,
            lqg_square_stop_time_ps=lqg_square_stop_time_ps,
            lqg_square_g_low_level=lqg_square_g_low_level,
            # LQG system identification parameters
            lqg_enable_system_id=lqg_enable_system_id,
            lqg_equilibrium_duration_ps=lqg_equilibrium_duration_ps,
            lqg_prbs_n_unique_steps=lqg_prbs_n_unique_steps,
            lqg_prbs_amplitude=lqg_prbs_amplitude,
            lqg_bath_n_unique_steps=lqg_bath_n_unique_steps,
            lqg_bath_excitation_amplitude=lqg_bath_excitation_amplitude,
            lqg_system_id_output_file=lqg_system_id_output_file,
            # PID controller parameters
            enable_pid_controller=enable_pid_controller,
            pid_signal_choice=pid_signal_choice,
            pid_target_temperature=pid_target_temperature,
            pid_self_loop=pid_self_loop,
            pid_Kp=pid_Kp,
            pid_Ti=pid_Ti,
            pid_Td=pid_Td,
            pid_auto_tune=pid_auto_tune,
            pid_auto_tune_step_size=pid_auto_tune_step_size,
            pid_auto_tune_duration_ps=pid_auto_tune_duration_ps,
            pid_turn_on_time_ps=pid_turn_on_time_ps,
            pid_turn_off_time_ps=pid_turn_off_time_ps,
            pid_update_interval_ps=pid_update_interval_ps,
            pid_apply_to=pid_apply_to,
            pid_T_min=pid_T_min,
            pid_T_max=pid_T_max,
            pid_rate_limit_K_per_ps=pid_rate_limit_K_per_ps,
            pid_derivative_filter_N=pid_derivative_filter_N,
            pid_enable_anti_windup=pid_enable_anti_windup,
            pid_output_file=pid_output_file,
            pid_console_output_period_ps=pid_console_output_period_ps,
            # Temperature tracker parameters
            enable_temp_tracker=enable_temp_tracker,
            temp_tracker_output_period_ps=temp_tracker_output_period_ps,
            temp_tracker_empirical_data_file=temp_tracker_empirical_data_file,
            # Molecular temperature decomposition parameters
            enable_molecular_temps=enable_molecular_temps,
            molecular_temps_output_period_ps=molecular_temps_output_period_ps,
            # Dipole moment FDR parameters
            enable_dipole_fdr=enable_dipole_fdr,
            dipole_fdr_output_period_ps=dipole_fdr_output_period_ps,
            dipole_fdr_max_correlation_time_ps=dipole_fdr_max_correlation_time_ps,
            dipole_fdr_field_direction=dipole_fdr_field_direction,
            dipole_fdr_exclude_cavity=dipole_fdr_exclude_cavity,
            enable_dipole_response=enable_dipole_response,
            dipole_response_field_strength=dipole_response_field_strength,
            dipole_response_sign=dipole_response_sign,
            # Mechanical cavity modulation parameters  
            mech_periodic=mech_periodic,
            mech_frequency_cm1=mech_frequency_cm1,
            mech_magnitude=mech_magnitude,
            # Switch-time parameters
            switch_time_ps=switch_time_ps,
            decay_time_constant_ps=decay_time_constant_ps,
            # Thermostat parameters
            molecular_thermostat_tau=molecular_tau,
            cavity_thermostat_tau=cavity_tau,
            # Timestep parameters
            dt_fs=timestep_fs,
            error_tolerance=error_tolerance,
            initial_fraction=initial_fraction,
            time_constant_ps=time_constant_ps,
            # Dynamic coupling detection parameters
            enable_dynamic_coupling_detection=enable_dynamic_coupling_detection,
            coupling_change_threshold=coupling_change_threshold,
            # Output and tracking parameters
            enable_energy_tracking=enable_energy_tracker,
            energy_output_period_ps=energy_output_period_ps,
            fkt_output_period_ps=fkt_output_period_ps,
            gsd_output_period_ps=gsd_output_period_ps,
            console_output_period_ps=console_output_period_ps,
            enable_fkt=enable_fkt,
            fkt_kmag=fkt_kmag,
            fkt_num_wavevectors=fkt_wavevectors,
            fkt_reference_interval_ps=fkt_ref_interval,
            fkt_max_references=fkt_max_refs,
            enable_dipole_autocorr=enable_dipole_autocorr,
            dipole_reference_interval_ps=dipole_ref_interval,
            dipole_max_references=dipole_max_refs,
            dipole_output_period_ps=dipole_output_period_ps,
            max_energy_output_time_ps=max_energy_output_time,
            # Hardware parameters
            device=device,
            gpu_id=gpu_id,
            truncate_gsd=truncate_gsd,
            seed=seed,
            restart_velocities=restart_velocities,
            zero_momentum_enabled=zero_momentum_enabled,
            zero_momentum_period_ps=zero_momentum_period_ps,
            # Empirical feedback parameters
            enable_empirical_feedback=enable_empirical_feedback,
            empirical_data_file=empirical_data_file,
            feedback_output_period_ps=feedback_output_period_ps,
            feedback_energy_component=feedback_energy_component,
            feedback_apply_to=feedback_apply_to,
            feedback_averaging_window_ps=feedback_averaging_window_ps,
            feedback_update_interval_ps=feedback_update_interval_ps,
            feedback_T_min=feedback_T_min,
            feedback_T_max=feedback_T_max,
            feedback_turn_off_time_ps=feedback_turn_off_time_ps,
            feedback_use_direct_harmonic_calculation=feedback_use_direct_harmonic,
            feedback_auto_adjust_molecular_tau=feedback_auto_adjust_molecular_tau,
            enable_auto_stopping=enable_auto_stopping,
            auto_stop_min_time_ps=auto_stop_min_time_ps,
            auto_stop_smoothing_window=auto_stop_smoothing_window
        )
        
        # Harmonic reset parameters are passed to CavityMDSimulation
        # The updaters will be created inside the simulation class
        
        return sim.run() == 0  # Return True for success (exit code 0)
        
    except Exception as e:
        print(f" Error running experiment: {e}")
        import traceback
        traceback.print_exc()
        return False

# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================

def create_parser():
    """Create command line argument parser."""
    parser = argparse.ArgumentParser(
        description='Advanced Unified Cavity Dynamics MD Experiment Runner',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Periodic coupling with mechanical modulation
  python 18_unified_cavity_dynamics.py --periodic --coupling 1e-3 --period 1.0 --mech-periodic --mech-frequency 100.0 --mech-magnitude 1e-4
  
  # Laser drive with mechanical modulation  
  python 18_unified_cavity_dynamics.py --laser --laser-frequency 1560 --laser-amplitude 1e-5 --mech-periodic --mech-frequency 100.0 --mech-magnitude 1e-4
  
  # Switch-time coupling with feedback
  python 18_unified_cavity_dynamics.py --coupling 1e-3 --switch-time 10.0 --enable-empirical-feedback --empirical-data-file data.txt
  
  # Pure mechanical modulation (no coupling oscillation)
  python 18_unified_cavity_dynamics.py --coupling 1e-3 --mech-periodic --mech-frequency 100.0 --mech-magnitude 1e-4

 ENHANCED COUPLING VARIANTS:
  --coupling-type constant       : g(t) = constant (default)
  --coupling-type step           : g(t) = 0 → g_target at t_switch (with optional decay/turn-off)
  --coupling-type periodic       : g(t) = A*sin(2π*t/T + φ) (with start/stop times)
  --coupling-type exponential    : g(t) = A*exp(-t/τ) (with turn-on/off times)
  --coupling-type square         : g(t) = periodic on/off pulses (with duty cycle)
  --coupling-type decaying_square: g(t) = A_n*square_wave(t), A_n = A₀*(1-r)ⁿ (amplitude decays per period)
  --coupling-type adaptive_square: g(t) = g_next*square_wave(t), g_next = g_target*√(T_target/T_harmonic) (adaptive amplitude control)
  --coupling-type exponentialwave: g(t) = A*exp(-t_period/τ) (periodic exponential decay pulses)
  --coupling-type lqg_square     : g(t) = LQG-controlled square wave (pulse-to-pulse optimal control)

 ENHANCED PERIODIC COUPLING PARAMETERS:
  --periodic-period          : Period in ps (default: 1.0)
  --periodic-phase-offset    : Phase offset in radians (default: 0.0)
  --periodic-start-time      : Start time in ps (default: 0.0)
  --periodic-stop-time       : Stop time in ps (default: None)

 ENHANCED DECAYING SQUARE WAVE PARAMETERS:
  --decaying-square-period           : Period in ps (default: 2.0)
  --decaying-square-duty-cycle       : Duty cycle 0.0-1.0 (default: 0.5)
  --decaying-square-phase-offset     : Phase offset in radians (default: 0.0)
  --decaying-square-start-time       : Start time in ps (default: 0.0)
  --decaying-square-stop-time        : Stop time in ps (default: None)
  --decaying-square-decay-rate       : Decay rate per period 0.0-1.0 (default: 0.1 = 10%)
  --decaying-square-minimum-amplitude: Minimum amplitude threshold (default: 1e-6)

 ADAPTIVE SQUARE WAVE PARAMETERS:
  --adaptive-square-period          : Period in ps (default: 5.0)
  --adaptive-square-duty-cycle      : Duty cycle 0.0-1.0 (default: 0.5)  
  --adaptive-square-phase-offset    : Phase offset in radians (default: 0.0)
  --adaptive-square-start-time      : Start time in ps (default: 0.0)
  --adaptive-square-stop-time       : Stop time in ps (default: None)
  --adaptive-square-min-amplitude   : Minimum amplitude limit (default: 1e-8)
  --adaptive-square-max-amplitude   : Maximum amplitude limit (default: 1e-1)

 EXPONENTIAL WAVE PARAMETERS:
  --exp-period                      : Period in ps (default: 2.0)
  --exp-tau                         : Exponential decay time constant in ps (default: 0.5)
  --exp-start-time                  : Start time in ps (default: 0.0)
  --exp-stop-time                   : Stop time in ps (default: None)
  --exp-adaptive                    : Enable adaptive amplitude scaling based on T_bath/T_harmonic ratio

 AUTO-STOP COUPLING PARAMETERS:
  --enable-auto-stop                : Enable automatic coupling stop when temperatures converge
  --auto-stop-tol                   : Temperature tolerance for convergence in K (default: 1.0)
  --auto-stop-window                : Averaging window for convergence check in ps (default: 10.0)

 ENHANCED LASER DRIVE:
  --laser                    : Enable laser with timing control F_laser = F_L*cos(ω_L*t)
  --laser-start-time         : Laser start time in ps
  --laser-stop-time          : Laser stop time in ps

 QUENCH CONTROLLER:
  --enable-quench-controller : Enable instantaneous temperature quench
  --quench-initial-temperature : Initial temperature before quench in K (default: 100.0)
  --quench-target-temperature  : Target temperature after quench in K (default: 50.0)
  --quench-time               : Time when quench occurs in ps (default: 50.0)
  --quench-apply-to           : Apply to 'molecular', 'cavity', or 'both' (default: both)

 COMPREHENSIVE TEMPERATURE TRACKER:
  --enable-temp-tracker      : Track all temperatures (kinetic, fictive, bath)
  --temp-tracker-empirical-data-file : Empirical data for LJ+Coul fictive temp

 MOLECULAR TEMPERATURE DECOMPOSITION:
  --enable-molecular-temps   : Decompose kinetic temperature into trans/rot/vib for dimers
  --molecular-temps-output-period-ps : Output period in ps (default: 1.0)

 DIPOLE MOMENT FDR ANALYSIS:
  --enable-dipole-fdr        : Enable dipole moment autocorrelation tracking
  --enable-dipole-response   : Enable electric field response for fork-and-clone FDR
  --dipole-response-sign +1/-1 : Sign for plus/minus clone in FDR measurements

 LEGACY MODES (still supported):
  --periodic             : Legacy periodic coupling
  --switch-time          : Legacy step coupling
  --mech-periodic        : Mechanical cavity modulation L(t) = L₀ + A*sin(ω*t)
        """
    )
    
    # Basic simulation parameters
    parser.add_argument('--molecular-bath', choices=['bussi', 'langevin', 'none'], default='bussi',
                       help='Molecular thermostat type (default: bussi)')
    parser.add_argument('--cavity-bath', choices=['bussi', 'langevin', 'none'], default='langevin',
                       help='Cavity thermostat type (default: langevin)')
    parser.add_argument('--finite-q', action='store_true',
                       help='Use finite-q cavity mode (default: q=0 mode)')
    parser.add_argument('--coupling', type=float, default=None,
                       help='[DEPRECATED] Cavity coupling strength epsilon in a.u. (default: None, use --lambda-coupling)')
    parser.add_argument('--lambda-coupling', type=float, default=1e-3,
                       help='Dimensionless coupling parameter lambda (epsilon = lambda * omega_c) (default: 1e-3)')
    
    # Enhanced coupling variant selection
    parser.add_argument('--coupling-type', choices=['constant', 'step', 'periodic', 'exponential', 'square', 'decaying_square', 'adaptive_square', 'exponentialwave', 'composite', 'lqg_square', 'lqg_coupling'], 
                       default='constant', help='Coupling variant type (default: constant)')
    
    # Enhanced step coupling parameters 
    parser.add_argument('--step-turn-off-time', type=float, default=None,
                       help='Time in ps to turn off step coupling (default: None, no turn-off)')
    
    # Exponential decay coupling parameters
    parser.add_argument('--exponential-amplitude', type=float, default=None,
                       help='Exponential decay amplitude (default: None, uses --coupling value)')
    parser.add_argument('--exponential-decay-time', type=float, default=10.0,
                       help='Exponential decay time constant in ps (default: 10.0)')
    parser.add_argument('--exponential-turn-on-time', type=float, default=0.0,
                       help='Exponential turn-on time in ps (default: 0.0)')
    parser.add_argument('--exponential-turn-off-time', type=float, default=None,
                       help='Exponential turn-off time in ps (default: None, no turn-off)')
    
    # Square wave coupling parameters
    parser.add_argument('--square-period', type=float, default=2.0,
                       help='Square wave period in ps (default: 2.0)')
    parser.add_argument('--square-duty-cycle', type=float, default=0.5,
                       help='Square wave duty cycle (0.0-1.0) (default: 0.5)')
    parser.add_argument('--square-phase-offset', type=float, default=0.0,
                       help='Square wave phase offset in radians (default: 0.0)')
    parser.add_argument('--square-start-time', type=float, default=0.0,
                       help='Square wave start time in ps (default: 0.0)')
    parser.add_argument('--square-stop-time', type=float, default=None,
                       help='Square wave stop time in ps (default: None, no stop)')
    
    # Decaying square wave coupling parameters
    parser.add_argument('--decaying-square-period', type=float, default=2.0,
                       help='Decaying square wave period in ps (default: 2.0)')
    parser.add_argument('--decaying-square-duty-cycle', type=float, default=0.5,
                       help='Decaying square wave duty cycle (0.0-1.0) (default: 0.5)')
    parser.add_argument('--decaying-square-phase-offset', type=float, default=0.0,
                       help='Decaying square wave phase offset in radians (default: 0.0)')
    parser.add_argument('--decaying-square-start-time', type=float, default=0.0,
                       help='Decaying square wave start time in ps (default: 0.0)')
    parser.add_argument('--decaying-square-stop-time', type=float, default=None,
                       help='Decaying square wave stop time in ps (default: None, no stop)')
    parser.add_argument('--decaying-square-decay-rate', type=float, default=0.1,
                       help='Decay rate per period (0.0-1.0) (default: 0.1 = 10%% per period)')
    parser.add_argument('--decaying-square-minimum-amplitude', type=float, default=1e-6,
                       help='Minimum amplitude below which wave stops (default: 1e-6)')
    
    # Adaptive square wave coupling parameters
    parser.add_argument('--adaptive-square-period', type=float, default=5.0,
                       help='Adaptive square wave period in ps (default: 5.0)')
    parser.add_argument('--adaptive-square-duty-cycle', type=float, default=0.5,
                       help='Adaptive square wave duty cycle (0.0-1.0) (default: 0.5)')
    parser.add_argument('--adaptive-square-phase-offset', type=float, default=0.0,
                       help='Adaptive square wave phase offset in radians (default: 0.0)')
    parser.add_argument('--adaptive-square-start-time', type=float, default=0.0,
                       help='Adaptive square wave start time in ps (default: 0.0)')
    parser.add_argument('--adaptive-square-stop-time', type=float, default=None,
                       help='Adaptive square wave stop time in ps (default: None, no stop)')
    parser.add_argument('--adaptive-square-min-amplitude', type=float, default=1e-8,
                       help='Minimum allowed amplitude (safety limit) (default: 1e-8)')
    parser.add_argument('--adaptive-square-max-amplitude', type=float, default=1e-1,
                       help='Maximum allowed amplitude (safety limit) (default: 1e-1)')
    
    # Exponential wave coupling parameters
    parser.add_argument('--exp-period', type=float, default=2.0,
                       help='Exponential wave period in ps (default: 2.0)')
    parser.add_argument('--exp-tau', type=float, default=0.5,
                       help='Exponential decay time constant in ps (default: 0.5)')
    parser.add_argument('--exp-start-time', type=float, default=0.0,
                       help='Exponential wave start time in ps (default: 0.0)')
    parser.add_argument('--exp-stop-time', type=float, default=None,
                       help='Exponential wave stop time in ps (default: None, no stop)')
    parser.add_argument('--exp-adaptive', action='store_true',
                       help='Enable adaptive exponential wave coupling based on T_bath/T_harmonic scaling')
    
    # Composite coupling parameters
    parser.add_argument('--composite-sinusoid-amplitude', type=float, default=1e-4,
                       help='Amplitude for sinusoidal component (default: 1e-4)')
    parser.add_argument('--composite-sinusoid-period', type=float, default=1.0,
                       help='Period for sinusoidal component in ps (default: 1.0)')
    parser.add_argument('--composite-sinusoid-phase', type=float, default=0.0,
                       help='Phase offset for sinusoidal component in radians (default: 0.0)')
    parser.add_argument('--composite-sinusoid-start-time', type=float, default=0.0,
                       help='Start time for sinusoidal component in ps (default: 0.0)')
    parser.add_argument('--composite-sinusoid-stop-time', type=float, default=None,
                       help='Stop time for sinusoidal component in ps (default: None, runs indefinitely)')
    parser.add_argument('--composite-square-amplitude', type=float, default=2e-4,
                       help='Amplitude for adaptive square component (default: 2e-4)')
    parser.add_argument('--composite-square-period', type=float, default=50.0,
                       help='Period for adaptive square component in ps (default: 50.0)')
    parser.add_argument('--composite-square-duty-cycle', type=float, default=0.02,
                       help='Duty cycle for adaptive square component (default: 0.02)')
    parser.add_argument('--composite-square-start-time', type=float, default=10.0,
                       help='Start time for adaptive square component in ps (default: 10.0)')
    parser.add_argument('--composite-square-stop-time', type=float, default=1000.0,
                       help='Stop time for square wave component in ps (default: 1000.0)')
    parser.add_argument('--composite-square-adaptive', action='store_true',
                       help='Use adaptive square wave (amplitude adjusts to temperature) instead of fixed amplitude')
    parser.add_argument('--composite-max-amplitude', type=float, default=None,
                       help='Maximum total amplitude for composite (default: None, no limit)')
    
    # Auto-stop coupling parameters
    parser.add_argument('--enable-auto-stop', action='store_true',
                       help='Enable automatic coupling stop when temperatures converge (default: False)')
    parser.add_argument('--auto-stop-tol', type=float, default=1.0,
                       help='Temperature tolerance for auto-stop in K (default: 1.0)')
    parser.add_argument('--auto-stop-window', type=float, default=10.0,
                       help='Averaging window for auto-stop in ps (default: 10.0)')
    
    # Enhanced periodic coupling parameters 
    parser.add_argument('--periodic-period', type=float, default=1.0,
                       help='Periodic coupling period in ps (default: 1.0)')
    parser.add_argument('--periodic-phase-offset', type=float, default=0.0,
                       help='Periodic coupling phase offset in radians (default: 0.0)')
    parser.add_argument('--periodic-start-time', type=float, default=0.0,
                       help='Time to start periodic oscillation in ps (default: 0.0)')
    parser.add_argument('--periodic-stop-time', type=float, default=None,
                       help='Time to stop periodic oscillation in ps (default: None, no stop)')
    
    # Enhanced laser drive parameters
    parser.add_argument('--laser-start-time', type=float, default=0.0,
                       help='Time to start laser drive in ps (default: 0.0)')
    parser.add_argument('--laser-stop-time', type=float, default=None,
                       help='Time to stop laser drive in ps (default: None, runs indefinitely)')
    parser.add_argument('--laser-kvector', nargs=3, type=float, default=None,
                       help='Laser k-vector components [kx, ky, kz] (default: [1.0, 0.0, 0.0])')
    
    # PI feedback controller parameters
    
    # Gradient descent feedback controller parameters
    parser.add_argument('--enable-gd-feedback', action='store_true',
                       help='Enable gradient descent feedback temperature controller')
    parser.add_argument('--gd-target-temperature', type=float, default=100.0,
                       help='GD controller target temperature in K (default: 100.0)')
    parser.add_argument('--gd-turn-on-time', type=float, default=0.0,
                       help='GD controller turn-on time in ps (default: 0.0)')
    parser.add_argument('--gd-turn-off-time', type=float, default=None,
                       help='GD controller turn-off time in ps (default: None, never turn off)')
    parser.add_argument('--gd-method', type=str, default='kinetic',
                       choices=['kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition', 'lj_coulombic_kinetic', 'lj_coulombic_bath'],
                       help='GD controller temperature method (default: kinetic)')
    parser.add_argument('--gd-update-interval', type=float, default=0.1,
                       help='GD controller update interval in ps (default: 0.1)')
    parser.add_argument('--gd-time-constant', type=float, default=10.0,
                       help='GD controller time constant in ps (default: 10.0)')
    parser.add_argument('--gd-apply-to', type=str, default='both', 
                       choices=['molecular', 'cavity', 'both'],
                       help='Which thermostats to control: molecular, cavity, or both (default: both)')
    parser.add_argument('--gd-T-min', type=float, default=0.0,
                       help='GD controller minimum temperature in K (default: 0.0)')
    parser.add_argument('--gd-T-max', type=float, default=None,
                       help='GD controller maximum temperature in K (default: None)')
    parser.add_argument('--gd-no-effective-temp', action='store_true',
                       help='Use raw measured temperature instead of effective temperature mixing (default: False)')
    
    # Multi-signal error function parameters
    parser.add_argument('--gd-enable-multi-signal', action='store_true',
                       help='Enable multi-signal error function combining three error signals (default: False)')
    parser.add_argument('--gd-weight-system-target', type=float, default=1.0,
                       help='Weight for (Ts - Ttarget) error signal (default: 1.0)')
    parser.add_argument('--gd-weight-bath-target', type=float, default=0.0,
                       help='Weight for (Tbath - Ttarget) error signal (default: 0.0)')  
    parser.add_argument('--gd-weight-system-bath', type=float, default=0.0,
                       help='Weight for (Ts - Tbath) error signal (default: 0.0)')
    
    # Dual independent feedback controller parameters
    parser.add_argument('--enable-dual-feedback', action='store_true',
                       help='Enable dual independent feedback controller (default: False)')
    parser.add_argument('--dual-cavity-method', type=str, default='harmonic_equipartition',
                       choices=['kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition', 'lj_coulombic_kinetic', 'lj_coulombic_bath'],
                       help='Temperature method for cavity bath control (default: harmonic_equipartition)')
    parser.add_argument('--dual-molecular-method', type=str, default='lj_coulombic_kinetic',
                       choices=['kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition', 'lj_coulombic_kinetic', 'lj_coulombic_bath'],
                       help='Temperature method for molecular bath control (default: lj_coulombic_kinetic)')
    parser.add_argument('--dual-cavity-target', type=float, default=100.0,
                       help='Target temperature for cavity bath in K (default: 100.0)')
    parser.add_argument('--dual-molecular-target', type=float, default=100.0,
                       help='Target temperature for molecular bath in K (default: 100.0)')
    parser.add_argument('--dual-cavity-time-constant', type=float, default=5.0,
                       help='Time constant for cavity bath control in ps (default: 5.0)')
    parser.add_argument('--dual-molecular-time-constant', type=float, default=10.0,
                       help='Time constant for molecular bath control in ps (default: 10.0)')
    parser.add_argument('--dual-turn-on-time', type=float, default=0.0,
                       help='Turn-on time for dual controller in ps (default: 0.0)')
    parser.add_argument('--dual-turn-off-time', type=float, default=None,
                       help='Turn-off time for dual controller in ps (default: None, never turn off)')
    parser.add_argument('--dual-update-interval', type=float, default=0.1,
                       help='Update interval for dual controller in ps (default: 0.1)')
    parser.add_argument('--dual-cavity-T-min', type=float, default=0.0,
                       help='Minimum cavity bath temperature in K (default: 0.0)')
    parser.add_argument('--dual-cavity-T-max', type=float, default=None,
                       help='Maximum cavity bath temperature in K (default: None)')
    parser.add_argument('--dual-molecular-T-min', type=float, default=0.0,
                       help='Minimum molecular bath temperature in K (default: 0.0)')
    parser.add_argument('--dual-molecular-T-max', type=float, default=None,
                       help='Maximum molecular bath temperature in K (default: None)')
    parser.add_argument('--dual-cavity-dynamic-target', action='store_true',
                       help='Use dynamic target temperature for cavity (measured at turn-on time) instead of fixed target (default: False)')
    parser.add_argument('--dual-molecular-dynamic-target', action='store_true',
                       help='Use dynamic target temperature for molecular (measured at turn-on time) instead of fixed target (default: False)')
    parser.add_argument('--dual-cavity-integral-time-constant', type=float, default=None,
                       help='Integral time constant for cavity PI controller in ps (default: None, P-only). Smaller values = more aggressive integral action')
    parser.add_argument('--dual-molecular-integral-time-constant', type=float, default=None,
                       help='Integral time constant for molecular PI controller in ps (default: None, P-only). Smaller values = more aggressive integral action')
    
    # Asymmetric gain factors
    parser.add_argument('--dual-cavity-heating-gain-factor', type=float, default=1.0,
                       help='Proportional gain factor for cavity heating (error < 0) (default: 1.0)')
    parser.add_argument('--dual-cavity-cooling-gain-factor', type=float, default=1.0,
                       help='Proportional gain factor for cavity cooling (error > 0) (default: 1.0)')
    parser.add_argument('--dual-molecular-heating-gain-factor', type=float, default=1.0,
                       help='Proportional gain factor for molecular heating (error < 0) (default: 1.0)')
    parser.add_argument('--dual-molecular-cooling-gain-factor', type=float, default=1.0,
                       help='Proportional gain factor for molecular cooling (error > 0) (default: 1.0)')
    
    # Asymmetric integral time constants
    parser.add_argument('--dual-cavity-integral-heating-time-constant', type=float, default=None,
                       help='Integral time constant for cavity heating in ps (default: None, use base integral constant)')
    parser.add_argument('--dual-cavity-integral-cooling-time-constant', type=float, default=None,
                       help='Integral time constant for cavity cooling in ps (default: None, use base integral constant)')
    parser.add_argument('--dual-molecular-integral-heating-time-constant', type=float, default=None,
                       help='Integral time constant for molecular heating in ps (default: None, use base integral constant)')
    parser.add_argument('--dual-molecular-integral-cooling-time-constant', type=float, default=None,
                       help='Integral time constant for molecular cooling in ps (default: None, use base integral constant)')
    
    # Sinusoidal bath temperature controller parameters
    parser.add_argument('--enable-sinusoidal-bath', action='store_true',
                       help='Enable sinusoidal bath temperature controller (default: False)')
    parser.add_argument('--sinusoidal-bath-period', type=float, default=1.0,
                       help='Period of sinusoidal temperature oscillation in ps (default: 1.0)')
    parser.add_argument('--sinusoidal-bath-amplitude-scale', type=float, default=0.1,
                       help='Scaling factor for amplitude based on harmonic equipartition deviation (default: 0.1)')
    parser.add_argument('--sinusoidal-bath-phase-offset', type=float, default=0.0,
                       help='Phase offset for sinusoidal function in radians (default: 0.0)')
    parser.add_argument('--sinusoidal-bath-target', type=float, default=100.0,
                       help='Static target temperature in K (ignored if dynamic target enabled) (default: 100.0)')
    parser.add_argument('--sinusoidal-bath-dynamic-target', action='store_true',
                       help='Set mean temperature dynamically at turn-on time (default: False)')
    parser.add_argument('--sinusoidal-bath-turn-on-time', type=float, default=0.0,
                       help='Time to start control in ps (default: 0.0)')
    parser.add_argument('--sinusoidal-bath-turn-off-time', type=float, default=None,
                       help='Time to stop control in ps (default: None = never turn off)')
    parser.add_argument('--sinusoidal-bath-update-interval', type=float, default=0.1,
                       help='Interval between control updates in ps (default: 0.1)')
    parser.add_argument('--sinusoidal-bath-apply-to', type=str, default='both',
                       choices=['molecular', 'cavity', 'both'],
                       help='Which thermostats to control (default: both)')
    parser.add_argument('--sinusoidal-bath-T-min', type=float, default=0.1,
                       help='Minimum allowed bath temperature in K (default: 0.1)')
    parser.add_argument('--sinusoidal-bath-T-max', type=float, default=None,
                       help='Maximum allowed bath temperature in K (default: None)')
    parser.add_argument('--sinusoidal-bath-empirical-data-file', type=str, default=None,
                       help='Path to empirical data file for harmonic equipartition calculation')
    parser.add_argument('--sinusoidal-bath-amplitude-update-interval', type=float, default=1.0,
                       help='Interval for updating amplitude based on harmonic temperature in ps (default: 1.0)')
    parser.add_argument('--sinusoidal-bath-amplitude-temperature-method', type=str, default='harmonic_equipartition',
                       choices=['harmonic_equipartition', 'kinetic', 'lj_coulombic', 'harmonic'],
                       help='Temperature method for amplitude calculation (default: harmonic_equipartition)')
    parser.add_argument('--sinusoidal-bath-adaptive-range-mode', action='store_true',
                       help='Enable adaptive range mode: sinusoid bounds fixed at target, range adapts to |Ts-Ttarget| (default: False)')
    
    # Harmonic bond reset parameters
    parser.add_argument('--enable-harmonic-reset', action='store_true',
                       help='Enable one-time harmonic bond reset (thermal sampling of bond stretch DOF for ALL bond types)')
    parser.add_argument('--harmonic-reset-turn-on-time', type=float, default=0.0,
                       help='Time in ps to trigger harmonic bond reset (default: 0.0)')
    parser.add_argument('--harmonic-reset-bond-type', type=str, default='bond',
                       help='[DEPRECATED] Bond type name - now auto-detects all bond types')
    parser.add_argument('--harmonic-reset-spring-constant', type=float, default=None,
                       help='[DEPRECATED] Spring constant K - now auto-detected from force field')
    parser.add_argument('--harmonic-reset-equilibrium-length', type=float, default=None,
                       help='[DEPRECATED] Equilibrium length r0 - now auto-detected from force field')
    parser.add_argument('--harmonic-reset-temperature', type=float, default=None,
                       help='Target temperature for internal stretch DOF (default: None, use dynamic bath temperature at reset time)')
    parser.add_argument('--harmonic-reset-seed', type=int, default=42,
                       help='Random seed for harmonic reset (default: 42)')
    
    # Adaptive bath temperature controller parameters  
    parser.add_argument('--enable-adaptive-bath', action='store_true',
                       help='Enable adaptive bath temperature controller (default: False)')
    parser.add_argument('--adaptive-bath-amplitude-scale', type=float, default=1.0,
                       help='Amplitude scale for bath temperature adjustment (default: 1.0)')
    parser.add_argument('--adaptive-bath-time-constant', type=float, default=1.0,
                       help='Time constant tau for GD evolution in ps (default: 1.0)')
    parser.add_argument('--adaptive-bath-target-temperature', type=float, default=100.0,
                       help='Static target temperature in K (ignored if --adaptive-bath-dynamic-target) (default: 100.0)')
    parser.add_argument('--adaptive-bath-dynamic-target', action='store_true',
                       help='Set target temperature dynamically at turn-on time (default: False)')
    parser.add_argument('--adaptive-bath-turn-on-time', type=float, default=0.0,
                       help='Time to start adaptive bath control in ps (default: 0.0)')
    parser.add_argument('--adaptive-bath-turn-off-time', type=float, default=None,
                       help='Time to stop adaptive bath control in ps (default: None)')
    parser.add_argument('--adaptive-bath-update-interval', type=float, default=0.1,
                       help='Update interval for adaptive bath control in ps (default: 0.1)')
    parser.add_argument('--adaptive-bath-apply-to', choices=['molecular', 'cavity', 'both'], default='both',
                       help='Which thermostats to apply adaptive bath control to (default: both)')
    parser.add_argument('--adaptive-bath-T-min', type=float, default=0.1,
                       help='Minimum allowed bath temperature in K (default: 0.1)')
    parser.add_argument('--adaptive-bath-T-max', type=float, default=None,
                       help='Maximum allowed bath temperature in K (default: None)')
    parser.add_argument('--adaptive-bath-empirical-data-file', type=str, default=None,
                       help='Path to empirical data file for signal temperature calculation')
    parser.add_argument('--adaptive-bath-signal-temperature-method', 
                       choices=['harmonic_equipartition', 'kinetic', 'lj_coulombic', 'harmonic'], 
                       default='harmonic_equipartition',
                       help='Temperature method for control signal (default: harmonic_equipartition)')
    parser.add_argument('--adaptive-bath-dynamic-target-temperature-method', 
                       choices=['harmonic_equipartition', 'kinetic', 'lj_coulombic', 'harmonic'], 
                       default=None,
                       help='Temperature method for setting dynamic target (default: same as signal method)')
    
    # Quench controller parameters
    parser.add_argument('--enable-quench-controller', action='store_true',
                       help='Enable quench controller for instantaneous temperature changes (default: False)')
    parser.add_argument('--quench-initial-temperature', type=float, default=100.0,
                       help='Initial temperature before quench in K (default: 100.0)')
    parser.add_argument('--quench-target-temperature', type=float, default=50.0,
                       help='Target temperature after quench in K (default: 50.0)')
    parser.add_argument('--quench-time', type=float, default=50.0,
                       help='Time when quench occurs in ps (default: 50.0)')
    parser.add_argument('--quench-apply-to', choices=['molecular', 'cavity', 'both'], default='both',
                       help='Which thermostats to apply quench to (default: both)')
    
    # Offset temperature controller parameters
    parser.add_argument('--enable-offset-controller', action='store_true',
                       help='Enable offset temperature controller for fixed temperature offsets (default: False)')
    parser.add_argument('--offset-temperature-method', 
                       choices=['kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition', 'lj_coulombic_bath'], 
                       default='kinetic',
                       help='Temperature calculation method for offset controller (default: kinetic)')
    parser.add_argument('--offset-temperature-offset-K', type=float, default=-50.0,
                       help='Temperature offset in K (can be positive or negative) (default: -50.0)')
    parser.add_argument('--offset-turn-on-time', type=float, default=0.0,
                       help='Turn-on time for offset controller in ps (default: 0.0)')
    parser.add_argument('--offset-turn-off-time', type=float, default=None,
                       help='Turn-off time for offset controller in ps (default: None, never turn off)')
    parser.add_argument('--offset-update-interval', type=float, default=0.1,
                       help='Update interval for offset controller in ps (default: 0.1)')
    parser.add_argument('--offset-apply-to',
                       choices=['molecular', 'cavity', 'both'], default='both',
                       help='Apply offset to molecular, cavity, or both thermostats (default: both)')
    parser.add_argument('--offset-T-min', type=float, default=0.0,
                       help='Minimum allowed bath temperature in K (default: 0.0)')
    parser.add_argument('--offset-T-max', type=float, default=None,
                       help='Maximum allowed bath temperature in K (default: None)')
    
    # Differential equation controller parameters
    parser.add_argument('--enable-diffeq-controller', action='store_true',
                       help='Enable differential equation temperature controller (default: False)')
    parser.add_argument('--diffeq-temperature-method', 
                       choices=['kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition', 'lj_coulombic_bath'], 
                       default='kinetic',
                       help='Temperature calculation method for differential equation controller (default: kinetic)')
    parser.add_argument('--diffeq-time-constant', type=float, default=5.0,
                       help='Time constant τ in ps for differential equation dT/dt = -(T_bath - T_signal)/τ (default: 5.0)')
    parser.add_argument('--diffeq-time-constant-auto', action='store_true',
                       help='Use adaptive time constant τ = T[T_bath] from relaxation model instead of fixed value (default: False)')
    parser.add_argument('--diffeq-turn-on-time', type=float, default=0.0,
                       help='Turn-on time for differential equation controller in ps (default: 0.0)')
    parser.add_argument('--diffeq-turn-off-time', type=float, default=None,
                       help='Turn-off time for differential equation controller in ps (default: None, never turn off)')
    parser.add_argument('--diffeq-update-interval', type=float, default=0.1,
                       help='Update interval for differential equation controller in ps (default: 0.1)')
    parser.add_argument('--diffeq-apply-to',
                       choices=['molecular', 'cavity', 'both'], default='both',
                       help='Apply differential equation control to molecular, cavity, or both thermostats (default: both)')
    parser.add_argument('--diffeq-T-min', type=float, default=0.0,
                       help='Minimum allowed bath temperature in K (default: 0.0)')
    parser.add_argument('--diffeq-T-max', type=float, default=None,
                       help='Maximum allowed bath temperature in K (default: None)')
    parser.add_argument('--diffeq-disable-bias-estimation', action='store_true',
                       help='Disable Kalman filter bias estimation in DiffEq controller (default: False, bias estimation enabled)')
    parser.add_argument('--diffeq-rate-limit', type=float, default=None,
                       help='Maximum rate of temperature change in K/ps (default: None, no limit)')
    
    # Exact-cancellation + PI control parameters for DiffEq controller  
    parser.add_argument('--diffeq-enable-pi-control', action='store_true',
                       help='Enable exact-cancellation + PI control in DiffEq controller (default: False)')
    parser.add_argument('--diffeq-pi-rho', type=float, default=1.0,
                       help='Expected disturbance rate for PI control in K/ps (default: 1.0)')
    parser.add_argument('--diffeq-pi-epsilon', type=float, default=0.5, 
                       help='Target tolerance for PI control in K (default: 0.5)')
    parser.add_argument('--diffeq-pi-zeta', type=float, default=0.8,
                       help='Damping coefficient for PI control (default: 0.8)')
    parser.add_argument('--diffeq-relaxation-data-file', type=str, default=None,
                       help='Path to relaxation time data file for exact cancellation (default: None)')
    parser.add_argument('--diffeq-filter-window', type=float, default=0.0,
                       help='Moving window size for filtering T_signal measurements in ps (default: 0.0, no filtering)')
    
    # Adaptive bias cancellation parameters for DiffEq controller (mutually exclusive with PI control)
    parser.add_argument('--diffeq-enable-bias-cancellation', action='store_true',
                       help='Enable adaptive bias cancellation (mutually exclusive with PI control)')
    parser.add_argument('--diffeq-bias-tau-b', type=float, default=50.0,
                       help='Bath time constant tau_b in ps (should be > tau_s)')
    parser.add_argument('--diffeq-bias-tau-b-auto', action='store_true',
                       help='Compute tau_b automatically as 5 * tau_s from relaxation model')
    parser.add_argument('--diffeq-bias-kappa', type=float, default=0.01,
                       help='Bias estimation gain kappa (small for robustness)')
    parser.add_argument('--diffeq-bias-kappa-auto', action='store_true',
                       help='Compute kappa automatically as 1/(prefactor*tau_s) for proper time-scale hierarchy')
    parser.add_argument('--diffeq-bias-tau-b-prefactor', type=float, default=5.0,
                       help='Prefactor for auto tau_b calculation: tau_b = prefactor * tau_s (default: 5.0)')
    parser.add_argument('--diffeq-bias-kappa-prefactor', type=float, default=50.0,
                       help='Prefactor for auto kappa calculation: kappa = 1/(prefactor * tau_s) (default: 50.0)')
    parser.add_argument('--diffeq-bias-calibration-time', type=float, default=10.0,
                       help='Calibration time to estimate initial bias in ps')
    
    # Simple Setpoint Controller parameters
    parser.add_argument('--enable-simple-setpoint-controller', action='store_true',
                       help='Enable simple setpoint temperature controller (default: False)')
    parser.add_argument('--simple-setpoint-signal-method',
                       choices=['kinetic', 'lj_coulombic', 'harmonic'], default='kinetic',
                       help='Temperature method for setpoint capture (default: kinetic)')
    parser.add_argument('--simple-setpoint-time-constant-ps', type=float, default=5.0,
                       help='Time constant τ in ps for control law dT/dt = -(T_kinetic - T_setpoint)/τ (default: 5.0)')
    parser.add_argument('--simple-setpoint-apply-to',
                       choices=['molecular', 'cavity', 'both'], default='both',
                       help='Apply simple setpoint control to molecular, cavity, or both thermostats (default: both)')
    parser.add_argument('--simple-setpoint-turn-on-time-ps', type=float, default=0.0,
                       help='Turn-on time for simple setpoint controller in ps (default: 0.0)')
    parser.add_argument('--simple-setpoint-turn-off-time-ps', type=float, default=None,
                       help='Turn-off time for simple setpoint controller in ps (default: None, never turn off)')
    parser.add_argument('--simple-setpoint-update-interval-ps', type=float, default=0.1,
                       help='Update interval for simple setpoint controller in ps (default: 0.1)')
    parser.add_argument('--simple-setpoint-T-min', type=float, default=0.0,
                       help='Minimum allowed bath temperature in K (default: 0.0)')
    parser.add_argument('--simple-setpoint-T-max', type=float, default=None,
                       help='Maximum allowed bath temperature in K (default: None)')
    parser.add_argument('--simple-setpoint-output-file', type=str, default='simple_setpoint_control.csv',
                       help='Output file for simple setpoint control data (default: simple_setpoint_control.csv)')
    parser.add_argument('--simple-setpoint-console-output-period-ps', type=float, default=1.0,
                       help='Console output period for simple setpoint controller in ps (default: 1.0)')
    
    # Adaptive MPC controller parameters
    parser.add_argument('--enable-adaptive-mpc-controller', action='store_true',
                       help='Enable adaptive MPC controller (default: False)')
    parser.add_argument('--adaptive-mpc-target-temperature', type=float, default=100.0,
                       help='Target temperature T_star for MPC in K (default: 100.0)')
    parser.add_argument('--adaptive-mpc-turn-on-time-ps', type=float, default=0.0,
                       help='Time to start system ID phase in ps (default: 0.0)')
    parser.add_argument('--adaptive-mpc-turn-off-time-ps', type=float, default=None,
                       help='Time to stop MPC in ps (default: None, never stop)')
    parser.add_argument('--adaptive-mpc-system-id-duration-ps', type=float, default=50.0,
                       help='Duration of system ID phase in ps (default: 50.0)')
    parser.add_argument('--adaptive-mpc-system-id-step-duration-ps', type=float, default=5.0,
                       help='Duration of each step change in ps (default: 5.0)')
    parser.add_argument('--adaptive-mpc-system-id-seed', type=int, default=42,
                       help='Random seed for system ID step sequence (default: 42)')
    parser.add_argument('--adaptive-mpc-update-interval-ps', type=float, default=0.1,
                       help='MPC update period in ps (default: 0.1)')
    parser.add_argument('--adaptive-mpc-prediction-horizon', type=int, default=10,
                       help='MPC prediction horizon Np (default: 10)')
    parser.add_argument('--adaptive-mpc-control-horizon', type=int, default=5,
                       help='MPC control horizon Nc (default: 5)')
    parser.add_argument('--adaptive-mpc-output-weight', type=float, default=100.0,
                       help='MPC output tracking weight Q (default: 100.0)')
    parser.add_argument('--adaptive-mpc-control-effort-weight', type=float, nargs=2, default=[1.0, 0.1],
                       help='MPC control effort weights [R_lambda, R_Tbath] (default: [1.0, 0.1])')
    parser.add_argument('--adaptive-mpc-rate-penalty-weight', type=float, nargs=2, default=[10.0, 1.0],
                       help='MPC rate-of-change penalty [S_lambda, S_Tbath] (default: [10.0, 1.0])')
    parser.add_argument('--adaptive-mpc-lambda-min', type=float, default=0.0,
                       help='Minimum coupling strength (default: 0.0)')
    parser.add_argument('--adaptive-mpc-lambda-max', type=float, default=1e-2,
                       help='Maximum coupling strength (default: 1e-2)')
    parser.add_argument('--adaptive-mpc-T-bath-min', type=float, default=0.1,
                       help='Minimum bath temperature in K (default: 0.1)')
    parser.add_argument('--adaptive-mpc-T-bath-max', type=float, default=500.0,
                       help='Maximum bath temperature in K (default: 500.0)')
    parser.add_argument('--adaptive-mpc-delta-lambda-max', type=float, default=1e-4,
                       help='Maximum rate of lambda change per update (default: 1e-4)')
    parser.add_argument('--adaptive-mpc-delta-T-bath-max', type=float, default=10.0,
                       help='Maximum rate of T_bath change in K per update (default: 10.0)')
    parser.add_argument('--adaptive-mpc-apply-to', choices=['molecular', 'cavity', 'both'], default='both',
                       help='Apply MPC to molecular, cavity, or both baths (default: both)')
    parser.add_argument('--adaptive-mpc-rls-forgetting-factor', type=float, default=0.995,
                       help='RLS forgetting factor (default: 0.995)')
    parser.add_argument('--adaptive-mpc-rls-initial-covariance', type=float, default=100.0,
                       help='RLS initial covariance scaling (default: 100.0)')
    parser.add_argument('--adaptive-mpc-model-update-interval', type=int, default=10,
                       help='Update model every N steps during control phase (default: 10)')
    parser.add_argument('--adaptive-mpc-output-file', type=str, default='adaptive_mpc_control.csv',
                       help='Output CSV file for MPC data (default: adaptive_mpc_control.csv)')
    parser.add_argument('--adaptive-mpc-console-output-period-ps', type=float, default=1.0,
                       help='Console output period in ps (default: 1.0)')
    parser.add_argument('--adaptive-mpc-regularization-param', type=float, default=1e-3,
                       help='Tikhonov regularization parameter for model ID (default: 1e-3)')
    parser.add_argument('--adaptive-mpc-use-scaling', type=lambda x: x.lower() in ['true', '1', 'yes'], default=True,
                       help='Enable state/control scaling for better conditioning (default: True)')
    parser.add_argument('--adaptive-mpc-debug-mode', action='store_true',
                       help='Enable detailed MPC debug output (default: False)')
    
    # PID controller parameters
    parser.add_argument('--enable-pid-controller', action='store_true',
                       help='Enable PID temperature controller (default: False)')
    parser.add_argument('--pid-signal-choice', choices=['lj_coulombic', 'harmonic_fictive', 'kinetic'], 
                       default='lj_coulombic',
                       help='Temperature signal for PID control (default: lj_coulombic)')
    parser.add_argument('--pid-target-temperature', type=float, default=100.0,
                       help='Target temperature setpoint in K (ignored if --pid-self-loop) (default: 100.0)')
    parser.add_argument('--pid-self-loop', action='store_true',
                       help='Enable self-loop mode where setpoint = T_bath (default: False)')
    parser.add_argument('--pid-Kp', type=float, default=None,
                       help='Proportional gain (manual tuning, default: None = auto-tune)')
    parser.add_argument('--pid-Ti', type=float, default=None,
                       help='Integral time constant in ps (manual tuning, default: None = auto-tune)')
    parser.add_argument('--pid-Td', type=float, default=None,
                       help='Derivative time constant in ps (manual tuning, default: None = auto-tune)')
    parser.add_argument('--pid-auto-tune', type=lambda x: x.lower() in ['true', '1', 'yes'], default=True,
                       help='Enable auto-tuning via FOPDT step response (default: True)')
    parser.add_argument('--pid-auto-tune-step-size', type=float, default=20.0,
                       help='Step change size for auto-tuning in K (default: 20.0)')
    parser.add_argument('--pid-auto-tune-duration-ps', type=float, default=50.0,
                       help='Duration to observe step response in ps (default: 50.0)')
    parser.add_argument('--pid-turn-on-time-ps', type=float, default=0.0,
                       help='Time to activate PID controller in ps (default: 0.0)')
    parser.add_argument('--pid-turn-off-time-ps', type=float, default=None,
                       help='Time to deactivate PID controller in ps (default: None = never)')
    parser.add_argument('--pid-update-interval-ps', type=float, default=0.1,
                       help='PID update period in ps (default: 0.1)')
    parser.add_argument('--pid-apply-to', choices=['molecular', 'cavity', 'both'], default='both',
                       help='Apply PID to molecular, cavity, or both baths (default: both)')
    parser.add_argument('--pid-T-min', type=float, default=0.1,
                       help='Minimum bath temperature in K (default: 0.1)')
    parser.add_argument('--pid-T-max', type=float, default=None,
                       help='Maximum bath temperature in K (default: None)')
    parser.add_argument('--pid-rate-limit-K-per-ps', type=float, default=None,
                       help='Maximum rate of temperature change in K/ps (default: None)')
    parser.add_argument('--pid-derivative-filter-N', type=float, default=10.0,
                       help='Derivative filter parameter N (tau_f = Td/N) (default: 10.0)')
    parser.add_argument('--pid-enable-anti-windup', type=lambda x: x.lower() in ['true', '1', 'yes'], default=True,
                       help='Enable anti-windup for integral term (default: True)')
    parser.add_argument('--pid-output-file', type=str, default='pid_control.csv',
                       help='Output CSV file for PID data (default: pid_control.csv)')
    parser.add_argument('--pid-console-output-period-ps', type=float, default=1.0,
                       help='Console output period in ps (default: 1.0)')
    
    # BathPI flux-based controller parameters
    parser.add_argument('--enable-bath-pi-controller', action='store_true',
                       help='Enable BathPI flux-based temperature controller (default: False)')
    parser.add_argument('--bath-pi-apply-to',
                       choices=['molecular', 'cavity', 'both'], default='both',
                       help='Apply BathPI control to molecular, cavity, or both baths (default: both)')
    parser.add_argument('--bath-pi-K-p-molecular', type=float, default=0.1,
                       help='Proportional gain for molecular bath PI control (default: 0.1)')
    parser.add_argument('--bath-pi-K-i-molecular', type=float, default=0.01,
                       help='Integral gain for molecular bath PI control (default: 0.01)')
    parser.add_argument('--bath-pi-K-T-molecular', type=float, default=0.0,
                       help='Temperature error gain for molecular bath (chases LJ+Coul temperature) (default: 0.0)')
    parser.add_argument('--bath-pi-K-p-cavity', type=float, default=0.1,
                       help='Proportional gain for cavity bath PI control (default: 0.1)')
    parser.add_argument('--bath-pi-K-i-cavity', type=float, default=0.01,
                       help='Integral gain for cavity bath PI control (default: 0.01)')
    parser.add_argument('--bath-pi-K-T-cavity', type=float, default=0.0,
                       help='Temperature error gain for cavity bath (chases LJ+Coul temperature) (default: 0.0)')
    parser.add_argument('--bath-pi-filter-window-ps', type=str, default='5.0',
                       help='Exponential moving average filter window in ps or "auto" for auto-tuning (default: 5.0)')
    parser.add_argument('--bath-pi-flux-source', type=str, default='reservoir',
                       choices=['reservoir', 'lj_coul'],
                       help='Energy flux source: reservoir (thermostat) or lj_coul (system energy) (default: reservoir)')
    parser.add_argument('--bath-pi-anti-windup-alpha', type=float, default=0.01,
                       help='Anti-windup integral leakage coefficient (default: 0.01)')
    parser.add_argument('--bath-pi-enable-feedforward', action='store_true',
                       help='Enable optional slow feedforward recentering to T_nominal (default: False)')
    parser.add_argument('--bath-pi-T-nominal', type=float, default=None,
                       help='Nominal temperature for feedforward in K (default: None)')
    parser.add_argument('--bath-pi-feedforward-tau-ps', type=float, default=1000.0,
                       help='Feedforward time constant in ps (default: 1000.0)')
    parser.add_argument('--bath-pi-turn-on-time-ps', type=float, default=0.0,
                       help='Turn-on time for BathPI controller in ps (default: 0.0)')
    parser.add_argument('--bath-pi-turn-off-time-ps', type=float, default=None,
                       help='Turn-off time for BathPI controller in ps (default: None)')
    parser.add_argument('--bath-pi-update-interval-ps', type=float, default=0.1,
                       help='Update interval for BathPI controller in ps (default: 0.1)')
    parser.add_argument('--bath-pi-T-min', type=float, default=0.1,
                       help='Minimum allowed bath temperature in K (default: 0.1)')
    parser.add_argument('--bath-pi-T-max', type=float, default=None,
                       help='Maximum allowed bath temperature in K (default: None)')
    parser.add_argument('--bath-pi-rate-limit-K-per-ps', type=float, default=None,
                       help='Maximum rate of temperature change in K/ps (default: None)')
    parser.add_argument('--bath-pi-output-file', type=str, default='bath_pi_control.csv',
                       help='Output CSV file for BathPI controller data (default: bath_pi_control.csv)')
    parser.add_argument('--bath-pi-relaxation-data-file', type=str, default=None,
                       help='Path to relaxation time data file for auto-tuning filter window (default: None)')
    
    # LQR temperature controller parameters
    parser.add_argument('--enable-lqr-controller', action='store_true',
                       help='Enable LQR optimal temperature controller with Kalman filter (default: False)')
    parser.add_argument('--lqr-signal-method', 
                       choices=['kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition'], 
                       default='lj_coulombic',
                       help='Temperature method for regulated signal (default: lj_coulombic)')
    parser.add_argument('--lqr-hot-method', 
                       choices=['kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition'], 
                       default='harmonic_equipartition',
                       help='Temperature method for disturbance signal (default: harmonic_equipartition)')
    parser.add_argument('--lqr-target-temperature', type=float, default=300.0,
                       help='Fixed target temperature in K (used if not dynamic) (default: 300.0)')
    parser.add_argument('--lqr-dynamic-target', action='store_true',
                       help='Set target temperature dynamically from signal at turn-on (default: False)')
    parser.add_argument('--lqr-dynamic-target-method', 
                       choices=['kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition'], 
                       default=None,
                       help='Temperature method for dynamic target (default: same as signal method)')
    parser.add_argument('--lqr-weight-signal', type=float, default=100.0,
                       help='LQR weight for signal regulation (higher = tighter) (default: 100.0)')
    parser.add_argument('--lqr-weight-hot', type=float, default=1.0,
                       help='LQR weight for hot temperature (lower = allow more variation) (default: 1.0)')
    parser.add_argument('--lqr-weight-bath', type=float, default=0.1,
                       help='[DEPRECATED - LQG uses 2D state, bath is control input] Old 3D formulation parameter (default: 0.1)')
    parser.add_argument('--lqr-weight-integral', type=float, default=10.0,
                       help='LQR weight for integral action (default: 10.0)')
    parser.add_argument('--lqr-control-effort', type=float, default=1.0,
                       help='LQR penalty on control effort (higher = gentler) (default: 1.0)')
    parser.add_argument('--lqr-process-noise-signal', type=float, default=0.1,
                       help='Process noise std dev for signal (K) (default: 0.1)')
    parser.add_argument('--lqr-process-noise-hot', type=float, default=0.1,
                       help='Process noise std dev for hot (K) (default: 0.1)')
    parser.add_argument('--lqr-measurement-noise-signal', type=float, default=0.5,
                       help='Measurement noise std dev for signal (K) (default: 0.5)')
    parser.add_argument('--lqr-measurement-noise-hot', type=float, default=0.5,
                       help='Measurement noise std dev for hot (K) (default: 0.5)')
    parser.add_argument('--lqr-system-id-mode', 
                       choices=['step', 'multi_step', 'load', 'skip', 'none'], 
                       default='step',
                       help='System ID mode: step (single quench), multi_step (4 sequential excitations), load (from file), skip (manual), none (adaptive no-quench) (default: step)')
    parser.add_argument('--lqr-system-id-temp', type=float, default=5.0,
                       help='Quench temperature for system ID (K) (default: 5.0)')
    parser.add_argument('--lqr-system-id-duration', type=float, default=50.0,
                       help='Duration of system ID cold quench (ps) (default: 50.0)')
    parser.add_argument('--lqr-system-id-file', type=str, default='lqr_system_params.json',
                       help='File to save/load system parameters (default: lqr_system_params.json)')
    parser.add_argument('--lqr-periodic-system-id', action='store_true',
                       help='Enable periodic re-identification of system parameters')
    parser.add_argument('--lqr-periodic-system-id-interval', type=float, default=1000.0,
                       help='Time interval between periodic system IDs (ps) (default: 1000.0)')
    
    # EKF-based Adaptive LQG parameters
    parser.add_argument('--lqr-use-ekf-adaptation', action='store_true', default=True,
                       help='Use EKF for online parameter estimation (default: True)')
    parser.add_argument('--lqr-no-ekf-adaptation', dest='lqr_use_ekf_adaptation', action='store_false',
                       help='Disable EKF adaptation and use RLS instead')
    parser.add_argument('--lqr-ekf-update-interval', type=int, default=50,
                       help='EKF update interval (control steps) (default: 50)')
    parser.add_argument('--lqr-ekf-process-noise-param', type=float, default=0.001,
                       help='EKF parameter drift noise std (default: 0.001)')
    parser.add_argument('--lqr-ekf-initial-covariance-param', type=float, default=0.1,
                       help='EKF initial parameter uncertainty (default: 0.1)')
    parser.add_argument('--lqr-adaptive-lqr-threshold', type=float, default=0.05,
                       help='Relative parameter change threshold for LQR redesign (default: 0.05 = 5%%)')
    
    # Gain Scheduling parameters
    parser.add_argument('--lqr-enable-gain-scheduling', action='store_true', default=True,
                       help='Enable gain scheduling based on temperature spread (default: True)')
    parser.add_argument('--lqr-disable-gain-scheduling', dest='lqr_enable_gain_scheduling', action='store_false',
                       help='Disable gain scheduling')
    parser.add_argument('--lqr-gain-schedule-far-threshold', type=float, default=20.0,
                       help='Temperature spread for full gain (K) (default: 20.0)')
    parser.add_argument('--lqr-gain-schedule-near-threshold', type=float, default=10.0,
                       help='Temperature spread for reduced gain (K) (default: 10.0)')
    
    # T_h low-pass filter parameters
    parser.add_argument('--lqr-th-filter-enabled', action='store_true', default=True,
                       help='Enable low-pass filtering on T_h measurements (default: True)')
    parser.add_argument('--lqr-th-filter-disabled', dest='lqr_th_filter_enabled',
                       action='store_false',
                       help='Disable T_h low-pass filtering')
    parser.add_argument('--lqr-th-filter-time-constant', type=float, default=20.0,
                       help='T_h filter time constant in ps (default: 20.0, cutoff ~8Hz)')
    
    # Gentle startup parameters
    parser.add_argument('--lqr-gentle-startup-steps', type=int, default=10,
                       help='Number of steps to ramp up control authority (default: 10)')
    parser.add_argument('--lqr-gentle-startup-min-authority', type=float, default=0.1,
                       help='Starting control authority fraction (default: 0.1 = 10%%)')
    
    # Kinetic temperature tracking parameters (3D state augmentation)
    parser.add_argument('--lqr-track-kinetic-temp', action='store_true', default=False,
                       help='Enable kinetic temperature as 3rd state variable (default: False, uses 2D system)')
    parser.add_argument('--lqr-weight-kinetic', type=float, default=100.0,
                       help='LQR weight for kinetic temperature error (default: 100.0, same as signal temp)')
    parser.add_argument('--lqr-process-noise-kinetic', type=float, default=2.0,
                       help='Kalman filter process noise for kinetic temp (K) (default: 2.0)')
    parser.add_argument('--lqr-measurement-noise-kinetic', type=float, default=2.0,
                       help='Kalman filter measurement noise for kinetic temp (K) (default: 2.0)')
    
    # Cross-coupling weights for thermal equilibration (3D mode only)
    parser.add_argument('--lqr-cross-coupling-signal-kinetic', type=float, default=0.0,
                       help='Cross-coupling weight for (T_signal - T_kinetic)² penalty (default: 0.0, disabled)')
    parser.add_argument('--lqr-cross-coupling-signal-hot', type=float, default=0.0,
                       help='Cross-coupling weight for (T_signal - T_hot)² penalty (default: 0.0, disabled)')
    parser.add_argument('--lqr-cross-coupling-hot-kinetic', type=float, default=0.0,
                       help='Cross-coupling weight for (T_hot - T_kinetic)² penalty (default: 0.0, disabled)')

    # Adaptive LQI Controller parameters
    parser.add_argument('--lqr-controller-type', type=str, choices=['standard', 'adaptive_lqi'], default='standard',
                       help='Controller type: standard (basic LQR) or adaptive_lqi (full mode-based) (default: standard)')
    parser.add_argument('--lqr-tau-L-initial', type=float, default=200.0,
                       help='Initial guess for signal time constant (ps) for adaptive_lqi (default: 200.0)')
    parser.add_argument('--lqr-tau-H-initial', type=float, default=30.0,
                       help='Initial guess for hot time constant (ps) for adaptive_lqi (default: 30.0)')
    parser.add_argument('--lqr-k-initial', type=float, default=0.01,
                       help='Initial guess for coupling coefficient for adaptive_lqi (default: 0.01)')
    parser.add_argument('--lqr-tau-b', type=float, default=1.0,
                       help='Bath actuator time constant (ps) for adaptive_lqi (default: 1.0)')
    parser.add_argument('--lqr-q-common', type=float, default=1000.0,
                       help='Common mode state weight for adaptive_lqi (default: 1000.0)')
    parser.add_argument('--lqr-q-diff', type=float, default=100.0,
                       help='Difference mode state weight for adaptive_lqi (default: 100.0)')
    parser.add_argument('--lqr-q-eta-common', type=float, default=1e6,
                       help='Common mode integral weight for adaptive_lqi (default: 1e6)')
    parser.add_argument('--lqr-q-eta-diff', type=float, default=10.0,
                       help='Difference mode integral weight for adaptive_lqi (default: 10.0)')
    parser.add_argument('--lqr-process-noise-drift', type=float, default=0.01,
                       help='Drift process noise for adaptive_lqi (default: 0.01)')
    parser.add_argument('--lqr-rls-forgetting', type=float, default=0.998,
                       help='RLS forgetting factor for adaptive_lqi (default: 0.998)')
    parser.add_argument('--lqr-rls-regularization', type=float, default=1e-6,
                       help='RLS regularization for adaptive_lqi (default: 1e-6)')
    parser.add_argument('--lqr-rls-update-interval', type=int, default=50,
                       help='RLS update interval (samples) for adaptive_lqi (default: 50)')
    parser.add_argument('--lqr-max-control-rate', type=float, default=10.0,
                       help='Maximum control rate (K/ps) for adaptive_lqi (default: 10.0)')
    parser.add_argument('--lqr-integral-max-common', type=float, default=1000.0,
                       help='Maximum common mode integral for adaptive_lqi (default: 1000.0)')
    parser.add_argument('--lqr-integral-max-diff', type=float, default=100.0,
                       help='Maximum difference mode integral for adaptive_lqi (default: 100.0)')
    parser.add_argument('--lqr-theta-change-threshold', type=float, default=0.05,
                       help='Parameter change threshold for relinearization for adaptive_lqi (default: 0.05)')
    
    parser.add_argument('--lqr-turn-on-time', type=float, default=0.0,
                       help='Turn-on time for LQR controller in ps (default: 0.0)')
    parser.add_argument('--lqr-turn-off-time', type=float, default=None,
                       help='Turn-off time for LQR controller in ps (default: None, never turn off)')
    parser.add_argument('--lqr-update-interval', type=float, default=0.1,
                       help='Control update interval in ps (default: 0.1)')
    parser.add_argument('--lqr-T-min', type=float, default=0.1,
                       help='Minimum allowed bath temperature in K (default: 0.1)')
    parser.add_argument('--lqr-T-max', type=float, default=None,
                       help='Maximum allowed bath temperature in K (default: None)')
    parser.add_argument('--lqr-apply-to',
                       choices=['molecular', 'cavity', 'both'], default='both',
                       help='Apply LQR control to molecular, cavity, or both thermostats (default: both)')
    parser.add_argument('--lqr-output-file', type=str, default='lqr_controller.csv',
                       help='Output file for LQR control data (default: lqr_controller.csv)')
    parser.add_argument('--lqr-empirical-data-file', type=str, default=None,
                       help='Path to empirical data file for LJ+Coulombic temperature calculation')
    
    # LQG pulse controller parameters
    parser.add_argument('--enable-lqg-controller', action='store_true',
                       help='Enable LQG pulse-to-pulse optimal controller (default: False)')
    parser.add_argument('--lqg-A-matrix', type=str, default=None,
                       help='JSON string or file path for A matrix (state transition)')
    parser.add_argument('--lqg-B-matrix', type=str, default=None,
                       help='JSON string or file path for B matrix (control input)')
    parser.add_argument('--lqg-C-matrix', type=str, default=None,
                       help='JSON string or file path for C matrix (measurement)')
    parser.add_argument('--lqg-S-matrix', type=str, default=None,
                       help='JSON string or file path for S matrix (tracked outputs)')
    parser.add_argument('--lqg-r-target', type=str, default=None,
                       help='JSON string or file path for target reference vector')
    parser.add_argument('--lqg-Qz-weights', type=str, default=None,
                       help='JSON string or file path for output tracking penalty matrix')
    parser.add_argument('--lqg-Ru-weights', type=str, default=None,
                       help='JSON string or file path for control effort penalty matrix')
    parser.add_argument('--lqg-Qe-weights', type=str, default=None,
                       help='JSON string or file path for integral error penalty matrix')
    parser.add_argument('--lqg-W-noise', type=str, default=None,
                       help='JSON string or file path for process noise covariance matrix')
    parser.add_argument('--lqg-V-noise', type=str, default=None,
                       help='JSON string or file path for measurement noise covariance matrix')
    parser.add_argument('--lqg-g0-baseline', type=float, default=1e-4,
                       help='Baseline coupling strength g0 in a.u. (default: 1e-4)')
    parser.add_argument('--lqg-dt', type=float, default=5.0,
                       help='LQG pulse period in ps (default: 5.0)')
    parser.add_argument('--lqg-delta-g-min', type=float, default=-0.3,
                       help='Minimum coupling deviation Δg in a.u. (default: -0.3)')
    parser.add_argument('--lqg-delta-g-max', type=float, default=0.3,
                       help='Maximum coupling deviation Δg in a.u. (default: 0.3)')
    parser.add_argument('--lqg-T-bath-min', type=float, default=10.0,
                       help='Minimum bath temperature in K (default: 10.0)')
    parser.add_argument('--lqg-T-bath-max', type=float, default=500.0,
                       help='Maximum bath temperature in K (default: 500.0)')
    parser.add_argument('--lqg-rate-limit', type=float, default=1.0,
                       help='Maximum rate of bath temperature change in K/ps (default: 1.0, None to disable)')
    parser.add_argument('--lqg-turn-on-time', type=float, default=0.0,
                       help='Turn-on time for LQG controller in ps (default: 0.0)')
    parser.add_argument('--lqg-turn-off-time', type=float, default=None,
                       help='Turn-off time for LQG controller in ps (default: None)')
    parser.add_argument('--lqg-update-interval', type=float, default=5.0,
                       help='LQG update interval in ps (default: 5.0)')
    parser.add_argument('--lqg-console-output-period', type=float, default=10.0,
                       help='Console output period for LQG in ps (default: 10.0)')
    parser.add_argument('--lqg-c-affine-term', type=str, default=None,
                       help='JSON string or file path for affine term c vector')
    parser.add_argument('--lqg-temperature-methods', type=str, default=None,
                       help='JSON string mapping state indices to temperature methods')
    
    # LQG square wave parameters
    parser.add_argument('--lqg-square-period', type=float, default=5.0,
                       help='LQG square wave period in ps (default: 5.0)')
    parser.add_argument('--lqg-square-duty-cycle', type=float, default=0.5,
                       help='LQG square wave duty cycle 0.0-1.0 (default: 0.5)')
    parser.add_argument('--lqg-square-phase-offset', type=float, default=0.0,
                       help='LQG square wave phase offset in radians (default: 0.0)')
    parser.add_argument('--lqg-square-start-time', type=float, default=0.0,
                       help='LQG square wave start time in ps (default: 0.0)')
    parser.add_argument('--lqg-square-stop-time', type=float, default=None,
                       help='LQG square wave stop time in ps (default: None)')
    parser.add_argument('--lqg-square-g-low-level', type=float, default=1e-6,
                       help='LQG square wave low coupling level in a.u. (default: 1e-6)')
    
    # LQG system identification parameters
    parser.add_argument('--lqg-enable-system-id', action='store_true',
                       help='Enable automatic system identification (default behavior, can be overridden by --lqg-disable-system-id)')
    parser.add_argument('--lqg-disable-system-id', action='store_true',
                       help='Disable system identification (use provided matrices or seed model)')
    parser.add_argument('--lqg-equilibrium-duration', type=float, default=500.0,
                       help='Duration for equilibrium data collection in ps (default: 500.0)')
    parser.add_argument('--lqg-prbs-n-unique-steps', type=int, default=10,
                       help='Number of unique coupling levels for PRBS excitation (default: 10)')
    parser.add_argument('--lqg-prbs-amplitude', type=float, default=5e-5,
                       help='PRBS excitation amplitude in a.u. (default: 5e-5)')
    parser.add_argument('--lqg-bath-n-unique-steps', type=int, default=8,
                       help='Number of unique bath temperature levels for excitation (default: 8)')
    parser.add_argument('--lqg-bath-excitation-amplitude', type=float, default=10.0,
                       help='Bath temperature excitation amplitude in K (default: 10.0)')
    parser.add_argument('--lqg-system-id-output-file', type=str, default='lqg_system_id.json',
                       help='Output file for identified model (default: lqg_system_id.json)')
    
    # LQG coupling controller parameters (new coupling-type: lqg_coupling)
    parser.add_argument('--lqg-coupling-target-temperature', type=float, default=100.0,
                       help='Target temperature for T_s regulation in K (default: 100.0)')
    parser.add_argument('--lqg-coupling-update-interval', type=float, default=0.1,
                       help='LQG coupling controller update period in ps (default: 0.1)')
    parser.add_argument('--lqg-coupling-equilibrium-duration', type=float, default=5.0,
                       help='Initial equilibrium phase duration for system ID in ps (default: 5.0)')
    parser.add_argument('--lqg-coupling-step-duration', type=float, default=5.0,
                       help='Duration for each step in system ID excitation in ps (default: 5.0)')
    parser.add_argument('--lqg-coupling-n-steps', type=int, default=2,
                       help='Number of step up/down cycles for system ID (default: 2)')
    parser.add_argument('--lqg-coupling-lambda-min', type=float, default=0.0,
                       help='Minimum allowed λ value (default: 0.0)')
    parser.add_argument('--lqg-coupling-lambda-max', type=float, default=1e-2,
                       help='Maximum allowed λ value (default: 1e-2)')
    parser.add_argument('--lqg-coupling-process-noise', type=float, default=0.1,
                       help='Process noise standard deviation for Kalman filter in K (default: 0.1)')
    parser.add_argument('--lqg-coupling-measurement-noise', type=float, default=0.5,
                       help='Measurement noise standard deviation for Kalman filter in K (default: 0.5)')
    parser.add_argument('--lqg-coupling-weight-signal', type=float, default=100.0,
                       help='LQR weight for T_s regulation (higher = tighter control) (default: 100.0)')
    parser.add_argument('--lqg-coupling-weight-harmonic', type=float, default=1.0,
                       help='LQR weight for T_harmonic observation (default: 1.0)')
    parser.add_argument('--lqg-coupling-weight-kinetic', type=float, default=1.0,
                       help='LQR weight for T_kinetic observation (default: 1.0)')
    parser.add_argument('--lqg-coupling-weight-bath', type=float, default=0.1,
                       help='LQR weight for T_bath observation (default: 0.1)')
    parser.add_argument('--lqg-coupling-control-effort', type=float, default=1.0,
                       help='LQR penalty on control effort (higher = gentler) (default: 1.0)')
    parser.add_argument('--lqg-coupling-integral-gain', type=float, default=2.0,
                       help='Integral action gain for zero steady-state error (default: 2.0)')
    parser.add_argument('--lqg-coupling-system-id-file', type=str, default='lqg_coupling_system_id.json',
                       help='JSON file to save identified system parameters (default: lqg_coupling_system_id.json)')
    parser.add_argument('--lqg-coupling-control-file', type=str, default='lqg_coupling_control.csv',
                       help='CSV file to save control trajectory (default: lqg_coupling_control.csv)')
    parser.add_argument('--lqg-coupling-temperature-methods', type=str, 
                       default='["lj_coulombic", "harmonic_equipartition", "kinetic"]',
                       help='JSON list of temperature methods for [T_s, T_harmonic, T_kinetic] (default: ["lj_coulombic", "harmonic_equipartition", "kinetic"])')
    parser.add_argument('--lqg-coupling-bath-temperature-method', type=str, default='kinetic',
                       help='Temperature method for T_bath observation (default: kinetic)')
    
    # Legacy periodic coupling parameters (for backward compatibility)
    parser.add_argument('--periodic', action='store_true',
                       help='Enable periodic coupling oscillation: g(t) = A * sin(2π*t/T + φ)')
    parser.add_argument('--period', type=float, default=1.0,
                       help='Oscillation period in ps (default: 1.0 ps = 1 THz)')
    parser.add_argument('--phase-offset', type=float, default=0.0,
                       help='Phase offset in radians (default: 0.0)')
    parser.add_argument('--start-time', type=float, default=0.0,
                       help='Time to start periodic oscillation in ps (default: 0.0)')
    
    # Laser drive parameters
    parser.add_argument('--laser', action='store_true',
                       help='Enable laser drive forcing: F_laser = F_L * cos(ω_L * t)')
    parser.add_argument('--laser-frequency', type=float, default=1560.0,
                       help='Laser frequency in cm⁻¹ (default: 1560.0)')
    parser.add_argument('--laser-amplitude', type=float, default=1e-5,
                       help='Laser amplitude in a.u. (default: 1e-5)')
    
    
    # Mechanical cavity modulation parameters
    parser.add_argument('--mech-periodic', action='store_true',
                       help='Enable mechanical cavity length modulation: L(t) = L₀ + A*sin(ω*t)')
    parser.add_argument('--mech-frequency', type=float, default=100.0,
                       help='Mechanical modulation frequency in cm⁻¹ (default: 100.0)')
    parser.add_argument('--mech-magnitude', type=float, default=1e-4,
                       help='Mechanical modulation magnitude in a.u. (default: 1e-4)')
    
    # Switch-time parameters
    parser.add_argument('--switch-time', type=float, default=None,
                       help='Time in ps when coupling and dissipation turn on (default: on from start)')
    parser.add_argument('--decay-time-constant', type=float, default=None,
                       help='Exponential decay time constant in ps for time-varying coupling (default: None)')
    parser.add_argument('--damping-ratio', type=float, default=0.0,
                       help='Damping ratio zeta for dissipation calculation (default: 0.0)')
    
    # Physical parameters
    parser.add_argument('--temperature', type=float, default=100.0,
                       help='System temperature in Kelvin (default: 100.0)')
    parser.add_argument('--frequency', type=float, default=1560.0,
                       help='Cavity frequency in cm⁻¹ (default: 1560.0)')
    parser.add_argument('--runtime', type=float, default=500.0,
                       help='Simulation runtime in ps (default: 500.0)')
    parser.add_argument('--no-cavity', action='store_true',
                       help='Run without cavity coupling (molecular-only simulation)')
    
    # Input/Output parameters
    parser.add_argument('--input-gsd', type=str, default='molecular-0.gsd',
                       help='Input GSD file path (default: molecular-0.gsd)')
    parser.add_argument('--frame', type=int, default=-1,
                       help='Frame number to read from GSD file (default: -1, last frame)')
    parser.add_argument('--replicas', type=str, default='0',
                       help='Replica specification: single number (e.g., "5") or range (e.g., "0-4") (default: "0")')
    
    # Thermostat parameters
    parser.add_argument('--molecular-tau', type=float, default=5.0,
                       help='Molecular thermostat time constant in ps (default: 5.0)')
    parser.add_argument('--cavity-tau', type=float, default=5.0,
                       help='Cavity thermostat time constant in ps (default: 5.0)')
    
    # Timestep parameters
    parser.add_argument('--fixed-timestep', action='store_true',
                       help='Use fixed timestep instead of adaptive timestep')
    parser.add_argument('--timestep', type=float, default=0.1,
                       help='Fixed timestep in fs (only used with --fixed-timestep) (default: 0.1)')
    parser.add_argument('--error-tolerance', type=float, default=0.01,
                       help='Error tolerance for adaptive timestep (ignored with --fixed-timestep) (default: 0.01)')
    parser.add_argument('--initial-fraction', type=float, default=1e-5,
                       help='Initial fraction for shock dampening (default: 1e-5)')
    parser.add_argument('--time-constant-ps', type=float, default=50.0,
                       help='Time constant for error tolerance ramping in ps (default: 50.0)')
    
    # Dynamic coupling detection parameters
    parser.add_argument('--enable-dynamic-coupling-detection', action='store_true', default=True,
                       help='Enable dynamic coupling change detection for shock dampening (default: True)')
    parser.add_argument('--disable-dynamic-coupling-detection', action='store_true',
                       help='Disable dynamic coupling change detection, use only legacy switch_time')
    parser.add_argument('--coupling-change-threshold', type=float, default=1e-5,
                       help='Threshold for detecting large coupling changes in a.u. (default: 1e-5)')
    
    # Energy tracking parameters
    parser.add_argument('--enable-energy-tracker', action='store_true', default=True,
                       help='Enable comprehensive energy tracking (default: True)')
    parser.add_argument('--disable-energy-tracker', dest='enable_energy_tracker', action='store_false',
                       help='Disable comprehensive energy tracking')
    parser.add_argument('--energy-output-period-ps', type=float, default=0.1,
                       help='Energy output period in ps (default: 0.1)')
    parser.add_argument('--fkt-output-period-ps', type=float, default=1.0,
                       help='F(k,t) output period in ps (default: 1.0)')
    parser.add_argument('--gsd-output-period-ps', type=float, default=50.0,
                       help='GSD trajectory output period in ps (default: 50.0)')
    parser.add_argument('--console-output-period-ps', type=float, default=1.0,
                       help='Console output period in ps (default: 1.0)')
    
    # F(k,t) tracking parameters
    parser.add_argument('--enable-fkt', action='store_true',
                       help='Enable F(k,t) density correlation tracking')
    parser.add_argument('--fkt-kmag', type=float, default=1.0,
                       help='F(k,t) k-magnitude (default: 1.0)')
    parser.add_argument('--fkt-wavevectors', type=int, default=50,
                       help='Number of F(k,t) wavevectors (default: 50)')
    parser.add_argument('--fkt-ref-interval', type=float, default=1.0,
                       help='F(k,t) reference interval in ps (default: 1.0)')
    parser.add_argument('--fkt-max-refs', type=int, default=10,
                       help='Maximum F(k,t) reference states (default: 10)')
    
    # Dipole autocorrelation parameters
    parser.add_argument('--enable-dipole-autocorr', action='store_true',
                       help='Enable dipole autocorrelation tracking')
    parser.add_argument('--dipole-ref-interval', type=float, default=1.0,
                       help='Dipole reference interval in ps (default: 1.0)')
    parser.add_argument('--dipole-max-refs', type=int, default=10,
                       help='Maximum dipole reference states (default: 10)')
    parser.add_argument('--dipole-output-period-ps', type=float, default=1.0,
                       help='Dipole output period in ps (default: 1.0)')
    
    # Limit energy output time
    parser.add_argument('--max-energy-output-time', type=float, default=None,
                       help='Maximum time for energy output in ps (default: None, no limit)')
    
    # Hardware parameters
    parser.add_argument('--device', choices=['CPU', 'GPU'], default='CPU',
                       help='Compute device (default: CPU)')
    parser.add_argument('--gpu-id', type=int, default=0,
                       help='GPU device ID (default: 0)')
    parser.add_argument('--truncate-gsd', action='store_true',
                       help='Truncate existing GSD file (default: append)')
    parser.add_argument('--seed', type=int, default=None,
                       help='Random seed for simulation (default: replica-based deterministic seed)')
    parser.add_argument('--no-restart-velocities', action='store_true',
                       help='Do not restart velocities - use existing velocities from GSD file (default: restart velocities)')
    
    # Momentum zeroing control
    parser.add_argument('--zero-momentum', action='store_true', 
                       help='Enable periodic momentum zeroing to prevent center-of-mass drift (default: disabled)')
    parser.add_argument('--zero-momentum-period-ps', type=float, default=1.0, 
                       help='Period for momentum zeroing in ps (default: 1.0)')
    
    # Empirical feedback parameters
    parser.add_argument('--enable-empirical-feedback', action='store_true',
                       help='Enable empirical temperature feedback control based on systemic temperature')
    parser.add_argument('--empirical-data-file', type=str, default=None,
                       help='Path to empirical energy-temperature data file (required if --enable-empirical-feedback)')
    parser.add_argument('--feedback-output-period-ps', type=float, default=0.1,
                       help='Empirical feedback CSV output period in ps (default: 0.1)')
    parser.add_argument('--feedback-energy-component', choices=['total_PE', 'lj_coulombic', 'harmonic'], default='lj_coulombic',
                       help='Energy component for empirical feedback (default: lj_coulombic)')
    parser.add_argument('--feedback-apply-to', choices=['molecular', 'cavity', 'both', 'none'], default='both',
                       help='Which thermostats to apply empirical feedback to (default: both)')
    parser.add_argument('--feedback-averaging-window-ps', type=float, default=5.0,
                       help='Time window for averaging fictive temperature in ps (default: 5.0)')
    parser.add_argument('--feedback-update-interval-ps', type=float, default=5.0,
                       help='Interval between temperature updates in ps (default: 5.0)')
    parser.add_argument('--feedback-T-min', type=float, default=50.0,
                       help='Minimum allowed temperature in K (default: 50.0)')
    parser.add_argument('--feedback-T-max', type=float, default=300.0,
                       help='Maximum allowed temperature in K (default: 300.0)')
    parser.add_argument('--feedback-turn-off-time-ps', type=float, default=None,
                       help='Time in ps to turn off empirical feedback (default: None, feedback runs until end)')
    parser.add_argument('--feedback-use-direct-harmonic', action='store_true',
                       help='Use direct harmonic calculation for fictive temperature (default: empirical fitting)')
    parser.add_argument('--feedback-auto-adjust-molecular-tau', action='store_true',
                       help='Auto-adjust molecular thermostat tau when feedback turns off (default: disabled)')
    parser.add_argument('--feedback-auto-stopping', action='store_true',
                       help='Enable automatic stopping when fictive and harmonic temperatures converge (default: disabled)')
    parser.add_argument('--auto-stop-min-time-ps', type=float, default=10.0,
                       help='Minimum time before auto-stop can trigger in ps (default: 10.0)')
    parser.add_argument('--auto-stop-smoothing-window', type=int, default=5,
                       help='Window size for temperature smoothing in auto-stop detection (default: 5)')
    
    # Comprehensive temperature tracker parameters
    parser.add_argument('--enable-temp-tracker', action='store_true', default=True,
                       help='Enable comprehensive temperature tracker (default: True)')
    parser.add_argument('--disable-temp-tracker', dest='enable_temp_tracker', action='store_false',
                       help='Disable comprehensive temperature tracker')
    parser.add_argument('--temp-tracker-output-period-ps', type=float, default=0.1,
                       help='Temperature tracker output period in ps (default: 0.1)')
    parser.add_argument('--temp-tracker-empirical-data-file', type=str, default=None,
                       help='Empirical data file for LJ+Coul fictive temperature (default: None)')

    # HDF5 observable output parameters
    parser.add_argument('--enable-hdf5-output', action='store_true', default=True,
                       help='Enable HDF5 observable output (default: True)')
    parser.add_argument('--disable-hdf5-output', dest='enable_hdf5_output', action='store_false',
                       help='Disable HDF5 observable output')
    parser.add_argument('--hdf5-output-file', type=str, default=None,
                       help='HDF5 output file path (default: observables_replica_X.h5)')
    parser.add_argument('--hdf5-output-period-ps', type=float, default=0.01,
                       help='HDF5 output period in ps (default: 0.01)')

    # Molecular temperature decomposition parameters
    parser.add_argument('--enable-molecular-temps', action='store_true',
                       help='Enable molecular temperature decomposition (translational/rotational/vibrational)')
    parser.add_argument('--molecular-temps-output-period-ps', type=float, default=1.0,
                       help='Molecular temperature output period in ps (default: 1.0)')
    
    # Dipole moment FDR parameters
    parser.add_argument('--enable-dipole-fdr', action='store_true',
                       help='Enable dipole moment FDR analysis (autocorrelation tracking)')
    parser.add_argument('--dipole-fdr-output-period-ps', type=float, default=0.1,
                       help='Dipole FDR output period in ps (default: 0.1)')
    parser.add_argument('--dipole-fdr-max-correlation-time-ps', type=float, default=100.0,
                       help='Maximum correlation time for dipole FDR in ps (default: 100.0)')
    parser.add_argument('--dipole-fdr-field-direction', type=float, nargs=3, default=[0, 0, 1],
                       help='Electric field direction for FDR [x, y, z] (default: [0, 0, 1])')
    parser.add_argument('--dipole-fdr-exclude-cavity', action='store_true', default=True,
                       help='Exclude cavity particles from dipole calculation (default: True)')
    parser.add_argument('--enable-dipole-response', action='store_true',
                       help='Enable dipole response force for fork-and-clone FDR measurements')
    parser.add_argument('--dipole-response-field-strength', type=float, default=1e-5,
                       help='Electric field strength for response measurements (default: 1e-5)')
    parser.add_argument('--dipole-response-sign', type=float, default=1.0,
                       help='Sign of electric field (+1 for plus clone, -1 for minus clone) (default: +1.0)')
    
    
    return parser

def parse_replica_specification(replica_spec: str) -> List[int]:
    """Parse replica specification string into list of replica numbers."""
    try:
        if '-' in replica_spec:
            # Range specification
            start, end = map(int, replica_spec.split('-'))
            return list(range(start, end + 1))
        else:
            # Single replica
            return [int(replica_spec)]
    except ValueError:
        raise ValueError(f"Invalid replica specification: {replica_spec}")

def parse_json_or_file(json_str_or_path):
    """Parse JSON string or load from file path."""
    import json
    import os
    
    if json_str_or_path is None:
        return None
    
    # If it's already a list or dict, return it as-is
    if isinstance(json_str_or_path, (list, dict)):
        return json_str_or_path
    
    # Check if it's a file path
    if isinstance(json_str_or_path, str) and os.path.isfile(json_str_or_path):
        with open(json_str_or_path, 'r') as f:
            return json.load(f)
    elif isinstance(json_str_or_path, str):
        # Assume it's a JSON string
        try:
            return json.loads(json_str_or_path)
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON string or file path: {json_str_or_path}. Error: {e}")
    else:
        raise ValueError(f"Invalid type for JSON parameter: {type(json_str_or_path)}")

def main():
    """Main function."""
    parser = create_parser()
    args = parser.parse_args()
    
    # DEBUG: Print parsed LQG system ID parameters
    print(f"\n>>> DEBUG: After argparse:")
    print(f"    args.lqg_equilibrium_duration = {args.lqg_equilibrium_duration}")
    print(f"    args.lqg_prbs_n_unique_steps = {args.lqg_prbs_n_unique_steps}")
    print(f"    args.enable_pid_controller = {args.enable_pid_controller}")
    print(f"    args.lqg_bath_n_unique_steps = {args.lqg_bath_n_unique_steps}")
    print(f"    args.lqg_coupling_step_duration = {args.lqg_coupling_step_duration}")
    print(f"    args.lqg_coupling_equilibrium_duration = {args.lqg_coupling_equilibrium_duration}\n")
    
    print(" Universal Cavity MD Experiment Runner - Phase 3 Complete")
    print("============================================================")
    print(" Enhanced Coupling Variants: constant, step, periodic, exponential, square")
    print(" Enhanced Laser Drive: timing control, arbitrary k-vectors")
    print(" PI Feedback Controller: kinetic, LJ+Coulombic, harmonic temperature")
    print(" Full GPU/CPU Compatibility: identical results on both devices")
    print(" Backward Compatibility: all legacy flags still supported")
    print("============================================================")
    
    # Parse replica specification
    try:
        replicas = parse_replica_specification(args.replicas)
    except ValueError as e:
        print(f"Error: {e}")
        return 1
    
    print(f"Local execution: Replicas {replicas}")
    print()
    
    # Coupling mode validation and information
    coupling_modes = sum([args.periodic, args.laser, args.switch_time is not None])
    if coupling_modes > 1:
        print("WARNING: Multiple coupling modes specified. This may create complex dynamics.")
        print("  Periodic coupling + Laser drive + Switch-time = Very complex system")
        print("  Consider using only one primary coupling mode for cleaner analysis.")
        print()
    
    # Show enabled modes
    if args.periodic and not args.no_cavity:
        print(f"*** PERIODIC COUPLING ENABLED ***")
        print(f"Amplitude: {args.coupling:.6e} a.u.")
        print(f"Period: {args.period:.1f} ps ({1000.0/args.period:.1f} THz)")
        print(f"Phase offset: {args.phase_offset:.3f} rad")
        print(f"Start time: {args.start_time:.1f} ps")
        print(f"Formula: g(t) = {args.coupling:.3e} * sin(2π*t/{args.period} + {args.phase_offset:.3f})")
        print()
    
    if args.laser:
        print(f"*** LASER DRIVE ENABLED ***")
        print(f"Laser frequency: {args.laser_frequency:.1f} cm⁻¹")
        print(f"Laser amplitude: {args.laser_amplitude:.6e} a.u.")
        laser_freq_au = args.laser_frequency * 4.556335e-6  # Convert to a.u.
        print(f"Laser frequency (a.u.): {laser_freq_au:.6e}")
        print(f"Start time: {args.laser_start_time:.1f} ps")
        print(f"Formula: F_laser = {args.laser_amplitude:.3e} * cos({args.laser_frequency:.1f} cm⁻¹ * t)")
        print()
    
    
    if args.mech_periodic:
        print(f"*** MECHANICAL CAVITY MODULATION ENABLED ***")
        print(f"Mechanical frequency: {args.mech_frequency:.1f} cm⁻¹")
        print(f"Mechanical amplitude: {args.mech_magnitude:.6e} a.u.")
        mech_freq_au = args.mech_frequency * 4.556335e-6  # Convert to a.u.
        print(f"Mechanical frequency (a.u.): {mech_freq_au:.6e}")
        print(f"Cavity length formula: L(t) = L₀ + {args.mech_magnitude:.3e} * sin(2π × {args.mech_frequency:.1f} cm⁻¹ × t)")
        print(f"This modulates both cavity frequency ω(t) ∝ 1/L(t) and coupling g(t) ∝ 1/L(t)")
        print()
    
    if args.switch_time is not None and not args.no_cavity:
        print(f"*** SWITCH-TIME COUPLING ENABLED ***")
        print(f"Switch time: {args.switch_time:.1f} ps")
        if args.decay_time_constant is not None:
            print(f"Decay time constant: {args.decay_time_constant:.1f} ps")
            print(f"Formula: g(t) = 0 (t < {args.switch_time}), g(t) = {args.coupling:.3e} * exp(-(t-{args.switch_time})/{args.decay_time_constant}) (t ≥ {args.switch_time})")
        else:
            print(f"Formula: g(t) = 0 (t < {args.switch_time}), g(t) = {args.coupling:.3e} (t ≥ {args.switch_time})")
        print()
    
    # Enhanced validation for Phase 3 features
    if args.enable_empirical_feedback and args.empirical_data_file is None:
        parser.error("--empirical-data-file is required when --enable-empirical-feedback is used")
    
    # Exponential decay validation
    if args.coupling_type == 'exponential':
        if args.exponential_decay_time <= 0:
            parser.error("Exponential decay time must be positive")
        if args.exponential_turn_on_time < 0:
            parser.error("Exponential turn-on time must be non-negative")
    
    # Square wave validation
    if args.coupling_type == 'square':
        if args.square_period <= 0:
            parser.error("Square wave period must be positive")
        if not (0.0 <= args.square_duty_cycle <= 1.0):
            parser.error("Square wave duty cycle must be between 0.0 and 1.0")
        if args.square_start_time < 0:
            parser.error("Square wave start time must be non-negative")
    
    # Legacy coupling mode parameters
    if args.periodic and args.period <= 0:
        parser.error("Period must be positive")
    
    if args.laser and args.laser_frequency <= 0:
        parser.error("Laser frequency must be positive")
    
    if args.laser and args.laser_amplitude <= 0:
        parser.error("Laser amplitude must be positive")
    
    # Enhanced laser validation
    if args.laser:
        if args.laser_start_time < 0:
            parser.error("Laser start time must be non-negative")
        if args.laser_stop_time is not None and args.laser_stop_time <= args.laser_start_time:
            parser.error("Laser stop time must be greater than start time")
    
    # Validate mechanical cavity parameters
    if args.mech_periodic and args.mech_frequency <= 0:
        parser.error("Mechanical frequency must be positive")
    
    if args.mech_periodic and args.mech_magnitude <= 0:
        parser.error("Mechanical magnitude must be positive")
    
    # Validate compatibility combinations
    if args.laser and args.mech_periodic:
        print("INFO: Both laser drive and mechanical modulation enabled - this creates a complex multi-frequency system")
    
    if args.periodic and args.mech_periodic:
        print("INFO: Both periodic coupling and mechanical modulation enabled - this creates complex coupled oscillations")
    
    if args.periodic and args.laser:
        print("INFO: Both periodic coupling and laser drive enabled - this creates a complex multi-mode forcing system")
    
    # Configuration summary
    print(f"Simulation Configuration:")
    print(f"  Cavity coupling: {'Disabled' if args.no_cavity else 'Enabled'}")
    if not args.no_cavity:
        # Determine effective coupling type (new system or legacy)
        effective_coupling_type = args.coupling_type
        
        # Override with legacy flags for backward compatibility
        if args.periodic:
            effective_coupling_type = 'periodic'
        elif args.switch_time is not None:
            effective_coupling_type = 'step'
        
        print(f"     Coupling type: {effective_coupling_type.upper()}")
        
        if effective_coupling_type == 'constant':
            print(f"    Coupling strength: {args.lambda_coupling:.6e} (dimensionless)")
            print(f"    Coupling strength: {args.lambda_coupling * args.frequency:.6e} a.u.")
        elif effective_coupling_type == 'step':
            print(f"    Coupling strength: {args.lambda_coupling:.6e} (dimensionless)")
            print(f"    Coupling strength: {args.lambda_coupling * args.frequency:.6e} a.u.")
            switch_time = args.switch_time if args.switch_time is not None else 0.0
            print(f"    Switch time: {switch_time:.1f} ps")
            if args.decay_time_constant is not None:
                print(f"    Decay time: {args.decay_time_constant:.1f} ps")
            if args.step_turn_off_time is not None:
                print(f"    Turn-off time: {args.step_turn_off_time:.1f} ps")
                
        elif effective_coupling_type == 'periodic':
            print(f"    Lambda coupling: {args.lambda_coupling:.6e} (dimensionless)")
            print(f"    Amplitude: {args.lambda_coupling * args.frequency:.6e} a.u.")
            period = args.periodic_period if args.coupling_type == 'periodic' else args.period
            phase_offset = args.periodic_phase_offset if args.coupling_type == 'periodic' else args.phase_offset
            start_time = args.periodic_start_time if args.coupling_type == 'periodic' else args.start_time
            print(f"    Period: {period:.1f} ps ({1000.0/period:.1f} THz)")
            print(f"    Phase offset: {phase_offset:.3f} rad")
            print(f"    Start time: {start_time:.1f} ps")
            if args.periodic_stop_time is not None:
                print(f"    Stop time: {args.periodic_stop_time:.1f} ps")
                
        elif effective_coupling_type == 'exponential':
            amplitude = args.exponential_amplitude if args.exponential_amplitude is not None else args.lambda_coupling
            print(f"    Lambda coupling: {args.lambda_coupling:.6e} (dimensionless)")
            print(f"    Lambda coupling: {args.lambda_coupling * args.frequency:.6e} a.u.")
            print(f"    Amplitude: {amplitude:.6e} a.u.")
            print(f"    Decay time: {args.exponential_decay_time:.1f} ps")
            print(f"    Turn-on time: {args.exponential_turn_on_time:.1f} ps")
            if args.exponential_turn_off_time is not None:
                print(f"    Turn-off time: {args.exponential_turn_off_time:.1f} ps")
                
        elif effective_coupling_type == 'square':
            print(f"    Lambda coupling: {args.lambda_coupling:.6e} (dimensionless)")
            print(f"    Lambda coupling: {args.lambda_coupling * args.frequency:.6e} a.u.")
            print(f"    Period: {args.square_period:.1f} ps")
            print(f"    Duty cycle: {args.square_duty_cycle:.1%}")
            print(f"    Phase offset: {args.square_phase_offset:.3f} rad")
            print(f"    Start time: {args.square_start_time:.1f} ps")
            if args.square_stop_time is not None:
                print(f"    Stop time: {args.square_stop_time:.1f} ps")
                
        elif effective_coupling_type == 'decaying_square':
            print(f"    Lambda coupling: {args.lambda_coupling:.6e} (dimensionless)")
            print(f"    Lambda coupling: {args.lambda_coupling * args.frequency:.6e} a.u.")
            print(f"    Period: {args.decaying_square_period:.1f} ps")
            print(f"    Duty cycle: {args.decaying_square_duty_cycle:.1%}")
            print(f"    Phase offset: {args.decaying_square_phase_offset:.3f} rad")
            print(f"    Start time: {args.decaying_square_start_time:.1f} ps")
            if args.decaying_square_stop_time is not None:
                print(f"    Stop time: {args.decaying_square_stop_time:.1f} ps")
            print(f"    Decay rate: {args.decaying_square_decay_rate:.1%} per period")
            print(f"    Minimum amplitude: {args.decaying_square_minimum_amplitude:.2e} a.u.")
            # Calculate amplitude after a few periods for preview
            periods_preview = 5
            amplitude_after_periods = args.lambda_coupling * (1 - args.decaying_square_decay_rate) ** periods_preview
            print(f"    Amplitude after {periods_preview} periods: {amplitude_after_periods:.2e} a.u. ({amplitude_after_periods/args.coupling:.1%} of initial)")
        
        print(f"    Frequency: {args.frequency:.1f} cm⁻¹")
        print(f"    Finite-q mode: {args.finite_q}")
        print(f"    Cavity thermostat: {args.cavity_bath}")
        
        # Enhanced laser information
        if args.laser:
            print(f"     Enhanced laser drive:")
            print(f"      Frequency: {args.laser_frequency:.1f} cm⁻¹")
            print(f"      Amplitude: {args.laser_amplitude:.6e} a.u.")
            print(f"      Start time: {args.laser_start_time:.1f} ps")
            if args.laser_stop_time is not None:
                print(f"      Stop time: {args.laser_stop_time:.1f} ps")
            else:
                print(f"      Stop time: None (runs indefinitely)")
            if args.laser_kvector is not None:
                print(f"      k-vector: {args.laser_kvector}")
        
        # Legacy mechanical modulation
        if args.mech_periodic:
            print(f"     Mechanical modulation: {args.mech_frequency:.1f} cm⁻¹, amplitude {args.mech_magnitude:.6e} a.u.")
            
    # Gradient descent feedback information
    if args.enable_gd_feedback:
        print(f"     Gradient descent feedback controller:")
        print(f"      Target temperature: {args.gd_target_temperature:.1f} K")
        print(f"      Method: {args.gd_method}")
        print(f"      Turn-on time: {args.gd_turn_on_time:.1f} ps")
        if args.gd_turn_off_time is not None:
            print(f"      Turn-off time: {args.gd_turn_off_time:.1f} ps")
        else:
            print(f"      Turn-off time: None (runs until end)")
        print(f"      Time constant: {args.gd_time_constant:.1f} ps")
        print(f"      Apply to: {args.gd_apply_to} thermostat(s)")
        print(f"      Update interval: {args.gd_update_interval:.1f} ps")
        alpha = 0.0001 / args.gd_time_constant  # MD timestep / time constant
        print(f"      Learning rate α: {alpha:.6f}")
        T_max_str = f"{args.gd_T_max:.1f}" if args.gd_T_max is not None else "∞"
        print(f"      Temperature limits: [{args.gd_T_min:.1f}, {T_max_str}] K")
    
    # Quench controller information
    if args.enable_quench_controller:
        print(f"     Quench controller:")
        print(f"      Initial temperature: {args.quench_initial_temperature:.1f} K")
        print(f"      Target temperature: {args.quench_target_temperature:.1f} K")
        print(f"      ΔT = {args.quench_target_temperature - args.quench_initial_temperature:+.1f} K")
        print(f"      Quench time: {args.quench_time:.1f} ps")
        print(f"      Apply to: {args.quench_apply_to} thermostat(s)")
        quench_type = "Cooling" if args.quench_target_temperature < args.quench_initial_temperature else "Heating"
        print(f"      Quench type: {quench_type}")
    
    # Temperature tracker information
    if args.enable_temp_tracker:
        print(f"     Comprehensive temperature tracker:")
        print(f"      Output period: {args.temp_tracker_output_period_ps:.1f} ps")
        print(f"      Tracking: kinetic, harmonic fictive, LJ+Coul fictive, cavity bath, molecular bath")
        if args.temp_tracker_empirical_data_file:
            print(f"      Empirical data: {args.temp_tracker_empirical_data_file}")
        else:
            print(f"      Empirical data: None (LJ+Coul fictive temp will be 0)")
    
    # Molecular temperature decomposition information
    if args.enable_molecular_temps:
        print(f"     Molecular temperature decomposition:")
        print(f"      Output period: {args.molecular_temps_output_period_ps:.1f} ps")
        print(f"      Tracking: T_trans, T_rot, T_vib (for O-O and N-N dimers)")
    
    # Dipole moment FDR information
    if args.enable_dipole_fdr:
        print(f"     Dipole moment FDR analysis:")
        print(f"      Output period: {args.dipole_fdr_output_period_ps:.1f} ps")
        print(f"      Max correlation time: {args.dipole_fdr_max_correlation_time_ps:.1f} ps")
        print(f"      Field direction: {args.dipole_fdr_field_direction}")
        print(f"      Exclude cavity: {args.dipole_fdr_exclude_cavity}")
        if args.enable_dipole_response:
            print(f"      Response mode: fork-and-clone FDR")
            print(f"      Field strength: {args.dipole_response_field_strength:.2e}")
            print(f"      Clone sign: {args.dipole_response_sign:+.0f}")
        else:
            print(f"      Response mode: autocorrelation only")
            
    print(f"  Molecular thermostat: {args.molecular_bath}")
    print(f"  Temperature: {args.temperature:.1f} K")
    print(f"  Runtime: {args.runtime:.1f} ps")
    print(f"  Device: {args.device}")
    if args.device == 'GPU':
        print(f"    GPU ID: {args.gpu_id}")
    print(f"  Random seed: {args.seed}")
    print(f"  Velocity restart: {'Disabled' if args.no_restart_velocities else 'Enabled'} - {'using existing velocities' if args.no_restart_velocities else 'thermalizing velocities'}")
    print(f"  Momentum zeroing: {'Enabled' if args.zero_momentum else 'Disabled'}")
    if args.enable_empirical_feedback:
        print(f"  Empirical feedback: Enabled")
        print(f"    Data file: {args.empirical_data_file}")
        print(f"    Energy component: {args.feedback_energy_component}")
        print(f"    Apply to: {args.feedback_apply_to}")
        print(f"    Averaging window: {args.feedback_averaging_window_ps:.1f} ps")
        print(f"    Update interval: {args.feedback_update_interval_ps:.1f} ps")
    
    print()
    print(f"Starting execution for {len(replicas)} replica(s)...")
    print("============================================================")
    print()
    
    # Run simulations for all replicas
    success_count = 0
    for replica in replicas:
        print(f"Running replica {replica}...")
        
        success = run_single_experiment(
            molecular_thermo=args.molecular_bath,
            cavity_thermo=args.cavity_bath,
            finite_q=args.finite_q,
            coupling=args.coupling,
            lambda_coupling=args.lambda_coupling,
            temperature=args.temperature,
            frequency=args.frequency,
            runtime_ps=args.runtime,
            input_gsd=args.input_gsd,
            frame=args.frame,
            
            # Enhanced coupling variant parameters
            coupling_variant_type=args.coupling_type,
            # Enhanced step coupling parameters
            step_turn_off_time_ps=args.step_turn_off_time,
            # Enhanced periodic coupling parameters
            periodic_coupling=args.periodic,
            periodic_period_ps=args.periodic_period if args.coupling_type == 'periodic' else args.period,
            periodic_phase_offset=args.periodic_phase_offset if args.coupling_type == 'periodic' else args.phase_offset,
            periodic_start_time_ps=args.periodic_start_time if args.coupling_type == 'periodic' else args.start_time,
            periodic_stop_time_ps=args.periodic_stop_time,
            # Exponential decay parameters
            exponential_amplitude=args.exponential_amplitude,
            exponential_decay_time_ps=args.exponential_decay_time,
            exponential_turn_on_time_ps=args.exponential_turn_on_time,
            exponential_turn_off_time_ps=args.exponential_turn_off_time,
            # Square wave parameters
            square_period_ps=args.square_period,
            square_duty_cycle=args.square_duty_cycle,
            square_phase_offset=args.square_phase_offset,
            square_start_time_ps=args.square_start_time,
            square_stop_time_ps=args.square_stop_time,
            # Decaying square wave parameters
            decaying_square_period_ps=args.decaying_square_period,
            decaying_square_duty_cycle=args.decaying_square_duty_cycle,
            decaying_square_phase_offset=args.decaying_square_phase_offset,
            decaying_square_start_time_ps=args.decaying_square_start_time,
            decaying_square_stop_time_ps=args.decaying_square_stop_time,
            decaying_square_decay_rate=args.decaying_square_decay_rate,
            decaying_square_minimum_amplitude=args.decaying_square_minimum_amplitude,
            # Adaptive square wave parameters
            adaptive_square_period_ps=args.adaptive_square_period,
            adaptive_square_duty_cycle=args.adaptive_square_duty_cycle,
            adaptive_square_phase_offset=args.adaptive_square_phase_offset,
            adaptive_square_start_time_ps=args.adaptive_square_start_time,
            adaptive_square_stop_time_ps=args.adaptive_square_stop_time,
            adaptive_square_min_amplitude=args.adaptive_square_min_amplitude,
            adaptive_square_max_amplitude=args.adaptive_square_max_amplitude,
            # Exponential wave parameters
            exp_period_ps=args.exp_period,
            exp_tau_ps=args.exp_tau,
            exp_start_time_ps=args.exp_start_time,
            exp_stop_time_ps=args.exp_stop_time,
            exp_adaptive=args.exp_adaptive,
            # Composite coupling parameters
            composite_sinusoid_amplitude=args.composite_sinusoid_amplitude,
            composite_sinusoid_period=args.composite_sinusoid_period,
            composite_sinusoid_phase=args.composite_sinusoid_phase,
            composite_sinusoid_start_time=args.composite_sinusoid_start_time,
            composite_sinusoid_stop_time=args.composite_sinusoid_stop_time,
            composite_square_amplitude=args.composite_square_amplitude,
            composite_square_period=args.composite_square_period,
            composite_square_duty_cycle=args.composite_square_duty_cycle,
            composite_square_start_time=args.composite_square_start_time,
            composite_square_stop_time=args.composite_square_stop_time,
            composite_square_adaptive=args.composite_square_adaptive,
            composite_max_amplitude=args.composite_max_amplitude,
            # Auto-stop coupling parameters
            enable_auto_stop=args.enable_auto_stop,
            auto_stop_tol=args.auto_stop_tol,
            auto_stop_window=args.auto_stop_window,
            # Enhanced laser drive parameters
            laser_enabled=args.laser,
            laser_frequency_cm1=args.laser_frequency,
            laser_amplitude=args.laser_amplitude,
            laser_start_time_ps=args.laser_start_time,
            laser_stop_time_ps=args.laser_stop_time,
            laser_kvector=args.laser_kvector,
            # Gradient descent feedback controller parameters
            enable_gd_feedback=args.enable_gd_feedback,
            gd_target_temperature=args.gd_target_temperature,
            gd_turn_on_time_ps=args.gd_turn_on_time,
            gd_turn_off_time_ps=args.gd_turn_off_time,
            gd_temperature_method=args.gd_method,
            gd_update_interval_ps=args.gd_update_interval,
            gd_time_constant_ps=args.gd_time_constant,
            gd_apply_to=args.gd_apply_to,
            gd_T_min=args.gd_T_min,
            gd_T_max=args.gd_T_max,
            gd_disable_effective_temp=args.gd_no_effective_temp,
            # Multi-signal error function parameters
            gd_enable_multi_signal=args.gd_enable_multi_signal,
            gd_weight_system_target=args.gd_weight_system_target,
            gd_weight_bath_target=args.gd_weight_bath_target,
            gd_weight_system_bath=args.gd_weight_system_bath,
            # Dual independent feedback controller parameters
            enable_dual_feedback=args.enable_dual_feedback,
            dual_cavity_method=args.dual_cavity_method,
            dual_molecular_method=args.dual_molecular_method,
            dual_cavity_target_temperature=args.dual_cavity_target,
            dual_molecular_target_temperature=args.dual_molecular_target,
            dual_cavity_time_constant_ps=args.dual_cavity_time_constant,
            dual_molecular_time_constant_ps=args.dual_molecular_time_constant,
            dual_turn_on_time_ps=args.dual_turn_on_time,
            dual_turn_off_time_ps=args.dual_turn_off_time,
            dual_update_interval_ps=args.dual_update_interval,
            dual_cavity_T_min=args.dual_cavity_T_min,
            dual_cavity_T_max=args.dual_cavity_T_max,
            dual_molecular_T_min=args.dual_molecular_T_min,
            dual_molecular_T_max=args.dual_molecular_T_max,
            dual_cavity_dynamic_target=args.dual_cavity_dynamic_target,
            dual_molecular_dynamic_target=args.dual_molecular_dynamic_target,
            dual_cavity_integral_time_constant_ps=args.dual_cavity_integral_time_constant,
            dual_molecular_integral_time_constant_ps=args.dual_molecular_integral_time_constant,
            # Sinusoidal bath temperature controller parameters
            enable_sinusoidal_bath=args.enable_sinusoidal_bath,
            sinusoidal_bath_period_ps=args.sinusoidal_bath_period,
            sinusoidal_bath_amplitude_scale=args.sinusoidal_bath_amplitude_scale,
            sinusoidal_bath_phase_offset=args.sinusoidal_bath_phase_offset,
            sinusoidal_bath_target_temperature=args.sinusoidal_bath_target,
            sinusoidal_bath_dynamic_target=args.sinusoidal_bath_dynamic_target,
            sinusoidal_bath_turn_on_time_ps=args.sinusoidal_bath_turn_on_time,
            sinusoidal_bath_turn_off_time_ps=args.sinusoidal_bath_turn_off_time,
            sinusoidal_bath_update_interval_ps=args.sinusoidal_bath_update_interval,
            sinusoidal_bath_apply_to=args.sinusoidal_bath_apply_to,
            sinusoidal_bath_T_min=args.sinusoidal_bath_T_min,
            sinusoidal_bath_T_max=args.sinusoidal_bath_T_max,
            sinusoidal_bath_empirical_data_file=args.sinusoidal_bath_empirical_data_file,
            sinusoidal_bath_amplitude_update_interval_ps=args.sinusoidal_bath_amplitude_update_interval,
            sinusoidal_bath_amplitude_temperature_method=args.sinusoidal_bath_amplitude_temperature_method,
            sinusoidal_bath_adaptive_range_mode=args.sinusoidal_bath_adaptive_range_mode,
            # Adaptive bath temperature controller parameters
            enable_adaptive_bath=args.enable_adaptive_bath,
            adaptive_bath_amplitude_scale=args.adaptive_bath_amplitude_scale,
            adaptive_bath_time_constant_ps=args.adaptive_bath_time_constant,
            adaptive_bath_target_temperature=args.adaptive_bath_target_temperature,
            adaptive_bath_dynamic_target=args.adaptive_bath_dynamic_target,
            adaptive_bath_turn_on_time_ps=args.adaptive_bath_turn_on_time,
            adaptive_bath_turn_off_time_ps=args.adaptive_bath_turn_off_time,
            adaptive_bath_update_interval_ps=args.adaptive_bath_update_interval,
            adaptive_bath_apply_to=args.adaptive_bath_apply_to,
            adaptive_bath_T_min=args.adaptive_bath_T_min,
            adaptive_bath_T_max=args.adaptive_bath_T_max,
            adaptive_bath_empirical_data_file=args.adaptive_bath_empirical_data_file,
            adaptive_bath_signal_temperature_method=args.adaptive_bath_signal_temperature_method,
            adaptive_bath_dynamic_target_temperature_method=args.adaptive_bath_dynamic_target_temperature_method,
            # Quench controller parameters
            enable_quench_controller=args.enable_quench_controller,
            quench_initial_temperature=args.quench_initial_temperature,
            quench_target_temperature=args.quench_target_temperature,
            quench_time_ps=args.quench_time,
            quench_apply_to=args.quench_apply_to,
            # Offset temperature controller parameters
            enable_offset_controller=args.enable_offset_controller,
            offset_temperature_method=args.offset_temperature_method,
            offset_temperature_offset_K=args.offset_temperature_offset_K,
            offset_turn_on_time_ps=args.offset_turn_on_time,
            offset_turn_off_time_ps=args.offset_turn_off_time,
            offset_update_interval_ps=args.offset_update_interval,
            offset_apply_to=args.offset_apply_to,
            offset_T_min=args.offset_T_min,
            offset_T_max=args.offset_T_max,
            # Harmonic bond reset parameters
            enable_harmonic_reset=args.enable_harmonic_reset,
            harmonic_reset_turn_on_time_ps=args.harmonic_reset_turn_on_time,
            harmonic_reset_temperature=args.harmonic_reset_temperature,
            harmonic_reset_seed=args.harmonic_reset_seed,
            # Differential equation controller parameters
            enable_diffeq_controller=args.enable_diffeq_controller,
            diffeq_temperature_method=args.diffeq_temperature_method,
            diffeq_time_constant_ps=args.diffeq_time_constant,
            diffeq_time_constant_auto=args.diffeq_time_constant_auto,
            diffeq_turn_on_time_ps=args.diffeq_turn_on_time,
            diffeq_turn_off_time_ps=args.diffeq_turn_off_time,
            diffeq_update_interval_ps=args.diffeq_update_interval,
            diffeq_apply_to=args.diffeq_apply_to,
            diffeq_T_min=args.diffeq_T_min,
            diffeq_T_max=args.diffeq_T_max,
            diffeq_rate_limit_K_per_ps=args.diffeq_rate_limit,
            diffeq_disable_bias_estimation=args.diffeq_disable_bias_estimation,
            # Exact-cancellation + PI control parameters
            diffeq_enable_pi_control=args.diffeq_enable_pi_control,
            diffeq_pi_rho=args.diffeq_pi_rho,
            diffeq_pi_epsilon=args.diffeq_pi_epsilon,  
            diffeq_pi_zeta=args.diffeq_pi_zeta,
            diffeq_relaxation_data_file=args.diffeq_relaxation_data_file,
            diffeq_filter_window=args.diffeq_filter_window,
            # Adaptive bias cancellation parameters for DiffEq controller
            diffeq_enable_bias_cancellation=args.diffeq_enable_bias_cancellation,
            diffeq_bias_tau_b_ps=args.diffeq_bias_tau_b,
            diffeq_bias_tau_b_auto=args.diffeq_bias_tau_b_auto,
            diffeq_bias_kappa=args.diffeq_bias_kappa,
            diffeq_bias_kappa_auto=args.diffeq_bias_kappa_auto,
            diffeq_bias_tau_b_prefactor=args.diffeq_bias_tau_b_prefactor,
            diffeq_bias_kappa_prefactor=args.diffeq_bias_kappa_prefactor,
            diffeq_bias_calibration_time_ps=args.diffeq_bias_calibration_time,
            # Simple Setpoint Controller parameters
            enable_simple_setpoint_controller=args.enable_simple_setpoint_controller,
            simple_setpoint_signal_method=args.simple_setpoint_signal_method,
            simple_setpoint_time_constant_ps=args.simple_setpoint_time_constant_ps,
            simple_setpoint_apply_to=args.simple_setpoint_apply_to,
            simple_setpoint_turn_on_time_ps=args.simple_setpoint_turn_on_time_ps,
            simple_setpoint_turn_off_time_ps=args.simple_setpoint_turn_off_time_ps,
            simple_setpoint_update_interval_ps=args.simple_setpoint_update_interval_ps,
            simple_setpoint_T_min=args.simple_setpoint_T_min,
            simple_setpoint_T_max=args.simple_setpoint_T_max,
            simple_setpoint_output_file=args.simple_setpoint_output_file,
            simple_setpoint_console_output_period_ps=args.simple_setpoint_console_output_period_ps,
            # Adaptive MPC controller parameters
            enable_adaptive_mpc_controller=args.enable_adaptive_mpc_controller,
            adaptive_mpc_target_temperature=args.adaptive_mpc_target_temperature,
            adaptive_mpc_turn_on_time_ps=args.adaptive_mpc_turn_on_time_ps,
            adaptive_mpc_turn_off_time_ps=args.adaptive_mpc_turn_off_time_ps,
            adaptive_mpc_system_id_duration_ps=args.adaptive_mpc_system_id_duration_ps,
            adaptive_mpc_system_id_step_duration_ps=args.adaptive_mpc_system_id_step_duration_ps,
            adaptive_mpc_system_id_seed=args.adaptive_mpc_system_id_seed,
            adaptive_mpc_update_interval_ps=args.adaptive_mpc_update_interval_ps,
            adaptive_mpc_prediction_horizon=args.adaptive_mpc_prediction_horizon,
            adaptive_mpc_control_horizon=args.adaptive_mpc_control_horizon,
            adaptive_mpc_output_weight=args.adaptive_mpc_output_weight,
            adaptive_mpc_control_effort_weight=args.adaptive_mpc_control_effort_weight,
            adaptive_mpc_rate_penalty_weight=args.adaptive_mpc_rate_penalty_weight,
            adaptive_mpc_lambda_min=args.adaptive_mpc_lambda_min,
            adaptive_mpc_lambda_max=args.adaptive_mpc_lambda_max,
            adaptive_mpc_T_bath_min=args.adaptive_mpc_T_bath_min,
            adaptive_mpc_T_bath_max=args.adaptive_mpc_T_bath_max,
            adaptive_mpc_delta_lambda_max=args.adaptive_mpc_delta_lambda_max,
            adaptive_mpc_delta_T_bath_max=args.adaptive_mpc_delta_T_bath_max,
            adaptive_mpc_apply_to=args.adaptive_mpc_apply_to,
            adaptive_mpc_rls_forgetting_factor=args.adaptive_mpc_rls_forgetting_factor,
            adaptive_mpc_rls_initial_covariance=args.adaptive_mpc_rls_initial_covariance,
            adaptive_mpc_model_update_interval=args.adaptive_mpc_model_update_interval,
            adaptive_mpc_output_file=args.adaptive_mpc_output_file,
            adaptive_mpc_console_output_period_ps=args.adaptive_mpc_console_output_period_ps,
            adaptive_mpc_regularization_param=args.adaptive_mpc_regularization_param,
            adaptive_mpc_use_scaling=args.adaptive_mpc_use_scaling,
            adaptive_mpc_debug_mode=args.adaptive_mpc_debug_mode,
            # PID controller parameters
            enable_pid_controller=args.enable_pid_controller,
            pid_signal_choice=args.pid_signal_choice,
            pid_target_temperature=args.pid_target_temperature,
            pid_self_loop=args.pid_self_loop,
            pid_Kp=args.pid_Kp,
            pid_Ti=args.pid_Ti,
            pid_Td=args.pid_Td,
            pid_auto_tune=args.pid_auto_tune,
            pid_auto_tune_step_size=args.pid_auto_tune_step_size,
            pid_auto_tune_duration_ps=args.pid_auto_tune_duration_ps,
            pid_turn_on_time_ps=args.pid_turn_on_time_ps,
            pid_turn_off_time_ps=args.pid_turn_off_time_ps,
            pid_update_interval_ps=args.pid_update_interval_ps,
            pid_apply_to=args.pid_apply_to,
            pid_T_min=args.pid_T_min,
            pid_T_max=args.pid_T_max,
            pid_rate_limit_K_per_ps=args.pid_rate_limit_K_per_ps,
            pid_derivative_filter_N=args.pid_derivative_filter_N,
            pid_enable_anti_windup=args.pid_enable_anti_windup,
            pid_output_file=args.pid_output_file,
            pid_console_output_period_ps=args.pid_console_output_period_ps,
            # BathPI controller parameters
            enable_bath_pi_controller=args.enable_bath_pi_controller,
            bath_pi_apply_to=args.bath_pi_apply_to,
            bath_pi_K_p_molecular=args.bath_pi_K_p_molecular,
            bath_pi_K_i_molecular=args.bath_pi_K_i_molecular,
            bath_pi_K_T_molecular=args.bath_pi_K_T_molecular,
            bath_pi_K_p_cavity=args.bath_pi_K_p_cavity,
            bath_pi_K_i_cavity=args.bath_pi_K_i_cavity,
            bath_pi_K_T_cavity=args.bath_pi_K_T_cavity,
            bath_pi_filter_window_ps=args.bath_pi_filter_window_ps,
            bath_pi_flux_source=args.bath_pi_flux_source,
            bath_pi_anti_windup_alpha=args.bath_pi_anti_windup_alpha,
            bath_pi_enable_feedforward=args.bath_pi_enable_feedforward,
            bath_pi_T_nominal=args.bath_pi_T_nominal,
            bath_pi_feedforward_tau_ps=args.bath_pi_feedforward_tau_ps,
            bath_pi_turn_on_time_ps=args.bath_pi_turn_on_time_ps,
            bath_pi_turn_off_time_ps=args.bath_pi_turn_off_time_ps,
            bath_pi_update_interval_ps=args.bath_pi_update_interval_ps,
            bath_pi_T_min=args.bath_pi_T_min,
            bath_pi_T_max=args.bath_pi_T_max,
            bath_pi_rate_limit_K_per_ps=args.bath_pi_rate_limit_K_per_ps,
            bath_pi_output_file=args.bath_pi_output_file,
            bath_pi_relaxation_data_file=args.bath_pi_relaxation_data_file,
            # LQR controller parameters
            enable_lqr_controller=args.enable_lqr_controller,
            lqr_signal_method=args.lqr_signal_method,
            lqr_hot_method=args.lqr_hot_method,
            lqr_target_temperature=args.lqr_target_temperature,
            lqr_dynamic_target=args.lqr_dynamic_target,
            lqr_dynamic_target_method=args.lqr_dynamic_target_method,
            lqr_weight_signal=args.lqr_weight_signal,
            lqr_weight_hot=args.lqr_weight_hot,
            lqr_weight_bath=args.lqr_weight_bath,  # Deprecated: not used in LQG 2D formulation
            lqr_weight_integral=args.lqr_weight_integral,
            lqr_control_effort=args.lqr_control_effort,
            lqr_process_noise_signal=args.lqr_process_noise_signal,
            lqr_process_noise_hot=args.lqr_process_noise_hot,
            lqr_measurement_noise_signal=args.lqr_measurement_noise_signal,
            lqr_measurement_noise_hot=args.lqr_measurement_noise_hot,
            lqr_system_id_mode=args.lqr_system_id_mode,
            lqr_system_id_temp_K=args.lqr_system_id_temp,
            lqr_system_id_duration_ps=args.lqr_system_id_duration,
            lqr_system_id_file=args.lqr_system_id_file,
            lqr_periodic_system_id=args.lqr_periodic_system_id,
            lqr_periodic_system_id_interval_ps=args.lqr_periodic_system_id_interval,
            # EKF adaptation
            lqr_use_ekf_adaptation=args.lqr_use_ekf_adaptation,
            lqr_ekf_update_interval=args.lqr_ekf_update_interval,
            lqr_ekf_process_noise_param=args.lqr_ekf_process_noise_param,
            lqr_ekf_initial_covariance_param=args.lqr_ekf_initial_covariance_param,
            lqr_adaptive_lqr_threshold=args.lqr_adaptive_lqr_threshold,
            # Gain scheduling
            lqr_enable_gain_scheduling=args.lqr_enable_gain_scheduling,
            lqr_gain_schedule_far_threshold=args.lqr_gain_schedule_far_threshold,
            lqr_gain_schedule_near_threshold=args.lqr_gain_schedule_near_threshold,
            # T_h low-pass filter
            lqr_th_filter_enabled=args.lqr_th_filter_enabled,
            lqr_th_filter_time_constant=args.lqr_th_filter_time_constant,
            # Gentle startup
            lqr_gentle_startup_steps=args.lqr_gentle_startup_steps,
            lqr_gentle_startup_min_authority=args.lqr_gentle_startup_min_authority,
            # Kinetic temperature tracking (3D state)
            lqr_track_kinetic_temp=args.lqr_track_kinetic_temp,
            lqr_weight_kinetic=args.lqr_weight_kinetic,
            lqr_process_noise_kinetic=args.lqr_process_noise_kinetic,
            lqr_measurement_noise_kinetic=args.lqr_measurement_noise_kinetic,
            # Cross-coupling weights for thermal equilibration
            lqr_cross_coupling_signal_kinetic=args.lqr_cross_coupling_signal_kinetic,
            lqr_cross_coupling_signal_hot=args.lqr_cross_coupling_signal_hot,
            lqr_cross_coupling_hot_kinetic=args.lqr_cross_coupling_hot_kinetic,
            # Timing
            lqr_turn_on_time_ps=args.lqr_turn_on_time,
            lqr_turn_off_time_ps=args.lqr_turn_off_time,
            lqr_update_interval_ps=args.lqr_update_interval,
            lqr_T_min=args.lqr_T_min,
            lqr_T_max=args.lqr_T_max,
            lqr_apply_to=args.lqr_apply_to,
            lqr_output_file=args.lqr_output_file,
            lqr_empirical_data_file=args.lqr_empirical_data_file,
            # Adaptive LQI controller parameters
            lqr_controller_type=args.lqr_controller_type,
            lqr_tau_L_initial=args.lqr_tau_L_initial,
            lqr_tau_H_initial=args.lqr_tau_H_initial,
            lqr_k_initial=args.lqr_k_initial,
            lqr_tau_b=args.lqr_tau_b,
            lqr_q_common=args.lqr_q_common,
            lqr_q_diff=args.lqr_q_diff,
            lqr_q_eta_common=args.lqr_q_eta_common,
            lqr_q_eta_diff=args.lqr_q_eta_diff,
            lqr_process_noise_drift=args.lqr_process_noise_drift,
            lqr_rls_forgetting=args.lqr_rls_forgetting,
            lqr_rls_regularization=args.lqr_rls_regularization,
            lqr_rls_update_interval=args.lqr_rls_update_interval,
            lqr_max_control_rate=args.lqr_max_control_rate,
            lqr_integral_max_common=args.lqr_integral_max_common,
            lqr_integral_max_diff=args.lqr_integral_max_diff,
            lqr_theta_change_threshold=args.lqr_theta_change_threshold,
            # Legacy LQG pulse controller parameters (for backward compatibility)
            enable_lqg_controller=args.enable_lqg_controller,
            lqg_A_matrix=parse_json_or_file(args.lqg_A_matrix) if args.lqg_A_matrix else None,
            lqg_B_matrix=parse_json_or_file(args.lqg_B_matrix) if args.lqg_B_matrix else None,
            lqg_C_matrix=parse_json_or_file(args.lqg_C_matrix) if args.lqg_C_matrix else None,
            lqg_S_matrix=parse_json_or_file(args.lqg_S_matrix) if args.lqg_S_matrix else None,
            lqg_r_target=parse_json_or_file(args.lqg_r_target) if args.lqg_r_target else None,
            lqg_Qz_weights=parse_json_or_file(args.lqg_Qz_weights) if args.lqg_Qz_weights else None,
            lqg_Ru_weights=parse_json_or_file(args.lqg_Ru_weights) if args.lqg_Ru_weights else None,
            lqg_Qe_weights=parse_json_or_file(args.lqg_Qe_weights) if args.lqg_Qe_weights else None,
            lqg_W_noise=parse_json_or_file(args.lqg_W_noise) if args.lqg_W_noise else None,
            lqg_V_noise=parse_json_or_file(args.lqg_V_noise) if args.lqg_V_noise else None,
            lqg_g0_baseline=args.lqg_g0_baseline,
            lqg_dt=args.lqg_dt,
            lqg_delta_g_min=args.lqg_delta_g_min,
            lqg_delta_g_max=args.lqg_delta_g_max,
            lqg_T_bath_min=args.lqg_T_bath_min,
            lqg_T_bath_max=args.lqg_T_bath_max,
            lqg_rate_limit_K_per_ps=args.lqg_rate_limit,
            lqg_turn_on_time_ps=args.lqg_turn_on_time,
            lqg_turn_off_time_ps=args.lqg_turn_off_time,
            lqg_update_interval_ps=args.lqg_update_interval,
            lqg_console_output_period_ps=args.lqg_console_output_period,
            lqg_c_affine_term=parse_json_or_file(args.lqg_c_affine_term) if args.lqg_c_affine_term else None,
            lqg_temperature_methods=parse_json_or_file(args.lqg_temperature_methods) if args.lqg_temperature_methods else None,
            lqg_square_period_ps=args.lqg_square_period,
            lqg_square_duty_cycle=args.lqg_square_duty_cycle,
            lqg_square_phase_offset=args.lqg_square_phase_offset,
            lqg_square_start_time_ps=args.lqg_square_start_time,
            lqg_square_stop_time_ps=args.lqg_square_stop_time,
            lqg_square_g_low_level=args.lqg_square_g_low_level,
            # LQG system identification parameters
            # Note: When coupling_type='lqg_coupling', disable legacy LQG system and use LQG coupling controller parameters
            lqg_enable_system_id=(args.lqg_enable_system_id or not args.lqg_disable_system_id) and args.coupling_type != 'lqg_coupling',
            lqg_equilibrium_duration_ps=args.lqg_coupling_equilibrium_duration if args.coupling_type == 'lqg_coupling' else args.lqg_equilibrium_duration,
            lqg_prbs_n_unique_steps=args.lqg_prbs_n_unique_steps,
            lqg_prbs_amplitude=args.lqg_prbs_amplitude,
            lqg_bath_n_unique_steps=args.lqg_bath_n_unique_steps,
            lqg_bath_excitation_amplitude=args.lqg_bath_excitation_amplitude,
            lqg_system_id_output_file=args.lqg_system_id_output_file,
            # LQG coupling controller parameters (new coupling-type: lqg_coupling)
            lqg_coupling_target_temperature=args.lqg_coupling_target_temperature,
            lqg_coupling_update_interval_ps=args.lqg_coupling_update_interval,
            lqg_coupling_equilibrium_duration_ps=args.lqg_coupling_equilibrium_duration,
            lqg_coupling_step_duration_ps=args.lqg_coupling_step_duration,
            lqg_coupling_n_steps=args.lqg_coupling_n_steps,
            lqg_coupling_lambda_min=args.lqg_coupling_lambda_min,
            lqg_coupling_lambda_max=args.lqg_coupling_lambda_max,
            lqg_coupling_process_noise_std=args.lqg_coupling_process_noise,
            lqg_coupling_measurement_noise_std=args.lqg_coupling_measurement_noise,
            lqg_coupling_weight_signal=args.lqg_coupling_weight_signal,
            lqg_coupling_weight_harmonic=args.lqg_coupling_weight_harmonic,
            lqg_coupling_weight_kinetic=args.lqg_coupling_weight_kinetic,
            lqg_coupling_weight_bath=args.lqg_coupling_weight_bath,
            lqg_coupling_control_effort=args.lqg_coupling_control_effort,
            lqg_coupling_integral_gain=args.lqg_coupling_integral_gain,
            lqg_coupling_system_id_file=args.lqg_coupling_system_id_file,
            lqg_coupling_control_file=args.lqg_coupling_control_file,
            lqg_coupling_temperature_methods=parse_json_or_file(args.lqg_coupling_temperature_methods) if isinstance(args.lqg_coupling_temperature_methods, str) else args.lqg_coupling_temperature_methods,
            lqg_coupling_bath_temperature_method=args.lqg_coupling_bath_temperature_method,
            # Temperature tracker parameters
            enable_temp_tracker=args.enable_temp_tracker,
            temp_tracker_output_period_ps=args.temp_tracker_output_period_ps,
            temp_tracker_empirical_data_file=args.temp_tracker_empirical_data_file,
            # Molecular temperature decomposition parameters
            enable_molecular_temps=args.enable_molecular_temps,
            molecular_temps_output_period_ps=args.molecular_temps_output_period_ps,
            # Dipole moment FDR parameters
            enable_dipole_fdr=args.enable_dipole_fdr,
            dipole_fdr_output_period_ps=args.dipole_fdr_output_period_ps,
            dipole_fdr_max_correlation_time_ps=args.dipole_fdr_max_correlation_time_ps,
            dipole_fdr_field_direction=args.dipole_fdr_field_direction,
            dipole_fdr_exclude_cavity=args.dipole_fdr_exclude_cavity,
            enable_dipole_response=args.enable_dipole_response,
            dipole_response_field_strength=args.dipole_response_field_strength,
            dipole_response_sign=args.dipole_response_sign,
            # Mechanical cavity modulation parameters  
            mech_periodic=args.mech_periodic,
            mech_frequency_cm1=args.mech_frequency,
            mech_magnitude=args.mech_magnitude,
            # Switch-time parameters
            switch_time_ps=args.switch_time,
            decay_time_constant_ps=args.decay_time_constant,
            damping_ratio=args.damping_ratio,
            # Thermostat parameters
            molecular_tau=args.molecular_tau,
            cavity_tau=args.cavity_tau,
            # Timestep parameters
            fixed_timestep=args.fixed_timestep,
            timestep_fs=args.timestep,
            error_tolerance=args.error_tolerance,
            initial_fraction=args.initial_fraction,
            time_constant_ps=args.time_constant_ps,
            # Dynamic coupling detection parameters
            enable_dynamic_coupling_detection=not args.disable_dynamic_coupling_detection,
            coupling_change_threshold=args.coupling_change_threshold,
            # Output and tracking parameters
            enable_energy_tracker=args.enable_energy_tracker,
            energy_output_period_ps=args.energy_output_period_ps,
            fkt_output_period_ps=args.fkt_output_period_ps,
            gsd_output_period_ps=args.gsd_output_period_ps,
            console_output_period_ps=args.console_output_period_ps,
            enable_fkt=args.enable_fkt,
            fkt_kmag=args.fkt_kmag,
            fkt_wavevectors=args.fkt_wavevectors,
            fkt_ref_interval=args.fkt_ref_interval,
            fkt_max_refs=args.fkt_max_refs,
            enable_dipole_autocorr=args.enable_dipole_autocorr,
            dipole_ref_interval=args.dipole_ref_interval,
            dipole_max_refs=args.dipole_max_refs,
            dipole_output_period_ps=args.dipole_output_period_ps,
            max_energy_output_time=args.max_energy_output_time,
            # Hardware parameters
            device=args.device,
            gpu_id=args.gpu_id,
            truncate_gsd=args.truncate_gsd,
            seed=args.seed,
            restart_velocities=not args.no_restart_velocities,
            zero_momentum_enabled=args.zero_momentum,
            zero_momentum_period_ps=args.zero_momentum_period_ps,
            # Empirical feedback parameters
            enable_empirical_feedback=args.enable_empirical_feedback,
            empirical_data_file=args.empirical_data_file,
            feedback_output_period_ps=args.feedback_output_period_ps,
            feedback_energy_component=args.feedback_energy_component,
            feedback_apply_to=args.feedback_apply_to,
            feedback_averaging_window_ps=args.feedback_averaging_window_ps,
            feedback_update_interval_ps=args.feedback_update_interval_ps,
            feedback_T_min=args.feedback_T_min,
            feedback_T_max=args.feedback_T_max,
            feedback_turn_off_time_ps=args.feedback_turn_off_time_ps,
            feedback_use_direct_harmonic=args.feedback_use_direct_harmonic,
            feedback_auto_adjust_molecular_tau=args.feedback_auto_adjust_molecular_tau,
            enable_auto_stopping=args.feedback_auto_stopping,
            auto_stop_min_time_ps=args.auto_stop_min_time_ps,
            auto_stop_smoothing_window=args.auto_stop_smoothing_window,
            # Experiment metadata
            replica=replica,
            incavity=not args.no_cavity
        )
        
        if success:
            print(f" Replica {replica} completed successfully")
            success_count += 1
        else:
            print(f" Replica {replica} failed")
        
        print()
    
    print("============================================================")
    print(f"Execution completed: {success_count}/{len(replicas)} replicas successful")
    
    if success_count == len(replicas):
        print(" All simulations completed successfully!")
        return 0
    else:
        print(f"  {len(replicas) - success_count} simulations failed")
        return 1

if __name__ == '__main__':
    sys.exit(main())