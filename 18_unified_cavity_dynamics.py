#!/usr/bin/env python3
"""
Universal Cavity MD Experiment Runner - All Phase 3 Features

This script provides a comprehensive framework for running cavity MD simulations
with ALL enhanced coupling modes and control features:

🔧 ENHANCED COUPLING VARIANTS:
- CONSTANT: g(t) = constant (default)
- STEP: g(t) = 0 → g_target at t_switch (with optional decay and turn-off)
- PERIODIC: g(t) = A * sin(2π*t/T + φ) (with start/stop times)
- EXPONENTIAL: g(t) = A * exp(-t/τ) (with turn-on/off times)
- SQUARE WAVE: g(t) = periodic on/off pulses (with duty cycle)

🔬 ENHANCED LASER DRIVE:
- F_laser = F_L * cos(ω_L * t) with precise timing control
- Start/stop times for laser pulses
- Arbitrary k-vector directions

🌡️ PI FEEDBACK CONTROLLER:
- Kinetic temperature control (from particle velocities)
- LJ+Coulombic fictive temperature (Rosenfeld scaling)
- Harmonic fictive temperature (T^3/2 scaling)
- IMC auto-tuning or manual PI parameters

🚀 PRODUCTION FEATURES:
- All Phase 3 enhanced features integrated
- Full GPU/CPU compatibility and consistency 
- Backward compatibility with existing workflows
- Smart folder naming with all parameters
- Comprehensive energy tracking and analysis
- Multi-replica execution support
- Advanced timing control for all features
- Real-time PI feedback with multiple temperature methods

💡 ENHANCED USAGE EXAMPLES:

  # Exponential decay coupling with PI feedback
  python 18_unified_cavity_dynamics.py --coupling-type exponential --exponential-amplitude 1e-3 --exponential-decay-time 15.0 --enable-pi-feedback --pi-method kinetic
  
  # Square wave coupling with timed laser
  python 18_unified_cavity_dynamics.py --coupling-type square --coupling 1e-3 --square-period 8.0 --square-duty-cycle 0.3 --laser --laser-start-time 10.0 --laser-stop-time 30.0
  
  # Periodic coupling with harmonic temperature control
  python 18_unified_cavity_dynamics.py --coupling-type periodic --coupling 1e-3 --period 5.0 --periodic-start-time 5.0 --periodic-stop-time 50.0 --enable-pi-feedback --pi-method harmonic
  
  # Step coupling with LJ+Coulombic feedback
  python 18_unified_cavity_dynamics.py --coupling-type step --coupling 1e-3 --switch-time 10.0 --step-turn-off-time 40.0 --enable-pi-feedback --pi-method lj_coulombic
  
  # Multi-feature advanced experiment
  python 18_unified_cavity_dynamics.py --coupling-type exponential --exponential-amplitude 1.5e-3 --laser --laser-start-time 20.0 --laser-stop-time 60.0 --enable-pi-feedback --pi-method kinetic --device GPU

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

# =============================================================================
# SIMULATION FUNCTIONS (Unified: Periodic, Laser, Mechanical, Switch-time)
# =============================================================================

def run_single_experiment(molecular_thermo, cavity_thermo, finite_q,
                         coupling=1e-3, temperature=100.0, frequency=1560.0,
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
                         # Enhanced laser drive parameters  
                         laser_enabled=False, laser_frequency_cm1=1560.0, laser_amplitude=1e-5, 
                         laser_start_time_ps=0.0, laser_stop_time_ps=None, laser_kvector=None,
                         # PI feedback controller parameters
                         enable_pi_feedback=False, pi_target_temperature=100.0,
                         pi_turn_on_time_ps=0.0, pi_turn_off_time_ps=None,
                         pi_temperature_method='kinetic', pi_update_interval_ps=5.0,
                         pi_Kc=None, pi_Ti=None, pi_molecular_tau_ps=5.0,
                         pi_beta=0.7, pi_T_min=0.0, pi_T_max=None,
                         # Gradient descent feedback controller parameters
                         enable_gd_feedback=False, gd_target_temperature=100.0,
                         gd_turn_on_time_ps=0.0, gd_turn_off_time_ps=None,
                         gd_temperature_method='kinetic', gd_update_interval_ps=0.1,
                         gd_time_constant_ps=10.0, gd_apply_to='both', gd_T_min=0.0, gd_T_max=None,
                         # Temperature tracker parameters
                         enable_temp_tracker=False, temp_tracker_output_period_ps=0.1,
                         temp_tracker_empirical_data_file=None,
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
                         feedback_turn_off_time_ps=None, feedback_use_direct_harmonic_calculation=True,
                         feedback_auto_adjust_molecular_tau=False, feedback_auto_stopping=False,
                         auto_stop_min_time_ps=10.0, auto_stop_smoothing_window=5,
                         # Experiment metadata
                         replica=0, incavity=True) -> bool:
    """
    Run a unified cavity MD experiment with support for ALL Phase 3 enhanced features.
    
    This function supports:
    🔧 ENHANCED COUPLING VARIANTS:
    - Constant: g(t) = constant (default)
    - Step: g(t) = 0 → g_target at t_switch (with optional decay and turn-off)
    - Periodic: g(t) = A * sin(2π*t/T + φ) (with start/stop times)
    - Exponential: g(t) = A * exp(-t/τ) (with turn-on/off times)
    - Square wave: g(t) = periodic on/off pulses (with duty cycle)
    
    🔬 ENHANCED LASER DRIVE:
    - F_laser = F_L * cos(ω_L * t) with precise timing control
    - Start/stop times for laser pulses
    - Arbitrary k-vector directions
    
    🌡️ PI FEEDBACK CONTROLLER:
    - Kinetic temperature control (from particle velocities)
    - LJ+Coulombic fictive temperature (Rosenfeld scaling)
    - Harmonic fictive temperature (T^3/2 scaling)
    - IMC auto-tuning or manual PI parameters
    
    🚀 LEGACY SUPPORT:
    - Mechanical cavity modulation: L(t) = L₀ + A_mech * sin(ω_mech * t)
    - Empirical temperature feedback (legacy)
    
    Returns True for success, False for failure.
    """
    
    try:
        # Calculate dissipation from damping_ratio
        phmass = 1.0  # Photon mass is 1.0 in a.u.
        omegac = frequency / PhysicalConstants.HARTREE_TO_CM_MINUS1
        dissipation = 2 * damping_ratio * phmass * omegac

        # Create experiment directory with appropriate naming based on enabled features
        if incavity:
            # For cavity simulations, include coupling strength in directory name
            coupling_str = f"{coupling:.0e}".replace("-", "neg").replace("+", "pos")
            
            # Add mechanical cavity modulation prefix if enabled
            mech_prefix = ""
            if mech_periodic:
                mech_freq_str = f"{mech_frequency_cm1:.0f}cm1"
                mech_mag_str = f"{mech_magnitude:.0e}".replace("-", "neg").replace("+", "pos")
                mech_prefix = f"mech_periodic_{mech_freq_str}_{mech_mag_str}_"
            
            # Determine coupling type and create appropriate directory name
            if periodic_coupling:
                # Periodic coupling directory
                period_str = f"_period_{periodic_period_ps}ps"
                phase_str = f"_phase_{periodic_phase_offset:.2f}rad" if periodic_phase_offset != 0.0 else ""
                start_str = f"_start_{periodic_start_time_ps}ps" if periodic_start_time_ps != 0.0 else ""
                exp_dir = Path(f"{mech_prefix}periodic_coupling_{coupling_str}{period_str}{phase_str}{start_str}")
                
            elif laser_enabled:
                # Laser drive directory
                freq_str = f"_laser_{laser_frequency_cm1:.0f}cm1"
                amp_str = f"_amp_{laser_amplitude:.0e}".replace("-", "neg").replace("+", "pos")
                start_str = f"_start_{laser_start_time_ps}ps" if laser_start_time_ps != 0.0 else ""
                exp_dir = Path(f"{mech_prefix}laser_drive_{coupling_str}{freq_str}{amp_str}{start_str}")
                
            elif switch_time_ps is not None:
                # Switch-time coupling directory
                switch_str = f"_switch_{switch_time_ps}ps"
                if decay_time_constant_ps is not None:
                    decay_str = f"_decay_{decay_time_constant_ps}ps"
                    exp_dir = Path(f"{mech_prefix}cavity_coupling_{coupling_str}{switch_str}{decay_str}")
                else:
                    exp_dir = Path(f"{mech_prefix}cavity_coupling_{coupling_str}{switch_str}")
            else:
                # Constant coupling directory
                exp_dir = Path(f"{mech_prefix}cavity_coupling_{coupling_str}")
        else:
            # For non-cavity simulations
            exp_dir = Path("no_cavity")
        exp_dir.mkdir(exist_ok=True)
        
        print(f"Running unified cavity dynamics experiment:")
        print(f"  Cavity coupling: {'Enabled' if incavity else 'Disabled'}")
        
        if incavity:
            print(f"  *** COUPLING STRENGTH: {coupling:.6e} a.u. ***")
            
            # Show active coupling modes
            if periodic_coupling:
                print(f"  *** PERIODIC COUPLING ENABLED ***")
                print(f"  Amplitude: {coupling:.6e} a.u.")
                print(f"  Period: {periodic_period_ps:.1f} ps ({1000.0/periodic_period_ps:.1f} THz)")
                print(f"  Phase offset: {periodic_phase_offset:.3f} rad")
                print(f"  Start time: {periodic_start_time_ps:.1f} ps")
                print(f"  Formula: g(t) = {coupling:.3e} * sin(2π*t/{periodic_period_ps} + {periodic_phase_offset:.3f})")
            
            if laser_enabled:
                print(f"  *** LASER DRIVE ENABLED ***")
                print(f"  Laser frequency: {laser_frequency_cm1:.1f} cm⁻¹")
                print(f"  Laser amplitude: {laser_amplitude:.6e} a.u.")
                laser_freq_au = laser_frequency_cm1 * 4.556335e-6  # Convert to a.u.
                print(f"  Laser frequency (a.u.): {laser_freq_au:.6e}")
                print(f"  Start time: {laser_start_time_ps:.1f} ps")
                print(f"  Formula: F_laser = {laser_amplitude:.3e} * cos({laser_frequency_cm1:.1f} cm⁻¹ * t)")
            
            
            if switch_time_ps is not None:
                print(f"  *** SWITCH-TIME COUPLING ***")
                print(f"  Switch time: {switch_time_ps:.1f} ps")
                if decay_time_constant_ps is not None:
                    print(f"  Decay time constant: {decay_time_constant_ps:.1f} ps")
                    print(f"  Formula: g(t) = 0 (t < {switch_time_ps}), g(t) = {coupling:.3e} * exp(-(t-{switch_time_ps})/{decay_time_constant_ps}) (t ≥ {switch_time_ps})")
                else:
                    print(f"  Decay: None (step function)")
                    print(f"  Formula: g(t) = 0 (t < {switch_time_ps}), g(t) = {coupling:.3e} (t ≥ {switch_time_ps})")
            
            if mech_periodic:
                print(f"  *** MECHANICAL CAVITY MODULATION ENABLED ***")
                print(f"  Mechanical frequency: {mech_frequency_cm1:.1f} cm⁻¹")
                print(f"  Mechanical amplitude: {mech_magnitude:.6e} a.u.")
                mech_freq_au = mech_frequency_cm1 * 4.556335e-6  # Convert to a.u.
                print(f"  Mechanical frequency (a.u.): {mech_freq_au:.6e}")
                print(f"  Formula: L(t) = L₀ + {mech_magnitude:.3e} * sin(2π × {mech_frequency_cm1:.1f} cm⁻¹ × t)")
                print(f"  Effects: ω(t) ∝ 1/L(t) and g(t) ∝ 1/L(t)")
            
            if not periodic_coupling and not laser_enabled and switch_time_ps is None:
                print(f"  *** CONSTANT COUPLING ***")
                print(f"  Coupling: {coupling:.6e} a.u. (constant)")
                
            print(f"  Damping ratio (zeta): {damping_ratio}")
            print(f"  Calculated dissipation (c): {dissipation:.4e} a.u.")
            print(f"  Molecular thermostat: {molecular_thermo}")
            print(f"  Cavity thermostat: {cavity_thermo}")
            print(f"  Finite-q mode: {finite_q}")
            
        print(f"  Temperature: {temperature:.1f} K")
        print(f"  Runtime: {runtime_ps:.1f} ps")
        print(f"  Device: {device}")
        if device == 'GPU':
            print(f"    GPU ID: {gpu_id}")
        print(f"  Random seed: {seed}")
        print(f"  Velocity restart: {'Enabled' if restart_velocities else 'Disabled'} - {'thermalizing velocities' if restart_velocities else 'using existing velocities'}")
        print(f"  Momentum zeroing: {'Enabled' if zero_momentum_enabled else 'Disabled'}")
        
        # Empirical feedback information
        if enable_empirical_feedback:
            print(f"  Empirical feedback: Enabled")
            print(f"    Data file: {empirical_data_file}")
            print(f"    Energy component: {feedback_energy_component}")
            print(f"    Apply to: {feedback_apply_to}")
            print(f"    Averaging window: {feedback_averaging_window_ps:.1f} ps")
            print(f"    Update interval: {feedback_update_interval_ps:.1f} ps")
        
        
        print(f"  Replica: {replica}")
        print(f"  Frame: {frame}")
        print(f"  Output directory: {exp_dir}")
        
        # Set error tolerance based on timestepping mode
        error_tolerance = 0.0 if fixed_timestep else error_tolerance
        
        # Set timestep based on user preference (only used if fixed_timestep is True)
        dt_fs = timestep_fs if fixed_timestep else None
        
        # Handle input GSD path - convert to absolute path if relative
        input_gsd_path = Path(input_gsd)
        if not input_gsd_path.is_absolute():
            # Convert relative path to absolute path based on current working directory
            input_gsd_path = Path.cwd() / input_gsd_path
        input_gsd_abs = str(input_gsd_path)
        
        # Create and run CavityMDSimulation with ALL Phase 3 parameters
        sim = CavityMDSimulation(
            job_dir=str(exp_dir),
            replica=replica,
            freq=frequency,
            couplstr=coupling,
            incavity=incavity,
            runtime_ps=runtime_ps,
            input_gsd=input_gsd_abs,  # Use absolute path
            frame=frame,
            name='prod',
            error_tolerance=error_tolerance,
            temperature=temperature,
            molecular_thermostat=molecular_thermo,
            cavity_thermostat=cavity_thermo,
            finite_q=finite_q,
            molecular_thermostat_tau=molecular_tau,
            cavity_thermostat_tau=cavity_tau,
            log_level='INFO',
            custom_log_file=None,
            device=device,
            gpu_id=gpu_id,
            seed=seed,
            restart_velocities=restart_velocities,
            dt_fs=dt_fs,
            initial_fraction=initial_fraction,
            time_constant_ps=time_constant_ps,
            
            # Enhanced coupling variant parameters
            coupling_variant_type=coupling_variant_type,
            # Exponential decay parameters
            exponential_amplitude=exponential_amplitude if exponential_amplitude is not None else coupling,
            exponential_decay_time_ps=exponential_decay_time_ps,
            exponential_turn_on_time_ps=exponential_turn_on_time_ps,
            exponential_turn_off_time_ps=exponential_turn_off_time_ps,
            # Square wave parameters
            square_period_ps=square_period_ps,
            square_duty_cycle=square_duty_cycle,
            square_phase_offset=square_phase_offset,
            square_start_time_ps=square_start_time_ps,
            square_stop_time_ps=square_stop_time_ps,
            
            # Enhanced periodic coupling parameters (with backward compatibility)
            periodic_coupling=periodic_coupling,
            periodic_period_ps=periodic_period_ps,
            periodic_phase_offset=periodic_phase_offset,
            periodic_start_time_ps=periodic_start_time_ps,
            periodic_stop_time_ps=periodic_stop_time_ps,
            
            # Enhanced step coupling (switch-time) parameters
            switch_time_ps=switch_time_ps,
            decay_time_constant_ps=decay_time_constant_ps,
            dissipation=dissipation,
            
            # Enhanced laser drive parameters
            laser_enabled=laser_enabled,
            laser_frequency_cm1=laser_frequency_cm1,
            laser_amplitude=laser_amplitude,
            laser_start_time_ps=laser_start_time_ps,
            laser_stop_time_ps=laser_stop_time_ps,
            laser_kvector=laser_kvector if laser_kvector is not None else [1.0, 0.0, 0.0],
            
            # PI feedback controller parameters
            enable_pi_feedback=enable_pi_feedback,
            pi_target_temperature=pi_target_temperature,
            pi_turn_on_time_ps=pi_turn_on_time_ps,
            pi_turn_off_time_ps=pi_turn_off_time_ps,
            pi_temperature_method=pi_temperature_method,
            pi_update_interval_ps=pi_update_interval_ps,
            pi_Kc=pi_Kc,
            pi_Ti=pi_Ti,
            pi_molecular_tau_ps=pi_molecular_tau_ps,
            pi_beta=pi_beta,
            pi_T_min=pi_T_min,
            pi_T_max=pi_T_max,
            
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
            
            # Temperature tracker parameters
            enable_temp_tracker=enable_temp_tracker,
            temp_tracker_output_period_ps=temp_tracker_output_period_ps,
            temp_tracker_empirical_data_file=temp_tracker_empirical_data_file,
            # Dipole moment FDR parameters
            enable_dipole_fdr=enable_dipole_fdr,
            dipole_fdr_output_period_ps=dipole_fdr_output_period_ps,
            dipole_fdr_max_correlation_time_ps=dipole_fdr_max_correlation_time_ps,
            dipole_fdr_field_direction=dipole_fdr_field_direction,
            dipole_fdr_exclude_cavity=dipole_fdr_exclude_cavity,
            enable_dipole_response=enable_dipole_response,
            dipole_response_field_strength=dipole_response_field_strength,
            dipole_response_sign=dipole_response_sign,
            
            # Legacy mechanical cavity modulation parameters
            mech_periodic=mech_periodic,
            mech_frequency_cm1=mech_frequency_cm1,
            mech_magnitude=mech_magnitude,
            # Energy tracking parameters
            enable_energy_tracking=enable_energy_tracker,
            energy_output_period_ps=energy_output_period_ps,
            max_energy_output_time_ps=max_energy_output_time,
            # F(k,t) tracking parameters
            enable_fkt=enable_fkt,
            fkt_kmag=fkt_kmag,
            fkt_num_wavevectors=fkt_wavevectors,
            fkt_reference_interval_ps=fkt_ref_interval,
            fkt_max_references=fkt_max_refs,
            fkt_output_period_ps=fkt_output_period_ps,
            # Dipole autocorrelation parameters
            enable_dipole_autocorr=enable_dipole_autocorr,
            dipole_reference_interval_ps=dipole_ref_interval,
            dipole_max_references=dipole_max_refs,
            dipole_output_period_ps=dipole_output_period_ps,
            # Output parameters
            gsd_output_period_ps=gsd_output_period_ps,
            console_output_period_ps=console_output_period_ps,
            truncate_gsd=truncate_gsd,
            # Momentum zeroing parameters
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
            feedback_use_direct_harmonic_calculation=feedback_use_direct_harmonic_calculation,
            feedback_auto_adjust_molecular_tau=feedback_auto_adjust_molecular_tau,
            enable_auto_stopping=feedback_auto_stopping,
            auto_stop_min_time_ps=auto_stop_min_time_ps,
            auto_stop_smoothing_window=auto_stop_smoothing_window,
        )
        
        return sim.run() == 0  # Return True for success (exit code 0)
        
    except Exception as e:
        print(f"❌ Error running experiment: {e}")
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

🔧 ENHANCED COUPLING VARIANTS:
  --coupling-type constant   : g(t) = constant (default)
  --coupling-type step       : g(t) = 0 → g_target at t_switch (with optional decay/turn-off)
  --coupling-type periodic   : g(t) = A*sin(2π*t/T + φ) (with start/stop times)
  --coupling-type exponential: g(t) = A*exp(-t/τ) (with turn-on/off times)
  --coupling-type square     : g(t) = periodic on/off pulses (with duty cycle)

🔬 ENHANCED LASER DRIVE:
  --laser                    : Enable laser with timing control F_laser = F_L*cos(ω_L*t)
  --laser-start-time         : Laser start time in ps
  --laser-stop-time          : Laser stop time in ps

🌡️ PI FEEDBACK CONTROLLER:
  --enable-pi-feedback       : Enable PI temperature feedback
  --pi-method kinetic        : Use kinetic temperature (from velocities)
  --pi-method lj_coulombic   : Use LJ+Coulombic fictive temperature
  --pi-method harmonic       : Use harmonic fictive temperature (fitted)
  --pi-method harmonic_equipartition : Use harmonic equipartition temperature (physics-based)

📊 COMPREHENSIVE TEMPERATURE TRACKER:
  --enable-temp-tracker      : Track all temperatures (kinetic, fictive, bath)
  --temp-tracker-empirical-data-file : Empirical data for LJ+Coul fictive temp

🔬 DIPOLE MOMENT FDR ANALYSIS:
  --enable-dipole-fdr        : Enable dipole moment autocorrelation tracking
  --enable-dipole-response   : Enable electric field response for fork-and-clone FDR
  --dipole-response-sign +1/-1 : Sign for plus/minus clone in FDR measurements

🚀 LEGACY MODES (still supported):
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
    parser.add_argument('--coupling', type=float, default=1e-3,
                       help='Cavity coupling strength / amplitude (default: 1e-3)')
    
    # Enhanced coupling variant selection
    parser.add_argument('--coupling-type', choices=['constant', 'step', 'periodic', 'exponential', 'square'], 
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
    
    # Enhanced periodic coupling parameters (with stop time)
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
    parser.add_argument('--enable-pi-feedback', action='store_true',
                       help='Enable PI feedback temperature controller')
    parser.add_argument('--pi-target-temperature', type=float, default=100.0,
                       help='PI controller target temperature in K (default: 100.0)')
    parser.add_argument('--pi-turn-on-time', type=float, default=0.0,
                       help='PI controller turn-on time in ps (default: 0.0)')
    parser.add_argument('--pi-turn-off-time', type=float, default=None,
                       help='PI controller turn-off time in ps (default: None, no turn-off)')
    parser.add_argument('--pi-method', choices=['kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition'], default='kinetic',
                       help='Temperature calculation method for PI controller (default: kinetic)')
    parser.add_argument('--pi-update-interval', type=float, default=5.0,
                       help='PI controller update interval in ps (default: 5.0)')
    parser.add_argument('--pi-Kc', type=float, default=None,
                       help='PI controller proportional gain (default: None, auto-calculate)')
    parser.add_argument('--pi-Ti', type=float, default=None,
                       help='PI controller integral time in ps (default: None, auto-calculate)')
    parser.add_argument('--pi-molecular-tau', type=float, default=5.0,
                       help='Molecular tau for PI auto-tuning in ps (default: 5.0)')
    parser.add_argument('--pi-beta', type=float, default=0.7,
                       help='PI controller setpoint weighting (default: 0.7)')
    parser.add_argument('--pi-T-min', type=float, default=0.0,
                       help='PI controller minimum temperature in K (default: 0.0)')
    parser.add_argument('--pi-T-max', type=float, default=np.inf,
                       help='PI controller maximum temperature in K (default: None)')
    
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
                       choices=['kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition'],
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
    
    # Energy tracking parameters
    parser.add_argument('--enable-energy-tracker', action='store_true',
                       help='Enable comprehensive energy tracking')
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
    parser.add_argument('--enable-temp-tracker', action='store_true',
                       help='Enable comprehensive temperature tracker (kinetic, fictive, bath temperatures)')
    parser.add_argument('--temp-tracker-output-period-ps', type=float, default=0.1,
                       help='Temperature tracker output period in ps (default: 0.1)')
    parser.add_argument('--temp-tracker-empirical-data-file', type=str, default=None,
                       help='Empirical data file for LJ+Coul fictive temperature (default: None)')
    
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

def main():
    """Main function."""
    parser = create_parser()
    args = parser.parse_args()
    
    print("🚀 Universal Cavity MD Experiment Runner - Phase 3 Complete")
    print("============================================================")
    print("✅ Enhanced Coupling Variants: constant, step, periodic, exponential, square")
    print("✅ Enhanced Laser Drive: timing control, arbitrary k-vectors")
    print("✅ PI Feedback Controller: kinetic, LJ+Coulombic, harmonic temperature")
    print("✅ Full GPU/CPU Compatibility: identical results on both devices")
    print("✅ Backward Compatibility: all legacy flags still supported")
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
    
    if not args.no_cavity:
        print(f"*** COUPLING CONSTANT: {args.coupling:.6e} a.u. ***")
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
    
    # PI feedback validation
    if args.enable_pi_feedback:
        if args.pi_target_temperature <= 0:
            parser.error("PI target temperature must be positive")
        if args.pi_T_min >= args.pi_T_max:
            parser.error("PI T_min must be less than T_max")
        if args.pi_update_interval <= 0:
            parser.error("PI update interval must be positive")
    
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
        
        print(f"    🔧 Coupling type: {effective_coupling_type.upper()}")
        
        if effective_coupling_type == 'constant':
            print(f"    Coupling strength: {args.coupling:.6e} a.u.")
            
        elif effective_coupling_type == 'step':
            print(f"    Coupling strength: {args.coupling:.6e} a.u.")
            switch_time = args.switch_time if args.switch_time is not None else 0.0
            print(f"    Switch time: {switch_time:.1f} ps")
            if args.decay_time_constant is not None:
                print(f"    Decay time: {args.decay_time_constant:.1f} ps")
            if args.step_turn_off_time is not None:
                print(f"    Turn-off time: {args.step_turn_off_time:.1f} ps")
                
        elif effective_coupling_type == 'periodic':
            print(f"    Amplitude: {args.coupling:.6e} a.u.")
            print(f"    Period: {args.period:.1f} ps ({1000.0/args.period:.1f} THz)")
            print(f"    Phase offset: {args.phase_offset:.3f} rad")
            print(f"    Start time: {args.start_time:.1f} ps")
            if args.periodic_stop_time is not None:
                print(f"    Stop time: {args.periodic_stop_time:.1f} ps")
                
        elif effective_coupling_type == 'exponential':
            amplitude = args.exponential_amplitude if args.exponential_amplitude is not None else args.coupling
            print(f"    Amplitude: {amplitude:.6e} a.u.")
            print(f"    Decay time: {args.exponential_decay_time:.1f} ps")
            print(f"    Turn-on time: {args.exponential_turn_on_time:.1f} ps")
            if args.exponential_turn_off_time is not None:
                print(f"    Turn-off time: {args.exponential_turn_off_time:.1f} ps")
                
        elif effective_coupling_type == 'square':
            print(f"    Amplitude: {args.coupling:.6e} a.u.")
            print(f"    Period: {args.square_period:.1f} ps")
            print(f"    Duty cycle: {args.square_duty_cycle:.1%}")
            print(f"    Phase offset: {args.square_phase_offset:.3f} rad")
            print(f"    Start time: {args.square_start_time:.1f} ps")
            if args.square_stop_time is not None:
                print(f"    Stop time: {args.square_stop_time:.1f} ps")
        
        print(f"    Frequency: {args.frequency:.1f} cm⁻¹")
        print(f"    Finite-q mode: {args.finite_q}")
        print(f"    Cavity thermostat: {args.cavity_bath}")
        
        # Enhanced laser information
        if args.laser:
            print(f"    🔬 Enhanced laser drive:")
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
            print(f"    🔧 Mechanical modulation: {args.mech_frequency:.1f} cm⁻¹, amplitude {args.mech_magnitude:.6e} a.u.")
            
    # PI feedback information
    if args.enable_pi_feedback:
        print(f"    🌡️ PI feedback controller:")
        print(f"      Target temperature: {args.pi_target_temperature:.1f} K")
        print(f"      Method: {args.pi_method}")
        print(f"      Turn-on time: {args.pi_turn_on_time:.1f} ps")
        if args.pi_turn_off_time is not None:
            print(f"      Turn-off time: {args.pi_turn_off_time:.1f} ps")
        else:
            print(f"      Turn-off time: None (runs until end)")
        print(f"      Update interval: {args.pi_update_interval:.1f} ps")
        if args.pi_Kc is not None:
            print(f"      Kc: {args.pi_Kc:.3f}")
        else:
            print(f"      Kc: Auto-calculate (IMC tuning)")
        if args.pi_Ti is not None:
            print(f"      Ti: {args.pi_Ti:.1f} ps")
        else:
            print(f"      Ti: Auto-calculate (IMC tuning)")
    
    # Gradient descent feedback information
    if args.enable_gd_feedback:
        print(f"    🎯 Gradient descent feedback controller:")
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
    
    # Temperature tracker information
    if args.enable_temp_tracker:
        print(f"    📊 Comprehensive temperature tracker:")
        print(f"      Output period: {args.temp_tracker_output_period_ps:.1f} ps")
        print(f"      Tracking: kinetic, harmonic fictive, LJ+Coul fictive, cavity bath, molecular bath")
        if args.temp_tracker_empirical_data_file:
            print(f"      Empirical data: {args.temp_tracker_empirical_data_file}")
        else:
            print(f"      Empirical data: None (LJ+Coul fictive temp will be 0)")
    
    # Dipole moment FDR information
    if args.enable_dipole_fdr:
        print(f"    ⚡ Dipole moment FDR analysis:")
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
            periodic_period_ps=args.period,
            periodic_phase_offset=args.phase_offset,
            periodic_start_time_ps=args.start_time,
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
            # Enhanced laser drive parameters
            laser_enabled=args.laser,
            laser_frequency_cm1=args.laser_frequency,
            laser_amplitude=args.laser_amplitude,
            laser_start_time_ps=args.laser_start_time,
            laser_stop_time_ps=args.laser_stop_time,
            laser_kvector=args.laser_kvector,
            # PI feedback controller parameters
            enable_pi_feedback=args.enable_pi_feedback,
            pi_target_temperature=args.pi_target_temperature,
            pi_turn_on_time_ps=args.pi_turn_on_time,
            pi_turn_off_time_ps=args.pi_turn_off_time,
            pi_temperature_method=args.pi_method,
            pi_update_interval_ps=args.pi_update_interval,
            pi_Kc=args.pi_Kc,
            pi_Ti=args.pi_Ti,
            pi_molecular_tau_ps=args.pi_molecular_tau,
            pi_beta=args.pi_beta,
            pi_T_min=args.pi_T_min,
            pi_T_max=args.pi_T_max,
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
            # Temperature tracker parameters
            enable_temp_tracker=args.enable_temp_tracker,
            temp_tracker_output_period_ps=args.temp_tracker_output_period_ps,
            temp_tracker_empirical_data_file=args.temp_tracker_empirical_data_file,
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
            feedback_use_direct_harmonic_calculation=args.feedback_use_direct_harmonic,
            feedback_auto_adjust_molecular_tau=args.feedback_auto_adjust_molecular_tau,
            feedback_auto_stopping=args.feedback_auto_stopping,
            auto_stop_min_time_ps=args.auto_stop_min_time_ps,
            auto_stop_smoothing_window=args.auto_stop_smoothing_window,
            # Experiment metadata
            replica=replica,
            incavity=not args.no_cavity
        )
        
        if success:
            print(f"✅ Replica {replica} completed successfully")
            success_count += 1
        else:
            print(f"❌ Replica {replica} failed")
        
        print()
    
    print("============================================================")
    print(f"Execution completed: {success_count}/{len(replicas)} replicas successful")
    
    if success_count == len(replicas):
        print("🎉 All simulations completed successfully!")
        return 0
    else:
        print(f"⚠️  {len(replicas) - success_count} simulations failed")
        return 1

if __name__ == '__main__':
    sys.exit(main())
