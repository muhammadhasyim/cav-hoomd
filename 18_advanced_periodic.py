#!/usr/bin/env python3
# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""
Advanced Periodic Cavity Molecular Dynamics Experiment Runner

This script provides a comprehensive framework for running cavity MD simulations
with PERIODIC COUPLING CONSTANTS using the CavityMDSimulation class from the 
hoomd.cavitymd plugin. It extends the standard advanced runner with support for
sinusoidal coupling oscillations.

NEW PERIODIC COUPLING FEATURES:
- Sinusoidal coupling oscillation: g(t) = A * sin(2π*t/T + φ)
- Customizable period (default: 1.0 ps = 1 THz frequency)
- Phase offset control (default: 0.0 rad)
- Delayed start times for equilibration
- Full integration with existing analysis features

FEATURES:
- Uses the plugin's CavityMDSimulation class directly (no code duplication)
- Multiple thermostat combinations (Bussi, Langevin)
- Advanced analysis trackers (energy, F(k,t), cavity modes)
- Adaptive or fixed timestep control
- GPU and CPU support
- Comprehensive logging and output control
- SLURM array job support and local multi-replica execution
- Smart cavity particle handling (auto-detects existing cavity particles)

BASIC USAGE:
   # Run a single experiment with periodic coupling (1 ps period)
   python 18_advanced_periodic.py --molecular-bath bussi --cavity-bath langevin --coupling 1e-3 --periodic --runtime 1000
   
   # Run with custom period and phase
   python 18_advanced_periodic.py --molecular-bath bussi --cavity-bath langevin --coupling 1e-3 --periodic --period 0.5 --phase-offset 1.57 --runtime 1000
   
   # Run with delayed periodic start (equilibrate first)
   python 18_advanced_periodic.py --molecular-bath bussi --cavity-bath langevin --coupling 1e-3 --periodic --start-time 50 --runtime 1000
   
   # Run without cavity (molecular-only simulation)
   python 18_advanced_periodic.py --molecular-bath bussi --no-cavity --runtime 1000
   
   # Run multiple replicas locally (replica 0 supported)
   python 18_advanced_periodic.py --molecular-bath bussi --cavity-bath langevin --coupling 1e-3 --periodic --replicas 0-4 --runtime 500
   
   # Run with finite-q cavity mode
   python 18_advanced_periodic.py --molecular-bath bussi --cavity-bath langevin --finite-q --coupling 1e-3 --periodic --runtime 1000

   # Standard constant coupling (fallback to regular behavior)
   python 18_advanced_periodic.py --molecular-bath bussi --cavity-bath langevin --coupling 1e-3 --runtime 1000
   
   # Run with fictive temperature feedback control
   python 18_advanced_periodic.py --molecular-bath bussi --cavity-bath langevin --coupling 1e-3 --runtime 1000 --enable-empirical-feedback --empirical-data-file equilibrium_data.txt

PERIODIC COUPLING PARAMETERS:
- --periodic: Enable periodic coupling oscillation
- --period: Oscillation period in ps (default: 1.0 ps = 1 THz)
- --phase-offset: Phase offset in radians (default: 0.0)
- --start-time: Time to start oscillation in ps (default: 0.0)

CAVITY PARTICLE HANDLING:
   The script automatically detects whether cavity particles already exist in your GSD file:
   - If cavity particles exist: Uses them (validates count and properties)
   - If no cavity particles exist: Adds new ones automatically (when cavity coupling enabled)
   - Clear error messages if configuration is invalid

THERMOSTAT OPTIONS:
- --molecular-bath: bussi, langevin, none (thermostat for molecular particles)
- --cavity-bath: bussi, langevin, none (thermostat for cavity particle)
- --finite-q: Enable finite-q cavity mode (default: q=0 mode)

REPLICA EXECUTION:
- Local execution: --replicas 0-4 or --replicas 0,1,2,3,4 (replica 0 supported)
- SLURM array jobs: Automatically detects SLURM_ARRAY_TASK_ID

ADVANCED FEATURES:
- --enable-energy-tracker: Detailed energy component tracking
- --enable-fkt: F(k,t) density correlation functions
- --fixed-timestep: Use fixed timestep instead of adaptive
- --device GPU: Run on GPU instead of CPU
- --seed: Control random seed for reproducibility
- --zero-momentum: Periodic momentum zeroing to prevent center-of-mass drift
- --switch-time: Time-varying coupling activation (alternative to periodic)
- --decay-time-constant: Exponential decay of coupling after switching
- --enable-empirical-feedback: Enable fictive temperature tracking and feedback
- --empirical-data-file: File with empirical energy-temperature calibration data

FICTIVE TEMPERATURE TRACKING:
- --enable-empirical-feedback: Enable real-time fictive temperature calculation and feedback
- --empirical-data-file: Path to energy-temperature calibration data (required for feedback)
- --feedback-energy-component: Energy component for feedback ('lj_coulombic', 'total_PE', 'harmonic')
- --feedback-apply-to: Apply feedback to ('molecular', 'cavity', 'both', 'none') thermostats
- --feedback-averaging-window-ps: Time window for temperature averaging
- --feedback-update-interval-ps: Interval between thermostat updates
- --feedback-turn-off-time-ps: Time to disable feedback (None = never)
- --feedback-auto-stopping: Enable automatic coupling stop when fictive and harmonic temperatures converge
  NOTE: This stops the PERIODIC COUPLING (sets to zero), NOT the feedback control.
  Allows reaching equilibrium supercooled states after temperature convergence.
- --auto-stop-min-time-ps: Minimum time before auto-stop triggers
- --auto-stop-smoothing-window: Window size for temperature smoothing

OUTPUT CONTROL:
- Separate output periods for different observables (energy, F(k,t), trajectories, console)
- Organized directory structure with periodic coupling parameters in names
- Comprehensive logging with timestamps
- Fictive temperature CSV output when empirical feedback is enabled

See --help for all available options.
"""

import sys
import argparse
import time
from pathlib import Path
import numpy as np

# Import the CavityForce and utilities from the plugin
from hoomd.cavitymd.utils import PhysicalConstants, get_slurm_info, parse_replicas
from hoomd.cavitymd.simulation import CavityMDSimulation

# =============================================================================
# SIMULATION FUNCTIONS (No need for local class - using plugin directly)
# =============================================================================

def run_single_experiment(molecular_thermo, cavity_thermo, finite_q, 
                         coupling, temperature, frequency, replica, frame, 
                         runtime_ps, molecular_tau, cavity_tau, enable_fkt, fkt_kmag, fkt_wavevectors, 
                         fkt_ref_interval, fkt_max_refs, max_energy_output_time=None, 
                         device='CPU', gpu_id=0, incavity=True, fixed_timestep=False, 
                         timestep_fs=1.0, enable_energy_tracking=False, 
                         energy_output_period_ps=0.1, fkt_output_period_ps=1.0, 
                         gsd_output_period_ps=50.0, console_output_period_ps=1.0, 
                         truncate_gsd=False, seed=None, restart_velocities=True,
                         switch_time_ps=None, decay_time_constant_ps=None, damping_ratio=0.0,
                         enable_dipole_autocorr=False, dipole_ref_interval=1.0, dipole_max_refs=10, 
                         dipole_output_period_ps=1.0, error_tolerance=5.0, initial_fraction=1e-5, 
                         time_constant_ps=50.0, zero_momentum_enabled=False, zero_momentum_period_ps=1.0,
                         # NEW: Periodic coupling parameters
                         periodic_coupling=False, periodic_period_ps=1.0, periodic_phase_offset=0.0, 
                         periodic_start_time_ps=0.0,
                         # NEW: Empirical feedback parameters for fictive temperature tracking
                         enable_empirical_feedback=False, empirical_data_file=None,
                         feedback_output_period_ps=0.1, feedback_energy_component='lj_coulombic',
                         feedback_apply_to='both', feedback_averaging_window_ps=5.0,
                         feedback_update_interval_ps=5.0, feedback_T_min=50.0, feedback_T_max=300.0,
                         feedback_turn_off_time_ps=None, feedback_use_direct_harmonic_calculation=True,
                         feedback_auto_adjust_molecular_tau=False,
                         feedback_auto_stopping=False, auto_stop_min_time_ps=10.0, auto_stop_smoothing_window=5,
                         # NEW: Simple temperature drop protocol (alternative to empirical feedback)
                         enable_simple_temp_drop=False, temp_drop_target=50.0, temp_drop_time=10.0, temp_drop_apply_to='both',
                         input_gsd='init-0.gsd'):
    """
    Run a single experiment using the CavityMDSimulation class from the plugin.
    Extended with periodic coupling and fictive temperature tracking support.
    
    NEW FEATURES:
    - Periodic coupling oscillation: g(t) = A * sin(2π*t/T + φ)
    - Real-time fictive temperature calculation and feedback
    - Empirical energy-temperature calibration
    - Thermostat temperature adjustment based on systemic temperature
    """
    
    try:
        # Calculate dissipation from damping_ratio
        phmass = 1.0  # Photon mass is 1.0 in a.u.
        omegac = frequency / PhysicalConstants.HARTREE_TO_CM_MINUS1
        dissipation = 2 * damping_ratio * phmass * omegac

        # Create experiment directory with appropriate naming
        if incavity:
            # For cavity simulations, include coupling strength in directory name
            coupling_str = f"{coupling:.0e}".replace("-", "neg").replace("+", "pos")
            
            if periodic_coupling:
                # Include periodic parameters in directory name
                period_str = f"_period_{periodic_period_ps}ps"
                if periodic_phase_offset != 0.0:
                    phase_str = f"_phase_{periodic_phase_offset:.2f}rad"
                else:
                    phase_str = ""
                if periodic_start_time_ps != 0.0:
                    start_str = f"_start_{periodic_start_time_ps}ps"
                else:
                    start_str = ""
                exp_dir = Path(f"periodic_coupling_{coupling_str}{period_str}{phase_str}{start_str}")
                
            elif switch_time_ps is not None:
                # Include switch time in directory name for time-varying simulations
                switch_str = f"_switch_{switch_time_ps}ps"
                if decay_time_constant_ps is not None:
                    # Include decay time constant for exponential decay simulations
                    decay_str = f"_decay_{decay_time_constant_ps}ps"
                    exp_dir = Path(f"cavity_coupling_{coupling_str}{switch_str}{decay_str}")
                else:
                    exp_dir = Path(f"cavity_coupling_{coupling_str}{switch_str}")
            else:
                exp_dir = Path(f"cavity_coupling_{coupling_str}")
        else:
            # For non-cavity simulations
            exp_dir = Path("no_cavity")
        exp_dir.mkdir(exist_ok=True)
        
        print(f"Running experiment:")
        print(f"  Cavity coupling: {'Enabled' if incavity else 'Disabled'}")
        if incavity:
            if periodic_coupling:
                print(f"  *** PERIODIC COUPLING ENABLED ***")
                print(f"  Amplitude: {coupling:.6e} a.u.")
                print(f"  Period: {periodic_period_ps} ps ({1000/periodic_period_ps:.1f} THz)")
                print(f"  Phase offset: {periodic_phase_offset:.3f} rad")
                print(f"  Start time: {periodic_start_time_ps} ps")
                print(f"  Formula: g(t) = {coupling:.3e} * sin(2π*t/{periodic_period_ps} + {periodic_phase_offset:.3f})")
            else:
                print(f"  *** COUPLING STRENGTH: {coupling:.6e} a.u. ***")  # Enhanced prominence
                if switch_time_ps is not None:
                    print(f"  Switch time: {switch_time_ps} ps")
                    if decay_time_constant_ps is not None:
                        print(f"  Decay time constant: {decay_time_constant_ps} ps (exponential decay)")
                    else:
                        print(f"  Decay: None (step function)")
            
            print(f"  Damping ratio (zeta): {damping_ratio}")
            print(f"  Calculated dissipation (c): {dissipation:.4e} a.u.")
            print(f"  Molecular thermostat: {molecular_thermo}")
            print(f"  Cavity thermostat: {cavity_thermo}")
            print(f"  Finite-q mode: {finite_q}")
        else:
            print(f"  Molecular thermostat: {molecular_thermo}")
        
        # Fictive temperature feedback information
        if enable_empirical_feedback:
            print(f"  *** FICTIVE TEMPERATURE FEEDBACK ENABLED ***")
            print(f"  Empirical data file: {empirical_data_file}")
            print(f"  Energy component: {feedback_energy_component}")
            print(f"  Apply to: {feedback_apply_to}")
            print(f"  Averaging window: {feedback_averaging_window_ps} ps")
            print(f"  Update interval: {feedback_update_interval_ps} ps")
            if feedback_turn_off_time_ps is not None:
                print(f"  Turn off time: {feedback_turn_off_time_ps} ps")
            if feedback_auto_stopping:
                print(f"  🤖 Auto-stopping: Enabled (min time: {auto_stop_min_time_ps} ps, smoothing: {auto_stop_smoothing_window})")
            print(f"  Temperature range: [{feedback_T_min:.1f}, {feedback_T_max:.1f}] K")
        
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
        
        # Create and run CavityMDSimulation
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
            enable_fkt=enable_fkt,
            fkt_kmag=fkt_kmag,
            fkt_num_wavevectors=fkt_wavevectors,
            fkt_reference_interval_ps=fkt_ref_interval,
            fkt_max_references=fkt_max_refs,
            max_energy_output_time_ps=max_energy_output_time,
            enable_energy_tracking=enable_energy_tracking,
            dt_fs=dt_fs,
            device=device,
            gpu_id=gpu_id,
            energy_output_period_ps=energy_output_period_ps,
            fkt_output_period_ps=fkt_output_period_ps,
            gsd_output_period_ps=gsd_output_period_ps,
            console_output_period_ps=console_output_period_ps,
            enable_text_output=False,
            text_output_file=None,
            truncate_gsd=truncate_gsd,
            seed=seed,
            restart_velocities=restart_velocities,
            switch_time_ps=switch_time_ps,
            decay_time_constant_ps=decay_time_constant_ps,
            dissipation=dissipation,
            enable_dipole_autocorr=enable_dipole_autocorr,
            dipole_reference_interval_ps=dipole_ref_interval,
            dipole_max_references=dipole_max_refs,
            dipole_output_period_ps=dipole_output_period_ps,
            initial_fraction=initial_fraction,
            time_constant_ps=time_constant_ps,
            zero_momentum_enabled=zero_momentum_enabled,
            zero_momentum_period_ps=zero_momentum_period_ps,
            # NEW: Periodic coupling parameters
            periodic_coupling=periodic_coupling,
            periodic_period_ps=periodic_period_ps,
            periodic_phase_offset=periodic_phase_offset,
            periodic_start_time_ps=periodic_start_time_ps,
            # NEW: Empirical feedback parameters for fictive temperature tracking
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
            # NEW: Simple temperature drop parameters
            enable_simple_temp_drop=enable_simple_temp_drop,
            temp_drop_target=temp_drop_target,
            temp_drop_time=temp_drop_time,
            temp_drop_apply_to=temp_drop_apply_to
        )
        
        # Run the simulation
        return sim.run() == 0  # Return True for success (exit code 0)
        
    except Exception as e:
        print(f"ERROR: Experiment failed: {e}")
        return False

def main():
    """Simplified main function for cavity MD experiments with periodic coupling."""
    parser = argparse.ArgumentParser(
        description='Advanced Periodic Cavity MD Experiment Runner',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    # Basic simulation parameters
    parser.add_argument('--molecular-bath', type=str, default='bussi', choices=['bussi', 'langevin', 'none'], 
                       help='Molecular thermostat type (default: bussi)')
    parser.add_argument('--cavity-bath', type=str, default='langevin', choices=['bussi', 'langevin', 'none'], 
                       help='Cavity thermostat type (default: langevin)')
    parser.add_argument('--finite-q', action='store_true', 
                       help='Use finite-q cavity mode (default: q=0 mode)')
    parser.add_argument('--coupling', type=float, default=1e-3, 
                       help='Cavity coupling strength / amplitude (default: 1e-3)')
    
    # NEW: Periodic coupling parameters
    parser.add_argument('--periodic', action='store_true',
                       help='Enable periodic coupling oscillation: g(t) = A * sin(2π*t/T + φ)')
    parser.add_argument('--period', type=float, default=1.0,
                       help='Oscillation period in ps (default: 1.0 ps = 1 THz)')
    parser.add_argument('--phase-offset', type=float, default=0.0,
                       help='Phase offset in radians (default: 0.0)')
    parser.add_argument('--start-time', type=float, default=0.0,
                       help='Time to start periodic oscillation in ps (default: 0.0)')
    
    # Existing time-varying coupling parameters (alternative to periodic)
    parser.add_argument('--switch-time', type=float, 
                       help='Time in ps when coupling and dissipation turn on (default: on from start)')
    parser.add_argument('--decay-time-constant', type=float,
                       help='Exponential decay time constant in ps for coupling after switch (default: no decay)')
    
    parser.add_argument('--damping-ratio', type=float, default=0.0,
                       help='Damping ratio (zeta) for the cavity mode (default: 0.0)')
    parser.add_argument('--temperature', type=float, default=100.0, 
                       help='Temperature in K (default: 100.0)')
    parser.add_argument('--frequency', type=float, default=2000.0, 
                       help='Cavity frequency in cm⁻¹ (default: 2000.0)')
    parser.add_argument('--runtime', type=float, default=500.0, 
                       help='Runtime in ps (default: 500.0)')
    parser.add_argument('--no-cavity', action='store_true', 
                       help='Disable cavity coupling (molecular-only simulation)')

    # Input file specification
    parser.add_argument('--input-gsd', type=str, default='init-0.gsd',
                       help='Path to input GSD file (default: init-0.gsd)')
    parser.add_argument('--frame', type=int, default=-1,
                       help='Frame number to read from GSD file (default: -1, last frame)')

    # Replica control
    parser.add_argument('--replicas', type=str, 
                       help='Replica specification (e.g., "0,1,2" or "0-5"). Replica 0 is supported.')
    
    # Thermostat parameters
    parser.add_argument('--molecular-tau', type=float, default=5.0, 
                       help='Molecular thermostat tau in ps (default: 5.0)')
    parser.add_argument('--cavity-tau', type=float, default=5.0, 
                       help='Cavity thermostat tau in ps (default: 5.0)')
    
    # Timestep control
    parser.add_argument('--fixed-timestep', action='store_true', 
                       help='Use fixed timestep instead of adaptive')
    parser.add_argument('--timestep', type=float, default=1.0, 
                       help='Fixed timestep in fs (default: 1.0)')
    
    # Adaptive timestep control
    parser.add_argument('--error-tolerance', type=float, default=5.0,
                       help='Target error tolerance for adaptive timestep (default: 5.0)')
    parser.add_argument('--initial-fraction', type=float, default=1e-5,
                       help='Initial fraction for shock dampening (ratio of initial to target error tolerance, default: 1e-5)')
    parser.add_argument('--time-constant-ps', type=float, default=50.0,
                       help='Time constant for error tolerance ramping in ps (default: 50.0)')
    
    # Energy tracking
    parser.add_argument('--enable-energy-tracker', action='store_true', 
                       help='Enable comprehensive energy tracking')
    
    # Output control options - separate periods for different observables
    parser.add_argument('--energy-output-period-ps', type=float, default=0.1, 
                       help='Energy tracker output period in ps (default: 0.1)')
    parser.add_argument('--fkt-output-period-ps', type=float, default=1.0, 
                       help='F(k,t) tracker output period in ps (default: 1.0)')
    parser.add_argument('--gsd-output-period-ps', type=float, default=50.0, 
                       help='GSD trajectory output period in ps (default: 50.0)')
    parser.add_argument('--console-output-period-ps', type=float, default=1.0, 
                       help='Console output period in ps (default: 1.0)')
    
    # F(k,t) options
    parser.add_argument('--enable-fkt', action='store_true', 
                       help='Enable F(k,t) density correlation calculation')
    parser.add_argument('--fkt-kmag', type=float, default=1.0, 
                       help='F(k,t) k magnitude (default: 1.0)')
    parser.add_argument('--fkt-wavevectors', type=int, default=50, 
                       help='F(k,t) number of wavevectors (default: 50)')
    parser.add_argument('--fkt-ref-interval', type=float, default=1.0, 
                       help='F(k,t) reference interval in ps (default: 1.0)')
    parser.add_argument('--fkt-max-refs', type=int, default=10, 
                       help='F(k,t) maximum references (default: 10)')
    
    # Dipole autocorrelation options
    parser.add_argument('--enable-dipole-autocorr', action='store_true', 
                       help='Enable dipole autocorrelation calculation')
    parser.add_argument('--dipole-ref-interval', type=float, default=1.0, 
                       help='Dipole autocorrelation reference interval in ps (default: 1.0)')
    parser.add_argument('--dipole-max-refs', type=int, default=10, 
                       help='Dipole autocorrelation maximum references (default: 10)')
    parser.add_argument('--dipole-output-period-ps', type=float, default=1.0, 
                       help='Dipole autocorrelation output period in ps (default: 1.0)')
    
    parser.add_argument('--max-energy-output-time', type=float, 
                       help='Maximum energy output time in ps (default: no limit)')
    
    # Device options
    parser.add_argument('--device', type=str, default='CPU', choices=['CPU', 'GPU'], 
                       help='Compute device (default: CPU)')
    parser.add_argument('--gpu-id', type=int, default=0, 
                       help='GPU ID when using GPU device (default: 0)')
    
    # GSD output control
    parser.add_argument('--truncate-gsd', action='store_true', 
                       help='Truncate GSD output file if it exists (default: append)')
    
    # Seed control
    parser.add_argument('--seed', type=int, 
                       help='Random seed for simulation (default: replica-based deterministic seed)')
    
    # Velocity control
    parser.add_argument('--no-restart-velocities', action='store_true', 
                       help='Do not restart velocities - use existing velocities from GSD file (default: restart velocities)')
    
    # Momentum zeroing control
    parser.add_argument('--zero-momentum', action='store_true', 
                       help='Enable periodic momentum zeroing to prevent center-of-mass drift (default: disabled)')
    parser.add_argument('--zero-momentum-period-ps', type=float, default=1.0, 
                       help='Period for momentum zeroing in ps (default: 1.0)')
    
    # Empirical feedback and fictive temperature tracking
    parser.add_argument('--enable-empirical-feedback', action='store_true',
                       help='Enable empirical temperature feedback control based on fictive temperature')
    parser.add_argument('--empirical-data-file', type=str,
                       help='Path to empirical energy-temperature calibration data file (required for feedback)')
    parser.add_argument('--feedback-output-period-ps', type=float, default=0.1,
                       help='Output period for empirical feedback CSV data in ps (default: 0.1)')
    parser.add_argument('--feedback-energy-component', type=str, default='lj_coulombic',
                       choices=['total_PE', 'lj_coulombic', 'harmonic'],
                       help='Energy component for feedback calculation (default: lj_coulombic)')
    parser.add_argument('--feedback-apply-to', type=str, default='both',
                       choices=['molecular', 'cavity', 'both', 'none'],
                       help='Which thermostats to apply feedback to: molecular, cavity, both, or none (default: both)')
    parser.add_argument('--feedback-averaging-window-ps', type=float, default=5.0,
                       help='Time window for averaging fictive temperature in ps (default: 5.0)')
    parser.add_argument('--feedback-update-interval-ps', type=float, default=5.0,
                       help='Interval between thermostat temperature updates in ps (default: 5.0)')
    parser.add_argument('--feedback-T-min', type=float, default=50.0,
                       help='Minimum allowed temperature in K (default: 50.0)')
    parser.add_argument('--feedback-T-max', type=float, default=300.0,
                       help='Maximum allowed temperature in K (default: 300.0)')
    parser.add_argument('--feedback-turn-off-time-ps', type=float,
                       help='Time in ps to turn off empirical feedback (default: None, never turn off)')
    parser.add_argument('--feedback-use-direct-harmonic', action='store_true',
                       help='Use direct calculation T=4*E/(N*kB) for harmonic component (default: False)')
    parser.add_argument('--feedback-auto-adjust-molecular-tau', action='store_true',
                       help='Automatically adjust molecular thermostat tau when feedback turns off (default: False)')
    
    # Auto-stopping functionality for empirical feedback
    parser.add_argument('--feedback-auto-stopping', action='store_true',
                       help='Enable automatic coupling stopping when fictive and harmonic temperatures converge (default: False)')
    parser.add_argument('--auto-stop-min-time-ps', type=float, default=10.0,
                       help='Minimum time before auto-stop can trigger in ps (default: 10.0)')
    parser.add_argument('--auto-stop-smoothing-window', type=int, default=5,
                       help='Window size for temperature smoothing in auto-stop detection (default: 5)')
    
    # Simple temperature drop protocol (alternative to empirical feedback)
    parser.add_argument('--enable-simple-temp-drop', action='store_true',
                       help='Enable simple temperature drop protocol instead of empirical feedback (default: False)')
    parser.add_argument('--temp-drop-target', type=float, default=50.0,
                       help='Target temperature for temperature drop in K (default: 50.0)')
    parser.add_argument('--temp-drop-time', type=float, default=10.0,
                       help='Time in ps when temperature drop occurs (default: 10.0)')
    parser.add_argument('--temp-drop-apply-to', type=str, default='both',
                       choices=['molecular', 'cavity', 'both'],
                       help='Which thermostats to apply temperature drop to (default: both)')
    
    args = parser.parse_args()
    
    print("Advanced Periodic Cavity MD Experiment Runner")
    print("="*60)
    
    # Validate periodic coupling parameters
    if args.periodic and args.switch_time and not args.enable_empirical_feedback and not args.enable_simple_temp_drop:
        parser.error("Cannot use both --periodic and --switch-time unless empirical feedback or simple temperature drop is enabled. "
                    "For these protocols, --switch-time controls activation timing.")
    
    if args.periodic and args.period <= 0:
        parser.error("Period must be positive")
    
    # Validate temperature control protocols (mutually exclusive)
    if args.enable_empirical_feedback and args.enable_simple_temp_drop:
        parser.error("Cannot use both --enable-empirical-feedback and --enable-simple-temp-drop. Choose one.")
    
    # Validate empirical feedback parameters
    if args.enable_empirical_feedback:
        if not args.empirical_data_file:
            parser.error("--empirical-data-file is required when --enable-empirical-feedback is used")
        
        # Check if the empirical data file exists
        empirical_path = Path(args.empirical_data_file)
        if not empirical_path.exists():
            parser.error(f"Empirical data file not found: {args.empirical_data_file}")
        
        # Automatically enable energy tracking when empirical feedback is enabled
        if not args.enable_energy_tracker:
            print("INFO: Automatically enabling energy tracking for empirical feedback")
            args.enable_energy_tracker = True
    
    # Validate simple temperature drop parameters
    if args.enable_simple_temp_drop:
        if not args.empirical_data_file:
            parser.error("--empirical-data-file is required when --enable-simple-temp-drop is used (needed for fictive temperature calculation)")
        if args.temp_drop_target <= 0:
            parser.error("--temp-drop-target must be positive")
        if args.temp_drop_time < 0:
            parser.error("--temp-drop-time must be non-negative")
        
        # Check if the empirical data file exists
        empirical_path = Path(args.empirical_data_file)
        if not empirical_path.exists():
            parser.error(f"Empirical data file not found: {args.empirical_data_file}")
        
        # Automatically enable energy tracking for auto-stopping functionality
        if not args.enable_energy_tracker:
            print("INFO: Automatically enabling energy tracking for temperature drop auto-stopping")
            args.enable_energy_tracker = True
    
    # Determine replica list
    task_id, job_id = get_slurm_info()
    
    if task_id is not None:
        # Running under SLURM array job
        replica_list = [task_id]
        print(f"SLURM array job detected: Task {task_id} (Job {job_id})")
    else:
        # Local execution - parse replicas
        replica_list = parse_replicas(args.replicas)
        print(f"Local execution: Replicas {replica_list}")
    
    # Set up simulation parameters
    incavity = not args.no_cavity
    molecular_thermo = args.molecular_bath
    cavity_thermo = args.cavity_bath if incavity else 'none'
    finite_q = args.finite_q
    
    # Display coupling information prominently
    if args.periodic:
        print(f"\n*** PERIODIC COUPLING ENABLED ***")
        print(f"Amplitude: {args.coupling:.6e} a.u.")
        print(f"Period: {args.period} ps ({1000/args.period:.1f} THz)")
        print(f"Phase offset: {args.phase_offset:.3f} rad")
        print(f"Start time: {args.start_time} ps")
        print(f"Formula: g(t) = {args.coupling:.3e} * sin(2π*t/{args.period} + {args.phase_offset:.3f})")
    else:
        print(f"\n*** COUPLING CONSTANT: {args.coupling:.6e} a.u. ***")
    
    print(f"\nSimulation Configuration:")
    print(f"  Cavity coupling: {'Enabled' if incavity else 'Disabled'}")
    if incavity:
        if args.periodic:
            print(f"    Coupling type: PERIODIC")
            print(f"    Amplitude: {args.coupling:.6e} a.u.")
            print(f"    Period: {args.period} ps ({1000/args.period:.1f} THz)")
            print(f"    Phase offset: {args.phase_offset:.3f} rad")
            print(f"    Start time: {args.start_time} ps")
        else:
            print(f"    Coupling type: CONSTANT")
            print(f"    Coupling strength: {args.coupling:.6e} a.u.")  # Enhanced format
        print(f"    Frequency: {args.frequency} cm⁻¹")
        print(f"    Finite-q mode: {finite_q}")
        print(f"    Cavity thermostat: {cavity_thermo}")
    print(f"  Molecular thermostat: {molecular_thermo}")
    print(f"  Temperature: {args.temperature} K")
    print(f"  Runtime: {args.runtime} ps")
    print(f"  Device: {args.device}")
    if args.device == 'GPU':
        print(f"    GPU ID: {args.gpu_id}")
    print(f"  Random seed: {args.seed if args.seed is not None else 'replica-based (deterministic)'}")
    print(f"  Velocity restart: {'Disabled - using GSD velocities' if args.no_restart_velocities else 'Enabled - thermalizing velocities'}")
    print(f"  Momentum zeroing: {'Enabled' if args.zero_momentum else 'Disabled'}")
    if args.zero_momentum:
        print(f"    Period: {args.zero_momentum_period_ps} ps")
    print(f"  Empirical feedback: {'Enabled' if args.enable_empirical_feedback else 'Disabled'}")
    if args.enable_empirical_feedback:
        print(f"    Data file: {args.empirical_data_file}")
        print(f"    Energy component: {args.feedback_energy_component}")
        print(f"    Apply to: {args.feedback_apply_to}")
        print(f"    Averaging window: {args.feedback_averaging_window_ps} ps")
        print(f"    Update interval: {args.feedback_update_interval_ps} ps")
        if args.feedback_turn_off_time_ps:
            print(f"    Turn off time: {args.feedback_turn_off_time_ps} ps")
    
    # Set up device configuration
    device = args.device.upper()
    
    # Performance tracking
    start_time = time.time()
    successful_experiments = 0
    failed_experiments = 0
    
    print(f"\nStarting execution for {len(replica_list)} replica(s)...")
    print("="*60)
    
    # Run replicas
    for replica in replica_list:
        frame = args.frame  # Use frame from command line argument
        
        print(f"\nRunning replica {replica}...")
        
        # Run experiment
        success = run_single_experiment(
            molecular_thermo=molecular_thermo,
            cavity_thermo=cavity_thermo,
            finite_q=finite_q,
            coupling=args.coupling,
            temperature=args.temperature,
            frequency=args.frequency,
            replica=replica,
            frame=frame,
            runtime_ps=args.runtime,
            molecular_tau=args.molecular_tau,
            cavity_tau=args.cavity_tau,
            enable_fkt=args.enable_fkt,
            fkt_kmag=args.fkt_kmag,
            fkt_wavevectors=args.fkt_wavevectors,
            fkt_ref_interval=args.fkt_ref_interval,
            fkt_max_refs=args.fkt_max_refs,
            max_energy_output_time=args.max_energy_output_time,
            device=device,
            gpu_id=args.gpu_id,
            incavity=incavity,
            fixed_timestep=args.fixed_timestep,
            timestep_fs=args.timestep,
            enable_energy_tracking=args.enable_energy_tracker,
            energy_output_period_ps=args.energy_output_period_ps,
            fkt_output_period_ps=args.fkt_output_period_ps,
            gsd_output_period_ps=args.gsd_output_period_ps,
            console_output_period_ps=args.console_output_period_ps,
            truncate_gsd=args.truncate_gsd,
            seed=args.seed,
            restart_velocities=not args.no_restart_velocities,
            switch_time_ps=args.switch_time,
            damping_ratio=args.damping_ratio,
            enable_dipole_autocorr=args.enable_dipole_autocorr,
            dipole_ref_interval=args.dipole_ref_interval,
            dipole_max_refs=args.dipole_max_refs,
            dipole_output_period_ps=args.dipole_output_period_ps,
            error_tolerance=args.error_tolerance,
            initial_fraction=args.initial_fraction,
            time_constant_ps=args.time_constant_ps,
            zero_momentum_enabled=args.zero_momentum,
            zero_momentum_period_ps=args.zero_momentum_period_ps,
            decay_time_constant_ps=args.decay_time_constant,
            # NEW: Periodic coupling parameters
            periodic_coupling=args.periodic,
            periodic_period_ps=args.period,
            periodic_phase_offset=args.phase_offset,
            periodic_start_time_ps=args.start_time,
            # NEW: Empirical feedback parameters for fictive temperature tracking
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
            # NEW: Simple temperature drop parameters
            enable_simple_temp_drop=args.enable_simple_temp_drop,
            temp_drop_target=args.temp_drop_target,
            temp_drop_time=args.temp_drop_time,
            temp_drop_apply_to=args.temp_drop_apply_to,
            input_gsd=args.input_gsd
        )
        
        if success:
            successful_experiments += 1
            if incavity:
                if args.periodic:
                    print(f"SUCCESS: Replica {replica} completed successfully (periodic coupling: {args.coupling:.6e} a.u., {args.period} ps period)")
                else:
                    print(f"SUCCESS: Replica {replica} completed successfully (coupling: {args.coupling:.6e} a.u.)")
            else:
                print(f"SUCCESS: Replica {replica} completed successfully (no cavity)")
        else:
            failed_experiments += 1
            if incavity:
                if args.periodic:
                    print(f"ERROR: Replica {replica} failed (periodic coupling: {args.coupling:.6e} a.u., {args.period} ps period)")
                else:
                    print(f"ERROR: Replica {replica} failed (coupling: {args.coupling:.6e} a.u.)")
            else:
                print(f"ERROR: Replica {replica} failed (no cavity)")
    
    # Final summary
    end_time = time.time()
    total_wall_time = end_time - start_time
    total_experiments = len(replica_list)
    
    print("\n" + "="*60)
    print("Execution Summary")
    print("="*60)
    if incavity:
        if args.periodic:
            print(f"Periodic coupling used: amplitude {args.coupling:.6e} a.u., period {args.period} ps")
        else:
            print(f"Constant coupling used: {args.coupling:.6e} a.u.")
    else:
        print("No cavity coupling (molecular-only simulation)")
    print(f"Total replicas: {total_experiments}")
    print(f"Successful: {successful_experiments}")
    print(f"Failed: {failed_experiments}")
    print(f"Wall time: {total_wall_time:.2f} seconds")
    
    if failed_experiments > 0:
        print(f"\nWARNING: {failed_experiments} replicas failed - check individual logs for details")
        return 1
    else:
        print("\nAll replicas completed successfully!")
        return 0

if __name__ == '__main__':
    sys.exit(main())
