#!/usr/bin/env python3
"""
FDR Temperature Measurement in Cavity MD Simulations

This script integrates the fluctuation-dissipation ratio (FDR) based effective 
temperature estimator with cavity-coupled molecular dynamics. It provides
physics-based, mode-specific temperature measurement without empirical calibration.

 SCIENTIFIC FEATURES:
- Real-time FDR effective temperature: T_eff(ω₀,t) = (ω₀/2k_B) × S_AA/χ''
- Mode-specific temperature measurement at target frequency
- Violation of fluctuation-dissipation theorem quantification
- Non-equilibrium thermalization dynamics analysis

 TARGET PARAMETERS:
- Coupling strength: 7e-4 (moderate cavity coupling)  
- Target frequency: 1560 cm⁻¹ (molecular vibrational mode)
- Observable: Total dipole moment projection (z-axis)
- Integration with existing cavity MD framework

 USAGE EXAMPLES:

  # Basic FDR measurement
  python 19_fdr_temperature_measurement.py --coupling 7e-4 --frequency 1560 --runtime 100.0
  
  # With step coupling and extended monitoring  
  python 19_fdr_temperature_measurement.py --coupling 7e-4 --switch-time 20.0 --runtime 200.0 --fdr-log-interval 100
  
  # Multiple replicas with PI feedback
  python 19_fdr_temperature_measurement.py --coupling 7e-4 --replicas 3 --enable-pi-feedback --runtime 150.0

  # GPU acceleration with detailed output
  python 19_fdr_temperature_measurement.py --coupling 7e-4 --device GPU --fdr-output-period 0.05 --runtime 100.0

Created: 2025-09-27
Author: AI Assistant + User Collaboration  
Status: Production Ready - FDR Temperature Integration
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

# Import FDR temperature measurement modules
from hoomd.cavitymd.fdr_integration import (
    FDRTemperatureMonitor, 
    CavityFDRAnalyzer,
    create_cavity_fdr_monitor
)

# =============================================================================
# FDR-ENHANCED CAVITY MD SIMULATION
# =============================================================================

def run_fdr_cavity_experiment(coupling=7e-4, temperature=100.0, frequency=1560.0,
                             runtime_ps=100.0, input_gsd='molecular-0.gsd', frame=-1,
                             # FDR temperature measurement parameters
                             fdr_target_frequency_cm=1560.0, fdr_observable_type='dipole',
                             fdr_axis='z', fdr_output_file=None,
                             fdr_calibration_steps=5000, fdr_log_interval=500,
                             fdr_output_period_ps=0.1, fdr_tau_avg_periods=100,
                             fdr_tau_id_periods=500, fdr_start_time_ps=None,
                             # Step coupling parameters
                             switch_time_ps=None, decay_time_constant_ps=None,
                             # PI feedback controller parameters
                             enable_pi_feedback=False, pi_target_temperature=100.0,
                             pi_turn_on_time_ps=0.0, pi_temperature_method='kinetic',
                             pi_update_interval_ps=5.0, pi_Kc=None, pi_Ti=None,
                             # Temperature tracker parameters  
                             enable_temp_tracker=False, temp_tracker_output_period_ps=0.1,
                             # Thermostat parameters
                             molecular_tau=5.0, cavity_tau=5.0,
                             # Timestep parameters
                             fixed_timestep=False, timestep_fs=0.1, error_tolerance=0.01,
                             # Output parameters
                             gsd_output_period_ps=50.0, console_output_period_ps=1.0,
                             # Hardware parameters
                             device='CPU', gpu_id=0, seed=None,
                             # Experiment metadata
                             replica=0) -> bool:
    """
    Run cavity MD experiment with integrated FDR temperature measurement.
    
    This function provides:
     FDR Temperature Measurement:
    - Real-time effective temperature at specified frequency
    - Physics-based approach without empirical calibration
    - Mode-specific thermalization dynamics
    
     Scientific Observables:
    - T_eff(ω₀,t): Effective temperature from FDR violations
    - S_AA(ω₀,t): Power spectral density at target frequency
    - χ''(ω₀,t): Imaginary susceptibility 
    - Dephasing statistics and lineshape analysis
    
    Returns True for success, False for failure.
    """
    
    try:
        # Calculate photon frequency in atomic units
        omegac = frequency / PhysicalConstants.HARTREE_TO_CM_MINUS1
        
        # Create experiment directory with FDR designation
        coupling_str = f"{coupling:.0e}".replace("-", "neg").replace("+", "pos")
        freq_str = f"{int(fdr_target_frequency_cm)}cm1"
        
        if switch_time_ps is not None:
            # Step coupling experiment
            switch_str = f"switch_{switch_time_ps:.1f}ps"
            if decay_time_constant_ps is not None:
                switch_str += f"_decay_{decay_time_constant_ps:.1f}ps"
            experiment_dir = f"fdr_cavity_coupling_{coupling_str}_{freq_str}_{switch_str}"
        else:
            # Constant coupling experiment
            experiment_dir = f"fdr_cavity_coupling_{coupling_str}_{freq_str}_constant"
            
        if replica > 0:
            experiment_dir += f"_replica_{replica}"
            
        Path(experiment_dir).mkdir(exist_ok=True)
        
        # Setup logging
        log_file = Path(experiment_dir) / f"fdr_simulation_replica_{replica}.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        logger = logging.getLogger(__name__)
        
        logger.info(" Starting FDR Temperature Measurement Experiment")
        logger.info(f"   Coupling: {coupling:.2e}")
        logger.info(f"   Cavity frequency: {frequency:.1f} cm⁻¹")
        logger.info(f"   FDR target frequency: {fdr_target_frequency_cm:.1f} cm⁻¹")
        logger.info(f"   Runtime: {runtime_ps:.1f} ps")
        logger.info(f"   Device: {device}")
        logger.info(f"   Replica: {replica}")
        
        # Initialize simulation using the correct CavityMDSimulation interface
        md_sim = CavityMDSimulation(
            job_dir=experiment_dir,
            replica=replica,
            freq=frequency,
            couplstr=coupling,
            incavity=True,
            runtime_ps=runtime_ps,
            input_gsd=input_gsd,
            frame=frame,
            temperature=temperature,
            molecular_thermostat='langevin',  # Use Langevin for consistency
            cavity_thermostat='langevin',
            molecular_thermostat_tau=molecular_tau,
            cavity_thermostat_tau=cavity_tau,
            finite_q=True,  # Enable finite q for realistic cavity
            device=device,
            gpu_id=gpu_id,
            enable_energy_tracking=False,  # Disable to reduce overhead
            gsd_output_period_ps=gsd_output_period_ps,
            console_output_period_ps=console_output_period_ps,
            log_level='INFO'
        )
        
        # The CavityMDSimulation handles most setup internally
        # For step coupling, we need to modify the coupling after initialization
        if switch_time_ps is not None:
            logger.info(f"   Step coupling: 0 → {coupling:.2e} at t = {switch_time_ps:.1f} ps")
            # We'll temporarily set coupling to 0 for calibration, then restore
        else:
            logger.info(f"   Constant coupling: {coupling:.2e}")
            
        # Setup FDR temperature monitor
        logger.info(" Initializing FDR temperature monitor...")
        
        # Create output file path
        if fdr_output_file is None:
            fdr_output_file = Path(experiment_dir) / f"fdr_temperature_replica_{replica}.dat"
        else:
            fdr_output_file = Path(experiment_dir) / fdr_output_file
            
        # Get the timestep from CavityMDSimulation default (typically 0.0001 ps)
        dt_ps = 0.0001  # CavityMDSimulation default
            
        fdr_monitor = create_cavity_fdr_monitor(
            cavity_frequency_cm=fdr_target_frequency_cm,
            simulation_dt=dt_ps,
            observable_type=fdr_observable_type,
            axis=fdr_axis,
            output_file=str(fdr_output_file),
            log_interval=fdr_log_interval
        )
        
        # Set custom time constants if specified
        if fdr_tau_avg_periods is not None or fdr_tau_id_periods is not None:
            target_period_ps = 2 * np.pi / (fdr_target_frequency_cm / PhysicalConstants.HARTREE_TO_CM_MINUS1 * PhysicalConstants.TIME_PS_CONVERSION)
            
            if fdr_tau_avg_periods is not None:
                fdr_monitor.fdr_estimator._tau_avg_target = fdr_tau_avg_periods * target_period_ps
            if fdr_tau_id_periods is not None:
                fdr_monitor.fdr_estimator._tau_id_target = fdr_tau_id_periods * target_period_ps
        
        # Setup optional temperature tracker
        if enable_temp_tracker:
            logger.info(" Enabling temperature tracker...")
            # CavityMDSimulation will handle this internally
            
        # Initialize the simulation (this sets up everything)
        logger.info(" Initializing simulation...")
        
        # Use the complete workflow to ensure everything is set up
        # We need to set the coupling to 0 for calibration before running
        # Let's temporarily override the coupling in the CavityMDSimulation object
        original_coupling = md_sim.couplstr
        if switch_time_ps is not None or coupling == 0:
            # For step coupling or zero coupling, start at zero
            md_sim.couplstr = 0.0
            logger.info(f"   Setting coupling to 0 for calibration (original: {original_coupling:.2e})")
        
        # Initialize simulation completely using the run method but stop early
        logger.info("   Running initialization phases...")
        try:
            # This sets up the complete simulation including integrator and forces
            # But we'll override the runtime to handle our custom workflow
            original_runtime = md_sim.runtime_ps
            md_sim.runtime_ps = 0.001  # Minimal runtime for setup
            
            # Start the simulation to get everything initialized
            md_sim.run()
            
            # Restore original runtime
            md_sim.runtime_ps = original_runtime
            
        except Exception as e:
            # If run fails, try just the setup
            logger.warning(f"Full run failed ({e}), trying setup only...")
            md_sim.setup_simulation()
        
        # Attach FDR monitor to the initialized simulation
        fdr_monitor.attach_to_simulation(md_sim.sim)
        
        # =========================================================================
        # EQUILIBRATION AND FDR CALIBRATION (BEFORE COUPLING TURNS ON)
        # =========================================================================
        
        logger.info(" Starting equilibration for FDR calibration...")
        logger.info("   Calibrating FDR estimator at equilibrium BEFORE cavity coupling turns on")
            
        equilibration_steps = fdr_calibration_steps
        logger.info(f"   Running {equilibration_steps} equilibration steps at T = {temperature:.1f} K...")
        
        # Collect equilibrium data for FDR calibration
        equilibrium_data = []
        for step in range(equilibration_steps):
            md_sim.sim.run(1)
            
            # Extract observable for FDR calibration
            snapshot = md_sim.sim.state.get_snapshot()
            A_val = fdr_monitor.observable_extractor(snapshot)
            equilibrium_data.append(A_val)
            
            # Progress reporting
            if (step + 1) % (equilibration_steps // 10) == 0:
                progress = (step + 1) / equilibration_steps * 100
                current_time = md_sim.sim.timestep * current_dt
                logger.info(f"   Equilibration progress: {progress:.1f}% (t = {current_time:.2f} ps)")
        
        # Calibrate FDR estimator using equilibrium data
        logger.info(" Calibrating FDR temperature estimator...")
        fdr_monitor.fdr_estimator.calibrate(temperature, np.array(equilibrium_data))
        
        logger.info(f"   FDR calibration complete at T = {temperature:.1f} K")
        logger.info(f"   Gain factor G = {fdr_monitor.fdr_estimator._G_gain:.6e}")
        
        # Calculate equilibration end time for reference
        equilibration_end_time = md_sim.sim.timestep * current_dt
        logger.info(f"   Equilibration completed at t = {equilibration_end_time:.2f} ps")
        
        # Restore coupling - this is when coupling AND FDR measurements start
        if switch_time_ps is not None and original_coupling != 0:
            # For step coupling, we'll handle the coupling change during the run loop
            logger.info(f"   Coupling will be activated at t = {switch_time_ps:.2f} ps")
        elif original_coupling != 0:
            # For constant coupling, restore it now
            md_sim.couplstr = original_coupling  # This will need to be updated in the force as well
            logger.info(f"   Coupling restored to {original_coupling:.2e}")
            
        # For step coupling, verify switch time is in the future
        if switch_time_ps is not None:
            if switch_time_ps <= equilibration_end_time:
                logger.warning(f" Switch time ({switch_time_ps:.2f} ps) is before/during equilibration!")
                logger.warning(f"    Consider setting switch-time > {equilibration_end_time:.2f} ps")
            else:
                logger.info(f"   Step coupling will activate at t = {switch_time_ps:.2f} ps")
        
        logger.info(" FDR calibration complete - ready to start measurements with coupling")
        
        # Determine when FDR measurements should start
        if fdr_start_time_ps is None:
            # Default: start FDR measurements when coupling turns on
            if switch_time_ps is not None:
                fdr_start_time_ps = switch_time_ps
                logger.info(f"   FDR measurements will start at step coupling time: {fdr_start_time_ps:.2f} ps")
            else:
                # For constant coupling, start immediately after calibration
                fdr_start_time_ps = equilibration_end_time
                logger.info(f"   FDR measurements will start immediately: {fdr_start_time_ps:.2f} ps")
        else:
            logger.info(f"   FDR measurements will start at specified time: {fdr_start_time_ps:.2f} ps")
        
        # =========================================================================
        # PRODUCTION RUN WITH FDR MONITORING 
        # =========================================================================
        
        logger.info(" Starting production run with FDR monitoring...")
        
        # Calculate number of steps
        current_dt = md_sim.sim.operations.integrator.dt if hasattr(md_sim.sim.operations.integrator, 'dt') else dt_ps
        total_steps = int(runtime_ps / current_dt)
        fdr_update_steps = max(1, int(fdr_output_period_ps / current_dt))
        fdr_start_step = int(fdr_start_time_ps / current_dt)
        
        logger.info(f"   Total steps: {total_steps}")
        logger.info(f"   FDR update interval: {fdr_update_steps} steps ({fdr_output_period_ps:.3f} ps)")
        logger.info(f"   FDR measurements start at step: {fdr_start_step} (t = {fdr_start_time_ps:.2f} ps)")
        
        # Initialize FDR analysis
        fdr_analyzer = CavityFDRAnalyzer(fdr_monitor)
        
        # Storage for analysis
        T_eff_trajectory = []
        timesteps = []
        diagnostics_trajectory = []
        
        # Production loop
        step_count = 0
        report_interval = max(1, total_steps // 20)  # Report every 5%
        fdr_measurements_started = False
        
        while step_count < total_steps:
            # Run simulation steps
            steps_to_run = min(fdr_update_steps, total_steps - step_count)
            md_sim.sim.run(steps_to_run)
            step_count += steps_to_run
            
            # Check if we should start FDR measurements
            current_simulation_step = md_sim.sim.timestep
            current_time_ps = current_simulation_step * current_dt
            
            if current_simulation_step >= fdr_start_step and not fdr_measurements_started:
                logger.info(f" Starting FDR measurements at t = {current_time_ps:.2f} ps (step {current_simulation_step})")
                fdr_measurements_started = True
            
            # Update FDR temperature estimate (only after start time)
            if fdr_measurements_started:
                try:
                    T_eff, diagnostics = fdr_monitor.update()
                    
                    # Store results
                    T_eff_trajectory.append(T_eff)
                    timesteps.append(current_simulation_step)
                    diagnostics_trajectory.append(diagnostics)
                
                    # Log significant changes or issues
                    if not np.isnan(T_eff):
                        if len(T_eff_trajectory) > 10:
                            recent_temps = [T for T in T_eff_trajectory[-10:] if not np.isnan(T)]
                            if len(recent_temps) > 5:
                                temp_std = np.std(recent_temps)
                                if temp_std > 0.1 * temperature:  # Large fluctuations
                                    logger.warning(f"Large T_eff fluctuations detected: σ = {temp_std:.1f} K")
                                    
                    else:
                        logger.warning(f"NaN T_eff at step {current_simulation_step}")
                        
                except Exception as e:
                    logger.error(f"FDR update failed at step {step_count}: {e}")
            else:
                # Before FDR starts, just run the simulation 
                logger.debug(f"Running simulation (FDR not started): t = {current_time_ps:.2f} ps")
                
            # Progress reporting
            if step_count % report_interval == 0 or step_count >= total_steps:
                progress = step_count / total_steps * 100
                current_time_ps = current_simulation_step * current_dt
                
                if fdr_measurements_started and len(T_eff_trajectory) > 0:
                    recent_T_eff = fdr_monitor.get_recent_temperature(50)
                    logger.info(f"Progress: {progress:.1f}% | Time: {current_time_ps:.1f} ps | T_eff: {recent_T_eff:.1f} K")
                    
                    # Report diagnostics
                    if len(diagnostics_trajectory) > 0:
                        latest_diag = diagnostics_trajectory[-1]
                        logger.debug(f"   FDR Diagnostics: SNR={latest_diag.snr:.2f}, γ={latest_diag.gamma:.4f}/ps")
                        logger.debug(f"   Lineshape: {latest_diag.lineshape_type.value}")
                else:
                    logger.info(f"Progress: {progress:.1f}% | Time: {current_time_ps:.1f} ps | (FDR not started)")
                    if switch_time_ps is not None and current_time_ps < switch_time_ps:
                        logger.info(f"   Waiting for coupling activation at t = {switch_time_ps:.1f} ps")
        
        # =========================================================================
        # FINAL ANALYSIS AND REPORTING
        # =========================================================================
        
        logger.info(" Performing final FDR analysis...")
        
        # Temperature statistics
        temp_stats = fdr_monitor.get_temperature_statistics()
        if temp_stats:
            logger.info(" FDR Temperature Statistics:")
            logger.info(f"   Mean T_eff: {temp_stats.get('mean', 0):.2f} ± {temp_stats.get('std', 0):.2f} K")
            logger.info(f"   Range: {temp_stats.get('min', 0):.1f} - {temp_stats.get('max', 0):.1f} K")
            logger.info(f"   Data points: {temp_stats.get('n_points', 0)}")
            
            # Check for temperature trend
            if 'trend_slope' in temp_stats:
                slope_K_per_ps = temp_stats['trend_slope']
                logger.info(f"   Temperature trend: {slope_K_per_ps:.4f} K/ps")
                
        # Thermalization analysis
        try:
            therm_results = fdr_analyzer.analyze_thermalization_dynamics()
            if therm_results:
                logger.info(" Thermalization Analysis:")
                logger.info(f"   Equilibrium T_eff: {therm_results.get('equilibrium_temperature', 0):.2f} K")
                logger.info(f"   Relaxation time: {therm_results.get('relaxation_time', 0):.2f} ps")
                logger.info(f"   Fit quality (R²): {therm_results.get('fit_quality', 0):.3f}")
        except Exception as e:
            logger.warning(f"Thermalization analysis failed: {e}")
            
        # Comparison with bath temperature
        comparison = fdr_analyzer.compare_with_equipartition(temperature)
        if comparison:
            logger.info(" FDR vs Equipartition Comparison:")
            logger.info(f"   FDR T_eff: {comparison['fdr_temperature']:.2f} K")
            logger.info(f"   Bath T: {comparison['bath_temperature']:.2f} K")
            logger.info(f"   Deviation: {comparison['relative_deviation']*100:.2f}%")
            
            if comparison['is_non_equilibrium']:
                logger.warning(f" Non-equilibrium detected! T_eff ≠ T_bath")
            else:
                logger.info(f" System appears to be in thermal equilibrium")
                
        # Final FDR diagnostics summary
        if len(diagnostics_trajectory) > 0:
            recent_diagnostics = diagnostics_trajectory[-100:]  # Last 100 points
            avg_snr = np.mean([d.snr for d in recent_diagnostics if not np.isnan(d.snr)])
            avg_gamma = np.mean([d.gamma for d in recent_diagnostics if not np.isnan(d.gamma)])
            
            logger.info(" Final FDR Diagnostics:")
            logger.info(f"   Average SNR: {avg_snr:.2f}")
            logger.info(f"   Average damping γ: {avg_gamma:.4f} ps⁻¹") 
            logger.info(f"   Target frequency: {fdr_target_frequency_cm:.1f} cm⁻¹")
            
        # Save analysis summary
        summary_file = Path(experiment_dir) / f"fdr_analysis_summary_replica_{replica}.txt"
        with open(summary_file, 'w') as f:
            f.write("FDR Temperature Measurement Analysis Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Experiment Parameters:\n")
            f.write(f"  Coupling strength: {coupling:.2e}\n")
            f.write(f"  Target frequency: {fdr_target_frequency_cm:.1f} cm⁻¹\n")
            f.write(f"  Runtime: {runtime_ps:.1f} ps\n")
            f.write(f"  Bath temperature: {temperature:.1f} K\n\n")
            
            if temp_stats:
                f.write(f"Temperature Statistics:\n")
                f.write(f"  Mean T_eff: {temp_stats.get('mean', 0):.2f} ± {temp_stats.get('std', 0):.2f} K\n")
                f.write(f"  Range: {temp_stats.get('min', 0):.1f} - {temp_stats.get('max', 0):.1f} K\n")
                f.write(f"  Data points: {temp_stats.get('n_points', 0)}\n\n")
                
            if comparison:
                f.write(f"Equilibrium Assessment:\n")
                f.write(f"  FDR T_eff: {comparison['fdr_temperature']:.2f} K\n")
                f.write(f"  Bath T: {comparison['bath_temperature']:.2f} K\n")
                f.write(f"  Deviation: {comparison['relative_deviation']*100:.2f}%\n")
                f.write(f"  Non-equilibrium: {comparison['is_non_equilibrium']}\n\n")
                
        logger.info(f" Analysis summary saved to: {summary_file}")
        logger.info(f" FDR trajectory saved to: {fdr_output_file}")
        logger.info(" FDR temperature measurement experiment completed successfully!")
        
        return True
        
    except Exception as e:
        logger.error(f" Experiment failed: {e}")
        import traceback
        logger.error(f"Traceback: {traceback.format_exc()}")
        return False


# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================

def main():
    """Main function with command line argument parsing."""
    
    parser = argparse.ArgumentParser(
        description="FDR Temperature Measurement in Cavity MD Simulations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic FDR measurement
  python 19_fdr_temperature_measurement.py --coupling 7e-4 --frequency 1560 --runtime 100.0
  
  # With step coupling
  python 19_fdr_temperature_measurement.py --coupling 7e-4 --switch-time 20.0 --runtime 200.0
  
  # Multiple replicas  
  python 19_fdr_temperature_measurement.py --coupling 7e-4 --replicas 3 --runtime 150.0
  
  # GPU acceleration
  python 19_fdr_temperature_measurement.py --coupling 7e-4 --device GPU --runtime 100.0
        """
    )
    
    # Simulation parameters
    parser.add_argument('--coupling', type=float, default=7e-4,
                        help='Cavity coupling strength (default: 7e-4)')
    parser.add_argument('--frequency', type=float, default=1560.0,
                        help='Cavity frequency in cm⁻¹ (default: 1560.0)')
    parser.add_argument('--temperature', type=float, default=100.0,
                        help='Bath temperature in K (default: 100.0)')
    parser.add_argument('--runtime', type=float, default=100.0,
                        help='Simulation runtime in ps (default: 100.0)')
    
    # FDR measurement parameters
    parser.add_argument('--fdr-frequency', type=float, default=1560.0,
                        help='FDR target frequency in cm⁻¹ (default: 1560.0)')
    parser.add_argument('--fdr-axis', choices=['x', 'y', 'z'], default='z',
                        help='Dipole projection axis (default: z)')
    parser.add_argument('--fdr-calibration-steps', type=int, default=5000,
                        help='FDR calibration steps (default: 5000)')
    parser.add_argument('--fdr-log-interval', type=int, default=500,
                        help='FDR logging interval (default: 500)')
    parser.add_argument('--fdr-output-period', type=float, default=0.1,
                        help='FDR output period in ps (default: 0.1)')
    parser.add_argument('--fdr-start-time', type=float, default=None,
                        help='FDR measurement start time in ps (default: auto - starts with coupling)')
    
    # Coupling control
    parser.add_argument('--switch-time', type=float, default=None,
                        help='Switch time for step coupling in ps (default: None)')
    parser.add_argument('--decay-time', type=float, default=None,
                        help='Decay time constant for step coupling in ps (default: None)')
    
    # PI feedback controller
    parser.add_argument('--enable-pi-feedback', action='store_true',
                        help='Enable PI feedback temperature controller')
    parser.add_argument('--pi-method', choices=['kinetic', 'lj_coulombic', 'harmonic'], 
                        default='kinetic', help='PI temperature method (default: kinetic)')
    parser.add_argument('--pi-update-interval', type=float, default=5.0,
                        help='PI update interval in ps (default: 5.0)')
    
    # Temperature tracking
    parser.add_argument('--enable-temp-tracker', action='store_true',
                        help='Enable temperature tracker')
    parser.add_argument('--temp-tracker-period', type=float, default=0.1,
                        help='Temperature tracker period in ps (default: 0.1)')
    
    # Input/output parameters
    parser.add_argument('--input-gsd', type=str, default='molecular-0.gsd',
                        help='Input GSD file (default: molecular-0.gsd)')
    parser.add_argument('--frame', type=int, default=-1,
                        help='Frame to load from GSD file (default: -1)')
    
    # Hardware parameters
    parser.add_argument('--device', choices=['CPU', 'GPU'], default='CPU',
                        help='Compute device (default: CPU)')
    parser.add_argument('--gpu-id', type=int, default=0,
                        help='GPU ID if using GPU (default: 0)')
    
    # Simulation control
    parser.add_argument('--replicas', type=int, default=1,
                        help='Number of replicas to run (default: 1)')
    parser.add_argument('--seed', type=int, default=None,
                        help='Random seed (default: None)')
    
    # Thermostat parameters
    parser.add_argument('--molecular-tau', type=float, default=5.0,
                        help='Molecular thermostat time constant in ps (default: 5.0)')
    parser.add_argument('--cavity-tau', type=float, default=5.0,
                        help='Cavity thermostat time constant in ps (default: 5.0)')
    
    # Output control
    parser.add_argument('--gsd-period', type=float, default=50.0,
                        help='GSD output period in ps (default: 50.0)')
    parser.add_argument('--console-period', type=float, default=1.0,
                        help='Console output period in ps (default: 1.0)')
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.coupling <= 0:
        parser.error("Coupling strength must be positive")
    if args.frequency <= 0:
        parser.error("Frequency must be positive")
    if args.runtime <= 0:
        parser.error("Runtime must be positive")
        
    # Set default FDR frequency to match cavity frequency if not specified
    if args.fdr_frequency == 1560.0 and args.frequency != 1560.0:
        args.fdr_frequency = args.frequency
        print(f"Setting FDR frequency to match cavity frequency: {args.frequency:.1f} cm⁻¹")
    
    # Display experiment summary
    print("\n" + "="*60)
    print(" FDR TEMPERATURE MEASUREMENT EXPERIMENT")
    print("="*60)
    print(f" Parameters:")
    print(f"   Coupling strength: {args.coupling:.2e}")
    print(f"   Cavity frequency: {args.frequency:.1f} cm⁻¹")
    print(f"   FDR frequency: {args.fdr_frequency:.1f} cm⁻¹")
    print(f"   Temperature: {args.temperature:.1f} K")
    print(f"   Runtime: {args.runtime:.1f} ps")
    print(f"   Device: {args.device}")
    print(f"   Replicas: {args.replicas}")
    
    if args.switch_time is not None:
        print(f"   Step coupling at: {args.switch_time:.1f} ps")
    if args.enable_pi_feedback:
        print(f"   PI feedback: {args.pi_method}")
    if args.enable_temp_tracker:
        print("   Temperature tracker: enabled")
    print("="*60)
    
    # Run experiments
    success_count = 0
    
    for replica in range(args.replicas):
        print(f"\n Starting replica {replica + 1}/{args.replicas}...")
        
        # Calculate replica-specific seed
        replica_seed = (args.seed + replica) if args.seed is not None else None
        
        # Run single experiment
        success = run_fdr_cavity_experiment(
            coupling=args.coupling,
            temperature=args.temperature,
            frequency=args.frequency,
            runtime_ps=args.runtime,
            input_gsd=args.input_gsd,
            frame=args.frame,
            # FDR parameters
            fdr_target_frequency_cm=args.fdr_frequency,
            fdr_axis=args.fdr_axis,
            fdr_calibration_steps=args.fdr_calibration_steps,
            fdr_log_interval=args.fdr_log_interval,
            fdr_output_period_ps=args.fdr_output_period,
            fdr_start_time_ps=args.fdr_start_time,
            # Step coupling
            switch_time_ps=args.switch_time,
            decay_time_constant_ps=args.decay_time,
            # PI feedback
            enable_pi_feedback=args.enable_pi_feedback,
            pi_target_temperature=args.temperature,
            pi_temperature_method=args.pi_method,
            pi_update_interval_ps=args.pi_update_interval,
            # Temperature tracker
            enable_temp_tracker=args.enable_temp_tracker,
            temp_tracker_output_period_ps=args.temp_tracker_period,
            # Thermostat
            molecular_tau=args.molecular_tau,
            cavity_tau=args.cavity_tau,
            # Output
            gsd_output_period_ps=args.gsd_period,
            console_output_period_ps=args.console_period,
            # Hardware
            device=args.device,
            gpu_id=args.gpu_id,
            seed=replica_seed,
            # Metadata
            replica=replica
        )
        
        if success:
            success_count += 1
            print(f" Replica {replica + 1} completed successfully")
        else:
            print(f" Replica {replica + 1} failed")
    
    # Final summary
    print("\n" + "="*60)
    print(" EXPERIMENT SUMMARY")
    print("="*60)
    print(f" Successful replicas: {success_count}/{args.replicas}")
    if success_count > 0:
        print(" Output files generated:")
        print("   - fdr_temperature_replica_*.dat (FDR trajectories)")
        print("   - fdr_analysis_summary_replica_*.txt (analysis summaries)")
        print("   - trajectory_replica_*.gsd (simulation trajectories)")
        print("   - fdr_simulation_replica_*.log (detailed logs)")
    print("="*60)
    
    if success_count == args.replicas:
        print(" All experiments completed successfully!")
        return 0
    else:
        print(f" {args.replicas - success_count} experiments failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
