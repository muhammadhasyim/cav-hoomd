#!/usr/bin/env python3
"""
LQG Controller Example Script

This script demonstrates how to use the LQG pulse-to-pulse controller
with the cavity-matter system. It shows the complete workflow from
setup to execution.

Usage:
    python lqg_controller_example.py --mode basic
    python lqg_controller_example.py --mode advanced --system-id
    python lqg_controller_example.py --mode validation --params my_params.json
"""

import argparse
import numpy as np
from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

# Import HOOMD and cavity MD components
import hoomd
from cavitymd.simulation import CavityMDSimulation
from cavitymd.lqg_controller import LQGParameters
from cavitymd.lqg_pulse_controller import LQGPulseController
from cavitymd.lqg_square_wave import LQGSquareWaveVariant


def create_basic_lqg_simulation():
    """
    Create a basic LQG-controlled cavity MD simulation.
    
    This example shows the minimal setup required to use the LQG controller
    with default parameters.
    """
    print("Creating basic LQG-controlled simulation...")
    
    # Create cavity MD simulation
    sim = CavityMDSimulation(
        job_dir="lqg_basic_example",
        replica=1,
        freq=2000.0,  # 2000 cm⁻¹ cavity frequency
        couplstr=1e-4,  # This will be overridden by LQG controller
        incavity=True,
        runtime_ps=1000.0,
        temperature=300.0,
        molecular_thermostat='bussi',
        cavity_thermostat='langevin',
        enable_energy_tracking=True,
        enable_temp_tracker=True,
        device='CPU'
    )
    
    # Set up simulation (this creates the HOOMD simulation object)
    snapshot = sim.setup_simulation()
    
    # Create LQG controller
    lqg_controller = LQGPulseController(
        time_tracker=sim.time_tracker,
        energy_tracker=sim.energy_tracker,
        simulation=sim.sim,
        molecular_thermostat=sim.molecular_thermostat_obj,
        cavity_thermostat=sim.cavity_thermostat_obj,
        temperature_tracker=getattr(sim, 'temperature_tracker', None),
        
        # Target temperatures
        target_temperatures={
            'structure': 280.0,    # Target structural temperature
            'kinetic': 300.0,      # Target kinetic temperature
            'vibrational': 320.0   # Target vibrational temperature
        },
        
        # Pulse parameters
        pulse_period_ps=5.0,
        pulse_duty_cycle=0.5,
        baseline_coupling=1e-5,
        
        # Temperature measurement methods
        structure_temp_method='lj_coulombic',
        kinetic_temp_method='kinetic',
        vibrational_temp_method='molecular_vibrational',
        
        # Control parameters
        turn_on_time_ps=50.0,  # Allow equilibration first
        apply_to='both',
        output_file='lqg_basic_output.csv'
    )
    
    # Create LQG-controlled square wave
    lqg_square_wave = LQGSquareWaveVariant(
        lqg_controller=lqg_controller,
        period_ps=5.0,
        time_tracker=sim.time_tracker,
        duty_cycle=0.5,
        start_time_ps=50.0
    )
    
    # Update cavity force to use LQG-controlled coupling
    if hasattr(sim, 'cavity_force'):
        sim.cavity_force.couplstr = lqg_square_wave
        print("Updated cavity force with LQG-controlled coupling")
    
    # Add LQG controller as updater
    sim.sim.operations.updaters.append(
        hoomd.update.CustomUpdater(
            action=lqg_controller,
            trigger=hoomd.trigger.Periodic(1)  # Check every timestep
        )
    )
    
    print("Basic LQG simulation setup complete!")
    return sim, lqg_controller, lqg_square_wave


def create_advanced_lqg_simulation(enable_system_id=True, custom_params=None):
    """
    Create an advanced LQG simulation with custom parameters and system identification.
    
    Parameters
    ----------
    enable_system_id : bool
        Whether to enable system identification
    custom_params : LQGParameters, optional
        Custom LQG parameters
    """
    print("Creating advanced LQG-controlled simulation...")
    
    # Create custom LQG parameters if not provided
    if custom_params is None:
        custom_params = LQGParameters()
        
        # Customize parameters for this example
        custom_params.Q_structure = 5.0      # Higher structure tracking weight
        custom_params.Q_kinetic = 1.0        # Standard kinetic tracking
        custom_params.Q_vibrational = 0.5    # Lower vibrational tracking
        
        custom_params.R_coupling = 0.05      # Lower coupling penalty (more aggressive)
        custom_params.R_bath = 0.01          # Lower bath penalty
        
        # Adjust model parameters based on expected physics
        custom_params.a_sk = 0.18            # Structure-kinetic exchange
        custom_params.a_sv = 0.08            # Structure-vibrational exchange
        custom_params.b_sg = 0.12            # Coupling cooling effect
        custom_params.b_vg = -0.15           # Coupling heating effect on vibrations
        
        # Process noise (model uncertainty)
        custom_params.W_structure = 2e-4
        custom_params.W_kinetic = 1e-4
        custom_params.W_vibrational = 3e-4
        
        # Measurement noise
        custom_params.V_structure = 1e-3     # LJ+Coulombic has more noise
        custom_params.V_kinetic = 5e-4       # Kinetic temp is cleaner
        custom_params.V_vibrational = 2e-3   # Vibrational temp is noisy
    
    # Create cavity MD simulation
    sim = CavityMDSimulation(
        job_dir="lqg_advanced_example",
        replica=1,
        freq=1800.0,  # Different frequency
        couplstr=2e-4,  # Will be overridden
        incavity=True,
        runtime_ps=2000.0,
        temperature=300.0,
        molecular_thermostat='bussi',
        cavity_thermostat='langevin',
        enable_energy_tracking=True,
        enable_temp_tracker=True,
        enable_molecular_temps=True,  # Enable molecular temperature decomposition
        device='CPU'
    )
    
    # Set up simulation
    snapshot = sim.setup_simulation()
    
    # Create advanced LQG controller
    lqg_controller = LQGPulseController(
        time_tracker=sim.time_tracker,
        energy_tracker=sim.energy_tracker,
        simulation=sim.sim,
        molecular_thermostat=sim.molecular_thermostat_obj,
        cavity_thermostat=sim.cavity_thermostat_obj,
        temperature_tracker=getattr(sim, 'molecular_temp_tracker', None),
        
        # Target temperatures with different setpoints
        target_temperatures={
            'structure': 275.0,
            'kinetic': 300.0,
            'vibrational': 325.0
        },
        
        # Pulse parameters
        pulse_period_ps=8.0,  # Longer period
        pulse_duty_cycle=0.3,  # Lower duty cycle
        baseline_coupling=2e-5,
        
        # Temperature methods
        structure_temp_method='lj_coulombic',
        kinetic_temp_method='molecular_kinetic',
        vibrational_temp_method='molecular_vibrational',
        
        # Custom parameters
        lqg_params=custom_params,
        
        # System identification
        enable_system_id=enable_system_id,
        system_id_duration_ps=200.0,
        system_id_file='advanced_system_id.json',
        
        # Online adaptation
        enable_online_adaptation=True,
        adaptation_update_interval=25,
        adaptation_forgetting_factor=0.995,
        
        # Control limits
        max_coupling_deviation=0.5,
        min_bath_temperature=50.0,
        max_bath_temperature=400.0,
        
        # Timing and output
        turn_on_time_ps=100.0,
        apply_to='both',
        output_file='lqg_advanced_output.csv',
        console_output_period_ps=20.0,
        empirical_data_file=getattr(sim, 'empirical_data_file', None)
    )
    
    # Create adaptive LQG square wave
    from cavitymd.lqg_square_wave import LQGAdaptiveSquareWaveVariant
    
    lqg_square_wave = LQGAdaptiveSquareWaveVariant(
        lqg_controller=lqg_controller,
        period_ps=8.0,
        time_tracker=sim.time_tracker,
        duty_cycle=0.3,
        start_time_ps=100.0,
        enable_period_adaptation=True,
        enable_duty_cycle_adaptation=False,
        min_period_ps=3.0,
        max_period_ps=15.0,
        adaptation_interval_pulses=50
    )
    
    # Update cavity force
    if hasattr(sim, 'cavity_force'):
        sim.cavity_force.couplstr = lqg_square_wave
        print("Updated cavity force with adaptive LQG-controlled coupling")
    
    # Add LQG controller as updater
    sim.sim.operations.updaters.append(
        hoomd.update.CustomUpdater(
            action=lqg_controller,
            trigger=hoomd.trigger.Periodic(1)
        )
    )
    
    print("Advanced LQG simulation setup complete!")
    print(f"System identification: {'Enabled' if enable_system_id else 'Disabled'}")
    print(f"Online adaptation: Enabled")
    print(f"Pulse period: {lqg_controller.pulse_period_ps} ps")
    print(f"Target temperatures: {lqg_controller.target_temperatures}")
    
    return sim, lqg_controller, lqg_square_wave


def run_validation_simulation(params_file):
    """
    Run a validation simulation with pre-trained parameters.
    
    Parameters
    ----------
    params_file : str
        Path to JSON file with LQG parameters
    """
    print(f"Running validation simulation with parameters from {params_file}")
    
    # Load parameters
    try:
        custom_params = LQGParameters.load(params_file)
        print("Successfully loaded custom parameters")
    except Exception as e:
        print(f"Warning: Could not load parameters ({e}), using defaults")
        custom_params = LQGParameters()
    
    # Create simulation with loaded parameters
    sim, lqg_controller, lqg_square_wave = create_advanced_lqg_simulation(
        enable_system_id=False,  # Skip system ID for validation
        custom_params=custom_params
    )
    
    # Set up validation-specific settings
    lqg_controller.turn_on_time_ps = 20.0  # Start control earlier
    lqg_controller.output_file = 'lqg_validation_output.csv'
    
    print("Validation simulation setup complete!")
    return sim, lqg_controller, lqg_square_wave


def demonstrate_controller_features(lqg_controller):
    """
    Demonstrate key features of the LQG controller.
    
    Parameters
    ----------
    lqg_controller : LQGPulseController
        The LQG controller instance
    """
    print("\n" + "="*50)
    print("LQG Controller Features Demonstration")
    print("="*50)
    
    # Show current parameters
    print("\nCurrent LQG Parameters:")
    print(f"  Structure tracking weight: {lqg_controller.lqg_params.Q_structure}")
    print(f"  Kinetic tracking weight: {lqg_controller.lqg_params.Q_kinetic}")
    print(f"  Vibrational tracking weight: {lqg_controller.lqg_params.Q_vibrational}")
    print(f"  Coupling penalty: {lqg_controller.lqg_params.R_coupling}")
    print(f"  Bath penalty: {lqg_controller.lqg_params.R_bath}")
    
    # Show control settings
    print("\nControl Settings:")
    print(f"  Pulse period: {lqg_controller.pulse_period_ps} ps")
    print(f"  Duty cycle: {lqg_controller.pulse_duty_cycle}")
    print(f"  Baseline coupling: {lqg_controller.lqg_params.g_baseline:.2e}")
    print(f"  Max coupling deviation: ±{lqg_controller.lqg_params.delta_g_max}")
    print(f"  Bath temperature range: {lqg_controller.lqg_params.T_bath_min}-{lqg_controller.lqg_params.T_bath_max} K")
    
    # Show target temperatures
    print("\nTarget Temperatures:")
    for temp_type, target in lqg_controller.target_temperatures.items():
        print(f"  {temp_type.capitalize()}: {target} K")
    
    # Show measurement methods
    print("\nTemperature Measurement Methods:")
    print(f"  Structure: {lqg_controller.structure_temp_method}")
    print(f"  Kinetic: {lqg_controller.kinetic_temp_method}")
    print(f"  Vibrational: {lqg_controller.vibrational_temp_method}")
    
    # Show system identification settings
    print("\nSystem Identification:")
    print(f"  Enabled: {lqg_controller.enable_system_id}")
    if lqg_controller.enable_system_id:
        print(f"  Duration: {lqg_controller.system_id_duration_ps} ps")
        print(f"  Output file: {lqg_controller.system_id_file}")
    
    # Show online adaptation settings
    print("\nOnline Adaptation:")
    print(f"  Enabled: {lqg_controller.enable_online_adaptation}")
    if lqg_controller.enable_online_adaptation:
        print(f"  Update interval: {lqg_controller.adaptation_update_interval} pulses")
        print(f"  Forgetting factor: {lqg_controller.adaptation_forgetting_factor}")
    
    print("\n" + "="*50)


def run_simulation_example(sim, lqg_controller, lqg_square_wave, steps=10000):
    """
    Run a short simulation example to demonstrate the controller.
    
    Parameters
    ----------
    sim : CavityMDSimulation
        The simulation object
    lqg_controller : LQGPulseController
        The LQG controller
    lqg_square_wave : LQGSquareWaveVariant
        The square wave variant
    steps : int
        Number of simulation steps to run
    """
    print(f"\nRunning simulation for {steps} steps...")
    
    try:
        # Run simulation
        sim.sim.run(steps)
        
        # Show final controller state
        print("\nFinal Controller State:")
        print(f"  Pulse count: {lqg_controller.pulse_count}")
        print(f"  Current coupling: {lqg_controller.get_current_coupling():.2e}")
        print(f"  Current bath temperature: {lqg_controller.get_current_bath_temperature():.1f} K")
        print(f"  System ID complete: {lqg_controller.system_id_complete}")
        
        # Show square wave state
        state_info = lqg_square_wave.get_state_info()
        print(f"  Square wave active: {state_info['is_active']}")
        print(f"  Currently high: {state_info['is_high']}")
        
        print(f"\nSimulation completed successfully!")
        print(f"Output saved to: {lqg_controller.output_file}")
        
    except Exception as e:
        print(f"Simulation failed: {e}")
        import traceback
        traceback.print_exc()


def main():
    """Main function for command-line interface."""
    parser = argparse.ArgumentParser(description="LQG Controller Example")
    parser.add_argument('--mode', choices=['basic', 'advanced', 'validation'], 
                       default='basic', help='Example mode')
    parser.add_argument('--params', help='Parameters file for validation mode')
    parser.add_argument('--system-id', action='store_true', 
                       help='Enable system identification (advanced mode)')
    parser.add_argument('--steps', type=int, default=5000, 
                       help='Number of simulation steps')
    parser.add_argument('--no-run', action='store_true',
                       help='Set up simulation but do not run')
    
    args = parser.parse_args()
    
    print("LQG Controller Example")
    print("=" * 30)
    
    # Create simulation based on mode
    if args.mode == 'basic':
        sim, lqg_controller, lqg_square_wave = create_basic_lqg_simulation()
        
    elif args.mode == 'advanced':
        sim, lqg_controller, lqg_square_wave = create_advanced_lqg_simulation(
            enable_system_id=args.system_id
        )
        
    elif args.mode == 'validation':
        if not args.params:
            print("Error: --params required for validation mode")
            return
        
        sim, lqg_controller, lqg_square_wave = run_validation_simulation(args.params)
    
    # Demonstrate controller features
    demonstrate_controller_features(lqg_controller)
    
    # Run simulation if requested
    if not args.no_run:
        run_simulation_example(sim, lqg_controller, lqg_square_wave, args.steps)
    else:
        print("\nSimulation setup complete (not running due to --no-run flag)")
        print("To run manually:")
        print("  sim.sim.run(10000)")


if __name__ == '__main__':
    main()
