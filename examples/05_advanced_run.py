#!/usr/bin/env python3
# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""
Advanced Cavity Molecular Dynamics Experiment Runner

This script provides a comprehensive framework for running cavity MD simulations
using the CavityMDSimulation class from the hoomd.cavitymd plugin. It supports
various thermostat combinations and advanced analysis features.

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
   # Run a single experiment with cavity coupling
   python 05_advanced_run.py --molecular-bath bussi --cavity-bath langevin --coupling 1e-3 --runtime 1000
   
   # Run without cavity (molecular-only simulation)
   python 05_advanced_run.py --molecular-bath bussi --no-cavity --runtime 1000
   
   # Run multiple replicas locally
   python 05_advanced_run.py --molecular-bath bussi --cavity-bath langevin --coupling 1e-3 --replicas 1-5 --runtime 500
   
   # Run with finite-q cavity mode
   python 05_advanced_run.py --molecular-bath bussi --cavity-bath langevin --finite-q --coupling 1e-3 --runtime 1000

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
- Local execution: --replicas 1-5 or --replicas 1,2,3,4,5
- SLURM array jobs: Automatically detects SLURM_ARRAY_TASK_ID

ADVANCED FEATURES:
- --enable-energy-tracker: Detailed energy component tracking
- --enable-fkt: F(k,t) density correlation functions
- --fixed-timestep: Use fixed timestep instead of adaptive
- --device GPU: Run on GPU instead of CPU
- --seed: Control random seed for reproducibility

OUTPUT CONTROL:
- Separate output periods for different observables (energy, F(k,t), trajectories, console)
- Organized directory structure
- Comprehensive logging with timestamps

See --help for all available options.
"""

import sys
import argparse
import time
from pathlib import Path

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
                         switch_time_ps=None, damping_ratio=0.0):
    """
    Run a single experiment using the CavityMDSimulation class from the plugin.
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
            if switch_time_ps is not None:
                # Include switch time in directory name for time-varying simulations
                switch_str = f"_switch_{switch_time_ps}ps"
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
            print(f"  Coupling strength: {coupling}")
            if switch_time_ps is not None:
                print(f"  Switch time: {switch_time_ps} ps")
                print(f"  Damping ratio (zeta): {damping_ratio}")
                print(f"  Calculated dissipation (c): {dissipation:.4e} a.u.")
            else:
                print(f"  Damping ratio (zeta): {damping_ratio}")
                print(f"  Calculated dissipation (c): {dissipation:.4e} a.u.")
            print(f"  Molecular thermostat: {molecular_thermo}")
            print(f"  Cavity thermostat: {cavity_thermo}")
            print(f"  Finite-q mode: {finite_q}")
        else:
            print(f"  Molecular thermostat: {molecular_thermo}")
        print(f"  Replica: {replica}")
        print(f"  Frame: {frame}")
        print(f"  Output directory: {exp_dir}")
        
        # Set error tolerance based on timestepping mode
        error_tolerance = 0.0 if fixed_timestep else 0.5
        
        # Set timestep based on user preference (only used if fixed_timestep is True)
        dt_fs = timestep_fs if fixed_timestep else None
        
        # Create and run CavityMDSimulation
        sim = CavityMDSimulation(
            job_dir=str(exp_dir),
            replica=replica,
            freq=frequency,
            couplstr=coupling,
            incavity=incavity,
            runtime_ps=runtime_ps,
            input_gsd='../init-0.gsd',  # Use relative path from job directory back to parent
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
            dissipation=dissipation
        )
        
        # Run the simulation
        return sim.run() == 0  # Return True for success (exit code 0)
        
    except Exception as e:
        print(f"ERROR: Experiment failed: {e}")
        return False

def main():
    """Simplified main function for cavity MD experiments."""
    parser = argparse.ArgumentParser(
        description='Advanced Cavity MD Experiment Runner',
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
                       help='Cavity coupling strength (default: 1e-3)')
    parser.add_argument('--switch-time', type=float, 
                       help='Time in ps when coupling and dissipation turn on (default: on from start)')
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
    
    # Replica control
    parser.add_argument('--replicas', type=str, 
                       help='Replica specification (e.g., "1,2,3" or "1-5")')
    
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
    
    args = parser.parse_args()
    
    print("Advanced Cavity MD Experiment Runner")
    print("="*50)
    
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
    
    print(f"\nSimulation Configuration:")
    print(f"  Cavity coupling: {'Enabled' if incavity else 'Disabled'}")
    if incavity:
        print(f"    Coupling strength: {args.coupling}")
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
    
    # Set up device configuration
    device = args.device.upper()
    
    # Performance tracking
    start_time = time.time()
    successful_experiments = 0
    failed_experiments = 0
    
    print(f"\nStarting execution for {len(replica_list)} replica(s)...")
    print("="*50)
    
    # Run replicas
    for replica in replica_list:
        frame = replica  # Use replica as frame number
        
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
            damping_ratio=args.damping_ratio
        )
        
        if success:
            successful_experiments += 1
            print(f"SUCCESS: Replica {replica} completed successfully")
        else:
            failed_experiments += 1
            print(f"ERROR: Replica {replica} failed")
    
    # Final summary
    end_time = time.time()
    total_wall_time = end_time - start_time
    total_experiments = len(replica_list)
    
    print("\n" + "="*50)
    print("Execution Summary")
    print("="*50)
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