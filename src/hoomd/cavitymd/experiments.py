# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Experiment runner and utilities for cavity molecular dynamics simulations."""

import sys
import os
import argparse
from pathlib import Path

from .simulation import CavityMDSimulation


# Available experiments: (name, molecular_thermostat, cavity_thermostat, finite_q)
BUSSI_LANGEVIN_EXPERIMENTS = [
    ("bussi_langevin_finiteq", "bussi", "langevin", True),
    ("bussi_langevin_no_finiteq", "bussi", "langevin", False),
]


def parse_replicas(replicas_str):
    """Parse replica specification string into list of integers.
    
    Args:
        replicas_str: String like "1-10", "1,3,5,7", etc.
        
    Returns:
        List of replica numbers (integers)
    """
    replicas = []
    
    # Handle comma-separated values
    parts = replicas_str.split(',')
    for part in parts:
        part = part.strip()
        if '-' in part:
            # Handle range like "1-10"
            start, end = part.split('-', 1)
            start, end = int(start.strip()), int(end.strip())
            replicas.extend(range(start, end + 1))
        else:
            # Handle single number
            replicas.append(int(part))
    
    return sorted(list(set(replicas)))  # Remove duplicates and sort


def get_slurm_info():
    """Get SLURM job information from environment variables."""
    task_id = os.environ.get('SLURM_ARRAY_TASK_ID')
    job_id = os.environ.get('SLURM_JOB_ID', 'unknown')
    
    if task_id is None:
        return None, None, job_id
    else:
        replica = int(task_id)
        frame = int(task_id)
        print(f"SLURM Array Job ID: {job_id}")
        print(f"SLURM Array Task ID: {task_id} (using as replica and frame number)")
        return replica, frame, job_id


def run_single_experiment(exp_name, molecular_thermo, cavity_thermo, finite_q, coupling, 
                         replica, frame, runtime_ps, molecular_tau, cavity_tau, log_to_file, 
                         log_to_console, enable_fkt, fkt_kmag, fkt_wavevectors, fkt_ref_interval, 
                         fkt_max_refs, max_energy_output_time=None, device='CPU', gpu_id=0, 
                         incavity=True, fixed_timestep=False, timestep_fs=1.0, 
                         enable_energy_tracking=True):
    """Run a single experiment directly (no subprocess)."""
    
    try:
        # Create experiment directory with appropriate naming
        if incavity:
            # For cavity simulations, include coupling strength in directory name
            coupling_str = f"{coupling:.0e}".replace("-", "neg").replace("+", "pos")
            exp_dir = Path(f"{exp_name}_coupling_{coupling_str}")
        else:
            # For non-cavity simulations, coupling doesn't matter - use simple naming
            exp_dir = Path(f"{exp_name}_no_cavity")
        exp_dir.mkdir(exist_ok=True)
        
        print(f"Running experiment: {exp_name}")
        print(f"Cavity coupling: {'Enabled' if incavity else 'Disabled'}")
        if incavity:
            print(f"Coupling strength: {coupling}")
        print(f"Replica: {replica}")
        print(f"Frame: {frame}")
        print(f"Output directory: {exp_dir}")
        
        # Set error tolerance based on timestepping mode
        error_tolerance = 0.0 if fixed_timestep else 1.0
        
        # Set timestep based on user preference (only used if fixed_timestep is True)
        dt_fs = timestep_fs if fixed_timestep else None
        
        # Run simulation
        sim = CavityMDSimulation(
            job_dir=str(exp_dir),
            replica=replica,
            freq=1560.0,
            couplstr=coupling,
            incavity=incavity,
            runtime_ps=runtime_ps,
            input_gsd='../init-0.gsd',
            frame=frame,
            name='prod',
            error_tolerance=error_tolerance,
            temperature=100.0,
            molecular_thermostat=molecular_thermo,
            cavity_thermostat=cavity_thermo,
            cavity_damping_factor=1.0,
            use_brownian_overdamped=False,
            add_cavity_particle=True,
            finite_q=finite_q,
            molecular_thermostat_tau=molecular_tau,
            cavity_thermostat_tau=cavity_tau,
            log_to_file=log_to_file,
            log_to_console=log_to_console,
            log_level='INFO',
            enable_fkt=enable_fkt,
            fkt_kmag=fkt_kmag,
            fkt_num_wavevectors=fkt_wavevectors,
            fkt_reference_interval_ps=fkt_ref_interval,
            fkt_max_references=fkt_max_refs,
            max_energy_output_time_ps=max_energy_output_time,
            enable_energy_tracking=enable_energy_tracking,
            dt_fs=dt_fs,
            device=device,
            gpu_id=gpu_id
        )
        
        print(f"🔄 Starting simulation for replica {replica}...")
        exit_code = sim.run()
        print(f"🏁 Simulation completed for replica {replica} with exit code {exit_code}")
        
        if exit_code == 0:
            print(f"✅ {exp_name} completed successfully (replica {replica}, frame {frame})")
            return True
        else:
            print(f"❌ {exp_name} failed with exit code {exit_code} (replica {replica}, frame {frame})")
            return False
            
    except Exception as e:
        print(f"❌ {exp_name} - error: {e}")
        import traceback
        traceback.print_exc()
        return False


def run_cavity_experiments():
    """Main function for running cavity experiments - command line interface."""
    parser = argparse.ArgumentParser(description="Run cavity MD experiments with SLURM or local loops")
    parser.add_argument('--coupling', type=float, required=True, help='Coupling strength (required)')
    parser.add_argument('--runtime', type=float, default=500.0, help='Runtime in ps (default: 500)')
    parser.add_argument('--experiment', type=str, default='bussi_langevin_finiteq', 
                       choices=[name for name, _, _, _ in BUSSI_LANGEVIN_EXPERIMENTS],
                       help='Specific experiment to run (default: bussi_langevin_finiteq)')
    parser.add_argument('--molecular-tau', type=float, default=5.0, help='Molecular thermostat time constant in ps (default: 5.0)')
    parser.add_argument('--cavity-tau', type=float, default=5.0, help='Cavity thermostat time constant in ps (default: 5.0)')
     
    # Cavity coupling option
    parser.add_argument('--no-cavity', action='store_true', help='Disable cavity coupling (default: False, cavity enabled)')
    
    # Timestepping options
    parser.add_argument('--fixed-timestep', action='store_true', help='Use fixed timestep instead of adaptive timestepping (default: False, adaptive enabled)')
    parser.add_argument('--timestep', type=float, default=1.0, help='Fixed timestep in femtoseconds (default: 1.0 fs, only used with --fixed-timestep)')
    
    # Energy tracking option
    parser.add_argument('--no-energy-tracking', action='store_true', help='Disable energy tracking to improve performance (default: False, energy tracking enabled)')
    
    # Replica specification options (for local loop mode)
    replica_group = parser.add_mutually_exclusive_group()
    replica_group.add_argument('--replicas', type=str, help='Replicas to run (e.g., "1-10" or "1,3,5,7")')
    replica_group.add_argument('--start-replica', type=int, help='Starting replica number (use with --end-replica)')
    parser.add_argument('--end-replica', type=int, help='Ending replica number (use with --start-replica)')
    
    # Logging options
    parser.add_argument('--log-to-file', action='store_true', help='Log output to files (default: False)')
    parser.add_argument('--log-to-console', action='store_true', help='Log output to console (default: False)')
    
    # F(k,t) calculation options
    parser.add_argument('--enable-fkt', action='store_true', help='Enable F(k,t) density correlation calculation (default: False)')
    parser.add_argument('--fkt-kmag', type=float, default=1.0, help='Wavevector magnitude for F(k,t) (default: 1.0)')
    parser.add_argument('--fkt-wavevectors', type=int, default=50, help='Number of wavevectors for F(k,t) (default: 50)')
    parser.add_argument('--fkt-ref-interval', type=float, default=1.0, help='Reference interval for F(k,t) in ps (default: 1.0)')
    parser.add_argument('--fkt-max-refs', type=int, default=10, help='Max reference frames for F(k,t) (default: 10)')
    parser.add_argument('--max-energy-output-time', type=float, default=None, help='Maximum time in ps to output energy data (default: no limit)')
    
    # Device options
    parser.add_argument('--device', type=str, default='CPU', choices=['CPU', 'GPU'], help='Device to use (default: CPU)')
    parser.add_argument('--gpu-id', type=int, default=0, help='GPU ID to use when device=GPU (default: 0)')
    
    args = parser.parse_args()
    
    # Handle cavity coupling logic (default: cavity enabled)
    incavity = not args.no_cavity
    
    # Validate replica arguments
    if args.start_replica is not None and args.end_replica is None:
        print("Error: --end-replica is required when using --start-replica")
        return 1
    if args.end_replica is not None and args.start_replica is None:
        print("Error: --start-replica is required when using --start-replica")
        return 1
    
    # Determine replicas to run
    replica_list = []
    mode = "SLURM"
    
    if args.replicas is not None:
        # Local loop mode with replicas string
        replica_list = parse_replicas(args.replicas)
        mode = "LOCAL_LOOP"
    elif args.start_replica is not None and args.end_replica is not None:
        # Local loop mode with start/end range
        replica_list = list(range(args.start_replica, args.end_replica + 1))
        mode = "LOCAL_LOOP"
    else:
        # SLURM mode - get single replica from environment
        replica, frame, job_id = get_slurm_info()
        if replica is None:
            print("Warning: SLURM_ARRAY_TASK_ID not found and no replica range specified.")
            print("Using default replica=0, frame=0")
            print("Use --replicas or --start-replica/--end-replica for local loop mode.")
            replica_list = [0]
            mode = "FALLBACK"
        else:
            replica_list = [replica]
            mode = "SLURM"
    
    # Find the selected experiment
    exp_dict = {name: (name, mol, cav, fq) for name, mol, cav, fq in BUSSI_LANGEVIN_EXPERIMENTS}
    if args.experiment not in exp_dict:
        print(f"Error: Unknown experiment '{args.experiment}'")
        print(f"Available experiments: {list(exp_dict.keys())}")
        return 1
    
    exp_name, molecular_thermo, cavity_thermo, finite_q = exp_dict[args.experiment]
    
    print("="*60)
    print(f"{'CAVITY' if incavity else 'STANDARD'} MD EXPERIMENT ({mode} MODE)")
    print("="*60)
    print(f"Experiment: {exp_name}")
    print(f"Cavity coupling: {'Enabled' if incavity else 'Disabled'}")
    if incavity:
        print(f"Coupling strength: {args.coupling}")
    print(f"Replicas to run: {replica_list}")
    print(f"Device: {args.device}" + (f" (GPU {args.gpu_id})" if args.device == 'GPU' else ""))
    print(f"Runtime: {args.runtime} ps")
    print(f"Timestepping: {'Fixed ' + str(args.timestep) + ' fs' if args.fixed_timestep else 'Adaptive'}")
    print(f"Molecular thermostat: {molecular_thermo} (tau={args.molecular_tau} ps)")
    if incavity:
        print(f"Cavity thermostat: {cavity_thermo} (tau={args.cavity_tau} ps)")
        print(f"Finite Q: {finite_q}")
    print(f"Logging: {'file ' if args.log_to_file else ''}{'console' if args.log_to_console else ''}{'none' if not args.log_to_file and not args.log_to_console else ''}".strip())
    print(f"F(k,t) calculation: {'Enabled' if args.enable_fkt else 'Disabled'}")
    if args.enable_fkt:
        print(f"  k magnitude: {args.fkt_kmag}")
        print(f"  Wavevectors: {args.fkt_wavevectors}")
        print(f"  Reference interval: {args.fkt_ref_interval:.1f} ps")
        print(f"  Max references: {args.fkt_max_refs}")
    print(f"Energy tracking: {'Disabled' if args.no_energy_tracking else 'Enabled'}")
    print("="*60)
    
    # Run experiments for all replicas
    total_experiments = len(replica_list)
    successful_experiments = 0
    failed_experiments = 0
    
    for i, replica in enumerate(replica_list, 1):
        frame = replica  # Frame number equals replica number as requested
        
        print(f"\n[{i}/{total_experiments}] Running replica {replica} (frame {frame})...")
        
        success = run_single_experiment(
            exp_name, molecular_thermo, cavity_thermo, finite_q, 
            args.coupling, replica, frame, args.runtime, args.molecular_tau, args.cavity_tau, 
            args.log_to_file, args.log_to_console, args.enable_fkt, 
            args.fkt_kmag, args.fkt_wavevectors, args.fkt_ref_interval, args.fkt_max_refs,
            args.max_energy_output_time, args.device, args.gpu_id, incavity, args.fixed_timestep, 
            args.timestep, not args.no_energy_tracking
        )
        
        if success:
            successful_experiments += 1
        else:
            failed_experiments += 1
    
    # Final summary
    print("\n" + "="*60)
    print("EXPERIMENT SUMMARY")
    print("="*60)
    print(f"Total experiments: {total_experiments}")
    print(f"Successful: {successful_experiments}")
    print(f"Failed: {failed_experiments}")
    print(f"Success rate: {100.0 * successful_experiments / total_experiments:.1f}%")
    
    if failed_experiments == 0:
        print("\n🎉 All experiments completed successfully!")
        return 0
    else:
        print(f"\n⚠️  {failed_experiments} experiment(s) failed!")
        return 1


if __name__ == '__main__':
    sys.exit(run_cavity_experiments()) 