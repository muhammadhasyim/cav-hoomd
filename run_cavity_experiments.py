#!/usr/bin/env python3
"""
Simple Cavity MD Experiment Runner for SLURM Array Jobs

Usage:
    python run_cavity_experiments.py --coupling COUPLSTR [--runtime PS] [--experiment EXPERIMENT] [--molecular-tau TAU] [--cavity-tau TAU] [--log-to-file] [--log-to-console] [--enable-fkt] [--fkt-kmag K]

Examples:
    python run_cavity_experiments.py --coupling 1e-3
    python run_cavity_experiments.py --coupling 5e-4 --experiment bussi_langevin_finiteq
    python run_cavity_experiments.py --coupling 1e-2 --molecular-tau 10.0 --cavity-tau 1.0
    python run_cavity_experiments.py --coupling 2e-3 --runtime 200 --molecular-tau 2.0 --cavity-tau 0.5
    
    # Logging examples:
    python run_cavity_experiments.py --coupling 1e-3 --log-to-console  # Log to console only
    python run_cavity_experiments.py --coupling 1e-3 --log-to-file  # Log to file only
    python run_cavity_experiments.py --coupling 1e-3 --log-to-file --log-to-console  # Log to both
    
    # F(k,t) calculation examples:
    python run_cavity_experiments.py --coupling 1e-3 --enable-fkt  # Enable F(k,t) calculation
    python run_cavity_experiments.py --coupling 1e-3 --enable-fkt --fkt-kmag 2.0 --fkt-wavevectors 100  # Custom F(k,t) settings

SLURM Usage:
    # Submit array job with different coupling values
    sbatch --array=1-100 run_job.sh
    
    # In run_job.sh:
    python run_cavity_experiments.py --coupling 1e-3
    
    # Both replica and frame will be automatically set to $SLURM_ARRAY_TASK_ID

Note on tau values:
    - For Langevin thermostats, tau must be > 0 (since gamma = 1/tau)
    - For overdamped dynamics (tau → 0), use molecular_thermostat='brownian' or cavity_thermostat='brownian' instead
    - Brownian dynamics is automatically used when use_brownian_overdamped=True and cavity_damping_factor > 5.0
    - Brownian dynamics explicitly models the overdamped limit where inertia is neglected
"""

import sys
import os
import argparse
from pathlib import Path
from cavitymd import CavityMDSimulation

# Available bussi_langevin experiments: (name, molecular_thermostat, cavity_thermostat, finite_q)
BUSSI_LANGEVIN_EXPERIMENTS = [
    ("bussi_langevin_finiteq", "bussi", "langevin", True),
    ("bussi_langevin_no_finiteq", "bussi", "langevin", False),
]

def get_slurm_info():
    """Get SLURM job information from environment variables."""
    task_id = os.environ.get('SLURM_ARRAY_TASK_ID')
    job_id = os.environ.get('SLURM_JOB_ID', 'unknown')
    
    if task_id is None:
        print("Warning: SLURM_ARRAY_TASK_ID not found. Using default replica=0, frame=1")
        print("This script is designed for SLURM array jobs.")
        replica = 0
        frame = -1
    else:
        replica = int(task_id)
        frame = int(task_id)
        print(f"SLURM Array Job ID: {job_id}")
        print(f"SLURM Array Task ID: {task_id} (using as replica and frame number)")
    
    return replica, frame

def run_single_experiment(exp_name, molecular_thermo, cavity_thermo, finite_q, coupling, replica, frame, runtime_ps, molecular_tau, cavity_tau, log_to_file, log_to_console, enable_fkt, fkt_kmag, fkt_wavevectors, fkt_ref_interval, fkt_max_refs, max_energy_output_time=None, device='CPU', gpu_id=0):
    """Run a single experiment."""
    
    try:
        # Create experiment directory with coupling info
        coupling_str = f"{coupling:.0e}".replace("-", "neg").replace("+", "pos")
        exp_dir = Path(f"{exp_name}_coupling_{coupling_str}")
        exp_dir.mkdir(exist_ok=True)
        
        print(f"Running experiment: {exp_name}")
        print(f"Coupling strength: {coupling}")
        print(f"Replica: {replica}")
        print(f"Frame: {frame}")
        print(f"Output directory: {exp_dir}")
        
        # Run simulation
        sim = CavityMDSimulation(
            job_dir=str(exp_dir),
            replica=replica,  # Use SLURM task ID as replica
            freq=1560.0,
            couplstr=coupling,  # Use user-specified coupling
            incavity=True,
            runtime_ps=runtime_ps,
            input_gsd='../init-0.gsd',
            frame=frame,  # Use SLURM task ID as frame
            name='prod',
            error_tolerance=1.0,
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
            device=device,
            gpu_id=gpu_id
        )
        
        exit_code = sim.run()
        
        if exit_code == 0:
            print(f"✅ {exp_name} completed successfully (replica {replica}, frame {frame})")
            return True
        else:
            print(f"❌ {exp_name} failed with exit code {exit_code} (replica {replica}, frame {frame})")
            return False
            
    except Exception as e:
        print(f"❌ {exp_name} - error: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Run cavity MD experiments with SLURM")
    parser.add_argument('--coupling', type=float, required=True, help='Coupling strength (required)')
    parser.add_argument('--runtime', type=float, default=500.0, help='Runtime in ps (default: 500)')
    parser.add_argument('--experiment', type=str, default='bussi_langevin_finiteq', 
                       choices=[name for name, _, _, _ in BUSSI_LANGEVIN_EXPERIMENTS],
                       help='Specific experiment to run (default: bussi_langevin_finiteq)')
    parser.add_argument('--molecular-tau', type=float, default=5.0, help='Molecular thermostat time constant in ps (default: 5.0)')
    parser.add_argument('--cavity-tau', type=float, default=5.0, help='Cavity thermostat time constant in ps (default: 5.0)')
     
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
    
    # Get SLURM task information
    replica, frame = get_slurm_info()
    
    # Find the selected experiment
    exp_dict = {name: (name, mol, cav, fq) for name, mol, cav, fq in BUSSI_LANGEVIN_EXPERIMENTS}
    if args.experiment not in exp_dict:
        print(f"Error: Unknown experiment '{args.experiment}'")
        print(f"Available experiments: {list(exp_dict.keys())}")
        return 1
    
    exp_name, molecular_thermo, cavity_thermo, finite_q = exp_dict[args.experiment]
    
    print("="*60)
    print("CAVITY MD EXPERIMENT (SLURM MODE)")
    print("="*60)
    print(f"Experiment: {exp_name}")
    print(f"Coupling strength: {args.coupling}")
    print(f"Replica: {replica}")
    print(f"Frame: {frame}")
    print(f"Device: {args.device}" + (f" (GPU {args.gpu_id})" if args.device == 'GPU' else ""))
    print(f"Runtime: {args.runtime} ps")
    print(f"Molecular thermostat: {molecular_thermo} (tau={args.molecular_tau} ps)")
    print(f"Cavity thermostat: {cavity_thermo} (tau={args.cavity_tau} ps)")
    print(f"Finite Q: {finite_q}")
    print(f"Logging: {'file ' if args.log_to_file else ''}{'console' if args.log_to_console else ''}{'none' if not args.log_to_file and not args.log_to_console else ''}".strip())
    print(f"F(k,t) calculation: {'Enabled' if args.enable_fkt else 'Disabled'}")
    if args.enable_fkt:
        print(f"  k magnitude: {args.fkt_kmag}")
        print(f"  Wavevectors: {args.fkt_wavevectors}")
        print(f"  Reference interval: {args.fkt_ref_interval:.1f} ps")
        print(f"  Max references: {args.fkt_max_refs}")
    print("="*60)
    
    # Run the experiment
    success = run_single_experiment(
        exp_name, molecular_thermo, cavity_thermo, finite_q, 
        args.coupling, replica, frame, args.runtime, args.molecular_tau, args.cavity_tau, 
        args.log_to_file, args.log_to_console, args.enable_fkt, 
        args.fkt_kmag, args.fkt_wavevectors, args.fkt_ref_interval, args.fkt_max_refs,
        args.max_energy_output_time, args.device, args.gpu_id
    )
    
    if success:
        print("\n🎉 Experiment completed successfully!")
        return 0
    else:
        print("\n💥 Experiment failed!")
        return 1

if __name__ == '__main__':
    sys.exit(main()) 
