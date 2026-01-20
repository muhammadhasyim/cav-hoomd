#!/usr/bin/env python3
"""
Finite-size scaling study with constant collective coupling on GPU.

This script runs a series of cavity MD simulations where the single-molecule
coupling lambda is scaled as lambda = g / sqrt(N) so that the collective
coupling g is constant across system size.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from hoomd.cavitymd.analysis.velocity_validation import compute_single_molecule_coupling
from hoomd.cavitymd.studies import compute_replica_count

# Use the existing lattice generator in examples
from initlattice_equilibrium import main as generate_equilibrium_gsd


def build_parser() -> argparse.ArgumentParser:
    """Build command-line parser for the finite-size scaling study."""
    parser = argparse.ArgumentParser(
        description="Finite-size scaling study with constant collective coupling."
    )
    parser.add_argument(
        "--collective-coupling",
        type=float,
        default=0.02,
        help="Target collective coupling g in atomic units (default: 0.02).",
    )
    parser.add_argument(
        "--n-values",
        type=str,
        default="100,250,500,1000,2500,5000,10000",
        help="Comma-separated list of molecule counts N (default: logarithmic set).",
    )
    parser.add_argument(
        "--base-molecules",
        type=int,
        default=250,
        help="Reference system size for replica scaling (default: 250).",
    )
    parser.add_argument(
        "--base-replicas",
        type=int,
        default=500,
        help="Reference replicas for base system size (default: 500).",
    )
    parser.add_argument(
        "--min-replicas",
        type=int,
        default=1,
        help="Minimum replicas per N (default: 1).",
    )
    parser.add_argument(
        "--max-replicas",
        type=int,
        default=None,
        help="Optional maximum replicas per N.",
    )
    parser.add_argument(
        "--replica-start",
        type=int,
        default=0,
        help="Starting replica index (default: 0).",
    )
    parser.add_argument(
        "--frequency",
        type=float,
        default=2000.0,
        help="Cavity frequency in cm^-1 (default: 2000).",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=100.0,
        help="Temperature in K (default: 100).",
    )
    parser.add_argument(
        "--runtime-ps",
        type=float,
        default=1000.0,
        help="Simulation runtime in ps (default: 1000).",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="finite_size_scaling",
        help="Output directory for the study.",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="GPU",
        choices=["CPU", "GPU"],
        help="Compute device (default: GPU).",
    )
    parser.add_argument(
        "--gpu-id",
        type=int,
        default=0,
        help="GPU ID to use when device is GPU (default: 0).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility (default: replica-based).",
    )
    parser.add_argument(
        "--molecular-thermostat",
        type=str,
        default="bussi",
        choices=["bussi", "langevin", "none"],
        help="Thermostat for molecular particles (default: bussi).",
    )
    parser.add_argument(
        "--cavity-thermostat",
        type=str,
        default="langevin",
        choices=["bussi", "langevin", "none"],
        help="Thermostat for cavity particle (default: langevin).",
    )
    parser.add_argument(
        "--molecular-thermostat-tau",
        type=float,
        default=5.0,
        help="Molecular thermostat time constant in ps (default: 5.0).",
    )
    parser.add_argument(
        "--cavity-thermostat-tau",
        type=float,
        default=5.0,
        help="Cavity thermostat time constant in ps (default: 5.0).",
    )
    parser.add_argument(
        "--gsd-output-period-ps",
        type=float,
        default=1.0,
        help="GSD output period in ps (default: 1.0).",
    )
    parser.add_argument(
        "--energy-output-period-ps",
        type=float,
        default=0.1,
        help="Energy output period in ps (default: 0.1).",
    )
    parser.add_argument(
        "--enable-hdf5-output",
        action="store_true",
        default=True,
        help="Enable HDF5 observable output (default: enabled).",
    )
    parser.add_argument(
        "--enable-dipole-autocorr",
        action="store_true",
        default=True,
        help="Enable dipole autocorrelation tracking (default: enabled).",
    )
    parser.add_argument(
        "--switch-time-ps",
        type=float,
        default=None,
        help="Time in ps when coupling turns on (default: on from start). Use for equilibration without coupling.",
    )
    parser.add_argument(
        "--damping-ratio",
        type=float,
        default=0.0,
        help="Damping ratio (zeta) for the cavity mode (default: 0.0).",
    )
    parser.add_argument(
        "--error-tolerance",
        type=float,
        default=5.0,
        help="Target error tolerance for adaptive timestep (default: 5.0).",
    )
    parser.add_argument(
        "--initial-fraction",
        type=float,
        default=1e-5,
        help="Initial fraction for shock dampening (default: 1e-5).",
    )
    parser.add_argument(
        "--time-constant-ps",
        type=float,
        default=50.0,
        help="Time constant for error tolerance ramping in ps (default: 50.0).",
    )
    # F(k,t) correlation tracking
    parser.add_argument(
        "--enable-fkt",
        action="store_true",
        default=False,
        help="Enable F(k,t) density correlation tracking (default: disabled).",
    )
    parser.add_argument(
        "--fkt-output-period-ps",
        type=float,
        default=0.1,
        help="F(k,t) output period in ps (default: 0.1).",
    )
    parser.add_argument(
        "--fkt-ref-interval-ps",
        type=float,
        default=200.0,
        help="F(k,t) reference interval in ps (default: 200.0).",
    )
    parser.add_argument(
        "--fkt-kmag",
        type=float,
        default=1.0,
        help="F(k,t) k-magnitude for density correlations (default: 1.0).",
    )
    parser.add_argument(
        "--fkt-num-wavevectors",
        type=int,
        default=50,
        help="Number of wavevectors for F(k,t) calculation (default: 50).",
    )
    # GSD output control
    parser.add_argument(
        "--disable-gsd",
        action="store_true",
        default=False,
        help="Disable GSD trajectory output (default: enabled).",
    )
    # SLURM single-replica mode arguments
    parser.add_argument(
        "--single-replica-mode",
        action="store_true",
        help="Run single replica only (for SLURM array tasks).",
    )
    parser.add_argument(
        "--replica-id",
        type=int,
        default=None,
        help="Specific replica ID to run (for SLURM tasks).",
    )
    parser.add_argument(
        "--n-molecules-single",
        type=int,
        default=None,
        help="Single system size to run (for SLURM tasks).",
    )
    return parser


def run_single_replica(args, n_molecules: int, replica: int) -> int:
    """Run a single replica for SLURM array tasks."""
    from hoomd.cavitymd.simulation import CavityMDSimulation
    from hoomd.cavitymd.utils import PhysicalConstants
    
    lambda_single = compute_single_molecule_coupling(
        args.collective_coupling, n_molecules
    )
    collective_check = lambda_single * (n_molecules ** 0.5)
    
    output_dir = Path(args.output_dir).resolve()
    job_dir = (output_dir / f"N_{n_molecules}" / f"replica_{replica}").resolve()
    job_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print(f"Single Replica Mode")
    print("=" * 60)
    print(f"  N_molecules = {n_molecules}")
    print(f"  replica = {replica}")
    print(f"  lambda_single = {lambda_single:.6e} a.u.")
    print(f"  g = {collective_check:.6e} a.u.")
    print(f"  output = {job_dir}")
    if args.switch_time_ps is not None:
        print(f"  Protocol: Equilibrate {args.switch_time_ps} ps without coupling, then {args.runtime_ps - args.switch_time_ps} ps with coupling")
    else:
        print(f"  Protocol: Coupling on from start for {args.runtime_ps} ps")
    print("=" * 60)
    
    # Generate equilibrium GSD
    generate_equilibrium_gsd(
        job_dir=str(job_dir),
        replica=replica,
        incavity=True,
        couplstr=lambda_single,
        frequency=args.frequency,
        Nmol=n_molecules,
        temperature_K=args.temperature,
        seed=args.seed,
    )
    
    input_gsd = job_dir / f"molecular-{replica}.gsd"
    
    # Calculate dissipation from damping_ratio
    phmass = 1.0  # Photon mass is 1.0 in a.u.
    omegac = args.frequency / PhysicalConstants.HARTREE_TO_CM_MINUS1
    dissipation = 2 * args.damping_ratio * phmass * omegac
    
    # Determine GSD output period (use very large value to effectively disable)
    gsd_period = 1e9 if args.disable_gsd else args.gsd_output_period_ps
    
    simulation = CavityMDSimulation(
        job_dir=str(job_dir),
        replica=replica,
        freq=args.frequency,
        incavity=True,
        couplstr=lambda_single,
        runtime_ps=args.runtime_ps,
        input_gsd=str(input_gsd),
        temperature=args.temperature,
        device=args.device,
        gpu_id=args.gpu_id,
        # F(k,t) configuration
        enable_fkt=args.enable_fkt,
        fkt_output_period_ps=args.fkt_output_period_ps,
        fkt_reference_interval_ps=args.fkt_ref_interval_ps,
        fkt_kmag=args.fkt_kmag,
        fkt_num_wavevectors=args.fkt_num_wavevectors,
        # Thermostat configuration
        molecular_thermostat=args.molecular_thermostat,
        cavity_thermostat=args.cavity_thermostat,
        molecular_thermostat_tau=args.molecular_thermostat_tau,
        cavity_thermostat_tau=args.cavity_thermostat_tau,
        # Output configuration
        enable_energy_tracking=True,
        enable_temp_tracker=True,
        enable_hdf5_output=args.enable_hdf5_output,
        enable_dipole_autocorr=args.enable_dipole_autocorr,
        energy_output_period_ps=args.energy_output_period_ps,
        gsd_output_period_ps=gsd_period,
        seed=args.seed,
        # Two-phase protocol
        switch_time_ps=args.switch_time_ps,
        dissipation=dissipation,
        # Adaptive timestepping
        error_tolerance=args.error_tolerance,
        initial_fraction=args.initial_fraction,
        time_constant_ps=args.time_constant_ps,
    )
    
    try:
        exit_code = simulation.run()
    except Exception as e:
        print(f"ERROR: simulation.run() raised exception: {e}", flush=True)
        exit_code = 1
    
    status = "success" if exit_code == 0 else "failure"
    
    # Write individual replica summary CSV
    summary_path = job_dir / "finite_size_scaling_summary.csv"
    with summary_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "N_molecules",
            "lambda_single",
            "collective_coupling",
            "replica",
            "replica_count",
            "device",
            "runtime_ps",
            "status",
        ])
        writer.writerow([
            n_molecules,
            f"{lambda_single:.8e}",
            f"{collective_check:.8e}",
            replica,
            1,  # replica_count not meaningful for single replica
            args.device,
            args.runtime_ps,
            status,
        ])
    
    print(f"\nReplica summary written to: {summary_path}")
    print(f"Status: {status}")
    
    return exit_code


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    # Handle single-replica mode (for SLURM array tasks)
    if args.single_replica_mode:
        if args.n_molecules_single is None or args.replica_id is None:
            print("ERROR: --single-replica-mode requires both --n-molecules-single and --replica-id")
            return 1
        return run_single_replica(args, args.n_molecules_single, args.replica_id)

    # Standard mode: run all replicas for all system sizes
    n_values = [int(value.strip()) for value in args.n_values.split(",") if value.strip()]
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_path = output_dir / "finite_size_scaling_summary.csv"
    with summary_path.open("w", newline="") as summary_file:
        writer = csv.writer(summary_file)
        writer.writerow(
            [
                "N_molecules",
                "lambda_single",
                "collective_coupling",
                "replica",
                "replica_count",
                "device",
                "runtime_ps",
                "status",
            ]
        )

        for n_molecules in n_values:
            lambda_single = compute_single_molecule_coupling(
                args.collective_coupling, n_molecules
            )
            collective_check = lambda_single * (n_molecules ** 0.5)
            replica_count = compute_replica_count(
                base_molecules=args.base_molecules,
                base_replicas=args.base_replicas,
                target_molecules=n_molecules,
                min_replicas=args.min_replicas,
                max_replicas=args.max_replicas,
            )

            print("=" * 60)
            print(f"Running N={n_molecules}")
            print(f"  lambda_single = {lambda_single:.6e} a.u.")
            print(f"  g = {collective_check:.6e} a.u.")
            print(f"  replicas = {replica_count}")
            if args.switch_time_ps is not None:
                print(f"  Protocol: Equilibrate {args.switch_time_ps} ps without coupling, then {args.runtime_ps - args.switch_time_ps} ps with coupling")
            else:
                print(f"  Protocol: Coupling on from start for {args.runtime_ps} ps")

            for replica_offset in range(replica_count):
                replica = args.replica_start + replica_offset
                job_dir = (output_dir / f"N_{n_molecules}" / f"replica_{replica}").resolve()
                job_dir.mkdir(parents=True, exist_ok=True)
                print(f"  replica {replica} -> {job_dir}")

                generate_equilibrium_gsd(
                    job_dir=str(job_dir),
                    replica=replica,
                    incavity=True,
                    coupling=lambda_single,
                    frequency=args.frequency,
                    Nmol=n_molecules,
                    temperature_K=args.temperature,
                    seed=args.seed,
                )

                input_gsd = job_dir / f"molecular-{replica}.gsd"

                from hoomd.cavitymd.simulation import CavityMDSimulation

                print(f"DEBUG: About to create CavityMDSimulation", flush=True)
                
                # Calculate dissipation from damping_ratio (matching 05_advanced_run.py)
                from hoomd.cavitymd.utils import PhysicalConstants
                phmass = 1.0  # Photon mass is 1.0 in a.u.
                omegac = args.frequency / PhysicalConstants.HARTREE_TO_CM_MINUS1
                dissipation = 2 * args.damping_ratio * phmass * omegac
                
                simulation = CavityMDSimulation(
                    job_dir=str(job_dir),
                    replica=replica,
                    freq=args.frequency,
                    incavity=True,
                    couplstr=lambda_single,
                    runtime_ps=args.runtime_ps,
                    input_gsd=str(input_gsd),
                    temperature=args.temperature,
                    device=args.device,
                    gpu_id=args.gpu_id,
                    enable_fkt=False,
                    molecular_thermostat=args.molecular_thermostat,
                    cavity_thermostat=args.cavity_thermostat,
                    molecular_thermostat_tau=args.molecular_thermostat_tau,
                    cavity_thermostat_tau=args.cavity_thermostat_tau,
                    enable_energy_tracking=True,
                    enable_temp_tracker=True,
                    enable_hdf5_output=args.enable_hdf5_output,
                    enable_dipole_autocorr=args.enable_dipole_autocorr,
                    energy_output_period_ps=args.energy_output_period_ps,
                    gsd_output_period_ps=args.gsd_output_period_ps,
                    seed=args.seed,
                    # Time-varying coupling parameters
                    switch_time_ps=args.switch_time_ps,
                    dissipation=dissipation,
                    # Adaptive timestepping parameters
                    error_tolerance=args.error_tolerance,
                    initial_fraction=args.initial_fraction,
                    time_constant_ps=args.time_constant_ps,
                )

                print(f"DEBUG: About to call simulation.run()", flush=True)
                try:
                    exit_code = simulation.run()
                    print(f"DEBUG: simulation.run() returned: {exit_code}", flush=True)
                except Exception as e:
                    print(f"DEBUG: simulation.run() raised exception: {e}", flush=True)
                    exit_code = 1
                
                status = "success" if exit_code == 0 else "failure"
                print(f"DEBUG: Writing to CSV - status={status}", flush=True)
                writer.writerow(
                    [
                        n_molecules,
                        f"{lambda_single:.8e}",
                        f"{collective_check:.8e}",
                        replica,
                        replica_count,
                        args.device,
                        args.runtime_ps,
                        status,
                    ]
                )
                summary_file.flush()

    print(f"\nSummary written to: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
