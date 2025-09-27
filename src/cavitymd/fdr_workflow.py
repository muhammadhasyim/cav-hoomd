# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""FDR workflow management for fork-and-clone susceptibility measurements."""

from typing import Optional, List, Dict, Any, Tuple, Union
import hoomd
import numpy as np
import logging
import time
from pathlib import Path
import threading
import queue
from dataclasses import dataclass

from .state_manager import StateSnapshotManager
from .fdr_forces import PerturbationForce
from .analysis import ElapsedTimeTracker, FieldAutocorrelationTracker
from .utils import PhysicalConstants


@dataclass
class FDRMeasurement:
    """Container for FDR measurement data."""
    waiting_time_ps: float
    kvector: np.ndarray
    k_magnitude: float
    times_ps: np.ndarray
    chi_k: np.ndarray  # Susceptibility χ_k(t,t_w)
    C_k: np.ndarray    # Correlation function C_k(t,t_w) 
    rho_k_plus: np.ndarray   # Density correlation for +h₀ clone
    rho_k_minus: np.ndarray  # Density correlation for -h₀ clone
    metadata: Dict[str, Any]


class FDRSusceptibilityTracker:
    """Tracker for measuring integrated response χ_k(t,t_w) in FDR experiments.
    
    This class handles the measurement of density correlations ρ_k(t) in 
    perturbed simulations and calculates the susceptibility:
    
    χ_k(t,t_w) = (ρ_k_plus(t) - ρ_k_minus(t)) / (2 * N * h₀)
    
    where ρ_k_plus/minus are measured in clones with ±h₀ perturbations.
    """
    
    def __init__(self,
                 kvector: np.ndarray,
                 perturbation_amplitude: float,
                 output_times_ps: np.ndarray,
                 output_file: str,
                 logger: Optional[logging.Logger] = None):
        """Initialize susceptibility tracker.
        
        Parameters
        ----------
        kvector : np.ndarray
            k-vector for density measurements
        perturbation_amplitude : float
            Perturbation amplitude h₀
        output_times_ps : np.ndarray
            Times at which to output results
        output_file : str
            Output file path
        logger : logging.Logger, optional
            Logger for status messages
        """
        self.kvector = np.array(kvector)
        self.k_magnitude = np.linalg.norm(self.kvector)
        self.perturbation_amplitude = perturbation_amplitude
        self.output_times_ps = output_times_ps
        self.output_file = output_file
        self.logger = logger or logging.getLogger(__name__)
        
        # Storage for density correlations
        self.rho_k_data_plus: List[complex] = []
        self.rho_k_data_minus: List[complex] = []
        self.measurement_times_ps: List[float] = []
        
        # Metadata
        self.waiting_time_ps: Optional[float] = None
        self.start_time = time.time()
        
    def measure_density_correlation(self, 
                                  simulation: hoomd.Simulation,
                                  time_tracker: ElapsedTimeTracker,
                                  sign: float) -> complex:
        """Measure density correlation ρ_k(t) = Σ_j exp(i k·r_j).
        
        Parameters
        ----------
        simulation : hoomd.Simulation
            HOOMD simulation object
        time_tracker : ElapsedTimeTracker
            Time tracking object
        sign : float
            Sign of perturbation (+1 or -1)
            
        Returns
        -------
        complex
            Density correlation ρ_k(t)
        """
        current_time_ps = time_tracker.elapsed_time()
        
        with simulation.state.cpu_local_snapshot as snap:
            positions = np.array(snap.particles.position)
            
            # Calculate ρ_k = Σ_j exp(i k·r_j)
            k_dot_r = np.dot(positions, self.kvector)
            rho_k = np.sum(np.exp(1j * k_dot_r))
            
        # Store data based on sign
        if sign > 0:
            self.rho_k_data_plus.append(rho_k)
        else:
            self.rho_k_data_minus.append(rho_k)
            
        self.measurement_times_ps.append(current_time_ps)
        
        return rho_k
    
    def calculate_susceptibility(self) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate susceptibility χ_k(t,t_w) from collected data.
        
        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            Times and susceptibility values
        """
        if len(self.rho_k_data_plus) == 0 or len(self.rho_k_data_minus) == 0:
            return np.array([]), np.array([])
        
        # Ensure both arrays have same length
        min_length = min(len(self.rho_k_data_plus), len(self.rho_k_data_minus))
        
        rho_plus = np.array(self.rho_k_data_plus[:min_length])
        rho_minus = np.array(self.rho_k_data_minus[:min_length])
        times = np.array(self.measurement_times_ps[:min_length])
        
        # Calculate susceptibility: χ_k = (ρ_k⁺ - ρ_k⁻) / (2 * N * h₀)
        # Note: We take the real part for the susceptibility
        N_particles = len(rho_plus)  # This should be number of particles, will be corrected
        chi_k = np.real(rho_plus - rho_minus) / (2.0 * N_particles * self.perturbation_amplitude)
        
        return times, chi_k
    
    def save_results(self, 
                    times_ps: np.ndarray,
                    chi_k: np.ndarray,
                    C_k: Optional[np.ndarray] = None) -> None:
        """Save susceptibility results to file.
        
        Parameters
        ----------
        times_ps : np.ndarray
            Time points in picoseconds
        chi_k : np.ndarray
            Susceptibility values
        C_k : np.ndarray, optional
            Correlation function values (from unperturbed simulation)
        """
        header = f"# FDR Susceptibility Data\n"
        header += f"# k-vector: [{self.kvector[0]:.6f}, {self.kvector[1]:.6f}, {self.kvector[2]:.6f}]\n"
        header += f"# |k| = {self.k_magnitude:.6f}\n"
        header += f"# h₀ = {self.perturbation_amplitude:.6e}\n"
        header += f"# t_w = {self.waiting_time_ps:.3f} ps\n"
        header += f"# Columns: time_ps, chi_k"
        if C_k is not None:
            header += ", C_k"
        header += "\n"
        
        # Prepare data
        if C_k is not None and len(C_k) == len(times_ps):
            data = np.column_stack([times_ps, chi_k, C_k])
        else:
            data = np.column_stack([times_ps, chi_k])
        
        # Save to file
        np.savetxt(self.output_file, data, header=header, fmt='%.6e')
        self.logger.info(f"Susceptibility data saved to {self.output_file}")


class FDRWorkflowManager:
    """Manager for FDR fork-and-clone workflow.
    
    This class orchestrates the complete FDR measurement process:
    1. Run unperturbed simulation with periodic snapshots at waiting times
    2. At each waiting time, create two clones with ±h₀ perturbations
    3. Evolve clones and measure density correlations ρ_k(t)
    4. Calculate susceptibility χ_k(t,t_w) and combine with unperturbed C_k(t,t_w)
    5. Output results for FDR analysis
    """
    
    def __init__(self,
                 waiting_times_ps: List[float],
                 kvector: np.ndarray,
                 perturbation_amplitude: float = 1e-6,
                 clone_runtime_ps: float = 100.0,
                 output_dir: str = "fdr_results",
                 enable_shell_averaging: bool = True,
                 num_shell_vectors: int = 50,
                 logger: Optional[logging.Logger] = None):
        """Initialize FDR workflow manager.
        
        Parameters
        ----------
        waiting_times_ps : List[float]
            Waiting times t_w at which to fork simulations
        kvector : np.ndarray
            k-vector for FDR measurements
        perturbation_amplitude : float, optional
            Perturbation amplitude h₀ (default: 1e-6)
        clone_runtime_ps : float, optional
            Runtime for each clone simulation (default: 100.0 ps)
        output_dir : str, optional
            Output directory for results (default: "fdr_results")
        enable_shell_averaging : bool, optional
            Whether to average over k-vectors of same magnitude (default: True)
        num_shell_vectors : int, optional
            Number of k-vectors for shell averaging (default: 50)
        logger : logging.Logger, optional
            Logger for status messages
        """
        self.waiting_times_ps = sorted(waiting_times_ps)
        self.kvector = np.array(kvector)
        self.k_magnitude = np.linalg.norm(self.kvector)
        self.perturbation_amplitude = perturbation_amplitude
        self.clone_runtime_ps = clone_runtime_ps
        self.output_dir = Path(output_dir)
        self.enable_shell_averaging = enable_shell_averaging
        self.num_shell_vectors = num_shell_vectors
        self.logger = logger or logging.getLogger(__name__)
        
        # Initialize components
        self.state_manager = StateSnapshotManager(logger=self.logger)
        self.measurements: List[FDRMeasurement] = []
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate shell k-vectors if needed
        if self.enable_shell_averaging:
            self.shell_kvectors = self._generate_shell_kvectors()
        else:
            self.shell_kvectors = [self.kvector]
        
        self.logger.info(f"FDR workflow initialized:")
        self.logger.info(f"  Waiting times: {self.waiting_times_ps}")
        self.logger.info(f"  k-vector: {self.kvector}")
        self.logger.info(f"  |k| = {self.k_magnitude:.6f}")
        self.logger.info(f"  h₀ = {self.perturbation_amplitude:.6e}")
        self.logger.info(f"  Shell averaging: {self.enable_shell_averaging}")
        if self.enable_shell_averaging:
            self.logger.info(f"  Shell vectors: {len(self.shell_kvectors)}")
        
    def run_fdr_measurement(self,
                           base_simulation: hoomd.Simulation,
                           base_time_tracker: ElapsedTimeTracker,
                           base_integrator_setup_func: callable,
                           waiting_time_ps: float) -> FDRMeasurement:
        """Run FDR measurement at a specific waiting time.
        
        Parameters
        ----------
        base_simulation : hoomd.Simulation
            Base simulation to fork from
        base_time_tracker : ElapsedTimeTracker
            Base time tracker
        base_integrator_setup_func : callable
            Function to set up integrator for clones
        waiting_time_ps : float
            Waiting time for this measurement
            
        Returns
        -------
        FDRMeasurement
            Complete FDR measurement data
        """
        self.logger.info(f"Starting FDR measurement at t_w = {waiting_time_ps:.3f} ps")
        
        # Create snapshot ID
        snapshot_id = f"tw_{waiting_time_ps:.3f}ps"
        
        # Capture current state
        snapshot = self.state_manager.capture_state(
            base_simulation, base_time_tracker, snapshot_id
        )
        
        # Create output times for this measurement
        output_times_ps = np.linspace(0, self.clone_runtime_ps, 100)  # Could be logarithmic
        
        # Run clone simulations for each k-vector in shell
        all_chi_k = []
        all_C_k = []
        
        for i, kvec in enumerate(self.shell_kvectors):
            self.logger.info(f"  Shell vector {i+1}/{len(self.shell_kvectors)}: {kvec}")
            
            # Run clone pair for this k-vector
            chi_k, C_k = self._run_clone_pair(
                base_simulation, base_time_tracker,
                base_integrator_setup_func, snapshot_id,
                kvec, output_times_ps, waiting_time_ps
            )
            
            all_chi_k.append(chi_k)
            all_C_k.append(C_k)
        
        # Average over shell if enabled
        if self.enable_shell_averaging and len(all_chi_k) > 1:
            chi_k_avg = np.mean(all_chi_k, axis=0)
            C_k_avg = np.mean(all_C_k, axis=0)
            self.logger.info(f"  Shell averaged over {len(all_chi_k)} k-vectors")
        else:
            chi_k_avg = all_chi_k[0]
            C_k_avg = all_C_k[0]
        
        # Create measurement object
        measurement = FDRMeasurement(
            waiting_time_ps=waiting_time_ps,
            kvector=self.kvector,
            k_magnitude=self.k_magnitude,
            times_ps=output_times_ps,
            chi_k=chi_k_avg,
            C_k=C_k_avg,
            rho_k_plus=np.array([]),  # These would be filled by detailed tracking
            rho_k_minus=np.array([]),
            metadata={
                'perturbation_amplitude': self.perturbation_amplitude,
                'clone_runtime_ps': self.clone_runtime_ps,
                'shell_averaging': self.enable_shell_averaging,
                'num_shell_vectors': len(self.shell_kvectors),
                'measurement_time': time.time()
            }
        )
        
        # Save individual measurement
        output_file = self.output_dir / f"fdr_tw_{waiting_time_ps:.3f}ps_k{self.k_magnitude:.3f}.dat"
        self._save_measurement(measurement, output_file)
        
        self.measurements.append(measurement)
        self.logger.info(f"FDR measurement completed at t_w = {waiting_time_ps:.3f} ps")
        
        return measurement
    
    def _run_clone_pair(self,
                       base_simulation: hoomd.Simulation,
                       base_time_tracker: ElapsedTimeTracker,
                       integrator_setup_func: callable,
                       snapshot_id: str,
                       kvector: np.ndarray,
                       output_times_ps: np.ndarray,
                       waiting_time_ps: float) -> Tuple[np.ndarray, np.ndarray]:
        """Run a pair of clones with ±h₀ perturbations.
        
        Parameters
        ----------
        base_simulation : hoomd.Simulation
            Base simulation
        base_time_tracker : ElapsedTimeTracker
            Base time tracker
        integrator_setup_func : callable
            Function to set up integrator
        snapshot_id : str
            Snapshot identifier
        kvector : np.ndarray
            k-vector for this measurement
        output_times_ps : np.ndarray
            Output time points
        waiting_time_ps : float
            Waiting time
            
        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            Susceptibility and correlation function
        """
        # Create susceptibility tracker
        output_file = self.output_dir / f"chi_k_tw_{waiting_time_ps:.3f}ps.dat"
        tracker = FDRSusceptibilityTracker(
            kvector=kvector,
            perturbation_amplitude=self.perturbation_amplitude,
            output_times_ps=output_times_ps,
            output_file=str(output_file),
            logger=self.logger
        )
        tracker.waiting_time_ps = waiting_time_ps
        
        # Run +h₀ clone
        self.logger.debug("  Running +h₀ clone")
        plus_clone, plus_tracker = self.state_manager.create_clone(
            base_simulation, base_time_tracker, snapshot_id
        )
        
        # Add +h₀ perturbation force
        plus_perturbation = PerturbationForce(
            kvector=kvector,
            amplitude=self.perturbation_amplitude,
            sign=+1.0
        )
        
        # Set up integrator for +h₀ clone
        integrator_setup_func(plus_clone, [plus_perturbation])
        
        # Run +h₀ clone simulation
        self._run_clone_simulation(plus_clone, plus_tracker, tracker, +1.0, output_times_ps)
        
        # Run -h₀ clone
        self.logger.debug("  Running -h₀ clone")
        minus_clone, minus_tracker = self.state_manager.create_clone(
            base_simulation, base_time_tracker, snapshot_id
        )
        
        # Add -h₀ perturbation force
        minus_perturbation = PerturbationForce(
            kvector=kvector,
            amplitude=self.perturbation_amplitude,
            sign=-1.0
        )
        
        # Set up integrator for -h₀ clone
        integrator_setup_func(minus_clone, [minus_perturbation])
        
        # Run -h₀ clone simulation
        self._run_clone_simulation(minus_clone, minus_tracker, tracker, -1.0, output_times_ps)
        
        # Calculate susceptibility
        times, chi_k = tracker.calculate_susceptibility()
        
        # For now, C_k is placeholder - would come from unperturbed simulation
        C_k = np.ones_like(chi_k)  # Placeholder
        
        return chi_k, C_k
    
    def _run_clone_simulation(self,
                             clone_sim: hoomd.Simulation,
                             clone_tracker: ElapsedTimeTracker,
                             susceptibility_tracker: FDRSusceptibilityTracker,
                             sign: float,
                             output_times_ps: np.ndarray) -> None:
        """Run a single clone simulation.
        
        Parameters
        ----------
        clone_sim : hoomd.Simulation
            Clone simulation
        clone_tracker : ElapsedTimeTracker
            Clone time tracker
        susceptibility_tracker : FDRSusceptibilityTracker
            Susceptibility measurement tracker
        sign : float
            Sign of perturbation (+1 or -1)
        output_times_ps : np.ndarray
            Output time points
        """
        # Add measurement updater
        def measurement_action(timestep):
            susceptibility_tracker.measure_density_correlation(
                clone_sim, clone_tracker, sign
            )
        
        measurement_updater = hoomd.update.CustomUpdater(
            action=hoomd.custom.Action(measurement_action),
            trigger=hoomd.trigger.Periodic(100)  # Measure every 100 steps
        )
        clone_sim.operations.updaters.append(measurement_updater)
        
        # Run clone simulation
        target_steps = int(self.clone_runtime_ps * 10000)  # Assuming 0.1 fs timestep
        clone_sim.run(target_steps)
    
    def _generate_shell_kvectors(self) -> List[np.ndarray]:
        """Generate k-vectors on shell of same magnitude for averaging.
        
        Returns
        -------
        List[np.ndarray]
            List of k-vectors with same magnitude
        """
        kvectors = []
        k_mag = self.k_magnitude
        
        # Generate uniformly distributed points on sphere
        # Using spherical coordinates
        for i in range(self.num_shell_vectors):
            # Random spherical angles
            theta = np.arccos(2 * np.random.random() - 1)  # [0, π]
            phi = 2 * np.pi * np.random.random()           # [0, 2π]
            
            # Convert to Cartesian
            kx = k_mag * np.sin(theta) * np.cos(phi)
            ky = k_mag * np.sin(theta) * np.sin(phi)
            kz = k_mag * np.cos(theta)
            
            kvectors.append(np.array([kx, ky, kz]))
        
        return kvectors
    
    def _save_measurement(self, measurement: FDRMeasurement, filepath: Path) -> None:
        """Save FDR measurement to file.
        
        Parameters
        ----------
        measurement : FDRMeasurement
            Measurement data
        filepath : Path
            Output file path
        """
        header = f"# FDR Measurement Data\n"
        header += f"# t_w = {measurement.waiting_time_ps:.3f} ps\n"
        header += f"# k-vector: [{measurement.kvector[0]:.6f}, {measurement.kvector[1]:.6f}, {measurement.kvector[2]:.6f}]\n"
        header += f"# |k| = {measurement.k_magnitude:.6f}\n"
        header += f"# h₀ = {measurement.metadata['perturbation_amplitude']:.6e}\n"
        header += f"# Shell averaging: {measurement.metadata['shell_averaging']}\n"
        header += f"# Columns: time_ps, chi_k, C_k\n"
        
        data = np.column_stack([
            measurement.times_ps,
            measurement.chi_k,
            measurement.C_k
        ])
        
        np.savetxt(filepath, data, header=header, fmt='%.6e')
        self.logger.info(f"FDR measurement saved to {filepath}")
    
    def save_summary(self, filepath: Optional[str] = None) -> None:
        """Save summary of all FDR measurements.
        
        Parameters
        ----------
        filepath : str, optional
            Output file path (default: auto-generated)
        """
        if filepath is None:
            filepath = self.output_dir / f"fdr_summary_k{self.k_magnitude:.3f}.dat"
        
        if not self.measurements:
            self.logger.warning("No measurements to save")
            return
        
        header = f"# FDR Summary\n"
        header += f"# k-vector: [{self.kvector[0]:.6f}, {self.kvector[1]:.6f}, {self.kvector[2]:.6f}]\n"
        header += f"# |k| = {self.k_magnitude:.6f}\n"
        header += f"# h₀ = {self.perturbation_amplitude:.6e}\n"
        header += f"# Number of waiting times: {len(self.measurements)}\n"
        header += f"# Waiting times (ps): {[m.waiting_time_ps for m in self.measurements]}\n"
        header += f"# Columns: waiting_time_ps, [time_dependent_data would go here]\n"
        
        # For now, just save waiting times
        # A full implementation would include comprehensive FDR data
        waiting_times = np.array([m.waiting_time_ps for m in self.measurements])
        np.savetxt(filepath, waiting_times, header=header, fmt='%.6f')
        
        self.logger.info(f"FDR summary saved to {filepath}")
        self.logger.info(f"Completed {len(self.measurements)} FDR measurements")
