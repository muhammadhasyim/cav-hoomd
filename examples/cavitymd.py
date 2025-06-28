import hoomd
import datetime 
import numpy as np
import gsd.hoomd
from numba import njit
from hoomd.bussi_reservoir.thermostats import BussiReservoir as Bussi
import logging
import sys
import os
import importlib

# Import the new CavityForce from the installed plugin
from hoomd.cavitymd import CavityForce

# Physical constants and unit conversions
class PhysicalConstants:
    HARTREE_TO_CM_MINUS1 = 219474.63
    KB_HARTREE_PER_K = 3.167e-6  # Boltzmann constant in Hartree/K
    ENERGY_JOULES = 4.35974e-18  # Hartree to Joules
    LENGTH_METERS = 5.29177210544e-11  # Bohr to meters
    MASS_KG = 9.1093837139e-31  # Electron mass in kg
    TIME_SECONDS = 2.418884e-17  # Atomic time unit to seconds
    TIME_PS_CONVERSION = 2.418884e-5  # a.u. to picoseconds (corrected from 0.02418884)
    
    @classmethod
    def ps_to_atomic_units(cls, time_ps):
        """
        Convert time from picoseconds to atomic units.
        
        Args:
            time_ps: Time in picoseconds
            
        Returns:
            Time in atomic units
        """
        return time_ps / cls.TIME_PS_CONVERSION
    
    @classmethod
    def atomic_units_to_ps(cls, time_au):
        """
        Convert time from atomic units to picoseconds.
        
        Args:
            time_au: Time in atomic units
            
        Returns:
            Time in picoseconds
        """
        return time_au * cls.TIME_PS_CONVERSION
    
    @classmethod
    def gamma_from_tau_ps(cls, tau_ps):
        """
        Calculate gamma (damping coefficient) from time constant in picoseconds.
        For Langevin dynamics: gamma = 1/tau
        
        Args:
            tau_ps: Time constant in picoseconds
            
        Returns:
            Gamma in atomic units (inverse time)
        """
        if tau_ps <= 0.0:
            raise ValueError(
                f"ERROR: tau_ps must be positive, got {tau_ps} ps.\n"
                f"For Langevin dynamics, gamma = 1/tau, so tau must be > 0.\n"
                f"For overdamped dynamics (tau → 0), use Brownian dynamics instead."
            )
        tau_au = cls.ps_to_atomic_units(tau_ps)
        return 1.0 / tau_au

def unwrap_positions(positions, images, box_lengths):
    """
    Unwrap particle positions across periodic boundaries.
    
    Args:
        positions: Array of wrapped positions (N x 3)
        images: Array of image flags (N x 3)
        box_lengths: Array of box dimensions (3,)
        
    Returns:
        Array of unwrapped positions (N x 3)
    """
    # Convert inputs to cupy arrays if they aren't already
    pos = np.asarray(positions)
    img = np.asarray(images)
    box = np.asarray(box_lengths)
    
    # Unwrap by adding box lengths multiplied by image flags
    return pos + img * box[None, :]

class Status:
    def __init__(self, simulation, chartime, time_tracker=None, runtime_ps=None):
        self.simulation = simulation
        self.chartime = chartime
        self.starttime = datetime.datetime.now()
        self.time_tracker = time_tracker
        self.runtime_ps = runtime_ps  # Total runtime in picoseconds
        self.last_timestep = 0
        self.last_wall_time = datetime.datetime.now()

    @property
    def seconds_remaining(self):
        """Calculate remaining time based on simulation progress and speed."""
        if self.runtime_ps is None or self.time_tracker is None:
            # Fallback to old calculation if runtime not available
            try:
                return (
                    self.simulation.final_timestep - self.simulation.timestep
                ) / self.simulation.tps
            except ZeroDivisionError:
                return 0
        
        # Calculate elapsed simulation time
        elapsed_sim_time_ps = self.time_tracker.elapsed_time
        
        # Calculate remaining simulation time
        remaining_sim_time_ps = self.runtime_ps - elapsed_sim_time_ps
        
        if remaining_sim_time_ps <= 0:
            return 0
        
        # Calculate wall time elapsed
        current_wall_time = datetime.datetime.now()
        wall_time_elapsed_seconds = (current_wall_time - self.starttime).total_seconds()
        
        if wall_time_elapsed_seconds <= 0 or elapsed_sim_time_ps <= 0:
            return 0
        
        # Calculate simulation rate: ps per second of wall time
        ps_per_second = elapsed_sim_time_ps / wall_time_elapsed_seconds
        
        if ps_per_second <= 0:
            return 0
        
        # Calculate remaining wall time
        remaining_wall_time_seconds = remaining_sim_time_ps / ps_per_second
        
        return remaining_wall_time_seconds

    @property
    def etr(self):
        """Estimated time remaining based on simulation progress and speed."""
        return str(datetime.timedelta(seconds=self.seconds_remaining))
    
    @property
    def nsd(self):
        # Calculate nanoseconds per day based on actual simulation progress
        current_timestep = self.simulation.timestep
        if current_timestep <= 0:
            return "0.0"
        
        # Use time_tracker if available for more accurate timing
        if self.time_tracker is not None:
            simulation_time_ps = self.time_tracker.elapsed_time
        else:
            # Fallback to dt * timestep calculation
            dt = float(self.simulation.operations.integrator.dt)
            simulation_time_ps = PhysicalConstants.atomic_units_to_ps(dt * current_timestep)
        
        # Calculate wall time elapsed
        current_wall_time = datetime.datetime.now()
        wall_time_elapsed = (current_wall_time - self.starttime).total_seconds()
        
        if wall_time_elapsed <= 0:
            return "0.0"
        
        # Calculate simulation rate: ps per second of wall time
        ps_per_second = simulation_time_ps / wall_time_elapsed
        
        # Convert to nanoseconds per day
        ns_per_second = ps_per_second / 1000.0  # ps to ns conversion
        ns_per_day = ns_per_second * 86400  # seconds to day conversion
        
        return str(np.round(ns_per_day, 6))
    
    @property
    def Dt(self):
        return str(np.round(float(self.simulation.operations.integrator.dt*self.chartime*1000000),6))
    
    @property
    def elapsed(self):
        curtime = datetime.datetime.now()
        return str(curtime-self.starttime)

class DipoleAutocorrelation(hoomd.custom.Action):
    def __init__(self, simulation, time_tracker=None, output_prefix='dipole_autocorr', output_period=1000):
        super().__init__()
        self.sim = simulation
        self.time_tracker = time_tracker
        self.output_prefix = output_prefix
        self.output_period = output_period
        self.output_file_path = f'{self.output_prefix}_dipole_autocorr.txt'
        self.reference_time = 0.0
        self.reference_dipole = None
        self.last_autocorr_value = None
        self.last_output_step = 0
        self.current_computed_value = None
        self.current_timestep = -1

        # Initialize reference dipole (t=0) and output file
        with self.sim.state.cpu_local_snapshot as snap:
            box_lengths = np.array([
                snap.global_box.L[0],
                snap.global_box.L[1],
                snap.global_box.L[2]
            ])
            uwrap_pos = unwrap_positions(
                snap.particles.position,
                snap.particles.image,
                box_lengths
            )
            self.reference_dipole = np.matmul(snap.particles.charge, np.array(uwrap_pos))
            self.initial_autocorr_value = np.dot(self.reference_dipole, self.reference_dipole)
            self.last_autocorr_value = self.initial_autocorr_value
            self.current_computed_value = self.initial_autocorr_value
            # Write header and t=0 value
            with open(self.output_file_path, 'w') as f:
                f.write(f'# Dipole autocorrelation data\n')
                f.write('# t0(ps) t(ps) C(t)\n')
                f.write(f'# Started at time 0.0 ps\n')
                f.write(f'{0.0:.6f} {0.0:.6f} {self.initial_autocorr_value:.6f}\n')
        print(f"Dipole autocorrelation tracker initialized. Initial C(0) = {self.initial_autocorr_value:.6e}")

    def compute_autocorr(self, current_dipole):
        # C(t) = d(0)·d(t)
        return np.dot(self.reference_dipole, current_dipole)

    def act(self, timestep):
        if timestep == 0:
            return
        with self.sim.state.cpu_local_snapshot as snap:
            box_lengths = np.array([
                snap.global_box.L[0],
                snap.global_box.L[1],
                snap.global_box.L[2]
            ])
            uwrap_pos = unwrap_positions(
                snap.particles.position,
                snap.particles.image,
                box_lengths
            )
            current_dipole = np.matmul(snap.particles.charge, np.array(uwrap_pos))
            if self.time_tracker is not None:
                current_time = self.time_tracker.elapsed_time
            else:
                current_time = timestep  # fallback
        autocorr_value = self.compute_autocorr(current_dipole)
        self.last_autocorr_value = autocorr_value
        self.current_computed_value = autocorr_value
        self.current_timestep = timestep
        self.output_latest_autocorr(current_time, autocorr_value)
        if timestep - self.last_output_step >= self.output_period:
            # print(f"C(t) at time {current_time:.4f} ps: {autocorr_value:.6e}")
            self.last_output_step = timestep

    def output_latest_autocorr(self, current_time, autocorr_value):
        dt = current_time - self.reference_time
        with open(self.output_file_path, 'a') as f:
            f.write(f'{self.reference_time:.6f} {dt:.6f} {autocorr_value:.6f}\n')

    @hoomd.logging.log
    def current_autocorr(self):
        timestep = self.sim.timestep
        if timestep == 0:
            return self.initial_autocorr_value
        if timestep != self.current_timestep:
            with self.sim.state.cpu_local_snapshot as snap:
                box_lengths = np.array([
                    snap.global_box.L[0],
                    snap.global_box.L[1],
                    snap.global_box.L[2]
                ])
                uwrap_pos = unwrap_positions(
                    snap.particles.position,
                    snap.particles.image,
                    box_lengths
                )
                current_dipole = np.matmul(snap.particles.charge, np.array(uwrap_pos))
                autocorr_value = self.compute_autocorr(current_dipole)
                self.current_computed_value = autocorr_value
                self.current_timestep = timestep
                return autocorr_value
        return self.current_computed_value

class EnergyContributionTracker(hoomd.custom.Action):
    def __init__(self, simulation, harmonic, lj, short, long, cavityforce=None, langevin_method=None, molecular_langevin_method=None, cavity_langevin_method=None, mttk_thermostat=None, cavity_mttk_thermostat=None, molecular_bussi_thermostat=None, cavity_bussi_thermostat=None, kinetic_tracker=None, cavity_mode_tracker=None, time_tracker=None, output_prefix='energy', output_period=1, max_time_ps=None):
        """
        Track individual potential and kinetic energy contributions and output to a text file.
        
        Args:
            simulation: HOOMD simulation object
            harmonic: Harmonic bond force object
            lj: Lennard-Jones force object  
            short: Short-range Coulomb force object
            long: Long-range Coulomb force object
            cavityforce: Cavity force object (optional)
            langevin_method: Langevin method for reservoir energy tracking (deprecated, use molecular_langevin_method)
            molecular_langevin_method: Molecular Langevin method for molecular reservoir energy tracking
            cavity_langevin_method: Cavity Langevin method for cavity reservoir energy tracking (optional)
            mttk_thermostat: Molecular MTTK thermostat for thermostat energy tracking (optional)
            cavity_mttk_thermostat: Cavity MTTK thermostat for thermostat energy tracking (optional)
            molecular_bussi_thermostat: Molecular BussiReservoir thermostat for reservoir energy tracking (optional)
            cavity_bussi_thermostat: Cavity BussiReservoir thermostat for reservoir energy tracking (optional)
            kinetic_tracker: KineticEnergyTracker for molecular kinetic energy (optional)
            cavity_mode_tracker: CavityModeTracker for cavity kinetic energy (optional)
            time_tracker: ElapsedTimeTracker for accurate timing
            output_prefix: Prefix for output file name
            output_period: How often to output (every N timesteps)
            max_time_ps: Maximum time in picoseconds to output energy data (None for no limit)
        """
        super().__init__()
        self.sim = simulation
        self.harmonic = harmonic
        self.lj = lj
        self.short = short
        self.long = long
        self.cavityforce = cavityforce
        
        # Handle backward compatibility for langevin_method
        if molecular_langevin_method is None and langevin_method is not None:
            molecular_langevin_method = langevin_method
        
        self.molecular_langevin_method = molecular_langevin_method
        self.cavity_langevin_method = cavity_langevin_method
        self.mttk_thermostat = mttk_thermostat
        self.cavity_mttk_thermostat = cavity_mttk_thermostat
        self.molecular_bussi_thermostat = molecular_bussi_thermostat
        self.cavity_bussi_thermostat = cavity_bussi_thermostat
        self.kinetic_tracker = kinetic_tracker
        self.cavity_mode_tracker = cavity_mode_tracker
        self.time_tracker = time_tracker
        self.output_prefix = output_prefix
        self.output_period = output_period
        self.max_time_ps = max_time_ps
        self.output_file_path = f'{self.output_prefix}_energy_contributions.txt'
        self.last_output_step = 0
        self.output_stopped = False  # Flag to track if output has been stopped
        
        # Initialize output file with header
        with open(self.output_file_path, 'w') as f:
            f.write('# Combined potential and kinetic energy contributions (all total/extensive quantities)\n')
            if max_time_ps is not None:
                f.write(f'# Energy output limited to first {max_time_ps:.2f} ps\n')
            
            # Create comprehensive header based on available components
            header_parts = ['time(ps)', 'timestep']
            
            # Potential energy components
            header_parts.extend(['harmonic_energy', 'lj_energy', 'coulomb_short_energy', 'coulomb_long_energy'])
            
            # Cavity potential energy components (if available)
            if cavityforce is not None:
                header_parts.extend(['cavity_harmonic_energy', 'cavity_coupling_energy', 'cavity_dipole_self_energy', 'cavity_total_potential_energy'])
            
            # Kinetic energy components
            if kinetic_tracker is not None:
                header_parts.append('molecular_kinetic_energy')
            
            if cavity_mode_tracker is not None:
                header_parts.append('cavity_mode_kinetic_energy')
            
            # Separate reservoir energies - always use standardized names
            header_parts.append('molecular_reservoir_energy')
            header_parts.append('cavity_reservoir_energy')
            
            # MTTK thermostat energies
            if mttk_thermostat is not None:
                header_parts.append('molecular_mttk_thermostat_energy')
            
            if cavity_mttk_thermostat is not None:
                header_parts.append('cavity_mttk_thermostat_energy')
            
            # BussiReservoir thermostat energies (additional detailed components)
            if molecular_bussi_thermostat is not None:
                header_parts.extend(['molecular_bussi_translational_energy', 'molecular_bussi_rotational_energy'])
            
            if cavity_bussi_thermostat is not None:
                header_parts.extend(['cavity_bussi_translational_energy', 'cavity_bussi_rotational_energy'])
            
            # Total energies
            header_parts.extend(['total_potential_energy', 'total_kinetic_energy', 'universe_total_energy'])
            
            # Write header
            f.write('# ' + ' '.join(header_parts) + '\n')
        
        print(f"Combined energy tracker initialized. Output file: {self.output_file_path}")
        if max_time_ps is not None:
            print(f"Energy output will stop after {max_time_ps:.2f} ps")
        if molecular_langevin_method is not None:
            print(f"Molecular reservoir energy tracking enabled - method: {type(molecular_langevin_method).__name__}")
            print(f"  Molecular Langevin tally_reservoir_energy: {getattr(molecular_langevin_method, 'tally_reservoir_energy', 'N/A')}")
        if cavity_langevin_method is not None:
            print(f"Cavity reservoir energy tracking enabled - method: {type(cavity_langevin_method).__name__}")
            print(f"  Cavity Langevin tally_reservoir_energy: {getattr(cavity_langevin_method, 'tally_reservoir_energy', 'N/A')}")
        if mttk_thermostat is not None:
            print(f"Molecular MTTK thermostat energy tracking enabled")
        if cavity_mttk_thermostat is not None:
            print(f"Cavity MTTK thermostat energy tracking enabled")
        if molecular_bussi_thermostat is not None:
            print(f"Molecular BussiReservoir thermostat energy tracking enabled")
        if cavity_bussi_thermostat is not None:
            print(f"Cavity BussiReservoir thermostat energy tracking enabled")
        if kinetic_tracker is not None:
            print(f"Molecular kinetic energy tracking enabled")
        if cavity_mode_tracker is not None:
            print(f"Cavity mode kinetic energy tracking enabled")

    def act(self, timestep):
        # Check if output has been stopped due to time limit
        if self.output_stopped:
            return
            
        if timestep - self.last_output_step >= self.output_period:
            # Get current time
            if self.time_tracker is not None:
                current_time = self.time_tracker.elapsed_time
            else:
                current_time = timestep * self.sim.operations.integrator.dt * 0.02418884  # Convert to ps
            
            # Check if we've exceeded the time limit
            if self.max_time_ps is not None and current_time > self.max_time_ps:
                if not self.output_stopped:
                    print(f"Energy output stopped: reached time limit of {self.max_time_ps:.2f} ps at t={current_time:.4f} ps")
                    self.output_stopped = True
                return
            
            # Get individual potential energy contributions
            harmonic_energy = self.harmonic.energy
            lj_energy = self.lj.energy
            short_energy = self.short.energy
            long_energy = self.long.energy
            
            # Calculate total potential energy (without cavity)
            total_potential_energy = harmonic_energy + lj_energy + short_energy + long_energy
            
            # Get cavity potential energy components if present
            cavity_harmonic_energy = 0.0
            cavity_coupling_energy = 0.0
            cavity_dipole_self_energy = 0.0
            cavity_total_potential_energy = 0.0
            
            if self.cavityforce is not None:
                # Use local CavityForce implementation (direct attributes)
                try:
                    cavity_harmonic_energy = getattr(self.cavityforce, 'harmonic_energy', 0.0)
                    cavity_coupling_energy = getattr(self.cavityforce, 'coupling_energy', 0.0)
                    cavity_dipole_self_energy = getattr(self.cavityforce, 'dipole_self_energy', 0.0)
                    # For total energy, try .energy property first, then sum components
                    if hasattr(self.cavityforce, 'energy'):
                        cavity_total_potential_energy = self.cavityforce.energy
                    else:
                        cavity_total_potential_energy = cavity_harmonic_energy + cavity_coupling_energy + cavity_dipole_self_energy
                except Exception as e:
                    # If any error occurs, set all to zero and continue
                    cavity_harmonic_energy = 0.0
                    cavity_coupling_energy = 0.0
                    cavity_dipole_self_energy = 0.0
                    cavity_total_potential_energy = 0.0
                    # Only print error for first few timesteps to avoid spam
                    if timestep < 50:
                        print(f"Warning: Error accessing cavity energy at timestep {timestep}: {e}")
                
                total_potential_energy += cavity_total_potential_energy
            
            # Get kinetic energy components
            molecular_kinetic_energy = 0.0
            if self.kinetic_tracker is not None:
                molecular_kinetic_energy = self.kinetic_tracker.kinetic_energy
            
            cavity_mode_kinetic_energy = 0.0
            if self.cavity_mode_tracker is not None:
                cavity_mode_kinetic_energy = self.cavity_mode_tracker.cavity_kinetic_energy
            
            total_kinetic_energy = molecular_kinetic_energy + cavity_mode_kinetic_energy
            
            # Get reservoir energy if Langevin method is available
            molecular_reservoir_energy = 0.0
            if self.molecular_langevin_method is not None:
                try:
                    molecular_reservoir_energy = self.molecular_langevin_method.reservoir_energy
                except AttributeError:
                    # molecular_reservoir_energy not available yet (before first simulation step)
                    molecular_reservoir_energy = 0.0
            
            cavity_reservoir_energy = 0.0
            if self.cavity_langevin_method is not None:
                try:
                    cavity_reservoir_energy = self.cavity_langevin_method.reservoir_energy
                except AttributeError:
                    # cavity_reservoir_energy not available yet (before first simulation step)
                    cavity_reservoir_energy = 0.0
            
            # Get MTTK thermostat energies if available
            molecular_mttk_thermostat_energy = 0.0
            if self.mttk_thermostat is not None:
                try:
                    molecular_mttk_thermostat_energy = self.mttk_thermostat.energy
                except AttributeError:
                    # MTTK energy not available yet (before first simulation step)
                    molecular_mttk_thermostat_energy = 0.0
            
            cavity_mttk_thermostat_energy = 0.0
            if self.cavity_mttk_thermostat is not None:
                try:
                    cavity_mttk_thermostat_energy = self.cavity_mttk_thermostat.energy
                except AttributeError:
                    # Cavity MTTK energy not available yet (before first simulation step)
                    cavity_mttk_thermostat_energy = 0.0
            
            # Get BussiReservoir thermostat energies if available
            molecular_bussi_total_reservoir_energy = 0.0
            molecular_bussi_translational_energy = 0.0
            molecular_bussi_rotational_energy = 0.0
            if self.molecular_bussi_thermostat is not None:
                try:
                    molecular_bussi_total_reservoir_energy = self.molecular_bussi_thermostat.total_reservoir_energy
                    molecular_bussi_translational_energy = self.molecular_bussi_thermostat.reservoir_energy_translational
                    molecular_bussi_rotational_energy = self.molecular_bussi_thermostat.reservoir_energy_rotational
                except (AttributeError, hoomd.error.DataAccessError):
                    # BussiReservoir energy not available yet (before first simulation step)
                    molecular_bussi_total_reservoir_energy = 0.0
                    molecular_bussi_translational_energy = 0.0
                    molecular_bussi_rotational_energy = 0.0
            
            cavity_bussi_total_reservoir_energy = 0.0
            cavity_bussi_translational_energy = 0.0
            cavity_bussi_rotational_energy = 0.0
            if self.cavity_bussi_thermostat is not None:
                try:
                    cavity_bussi_total_reservoir_energy = self.cavity_bussi_thermostat.total_reservoir_energy
                    cavity_bussi_translational_energy = self.cavity_bussi_thermostat.reservoir_energy_translational
                    cavity_bussi_rotational_energy = self.cavity_bussi_thermostat.reservoir_energy_rotational
                except (AttributeError, hoomd.error.DataAccessError):
                    # BussiReservoir energy not available yet (before first simulation step)
                    cavity_bussi_total_reservoir_energy = 0.0
                    cavity_bussi_translational_energy = 0.0
                    cavity_bussi_rotational_energy = 0.0
            
            # Calculate total energy (system only, without reservoir)
            system_total_energy = total_potential_energy + total_kinetic_energy
            
            # Calculate total thermostat/reservoir energy from all sources
            total_thermostat_energy = (molecular_reservoir_energy + 
                                     cavity_reservoir_energy + 
                                     molecular_mttk_thermostat_energy + 
                                     cavity_mttk_thermostat_energy + 
                                     molecular_bussi_total_reservoir_energy + 
                                     cavity_bussi_total_reservoir_energy)
            
            # Calculate universe total energy (system + all thermostats)
            universe_total_energy = system_total_energy + total_thermostat_energy
            
            # Build output line based on available components
            output_values = [current_time, timestep]
            
            # Potential energy components
            output_values.extend([harmonic_energy, lj_energy, short_energy, long_energy])
            
            # Cavity potential energy components (if available)
            if self.cavityforce is not None:
                output_values.extend([cavity_harmonic_energy, cavity_coupling_energy, 
                                    cavity_dipole_self_energy, cavity_total_potential_energy])
            
            # Kinetic energy components
            if self.kinetic_tracker is not None:
                output_values.append(molecular_kinetic_energy)
            
            if self.cavity_mode_tracker is not None:
                output_values.append(cavity_mode_kinetic_energy)
            
            # Combined reservoir energies - always output standardized columns
            # For molecular: combine Langevin + Bussi total reservoir energy
            combined_molecular_reservoir = molecular_reservoir_energy + molecular_bussi_total_reservoir_energy
            output_values.append(combined_molecular_reservoir)
            
            # For cavity: combine Langevin + Bussi total reservoir energy  
            combined_cavity_reservoir = cavity_reservoir_energy + cavity_bussi_total_reservoir_energy
            output_values.append(combined_cavity_reservoir)
            
            # MTTK thermostat energies
            if self.mttk_thermostat is not None:
                output_values.append(molecular_mttk_thermostat_energy)
            
            if self.cavity_mttk_thermostat is not None:
                output_values.append(cavity_mttk_thermostat_energy)
            
            # BussiReservoir detailed components (if available)
            if self.molecular_bussi_thermostat is not None:
                output_values.extend([molecular_bussi_translational_energy,
                                     molecular_bussi_rotational_energy])
            
            if self.cavity_bussi_thermostat is not None:
                output_values.extend([cavity_bussi_translational_energy,
                                     cavity_bussi_rotational_energy])
            
            # Total energies
            output_values.extend([total_potential_energy, total_kinetic_energy, universe_total_energy])
            
            # Write to file with better formatting and error handling
            try:
                with open(self.output_file_path, 'a') as f:
                    formatted_values = [f'{val:.6f}' if isinstance(val, float) else str(val) for val in output_values]
                    f.write(' '.join(formatted_values) + '\n')
            except IOError as e:
                print(f"Warning: Failed to write energy data at timestep {timestep}: {e}")
            
            self.last_output_step = timestep

    @hoomd.logging.log
    def total_potential_energy(self):
        """Calculate total potential energy from individual contributions."""
        try:
            total = self.harmonic.energy + self.lj.energy + self.short.energy + self.long.energy
            if self.cavityforce is not None:
                # Use local CavityForce implementation (direct attributes)
                try:
                    cavity_harmonic_energy = getattr(self.cavityforce, 'harmonic_energy', 0.0)
                    cavity_coupling_energy = getattr(self.cavityforce, 'coupling_energy', 0.0)
                    cavity_dipole_self_energy = getattr(self.cavityforce, 'dipole_self_energy', 0.0)
                    # For total energy, try .energy property first, then sum components
                    if hasattr(self.cavityforce, 'energy'):
                        cavity_total_potential_energy = self.cavityforce.energy
                    else:
                        cavity_total_potential_energy = cavity_harmonic_energy + cavity_coupling_energy + cavity_dipole_self_energy
                    total += cavity_total_potential_energy
                except Exception:
                    # If cavity energy access fails, just continue without it
                    pass
            return total
        except hoomd.error.DataAccessError:
            # Energy not available yet (before first simulation step)
            return 0.0

class KineticEnergyTracker(hoomd.custom.Action):
    """
    Tracks the total kinetic energy of molecular particles (excluding photon).
    Computes KE = (1/2) * Σ_i m_i * v_i² for all molecular particles.
    """
    def __init__(self, simulation, time_tracker=None, output_prefix='kinetic_energy', output_period=1):
        super().__init__()
        self.simulation = simulation
        self.time_tracker = time_tracker
        self.output_prefix = output_prefix
        self.output_period = output_period
        self.output_file_path = f'{self.output_prefix}_kinetic_energy.txt'
        self.last_output_step = 0
        
        # For caching computed values (fixes logging delay issue)
        self.current_kinetic_energy = 0.0
        self.current_temperature = 0.0
        self.current_timestep = -1
        
        # Initialize output file with header
        with open(self.output_file_path, 'w') as f:
            f.write('# Molecular kinetic energy data (excluding photon particle)\n')
            f.write('# time(ps) timestep kinetic_energy(a.u.) temperature(K) avg_vel_magnitude(a.u.)\n')
        
        print(f"KineticEnergyTracker initialized, output file: {self.output_file_path}")

    def compute_kinetic_energy(self):
        """
        Compute total kinetic energy and related properties for molecular particles.
        Returns kinetic energy, temperature, and average velocity magnitude.
        """
        with self.simulation.state.cpu_local_snapshot as snap:
            # Filter out photon particle (typeid == 2)
            molecular_mask = snap.particles.typeid != 2
            
            if not np.any(molecular_mask):
                return 0.0, 0.0, 0.0
            
            # Get molecular particle properties
            masses = snap.particles.mass[molecular_mask]
            velocities = snap.particles.velocity[molecular_mask]
            
            # Compute kinetic energy: KE = (1/2) * m * v²
            vel_squared = np.sum(velocities**2, axis=1)  # |v|² for each particle
            kinetic_energies = 0.5 * masses * vel_squared
            total_kinetic_energy = np.sum(kinetic_energies)
            
            # Compute instantaneous temperature from kinetic energy
            # T = (2/3) * KE_total / (N * kB) for 3D system
            n_molecular_particles = np.sum(molecular_mask)
            kB = 3.167e-6  # Hartree/K (same as in main simulation)
            
            if n_molecular_particles > 0:
                temperature = (2.0/3.0) * total_kinetic_energy / (n_molecular_particles * kB)
            else:
                temperature = 0.0
            
            # Compute average velocity magnitude
            vel_magnitudes = np.sqrt(vel_squared)
            avg_vel_magnitude = np.mean(vel_magnitudes) if len(vel_magnitudes) > 0 else 0.0
            
            return total_kinetic_energy, temperature, avg_vel_magnitude

    def act(self, timestep):
        if timestep - self.last_output_step >= self.output_period:
            # Get current time
            if self.time_tracker is not None:
                current_time = self.time_tracker.elapsed_time
            else:
                current_time = timestep * self.sim.operations.integrator.dt * 0.02418884  # Convert to ps
            
            # Compute kinetic energy and related properties
            kinetic_energy, temperature, avg_vel_magnitude = self.compute_kinetic_energy()
            
            # Cache values for logging
            self.current_kinetic_energy = kinetic_energy
            self.current_temperature = temperature
            self.current_timestep = timestep
            
            # Write to file
            with open(self.output_file_path, 'a') as f:
                f.write(f'{current_time:.6f} {timestep} {kinetic_energy:.6f} {temperature:.6f} {avg_vel_magnitude:.6f}\n')
            
            self.last_output_step = timestep

    @hoomd.logging.log
    def kinetic_energy(self):
        """Return total kinetic energy of molecular particles for logging."""
        timestep = self.simulation.timestep
        if timestep != self.current_timestep:
            # Need to compute current value
            kinetic_energy, temperature, avg_vel_magnitude = self.compute_kinetic_energy()
            self.current_kinetic_energy = kinetic_energy
            self.current_temperature = temperature
            self.current_timestep = timestep
        
        return self.current_kinetic_energy

    @hoomd.logging.log
    def temperature(self):
        """Return instantaneous temperature from kinetic energy for logging."""
        timestep = self.simulation.timestep
        if timestep != self.current_timestep:
            # Compute if not already cached
            kinetic_energy, temperature, avg_vel_magnitude = self.compute_kinetic_energy()
            self.current_kinetic_energy = kinetic_energy
            self.current_temperature = temperature
            self.current_timestep = timestep
        return self.current_temperature

class CavityModeTracker(hoomd.custom.Action):
    """
    Tracks the kinetic energy and other properties of the cavity mode (photon particle).
    Computes cavity mode KE, PE (harmonic only), and total cavity oscillator energy.
    """
    def __init__(self, simulation, cavityforce, time_tracker=None, output_prefix='cavity_mode', output_period=1):
        super().__init__()
        self.simulation = simulation
        self.cavityforce = cavityforce
        self.time_tracker = time_tracker
        self.output_prefix = output_prefix
        self.output_period = output_period
        self.output_file_path = f'{self.output_prefix}_cavity_mode.txt'
        self.last_output_step = 0
        
        # For caching computed values (fixes logging delay issue)
        self.current_cavity_ke = 0.0
        self.current_cavity_pe_harmonic = 0.0
        self.current_cavity_total_energy = 0.0
        self.current_cavity_position = np.zeros(3)
        self.current_cavity_velocity = np.zeros(3)
        self.current_timestep = -1
        
        # Initialize output file with header
        with open(self.output_file_path, 'w') as f:
            f.write('# Cavity mode energy data\n')
            f.write('# time(ps) timestep cavity_ke(a.u.) cavity_pe_harmonic(a.u.) cavity_total_energy(a.u.) cavity_temperature(K) cavity_pos_x cavity_pos_y cavity_pos_z cavity_vel_x cavity_vel_y cavity_vel_z\n')
        
        print(f"CavityModeTracker initialized, output file: {self.output_file_path}")

    def compute_cavity_properties(self):
        """
        Compute cavity mode kinetic energy, harmonic potential energy, and other properties.
        Returns cavity KE, harmonic PE, total energy, position, and velocity.
        """
        with self.simulation.state.cpu_local_snapshot as snap:
            # Find photon particle (typeid == 2)
            photon_mask = snap.particles.typeid == 2
            
            if not np.any(photon_mask):
                return 0.0, 0.0, 0.0, np.zeros(3), np.zeros(3)
            
            # Get photon properties with unwrapped positions
            photon_mass = snap.particles.mass[photon_mask][0]
            photon_velocity = snap.particles.velocity[photon_mask][0]
            
            # Get unwrapped position for the photon
            box_lengths = np.array([
                snap.global_box.L[0],
                snap.global_box.L[1],
                snap.global_box.L[2]
            ])
            unwrapped_positions = unwrap_positions(
                snap.particles.position,
                snap.particles.image,
                box_lengths
            )
            photon_position = unwrapped_positions[photon_mask][0]
            
            # Compute cavity kinetic energy: KE = (1/2) * m * v²
            vel_squared = np.sum(photon_velocity**2)
            cavity_ke = 0.5 * photon_mass * vel_squared
            
            # Get harmonic potential energy from cavity force (without coupling and self-energy terms)
            # Use local CavityForce with direct attribute access
            if hasattr(self.cavityforce, 'harmonic_energy'):
                cavity_pe_harmonic = self.cavityforce.harmonic_energy
            else:
                cavity_pe_harmonic = 0.0
            
            # Total cavity oscillator energy (KE + harmonic PE only)
            cavity_total_energy = cavity_ke + cavity_pe_harmonic
            
            return cavity_ke, cavity_pe_harmonic, cavity_total_energy, photon_position, photon_velocity

    def act(self, timestep):
        if timestep - self.last_output_step >= self.output_period:
            # Get current time
            if self.time_tracker is not None:
                current_time = self.time_tracker.elapsed_time
            else:
                current_time = timestep * self.sim.operations.integrator.dt * 0.02418884  # Convert to ps
            
            # Compute cavity properties
            cavity_ke, cavity_pe_harmonic, cavity_total_energy, cavity_position, cavity_velocity = self.compute_cavity_properties()
            
            # Calculate cavity temperature
            kB = 3.167e-6  # Hartree/K (same as in main simulation)
            if cavity_ke > 0:
                # For 3D motion: <KE> = (3/2) * kB * T, so T = (2/3) * KE / kB
                # The cavity particle now has full 3D motion capability
                cavity_temperature = (2.0/3.0) * cavity_ke / kB
            else:
                cavity_temperature = 0.0
            
            # Cache values for logging
            self.current_cavity_ke = cavity_ke
            self.current_cavity_pe_harmonic = cavity_pe_harmonic
            self.current_cavity_total_energy = cavity_total_energy
            self.current_cavity_position = cavity_position.copy()
            self.current_cavity_velocity = cavity_velocity.copy()
            self.current_timestep = timestep
            
            # Write to file
            with open(self.output_file_path, 'a') as f:
                f.write(f'{current_time:.6f} {timestep} {cavity_ke:.6f} {cavity_pe_harmonic:.6f} {cavity_total_energy:.6f} {cavity_temperature:.6f} '
                       f'{cavity_position[0]:.6f} {cavity_position[1]:.6f} {cavity_position[2]:.6f} '
                       f'{cavity_velocity[0]:.6f} {cavity_velocity[1]:.6f} {cavity_velocity[2]:.6f}\n')
            
            self.last_output_step = timestep

    @hoomd.logging.log
    def cavity_kinetic_energy(self):
        """Return cavity mode kinetic energy for logging."""
        timestep = self.simulation.timestep
        if timestep != self.current_timestep:
            # Need to compute current value
            cavity_ke, cavity_pe_harmonic, cavity_total_energy, cavity_position, cavity_velocity = self.compute_cavity_properties()
            self.current_cavity_ke = cavity_ke
            self.current_cavity_pe_harmonic = cavity_pe_harmonic
            self.current_cavity_total_energy = cavity_total_energy
            self.current_cavity_position = cavity_position.copy()
            self.current_cavity_velocity = cavity_velocity.copy()
            self.current_timestep = timestep
        
        return self.current_cavity_ke

    @hoomd.logging.log
    def cavity_potential_energy_harmonic(self):
        """Return cavity mode harmonic potential energy for logging."""
        timestep = self.simulation.timestep
        if timestep != self.current_timestep:
            # Compute if not already cached
            cavity_ke, cavity_pe_harmonic, cavity_total_energy, cavity_position, cavity_velocity = self.compute_cavity_properties()
            self.current_cavity_ke = cavity_ke
            self.current_cavity_pe_harmonic = cavity_pe_harmonic
            self.current_cavity_total_energy = cavity_total_energy
            self.current_cavity_position = cavity_position.copy()
            self.current_cavity_velocity = cavity_velocity.copy()
            self.current_timestep = timestep
        return self.current_cavity_pe_harmonic

    @hoomd.logging.log
    def cavity_total_energy(self):
        """Return total cavity oscillator energy (KE + harmonic PE only) for logging."""
        timestep = self.simulation.timestep
        if timestep != self.current_timestep:
            # Compute if not already cached
            cavity_ke, cavity_pe_harmonic, cavity_total_energy, cavity_position, cavity_velocity = self.compute_cavity_properties()
            self.current_cavity_ke = cavity_ke
            self.current_cavity_pe_harmonic = cavity_pe_harmonic
            self.current_cavity_total_energy = cavity_total_energy
            self.current_cavity_position = cavity_position.copy()
            self.current_cavity_velocity = cavity_velocity.copy()
            self.current_timestep = timestep
        return self.current_cavity_total_energy

    @hoomd.logging.log
    def cavity_temperature(self):
        """Return cavity mode temperature in Kelvin for logging."""
        timestep = self.simulation.timestep
        if timestep != self.current_timestep:
            # Compute if not already cached
            cavity_ke, cavity_pe_harmonic, cavity_total_energy, cavity_position, cavity_velocity = self.compute_cavity_properties()
            self.current_cavity_ke = cavity_ke
            self.current_cavity_pe_harmonic = cavity_pe_harmonic
            self.current_cavity_total_energy = cavity_total_energy
            self.current_cavity_position = cavity_position.copy()
            self.current_cavity_velocity = cavity_velocity.copy()
            self.current_timestep = timestep
        
        # Calculate temperature from kinetic energy
        # For 3D motion: <KE> = (3/2) * kB * T, so T = (2/3) * KE / kB
        # The cavity particle now has full 3D motion capability
        kB = 3.167e-6  # Hartree/K (same as in main simulation)
        
        if self.current_cavity_ke > 0:
            temperature = (2.0/3.0) * self.current_cavity_ke / kB
            
            # Debug output for first few timesteps
            timestep = self.simulation.timestep
            if timestep < 100 and timestep % 10 == 0:
                print(f"Debug cavity temp at timestep {timestep}: KE={self.current_cavity_ke:.6f} a.u., T={temperature:.1f} K")
        else:
            temperature = 0.0
        
        return temperature

class IndividualBondTracker(hoomd.custom.Action):
    """
    Tracks individual harmonic bond energies to capture fundamental vibrational frequencies.
    This allows us to see O-O and N-N bond frequencies separately rather than the total sum.
    """
    def __init__(self, simulation, harmonic, time_tracker=None, output_prefix='individual_bonds', output_period=1, max_time_ps=None):
        super().__init__()
        self.simulation = simulation
        self.harmonic = harmonic
        self.time_tracker = time_tracker
        self.output_prefix = output_prefix
        self.output_period = output_period
        self.max_time_ps = max_time_ps
        self.output_file_path = f'{self.output_prefix}_individual_bonds.txt'
        self.last_output_step = 0
        self.output_stopped = False
        
        # For caching computed values
        self.current_bond_energies = {}
        self.current_timestep = -1
        
        # Initialize output file with header
        with open(self.output_file_path, 'w') as f:
            f.write('# Individual harmonic bond energies (extensive quantities)\n')
            if max_time_ps is not None:
                f.write(f'# Energy output limited to first {max_time_ps:.2f} ps\n')
            f.write('# time(ps) timestep total_bonds oo_bonds_energy nn_bonds_energy oo_bonds_count nn_bonds_count avg_oo_bond_energy avg_nn_bond_energy\n')
        
        print(f"IndividualBondTracker initialized. Output file: {self.output_file_path}")
        if max_time_ps is not None:
            print(f"Individual bond energy output will stop after {max_time_ps:.2f} ps")

    def compute_individual_bond_energies(self):
        """
        Compute individual bond energies by bond type.
        Returns total bond count, O-O energy, N-N energy, O-O count, N-N count.
        """
        with self.simulation.state.cpu_local_snapshot as snap:
            # Get bond information
            if not hasattr(snap.bonds, 'typeid') or len(snap.bonds.typeid) == 0:
                return 0, 0.0, 0.0, 0, 0, 0.0, 0.0
            
            # Get particle positions and types
            box_lengths = np.array([
                snap.global_box.L[0],
                snap.global_box.L[1],
                snap.global_box.L[2]
            ])
            positions = unwrap_positions(
                snap.particles.position,
                snap.particles.image,
                box_lengths
            )
            
            # Get bond parameters from the harmonic force
            # O-O bonds: type 0, N-N bonds: type 1 (typically)
            bond_types = snap.bonds.typeid
            bond_groups = snap.bonds.group
            
            oo_energy = 0.0
            nn_energy = 0.0
            oo_count = 0
            nn_count = 0
            
            # Iterate through all bonds
            for i, bond_type in enumerate(bond_types):
                # Get the two particles in this bond
                particle_i = bond_groups[i][0]
                particle_j = bond_groups[i][1]
                
                # Get positions
                pos_i = positions[particle_i]
                pos_j = positions[particle_j]
                
                # Calculate bond distance
                dr = pos_j - pos_i
                # Handle periodic boundary conditions
                dr = dr - np.round(dr / box_lengths) * box_lengths
                r = np.linalg.norm(dr)
                
                # Get bond parameters based on type
                if bond_type == 0:  # O-O bonds
                    k = 2*0.36602  # from your harmonic.params['O-O']
                    r0 = 2.281655158
                    bond_energy = 0.5 * k * (r - r0)**2
                    oo_energy += bond_energy
                    oo_count += 1
                elif bond_type == 1:  # N-N bonds
                    k = 2*0.71625  # from your harmonic.params['N-N']
                    r0 = 2.0743522177
                    bond_energy = 0.5 * k * (r - r0)**2
                    nn_energy += bond_energy
                    nn_count += 1
            
            total_bonds = oo_count + nn_count
            avg_oo_energy = oo_energy / oo_count if oo_count > 0 else 0.0
            avg_nn_energy = nn_energy / nn_count if nn_count > 0 else 0.0
            
            return total_bonds, oo_energy, nn_energy, oo_count, nn_count, avg_oo_energy, avg_nn_energy

    def act(self, timestep):
        # Check if output has been stopped due to time limit
        if self.output_stopped:
            return
            
        if timestep - self.last_output_step >= self.output_period:
            # Get current time
            if self.time_tracker is not None:
                current_time = self.time_tracker.elapsed_time
            else:
                current_time = timestep * self.sim.operations.integrator.dt * 0.02418884  # Convert to ps
            
            # Check if we've exceeded the time limit
            if self.max_time_ps is not None and current_time > self.max_time_ps:
                if not self.output_stopped:
                    print(f"Individual bond energy output stopped: reached time limit of {self.max_time_ps:.2f} ps at t={current_time:.4f} ps")
                    self.output_stopped = True
                return
            
            # Compute individual bond energies
            total_bonds, oo_energy, nn_energy, oo_count, nn_count, avg_oo_energy, avg_nn_energy = self.compute_individual_bond_energies()
            
            # Cache values for logging
            self.current_bond_energies = {
                'total_bonds': total_bonds,
                'oo_energy': oo_energy,
                'nn_energy': nn_energy,
                'oo_count': oo_count,
                'nn_count': nn_count,
                'avg_oo_energy': avg_oo_energy,
                'avg_nn_energy': avg_nn_energy
            }
            self.current_timestep = timestep
            
            # Write to file
            with open(self.output_file_path, 'a') as f:
                f.write(f'{current_time:.6f} {timestep} {total_bonds} {oo_energy:.6f} {nn_energy:.6f} {oo_count} {nn_count} {avg_oo_energy:.6f} {avg_nn_energy:.6f}\n')
            
            self.last_output_step = timestep

    @hoomd.logging.log
    def oo_bonds_energy(self):
        """Return total O-O bond energy for logging."""
        timestep = self.simulation.timestep
        if timestep != self.current_timestep:
            # Need to compute current values
            total_bonds, oo_energy, nn_energy, oo_count, nn_count, avg_oo_energy, avg_nn_energy = self.compute_individual_bond_energies()
            self.current_bond_energies = {
                'total_bonds': total_bonds,
                'oo_energy': oo_energy,
                'nn_energy': nn_energy,
                'oo_count': oo_count,
                'nn_count': nn_count,
                'avg_oo_energy': avg_oo_energy,
                'avg_nn_energy': avg_nn_energy
            }
            self.current_timestep = timestep
        return self.current_bond_energies['oo_energy']

    @hoomd.logging.log
    def nn_bonds_energy(self):
        """Return total N-N bond energy for logging."""
        timestep = self.simulation.timestep
        if timestep != self.current_timestep:
            self.oo_bonds_energy()  # This will compute and cache all values
        return self.current_bond_energies['nn_energy']

    @hoomd.logging.log
    def avg_oo_bond_energy(self):
        """Return average O-O bond energy for logging."""
        timestep = self.simulation.timestep
        if timestep != self.current_timestep:
            self.oo_bonds_energy()  # This will compute and cache all values
        return self.current_bond_energies['avg_oo_energy']

    @hoomd.logging.log
    def avg_nn_bond_energy(self):
        """Return average N-N bond energy for logging."""
        timestep = self.simulation.timestep
        if timestep != self.current_timestep:
            self.oo_bonds_energy()  # This will compute and cache all values
        return self.current_bond_energies['avg_nn_energy']

class ElapsedTimeTracker(hoomd.custom.Action):
    """Tracks the total elapsed time in a simulation with variable timesteps."""
    def __init__(self, simulation, runtime):
        super().__init__()
        self.simulation = simulation
        self.total_time = 0.0
        self.runtime = runtime
        self.last_timestep = 0  # Start from 0, not simulation.timestep
        self.initial_timestep = None  # Track the starting timestep to handle inheritance
        
    def act(self, timestep):
        """Update the total elapsed time by accumulating time increments."""
        # Get current timestep size
        dt = self.simulation.operations.integrator.dt
        
        # For the first call, handle initialization
        if self.last_timestep == 0:
            # Initialize - record the starting timestep but don't add its time
            self.initial_timestep = timestep
            self.last_timestep = timestep
            self.total_time = 0.0  # Always start elapsed time from 0, regardless of inherited timestep
            if timestep > 0:
                print(f"NOTICE: Starting from inherited timestep {timestep}")
                print(f"  Elapsed time will start from 0, not from inherited simulation time")
            return
        
        # Calculate time increment since last update
        if timestep > self.last_timestep:
            timestep_increment = timestep - self.last_timestep
            time_increment = timestep_increment * dt
            self.total_time += time_increment
            

        
        # Update last timestep for next iteration
        self.last_timestep = timestep
        
        # Check if we've reached the runtime and exit if so
        if PhysicalConstants.atomic_units_to_ps(self.total_time) >= self.runtime:
            print(f"Runtime {self.runtime} ps reached. Exiting simulation.")
            import sys
            sys.exit(0)

    @hoomd.logging.log
    def elapsed_time(self):
        """Expose the total elapsed time as a property in (ps)."""
        return PhysicalConstants.atomic_units_to_ps(self.total_time)

class DensityCorrelationTracker(hoomd.custom.Action):
    """
    Tracks the density-density correlation function F(k,t) during simulation.
    Computes F(k,t) = <ρₖ(t)·ρₖ*(t₀)> for multiple reference frames.
    Each reference frame gets its own output file, with lag time as the x-axis.
    """
    def __init__(self, simulation, time_tracker, kmag=1.0, num_wavevectors=50, output_period=1000, output_prefix='fskt', reference_interval_ps=1.0, max_references=10, file_buffer_size=1000):
        super().__init__()
        self.simulation = simulation
        self.time_tracker = time_tracker  # Reference to ElapsedTimeTracker for accurate timing
        self.kmag = kmag  # Wavevector magnitude
        self.num_wavevectors = num_wavevectors  # Number of wavevectors to sample on sphere
        self.output_period = output_period  # How often to output summary results
        self.output_prefix = output_prefix  # Prefix for output files
        self.reference_interval_ps = reference_interval_ps
        self.max_references = max_references
        self.file_buffer_size = file_buffer_size  # Buffer size for file operations

        # Generate wavevectors
        self.wavevectors = self._fibonacci_sphere(samples=num_wavevectors) * kmag

        # List of reference frames: each is a dict with keys: 'timestep', 'time', 'rhok_real', 'file_path', 'buffer'
        self.references = []
        self.last_reference_time = None

        # Counter for tracking when to output results
        self.last_output_step = 0

        # For caching computed values in current timestep (fixes logging delay issue)
        self.current_computed_value = None
        self.current_timestep = -1
        
        # Cache for current rhok to avoid recomputation
        self.cached_rhok_real = None
        self.cached_rhok_imag = None
        self.cached_rhok_timestep = -1

        # Initialize with the first reference at t=0
        with self.simulation.state.cpu_local_snapshot as snap:
            box_lengths = np.array([
                snap.global_box.L[0],
                snap.global_box.L[1],
                snap.global_box.L[2]
            ])
            positions = unwrap_positions(
                snap.particles.position,
                snap.particles.image,
                box_lengths
            )
            particle_mask = snap.particles.typeid != 2
            initial_positions = positions[particle_mask]
            rhok_real, rhok_imag = self.compute_rhok(initial_positions)
            t0 = 0.0
            file_path = f'{self.output_prefix}_fskt_k{self.kmag:.2f}_ref0.txt'
            
            # Initialize buffer for this reference
            buffer_data = []
            
            # Write header to file
            try:
                with open(file_path, 'w') as f:
                    f.write(f'# F(k,t) data for k={self.kmag:.4f}, t0={t0:.6f} ps\n')
                    f.write('# lag_time(ps) F(k,t)\n')
            except IOError as e:
                print(f"Warning: Failed to initialize F(k,t) file {file_path}: {e}")
            
            self.references.append({
                'timestep': 0,
                'time': t0,
                'rhok_real': rhok_real.copy(),
                'rhok_imag': rhok_imag.copy(),
                'file_path': file_path,
                'buffer': buffer_data
            })
            self.last_reference_time = t0
        print(f"DensityCorrelationTracker initialized with k={kmag:.2f}, {num_wavevectors} wavevectors, reference_interval_ps={reference_interval_ps}, max_references={max_references}")

    def _fibonacci_sphere(self, samples=100):
        """
        Generates points that uniformly sample the surface of a sphere using the Fibonacci lattice method.
        """
        points = np.zeros((samples, 3))
        phi = np.pi * (3. - np.sqrt(5.))  # golden angle in radians
        
        for i in range(samples):
            y = 1 - (i / float(samples - 1)) * 2
            radius = np.sqrt(1 - y * y)
            theta = phi * i
            x = np.cos(theta) * radius
            z = np.sin(theta) * radius
            points[i] = np.array([x, y, z])

        return points

    def compute_rhok(self, positions):
        n_particles = positions.shape[0]
        n_wavevectors = len(self.wavevectors)
        rhok_real = np.zeros(n_wavevectors)
        rhok_imag = np.zeros(n_wavevectors)
        for i in range(n_wavevectors):
            kvec = self.wavevectors[i]
            fac = np.sum(kvec * positions, axis=1)  # k·r for each particle
            rhok_real[i] = np.sum(np.cos(fac))
            rhok_imag[i] = np.sum(np.sin(fac))
        return rhok_real, rhok_imag

    def compute_fskt(self, rhok0_real, rhok0_imag, rhok_real, rhok_imag):
        n_wavevectors = len(self.wavevectors)
        correlations = np.zeros(n_wavevectors)
        for k in range(n_wavevectors):
            # ρₖ(t) · ρₖ*(t₀) = [real(t) + i·imag(t)] · [real(t₀) - i·imag(t₀)]
            # = real(t)·real(t₀) + imag(t)·imag(t₀) + i·[imag(t)·real(t₀) - real(t)·imag(t₀)]
            # We take the real part: real(t)·real(t₀) + imag(t)·imag(t₀)
            correlation = (rhok_real[k] * rhok0_real[k] + 
                          rhok_imag[k] * rhok0_imag[k])
            correlations[k] = correlation
        return np.mean(correlations)

    def _get_current_rhok(self, timestep):
        """Get rhok for current timestep with caching to avoid recomputation."""
        if timestep != self.cached_rhok_timestep:
            # Need to recompute
            with self.simulation.state.cpu_local_snapshot as snap:
                box_lengths = np.array([
                    snap.global_box.L[0],
                    snap.global_box.L[1],
                    snap.global_box.L[2]
                ])
                positions = unwrap_positions(
                    snap.particles.position,
                    snap.particles.image,
                    box_lengths
                )
                particle_mask = snap.particles.typeid != 2
                positions = positions[particle_mask]
            
            self.cached_rhok_real, self.cached_rhok_imag = self.compute_rhok(positions)
            self.cached_rhok_timestep = timestep
        
        return self.cached_rhok_real, self.cached_rhok_imag

    def _flush_buffers(self):
        """Flush all buffered data to files."""
        for ref in self.references:
            if ref['buffer']:
                try:
                    with open(ref['file_path'], 'a') as f:
                        for line in ref['buffer']:
                            f.write(line)
                    ref['buffer'].clear()
                except IOError as e:
                    print(f"Warning: Failed to flush F(k,t) buffer to {ref['file_path']}: {e}")

    def act(self, timestep):
        # Get current time
        current_time = self.time_tracker.elapsed_time
        
        # Get rhok for current frame (with caching)
        current_rhok_real, current_rhok_imag = self._get_current_rhok(timestep)

        # Add a new reference if enough time has passed and we haven't hit max_references
        if (len(self.references) < self.max_references and
            (self.last_reference_time is None or current_time - self.last_reference_time >= self.reference_interval_ps)):
            file_path = f'{self.output_prefix}_fskt_k{self.kmag:.2f}_ref{len(self.references)}.txt'
            
            # Initialize buffer for new reference
            buffer_data = []
            
            # Write header to file
            try:
                with open(file_path, 'w') as f:
                    f.write(f'# F(k,t) data for k={self.kmag:.4f}, t0={current_time:.6f} ps\n')
                    f.write('# lag_time(ps) F(k,t)\n')
            except IOError as e:
                print(f"Warning: Failed to create F(k,t) file {file_path}: {e}")
            
            self.references.append({
                'timestep': timestep,
                'time': current_time,
                'rhok_real': current_rhok_real.copy(),
                'rhok_imag': current_rhok_imag.copy(),
                'file_path': file_path,
                'buffer': buffer_data
            })
            self.last_reference_time = current_time

        # For each reference, compute and buffer F(k, t-t0)
        for ref in self.references:
            lag_time = current_time - ref['time']
            if lag_time < 0:  # Should not happen, but skip if so
                continue
            fskt_value = self.compute_fskt(ref['rhok_real'], ref['rhok_imag'], 
                                         current_rhok_real, current_rhok_imag)
            
            # Add to buffer instead of writing immediately
            ref['buffer'].append(f'{lag_time:.6f} {fskt_value:.6f}\n')
            
            # Flush buffer if it gets too large
            if len(ref['buffer']) >= self.file_buffer_size:
                try:
                    with open(ref['file_path'], 'a') as f:
                        for line in ref['buffer']:
                            f.write(line)
                    ref['buffer'].clear()
                except IOError as e:
                    print(f"Warning: Failed to write F(k,t) buffer to {ref['file_path']}: {e}")

        # Optionally, print debug info and flush buffers periodically
        if timestep - self.last_output_step >= self.output_period:
            print(f"F(k,t) computed for {len(self.references)} references at t={current_time:.4f} ps")
            # Flush all buffers periodically
            self._flush_buffers()
            self.last_output_step = timestep

    def __del__(self):
        """Ensure buffers are flushed when object is destroyed."""
        try:
            self._flush_buffers()
        except:
            pass  # Ignore errors during cleanup

    @hoomd.logging.log
    def current_fskt(self):
        # For logging, just return the value for the first reference (if any)
        if not self.references:
            return 0.0
        
        # Use cached rhok if available
        timestep = self.simulation.timestep
        current_rhok_real, current_rhok_imag = self._get_current_rhok(timestep)
        
        ref = self.references[0]
        return self.compute_fskt(ref['rhok_real'], ref['rhok_imag'], 
                               current_rhok_real, current_rhok_imag)

class AdaptiveTimestepUpdater(hoomd.custom.Action):
    def __init__(self, state, integrator, error_tolerance, time_constant_ps=50.0, initial_fraction=0.01, adaptiveerror=True, cavity_damping_factor=1.0, molecular_thermostat_tau=5.0, cavity_thermostat_tau=5.0):
        """
        Initialize the adaptive timestep updater.
        
        Args:
            state: Simulation state
            integrator: Integrator to update timestep for
            error_tolerance: Target error tolerance for timestepping
            time_constant_ps: Time constant for error tolerance approach (ps)
            initial_fraction: Initial fraction of target error tolerance
            adaptiveerror: Whether to adaptively change error tolerance
            cavity_damping_factor: Multiplier for cavity thermostat damping
            molecular_thermostat_tau: Time constant for molecular thermostat in ps
            cavity_thermostat_tau: Time constant for cavity thermostat in ps
        """
        super().__init__()
        print("Performing error tolerance ramping with time constant", time_constant_ps, "ps")
        self.state = state
        self.integrator = integrator
        self.target_error_tolerance = error_tolerance
        self.initial_error_tolerance = error_tolerance * initial_fraction
        self.current_error_tolerance = self.initial_error_tolerance
        self.time_constant_ps = time_constant_ps
        self.accumulated_time_ps = 0.0
        self.last_timestep = 0  # Will be set correctly in first act() call
        self.adaptiveerror = adaptiveerror
        self.cavity_damping_factor = cavity_damping_factor
        self.molecular_thermostat_tau = molecular_thermostat_tau
        self.cavity_thermostat_tau = cavity_thermostat_tau
        
    def act(self, timestep):
        """
        Custom action to update the timestep size.

        Parameters:
        - timestep: Current simulation timestep.
        """
        # Initialize last_timestep on first call to handle resuming from checkpoints
        if self.last_timestep == 0:
            self.last_timestep = timestep
        
        # Update accumulated simulation time
        if timestep > self.last_timestep:
            # Convert dt to picoseconds using correct conversion
            dt_ps = PhysicalConstants.atomic_units_to_ps(self.integrator.dt)
            self.accumulated_time_ps += (timestep - self.last_timestep) * dt_ps
        self.last_timestep = timestep
        
        # Update error tolerance based on exponential approach
        # formula: current = target - (target - initial) * exp(-t/tau)
        if self.adaptiveerror:
            exp_factor = np.exp(-self.accumulated_time_ps / self.time_constant_ps)
            self.current_error_tolerance = self.target_error_tolerance - \
                                          (self.target_error_tolerance - self.initial_error_tolerance) * exp_factor
        else:
            self.current_error_tolerance = self.target_error_tolerance
        
        # Collect forces and masses
        particle_data = self.state.get_snapshot().particles
        masses = np.array(particle_data.mass)
        n_particles = len(masses)
        
        # Initialize total forces array
        total_forces = np.zeros((n_particles, 3))
        
        # Sum forces from all force objects
        for force in self.integrator.forces:
            try:
                particle_forces = force.forces
                if particle_forces is not None and len(particle_forces) == n_particles:
                    total_forces += np.asarray(particle_forces)
            except (AttributeError, RuntimeError) as e:
                # Skip forces that can't be accessed (e.g., not yet computed)
                continue
        
        # Calculate sum |f_i| / m_i
        force_norm = np.array([np.linalg.norm(f) for f in total_forces])
        force_mass_sum = np.sum(force_norm / masses)
        
        # Update timestep using current error tolerance
        if force_mass_sum > 0:
            new_dt = np.sqrt(self.current_error_tolerance / force_mass_sum)
            
            self.integrator.dt = new_dt
            
            # Convert thermostat time constants from ps to atomic units
            molecular_tau_au = PhysicalConstants.ps_to_atomic_units(self.molecular_thermostat_tau)
            cavity_tau_au = PhysicalConstants.ps_to_atomic_units(self.cavity_thermostat_tau)
            
            # Update thermostat parameters for both molecular and cavity methods
            # methods[0] is always molecular method, methods[1] is cavity method (if present)
            
            # Update molecular method thermostat parameters
            molecular_method = self.integrator.methods[0]
            if hasattr(molecular_method, 'default_gamma'):  # Langevin thermostat
                molecular_method.default_gamma = PhysicalConstants.gamma_from_tau_ps(self.molecular_thermostat_tau)
            elif hasattr(molecular_method, 'thermostat'):  # Bussi or MTTK thermostat
                if hasattr(molecular_method.thermostat, 'tau'):
                    if hasattr(molecular_method.thermostat, 'translational_dof'):  # MTTK thermostat
                        molecular_method.thermostat.tau = molecular_tau_au  # Use configurable tau
                    else:  # Bussi thermostat
                        molecular_method.thermostat.tau = molecular_tau_au  # Use configurable tau
            
            # Update cavity method thermostat parameters if present
            if len(self.integrator.methods) > 1:  # Cavity particle present
                cavity_method = self.integrator.methods[1]
                if hasattr(cavity_method, 'default_gamma'):  # Langevin or Brownian thermostat
                    if hasattr(cavity_method, 'gamma'):  # Brownian dynamics
                        # Brownian dynamics has per-type gamma parameters - apply damping factor
                        base_gamma = PhysicalConstants.gamma_from_tau_ps(self.cavity_thermostat_tau)
                        cavity_method.default_gamma = self.cavity_damping_factor * base_gamma
                    else:  # Langevin thermostat
                        # Apply damping factor to base gamma calculated from tau
                        base_gamma = PhysicalConstants.gamma_from_tau_ps(self.cavity_thermostat_tau)
                        cavity_method.default_gamma = self.cavity_damping_factor * base_gamma
                elif hasattr(cavity_method, 'thermostat'):  # Bussi or MTTK thermostat
                    if hasattr(cavity_method.thermostat, 'tau'):
                        if hasattr(cavity_method.thermostat, 'translational_dof'):  # MTTK thermostat
                            cavity_method.thermostat.tau = cavity_tau_au  # Use configurable tau
                        else:  # Bussi thermostat
                            cavity_method.thermostat.tau = cavity_tau_au  # Use configurable tau
            
    @hoomd.logging.log
    def error_tolerance(self):
        """Log the current effective error tolerance."""
        return self.current_error_tolerance
        
    @hoomd.logging.log
    def elapsed_time_ps(self):
        """Log the elapsed simulation time in picoseconds."""
        return self.accumulated_time_ps

class CavityMDSimulation:
    """
    A class to encapsulate cavity MD simulation setup and execution.
    Handles cavity particle creation, thermostat configuration, and force setup.
    """
    
    def __init__(self, job_dir, replica, freq, couplstr, incavity, runtime_ps=500.0, 
                 input_gsd='molecular-0.gsd', frame=-1, name='prod', error_tolerance=0.01,
                 temperature=100.0, molecular_thermostat='bussi', cavity_thermostat='langevin',
                 cavity_damping_factor=1.0, use_brownian_overdamped=True, add_cavity_particle=True,
                 finite_q=False, molecular_thermostat_tau=5.0, cavity_thermostat_tau=5.0,
                 log_to_file=True, log_to_console=True, log_level='INFO', custom_log_file=None,
                 enable_fkt=True, fkt_kmag=1.0, fkt_num_wavevectors=50, fkt_reference_interval_ps=1.0, fkt_max_references=10,
                 max_energy_output_time_ps=None, enable_energy_tracking=True, dt_fs=None, device='CPU', gpu_id=0):
        """
        Initialize the CavityMDSimulation with simulation parameters.
        
        Args:
            job_dir: Directory for simulation input/output
            replica: Replica index for this simulation
            freq: Cavity frequency in cm^-1
            couplstr: Coupling strength
            incavity: Whether to include cavity coupling
            runtime_ps: Total simulation time in picoseconds
            input_gsd: Input GSD file to initialize the simulation from
            frame: Frame number to use from the input GSD file (-1 for last frame)
            name: Prefix for output files
            error_tolerance: Error tolerance for adaptive timestepping
            temperature: Target temperature in Kelvin
            molecular_thermostat: Type of thermostat for molecular system
            cavity_thermostat: Type of thermostat for cavity particle
            cavity_damping_factor: Multiplier for cavity thermostat damping
            use_brownian_overdamped: Use Brownian dynamics for highly overdamped cavity
            add_cavity_particle: Whether to add a cavity particle to the system
            finite_q: Allow finite-q photon displacement based on dipole-coupling interaction
            molecular_thermostat_tau: Time constant for molecular thermostat in ps
            cavity_thermostat_tau: Time constant for cavity thermostat in ps
            log_to_file: Whether to log output to a file (default: True)
            log_to_console: Whether to log output to console (default: True)
            log_level: Logging level ('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL')
            custom_log_file: Custom log file name (if None, auto-generated)
            enable_fkt: Whether to enable F(k,t) density correlation calculation (default: True)
            fkt_kmag: Wavevector magnitude for F(k,t) calculation (default: 1.0)
            fkt_num_wavevectors: Number of wavevectors to sample on sphere (default: 50)
            fkt_reference_interval_ps: Time interval between reference frames in ps (default: 1.0)
            fkt_max_references: Maximum number of reference frames to track (default: 10)
            max_energy_output_time_ps: Maximum time in ps to output energy data (None for no limit)
            enable_energy_tracking: Whether to enable energy tracking (default: True)
            dt_fs: Fixed timestep in femtoseconds (None for adaptive/default, only used when error_tolerance <= 0)
            device: Device type to use ('CPU' or 'GPU', default: 'CPU')
            gpu_id: GPU ID to use when device='GPU' (default: 0)
        """
        self.job_dir = job_dir
        self.replica = replica
        self.freq = freq
        self.couplstr = couplstr
        self.incavity = incavity
        self.runtime_ps = runtime_ps
        self.input_gsd = input_gsd
        self.frame = frame
        print(f"Frame: {frame}")
        self.name = name
        self.error_tolerance = error_tolerance
        self.temperature = temperature
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        self.cavity_damping_factor = cavity_damping_factor
        self.use_brownian_overdamped = use_brownian_overdamped
        self.add_cavity_particle = add_cavity_particle
        self.finite_q = finite_q
        self.molecular_thermostat_tau = molecular_thermostat_tau
        self.cavity_thermostat_tau = cavity_thermostat_tau
        
        # Logging parameters
        self.log_to_file = log_to_file
        self.log_to_console = log_to_console
        self.log_level = log_level
        self.custom_log_file = custom_log_file
        
        # F(k,t) calculation parameters
        self.enable_fkt = enable_fkt
        self.fkt_kmag = fkt_kmag
        self.fkt_num_wavevectors = fkt_num_wavevectors
        self.fkt_reference_interval_ps = fkt_reference_interval_ps
        self.fkt_max_references = fkt_max_references
        
        # Energy output limit parameter
        self.max_energy_output_time_ps = max_energy_output_time_ps
        
        # Energy tracking parameter
        self.enable_energy_tracking = enable_energy_tracking
        
        # Manual timestep parameter (in femtoseconds)
        self.dt_fs = dt_fs
        
        # Device configuration
        self.device = device.upper()
        self.gpu_id = gpu_id
        
        # Physical constants
        self.kB = PhysicalConstants.KB_HARTREE_PER_K
        
        # Initialize simulation components (will be set during setup)
        self.sim = None
        self.integrator = None
        self.forces = []
        self.methods = []
        self.trackers = {}
        self.logger = None
        
    def get_snapshot(self):
        """
        Get a snapshot from the simulation state, device-agnostic.
        
        Returns:
            Snapshot context manager that works with both CPU and GPU
        """
        if self.device == 'GPU':
            return self.sim.state.cpu_local_snapshot
        else:
            return self.sim.state.cpu_local_snapshot
    
    def setup_logging(self):
        """
        Set up logging configuration for the simulation.
        This method configures both file and console logging based on the user's preferences.
        """
        # Create a custom logger for this simulation
        logger_name = f"CavityMD_{self.name}_{self.replica}"
        self.logger = logging.getLogger(logger_name)
        self.logger.setLevel(getattr(logging, self.log_level.upper()))
        
        # Clear any existing handlers to avoid duplication
        self.logger.handlers.clear()
        
        # Create formatter for log messages
        formatter = logging.Formatter(
            '%(asctime)s | %(levelname)s | %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        # Set up file logging if requested
        if self.log_to_file:
            if self.custom_log_file:
                log_filename = self.custom_log_file
            else:
                # Auto-generate log filename based on simulation parameters
                cavity_suffix = "_cavity" if self.incavity else "_no_cavity"
                thermostat_suffix = f"_{self.molecular_thermostat}_{self.cavity_thermostat}" if self.incavity else f"_{self.molecular_thermostat}"
                device_suffix = f"_{self.device.lower()}"
                if self.device == 'GPU':
                    device_suffix += f"{self.gpu_id}"
                log_filename = f"{self.name}-{self.replica}{cavity_suffix}{thermostat_suffix}_{self.temperature}K{device_suffix}.log"
            
            # Create file handler
            file_handler = logging.FileHandler(log_filename, mode='w')
            file_handler.setLevel(getattr(logging, self.log_level.upper()))
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)
            
            # Store log filename for reference
            self.log_filename = log_filename
        
        # Set up console logging if requested
        if self.log_to_console:
            console_handler = logging.StreamHandler(sys.stdout)
            console_handler.setLevel(getattr(logging, self.log_level.upper()))
            console_handler.setFormatter(formatter)
            self.logger.addHandler(console_handler)
        
        # If not logging to console, redirect stdout to capture print statements
        if not self.log_to_console and self.log_to_file:
            # Create a custom stream that writes to the logger
            class LoggerWriter:
                def __init__(self, logger, level):
                    self.logger = logger
                    self.level = level
                    self.linebuf = ''
                
                def write(self, buf):
                    for line in buf.rstrip().splitlines():
                        if line.strip():  # Only log non-empty lines
                            self.logger.log(self.level, line.rstrip())
                
                def flush(self):
                    pass
            
            # Redirect stdout and stderr to logger
            sys.stdout = LoggerWriter(self.logger, logging.INFO)
            sys.stderr = LoggerWriter(self.logger, logging.ERROR)
        
        # Log initial setup information
        self.log_info("="*60)
        self.log_info("CAVITY MD SIMULATION STARTED")
        self.log_info("="*60)
        self.log_info(f"Simulation: {self.name}-{self.replica}")
        self.log_info(f"Device: {self.device}" + (f" (GPU {self.gpu_id})" if self.device == 'GPU' else ""))
        self.log_info(f"Runtime: {self.runtime_ps} ps")
        self.log_info(f"Temperature: {self.temperature} K")
        self.log_info(f"Cavity coupling: {'Enabled' if self.incavity else 'Disabled'}")
        if self.incavity:
            self.log_info(f"  Frequency: {self.freq} cm^-1")
            self.log_info(f"  Coupling strength: {self.couplstr}")
            self.log_info(f"  Finite-q mode: {self.finite_q}")
        self.log_info(f"Molecular thermostat: {self.molecular_thermostat} (tau={self.molecular_thermostat_tau} ps)")
        if self.incavity:
            self.log_info(f"Cavity thermostat: {self.cavity_thermostat} (tau={self.cavity_thermostat_tau} ps)")
        self.log_info(f"F(k,t) calculation: {'Enabled' if self.enable_fkt else 'Disabled'}")
        if self.enable_fkt:
            self.log_info(f"  k magnitude: {self.fkt_kmag}")
            self.log_info(f"  Wavevectors: {self.fkt_num_wavevectors}")
        self.log_info(f"Energy tracking: {'Enabled' if self.enable_energy_tracking else 'Disabled'}")
        if self.log_to_file:
            self.log_info(f"Logging to file: {getattr(self, 'log_filename', 'Unknown')}")
        self.log_info("="*60)
    
    def log_info(self, message):
        """Log an info message."""
        if self.logger:
            self.logger.info(message)
        else:
            print(message)
    
    def log_warning(self, message):
        """Log a warning message."""
        if self.logger:
            self.logger.warning(message)
        else:
            print(f"WARNING: {message}")
    
    def log_error(self, message):
        """Log an error message."""
        if self.logger:
            self.logger.error(message)
        else:
            print(f"ERROR: {message}")
    
    def log_debug(self, message):
        """Log a debug message."""
        if self.logger:
            self.logger.debug(message)
        else:
            print(f"DEBUG: {message}")

    def calculate_physical_parameters(self):
        """Calculate physical parameters and unit conversions."""
        # Time stepping parameters - keep ps and atomic unit values separate
        dt_ps = 0.0001  # timestep in ps (1 fs - reasonable for MD simulations)
        runtime_real = self.runtime_ps  # runtime in ps
        
        period_ps = 0.1  # ps
        table_period_ps = 0.1  # ps
        
        # Calculate periods in timesteps using ps values
        period = int(period_ps / dt_ps)
        table_period = int(table_period_ps / dt_ps)
        
        if period < 1:
            period = 1
        if table_period < 1:
            table_period = 1
        
        print(f"Period: {period}, Table period: {table_period}")
        
        # Convert timestep from ps to atomic units using helper method
        dt_au = PhysicalConstants.ps_to_atomic_units(dt_ps)
        
        # Calculate total steps needed for the specified runtime
        # Note: This is only used for fixed timestep mode
        # For adaptive timestep, ElapsedTimeTracker handles runtime termination
        total_steps_needed = int(runtime_real / dt_ps)
        
        # Store converted values as instance variables for later use
        self.dt = dt_au  # timestep in atomic units (for HOOMD)
        self.dt_ps = dt_ps  # timestep in ps (for calculations)
        self.runtime = total_steps_needed  # total steps needed
        self.period = period
        self.table_period = table_period
        
        self.log_info(f"Period: {period}, Table period: {table_period}")
        self.log_info(f"Time conversions:")
        self.log_info(f"  Timestep: {dt_ps} ps = {dt_au:.6f} a.u.")
        self.log_info(f"  Runtime: {self.runtime_ps:.1f} ps = {total_steps_needed} steps")
        self.log_info(f"  Steps per ps: {1.0/dt_ps:.1f}")
        
        return dt_au, total_steps_needed, period, table_period
        
    def create_cavity_particle(self, snapshot):
        """
        Add a cavity particle to the simulation snapshot.
        
        Args:
            snapshot: GSD snapshot to modify
        
        Returns:
            Modified snapshot with cavity particle added
        """
        self.log_info("Adding cavity particle to system...")
        
        # Calculate dipole moment and photon position
        positions = unwrap_positions(snapshot.particles.position, snapshot.particles.image, 
                                   snapshot.configuration.box[:3])
        dipmom = np.einsum('i,ij->j', snapshot.particles.charge, positions)

        omegac = self.freq / PhysicalConstants.HARTREE_TO_CM_MINUS1
        
        if self.finite_q:
            # Allow finite-q photon displacement based on dipole-coupling interaction
            newpos = -dipmom * self.couplstr / omegac**2
            newpos[-1] = 0.0
            # Only add thermal fluctuations if coupling is non-zero
            if self.couplstr != 0.0:
                sigma = np.sqrt(self.kB * self.temperature / omegac**2)
                newpos = np.random.normal(loc=newpos, scale=sigma, size=3)
                self.log_info(f"Finite-q mode: Photon displaced by dipole interaction to {newpos} (with thermal fluctuations)")
            else:
                self.log_info(f"Finite-q mode: Photon at equilibrium position {newpos} (no thermal fluctuations due to zero coupling)")
        else:
            # Start photon at origin (q=0 limit)
            newpos = np.array([0.0, 0.0, 0.0])
            # Only add thermal fluctuations if coupling is non-zero
            if self.couplstr != 0.0:
                sigma = np.sqrt(self.kB * self.temperature / omegac**2)
                newpos = np.random.normal(loc=newpos, scale=sigma, size=3)
                self.log_info("q=0 mode: Photon positioned at origin + thermal fluctuations")
            else:
                self.log_info("q=0 mode: Photon positioned exactly at origin (no thermal fluctuations due to zero coupling)")
        
        # Wrap position and get image flags
        def wrap_position(x, L):
            # Compute the image flags (how many box lengths away from the primary box)
            image_flags = np.floor((x + L/2) / L)
            # Compute the wrapped position inside the primary box
            wrapped_position = x - image_flags * L
            return wrapped_position, image_flags.astype(int)
            
        newpos, image_flags = wrap_position(newpos, np.array(snapshot.configuration.box[:3]))
        
        
        # Add photon particle
        if 'L' not in snapshot.particles.types:
            snapshot.particles.types.append('L')
        snapshot.particles.N += 1
        snapshot.particles.typeid = np.append(snapshot.particles.typeid, [2])
        snapshot.particles.position = np.append(snapshot.particles.position, [newpos], axis=0)
        snapshot.particles.charge = np.append(snapshot.particles.charge, [0.0])
        snapshot.particles.mass = np.append(snapshot.particles.mass, [1.0])
        snapshot.particles.diameter = np.append(snapshot.particles.diameter, [1.0])
        snapshot.particles.image = np.vstack([snapshot.particles.image, image_flags])

        # Set additional particle properties for the photon
        if hasattr(snapshot.particles, "body"):
            snapshot.particles.body = np.append(snapshot.particles.body, [-1], axis=0)

        if hasattr(snapshot.particles, "orientation"):
            snapshot.particles.orientation = np.append(
                snapshot.particles.orientation,
                [[1.0, 0.0, 0.0, 0.0]],
                axis=0
            )

        if hasattr(snapshot.particles, "moment_inertia"):
            snapshot.particles.moment_inertia = np.vstack([
                snapshot.particles.moment_inertia,
                np.zeros((1, 3))
            ])

        if hasattr(snapshot.particles, "velocity"):
            snapshot.particles.velocity = np.vstack([
                snapshot.particles.velocity,
                np.zeros((1, 3))
            ])

        if hasattr(snapshot.particles, "angmom"):
            snapshot.particles.angmom = np.vstack([
                snapshot.particles.angmom,
                np.zeros((1, 4))
            ])
            
        self.log_info(f"Cavity particle added at position {newpos}")
        return snapshot
        
    def setup_force_parameters(self, dt, rcut=15):
        """
        Set up force parameters for the simulation.
        
        Args:
            dt: Timestep in atomic units
            rcut: Cutoff radius in Bohr units
            
        Returns:
            List of configured forces
        """
        forces = []
        
        # Setup cavity force if requested
        if self.incavity:
            # Calculate omegac here for consistency with cavity particle positioning
            omegac = self.freq / PhysicalConstants.HARTREE_TO_CM_MINUS1
            cavityforce = CavityForce(kvector=np.array([0,0,1]), couplstr=self.couplstr, omegac=omegac)
            forces.append(cavityforce)
            self.trackers['cavityforce'] = cavityforce
        
        # Setup harmonic bonds
        harmonic = hoomd.md.bond.Harmonic()
        harmonic.params['O-O'] = dict(k=2*0.36602, r0=2.281655158)
        harmonic.params['N-N'] = dict(k=2*0.71625, r0=2.0743522177)
        forces.append(harmonic)
        self.trackers['harmonic'] = harmonic
        
        # Setup neighbor list
        cell = hoomd.md.nlist.Cell(buffer=1.0, exclusions=('bond',))
        
        # Setup Lennard-Jones interactions
        lj = hoomd.md.pair.LJ(nlist=cell, mode='shift')
        lj.params[('O', 'O')] = dict(epsilon=0.00016685201, sigma=6.230426584)
        lj.r_cut[('O', 'O')] = rcut
        lj.params[('N', 'N')] = dict(epsilon=0.000083426, sigma=5.48277488)
        lj.r_cut[('N', 'N')] = rcut
        lj.params[('N', 'O')] = dict(epsilon=0.00025027802, sigma=4.9832074319)
        lj.r_cut[('N', 'O')] = rcut

        # Disable pair interaction with 'L' particle (photon)
        if self.incavity:
            lj.params[('L', 'N')] = dict(epsilon=0.0, sigma=1.0)
            lj.r_cut[('L', 'N')] = 0.0
            lj.params[('N', 'L')] = dict(epsilon=0.0, sigma=1.0)
            lj.r_cut[('N', 'L')] = 0.0
            lj.params[('O', 'L')] = dict(epsilon=0.0, sigma=1.0)
            lj.r_cut[('O', 'L')] = 0.0
            lj.params[('L', 'O')] = dict(epsilon=0.0, sigma=1.0)
            lj.r_cut[('L', 'O')] = 0.0
            lj.params[('L', 'L')] = dict(epsilon=0.0, sigma=1.0)
            lj.r_cut[('L', 'L')] = 0.0
        forces.append(lj)
        self.trackers['lj'] = lj

        # Setup long-range Coulomb interactions using PPPM method
        numpoints = 32
        order = 6
        short, long = hoomd.md.long_range.pppm.make_pppm_coulomb_forces(
            nlist=cell, resolution=[numpoints, numpoints, numpoints], 
            order=order, r_cut=rcut, alpha=0.0
        )
        forces.append(short)
        forces.append(long)
        self.trackers['short'] = short
        self.trackers['long'] = long
        
        return forces
        
    def setup_thermostat_parameters(self, dt):
        """
        Set up thermostat parameters for molecular and cavity systems.
        
        Args:
            dt: Timestep in atomic units
            
        Returns:
            Tuple of (molecular_method, cavity_method, thermostat_references)
        """
        kT = self.kB * self.temperature
        molecular_filter = hoomd.filter.Type(['O', 'N'])  # Molecular particles only
        
        # Convert thermostat time constants from ps to atomic units using helper methods
        molecular_tau_au = PhysicalConstants.ps_to_atomic_units(self.molecular_thermostat_tau)
        cavity_tau_au = PhysicalConstants.ps_to_atomic_units(self.cavity_thermostat_tau)
        
        print(f"Thermostat time constant conversions:")
        print(f"  Molecular: {self.molecular_thermostat_tau:.3f} ps = {molecular_tau_au:.6f} a.u.")
        print(f"  Cavity: {self.cavity_thermostat_tau:.3f} ps = {cavity_tau_au:.6f} a.u.")
        
        # Validate tau=0.0 with Langevin thermostats
        if self.molecular_thermostat.lower() == 'langevin' and self.molecular_thermostat_tau <= 0.0:
            raise ValueError(
                f"ERROR: Cannot use Langevin thermostat with molecular_thermostat_tau={self.molecular_thermostat_tau} ps.\n"
                f"Langevin dynamics requires tau > 0 since gamma = 1/tau.\n"
                f"For overdamped dynamics (tau → 0), use:\n"
                f"  - molecular_thermostat='brownian' (Brownian dynamics), or\n"
                f"  - molecular_thermostat='none' (no molecular thermostat), or\n"
                f"  - A small positive tau value (e.g., 0.01 ps), or\n"
                f"  - Consider if you need thermostatting for the molecular system at all."
            )
        
        if self.incavity and self.cavity_thermostat.lower() == 'langevin' and self.cavity_thermostat_tau <= 0.0:
            raise ValueError(
                f"ERROR: Cannot use Langevin thermostat with cavity_thermostat_tau={self.cavity_thermostat_tau} ps.\n"
                f"Langevin dynamics requires tau > 0 since gamma = 1/tau.\n"
                f"For overdamped dynamics (tau → 0), use:\n"
                f"  - cavity_thermostat='brownian' (Brownian dynamics), or\n"
                f"  - cavity_thermostat='none' (no cavity thermostat), or\n"
                f"  - A small positive tau value (e.g., 0.01 ps)\n"
                f"Brownian dynamics is the correct choice for highly overdamped cavity modes."
            )
        
        # Store references to thermostats for energy tracking
        thermostat_refs = {
            'molecular_langevin': None,
            'cavity_langevin': None,
            'molecular_mttk': None,
            'cavity_mttk': None,
            'molecular_bussi': None,
            'cavity_bussi': None
        }
        
        # Configure molecular thermostat
        if self.molecular_thermostat.lower() == 'bussi':
            print("Running molecular system with Bussi thermostat (NVT ensemble)")
            molecular_bussi = Bussi(kT=kT, tau=molecular_tau_au)
            molecular_method = hoomd.md.methods.ConstantVolume(filter=molecular_filter, thermostat=molecular_bussi)
            thermostat_refs['molecular_bussi'] = molecular_bussi
            print(f"Molecular Bussi thermostat configured: T = {self.temperature:.1f} K, kT = {kT:.6f} a.u., tau = {self.molecular_thermostat_tau:.3f} ps")
            print("  Reservoir energy tracking enabled for molecular system")
        elif self.molecular_thermostat.lower() == 'langevin':
            print("Running molecular system with Langevin thermostat (NVT ensemble)")
            molecular_gamma = PhysicalConstants.gamma_from_tau_ps(self.molecular_thermostat_tau)
            molecular_method = hoomd.md.methods.Langevin(filter=molecular_filter, kT=kT, default_gamma=molecular_gamma, tally_reservoir_energy=True)
            thermostat_refs['molecular_langevin'] = molecular_method
            print(f"Molecular Langevin thermostat configured: T = {self.temperature:.1f} K, kT = {kT:.6f} a.u.")
            print(f"  gamma = {molecular_gamma:.6f} a.u.^-1 (tau = {self.molecular_thermostat_tau:.3f} ps)")
        elif self.molecular_thermostat.lower() == 'brownian':
            print("Running molecular system with Brownian dynamics (overdamped limit)")
            molecular_gamma = PhysicalConstants.gamma_from_tau_ps(self.molecular_thermostat_tau)
            molecular_method = hoomd.md.methods.Brownian(filter=molecular_filter, kT=kT, default_gamma=molecular_gamma, tally_reservoir_energy=True)
            print(f"Molecular Brownian dynamics configured: T = {self.temperature:.1f} K, kT = {kT:.6f} a.u.")
            print(f"  gamma = {molecular_gamma:.6f} a.u.^-1 (tau = {self.molecular_thermostat_tau:.3f} ps)")
            print("  OVERDAMPED - using Brownian dynamics")
            print("  Note: Velocities are drawn from Maxwell-Boltzmann distribution each timestep")
        elif self.molecular_thermostat.lower() == 'mttk':
            print("Running molecular system with MTTK (Nosé-Hoover) thermostat (NVT ensemble)")
            mttk_thermostat = hoomd.md.methods.thermostats.MTTK(kT=kT, tau=molecular_tau_au)
            molecular_method = hoomd.md.methods.ConstantVolume(filter=molecular_filter, thermostat=mttk_thermostat)
            thermostat_refs['molecular_mttk'] = mttk_thermostat
            print(f"Molecular MTTK thermostat configured: T = {self.temperature:.1f} K, kT = {kT:.6f} a.u., tau = {self.molecular_thermostat_tau:.3f} ps")
        elif self.molecular_thermostat.lower() == 'none':
            print("Running molecular system without thermostat (NVE ensemble)")
            molecular_method = hoomd.md.methods.ConstantVolume(filter=molecular_filter)
            print("No molecular thermostat - microcanonical (NVE) ensemble")
        else:
            raise ValueError(f"Invalid molecular_thermostat option: {self.molecular_thermostat}. Must be 'bussi', 'langevin', 'mttk', 'brownian', or 'none'")
        
        # Set up thermostat for cavity particle if present
        cavity_method = None
        if self.incavity:
            cavity_filter = hoomd.filter.Type(['L'])  # Cavity particle only
            
            # Configure cavity thermostat based on cavity_thermostat parameter
            if self.cavity_thermostat.lower() == 'langevin':
                # Determine if we should use Brownian dynamics (overdamped limit)
                use_brownian = self.use_brownian_overdamped and self.cavity_damping_factor > 5.0
                
                if use_brownian:
                    print("Running cavity with Brownian dynamics (overdamped limit)")
                    # Calculate base gamma from tau, then apply damping factor
                    base_gamma = PhysicalConstants.gamma_from_tau_ps(self.cavity_thermostat_tau)
                    cavity_gamma_brownian = self.cavity_damping_factor * base_gamma
                    cavity_method = hoomd.md.methods.Brownian(filter=cavity_filter, kT=kT, default_gamma=cavity_gamma_brownian, tally_reservoir_energy=True)
                    print(f"Cavity Brownian dynamics configured: T = {self.temperature:.1f} K, kT = {kT:.6f} a.u.")
                    print(f"  base_gamma = {base_gamma:.6f} a.u.^-1 (tau = {self.cavity_thermostat_tau:.3f} ps)")
                    print(f"  effective_gamma = {cavity_gamma_brownian:.6f} a.u.^-1 (damping_factor = {self.cavity_damping_factor:.1f}x)")
                    print(f"  effective_tau = {PhysicalConstants.atomic_units_to_ps(1.0/cavity_gamma_brownian):.6f} ps")
                    print(f"  OVERDAMPED - using Brownian dynamics")
                    print("  Note: Velocities are drawn from Maxwell-Boltzmann distribution each timestep")
                else:
                    print("Running cavity with Langevin thermostat")
                    # Calculate base gamma from tau, then apply damping factor
                    base_gamma = PhysicalConstants.gamma_from_tau_ps(self.cavity_thermostat_tau)
                    cavity_gamma = self.cavity_damping_factor * base_gamma
                    cavity_method = hoomd.md.methods.Langevin(filter=cavity_filter, 
                                                             kT=kT, default_gamma=cavity_gamma, tally_reservoir_energy=True)
                    thermostat_refs['cavity_langevin'] = cavity_method
                    print(f"Cavity Langevin thermostat configured: T = {self.temperature:.1f} K, kT = {kT:.6f} a.u.")
                    print(f"  base_gamma = {base_gamma:.6f} a.u.^-1 (tau = {self.cavity_thermostat_tau:.3f} ps)")
                    print(f"  effective_gamma = {cavity_gamma:.6f} a.u.^-1 (damping_factor = {self.cavity_damping_factor:.1f}x)")
                    print(f"  effective_tau = {PhysicalConstants.atomic_units_to_ps(1.0/cavity_gamma):.6f} ps")
                    print(f"  {'(OVERDAMPED)' if self.cavity_damping_factor > 1.0 else '(normal)'}")
                    if self.cavity_damping_factor > 5.0:
                        print(f"  Note: Consider using Brownian dynamics (--use_brownian_overdamped) for damping_factor > 5.0")
                        
            elif self.cavity_thermostat.lower() == 'brownian':
                print("Running cavity with Brownian dynamics (explicit)")
                # Calculate base gamma from tau, then apply damping factor
                base_gamma = PhysicalConstants.gamma_from_tau_ps(self.cavity_thermostat_tau)
                cavity_gamma_brownian = self.cavity_damping_factor * base_gamma
                cavity_method = hoomd.md.methods.Brownian(filter=cavity_filter, kT=kT, default_gamma=cavity_gamma_brownian, tally_reservoir_energy=True)
                print(f"Cavity Brownian dynamics configured: T = {self.temperature:.1f} K, kT = {kT:.6f} a.u.")
                print(f"  base_gamma = {base_gamma:.6f} a.u.^-1 (tau = {self.cavity_thermostat_tau:.3f} ps)")
                print(f"  effective_gamma = {cavity_gamma_brownian:.6f} a.u.^-1 (damping_factor = {self.cavity_damping_factor:.1f}x)")
                print(f"  effective_tau = {PhysicalConstants.atomic_units_to_ps(1.0/cavity_gamma_brownian):.6f} ps")
                
            elif self.cavity_thermostat.lower() == 'bussi':
                print("Running cavity with Bussi thermostat")
                cavity_bussi = Bussi(kT=kT, tau=cavity_tau_au)
                cavity_method = hoomd.md.methods.ConstantVolume(filter=cavity_filter, thermostat=cavity_bussi)
                thermostat_refs['cavity_bussi'] = cavity_bussi
                print(f"Cavity Bussi thermostat configured: kT = {kT:.6f} a.u., tau = {self.cavity_thermostat_tau:.3f} ps")
                print("  Reservoir energy tracking enabled for cavity system")
                if self.cavity_damping_factor != 1.0:
                    print(f"  WARNING: cavity_damping_factor={self.cavity_damping_factor} has no effect on Bussi thermostat")
            elif self.cavity_thermostat.lower() == 'mttk':
                print("Running cavity with MTTK (Nosé-Hoover) thermostat")
                cavity_mttk_thermostat = hoomd.md.methods.thermostats.MTTK(kT=kT, tau=cavity_tau_au)
                cavity_method = hoomd.md.methods.ConstantVolume(filter=cavity_filter, thermostat=cavity_mttk_thermostat)
                thermostat_refs['cavity_mttk'] = cavity_mttk_thermostat
                print(f"Cavity MTTK thermostat configured: kT = {kT:.6f} a.u., tau = {self.cavity_thermostat_tau:.3f} ps")
                if self.cavity_damping_factor != 1.0:
                    print(f"  WARNING: cavity_damping_factor={self.cavity_damping_factor} has no effect on MTTK thermostat")
            elif self.cavity_thermostat.lower() == 'none':
                print("Running cavity without thermostat (NVE ensemble)")
                cavity_method = hoomd.md.methods.ConstantVolume(filter=cavity_filter)
                print("No cavity thermostat - microcanonical (NVE) ensemble")
                if self.cavity_damping_factor != 1.0:
                    print(f"  WARNING: cavity_damping_factor={self.cavity_damping_factor} has no effect when no thermostat is used")
            else:
                raise ValueError(f"Invalid cavity_thermostat option: {self.cavity_thermostat}. Must be 'langevin', 'brownian', 'bussi', 'mttk', or 'none'")
        
        # Store thermostat references in the class for later use
        self.trackers.update(thermostat_refs)
        
        return molecular_method, cavity_method, thermostat_refs

    def setup_simulation(self):
        """
        Create HOOMD simulation object and initialize state from GSD file.
        Handles cavity particle addition/validation and snapshot modification.
        
        Returns:
            Modified snapshot if cavity particle was added, None otherwise
        """
        import os
        import gsd.hoomd
        
        # Save current directory and change to job directory
        self.original_cwd = os.getcwd()
        os.chdir(self.job_dir)
        
        # Load GSD file
        with gsd.hoomd.open(self.input_gsd, 'r') as f:
            if self.frame < 0:
                self.frame = len(f) + self.frame  # Convert negative index to positive
                if self.frame < 0:  # Handle case where abs(frame) > len(f)
                    self.frame = 0
            snapshot = f[self.frame]
            
            if self.incavity:
                if self.add_cavity_particle:
                    # Add new cavity particle
                    print("Adding cavity particle to system...")
                    snapshot = self.create_cavity_particle(snapshot)
                else:
                    # Validate that cavity particle already exists
                    print("Validating that cavity particle already exists in GSD file...")
                    
                    if 'L' not in snapshot.particles.types:
                        raise ValueError("ERROR: Cavity simulation requested but no cavity particle type 'L' found in GSD file. "
                                       "Make sure to run equilibration phase first with --add_cavity_particle flag.")
                    
                    if 2 not in snapshot.particles.typeid:
                        raise ValueError("ERROR: Cavity simulation requested but no cavity particles found in GSD file. "
                                       "Make sure to run equilibration phase first with --add_cavity_particle flag.")
                    
                    cavity_count = np.sum(snapshot.particles.typeid == 2)
                    if cavity_count != 1:
                        raise ValueError(f"ERROR: Expected exactly 1 cavity particle but found {cavity_count} in GSD file.")
                    
                    cavity_index = np.where(snapshot.particles.typeid == 2)[0][0]
                    cavity_position = snapshot.particles.position[cavity_index]
                    print(f"Cavity particle validated at index {cavity_index}, position {cavity_position}")
                    print("Continuing with existing cavity particle...")
        
        # Create simulation object with appropriate device
        if self.device == 'GPU':
            try:
                # Try different GPU initialization methods for different HOOMD versions
                try:
                    # First try with gpu_ids parameter (newer HOOMD versions)
                    device = hoomd.device.GPU(gpu_ids=[self.gpu_id])
                except TypeError:
                    # Fall back to older parameter name
                    device = hoomd.device.GPU(gpu_id=self.gpu_id)
                self.log_info(f"Initializing simulation on GPU {self.gpu_id}")
            except Exception as e:
                self.log_warning(f"Failed to initialize GPU {self.gpu_id}: {str(e)}")
                self.log_warning("Falling back to CPU")
                device = hoomd.device.CPU()
                self.device = 'CPU'  # Update device setting
        elif self.device == 'CPU':
            device = hoomd.device.CPU()
            self.log_info("Initializing simulation on CPU")
        else:
            raise ValueError(f"Invalid device '{self.device}'. Must be 'CPU' or 'GPU'")
        
        self.sim = hoomd.Simulation(device=device, seed=np.random.randint(int(10**4)))
        
        # Create state from snapshot or GSD file
        if self.incavity and self.add_cavity_particle:
            # Use modified snapshot
            self.sim.create_state_from_snapshot(snapshot)
            print(f"Simulation state created from modified snapshot with cavity particle")
            
            # Verify cavity particle position after state creation
            with self.sim.state.cpu_local_snapshot as snap:
                cavity_indices = np.where(snap.particles.typeid == 2)[0]
                if len(cavity_indices) > 0:
                    cavity_pos = snap.particles.position[cavity_indices[0]]
                    print(f"VERIFICATION: Cavity particle position in simulation state: {cavity_pos}")
                else:
                    print("WARNING: No cavity particle found in simulation state!")
            return snapshot
        else:
            # Use original GSD file
            self.sim.create_state_from_gsd(filename=self.input_gsd, frame=self.frame)
            print(f"Simulation state created from original GSD file frame {self.frame}")
            
            # Check if cavity particle exists in loaded state
            if self.incavity:
                with self.sim.state.cpu_local_snapshot as snap:
                    cavity_indices = np.where(snap.particles.typeid == 2)[0]
                    if len(cavity_indices) > 0:
                        cavity_pos = snap.particles.position[cavity_indices[0]]
                        print(f"VERIFICATION: Existing cavity particle position: {cavity_pos}")
                    else:
                        print("WARNING: No cavity particle found in loaded GSD file!")
            return None

    def setup_integrator(self, forces, methods):
        """
        Configure the integrator with forces and integration methods.
        
        Args:
            forces: List of force objects
            methods: List of integration methods [molecular_method, cavity_method]
        """
        # Setup integrator with initial dt
        integrator = hoomd.md.Integrator(dt=self.dt, forces=forces)
        self.sim.operations.integrator = integrator
        
        # Set integration methods (filter out None methods)
        valid_methods = [method for method in methods if method is not None]
        self.sim.operations.integrator.methods = valid_methods
        
        print(f"Integrator configured with initial dt = {self.dt:.6f} a.u. ({self.dt_ps:.6f} ps)")
        print(f"Number of integration methods: {len(valid_methods)}")

    def compute_and_set_optimal_timestep(self):
        """
        Compute and set the optimal timestep after running one step to initialize forces.
        This should be called after thermalization but before setting up trackers/loggers.
        """
        if self.error_tolerance <= 0:
            # Fixed timestep mode
            if self.dt_fs is not None:
                # Use user-specified timestep
                dt_au = PhysicalConstants.ps_to_atomic_units(self.dt_fs / 1000.0)  # Convert fs to ps, then to a.u.
                self.sim.operations.integrator.dt = dt_au
                self.dt = dt_au
                print(f"Using user-specified fixed timestep: {dt_au:.6f} a.u. ({self.dt_fs:.3f} fs)")
            else:
                # Keep current dt (HOOMD default)
                print(f"Using default fixed timestep: {self.sim.operations.integrator.dt:.6f} a.u. ({PhysicalConstants.atomic_units_to_ps(self.sim.operations.integrator.dt) * 1000:.3f} fs)")
            return
        
        try:
            print("Computing optimal timestep...")
            
            # Run one step to initialize forces (required by HOOMD)
            self.sim.run(1)
            
            # Use initial error tolerance (same as AdaptiveTimestepUpdater)
            initial_error_tolerance = self.error_tolerance * 1e-3  # initial_fraction = 1e-3
            
            # Collect forces and masses
            particle_data = self.sim.state.get_snapshot().particles
            masses = np.array(particle_data.mass)
            n_particles = len(masses)
            
            # Initialize total forces array
            total_forces = np.zeros((n_particles, 3))
            
            # Sum forces from all force objects
            for force in self.sim.operations.integrator.forces:
                try:
                    particle_forces = force.forces
                    if particle_forces is not None and len(particle_forces) == n_particles:
                        total_forces += np.asarray(particle_forces)
                except (AttributeError, RuntimeError) as e:
                    # Skip forces that can't be accessed (e.g., not yet computed)
                    continue
            
            # Calculate sum |f_i| / m_i (same logic as AdaptiveTimestepUpdater)
            force_norm = np.array([np.linalg.norm(f) for f in total_forces])
            force_mass_sum = np.sum(force_norm / masses)
            
            # Compute optimal timestep
            if force_mass_sum > 0:
                optimal_dt = np.sqrt(initial_error_tolerance / force_mass_sum)
                
                # Update integrator timestep
                self.sim.operations.integrator.dt = optimal_dt
                self.dt = optimal_dt  # Update stored dt
                
                print(f"Optimal timestep computed and set:")
                print(f"  Initial error tolerance: {initial_error_tolerance:.2e}")
                print(f"  Force/mass sum: {force_mass_sum:.6e}")
                print(f"  Optimal dt: {optimal_dt:.6f} a.u. ({PhysicalConstants.atomic_units_to_ps(optimal_dt) * 1000:.3f} fs)")
                
                # Reset simulation timestep counter to 0 so output starts clean
                # Note: We can't directly reset the timestep, but we can note this in output
                print(f"  Note: One initialization step was used for force computation")
            else:
                print("WARNING: Force/mass sum is zero - keeping initial timestep")
                
        except Exception as e:
            print(f"WARNING: Failed to compute optimal timestep: {str(e)}")
            print(f"Using initial timestep: {self.sim.operations.integrator.dt:.6f} a.u.")

    def thermalize_system(self):
        """
        Initialize particle velocities and thermostat degrees of freedom.
        Handles molecular and cavity particle thermalization separately.
        """
        kT = self.kB * self.temperature
        
        # Thermalize particle momenta
        if self.incavity:
            # Only thermalize molecular particles, not cavity particle
            molecular_filter = hoomd.filter.Type(['O', 'N'])
            self.sim.state.thermalize_particle_momenta(kT=kT, filter=molecular_filter)
            print("Thermalized molecular particles only (cavity particle excluded)")
            
            # Initialize cavity particle velocity based on thermostat type
            with self.get_snapshot() as snap:
                cavity_indices = np.where(snap.particles.typeid == 2)[0]
                if len(cavity_indices) > 0:
                    cavity_idx = cavity_indices[0]  # Get first cavity particle
                    
                    # For 3D Maxwell-Boltzmann: each component has variance kT/m
                    # With mass = 1.0 a.u., std dev per component = sqrt(kT)
                    cavity_velocity = np.random.normal(0.0, np.sqrt(kT), size=3)
                    
                    # Calculate expected kinetic energy and temperature
                    expected_ke = 0.5 * 1.0 * np.sum(cavity_velocity**2)  # KE = (1/2) * m * v²
                    expected_temp = (2.0/3.0) * expected_ke / self.kB  # T = (2/3) * KE / kB for 3D
                
                    print(f"Cavity particle thermalization:")
                    print(f"  Target temperature: {self.temperature:.1f} K")
                    print(f"  kT = {kT:.6f} a.u.")
                    print(f"  Initial velocity: {cavity_velocity}")
                    print(f"  Initial KE: {expected_ke:.6f} a.u.")
                    print(f"  Expected temperature: {expected_temp:.1f} K")
                    print(f"  Thermostat: {self.cavity_thermostat}")
                    
                    snap.particles.velocity[cavity_idx] = cavity_velocity
                    
                    # Verify the velocity was set correctly
                    set_velocity = snap.particles.velocity[cavity_idx]
                    print(f"  Verified velocity: {set_velocity}")
                    
                else:
                    print("WARNING: No cavity particle found for thermalization!")
        else:
            # No cavity particle, thermalize all particles
            self.sim.state.thermalize_particle_momenta(kT=kT, filter=hoomd.filter.All())
            print("Thermalized all molecular particles")
        
        # Initialize MTTK thermostat degrees of freedom if using MTTK
        molecular_mttk = self.trackers.get('molecular_mttk')
        cavity_mttk = self.trackers.get('cavity_mttk')
        
        if molecular_mttk is not None:
            print("Thermalizing molecular MTTK thermostat degrees of freedom...")
            self.sim.run(1)
            molecular_mttk.thermalize_dof()
            print("Molecular MTTK thermostat initialized")
        
        if cavity_mttk is not None:
            print("Thermalizing cavity MTTK thermostat degrees of freedom...")
            if molecular_mttk is None:  # Only run if we haven't already
                self.sim.run(1)
            cavity_mttk.thermalize_dof()
            print("Cavity MTTK thermostat initialized")
        
        # Initialize reservoir energy logging (requires at least one simulation step)
        if molecular_mttk is None and cavity_mttk is None:
            print("Initializing reservoir energy tracking...")
            self.sim.run(1)

    def setup_trackers_and_loggers(self):
        """
        Set up all tracking and logging objects for the simulation.
        Creates trackers for energy, time, cavity mode, density correlation, etc.
        """
        # Initialize time tracker
        self.time_tracker = ElapsedTimeTracker(self.sim, self.runtime_ps)
        
        # Setup status update logging
        self.status = Status(self.sim, PhysicalConstants.TIME_SECONDS * 1e9, self.time_tracker, self.runtime_ps)
        
        # Create timestep formatter
        self.timestep_formatter = TimestepFormatter(self.sim.operations.integrator)
        
        # Create console logger
        self.console_logger = hoomd.logging.Logger(categories=['scalar', 'string'])
        self.console_logger.add(self.sim, quantities=['timestep'])
        self.console_logger[('Status', 'nsd')] = (self.status, 'nsd', 'string')
        self.console_logger[('Status', 'etr')] = (self.status, 'etr', 'string')
        self.console_logger[('Status', 'elapsed')] = (self.status, 'elapsed', 'string')
        self.console_logger.add(self.time_tracker, quantities=['elapsed_time'])
        
        # Add timestep in femtoseconds to console logger
        self.console_logger[('Integrator', 'dt_fs')] = (self.timestep_formatter, 'dt_fs', 'scalar')
        
        # Set up kinetic energy tracking
        self.kinetic_energy_tracker = KineticEnergyTracker(
            simulation=self.sim,
            time_tracker=self.time_tracker,
            output_prefix=f'{self.name}-{self.replica}',
            output_period=1
        )
        
        # Add kinetic energy tracker as updater
        kinetic_energy_updater = hoomd.update.CustomUpdater(
            action=self.kinetic_energy_tracker,
            trigger=hoomd.trigger.Periodic(1)
        )
        self.sim.operations.updaters.append(kinetic_energy_updater)
        
        # Add to console logger
        self.console_logger[('MolecularTemp', 'temperature_K')] = (self.kinetic_energy_tracker, 'temperature', 'scalar')
        
        # Set up cavity mode tracking if in cavity simulation
        if self.incavity:
            cavityforce = self.trackers.get('cavityforce')
            self.cavity_mode_tracker = CavityModeTracker(
                simulation=self.sim,
                cavityforce=cavityforce,
                time_tracker=self.time_tracker,
                output_prefix=f'{self.name}-{self.replica}',
                output_period=1
            )
            
            # Add cavity mode tracker as updater
            cavity_mode_updater = hoomd.update.CustomUpdater(
                action=self.cavity_mode_tracker,
                trigger=hoomd.trigger.Periodic(1)
            )
            self.sim.operations.updaters.append(cavity_mode_updater)
            
            # Add to console logger
            self.console_logger[('CavityTemp', 'temperature_K')] = (self.cavity_mode_tracker, 'cavity_temperature', 'scalar')
        else:
            self.cavity_mode_tracker = None
        
        # Set up energy contribution tracker if enabled
        if self.enable_energy_tracking:
            self.energy_tracker = EnergyContributionTracker(
                simulation=self.sim,
                harmonic=self.trackers.get('harmonic'),
                lj=self.trackers.get('lj'),
                short=self.trackers.get('short'),
                long=self.trackers.get('long'),
                cavityforce=self.trackers.get('cavityforce'),
                cavity_langevin_method=self.trackers.get('cavity_langevin'),
                molecular_langevin_method=self.trackers.get('molecular_langevin'),
                mttk_thermostat=self.trackers.get('molecular_mttk'),
                cavity_mttk_thermostat=self.trackers.get('cavity_mttk'),
                molecular_bussi_thermostat=self.trackers.get('molecular_bussi'),
                cavity_bussi_thermostat=self.trackers.get('cavity_bussi'),
                kinetic_tracker=self.kinetic_energy_tracker,
                cavity_mode_tracker=self.cavity_mode_tracker,
                time_tracker=self.time_tracker,
                output_prefix=f'{self.name}-{self.replica}',
                output_period=1,
                max_time_ps=self.max_energy_output_time_ps
            )
            
            # Add energy tracker as updater
            energy_tracker_updater = hoomd.update.CustomUpdater(
                action=self.energy_tracker,
                trigger=hoomd.trigger.Periodic(1)
            )
            self.sim.operations.updaters.append(energy_tracker_updater)
        else:
            self.energy_tracker = None
            self.log_info("Energy tracking disabled")
        
        # Add reservoir energy tracking to console logger
        molecular_bussi = self.trackers.get('molecular_bussi')
        cavity_bussi = self.trackers.get('cavity_bussi')
        cavity_langevin = self.trackers.get('cavity_langevin')
        molecular_langevin = self.trackers.get('molecular_langevin')
        
        if molecular_bussi is not None:
            self.console_logger[('MolecularReservoir', 'total_energy')] = (molecular_bussi, 'total_reservoir_energy', 'scalar')
        if cavity_bussi is not None:
            self.console_logger[('CavityReservoir', 'total_energy')] = (cavity_bussi, 'total_reservoir_energy', 'scalar')
        if cavity_langevin is not None:
            self.console_logger[('CavityLangevin', 'reservoir_energy')] = (cavity_langevin, 'reservoir_energy', 'scalar')
        if molecular_langevin is not None:
            self.console_logger[('MolecularLangevin', 'reservoir_energy')] = (molecular_langevin, 'reservoir_energy', 'scalar')
        
        # Add time tracker as updater
        time_updater = hoomd.update.CustomUpdater(
            action=self.time_tracker,
            trigger=hoomd.trigger.Periodic(1)
        )
        self.sim.operations.updaters.append(time_updater)
        
        # Set up adaptive timestep updater if error_tolerance is positive
        if self.error_tolerance > 0:
            self.adaptive_action = AdaptiveTimestepUpdater(
                state=self.sim.state,
                integrator=self.sim.operations.integrator,
                error_tolerance=self.error_tolerance,
                time_constant_ps=10.0,
                initial_fraction=1e-3,
                adaptiveerror=True,
                cavity_damping_factor=self.cavity_damping_factor,
                molecular_thermostat_tau=self.molecular_thermostat_tau,
                cavity_thermostat_tau=self.cavity_thermostat_tau
            )
            
            # Add adaptive updater
            adaptive_updater = hoomd.update.CustomUpdater(
                action=self.adaptive_action,
                trigger=hoomd.trigger.Periodic(self.period)
            )
            self.sim.operations.updaters.append(adaptive_updater)
        
        # Add dipole autocorrelation updater
        self.dipole_autocorr = DipoleAutocorrelation(
            self.sim, 
            time_tracker=self.time_tracker, 
            output_prefix=f'{self.name}-{self.replica}_dipole'
        )
        dipole_autocorr_updater = hoomd.update.CustomUpdater(
            action=self.dipole_autocorr,
            trigger=hoomd.trigger.Periodic(1)
        )
        self.sim.operations.updaters.append(dipole_autocorr_updater)
        
        # Add F(k,t) density correlation tracker if enabled
        if self.enable_fkt:
            self.log_info(f"Setting up F(k,t) density correlation calculation:")
            self.log_info(f"  k magnitude: {self.fkt_kmag:.2f}")
            self.log_info(f"  Number of wavevectors: {self.fkt_num_wavevectors}")
            self.log_info(f"  Reference interval: {self.fkt_reference_interval_ps:.1f} ps")
            self.log_info(f"  Max references: {self.fkt_max_references}")
            
            # Calculate optimal buffer size and output period based on simulation runtime
            # For long simulations, use larger buffers and less frequent flushing
            runtime_scale = max(1.0, self.runtime_ps / 100.0)  # Scale factor based on runtime
            buffer_size = min(5000, int(1000 * runtime_scale))  # Between 1000 and 5000
            output_period = min(10000, max(1000, int(1000 * runtime_scale)))  # Between 1000 and 10000
            
            self.log_info(f"  Buffer size: {buffer_size}")
            self.log_info(f"  Output period: {output_period} steps")
            
            self.density_corr_tracker = DensityCorrelationTracker(
                simulation=self.sim,
                time_tracker=self.time_tracker,
                kmag=self.fkt_kmag,
                num_wavevectors=self.fkt_num_wavevectors,
                output_period=output_period,
                output_prefix=f'{self.name}-{self.replica}',
                reference_interval_ps=self.fkt_reference_interval_ps,
                max_references=self.fkt_max_references,
                file_buffer_size=buffer_size
            )
            
            # Use periodic trigger that's optimized for file I/O performance
            # Update every step for computation, but buffer writes for efficiency
            density_corr_updater = hoomd.update.CustomUpdater(
                action=self.density_corr_tracker,
                trigger=hoomd.trigger.Periodic(1)
            )
            self.sim.operations.updaters.append(density_corr_updater)
            
            # Add to console logger
            self.console_logger[('F(k,t)', 'current_fskt')] = (self.density_corr_tracker, 'current_fskt', 'scalar')
        else:
            self.density_corr_tracker = None
            self.log_info("F(k,t) density correlation calculation disabled")

    def setup_output_writers(self):
        """
        Configure GSD writer and console table for simulation output.
        Sets up logging for all tracked quantities.
        """
        # Create GSD logger
        self.gsd_logger = hoomd.logging.Logger(categories=['scalar'])
        self.gsd_logger.add(self.sim, quantities=['timestep', 'tps'])
        self.gsd_logger.add(self.time_tracker)
        self.gsd_logger[('integrator', 'dt')] = (self.sim.operations.integrator, 'dt', 'scalar')
        
        # Add energy tracking to GSD logger if enabled
        if self.enable_energy_tracking:
            self.gsd_logger[('EnergyTracker', 'total_potential_energy')] = (self.energy_tracker, 'total_potential_energy', 'scalar')
        self.gsd_logger[('KineticEnergy', 'kinetic_energy')] = (self.kinetic_energy_tracker, 'kinetic_energy', 'scalar')
        self.gsd_logger[('MolecularTemp', 'temperature_K')] = (self.kinetic_energy_tracker, 'temperature', 'scalar')
        
        # Add cavity mode tracking to GSD logger if applicable
        if self.incavity:
            self.gsd_logger[('CavityMode', 'cavity_kinetic_energy')] = (self.cavity_mode_tracker, 'cavity_kinetic_energy', 'scalar')
            self.gsd_logger[('CavityMode', 'cavity_potential_energy_harmonic')] = (self.cavity_mode_tracker, 'cavity_potential_energy_harmonic', 'scalar')
            self.gsd_logger[('CavityMode', 'cavity_total_energy')] = (self.cavity_mode_tracker, 'cavity_total_energy', 'scalar')
            self.gsd_logger[('CavityTemp', 'temperature_K')] = (self.cavity_mode_tracker, 'cavity_temperature', 'scalar')
        
        # Add F(k,t) density correlation tracking to GSD logger if enabled
        if self.enable_fkt and hasattr(self, 'density_corr_tracker') and self.density_corr_tracker is not None:
            self.gsd_logger[('F(k,t)', 'current_fskt')] = (self.density_corr_tracker, 'current_fskt', 'scalar')
        
        # Add adaptive timestep logging if enabled
        if hasattr(self, 'adaptive_action'):
            self.gsd_logger.add(self.adaptive_action, quantities=['error_tolerance', 'elapsed_time_ps'])
        
        # Set up GSD writer
        self.gsd_writer = hoomd.write.GSD(
            filename=f'{self.name}-{self.replica}.gsd',
            trigger=hoomd.trigger.Periodic(300000),
            dynamic=['property', 'momentum', 'particles/diameter', 'topology'],
            mode='wb',
            truncate=False,#True,  # Could be made configurable
            filter=hoomd.filter.All()
        )
        self.gsd_writer.logger = self.gsd_logger
        
        # Write initial frame
        self.gsd_writer.write(self.sim.state, filename=f'{self.name}-{self.replica}.gsd',
                             mode='wb', filter=hoomd.filter.All(), logger=self.gsd_logger)
        
        # Add GSD writer to simulation
        self.sim.operations.writers.append(self.gsd_writer)
        print(f"GSD writer added for file: {self.name}-{self.replica}.gsd")
        
        # Set up console output table
        self.table = hoomd.write.Table(
            trigger=hoomd.trigger.Periodic(period=self.table_period),
            logger=self.console_logger
        )
        self.sim.operations.writers.append(self.table)

    def run_simulation(self):
        """
        Execute the main simulation loop.
        For adaptive timestep, uses large step count and relies on ElapsedTimeTracker for runtime control.
        """
        # For adaptive timestep, use a very large number of steps
        # The ElapsedTimeTracker will exit when runtime is reached
        if self.error_tolerance > 0:
            # Adaptive timestep mode - use effectively infinite steps
            total_steps = 999999999  # Very large number - will be stopped by ElapsedTimeTracker
            self.log_info(f"Starting adaptive timestep simulation for {self.runtime_ps:.1f} ps")
            self.log_info(f"Using max steps: {total_steps} (will exit when runtime reached)")
        else:
            # Fixed timestep mode - use calculated steps
            total_steps = self.runtime
            self.log_info(f"Starting fixed timestep simulation for {self.runtime_ps:.1f} ps ({total_steps} steps)")
        
        # Get the actual timestep being used (may have been updated by adaptive timestep computation)
        actual_dt = self.sim.operations.integrator.dt
        actual_dt_ps = PhysicalConstants.atomic_units_to_ps(actual_dt)
        
        self.log_info(f"Initial timestep: {actual_dt:.6f} a.u. ({actual_dt_ps * 1000:.3f} fs)")
        if self.error_tolerance > 0:
            self.log_info(f"Timestep will adapt dynamically (error_tolerance = {self.error_tolerance})")
        else:
            self.log_info(f"Fixed timestep mode - steps per ps: {1.0/actual_dt_ps:.1f}")
        
        # Run the simulation
        self.sim.run(total_steps, write_at_start=True)
        
        self.log_info(f"Simulation completed successfully")

    def cleanup(self):
        """
        Handle post-simulation cleanup and restore original directory.
        """
        import os
        if hasattr(self, 'original_cwd'):
            os.chdir(self.original_cwd)
            print(f"Returned to original directory: {self.original_cwd}")

    def run(self):
        """
        Main orchestrator method that runs the complete simulation workflow.
        Coordinates all setup, execution, and cleanup phases.
        
        Returns:
            int: Exit code (0 for success)
        """
        try:
            # Phase 0: Setup logging
            self.setup_logging()
            
            # Phase 1: Setup simulation and state
            self.log_info("=== Phase 1: Setting up simulation ===")
            self.calculate_physical_parameters()
            snapshot = self.setup_simulation()
            
            # Phase 2: Configure forces and thermostats
            self.log_info("=== Phase 2: Configuring forces and thermostats ===")
            forces = self.setup_force_parameters(self.dt)
            molecular_method, cavity_method, thermostat_refs = self.setup_thermostat_parameters(self.dt)
            
            # Phase 3: Setup integrator and thermalization
            self.log_info("=== Phase 3: Setting up integrator and thermalization ===")
            methods = [molecular_method]
            if cavity_method is not None:
                methods.append(cavity_method)
            self.setup_integrator(forces, methods)
            self.thermalize_system()
            
            # Phase 3.5: Compute and set optimal timestep (after thermalization, before logging)
            self.log_info("=== Phase 3.5: Computing optimal timestep ===")
            self.compute_and_set_optimal_timestep()
            
            # Phase 4: Setup trackers and loggers
            self.log_info("=== Phase 4: Setting up trackers and loggers ===")
            self.setup_trackers_and_loggers()
            
            # Phase 5: Setup output writers
            self.log_info("=== Phase 5: Setting up output writers ===")
            self.setup_output_writers()
            
            # Phase 6: Run simulation
            self.log_info("=== Phase 6: Running simulation ===")
            self.run_simulation()
            
            # Phase 7: Cleanup
            self.log_info("=== Phase 7: Cleanup ===")
            self.cleanup()
            
            self.log_info("=== SIMULATION COMPLETED SUCCESSFULLY ===")
            return 0  # Success
            
        except Exception as e:
            self.log_error(f"CRITICAL ERROR in simulation: {str(e)}")
            import traceback
            self.log_error("Full traceback:")
            for line in traceback.format_exc().split('\n'):
                if line.strip():
                    self.log_error(line)
            self.cleanup()  # Try to cleanup even on error
            return 1  # Failure

class TimestepFormatter(hoomd.custom.Action):
    """Custom formatter to display timestep in femtoseconds."""
    def __init__(self, integrator):
        super().__init__()
        self.integrator = integrator
    
    def act(self, timestep):
        pass  # No action needed, just for logging
    
    
    @hoomd.logging.log
    def dt_fs(self):
        """Return timestep in femtoseconds."""
        dt_au = self.integrator.dt
        dt_fs = PhysicalConstants.atomic_units_to_ps(dt_au) * 1000.0  # Convert to fs
        return dt_fs


