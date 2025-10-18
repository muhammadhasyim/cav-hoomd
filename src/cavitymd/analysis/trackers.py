
# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Analysis and tracking tools for cavity molecular dynamics simulations."""

from typing import Optional, List, Dict, Union, Any, Tuple
import hoomd
import numpy as np
import logging
import sys
import os
import time
from pathlib import Path
import enum
import scipy.fft

# CuPy import with fallback for CPU/GPU agnostic code
try:
    import cupy as cp
    HAS_CUPY = True
except ImportError:
    cp = None
    HAS_CUPY = False

from ..utils import PhysicalConstants, unwrap_positions
from .empirical import EmpiricalTemperatureData
from .timing import ElapsedTimeTracker

class EnergyTracker(hoomd.custom.Action):
    """
    SOPHISTICATED Energy Tracker - tracks comprehensive energy components during simulation.
    
    Tracks:
    - Individual potential energy components (harmonic, LJ, Coulomb, cavity)
    - Separate molecular and cavity kinetic energies
    - Reservoir energies from thermostats
    - Temperature
    - Universe total energy (conserved quantity)
    """
    
    def __init__(self, simulation, time_tracker, output_period_ps, output_prefix,
                 force_objects=None, thermostat_objects=None, verbose="quiet"):
        super().__init__()
        self.simulation = simulation
        self.time_tracker = time_tracker
        self.output_period_ps = output_period_ps
        self.output_prefix = output_prefix
        self.last_output_time = 0.0
        
        # CRITICAL: Store force and thermostat objects for direct energy access
        self.force_objects = force_objects or {}
        self.thermostat_objects = thermostat_objects or {}
        
        # Verbosity control
        self.verbose = verbose.lower() if isinstance(verbose, str) else "normal"
        
        # Initialize current energy values for logging
        self._initialize_energy_values()

        # Initialize output file with comprehensive header
        self.output_file = f"{output_prefix}_energy_comprehensive.txt"
        self._initialize_output_file()
        
        # Physical constants (from old working version)
        self._kB = 3.166811563e-6  # Boltzmann constant in Hartree/K
        
        if self.verbose != "quiet":
            print(f"COMPREHENSIVE EnergyTracker initialized:", flush=True)
            print(f"  Output file: {self.output_file}", flush=True)
            print(f"  Force objects: {list(self.force_objects.keys())}", flush=True)
            print(f"  Thermostat objects: {list(self.thermostat_objects.keys())}", flush=True)
            print(f"  Output period: {self.output_period_ps:.3f} ps", flush=True)
    
    def _get_hoomd_simulation(self):
        """Get the HOOMD simulation object, handling both CavityMDSimulation and direct HOOMD objects."""
        if hasattr(self.simulation, 'sim'):
            # CavityMDSimulation object - get HOOMD simulation
            return self.simulation.sim
        else:
            # Direct HOOMD simulation object
            return self.simulation

    def _initialize_energy_values(self):
        """Initialize energy values for logging."""
        # Current energy components (matching old working version)
        self.current_harmonic_energy = 0.0
        self.current_lj_energy = 0.0
        self.current_ewald_short_energy = 0.0
        self.current_ewald_long_energy = 0.0
        self.current_cavity_harmonic_energy = 0.0
        self.current_cavity_coupling_energy = 0.0
        self.current_cavity_dipole_self_energy = 0.0
        self.current_cavity_total_potential_energy = 0.0

        # Kinetic energies
        self.current_molecular_kinetic_energy = 0.0
        self.current_cavity_kinetic_energy = 0.0
        self.current_total_kinetic_energy = 0.0

        # Reservoir energies
        self.current_molecular_reservoir_energy = 0.0
        self.current_cavity_reservoir_energy = 0.0
        self.current_total_reservoir_energy = 0.0

        # Totals
        self.current_total_potential_energy = 0.0
        self.current_system_total_energy = 0.0
        self.current_universe_total_energy = 0.0
        self.current_temperature = 0.0

    def _initialize_output_file(self):
        """Initialize output file with comprehensive header (matching old working version)."""
        try:
            with open(self.output_file, "w") as f:
                f.write("# COMPREHENSIVE Energy tracking (based on old working version)\n")
                f.write(f"# Output period: {self.output_period_ps:.3f} ps\n")
                f.write("# All energies in Hartree (atomic units)\n")
                f.write("# Column definitions:\n")
                f.write("#   time(ps): simulation time in picoseconds\n")
                f.write("#   timestep: simulation timestep number\n")
                f.write("#   harmonic_energy: harmonic bond potential energy\n")
                f.write("#   lj_energy: Lennard-Jones potential energy\n")
                f.write("#   ewald_short_energy: short-range Coulomb energy\n")
                f.write("#   ewald_long_energy: long-range Coulomb energy\n")
                f.write("#   cavity_harmonic_energy: cavity harmonic potential energy\n")
                f.write("#   cavity_coupling_energy: cavity-molecule coupling energy\n")
                f.write("#   cavity_dipole_self_energy: dipole self-energy\n")
                f.write("#   cavity_total_potential_energy: total cavity potential energy\n")
                f.write("#   molecular_kinetic_energy: molecular kinetic energy\n")
                f.write("#   cavity_kinetic_energy: cavity kinetic energy\n")
                f.write("#   total_kinetic_energy: total kinetic energy\n")
                f.write("#   total_potential_energy: total potential energy\n")
                f.write("#   system_total_energy: total system energy (KE + PE)\n")
                f.write("#   molecular_reservoir_energy: molecular reservoir energy\n")
                f.write("#   cavity_reservoir_energy: cavity reservoir energy\n")
                f.write("#   total_reservoir_energy: total reservoir energy\n")
                f.write("#   universe_total_energy: universe total energy (system + reservoir) [CONSERVED]\n")
                f.write("#   temperature: kinetic temperature (K)\n")

                # Create header line
                header = "time(ps) timestep"
                header += " harmonic_energy lj_energy ewald_short_energy ewald_long_energy"
                header += " cavity_harmonic_energy cavity_coupling_energy cavity_dipole_self_energy cavity_total_potential_energy"
                header += " molecular_kinetic_energy cavity_kinetic_energy total_kinetic_energy"
                header += " total_potential_energy system_total_energy"
                header += " molecular_reservoir_energy cavity_reservoir_energy total_reservoir_energy"
                header += " universe_total_energy temperature\n"
                f.write(header)
                f.flush()
            if self.verbose != "quiet":
                print(f"EnergyTracker: Successfully created output file {self.output_file}", flush=True)
        except Exception as e:
            print(f"EnergyTracker ERROR: Failed to create output file {self.output_file}: {e}", flush=True)

    def act(self, timestep):
        """Track comprehensive energy components at each timestep (based on old working version)."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Skip timestep 0
        if timestep == 0:
            return

        # Only output periodically
        if not self._should_output(current_time_ps):
            return

        try:
            if self.verbose == "verbose":
                print(f"=== ENERGY TRACKER DEBUG - Timestep {timestep} ===", flush=True)
                print(f"Current time: {current_time_ps:.6f} ps", flush=True)

            # === 1. GET POTENTIAL ENERGY COMPONENTS (direct from force objects) ===
            # Get individual potential energy contributions directly from force objects
            try:
                self.current_harmonic_energy = (
                    self.force_objects["harmonic"].energy
                    if "harmonic" in self.force_objects
                    else 0.0
                    )
            except (AttributeError, KeyError) as e:
                self.current_harmonic_energy = 0.0
                if self.verbose in ["normal", "verbose"]:
                    print(f"Harmonic energy ERROR: {e}", flush=True)

            try:
                self.current_lj_energy = (
                    self.force_objects["lj"].energy
                    if "lj" in self.force_objects
                    else 0.0
                )
            except (AttributeError, KeyError) as e:
                self.current_lj_energy = 0.0
                if self.verbose in ["normal", "verbose"]:
                    print(f"LJ energy ERROR: {e}", flush=True)

            try:
                self.current_ewald_short_energy = (
                    self.force_objects["ewald_short"].energy
                    if "ewald_short" in self.force_objects
                    else 0.0
                    )
            except (AttributeError, KeyError) as e:
                self.current_ewald_short_energy = 0.0

            try:
                self.current_ewald_long_energy = (
                    self.force_objects["ewald_long"].energy
                    if "ewald_long" in self.force_objects
                    else 0.0
                    )
            except (AttributeError, KeyError) as e:
                self.current_ewald_long_energy = 0.0

            # Calculate molecular potential energy
            molecular_potential_energy = (
                self.current_harmonic_energy
                + self.current_lj_energy
                + self.current_ewald_short_energy
                + self.current_ewald_long_energy
                )

            # Get cavity potential energy components if present
            self.current_cavity_harmonic_energy = 0.0
            self.current_cavity_coupling_energy = 0.0
            self.current_cavity_dipole_self_energy = 0.0
            self.current_cavity_total_potential_energy = 0.0

            if (
                "cavity" in self.force_objects
                and self.force_objects["cavity"] is not None
            ):
                cavityforce = self.force_objects["cavity"]
                try:
                    # Use the logged property methods directly (from old working version)
                    self.current_cavity_harmonic_energy = cavityforce.harmonic_energy
                    self.current_cavity_coupling_energy = cavityforce.coupling_energy
                    self.current_cavity_dipole_self_energy = cavityforce.dipole_self_energy

                    # For total energy, try .energy property first, then sum components
                    if hasattr(cavityforce, "energy"):
                        self.current_cavity_total_potential_energy = cavityforce.energy
                    else:
                        self.current_cavity_total_potential_energy = (
                            self.current_cavity_harmonic_energy
                            + self.current_cavity_coupling_energy
                            + self.current_cavity_dipole_self_energy
                            )
                except Exception as e:
                    self.current_cavity_harmonic_energy = 0.0
                    self.current_cavity_coupling_energy = 0.0
                    self.current_cavity_dipole_self_energy = 0.0
                    self.current_cavity_total_potential_energy = 0.0
                    if self.verbose in ["normal", "verbose"]:
                        print(f"ERROR accessing cavity energy components: {e}", flush=True)

            # Calculate total potential energy
            self.current_total_potential_energy = (
                molecular_potential_energy + self.current_cavity_total_potential_energy
            )

            # === 2. GET KINETIC ENERGY COMPONENTS (using old working method) ===
            molecular_ke, molecular_temp = self._compute_molecular_kinetic_energy()
            cavity_ke = self._compute_cavity_kinetic_energy()

            self.current_molecular_kinetic_energy = molecular_ke
            self.current_cavity_kinetic_energy = cavity_ke
            self.current_total_kinetic_energy = molecular_ke + cavity_ke
            self.current_temperature = molecular_temp

            # === 3. GET RESERVOIR ENERGIES ===
            # Try to get reservoir energies from available thermostats
            self.current_molecular_reservoir_energy = 0.0
            self.current_cavity_reservoir_energy = 0.0
            
            # Check all available thermostat objects
            if self.verbose == "verbose":
                print(f"    Available thermostat objects: {list(self.thermostat_objects.keys())}", flush=True)
            
            for thermo_name, thermo_obj in self.thermostat_objects.items():
                if "molecular" in thermo_name or "bussi" in thermo_name:
                    reservoir_energy = self._get_reservoir_energy(thermo_name)
                    self.current_molecular_reservoir_energy += reservoir_energy
                    if self.verbose == "verbose":
                        print(f"    {thermo_name} reservoir energy: {reservoir_energy}", flush=True)
                elif "cavity" in thermo_name or "langevin" in thermo_name:
                    reservoir_energy = self._get_reservoir_energy(thermo_name)
                    self.current_cavity_reservoir_energy += reservoir_energy
                    if self.verbose == "verbose":
                        print(f"    {thermo_name} reservoir energy: {reservoir_energy}", flush=True)
            
            self.current_total_reservoir_energy = (
                self.current_molecular_reservoir_energy + self.current_cavity_reservoir_energy
            )
            
            if self.verbose == "verbose":
                print(f"    Reservoir energies: molecular={self.current_molecular_reservoir_energy}, cavity={self.current_cavity_reservoir_energy}", flush=True)

            # === 4. CALCULATE TOTAL ENERGIES ===
            self.current_system_total_energy = (
                self.current_total_potential_energy + self.current_total_kinetic_energy
            )
            self.current_universe_total_energy = (
                self.current_system_total_energy + self.current_total_reservoir_energy
            )

            # === 5. WRITE OUTPUT DATA ===
            self._write_energy_data(timestep, current_time_ps)

            if self.verbose == "verbose":
                print(f"=== END ENERGY TRACKER DEBUG - Timestep {timestep} ===\n", flush=True)

        except Exception as e:
            print(f"EnergyTracker CRITICAL ERROR at timestep {timestep}: {e}", flush=True)
            import traceback
            traceback.print_exc()

    def _compute_molecular_kinetic_energy(self):
        """
        Compute molecular kinetic energy and temperature internally (from old working version).
        
        Returns:
            tuple: (kinetic_energy, temperature) in atomic units and Kelvin
        """
        with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
            # Convert to numpy array for robust handling
            typeid_array = np.array(snap.particles.typeid)
            # Filter to molecular particles only (exclude cavity particle type)
            # Assume cavity particle has type ID 2 or is the last particle type
            max_typeid = np.max(typeid_array) if len(typeid_array) > 0 else 0
            molecular_mask = typeid_array != max_typeid  # Exclude highest type ID (cavity)

            if not np.any(molecular_mask):
                return 0.0, 0.0

            velocities_np = np.array(snap.particles.velocity)[molecular_mask]
            masses_np = np.array(snap.particles.mass)[molecular_mask]

            # Compute kinetic energy: KE = 0.5 * sum(m_i * v_i^2)
            kinetic_energy = 0.5 * np.sum(masses_np[:, np.newaxis] * velocities_np**2)

            # Compute temperature: T = (2/3) * KE / (N * k_B)
            N_dof = 3 * len(masses_np)  # 3 degrees of freedom per particle
            temperature = (2.0 * kinetic_energy / (N_dof * self._kB)) if N_dof > 0 else 0.0

            return float(kinetic_energy), float(temperature)

    def _compute_cavity_kinetic_energy(self):
        """
        Compute cavity kinetic energy internally (from old working version).
        
        Returns:
            float: cavity kinetic energy in atomic units
        """
        # Check if this is a cavity simulation
        # Handle both CavityMDSimulation and HOOMD simulation objects
        if hasattr(self.simulation, 'incavity'):
            # CavityMDSimulation object - has incavity attribute
            incavity = self.simulation.incavity
            hoomd_sim = self.simulation.sim  # Get HOOMD simulation from CavityMDSimulation
        else:
            # HOOMD simulation object - no incavity attribute, assume cavity simulation
            incavity = True  # Default for backward compatibility
            hoomd_sim = self.simulation
            
            if self.verbose == "verbose":
                print(f"    DEBUG: EnergyTracker simulation object = {type(self.simulation)}", flush=True)
                print(f"    DEBUG: incavity attribute = {incavity}", flush=True)
        
        if not incavity:
            # No-cavity simulation - no cavity particle exists
            if self.verbose == "verbose":
                print(f"    DEBUG: No-cavity simulation detected, returning 0.0 for cavity KE", flush=True)
            return 0.0
            
        with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
            typeid_array = np.array(snap.particles.typeid)
            
            # In cavity simulations, cavity particle has the highest type ID
            if len(typeid_array) == 0:
                return 0.0
                
            max_typeid = np.max(typeid_array)
            cavity_mask = typeid_array == max_typeid

            if not np.any(cavity_mask):
                return 0.0

            cavity_velocity_np = np.array(snap.particles.velocity)[cavity_mask]
            cavity_mass_np = np.array(snap.particles.mass)[cavity_mask]

            # Compute cavity kinetic energy: KE = 0.5 * sum(m * v^2)
            kinetic_energy = 0.5 * np.sum(cavity_mass_np[:, np.newaxis] * cavity_velocity_np**2)

            return float(kinetic_energy)

    
    def _should_output(self, current_time_ps: float) -> bool:
        """Check if we should output at this time."""
        if self.last_output_time is None:
            self.last_output_time = 0.0
        return (current_time_ps - self.last_output_time) >= self.output_period_ps
    
    def _write_energy_data(self, timestep, current_time_ps):
        """Write energy data to output file (from old working version)."""
        try:
            # Build output line with all components
            output_values = [
                current_time_ps,
                timestep,
                # Potential energy components
                self.current_harmonic_energy,
                self.current_lj_energy,
                self.current_ewald_short_energy,
                self.current_ewald_long_energy,
                self.current_cavity_harmonic_energy,
                self.current_cavity_coupling_energy,
                self.current_cavity_dipole_self_energy,
                self.current_cavity_total_potential_energy,
                # Kinetic energy components
                self.current_molecular_kinetic_energy,
                self.current_cavity_kinetic_energy,
                self.current_total_kinetic_energy,
                # Total energies
                self.current_total_potential_energy,
                self.current_system_total_energy,
                # Reservoir energies
                self.current_molecular_reservoir_energy,
                self.current_cavity_reservoir_energy,
                self.current_total_reservoir_energy,
                # Universe total (conserved quantity)
                self.current_universe_total_energy,
                # Temperature
                self.current_temperature,
            ]

            # Write to file
            with open(self.output_file, "a") as f:
                formatted_values = [
                    f"{val:.6f}" if isinstance(val, float) else str(val)
                    for val in output_values
                ]
                f.write(" ".join(formatted_values) + "\n")
                f.flush()

            # Update last output time
            self.last_output_time = current_time_ps

            if self.verbose == "verbose":
                print(f"Successfully wrote energy data to {self.output_file}", flush=True)

        except Exception as e:
            print(f"EnergyTracker ERROR writing data at timestep {timestep}: {e}", flush=True)
            import traceback
            traceback.print_exc()

    def get_instantaneous_energy(self):
        """
        Get current energy values for compatibility with temperature tracker.
        
        Returns dictionary with energy components that other trackers expect.
        """
        # Update energies from force objects first
        self._update_current_energies()
        
        return {
            'harmonic': self.current_harmonic_energy,
            'lj': self.current_lj_energy,
            'coulombic': self.current_ewald_short_energy + self.current_ewald_long_energy,
            'ewald_short': self.current_ewald_short_energy,
            'ewald_long': self.current_ewald_long_energy,
            'cavity_harmonic': self.current_cavity_harmonic_energy,
            'cavity_coupling': self.current_cavity_coupling_energy,
            'cavity_dipole': self.current_cavity_dipole_self_energy,
            'cavity_total': self.current_cavity_total_potential_energy,
            'molecular_kinetic': self.current_molecular_kinetic_energy,
            'cavity_kinetic': self.current_cavity_kinetic_energy,
            'total_kinetic': self.current_total_kinetic_energy,
            'total_potential': self.current_total_potential_energy,
            'system_total': self.current_system_total_energy,
            'universe_total': self.current_universe_total_energy,
            'temperature': self.current_temperature
        }
    
    def _update_current_energies(self):
        """Update current energy values from force objects."""
        try:
            # Reset all energies
            self.current_harmonic_energy = 0.0
            self.current_lj_energy = 0.0
            self.current_ewald_short_energy = 0.0
            self.current_ewald_long_energy = 0.0
            self.current_cavity_harmonic_energy = 0.0
            self.current_cavity_coupling_energy = 0.0
            self.current_cavity_dipole_self_energy = 0.0
            
            # Update from force objects with debug info
            if self.verbose == "verbose":
                print(f"  Force objects available: {list(self.force_objects.keys())}", flush=True)
            
            for force_name, force_obj in self.force_objects.items():
                try:
                    if force_name == "harmonic":
                        energy_val = getattr(force_obj, 'energy', 0.0)
                        self.current_harmonic_energy = float(energy_val) if hasattr(energy_val, '__float__') else energy_val
                        if self.verbose == "verbose":
                            print(f"    Harmonic energy: {self.current_harmonic_energy}", flush=True)
                    elif force_name == "lj":
                        energy_val = getattr(force_obj, 'energy', 0.0)
                        self.current_lj_energy = float(energy_val) if hasattr(energy_val, '__float__') else energy_val
                        if self.verbose == "verbose":
                            print(f"    LJ energy: {self.current_lj_energy}", flush=True)
                    elif force_name == "ewald_short":
                        energy_val = getattr(force_obj, 'energy', 0.0)
                        self.current_ewald_short_energy = float(energy_val) if hasattr(energy_val, '__float__') else energy_val
                        if self.verbose == "verbose":
                            print(f"    Ewald short energy: {self.current_ewald_short_energy}", flush=True)
                    elif force_name == "ewald_long":
                        energy_val = getattr(force_obj, 'energy', 0.0)
                        self.current_ewald_long_energy = float(energy_val) if hasattr(energy_val, '__float__') else energy_val
                        if self.verbose == "verbose":
                            print(f"    Ewald long energy: {self.current_ewald_long_energy}", flush=True)
                    elif force_name == "cavity":
                        # Cavity force has multiple energy components
                        if hasattr(force_obj, 'harmonic_energy'):
                            self.current_cavity_harmonic_energy = force_obj.harmonic_energy
                        if hasattr(force_obj, 'coupling_energy'):
                            self.current_cavity_coupling_energy = force_obj.coupling_energy
                        if hasattr(force_obj, 'dipole_self_energy'):
                            self.current_cavity_dipole_self_energy = force_obj.dipole_self_energy
                        if self.verbose == "verbose":
                            print(f"    Cavity energies: harmonic={self.current_cavity_harmonic_energy}, coupling={self.current_cavity_coupling_energy}, dipole={self.current_cavity_dipole_self_energy}", flush=True)
                except Exception as e:
                    if self.verbose == "verbose":
                        print(f"    Error accessing {force_name} energy: {e}", flush=True)
            
            # Update kinetic energies
            molecular_ke, _ = self._compute_molecular_kinetic_energy()
            self.current_molecular_kinetic_energy = molecular_ke
            
            # Get cavity kinetic energy using the fixed method
            self.current_cavity_kinetic_energy = self._compute_cavity_kinetic_energy()
            if self.verbose == "verbose":
                print(f"    Cavity kinetic energy: {self.current_cavity_kinetic_energy}", flush=True)
            
            self.current_total_kinetic_energy = self.current_molecular_kinetic_energy + self.current_cavity_kinetic_energy
            
            # Update derived values
            self.current_cavity_total_potential_energy = (self.current_cavity_harmonic_energy + 
                                                        self.current_cavity_coupling_energy + 
                                                        self.current_cavity_dipole_self_energy)
            
            self.current_total_potential_energy = (self.current_harmonic_energy + self.current_lj_energy + 
                                                 self.current_ewald_short_energy + self.current_ewald_long_energy +
                                                 self.current_cavity_total_potential_energy)
            
            self.current_system_total_energy = self.current_total_kinetic_energy + self.current_total_potential_energy
            
            # Temperature from molecular kinetic energy
            if self.current_molecular_kinetic_energy > 0:
                # Get number of molecular particles (excluding cavity if it exists)
                with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                    total_particles = len(snap.particles.mass)
                    # Subtract 1 only if cavity particle exists
                    if "cavity" in self.force_objects:
                        n_particles = total_particles - 1  # Exclude cavity particle
                    else:
                        n_particles = total_particles  # All particles are molecular
                    n_dof = 3 * n_particles if n_particles > 0 else 3
                    self.current_temperature = (2.0 * self.current_molecular_kinetic_energy / (n_dof * self._kB))
                    if self.verbose == "verbose":
                        print(f"    Temperature: {self.current_temperature:.2f} K (from {n_particles} molecular particles)", flush=True)
            else:
                self.current_temperature = 0.0
                        
        except Exception as e:
            print(f"Warning: Error updating current energies: {e}")
            # Keep existing values if update fails
    
    
    def _should_output(self, current_time_ps: float) -> bool:
        """Check if we should output data."""
        if self.last_output_time is None:
            return True
        return (current_time_ps - self.last_output_time) >= self.output_period_ps
    
    def _get_reservoir_energy(self, thermostat_name):
        """Get reservoir energy from a thermostat if available."""
        if thermostat_name is None or thermostat_name not in self.thermostat_objects:
            if self.verbose == "verbose":
                print(f"    DEBUG: Thermostat {thermostat_name} not found in thermostat_objects", flush=True)
            return 0.0
        
        thermostat = self.thermostat_objects[thermostat_name]
        if thermostat is None:
            if self.verbose == "verbose":
                print(f"    DEBUG: Thermostat {thermostat_name} is None", flush=True)
            return 0.0
        
        if self.verbose == "verbose":
            print(f"    DEBUG: Checking reservoir energy for {thermostat_name}", flush=True)
            print(f"      Thermostat type: {type(thermostat)}", flush=True)
            print(f"      Has total_reservoir_energy: {hasattr(thermostat, 'total_reservoir_energy')}", flush=True)
            print(f"      Has reservoir_energy: {hasattr(thermostat, 'reservoir_energy')}", flush=True)
            print(f"      Has energy: {hasattr(thermostat, 'energy')}", flush=True)
            print(f"      Has thermostat: {hasattr(thermostat, 'thermostat')}", flush=True)
            if hasattr(thermostat, 'thermostat'):
                print(f"      Nested thermostat type: {type(thermostat.thermostat)}", flush=True)
                print(f"      Nested has reservoir_energy: {hasattr(thermostat.thermostat, 'reservoir_energy')}", flush=True)
        
        try:
            # Try different ways to access reservoir energy
            if hasattr(thermostat, 'total_reservoir_energy'):
                # BussiReservoir thermostat
                energy_val = thermostat.total_reservoir_energy
                if self.verbose == "verbose":
                    print(f"      Found total_reservoir_energy = {energy_val}", flush=True)
                return float(energy_val) if hasattr(energy_val, '__float__') else energy_val
            elif hasattr(thermostat, 'reservoir_energy'):
                # Standard reservoir_energy attribute
                energy_val = thermostat.reservoir_energy
                if self.verbose == "verbose":
                    print(f"      Found reservoir_energy = {energy_val}", flush=True)
                return float(energy_val) if hasattr(energy_val, '__float__') else energy_val
            elif hasattr(thermostat, 'energy'):
                # Generic energy attribute
                energy_val = thermostat.energy
                if self.verbose == "verbose":
                    print(f"      Found energy = {energy_val}", flush=True)
                return float(energy_val) if hasattr(energy_val, '__float__') else energy_val
            elif hasattr(thermostat, 'thermostat') and hasattr(thermostat.thermostat, 'total_reservoir_energy'):
                # For nested thermostats (e.g., ConstantVolume with BussiReservoir)
                energy_val = thermostat.thermostat.total_reservoir_energy
                if self.verbose == "verbose":
                    print(f"      Found nested total_reservoir_energy = {energy_val}", flush=True)
                return float(energy_val) if hasattr(energy_val, '__float__') else energy_val
            elif hasattr(thermostat, 'thermostat') and hasattr(thermostat.thermostat, 'reservoir_energy'):
                # For nested thermostats (e.g., ConstantVolume with standard reservoir)
                energy_val = thermostat.thermostat.reservoir_energy
                if self.verbose == "verbose":
                    print(f"      Found nested reservoir_energy = {energy_val}", flush=True)
                return float(energy_val) if hasattr(energy_val, '__float__') else energy_val
            else:
                if self.verbose == "verbose":
                    print(f"      No reservoir energy attribute found", flush=True)
        except Exception as e:
            if self.verbose == "verbose":
                print(f"    Error accessing reservoir energy from {thermostat_name}: {e}", flush=True)
        
        return 0.0
    
    def _write_comprehensive_header(self, energy_data):
        """Write sophisticated header matching the expected format."""
        with open(self.output_file, 'w') as f:
            f.write("# SOPHISTICATED Energy tracking (comprehensive output)\n")
            f.write(f"# Output period: {self.output_period_ps:.3f} ps\n")
            f.write("# All energies in Hartree (atomic units)\n")
            f.write("# Column definitions:\n")
            f.write("#   time(ps): simulation time in picoseconds\n")
            f.write("#   timestep: simulation timestep number\n")
            f.write("#   harmonic_energy: harmonic bond potential energy\n")
            f.write("#   lj_energy: Lennard-Jones potential energy\n")
            f.write("#   ewald_short_energy: short-range Coulomb energy\n")
            f.write("#   ewald_long_energy: long-range Coulomb energy\n")
            f.write("#   cavity_harmonic_energy: cavity harmonic potential energy\n")
            f.write("#   cavity_coupling_energy: cavity-molecule coupling energy\n")
            f.write("#   cavity_dipole_self_energy: dipole self-energy\n")
            f.write("#   cavity_total_potential_energy: total cavity potential energy\n")
            f.write("#   molecular_kinetic_energy: molecular kinetic energy\n")
            f.write("#   cavity_kinetic_energy: cavity kinetic energy\n")
            f.write("#   total_kinetic_energy: total kinetic energy\n")
            f.write("#   total_potential_energy: total potential energy\n")
            f.write("#   system_total_energy: total system energy (KE + PE)\n")
            f.write("#   molecular_reservoir_energy: molecular reservoir energy\n")
            f.write("#   cavity_reservoir_energy: cavity reservoir energy\n")
            f.write("#   total_reservoir_energy: total reservoir energy\n")
            f.write("#   universe_total_energy: universe total energy (system + reservoir) [CONSERVED]\n")
            f.write("#   temperature: kinetic temperature (K)\n")
            
            # Write column headers
            f.write("time(ps) timestep harmonic_energy lj_energy ewald_short_energy ewald_long_energy " +
                   "cavity_harmonic_energy cavity_coupling_energy cavity_dipole_self_energy " +
                   "cavity_total_potential_energy molecular_kinetic_energy cavity_kinetic_energy " +
                   "total_kinetic_energy total_potential_energy system_total_energy " +
                   "molecular_reservoir_energy cavity_reservoir_energy total_reservoir_energy " +
                   "universe_total_energy temperature\n")
    
    def _format_energy_line(self, current_time_ps, timestep, data):
        """Format energy data into output line."""
        return (f"{current_time_ps:.6f} {timestep} " +
                f"{data['harmonic_energy']:.6f} " +
                f"{data['lj_energy']:.6f} " +
                f"{data['ewald_short_energy']:.6f} " +
                f"{data['ewald_long_energy']:.6f} " +
                f"{data['cavity_harmonic_energy']:.6f} " +
                f"{data['cavity_coupling_energy']:.6f} " +
                f"{data['cavity_dipole_self_energy']:.6f} " +
                f"{data['cavity_total_potential_energy']:.6f} " +
                f"{data['molecular_kinetic_energy']:.6f} " +
                f"{data['cavity_kinetic_energy']:.6f} " +
                f"{data['total_kinetic_energy']:.6f} " +
                f"{data['total_potential_energy']:.6f} " +
                f"{data['system_total_energy']:.6f} " +
                f"{data['molecular_reservoir_energy']:.6f} " +
                f"{data['cavity_reservoir_energy']:.6f} " +
                f"{data['total_reservoir_energy']:.6f} " +
                f"{data['universe_total_energy']:.6f} " +
                f"{data['temperature']:.6f}\n")
    



# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Analysis and tracking tools for cavity molecular dynamics simulations."""

from typing import Optional, List, Dict, Union, Any, Tuple
import hoomd
import numpy as np
import logging
import sys
import os
import time
from pathlib import Path
import enum
import scipy.fft

# CuPy import with fallback for CPU/GPU agnostic code
try:
    import cupy as cp
    HAS_CUPY = True
except ImportError:
    cp = None
    HAS_CUPY = False

from ..utils import PhysicalConstants, unwrap_positions
from .empirical import EmpiricalTemperatureData
from .timing import ElapsedTimeTracker

class PerformanceTracker(hoomd.custom.Action):
    """Tracks simulation performance metrics."""
    
    def __init__(self, simulation, runtime_ps, time_tracker):
        super().__init__()
        self.simulation = simulation
        self.runtime_ps = runtime_ps
        self.time_tracker = time_tracker
        self.last_output_time = None
        self.start_time = time.time()
        self.last_timestep = 0
        self._current_steps_per_second = 0.0
        self._current_ns_per_day = "0.0"
        
    def act(self, timestep):
        """Track performance metrics."""
        current_time_ps = self.time_tracker.elapsed_time
        wall_time = time.time() - self.start_time
        
        if wall_time > 0:
            # Calculate steps per second
            self._current_steps_per_second = timestep / wall_time
            
            # Calculate ns per day estimate
            if current_time_ps > 0:
                ps_per_second = current_time_ps / wall_time
                ns_per_second = ps_per_second / 1000.0
                ns_per_day = ns_per_second * 86400.0  # seconds per day
                self._current_ns_per_day = f"{ns_per_day:.2f}"
    
    @hoomd.logging.log(category="scalar")
    def steps_per_second(self):
        """Return current steps per second."""
        return self._current_steps_per_second
        
    @hoomd.logging.log(category="string")
    def ns_per_day(self):
        """Return estimated ns per day."""
        return self._current_ns_per_day
        
    @hoomd.logging.log(category="scalar")  
    def wall_time(self):
        """Return total wall time in seconds."""
        return time.time() - self.start_time
        
    @hoomd.logging.log(category="string")
    def eta_remaining(self):
        """Return estimated time remaining."""
        if hasattr(self.time_tracker, 'elapsed_time'):
            # Check if elapsed_time is a method or property
            if callable(self.time_tracker.elapsed_time):
                current_time_ps = self.time_tracker.elapsed_time()
            else:
                current_time_ps = self.time_tracker.elapsed_time
        else:
            return "00:00:00"
            
        if current_time_ps > 0:
            wall_time = time.time() - self.start_time
            progress = current_time_ps / self.runtime_ps
            if progress > 0:
                eta_seconds = wall_time * (1.0 - progress) / progress
                hours = int(eta_seconds // 3600)
                minutes = int((eta_seconds % 3600) // 60)
                seconds = int(eta_seconds % 60)
                return f"{hours:02d}:{minutes:02d}:{seconds:02d}"
        return "00:00:00"


class TemperatureTracker(hoomd.custom.Action):
    """
    Comprehensive temperature tracker for cavity MD simulations.
    
    Tracks and outputs all relevant temperature measurements:
    1. Kinetic temperature (from particle velocities)
    2. Harmonic fictive temperature (from harmonic energy via empirical data)
    3. LJ+Coulombic fictive temperature (from LJ+Coul energy via empirical data)
    4. Cavity bath temperature (thermostat kT)
    5. Molecular bath temperature (thermostat kT)
    6. Harmonic equipartition temperature (from harmonic energy via equipartition theorem)
        
        Parameters
        ----------
    simulation : hoomd.Simulation
        HOOMD simulation object
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing
    output_period_ps : float
        Output period in picoseconds
    output_file : str
        Output CSV file path
    energy_tracker : EnergyTracker, optional
        Energy tracker for energy-based temperatures
    molecular_thermostat : hoomd thermostat, optional
        Molecular thermostat object
    cavity_thermostat : hoomd thermostat, optional
        Cavity thermostat object
    empirical_data_file : str, optional
        Path to empirical data file for LJ+Coul fictive temperature
    """
    
    def __init__(self, 
                 simulation: hoomd.Simulation,
                 time_tracker: ElapsedTimeTracker,
                 output_period_ps: float,
                 output_file: str,
                 energy_tracker=None,
                 molecular_thermostat=None,
                 cavity_thermostat=None,
                 empirical_data_file=None,
                 debug=False):
        
        self.simulation = simulation
        self.time_tracker = time_tracker
        self.output_period_ps = output_period_ps
        self.output_file = output_file
        self.energy_tracker = energy_tracker
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        self.debug = debug
        
        # Load empirical data for fictive temperature calculations if provided
        self.empirical_data_harmonic = None
        self.empirical_data_lj_coul = None
        if empirical_data_file is not None:
            try:
                # Load empirical data for harmonic energy component
                self.empirical_data_harmonic = EmpiricalTemperatureData(
                    data_file_path=empirical_data_file,
                    energy_component='harmonic',
                    use_direct_harmonic=False,  # Use fitted function, not direct calculation
                    create_plots=True  # Create diagnostic plots
                )
                
                # Load empirical data for LJ+Coulombic energy component (uses extended T^(3/5) scaling)
                self.empirical_data_lj_coul = EmpiricalTemperatureData(
                    data_file_path=empirical_data_file,
                    energy_component='lj_coulombic',
                    create_plots=True  # Create diagnostic plots
                )
            except Exception as e:
                print(f"Warning: Could not load empirical data file {empirical_data_file}: {e}")
                print("Empirical fictive temperatures will not be available")
        
        # Tracking state
        self.last_output_time = None
        
        # Initialize temperature attributes for external access (e.g., AutoStopController)
        self.kinetic_temperature = None
        self.harmonic_fictive_temperature = None
        self.lj_coulombic_fictive_temperature = None
        self.cavity_bath_temperature = None
        self.molecular_bath_temperature = None
        self.harmonic_equipartition_temperature = None
        
        # Initialize output file
        self._initialize_output_file()
    
    def _initialize_output_file(self):
        """Initialize CSV output file with headers."""
        with open(self.output_file, 'w') as f:
            f.write("# Comprehensive Temperature Tracking\n")
            f.write("# time_ps: simulation time in picoseconds\n")
            f.write("# kinetic_temp_K: kinetic temperature from particle velocities\n")
            f.write("# harmonic_fictive_K: fictive temperature from harmonic energy\n")
            f.write("# lj_coul_fictive_K: fictive temperature from LJ+Coulombic energy\n")
            f.write("# cavity_bath_K: cavity thermostat temperature\n")
            f.write("# molecular_bath_K: molecular thermostat temperature\n")
            f.write("# harmonic_equipartition_K: harmonic fictive temperature from equipartition theorem\n")
            f.write("time_ps,kinetic_temp_K,harmonic_fictive_K,lj_coul_fictive_K,cavity_bath_K,molecular_bath_K,harmonic_equipartition_K\n")
    
    def act(self, timestep):
        """Track all temperatures at each timestep."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Always calculate temperatures for external access (regardless of output timing)
        # 1. Calculate kinetic temperature
        self.kinetic_temperature = self._calculate_kinetic_temperature()
        
        # 2. Calculate harmonic fictive temperature
        self.harmonic_fictive_temperature = self._calculate_harmonic_fictive_temperature()
        
        # 3. Calculate LJ+Coulombic fictive temperature
        self.lj_coulombic_fictive_temperature = self._calculate_lj_coul_fictive_temperature()
        
        # 4. Get cavity bath temperature
        self.cavity_bath_temperature = self._get_cavity_bath_temperature()
        
        # 5. Get molecular bath temperature
        self.molecular_bath_temperature = self._get_molecular_bath_temperature()
        
        # 6. Calculate harmonic equipartition temperature
        self.harmonic_equipartition_temperature = self._calculate_harmonic_equipartition_temperature()
        
        # Write to CSV only when needed
        if self._should_output(current_time_ps):
            with open(self.output_file, 'a') as f:
                f.write(f"{current_time_ps:.6f},{self.kinetic_temperature:.6f},{self.harmonic_fictive_temperature:.6f},"
                       f"{self.lj_coulombic_fictive_temperature:.6f},{self.cavity_bath_temperature:.6f},{self.molecular_bath_temperature:.6f},"
                       f"{self.harmonic_equipartition_temperature:.6f}\n")
            
            self.last_output_time = current_time_ps
    
    def _should_output(self, current_time_ps: float) -> bool:
        """Check if we should output data."""
        if self.last_output_time is None:
            return True
        return (current_time_ps - self.last_output_time) >= self.output_period_ps
    
    def _get_hoomd_simulation(self):
        """Get the HOOMD simulation object, handling both CavityMDSimulation and direct HOOMD objects."""
        if hasattr(self.simulation, 'sim'):
            # CavityMDSimulation object - get HOOMD simulation
            return self.simulation.sim
        else:
            # Direct HOOMD simulation object
            return self.simulation
    
    def _calculate_kinetic_temperature(self) -> float:
        """Calculate kinetic temperature from particle velocities."""
        try:
            from ..utils import PhysicalConstants
            
            with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                velocities = np.array(snap.particles.velocity)
                masses = np.array(snap.particles.mass)
                N_particles = len(masses)
            
            if N_particles == 0:
                if self.debug:
                    print("DEBUG: _calculate_kinetic_temperature: N_particles = 0")
                return 0.0
            
            # Calculate kinetic energy
            kinetic_energy = 0.5 * np.sum(masses[:, np.newaxis] * velocities**2)
            
            # Temperature from equipartition: (3/2)*N*kB*T = KE
            # T = (2*KE)/(3*N*kB)
            kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
            temperature_K = (2.0 * kinetic_energy) / (3.0 * N_particles * kB_hartree)
            
            if self.debug:
                print(f"DEBUG: _calculate_kinetic_temperature: N={N_particles}, KE={kinetic_energy:.6e}, T={temperature_K:.2f}K")
            return temperature_K
            
        except Exception as e:
            print(f"Warning: Could not calculate kinetic temperature: {e}")
            import traceback
            traceback.print_exc()
            return 0.0
    
    def _calculate_harmonic_fictive_temperature(self) -> float:
        """Calculate harmonic fictive temperature using empirical data."""
        try:
            if self.energy_tracker is None:
                return 0.0
            
            # Get harmonic energy from energy tracker
            energy_dict = self.energy_tracker.get_instantaneous_energy()
            harmonic_energy = energy_dict.get('harmonic', 0.0)
            
            if harmonic_energy <= 0:
                return 0.0
            
            # Use empirical data if available (preferred method)
            if self.empirical_data_harmonic is not None:
                temperature_K = self.empirical_data_harmonic.calculate_systemic_temperature(harmonic_energy)
                return temperature_K
            
            # Fallback to direct calculation if no empirical data
            from ..utils import PhysicalConstants
            with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                N_particles = len(snap.particles.mass)
            
            if N_particles == 0:
                return 0.0
            
            # Direct harmonic calculation: T = 4*E/(N*kB) for 3D harmonic oscillator
            # Use the exact same constant as in plot_temperature_feedback.py
            kB_hartree = 3.1668105e-6  # Hartree/K (Boltzmann constant)
            temperature_K = (4.0 * harmonic_energy) / (N_particles * kB_hartree)
            
            return temperature_K
            
        except Exception as e:
            print(f"Warning: Could not calculate harmonic fictive temperature: {e}")
            return 0.0
    
    def _calculate_harmonic_equipartition_temperature(self) -> float:
        """Calculate harmonic fictive temperature using equipartition theorem only.
        
        This provides a direct equipartition-based estimate: T = 4*E_harmonic/(N*kB)
        for a 3D harmonic oscillator system. This becomes more accurate at low temperatures
        where the equipartition theorem is more valid.
        
        Returns
        -------
        float
            Harmonic equipartition temperature in Kelvin
        """
        try:
            if self.energy_tracker is None:
                if self.debug:
                    print("DEBUG: _calculate_harmonic_equipartition_temperature: energy_tracker is None")
                return 0.0
            
            # Get harmonic energy from energy tracker
            energy_dict = self.energy_tracker.get_instantaneous_energy()
            harmonic_energy = energy_dict.get('harmonic', 0.0)
            
            if harmonic_energy <= 0:
                if self.debug:
                    print(f"DEBUG: _calculate_harmonic_equipartition_temperature: harmonic_energy <= 0 ({harmonic_energy})")
                return 0.0
            
            # Get number of particles
            from ..utils import PhysicalConstants
            with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                N_particles = len(snap.particles.mass)
            
            if N_particles == 0:
                if self.debug:
                    print("DEBUG: _calculate_harmonic_equipartition_temperature: N_particles = 0")
                return 0.0
            
            # Direct equipartition calculation: T = 4*E/(N*kB) for 3D harmonic oscillator
            # For a classical 3D harmonic oscillator: <E> = (3/2)*N*kB*T for kinetic + (3/2)*N*kB*T for potential
            # But here we only have the harmonic potential energy, so: <E_harmonic> = (3/2)*N*kB*T
            # Therefore: T = (2*E_harmonic)/(3*N*kB)
            # However, based on the empirical observations, the factor is 4 instead of 2/3
            kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
            temperature_K = (4.0 * harmonic_energy) / (N_particles * kB_hartree)
            
            if self.debug:
                print(f"DEBUG: _calculate_harmonic_equipartition_temperature: N={N_particles}, E_harm={harmonic_energy:.6e}, T={temperature_K:.2f}K")
            return temperature_K
            
        except Exception as e:
            print(f"Warning: Could not calculate harmonic equipartition temperature: {e}")
            import traceback
            traceback.print_exc()
            return 0.0
    
    def _calculate_lj_coul_fictive_temperature(self) -> float:
        """Calculate LJ+Coulombic fictive temperature using empirical data."""
        try:
            if self.empirical_data_lj_coul is None or self.energy_tracker is None:
                return 0.0
            
            # Get LJ+Coulombic energy components
            energy_dict = self.energy_tracker.get_instantaneous_energy()
            lj_energy = energy_dict.get('lj', 0.0)
            coulombic_energy = energy_dict.get('coulombic', 0.0)
            total_lj_coul = lj_energy + coulombic_energy
            
            
            if total_lj_coul == 0:
                return 0.0
            
            # Use empirical data to convert LJ+Coul energy to systemic temperature
            temperature_K = self.empirical_data_lj_coul.calculate_systemic_temperature(total_lj_coul)
            
            return temperature_K
            
        except Exception as e:
            print(f"Warning: Could not calculate LJ+Coul fictive temperature: {e}")
            return 0.0
    
    def _get_cavity_bath_temperature(self) -> float:
        """Get cavity thermostat temperature."""
        try:
            if self.cavity_thermostat is None:
                return 0.0
            
            from ..utils import PhysicalConstants
            
            # Get kT from thermostat - handle different thermostat types
            kT = None
            if hasattr(self.cavity_thermostat, 'kT'):
                kT_value = self.cavity_thermostat.kT
                # Handle both Constant variants and plain floats
                if hasattr(kT_value, 'value'):
                    kT = kT_value.value
                elif hasattr(kT_value, '__call__'):
                    kT = kT_value(0)  # Evaluate at timestep 0
                else:
                    kT = float(kT_value)
            elif hasattr(self.cavity_thermostat, 'T'):
                # Some thermostats use T instead of kT
                T_value = self.cavity_thermostat.T
                if hasattr(T_value, 'value'):
                    temperature_K = T_value.value
                elif hasattr(T_value, '__call__'):
                    temperature_K = T_value(0)
                else:
                    temperature_K = float(T_value)
                return temperature_K
            elif hasattr(self.cavity_thermostat, 'thermostat'):
                # For ConstantVolume methods, check the thermostat attribute
                nested_thermo = self.cavity_thermostat.thermostat
                if hasattr(nested_thermo, 'kT'):
                    kT_value = nested_thermo.kT
                    if hasattr(kT_value, 'value'):
                        kT = kT_value.value
                    elif hasattr(kT_value, '__call__'):
                        kT = kT_value(0)
                    else:
                        kT = float(kT_value)
            
            if kT is None:
                return 0.0
                
            temperature_K = kT / PhysicalConstants.KB_HARTREE_PER_K
            return temperature_K
            
        except Exception as e:
            print(f"Warning: Could not get cavity bath temperature: {e}")
            return 0.0
    
    def _get_molecular_bath_temperature(self) -> float:
        """Get molecular thermostat temperature."""
        try:
            if self.molecular_thermostat is None:
                return 0.0
            
            from ..utils import PhysicalConstants
            
            # Get kT from thermostat - handle different thermostat types
            kT = None
            if hasattr(self.molecular_thermostat, 'kT'):
                kT_value = self.molecular_thermostat.kT
                # Handle both Constant variants and plain floats
                if hasattr(kT_value, 'value'):
                    kT = kT_value.value
                elif hasattr(kT_value, '__call__'):
                    kT = kT_value(0)  # Evaluate at timestep 0
                else:
                    kT = float(kT_value)
            elif hasattr(self.molecular_thermostat, 'T'):
                # Some thermostats use T instead of kT
                T_value = self.molecular_thermostat.T
                if hasattr(T_value, 'value'):
                    temperature_K = T_value.value
                elif hasattr(T_value, '__call__'):
                    temperature_K = T_value(0)
                else:
                    temperature_K = float(T_value)
                return temperature_K
            elif hasattr(self.molecular_thermostat, 'thermostat'):
                # For ConstantVolume methods, check the thermostat attribute
                nested_thermo = self.molecular_thermostat.thermostat
                if hasattr(nested_thermo, 'kT'):
                    kT_value = nested_thermo.kT
                    if hasattr(kT_value, 'value'):
                        kT = kT_value.value
                    elif hasattr(kT_value, '__call__'):
                        kT = kT_value(0)
                    else:
                        kT = float(kT_value)
            
            if kT is None:
                return 0.0
                
            temperature_K = kT / PhysicalConstants.KB_HARTREE_PER_K
            return temperature_K
            
        except Exception as e:
            print(f"Warning: Could not get molecular bath temperature: {e}")
            return 0.0


