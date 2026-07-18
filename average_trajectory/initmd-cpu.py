"""
Main script for running molecular dynamics simulations with optional cavity coupling using HOOMD-blue.
Supports both regular NVT simulations and cavity-coupled simulations with adaptive timesteps.
"""

from cavitymd import *
import argparse
import os
import sys
import gsd.hoomd  # Add this import at the top of the file with other imports
import numpy as np
from numba import njit

@njit
def fibonacci_sphere(samples=100):
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

class DensityCorrelationTracker(hoomd.custom.Action):
    """
    Tracks the density-density correlation function F(k,t) during simulation.
    Computes F(k,t) = <ρₖ(t)·ρₖ*(t₀)> for multiple reference frames.
    Each reference frame gets its own output file, with lag time as the x-axis.
    """
    def __init__(self, simulation, time_tracker, kmag=1.0, num_wavevectors=50, output_period=1000, output_prefix='fskt', reference_interval_ps=1.0, max_references=10):
        super().__init__()
        self.simulation = simulation
        self.time_tracker = time_tracker  # Reference to ElapsedTimeTracker for accurate timing
        self.kmag = kmag  # Wavevector magnitude
        self.num_wavevectors = num_wavevectors  # Number of wavevectors to sample on sphere
        self.output_period = output_period  # How often to output summary results
        self.output_prefix = output_prefix  # Prefix for output files
        self.reference_interval_ps = reference_interval_ps
        self.max_references = max_references

        # Generate wavevectors
        self.wavevectors = fibonacci_sphere(samples=num_wavevectors) * kmag

        # List of reference frames: each is a dict with keys: 'timestep', 'time', 'rhok_real', 'file_path'
        self.references = []
        self.last_reference_time = None

        # Counter for tracking when to output results
        self.last_output_step = 0

        # For caching computed values in current timestep (fixes logging delay issue)
        self.current_computed_value = None
        self.current_timestep = -1

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
            rhok_real = self.compute_rhok(initial_positions)
            t0 = 0.0
            file_path = f'{self.output_prefix}_fskt_k{self.kmag:.2f}_ref0.txt'
            with open(file_path, 'w') as f:
                f.write(f'# F(k,t) data for k={self.kmag:.4f}, t0={t0:.6f} ps\n')
                f.write('# lag_time(ps) F(k,t)\n')
            self.references.append({
                'timestep': 0,
                'time': t0,
                'rhok_real': rhok_real.copy(),
                'file_path': file_path
            })
            self.last_reference_time = t0
        print(f"DensityCorrelationTracker initialized with k={kmag:.2f}, {num_wavevectors} wavevectors, reference_interval_ps={reference_interval_ps}, max_references={max_references}")

    def compute_rhok(self, positions):
        n_particles = positions.shape[0]
        n_wavevectors = len(self.wavevectors)
        rhok_real = np.zeros(n_wavevectors)
        for i in range(n_wavevectors):
            kvec = self.wavevectors[i]
            fac = np.sum(kvec * positions, axis=1)
            cos_fac = np.cos(fac)
            rhok_real[i] = np.sum(cos_fac)
        return rhok_real

    def compute_fskt(self, rhok0_real, rhok_real):
        n_wavevectors = len(self.wavevectors)
        correlations = np.zeros(n_wavevectors)
        for k in range(n_wavevectors):
            real_part = rhok0_real[k] * rhok_real[k]  # Only use real part
            correlations[k] = real_part
        return np.mean(correlations)

    def act(self, timestep):
        # Get current positions and timestamp
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
            current_time = self.time_tracker.elapsed_time

        # Compute rhok for current frame
        current_rhok_real = self.compute_rhok(positions)

        # Add a new reference if enough time has passed and we haven't hit max_references
        if (len(self.references) < self.max_references and
            (self.last_reference_time is None or current_time - self.last_reference_time >= self.reference_interval_ps)):
            file_path = f'{self.output_prefix}_fskt_k{self.kmag:.2f}_ref{len(self.references)}.txt'
            with open(file_path, 'w') as f:
                f.write(f'# F(k,t) data for k={self.kmag:.4f}, t0={current_time:.6f} ps\n')
                f.write('# lag_time(ps) F(k,t)\n')
            self.references.append({
                'timestep': timestep,
                'time': current_time,
                'rhok_real': current_rhok_real.copy(),
                'file_path': file_path
            })
            self.last_reference_time = current_time

        # For each reference, compute and save F(k, t-t0)
        for ref in self.references:
            lag_time = current_time - ref['time']
            if lag_time < 0:  # Should not happen, but skip if so
                continue
            fskt_value = self.compute_fskt(ref['rhok_real'], current_rhok_real)
            with open(ref['file_path'], 'a') as f:
                f.write(f'{lag_time:.6f} {fskt_value:.6f}\n')

        # Optionally, print debug info
        if timestep - self.last_output_step >= self.output_period:
            print(f"F(k,t) written for {len(self.references)} references at t={current_time:.4f} ps")
            self.last_output_step = timestep

    @hoomd.logging.log
    def current_fskt(self):
        # For logging, just return the value for the first reference (if any)
        if not self.references:
            return 0.0
        # Get current positions
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
        current_rhok_real = self.compute_rhok(positions)
        ref = self.references[0]
        return self.compute_fskt(ref['rhok_real'], current_rhok_real)

class ElapsedTimeTracker(hoomd.custom.Action):
    """Tracks the total elapsed time in a simulation with variable timesteps."""
    def __init__(self, simulation, runtime):
        super().__init__()
        self.simulation = simulation
        self.total_time = 0.0
        self.runtime = runtime
        self.last_timestep = simulation.timestep
        
    def act(self, timestep):
        """Update the total elapsed time."""
        current_timestep = self.simulation.timestep
        if current_timestep > self.last_timestep:
            self.total_time += (current_timestep - self.last_timestep) * self.simulation.operations.integrator.dt
        self.last_timestep = current_timestep
        
        # Check if we've reached the runtime and exit if so
        if self.total_time*0.02418884/1000.0 >= self.runtime:
            print(f"Runtime {self.runtime} ps reached. Exiting simulation.")
            sys.exit(0)

    @hoomd.logging.log #@property
    def elapsed_time(self):
        """Expose the total elapsed time as a property in (ps)."""
        return self.total_time*0.02418884/1000.0

class AdaptiveTimestepUpdater(hoomd.custom.Action):
    def __init__(self, state, integrator, error_tolerance, time_constant_ps=50.0, initial_fraction=0.01, adaptiveerror=True):
        """
        Parameters:
        - state: The simulation state.
        - integrator: The integrator in the simulation.
        - error_tolerance: Target error tolerance to reach.
        - time_constant_ps: Time constant in picoseconds for exponential approach.
        - initial_fraction: Initial error tolerance as a fraction of the target value.
        """
        super().__init__()
        print("Performing error tolerance ramping with time constant", time_constant_ps, "ps")
        self.integrator = integrator
        self.target_error_tolerance = error_tolerance
        self.state = state
        self.initial_fraction = initial_fraction
        self.time_constant_ps = time_constant_ps
        
        # Calculate initial values - use a more conservative starting point for safety
        self.initial_error_tolerance = error_tolerance * initial_fraction
        
        # For continuation runs, we want to be more conservative with the initial timestep
        # to prevent instability when restarting
        #particle_data = self.state.get_snapshot().particles
        # positions = np.array(particle_data.position)
        # velocities = np.array(particle_data.velocity)
        self.current_error_tolerance = self.initial_error_tolerance
        # # Check if this seems to be a continuation (non-zero velocities)
        # avg_vel_mag = np.mean(np.linalg.norm(velocities, axis=1))
        
        # if avg_vel_mag > 0.01:  # This is likely a continuation run
        #     print(f"Detected continuation run (avg velocity: {avg_vel_mag:.4f})")
        #     print("Using more conservative initial timestep for stability")
        #     # Start with a safety factor for continuation runs
        #     self.current_error_tolerance = self.initial_error_tolerance * 0.1
        # else:
        #     # For fresh runs, use the normal initial fraction
        #     self.current_error_tolerance = self.initial_error_tolerance
        
        # Keep track of accumulated simulation time
        self.accumulated_time_ps = 0.0
        self.last_timestep = 0
        self.adaptiveerror = adaptiveerror
        
    def act(self, timestep):
        """
        Custom action to update the timestep size.

        Parameters:
        - timestep: Current simulation timestep.
        """
        # Update accumulated simulation time
        if timestep > self.last_timestep:
            # Convert dt to picoseconds
            dt_ps = self.integrator.dt * 0.02418884  # Convert from a.u. to ps
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
        forces = []
        for force in self.integrator.forces:
            particle_forces = force.forces
            forces.append(particle_forces)
        
        particle_data = self.state.get_snapshot().particles
        masses = np.array(particle_data.mass)
        forces = np.array(forces)
        
        # Calculate sum |f_i| / m_i
        forces = np.sum(forces, axis=0)    
        force_norm = np.array([np.linalg.norm(f) for f in forces])
        force_mass_sum = np.sum(force_norm / masses)
        
        # Update timestep using current error tolerance
        if force_mass_sum > 0:
            new_dt = np.sqrt(self.current_error_tolerance / force_mass_sum)
            
            # Add safety limit for the timestep to prevent large jumps
            #current_dt = self.integrator.dt
            #max_change_factor = 1.2  # Limit timestep increases to 20% at a time
            #if new_dt > current_dt * max_change_factor:
            #    new_dt = current_dt * max_change_factor
                
            self.integrator.dt = new_dt
            # Update gamma for the common Langevin method
            self.integrator.methods[0].default_gamma = 1/(500*new_dt)
    
    @hoomd.logging.log
    def error_tolerance(self):
        """Log the current effective error tolerance."""
        return self.current_error_tolerance
        
    @hoomd.logging.log
    def elapsed_time_ps(self):
        """Log the elapsed simulation time in picoseconds."""
        return self.accumulated_time_ps

class BussiHeatExchangeTracker(hoomd.custom.Action):
    """
    Tracks heat exchange from the Bussi thermostat by monitoring kinetic energy changes
    of the molecular system (excluding cavity particles).
    """
    def __init__(self, simulation, time_tracker, output_prefix='bussi_heat', output_period=1000, max_time_ps=None):
        super().__init__()
        self.simulation = simulation
        self.time_tracker = time_tracker
        self.output_prefix = output_prefix
        self.output_period = output_period
        self.max_time_ps = max_time_ps
        self.last_output_step = 0
        
        # Previous step data
        self.previous_ke = None
        self.previous_time = None
        self.previous_velocities = None
        self.previous_forces = None
        
        # Accumulated heat exchange
        self.cumulative_heat_exchange = 0.0
        self.instantaneous_heat_rate = 0.0
        
        # File for detailed output
        self.heat_file = open(f'{output_prefix}.txt', 'w')
        self.heat_file.write('# Time(ps) InstantaneousKE KEChange ForceContribution ThermostatContribution CumulativeHeatExchange HeatRate(per_ps)\n')
        
    def compute_molecular_kinetic_energy(self, snap):
        """Compute kinetic energy of molecular particles (O and N types only)."""
        # Filter for molecular particles (exclude cavity particle type 'L' which is typeid=2)
        molecular_mask = snap.particles.typeid != 2
        velocities = snap.particles.velocity[molecular_mask]
        masses = snap.particles.mass[molecular_mask]
        
        # Calculate kinetic energy: KE = 0.5 * m * v^2
        ke = 0.5 * np.sum(masses * np.sum(velocities**2, axis=1))
        return ke, velocities, masses
    
    def estimate_force_contribution(self, snap, dt):
        """Estimate the kinetic energy change due to forces."""
        if self.previous_velocities is None or self.previous_forces is None:
            return 0.0
        
        # Get current forces on molecular particles
        molecular_mask = snap.particles.typeid != 2
        current_forces = np.zeros_like(snap.particles.position)
        
        # We can't directly access forces in HOOMD easily, so we'll use an approximation
        # based on acceleration from velocity changes
        current_velocities = snap.particles.velocity[molecular_mask]
        masses = snap.particles.mass[molecular_mask]
        
        # Estimate acceleration from velocity change
        dv = current_velocities - self.previous_velocities
        acceleration = dv / dt
        estimated_forces = masses[:, np.newaxis] * acceleration
        
        # Work done by forces: W = F·v·dt (average of initial and final velocities)
        avg_velocities = 0.5 * (self.previous_velocities + current_velocities)
        work_done = np.sum(estimated_forces * avg_velocities) * dt
        
        return work_done
    
    def act(self, timestep):
        current_time = self.time_tracker.elapsed_time
        
        # Check if we should stop outputting based on max_time_ps
        if self.max_time_ps is not None and current_time > self.max_time_ps:
            return
        
        with self.simulation.state.cpu_local_snapshot as snap:
            # Compute current kinetic energy of molecular system
            current_ke, current_velocities, masses = self.compute_molecular_kinetic_energy(snap)
            
            if self.previous_ke is not None and self.previous_time is not None:
                dt = current_time - self.previous_time
                if dt > 0:
                    # Calculate kinetic energy change
                    ke_change = current_ke - self.previous_ke
                    
                    # Estimate the contribution from forces
                    force_contribution = self.estimate_force_contribution(snap, dt)
                    
                    # The remainder is attributed to thermostat action
                    thermostat_contribution = ke_change - force_contribution
                    
                    # Accumulate heat exchange (negative means heat removed from system)
                    self.cumulative_heat_exchange += thermostat_contribution
                    self.instantaneous_heat_rate = thermostat_contribution / dt
                    
                    # Write to file
                    self.heat_file.write(f'{current_time:.6f} {current_ke:.6f} {ke_change:.6f} {force_contribution:.6f} {thermostat_contribution:.6f} {self.cumulative_heat_exchange:.6f} {self.instantaneous_heat_rate:.6f}\n')
                    self.heat_file.flush()
            
            # Store current data for next step
            self.previous_ke = current_ke
            self.previous_time = current_time
            self.previous_velocities = current_velocities.copy()
            
        # Debug output
        if timestep - self.last_output_step >= self.output_period:
            print(f"Heat exchange at t={current_time:.4f} ps: Cumulative={self.cumulative_heat_exchange:.6f}, Rate={self.instantaneous_heat_rate:.6f}")
            self.last_output_step = timestep
    
    @hoomd.logging.log
    def molecular_kinetic_energy(self):
        """Current kinetic energy of molecular system."""
        if self.previous_ke is not None:
            return self.previous_ke
        return 0.0
    
    @hoomd.logging.log
    def cumulative_heat_removed(self):
        """Cumulative heat removed by thermostat (negative of heat exchange)."""
        return -self.cumulative_heat_exchange
    
    @hoomd.logging.log
    def heat_removal_rate(self):
        """Instantaneous rate of heat removal by thermostat."""
        return -self.instantaneous_heat_rate
    
    @hoomd.logging.log
    def total_heat_exchange(self):
        """Total cumulative heat exchange."""
        return self.cumulative_heat_exchange
    
    def __del__(self):
        if hasattr(self, 'heat_file'):
            self.heat_file.close()

def wrap_position(x, L):
    # Compute the image flags (how many box lengths away from the primary box)
    image_flags = np.floor((x + L/2) / L)

    # Compute the wrapped position inside the primary box
    wrapped_position = x - image_flags * L

    return wrapped_position, image_flags.astype(int)

def main(job_dir: str, replica: int, freq: float, couplstr: float, incavity: bool, runtime_ps: float = 500.0, input_gsd: str = 'molecular-0.gsd', frame: int = -1, name: str = 'prod', error_tolerance: float = 0.01, compute_fskt: bool = False, kmag: float = 1.0, truncate_gsd: bool = True, compute_density_mode: bool = False, density_kx: float = 1.0):
    global stop_simulation
    """
    Main function to run molecular dynamics simulation.
    
    Sets up and runs either a standard NVT simulation or a cavity-coupled simulation
    with specified parameters.
    
    Args:
        job_dir: Directory for simulation input/output
        replica: Replica index for this simulation
        freq: Cavity frequency in cm^-1 (for cavity coupling)
        couplstr: Coupling strength (for cavity coupling)
        incavity: Whether to include cavity coupling
        runtime_ps: Total simulation time in picoseconds
        input_gsd: Input GSD file to initialize the simulation from
        frame: Frame number to use from the input GSD file (-1 for last frame)
        name: Prefix for output files (default: 'prod')
        error_tolerance: Error tolerance for adaptive timestepping (default: 0.01)
        compute_fskt: Whether to compute F(k,t) on the fly (default: False)
        kmag: Wavevector magnitude for F(k,t) calculation (default: 1.0)
        truncate_gsd: Whether to truncate existing GSD files (default: True)
        compute_density_mode: Whether to compute density mode at wavevector in x direction
        density_kx: Magnitude of wavevector in x direction for density mode (default: 1.0)
    """
    #try:
    cwd = os.getcwd()
    os.chdir(job_dir)
    

    # Define physical constants and unit conversions
    T = 100 #K
    
    #We're going to follow atomic units!
    energy = 4.35974e-18 #Joules
    length = 5.29177210544e-11 #m 
    mass = 9.1093837139e-31 #kg
    time = length*(mass/energy)**0.5 #s, same as a.u. units convention
    kB = 3.167e-6 #Hartree/K
    # Simulation parameters
    omegac = freq

    # Time stepping parameters
    dt = 0.0005 #ps 
    runtime_real = runtime_ps 
    runtime = int(runtime_real/dt)
    
    period = 0.1 #ps
    table_period = 0.1 #ps
    
    period = int(period/dt)
    table_period = int(table_period/dt)
    
    if period < 1:
        period = 1
    if table_period < 1:
        table_period = 1
    
    print(f"Period: {period}, Table period: {table_period}")
    timeps = time*1e+12
    dt = dt/timeps #in a.u. units
    
    rcut = 15 #Bohr units
    
    ## Setup the MD Simulation
    # First load the GSD file using gsd.hoomd
    with gsd.hoomd.open(input_gsd, 'r') as f:
        if frame < 0:
            frame = len(f) + frame  # Convert negative index to positive
            if frame < 0:  # Handle case where abs(frame) > len(f)
                frame = 0
        snapshot = f[frame]
        
        if incavity:
            # Calculate dipole moment and photon position
            positions = unwrap_positions(snapshot.particles.position, snapshot.particles.image, snapshot.configuration.box[:3])
            dipmom = np.einsum('i,ij->j', snapshot.particles.charge, positions)

            hartree_to_cm_minus1 = 219474.63
            omegac = freq/hartree_to_cm_minus1
            newpos = -0.0*dipmom*couplstr/omegac**2
            sigma = np.sqrt(kB*T / omegac**2)
            newpos = np.random.normal(loc=newpos, scale=sigma, size=3)
            newpos[-1] = 0.0
            newpos, image_flags = wrap_position(newpos, np.array(snapshot.configuration.box[:3]))
            newpos[-1] = 0.0
            image_flags[-1] = 0
            #image_flags = np.zeros(3, dtype=int)
            
            # Add photon particle
            snapshot.particles.types.append('L')
            snapshot.particles.N += 1
            snapshot.particles.typeid = np.append(snapshot.particles.typeid, [2])
            snapshot.particles.position = np.append(snapshot.particles.position, [newpos], axis=0)
            snapshot.particles.charge = np.append(snapshot.particles.charge, [0.0])
            snapshot.particles.mass = np.append(snapshot.particles.mass, [1.0])
            snapshot.particles.diameter = np.append(snapshot.particles.diameter, [1.0])
            snapshot.particles.image = np.vstack([snapshot.particles.image, image_flags])

            # Set body to -1 for the photon
            if hasattr(snapshot.particles, "body"):
                snapshot.particles.body = np.append(snapshot.particles.body, [-1], axis=0)

            # Set orientation to [1,0,0,0] for the photon
            if hasattr(snapshot.particles, "orientation"):
                snapshot.particles.orientation = np.append(
                    snapshot.particles.orientation,
                    [[1.0, 0.0, 0.0, 0.0]],
                    axis=0
                )

            # Set moment_inertia to [0,0,0] for the photon
            if hasattr(snapshot.particles, "moment_inertia"):
                snapshot.particles.moment_inertia = np.vstack([
                    snapshot.particles.moment_inertia,
                    np.zeros((1, 3))
                ])

            # Set velocity to [0,0,0] for the photon
            if hasattr(snapshot.particles, "velocity"):
                snapshot.particles.velocity = np.vstack([
                    snapshot.particles.velocity,
                    np.zeros((1, 3))
                ])

            # Set angmom to [0,0,0,0] for the photon
            if hasattr(snapshot.particles, "angmom"):
                snapshot.particles.angmom = np.vstack([
                    snapshot.particles.angmom,
                    np.zeros((1, 4))
                ])

    # Create a communicator instance
    comm = hoomd.communicator.Communicator()

    
    # Convert to HOOMD snapshot
    hoomd_snapshot = hoomd.Snapshot.from_gsd_snapshot(snapshot, communicator=comm)

    # Add a cavity particle to the snapshot if needed
    # Create simulation from the original snapshot
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=np.random.randint(int(10**4)))
    sim.create_state_from_snapshot(hoomd_snapshot)
    #sim.state.thermalize_particle_momenta(kT=kB*T,filter=hoomd.filter.All())
    
    cell = hoomd.md.nlist.Cell(buffer=1.0,exclusions=('bond',))
    
    #Setup baro-thermostat for the particle system
    list_forces = []

    #Setup cavity force if requested
    if incavity:
        # Create the cavity force
        cavityforce = CavityForce(kvector=np.array([0,0,1]), couplstr=couplstr, omegac=freq)#omegac)
        # Add it to the list of forces
        list_forces.append(cavityforce)

    # Set up a single Langevin integrator for all particle types with reservoir energy tracking
    common_langevin = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), 
                                               kT=kB*T, default_gamma=1/(500*dt), tally_reservoir_energy=True)
    
    sim.state.thermalize_particle_momenta(kT=kB*T,filter=hoomd.filter.All())

    #Apply the harmonic potential on the bonds.
    harmonic = hoomd.md.bond.Harmonic()
    harmonic.params['O-O'] = dict(k=2*0.36602, r0=2.281655158)
    harmonic.params['N-N'] = dict(k=2*0.71625, r0=2.0743522177)
    list_forces.append(harmonic)

    #Apply lennard jones interaction
    lj = hoomd.md.pair.LJ(nlist=cell,mode='shift')
    lj.params[('O', 'O')] = dict(epsilon=0.00016685201, sigma=6.230426584)
    lj.r_cut[('O', 'O')] = rcut
    lj.params[('N', 'N')] = dict(epsilon=0.000083426, sigma=5.48277488)
    lj.r_cut[('N', 'N')] = rcut
    lj.params[('N', 'O')] = dict(epsilon=0.00025027802, sigma=4.9832074319)
    lj.r_cut[('N', 'O')] = rcut

    #Disable pair interaction with 'L' particle (photon)
    if incavity:
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
    list_forces.append(lj)

    #Apply long range Coulomb interactions using PPPM method
    r_cut = rcut
    numpoints = 32
    order = 6
    short, long  = hoomd.md.long_range.pppm.make_pppm_coulomb_forces(nlist=cell, resolution=[numpoints,numpoints,numpoints], order=order, r_cut=r_cut, alpha=0.0)
    list_forces.append(short)
    list_forces.append(long)
    
    #Setup integrator and run simulation
    integrator = hoomd.md.Integrator(dt=dt,forces=list_forces)
    sim.operations.integrator = integrator
    sim.operations.integrator.methods = [common_langevin]
    
    #Setup status update logging
    status = Status(sim, time*1e9)
    
    # Initialize time tracker for variable timestep simulations
    time_tracker = ElapsedTimeTracker(sim, runtime_ps)
    
    # Create a logger for console output (can include strings)
    console_logger = hoomd.logging.Logger(categories=['scalar', 'string'])
    console_logger.add(sim, quantities=['timestep', 'tps'])
    console_logger[('Status', 'nsd')] = (status, 'nsd', 'string')
    console_logger[('Status', 'etr')] = (status, 'etr', 'string')
    console_logger[('Status', 'Dt')] = (status, 'Dt', 'scalar')
    console_logger[('Status', 'elapsed')] = (status, 'elapsed', 'string')
    console_logger.add(time_tracker)
    
    # Set up F(k,t) calculation if requested
    if compute_fskt:
        # Set up DensityCorrelationTracker
        fskt_tracker = DensityCorrelationTracker(
            simulation=sim,
            time_tracker=time_tracker,
            kmag=kmag,
            num_wavevectors=50,
            output_period=1,
            output_prefix=f'{name}-{replica}',
            reference_interval_ps=args.fskt_reference_interval_ps,
            max_references=args.fskt_max_references
        )
        
        # Add the tracker to loggers
        console_logger[('FsktTracker', 'current_fskt')] = (fskt_tracker, 'current_fskt', 'scalar')
        
        # Add the tracker as an updater
        fskt_updater = hoomd.update.CustomUpdater(
            action=fskt_tracker,
            trigger=hoomd.trigger.Periodic(period)  # Same period as GSD output
        )
        sim.operations.updaters.append(fskt_updater)
    
    # Set up density mode calculation if requested
    if compute_density_mode:
        # Set up DensityModeTracker for single wavevector in x direction
        density_mode_tracker = DensityModeTracker(
            simulation=sim,
            time_tracker=time_tracker,
            kx_magnitude=density_kx,
            output_prefix=f'{name}-{replica}',
            output_period=1  # Output at every timestep
        )
        
        # Add the tracker to console logger
        console_logger[('DensityMode', 'density_real')] = (density_mode_tracker, 'density_real', 'scalar')
        console_logger[('DensityMode', 'density_imag')] = (density_mode_tracker, 'density_imag', 'scalar')
        console_logger[('DensityMode', 'density_magnitude')] = (density_mode_tracker, 'density_magnitude', 'scalar')
        
        # Add the tracker as an updater
        density_mode_updater = hoomd.update.CustomUpdater(
            action=density_mode_tracker,
            trigger=hoomd.trigger.Periodic(1)  # Output at every timestep
        )
        sim.operations.updaters.append(density_mode_updater)
        
        # Add to GSD logger if available
        if 'gsd_logger' in locals():
            gsd_logger[('DensityMode', 'density_real')] = (density_mode_tracker, 'density_real', 'scalar')
            gsd_logger[('DensityMode', 'density_imag')] = (density_mode_tracker, 'density_imag', 'scalar')
            gsd_logger[('DensityMode', 'density_magnitude')] = (density_mode_tracker, 'density_magnitude', 'scalar')
    
    # Set up kinetic energy tracking for molecular particles
    kinetic_energy_tracker = KineticEnergyTracker(
        simulation=sim,
        time_tracker=time_tracker,
        output_prefix=f'{name}-{replica}',
        output_period=100000  # Disable separate file output (use large period)
    )
    
    # Add the kinetic energy tracker to console logger
    console_logger[('KineticEnergy', 'kinetic_energy')] = (kinetic_energy_tracker, 'kinetic_energy', 'scalar')
    console_logger[('KineticEnergy', 'temperature')] = (kinetic_energy_tracker, 'temperature', 'scalar')
    
    # Don't add kinetic energy tracker as a separate updater since the combined energy tracker will handle it
    
    # Set up cavity mode tracking if in cavity simulation
    if incavity:
        cavity_mode_tracker = CavityModeTracker(
            simulation=sim,
            cavityforce=cavityforce,
            time_tracker=time_tracker,
            output_prefix=f'{name}-{replica}',
            output_period=100000  # Disable separate file output (use large period)
        )
        
        # Add the cavity mode tracker to console logger
        console_logger[('CavityMode', 'cavity_kinetic_energy')] = (cavity_mode_tracker, 'cavity_kinetic_energy', 'scalar')
        console_logger[('CavityMode', 'cavity_potential_energy_harmonic')] = (cavity_mode_tracker, 'cavity_potential_energy_harmonic', 'scalar')
        console_logger[('CavityMode', 'cavity_total_energy')] = (cavity_mode_tracker, 'cavity_total_energy', 'scalar')
        
        # Don't add cavity mode tracker as a separate updater since the combined energy tracker will handle it
    else:
        cavity_mode_tracker = None
    
    # Set up energy contribution tracker (now combining both potential and kinetic energies)
    energy_tracker = EnergyContributionTracker(
        simulation=sim,
        harmonic=harmonic,
        lj=lj,
        short=short,
        long=long,
        cavityforce=cavityforce if incavity else None,
        langevin_method=common_langevin,  # Pass the common Langevin method for reservoir energy
        molecular_langevin_method=common_langevin,  # Pass the common Langevin method for reservoir energy
        kinetic_tracker=kinetic_energy_tracker,  # Pass molecular kinetic energy tracker
        cavity_mode_tracker=cavity_mode_tracker,  # Pass cavity mode tracker (None if not in cavity)
        time_tracker=time_tracker,
        output_prefix=f'{name}-{replica}',
        output_period=1,  # Output at same frequency as GSD
        max_time_ps=args.max_energy_output_time_ps
    )
    
    # Set up GSD output if requested
    print("Setting up GSD output")
    # Create a logger for GSD file (numeric data only)
    gsd_logger = hoomd.logging.Logger(categories=['scalar'])
    gsd_logger.add(sim, quantities=['timestep', 'tps'])
    # Only include numeric data in the GSD logger
    gsd_logger.add(time_tracker)
    # Optionally add other scalar quantities if needed
    gsd_logger[('integrator', 'dt')] = (integrator, 'dt', 'scalar')
    
    # Add total potential energy from energy tracker instead of individual contributions
    gsd_logger[('EnergyTracker', 'total_potential_energy')] = (energy_tracker, 'total_potential_energy', 'scalar')

    # Add kinetic energy tracking to GSD logger
    gsd_logger[('KineticEnergy', 'kinetic_energy')] = (kinetic_energy_tracker, 'kinetic_energy', 'scalar')
    gsd_logger[('KineticEnergy', 'temperature')] = (kinetic_energy_tracker, 'temperature', 'scalar')

    # Add cavity mode energy tracking to GSD logger if in cavity simulation
    if incavity:
        gsd_logger[('CavityMode', 'cavity_kinetic_energy')] = (cavity_mode_tracker, 'cavity_kinetic_energy', 'scalar')
        gsd_logger[('CavityMode', 'cavity_potential_energy_harmonic')] = (cavity_mode_tracker, 'cavity_potential_energy_harmonic', 'scalar')
        gsd_logger[('CavityMode', 'cavity_total_energy')] = (cavity_mode_tracker, 'cavity_total_energy', 'scalar')

    # Add F(k,t) to GSD logger if enabled
    if compute_fskt:
        gsd_logger[('FsktTracker', 'current_fskt')] = (fskt_tracker, 'current_fskt', 'scalar')
    
    # Set up GSD writer with the GSD logger to record numeric quantities
    gsd_writer = hoomd.write.GSD(
        filename=f'{name}-{replica}.gsd',
        trigger=hoomd.trigger.Periodic(100),
        dynamic=['property', 'momentum', 'particles/diameter', 'topology'],
        mode='wb',
        truncate=truncate_gsd,
        filter=hoomd.filter.All()
    )
    gsd_writer.logger = gsd_logger  # Attach the GSD logger to the GSD writer
    
    # Write the initial frame (frame 0)
    gsd_writer.write(sim.state, filename=f'{name}-{replica}.gsd',
        mode='wb',
        filter=hoomd.filter.All(), logger=gsd_logger)
    
    # Add the GSD writer to the simulation's writers
    sim.operations.writers.append(gsd_writer)
    print(f"GSD writer added for file: {name}-{replica}.gsd")

    # Add the energy tracker as an updater
    energy_tracker_updater = hoomd.update.CustomUpdater(
        action=energy_tracker,
        trigger=hoomd.trigger.Periodic(1)  # Output at every timestep
    )
    sim.operations.updaters.append(energy_tracker_updater)
    
    # Add total potential energy to console logger for monitoring
    console_logger[('EnergyTracker', 'total_potential_energy')] = (energy_tracker, 'total_potential_energy', 'scalar')

    # Set up the console output table with the console logger
    table = hoomd.write.Table(
        trigger=hoomd.trigger.Periodic(period=table_period),
        logger=console_logger
    )
    sim.operations.writers.append(table)

    # Add the time tracker as an updater so it records elapsed time
    time_updater = hoomd.update.CustomUpdater(
        action=time_tracker, 
        trigger=hoomd.trigger.Periodic(1)
    )
    sim.operations.updaters.append(time_updater)
    
    # Set up the adaptive timestep updater action if error_tolerance is positive
    if error_tolerance > 0:
        adaptive_action = AdaptiveTimestepUpdater(
            state=sim.state, 
            integrator=sim.operations.integrator, 
            error_tolerance=error_tolerance,
            time_constant_ps=15.0,
            initial_fraction=1e-3,
            adaptiveerror=True 
        )
        
        # Add error tolerance to GSD logger if available
        if 'gsd_logger' in locals():
            gsd_logger.add(adaptive_action, quantities=['error_tolerance', 'elapsed_time_ps'])
        
        # Wrap the action in a CustomUpdater
        adaptive_updater = hoomd.update.CustomUpdater(
            action=adaptive_action, 
            trigger=hoomd.trigger.Periodic(period)
        )
        sim.operations.updaters.append(adaptive_updater)

    # Add dipole autocorrelation updater
    dipole_autocorr = DipoleAutocorrelation(sim, time_tracker=time_tracker, output_prefix=f'{name}-{replica}_dipole')
    dipole_autocorr_updater = hoomd.update.CustomUpdater(
        action=dipole_autocorr,
        trigger=hoomd.trigger.Periodic(1)
    )
    sim.operations.updaters.append(dipole_autocorr_updater)
    
    # Add dipole autocorrelation to console logger
    console_logger[('DipoleAutocorr', 'current_autocorr')] = (dipole_autocorr, 'current_autocorr', 'scalar')

    # Run the simulation
    
    sim.run(int(100*runtime),write_at_start=True)
    os.chdir(cwd)
    return 0  # Success

def parse_args():
    """
    Parse command line arguments for simulation parameters.
    
    Returns:
        Namespace containing parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Parameters for Cavity MD Simulation"
    )
    parser.add_argument(
        "--job_dir",
        type=str,
        required=True,
        help="Directory containing inputs as well as where outputs should be stored",
    )
    parser.add_argument(
        "--replica", type=int, required=True,help="Run simulation for this replica index"
    )
    parser.add_argument(
        "--freq", type=float, required=True,help="Photon frequency (in cm-1)"
    )
    parser.add_argument(
        "--couplstr", type=float, required=True,help="Coupling strength without square root N factor (in a.u.)"
    )
    parser.add_argument(
        "--incavity", action="store_true", help="Do simulations inside a cavity or not"
    )
    parser.add_argument(
        "--runtime", type=float, default=10.0, help="Runtime in picoseconds (default: 10 ps)"
    )
    parser.add_argument(
        "--input_gsd", type=str, default="molecular-0.gsd", help="Input GSD file to initialize the simulation from"
    )
    parser.add_argument(
        "--frame", type=int, default=-1, help="Frame number to use from the input GSD file (-1 for last frame)"
    )
    parser.add_argument(
        "--name", type=str, default="prod", help="Prefix for output files (default: 'prod')"
    )
    parser.add_argument(
        "--error_tolerance", type=float, default=1.0, 
        help="Error tolerance for adaptive timestepping (default: 1.0, set to 0 to disable)"
    )
    parser.add_argument(
        "--compute_fskt", action="store_true", help="Compute F(k,t) on the fly"
    )
    parser.add_argument(
        "--kmag", type=float, default=1.0, help="Wavevector magnitude for F(k,t) calculation (default: 1.0)"
    )
    parser.add_argument(
        "--truncate_gsd", action="store_true", default=False, 
        help="Truncate existing GSD files (default: False)"
    )
    parser.add_argument(
        "--fskt_reference_interval_ps", type=float, default=1.0,
        help="Interval in ps between F(k,t) reference frames (default: 1.0 ps)"
    )
    parser.add_argument(
        "--fskt_max_references", type=int, default=10,
        help="Maximum number of F(k,t) reference frames (default: 10)"
    )
    parser.add_argument(
        "--max_energy_output_time_ps", type=float, default=None,
        help="Maximum time in ps to output energy contributions (default: None, no limit)"
    )
    parser.add_argument(
        "--compute_density_mode", action="store_true", 
        help="Compute density mode at wavevector in x direction"
    )
    parser.add_argument(
        "--density_kx", type=float, default=1.0,
        help="Magnitude of wavevector in x direction for density mode (default: 1.0)"
    )
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    try:
        exit_code = main(args.job_dir, args.replica, args.freq, args.couplstr, args.incavity, 
                         args.runtime, args.input_gsd, args.frame, args.name, args.error_tolerance,
                         args.compute_fskt, args.kmag, args.truncate_gsd, args.compute_density_mode, args.density_kx)
        sys.exit(exit_code)
    except Exception as e:
        import traceback
        print(f"CRITICAL ERROR: {str(e)}")
        traceback.print_exc()
        sys.exit(1)
