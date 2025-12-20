# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Main simulation framework for cavity molecular dynamics."""

from typing import Optional, List, Dict, Union, Any, Tuple
import hoomd
import numpy as np
import logging
import sys
import os
import time
from pathlib import Path
import gsd.hoomd
from hoomd.bussi_reservoir.thermostats import BussiReservoir as Bussi

# CuPy import with fallback for CPU/GPU agnostic code
try:
    import cupy as cp
    HAS_CUPY = True
except ImportError:
    cp = None
    HAS_CUPY = False

from .utils import PhysicalConstants, unwrap_positions
from .forces import CavityForce
from .analysis import (
    Status, ElapsedTimeTracker, TimestepFormatter, FieldAutocorrelationTracker,
    EnergyTracker, PerformanceTracker, AutocorrelationTracker
)
from .updaters import CavityParticleDisplacer
from .variants import StepVariant


class CavityMDSimulation:
    r"""
    Main class for cavity molecular dynamics simulation setup and execution.
    
    This class provides a comprehensive framework for running cavity MD simulations
    with smart cavity particle handling, parameter validation, and advanced analysis
    capabilities. It implements the single-mode cavity-molecule coupling theory
    with support for both constant and time-varying parameters.
    
    **Physics Implementation:**
    
    The simulation implements the equations of motion from the single-mode cavity theory:
    
    **Nuclear Motion:**
    
    .. math::
        
        M_{nj}\ddot{R}_{nj} = F_{nj}^{(0)} - g(t)\tilde{\varepsilon}_{0,\lambda}^{(0)}\tilde{q}_{0,\lambda} \frac{\partial d_{ng,\lambda}}{\partial R_{nj}} + 
                              \frac{g(t)^2(\tilde{\varepsilon}_{0,\lambda}^{(0)})^2}{K} \sum_{l=1}^{N_{\text{sub}}} d_{lg,\lambda} \frac{\partial d_{ng,\lambda}}{\partial R_{nj}}
    
    **Cavity Mode Dynamics:**
    
    .. math::
        
        m_{0,\lambda}\ddot{\tilde{q}}_{0,\lambda} = -K \tilde{q}_{0,\lambda} - g(t)\tilde{\varepsilon}_{0,\lambda}^{(0)} \sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda} - \gamma(t) \dot{\tilde{q}}_{0,\lambda}
    
    where :math:`g(t)` and :math:`\gamma(t)` can be time-varying for dynamic experiments.
    
    **Smart Cavity Particle Handling:**
    
    The simulation automatically:
    
    - Detects existing cavity particles in GSD files
    - Validates cavity particle properties and count
    - Adds new cavity particles only when needed
    - Provides clear error messages for misconfigurations
    - Supports both q=0 and finite-q cavity modes
    
    **Time-Varying Coupling Support:**
    
    For dynamic experiments, coupling can switch at time :math:`t_{\text{switch}}`:
    
    .. math::
        
        g(t) = \begin{cases}
        0 & \text{if } t < t_{\text{switch}} \\
        g_{\text{target}} & \text{if } t \geq t_{\text{switch}}
        \end{cases}
    
    This enables study of non-equilibrium cavity-coupling phenomena.
    
    Parameters
    ----------
    job_dir : str
        Output directory for simulation files
    replica : int
        Replica number for this simulation instance
    freq : float
        Cavity frequency in cm⁻¹ (wavenumbers)
    couplstr : float
        Cavity coupling strength in atomic units
    incavity : bool
        Whether to enable cavity coupling
    runtime_ps : float, optional
        Simulation runtime in picoseconds. Default: 500.0
    input_gsd : str, optional
        Path to input GSD file. Default: 'molecular-0.gsd'
    frame : int, optional
        Frame number to read from GSD file. Default: -1 (last frame)
    name : str, optional
        Simulation name prefix. Default: 'prod'
    error_tolerance : float, optional
        Error tolerance for adaptive timestep. Default: 0.01. Set to 0.0 for fixed timestep.
    temperature : float, optional
        System temperature in Kelvin. Default: 100.0
    molecular_thermostat : str, optional
        Thermostat for molecular particles ('bussi', 'langevin', 'none'). Default: 'bussi'
    cavity_thermostat : str, optional  
        Thermostat for cavity particle ('bussi', 'langevin', 'none'). Default: 'langevin'
    finite_q : bool, optional
        Whether to use finite-q cavity mode. Default: False (q=0 mode)
    molecular_thermostat_tau : float, optional
        Molecular thermostat time constant in ps. Default: 5.0
    cavity_thermostat_tau : float, optional
        Cavity thermostat time constant in ps. Default: 5.0
    switch_time_ps : float, optional
        Time in ps to switch coupling from 0 to target value. Default: None (constant coupling)
    dissipation : float, optional
        Cavity dissipation coefficient in atomic units. Default: 0.0
    decay_time_constant_ps : float, optional
        Exponential decay time constant in picoseconds for time-varying coupling.
        If provided with switch_time_ps, coupling decays as exp(-(t-t_switch)/τ_decay)
        after switching. Default: None (no decay, standard step function)
    enable_energy_tracking : bool, optional
        Whether to enable comprehensive energy tracking. Default: False
    enable_fkt : bool, optional
        Whether to enable F(k,t) density correlation analysis. Default: True
    device : str, optional
        Compute device ('CPU' or 'GPU'). Default: 'CPU'
    seed : int, optional
        Random seed for reproducibility. Default: None (replica-based seed)
    restart_velocities : bool, optional
        Whether to thermalize velocities at start. Default: True
    initial_fraction : float, optional
        Initial fraction for shock dampening (ratio of initial to target error tolerance). Default: 1e-5
    time_constant_ps : float, optional
        Time constant for error tolerance ramping in ps. Default: 50.0
    zero_momentum_enabled : bool, optional
        Whether to periodically zero out system momentum to prevent center-of-mass drift. Default: False
    zero_momentum_period_ps : float, optional
        Period for momentum zeroing in picoseconds. Default: 1.0
        
    **Advanced Parameters:**
        
    dt_fs : float, optional
        Fixed timestep in femtoseconds (only used if error_tolerance=0.0)
    fkt_kmag : float, optional
        F(k,t) k-magnitude for density correlations. Default: 1.0
    fkt_num_wavevectors : int, optional
        Number of wavevectors for F(k,t) calculation. Default: 50
    fkt_reference_interval_ps : float, optional
        Reference interval for F(k,t) in ps. Default: 1.0
    fkt_max_references : int, optional
        Maximum F(k,t) reference states. Default: 10
    energy_output_period_ps : float, optional
        Energy tracker output period in ps. Default: 0.1
    fkt_output_period_ps : float, optional
        F(k,t) output period in ps. Default: 1.0
    gsd_output_period_ps : float, optional
        Trajectory output period in ps. Default: 50.0
    console_output_period_ps : float, optional
        Console output period in ps. Default: 1.0
        
    Attributes
    ----------
    sim : hoomd.Simulation
        HOOMD simulation object (available after setup)
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing (available after setup)
    cavity_force : CavityForce
        Cavity force object (available after setup if incavity=True)
    device_config : str
        Current device configuration ('CPU' or 'GPU')
        
    Methods
    -------
    run() : int
        Execute the complete simulation workflow
    setup_simulation() : object
        Set up HOOMD simulation state
    setup_force_parameters(dt) : List[hoomd.md.force.Force]
        Configure all force field components
    setup_thermostat_parameters(dt) : Tuple
        Configure thermostats for molecular and cavity particles
        
    Examples
    --------
    **Basic cavity simulation:**
    
    >>> from hoomd.cavitymd.simulation import CavityMDSimulation
    >>> 
    >>> sim = CavityMDSimulation(
    ...     job_dir="cavity_run",
    ...     replica=1,
    ...     freq=2000.0,  # 2000 cm⁻¹
    ...     couplstr=0.001,  # 1e-3 a.u.
    ...     incavity=True,
    ...     runtime_ps=1000.0
    ... )
    >>> 
    >>> # Run simulation
    >>> exit_code = sim.run()
    
    **Time-varying coupling experiment:**
    
    >>> sim = CavityMDSimulation(
    ...     job_dir="switching_experiment", 
    ...     replica=1,
    ...     freq=2000.0,
    ...     couplstr=0.001,
    ...     incavity=True,
    ...     runtime_ps=500.0,
    ...     switch_time_ps=100.0,  # Switch at 100 ps
    ...     dissipation=0.0001,    # Add damping
    ...     enable_energy_tracking=True
    ... )
    >>> exit_code = sim.run()
    
    **No-cavity control simulation:**
    
    >>> control_sim = CavityMDSimulation(
    ...     job_dir="no_cavity_control",
    ...     replica=1, 
    ...     freq=2000.0,
    ...     couplstr=0.0,
    ...     incavity=False,  # Disable cavity
    ...     runtime_ps=1000.0
    ... )
    >>> exit_code = control_sim.run()
    
    **GPU-accelerated simulation:**
    
    >>> gpu_sim = CavityMDSimulation(
    ...     job_dir="gpu_cavity",
    ...     replica=1,
    ...     freq=2000.0,
    ...     couplstr=0.001, 
    ...     incavity=True,
    ...     runtime_ps=1000.0,
    ...     device='GPU',
    ...     gpu_id=0
    ... )
    >>> exit_code = gpu_sim.run()
    
    **Simulation with momentum zeroing:**
    
    >>> momentum_sim = CavityMDSimulation(
    ...     job_dir="momentum_control",
    ...     replica=1,
    ...     freq=2000.0,
    ...     couplstr=0.001,
    ...     incavity=True,
    ...     runtime_ps=1000.0,
    ...     zero_momentum_enabled=True,
    ...     zero_momentum_period_ps=2.0  # Zero momentum every 2 ps
    ... )
    >>> exit_code = momentum_sim.run()
    
    Notes
    -----
    - Automatically handles cavity particle setup and validation
    - Supports both adaptive and fixed timestep integration
    - Provides comprehensive energy conservation monitoring
    - Compatible with SLURM array jobs and local execution
    - Outputs organized trajectories, logs, and analysis data
    - Maintains energy conservation during time-varying coupling
    - Optional momentum zeroing prevents center-of-mass drift during long simulations
    
    See Also
    --------
    hoomd.cavitymd.forces.CavityForce : For cavity-molecule interactions
    hoomd.cavitymd.analysis.EnergyTracker : For energy monitoring
    hoomd.cavitymd.variants.StepVariant : For time-varying parameters
    hoomd.bussi_reservoir.BussiReservoir : For Bussi thermostat
    
    References
    ----------
    For the underlying theory, see the documentation sections on single-mode
    cavity theory and time-varying coupling dynamics.
    """
    
    def __init__(self, 
                 job_dir: str, 
                 replica: int, 
                 freq: float, 
                 couplstr: float, 
                 incavity: bool, 
                 runtime_ps: float = 500.0,
                 input_gsd: str = 'molecular-0.gsd', 
                 frame: int = -1, 
                 name: str = 'prod', 
                 error_tolerance: float = 0.01,
                 temperature: float = 100.0, 
                 molecular_thermostat: str = 'bussi', 
                 cavity_thermostat: str = 'langevin',
                 cavity_damping_factor: float = 1.0, 
                 use_brownian_overdamped: bool = True, 
                 add_cavity_particle: bool = True,
                 finite_q: bool = False, 
                 molecular_thermostat_tau: float = 5.0, 
                 cavity_thermostat_tau: float = 5.0,
                 log_level: str = 'INFO', 
                 custom_log_file: Optional[str] = None,
                 enable_fkt: bool = True, 
                 fkt_kmag: float = 1.0, 
                 fkt_num_wavevectors: int = 50, 
                 fkt_reference_interval_ps: float = 1.0, 
                 fkt_max_references: int = 10,
                 enable_dipole_autocorr: bool = False,
                 dipole_reference_interval_ps: float = 1.0,
                 dipole_max_references: int = 10,
                 max_energy_output_time_ps: Optional[float] = None, 
                 enable_energy_tracking: bool = False, 
                 dt_fs: Optional[float] = None, 
                 device: str = 'CPU', 
                 gpu_id: int = 0,
                 energy_output_period_ps: float = 0.1, 
                 fkt_output_period_ps: float = 1.0, 
                 dipole_output_period_ps: float = 1.0,
                 gsd_output_period_ps: float = 50.0, 
                 console_output_period_ps: float = 1.0,
                 enable_text_output: bool = False, 
                 text_output_file: Optional[str] = None, 
                 truncate_gsd: bool = False, 
                 seed: Optional[int] = None, 
                 restart_velocities: bool = True,
                 switch_time_ps: Optional[float] = None, 
                 dissipation: float = 0.0,
                 decay_time_constant_ps: Optional[float] = None,
                 initial_fraction: float = 1e-5,
                 time_constant_ps: float = 50.0,
                 transition_time_ps: float = 0.5,
                 use_smooth_switching: bool = True,
                 zero_momentum_enabled: bool = False,
                 zero_momentum_period_ps: float = 1.0) -> None:
        """Initialize the CavityMDSimulation with simulation parameters."""
        self.job_dir = job_dir
        self.replica = replica
        self.freq = freq
        self.couplstr = couplstr
        self.incavity = incavity
        self.runtime_ps = runtime_ps
        self.input_gsd = input_gsd
        self.frame = frame
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
        
        # Time-varying parameters
        self.switch_time_ps = switch_time_ps
        self.dissipation = dissipation
        self.decay_time_constant_ps = decay_time_constant_ps
        
        # Adaptive timestep parameters
        self.initial_fraction = initial_fraction
        self.time_constant_ps = time_constant_ps
        
        # Momentum zeroing parameters
        self.zero_momentum_enabled = zero_momentum_enabled
        self.zero_momentum_period_ps = zero_momentum_period_ps
        
        # Logging parameters - simplified to console only
        self.log_level = log_level
        self.custom_log_file = custom_log_file
        
        # F(k,t) calculation parameters
        self.enable_fkt = enable_fkt
        self.fkt_kmag = fkt_kmag
        self.fkt_num_wavevectors = fkt_num_wavevectors
        self.fkt_reference_interval_ps = fkt_reference_interval_ps
        self.fkt_max_references = fkt_max_references
        
        # Dipole autocorrelation parameters
        self.enable_dipole_autocorr = enable_dipole_autocorr
        self.dipole_reference_interval_ps = dipole_reference_interval_ps
        self.dipole_max_references = dipole_max_references
        
        # Energy output limit parameter
        self.max_energy_output_time_ps = max_energy_output_time_ps
        
        # Energy tracking parameter
        self.enable_energy_tracking = enable_energy_tracking
        
        # Manual timestep parameter (in femtoseconds)
        self.dt_fs = dt_fs
        
        # Device configuration
        self.device = device.upper()
        self.gpu_id = gpu_id
        
        # Seed configuration
        self.seed = seed
        
        # Velocity restart configuration
        self.restart_velocities = restart_velocities
        
        # Physical constants
        self.kB = PhysicalConstants.KB_HARTREE_PER_K
        
        # Output period parameters - detailed control for each observable
        self.energy_output_period_ps = energy_output_period_ps
        self.fkt_output_period_ps = fkt_output_period_ps
        self.dipole_output_period_ps = dipole_output_period_ps
        self.gsd_output_period_ps = gsd_output_period_ps
        self.console_output_period_ps = console_output_period_ps
        
        # Text output parameters (deprecated but kept for compatibility)
        self.enable_text_output = enable_text_output
        self.text_output_file = text_output_file
        
        # GSD file control
        self.truncate_gsd = truncate_gsd
        
        # Initialize simulation components (will be set during setup)
        self.sim = None
        self.logger = None

    def get_array_module(self, *arrays):
        """
        Get the appropriate array module (numpy or cupy) based on the input arrays.
        
        This follows CuPy's guidelines for writing CPU/GPU agnostic code.
        If any input array is on GPU, returns cupy module, otherwise numpy.
        
        Parameters
        ----------
        *arrays : array-like
            Input arrays to check
            
        Returns
        -------
        module
            Either numpy or cupy module
        """
        if not HAS_CUPY or self.device == 'CPU':
            return np
            
        # Check if any array is a CuPy array or if we're using GPU device
        for array in arrays:
            if hasattr(array, 'device') and hasattr(array.device, 'id'):
                # This is a CuPy array
                return cp
        
        # Default to cupy if we're configured for GPU
        if self.device == 'GPU':
            return cp
        else:
            return np

    def to_device_array(self, array_data, target_arrays=None):
        """
        Convert array to the appropriate device format (NumPy for CPU, CuPy for GPU).
        
        This follows CuPy's guidelines using cupy.asarray() for device conversion.
        
        Parameters
        ----------
        array_data : array-like
            Input array data
        target_arrays : array-like, optional
            Target arrays to match device type. If provided, uses their device.
            
        Returns
        -------
        array
            Array on the appropriate device (NumPy for CPU, CuPy for GPU)
        """
        if not HAS_CUPY or self.device == 'CPU':
            # Always use NumPy for CPU
            return np.asarray(array_data)
        
        if target_arrays is not None:
            # Use get_array_module to determine the appropriate module
            xp = self.get_array_module(target_arrays)
            return xp.asarray(array_data)
        else:
            # Use device configuration
            if self.device == 'GPU':
                return cp.asarray(array_data)
            else:
                return np.asarray(array_data)

    def to_numpy(self, array_data):
        """
        Convert array to NumPy array, following CuPy guidelines using cupy.asnumpy().
        
        Parameters
        ----------
        array_data : array-like
            Input array data (NumPy or CuPy)
            
        Returns
        -------
        np.ndarray
            NumPy array
        """
        if HAS_CUPY and hasattr(array_data, 'get'):
            # CuPy array - use .get() method
            return array_data.get()
        elif HAS_CUPY and hasattr(cp, 'asnumpy'):
            # Use cupy.asnumpy() for robust conversion
            return cp.asnumpy(array_data)
        else:
            # Fallback to numpy.asarray
            return np.asarray(array_data)

    def calculate_physical_parameters(self):
        """Calculate physical parameters and unit conversions."""
        # Determine the actual timestep to use for calculations
        if self.error_tolerance <= 0 and self.dt_fs is not None:
            # Fixed timestep mode with user-specified timestep
            dt_ps = self.dt_fs / 1000.0  # Convert fs to ps
            self.log_info(f"Using user-specified timestep for calculations: {self.dt_fs} fs = {dt_ps} ps")
        else:
            # Adaptive timestep mode or no user-specified timestep - use default
            dt_ps = 0.0001  # timestep in ps (0.1 fs - reasonable default for adaptive mode)
            self.log_info(f"Using default timestep for calculations: {dt_ps} ps")
        
        runtime_real = self.runtime_ps  # runtime in ps
        
        # Calculate different output periods in timesteps using the correct dt_ps
        energy_period = max(1, int(self.energy_output_period_ps / dt_ps))
        fkt_period = max(1, int(self.fkt_output_period_ps / dt_ps))
        gsd_period = max(1, int(self.gsd_output_period_ps / dt_ps))
        console_period = max(1, int(self.console_output_period_ps / dt_ps))
        
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
        self.energy_period = energy_period
        self.fkt_period = fkt_period
        self.gsd_period = gsd_period
        self.console_period = console_period
        
        self.log_info(f"Time conversions:")
        self.log_info(f"  Timestep: {dt_ps} ps = {dt_au:.6f} a.u.")
        self.log_info(f"  Runtime: {self.runtime_ps:.1f} ps = {total_steps_needed} steps")
        self.log_info(f"  Steps per ps: {1.0/dt_ps:.1f}")
        self.log_info(f"Output periods (steps):")
        self.log_info(f"  Energy: {energy_period} ({self.energy_output_period_ps:.3f} ps)")
        self.log_info(f"  F(k,t): {fkt_period} ({self.fkt_output_period_ps:.3f} ps)")
        self.log_info(f"  GSD: {gsd_period} ({self.gsd_output_period_ps:.3f} ps)")
        self.log_info(f"  Console: {console_period} ({self.console_output_period_ps:.3f} ps)")
        
        return dt_au, total_steps_needed, energy_period, fkt_period, gsd_period, console_period

    def setup_device(self):
        """Setup the HOOMD device (CPU or GPU)."""
        if self.device == 'GPU':
            # Try different GPU initialization methods for different HOOMD versions
            try:
                # First try with gpu_ids parameter (newer HOOMD versions)
                device = hoomd.device.GPU(gpu_ids=[self.gpu_id])
            except TypeError:
                # Fall back to older parameter name
                device = hoomd.device.GPU(gpu_id=self.gpu_id)
            self.log_info(f"Initializing simulation on GPU {self.gpu_id}")
        elif self.device == 'CPU':
            device = hoomd.device.CPU()
            self.log_info("Initializing simulation on CPU")
        else:
            raise ValueError(f"Invalid device '{self.device}'. Must be 'CPU' or 'GPU'")
        
        return device

    def setup_logging(self):
        """Set up simplified logging configuration for the simulation."""
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
        
        # Always set up console logging (simplified)
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(getattr(logging, self.log_level.upper()))
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)
        
        # Log initial setup information
        self.log_info("="*60)
        self.log_info("CAVITY MD SIMULATION STARTED")
        self.log_info("="*60)
        self.log_info(f"Simulation: {self.name}-{self.replica}")
        self.log_info(f"Device: {self.device}" + (f" (GPU {self.gpu_id})" if self.device == 'GPU' else ""))
        self.log_info(f"Runtime: {self.runtime_ps} ps")
        self.log_info(f"Temperature: {self.temperature} K")
        self.log_info(f"Random seed: {self.seed if self.seed is not None else 'replica-based'}")
        self.log_info(f"Cavity coupling: {'Enabled' if self.incavity else 'Disabled'}")
        if self.incavity:
            self.log_info(f"  Frequency: {self.freq} cm^-1")
            self.log_info(f"  Coupling strength: {self.couplstr}")
            self.log_info(f"  Finite-q mode: {self.finite_q}")
        self.log_info(f"Molecular thermostat: {self.molecular_thermostat} (tau={self.molecular_thermostat_tau} ps)")
        if self.incavity:
            self.log_info(f"Cavity thermostat: {self.cavity_thermostat} (tau={self.cavity_thermostat_tau} ps)")
        self.log_info("="*60)
    
    def log_info(self, message):
        """Log an info message."""
        if hasattr(self, 'logger') and self.logger:
            self.logger.info(message)
        else:
            print(f"INFO: {message}", flush=True)
    
    def log_warning(self, message):
        """Log a warning message."""
        if hasattr(self, 'logger') and self.logger:
            self.logger.warning(message)
        else:
            print(f"WARNING: {message}", flush=True)
    
    def log_error(self, message):
        """Log an error message."""
        if hasattr(self, 'logger') and self.logger:
            self.logger.error(message)
        else:
            print(f"ERROR: {message}", flush=True)

    def setup_force_parameters(self, dt, rcut=15):
        """Set up force parameters for the simulation."""
        forces = []
        
        # Initialize cavity force reference
        self.cavity_force = None
        
        # Setup cavity force if requested
        if self.incavity:
            from .utils import PhysicalConstants
            omegac = self.freq / PhysicalConstants.HARTREE_TO_CM_MINUS1
            
            # Create time-varying coupling and dissipation if switch_time is specified
            if self.switch_time_ps is not None:
                # Import variants from plugin
                from .variants import StepVariant
                
                # Create step variants for instantaneous switching
                coupling_variant = StepVariant(
                    target_value=self.couplstr,
                    switch_time_ps=self.switch_time_ps,
                    time_tracker=self.time_tracker,
                    decay_time_constant_ps=self.decay_time_constant_ps
                )
                
                dissipation_variant = StepVariant(
                    target_value=self.dissipation,
                    switch_time_ps=self.switch_time_ps,
                    time_tracker=self.time_tracker,
                    decay_time_constant_ps=self.decay_time_constant_ps
                )
                
                # Update log messages based on whether decay is enabled
                if self.decay_time_constant_ps is not None:
                    self.log_info(f"Using step switching with exponential decay:")
                    self.log_info(f"  Switch time: {self.switch_time_ps} ps")
                    self.log_info(f"  Decay time constant: {self.decay_time_constant_ps} ps")
                    self.log_info(f"  Coupling: 0 → {self.couplstr} a.u. (with decay)")
                    self.log_info(f"  Dissipation: 0 → {self.dissipation} a.u. (with decay)")
                else:
                    self.log_info(f"Using instantaneous switching:")
                    self.log_info(f"  Switch time: {self.switch_time_ps} ps")
                    self.log_info(f"  Coupling: 0 → {self.couplstr} a.u. (instantaneous)")
                    self.log_info(f"  Dissipation: 0 → {self.dissipation} a.u. (instantaneous)")
                
                # Create cavity force with variants
                from .forces import CavityForce
                cavityforce = CavityForce(
                    kvector=np.array([0,0,1]), 
                    couplstr=coupling_variant, 
                    omegac=omegac,
                    dissipation=dissipation_variant
                )

                # Store variants and cavity force for later use in finite-q setup and console output
                self.coupling_variant = coupling_variant
                self.omegac = omegac
                self.cavity_force = cavityforce  # Store for console output access

            else:
                # Use constant values (backward compatibility)
                self.log_info(f"Using constant coupling and dissipation:")
                self.log_info(f"  Coupling: {self.couplstr} a.u.")
                self.log_info(f"  Dissipation: {self.dissipation} a.u.")
                
                from .forces import CavityForce
                cavityforce = CavityForce(
                    kvector=np.array([0,0,1]), 
                    couplstr=self.couplstr, 
                    omegac=omegac,
                    dissipation=self.dissipation
                )
                
                # Clear variants for later use but store cavity force
                self.coupling_variant = None
                self.omegac = omegac
                self.cavity_force = cavityforce  # Store for console output access
            
            forces.append(cavityforce)
        
        # Setup harmonic bonds
        harmonic = hoomd.md.bond.Harmonic()
        harmonic.params['O-O'] = dict(k=2*0.36602, r0=2.281655158)
        harmonic.params['N-N'] = dict(k=2*0.71625, r0=2.0743522177)
        forces.append(harmonic)
        
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

        # Setup long-range Coulomb interactions using PPPM method
        numpoints = 32
        order = 6
        short, long = hoomd.md.long_range.pppm.make_pppm_coulomb_forces(
            nlist=cell, resolution=[numpoints, numpoints, numpoints], 
            order=order, r_cut=rcut, alpha=0.0
        )
        forces.append(short)
        forces.append(long)
        
        return forces

    def setup_thermostat_parameters(self, dt):
        """Set up thermostat parameters for molecular and cavity systems."""
        from .utils import PhysicalConstants
        from hoomd.bussi_reservoir.thermostats import BussiReservoir as Bussi
        
        kT = self.kB * self.temperature
        molecular_filter = hoomd.filter.Type(['O', 'N'])  # Molecular particles only
        
        # Convert thermostat time constants from ps to atomic units
        molecular_tau_au = PhysicalConstants.ps_to_atomic_units(self.molecular_thermostat_tau)
        cavity_tau_au = PhysicalConstants.ps_to_atomic_units(self.cavity_thermostat_tau)
        
        self.log_info(f"Thermostat time constant conversions:")
        self.log_info(f"  Molecular: {self.molecular_thermostat_tau:.3f} ps = {molecular_tau_au:.6f} a.u.")
        self.log_info(f"  Cavity: {self.cavity_thermostat_tau:.3f} ps = {cavity_tau_au:.6f} a.u.")
        
        # Validate tau=0.0 with Langevin thermostats
        if self.molecular_thermostat.lower() == 'langevin' and self.molecular_thermostat_tau <= 0.0:
            raise ValueError(
                f"ERROR: Cannot use Langevin thermostat with molecular_thermostat_tau={self.molecular_thermostat_tau} ps.\n"
                f"Langevin dynamics requires tau > 0 since gamma = 1/tau.\n"
                f"For overdamped dynamics (tau → 0), use Brownian dynamics instead."
            )
        
        if self.incavity and self.cavity_thermostat.lower() == 'langevin' and self.cavity_thermostat_tau <= 0.0:
            raise ValueError(
                f"ERROR: Cannot use Langevin thermostat with cavity_thermostat_tau={self.cavity_thermostat_tau} ps.\n"
                f"Langevin dynamics requires tau > 0 since gamma = 1/tau.\n"
                f"For overdamped dynamics (tau → 0), use Brownian dynamics instead."
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
            self.log_info("Running molecular system with Bussi thermostat (NVT ensemble)")
            molecular_bussi = Bussi(kT=kT, tau=molecular_tau_au)
            molecular_method = hoomd.md.methods.ConstantVolume(filter=molecular_filter, thermostat=molecular_bussi)
            thermostat_refs['molecular_bussi'] = molecular_bussi
            self.log_info(f"Molecular Bussi thermostat configured: T = {self.temperature:.1f} K, kT = {kT:.6f} a.u., tau = {self.molecular_thermostat_tau:.3f} ps")
        elif self.molecular_thermostat.lower() == 'langevin':
            self.log_info("Running molecular system with Langevin thermostat (NVT ensemble)")
            molecular_gamma = PhysicalConstants.gamma_from_tau_ps(self.molecular_thermostat_tau)
            molecular_method = hoomd.md.methods.Langevin(filter=molecular_filter, kT=kT, default_gamma=molecular_gamma, tally_reservoir_energy=True)
            thermostat_refs['molecular_langevin'] = molecular_method
            self.log_info(f"Molecular Langevin thermostat configured: T = {self.temperature:.1f} K, kT = {kT:.6f} a.u.")
            self.log_info(f"  gamma = {molecular_gamma:.6f} a.u.^-1 (tau = {self.molecular_thermostat_tau:.3f} ps)")
        elif self.molecular_thermostat.lower() == 'none':
            self.log_info("Running molecular system without thermostat (NVE ensemble)")
            molecular_method = hoomd.md.methods.ConstantVolume(filter=molecular_filter)
        else:
            raise ValueError(f"Invalid molecular_thermostat option: {self.molecular_thermostat}")
        
        # Set up thermostat for cavity particle if present
        cavity_method = None
        if self.incavity:
            cavity_filter = hoomd.filter.Type(['L'])  # Cavity particle only
            
            if self.cavity_thermostat.lower() == 'langevin':
                self.log_info("Running cavity with Langevin thermostat")
                base_gamma = PhysicalConstants.gamma_from_tau_ps(self.cavity_thermostat_tau)
                cavity_gamma = self.cavity_damping_factor * base_gamma
                cavity_method = hoomd.md.methods.Langevin(filter=cavity_filter, 
                                                         kT=kT, default_gamma=cavity_gamma, tally_reservoir_energy=True)
                thermostat_refs['cavity_langevin'] = cavity_method
                self.log_info(f"Cavity Langevin thermostat configured: T = {self.temperature:.1f} K, kT = {kT:.6f} a.u.")
                self.log_info(f"  base_gamma = {base_gamma:.6f} a.u.^-1 (tau = {self.cavity_thermostat_tau:.3f} ps)")
                self.log_info(f"  effective_gamma = {cavity_gamma:.6f} a.u.^-1 (damping_factor = {self.cavity_damping_factor:.1f}x)")
            elif self.cavity_thermostat.lower() == 'bussi':
                self.log_info("Running cavity with Bussi thermostat")
                cavity_bussi = Bussi(kT=kT, tau=cavity_tau_au)
                cavity_method = hoomd.md.methods.ConstantVolume(filter=cavity_filter, thermostat=cavity_bussi)
                thermostat_refs['cavity_bussi'] = cavity_bussi
                self.log_info(f"Cavity Bussi thermostat configured: kT = {kT:.6f} a.u., tau = {self.cavity_thermostat_tau:.3f} ps")
            elif self.cavity_thermostat.lower() == 'none':
                self.log_info("Running cavity without thermostat (NVE ensemble)")
                cavity_method = hoomd.md.methods.ConstantVolume(filter=cavity_filter)
            else:
                raise ValueError(f"Invalid cavity_thermostat option: {self.cavity_thermostat}")
        
        return molecular_method, cavity_method, thermostat_refs

    def setup_integrator(self, forces, methods):
        """Configure the integrator with forces and integration methods."""
        # Setup integrator with initial dt
        integrator = hoomd.md.Integrator(dt=self.dt, forces=forces)
        self.sim.operations.integrator = integrator
        
        # Set integration methods (filter out None methods)
        valid_methods = [method for method in methods if method is not None]
        self.sim.operations.integrator.methods = valid_methods
        
        self.log_info(f"Integrator configured with initial dt = {self.dt:.6f} a.u. ({self.dt_ps:.6f} ps)")
        self.log_info(f"Number of integration methods: {len(valid_methods)}")
        
        # Setup CavityParticleDisplacer if using finite-q mode with time-varying coupling
        if self.finite_q and self.switch_time_ps is not None:
            self.log_info("Setting up CavityParticleDisplacer for finite-q mode with time-varying coupling")
            
            # Create Python custom action for cavity particle displacement
            class CavityParticleDisplacer(hoomd.custom.Action):
                def __init__(self, coupling_variant, omegac, temperature, kB, parent_sim):
                    self.coupling_variant = coupling_variant
                    self.omegac = omegac
                    self.temperature = temperature
                    self.kB = kB
                    self.phmass = 1.0  # Photon mass
                    self.K = self.phmass * self.omegac * self.omegac
                    self.has_run = False
                    self.parent_sim = parent_sim  # Reference to parent simulation for CPU/GPU methods
                    
                def act(self, timestep):
                    if self.has_run:
                        return
                    
                    # Get current coupling strength
                    g_current = self.coupling_variant(timestep)
                    
                    # Check if coupling has switched from 0 to non-zero
                    if g_current == 0.0:
                        return  # Still zero, nothing to do
                    
                    # Coupling has switched ON! Perform displacement
                    print(f"CavityParticleDisplacer: Coupling switched ON at timestep {timestep}, g = {g_current}", flush=True)
                    
                    # Get system snapshot
                    snap = self._state.get_snapshot()
                    
                    if snap.communicator.rank == 0:
                        # Find cavity particle (type 'L', typeid = 2)
                        cavity_indices = [i for i in range(snap.particles.N) if snap.particles.typeid[i] == 2]
                        
                        if len(cavity_indices) == 0:
                            print("Warning: No cavity particle found for displacement", flush=True)
                            self.has_run = True
                            return
                        
                        cavity_idx = int(cavity_indices[0])  # Use first cavity particle
                        
                        # Get appropriate array module for this device
                        xp = self.parent_sim.get_array_module(snap.particles.position)
                        
                        # Ensure all arrays are on the same device before operations
                        positions = self.parent_sim.to_device_array(snap.particles.position, snap.particles.position)
                        images = self.parent_sim.to_device_array(snap.particles.image, snap.particles.position)
                        box_L = self.parent_sim.to_device_array(
                            np.array(snap.configuration.box[:3]), 
                            snap.particles.position
                        )
                        
                        unwrapped_pos = positions + images * box_L
                        
                        # Compute molecular dipole moment (excluding cavity particle)
                        # Convert to NumPy for calculations
                        unwrapped_pos_np = self.parent_sim.to_numpy(unwrapped_pos)
                        charges_np = self.parent_sim.to_numpy(snap.particles.charge)
                        
                        dipole = np.array([0.0, 0.0, 0.0])
                        for i in range(snap.particles.N):
                            if i != cavity_idx:
                                dipole += charges_np[i] * unwrapped_pos_np[i]
                        
                        # Only use x,y components for displacement
                        dipole_xy = np.array([dipole[0], dipole[1], 0.0])
                        
                        # Calculate new equilibrium position
                        q_eq = -(g_current / self.K) * dipole_xy
                        
                        # Get current cavity position
                        q_old = unwrapped_pos_np[cavity_idx]
                        
                        print(f"  Dipole moment (xy): {dipole_xy}", flush=True)
                        print(f"  Old cavity position: {q_old}", flush=True)
                        print(f"  New equilibrium position: {q_eq}", flush=True)
                        
                        # Create new position (preserve z coordinate)
                        new_pos_unwrapped = np.array([q_eq[0], q_eq[1], q_old[2]])
                        
                        # Convert box lengths to NumPy for calculations
                        box_L_np = self.parent_sim.to_numpy(box_L)
                        
                        # Wrap position back into box
                        new_image = np.floor((new_pos_unwrapped + box_L_np/2) / box_L_np).astype(int)
                        new_pos_wrapped = new_pos_unwrapped - new_image * box_L_np
                        
                        # Convert to device-appropriate arrays for assignment
                        new_pos_device = self.parent_sim.to_device_array(new_pos_wrapped, snap.particles.position)
                        new_image_device = self.parent_sim.to_device_array(new_image, snap.particles.image)
                        
                        # Update particle position and image
                        # CRITICAL: HOOMD snapshots ALWAYS expect NumPy arrays, regardless of device
                        # Convert all arrays to NumPy before assignment
                        snap.particles.position[cavity_idx] = self.parent_sim.to_numpy(new_pos_device)
                        snap.particles.image[cavity_idx] = self.parent_sim.to_numpy(new_image_device)
                        
                        displacement_xy = np.linalg.norm(new_pos_unwrapped[:2] - q_old[:2])
                        print(f"  Displacement magnitude (xy): {displacement_xy}", flush=True)
                        print(f"  Final position: {new_pos_wrapped}", flush=True)
                        print(f"  Using array module: {xp.__name__}", flush=True)
                    
                    # Set updated snapshot back to system
                    self._state.set_snapshot(snap)
                    
                    print("CavityParticleDisplacer: Displacement completed and system state updated", flush=True)
                    self.has_run = True
            
            # Create and add the custom updater
            cavity_displacer_action = CavityParticleDisplacer(
                coupling_variant=self.coupling_variant,
                omegac=self.omegac,
                temperature=self.temperature,
                kB=self.kB,
                parent_sim=self  # Pass reference to parent simulation for CPU/GPU methods
            )
            
            displacer_updater = hoomd.update.CustomUpdater(
                action=cavity_displacer_action,
                trigger=hoomd.trigger.Periodic(1)
            )
            self.sim.operations.updaters.append(displacer_updater)
            
            self.log_info("CavityParticleDisplacer configured as Python custom action")
        else:
            self.log_info("CavityParticleDisplacer not needed (constant coupling or q=0 mode)")

    def thermalize_system(self):
        """Initialize particle velocities using HOOMD State for molecular particles and manual Maxwell-Boltzmann for cavity particle."""
        
        if not self.restart_velocities:
            self.log_info("Velocity restart disabled - keeping existing velocities from GSD file")
            # Check if we need to verify existing velocities
            with self.get_local_snapshot() as snap:
                # Get masses and velocities
                masses = self.to_numpy(snap.particles.mass)
                velocities = self.to_numpy(snap.particles.velocity)
                typeid_array = self.to_numpy(snap.particles.typeid)
                
                # Filter for molecular particles only (typeid != 2)
                molecular_mask = typeid_array != 2
                
                if np.any(molecular_mask):
                    # Extract molecular particles data
                    mol_masses = masses[molecular_mask]
                    mol_velocities = velocities[molecular_mask]
                    
                    # Calculate kinetic energy: KE = 0.5 * m * v²
                    v_squared = np.sum(mol_velocities**2, axis=1)
                    kinetic_energy = 0.5 * mol_masses * v_squared
                    total_ke = np.sum(kinetic_energy)
                    
                    # Calculate temperature: T = (2/3) * KE / (N * kB)
                    n_mol_particles = np.sum(molecular_mask)
                    degrees_of_freedom = 3 * n_mol_particles
                    current_temp = (2.0) * total_ke / (degrees_of_freedom * self.kB)
                    
                    self.log_info(f"Initial kinetic temperature from existing velocities: {current_temp:.2f} K")
                    
                    # Exit if temperature is too high (> 10K)
                    if current_temp-100.0 > 10.0:
                        self.log_info(f"ERROR: Initial temperature too high ({current_temp:.2f} K > 10 K). Exiting.")
                        import sys
                        sys.exit(1)
                    
                    # FIX: Check for zero velocities when using Bussi thermostat
                    # Bussi thermostat requires non-zero initial momenta to function properly
                    velocity_magnitude_threshold = 1e-12  # Very small threshold for "zero" velocities
                    max_velocity_magnitude = np.max(np.sqrt(np.sum(mol_velocities**2, axis=1)))
                    
                    # Also check cavity particle velocities if cavity thermostat is Bussi
                    cavity_max_velocity = 0.0
                    if self.incavity:
                        cavity_mask = typeid_array == 2
                        if np.any(cavity_mask):
                            cavity_velocities = velocities[cavity_mask]
                            cavity_max_velocity = np.max(np.sqrt(np.sum(cavity_velocities**2, axis=1)))
                            max_velocity_magnitude = max(max_velocity_magnitude, cavity_max_velocity)
                    
                    is_bussi_molecular = (self.molecular_thermostat.lower() == 'bussi')
                    is_bussi_cavity = (self.cavity_thermostat.lower() == 'bussi') if self.incavity else False
                    
                    if max_velocity_magnitude < velocity_magnitude_threshold and (is_bussi_molecular or is_bussi_cavity):
                        self.log_warning("CRITICAL: Zero velocities detected with Bussi thermostat!")
                        self.log_warning("Bussi thermostat requires non-zero initial momenta.")
                        self.log_warning("Forcing velocity thermalization to prevent simulation failure...")
                        
                        # Force thermalization despite restart_velocities=False
                        self._force_thermalization_for_bussi()
                        return
                        
            return
        
        kT = self.kB * self.temperature
        
        # Set numpy seed for reproducible thermalization if seed is provided
        if self.seed is not None:
            np.random.seed(self.seed + 1)  # Use seed+1 to differentiate from HOOMD seed
            self.log_info(f"Using seed {self.seed + 1} for thermalization randomness")
        
        molecular_filter = hoomd.filter.Type(['O', 'N'])  # Molecular particles only
        
        if self.incavity:
            # Use HOOMD's State.thermalize_particle_momenta for molecular system
            self.sim.state.thermalize_particle_momenta(kT=kT, filter=molecular_filter)
            self.log_info("Thermalized molecular particles using HOOMD State object")
            
            # Manual thermalization for cavity particle
            self._thermalize_cavity_particle_manually(kT)
            
        else:
            # No cavity particle, thermalize all particles
            self.sim.state.thermalize_particle_momenta(kT=kT, filter=hoomd.filter.All())
            self.log_info("Thermalized all molecular particles")
        # Check temperature of the molecular system
        # Get current snapshot to analyze velocities
        with self.get_local_snapshot() as snap:
            # Get masses and velocities
            masses = self.to_numpy(snap.particles.mass)
            velocities = self.to_numpy(snap.particles.velocity)
            typeid_array = self.to_numpy(snap.particles.typeid)
            
            # Filter for molecular particles only (typeid != 2)
            molecular_mask = typeid_array != 2
            
            if np.any(molecular_mask):
                # Extract molecular particles data
                mol_masses = masses[molecular_mask]
                mol_velocities = velocities[molecular_mask]
                
                # Calculate kinetic energy: KE = 0.5 * m * v²
                v_squared = np.sum(mol_velocities**2, axis=1)
                kinetic_energy = 0.5 * mol_masses * v_squared
                total_ke = np.sum(kinetic_energy)
                
                # Calculate temperature: T = (2/3) * KE / (N * kB)
                n_mol_particles = np.sum(molecular_mask)
                degrees_of_freedom = 3 * n_mol_particles
                current_temp = (2.0 ) * total_ke / (degrees_of_freedom  * self.kB)
                
                self.log_info(f"Checking molecular system temperature after thermalization:")
                self.log_info(f"  Number of molecular particles: {n_mol_particles}")
                self.log_info(f"  Total kinetic energy: {total_ke:.6f} a.u.")
                self.log_info(f"  Current temperature: {current_temp:.1f} K")
                
                # Exit if temperature is too high
                if current_temp-100.0 > 10.0:
                    self.log_info(f"ERROR: Molecular system temperature ({current_temp:.1f} K) exceeds 10 K threshold!")
                    self.log_info("Terminating simulation for safety.")
                    import sys
                    sys.exit(1)
                else:
                    self.log_info(f"Temperature check passed: {current_temp:.1f} K is within acceptable range (≤ 10 K)")
            else:
                self.log_info("WARNING: No molecular particles found for temperature check!")
        self.log_info("Thermalization completed - velocities properly initialized")

    def _thermalize_cavity_particle_manually(self, kT):
        """Manually thermalize cavity particle using Maxwell-Boltzmann distribution."""
        with self.get_local_snapshot() as snap:
            # Get the appropriate array module using CuPy's guidelines
            xp = self.get_array_module(snap.particles.velocity)
            
            # Convert to numpy array for CPU/GPU agnostic operations
            typeid_array = self.to_numpy(snap.particles.typeid)
            cavity_mask = typeid_array == 2
            
            # Find cavity particle indices using proper array handling
            if np.any(cavity_mask):
                cavity_indices = np.where(cavity_mask)[0]
                cavity_idx = cavity_indices[0]  # Get first cavity particle
                
                # Maxwell-Boltzmann distribution: each component has variance kT/m
                # With mass = 1.0 a.u., std dev per component = sqrt(kT)
                # Generate velocities using NumPy (always works)
                cavity_velocity_np = np.random.normal(0.0, np.sqrt(kT), size=3)
                
                # Convert to appropriate array type for the current device
                if xp == cp and HAS_CUPY:
                    # GPU case: convert NumPy to CuPy for assignment
                    cavity_velocity_device = cp.asarray(cavity_velocity_np)
                else:
                    # CPU case: use NumPy array directly
                    cavity_velocity_device = cavity_velocity_np
                
                # Calculate expected kinetic energy and temperature for logging
                expected_ke = 0.5 * 1.0 * np.sum(cavity_velocity_np**2)  # KE = (1/2) * m * v²
                expected_temp = (2.0/3.0) * expected_ke / self.kB  # T = (2/3) * KE / kB for 3D

                self.log_info(f"Cavity particle manually thermalized:")
                self.log_info(f"  Target temperature: {self.temperature:.1f} K")
                self.log_info(f"  Initial velocity: {cavity_velocity_np}")
                self.log_info(f"  Expected KE: {expected_ke:.6f} a.u.")
                self.log_info(f"  Expected temperature: {expected_temp:.1f} K")
                self.log_info(f"  Using array module: {xp.__name__}")
                
                # Assign velocity using device-appropriate array type
                snap.particles.velocity[cavity_idx] = cavity_velocity_device
                
            else:
                self.log_info("WARNING: No cavity particle found for thermalization!")

    def setup_trackers_and_loggers(self):
        """Set up comprehensive tracking and logging objects for the simulation."""
        
        # Create custom performance tracker
        self.performance_tracker = PerformanceTracker(self.sim, self.runtime_ps, self.time_tracker)
        self.sim.operations.updaters.append(hoomd.update.CustomUpdater(
            action=self.performance_tracker, trigger=hoomd.trigger.Periodic(10)
        ))
        
        # Set up adaptive timestep updater if error_tolerance is positive
        if hasattr(self, 'error_tolerance') and self.error_tolerance > 0:
            self.log_info(f"Setting up adaptive timestep updater (error_tolerance = {self.error_tolerance})")
            
            # Get required parameters with defaults
            cavity_damping_factor = getattr(self, 'cavity_damping_factor', 1.0)
            molecular_thermostat_tau = getattr(self, 'molecular_thermostat_tau', 5.0)
            cavity_thermostat_tau = getattr(self, 'cavity_thermostat_tau', 5.0)
            switch_time_ps = getattr(self, 'switch_time_ps', None)
            
            self.adaptive_action = AdaptiveTimestepUpdater(
                state=self.sim.state,
                integrator=self.sim.operations.integrator,
                error_tolerance=self.error_tolerance,
                time_constant_ps=self.time_constant_ps,
                initial_fraction=self.initial_fraction,
                adaptiveerror=True,
                cavity_damping_factor=cavity_damping_factor,
                molecular_thermostat_tau=molecular_thermostat_tau,
                cavity_thermostat_tau=cavity_thermostat_tau,
                time_tracker=self.time_tracker,
                switch_time_ps=switch_time_ps,
                timestep_change_threshold=0.1,  # Only update if change is >10%
                max_timestep_change_factor=1.5,  # Limit maximum change to 10%
                shock_dampening_factor=self.initial_fraction,  # Use the actual initial_fraction parameter
                shock_dampening_enabled=(switch_time_ps is not None),  # Enable whenever there's a switch time
                shock_dampening_time_constant_ps=self.time_constant_ps  # Use the same time constant for shock dampening
            )
            
            # Add adaptive updater - use a less frequent trigger to reduce computational overhead
            # The AdaptiveTimestepUpdater now has internal logic to prevent frequent updates
            adaptive_timestep_trigger_period = 1  # Check every 100 steps instead of every energy_period
            adaptive_updater = hoomd.update.CustomUpdater(
                action=self.adaptive_action,
                trigger=hoomd.trigger.Periodic(adaptive_timestep_trigger_period)
            )
            self.sim.operations.updaters.append(adaptive_updater)
            self.log_info(f"Adaptive timestep updater enabled (check every {adaptive_timestep_trigger_period} steps)")
        else:
            self.adaptive_action = None
            self.log_info("Fixed timestep mode (error_tolerance = 0)")
        
        # Create status tracker for performance metrics (kept for compatibility)
        self.status = Status(self.sim, self.runtime_ps, self.time_tracker)
        
        # Create timestep formatter for fs display
        self.timestep_formatter = TimestepFormatter(self.sim.operations.integrator)
        
        # Create comprehensive logger
        logger = hoomd.logging.Logger(categories=['scalar', 'string'])
        
        # Basic simulation quantities
        logger.add(self.sim, quantities=['timestep', 'tps'])
        
        # Time and performance information
        logger[('Time', 'elapsed_ps')] = (self.time_tracker, 'elapsed_time', 'scalar')
        logger[('Performance', 'ns_per_day')] = (self.performance_tracker, 'ns_per_day', 'string')
        logger[('Performance', 'eta')] = (self.performance_tracker, 'eta_remaining', 'string')
        logger[('Timestep', 'dt_fs')] = (self.timestep_formatter, 'dt_fs', 'scalar')
        
        # Add adaptive timestep logging if enabled
        if hasattr(self, 'adaptive_action') and self.adaptive_action is not None:
            logger[('Adaptive', 'error_tolerance')] = (self.adaptive_action, 'error_tolerance', 'scalar')
        
        # Add thermodynamic quantities
        # Add integrator energies - DISABLED: integrator doesn't expose these quantities
        # logger.add(self.sim.operations.integrator, quantities=['kinetic_energy', 'potential_energy'])
        
        # Add system thermodynamics
        all_filter = hoomd.filter.All()
        self.thermo = hoomd.md.compute.ThermodynamicQuantities(filter=all_filter)
        self.sim.operations.computes.append(self.thermo)
        
        # Create molecular thermodynamics computer
        molecular_filter = hoomd.filter.Type(['O', 'N'])
        self.molecular_thermo = hoomd.md.compute.ThermodynamicQuantities(filter=molecular_filter)
        self.sim.operations.computes.append(self.molecular_thermo)
        
        if hasattr(self, 'incavity') and self.incavity:
            # Add cavity particle thermodynamics
            cavity_filter = hoomd.filter.Type(['L'])
            self.cavity_thermo = hoomd.md.compute.ThermodynamicQuantities(filter=cavity_filter)
            self.sim.operations.computes.append(self.cavity_thermo)
        
        # Run a single step to initialize all compute objects before setting up logger
        # This ensures all properties are accessible when we add them to the logger
        self.sim.run(1, write_at_start=False)
        
        # Now add thermodynamic quantities to logger after compute objects are initialized
        logger[('System', 'kinetic_energy')] = (self.thermo, 'kinetic_energy', 'scalar')
        logger[('System', 'potential_energy')] = (self.thermo, 'potential_energy', 'scalar')
        logger[('System', 'temperature')] = (self.thermo, 'kinetic_temperature', 'scalar')
        logger[('System', 'pressure')] = (self.thermo, 'pressure', 'scalar')
        logger[('System', 'volume')] = (self.thermo, 'volume', 'scalar')
        logger[('System', 'N_particles')] = (self.thermo, 'num_particles', 'scalar')
        
        logger[('Molecular', 'temperature')] = (self.molecular_thermo, 'kinetic_temperature', 'scalar')
        logger[('Molecular', 'kinetic_energy')] = (self.molecular_thermo, 'kinetic_energy', 'scalar')
        
        if hasattr(self, 'incavity') and self.incavity:
            logger[('Cavity', 'temperature')] = (self.cavity_thermo, 'kinetic_temperature', 'scalar')
            logger[('Cavity', 'kinetic_energy')] = (self.cavity_thermo, 'kinetic_energy', 'scalar')
        
        try:
            # Track thermostat energies
            method_count = 0
            for i, method in enumerate(self.sim.operations.integrator.methods):
                if hasattr(method, 'thermostat') and method.thermostat is not None:
                    if hasattr(method.thermostat, 'energy'):
                        try:
                            logger[('Thermostat', f'energy_{i}')] = (method.thermostat, 'energy', 'scalar')
                            method_count += 1
                            self.log_info(f"  Added thermostat energy tracking for method {i}")
                        except Exception as e:
                            self.log_warning(f"  Failed to add thermostat energy for method {i}: {e}")
                    
                    # Track detailed BussiReservoir energies if available
                    if hasattr(method.thermostat, 'total_reservoir_energy'):
                        try:
                            logger[('Thermostat', f'total_reservoir_energy_{i}')] = (method.thermostat, 'total_reservoir_energy', 'scalar')
                            logger[('Thermostat', f'translational_energy_{i}')] = (method.thermostat, 'reservoir_energy_translational', 'scalar')
                            logger[('Thermostat', f'rotational_energy_{i}')] = (method.thermostat, 'reservoir_energy_rotational', 'scalar')
                            self.log_info(f"  Added detailed BussiReservoir energy tracking for method {i}")
                        except Exception as e:
                            self.log_warning(f"  Failed to add detailed thermostat energies for method {i}: {e}")
            
            # Track individual force energies
            force_count = 0
            for i, force in enumerate(self.sim.operations.integrator.forces):
                force_name = type(force).__name__.lower()
                try:
                    if hasattr(force, 'energy'):
                        logger[('Force', f'{force_name}_energy')] = (force, 'energy', 'scalar')
                        force_count += 1
                        self.log_info(f"  Added energy tracking for {force_name}")
                    
                    # Special handling for CavityForce with component energies
                    if hasattr(force, 'harmonic_energy'):
                        logger[('Cavity', 'harmonic_energy')] = (force, 'harmonic_energy', 'scalar')
                        logger[('Cavity', 'coupling_energy')] = (force, 'coupling_energy', 'scalar')
                        logger[('Cavity', 'dipole_self_energy')] = (force, 'dipole_self_energy', 'scalar')
                        self.log_info(f"  Added detailed cavity energy components")
                except Exception as e:
                    self.log_warning(f"  Failed to add energy tracking for {force_name}: {e}")
            
            self.log_info(f"Energy tracking setup completed: {method_count} methods, {force_count} forces")
            
            # Set up energy tracker from plugin for proper reservoir energy tracking
            try:
                name = getattr(self, 'name', 'sim')
                replica = getattr(self, 'replica', 0)
                energy_output_period_ps = getattr(self, 'energy_output_period_ps', 0.1)
                max_energy_output_time_ps = getattr(self, 'max_energy_output_time_ps', None)
                
                # Fix duplicate "energy" in filename - EnergyTracker adds "_energy_tracker.txt" automatically
                output_prefix = f'{name}-{replica}'
                actual_energy_filename = f'{output_prefix}_energy_tracker.txt'
                self.log_info(f"Setting up energy tracker with output file: {actual_energy_filename}")
                
                self.log_info(f"Energy tracker configuration:")
                self.log_info(f"  Using time-based output: {energy_output_period_ps:.3f} ps")
                if max_energy_output_time_ps:
                    self.log_info(f"  Max output time: {max_energy_output_time_ps:.3f} ps")
                else:
                    self.log_info(f"  No max time limit")
                
                # Prepare individual force objects for EnergyTracker
                force_objects = {}
                thermostat_objects = {}
                
                for force in self.sim.operations.integrator.forces:
                    force_name = type(force).__name__.lower()
                    if 'cavity' in force_name:
                        force_objects['cavity'] = force
                    elif 'lj' in force_name or 'lennard' in force_name:
                        force_objects['lj'] = force
                    elif 'harmonic' in force_name or 'bond' in force_name:
                        force_objects['harmonic'] = force
                    elif 'ewald' in force_name:
                        force_objects['ewald_short'] = force  # Ewald is the short-range PPPM force
                    elif 'coulomb' in force_name:
                        force_objects['ewald_long'] = force   # Coulomb is the long-range PPPM force
                
                # Extract thermostat objects for EnergyTracker
                for i, method in enumerate(self.sim.operations.integrator.methods):
                    method_name = type(method).__name__.lower()
                    
                    # Check for Langevin methods (have reservoir_energy)
                    if 'langevin' in method_name and hasattr(method, 'reservoir_energy'):
                        # Determine if this is molecular or cavity based on filter
                        if hasattr(method, 'filter'):
                            filter_types = getattr(method.filter, '_types', [])
                            if 'L' in filter_types:
                                thermostat_objects['langevin_cavity'] = method
                            else:
                                thermostat_objects['langevin_molecular'] = method
                    
                    # Check for Bussi thermostats
                    if hasattr(method, 'thermostat') and method.thermostat is not None:
                        thermostat_type = type(method.thermostat).__name__.lower()
                        if 'bussi' in thermostat_type:
                            # Determine if this is molecular or cavity based on filter
                            if hasattr(method, 'filter'):
                                filter_types = getattr(method.filter, '_types', [])
                                if 'L' in filter_types:
                                    thermostat_objects['bussi_cavity'] = method.thermostat
                                else:
                                    thermostat_objects['bussi_molecular'] = method.thermostat
                
                # Get kinetic trackers 
                kinetic_tracker = None  # Use internal kinetic computation
                
                self.log_info(f"Energy tracker configuration:")
                self.log_info(f"  Force objects: {list(force_objects.keys())}")
                self.log_info(f"  Thermostat objects: {list(thermostat_objects.keys())}")
                self.log_info(f"  Using internal kinetic computation (no external tracker needed)")
                
                # Use time-based output period for accurate timing
                self.energy_tracker = EnergyTracker(
                    simulation=self.sim,
                    components=['kinetic', 'harmonic', 'lj', 'ewald_short', 'ewald_long', 'cavity'],
                    force_objects=force_objects,
                    thermostat_objects=thermostat_objects,
                    kinetic_tracker=kinetic_tracker,  # Use internal kinetic computation
                    time_tracker=self.time_tracker,
                    output_prefix=output_prefix,
                    output_period_ps=energy_output_period_ps,  # Use time-based output!
                    max_time_ps=max_energy_output_time_ps,  # Use time-based limit
                    compute_temperature=True,
                    track_reservoirs=True,
                    verbose='quiet'  # Reduce verbose output - use 'verbose' for full debug output
                )
                
                # Add energy tracker to simulation - trigger period doesn't matter since it uses internal timing
                energy_updater = hoomd.update.CustomUpdater(
                    action=self.energy_tracker,
                    trigger=hoomd.trigger.Periodic(1)  # Check every step but tracker handles timing internally
                )
                self.sim.operations.updaters.append(energy_updater)
                
                self.log_info(f"Energy tracker setup completed with time-based output:")
                self.log_info(f"  Output period: {energy_output_period_ps:.3f} ps (accurate timing)")
                self.log_info(f"  Tracker handles timing internally using ElapsedTimeTracker")
                if max_energy_output_time_ps:
                    self.log_info(f"  Output limited to {max_energy_output_time_ps:.3f} ps")
                
            except Exception as e:
                self.log_error(f"Failed to setup energy tracker: {e}")
                import traceback
                self.log_error("Full traceback:")
                for line in traceback.format_exc().split('\n'):
                    if line.strip():
                        self.log_error(line)
                self.energy_tracker = None
            
        except Exception as e:
            self.log_warning(f"Could not complete energy tracking setup: {str(e)}")
            self.log_warning(f"Error details: {type(e).__name__}: {str(e)}")
            self.log_warning(f"Full traceback:")
            import traceback
            traceback.print_exc()
        
        # Set up F(k,t) density correlation tracker if enabled
        enable_fkt = getattr(self, 'enable_fkt', False)
        if enable_fkt:
            try:
                fkt_kmag = getattr(self, 'fkt_kmag', 1.0)
                fkt_num_wavevectors = getattr(self, 'fkt_num_wavevectors', 50)
                fkt_reference_interval_ps = getattr(self, 'fkt_reference_interval_ps', 1.0)
                fkt_max_references = getattr(self, 'fkt_max_references', 10)
                fkt_output_period_ps = getattr(self, 'fkt_output_period_ps', 1.0)
                name = getattr(self, 'name', 'sim')
                replica = getattr(self, 'replica', 0)
                
                self.log_info("Setting up F(k,t) density correlation tracker")
                self.log_info(f"  k magnitude: {fkt_kmag:.2f}")
                self.log_info(f"  Number of wavevectors: {fkt_num_wavevectors}")
                self.log_info(f"  Reference interval: {fkt_reference_interval_ps:.3f} ps")
                self.log_info(f"  Max references: {fkt_max_references}")
                self.log_info(f"  Output period: {fkt_output_period_ps:.3f} ps")
                
                # Use time-based intervals for adaptive timestep compatibility
                fkt_trigger_period = 1  # Check every step for best timing
                
                self.log_info(f"  Using time-based reference intervals for adaptive timestep compatibility")
                self.log_info(f"  Reference interval: {fkt_reference_interval_ps:.3f} ps")
                self.log_info(f"  Trigger period: {fkt_trigger_period} step")
                
                # Create density correlation tracker with time-based reference interval
                self.density_corr_tracker = FieldAutocorrelationTracker(
                    simulation=self.sim,
                    observable="density_correlation",
                    time_tracker=self.time_tracker,
                    output_period_ps=fkt_output_period_ps,  # Use time-based output period!
                    output_prefix=f'{name}-{replica}',
                    reference_interval_ps=fkt_reference_interval_ps,  # Use time-based interval
                    max_references=fkt_max_references,
                    kmag=fkt_kmag,
                    num_wavevectors=fkt_num_wavevectors
                )
                
                # Add F(k,t) tracker to simulation - trigger period doesn't matter since it uses internal timing
                fkt_updater = hoomd.update.CustomUpdater(
                    action=self.density_corr_tracker,
                    trigger=hoomd.trigger.Periodic(1)  # Check every step but tracker handles timing internally
                )
                self.sim.operations.updaters.append(fkt_updater)
                
                # Add F(k,t) data to logger
                logger[('F(k,t)', 'current_autocorr')] = (self.density_corr_tracker, 'current_autocorr', 'scalar')
                
                self.log_info("F(k,t) tracker successfully enabled with time-based output:")
                self.log_info(f"  Output period: {fkt_output_period_ps:.3f} ps (accurate timing)")
                self.log_info(f"  Reference interval: {fkt_reference_interval_ps:.3f} ps")
                self.log_info(f"  Tracker handles timing internally using ElapsedTimeTracker")
                
            except Exception as e:
                self.log_warning(f"Could not set up F(k,t) tracker: {str(e)}")
                self.density_corr_tracker = None
        else:
            self.density_corr_tracker = None
            self.log_info("F(k,t) tracking disabled")
        
        # ===== DIPOLE AUTOCORRELATION TRACKER SETUP =====
        if self.enable_dipole_autocorr:
            try:
                self.log_info("Setting up dipole autocorrelation tracker...")
                self.log_info(f"  Reference interval: {self.dipole_reference_interval_ps:.3f} ps")
                self.log_info(f"  Max references: {self.dipole_max_references}")
                self.log_info(f"  Output period: {self.dipole_output_period_ps:.3f} ps")
                
                # Use time-based intervals for adaptive timestep compatibility
                dipole_trigger_period = 1  # Check every step for best timing
                
                self.log_info(f"  Using time-based reference intervals for adaptive timestep compatibility")
                self.log_info(f"  Reference interval: {self.dipole_reference_interval_ps:.3f} ps")
                self.log_info(f"  Trigger period: {dipole_trigger_period} step")
                
                # Create dipole autocorrelation tracker with time-based reference interval
                self.dipole_autocorr_tracker = AutocorrelationTracker(
                    simulation=self.sim,
                    observable="dipole",
                    time_tracker=self.time_tracker,
                    output_period_ps=self.dipole_output_period_ps,  # Use time-based output period!
                    output_prefix=f'{name}-{replica}_dipole_autocorr',
                    reference_interval_ps=self.dipole_reference_interval_ps,  # Use time-based interval
                    max_references=self.dipole_max_references
                )
                
                # Add dipole autocorrelation tracker to simulation - trigger period doesn't matter since it uses internal timing
                dipole_updater = hoomd.update.CustomUpdater(
                    action=self.dipole_autocorr_tracker,
                    trigger=hoomd.trigger.Periodic(1)  # Check every step but tracker handles timing internally
                )
                self.sim.operations.updaters.append(dipole_updater)
                
                # Add dipole autocorrelation data to logger
                logger[('Dipole', 'current_autocorr')] = (self.dipole_autocorr_tracker, 'current_autocorr', 'scalar')
                
                self.log_info("Dipole autocorrelation tracker successfully enabled with time-based output:")
                self.log_info(f"  Output period: {self.dipole_output_period_ps:.3f} ps (accurate timing)")
                self.log_info(f"  Reference interval: {self.dipole_reference_interval_ps:.3f} ps")
                self.log_info(f"  Tracker handles timing internally using ElapsedTimeTracker")
                
            except Exception as e:
                self.log_warning(f"Could not set up dipole autocorrelation tracker: {str(e)}")
                self.dipole_autocorr_tracker = None
        else:
            self.dipole_autocorr_tracker = None
            self.log_info("Dipole autocorrelation tracking disabled")
        
        # Store logger for later use
        self.logger_hoomd = logger
        
        # ===== MOMENTUM ZEROING SETUP =====
        if self.zero_momentum_enabled:
            try:
                self.log_info("Setting up momentum zeroing...")
                self.log_info(f"  Zero momentum period: {self.zero_momentum_period_ps:.3f} ps")
                
                # Calculate trigger period in timesteps for momentum zeroing
                # Use conservative step-based trigger for efficiency since exact timing is less critical
                if self.error_tolerance > 0:
                    # Adaptive timestep mode - use estimated steps based on expected timestep
                    estimated_dt_ps = 0.001  # Assume ~1 fs average timestep
                    zero_momentum_steps = max(1, int(self.zero_momentum_period_ps / estimated_dt_ps))
                    zero_momentum_steps = min(zero_momentum_steps, 10000)  # Cap at reasonable value
                    self.log_info(f"  Using estimated trigger: every {zero_momentum_steps} steps (adaptive mode)")
                else:
                    # Fixed timestep mode - use precise calculation
                    zero_momentum_steps = max(1, int(self.zero_momentum_period_ps / self.dt_ps))
                    self.log_info(f"  Using calculated trigger: every {zero_momentum_steps} steps (fixed mode)")
                
                # Create ZeroMomentum updater
                zero_momentum = hoomd.md.update.ZeroMomentum(
                    trigger=hoomd.trigger.Periodic(zero_momentum_steps)
                )
                
                # Add to simulation operations
                self.sim.operations.updaters.append(zero_momentum)
                
                self.log_info("Momentum zeroing successfully enabled:")
                self.log_info(f"  Trigger period: every {zero_momentum_steps} steps")
                self.log_info(f"  Target period: {self.zero_momentum_period_ps:.3f} ps")
                self.log_info("  Prevents center-of-mass drift during simulation")
                
            except Exception as e:
                self.log_warning(f"Could not set up momentum zeroing: {str(e)}")
                self.log_warning("Continuing simulation without momentum zeroing")
        else:
            self.log_info("Momentum zeroing disabled")
        
        # Create console output with only performance and time metrics
        console_items = ["timestep", "tps", "elapsed_time", "ns_per_day", "eta", "dt(fs)"]
        
        # Add adaptive timestep info (performance related)
        if hasattr(self, 'adaptive_action') and self.adaptive_action is not None:
            console_items.append("adaptive_error_tolerance")
        
        self.log_info("TRACKING AND LOGGING SETUP COMPLETED:")
        self.log_info("  All output systems now use precise time-based periods")
        
        # Log detailed summary of what's enabled
        enabled_features = []
        energy_output_period_ps = getattr(self, 'energy_output_period_ps', 0.1)
        fkt_output_period_ps = getattr(self, 'fkt_output_period_ps', 1.0)
        fkt_kmag = getattr(self, 'fkt_kmag', 1.0)
        dipole_output_period_ps = getattr(self, 'dipole_output_period_ps', 1.0)
        
        if self.enable_energy_tracking:
            enabled_features.append(f"detailed energy tracking ({energy_output_period_ps:.3f} ps)")
        if self.enable_fkt:
            enabled_features.append(f"F(k,t) density correlation ({fkt_output_period_ps:.3f} ps, k={fkt_kmag})")
        if self.enable_dipole_autocorr:
            enabled_features.append(f"dipole autocorrelation ({dipole_output_period_ps:.3f} ps)")
        if hasattr(self, 'adaptive_action') and self.adaptive_action is not None:
            enabled_features.append("adaptive timestep control")
        if self.zero_momentum_enabled:
            enabled_features.append(f"momentum zeroing ({self.zero_momentum_period_ps:.3f} ps)")
        
        if enabled_features:
            self.log_info(f"Advanced features enabled: {', '.join(enabled_features)}")
        else:
            self.log_info("Running with basic tracking only")
        
        console_output_period_ps = getattr(self, 'console_output_period_ps', 1.0)
        self.log_info(f"Console output: {console_output_period_ps:.3f} ps periods (time-based)")
        self.log_info("  Works accurately for both adaptive and fixed timestep modes")

    def setup_output_writers(self):
        """Configure GSD writer and console table for simulation output."""
        
        # For GSD output, we can still use step-based for efficiency in most cases
        # since GSD output is typically less frequent and timing precision is less critical
        if self.error_tolerance > 0:
            # Adaptive timestep mode - use smaller intervals for better timing
            gsd_trigger_steps = max(1, int(self.gsd_output_period_ps / 0.001))  # Assume ~1 fs effective timestep
            gsd_trigger_steps = min(gsd_trigger_steps, 10000)  # Cap at reasonable value
            gsd_trigger = hoomd.trigger.Periodic(gsd_trigger_steps)
            
            self.log_info("GSD output setup for adaptive timestep mode:")
            self.log_info(f"  GSD trigger: every {gsd_trigger_steps} steps (target: {self.gsd_output_period_ps:.3f} ps)")
            self.log_info("  Note: GSD timing may vary slightly due to adaptive timestep changes")
            
        else:
            # Fixed timestep mode - use calculated step-based triggers (most efficient)
            gsd_trigger = hoomd.trigger.Periodic(self.gsd_period)
            self.log_info("GSD output setup for fixed timestep mode:")
            self.log_info(f"  GSD trigger: every {self.gsd_period} steps ({self.gsd_output_period_ps:.3f} ps)")
        
        # Set up GSD writer
        gsd_writer = hoomd.write.GSD(
            filename=f'{self.name}-{self.replica}.gsd',
            trigger=gsd_trigger,
            dynamic=['property', 'momentum', 'particles/diameter', 'topology'],
            mode='wb',
            truncate=self.truncate_gsd,
            filter=hoomd.filter.All()
        )
        gsd_writer.logger = self.logger_hoomd
        
        # Write initial frame
        gsd_writer.write(self.sim.state, filename=f'{self.name}-{self.replica}.gsd',
                         mode='wb', filter=hoomd.filter.All(), logger=self.logger_hoomd)
        
        # Add GSD writer to simulation
        self.sim.operations.writers.append(gsd_writer)
        self.log_info(f"GSD writer added for file: {self.name}-{self.replica}.gsd")
        self.log_info(f"  GSD output period: {self.gsd_output_period_ps:.3f} ps")
        self.log_info(f"  GSD truncate mode: {self.truncate_gsd} ({'overwrite existing file' if self.truncate_gsd else 'append to existing file'})")
        
        # Create custom time-based console output tracker for accurate timing
        class ConsoleOutputTracker(hoomd.custom.Action):
            """Time-based console output tracker for accurate timing in both adaptive and fixed timestep modes."""
            
            def __init__(self, sim, time_tracker, performance_tracker, timestep_formatter, adaptive_action, output_period_ps, coupling_constant, incavity, cavity_force):
                super().__init__()
                self.sim = sim
                self.time_tracker = time_tracker
                self.performance_tracker = performance_tracker
                self.timestep_formatter = timestep_formatter
                self.adaptive_action = adaptive_action
                self.output_period_ps = output_period_ps
                self.coupling_constant = coupling_constant  # Target coupling (fallback)
                self.incavity = incavity
                self.cavity_force = cavity_force  # Reference to cavity force for current coupling
                self.last_output_time = 0.0
                self.header_printed = False
                
            def _get_current_coupling(self):
                """Get the current coupling value from the cavity force."""
                if not self.incavity or self.cavity_force is None:
                    return 0.0
                
                try:
                    # Check if the couplstr is a variant (time-varying) or constant
                    couplstr_param = self.cavity_force.couplstr
                    if hasattr(couplstr_param, 'current_value'):
                        # Time-varying coupling - get current value from variant
                        return couplstr_param.current_value
                    elif hasattr(couplstr_param, '__call__'):
                        # It's a variant but may not have current_value - call it with current timestep
                        return float(couplstr_param(self.sim.timestep))
                    else:
                        # Constant coupling
                        return float(couplstr_param)
                except Exception:
                    # Fallback to target coupling if we can't get current value
                    return self.coupling_constant
                
            def _get_current_time(self, timestep):
                """Get current simulation time."""
                if self.time_tracker is not None:
                    return self.time_tracker.elapsed_time
                else:
                    dt = float(self.sim.operations.integrator.dt)
                    return PhysicalConstants.atomic_units_to_ps(dt * timestep)
                
            def _should_output(self, timestep):
                """Check if we should output at this timestep using time-based logic."""
                current_time = self._get_current_time(timestep)
                time_since_last = current_time - self.last_output_time
                return time_since_last >= self.output_period_ps
                
            def act(self, timestep):
                if timestep == 0:
                    return
                    
                # Check if we should output using time-based logic
                if not self._should_output(timestep):
                    return
                    
                # Print header on first output
                if not self.header_printed:
                    header_parts = [
                        "Simulation.timestep",
                        "Simulation.tps", 
                        "Time.elapsed_ps",
                        "Performance.ns_per_day",
                        "Performance.eta",
                        "Timestep.dt_fs"
                    ]
                    if self.incavity:
                        header_parts.append("Cavity.coupling_au")
                    if self.adaptive_action is not None:
                        header_parts.append("Adaptive.error_tolerance")
                    print(" ".join(f"{h:>15s}" for h in header_parts), flush=True)
                    self.header_printed = True
                    
                # Collect current values
                current_time = self._get_current_time(timestep)
                tps = self.sim.tps if hasattr(self.sim, 'tps') else 0.0
                ns_per_day = self.performance_tracker.ns_per_day  # Property, not method
                eta = self.performance_tracker.eta_remaining      # Property, not method
                dt_fs = self.timestep_formatter.dt_fs             # Property, not method
                
                # Build output line
                output_parts = [
                    f"{timestep:15d}",
                    f"{tps:15.5f}",
                    f"{current_time:15.5f}",
                    f"{ns_per_day:>15s}",
                    f"{eta:>15s}",
                    f"{dt_fs:15.5f}"
                ]
                
                if self.incavity:
                    # Get current coupling value (changes with time for time-varying simulations)
                    current_coupling = self._get_current_coupling()
                    output_parts.append(f"{current_coupling:15.2e}")
                
                if self.adaptive_action is not None:
                    error_tol = self.adaptive_action.error_tolerance  # Property, not method
                    output_parts.append(f"{error_tol:15.2e}")
                    
                print(" ".join(output_parts), flush=True)
                
                # Update last output time
                self.last_output_time = current_time
        
        # Create and set up console output tracker
        console_tracker = ConsoleOutputTracker(
            sim=self.sim,
            time_tracker=self.time_tracker,
            performance_tracker=self.performance_tracker,
            timestep_formatter=self.timestep_formatter,
            adaptive_action=getattr(self, 'adaptive_action', None),
            output_period_ps=self.console_output_period_ps,
            coupling_constant=self.couplstr,
            incavity=self.incavity,
            cavity_force=self.cavity_force
        )
        
        # Add console tracker to simulation (check every step but handle timing internally)
        console_updater = hoomd.update.CustomUpdater(
            action=console_tracker,
            trigger=hoomd.trigger.Periodic(1)  # Check every step but tracker handles timing internally
        )
        self.sim.operations.updaters.append(console_updater)
        
        self.log_info("Console output setup completed:")
        self.log_info(f"  Output period: {self.console_output_period_ps:.3f} ps (accurate time-based)")
        if self.incavity:
            self.log_info(f"  Coupling constant column: {self.couplstr:.6e} a.u. (displayed in real-time)")
        self.log_info("  Console tracker handles timing internally using ElapsedTimeTracker")
        self.log_info("  Works accurately for both adaptive and fixed timestep modes")
        
        # Log final setup summary
        if self.error_tolerance > 0:
            self.log_info("FIXED: Console output now uses precise time-based logic")
            self.log_info("  Both console and energy tracker use accurate timing")
            self.log_info("  GSD output timing may vary slightly but is less critical")
        else:
            self.log_info("Time-based console output setup for fixed timestep mode")
            self.log_info("  Provides consistent behavior across all timestep modes")

    def setup_simulation(self):
        """Create HOOMD simulation object and initialize state from GSD file."""
        import os
        import gsd.hoomd
        
        # Save current directory and change to job directory
        self.original_cwd = os.getcwd()
        os.chdir(self.job_dir)
        
        # Setup device
        device = self.setup_device()
        
        # Create simulation object with seed control
        if self.seed is not None:
            simulation_seed = self.seed
            self.log_info(f"Using user-specified seed: {simulation_seed}")
        else:
            # Generate deterministic seed based on replica for reproducibility
            simulation_seed = hash(str(self.replica)) % (2**31)  # Use replica-based seed
            self.log_info(f"Using replica-based seed: {simulation_seed} (replica {self.replica})")
        
        self.sim = hoomd.Simulation(device=device, seed=simulation_seed)
        
        # Load GSD file and handle cavity particle
        with gsd.hoomd.open(self.input_gsd, 'r') as f:
            if self.frame < 0:
                self.frame = len(f) + self.frame  # Convert negative index to positive
                if self.frame < 0:  # Handle case where abs(frame) > len(f)
                    self.frame = 0
            snapshot = f[self.frame]
            
            if self.incavity:
                # Check if cavity particle already exists
                cavity_exists, cavity_count = self.check_cavity_particle_exists(snapshot)
                
                if cavity_exists:
                    # Cavity particle already exists - use original snapshot
                    self.log_info(f"Cavity particle already exists in GSD file ({cavity_count} found)")
                    if cavity_count > 1:
                        self.log_warning(f"Multiple cavity particles found ({cavity_count}). Expected exactly 1.")
                    self.sim.create_state_from_gsd(filename=self.input_gsd, frame=self.frame)
                    self.log_info(f"Simulation state created from original GSD file with existing cavity particle")
                    # Validate the existing cavity particle
                    self.validate_cavity_particle()
                    
                elif self.add_cavity_particle:
                    # No cavity particle exists, but we want to add one
                    self.log_info("No cavity particle found in GSD file - adding new cavity particle")
                    snapshot = self.create_cavity_particle(snapshot)
                    self.sim.create_state_from_snapshot(snapshot)
                    self.log_info(f"Simulation state created from modified snapshot with new cavity particle")
                    
                else:
                    # No cavity particle exists and we don't want to add one - this is an error
                    raise ValueError(
                        "ERROR: Cavity simulation requested but no cavity particle found in GSD file "
                        "and add_cavity_particle=False. Either:\n"
                        "  1. Set add_cavity_particle=True to automatically add a cavity particle, or\n"
                        "  2. Use a GSD file that already contains a cavity particle (type 'L' with typeid=2)"
                    )
            else:
                # No cavity simulation - use original GSD file
                self.sim.create_state_from_gsd(filename=self.input_gsd, frame=self.frame)
                self.log_info(f"Simulation state created from original GSD file (no cavity)")
        
        return snapshot

    def check_cavity_particle_exists(self, snapshot):
        """Check if a cavity particle already exists in the snapshot."""
        # Check if 'L' particle type exists
        if 'L' not in snapshot.particles.types:
            return False, 0
        
        # Check if any particles have typeid == 2 (cavity particles)
        cavity_count = np.sum(snapshot.particles.typeid == 2)
        return cavity_count > 0, cavity_count

    def create_cavity_particle(self, snapshot):
        """
        Add a new cavity particle to the simulation snapshot.
        
        Note: This method assumes no cavity particle already exists in the snapshot.
        Use check_cavity_particle_exists() to verify this before calling.
        """
        self.log_info("Creating new cavity particle and adding to system...")
        
        # Set numpy seed for reproducible cavity particle positioning if seed is provided
        if self.seed is not None:
            np.random.seed(self.seed + 2)  # Use seed+2 to differentiate from other random operations
            self.log_info(f"Using seed {self.seed + 2} for cavity particle positioning")
        
        # Calculate dipole moment and photon position
        positions = unwrap_positions(snapshot.particles.position, snapshot.particles.image, 
                                   snapshot.configuration.box[:3])
        dipmom = np.einsum('i,ij->j', snapshot.particles.charge, positions)

        omegac = self.freq / PhysicalConstants.HARTREE_TO_CM_MINUS1
        
        # Determine coupling strength for initial placement
        initial_couplstr = 0.0 if self.switch_time_ps is not None else self.couplstr

        if self.finite_q:
            # Allow finite-q photon displacement based on dipole-coupling interaction
            if self.switch_time_ps is not None:
                # For instantaneous switching: start at zero coupling position
                # Displacement will happen at switch time via CavityParticleDisplacer
                newpos = np.array([0.0, 0.0, 0.0])
                self.log_info(f"Finite-q + instantaneous switching: Photon starts at origin")
                self.log_info(f"  CavityParticleDisplacer will handle displacement at switch time")
            else:
                # Original finite-q behavior for constant coupling
                newpos = -dipmom * initial_couplstr / omegac**2
                newpos[-1] = 0.0
                # Only add thermal fluctuations if coupling is non-zero
                if initial_couplstr != 0.0:
                    sigma = np.sqrt(self.kB * self.temperature / omegac**2)
                    newpos = np.random.normal(loc=newpos, scale=sigma, size=3)
                    self.log_info(f"Finite-q mode: Photon displaced by dipole interaction to {newpos} (with thermal fluctuations)")
                else:
                    self.log_info(f"Finite-q mode: Photon at equilibrium position {newpos} (no thermal fluctuations due to zero coupling)")
        else:
            # Start photon at origin (q=0 limit)
            newpos = np.array([0.0, 0.0, 0.0])
            # Only add thermal fluctuations if coupling is non-zero
            if initial_couplstr != 0.0:
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

    def validate_cavity_particle(self):
        """Validate that cavity particle exists when required."""
        # Get particle types from simulation state, not local snapshot
        particle_types = self.sim.state.particle_types
        
        if 'L' not in particle_types:
            raise ValueError("ERROR: Cavity simulation requested but no cavity particle type 'L' found in GSD file.")
        
        with self.get_local_snapshot() as snap:
            # Convert to numpy array for robust handling across CPU/GPU
            typeid_array = self.to_numpy(snap.particles.typeid)
            
            if 2 not in typeid_array:
                raise ValueError("ERROR: Cavity simulation requested but no cavity particles found in GSD file.")
            
            cavity_count = np.sum(typeid_array == 2)
            if cavity_count != 1:
                raise ValueError(f"ERROR: Expected exactly 1 cavity particle but found {cavity_count} in GSD file.")
            
            cavity_indices = np.where(typeid_array == 2)[0]
            cavity_index = cavity_indices[0]
            cavity_position = snap.particles.position[cavity_index]
            self.log_info(f"Cavity particle validated at index {cavity_index}, position {cavity_position}")

    def run(self):
        """Main orchestrator method that runs the complete simulation workflow."""
        try:
            # Phase 0: Setup logging
            self.setup_logging()
            
            # Phase 1: Setup simulation and state
            self.log_info("=== Phase 1: Setting up simulation ===")
            self.calculate_physical_parameters()
            snapshot = self.setup_simulation()

            # Phase 1.5: Setup time tracker (must be after sim setup, before force setup)
            self.log_info("=== Phase 1.5: Setting up time tracker ===")
            self.time_tracker = ElapsedTimeTracker(self.sim, self.runtime_ps)
            self.sim.operations.updaters.append(hoomd.update.CustomUpdater(
                action=self.time_tracker, trigger=hoomd.trigger.Periodic(1)
            ))
            
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

    def run_simulation(self):
        """Execute the main simulation loop."""
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
        """Handle post-simulation cleanup and restore original directory."""
        # Note: Trackers from analysis.py write directly to files, no buffering needed
        self.log_info("Cleanup initiated...")
        
        # Restore original directory
        if hasattr(self, 'original_cwd'):
            os.chdir(self.original_cwd)
            self.log_info(f"Returned to original directory: {self.original_cwd}")
        
        self.log_info("Cleanup completed")

    def get_local_snapshot(self):
        """
        Get the appropriate local snapshot context manager based on device type.
        
        Returns
        -------
        Local snapshot context manager (CPU or GPU)
        """
        if self.device == 'GPU':
            return self.sim.state.gpu_local_snapshot
        else:
            return self.sim.state.cpu_local_snapshot

    def get_particle_count(self, snap):
        """
        Get the number of particles from a local snapshot in a device-agnostic way.
        
        Parameters
        ----------
        snap : Local snapshot object
            Either CPU or GPU local snapshot
            
        Returns
        -------
        int
            Number of particles
        """
        if self.device == 'GPU':
            # GPU snapshots don't have .N attribute, use length of an array
            return len(snap.particles.typeid)
        else:
            # CPU snapshots have .N attribute
            return snap.particles.N

    def convert_to_numpy(self, array_data):
        """
        Convert array data to NumPy array, handling both CPU (NumPy) and GPU (CuPy) arrays.
        
        DEPRECATED: Use to_numpy() instead for better CuPy compatibility.
        
        Parameters
        ----------
        array_data : array-like
            Array data from local snapshot (could be NumPy or CuPy)
            
        Returns
        -------
        np.ndarray
            NumPy array
        """
        return self.to_numpy(array_data)
    def _force_thermalization_for_bussi(self):
        """
        Force thermalization when zero velocities are detected with Bussi thermostat.
        
        This is a safety measure to prevent Bussi thermostat failures when the GSD
        file contains particles at rest.
        """
        self.log_info("Performing forced thermalization for Bussi thermostat compatibility...")
        
        kT = self.kB * self.temperature
        
        # Set numpy seed for reproducible thermalization if seed is provided
        if self.seed is not None:
            np.random.seed(self.seed + 1)  # Use seed+1 to differentiate from HOOMD seed
            self.log_info(f"Using seed {self.seed + 1} for forced thermalization")
        
        molecular_filter = hoomd.filter.Type(['O', 'N'])  # Molecular particles only
        
        if self.incavity:
            # Use HOOMD's State.thermalize_particle_momenta for molecular system
            self.sim.state.thermalize_particle_momenta(kT=kT, filter=molecular_filter)
            self.log_info("Force-thermalized molecular particles using HOOMD State object")
            
            # Manual thermalization for cavity particle
            self._thermalize_cavity_particle_manually(kT)
            
        else:
            # No cavity particle, thermalize all particles
            self.sim.state.thermalize_particle_momenta(kT=kT, filter=hoomd.filter.All())
            self.log_info("Force-thermalized all molecular particles")
            
        # Verification after forced thermalization
        with self.get_local_snapshot() as snap:
            masses = self.to_numpy(snap.particles.mass)
            velocities = self.to_numpy(snap.particles.velocity)
            typeid_array = self.to_numpy(snap.particles.typeid)
            
            # Filter for molecular particles only (typeid != 2)
            molecular_mask = typeid_array != 2
            
            if np.any(molecular_mask):
                mol_masses = masses[molecular_mask]
                mol_velocities = velocities[molecular_mask]
                
                # Calculate kinetic energy and temperature
                v_squared = np.sum(mol_velocities**2, axis=1)
                kinetic_energy = 0.5 * mol_masses * v_squared
                total_ke = np.sum(kinetic_energy)
                
                n_mol_particles = np.sum(molecular_mask)
                degrees_of_freedom = 3 * n_mol_particles
                current_temp = (2.0) * total_ke / (degrees_of_freedom * self.kB)
                
                self.log_info(f"Post-thermalization kinetic temperature: {current_temp:.2f} K")
                
                # Check maximum velocity magnitude to confirm non-zero velocities
                max_velocity_magnitude = np.max(np.sqrt(np.sum(mol_velocities**2, axis=1)))
                self.log_info(f"Maximum particle velocity magnitude: {max_velocity_magnitude:.6e} a.u.")
                
                if max_velocity_magnitude > 1e-12:
                    self.log_info("Forced thermalization successful - Bussi thermostat should work now")
                else:
                    self.log_error("ERROR: Forced thermalization failed - velocities are still zero!")
                    raise RuntimeError("Unable to initialize non-zero velocities for Bussi thermostat")



class AdaptiveTimestepUpdater(hoomd.custom.Action):
    """Update timestep adaptively based on energy conservation and stability."""
    
    def __init__(self, state, integrator, error_tolerance, time_constant_ps=50.0, 
                 initial_fraction=0.01, adaptiveerror=True, cavity_damping_factor=1.0, 
                 molecular_thermostat_tau=5.0, cavity_thermostat_tau=5.0, time_tracker=None,
                 switch_time_ps=None, timestep_change_threshold=0.1, max_timestep_change_factor=1.5,
                 shock_dampening_factor=1e-3, shock_dampening_enabled=True, shock_dampening_time_constant_ps=50.0):
        super().__init__()
        

        
        self.state = state
        self.integrator = integrator
        self.target_error_tolerance = error_tolerance
        self.initial_error_tolerance = error_tolerance * initial_fraction
        self.time_constant_ps = time_constant_ps
        self.accumulated_time_ps = 0.0  # Fallback for backward compatibility
        self.last_timestep = 0  # Will be set correctly in first act() call
        self.adaptiveerror = adaptiveerror
        self.cavity_damping_factor = cavity_damping_factor
        self.molecular_thermostat_tau = molecular_thermostat_tau
        self.cavity_thermostat_tau = cavity_thermostat_tau
        self.time_tracker = time_tracker  # Reference to ElapsedTimeTracker for accurate timing
        self.switch_time_ps = switch_time_ps #Time when ramping should start
        
        # NEW: Shock dampening parameters
        self.shock_dampening_factor = shock_dampening_factor
        self.shock_dampening_enabled = shock_dampening_enabled
        
        # Calculate shock dampening error tolerance
        if shock_dampening_enabled and switch_time_ps is not None:
            self.shock_error_tolerance = error_tolerance * shock_dampening_factor
        else:
            self.shock_error_tolerance = self.initial_error_tolerance
        
        # Initialize current_error_tolerance based on switch time behavior
        if switch_time_ps is not None:
            # Start with target tolerance before switch
            self.current_error_tolerance = self.target_error_tolerance
        else:
            # Original behavior: start with low tolerance for immediate ramping
            self.current_error_tolerance = self.initial_error_tolerance
        
        # Conservative timestep update parameters
        self.timestep_change_threshold = timestep_change_threshold  # Minimum relative change to trigger update
        self.max_timestep_change_factor = max_timestep_change_factor  # Maximum factor by which timestep can change
        self.last_timestep_update = 0  # Track when we last updated the timestep
        self.min_update_interval = 1000  # Minimum steps between timestep updates
        
        # NEW: Switch detection variables
        self.last_elapsed_time_ps = 0.0  # Track previous elapsed time to detect switch crossing
        self.switch_detected = False  # Flag to track if we've detected the switch
        self.switch_detection_tolerance = 0.001  # ps - tolerance for detecting switch crossing
        
        # Log the shock dampening behavior
        if shock_dampening_enabled and switch_time_ps is not None:
            print(f"Shock dampening enabled with time constant {time_constant_ps} ps", flush=True)
            print(f"Before switch: error_tolerance = {self.target_error_tolerance:.2e} (normal tolerance)", flush=True)
            print(f"At switch: error_tolerance drops to {self.shock_error_tolerance:.2e} (shock dampening factor: {shock_dampening_factor:.2e})", flush=True)
            print(f"After switch: error_tolerance ramps from {self.shock_error_tolerance:.2e} to {self.target_error_tolerance:.2e} with τ = {time_constant_ps} ps", flush=True)
            print(f"Switch detection: immediate timestep adjustment within {self.switch_detection_tolerance} ps of switch", flush=True)
        elif switch_time_ps is not None:
            print(f"Error tolerance ramping will start at t = {self.switch_time_ps} ps", flush=True)
            print(f"Before switch: error_tolerance = {self.target_error_tolerance:.2e} (final tolerance for efficiency)", flush=True)
            if self.shock_dampening_factor > 0:
                print(f"At switch: error_tolerance drops to {self.initial_error_tolerance:.2e} (high precision)", flush=True)
                print(f"After switch: error_tolerance ramps from {self.initial_error_tolerance:.2e} to {self.target_error_tolerance:.2e} with τ = {time_constant_ps} ps", flush=True)
            else:
                print(f"At switch: error_tolerance remains at {self.target_error_tolerance:.2e} (no ramping, initial_fraction = 0.0)", flush=True)
            print(f"Switch detection: immediate timestep adjustment within {self.switch_detection_tolerance} ps of switch", flush=True)
        else:
            print("Error tolerance ramping starts immediately from t = 0 ps", flush=True)
        
        # Log conservative timestep parameters
        print(f"Conservative timestep parameters:", flush=True)
        print(f"  Change threshold: {self.timestep_change_threshold:.1%}", flush=True)
        print(f"  Max change factor: {self.max_timestep_change_factor:.1f}", flush=True)
        print(f"  Min update interval: {self.min_update_interval} steps", flush=True)
        print(f"  Switch detection tolerance: {self.switch_detection_tolerance} ps", flush=True)

        # DEBUG: Show key variables for shock dampening
        #print(f"DEBUG: shock_dampening_factor = {shock_dampening_factor:.2e}", flush=True)
        #print(f"DEBUG: shock_dampening_enabled = {shock_dampening_enabled}", flush=True)
        #print(f"DEBUG: initial_error_tolerance = {self.initial_error_tolerance:.2e}", flush=True)
        #print(f"DEBUG: shock_error_tolerance = {self.shock_error_tolerance:.2e}", flush=True)

    def act(self, timestep):
        """
        Custom action to update the timestep size.

        Parameters:
        - timestep: Current simulation timestep.
        """
        # Initialize last_timestep on first call to handle resuming from checkpoints
        if self.last_timestep == 0:
            self.last_timestep = timestep
        
        # Update accumulated simulation time (fallback method for backward compatibility)
        if timestep > self.last_timestep:
            # Convert dt to picoseconds using correct conversion
            dt_ps = PhysicalConstants.atomic_units_to_ps(self.integrator.dt)
            self.accumulated_time_ps += (timestep - self.last_timestep) * dt_ps
        self.last_timestep = timestep
        
        # Get accurate elapsed time for error tolerance ramping
        if self.time_tracker is not None:
            current_elapsed_time_ps = self.time_tracker.elapsed_time
        else:
            current_elapsed_time_ps = self.accumulated_time_ps
        
        # NEW: Detect if we've just crossed the switch time
        force_immediate_update = False
        #print(f"DEBUG: self.switch_time_ps = {self.switch_time_ps:.6f} ps", flush=True)
        #print(f"DEBUG: current_elapsed_time_ps = {current_elapsed_time_ps:.6f} ps", flush=True)
        #print(f"DEBUG: self.current_error_tolerance = {self.current_error_tolerance:.6f}", flush=True)
        #print(f"DEBUG: self.shock_dampening_factor = {self.shock_dampening_factor:.6f}", flush=True)
        if (self.switch_time_ps is not None and not self.switch_detected and
            self.last_elapsed_time_ps < self.switch_time_ps <= current_elapsed_time_ps):
            
            # We've just crossed the switch time!
            self.switch_detected = True
            force_immediate_update = True
            print(f"SWITCH DETECTED at timestep {timestep}: t = {current_elapsed_time_ps:.6f} ps", flush=True)
            print(f"  Forcing immediate timestep adjustment for shock dampening", flush=True)
        
        # Update last elapsed time for next iteration
        self.last_elapsed_time_ps = current_elapsed_time_ps
        
        # Update error tolerance based on shock dampening and exponential ramping
        if self.adaptiveerror:
            if self.switch_time_ps is not None:
                if current_elapsed_time_ps < self.switch_time_ps:
                    # Before switch: use target error tolerance
                    self.current_error_tolerance = self.target_error_tolerance
                else:
                    if current_elapsed_time_ps < self.switch_time_ps:
                        self.current_error_tolerance = self.shock_error_tolerance
                        #_time_ps = {current_elapsed_time_ps} ps", flush=True)
                        #print(f"DEBUG: self.current_error_tolerance = {self.current_error_tolerance}", flush=True)
                    else:
                    # After switch: implement shock dampening and exponential recovery
                        ramping_time = current_elapsed_time_ps - self.switch_time_ps
                        
                        #if self.shock_dampening_enabled:
                            # NEW: Shock dampening mode - exponential recovery from shock tolerance
                            # SPECIAL CASE: If switch was just detected, start with exact shock tolerance
                        if self.shock_dampening_enabled:
                            if force_immediate_update:
                                # At the exact moment of switch detection, use shock tolerance
                                self.current_error_tolerance = self.shock_error_tolerance
                            else:
                                # Normal exponential recovery after switch
                                exp_factor = np.exp(-ramping_time / self.time_constant_ps)
                                self.current_error_tolerance = self.target_error_tolerance - \
                                                                (self.target_error_tolerance - self.shock_error_tolerance) * exp_factor
                        #else:
                        #    # No ramping when initial_fraction = 0.0 - just use target tolerance
                        #    self.current_error_tolerance = self.target_error_tolerance
                    #if timestep % 1000 == 0:  # Debug output
                    #    print(f"DEBUG: shock_dampening mode, ramping_time={ramping_time:.6f}, current_error_tolerance={self.current_error_tolerance:.2e}", flush=True)
                    #elif self.shock_dampening_factor > 0:
                        # Original mode - exponential recovery from initial tolerance (only if initial_fraction > 0)
                    #    exp_factor = np.exp(-ramping_time / self.time_constant_ps)
                    #    self.current_error_tolerance = self.target_error_tolerance - \
                    #                                  (self.target_error_tolerance - self.initial_error_tolerance) * exp_factor
                        #if timestep % 1000 == 0:  # Debug output
                        #    print(f"DEBUG: ramping mode, ramping_time={ramping_time:.6f}, exp_factor={exp_factor:.6f}, current_error_tolerance={self.current_error_tolerance:.2e}", flush=True)
                        #if timestep % 1000 == 0:  # Debug output
                        #    print(f"DEBUG: no ramping mode, current_error_tolerance={self.current_error_tolerance:.2e}", flush=True)
            else:
                # Original behavior: start ramping immediately from t=0 (only if initial_fraction > 0)
                if self.shock_dampening_factor > 0:
                    exp_factor = np.exp(-current_elapsed_time_ps / self.time_constant_ps)
                    self.current_error_tolerance = self.target_error_tolerance - \
                                                  (self.target_error_tolerance - self.initial_error_tolerance) * exp_factor
                else:
                    # No ramping when initial_fraction = 0.0 - just use target tolerance
                    self.current_error_tolerance = self.target_error_tolerance
        else:
            self.current_error_tolerance = self.target_error_tolerance
        
        # Check if we should update the timestep (normal interval OR forced immediate update)
        steps_since_last_update = timestep - self.last_timestep_update
        should_update_timestep = (force_immediate_update or 
                                 steps_since_last_update >= self.min_update_interval)
        
        if not should_update_timestep:
            # Too soon to update timestep again (and no forced update)
            return
        
        # Collect forces and masses
        particle_data = self.state.get_snapshot().particles
        masses = np.array(particle_data.mass)
        n_particles = len(masses)
        
        # Initialize total forces array
        total_forces = np.zeros((n_particles, 3))
        
        # Sum forces from all force objects - let errors propagate clearly
        for force in self.integrator.forces:
            particle_forces = force.forces
            if particle_forces is not None and len(particle_forces) == n_particles:
                total_forces += np.asarray(particle_forces)
        
        # Calculate sum |f_i| / m_i
        force_norm = np.array([np.linalg.norm(f) for f in total_forces])
        force_max_norm = np.max(force_norm / masses) * n_particles
        
        # DEBUG: Show force information when timestep becomes very large
        # if timestep % 1000 == 0 or force_immediate_update:
        #     max_force = np.max(force_norm) if len(force_norm) > 0 else 0.0
        #     mean_force = np.mean(force_norm) if len(force_norm) > 0 else 0.0
        #     min_mass = np.min(masses) if len(masses) > 0 else 0.0
        #     max_mass = np.max(masses) if len(masses) > 0 else 0.0
            #print(f"DEBUG FORCES: max_force={max_force:.2e} a.u., mean_force={mean_force:.2e} a.u., force_mass_sum={force_mass_sum:.2e} a.u.^-1", flush=True)
            #print(f"DEBUG MASSES: min_mass={min_mass:.2e} a.u., max_mass={max_mass:.2e} a.u., n_particles={len(masses)}", flush=True)
            
        # Compute optimal timestep using current error tolerance
        if force_max_norm > 0:
            # Use error tolerance directly (no minimum protection)
            # This allows observing dt=0.0 when error tolerance becomes zero during shock dampening
            effective_error_tolerance = self.current_error_tolerance
            #print(f"DEBUG: effective_error_tolerance = {effective_error_tolerance}", flush=True)

            optimal_dt = np.sqrt(effective_error_tolerance / force_max_norm)
            current_dt = self.integrator.dt
            # print(f"DEBUG: current_dt = {current_dt}", flush=True)
            # print(f"DEBUG: force_max_norm = {force_max_norm}", flush=True)
            # print(f"DEBUG: optimal_dt = {optimal_dt}", flush=True)
            # DEBUG: Show exact values when switch is detected or shock dampening is active
            #if force_immediate_update or (self.shock_dampening_enabled and self.switch_detected):
            #    print(f"DEBUG SHOCK DAMPENING: effective_error_tolerance = {effective_error_tolerance:.2e}, optimal_dt = {optimal_dt:.6e} a.u. ({PhysicalConstants.atomic_units_to_ps(optimal_dt):.1f} ps)", flush=True)
            
            # Apply conservative timestep update logic
            # dt_ratio = optimal_dt / current_dt
            
            # # For forced updates (switch detected), be more aggressive with changes
            # if force_immediate_update:
            #     # Allow larger changes when switch is detected for immediate stabilization
            #     effective_change_threshold = self.timestep_change_threshold * 0.1  # 10x more sensitive
            #     effective_max_change_factor = self.max_timestep_change_factor * 2.0  # Allow 2x larger changes
            # else:
            #     # Normal conservative parameters
            #     effective_change_threshold = self.timestep_change_threshold
            #     effective_max_change_factor = self.max_timestep_change_factor
            
            # Only update if the change is significant (or forced)
            #if force_immediate_update or abs(dt_ratio - 1.0) > effective_change_threshold:
            # SPECIAL CASE: Allow dt=0 when optimal_dt is exactly zero (shock dampening)
            # if optimal_dt == 0.0:
            #     new_dt = 0.0
            #     clamped = False
            #     print(f"  SPECIAL: Setting dt = 0.0 exactly (optimal_dt = 0.0 from shock dampening)", flush=True)
            # # Limit the maximum change to prevent dramatic jumps (normal cases)
            # elif dt_ratio > effective_max_change_factor:
            #     new_dt = current_dt * effective_max_change_factor
            #     clamped = True
            # elif dt_ratio < (1.0 / effective_max_change_factor):
            #     new_dt = current_dt / effective_max_change_factor
            #     clamped = True
            # else:
            #     new_dt = optimal_dt
            #     clamped = False
            
            # Apply the timestep change
            self.integrator.dt = optimal_dt
            self.last_timestep_update = timestep
                
        else:
            if timestep % 5000 == 0:
                print(f"WARNING: Zero force detected at step {timestep} - keeping current timestep", flush=True)

    @property
    def error_tolerance(self):
        """Access the current effective error tolerance."""
        return self.current_error_tolerance

    @hoomd.logging.log
    def elapsed_time_ps(self):
        """Return the elapsed time in picoseconds."""
        return self.accumulated_time_ps


