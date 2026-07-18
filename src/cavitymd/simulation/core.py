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

from ..utils import PhysicalConstants, unwrap_positions
from ..forces import CavityForce
from ..analysis import (
    Status, ElapsedTimeTracker, TimestepFormatter, FieldAutocorrelationTracker,
    EnergyTracker, PerformanceTracker, AutocorrelationTracker,
    TemperatureTracker
)
from ..controllers.empirical import EmpiricalTemperatureData
from .timestep import AdaptiveTimestepUpdater
from ..updaters import CavityParticleDisplacer
from ..variants import StepVariant, PeriodicVariant, ExponentialDecayVariant, SquareWaveVariant, DecayingSquareWaveVariant, AdaptiveSquareWaveVariant
from ..composite_variant import CompositeVariant
from ..controllers import DiffEqController, SimpleSetpointController
from ..data import ObservableWriter


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
    incavity : bool
        Whether to enable cavity coupling
    lambda_coupling : float
        Dimensionless coupling parameter :math:`\lambda`. The effective coupling strength
        is computed as :math:`\varepsilon = \lambda \cdot \omega_c` where :math:`\omega_c`
        is the cavity characteristic frequency. This parameter is independent of the
        cavity frequency, allowing the same :math:`\lambda` values to be used across
        different cavity frequencies.
    couplstr : float, optional
        **DEPRECATED**: Use `lambda_coupling` instead. For backward compatibility, if
        provided, it will be used as the effective coupling strength directly (not scaled
        by omega_c). This parameter will be removed in a future version.
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
    enable_empirical_feedback : bool, optional
        Whether to enable empirical temperature feedback control. Default: False
    empirical_data_file : str, optional
        Path to empirical energy-temperature data file. Required if enable_empirical_feedback=True
    feedback_output_period_ps : float, optional
        Output period for empirical feedback CSV data in ps. Default: 0.1
    feedback_energy_component : str, optional
        Energy component for feedback ('total_PE', 'lj_coulombic', etc.). Default: 'lj_coulombic'
    feedback_apply_to : str, optional
        Which thermostats to apply feedback to ('molecular', 'cavity', 'both'). Default: 'both'
    feedback_averaging_window_ps : float, optional
        Time window for averaging fictive temperature in ps. Default: 5.0
    feedback_update_interval_ps : float, optional
        Interval between temperature updates in ps. Default: 5.0
    feedback_T_min : float, optional
        Minimum allowed temperature in K. Default: 50.0
    feedback_T_max : float, optional
        Maximum allowed temperature in K. Default: 300.0
    feedback_turn_off_time_ps : float, optional
        Time in ps to turn off empirical feedback. If None, feedback runs until end of simulation. Default: None
    feedback_auto_adjust_molecular_tau : bool, optional
        Whether to automatically set molecular thermostat tau to 2×Δt when feedback turns off. Default: False
        
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
    fkt_log_time_spacing : bool, optional
        Whether to use logarithmic time spacing for F(k,t) output. Default: False
    fkt_min_log_time_ps : float, optional
        Minimum time for logarithmic spacing (required if fkt_log_time_spacing=True). Default: None
    fkt_max_log_time_ps : float, optional
        Maximum time for logarithmic spacing (required if fkt_log_time_spacing=True). Default: None
    fkt_log_num_points : int, optional
        Number of points in logarithmic time spacing. Default: 50
    energy_output_period_ps : float, optional
        Energy tracker output period in ps. Default: 0.1
    fkt_output_period_ps : float, optional
        F(k,t) output period in ps. Default: 1.0
    gsd_output_period_ps : float, optional
        Trajectory output period in ps. Default: 50.0
    enable_gsd_output : bool, optional
        Write production GSD trajectories. When False, no ``prod-*.gsd``
        writer is created (HDF5 and F(k,t) outputs are unaffected). Default: True
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
    
    **Empirical temperature feedback experiment:**
    
    >>> feedback_sim = CavityMDSimulation(
    ...     job_dir="empirical_feedback",
    ...     replica=1,
    ...     freq=2000.0,
    ...     couplstr=0.001,
    ...     incavity=True,
    ...     runtime_ps=1000.0,
    ...     enable_empirical_feedback=True,
    ...     empirical_data_file="equilibrium_data.txt",
    ...     feedback_energy_component="lj_coulombic",
    ...     feedback_averaging_window_ps=5.0,
    ...     feedback_update_interval_ps=5.0,
    ...     enable_energy_tracking=True  # Required for empirical feedback
    ... )
    >>> exit_code = feedback_sim.run()
    
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
    hoomd.cavitymd.coupling.StepVariant : For time-varying parameters
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
                 incavity: bool,
                 lambda_coupling: Optional[float] = None,
                 couplstr: Optional[float] = None,
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
                 # F(k,t) logarithmic time spacing parameters
                 fkt_log_time_spacing: bool = False,
                 fkt_min_log_time_ps: Optional[float] = None,
                 fkt_max_log_time_ps: Optional[float] = None,
                 fkt_log_num_points: int = 50,
                 enable_dipole_autocorr: bool = False,
                 dipole_reference_interval_ps: float = 1.0,
                 dipole_max_references: int = 10,
                 max_energy_output_time_ps: Optional[float] = None, 
                 enable_energy_tracking: bool = True, 
                 dt_fs: Optional[float] = None, 
                 device: str = 'CPU', 
                 gpu_id: int = 0,
                 energy_output_period_ps: float = 0.1, 
                 fkt_output_period_ps: float = 1.0, 
                 dipole_output_period_ps: float = 1.0,
                 gsd_output_period_ps: float = 50.0,
                 enable_gsd_output: bool = True,
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
                 # Dynamic coupling detection parameters
                 enable_dynamic_coupling_detection: bool = True,
                 coupling_change_threshold: float = 1e-5,
                 transition_time_ps: float = 0.5,
                 use_smooth_switching: bool = True,
                 zero_momentum_enabled: bool = False,
                 zero_momentum_period_ps: float = 1.0,
                 # Empirical temperature feedback parameters
                 enable_empirical_feedback: bool = False,
                 empirical_data_file: Optional[str] = None,
                 feedback_output_period_ps: float = 0.1,
                 feedback_energy_component: str = 'lj_coulombic',
                 feedback_apply_to: str = 'both',
                 feedback_averaging_window_ps: float = 5.0,
                 feedback_update_interval_ps: float = 5.0,
                 feedback_T_min: float = 50.0,
                 feedback_T_max: float = 300.0,
                 feedback_turn_off_time_ps: Optional[float] = None,
                 feedback_use_direct_harmonic_calculation: bool = True,
                 feedback_auto_adjust_molecular_tau: bool = False,
                 # Periodic coupling parameters
                 periodic_coupling: bool = False,
                 periodic_period_ps: float = 1.0,
                 periodic_phase_offset: float = 0.0,
                 periodic_start_time_ps: float = 0.0,
                 periodic_stop_time_ps: Optional[float] = None,
                 # Enhanced coupling variant parameters
                 coupling_variant_type: str = 'constant',  # 'constant', 'step', 'periodic', 'exponential', 'square', 'exponentialwave'
                 # Exponential decay parameters
                 exponential_amplitude: Optional[float] = None,
                 exponential_decay_time_ps: float = 10.0,
                 exponential_turn_on_time_ps: float = 0.0,
                 exponential_turn_off_time_ps: Optional[float] = None,
                 # Square wave parameters  
                 square_period_ps: float = 2.0,
                 square_duty_cycle: float = 0.5,
                 square_phase_offset: float = 0.0,
                 square_start_time_ps: float = 0.0,
                 square_stop_time_ps: Optional[float] = None,
                 # Decaying square wave parameters
                 decaying_square_period_ps: float = 2.0,
                 decaying_square_duty_cycle: float = 0.5,
                 decaying_square_phase_offset: float = 0.0,
                 decaying_square_start_time_ps: float = 0.0,
                 decaying_square_stop_time_ps: Optional[float] = None,
                 decaying_square_decay_rate: float = 0.1,
                 decaying_square_minimum_amplitude: float = 1e-6,
                 # Adaptive square wave parameters
                 adaptive_square_period_ps: float = 5.0,
                 adaptive_square_duty_cycle: float = 0.5,
                 adaptive_square_phase_offset: float = 0.0,
                 adaptive_square_start_time_ps: float = 0.0,
                 adaptive_square_stop_time_ps: Optional[float] = None,
                 adaptive_square_min_amplitude: float = 1e-8,
                 adaptive_square_max_amplitude: float = 1e-1,
                 # Exponential wave parameters
                 exp_period_ps: float = 2.0,
                 exp_tau_ps: float = 0.5,
                 exp_start_time_ps: float = 0.0,
                 exp_stop_time_ps: Optional[float] = None,
                 exp_adaptive: bool = False,
                 # Composite coupling parameters
                 composite_sinusoid_amplitude: float = 1e-4,
                 composite_sinusoid_period: float = 1.0,
                 composite_sinusoid_phase: float = 0.0,
                 composite_sinusoid_start_time: float = 0.0,
                 composite_sinusoid_stop_time: Optional[float] = None,
                 composite_square_amplitude: float = 2e-4,
                 composite_square_period: float = 50.0,
                 composite_square_duty_cycle: float = 0.02,
                 composite_square_start_time: float = 10.0,
                 composite_square_stop_time: float = 1000.0,
                 composite_square_adaptive: bool = False,
                 composite_max_amplitude: Optional[float] = None,
                 # Auto-stop coupling parameters
                 enable_auto_stop: bool = False,
                 auto_stop_tol: float = 1.0,
                 auto_stop_window: float = 10.0,
                 # Enhanced laser drive parameters  
                 laser_enabled: bool = False,
                 laser_frequency_cm1: float = 0.0,
                 laser_amplitude: float = 0.0,
                 laser_start_time_ps: float = 0.0,
                 laser_stop_time_ps: Optional[float] = None,
                 laser_kvector: Optional[List[float]] = None,
                 # PI feedback controller parameters
                 enable_pi_feedback: bool = False,
                 pi_target_temperature: float = 100.0,
                 pi_turn_on_time_ps: float = 0.0,
                 pi_turn_off_time_ps: Optional[float] = None,
                 pi_temperature_method: str = 'kinetic',  # 'kinetic', 'lj_coulombic', 'harmonic'
                 pi_update_interval_ps: float = 5.0,
                 pi_Kc: Optional[float] = None,  # Auto-calculate if None
                 pi_Ti: Optional[float] = None,  # Auto-calculate if None
                 pi_molecular_tau_ps: float = 5.0,  # For IMC auto-tuning
                 pi_beta: float = 0.7,  # Setpoint weighting
                 pi_T_min: float = 0.0,
                 pi_T_max: Optional[float] = None,
                 # Gradient descent feedback controller parameters
                 enable_gd_feedback: bool = False,
                 gd_target_temperature: float = 100.0,
                 gd_turn_on_time_ps: float = 0.0,
                 gd_turn_off_time_ps: Optional[float] = None,
                 gd_temperature_method: str = 'kinetic',  # 'kinetic', 'lj_coulombic', 'harmonic'
                 gd_update_interval_ps: float = 0.1,
                 gd_time_constant_ps: float = 10.0,
                 gd_apply_to: str = 'both',
                 gd_T_min: float = 0.0,
                 gd_T_max: Optional[float] = None,
                 gd_disable_effective_temp: bool = False,
                 # Multi-signal error function parameters
                 gd_enable_multi_signal: bool = False,
                 gd_weight_system_target: float = 1.0,
                 gd_weight_bath_target: float = 0.0,
                 gd_weight_system_bath: float = 0.0,
                 # Dual independent feedback controller parameters
                 enable_dual_feedback: bool = False,
                 dual_cavity_method: str = 'harmonic_equipartition',
                 dual_molecular_method: str = 'lj_coulombic_kinetic',
                 dual_cavity_target_temperature: float = 100.0,
                 dual_molecular_target_temperature: float = 100.0,
                 dual_cavity_time_constant_ps: float = 5.0,
                 dual_molecular_time_constant_ps: float = 10.0,
                 dual_turn_on_time_ps: float = 0.0,
                 dual_turn_off_time_ps: Optional[float] = None,
                 dual_update_interval_ps: float = 0.1,
                 dual_cavity_T_min: float = 0.0,
                 dual_cavity_T_max: Optional[float] = None,
                 dual_molecular_T_min: float = 0.0,
                 dual_molecular_T_max: Optional[float] = None,
                 dual_cavity_dynamic_target: bool = False,
                 dual_molecular_dynamic_target: bool = False,
                 dual_cavity_integral_time_constant_ps: Optional[float] = None,
                 dual_molecular_integral_time_constant_ps: Optional[float] = None,
                 # Sinusoidal bath temperature controller parameters
                 enable_sinusoidal_bath: bool = False,
                 sinusoidal_bath_period_ps: float = 1.0,
                 sinusoidal_bath_amplitude_scale: float = 0.1,
                 sinusoidal_bath_phase_offset: float = 0.0,
                 sinusoidal_bath_target_temperature: float = 100.0,
                 sinusoidal_bath_dynamic_target: bool = False,
                 sinusoidal_bath_turn_on_time_ps: float = 0.0,
                 sinusoidal_bath_turn_off_time_ps: Optional[float] = None,
                 sinusoidal_bath_update_interval_ps: float = 0.1,
                 sinusoidal_bath_apply_to: str = 'both',
                 sinusoidal_bath_T_min: float = 0.1,
                 sinusoidal_bath_T_max: Optional[float] = None,
                 sinusoidal_bath_empirical_data_file: Optional[str] = None,
                 sinusoidal_bath_amplitude_update_interval_ps: float = 1.0,
                 sinusoidal_bath_amplitude_temperature_method: str = 'harmonic_equipartition',
                 sinusoidal_bath_adaptive_range_mode: bool = False,
                 # Adaptive bath temperature controller parameters
                 enable_adaptive_bath: bool = False,
                 adaptive_bath_amplitude_scale: float = 1.0,
                 adaptive_bath_time_constant_ps: float = 1.0,
                 adaptive_bath_target_temperature: float = 100.0,
                 adaptive_bath_dynamic_target: bool = True,
                 adaptive_bath_turn_on_time_ps: float = 0.0,
                 adaptive_bath_turn_off_time_ps: Optional[float] = None,
                 adaptive_bath_update_interval_ps: float = 0.1,
                 adaptive_bath_apply_to: str = 'both',
                 adaptive_bath_T_min: float = 0.1,
                adaptive_bath_T_max: Optional[float] = None,
                adaptive_bath_empirical_data_file: Optional[str] = None,
                adaptive_bath_signal_temperature_method: str = 'harmonic_equipartition',
                adaptive_bath_dynamic_target_temperature_method: Optional[str] = None,
                # Quench controller parameters
                 enable_quench_controller: bool = False,
                 quench_initial_temperature: float = 100.0,
                 quench_target_temperature: float = 50.0,
                 quench_time_ps: float = 50.0,
                 quench_apply_to: str = 'both',
                 # Offset temperature controller parameters
                 enable_offset_controller: bool = False,
                 offset_temperature_method: str = 'kinetic',
                 offset_temperature_offset_K: float = -50.0,
                 offset_turn_on_time_ps: float = 0.0,
                 offset_turn_off_time_ps: Optional[float] = None,
                 offset_update_interval_ps: float = 0.1,
                 offset_apply_to: str = 'both',
                 offset_T_min: float = 0.0,
                 offset_T_max: Optional[float] = None,
                 # Harmonic bond reset parameters
                 enable_harmonic_reset: bool = False,
                 harmonic_reset_turn_on_time_ps: float = 0.0,
                 harmonic_reset_temperature: Optional[float] = None,
                 harmonic_reset_seed: int = 42,
                 # Differential equation controller parameters
                 enable_diffeq_controller: bool = False,
                 diffeq_temperature_method: str = 'kinetic',
                 diffeq_time_constant_ps: float = 5.0,
                 diffeq_time_constant_auto: bool = False,
                 diffeq_turn_on_time_ps: float = 0.0,
                 diffeq_turn_off_time_ps: Optional[float] = None,
                 diffeq_update_interval_ps: float = 0.1,
                 diffeq_apply_to: str = 'both',
                 diffeq_T_min: float = 0.0,
                 diffeq_T_max: Optional[float] = None,
                 diffeq_rate_limit_K_per_ps: Optional[float] = None,
                 diffeq_disable_bias_estimation: bool = False,
                 # DiffEq Kalman filter bias estimation parameters
                 diffeq_enable_bias_estimation: bool = True,
                 diffeq_bias_process_noise: float = 1e-6,
                 diffeq_bias_initial_covariance: float = 100.0,
                 # DiffEq auto-tuning parameters
                 diffeq_enable_autotuning: bool = False,
                 diffeq_autotune_duration_ps: float = 50.0,
                 diffeq_autotune_perturbation_amplitude_K: float = 2.0,
                 diffeq_autotune_measurement_noise_K: float = 1.0,
                 diffeq_autotune_min_bias_timescale_factor: float = 10.0,
                 diffeq_periodic_coupling_period_ps: Optional[float] = None,
                 # DiffEq bias correction parameters (new clean interface)
                 diffeq_bias_correction_mode: str = 'none',
                 diffeq_equilibrium_temperature: Optional[float] = None,
                 # Exact-cancellation + PI control parameters
                 diffeq_enable_pi_control: bool = False,
                 diffeq_pi_rho: float = 1.0,
                 diffeq_pi_epsilon: float = 0.5,
                 diffeq_pi_zeta: float = 0.8,
                 diffeq_relaxation_data_file: Optional[str] = None,
                 diffeq_filter_window: float = 0.0,
                 # Adaptive bias cancellation parameters (mutually exclusive with PI control)
                 diffeq_enable_bias_cancellation: bool = False,
                 diffeq_bias_tau_b_ps: float = 50.0,
                 diffeq_bias_tau_b_auto: bool = False,
                 diffeq_bias_kappa: float = 0.01,
                 diffeq_bias_kappa_auto: bool = False,
                 diffeq_bias_tau_b_prefactor: float = 5.0,
                 diffeq_bias_kappa_prefactor: float = 50.0,
                 diffeq_bias_calibration_time_ps: float = 10.0,
                 # BathPI controller parameters
                 enable_bath_pi_controller: bool = False,
                 bath_pi_apply_to: str = 'both',
                 bath_pi_K_p_molecular: float = 0.1,
                 bath_pi_K_i_molecular: float = 0.01,
                 bath_pi_K_T_molecular: float = 0.0,
                 bath_pi_K_p_cavity: float = 0.1,
                 bath_pi_K_i_cavity: float = 0.01,
                 bath_pi_K_T_cavity: float = 0.0,
                 bath_pi_filter_window_ps: Union[float, str] = 5.0,
                 bath_pi_flux_source: str = 'reservoir',
                 bath_pi_anti_windup_alpha: float = 0.01,
                 bath_pi_enable_feedforward: bool = False,
                 bath_pi_T_nominal: Optional[float] = None,
                 bath_pi_feedforward_tau_ps: float = 1000.0,
                 bath_pi_turn_on_time_ps: float = 0.0,
                 bath_pi_turn_off_time_ps: Optional[float] = None,
                 bath_pi_update_interval_ps: float = 0.1,
                 bath_pi_T_min: float = 0.1,
                 bath_pi_T_max: Optional[float] = None,
                 bath_pi_rate_limit_K_per_ps: Optional[float] = None,
                 bath_pi_output_file: str = 'bath_pi_control.csv',
                 bath_pi_relaxation_data_file: Optional[str] = None,
                 # Simple Setpoint Controller parameters
                 enable_simple_setpoint_controller: bool = False,
                 simple_setpoint_signal_method: str = 'kinetic',
                 simple_setpoint_time_constant_ps: float = 5.0,
                 simple_setpoint_apply_to: str = 'both',
                 simple_setpoint_turn_on_time_ps: float = 0.0,
                 simple_setpoint_turn_off_time_ps: Optional[float] = None,
                 simple_setpoint_update_interval_ps: float = 0.1,
                 simple_setpoint_T_min: float = 0.0,
                 simple_setpoint_T_max: Optional[float] = None,
                 simple_setpoint_output_file: str = 'simple_setpoint_control.csv',
                 simple_setpoint_console_output_period_ps: float = 1.0,
                 # Adaptive MPC controller parameters
                 enable_adaptive_mpc_controller: bool = False,
                 adaptive_mpc_target_temperature: float = 100.0,
                 adaptive_mpc_turn_on_time_ps: float = 0.0,
                 adaptive_mpc_turn_off_time_ps: Optional[float] = None,
                 adaptive_mpc_system_id_duration_ps: float = 50.0,
                 adaptive_mpc_system_id_step_duration_ps: float = 5.0,
                 adaptive_mpc_system_id_seed: int = 42,
                 adaptive_mpc_update_interval_ps: float = 0.1,
                 adaptive_mpc_prediction_horizon: int = 10,
                 adaptive_mpc_control_horizon: int = 5,
                 adaptive_mpc_output_weight: float = 100.0,
                 adaptive_mpc_control_effort_weight: List[float] = None,
                 adaptive_mpc_rate_penalty_weight: List[float] = None,
                 adaptive_mpc_lambda_min: float = 0.0,
                 adaptive_mpc_lambda_max: float = 1e-2,
                 adaptive_mpc_T_bath_min: float = 0.1,
                 adaptive_mpc_T_bath_max: float = 500.0,
                 adaptive_mpc_delta_lambda_max: float = 1e-4,
                 adaptive_mpc_delta_T_bath_max: float = 10.0,
                 adaptive_mpc_apply_to: str = 'both',
                 adaptive_mpc_rls_forgetting_factor: float = 0.995,
                 adaptive_mpc_rls_initial_covariance: float = 100.0,
                 adaptive_mpc_model_update_interval: int = 10,
                 adaptive_mpc_output_file: str = 'adaptive_mpc_control.csv',
                 adaptive_mpc_console_output_period_ps: float = 1.0,
                 adaptive_mpc_regularization_param: float = 1e-3,
                 adaptive_mpc_use_scaling: bool = True,
                 adaptive_mpc_debug_mode: bool = False,
                 # PID controller parameters
                 enable_pid_controller: bool = False,
                 pid_signal_choice: str = 'lj_coulombic',
                 pid_target_temperature: float = 100.0,
                 pid_self_loop: bool = False,
                 pid_Kp: Optional[float] = None,
                 pid_Ti: Optional[float] = None,
                 pid_Td: Optional[float] = None,
                 pid_auto_tune: bool = True,
                 pid_auto_tune_step_size: float = 20.0,
                 pid_auto_tune_duration_ps: float = 50.0,
                 pid_turn_on_time_ps: float = 0.0,
                 pid_turn_off_time_ps: Optional[float] = None,
                 pid_update_interval_ps: float = 0.1,
                 pid_apply_to: str = 'both',
                 pid_T_min: float = 0.1,
                 pid_T_max: Optional[float] = None,
                 pid_rate_limit_K_per_ps: Optional[float] = None,
                 pid_derivative_filter_N: float = 10.0,
                 pid_enable_anti_windup: bool = True,
                 pid_output_file: str = 'pid_control.csv',
                 pid_console_output_period_ps: float = 1.0,
                 # LQR controller parameters
                 enable_lqr_controller: bool = False,
                 lqr_signal_method: str = 'lj_coulombic',
                 lqr_hot_method: str = 'harmonic_equipartition',
                 lqr_target_temperature: float = 300.0,
                 lqr_dynamic_target: bool = False,
                 lqr_dynamic_target_method: Optional[str] = None,
                 lqr_weight_signal: float = 100.0,
                 lqr_weight_hot: float = 1.0,
                 lqr_weight_bath: float = 0.1,
                 lqr_weight_integral: float = 10.0,
                 lqr_control_effort: float = 1.0,
                 lqr_process_noise_signal: float = 0.1,
                 lqr_process_noise_hot: float = 0.1,
                 lqr_measurement_noise_signal: float = 0.5,
                 lqr_measurement_noise_hot: float = 0.5,
                 lqr_system_id_mode: str = 'step',
                 lqr_system_id_temp_K: float = 5.0,
                 lqr_system_id_duration_ps: float = 50.0,
                 lqr_system_id_file: str = 'lqr_system_params.json',
                lqr_periodic_system_id: bool = False,
                lqr_periodic_system_id_interval_ps: float = 1000.0,
                # EKF-based adaptation parameters
                lqr_use_ekf_adaptation: bool = True,
                lqr_ekf_update_interval: int = 50,
                lqr_ekf_process_noise_param: float = 0.001,
                lqr_ekf_initial_covariance_param: float = 0.1,
                lqr_adaptive_lqr_threshold: float = 0.05,
                # Gain scheduling parameters
                lqr_enable_gain_scheduling: bool = True,
                lqr_gain_schedule_far_threshold: float = 20.0,
                lqr_gain_schedule_near_threshold: float = 10.0,
                # T_h low-pass filter parameters
                lqr_th_filter_enabled: bool = True,
                lqr_th_filter_time_constant: float = 20.0,
                # Gentle startup parameters
                lqr_gentle_startup_steps: int = 10,
                lqr_gentle_startup_min_authority: float = 0.1,
                # Kinetic temperature tracking (3D state augmentation)
                lqr_track_kinetic_temp: bool = False,
                lqr_weight_kinetic: float = 100.0,
                lqr_process_noise_kinetic: float = 2.0,
                lqr_measurement_noise_kinetic: float = 2.0,
                # Cross-coupling weights for thermal equilibration (penalize temperature differences)
                lqr_cross_coupling_signal_kinetic: float = 0.0,
                lqr_cross_coupling_signal_hot: float = 0.0,
                lqr_cross_coupling_hot_kinetic: float = 0.0,
                # Timing and limits
                lqr_turn_on_time_ps: float = 0.0,
                 lqr_turn_off_time_ps: Optional[float] = None,
                 lqr_update_interval_ps: float = 0.1,
                 lqr_T_min: float = 0.1,
                 lqr_T_max: Optional[float] = None,
                 lqr_apply_to: str = 'both',
                 lqr_output_file: str = 'lqr_controller.csv',
                 lqr_empirical_data_file: Optional[str] = None,
                 # Adaptive LQI controller parameters
                 lqr_controller_type: str = 'standard',
                 lqr_tau_L_initial: float = 200.0,
                 lqr_tau_H_initial: float = 30.0,
                 lqr_k_initial: float = 0.01,
                 lqr_tau_b: float = 1.0,
                 lqr_q_common: float = 1000.0,
                 lqr_q_diff: float = 100.0,
                 lqr_q_eta_common: float = 1e6,
                 lqr_q_eta_diff: float = 10.0,
                 lqr_process_noise_drift: float = 0.01,
                 lqr_rls_forgetting: float = 0.998,
                 lqr_rls_regularization: float = 1e-6,
                 lqr_rls_update_interval: int = 50,
                 lqr_max_control_rate: float = 10.0,
                 lqr_integral_max_common: float = 1000.0,
                 lqr_integral_max_diff: float = 100.0,
                 lqr_theta_change_threshold: float = 0.05,
                 # LQG coupling controller parameters (new coupling-type: lqg_coupling)
                 lqg_coupling_target_temperature: float = 100.0,
                 lqg_coupling_update_interval_ps: float = 0.1,
                 lqg_coupling_equilibrium_duration_ps: float = 5.0,
                 lqg_coupling_step_duration_ps: float = 5.0,
                 lqg_coupling_n_steps: int = 2,
                 lqg_coupling_lambda_min: float = 0.0,
                 lqg_coupling_lambda_max: float = 1e-2,
                 lqg_coupling_process_noise_std: float = 0.1,
                 lqg_coupling_measurement_noise_std: float = 0.5,
                 lqg_coupling_weight_signal: float = 100.0,
                 lqg_coupling_weight_harmonic: float = 1.0,
                 lqg_coupling_weight_kinetic: float = 1.0,
                 lqg_coupling_weight_bath: float = 0.1,
                 lqg_coupling_control_effort: float = 1.0,
                 lqg_coupling_integral_gain: float = 2.0,
                 lqg_coupling_system_id_file: str = 'lqg_coupling_system_id.json',
                 lqg_coupling_control_file: str = 'lqg_coupling_control.csv',
                 lqg_coupling_temperature_methods: List[str] = None,
                 lqg_coupling_bath_temperature_method: str = 'kinetic',
                 # Legacy LQG pulse controller parameters (for backward compatibility)
                 enable_lqg_controller: bool = False,
                 lqg_A_matrix=None, lqg_B_matrix=None, lqg_C_matrix=None,
                 lqg_S_matrix=None, lqg_r_target=None,
                 lqg_Qz_weights=None, lqg_Ru_weights=None, lqg_Qe_weights=None,
                 lqg_W_noise=None, lqg_V_noise=None,
                 lqg_g0_baseline: float = 1e-4, lqg_dt: float = 5.0,
                 lqg_delta_g_min: float = -0.3, lqg_delta_g_max: float = 0.3,
                 lqg_T_bath_min: float = 10.0, lqg_T_bath_max: float = 500.0,
                 lqg_rate_limit_K_per_ps: float = 1.0,
                 lqg_turn_on_time_ps: float = 0.0, lqg_turn_off_time_ps: Optional[float] = None,
                 lqg_update_interval_ps: float = 5.0, lqg_console_output_period_ps: float = 10.0,
                 lqg_c_affine_term=None, lqg_temperature_methods=None,
                 lqg_square_period_ps: float = 5.0, lqg_square_duty_cycle: float = 0.5,
                 lqg_square_phase_offset: float = 0.0, lqg_square_start_time_ps: float = 0.0,
                 lqg_square_stop_time_ps: Optional[float] = None, lqg_square_g_low_level: float = 1e-6,
                 lqg_enable_system_id: bool = True, lqg_equilibrium_duration_ps: float = 500.0,
                 lqg_prbs_n_unique_steps: int = 10, lqg_prbs_amplitude: float = 5e-5,
                 lqg_bath_n_unique_steps: int = 8, lqg_bath_excitation_amplitude: float = 10.0,
                lqg_system_id_output_file: str = 'lqg_system_id.json',
                enable_temp_tracker: bool = True,
                temp_tracker_output_period_ps: float = 0.1,
                temp_tracker_empirical_data_file: Optional[str] = None,
                # HDF5 observable output parameters
                enable_hdf5_output: bool = True,
                 hdf5_output_file: Optional[str] = None,
                 hdf5_output_period_ps: float = 0.01,
                 # Molecular temperature decomposition parameters
                 enable_molecular_temps: bool = False,
                 molecular_temps_output_period_ps: float = 1.0,
                 # Dipole moment FDR parameters
                 enable_dipole_fdr: bool = False,
                 dipole_fdr_output_period_ps: float = 0.1,
                 dipole_fdr_max_correlation_time_ps: float = 100.0,
                 dipole_fdr_field_direction: List[float] = [0, 0, 1],
                 dipole_fdr_exclude_cavity: bool = True,
                 enable_dipole_response: bool = False,
                 dipole_response_field_strength: float = 1e-5,
                 dipole_response_sign: float = 1.0,
                 # Mechanical cavity modulation parameters
                 mech_periodic: bool = False,
                 mech_frequency_cm1: float = 100.0,
                 mech_magnitude: float = 1e-4,
                 # Auto-stopping parameters
                 enable_auto_stopping: bool = False,
                 auto_stop_min_time_ps: float = 10.0,
                 auto_stop_smoothing_window: int = 5,
) -> None:
        """Initialize the CavityMDSimulation with simulation parameters."""
        self.job_dir = job_dir
        self.replica = replica
        self.freq = freq
        
        # Calculate omegac for lambda-to-epsilon conversion
        from ..utils import PhysicalConstants
        self.omegac = self.freq / PhysicalConstants.HARTREE_TO_CM_MINUS1

        # Handle coupling constant: lambda_coupling is now the primary parameter
        if couplstr is not None:
            import warnings
            warnings.warn(
                "The 'couplstr' parameter is DEPRECATED and will be removed in a future version. "
                "Please use 'lambda_coupling' instead. "
                "The relationship is: epsilon = lambda_coupling * omega_c, "
                "where omega_c is the cavity characteristic frequency.",
                DeprecationWarning,
                stacklevel=2
            )
            if lambda_coupling is not None:
                raise ValueError(
                    "Cannot specify both 'couplstr' (deprecated) and 'lambda_coupling'. "
                    "Please use only 'lambda_coupling' going forward."
                )
            # For backward compatibility, treat couplstr as the effective coupling (not scaled)
            self.lambda_coupling = None
            self.couplstr = couplstr
            self.use_deprecated_couplstr = True
        else:
            if lambda_coupling is None:
                raise ValueError(
                    "Must specify 'lambda_coupling' (dimensionless coupling parameter). "
                    "The effective coupling strength will be computed as: epsilon = lambda_coupling * omega_c"
                )
            self.lambda_coupling = lambda_coupling
            self.couplstr = None
            self.use_deprecated_couplstr = False

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
        
        # Dynamic coupling detection parameters
        self.enable_dynamic_coupling_detection = enable_dynamic_coupling_detection
        self.coupling_change_threshold = coupling_change_threshold
        
        # Momentum zeroing parameters
        self.zero_momentum_enabled = zero_momentum_enabled
        self.zero_momentum_period_ps = zero_momentum_period_ps
        
        # Empirical temperature feedback parameters
        self.enable_empirical_feedback = enable_empirical_feedback
        self.empirical_data_file = empirical_data_file
        self.feedback_output_period_ps = feedback_output_period_ps
        self.feedback_energy_component = feedback_energy_component
        self.feedback_apply_to = feedback_apply_to
        self.feedback_averaging_window_ps = feedback_averaging_window_ps
        self.feedback_update_interval_ps = feedback_update_interval_ps
        self.feedback_T_min = feedback_T_min
        self.feedback_T_max = feedback_T_max
        self.feedback_turn_off_time_ps = feedback_turn_off_time_ps
        self.feedback_use_direct_harmonic_calculation = feedback_use_direct_harmonic_calculation
        self.feedback_auto_adjust_molecular_tau = feedback_auto_adjust_molecular_tau
        
        # Periodic coupling parameters  
        self.periodic_coupling = periodic_coupling
        self.periodic_period_ps = periodic_period_ps
        self.periodic_phase_offset = periodic_phase_offset
        self.periodic_start_time_ps = periodic_start_time_ps
        self.periodic_stop_time_ps = periodic_stop_time_ps
        
        # Enhanced coupling variant parameters
        self.coupling_variant_type = coupling_variant_type
        self.exponential_amplitude = exponential_amplitude if exponential_amplitude is not None else couplstr
        self.exponential_decay_time_ps = exponential_decay_time_ps
        self.exponential_turn_on_time_ps = exponential_turn_on_time_ps
        self.exponential_turn_off_time_ps = exponential_turn_off_time_ps
        self.square_period_ps = square_period_ps
        self.square_duty_cycle = square_duty_cycle
        self.square_phase_offset = square_phase_offset
        self.square_start_time_ps = square_start_time_ps
        self.square_stop_time_ps = square_stop_time_ps
        
        # Decaying square wave parameters
        self.decaying_square_period_ps = decaying_square_period_ps
        self.decaying_square_duty_cycle = decaying_square_duty_cycle
        self.decaying_square_phase_offset = decaying_square_phase_offset
        self.decaying_square_start_time_ps = decaying_square_start_time_ps
        self.decaying_square_stop_time_ps = decaying_square_stop_time_ps
        self.decaying_square_decay_rate = decaying_square_decay_rate
        self.decaying_square_minimum_amplitude = decaying_square_minimum_amplitude
        
        # Adaptive square wave parameters
        self.adaptive_square_period_ps = adaptive_square_period_ps
        self.adaptive_square_duty_cycle = adaptive_square_duty_cycle
        self.adaptive_square_phase_offset = adaptive_square_phase_offset
        self.adaptive_square_start_time_ps = adaptive_square_start_time_ps
        self.adaptive_square_stop_time_ps = adaptive_square_stop_time_ps
        self.adaptive_square_min_amplitude = adaptive_square_min_amplitude
        self.adaptive_square_max_amplitude = adaptive_square_max_amplitude
        
        # Exponential wave parameters
        self.exp_period_ps = exp_period_ps
        self.exp_tau_ps = exp_tau_ps
        self.exp_start_time_ps = exp_start_time_ps
        self.exp_stop_time_ps = exp_stop_time_ps
        self.exp_adaptive = exp_adaptive
        
        # Composite coupling parameters
        self.composite_sinusoid_amplitude = composite_sinusoid_amplitude
        self.composite_sinusoid_period = composite_sinusoid_period
        self.composite_sinusoid_phase = composite_sinusoid_phase
        self.composite_sinusoid_start_time = composite_sinusoid_start_time
        self.composite_sinusoid_stop_time = composite_sinusoid_stop_time
        self.composite_square_amplitude = composite_square_amplitude
        self.composite_square_period = composite_square_period
        self.composite_square_duty_cycle = composite_square_duty_cycle
        self.composite_square_start_time = composite_square_start_time
        self.composite_square_stop_time = composite_square_stop_time
        self.composite_square_adaptive = composite_square_adaptive
        self.composite_max_amplitude = composite_max_amplitude
        
        # Auto-stop coupling parameters
        self.enable_auto_stop = enable_auto_stop
        self.auto_stop_tol = auto_stop_tol
        self.auto_stop_window = auto_stop_window
        
        # Enhanced laser drive parameters
        self.laser_enabled = laser_enabled
        self.laser_frequency_cm1 = laser_frequency_cm1
        self.laser_amplitude = laser_amplitude
        self.laser_start_time_ps = laser_start_time_ps
        self.laser_stop_time_ps = laser_stop_time_ps
        self.laser_kvector = laser_kvector if laser_kvector is not None else [1.0, 0.0, 0.0]
        
        # PI feedback controller parameters
        self.enable_pi_feedback = enable_pi_feedback
        self.pi_target_temperature = pi_target_temperature
        self.pi_turn_on_time_ps = pi_turn_on_time_ps
        self.pi_turn_off_time_ps = pi_turn_off_time_ps
        self.pi_temperature_method = pi_temperature_method
        self.pi_update_interval_ps = pi_update_interval_ps
        self.pi_Kc = pi_Kc
        self.pi_Ti = pi_Ti
        self.pi_molecular_tau_ps = pi_molecular_tau_ps
        self.pi_beta = pi_beta
        self.pi_T_min = pi_T_min
        self.pi_T_max = pi_T_max
        
        # Gradient descent feedback controller parameters
        self.enable_gd_feedback = enable_gd_feedback
        self.gd_target_temperature = gd_target_temperature
        self.gd_turn_on_time_ps = gd_turn_on_time_ps
        self.gd_turn_off_time_ps = gd_turn_off_time_ps
        self.gd_temperature_method = gd_temperature_method
        self.gd_update_interval_ps = gd_update_interval_ps
        self.gd_time_constant_ps = gd_time_constant_ps
        self.gd_apply_to = gd_apply_to
        self.gd_T_min = gd_T_min
        self.gd_T_max = gd_T_max
        self.gd_disable_effective_temp = gd_disable_effective_temp
        
        # Multi-signal error function parameters
        self.gd_enable_multi_signal = gd_enable_multi_signal
        self.gd_weight_system_target = gd_weight_system_target
        self.gd_weight_bath_target = gd_weight_bath_target
        self.gd_weight_system_bath = gd_weight_system_bath
        
        # Dual independent feedback controller parameters
        self.enable_dual_feedback = enable_dual_feedback
        self.dual_cavity_method = dual_cavity_method
        self.dual_molecular_method = dual_molecular_method
        self.dual_cavity_target_temperature = dual_cavity_target_temperature
        self.dual_molecular_target_temperature = dual_molecular_target_temperature
        self.dual_cavity_time_constant_ps = dual_cavity_time_constant_ps
        self.dual_molecular_time_constant_ps = dual_molecular_time_constant_ps
        self.dual_turn_on_time_ps = dual_turn_on_time_ps
        self.dual_turn_off_time_ps = dual_turn_off_time_ps
        self.dual_update_interval_ps = dual_update_interval_ps
        self.dual_cavity_T_min = dual_cavity_T_min
        self.dual_cavity_T_max = dual_cavity_T_max
        self.dual_molecular_T_min = dual_molecular_T_min
        self.dual_molecular_T_max = dual_molecular_T_max
        self.dual_cavity_dynamic_target = dual_cavity_dynamic_target
        self.dual_molecular_dynamic_target = dual_molecular_dynamic_target
        self.dual_cavity_integral_time_constant_ps = dual_cavity_integral_time_constant_ps
        self.dual_molecular_integral_time_constant_ps = dual_molecular_integral_time_constant_ps
        
        # Sinusoidal bath temperature controller parameters
        self.enable_sinusoidal_bath = enable_sinusoidal_bath
        self.sinusoidal_bath_period_ps = sinusoidal_bath_period_ps
        self.sinusoidal_bath_amplitude_scale = sinusoidal_bath_amplitude_scale
        self.sinusoidal_bath_phase_offset = sinusoidal_bath_phase_offset
        self.sinusoidal_bath_target_temperature = sinusoidal_bath_target_temperature
        self.sinusoidal_bath_dynamic_target = sinusoidal_bath_dynamic_target
        self.sinusoidal_bath_turn_on_time_ps = sinusoidal_bath_turn_on_time_ps
        self.sinusoidal_bath_turn_off_time_ps = sinusoidal_bath_turn_off_time_ps
        self.sinusoidal_bath_update_interval_ps = sinusoidal_bath_update_interval_ps
        self.sinusoidal_bath_apply_to = sinusoidal_bath_apply_to
        self.sinusoidal_bath_T_min = sinusoidal_bath_T_min
        self.sinusoidal_bath_T_max = sinusoidal_bath_T_max
        self.sinusoidal_bath_empirical_data_file = sinusoidal_bath_empirical_data_file
        self.sinusoidal_bath_amplitude_update_interval_ps = sinusoidal_bath_amplitude_update_interval_ps
        self.sinusoidal_bath_amplitude_temperature_method = sinusoidal_bath_amplitude_temperature_method
        self.sinusoidal_bath_adaptive_range_mode = sinusoidal_bath_adaptive_range_mode
        
        # Adaptive bath temperature controller parameters  
        self.enable_adaptive_bath = enable_adaptive_bath
        self.adaptive_bath_amplitude_scale = adaptive_bath_amplitude_scale
        self.adaptive_bath_time_constant_ps = adaptive_bath_time_constant_ps
        self.adaptive_bath_target_temperature = adaptive_bath_target_temperature
        self.adaptive_bath_dynamic_target = adaptive_bath_dynamic_target
        self.adaptive_bath_turn_on_time_ps = adaptive_bath_turn_on_time_ps
        self.adaptive_bath_turn_off_time_ps = adaptive_bath_turn_off_time_ps
        self.adaptive_bath_update_interval_ps = adaptive_bath_update_interval_ps
        self.adaptive_bath_T_min = adaptive_bath_T_min
        self.adaptive_bath_T_max = adaptive_bath_T_max
        self.adaptive_bath_apply_to = adaptive_bath_apply_to
        self.adaptive_bath_empirical_data_file = adaptive_bath_empirical_data_file
        self.adaptive_bath_signal_temperature_method = adaptive_bath_signal_temperature_method
        self.adaptive_bath_dynamic_target_temperature_method = adaptive_bath_dynamic_target_temperature_method
        
        # Quench controller parameters
        self.enable_quench_controller = enable_quench_controller
        self.quench_initial_temperature = quench_initial_temperature
        self.quench_target_temperature = quench_target_temperature
        self.quench_time_ps = quench_time_ps
        self.quench_apply_to = quench_apply_to
        
        # Offset temperature controller parameters
        self.enable_offset_controller = enable_offset_controller
        self.offset_temperature_method = offset_temperature_method
        self.offset_temperature_offset_K = offset_temperature_offset_K
        self.offset_turn_on_time_ps = offset_turn_on_time_ps
        self.offset_turn_off_time_ps = offset_turn_off_time_ps
        self.offset_update_interval_ps = offset_update_interval_ps
        self.offset_apply_to = offset_apply_to
        self.offset_T_min = offset_T_min
        self.offset_T_max = offset_T_max
        
        # Harmonic bond reset parameters
        self.enable_harmonic_reset = enable_harmonic_reset
        self.harmonic_reset_turn_on_time_ps = harmonic_reset_turn_on_time_ps
        self.harmonic_reset_temperature = harmonic_reset_temperature
        self.harmonic_reset_seed = harmonic_reset_seed
        
        # Differential equation controller parameters
        self.enable_diffeq_controller = enable_diffeq_controller
        self.diffeq_temperature_method = diffeq_temperature_method
        self.diffeq_time_constant_ps = diffeq_time_constant_ps
        self.diffeq_time_constant_auto = diffeq_time_constant_auto
        self.diffeq_turn_on_time_ps = diffeq_turn_on_time_ps
        self.diffeq_turn_off_time_ps = diffeq_turn_off_time_ps
        self.diffeq_update_interval_ps = diffeq_update_interval_ps
        self.diffeq_apply_to = diffeq_apply_to
        self.diffeq_T_min = diffeq_T_min
        self.diffeq_T_max = diffeq_T_max
        self.diffeq_rate_limit_K_per_ps = diffeq_rate_limit_K_per_ps
        self.diffeq_disable_bias_estimation = diffeq_disable_bias_estimation
        # Kalman filter bias estimation
        self.diffeq_enable_bias_estimation = diffeq_enable_bias_estimation
        self.diffeq_bias_process_noise = diffeq_bias_process_noise
        self.diffeq_bias_initial_covariance = diffeq_bias_initial_covariance
        # Auto-tuning
        self.diffeq_enable_autotuning = diffeq_enable_autotuning
        self.diffeq_autotune_duration_ps = diffeq_autotune_duration_ps
        self.diffeq_autotune_perturbation_amplitude_K = diffeq_autotune_perturbation_amplitude_K
        self.diffeq_autotune_measurement_noise_K = diffeq_autotune_measurement_noise_K
        self.diffeq_autotune_min_bias_timescale_factor = diffeq_autotune_min_bias_timescale_factor
        self.diffeq_periodic_coupling_period_ps = diffeq_periodic_coupling_period_ps
        # DiffEq bias correction parameters
        self.diffeq_bias_correction_mode = diffeq_bias_correction_mode
        self.diffeq_equilibrium_temperature = diffeq_equilibrium_temperature
        # Exact-cancellation + PI control parameters
        self.diffeq_enable_pi_control = diffeq_enable_pi_control
        self.diffeq_pi_rho = diffeq_pi_rho
        self.diffeq_pi_epsilon = diffeq_pi_epsilon
        self.diffeq_pi_zeta = diffeq_pi_zeta
        self.diffeq_relaxation_data_file = diffeq_relaxation_data_file
        self.diffeq_filter_window = diffeq_filter_window
        # Adaptive bias cancellation parameters
        self.diffeq_enable_bias_cancellation = diffeq_enable_bias_cancellation
        self.diffeq_bias_tau_b_ps = diffeq_bias_tau_b_ps
        self.diffeq_bias_tau_b_auto = diffeq_bias_tau_b_auto
        self.diffeq_bias_kappa = diffeq_bias_kappa
        self.diffeq_bias_kappa_auto = diffeq_bias_kappa_auto
        self.diffeq_bias_tau_b_prefactor = diffeq_bias_tau_b_prefactor
        self.diffeq_bias_kappa_prefactor = diffeq_bias_kappa_prefactor
        self.diffeq_bias_calibration_time_ps = diffeq_bias_calibration_time_ps
        
        # BathPI controller parameters
        self.enable_bath_pi_controller = enable_bath_pi_controller
        self.bath_pi_apply_to = bath_pi_apply_to
        self.bath_pi_K_p_molecular = bath_pi_K_p_molecular
        self.bath_pi_K_i_molecular = bath_pi_K_i_molecular
        self.bath_pi_K_T_molecular = bath_pi_K_T_molecular
        self.bath_pi_K_p_cavity = bath_pi_K_p_cavity
        self.bath_pi_K_i_cavity = bath_pi_K_i_cavity
        self.bath_pi_K_T_cavity = bath_pi_K_T_cavity
        self.bath_pi_filter_window_ps = bath_pi_filter_window_ps
        self.bath_pi_flux_source = bath_pi_flux_source
        self.bath_pi_anti_windup_alpha = bath_pi_anti_windup_alpha
        self.bath_pi_enable_feedforward = bath_pi_enable_feedforward
        self.bath_pi_T_nominal = bath_pi_T_nominal
        self.bath_pi_feedforward_tau_ps = bath_pi_feedforward_tau_ps
        self.bath_pi_turn_on_time_ps = bath_pi_turn_on_time_ps
        self.bath_pi_turn_off_time_ps = bath_pi_turn_off_time_ps
        self.bath_pi_update_interval_ps = bath_pi_update_interval_ps
        self.bath_pi_T_min = bath_pi_T_min
        self.bath_pi_T_max = bath_pi_T_max
        self.bath_pi_rate_limit_K_per_ps = bath_pi_rate_limit_K_per_ps
        self.bath_pi_output_file = bath_pi_output_file
        self.bath_pi_relaxation_data_file = bath_pi_relaxation_data_file
        
        # Simple Setpoint Controller parameters
        self.enable_simple_setpoint_controller = enable_simple_setpoint_controller
        self.simple_setpoint_signal_method = simple_setpoint_signal_method
        self.simple_setpoint_time_constant_ps = simple_setpoint_time_constant_ps
        self.simple_setpoint_apply_to = simple_setpoint_apply_to
        self.simple_setpoint_turn_on_time_ps = simple_setpoint_turn_on_time_ps
        self.simple_setpoint_turn_off_time_ps = simple_setpoint_turn_off_time_ps
        self.simple_setpoint_update_interval_ps = simple_setpoint_update_interval_ps
        self.simple_setpoint_T_min = simple_setpoint_T_min
        self.simple_setpoint_T_max = simple_setpoint_T_max
        self.simple_setpoint_output_file = simple_setpoint_output_file
        self.simple_setpoint_console_output_period_ps = simple_setpoint_console_output_period_ps
        
        # Adaptive MPC controller parameters
        self.enable_adaptive_mpc_controller = enable_adaptive_mpc_controller
        self.adaptive_mpc_target_temperature = adaptive_mpc_target_temperature
        self.adaptive_mpc_turn_on_time_ps = adaptive_mpc_turn_on_time_ps
        self.adaptive_mpc_turn_off_time_ps = adaptive_mpc_turn_off_time_ps
        self.adaptive_mpc_system_id_duration_ps = adaptive_mpc_system_id_duration_ps
        self.adaptive_mpc_system_id_step_duration_ps = adaptive_mpc_system_id_step_duration_ps
        self.adaptive_mpc_system_id_seed = adaptive_mpc_system_id_seed
        self.adaptive_mpc_update_interval_ps = adaptive_mpc_update_interval_ps
        self.adaptive_mpc_prediction_horizon = adaptive_mpc_prediction_horizon
        self.adaptive_mpc_control_horizon = adaptive_mpc_control_horizon
        self.adaptive_mpc_output_weight = adaptive_mpc_output_weight
        self.adaptive_mpc_control_effort_weight = adaptive_mpc_control_effort_weight if adaptive_mpc_control_effort_weight else [1.0, 0.1]
        self.adaptive_mpc_rate_penalty_weight = adaptive_mpc_rate_penalty_weight if adaptive_mpc_rate_penalty_weight else [10.0, 1.0]
        self.adaptive_mpc_lambda_min = adaptive_mpc_lambda_min
        self.adaptive_mpc_lambda_max = adaptive_mpc_lambda_max
        self.adaptive_mpc_T_bath_min = adaptive_mpc_T_bath_min
        self.adaptive_mpc_T_bath_max = adaptive_mpc_T_bath_max
        self.adaptive_mpc_delta_lambda_max = adaptive_mpc_delta_lambda_max
        self.adaptive_mpc_delta_T_bath_max = adaptive_mpc_delta_T_bath_max
        self.adaptive_mpc_apply_to = adaptive_mpc_apply_to
        self.adaptive_mpc_rls_forgetting_factor = adaptive_mpc_rls_forgetting_factor
        self.adaptive_mpc_rls_initial_covariance = adaptive_mpc_rls_initial_covariance
        self.adaptive_mpc_model_update_interval = adaptive_mpc_model_update_interval
        self.adaptive_mpc_output_file = adaptive_mpc_output_file
        self.adaptive_mpc_console_output_period_ps = adaptive_mpc_console_output_period_ps
        self.adaptive_mpc_regularization_param = adaptive_mpc_regularization_param
        self.adaptive_mpc_use_scaling = adaptive_mpc_use_scaling
        self.adaptive_mpc_debug_mode = adaptive_mpc_debug_mode
        
        # PID controller parameters
        self.enable_pid_controller = enable_pid_controller
        self.pid_signal_choice = pid_signal_choice
        self.pid_target_temperature = pid_target_temperature
        self.pid_self_loop = pid_self_loop
        self.pid_Kp = pid_Kp
        self.pid_Ti = pid_Ti
        self.pid_Td = pid_Td
        self.pid_auto_tune = pid_auto_tune
        self.pid_auto_tune_step_size = pid_auto_tune_step_size
        self.pid_auto_tune_duration_ps = pid_auto_tune_duration_ps
        self.pid_turn_on_time_ps = pid_turn_on_time_ps
        self.pid_turn_off_time_ps = pid_turn_off_time_ps
        self.pid_update_interval_ps = pid_update_interval_ps
        self.pid_apply_to = pid_apply_to
        self.pid_T_min = pid_T_min
        self.pid_T_max = pid_T_max
        self.pid_rate_limit_K_per_ps = pid_rate_limit_K_per_ps
        self.pid_derivative_filter_N = pid_derivative_filter_N
        self.pid_enable_anti_windup = pid_enable_anti_windup
        self.pid_output_file = pid_output_file
        self.pid_console_output_period_ps = pid_console_output_period_ps
        
        # LQR controller parameters
        self.enable_lqr_controller = enable_lqr_controller
        self.lqr_signal_method = lqr_signal_method
        self.lqr_hot_method = lqr_hot_method
        self.lqr_target_temperature = lqr_target_temperature
        self.lqr_dynamic_target = lqr_dynamic_target
        self.lqr_dynamic_target_method = lqr_dynamic_target_method
        self.lqr_weight_signal = lqr_weight_signal
        self.lqr_weight_hot = lqr_weight_hot
        self.lqr_weight_bath = lqr_weight_bath
        self.lqr_weight_integral = lqr_weight_integral
        self.lqr_control_effort = lqr_control_effort
        self.lqr_process_noise_signal = lqr_process_noise_signal
        self.lqr_process_noise_hot = lqr_process_noise_hot
        self.lqr_measurement_noise_signal = lqr_measurement_noise_signal
        self.lqr_measurement_noise_hot = lqr_measurement_noise_hot
        self.lqr_system_id_mode = lqr_system_id_mode
        self.lqr_system_id_temp_K = lqr_system_id_temp_K
        self.lqr_system_id_duration_ps = lqr_system_id_duration_ps
        self.lqr_system_id_file = lqr_system_id_file
        self.lqr_periodic_system_id = lqr_periodic_system_id
        self.lqr_periodic_system_id_interval_ps = lqr_periodic_system_id_interval_ps
        # EKF-based adaptation parameters
        self.lqr_use_ekf_adaptation = lqr_use_ekf_adaptation
        self.lqr_ekf_update_interval = lqr_ekf_update_interval
        self.lqr_ekf_process_noise_param = lqr_ekf_process_noise_param
        self.lqr_ekf_initial_covariance_param = lqr_ekf_initial_covariance_param
        self.lqr_adaptive_lqr_threshold = lqr_adaptive_lqr_threshold
        # Gain scheduling parameters
        self.lqr_enable_gain_scheduling = lqr_enable_gain_scheduling
        self.lqr_gain_schedule_far_threshold = lqr_gain_schedule_far_threshold
        self.lqr_gain_schedule_near_threshold = lqr_gain_schedule_near_threshold
        # T_h low-pass filter parameters
        self.lqr_th_filter_enabled = lqr_th_filter_enabled
        self.lqr_th_filter_time_constant = lqr_th_filter_time_constant
        # Gentle startup parameters
        self.lqr_gentle_startup_steps = lqr_gentle_startup_steps
        self.lqr_gentle_startup_min_authority = lqr_gentle_startup_min_authority
        # Kinetic temperature tracking (3D state augmentation)
        self.lqr_track_kinetic_temp = lqr_track_kinetic_temp
        self.lqr_weight_kinetic = lqr_weight_kinetic
        self.lqr_process_noise_kinetic = lqr_process_noise_kinetic
        self.lqr_measurement_noise_kinetic = lqr_measurement_noise_kinetic
        # Cross-coupling weights
        self.lqr_cross_coupling_signal_kinetic = lqr_cross_coupling_signal_kinetic
        self.lqr_cross_coupling_signal_hot = lqr_cross_coupling_signal_hot
        self.lqr_cross_coupling_hot_kinetic = lqr_cross_coupling_hot_kinetic
        # Timing and limits
        self.lqr_turn_on_time_ps = lqr_turn_on_time_ps
        self.lqr_turn_off_time_ps = lqr_turn_off_time_ps
        self.lqr_update_interval_ps = lqr_update_interval_ps
        self.lqr_T_min = lqr_T_min
        self.lqr_T_max = lqr_T_max
        self.lqr_apply_to = lqr_apply_to
        self.lqr_output_file = lqr_output_file
        self.lqr_empirical_data_file = lqr_empirical_data_file
        # Adaptive LQI controller parameters
        self.lqr_controller_type = lqr_controller_type
        self.lqr_tau_L_initial = lqr_tau_L_initial
        self.lqr_tau_H_initial = lqr_tau_H_initial
        self.lqr_k_initial = lqr_k_initial
        self.lqr_tau_b = lqr_tau_b
        self.lqr_q_common = lqr_q_common
        self.lqr_q_diff = lqr_q_diff
        self.lqr_q_eta_common = lqr_q_eta_common
        self.lqr_q_eta_diff = lqr_q_eta_diff
        self.lqr_process_noise_drift = lqr_process_noise_drift
        self.lqr_rls_forgetting = lqr_rls_forgetting
        self.lqr_rls_regularization = lqr_rls_regularization
        self.lqr_rls_update_interval = lqr_rls_update_interval
        self.lqr_max_control_rate = lqr_max_control_rate
        self.lqr_integral_max_common = lqr_integral_max_common
        self.lqr_integral_max_diff = lqr_integral_max_diff
        self.lqr_theta_change_threshold = lqr_theta_change_threshold
        
        # LQG coupling controller parameters
        self.lqg_coupling_target_temperature = lqg_coupling_target_temperature
        self.lqg_coupling_update_interval_ps = lqg_coupling_update_interval_ps
        self.lqg_coupling_equilibrium_duration_ps = lqg_coupling_equilibrium_duration_ps
        self.lqg_coupling_step_duration_ps = lqg_coupling_step_duration_ps
        self.lqg_coupling_n_steps = lqg_coupling_n_steps
        self.lqg_coupling_lambda_min = lqg_coupling_lambda_min
        self.lqg_coupling_lambda_max = lqg_coupling_lambda_max
        self.lqg_coupling_process_noise_std = lqg_coupling_process_noise_std
        self.lqg_coupling_measurement_noise_std = lqg_coupling_measurement_noise_std
        self.lqg_coupling_weight_signal = lqg_coupling_weight_signal
        self.lqg_coupling_weight_harmonic = lqg_coupling_weight_harmonic
        self.lqg_coupling_weight_kinetic = lqg_coupling_weight_kinetic
        self.lqg_coupling_weight_bath = lqg_coupling_weight_bath
        self.lqg_coupling_control_effort = lqg_coupling_control_effort
        self.lqg_coupling_integral_gain = lqg_coupling_integral_gain
        self.lqg_coupling_system_id_file = lqg_coupling_system_id_file
        self.lqg_coupling_control_file = lqg_coupling_control_file
        self.lqg_coupling_temperature_methods = lqg_coupling_temperature_methods or ['lj_coulombic', 'harmonic_equipartition', 'kinetic']
        self.lqg_coupling_bath_temperature_method = lqg_coupling_bath_temperature_method
        
        # Store legacy LQG pulse controller parameters
        self.enable_lqg_controller = enable_lqg_controller
        self.lqg_A_matrix = lqg_A_matrix
        self.lqg_B_matrix = lqg_B_matrix
        self.lqg_C_matrix = lqg_C_matrix
        self.lqg_S_matrix = lqg_S_matrix
        self.lqg_r_target = lqg_r_target
        self.lqg_Qz_weights = lqg_Qz_weights
        self.lqg_Ru_weights = lqg_Ru_weights
        self.lqg_Qe_weights = lqg_Qe_weights
        self.lqg_W_noise = lqg_W_noise
        self.lqg_V_noise = lqg_V_noise
        self.lqg_g0_baseline = lqg_g0_baseline
        self.lqg_dt = lqg_dt
        self.lqg_delta_g_min = lqg_delta_g_min
        self.lqg_delta_g_max = lqg_delta_g_max
        self.lqg_T_bath_min = lqg_T_bath_min
        self.lqg_T_bath_max = lqg_T_bath_max
        self.lqg_rate_limit_K_per_ps = lqg_rate_limit_K_per_ps
        self.lqg_turn_on_time_ps = lqg_turn_on_time_ps
        self.lqg_turn_off_time_ps = lqg_turn_off_time_ps
        self.lqg_update_interval_ps = lqg_update_interval_ps
        self.lqg_console_output_period_ps = lqg_console_output_period_ps
        self.lqg_c_affine_term = lqg_c_affine_term
        self.lqg_temperature_methods = lqg_temperature_methods
        self.lqg_square_period_ps = lqg_square_period_ps
        self.lqg_square_duty_cycle = lqg_square_duty_cycle
        self.lqg_square_phase_offset = lqg_square_phase_offset
        self.lqg_square_start_time_ps = lqg_square_start_time_ps
        self.lqg_square_stop_time_ps = lqg_square_stop_time_ps
        self.lqg_square_g_low_level = lqg_square_g_low_level
        self.lqg_enable_system_id = lqg_enable_system_id
        self.lqg_equilibrium_duration_ps = lqg_equilibrium_duration_ps
        self.lqg_prbs_n_unique_steps = lqg_prbs_n_unique_steps
        self.lqg_prbs_amplitude = lqg_prbs_amplitude
        self.lqg_bath_n_unique_steps = lqg_bath_n_unique_steps
        self.lqg_bath_excitation_amplitude = lqg_bath_excitation_amplitude
        self.lqg_system_id_output_file = lqg_system_id_output_file
        
        self.enable_temp_tracker = enable_temp_tracker
        self.temp_tracker_output_period_ps = temp_tracker_output_period_ps
        self.temp_tracker_empirical_data_file = temp_tracker_empirical_data_file
        
        # HDF5 output parameters
        self.enable_hdf5_output = enable_hdf5_output
        self.hdf5_output_file = hdf5_output_file
        self.hdf5_output_period_ps = hdf5_output_period_ps
        self.hdf5_writer = None  # Will be initialized during setup
        
        # Molecular temperature decomposition parameters
        self.enable_molecular_temps = enable_molecular_temps
        self.molecular_temps_output_period_ps = molecular_temps_output_period_ps
        
        # Dipole moment FDR parameters
        self.enable_dipole_fdr = enable_dipole_fdr
        self.dipole_fdr_output_period_ps = dipole_fdr_output_period_ps
        self.dipole_fdr_max_correlation_time_ps = dipole_fdr_max_correlation_time_ps
        self.dipole_fdr_field_direction = dipole_fdr_field_direction
        self.dipole_fdr_exclude_cavity = dipole_fdr_exclude_cavity
        self.enable_dipole_response = enable_dipole_response
        self.dipole_response_field_strength = dipole_response_field_strength
        self.dipole_response_sign = dipole_response_sign
        
        # Mechanical cavity modulation parameters
        self.mech_periodic = mech_periodic
        self.mech_frequency_cm1 = mech_frequency_cm1
        self.mech_magnitude = mech_magnitude
        
        # Auto-stopping parameters
        self.enable_auto_stopping = enable_auto_stopping
        self.auto_stop_min_time_ps = auto_stop_min_time_ps
        self.auto_stop_smoothing_window = auto_stop_smoothing_window
        
        
        
        # Logging parameters - simplified to console only
        self.log_level = log_level
        self.custom_log_file = custom_log_file
        
        # F(k,t) calculation parameters
        self.enable_fkt = enable_fkt
        self.fkt_kmag = fkt_kmag
        self.fkt_num_wavevectors = fkt_num_wavevectors
        self.fkt_reference_interval_ps = fkt_reference_interval_ps
        self.fkt_max_references = fkt_max_references
        
        # F(k,t) logarithmic time spacing parameters
        self.fkt_log_time_spacing = fkt_log_time_spacing
        self.fkt_min_log_time_ps = fkt_min_log_time_ps
        self.fkt_max_log_time_ps = fkt_max_log_time_ps
        self.fkt_log_num_points = fkt_log_num_points
        
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
        self.enable_gsd_output = enable_gsd_output
        self.console_output_period_ps = console_output_period_ps
        
        # Text output parameters (deprecated but kept for compatibility)
        self.enable_text_output = enable_text_output
        self.text_output_file = text_output_file
        
        # GSD file control
        self.truncate_gsd = truncate_gsd
        
        # Initialize simulation components (will be set during setup)
        self.sim = None
        self.logger = None
        
        # Delayed controller management
        self._delayed_controllers = []  # List of controllers to add during simulation
        self._delayed_controller_monitor = None
        
    def _add_delayed_controller(self, controller_setup_func, activation_time_ps, controller_name):
        """Register a controller to be added during simulation at a specific time."""
        self._delayed_controllers.append({
            'setup_func': controller_setup_func,
            'activation_time_ps': activation_time_ps,
            'controller_name': controller_name,
            'activated': False
        })
        self.log_info(f"Delayed controller registered: {controller_name} (activates at {activation_time_ps:.1f} ps)")
    
    def _create_delayed_controller_monitor(self):
        """Create a monitor that activates delayed controllers during simulation."""
        if not self._delayed_controllers:
            return None
            
        class DelayedControllerMonitor(hoomd.custom.Action):
            def __init__(self, simulation_obj):
                self.simulation = simulation_obj
                self.time_tracker = simulation_obj.time_tracker
                
            def act(self, timestep):
                current_time_ps = self.time_tracker.elapsed_time
                
                for controller_data in self.simulation._delayed_controllers:
                    if not controller_data['activated'] and current_time_ps >= controller_data['activation_time_ps']:
                        try:
                            # Activate the controller
                            controller_data['setup_func']()
                            controller_data['activated'] = True
                            self.simulation.log_info(f"Delayed controller activated: {controller_data['controller_name']} at t = {current_time_ps:.2f} ps")
                        except Exception as e:
                            self.simulation.log_error(f"Failed to activate delayed controller {controller_data['controller_name']}: {e}")
        
        return DelayedControllerMonitor(self)

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
            if self.use_deprecated_couplstr:
                self.log_info(f"  DEPRECATED: Coupling strength (couplstr): {self.couplstr} a.u.")
            else:
                self.log_info(f"  Lambda coupling: {self.lambda_coupling}")
                self.log_info(f"  Epsilon (lambda * omega_c): {self.lambda_coupling * self.omegac:.6e} a.u.")
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
    def _create_coupling_variant(self):
        """Create coupling variant based on coupling_variant_type."""
        if not hasattr(self, 'time_tracker') or self.time_tracker is None:
            raise ValueError("Time tracker must be set up before creating coupling variants")

        # Determine the base coupling value to use
        if self.use_deprecated_couplstr:
            # Using deprecated couplstr - use it directly as epsilon
            base_coupling = self.couplstr
            self.log_info(f"DEPRECATED: Using couplstr={base_coupling} a.u. directly (not scaled by omega_c)")
            self.log_info(f"    Please migrate to lambda_coupling parameter for future compatibility")
        else:
            # Using new lambda_coupling - will be scaled by omega_c
            base_coupling = self.lambda_coupling
            self.log_info(f"Using lambda_coupling={base_coupling} (will be scaled by omega_c={self.omegac:.6f} a.u.)")

        variant_type = self.coupling_variant_type.lower()

        if variant_type == 'constant':
            # Constant coupling (default)
            from hoomd.variant import Constant

            if self.use_deprecated_couplstr:
                # For backward compatibility, use couplstr directly
                coupling_variant = Constant(base_coupling)
                self.log_info(f"Using constant coupling: {base_coupling} a.u. (deprecated mode)")
            else:
                # Pass lambda (dimensionless) directly - C++ will multiply by omegac
                coupling_variant = Constant(base_coupling)
                epsilon = base_coupling * self.omegac
                self.log_info(f"Using constant coupling: lambda={base_coupling}, epsilon={epsilon:.6e} a.u.")
        
        elif variant_type == 'step':
            # Step coupling (with optional turn-off and decay)
            step_variant = StepVariant(
                target_value=base_coupling,
                switch_time_ps=self.switch_time_ps if self.switch_time_ps is not None else 0.0,
                time_tracker=self.time_tracker,
                decay_time_constant_ps=self.decay_time_constant_ps,
                turn_off_time_ps=getattr(self, 'step_turn_off_time_ps', None)
            )

            if self.use_deprecated_couplstr:
                coupling_variant = step_variant
                self.log_info(f"Using step coupling (deprecated mode):")
                self.log_info(f"  Target value: {base_coupling} a.u.")
            else:
                # Pass lambda (dimensionless) directly - C++ will multiply by omegac
                coupling_variant = step_variant
                epsilon = base_coupling * self.omegac
                self.log_info(f"Using step coupling:")
                self.log_info(f"  Target lambda: {base_coupling}")
                self.log_info(f"  Target epsilon: {epsilon:.6e} a.u.")

            self.log_info(f"  Switch time: {self.switch_time_ps} ps")
            if self.decay_time_constant_ps:
                self.log_info(f"  Decay time constant: {self.decay_time_constant_ps} ps")
        
        elif variant_type == 'periodic':
            # Periodic coupling
            periodic_variant = PeriodicVariant(
                amplitude=base_coupling,
                time_tracker=self.time_tracker,
                period_ps=getattr(self, 'periodic_period_ps', 1.0),
                phase_offset=getattr(self, 'periodic_phase_offset', 0.0),
                start_time_ps=getattr(self, 'periodic_start_time_ps', 0.0),
                stop_time_ps=getattr(self, 'periodic_stop_time_ps', None)
            )

            if self.use_deprecated_couplstr:
                coupling_variant = periodic_variant
                self.log_info(f"Using periodic coupling (deprecated mode):")
                self.log_info(f"  Amplitude: {base_coupling} a.u.")
            else:
                # Pass lambda (dimensionless) directly - C++ will multiply by omegac
                coupling_variant = periodic_variant
                epsilon = base_coupling * self.omegac
                self.log_info(f"Using periodic coupling:")
                self.log_info(f"  Amplitude lambda: {base_coupling}")
                self.log_info(f"  Amplitude epsilon: {epsilon:.6e} a.u.")

            self.log_info(f"  Period: {getattr(self, 'periodic_period_ps', 1.0)} ps")
            self.log_info(f"  Phase offset: {getattr(self, 'periodic_phase_offset', 0.0):.3f} rad")
        
        elif variant_type == 'exponential':
            # Exponential decay coupling
            # For exponential, use exponential_amplitude if provided, else use base_coupling
            exp_amplitude = self.exponential_amplitude if self.exponential_amplitude is not None else base_coupling

            exp_variant = ExponentialDecayVariant(
                amplitude=exp_amplitude,
                time_tracker=self.time_tracker,
                decay_time_constant_ps=self.exponential_decay_time_ps,
                turn_on_time_ps=self.exponential_turn_on_time_ps,
                turn_off_time_ps=self.exponential_turn_off_time_ps
            )

            if self.use_deprecated_couplstr:
                coupling_variant = exp_variant
                self.log_info(f"Using exponential decay coupling (deprecated mode):")
                self.log_info(f"  Amplitude: {exp_amplitude} a.u.")
            else:
                # Pass lambda (dimensionless) directly - C++ will multiply by omegac
                coupling_variant = exp_variant
                epsilon = exp_amplitude * self.omegac
                self.log_info(f"Using exponential decay coupling:")
                self.log_info(f"  Amplitude lambda: {exp_amplitude}")
                self.log_info(f"  Amplitude epsilon: {epsilon:.6e} a.u.")

            self.log_info(f"  Decay time constant: {self.exponential_decay_time_ps} ps")
            self.log_info(f"  Turn-on time: {self.exponential_turn_on_time_ps} ps")
        
        elif variant_type == 'square':
            # Square wave coupling
            square_variant = SquareWaveVariant(
                amplitude=base_coupling,
                period_ps=self.square_period_ps,
                time_tracker=self.time_tracker,
                duty_cycle=self.square_duty_cycle,
                phase_offset=self.square_phase_offset,
                start_time_ps=self.square_start_time_ps,
                stop_time_ps=self.square_stop_time_ps
            )

            if self.use_deprecated_couplstr:
                coupling_variant = square_variant
                self.log_info(f"Using square wave coupling (deprecated mode):")
                self.log_info(f"  Amplitude: {base_coupling} a.u.")
            else:
                # Pass lambda (dimensionless) directly - C++ will multiply by omegac
                coupling_variant = square_variant
                epsilon = base_coupling * self.omegac
                self.log_info(f"Using square wave coupling:")
                self.log_info(f"  Amplitude lambda: {base_coupling}")
                self.log_info(f"  Amplitude epsilon: {epsilon:.6e} a.u.")

            self.log_info(f"  Period: {self.square_period_ps} ps")
            self.log_info(f"  Duty cycle: {self.square_duty_cycle:.1%}")
        
        elif variant_type == 'decaying_square':
            # Decaying square wave coupling
            decay_square_variant = DecayingSquareWaveVariant(
                initial_amplitude=base_coupling,
                period_ps=self.decaying_square_period_ps,
                time_tracker=self.time_tracker,
                decay_rate_per_period=self.decaying_square_decay_rate,
                duty_cycle=self.decaying_square_duty_cycle,
                phase_offset=self.decaying_square_phase_offset,
                start_time_ps=self.decaying_square_start_time_ps,
                stop_time_ps=self.decaying_square_stop_time_ps,
                minimum_amplitude=self.decaying_square_minimum_amplitude
            )

            if self.use_deprecated_couplstr:
                coupling_variant = decay_square_variant
                self.log_info(f"Using decaying square wave coupling (deprecated mode):")
                self.log_info(f"  Initial amplitude: {base_coupling} a.u.")
            else:
                # Pass lambda (dimensionless) directly - C++ will multiply by omegac
                coupling_variant = decay_square_variant
                epsilon = base_coupling * self.omegac
                self.log_info(f"Using decaying square wave coupling:")
                self.log_info(f"  Initial amplitude lambda: {base_coupling}")
                self.log_info(f"  Initial amplitude epsilon: {epsilon:.6e} a.u.")

            self.log_info(f"  Period: {self.decaying_square_period_ps} ps")
            self.log_info(f"  Duty cycle: {self.decaying_square_duty_cycle:.1%}")
            self.log_info(f"  Decay rate: {self.decaying_square_decay_rate:.1%} per period")
            self.log_info(f"  Minimum amplitude: {self.decaying_square_minimum_amplitude:.2e} a.u.")
        
        elif variant_type == 'adaptive_square':
            # Adaptive square wave coupling - needs target temperature and temperature tracker
            
            # Determine target temperature from controllers
            target_temperature = None
            
            # Debug: Check what controllers are enabled
            print(f"DEBUG: Checking temperature controllers for adaptive square wave:")
            print(f"  enable_quench_controller: {getattr(self, 'enable_quench_controller', 'NOT_SET')}")
            print(f"  enable_gd_feedback: {getattr(self, 'enable_gd_feedback', 'NOT_SET')}")
            print(f"  enable_dual_feedback: {getattr(self, 'enable_dual_feedback', 'NOT_SET')}")
            print(f"  enable_diffeq_controller: {getattr(self, 'enable_diffeq_controller', 'NOT_SET')}")
            
            if self.enable_quench_controller:
                target_temperature = self.quench_target_temperature
                temp_source = "quench controller"
            elif self.enable_gd_feedback:
                target_temperature = self.gd_target_temperature
                temp_source = "gradient descent controller"
            elif (hasattr(self, 'enable_dual_feedback') and self.enable_dual_feedback):
                # Use cavity target temperature from dual controller
                target_temperature = getattr(self, 'dual_cavity_target_temperature', 100.0)
                temp_source = "dual independent controller (cavity target)"
                print(f"DEBUG: Using dual controller, target_temperature = {target_temperature}")
            elif (hasattr(self, 'enable_diffeq_controller') and self.enable_diffeq_controller):
                # Use initial temperature from DiffEq controller (T_initial in the formula)
                target_temperature = self.temperature  # This is the --temperature parameter (T_initial)
                temp_source = "differential equation controller (initial temperature)"
                print(f"DEBUG: Using DiffEq controller, target_temperature = {target_temperature}")
            else:
                raise ValueError("Adaptive square wave coupling requires either --enable-quench-controller, "
                    "--enable-gd-feedback, --enable-dual-feedback, or --enable-diffeq-controller to provide target temperature")
            
            # Ensure temperature tracker is enabled
            if not self.enable_temp_tracker:
                raise ValueError("Adaptive square wave coupling requires --enable-temp-tracker "
                    "to measure harmonic equipartition temperature")
            
            # Create adaptive square wave variant (temperature tracker will be set later)
            coupling_variant = AdaptiveSquareWaveVariant(
                target_coupling=self.couplstr,
                target_temperature=target_temperature,
                period_ps=self.adaptive_square_period_ps,
                time_tracker=self.time_tracker,
                temperature_tracker=None,  # Will be set after temperature tracker is created
                duty_cycle=self.adaptive_square_duty_cycle,
                phase_offset=self.adaptive_square_phase_offset,
                start_time_ps=self.adaptive_square_start_time_ps,
                stop_time_ps=self.adaptive_square_stop_time_ps,
                min_amplitude=self.adaptive_square_min_amplitude,
                max_amplitude=self.adaptive_square_max_amplitude,
                simulation=self  # Pass simulation reference for auto-stop signal
            )
            
            # Store reference for later temperature tracker assignment
            self._adaptive_coupling_variant = coupling_variant
            
            self.log_info(f"Using adaptive square wave coupling:")
            self.log_info(f"  Target coupling (g_target): {self.couplstr} a.u.")
            self.log_info(f"  Target temperature (T_target): {target_temperature:.1f} K (from {temp_source})")
            self.log_info(f"  Period: {self.adaptive_square_period_ps} ps")
            self.log_info(f"  Duty cycle: {self.adaptive_square_duty_cycle:.1%}")
            self.log_info(f"  Amplitude limits: [{self.adaptive_square_min_amplitude:.2e}, {self.adaptive_square_max_amplitude:.2e}] a.u.")
            self.log_info(f"  Algorithm: g_next = g_target × √(T_target / T_harmonic)")
        
        elif variant_type == 'composite':
            # Create composite coupling from sinusoidal + adaptive square wave
            
            # First, create sinusoidal component
            sinusoid_variant = PeriodicVariant(
                amplitude=self.composite_sinusoid_amplitude,
                period_ps=self.composite_sinusoid_period,
                time_tracker=self.time_tracker,
                phase_offset=self.composite_sinusoid_phase,
                start_time_ps=self.composite_sinusoid_start_time,
                stop_time_ps=self.composite_sinusoid_stop_time
            )
            
            # Create square wave component (adaptive or fixed)
            if self.composite_square_adaptive:
                # Determine target temperature for adaptive square component
                if hasattr(self, 'enable_quench_controller') and self.enable_quench_controller:
                    target_temperature = getattr(self, 'quench_target_temperature', 100.0)
                    temp_source = "quench controller"
                elif hasattr(self, 'enable_gd_feedback') and self.enable_gd_feedback:
                    target_temperature = getattr(self, 'gd_target_temperature', 100.0)
                    temp_source = "gradient descent controller"
                elif hasattr(self, 'enable_dual_feedback') and self.enable_dual_feedback:
                    target_temperature = getattr(self, 'dual_cavity_target_temperature', 100.0)
                    temp_source = "dual independent controller (cavity target)"
            else:
                raise ValueError("Adaptive square wave in composite coupling requires a controller (quench, gd, or dual feedback) to provide target temperature")
            
            # Check if temperature tracker is available yet
            temp_tracker = getattr(self, 'temperature_tracker', None)
            if temp_tracker is None:
                # Temperature tracker not set up yet, fall back to fixed square wave
                self.log_warning("Temperature tracker not available yet for adaptive square wave, using fixed amplitude instead")
                square_wave_variant = SquareWaveVariant(
                    amplitude=self.composite_square_amplitude,
                    period_ps=self.composite_square_period,
                    time_tracker=self.time_tracker,
                    duty_cycle=self.composite_square_duty_cycle,
                    phase_offset=0.0,
                    start_time_ps=self.composite_square_start_time,
                    stop_time_ps=self.composite_square_stop_time
                )
            else:
                square_wave_variant = AdaptiveSquareWaveVariant(
                    target_amplitude=self.composite_square_amplitude,
                    period_ps=self.composite_square_period,
                    duty_cycle=self.composite_square_duty_cycle,
                    time_tracker=self.time_tracker,
                    temperature_tracker=temp_tracker,
                    target_temperature=target_temperature,
                    start_time_ps=self.composite_square_start_time,
                    stop_time_ps=self.composite_square_stop_time,
                    target_coupling=0.0  # This should be set to the desired constant value post-auto-stop
                )
            
            # Create composite variant
            coupling_variant = CompositeVariant(
                variants=[sinusoid_variant, square_wave_variant],
                max_amplitude=self.composite_max_amplitude
            )
            
            # Store reference for later temperature tracker assignment
            self._composite_coupling_variant = coupling_variant
            
            self.log_info(f"Using composite coupling:")
            self.log_info(f"  Sinusoidal component:")
            self.log_info(f"    Amplitude: {self.composite_sinusoid_amplitude:.2e} a.u.")
            self.log_info(f"    Period: {self.composite_sinusoid_period} ps")
            self.log_info(f"    Phase: {self.composite_sinusoid_phase:.2f} rad")
            self.log_info(f"    Start time: {self.composite_sinusoid_start_time} ps")
            self.log_info(f"    Stop time: {self.composite_sinusoid_stop_time} ps")
            
            if self.composite_square_adaptive:
                if getattr(self, 'temperature_tracker', None) is not None:
                    self.log_info(f"  Adaptive square wave component:")
                    self.log_info(f"    Target amplitude: {self.composite_square_amplitude:.2e} a.u.")
                    self.log_info(f"    Target temperature: {target_temperature:.1f} K (from {temp_source})")
                    self.log_info(f"    Period: {self.composite_square_period} ps")
                    self.log_info(f"    Duty cycle: {self.composite_square_duty_cycle:.1%}")
                    self.log_info(f"    Start time: {self.composite_square_start_time} ps")
                    self.log_info(f"    Stop time: {self.composite_square_stop_time} ps")
                else:
                    self.log_info(f"  Fixed square wave component (adaptive requested but temperature tracker not available):")
                    self.log_info(f"    Amplitude: {self.composite_square_amplitude:.2e} a.u.")
                    self.log_info(f"    Period: {self.composite_square_period} ps")
                    self.log_info(f"    Duty cycle: {self.composite_square_duty_cycle:.1%}")
                    self.log_info(f"    Start time: {self.composite_square_start_time} ps")
                    self.log_info(f"    Stop time: {self.composite_square_stop_time} ps")
            else:
                self.log_info(f"  Fixed square wave component:")
                self.log_info(f"    Amplitude: {self.composite_square_amplitude:.2e} a.u.")
                self.log_info(f"    Period: {self.composite_square_period} ps")
                self.log_info(f"    Duty cycle: {self.composite_square_duty_cycle:.1%}")
                self.log_info(f"    Start time: {self.composite_square_start_time} ps")
                self.log_info(f"    Stop time: {self.composite_square_stop_time} ps")
            if self.composite_max_amplitude is not None:
                self.log_info(f"  Maximum total amplitude: {self.composite_max_amplitude:.2e} a.u.")
            else:
                self.log_info(f"  Maximum total amplitude: No limit")
        
        elif variant_type == 'exponentialwave':
            # Exponential wave coupling
            from ..variants import ExponentialWaveVariant
            
            # Check if adaptive mode is requested
            if self.exp_adaptive:
                # Adaptive exponential wave requires temperature tracker
                temp_tracker = getattr(self, 'temperature_tracker', None)
                if temp_tracker is None:
                    self.log_warning("Temperature tracker not available for adaptive exponential wave, using fixed amplitude instead")
                    adaptive_mode = False
                else:
                    adaptive_mode = True
            else:
                adaptive_mode = False
            
            coupling_variant = ExponentialWaveVariant(
                amplitude=self.couplstr,
                period_ps=self.exp_period_ps,
                tau_ps=self.exp_tau_ps,
                time_tracker=self.time_tracker,
                start_time_ps=self.exp_start_time_ps,
                stop_time_ps=self.exp_stop_time_ps,
                adaptive=adaptive_mode,
                temperature_tracker=getattr(self, 'temperature_tracker', None),
                simulation=self  # Pass simulation reference for auto-stop signal
            )
            
            self.log_info(f"Using exponential wave coupling:")
            self.log_info(f"  Target amplitude: {self.couplstr} a.u.")
            self.log_info(f"  Period: {self.exp_period_ps} ps")
            self.log_info(f"  Decay time constant: {self.exp_tau_ps} ps")
            self.log_info(f"  Start time: {self.exp_start_time_ps} ps")
            if self.exp_stop_time_ps is not None:
                self.log_info(f"  Stop time: {self.exp_stop_time_ps} ps")
            else:
                self.log_info(f"  No stop time (runs indefinitely)")
            if adaptive_mode:
                self.log_info(f"  Adaptive mode: ENABLED (T_bath/T_harmonic scaling)")
                self.log_info(f"  Algorithm: coupling_new = coupling_target × √(T_bath / T_harmonic)")
            else:
                self.log_info(f"  Adaptive mode: DISABLED (fixed amplitude)")

        elif variant_type == 'lqg_coupling':
            # LQG coupling controller has been removed from codebase
            self.log_error("LQG coupling controller has been removed. Please use a different coupling variant.")
            raise ValueError("LQG coupling controller has been removed from codebase.")
        else:
            raise ValueError(f"Unknown coupling_variant_type: {variant_type}")
        
        return coupling_variant
    
    def _create_laser_force(self):
        """Create enhanced laser force with timing control."""
        from .fdr_forces import PerturbationForce, PerturbationTimingUpdater
        
        # Create laser force with timing
        laser_force = PerturbationForce(
            kvector=self.laser_kvector,
            amplitude=self.laser_amplitude,
            start_time_ps=self.laser_start_time_ps,
            stop_time_ps=self.laser_stop_time_ps,
            time_tracker=self.time_tracker
        )
        
        # Store for timing updates (will be added to simulation operations)
        self.laser_force = laser_force
        
        self.log_info(f"Enhanced laser drive configured:")
        self.log_info(f"  k-vector: {self.laser_kvector}")
        self.log_info(f"  Amplitude: {self.laser_amplitude}")
        self.log_info(f"  Start time: {self.laser_start_time_ps} ps")
        if self.laser_stop_time_ps is not None:
            self.log_info(f"  Stop time: {self.laser_stop_time_ps} ps")
        else:
            self.log_info(f"  Stop time: None (runs indefinitely)")
        
        return laser_force

    def setup_force_parameters(self, dt, rcut=15):
        """Set up force parameters for the simulation."""
        forces = []
        
        # Initialize cavity force reference
        self.cavity_force = None
        
        # Setup cavity force if requested
        if self.incavity:
            from ..utils import PhysicalConstants
            omegac = self.freq / PhysicalConstants.HARTREE_TO_CM_MINUS1
            
            # Create coupling variant based on type selection
            coupling_variant = self._create_coupling_variant()
            
            # Create cavity force with variant
            from ..forces import CavityForce
            cavityforce = CavityForce(
                kvector=np.array([0,0,1]), 
                lambda_coupling=coupling_variant, 
                omegac=omegac,
            )
            
            # Store variants and cavity force for later use in finite-q setup and console output
            self.coupling_variant = coupling_variant
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

        # Always disable pair interaction with 'L' particle (photon) if it exists
        # This ensures compatibility whether incavity=True or False
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
        
        # Setup enhanced laser driving if enabled
        if self.laser_enabled:
            laser_force = self._create_laser_force()
            forces.append(laser_force)
        
        return forces

    def setup_thermostat_parameters(self, dt):
        """Set up thermostat parameters for molecular and cavity systems."""
        from ..utils import PhysicalConstants
        # Stock HOOMD Bussi: BussiReservoir currently leaves GPU simulations in
        # an "Invalid data location state" on HOOMD 5.4 + CUDA 13 before the
        # first sim.run (reproduced with and without CavityForce).
        Bussi = hoomd.md.methods.thermostats.Bussi
        
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
            self.log_info(
                "Running molecular system with stock HOOMD Bussi thermostat "
                "(NVT ensemble; not BussiReservoir)"
            )
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
                self.log_info("Running cavity with stock HOOMD Bussi thermostat")
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
        # Defer host-snapshot T-check until after the first sim.run: reading a
        # full Snapshot before CavityForceComputeGPU initializes has been
        # observed to leave an Invalid data location state on GPU.
        self.log_info("Thermalization completed - velocities properly initialized")
    def _thermalize_cavity_particle_manually(self, kT):
        """Thermalize the cavity particle with HOOMD's State API.

        Manual velocity writes through ``gpu_local_snapshot`` are unreliable
        when CuPy is unavailable (``np.asarray`` does not materialize
        ``HOOMDGPUArray`` buffers). ``thermalize_particle_momenta`` keeps
        data-location tracking consistent before the first ``sim.run``.
        """
        particle_types = list(self.sim.state.particle_types)
        if "L" not in particle_types:
            self.log_info("WARNING: No cavity particle type 'L' for thermalization!")
            return
        cavity_filter = hoomd.filter.Type(["L"])
        self.sim.state.thermalize_particle_momenta(kT=kT, filter=cavity_filter)
        self.log_info("Cavity particle thermalized via HOOMD State API:")
        self.log_info(f"  Target temperature: {self.temperature:.1f} K")
        self.log_info(f"  kT: {kT:.6e} a.u.")

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
            
            # Enable dynamic coupling detection by default, but still support legacy switch_time
            dynamic_coupling_detection = getattr(self, 'enable_dynamic_coupling_detection', True)
            coupling_change_threshold = getattr(self, 'coupling_change_threshold', 1e-5)
            
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
                shock_dampening_time_constant_ps=self.time_constant_ps,  # Use the same time constant for shock dampening
                # New dynamic coupling detection parameters
                dynamic_coupling_detection=dynamic_coupling_detection,
                coupling_change_threshold=coupling_change_threshold
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
        # logger[('Timestep', 'dt_fs')] = (self.timestep_formatter, 'dt_fs', 'scalar')  # Temporarily disabled for testing
        
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
                    simulation=self,  # ← Pass CavityMDSimulation object which has incavity attribute
                    time_tracker=self.time_tracker,
                    output_period_ps=energy_output_period_ps,
                    output_prefix=output_prefix,
                    force_objects=force_objects,        # CRITICAL: Pass force objects
                    thermostat_objects=thermostat_objects,  # CRITICAL: Pass thermostat objects
                    verbose="quiet",  # Suppress debug output by default
                    enable_csv_output=False  # Use HDF5 output via ObservableWriter
                )
                
                # Add energy tracker to simulation - trigger period doesn't matter since it uses internal timing
                energy_updater = hoomd.update.CustomUpdater(
                    action=self.energy_tracker,
                    trigger=hoomd.trigger.Periodic(1)  # Check every step but tracker handles timing internally
                )
                self.sim.operations.updaters.append(energy_updater)
                
                self.log_info(f" Energy tracker setup completed with time-based output:")
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
                
                # Get logarithmic time spacing parameters
                fkt_log_time_spacing = getattr(self, 'fkt_log_time_spacing', False)
                fkt_min_log_time_ps = getattr(self, 'fkt_min_log_time_ps', None)
                fkt_max_log_time_ps = getattr(self, 'fkt_max_log_time_ps', None)
                fkt_log_num_points = getattr(self, 'fkt_log_num_points', 50)
                
                # Log logarithmic time spacing configuration
                if fkt_log_time_spacing:
                    self.log_info(f"  Logarithmic time spacing: ENABLED")
                    self.log_info(f"  Log time range: {fkt_min_log_time_ps} - {fkt_max_log_time_ps} ps")
                    self.log_info(f"  Log time points: {fkt_log_num_points}")
                else:
                    self.log_info(f"  Logarithmic time spacing: DISABLED (regular spacing)")
                
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
                    num_wavevectors=fkt_num_wavevectors,
                    # Logarithmic time spacing parameters
                    log_time_spacing=fkt_log_time_spacing,
                    min_log_time_ps=fkt_min_log_time_ps,
                    max_log_time_ps=fkt_max_log_time_ps,
                    log_num_points=fkt_log_num_points
                )
                
                # Add F(k,t) tracker to simulation - trigger period doesn't matter since it uses internal timing
                fkt_updater = hoomd.update.CustomUpdater(
                    action=self.density_corr_tracker,
                    trigger=hoomd.trigger.Periodic(1)  # Check every step but tracker handles timing internally
                )
                self.sim.operations.updaters.append(fkt_updater)
                
                # Add F(k,t) data to logger
                logger[('F(k,t)', 'current_autocorr')] = (self.density_corr_tracker, 'current_autocorr', 'scalar')
                
                self.log_info(" F(k,t) tracker successfully enabled with time-based output:")
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
                
                self.log_info(" Dipole autocorrelation tracker successfully enabled with time-based output:")
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
                
                self.log_info(" Momentum zeroing successfully enabled:")
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
        
        self.log_info(" TRACKING AND LOGGING SETUP COMPLETED:")
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
        
        # Set up empirical temperature feedback if enabled
        if getattr(self, 'enable_empirical_feedback', False):
            self._setup_empirical_feedback()
            enabled_features.append(f"empirical temperature feedback ({self.feedback_update_interval_ps:.1f} ps)")
        
        # Set up enhanced Phase 3 features
        self._setup_enhanced_features(enabled_features)
        
        if enabled_features:
            self.log_info(f"Advanced features enabled: {', '.join(enabled_features)}")
        else:
            self.log_info("Running with basic tracking only")
        
        console_output_period_ps = getattr(self, 'console_output_period_ps', 1.0)
        self.log_info(f"Console output: {console_output_period_ps:.3f} ps periods (time-based)")
        self.log_info("  Works accurately for both adaptive and fixed timestep modes")

    def _setup_empirical_feedback(self):
        """Set up empirical temperature feedback system."""
        try:
            from ..controllers.empirical import EmpiricalTemperatureData
            from ..controllers.feedback import EmpiricalTemperatureFeedback
            
            # Validate parameters
            if not self.empirical_data_file:
                raise ValueError("empirical_data_file must be provided when enable_empirical_feedback=True")
            
            # Check that energy tracking is enabled
            if not self.enable_energy_tracking:
                self.log_warning("Energy tracking is not enabled - empirical feedback may not work correctly")
                self.log_warning("Consider setting enable_energy_tracking=True for proper empirical feedback")
            
            self.log_info("Setting up empirical temperature feedback system:")
            self.log_info(f"  Data file: {self.empirical_data_file}")
            self.log_info(f"  Energy component: {self.feedback_energy_component}")
            self.log_info(f"  Apply to: {self.feedback_apply_to}")
            self.log_info(f"  Averaging window: {self.feedback_averaging_window_ps:.1f} ps")
            self.log_info(f"  Update interval: {self.feedback_update_interval_ps:.1f} ps")
            self.log_info(f"  Temperature range: [{self.feedback_T_min:.1f}, {self.feedback_T_max:.1f}] K")
            if self.switch_time_ps is not None:
                self.log_info(f"  Switch time: {self.switch_time_ps:.1f} ps (feedback starts after equilibration)")
            else:
                self.log_info("  Switch time: None (feedback active from start)")
            if self.feedback_turn_off_time_ps is not None:
                self.log_info(f"  Turn-off time: {self.feedback_turn_off_time_ps:.1f} ps (feedback stops here)")
            if self.feedback_energy_component == 'harmonic':
                if self.feedback_use_direct_harmonic_calculation:
                    self.log_info("  Harmonic method: Direct calculation T = 4*E/(N*kB)")
                else:
                    self.log_info("  Harmonic method: Empirical interpolation")
            
            # Load empirical data
            self.empirical_data = EmpiricalTemperatureData(
                data_file_path=self.empirical_data_file,
                energy_component=self.feedback_energy_component
            )
            
            # Create output file path
            output_prefix = f"empirical_feedback_replica_{self.replica}"
            
            # Create feedback tracker
            self.empirical_feedback = EmpiricalTemperatureFeedback(
                empirical_data=self.empirical_data,
                energy_tracker=self.energy_tracker,
                molecular_thermostat=self.molecular_thermostat_obj,
                cavity_thermostat=self.cavity_thermostat_obj,
                apply_to=self.feedback_apply_to,
                output_period_ps=self.feedback_output_period_ps,
                averaging_window_ps=self.feedback_averaging_window_ps,
                update_interval_ps=self.feedback_update_interval_ps,
                T_min=self.feedback_T_min,
                T_max=self.feedback_T_max,
                turn_off_time_ps=self.feedback_turn_off_time_ps,
                switch_time_ps=self.switch_time_ps,
                output_file=f"{output_prefix}_feedback.csv",
                time_tracker=self.time_tracker,  # Add time tracker for accurate timing
                initial_temperature=self.temperature  # Pass initial bath temperature
            )
            
            # Calculate trigger period for CSV output (convert ps to steps)
            # Use conservative estimate: assume ~1000 steps per ps for trigger calculation
            feedback_trigger_steps = max(1, int(self.feedback_output_period_ps * 1000))
            
            # Add to simulation with configurable output frequency
            feedback_updater = hoomd.update.CustomUpdater(
                action=self.empirical_feedback,
                trigger=hoomd.trigger.Periodic(feedback_trigger_steps)
            )
            self.sim.operations.updaters.append(feedback_updater)
            
            self.log_info(" Empirical temperature feedback system enabled")
            
            # Print basic empirical data info
            self.log_info(f"  Energy component: {self.empirical_data.energy_component}")
            self.log_info(f"  Temperature range: {self.empirical_data.temperatures.min():.1f} - {self.empirical_data.temperatures.max():.1f} K")
            self.log_info(f"  Energy range: {self.empirical_data.energies.min():.6f} - {self.empirical_data.energies.max():.6f} Hartree")
            self.log_info(f"  Data points: {len(self.empirical_data.temperatures)}")
            
        except Exception as e:
            self.log_error(f"Failed to setup empirical temperature feedback: {e}")
            import traceback
            self.log_error("Full traceback:")
            for line in traceback.format_exc().split('\n'):
                if line.strip():
                    self.log_error(line)
            self.empirical_feedback = None
    
    def _setup_enhanced_features(self, enabled_features):
        """Set up enhanced Phase 3 features: laser timing and PI feedback."""
        
        # Set up laser timing updater if laser is enabled
        if getattr(self, 'laser_enabled', False) and hasattr(self, 'laser_force'):
            from .fdr_forces import PerturbationTimingUpdater
            
            laser_timing_updater = PerturbationTimingUpdater([self.laser_force])
            laser_updater = hoomd.update.CustomUpdater(
                action=laser_timing_updater,
                trigger=hoomd.trigger.Periodic(1)  # Check every timestep
            )
            self.sim.operations.updaters.append(laser_updater)
            enabled_features.append(f"laser drive timing ({self.laser_start_time_ps:.1f}-{self.laser_stop_time_ps or '∞'} ps)")
        
        # Set up PI feedback controller if enabled
        if getattr(self, 'enable_pi_feedback', False):
            self._setup_pi_feedback()
            enabled_features.append(f"PI feedback controller ({self.pi_temperature_method})")
        
        # Set up gradient descent feedback controller if enabled
        # NOTE: GD and quench controllers are mutually exclusive for temperature control
        if getattr(self, 'enable_gd_feedback', False) and not getattr(self, 'enable_quench_controller', False):
            self._setup_gd_feedback()
            enabled_features.append(f"gradient descent feedback controller ({self.gd_temperature_method}, τ={self.gd_time_constant_ps:.1f}ps)")
        elif getattr(self, 'enable_gd_feedback', False) and getattr(self, 'enable_quench_controller', False):
            self.log_info("  WARNING: Both GD feedback and quench controller enabled - disabling GD to avoid conflicts")
            self.gd_feedback = None
        
        # Set up dual independent feedback controller if enabled
        # NOTE: Dual controller can be immediate or delayed
        if getattr(self, 'enable_dual_feedback', False):
            # Check if this should be a delayed controller (turn-on time > current time)
            dual_turn_on_time = getattr(self, 'dual_turn_on_time_ps', 0.0)
            current_time = self.time_tracker.elapsed_time if self.time_tracker else 0.0
            
            if dual_turn_on_time > current_time + 1.0:  # Allow 1ps buffer for immediate activation
                # Register as delayed controller
                self._add_delayed_controller(
                    controller_setup_func=self._setup_dual_feedback,
                    activation_time_ps=dual_turn_on_time,
                    controller_name="DualIndependentTemperatureFeedback"
                )
                enabled_features.append(f"dual independent feedback controller (DELAYED to {dual_turn_on_time:.1f}ps)")
            else:
                # Check for conflicts with other immediate controllers
                if (not getattr(self, 'enable_gd_feedback', False) and 
                    not getattr(self, 'enable_quench_controller', False) and
                    not getattr(self, 'enable_diffeq_controller', False)):
                    # Setup immediately
                    self._setup_dual_feedback()
                    enabled_features.append(f"dual independent feedback controller (cavity:{self.dual_cavity_method}, molecular:{self.dual_molecular_method})")
                else:
                    self.log_info("  WARNING: Dual feedback enabled with other immediate controllers - setting up as delayed controller")
                    # Force it to be delayed to avoid conflicts
                    self._add_delayed_controller(
                        controller_setup_func=self._setup_dual_feedback,
                        activation_time_ps=max(dual_turn_on_time, current_time + 5.0),  # At least 5ps delay
                        controller_name="DualIndependentTemperatureFeedback"
                    )
                    enabled_features.append(f"dual independent feedback controller (DELAYED due to conflicts)")
        else:
            self.dual_feedback = None
        
        # Set up differential equation controller if enabled
        # NOTE: DiffEq controller can coexist with delayed controllers but not immediate ones
        if getattr(self, 'enable_diffeq_controller', False):
            # Check for immediate conflicts (controllers that start right away)
            immediate_conflicts = []
            if getattr(self, 'enable_gd_feedback', False):
                immediate_conflicts.append("gradient descent")
            if getattr(self, 'enable_quench_controller', False):
                immediate_conflicts.append("quench")
            
            # Check if dual controller is immediate (not delayed)
            dual_turn_on_time = getattr(self, 'dual_turn_on_time_ps', 0.0)
            current_time = self.time_tracker.elapsed_time if self.time_tracker else 0.0
            if (getattr(self, 'enable_dual_feedback', False) and 
                dual_turn_on_time <= current_time + 1.0):  # Immediate activation
                immediate_conflicts.append(f"dual independent (immediate, turn_on={dual_turn_on_time:.1f}ps)")
            
            if not immediate_conflicts:
                self._setup_diffeq_controller()
                if self.diffeq_time_constant_auto:
                    enabled_features.append(f"differential equation controller ({self.diffeq_temperature_method}, adaptive τ)")
                else:
                    enabled_features.append(f"differential equation controller ({self.diffeq_temperature_method}, τ={self.diffeq_time_constant_ps:.1f}ps)")
            else:
                self.log_info(f"  WARNING: DiffEq controller conflicts with immediate controllers: {', '.join(immediate_conflicts)} - disabling diffeq")
                self.diffeq_controller = None
        
        # Set up comprehensive temperature tracker if enabled
        # IMPORTANT: Must be set up BEFORE BathPI and LQR controllers (which need temperature measurements)
        if getattr(self, 'enable_temp_tracker', False):
            self._setup_temperature_tracker()
            enabled_features.append(f"comprehensive temperature tracker ({self.temp_tracker_output_period_ps:.1f} ps)")
        
        # Set up BathPI controller if enabled (AFTER temperature tracker!)
        if getattr(self, 'enable_bath_pi_controller', False):
            # Check for conflicts with immediate controllers
            immediate_conflicts = []
            if getattr(self, 'enable_empirical_feedback', False):
                immediate_conflicts.append("empirical feedback")
            if getattr(self, 'enable_lqr_controller', False):
                immediate_conflicts.append("LQR")
            if getattr(self, 'enable_quench_controller', False):
                immediate_conflicts.append("quench")
            
            # Check if dual controller is immediate (not delayed)
            dual_turn_on_time = getattr(self, 'dual_turn_on_time_ps', 0.0)
            current_time = self.time_tracker.elapsed_time if self.time_tracker else 0.0
            if (getattr(self, 'enable_dual_feedback', False) and 
                dual_turn_on_time <= current_time + 1.0):
                immediate_conflicts.append(f"dual independent (immediate, turn_on={dual_turn_on_time:.1f}ps)")
            
            if not immediate_conflicts:
                self.log_info(f"  WARNING: BathPI controller has been removed from codebase")
                self.bath_pi_controller = None
            else:
                self.log_info(f"  WARNING: BathPI controller has been removed from codebase")
                self.bath_pi_controller = None
        
        # Set up SimpleSetpointController if enabled
        if getattr(self, 'enable_simple_setpoint_controller', False):
            self._setup_simple_setpoint_controller()
            if self.simple_setpoint_controller is not None:
                enabled_features.append(f"SimpleSetpointController ({self.simple_setpoint_apply_to}, signal={self.simple_setpoint_signal_method})")
        
        # Set up Adaptive MPC Controller if enabled
        if getattr(self, 'enable_adaptive_mpc_controller', False):
            self._setup_adaptive_mpc_controller()
            if self.adaptive_mpc_controller is not None:
                enabled_features.append(f"Adaptive MPC Controller (target={getattr(self, 'adaptive_mpc_target_temperature', 100.0):.1f}K, Np={getattr(self, 'adaptive_mpc_prediction_horizon', 10)})")
        
        # Note: PID Controller setup moved to run() method after thermostat creation
        # to ensure proper thermostat object references
        if getattr(self, 'enable_pid_controller', False):
            mode = "self-loop" if getattr(self, 'pid_self_loop', False) else f"setpoint={getattr(self, 'pid_target_temperature', 100.0):.1f}K"
            enabled_features.append(f"PID Controller ({mode}, signal={getattr(self, 'pid_signal_choice', 'lj_coulombic')})")
        
        # Set up harmonic bond reset if enabled
        if getattr(self, 'enable_harmonic_reset', False):
            self._setup_harmonic_reset()
            enabled_features.append(f"harmonic bond reset (t={self.harmonic_reset_turn_on_time_ps:.1f}ps)")
        
        # Set up HDF5 observable output if enabled
        # NOTE: Must be set up AFTER trackers (energy, temperature) so they can be registered
        if getattr(self, 'enable_hdf5_output', False):
            self._setup_hdf5_output()
            enabled_features.append(f"HDF5 observable output ({self.hdf5_output_period_ps:.3f} ps)")
        
        # Note: Quench controller setup moved to after thermostat creation
        # to ensure proper thermostat object references
        if getattr(self, 'enable_quench_controller', False):
            enabled_features.append(f"quench controller ({self.quench_initial_temperature:.1f}K→{self.quench_target_temperature:.1f}K at {self.quench_time_ps:.1f}ps)")
        
        # Set up molecular temperature decomposition if enabled
        # NOTE: Molecular temperatures are now tracked by TemperatureTracker with track_molecular=True
        if getattr(self, 'enable_molecular_temps', False):
            # self._setup_molecular_temperature_tracker()  # DEPRECATED - now handled by TemperatureTracker
            enabled_features.append(f"molecular temperature decomposition (integrated in TemperatureTracker)")
            
        if getattr(self, 'enable_temp_tracker', False):
            # Set up auto-stop controller if enabled (requires temperature tracker)
            if getattr(self, 'enable_auto_stop', False):
                self._setup_auto_stop_controller()
                # Add auto-stop controller to operations if successfully created
                if hasattr(self, 'auto_stop_controller') and self.auto_stop_controller is not None:
                    auto_stop_updater = hoomd.update.CustomUpdater(
                        action=self.auto_stop_controller,
                        trigger=hoomd.trigger.Periodic(10)  # Check every 10 steps
                    )
                    self.sim.operations.updaters.append(auto_stop_updater)
                    enabled_features.append(f"auto-stop coupling (tol={self.auto_stop_tol:.1f}K, window={self.auto_stop_window:.1f}ps)")
                    self.log_info("Auto-stop controller added to simulation operations")
        
        # Check for auto-stop without temperature tracker
        elif getattr(self, 'enable_auto_stop', False):
            self.log_warning("Auto-stop controller requires --enable-temp-tracker to be enabled")
            self.log_warning("Auto-stop feature will be disabled for this simulation")
        
        # Set up sinusoidal bath temperature controller if enabled
        if getattr(self, 'enable_sinusoidal_bath', False):
            self._setup_sinusoidal_bath_controller()
            enabled_features.append(f"sinusoidal bath controller (period={self.sinusoidal_bath_period_ps:.1f}ps, scale={self.sinusoidal_bath_amplitude_scale:.2f})")
        
        # Set up adaptive bath temperature controller if enabled
        if getattr(self, 'enable_adaptive_bath', False):
            self._setup_adaptive_bath_controller()
            enabled_features.append(f"adaptive bath controller (scale={self.adaptive_bath_amplitude_scale:.2f}, method={self.adaptive_bath_signal_temperature_method})")
        
        # Set up dipole moment FDR tracking if enabled
        if getattr(self, 'enable_dipole_fdr', False):
            self._setup_dipole_fdr()
            enabled_features.append(f"dipole moment FDR tracker ({self.dipole_fdr_output_period_ps:.1f} ps, τ_max={self.dipole_fdr_max_correlation_time_ps:.1f} ps)")
        
        # Set up dipole response force if enabled
        if getattr(self, 'enable_dipole_response', False):
            self._setup_dipole_response()
            enabled_features.append(f"dipole response force (E₀={self.dipole_response_field_strength:.2e}, sign={self.dipole_response_sign:+.0f})")
    def _setup_gd_feedback(self):
        """Set up gradient descent feedback temperature controller."""
        if not self.enable_gd_feedback:
            self.gd_feedback = None
            return
        
        try:
            from ..controllers.feedback import GradientDescentTemperatureFeedback
            
            # Create output file path
            output_file = f"gd_feedback_replica_{self.replica}.csv"
            
            # Create gradient descent controller
            self.gd_feedback = GradientDescentTemperatureFeedback(
                temperature_method=self.gd_temperature_method,
                time_constant_ps=self.gd_time_constant_ps,
                time_tracker=self.time_tracker,
                energy_tracker=getattr(self, 'energy_tracker', None),
                simulation=self.sim,
                molecular_thermostat=getattr(self, 'molecular_thermostat_obj', None),
                cavity_thermostat=getattr(self, 'cavity_thermostat_obj', None),
                target_temperature=self.gd_target_temperature,
                apply_to=self.gd_apply_to,
                turn_on_time_ps=self.gd_turn_on_time_ps,
                turn_off_time_ps=self.gd_turn_off_time_ps,
                update_interval_ps=self.gd_update_interval_ps,
                T_min=self.gd_T_min,
                T_max=self.gd_T_max,
                output_file=output_file,
                empirical_data_file=getattr(self, 'temp_tracker_empirical_data_file', None),
                console_output_period_ps=self.console_output_period_ps
            )
            
            # Add to simulation 
            gd_trigger_steps = max(1, int(self.gd_update_interval_ps * 1000))
            gd_updater = hoomd.update.CustomUpdater(
                action=self.gd_feedback,
                trigger=hoomd.trigger.Periodic(gd_trigger_steps)
            )
            self.sim.operations.updaters.append(gd_updater)
            
            self.log_info(" Gradient descent feedback controller enabled")
            self.log_info(f"  Target temperature: {self.gd_target_temperature:.1f} K")
            self.log_info(f"  Method: {self.gd_temperature_method}")
            self.log_info(f"  Time constant: {self.gd_time_constant_ps:.1f} ps")
            self.log_info(f"  Learning rate α: {0.0001/self.gd_time_constant_ps:.6f}")
            self.log_info(f"  Update interval: {self.gd_update_interval_ps:.1f} ps")
            
        except Exception as e:
            self.log_error(f"Failed to setup gradient descent feedback controller: {e}")
            self.gd_feedback = None
    
    def _setup_dual_feedback(self):
        """Set up dual independent temperature feedback controller."""
        if not self.enable_dual_feedback:
            self.dual_feedback = None
            return
        
        try:
            from ..controllers.dual_feedback import DualIndependentTemperatureFeedback
            
            # Create output file path
            output_file = f"dual_feedback_replica_{self.replica}.csv"
            
            # Create dual independent controller
            self.dual_feedback = DualIndependentTemperatureFeedback(
                cavity_method=self.dual_cavity_method,
                molecular_method=self.dual_molecular_method,
                cavity_time_constant_ps=self.dual_cavity_time_constant_ps,
                molecular_time_constant_ps=self.dual_molecular_time_constant_ps,
                time_tracker=self.time_tracker,
                energy_tracker=getattr(self, 'energy_tracker', None),
                simulation=self.sim,
                molecular_thermostat=getattr(self, 'molecular_thermostat_obj', None),
                cavity_thermostat=getattr(self, 'cavity_thermostat_obj', None),
                cavity_target_temperature=self.dual_cavity_target_temperature,
                molecular_target_temperature=self.dual_molecular_target_temperature,
                turn_on_time_ps=self.dual_turn_on_time_ps,
                turn_off_time_ps=self.dual_turn_off_time_ps,
                update_interval_ps=self.dual_update_interval_ps,
                cavity_T_min=self.dual_cavity_T_min,
                cavity_T_max=self.dual_cavity_T_max,
                molecular_T_min=self.dual_molecular_T_min,
                molecular_T_max=self.dual_molecular_T_max,
                cavity_dynamic_target=self.dual_cavity_dynamic_target,
                molecular_dynamic_target=self.dual_molecular_dynamic_target,
                cavity_integral_time_constant_ps=getattr(self, 'dual_cavity_integral_time_constant_ps', None),
                molecular_integral_time_constant_ps=getattr(self, 'dual_molecular_integral_time_constant_ps', None),
                cavity_heating_gain_factor=getattr(self, 'dual_cavity_heating_gain_factor', 1.0),
                cavity_cooling_gain_factor=getattr(self, 'dual_cavity_cooling_gain_factor', 1.0),
                molecular_heating_gain_factor=getattr(self, 'dual_molecular_heating_gain_factor', 1.0),
                molecular_cooling_gain_factor=getattr(self, 'dual_molecular_cooling_gain_factor', 1.0),
                cavity_integral_heating_time_constant_ps=getattr(self, 'dual_cavity_integral_heating_time_constant_ps', None),
                cavity_integral_cooling_time_constant_ps=getattr(self, 'dual_cavity_integral_cooling_time_constant_ps', None),
                molecular_integral_heating_time_constant_ps=getattr(self, 'dual_molecular_integral_heating_time_constant_ps', None),
                molecular_integral_cooling_time_constant_ps=getattr(self, 'dual_molecular_integral_cooling_time_constant_ps', None),
                output_file=output_file,
                empirical_data_file=getattr(self, 'temp_tracker_empirical_data_file', None),
                console_output_period_ps=self.console_output_period_ps
            )
            
            # Add to simulation 
            dual_trigger_steps = max(1, int(self.dual_update_interval_ps * 1000))
            dual_updater = hoomd.update.CustomUpdater(
                action=self.dual_feedback,
                trigger=hoomd.trigger.Periodic(dual_trigger_steps)
            )
            self.sim.operations.updaters.append(dual_updater)
            
            self.log_info(" Dual independent feedback controller enabled")
            self.log_info(f"  Cavity: {self.dual_cavity_method} → {self.dual_cavity_target_temperature:.1f}K (τ={self.dual_cavity_time_constant_ps:.1f}ps)")
            self.log_info(f"  Molecular: {self.dual_molecular_method} → {self.dual_molecular_target_temperature:.1f}K (τ={self.dual_molecular_time_constant_ps:.1f}ps)")
            self.log_info(f"  Update interval: {self.dual_update_interval_ps:.1f} ps")
            
        except Exception as e:
            self.log_error(f"Failed to setup dual independent feedback controller: {e}")
            self.dual_feedback = None
    
    def _setup_diffeq_controller(self):
        """Set up differential equation temperature controller."""
        if not self.enable_diffeq_controller:
            self.diffeq_controller = None
            return
        
        try:
            from ..controllers import DiffEqController
            
            # Create output file path
            output_file = f"diffeq_control_replica_{self.replica}.csv"
            
            # Create differential equation controller with full Kalman filter + auto-tuning support
            self.diffeq_controller = DiffEqController(
                temperature_method=self.diffeq_temperature_method,
                time_constant_ps=self.diffeq_time_constant_ps,
                time_tracker=self.time_tracker,
                energy_tracker=getattr(self, 'energy_tracker', None),
                simulation=self.sim,
                time_constant_auto=self.diffeq_time_constant_auto,
                molecular_thermostat=getattr(self, 'molecular_thermostat_obj', None),
                cavity_thermostat=getattr(self, 'cavity_thermostat_obj', None),
                apply_to=self.diffeq_apply_to,
                turn_on_time_ps=self.diffeq_turn_on_time_ps,
                update_interval_ps=self.diffeq_update_interval_ps,
                T_min=self.diffeq_T_min,
                T_max=self.diffeq_T_max,
                output_file=output_file,
                empirical_data_file=getattr(self, 'temp_tracker_empirical_data_file', None),
                # Bias correction parameters (using actual DiffEqController interface)
                enable_bias_estimation=not self.diffeq_disable_bias_estimation,  # Controlled by --diffeq-disable-bias-estimation flag
                console_output_period_ps=self.console_output_period_ps,
                # Exact-cancellation + PI control parameters
                enable_pi_control=self.diffeq_enable_pi_control,
                pi_rho=self.diffeq_pi_rho,
                pi_epsilon=self.diffeq_pi_epsilon,
                pi_zeta=self.diffeq_pi_zeta,
                relaxation_data_file=self.diffeq_relaxation_data_file,
                filter_window_ps=self.diffeq_filter_window,
                # Adaptive bias cancellation parameters
                enable_bias_cancellation=self.diffeq_enable_bias_cancellation,
                bias_tau_b_ps=self.diffeq_bias_tau_b_ps,
                bias_tau_b_auto=self.diffeq_bias_tau_b_auto,
                bias_kappa=self.diffeq_bias_kappa,
                bias_kappa_auto=self.diffeq_bias_kappa_auto,
                bias_tau_b_prefactor=self.diffeq_bias_tau_b_prefactor,
                bias_kappa_prefactor=self.diffeq_bias_kappa_prefactor,
                bias_calibration_time_ps=self.diffeq_bias_calibration_time_ps,
                enable_csv_output=False  # Use HDF5 output via ObservableWriter
            )
            
            # Add to simulation with appropriate trigger frequency
            # Use update interval to determine trigger frequency
            dt_ps = PhysicalConstants.atomic_units_to_ps(self.sim.operations.integrator.dt)
            steps_per_ps = 1.0 / dt_ps
            trigger_steps = max(1, int(self.diffeq_update_interval_ps * steps_per_ps))
            
            diffeq_updater = hoomd.update.CustomUpdater(
                action=self.diffeq_controller,
                trigger=hoomd.trigger.Periodic(trigger_steps)
            )
            self.sim.operations.updaters.append(diffeq_updater)
            
            self.log_info(f" Differential equation controller enabled")
            self.log_info(f"  Method: {self.diffeq_temperature_method}")
            if self.diffeq_time_constant_auto:
                self.log_info(f"  Time constant: ADAPTIVE τ = T[T_bath] (fallback: {self.diffeq_time_constant_ps:.2f} ps)")
            else:
                self.log_info(f"  Time constant: {self.diffeq_time_constant_ps:.2f} ps (fixed)")
            self.log_info(f"  Turn on time: {self.diffeq_turn_on_time_ps:.1f} ps")
            self.log_info(f"  Update interval: {self.diffeq_update_interval_ps:.3f} ps")
            self.log_info(f"  Apply to: {self.diffeq_apply_to}")
            if self.diffeq_T_min is not None:
                self.log_info(f"  Minimum temperature: {self.diffeq_T_min:.1f} K")
            if self.diffeq_T_max is not None:
                self.log_info(f"  Maximum temperature: {self.diffeq_T_max:.1f} K")
            if self.diffeq_turn_off_time_ps is not None:
                self.log_info(f"  Turn off time: {self.diffeq_turn_off_time_ps:.1f} ps")
            if self.diffeq_rate_limit_K_per_ps is not None:
                self.log_info(f"  Rate limit: {self.diffeq_rate_limit_K_per_ps:.2f} K/ps")
            
            # Log Kalman filter and auto-tuning status
            enable_bias_estimation_actual = not self.diffeq_disable_bias_estimation
            self.log_info(f"  Kalman filter bias estimation: {'ENABLED' if enable_bias_estimation_actual else 'DISABLED'}")
            if enable_bias_estimation_actual:
                self.log_info(f"    Process noise: {self.diffeq_bias_process_noise:.2e} K²/ps")
                self.log_info(f"    Initial covariance: {self.diffeq_bias_initial_covariance:.1f} K²")
            
            self.log_info(f"  Auto-tuning: {'ENABLED' if self.diffeq_enable_autotuning else 'DISABLED'}")
            if self.diffeq_enable_autotuning:
                self.log_info(f"    Duration: {self.diffeq_autotune_duration_ps:.1f} ps")
                self.log_info(f"    PRBS perturbation: ±{self.diffeq_autotune_perturbation_amplitude_K:.1f} K")
                self.log_info(f"    Measurement noise: {self.diffeq_autotune_measurement_noise_K:.1f} K")
                self.log_info(f"    Bias timescale factor: {self.diffeq_autotune_min_bias_timescale_factor:.1f}×")
            
            if self.diffeq_periodic_coupling_period_ps is not None:
                self.log_info(f"  Periodic coupling period: {self.diffeq_periodic_coupling_period_ps:.2f} ps (Floquet/Lifted KF mode)")
            
        except Exception as e:
            import traceback
            self.log_error(f"Failed to setup differential equation controller: {e}")
            self.log_error(f"Full traceback: {traceback.format_exc()}")
            self.diffeq_controller = None
    
    def _setup_simple_setpoint_controller(self):
        """Set up SimpleSetpointController for temperature control."""
        if not self.enable_simple_setpoint_controller:
            self.simple_setpoint_controller = None
            return
        
        try:
            # Create output file path
            output_file = f"{self.simple_setpoint_output_file.replace('.csv', '')}_replica_{self.replica}.csv"
            
            # Get empirical data file path if needed
            empirical_data_file = None
            if self.simple_setpoint_signal_method in ['lj_coulombic', 'harmonic']:
                empirical_data_file = getattr(self, 'temp_tracker_empirical_data_file', None)
            
            # Create SimpleSetpointController
            self.simple_setpoint_controller = SimpleSetpointController(
                signal_method=self.simple_setpoint_signal_method,
                time_constant_ps=self.simple_setpoint_time_constant_ps,
                time_tracker=self.time_tracker,
                energy_tracker=getattr(self, 'energy_tracker', None),
                molecular_thermostat=getattr(self, 'molecular_thermostat_obj', None),
                cavity_thermostat=getattr(self, 'cavity_thermostat_obj', None),
                apply_to=self.simple_setpoint_apply_to,
                turn_on_time_ps=self.simple_setpoint_turn_on_time_ps,
                turn_off_time_ps=self.simple_setpoint_turn_off_time_ps,
                update_interval_ps=self.simple_setpoint_update_interval_ps,
                T_min=self.simple_setpoint_T_min,
                T_max=self.simple_setpoint_T_max,
                output_file=output_file,
                empirical_data_file=empirical_data_file,
                console_output_period_ps=self.simple_setpoint_console_output_period_ps,
                enable_csv_output=False  # Use HDF5 output via ObservableWriter
            )
            
            # Add to simulation operations
            steps_per_ps = 1.0 / (self.dt_fs * 1e-3)  # Convert fs to ps
            trigger_steps = max(1, int(self.simple_setpoint_update_interval_ps * steps_per_ps))
            
            simple_setpoint_updater = hoomd.update.CustomUpdater(
                action=self.simple_setpoint_controller,
                trigger=hoomd.trigger.Periodic(trigger_steps)
            )
            self.sim.operations.updaters.append(simple_setpoint_updater)
            
            self.log_info(f"\n SimpleSetpointController enabled")
            self.log_info(f"   Signal method: {self.simple_setpoint_signal_method}")
            self.log_info(f"   Time constant: {self.simple_setpoint_time_constant_ps:.2f} ps")
            self.log_info(f"   Apply to: {self.simple_setpoint_apply_to}")
            self.log_info(f"   Turn on time: {self.simple_setpoint_turn_on_time_ps:.2f} ps")
            if self.simple_setpoint_turn_off_time_ps is not None:
                self.log_info(f"   Turn off time: {self.simple_setpoint_turn_off_time_ps:.2f} ps")
            self.log_info(f"   Output file: {output_file}")
            if empirical_data_file is not None:
                self.log_info(f"   Empirical data: {empirical_data_file}")
            
        except Exception as e:
            import traceback
            self.log_error(f"Failed to setup SimpleSetpointController: {e}")
            self.log_error(f"Full traceback: {traceback.format_exc()}")
            self.simple_setpoint_controller = None
    
    def _setup_adaptive_mpc_controller(self):
        """Set up Adaptive MPC Controller for temperature regulation."""
        if not getattr(self, 'enable_adaptive_mpc_controller', False):
            self.adaptive_mpc_controller = None
            return
        
        try:
            from ..controllers.adaptive_mpc import AdaptiveMPCController
            
            # Create output file path
            output_file = f"{getattr(self, 'adaptive_mpc_output_file', 'adaptive_mpc_control.csv')}"
            output_file = output_file.replace('.csv', f'_replica_{self.replica}.csv')
            
            # Create Adaptive MPC Controller
            self.adaptive_mpc_controller = AdaptiveMPCController(
                simulation=self,
                time_tracker=self.time_tracker,
                temperature_tracker=self.temperature_tracker,
                target_temperature=getattr(self, 'adaptive_mpc_target_temperature', 100.0),
                turn_on_time_ps=getattr(self, 'adaptive_mpc_turn_on_time_ps', 0.0),
                turn_off_time_ps=getattr(self, 'adaptive_mpc_turn_off_time_ps', None),
                system_id_duration_ps=getattr(self, 'adaptive_mpc_system_id_duration_ps', 50.0),
                update_interval_ps=getattr(self, 'adaptive_mpc_update_interval_ps', 0.1),
                prediction_horizon=getattr(self, 'adaptive_mpc_prediction_horizon', 10),
                control_horizon=getattr(self, 'adaptive_mpc_control_horizon', 5),
                output_weight=getattr(self, 'adaptive_mpc_output_weight', 100.0),
                control_effort_weight=getattr(self, 'adaptive_mpc_control_effort_weight', [1.0, 0.1]),
                rate_penalty_weight=getattr(self, 'adaptive_mpc_rate_penalty_weight', [10.0, 1.0]),
                lambda_min=getattr(self, 'adaptive_mpc_lambda_min', 0.0),
                lambda_max=getattr(self, 'adaptive_mpc_lambda_max', 1e-2),
                T_bath_min=getattr(self, 'adaptive_mpc_T_bath_min', 0.1),
                T_bath_max=getattr(self, 'adaptive_mpc_T_bath_max', 500.0),
                delta_lambda_max=getattr(self, 'adaptive_mpc_delta_lambda_max', 1e-4),
                delta_T_bath_max=getattr(self, 'adaptive_mpc_delta_T_bath_max', 10.0),
                apply_to=getattr(self, 'adaptive_mpc_apply_to', 'both'),
                system_id_step_duration_ps=getattr(self, 'adaptive_mpc_system_id_step_duration_ps', 5.0),
                system_id_seed=getattr(self, 'adaptive_mpc_system_id_seed', 42),
                rls_forgetting_factor=getattr(self, 'adaptive_mpc_rls_forgetting_factor', 0.995),
                rls_initial_covariance=getattr(self, 'adaptive_mpc_rls_initial_covariance', 100.0),
                model_update_interval=getattr(self, 'adaptive_mpc_model_update_interval', 10),
                empirical_data_file=getattr(self, 'temp_tracker_empirical_data_file', None),
                output_file=output_file,
                console_output_period_ps=getattr(self, 'adaptive_mpc_console_output_period_ps', 1.0),
                regularization_param=getattr(self, 'adaptive_mpc_regularization_param', 1e-3),
                use_scaling=getattr(self, 'adaptive_mpc_use_scaling', True),
                debug_mode=getattr(self, 'adaptive_mpc_debug_mode', False)
            )
            
            # Add to simulation operations
            steps_per_ps = 1.0 / (self.dt_fs * 1e-3)  # Convert fs to ps
            trigger_steps = max(1, int(getattr(self, 'adaptive_mpc_update_interval_ps', 0.1) * steps_per_ps))
            
            adaptive_mpc_updater = hoomd.update.CustomUpdater(
                action=self.adaptive_mpc_controller,
                trigger=hoomd.trigger.Periodic(trigger_steps)
            )
            self.sim.operations.updaters.append(adaptive_mpc_updater)
            
            self.log_info(f"\n Adaptive MPC Controller enabled")
            self.log_info(f"   Target temperature: {getattr(self, 'adaptive_mpc_target_temperature', 100.0):.2f} K")
            self.log_info(f"   Turn on time: {getattr(self, 'adaptive_mpc_turn_on_time_ps', 0.0):.2f} ps")
            self.log_info(f"   System ID duration: {getattr(self, 'adaptive_mpc_system_id_duration_ps', 50.0):.2f} ps")
            self.log_info(f"   Prediction horizon: {getattr(self, 'adaptive_mpc_prediction_horizon', 10)}")
            self.log_info(f"   Control horizon: {getattr(self, 'adaptive_mpc_control_horizon', 5)}")
            self.log_info(f"   Lambda range: [{getattr(self, 'adaptive_mpc_lambda_min', 0.0):.2e}, {getattr(self, 'adaptive_mpc_lambda_max', 1e-2):.2e}]")
            self.log_info(f"   T_bath range: [{getattr(self, 'adaptive_mpc_T_bath_min', 0.1):.2f}, {getattr(self, 'adaptive_mpc_T_bath_max', 500.0):.2f}] K")
            self.log_info(f"   Apply to: {getattr(self, 'adaptive_mpc_apply_to', 'both')}")
            self.log_info(f"   Output file: {output_file}")
            
        except Exception as e:
            import traceback
            self.log_error(f"Failed to setup Adaptive MPC Controller: {e}")
            self.log_error(f"Full traceback: {traceback.format_exc()}")
            self.adaptive_mpc_controller = None
    
    def _setup_pid_controller(self):
        """Set up PID Controller for temperature regulation."""
        if not getattr(self, 'enable_pid_controller', False):
            self.pid_controller = None
            return
        
        try:
            from ..controllers.pid_control import PIDControl
            
            # Verify prerequisites
            if not hasattr(self, 'temperature_tracker') or self.temperature_tracker is None:
                raise RuntimeError("FATAL: PID controller requires temperature_tracker to be initialized first")
            if not hasattr(self, 'time_tracker') or self.time_tracker is None:
                raise RuntimeError("FATAL: PID controller requires time_tracker to be initialized first")
            if not hasattr(self, 'molecular_thermostat_obj') or self.molecular_thermostat_obj is None:
                raise RuntimeError("FATAL: PID controller requires molecular_thermostat_obj to be initialized first")
            if not hasattr(self, 'cavity_thermostat_obj') or self.cavity_thermostat_obj is None:
                raise RuntimeError("FATAL: PID controller requires cavity_thermostat_obj to be initialized first")
            
            # Create output file path
            output_file = f"{getattr(self, 'pid_output_file', 'pid_control.csv')}"
            output_file = output_file.replace('.csv', f'_replica_{self.replica}.csv')
            
            # Determine if manual gains are provided
            Kp = getattr(self, 'pid_Kp', None)
            Ti = getattr(self, 'pid_Ti', None)
            Td = getattr(self, 'pid_Td', None)
            
            # Create PID Controller
            self.pid_controller = PIDControl(
                temperature_tracker=self.temperature_tracker,
                time_tracker=self.time_tracker,
                simulation=self,
                molecular_thermostat=self.molecular_thermostat_obj,
                cavity_thermostat=self.cavity_thermostat_obj,
                signal_choice=getattr(self, 'pid_signal_choice', 'lj_coulombic'),
                target_temperature=getattr(self, 'pid_target_temperature', 100.0),
                self_loop=getattr(self, 'pid_self_loop', False),
                Kp=Kp,
                Ti=Ti,
                Td=Td,
                auto_tune=getattr(self, 'pid_auto_tune', True),
                auto_tune_step_size=getattr(self, 'pid_auto_tune_step_size', 20.0),
                auto_tune_duration_ps=getattr(self, 'pid_auto_tune_duration_ps', 50.0),
                turn_on_time_ps=getattr(self, 'pid_turn_on_time_ps', 0.0),
                turn_off_time_ps=getattr(self, 'pid_turn_off_time_ps', None),
                update_interval_ps=getattr(self, 'pid_update_interval_ps', 0.1),
                apply_to=getattr(self, 'pid_apply_to', 'both'),
                T_min=getattr(self, 'pid_T_min', 0.1),
                T_max=getattr(self, 'pid_T_max', None),
                rate_limit_K_per_ps=getattr(self, 'pid_rate_limit_K_per_ps', None),
                derivative_filter_N=getattr(self, 'pid_derivative_filter_N', 10.0),
                enable_anti_windup=getattr(self, 'pid_enable_anti_windup', True),
                output_file=output_file,
                console_output_period_ps=getattr(self, 'pid_console_output_period_ps', 1.0),
                empirical_data_file=getattr(self, 'temp_tracker_empirical_data_file', None)
            )
            
            # Add to simulation operations
            steps_per_ps = 1.0 / (self.dt_fs * 1e-3)  # Convert fs to ps
            trigger_steps = max(1, int(getattr(self, 'pid_update_interval_ps', 0.1) * steps_per_ps))
            
            pid_updater = hoomd.update.CustomUpdater(
                action=self.pid_controller,
                trigger=hoomd.trigger.Periodic(trigger_steps)
            )
            self.sim.operations.updaters.append(pid_updater)
            
            self.log_info(f"\n PID Controller enabled")
            self.log_info(f"   Signal: {getattr(self, 'pid_signal_choice', 'lj_coulombic')}")
            if getattr(self, 'pid_self_loop', False):
                self.log_info(f"   Mode: Self-loop (T_setpoint = T_bath)")
            else:
                self.log_info(f"   Mode: Setpoint tracking (T* = {getattr(self, 'pid_target_temperature', 100.0):.2f} K)")
            self.log_info(f"   Turn on time: {getattr(self, 'pid_turn_on_time_ps', 0.0):.2f} ps")
            if Kp is not None or Ti is not None or Td is not None:
                self.log_info(f"   Tuning: Manual (Kp={Kp}, Ti={Ti}, Td={Td})")
            else:
                self.log_info(f"   Tuning: Auto (step={getattr(self, 'pid_auto_tune_step_size', 20.0):.1f} K, duration={getattr(self, 'pid_auto_tune_duration_ps', 50.0):.1f} ps)")
            self.log_info(f"   Apply to: {getattr(self, 'pid_apply_to', 'both')}")
            self.log_info(f"   Output file: {output_file}")
            
        except Exception as e:
            import traceback
            self.log_error(f"FATAL: Failed to setup PID Controller: {e}")
            self.log_error(f"Full traceback: {traceback.format_exc()}")
            raise RuntimeError(f"FATAL: PID controller setup failed: {e}") from e
    
    def _setup_harmonic_reset(self):
        """Set up harmonic bond reset."""
        if not self.enable_harmonic_reset:
            return
        
        try:
            from ..updaters import HarmonicBondReset
            
            # Use dynamic bath temperature if reset_temperature is None
            use_dynamic_temperature = (self.harmonic_reset_temperature is None)
            initial_temp = self.harmonic_reset_temperature if self.harmonic_reset_temperature is not None else self.temperature
            
            self.log_info(f"\n{'='*60}")
            self.log_info("HARMONIC BOND RESET CONFIGURATION")
            self.log_info(f"{'='*60}")
            self.log_info(f"Bond parameters: Will auto-detect from simulation at reset time")
            if use_dynamic_temperature:
                self.log_info(f"Reset temperature: DYNAMIC (will use bath temperature at reset time)")
            else:
                self.log_info(f"Reset temperature: {initial_temp} K (fixed)")
            self.log_info(f"Turn-on time: {self.harmonic_reset_turn_on_time_ps} ps")
            self.log_info(f"Random seed: {self.harmonic_reset_seed}")
            self.log_info(f"{'='*60}\n")
            
            # Create a placeholder bond reset action
            # Parameters will be populated at reset time
            bond_reset = HarmonicBondReset(
                bond_params={'placeholder': {'K': 1.0, 'r0': 1.0}},  # Temporary
                temperature=initial_temp,  # Initial value
                kB=PhysicalConstants.KB_HARTREE_PER_K,  # Correct atomic units!
                seed=self.harmonic_reset_seed
            )
            bond_reset._sim_obj = self  # Store reference to simulation object
            
            # Load empirical data for harmonic fictive temperature calculation
            # Use the same empirical data file as the temperature tracker
            if hasattr(self, 'empirical_data_file') and self.empirical_data_file is not None:
                from pathlib import Path
                from ..controllers.empirical import EmpiricalTemperatureData
                empirical_path = Path(self.empirical_data_file)
                if empirical_path.exists():
                    try:
                        harmonic_empirical_data = EmpiricalTemperatureData(
                            str(empirical_path),
                            energy_component='harmonic',
                            use_direct_harmonic=False,  # Use extended fits
                            create_plots=False
                        )
                        bond_reset.set_empirical_data(harmonic_empirical_data)
                        print(f"✓ Loaded harmonic empirical data for bond reset from {empirical_path}")
                    except Exception as e:
                        print(f"Warning: Could not load empirical data for bond reset: {e}")
            
            # Add as custom updater that checks periodically
            reset_updater = hoomd.update.CustomUpdater(
                action=bond_reset,
                trigger=hoomd.trigger.Periodic(period=100)  # Check every 100 steps
            )
            self.sim.operations.updaters.append(reset_updater)
            
            # Create a trigger callback to enable the reset at the specified time
            class HarmonicResetTrigger(hoomd.custom.Action):
                def __init__(self, bond_reset_action, trigger_time_ps, use_dynamic_temp, sim_obj):
                    self.bond_reset_action = bond_reset_action
                    self.trigger_time_ps = trigger_time_ps
                    self.triggered = False
                    self.use_dynamic_temp = use_dynamic_temp
                    self.sim_obj = sim_obj  # Reference to CavityMDSimulation object
                
                def act(self, timestep):
                    # Use actual simulation time instead of timestep count (for adaptive timestep compatibility)
                    current_time_ps = self.sim_obj.time_tracker.elapsed_time
                    
                    if not self.triggered and current_time_ps >= self.trigger_time_ps:
                        # Update temperature dynamically from bath if requested
                        if self.use_dynamic_temp:
                            try:
                                # Get actual kinetic temperature from particle velocities, not thermostat setpoint!
                                # CRITICAL: Must exclude cavity particle (type 'X') from molecular bath calculation
                                
                                # Get type names from state (not from snapshot, which doesn't have it)
                                type_names = self.sim_obj.sim.state.particle_types
                                cavity_type_idx = None
                                if 'X' in type_names:
                                    cavity_type_idx = type_names.index('X')
                                
                                with self.sim_obj.sim.state.cpu_local_snapshot as snap:
                                    velocities = snap.particles.velocity
                                    masses = snap.particles.mass
                                    types = snap.particles.typeid
                                    
                                    # Exclude cavity particle from kinetic temperature calculation
                                    if cavity_type_idx is not None:
                                        molecular_mask = types != cavity_type_idx
                                        mol_velocities = velocities[molecular_mask]
                                        mol_masses = masses[molecular_mask]
                                    else:
                                        mol_velocities = velocities
                                        mol_masses = masses
                                    
                                    # Calculate kinetic energy: KE = 0.5 * sum(m * v^2)
                                    KE_total = 0.5 * np.sum(mol_masses[:, np.newaxis] * mol_velocities**2)
                                    
                                    # Number of degrees of freedom (3N for N molecular particles, excluding cavity)
                                    N_mol_dof = 3 * len(mol_masses)
                                    
                                    # Equipartition: KE = 0.5 * N_dof * kB * T
                                    # T = 2 * KE / (N_dof * kB)
                                    molecular_bath_T = 2.0 * KE_total / (N_mol_dof * PhysicalConstants.KB_HARTREE_PER_K)
                                    
                                    self.bond_reset_action.T = molecular_bath_T
                                    print(f"\n>>> Harmonic bond reset: Using dynamic kinetic temperature T = {molecular_bath_T:.2f} K <<<")
                                    print(f"    (Calculated from {len(mol_masses)} molecular particles, excluding cavity)")
                                    print(f"    (Actual instantaneous kinetic energy, not thermostat setpoint)")
                                    
                            except Exception as e:
                                    print(f"\n>>> Warning: Could not calculate kinetic temperature ({e}), using default T = {float(self.bond_reset_action.T):.2f} K <<<")
                                    import traceback
                                    traceback.print_exc()
                        
                        print(f"\n>>> Enabling harmonic bond reset at t={current_time_ps:.2f} ps (timestep {timestep}) <<<\n")
                        self.bond_reset_action.enabled = True
                        self.triggered = True
            
            # Add trigger updater (no dt_fs needed - uses simulation time)
            trigger_action = HarmonicResetTrigger(
                bond_reset, 
                self.harmonic_reset_turn_on_time_ps, 
                use_dynamic_temperature,
                self  # Pass CavityMDSimulation object
            )
            trigger_updater = hoomd.update.CustomUpdater(
                action=trigger_action,
                trigger=hoomd.trigger.Periodic(period=10)  # Check every 10 steps
            )
            self.sim.operations.updaters.append(trigger_updater)
            
        except Exception as e:
            import traceback
            self.log_error(f"Failed to setup harmonic bond reset: {e}")
            self.log_error(f"Full traceback: {traceback.format_exc()}")
    
    def _setup_hdf5_output(self):
        """Set up HDF5-based observable output."""
        if not self.enable_hdf5_output:
            self.hdf5_writer = None
            return
        
        try:
            from ..data import ObservableWriter
            
            # Determine output file name
            if self.hdf5_output_file is None:
                self.hdf5_output_file = f"observables_replica_{self.replica}.h5"
            
            self.log_info(f"\n{'='*60}")
            self.log_info("HDF5 OBSERVABLE OUTPUT CONFIGURATION")
            self.log_info(f"{'='*60}")
            self.log_info(f"Output file: {self.hdf5_output_file}")
            self.log_info(f"Output period: {self.hdf5_output_period_ps:.3f} ps")
            self.log_info(f"SWMR mode: Enabled (concurrent read access)")
            self.log_info(f"{'='*60}\n")
            
            # Create HDF5 writer with runtime for SWMR pre-allocation
            self.hdf5_writer = ObservableWriter(
                output_file=self.hdf5_output_file,
                time_tracker=self.time_tracker,
                output_period_ps=self.hdf5_output_period_ps,
                enable_swmr=True,
                runtime_ps=self.runtime_ps
            )
            
            # Register energy tracker if available
            if hasattr(self, 'energy_tracker') and self.energy_tracker is not None:
                self.hdf5_writer.add_energy_tracker(self.energy_tracker)
                self.log_info("  Registered energy tracker for HDF5 output")
            
            # Register temperature tracker if available
            if hasattr(self, 'temperature_tracker') and self.temperature_tracker is not None:
                self.hdf5_writer.add_temperature_tracker(self.temperature_tracker)
                self.log_info("  Registered temperature tracker for HDF5 output")
            
            # Register controllers if available
            if hasattr(self, 'diffeq_controller') and self.diffeq_controller is not None:
                self.hdf5_writer.add_controller('diffeq', self.diffeq_controller)
                self.log_info("  Registered DiffEq controller for HDF5 output")
            
            if hasattr(self, 'simple_setpoint_controller') and self.simple_setpoint_controller is not None:
                self.hdf5_writer.add_controller('simple_setpoint', self.simple_setpoint_controller)
                self.log_info("  Registered SimpleSetpoint controller for HDF5 output")
            
            # Register dipole FDR tracker if available
            if hasattr(self, 'dipole_fdr_tracker') and self.dipole_fdr_tracker is not None:
                self.hdf5_writer.add_dipole_tracker(self.dipole_fdr_tracker)
                self.log_info("  Registered dipole FDR tracker for HDF5 output")
            
            # Add HDF5 writer as custom updater
            hdf5_updater = hoomd.update.CustomUpdater(
                action=self.hdf5_writer,
                trigger=hoomd.trigger.Periodic(period=1)  # Check every step, writer handles timing
            )
            self.sim.operations.updaters.append(hdf5_updater)
            
            self.log_info("✓ HDF5 observable output enabled\n")
            
        except Exception as e:
            import traceback
            self.log_error(f"Failed to setup HDF5 output: {e}")
            self.log_error(f"Full traceback: {traceback.format_exc()}")
            self.hdf5_writer = None
            # CRITICAL: If HDF5 output was explicitly requested, abort simulation
            # Don't continue silently when a requested feature fails
            self.log_error("\n" + "="*60)
            self.log_error("SIMULATION ABORTED: HDF5 output setup failed")
            self.log_error("="*60)
            self.log_error("HDF5 output was explicitly enabled but could not be initialized.")
            self.log_error("The simulation cannot continue without this requested feature.")
            self.log_error("Please fix the error above and try again.")
            raise RuntimeError(f"HDF5 output setup failed: {e}") from e
    def _setup_temperature_tracker(self):
        """Set up comprehensive temperature tracker."""
        try:
            from ..analysis import TemperatureTracker
            
            # Create output file path
            output_file = f"temperature_tracker_replica_{self.replica}.csv"
            
            # Create temperature tracker
            self.temperature_tracker = TemperatureTracker(
                simulation=self,  # Pass CavityMDSimulation object, not HOOMD sim
                time_tracker=self.time_tracker,
                output_period_ps=self.temp_tracker_output_period_ps,
                output_file=output_file,
                energy_tracker=getattr(self, 'energy_tracker', None),
                molecular_thermostat=getattr(self, 'molecular_thermostat_obj', None),
                cavity_thermostat=getattr(self, 'cavity_thermostat_obj', None),
                empirical_data_file=self.temp_tracker_empirical_data_file,
                debug=False,  # Suppress debug output by default
                enable_csv_output=False,  # Use HDF5 output via ObservableWriter
                track_molecular=getattr(self, 'enable_molecular_temps', False)  # Enable molecular tracking if requested
            )
            
            # Add to simulation
            temp_trigger_steps = max(1, int(self.temp_tracker_output_period_ps * 1000))
            temp_updater = hoomd.update.CustomUpdater(
                action=self.temperature_tracker,
                trigger=hoomd.trigger.Periodic(temp_trigger_steps)
            )
            self.sim.operations.updaters.append(temp_updater)
            
            self.log_info(" Comprehensive temperature tracker enabled")
            self.log_info(f"  Output file: {output_file}")
            self.log_info(f"  Output period: {self.temp_tracker_output_period_ps:.1f} ps")
            self.log_info(f"  Tracking: kinetic, harmonic fictive, LJ+Coul fictive, cavity bath, molecular bath")
            if self.temp_tracker_empirical_data_file:
                self.log_info(f"  Empirical data: {self.temp_tracker_empirical_data_file}")
            else:
                self.log_info(f"  Empirical data: None (LJ+Coul fictive temp will be 0)")
            
            # Connect temperature tracker to adaptive coupling variant if needed
            if hasattr(self, '_adaptive_coupling_variant') and self._adaptive_coupling_variant is not None:
                self._adaptive_coupling_variant.temperature_tracker = self.temperature_tracker
                self.log_info(f"Connected temperature tracker to adaptive square wave coupling")
            
            # Update composite variant if it needs temperature tracker for adaptive components
            self._update_composite_adaptive_components()
            
            # Connect temperature tracker to GD controller if needed
            if hasattr(self, 'gd_feedback') and self.gd_feedback is not None:
                self.gd_feedback.temperature_tracker = self.temperature_tracker
                self.log_info(f"Connected temperature tracker to gradient descent controller")
            
        except Exception as e:
            self.log_error(f"Failed to setup temperature tracker: {e}")
            self.temperature_tracker = None
    
    def _setup_molecular_temperature_tracker(self):
        """
        DEPRECATED: Set up molecular temperature decomposition tracker.
        
        This method is deprecated. Molecular temperatures are now tracked by
        TemperatureTracker with track_molecular=True parameter.
        """
        self.log_warning("_setup_molecular_temperature_tracker is deprecated. Use TemperatureTracker with track_molecular=True instead.")
        return
    
    def _update_composite_adaptive_components(self):
        """Update composite variant to use adaptive square wave if requested and temperature tracker is available."""
        if not hasattr(self, 'coupling_variant_type') or self.coupling_variant_type != 'composite':
            return
            
        if not self.composite_square_adaptive:
            return
            
        if not hasattr(self, 'temperature_tracker') or self.temperature_tracker is None:
            return
            
        # Find the composite variant stored reference
        if not hasattr(self, '_composite_coupling_variant') or self._composite_coupling_variant is None:
            return
            
        composite_variant = self._composite_coupling_variant
        if not hasattr(composite_variant, 'variants'):
            return
            
        # Get target temperature for adaptive square component
        if hasattr(self, 'enable_quench_controller') and self.enable_quench_controller:
            target_temperature = getattr(self, 'quench_target_temperature', 100.0)
            temp_source = "quench controller"
        elif hasattr(self, 'enable_gd_feedback') and self.enable_gd_feedback:
            target_temperature = getattr(self, 'gd_target_temperature', 100.0)
            temp_source = "gradient descent controller"
        elif hasattr(self, 'enable_dual_feedback') and self.enable_dual_feedback:
            target_temperature = getattr(self, 'dual_cavity_target_temperature', 100.0)
            temp_source = "dual independent controller (cavity target)"
        else:
            self.log_warning("Cannot update to adaptive square wave: no controller available for target temperature")
            return
            
        # Replace the fixed square wave component with adaptive
        for i, variant in enumerate(composite_variant.variants):
            if hasattr(variant, 'amplitude') and hasattr(variant, 'duty_cycle'):  # This is the square wave
                from ..variants import AdaptiveSquareWaveVariant
                
                # Create adaptive square wave variant
                adaptive_variant = AdaptiveSquareWaveVariant(
                    target_coupling=self.composite_square_amplitude,
                    target_temperature=target_temperature,
                    period_ps=self.composite_square_period,
                    time_tracker=self.time_tracker,
                    temperature_tracker=self.temperature_tracker,
                    duty_cycle=self.composite_square_duty_cycle,
                    phase_offset=0.0,
                    start_time_ps=self.composite_square_start_time,
                    stop_time_ps=self.composite_square_stop_time
                )
                
                # Replace the fixed square wave with adaptive
                composite_variant.variants[i] = adaptive_variant
                
                self.log_info("Updated composite coupling to use adaptive square wave:")
                self.log_info(f"   Target amplitude: {self.composite_square_amplitude:.2e} a.u.")
                self.log_info(f"   Target temperature: {target_temperature:.1f} K (from {temp_source})")
                self.log_info(f"   Temperature tracker connected successfully")
                break
        
    def _setup_dipole_fdr(self):
        """Set up dipole moment FDR tracker."""
        try:
            from ..analysis import DipoleMomentFDRTracker
            
            # Create output file path
            output_file = f"dipole_fdr_replica_{self.replica}.csv"
            
            # Create dipole FDR tracker
            self.dipole_fdr_tracker = DipoleMomentFDRTracker(
                time_tracker=self.time_tracker,
                output_file=output_file,
                max_correlation_time_ps=self.dipole_fdr_max_correlation_time_ps,
                correlation_output_interval_ps=self.dipole_fdr_output_period_ps,
                exclude_cavity=self.dipole_fdr_exclude_cavity,
                field_direction=self.dipole_fdr_field_direction,
                enable_response_measurement=self.enable_dipole_response,
                output_period_ps=self.dipole_fdr_output_period_ps,
                enable_csv_output=False  # Use HDF5 output via ObservableWriter
            )
            
            # Add to simulation
            fdr_trigger_steps = max(1, int(self.dipole_fdr_output_period_ps * 1000))
            fdr_updater = hoomd.update.CustomUpdater(
                action=self.dipole_fdr_tracker,
                trigger=hoomd.trigger.Periodic(fdr_trigger_steps)
            )
            self.sim.operations.updaters.append(fdr_updater)
            
            self.log_info(" Dipole moment FDR tracker enabled")
            self.log_info(f"  Output file: {output_file}")
            self.log_info(f"  Output period: {self.dipole_fdr_output_period_ps:.1f} ps")
            self.log_info(f"  Max correlation time: {self.dipole_fdr_max_correlation_time_ps:.1f} ps")
            self.log_info(f"  Field direction: {self.dipole_fdr_field_direction}")
            self.log_info(f"  Exclude cavity: {self.dipole_fdr_exclude_cavity}")
            if self.enable_dipole_response:
                self.log_info(f"  Response measurement: Enabled (fork-and-clone)")
            else:
                self.log_info(f"  Response measurement: Disabled (autocorrelation only)")
            
        except Exception as e:
            self.log_error(f"Failed to setup dipole FDR tracker: {e}")
            self.dipole_fdr_tracker = None
    
    def _setup_dipole_response(self):
        """Set up dipole response force for FDR measurements."""
        try:
            from .fdr_forces import DipoleResponseForce
            
            # Create dipole response force
            self.dipole_response_force = DipoleResponseForce(
                field_vector=self.dipole_fdr_field_direction,
                field_strength=self.dipole_response_field_strength,
                sign=self.dipole_response_sign,
                exclude_cavity=self.dipole_fdr_exclude_cavity
            )
            
            # Add to simulation forces
            self.sim.operations.integrator.forces.append(self.dipole_response_force)
            
            self.log_info(" Dipole response force enabled")
            self.log_info(f"  Field strength: {self.dipole_response_field_strength:.2e}")
            self.log_info(f"  Field direction: {self.dipole_fdr_field_direction}")
            self.log_info(f"  Sign: {self.dipole_response_sign:+.0f}")
            self.log_info(f"  Exclude cavity: {self.dipole_fdr_exclude_cavity}")
            
        except Exception as e:
            self.log_error(f"Failed to setup dipole response force: {e}")
            self.dipole_response_force = None
    def setup_output_writers(self):
        """Configure GSD writer and console table for simulation output."""
        
        if not self.enable_gsd_output:
            self.log_info("GSD trajectory output disabled (enable_gsd_output=False)")
        else:
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
            
            def __init__(self, sim, time_tracker, performance_tracker, timestep_formatter, adaptive_action, output_period_ps, coupling_constant, incavity, cavity_force, omegac):
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
                self.omegac = omegac  # Cavity frequency for epsilon -> lambda conversion
                self.last_output_time = 0.0
                self.header_printed = False
                
            def _get_current_coupling(self):
                """Get the current coupling value (lambda, dimensionless) from the cavity force."""
                if not self.incavity or self.cavity_force is None:
                    return 0.0
                
                try:
                    # Get lambda value from the variant (dimensionless)
                    lambda_param = self.cavity_force.lambda_coupling_variant
                    lambda_value = None
                    
                    if hasattr(lambda_param, 'current_value'):
                        # Time-varying coupling - get current value from variant
                        lambda_value = lambda_param.current_value
                    elif hasattr(lambda_param, '__call__'):
                        # It's a variant but may not have current_value - call it with current timestep
                        lambda_value = float(lambda_param(self.sim.timestep))
                    else:
                        # Constant coupling
                        lambda_value = float(lambda_param)
                    
                    # Return lambda directly (dimensionless)
                    # Note: Variants now return lambda directly, NOT epsilon
                    return lambda_value if lambda_value is not None else 0.0
                    
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
                        # "Timestep.dt_fs"  # Temporarily disabled for testing
                    ]
                    if self.incavity:
                        header_parts.append("Cavity.lambda")
                    if self.adaptive_action is not None:
                        header_parts.append("Adaptive.error_tolerance")
                    print(" ".join(f"{h:>15s}" for h in header_parts), flush=True)
                    self.header_printed = True
                    
                # Collect current values
                current_time = self._get_current_time(timestep)
                tps = self.sim.tps if hasattr(self.sim, 'tps') else 0.0
                ns_per_day = self.performance_tracker.ns_per_day  # Property, not method
                eta = self.performance_tracker.eta_remaining      # Property, not method
                # dt_fs = self.timestep_formatter.dt_fs             # Property, not method - temporarily disabled
                
                # Build output line
                output_parts = [
                    f"{timestep:15d}",
                    f"{tps:15.5f}",
                    f"{current_time:15.5f}",
                    f"{ns_per_day:>15s}",
                    f"{eta:>15s}",
                    # f"{dt_fs:15.5f}"  # Temporarily disabled
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
            coupling_constant=self.lambda_coupling if self.lambda_coupling is not None else self.couplstr,
            incavity=self.incavity,
            cavity_force=self.cavity_force,
            omegac=self.omegac if hasattr(self, 'omegac') else 1.0
        )
        
        # Add console tracker to simulation (check every step but handle timing internally)
        console_updater = hoomd.update.CustomUpdater(
            action=console_tracker,
            trigger=hoomd.trigger.Periodic(1)  # Check every step but tracker handles timing internally
        )
        self.sim.operations.updaters.append(console_updater)
        
        self.log_info(" Console output setup completed:")
        self.log_info(f"  Output period: {self.console_output_period_ps:.3f} ps (accurate time-based)")
        if self.incavity:
            if hasattr(self, 'lambda_coupling') and self.lambda_coupling is not None:
                self.log_info(f"  Lambda coupling column: {self.lambda_coupling:.6e} (dimensionless) (displayed in real-time)")
            elif self.couplstr is not None:
                self.log_info(f"  Coupling constant column: {self.couplstr:.6e} a.u. (displayed in real-time)")
        self.log_info("  Console tracker handles timing internally using ElapsedTimeTracker")
        self.log_info("  Works accurately for both adaptive and fixed timestep modes")
        
        # Log final setup summary
        if self.error_tolerance > 0:
            self.log_info(" FIXED: Console output now uses precise time-based logic")
            self.log_info("  Both console and energy tracker use accurate timing")
            self.log_info("  GSD output timing may vary slightly but is less critical")
        else:
            self.log_info(" Time-based console output setup for fixed timestep mode")
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
        # Use the appropriate parameter depending on whether deprecated couplstr or new lambda_coupling was used
        if self.use_deprecated_couplstr:
            initial_lambda = 0.0 if self.switch_time_ps is not None else self.couplstr
        else:
            initial_lambda = 0.0 if self.switch_time_ps is not None else self.lambda_coupling

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
                # Formula: q_eq = -(lambda / omegac) * dipole, where lambda is dimensionless
                newpos = -dipmom * initial_lambda / omegac
                newpos[-1] = 0.0
                # Only add thermal fluctuations if coupling is non-zero
                if initial_lambda != 0.0:
                    sigma = np.sqrt(self.kB * self.temperature / omegac**2)
                    newpos = np.random.normal(loc=newpos, scale=sigma, size=3)
                    self.log_info(f"Finite-q mode: Photon displaced by dipole interaction to {newpos} (with thermal fluctuations)")
                else:
                    self.log_info(f"Finite-q mode: Photon at equilibrium position {newpos} (no thermal fluctuations due to zero coupling)")
        else:
            # Start photon at origin (q=0 limit)
            newpos = np.array([0.0, 0.0, 0.0])
            # Only add thermal fluctuations if coupling is non-zero
            if initial_lambda != 0.0:
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
        """Validate that cavity particle exists when required.

        Uses ``State.get_snapshot()`` so typeid is a real host ndarray.
        ``np.asarray`` on a ``HOOMDGPUArray`` from ``gpu_local_snapshot`` yields
        a 0-d object wrapper (empty logical particle buffer), which falsely
        reports zero cavity particles on GPU devices without CuPy.
        """
        particle_types = self.sim.state.particle_types

        if "L" not in particle_types:
            raise ValueError(
                "ERROR: Cavity simulation requested but no cavity particle "
                "type 'L' found in GSD file."
            )

        snap = self.sim.state.get_snapshot()
        typeid_array = np.asarray(snap.particles.typeid)
        positions = np.asarray(snap.particles.position)
        cavity_typeid = int(list(particle_types).index("L"))
        cavity_count = int(np.sum(typeid_array == cavity_typeid))
        if cavity_count == 0:
            raise ValueError(
                "ERROR: Cavity simulation requested but no cavity "
                "particles found in GSD file."
            )
        if cavity_count != 1:
            raise ValueError(
                f"ERROR: Expected exactly 1 cavity particle but found "
                f"{cavity_count} in GSD file."
            )
        cavity_index = int(np.where(typeid_array == cavity_typeid)[0][0])
        cavity_position = positions[cavity_index]
        self.log_info(
            f"Cavity particle validated at index {cavity_index}, "
            f"position {cavity_position}"
        )

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
            
            # Store thermostat objects (extract actual thermostats from methods or refs)
            # For controllers, we need the actual thermostat objects, not the integration methods
            self.molecular_thermostat_obj = (thermostat_refs.get('molecular_bussi') or 
                thermostat_refs.get('molecular_langevin') or 
                    thermostat_refs.get('molecular_mttk') or 
                        molecular_method)
            self.cavity_thermostat_obj = (thermostat_refs.get('cavity_bussi') or 
                thermostat_refs.get('cavity_langevin') or 
                    thermostat_refs.get('cavity_mttk') or 
                        cavity_method)
            self.thermostat_refs = thermostat_refs
            
            # Set up quench controller now that thermostat objects exist
            if getattr(self, 'enable_quench_controller', False):
                self._setup_quench_controller()
            
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
            
            # Set up offset controller now that energy tracker exists
            if getattr(self, 'enable_offset_controller', False):
                self._setup_offset_controller()
            
            # Set up PID controller now that temperature tracker exists
            if getattr(self, 'enable_pid_controller', False):
                self.log_info(f"DEBUG: Attempting to set up PID controller (enable_pid_controller={getattr(self, 'enable_pid_controller', False)})")
                self._setup_pid_controller()
                self.log_info("DEBUG: PID controller setup completed")
            else:
                self.log_info(f"DEBUG: PID controller NOT enabled (enable_pid_controller={getattr(self, 'enable_pid_controller', False)})")
            
            # Phase 5: Setup output writers
            self.log_info("=== Phase 5: Setting up output writers ===")
            self.setup_output_writers()
            
            # Phase 6: Run simulation
            self.log_info("=== Phase 6: Running simulation ===")
            try:
                self.run_simulation()
            except StopIteration:
                # Normal exit when runtime is reached
                self.log_info("Simulation runtime reached - normal exit")
            
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
        
        # Set up delayed controller monitor if needed
        if self._delayed_controllers:
            self._delayed_controller_monitor = self._create_delayed_controller_monitor()
            monitor_updater = hoomd.update.CustomUpdater(
                action=self._delayed_controller_monitor,
                trigger=hoomd.trigger.Periodic(100)  # Check every 100 timesteps
            )
            self.sim.operations.updaters.append(monitor_updater)
            self.log_info(f"Delayed controller monitor enabled (checking every 100 steps)")
        
        # Run the simulation
        self.sim.run(total_steps, write_at_start=True)
        
        self.log_info(f"Simulation completed successfully")

    def cleanup(self):
        """Handle post-simulation cleanup and restore original directory."""
        self.log_info("Cleanup initiated...")
        
        # Close HDF5 writer if it exists
        if hasattr(self, 'hdf5_writer') and self.hdf5_writer is not None:
            try:
                self.log_info("Closing HDF5 observable writer...")
                self.hdf5_writer.close()
            except Exception as e:
                self.log_error(f"Error closing HDF5 writer: {e}")
        
        # Finalize F(k,t) logarithmic output if enabled
        print("DEBUG: Cleanup called - checking for fkt_tracker")
        if hasattr(self, 'fkt_tracker') and self.fkt_tracker is not None:
            print(f"DEBUG: fkt_tracker found: {type(self.fkt_tracker)}")
            if hasattr(self.fkt_tracker, 'finalize_output'):
                print("DEBUG: Calling fkt_tracker.finalize_output()")
                self.log_info("Finalizing F(k,t) logarithmic output...")
                self.fkt_tracker.finalize_output()
            else:
                print("DEBUG: fkt_tracker has no finalize_output method")
        else:
            print("DEBUG: No fkt_tracker found or it's None")
        
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
