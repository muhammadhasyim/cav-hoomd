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
from .variants import StepVariant, PeriodicVariant, ExponentialDecayVariant, SquareWaveVariant, DecayingSquareWaveVariant, AdaptiveSquareWaveVariant
from .composite_variant import CompositeVariant


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
                 # F(k,t) logarithmic time spacing parameters
                 fkt_log_time_spacing: bool = False,
                 fkt_min_log_time_ps: Optional[float] = None,
                 fkt_max_log_time_ps: Optional[float] = None,
                 fkt_log_num_points: int = 50,
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
                 # Differential equation controller parameters
                 enable_diffeq_controller: bool = False,
                 diffeq_temperature_method: str = 'kinetic',
                 diffeq_time_constant_ps: float = 5.0,
                 diffeq_turn_on_time_ps: float = 0.0,
                 diffeq_turn_off_time_ps: Optional[float] = None,
                 diffeq_update_interval_ps: float = 0.1,
                 diffeq_apply_to: str = 'both',
                 diffeq_T_min: float = 0.0,
                 diffeq_T_max: Optional[float] = None,
                 diffeq_rate_limit_K_per_ps: Optional[float] = None,
                 # Temperature tracker parameters
                 enable_temp_tracker: bool = False,
                 temp_tracker_output_period_ps: float = 0.1,
                 temp_tracker_empirical_data_file: Optional[str] = None,
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
        
        # Differential equation controller parameters
        self.enable_diffeq_controller = enable_diffeq_controller
        self.diffeq_temperature_method = diffeq_temperature_method
        self.diffeq_time_constant_ps = diffeq_time_constant_ps
        self.diffeq_turn_on_time_ps = diffeq_turn_on_time_ps
        self.diffeq_turn_off_time_ps = diffeq_turn_off_time_ps
        self.diffeq_update_interval_ps = diffeq_update_interval_ps
        self.diffeq_apply_to = diffeq_apply_to
        self.diffeq_T_min = diffeq_T_min
        self.diffeq_T_max = diffeq_T_max
        self.diffeq_rate_limit_K_per_ps = diffeq_rate_limit_K_per_ps
        
        # Temperature tracker parameters
        self.enable_temp_tracker = enable_temp_tracker
        self.temp_tracker_output_period_ps = temp_tracker_output_period_ps
        self.temp_tracker_empirical_data_file = temp_tracker_empirical_data_file
        
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
                            self.simulation.log_info(f"✅ Delayed controller activated: {controller_data['controller_name']} at t = {current_time_ps:.2f} ps")
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

    def _create_coupling_variant(self):
        """Create coupling variant based on coupling_variant_type."""
        if not hasattr(self, 'time_tracker') or self.time_tracker is None:
            raise ValueError("Time tracker must be set up before creating coupling variants")
        
        variant_type = self.coupling_variant_type.lower()
        
        if variant_type == 'constant':
            # Constant coupling (default)
            from hoomd.variant import Constant
            coupling_variant = Constant(self.couplstr)
            self.log_info(f"Using constant coupling: {self.couplstr} a.u.")
        
        elif variant_type == 'step':
            # Step coupling (with optional turn-off and decay)  
            coupling_variant = StepVariant(
                target_value=self.couplstr,
                switch_time_ps=self.switch_time_ps if self.switch_time_ps is not None else 0.0,
                time_tracker=self.time_tracker,
                decay_time_constant_ps=self.decay_time_constant_ps,
                turn_off_time_ps=getattr(self, 'step_turn_off_time_ps', None)
            )
            self.log_info(f"Using step coupling:")
            self.log_info(f"  Target value: {self.couplstr} a.u.")
            self.log_info(f"  Switch time: {self.switch_time_ps} ps")
            if self.decay_time_constant_ps:
                self.log_info(f"  Decay time constant: {self.decay_time_constant_ps} ps")
        
        elif variant_type == 'periodic':
            # Periodic coupling  
            coupling_variant = PeriodicVariant(
                amplitude=self.couplstr,
                time_tracker=self.time_tracker,
                period_ps=getattr(self, 'periodic_period_ps', 1.0),
                phase_offset=getattr(self, 'periodic_phase_offset', 0.0),
                start_time_ps=getattr(self, 'periodic_start_time_ps', 0.0),
                stop_time_ps=getattr(self, 'periodic_stop_time_ps', None)
            )
            self.log_info(f"Using periodic coupling:")
            self.log_info(f"  Amplitude: {self.couplstr} a.u.")
            self.log_info(f"  Period: {getattr(self, 'periodic_period_ps', 1.0)} ps")
            self.log_info(f"  Phase offset: {getattr(self, 'periodic_phase_offset', 0.0):.3f} rad")
        
        elif variant_type == 'exponential':
            # Exponential decay coupling
            coupling_variant = ExponentialDecayVariant(
                amplitude=self.exponential_amplitude,
                time_tracker=self.time_tracker,
                decay_time_constant_ps=self.exponential_decay_time_ps,
                turn_on_time_ps=self.exponential_turn_on_time_ps,
                turn_off_time_ps=self.exponential_turn_off_time_ps
            )
            self.log_info(f"Using exponential decay coupling:")
            self.log_info(f"  Amplitude: {self.exponential_amplitude} a.u.")
            self.log_info(f"  Decay time constant: {self.exponential_decay_time_ps} ps")
            self.log_info(f"  Turn-on time: {self.exponential_turn_on_time_ps} ps")
        
        elif variant_type == 'square':
            # Square wave coupling
            coupling_variant = SquareWaveVariant(
                amplitude=self.couplstr,
                period_ps=self.square_period_ps,
                time_tracker=self.time_tracker,
                duty_cycle=self.square_duty_cycle,
                phase_offset=self.square_phase_offset,
                start_time_ps=self.square_start_time_ps,
                stop_time_ps=self.square_stop_time_ps
            )
            self.log_info(f"Using square wave coupling:")
            self.log_info(f"  Amplitude: {self.couplstr} a.u.")
            self.log_info(f"  Period: {self.square_period_ps} ps")
            self.log_info(f"  Duty cycle: {self.square_duty_cycle:.1%}")
        
        elif variant_type == 'decaying_square':
            # Decaying square wave coupling
            coupling_variant = DecayingSquareWaveVariant(
                initial_amplitude=self.couplstr,
                period_ps=self.decaying_square_period_ps,
                time_tracker=self.time_tracker,
                decay_rate_per_period=self.decaying_square_decay_rate,
                duty_cycle=self.decaying_square_duty_cycle,
                phase_offset=self.decaying_square_phase_offset,
                start_time_ps=self.decaying_square_start_time_ps,
                stop_time_ps=self.decaying_square_stop_time_ps,
                minimum_amplitude=self.decaying_square_minimum_amplitude
            )
            self.log_info(f"Using decaying square wave coupling:")
            self.log_info(f"  Initial amplitude: {self.couplstr} a.u.")
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
            from .variants import ExponentialWaveVariant
            
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
            from .utils import PhysicalConstants
            omegac = self.freq / PhysicalConstants.HARTREE_TO_CM_MINUS1
            
            # Create coupling variant based on type selection
            coupling_variant = self._create_coupling_variant()
            
            # Create cavity force with variant
            from .forces import CavityForce
            cavityforce = CavityForce(
                kvector=np.array([0,0,1]), 
                couplstr=coupling_variant, 
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
        
        # Setup enhanced laser driving if enabled
        if self.laser_enabled:
            laser_force = self._create_laser_force()
            forces.append(laser_force)
        
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
                    verbose="quiet"  # Suppress debug output by default
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
            from .analysis import EmpiricalTemperatureData, EmpiricalTemperatureFeedback
            
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
                enabled_features.append(f"differential equation controller ({self.diffeq_temperature_method}, τ={self.diffeq_time_constant_ps:.1f}ps)")
            else:
                self.log_info(f"  WARNING: DiffEq controller conflicts with immediate controllers: {', '.join(immediate_conflicts)} - disabling diffeq")
                self.diffeq_controller = None
        
        # Note: Quench controller setup moved to after thermostat creation
        # to ensure proper thermostat object references
        if getattr(self, 'enable_quench_controller', False):
            enabled_features.append(f"quench controller ({self.quench_initial_temperature:.1f}K→{self.quench_target_temperature:.1f}K at {self.quench_time_ps:.1f}ps)")
        
        # Set up comprehensive temperature tracker if enabled
        if getattr(self, 'enable_temp_tracker', False):
            self._setup_temperature_tracker()
            enabled_features.append(f"comprehensive temperature tracker ({self.temp_tracker_output_period_ps:.1f} ps)")
        
        # Set up molecular temperature decomposition if enabled
        if getattr(self, 'enable_molecular_temps', False):
            self._setup_molecular_temperature_tracker()
            enabled_features.append(f"molecular temperature decomposition ({self.molecular_temps_output_period_ps:.1f} ps)")
            
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
            from .analysis import GradientDescentTemperatureFeedback
            
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
            from .analysis import DualIndependentTemperatureFeedback
            
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
    
    def _setup_sinusoidal_bath_controller(self):
        """Set up sinusoidal bath temperature controller."""
        if not self.enable_sinusoidal_bath:
            self.sinusoidal_bath_controller = None
            return
        
        try:
            from .analysis import SinusoidalBathTemperatureController
            
            # Create output file path
            output_file = f"sinusoidal_bath_controller_replica_{self.replica}.csv"
            
            # Create sinusoidal bath controller
            self.sinusoidal_bath_controller = SinusoidalBathTemperatureController(
                period_ps=self.sinusoidal_bath_period_ps,
                amplitude_scale=self.sinusoidal_bath_amplitude_scale,
                phase_offset=self.sinusoidal_bath_phase_offset,
                time_tracker=self.time_tracker,
                energy_tracker=getattr(self, 'energy_tracker', None),
                simulation=self.sim,
                molecular_thermostat=getattr(self, 'molecular_thermostat_obj', None),
                cavity_thermostat=getattr(self, 'cavity_thermostat_obj', None),
                apply_to=self.sinusoidal_bath_apply_to,
                target_temperature=self.sinusoidal_bath_target_temperature,
                dynamic_target=self.sinusoidal_bath_dynamic_target,
                turn_on_time_ps=self.sinusoidal_bath_turn_on_time_ps,
                turn_off_time_ps=self.sinusoidal_bath_turn_off_time_ps,
                update_interval_ps=self.sinusoidal_bath_update_interval_ps,
                T_min=self.sinusoidal_bath_T_min,
                T_max=self.sinusoidal_bath_T_max,
                output_file=output_file,
                empirical_data_file=self.sinusoidal_bath_empirical_data_file,
                console_output_period_ps=self.console_output_period_ps,
                amplitude_update_interval_ps=self.sinusoidal_bath_amplitude_update_interval_ps,
                amplitude_temperature_method=self.sinusoidal_bath_amplitude_temperature_method,
                adaptive_range_mode=self.sinusoidal_bath_adaptive_range_mode
            )
            
            print(f"Sinusoidal bath temperature controller enabled:")
            print(f"  Period: {self.sinusoidal_bath_period_ps:.2f} ps")
            print(f"  Amplitude scale: {self.sinusoidal_bath_amplitude_scale:.3f}")
            print(f"  Mode: {'Adaptive Range' if self.sinusoidal_bath_adaptive_range_mode else 'Fixed Mean'}")
            print(f"  Amplitude method: {self.sinusoidal_bath_amplitude_temperature_method}")
            print(f"  Phase offset: {self.sinusoidal_bath_phase_offset:.2f} rad")
            print(f"  Apply to: {self.sinusoidal_bath_apply_to}")
            print(f"  Turn on time: {self.sinusoidal_bath_turn_on_time_ps:.1f} ps")
            print(f"  Dynamic target: {self.sinusoidal_bath_dynamic_target}")
            if not self.sinusoidal_bath_dynamic_target:
                print(f"  Static target: {self.sinusoidal_bath_target_temperature:.1f} K")
            print(f"  Temperature range: {self.sinusoidal_bath_T_min:.1f} - {self.sinusoidal_bath_T_max or 'inf'} K")
            print(f"  Output file: {output_file}")
            
            # Add to simulation with appropriate trigger frequency
            # Use update interval to determine trigger frequency
            from .utils import PhysicalConstants
            dt_ps = PhysicalConstants.atomic_units_to_ps(self.sim.operations.integrator.dt)
            steps_per_ps = 1.0 / dt_ps
            trigger_steps = max(1, int(self.sinusoidal_bath_update_interval_ps * steps_per_ps))
            
            sinusoidal_updater = hoomd.update.CustomUpdater(
                action=self.sinusoidal_bath_controller,
                trigger=hoomd.trigger.Periodic(trigger_steps)
            )
            self.sim.operations.updaters.append(sinusoidal_updater)
            
            self.log_info(" Sinusoidal bath temperature controller enabled")
            self.log_info(f"  Update interval: {self.sinusoidal_bath_update_interval_ps:.3f} ps")
            self.log_info(f"  Trigger steps: {trigger_steps}")
            
        except Exception as e:
            print(f"Warning: Failed to set up sinusoidal bath controller: {e}")
            import traceback
            traceback.print_exc()
            self.sinusoidal_bath_controller = None

    def _setup_adaptive_bath_controller(self):
        """Set up adaptive bath temperature controller."""
        if not self.enable_adaptive_bath:
            self.adaptive_bath_controller = None
            return
        
        try:
            from .analysis import AdaptiveBathTemperatureController
            
            # Create output file path
            output_file = f"adaptive_bath_controller_replica_{self.replica}.csv"
            
            # Create adaptive bath controller
            self.adaptive_bath_controller = AdaptiveBathTemperatureController(
                amplitude_scale=self.adaptive_bath_amplitude_scale,
                time_constant_ps=self.adaptive_bath_time_constant_ps,
                time_tracker=self.time_tracker,
                energy_tracker=getattr(self, 'energy_tracker', None),
                simulation=self.sim,
                molecular_thermostat=getattr(self, 'molecular_thermostat_obj', None),
                cavity_thermostat=getattr(self, 'cavity_thermostat_obj', None),
                apply_to=self.adaptive_bath_apply_to,
                target_temperature=self.adaptive_bath_target_temperature,
                dynamic_target=self.adaptive_bath_dynamic_target,
                turn_on_time_ps=self.adaptive_bath_turn_on_time_ps,
                turn_off_time_ps=self.adaptive_bath_turn_off_time_ps,
                update_interval_ps=self.adaptive_bath_update_interval_ps,
                T_min=self.adaptive_bath_T_min,
                T_max=self.adaptive_bath_T_max,
                output_file=output_file,
                empirical_data_file=self.adaptive_bath_empirical_data_file,
                console_output_period_ps=self.console_output_period_ps,
                signal_temperature_method=self.adaptive_bath_signal_temperature_method
            )
            
            print(f"Adaptive bath temperature controller enabled (GD Framework):")
            print(f"  Amplitude scale: {self.adaptive_bath_amplitude_scale:.3f}")
            print(f"  Time constant: {self.adaptive_bath_time_constant_ps:.3f} ps")
            print(f"  Signal method: {self.adaptive_bath_signal_temperature_method}")
            print(f"  Apply to: {self.adaptive_bath_apply_to}")
            print(f"  Turn on time: {self.adaptive_bath_turn_on_time_ps:.1f} ps")
            print(f"  Dynamic target: {self.adaptive_bath_dynamic_target}")
            if not self.adaptive_bath_dynamic_target:
                print(f"  Static target: {self.adaptive_bath_target_temperature:.1f} K")
            print(f"  Temperature range: {self.adaptive_bath_T_min:.1f} - {self.adaptive_bath_T_max or 'inf'} K")
            print(f"  Output file: {output_file}")
            
            # Add to simulation with appropriate trigger frequency
            from .utils import PhysicalConstants
            dt_ps = PhysicalConstants.atomic_units_to_ps(self.sim.operations.integrator.dt)
            trigger_steps = max(1, int(self.adaptive_bath_update_interval_ps / dt_ps))
            
            adaptive_bath_updater = hoomd.update.CustomUpdater(
                action=self.adaptive_bath_controller,
                trigger=hoomd.trigger.Periodic(trigger_steps)
            )
            self.sim.operations.updaters.append(adaptive_bath_updater)
            
            self.log_info(f"  Update interval: {self.adaptive_bath_update_interval_ps:.3f} ps")
            self.log_info(f"  Trigger steps: {trigger_steps}")
            
        except Exception as e:
            print(f"Warning: Failed to set up adaptive bath controller: {e}")
            import traceback
            traceback.print_exc()
            self.adaptive_bath_controller = None
    
    def _setup_quench_controller(self):
        """Set up quench controller for instantaneous temperature changes."""
        if not self.enable_quench_controller:
            self.quench_controller = None
            return
        
        try:
            from .analysis import QuenchController
            
            # Create output file path
            output_file = f"quench_control_replica_{self.replica}.csv"
            
            # Create quench controller
            self.quench_controller = QuenchController(
                initial_temperature=self.quench_initial_temperature,
                target_temperature=self.quench_target_temperature,
                quench_time_ps=self.quench_time_ps,
                time_tracker=self.time_tracker,
                molecular_thermostat=getattr(self, 'molecular_thermostat_obj', None),
                cavity_thermostat=getattr(self, 'cavity_thermostat_obj', None),
                apply_to=self.quench_apply_to,
                output_file=output_file,
                console_output_period_ps=self.console_output_period_ps
            )
            
            # Add to simulation with high frequency trigger (every timestep)
            quench_trigger_steps = 1  # Check every timestep for precise timing
            quench_updater = hoomd.update.CustomUpdater(
                action=self.quench_controller,
                trigger=hoomd.trigger.Periodic(quench_trigger_steps)
            )
            self.sim.operations.updaters.append(quench_updater)
            
            self.log_info(" Quench controller enabled")
            self.log_info(f"  Initial temperature: {self.quench_initial_temperature:.1f} K")
            self.log_info(f"  Target temperature: {self.quench_target_temperature:.1f} K")
            self.log_info(f"  Quench time: {self.quench_time_ps:.1f} ps")
            self.log_info(f"  Applied to: {self.quench_apply_to}")
            self.log_info(f"  ΔT = {self.quench_target_temperature - self.quench_initial_temperature:+.1f} K")
            
        except Exception as e:
            import traceback
            self.log_error(f"Failed to setup quench controller: {e}")
            self.log_error(f"Full traceback: {traceback.format_exc()}")
            self.quench_controller = None
    
    def _setup_offset_controller(self):
        """Set up offset temperature controller for fixed temperature offsets."""
        if not self.enable_offset_controller:
            self.offset_controller = None
            return
        
        try:
            from .analysis import OffsetTemperatureController
            
            # Create output file path
            output_file = f"offset_control_replica_{self.replica}.csv"
            
            # Create offset controller
            self.offset_controller = OffsetTemperatureController(
                temperature_method=self.offset_temperature_method,
                temperature_offset_K=self.offset_temperature_offset_K,
                time_tracker=self.time_tracker,
                energy_tracker=getattr(self, 'energy_tracker', None),
                simulation=self.sim,
                molecular_thermostat=getattr(self, 'molecular_thermostat_obj', None),
                cavity_thermostat=getattr(self, 'cavity_thermostat_obj', None),
                apply_to=self.offset_apply_to,
                turn_on_time_ps=self.offset_turn_on_time_ps,
                turn_off_time_ps=self.offset_turn_off_time_ps,
                update_interval_ps=self.offset_update_interval_ps,
                T_min=self.offset_T_min,
                T_max=self.offset_T_max,
                output_file=output_file,
                empirical_data_file=getattr(self, 'temp_tracker_empirical_data_file', None),
                console_output_period_ps=self.console_output_period_ps
            )
            
            # Add to simulation with appropriate trigger frequency
            # Calculate steps per ps from current timestep
            dt_ps = PhysicalConstants.atomic_units_to_ps(self.sim.operations.integrator.dt)
            steps_per_ps = 1.0 / dt_ps
            
            # Use update interval to determine trigger frequency
            trigger_steps = max(1, int(self.offset_update_interval_ps * steps_per_ps))
            offset_updater = hoomd.update.CustomUpdater(
                action=self.offset_controller,
                trigger=hoomd.trigger.Periodic(trigger_steps)
            )
            self.sim.operations.updaters.append(offset_updater)
            
            self.log_info(" Offset temperature controller enabled")
            self.log_info(f"  Method: {self.offset_temperature_method}")
            self.log_info(f"  Temperature offset: {self.offset_temperature_offset_K:+.1f} K")
            self.log_info(f"  Turn on time: {self.offset_turn_on_time_ps:.1f} ps")
            self.log_info(f"  Update interval: {self.offset_update_interval_ps:.3f} ps")
            self.log_info(f"  Apply to: {self.offset_apply_to}")
            
        except Exception as e:
            import traceback
            self.log_error(f"Failed to setup offset controller: {e}")
            self.log_error(f"Full traceback: {traceback.format_exc()}")
            self.offset_controller = None
    
    def _setup_diffeq_controller(self):
        """Set up differential equation temperature controller."""
        if not self.enable_diffeq_controller:
            self.diffeq_controller = None
            return
        
        try:
            from .analysis import DiffEqController
            
            # Create output file path
            output_file = f"diffeq_control_replica_{self.replica}.csv"
            
            # Create differential equation controller
            self.diffeq_controller = DiffEqController(
                temperature_method=self.diffeq_temperature_method,
                time_constant_ps=self.diffeq_time_constant_ps,
                time_tracker=self.time_tracker,
                energy_tracker=getattr(self, 'energy_tracker', None),
                simulation=self.sim,
                molecular_thermostat=getattr(self, 'molecular_thermostat_obj', None),
                cavity_thermostat=getattr(self, 'cavity_thermostat_obj', None),
                apply_to=self.diffeq_apply_to,
                turn_on_time_ps=self.diffeq_turn_on_time_ps,
                turn_off_time_ps=self.diffeq_turn_off_time_ps,
                update_interval_ps=self.diffeq_update_interval_ps,
                T_min=self.diffeq_T_min,
                T_max=self.diffeq_T_max,
                rate_limit_K_per_ps=self.diffeq_rate_limit_K_per_ps,
                output_file=output_file,
                empirical_data_file=getattr(self, 'temp_tracker_empirical_data_file', None),
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
            self.log_info(f"  Time constant: {self.diffeq_time_constant_ps:.2f} ps")
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
            
        except Exception as e:
            import traceback
            self.log_error(f"Failed to setup differential equation controller: {e}")
            self.log_error(f"Full traceback: {traceback.format_exc()}")
            self.diffeq_controller = None
    
    def _setup_temperature_tracker(self):
        """Set up comprehensive temperature tracker."""
        try:
            from .analysis import TemperatureTracker
            
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
                debug=False  # Suppress debug output by default
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
        """Set up molecular temperature decomposition tracker."""
        try:
            from .molecular_temperatures import DiatomicMolecularTemperatures
            
            # Create output file path
            output_file = f"molecular_temperatures.csv"
            
            # Create molecular temperature tracker
            self.molecular_temp_tracker = DiatomicMolecularTemperatures(
                simulation=self.sim,  # Pass HOOMD sim object
                time_tracker=self.time_tracker,
                output_period_ps=self.molecular_temps_output_period_ps,
                output_file=output_file,
                debug=False  # Suppress debug output by default
            )
            
            # Add to simulation
            # Use time-based trigger that works for both adaptive and fixed timestep
            mol_temp_trigger_steps = max(1, int(self.molecular_temps_output_period_ps * 1000))
            mol_temp_writer = hoomd.write.CustomWriter(
                action=self.molecular_temp_tracker,
                trigger=hoomd.trigger.Periodic(mol_temp_trigger_steps)
            )
            self.sim.operations.writers.append(mol_temp_writer)
            
            self.log_info(" Molecular temperature decomposition enabled")
            self.log_info(f"  Output file: {output_file}")
            self.log_info(f"  Output period: {self.molecular_temps_output_period_ps:.1f} ps")
            self.log_info(f"  Tracking: T_trans, T_rot, T_vib for O-O and N-N dimers")
            
        except Exception as e:
            import traceback
            self.log_error(f"Failed to setup molecular temperature tracker: {e}")
            self.log_error(f"Full traceback: {traceback.format_exc()}")
            self.molecular_temp_tracker = None
    
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
                from .variants import AdaptiveSquareWaveVariant
                
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
                
                self.log_info("✅ Updated composite coupling to use adaptive square wave:")
                self.log_info(f"   Target amplitude: {self.composite_square_amplitude:.2e} a.u.")
                self.log_info(f"   Target temperature: {target_temperature:.1f} K (from {temp_source})")
                self.log_info(f"   Temperature tracker connected successfully")
                break
        
    def _setup_auto_stop_controller(self):
        """Set up auto-stop controller for coupling convergence detection."""
        try:
            from .analysis import AutoStopController
            
            # Ensure temperature tracker is available
            if not hasattr(self, 'temperature_tracker') or self.temperature_tracker is None:
                self.log_error("Auto-stop controller requires temperature tracker to be set up first")
                self.auto_stop_controller = None
                return
            
            # Determine coupling start time based on coupling variant
            coupling_start_time = 0.0
            if self.coupling_variant_type == 'step' and self.switch_time_ps is not None:
                coupling_start_time = self.switch_time_ps
            elif self.coupling_variant_type == 'square':
                coupling_start_time = self.square_start_time_ps
            elif self.coupling_variant_type == 'decaying_square':
                coupling_start_time = self.decaying_square_start_time_ps
            elif self.coupling_variant_type == 'adaptive_square':
                coupling_start_time = self.adaptive_square_start_time_ps
            elif self.coupling_variant_type == 'periodic':
                coupling_start_time = self.periodic_start_time_ps
            elif self.coupling_variant_type == 'exponential':
                coupling_start_time = self.exponential_turn_on_time_ps
            elif self.coupling_variant_type == 'exponentialwave':
                coupling_start_time = self.exp_start_time_ps
            
            # Determine target temperature for auto-stop (from active controller)
            target_temp = None
            if hasattr(self, 'gd_target_temperature') and self.gd_target_temperature is not None:
                target_temp = self.gd_target_temperature
            elif hasattr(self, 'enable_dual_feedback') and self.enable_dual_feedback and hasattr(self, 'dual_cavity_target_temperature'):
                # For dual controller, use cavity target (as it drives the adaptive square wave)
                target_temp = self.dual_cavity_target_temperature
            elif hasattr(self, 'quench_target_temperature') and self.quench_target_temperature is not None:
                target_temp = self.quench_target_temperature
            else:
                target_temp = self.temperature  # Fall back to system temperature
            
            # Create auto-stop controller
            self.auto_stop_controller = AutoStopController(
                temperature_tracker=self.temperature_tracker,
                time_tracker=self.time_tracker,
                tolerance=self.auto_stop_tol,
                window_ps=self.auto_stop_window,
                coupling_start_time_ps=coupling_start_time,
                target_temperature=target_temp
            )
            
            # Set simulation reference for auto-stop signaling
            self.auto_stop_controller.set_simulation(self)
            
            # Initialize auto-stop signal
            self.auto_stop_coupling_signal = False
            
            self.log_info(f"Auto-stop controller setup complete:")
            self.log_info(f"  Tolerance: {self.auto_stop_tol:.1f} K")
            self.log_info(f"  Averaging window: {self.auto_stop_window:.1f} ps")
            self.log_info(f"  Coupling start time: {coupling_start_time:.1f} ps")
            self.log_info(f"  Monitoring: |T_fictive_avg - T_kinetic_avg| < {self.auto_stop_tol:.1f} K")
            
        except Exception as e:
            self.log_error(f"Failed to setup auto-stop controller: {e}")
            import traceback
            self.log_error(f"Full traceback: {traceback.format_exc()}")
            self.auto_stop_controller = None
    
    def _setup_dipole_fdr(self):
        """Set up dipole moment FDR tracker."""
        try:
            from .analysis import DipoleMomentFDRTracker
            
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
                enable_response_measurement=self.enable_dipole_response
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
                        # "Timestep.dt_fs"  # Temporarily disabled for testing
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
        
        self.log_info(" Console output setup completed:")
        self.log_info(f"  Output period: {self.console_output_period_ps:.3f} ps (accurate time-based)")
        if self.incavity:
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
            
            # Store thermostat objects for empirical feedback
            self.molecular_thermostat_obj = molecular_method
            self.cavity_thermostat_obj = cavity_method
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
                    self.log_info(" Forced thermalization successful - Bussi thermostat should work now")
                else:
                    self.log_error(" Forced thermalization failed - velocities are still zero!")
                    raise RuntimeError("Unable to initialize non-zero velocities for Bussi thermostat")



class AdaptiveTimestepUpdater(hoomd.custom.Action):
    """Update timestep adaptively based on energy conservation and stability."""
    
    def __init__(self, state, integrator, error_tolerance, time_constant_ps=50.0, 
                 initial_fraction=0.01, adaptiveerror=True, cavity_damping_factor=1.0, 
                 molecular_thermostat_tau=5.0, cavity_thermostat_tau=5.0, time_tracker=None,
                 switch_time_ps=None, timestep_change_threshold=0.1, max_timestep_change_factor=1.5,
                 shock_dampening_factor=1e-3, shock_dampening_enabled=True, shock_dampening_time_constant_ps=50.0,
                 # New dynamic coupling change detection parameters
                 dynamic_coupling_detection=True, coupling_change_threshold=1e-5):
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
        
        # Initialize current_error_tolerance based on shock dampening mode
        if dynamic_coupling_detection:
            # Dynamic mode: start with target tolerance, only drop on actual shocks
            self.current_error_tolerance = self.target_error_tolerance
        elif switch_time_ps is not None:
            # Legacy switch time mode: start with target tolerance before switch
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
        
        # NEW: Dynamic coupling change detection
        self.dynamic_coupling_detection = dynamic_coupling_detection
        self.coupling_change_threshold = coupling_change_threshold  # a.u. - threshold for large coupling changes
        self.last_coupling_value = None  # Track previous coupling strength
        self.coupling_shock_detected_time = None  # When we last detected a coupling shock
        self.coupling_recovery_active = False  # Whether we're currently in recovery mode
        self.first_coupling_check = True  # Flag to initialize without triggering shock
        
        # Log the shock dampening behavior
        if dynamic_coupling_detection:
            print(f"Dynamic coupling change detection enabled", flush=True)
            print(f"Coupling change threshold: ±{coupling_change_threshold:.2e} a.u.", flush=True)
            print(f"Shock dampening: error_tolerance drops to {self.shock_error_tolerance:.2e} (factor: {shock_dampening_factor:.2e})", flush=True)
            print(f"Recovery: exponential ramp to {self.target_error_tolerance:.2e} with τ = {time_constant_ps} ps", flush=True)
        elif shock_dampening_enabled and switch_time_ps is not None:
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
    
    def _get_current_coupling_strength(self):
        """Get the current coupling strength from cavity force."""
        try:
            # Look for cavity force in the integrator's forces
            for force in self.integrator.forces:
                if hasattr(force, 'couplstr_variant'):
                    # This is the cavity force - get the variant
                    coupling_variant = force.couplstr_variant
                    if hasattr(coupling_variant, '__call__'):
                        # It's a variant - call it with current timestep
                        # Access timestep through the integrator's simulation
                        if hasattr(self.integrator, 'simulation'):
                            current_timestep = self.integrator.simulation.timestep
                        else:
                            # Fallback - use a reasonable default
                            current_timestep = 0
                        return coupling_variant(current_timestep)
                    else:
                        # It's a constant value
                        return float(coupling_variant.value) if hasattr(coupling_variant, 'value') else float(coupling_variant)
            return None
        except Exception as e:
            # If we can't get coupling strength, return None
            return None

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
        
        # NEW: Detect sudden coupling changes OR switch time crossing
        force_immediate_update = False
        coupling_shock_detected = False
        
        # Dynamic coupling change detection
        if self.dynamic_coupling_detection:
            # Get current coupling value from cavity force
            current_coupling = self._get_current_coupling_strength()
            
            if self.first_coupling_check:
                # First time - just store the value without triggering shock
                self.last_coupling_value = current_coupling
                self.first_coupling_check = False
            elif self.last_coupling_value is not None and current_coupling is not None:
                coupling_change = abs(current_coupling - self.last_coupling_value)
                
                
                if coupling_change >= self.coupling_change_threshold:
                    # Large coupling change detected!
                    coupling_shock_detected = True
                    force_immediate_update = True
                    self.coupling_shock_detected_time = current_elapsed_time_ps
                    self.coupling_recovery_active = True
                    
                    print(f"COUPLING SHOCK DETECTED at timestep {timestep}: t = {current_elapsed_time_ps:.6f} ps", flush=True)
                    print(f"  Coupling change: {self.last_coupling_value:.6e} → {current_coupling:.6e} a.u. (Δ = {coupling_change:.6e})", flush=True)
                    print(f"  Triggering shock dampening (threshold: ±{self.coupling_change_threshold:.2e} a.u.)", flush=True)
                
                # Update last coupling value for next iteration
                self.last_coupling_value = current_coupling
        
        # Legacy switch time detection (still supported)
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
            # NEW: Dynamic coupling shock recovery
            if self.dynamic_coupling_detection:
                if coupling_shock_detected:
                    # Just detected a coupling shock - immediately apply shock tolerance
                    self.current_error_tolerance = self.shock_error_tolerance
                elif self.coupling_recovery_active:
                    # In recovery mode - exponential recovery from shock tolerance
                    recovery_time = current_elapsed_time_ps - self.coupling_shock_detected_time
                    if recovery_time >= 0:
                        exp_factor = np.exp(-recovery_time / self.time_constant_ps)
                        self.current_error_tolerance = self.target_error_tolerance - \
                                                      (self.target_error_tolerance - self.shock_error_tolerance) * exp_factor
                        
                        # Stop recovery when we're close to target tolerance
                        if exp_factor < 0.01:  # 99% recovered
                            self.coupling_recovery_active = False
                    else:
                        self.current_error_tolerance = self.shock_error_tolerance
                else:
                    # No shock detected and not in recovery - use target tolerance
                    self.current_error_tolerance = self.target_error_tolerance
            
            # Legacy switch time logic (still supported)
            elif self.switch_time_ps is not None:
                if current_elapsed_time_ps < self.switch_time_ps:
                    # Before switch: use target error tolerance
                    self.current_error_tolerance = self.target_error_tolerance
                else:
                    # After switch: implement shock dampening and exponential recovery
                    ramping_time = current_elapsed_time_ps - self.switch_time_ps
                    
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
            # Add minimum timestep protection to prevent numerical instability
            # Even during shock dampening, we need a non-zero timestep for the integrator
            min_dt = PhysicalConstants.fs_to_atomic_units(0.001)  # 0.001 fs minimum
            effective_error_tolerance = max(self.current_error_tolerance, 
                                          self.target_error_tolerance * 1e-10)  # Never below 1e-10 * target
            #print(f"DEBUG: effective_error_tolerance = {effective_error_tolerance}", flush=True)

            optimal_dt = np.sqrt(effective_error_tolerance / force_max_norm)
            current_dt = self.integrator.dt
            
            # Apply minimum timestep protection
            if optimal_dt < min_dt:
                optimal_dt = min_dt
                if force_immediate_update:
                    print(f"  WARNING: Computed timestep {np.sqrt(self.current_error_tolerance / force_max_norm):.2e} a.u. below minimum", flush=True)
                    print(f"  Clamping to minimum timestep: {min_dt:.2e} a.u. ({min_dt * PhysicalConstants.TIME_PS_CONVERSION * 1000:.3f} fs)", flush=True)
            
            # Apply maximum timestep protection  
            max_dt = PhysicalConstants.fs_to_atomic_units(10.0)  # 10 fs maximum
            if optimal_dt > max_dt:
                optimal_dt = max_dt
                if force_immediate_update:
                    print(f"  WARNING: Computed timestep above maximum, clamping to {max_dt * PhysicalConstants.TIME_PS_CONVERSION * 1000:.1f} fs", flush=True)
            
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











