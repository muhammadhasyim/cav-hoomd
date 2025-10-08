"""
Advanced control systems for cavity molecular dynamics simulations.

This module provides optimal and robust control algorithms for temperature regulation
in coupled thermal systems, including:
- LQR (Linear-Quadratic Regulator) with integral action
- Kalman filtering for state estimation
- In-situ system identification
"""

from typing import Optional, Dict, Any, Tuple
import hoomd
import numpy as np
import json
from datetime import datetime
from pathlib import Path

try:
    from scipy.linalg import solve_discrete_are
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("Warning: scipy not available. LQR controller will not be functional.")

from .utils import PhysicalConstants


class ParameterEKF:
    """
    Extended Kalman Filter for online system parameter estimation.
    
    This class implements an EKF that treats system parameters as augmented states,
    allowing for robust online identification of time-varying dynamics.
    
    Augmented State Vector:
    -----------------------
    z = [x_s, x_h, a11, a12, a21, a22, b1, b2]^T  (8-dimensional)
    
    Where:
    - x_s, x_h: Physical states (signal and hot temperatures)
    - a11, a12, a21, a22: Elements of discrete-time A matrix (2x2)
    - b1, b2: Elements of discrete-time B vector (2x1)
    
    Dynamics:
    ---------
    Physical states evolve according to identified parameters:
        x[k+1] = A(θ[k]) @ x[k] + B(θ[k]) @ u[k] + w_x
    
    Parameters drift slowly (random walk model):
        θ[k+1] = θ[k] + w_θ
    
    Measurement model (only physical states are measured):
        y[k] = C @ x[k] + v
        where C = [I_2×2 | 0_2×6] extracts [x_s, x_h] from augmented state
    
    Parameters
    ----------
    dt : float
        Sampling time (ps)
    process_noise_state : float
        Process noise std for physical states (K)
    process_noise_param : float
        Process noise std for parameters (unitless)
    measurement_noise : float
        Measurement noise std (K)
    initial_covariance_state : float
        Initial uncertainty in states (K²)
    initial_covariance_param : float
        Initial uncertainty in parameters (unitless²)
    """
    
    def __init__(
        self,
        dt: float,
        process_noise_state: float = 0.1,
        process_noise_param: float = 0.01,
        measurement_noise_signal: float = 0.5,
        measurement_noise_hot: float = 0.5,
        initial_covariance_state: float = 1.0,
        initial_covariance_param: float = 0.1
    ):
        """Initialize EKF for parameter estimation."""
        self.dt = dt
        
        # Augmented state: z = [x_s, x_h, a11, a12, a21, a22, b1, b2]
        self.z_hat = np.zeros((8, 1))  # State estimate
        
        # Covariance matrix P (8x8)
        self.P = np.eye(8)
        self.P[:2, :2] *= initial_covariance_state  # State uncertainty
        self.P[2:, 2:] *= initial_covariance_param  # Parameter uncertainty
        
        # Process noise covariance Q (8x8)
        self.Q = np.eye(8)
        self.Q[:2, :2] *= process_noise_state**2  # State process noise
        self.Q[2:, 2:] *= process_noise_param**2  # Parameter drift noise
        
        # Measurement noise covariance R (2x2)
        self.R = np.diag([measurement_noise_signal**2, measurement_noise_hot**2])
        
        # Measurement matrix H (2x8): extracts [x_s, x_h] from augmented state
        self.H = np.zeros((2, 8))
        self.H[0, 0] = 1.0  # Measure x_s
        self.H[1, 1] = 1.0  # Measure x_h
        
        # Statistics
        self.num_updates = 0
        self.innovation_history = []
    
    def initialize(self, x_initial: np.ndarray, A_initial: np.ndarray, B_initial: np.ndarray):
        """
        Initialize the augmented state with initial estimates.
        
        Parameters
        ----------
        x_initial : np.ndarray, shape (2, 1)
            Initial state [x_s, x_h]
        A_initial : np.ndarray, shape (2, 2)
            Initial A matrix estimate
        B_initial : np.ndarray, shape (2, 1)
            Initial B vector estimate
        """
        self.z_hat[0:2, 0] = x_initial.flatten()
        self.z_hat[2:6, 0] = A_initial.flatten()  # [a11, a12, a21, a22]
        self.z_hat[6:8, 0] = B_initial.flatten()  # [b1, b2]
    
    def predict(self, u_current: float):
        """
        EKF prediction step: propagate state and covariance forward.
        
        Nonlinear state dynamics:
            x[k+1] = A(θ[k]) @ x[k] + B(θ[k]) @ u[k]
            θ[k+1] = θ[k]  (parameters constant between measurements)
        
        Parameters
        ----------
        u_current : float
            Current control input (bath temperature deviation)
        """
        # Extract current estimates
        x_hat = self.z_hat[0:2, 0].reshape(2, 1)
        A_hat = self.z_hat[2:6, 0].reshape(2, 2)
        B_hat = self.z_hat[6:8, 0].reshape(2, 1)
        
        # Predict next physical state (nonlinear!)
        x_hat_next = A_hat @ x_hat + B_hat * u_current
        
        # Parameters remain constant (random walk)
        # z_hat_next = [x_hat_next; A_hat; B_hat]
        self.z_hat[0:2, 0] = x_hat_next.flatten()
        # Parameters unchanged: z_hat[2:8] stays same
        
        # Compute Jacobian F of state transition w.r.t. augmented state
        F = self._compute_jacobian(x_hat, A_hat, B_hat, u_current)
        
        # Propagate covariance: P[k+1|k] = F @ P @ F^T + Q
        self.P = F @ self.P @ F.T + self.Q
    
    def update(self, y_measured: np.ndarray):
        """
        EKF measurement update step: correct state estimate with measurement.
        
        Parameters
        ----------
        y_measured : np.ndarray, shape (2, 1)
            Measured physical states [x_s_meas, x_h_meas]
        """
        # Innovation: e = y - H @ z_hat
        y_pred = self.H @ self.z_hat
        innovation = y_measured - y_pred
        
        # Innovation covariance: S = H @ P @ H^T + R
        S = self.H @ self.P @ self.H.T + self.R
        
        # Kalman gain: K = P @ H^T @ S^{-1}
        K_gain = self.P @ self.H.T @ np.linalg.inv(S)
        
        # State update: z_hat = z_hat + K @ e
        self.z_hat = self.z_hat + K_gain @ innovation
        
        # Covariance update: P = (I - K @ H) @ P
        I = np.eye(8)
        self.P = (I - K_gain @ self.H) @ self.P
        
        # Track innovation for diagnostics
        self.innovation_history.append(np.linalg.norm(innovation))
        if len(self.innovation_history) > 100:
            self.innovation_history.pop(0)
        
        self.num_updates += 1
    
    def _compute_jacobian(
        self,
        x: np.ndarray,
        A: np.ndarray,
        B: np.ndarray,
        u: float
    ) -> np.ndarray:
        """
        Compute Jacobian F of augmented state dynamics.
        
        F = ∂f/∂z where f(z, u) is the nonlinear state transition function.
        
        Structure (8x8):
            ┌────────────┬─────────────────────────────┬──────────┐
            │  ∂x'/∂x    │         ∂x'/∂A             │  ∂x'/∂B  │
            │   (2x2)    │         (2x4)               │  (2x2)   │
            ├────────────┼─────────────────────────────┼──────────┤
            │   0        │          I_4                │    0     │
            │  (4x2)     │         (4x4)               │  (4x2)   │
            ├────────────┼─────────────────────────────┼──────────┤
            │   0        │           0                 │   I_2    │
            │  (2x2)     │         (2x4)               │  (2x2)   │
            └────────────┴─────────────────────────────┴──────────┘
        
        Returns
        -------
        F : np.ndarray, shape (8, 8)
            Jacobian matrix
        """
        F = np.eye(8)  # Start with identity (parameters don't change)
        
        # Top-left block: ∂x'/∂x = A
        F[0:2, 0:2] = A
        
        # Top-middle block: ∂x'/∂A = [x^T ⊗ I_2]
        # For x' = A @ x + B @ u, we have ∂x'/∂A = [x^T, 0; 0, x^T]
        x_flat = x.flatten()
        F[0, 2:4] = x_flat  # ∂x'_s/∂[a11, a12]
        F[1, 4:6] = x_flat  # ∂x'_h/∂[a21, a22]
        
        # Top-right block: ∂x'/∂B = u * I_2
        F[0:2, 6:8] = u * np.eye(2)
        
        # Lower blocks remain identity (parameters unchanged)
        
        return F
    
    def get_state_estimate(self) -> np.ndarray:
        """
        Get estimated physical states [x_s, x_h].
        
        Returns
        -------
        x_hat : np.ndarray, shape (2, 1)
            State estimate
        """
        return self.z_hat[0:2, :].copy()
    
    def get_parameter_estimates(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get estimated system matrices A and B.
        
        Returns
        -------
        A_hat : np.ndarray, shape (2, 2)
            Estimated A matrix
        B_hat : np.ndarray, shape (2, 1)
            Estimated B vector
        """
        A_hat = self.z_hat[2:6, 0].reshape(2, 2)
        B_hat = self.z_hat[6:8, 0].reshape(2, 1)
        return A_hat.copy(), B_hat.copy()
    
    def get_parameter_uncertainty(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get parameter uncertainties (standard deviations).
        
        Returns
        -------
        A_std : np.ndarray, shape (2, 2)
            Standard deviations of A elements
        B_std : np.ndarray, shape (2, 1)
            Standard deviations of B elements
        """
        A_var = np.diag(self.P[2:6, 2:6])
        B_var = np.diag(self.P[6:8, 6:8])
        
        A_std = np.sqrt(np.maximum(A_var, 0)).reshape(2, 2)
        B_std = np.sqrt(np.maximum(B_var, 0)).reshape(2, 1)
        
        return A_std, B_std
    
    def is_stable(self) -> bool:
        """
        Check if estimated A matrix is stable (eigenvalues < 1).
        
        Returns
        -------
        stable : bool
            True if all eigenvalues have magnitude < 0.99
        """
        A_hat, _ = self.get_parameter_estimates()
        try:
            eigs = np.linalg.eigvals(A_hat)
            return np.max(np.abs(eigs)) < 0.99
        except:
            return False
    
    def get_innovation_statistics(self) -> Dict[str, float]:
        """
        Get statistics on innovation magnitude for diagnostics.
        
        Returns
        -------
        stats : dict
            Dictionary with 'mean' and 'std' of recent innovations
        """
        if len(self.innovation_history) < 5:
            return {'mean': 0.0, 'std': 0.0}
        
        innovations = np.array(self.innovation_history)
        return {
            'mean': float(np.mean(innovations)),
            'std': float(np.std(innovations))
        }


class LQRTemperatureController(hoomd.custom.Action):
    """
    LQR-based optimal temperature controller with Kalman filter and integral action.
    
    This controller implements Linear-Quadratic Regulator (LQR) optimal control
    for coupled thermal subsystems in molecular dynamics simulations. It provides:
    
    1. **Zero steady-state error** via integral action (Internal Model Principle)
    2. **Optimal variance reduction** for stochastic thermal noise
    3. **In-situ system identification** using cold quench experiments
    4. **Kalman filter** for robust state estimation from noisy measurements
    
    The controller models three coupled temperature deviations from target:
    - x_s: Signal temperature (e.g., LJ+Coulombic) - regulated output
    - x_h: Hot temperature (e.g., harmonic) - disturbance signal
    - x_c: Bath temperature - actuator
    - ξ: Integral of signal error - ensures zero drift
    
    Theory
    ------
    Linear state-space model (discrete-time):
        x[k+1] = A @ x[k] + B @ u[k] + w[k]
        y[k] = C @ x[k] + v[k]
    
    Control law:
        u[k] = -K @ x_hat[k] - k_I @ ξ_hat[k]
    
    Where:
    - A, B: System matrices (identified from cold quench)
    - K, k_I: LQR optimal gains
    - x_hat, ξ_hat: Kalman filter state estimates
    - w, v: Process and measurement noise
    
    Guarantees:
    - Zero DC gain from deterministic disturbances to regulated output
    - Minimal variance for stochastic noise (set by LQR/Kalman weights)
    - Bounded control effort (configurable)
    
    Parameters
    ----------
    signal_temperature_method : str
        Temperature method for regulated output ('lj_coulombic', 'harmonic_equipartition', etc.)
    hot_temperature_method : str
        Temperature method for disturbance signal ('harmonic', 'kinetic', etc.)
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing
    energy_tracker : EnergyTracker
        Energy tracker for temperature calculations
    simulation : hoomd.Simulation
        HOOMD simulation object
    molecular_thermostat : hoomd.md.methods.thermostats.Thermostat
        Molecular thermostat to control
    cavity_thermostat : hoomd.md.methods.thermostats.Thermostat
        Cavity thermostat to control
    target_temperature : float, optional
        Fixed target temperature in Kelvin (used if dynamic_target=False). Default: 300.0
    dynamic_target : bool, optional
        If True, set target from signal at turn_on_time. Default: True
    dynamic_target_method : str, optional
        Temperature method for dynamic target (None = use signal_temperature_method)
    lqr_weight_signal : float, optional
        LQR weight for signal regulation (higher = tighter). Default: 100.0
    lqr_weight_hot : float, optional
        LQR weight for hot temperature (higher = tighter control). Default: 80.0
    lqr_weight_bath : float, optional
        LQR weight for bath temperature. Default: 0.1
    lqr_weight_integral : float, optional
        LQR weight for integral action. Default: 10.0
    lqr_control_effort : float, optional
        LQR penalty on control effort (higher = gentler control). Default: 1.0
    process_noise_signal : float, optional
        Process noise standard deviation for signal (K). Default: 0.1
    process_noise_hot : float, optional
        Process noise standard deviation for hot (K). Default: 0.1
    measurement_noise_signal : float, optional
        Measurement noise standard deviation for signal (K). Default: 0.5
    measurement_noise_hot : float, optional
        Measurement noise standard deviation for hot (K). Default: 0.5
    system_id_mode : str, optional
        System identification mode: 'step' (cold quench), 'load' (from file), 'skip' (manual). Default: 'step'
    system_id_temp_K : float, optional
        Quench temperature for system ID (K). Default: 5.0
    system_id_duration_ps : float, optional
        Duration of system ID cold quench (ps). Default: 50.0
    system_id_file : str, optional
        File to save/load system parameters. Default: 'lqr_system_params.json'
    system_params : dict, optional
        Manual system parameters (for system_id_mode='skip')
    turn_on_time_ps : float, optional
        Time to activate controller (ps). Default: 0.0
    turn_off_time_ps : float, optional
        Time to deactivate controller (ps). Default: None (never)
    update_interval_ps : float, optional
        Control update interval (ps). Default: 0.1
    T_min : float, optional
        Minimum bath temperature (K). Default: 0.1
    T_max : float, optional
        Maximum bath temperature (K). Default: None
    apply_to : str, optional
        Which thermostats to control ('molecular', 'cavity', 'both'). Default: 'both'
    output_file : str, optional
        Output CSV file for control data. Default: 'lqr_controller.csv'
    empirical_data_file : str, optional
        Empirical temperature calibration data file
    console_output_period_ps : float, optional
        Console output period (ps). Default: 1.0
    
    Examples
    --------
    **Basic LQR control after DiffEq controller:**
    
    >>> lqr_controller = LQRTemperatureController(
    ...     signal_temperature_method='lj_coulombic',
    ...     hot_temperature_method='harmonic_equipartition',
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     simulation=sim,
    ...     molecular_thermostat=molecular_thermo,
    ...     cavity_thermostat=cavity_thermo,
    ...     dynamic_target=True,
    ...     turn_on_time_ps=510.0,
    ...     system_id_mode='step',
    ...     system_id_temp_K=5.0,
    ...     system_id_duration_ps=50.0,
    ...     lqr_weight_signal=100.0,
    ...     lqr_weight_hot=80.0,
    ...     apply_to='both'
    ... )
    
    **Load pre-identified parameters:**
    
    >>> lqr_controller = LQRTemperatureController(
    ...     signal_temperature_method='lj_coulombic',
    ...     hot_temperature_method='harmonic_equipartition',
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     simulation=sim,
    ...     molecular_thermostat=molecular_thermo,
    ...     turn_on_time_ps=510.0,
    ...     system_id_mode='load',
    ...     system_id_file='lqr_system_params_nocoupling.json',
    ...     apply_to='molecular'
    ... )
    
    Notes
    -----
    - Designed for sequential control after DiffEqController with cavity coupling OFF
    - System identification runs at turn_on_time before control begins
    - Requires scipy for LQR Riccati equation solver
    - Controller operates in discrete time with sample period = update_interval_ps
    - Integral action ensures zero mean tracking error despite constant disturbances
    - Kalman filter provides optimal state estimates from noisy measurements
    
    References
    ----------
    [1] Åström, K. J., & Murray, R. M. (2021). Feedback Systems (2nd ed.). Princeton University Press.
    [2] Anderson, B. D., & Moore, J. B. (2007). Optimal Control: Linear Quadratic Methods. Dover.
    """
    
    def __init__(self,
                 signal_temperature_method: str,
                 hot_temperature_method: str,
                 time_tracker,
                 energy_tracker,
                 simulation=None,
                 molecular_thermostat=None,
                 cavity_thermostat=None,
                 temperature_tracker=None,
                 target_temperature: float = 300.0,
                 dynamic_target: bool = True,
                 dynamic_target_method: Optional[str] = None,
                 lqr_weight_signal: float = 100.0,
                 lqr_weight_hot: float = 80.0,
                 lqr_weight_bath: float = 0.1,
                 lqr_weight_integral: float = 10.0,
                 lqr_control_effort: float = 1.0,
                 process_noise_signal: float = 0.1,
                 process_noise_hot: float = 0.1,
                 measurement_noise_signal: float = 0.5,
                 measurement_noise_hot: float = 0.5,
                 system_id_mode: str = 'step',
                 system_id_temp_K: float = 5.0,
                 system_id_duration_ps: float = 50.0,
                 system_id_file: str = 'lqr_system_params.json',
                 system_params: Optional[Dict] = None,
                 periodic_system_id: bool = False,
                 periodic_system_id_interval_ps: float = 1000.0,
                 turn_on_time_ps: float = 0.0,
                 turn_off_time_ps: Optional[float] = None,
                 update_interval_ps: float = 0.1,
                 T_min: float = 0.1,
                 T_max: Optional[float] = None,
                 apply_to: str = 'both',
                 output_file: str = 'lqr_controller.csv',
                 empirical_data_file: Optional[str] = None,
                 console_output_period_ps: float = 1.0,
                 use_ekf_adaptation: bool = True,
                 ekf_update_interval: int = 50,
                 ekf_process_noise_param: float = 0.001,
                 ekf_initial_covariance_param: float = 0.1,
                 adaptive_lqr_threshold: float = 0.05,
                 enable_gain_scheduling: bool = True,
                 gain_schedule_far_threshold: float = 20.0,
                 gain_schedule_near_threshold: float = 10.0,
                 th_filter_enabled: bool = True,
                 th_filter_time_constant: float = 20.0,
                 # Gentle startup parameters
                 gentle_startup_steps: int = 10,
                 gentle_startup_min_authority: float = 0.1,
                # Kinetic temperature tracking (3D state augmentation)
                track_kinetic_temp: bool = False,
                weight_kinetic: float = 100.0,
                process_noise_kinetic: float = 2.0,
                measurement_noise_kinetic: float = 2.0,
                # Cross-coupling weights for thermal equilibration
                cross_coupling_signal_kinetic: float = 0.0,
                cross_coupling_signal_hot: float = 0.0,
                cross_coupling_hot_kinetic: float = 0.0):
        
        if not HAS_SCIPY:
            raise ImportError("LQR controller requires scipy. Please install: pip install scipy")
        
        # Store references
        self.signal_temperature_method = signal_temperature_method
        self.hot_temperature_method = hot_temperature_method
        self.time_tracker = time_tracker
        self.energy_tracker = energy_tracker
        self.temperature_tracker = temperature_tracker
        self.simulation = simulation
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        
        # Target settings
        self.target_temperature = target_temperature
        self.dynamic_target = dynamic_target
        self.dynamic_target_method = dynamic_target_method or signal_temperature_method
        self.target_set = False
        
        # LQR weights (LQG formulation: 2D state [x_s, x_h], 1D control [u_bath])
        self.lqr_weight_signal = lqr_weight_signal
        self.lqr_weight_hot = lqr_weight_hot
        self.lqr_weight_bath = lqr_weight_bath  # DEPRECATED: not used in LQG 2D formulation
        self.lqr_weight_integral = lqr_weight_integral
        self.lqr_control_effort = lqr_control_effort
        
        # Warn if deprecated parameter is being used
        if lqr_weight_bath != 0.1:  # Check if user changed from default
            print(f"⚠ WARNING: lqr_weight_bath={lqr_weight_bath} is DEPRECATED in LQG 2D formulation")
            print(f"  Bath temperature is now the control input, not a state variable.")
            print(f"  This parameter is ignored. Use lqr_control_effort to penalize control action.")
        
        # Kalman filter noise parameters
        self.process_noise_signal = process_noise_signal
        self.process_noise_hot = process_noise_hot
        self.measurement_noise_signal = measurement_noise_signal
        self.measurement_noise_hot = measurement_noise_hot
        
        # System identification
        self.system_id_mode = system_id_mode
        self.system_id_temp_K = system_id_temp_K
        self.system_id_duration_ps = system_id_duration_ps
        self.system_id_file = system_id_file
        self.system_params = system_params
        
        # Periodic system ID for adaptive re-identification
        self.periodic_system_id = periodic_system_id
        self.periodic_system_id_interval_ps = periodic_system_id_interval_ps
        self.last_system_id_completion_time = 0.0
        self.system_id_count = 0  # Track number of system IDs performed
        
        # Timing
        self.turn_on_time_ps = turn_on_time_ps
        self.turn_off_time_ps = turn_off_time_ps
        self.update_interval_ps = update_interval_ps
        self.dt = update_interval_ps  # Sample time
        
        # Temperature limits
        self.T_min = T_min
        self.T_max = T_max
        
        # Control application
        self.apply_to = apply_to
        
        # Output
        self.output_file = output_file
        self.console_output_period_ps = console_output_period_ps
        
        # Empirical data for temperature calculations
        self.empirical_data = None
        if empirical_data_file:
            from .analysis import EmpiricalTemperatureData
            self.empirical_data = EmpiricalTemperatureData(empirical_data_file)
        
        # Controller state
        self.phase = 'dormant'  # 'dormant' -> 'system_id' -> 'control'
        self.is_active = False
        self.last_update_time = -float('inf')
        self.last_console_output_time = 0.0
        self.num_updates = 0  # Counter for control updates
        
        # System ID data collection
        self.system_id_started = False
        self.system_id_start_time = 0.0
        self.system_id_data = {
            'time': [],
            'x_s': [],
            'x_h': [],
            'x_kin': [],  # Kinetic temperature (only used if track_kinetic_temp=True)
            'x_c': [],
            'u': []
        }
        
        # System matrices (to be identified or loaded)
        # LQG formulation: x = [x_s, x_h] (2D state), u = x_c (1D control input)
        self.A_d = None  # Discrete-time state matrix (state_dim x state_dim)
        self.B_d = None  # Discrete-time input matrix (state_dim x 1)
        self.C = None  # Measurement matrix (state_dim x state_dim) - will be set after track_kinetic_temp is known
        
        # Control gains (to be computed)
        self.K = None  # LQR state feedback gain (1 x state_dim+1): acts on [x_s, x_h, ξ] or [x_s, x_h, x_kin, ξ]
        self.L = None  # Kalman filter gain (state_dim x state_dim): estimates [x_s, x_h] or [x_s, x_h, x_kin]
        
        # Kinetic temperature tracking (3D state augmentation) - must be set before state init
        self.track_kinetic_temp = track_kinetic_temp  # Enable kinetic temp as 3rd state
        self.weight_kinetic = weight_kinetic  # LQR weight for kinetic temperature
        self.process_noise_kinetic = process_noise_kinetic  # Kalman process noise for T_kin
        self.measurement_noise_kinetic = measurement_noise_kinetic  # Kalman measurement noise for T_kin
        self.T_kin_filtered = None  # Low-pass filtered kinetic temperature
        
        # Cross-coupling weights for thermal equilibration
        # Penalties on pairwise temperature differences (T_i - T_j)²
        self.cross_coupling_signal_kinetic = cross_coupling_signal_kinetic
        self.cross_coupling_signal_hot = cross_coupling_signal_hot
        self.cross_coupling_hot_kinetic = cross_coupling_hot_kinetic
        
        # Determine state dimension based on tracking options
        # 2D: [x_s, x_h], 3D: [x_s, x_h, x_kin]
        if track_kinetic_temp:
            self.state_dim = 3  # 3D: [x_s, x_h, x_kin]
        else:
            self.state_dim = 2  # 2D: [x_s, x_h]
        
        # Set measurement matrix based on state dimension
        if self.state_dim == 2:
            self.C = np.array([[1.0, 0.0],  # Measure x_s
                              [0.0, 1.0]])  # Measure x_h
        else:  # 3D state
            self.C = np.array([[1.0, 0.0, 0.0],  # Measure x_s
                              [0.0, 1.0, 0.0],  # Measure x_h
                              [0.0, 0.0, 1.0]])  # Measure x_kin
        
        # State estimates (size depends on state_dim)
        self.x_hat = np.zeros((self.state_dim, 1))  # [x_s, x_h] or [x_s, x_h, x_kin]
        self.xi_hat = np.zeros((1, 1))  # Integral state: ξ = ∫(x_s)dt
        
        # Low-pass filter for T_h measurement noise suppression
        self.T_h_filtered = None  # Exponential moving average of hot temperature
        self.T_h_filter_enabled = th_filter_enabled  # Enable low-pass filtering on T_h
        self.T_h_filter_time_constant = th_filter_time_constant  # ps (filters out oscillations faster than ~20ps)
        
        # Adaptive control: online parameter estimation with EKF or RLS
        self.use_ekf_adaptation = use_ekf_adaptation  # Use EKF (True) or RLS (False)
        self.enable_adaptation = True  # Enable adaptive LQR
        self.adaptation_window = 200  # Number of samples for RLS window
        self.adaptation_interval = 50  # Update parameters every N samples
        self.adaptation_counter = 0
        self.rls_data_buffer = {
            'x': [],  # State vectors [x_s, x_h]
            'u': [],  # Control inputs
            'x_next': []  # Next states
        }
        self.rls_regularization = 1e-6  # Ridge regression for numerical stability
        
        # EKF-based adaptation parameters
        self.ekf_update_interval = ekf_update_interval  # Update EKF every N control steps
        self.ekf_process_noise_param = ekf_process_noise_param  # Parameter drift noise
        self.ekf_initial_covariance_param = ekf_initial_covariance_param  # Initial parameter uncertainty
        self.parameter_ekf = None  # Will be initialized after system ID
        
        # Adaptive LQR redesign
        self.adaptive_lqr_threshold = adaptive_lqr_threshold  # Re-design LQR if ||θ_new - θ_old||/||θ_old|| > threshold
        self.last_lqr_redesign_params = None  # Track when LQR was last re-designed
        self.num_lqr_redesigns = 0  # Counter for LQR redesigns
        
        # Gain scheduling
        self.enable_gain_scheduling = enable_gain_scheduling  # Enable gain scheduling
        self.gain_schedule_far_threshold = gain_schedule_far_threshold  # |T_h - T_s| > this → full gain
        self.gain_schedule_near_threshold = gain_schedule_near_threshold  # |T_h - T_s| < this → reduced gain
        self.gain_schedule_alpha = 1.0  # Current gain scheduling factor
        
        # Gentle startup parameters
        self.gentle_startup_steps = gentle_startup_steps  # Number of steps to ramp up (default: 10)
        self.gentle_startup_min_authority = gentle_startup_min_authority  # Starting authority fraction (default: 0.1 = 10%)
        
        # Adaptive noise estimation
        self.noise_estimation_window = 100  # Samples for variance estimation
        self.noise_buffer_signal = []
        self.noise_buffer_hot = []
        self.estimated_noise_signal = measurement_noise_signal
        self.estimated_noise_hot = measurement_noise_hot
        self.last_noise_update = 0.0
        self.noise_update_interval = 10.0  # Update noise estimate every 10ps
        
        # Median filter for outlier rejection
        self.median_filter_window = 5  # Use median of last 5 measurements
        self.temp_buffer_signal = []
        self.temp_buffer_hot = []
        
        # Adaptation statistics
        self.num_adaptations = 0
        self.last_r_squared = [0.0, 0.0]  # Track model quality (2D state system)
        
        # Control rate limiting
        self.max_control_rate = 10.0  # Maximum K/ps change (safety)
        self.u_prev_temp = None  # Previous bath temperature for rate limiting
        
        # Integral anti-windup
        self.integral_max = 1000.0  # Maximum integral accumulation
        
        # Load parameters if requested
        if self.system_id_mode == 'load':
            self._load_system_parameters()
            self._design_lqr_controller()
            self.phase = 'control'
        elif self.system_id_mode == 'skip':
            if self.system_params is None:
                raise ValueError("system_id_mode='skip' requires system_params dict")
            self._set_manual_parameters()
            self._design_lqr_controller()
            self.phase = 'control'
        elif self.system_id_mode == 'none':
            # No quench: initialize with physics-based guesses and learn online
            print("\n" + "="*70)
            print("ADAPTIVE NO-QUENCH MODE: Initializing with physics-based guesses")
            print("="*70)
            self._initialize_from_physics()
            self._design_lqr_controller()
            self.phase = 'control'
            print("Starting control immediately - RLS will learn online from control actions")
            print("="*70 + "\n")
        else:  # 'step'
            self.phase = 'system_id'
        
        # Initialize output file
        self._initialize_output_file()
    
    def _initialize_output_file(self):
        """Initialize CSV output file with header."""
        try:
            with open(self.output_file, 'w', encoding='utf-8') as f:
                f.write("# LQR Temperature Controller Output\n")
                f.write(f"# Signal method: {self.signal_temperature_method}\n")
                f.write(f"# Hot method: {self.hot_temperature_method}\n")
                f.write(f"# Target temperature: {self.target_temperature} K\n")
                f.write(f"# System ID mode: {self.system_id_mode}\n")
                f.write("time_ps,phase,signal_temp_K,hot_temp_K,bath_temp_K,"
                       "x_s_hat,x_h_hat,x_c_hat,xi_hat,control_u,is_active\n")
        except Exception as e:
            print(f"Warning: Could not initialize LQR controller output file: {e}")
    
    def _get_hoomd_simulation(self):
        """Get HOOMD simulation object."""
        if self.simulation is not None:
            return self.simulation
        
        # Try to get from global state
        try:
            import hoomd
            return hoomd.simulation.Simulation.get_current()
        except:
            return None
    
    def act(self, timestep):
        """Execute LQR controller at each timestep."""
        try:
            current_time_ps = self.time_tracker.elapsed_time
            
            # Check if we should be active
            should_be_active = (current_time_ps >= self.turn_on_time_ps and 
                              (self.turn_off_time_ps is None or current_time_ps < self.turn_off_time_ps))
            
            if not should_be_active:
                # Controller is dormant
                if self.is_active:
                    self.is_active = False
                    print(f"LQR controller turned OFF at t = {current_time_ps:.2f} ps")
                return
            
            # Controller should be active
            if not self.is_active:
                self.is_active = True
                print(f"\n{'='*70}")
                print(f"LQR Temperature Controller activated at t = {current_time_ps:.2f} ps")
                print(f"{'='*70}")
                print(f"Signal method: {self.signal_temperature_method}")
                print(f"Hot method: {self.hot_temperature_method}")
                print(f"System ID mode: {self.system_id_mode}")
                
                # Set dynamic target if requested
                if self.dynamic_target and not self.target_set:
                    measured_temp = self._measure_temperature(self.dynamic_target_method, current_time_ps)
                    if measured_temp is not None:
                        self.target_temperature = measured_temp
                        self.target_set = True
                        print(f"Dynamic target set: T_target = {self.target_temperature:.2f} K")
                    else:
                        print(f"Warning: Could not measure dynamic target, using {self.target_temperature:.2f} K")
            
            # Execute appropriate phase
            if self.phase == 'system_id':
                self._execute_system_identification(timestep, current_time_ps)
            elif self.phase == 'control':
                # Check if periodic system ID is needed
                if self.periodic_system_id:
                    time_since_last_id = current_time_ps - self.last_system_id_completion_time
                    if time_since_last_id >= self.periodic_system_id_interval_ps:
                        self._trigger_periodic_system_id(current_time_ps)
                        return  # Will execute system ID on next call
                
                self._execute_lqr_control(timestep, current_time_ps)
        
        except Exception as e:
            print(f"Error in LQR controller: {e}")
            import traceback
            traceback.print_exc()
    
    def _execute_system_identification(self, timestep, current_time_ps):
        """Execute system identification phase with single or multi-step excitation.
        
        Modes:
        - 'step': Single cold quench (traditional)
        - 'multi_step': Multiple sequential steps for better frequency coverage
        
        IMPORTANT: For adaptive timestep, this is called at variable intervals.
        We collect data EVERY call (not time-filtered) to get all available samples.
        """
        # Initialize system ID on first call
        if not self.system_id_started:
            self.system_id_started = True
            self.system_id_start_time = current_time_ps
            self.system_id_sample_count = 0
            self.system_id_current_step = 0  # Track which excitation step we're on
            
            # SAVE current bath temperature before system ID
            # Excitations will be relative to THIS, not the target
            current_bath_temp = self._get_current_bath_temperature()
            if current_bath_temp is None:
                print(f"⚠ WARNING: Could not read current bath temperature, using target")
                current_bath_temp = self.target_temperature
            self.bath_temp_before_system_id = current_bath_temp
            
            print(f"\n{'='*70}")
            print(f"System Identification Phase Started")
            print(f"{'='*70}")
            print(f"True target temperature: {self.target_temperature:.1f} K (will be restored later)")
            print(f"Current bath temperature: {current_bath_temp:.1f} K (excitations relative to this)")
            print(f"Mode: {self.system_id_mode}")
            print(f"Duration: {self.system_id_duration_ps} ps")
            print(f"NOTE: Using adaptive timestep - collecting ALL samples regardless of dt")
            
            # Define excitation sequence based on mode
            if self.system_id_mode == 'multi_step':
                # Multi-step excitation: 4 steps for comprehensive frequency coverage
                # Each step is 1/4 of total duration
                step_duration = self.system_id_duration_ps / 4.0
                
                # Define 4 excitation levels (relative to CURRENT bath temp):
                # Step 1: Down (cooling)
                # Step 2: Up (heating)  
                # Step 3: Down again (different magnitude)
                # Step 4: Return to current
                excitation_magnitudes = [
                    -self.system_id_temp_K,      # -5K: cool down
                    +self.system_id_temp_K * 0.7, # +3.5K: heat up
                    -self.system_id_temp_K * 0.5, # -2.5K: cool slightly
                    0.0                           # 0K: return to current
                ]
                
                self.system_id_steps = [
                    {
                        'start_time': 0.0,
                        'duration': step_duration,
                        'excitation': excitation_magnitudes[0],
                        'temperature': current_bath_temp + excitation_magnitudes[0]
                    },
                    {
                        'start_time': step_duration,
                        'duration': step_duration,
                        'excitation': excitation_magnitudes[1],
                        'temperature': current_bath_temp + excitation_magnitudes[1]
                    },
                    {
                        'start_time': 2 * step_duration,
                        'duration': step_duration,
                        'excitation': excitation_magnitudes[2],
                        'temperature': current_bath_temp + excitation_magnitudes[2]
                    },
                    {
                        'start_time': 3 * step_duration,
                        'duration': step_duration,
                        'excitation': excitation_magnitudes[3],
                        'temperature': current_bath_temp + excitation_magnitudes[3]
                    }
                ]
                
                print(f"\nMulti-Step Excitation Sequence (relative to current bath {current_bath_temp:.1f}K):")
                print(f"  Step 1 (0-{step_duration:.1f}ps):     T = {self.system_id_steps[0]['temperature']:.1f}K (Δ={excitation_magnitudes[0]:+.1f}K)")
                print(f"  Step 2 ({step_duration:.1f}-{2*step_duration:.1f}ps): T = {self.system_id_steps[1]['temperature']:.1f}K (Δ={excitation_magnitudes[1]:+.1f}K)")
                print(f"  Step 3 ({2*step_duration:.1f}-{3*step_duration:.1f}ps): T = {self.system_id_steps[2]['temperature']:.1f}K (Δ={excitation_magnitudes[2]:+.1f}K)")
                print(f"  Step 4 ({3*step_duration:.1f}-{self.system_id_duration_ps:.1f}ps): T = {self.system_id_steps[3]['temperature']:.1f}K (Δ={excitation_magnitudes[3]:+.1f}K)")
                print(f"  → 4 excitations, both heating & cooling, returning to current bath temp")
                
            else:  # 'step' mode (traditional single quench)
                desired_quench_temp = current_bath_temp - self.system_id_temp_K
                actual_quench_temp = max(self.T_min * 2.0, desired_quench_temp)
                actual_quench_magnitude = current_bath_temp - actual_quench_temp
                
                print(f"Requested quench magnitude: {self.system_id_temp_K:.1f} K (relative to current)")
                if actual_quench_temp != desired_quench_temp:
                    print(f"⚠ WARNING: Quench limited to prevent negative temperature!")
                    print(f"  Desired: {desired_quench_temp:.1f} K → Limited to: {actual_quench_temp:.1f} K")
                print(f"Applying cold quench: T_bath → {actual_quench_temp:.1f} K (Δ = -{actual_quench_magnitude:.1f} K)")
                
                self.system_id_steps = [
                    {
                        'start_time': 0.0,
                        'duration': self.system_id_duration_ps,
                        'excitation': -(actual_quench_magnitude),
                        'temperature': actual_quench_temp
                    }
                ]
            
            # Apply initial excitation
            print(f"\nApplying initial excitation: T_bath = {self.system_id_steps[0]['temperature']:.1f} K")
            self._update_thermostats(self.system_id_steps[0]['temperature'])
        
        # Collect data EVERY call (critical for adaptive timestep!)
        elapsed_time = current_time_ps - self.system_id_start_time
        
        # Check if we need to transition to next excitation step (multi_step mode only)
        if self.system_id_mode == 'multi_step':
            for i in range(self.system_id_current_step + 1, len(self.system_id_steps)):
                step = self.system_id_steps[i]
                if elapsed_time >= step['start_time']:
                    # Transition to next step
                    self.system_id_current_step = i
                    print(f"\n[System ID] t={current_time_ps:.1f}ps: Transitioning to Step {i+1}")
                    print(f"  Setting T_bath = {step['temperature']:.1f}K (Δ={step['excitation']:+.1f}K)")
                    self._update_thermostats(step['temperature'])
                    break
        
        # Measure temperatures (raw)
        signal_temp_raw = self._measure_temperature(self.signal_temperature_method, current_time_ps)
        hot_temp_raw = self._measure_temperature(self.hot_temperature_method, current_time_ps)
        kinetic_temp_raw = self._measure_temperature('kinetic', current_time_ps) if self.track_kinetic_temp else None
        bath_temp = self._get_current_bath_temperature()
        
        # Check if we have all required measurements
        required_temps_valid = (signal_temp_raw is not None and hot_temp_raw is not None and bath_temp is not None)
        if self.track_kinetic_temp:
            required_temps_valid = required_temps_valid and (kinetic_temp_raw is not None)
        
        if required_temps_valid:
            # Apply median filter for outlier rejection during system ID
            signal_temp = self._apply_median_filter(signal_temp_raw, self.temp_buffer_signal, self.median_filter_window)
            hot_temp = self._apply_median_filter(hot_temp_raw, self.temp_buffer_hot, self.median_filter_window)
            if self.track_kinetic_temp:
                kinetic_temp = self._apply_median_filter(kinetic_temp_raw, self.temp_buffer_signal, self.median_filter_window)
            else:
                kinetic_temp = None
            
            # Store deviations from target
            # Get current excitation magnitude
            current_step = self.system_id_steps[self.system_id_current_step]
            control_deviation = current_step['excitation']
            
            self.system_id_data['time'].append(elapsed_time)
            self.system_id_data['x_s'].append(signal_temp - self.target_temperature)
            self.system_id_data['x_h'].append(hot_temp - self.target_temperature)
            if self.track_kinetic_temp:
                self.system_id_data['x_kin'].append(kinetic_temp - self.target_temperature)
            else:
                self.system_id_data['x_kin'].append(0.0)  # Placeholder for 2D mode
            self.system_id_data['x_c'].append(bath_temp - self.target_temperature)
            self.system_id_data['u'].append(control_deviation)
            self.system_id_sample_count += 1
            
            # Log to file
            try:
                with open(self.output_file, 'a', encoding='utf-8') as f:
                    f.write(f"{current_time_ps:.6f},system_id,{signal_temp:.6f},{hot_temp:.6f},"
                           f"{bath_temp:.6f},0.0,0.0,0.0,0.0,"
                           f"{control_deviation:.6f},1\n")
            except:
                pass
            
            # Progress update every 100 samples
            if self.system_id_sample_count % 100 == 0:
                avg_dt = elapsed_time / self.system_id_sample_count if self.system_id_sample_count > 0 else 0
                print(f"  System ID: {self.system_id_sample_count} samples collected, "
                      f"t={current_time_ps:.1f}ps, avg Δt={avg_dt:.3f}ps")
        
        # Check if identification complete
        if elapsed_time >= self.system_id_duration_ps:
            print(f"\nSystem ID data collection complete at t = {current_time_ps:.2f} ps")
            print(f"Collected {len(self.system_id_data['time'])} data points")
            if len(self.system_id_data['time']) > 1:
                times = np.array(self.system_id_data['time'])
                dts = np.diff(times)
                print(f"  Sampling statistics:")
                print(f"    Mean Δt: {np.mean(dts):.3f} ps")
                print(f"    Median Δt: {np.median(dts):.3f} ps")
                print(f"    Min Δt: {np.min(dts):.3f} ps")
                print(f"    Max Δt: {np.max(dts):.3f} ps")
            
            # Fit system parameters
            self._fit_system_parameters()
            
            # Design LQR controller
            self._design_lqr_controller()
            
            # Record system ID completion
            self.last_system_id_completion_time = current_time_ps
            self.system_id_count += 1
            
            # RESTORE bath temperature to pre-system-ID value
            # This ensures smooth transition back to control from wherever we left off
            if hasattr(self, 'bath_temp_before_system_id'):
                restore_temp = self.bath_temp_before_system_id
                print(f"\n{'='*70}")
                print(f"Restoring Bath Temperature")
                print(f"{'='*70}")
                print(f"Restoring bath to pre-system-ID value: {restore_temp:.1f} K")
                print(f"True target temperature (unchanged): {self.target_temperature:.1f} K")
                print(f"Controller will regulate toward target from this starting point.")
                self._update_thermostats(restore_temp)
            
            # Transition to control phase
            self.phase = 'control'
            self.num_updates = 0  # Reset gentle startup counter
            
            # Initialize state estimates from CURRENT measurements for smooth handoff
            # This prevents aggressive control from assuming we're at target when we're not
            signal_temp_now = self._measure_temperature(self.signal_temperature_method, current_time_ps)
            hot_temp_now = self._measure_temperature(self.hot_temperature_method, current_time_ps)
            kinetic_temp_now = self._measure_temperature('kinetic', current_time_ps) if self.track_kinetic_temp else None
            
            # Check if we have all required measurements
            required_temps_valid = (signal_temp_now is not None and hot_temp_now is not None)
            if self.track_kinetic_temp:
                required_temps_valid = required_temps_valid and (kinetic_temp_now is not None)
            
            if required_temps_valid:
                # Initialize state estimates as deviations from target
                x_s_init = signal_temp_now - self.target_temperature
                x_h_init = hot_temp_now - self.target_temperature
                
                if self.track_kinetic_temp:
                    x_kin_init = kinetic_temp_now - self.target_temperature
                    self.x_hat = np.array([[x_s_init], [x_h_init], [x_kin_init]])
                else:
                    self.x_hat = np.array([[x_s_init], [x_h_init]])
                
                # RESET INTEGRAL STATE: Simple and robust approach
                # After system ID, we have new system matrices and gains.
                # Reset integral to zero and let the controller naturally settle.
                # Gentle startup (10% → 100% ramp) provides damping to prevent aggressive transients.
                if hasattr(self, 'bath_temp_before_system_id'):
                    T_bath_current = self.bath_temp_before_system_id
                    
                    # Reset integral state to zero (fresh start with new controller)
                    self.xi_hat = np.zeros((1, 1))
                    
                    # Reset gain scheduling to initial state
                    self.gain_schedule_alpha = 1.0
                    
                    print(f"\n{'='*70}")
                    print(f"Controller Re-initialization After System ID")
                    print(f"{'='*70}")
                    print(f"Signal temp: {signal_temp_now:.2f}K → x_s = {x_s_init:+.2f}K")
                    print(f"Hot temp: {hot_temp_now:.2f}K → x_h = {x_h_init:+.2f}K")
                    if self.track_kinetic_temp:
                        print(f"Kinetic temp: {kinetic_temp_now:.2f}K → x_kin = {x_kin_init:+.2f}K")
                    print(f"Target: {self.target_temperature:.2f}K")
                    print(f"Bath restored to: {T_bath_current:.2f}K")
                    print(f"Integral state: ξ = 0 (reset for new controller)")
                    print(f"Gentle startup will ramp control authority: {self.gentle_startup_min_authority*100:.0f}% → 100% over {self.gentle_startup_steps} steps")
                    print(f"Small transient expected as controller re-learns from ξ=0")
                    print(f"{'='*70}\n")
                else:
                    # No saved bath temp, use zero integral
                    self.xi_hat = np.zeros((1, 1))
                    print(f"\n{'='*70}")
                    print(f"Initializing State Estimates from Current Measurements")
                    print(f"{'='*70}")
                    print(f"Signal temp: {signal_temp_now:.2f}K → x_s = {x_s_init:+.2f}K")
                    print(f"Hot temp: {hot_temp_now:.2f}K → x_h = {x_h_init:+.2f}K")
                    if self.track_kinetic_temp:
                        print(f"Kinetic temp: {kinetic_temp_now:.2f}K → x_kin = {x_kin_init:+.2f}K")
                    print(f"Target: {self.target_temperature:.2f}K")
                    print(f"Integral state: ξ = 0 (no pre-system-ID bath temp available)")
                    print(f"{'='*70}\n")
            else:
                # Fallback: if measurements fail, use zeros (previous behavior)
                self.x_hat = np.zeros((self.state_dim, 1))
                self.xi_hat = np.zeros((1, 1))
                print("⚠ WARNING: Could not measure temperatures, initializing state estimates to zeros")
            self.u_prev = 0.0  # Previous control input (for Kalman prediction)
            self.u_prev_temp = None  # Reset rate limiter memory
            self.T_h_filtered = None  # Reset low-pass filter
            
            print(f"{'='*70}")
            print(f"Transitioning to LQR Control Phase")
            print(f"{'='*70}")
            print(f"LQR gains computed. Starting optimal control...")
            print(f"Kalman filter initialized with current state. Bath temp maintained.")
            if self.periodic_system_id:
                next_id_time = current_time_ps + self.periodic_system_id_interval_ps
                print(f"Next periodic system ID scheduled at t = {next_id_time:.1f} ps")
            print(f"{'='*70}\n")
    
    def _fit_system_parameters(self):
        """Fit system parameters using least-squares from collected data.
        
        LQG formulation: 
        - State x = [x_s, x_h] (2D) or [x_s, x_h, x_kin] (3D): temperature deviations
        - Control u = x_c (1D): bath temperature deviation
        - Dynamics: x[k+1] = A_d @ x[k] + B_d @ u[k]
        """
        print(f"\nFitting system parameters...")
        
        # Convert lists to arrays
        N = len(self.system_id_data['time'])
        x_s = np.array(self.system_id_data['x_s'])
        x_h = np.array(self.system_id_data['x_h'])
        x_kin = np.array(self.system_id_data['x_kin']) if self.track_kinetic_temp else None
        x_c = np.array(self.system_id_data['x_c'])
        u = np.array(self.system_id_data['u'])
        
        # CRITICAL FIX: Compute actual sample time from collected data
        # Previously used self.update_interval_ps, which could be 0.0, causing division by zero
        times = np.array(self.system_id_data['time'])
        if len(times) > 1:
            dts = np.diff(times)
            self.dt = float(np.median(dts))  # Use median for robustness against outliers
            print(f"\n{'='*70}")
            print(f"Sample Time Computation (Critical for Parameter Extraction)")
            print(f"{'='*70}")
            print(f"  Median Δt: {self.dt:.6f} ps (used for discrete→continuous conversion)")
            print(f"  Mean Δt:   {np.mean(dts):.6f} ps")
            print(f"  Std Dev:   {np.std(dts):.6f} ps")
            print(f"  Range:     [{np.min(dts):.6f}, {np.max(dts):.6f}] ps")
            print(f"  Total samples: {len(times)}")
            
            # Warn if sampling is too irregular
            dt_std = np.std(dts)
            dt_cv = dt_std / self.dt if self.dt > 0 else float('inf')
            if dt_cv > 0.1:
                print(f"\n  ⚠ WARNING: Sampling is irregular (CV={dt_cv:.2%})")
                print(f"     This may affect parameter extraction accuracy.")
            
            # Safeguard against zero or negative dt
            if self.dt <= 0:
                print(f"\n  ⚠⚠⚠ CRITICAL ERROR: Computed dt = {self.dt} ps (invalid!) ⚠⚠⚠")
                print(f"      Using fallback dt = 0.01 ps")
                self.dt = 0.01
        else:
            print(f"\n  ⚠ ERROR: Insufficient data for dt estimation (only {len(times)} samples)!")
            # Fallback: use update_interval_ps if valid, otherwise 0.01 ps
            if self.update_interval_ps > 0:
                self.dt = self.update_interval_ps
                print(f"     Using fallback dt = {self.dt} ps (from update_interval_ps)")
            else:
                self.dt = 0.01
                print(f"     Using fallback dt = {self.dt} ps (hardcoded default)")
        
        print(f"{'='*70}\n")
        
        # Build state trajectory matrix (2D or 3D depending on configuration)
        if self.track_kinetic_temp:
            x = np.vstack([x_s, x_h, x_kin])  # (3, N) - signal, hot, and kinetic
        else:
            x = np.vstack([x_s, x_h])  # (2, N) - only signal and hot
        
        # Control input is bath temperature deviation
        u_input = x_c  # (N,) - bath is control, not state!
        
        # Least squares identification: x[k+1] = A_d @ x[k] + B_d @ u[k]
        n_states = x.shape[0]  # 2 or 3
        X_curr = x[:, :-1]  # (n_states, N-1)
        X_next = x[:, 1:]   # (n_states, N-1)
        U = u_input[:-1].reshape(1, -1)  # (1, N-1)
        
        # Stack [X_curr; U] as regression matrix
        Z = np.vstack([X_curr, U])  # (n_states+1, N-1)
        
        # Solve: X_next = [A_d, B_d] @ Z
        # [A_d, B_d] = X_next @ Z^T @ (Z @ Z^T)^{-1}
        AB = X_next @ Z.T @ np.linalg.pinv(Z @ Z.T)
        
        self.A_d = AB[:, :n_states]  # (n_states, n_states) state matrix
        self.B_d = AB[:, n_states:n_states+1]  # (n_states, 1) input matrix
        
        # Compute fit quality metrics (for all states)
        X_pred = AB @ Z  # (n_states, N-1)
        residual = X_next - X_pred
        
        # 1. R² (variance explained - sensitive to noise)
        ss_res = np.sum(residual**2, axis=1)
        ss_tot = np.sum((X_next - X_next.mean(axis=1, keepdims=True))**2, axis=1)
        r_squared = 1 - ss_res / (ss_tot + 1e-10)
        
        # Store R² for display (update self.last_r_squared for 2D or 3D)
        if n_states == 2:
            self.last_r_squared = [float(r_squared[0]), float(r_squared[1])]
        else:  # n_states == 3
            # For 3D, store all three but display will only show first two for compatibility
            if not hasattr(self, 'last_r_squared_3d'):
                self.last_r_squared_3d = [0.0, 0.0, 0.0]
            self.last_r_squared_3d = [float(r_squared[0]), float(r_squared[1]), float(r_squared[2])]
            self.last_r_squared = [float(r_squared[0]), float(r_squared[1])]  # Keep 2D for compatibility
        
        # 2. MAE (Mean Absolute Error - robust to outliers)
        mae = np.mean(np.abs(residual), axis=1)
        
        # 3. Smoothed R² (apply moving average filter to reduce noise impact)
        def smooth(data, window=5):
            """Simple moving average filter."""
            if len(data.shape) == 1:
                return np.convolve(data, np.ones(window)/window, mode='valid')
            else:
                return np.array([np.convolve(row, np.ones(window)/window, mode='valid') 
                                for row in data])
        
        X_next_smooth = smooth(X_next, window=5)
        X_pred_smooth = smooth(X_pred, window=5)
        residual_smooth = X_next_smooth - X_pred_smooth
        ss_res_smooth = np.sum(residual_smooth**2, axis=1)
        ss_tot_smooth = np.sum((X_next_smooth - X_next_smooth.mean(axis=1, keepdims=True))**2, axis=1)
        r_squared_smooth = 1 - ss_res_smooth / (ss_tot_smooth + 1e-10)
        
        # 4. Pearson correlation coefficient
        corr = np.array([np.corrcoef(X_next[i, :], X_pred[i, :])[0, 1] for i in range(n_states)])
        
        print(f"\nIdentified System Matrices:")
        print(f"A_d (discrete-time state matrix, {n_states}x{n_states}):")
        print(f"{self.A_d}")
        print(f"\nB_d (discrete-time input matrix, {n_states}x1):")
        print(f"{self.B_d.flatten()}")
        
        # Print fit quality metrics with appropriate state names
        print(f"\nFit Quality Metrics ({n_states}D State System):")
        if n_states == 2:
            headers = f"{'':20s} {'Signal':>12s} {'Hot':>12s}"
            separator = '-' * 46
            r2_line = f"{'R² (raw)':20s} {r_squared[0]:12.4f} {r_squared[1]:12.4f}"
            r2_smooth_line = f"{'R² (smoothed)':20s} {r_squared_smooth[0]:12.4f} {r_squared_smooth[1]:12.4f}"
            corr_line = f"{'Correlation':20s} {corr[0]:12.4f} {corr[1]:12.4f}"
            mae_line = f"{'MAE (K)':20s} {mae[0]:12.4f} {mae[1]:12.4f}"
        else:  # n_states == 3
            headers = f"{'':20s} {'Signal':>12s} {'Hot':>12s} {'Kinetic':>12s}"
            separator = '-' * 59
            r2_line = f"{'R² (raw)':20s} {r_squared[0]:12.4f} {r_squared[1]:12.4f} {r_squared[2]:12.4f}"
            r2_smooth_line = f"{'R² (smoothed)':20s} {r_squared_smooth[0]:12.4f} {r_squared_smooth[1]:12.4f} {r_squared_smooth[2]:12.4f}"
            corr_line = f"{'Correlation':20s} {corr[0]:12.4f} {corr[1]:12.4f} {corr[2]:12.4f}"
            mae_line = f"{'MAE (K)':20s} {mae[0]:12.4f} {mae[1]:12.4f} {mae[2]:12.4f}"
        
        print(headers)
        print(separator)
        print(r2_line)
        print(r2_smooth_line)
        print(corr_line)
        print(mae_line)
        print(f"\nNOTE: Bath is control input (u), NOT a dynamic state!")
        
        # Warn if fit quality is poor
        min_r_squared = np.min(r_squared)
        if min_r_squared < 0.5:
            print(f"\n⚠⚠⚠ WARNING: POOR SYSTEM ID FIT (min R²={min_r_squared:.4f}) ⚠⚠⚠")
            print(f"  This will cause poor LQR performance!")
            print(f"  Common causes:")
            print(f"    1. System ID duration too short (current: {self.system_id_duration_ps}ps)")
            print(f"    2. States don't respond to bath changes (check thermostat coupling)")
            print(f"    3. In 3D mode: T_kinetic thermalizes very slowly")
            print(f"  Suggestions:")
            print(f"    - Increase --lqr-system-id-duration (try 200-500ps)")
            print(f"    - Reduce --molecular-tau for stronger thermostat")
            print(f"    - For 3D mode: Expect T_kinetic R² to be low if tau >> duration")
        
        # Plot fit quality
        self._plot_system_id_fit(X_next, X_pred, u_input[:-1], np.array(self.system_id_data['time'])[1:])
        
        # Extract continuous-time parameters (approximate for small dt)
        # Use robust extraction with safeguards for very small dt
        try:
            # CRITICAL SAFEGUARD: Ensure dt is valid before division
            if self.dt <= 0:
                print(f"\n⚠⚠⚠ CRITICAL ERROR: dt = {self.dt} ps (invalid for parameter extraction!) ⚠⚠⚠")
                print(f"  Cannot extract continuous-time parameters with dt ≤ 0")
                print(f"  This will cause Infinity/-Infinity in continuous parameters!")
                print(f"  Using emergency fallback dt = 0.01 ps")
                self.dt = 0.01
            
            # For very small dt (< 0.01 ps), use matrix logarithm for better accuracy
            if self.dt < 0.01:
                print(f"\n⚠ Very small timestep (dt={self.dt:.4f}ps) detected!")
                print(f"  Using matrix logarithm method for continuous-time parameter extraction...")
                # A_c = (1/dt) * log(A_d), but use approximate A_c = (A_d - I)/dt for stability
                A_c_approx = (self.A_d - np.eye(n_states)) / self.dt
                B_c_approx = self.B_d / self.dt
            else:
                A_c_approx = (self.A_d - np.eye(n_states)) / self.dt  # Standard first-order approximation
                B_c_approx = self.B_d / self.dt
            
            k_s = -A_c_approx[0, 0]  # Signal decay rate
            k_h = -A_c_approx[1, 1]  # Hot decay rate
            g_sh = A_c_approx[0, 1]  # Signal-hot cross-coupling
            g_hs = A_c_approx[1, 0]  # Hot-signal cross-coupling
            b_s = B_c_approx[0, 0]   # Signal-bath coupling
            b_h = B_c_approx[1, 0]   # Hot-bath coupling
            
            # Clamp extreme values to prevent overflow in display
            def safe_print(val, name, unit=""):
                if np.abs(val) > 1e6:
                    return f"  {name} = >1e6 ({unit}) [CLIPPED - dt too small for accurate extraction]"
                elif np.isnan(val) or np.isinf(val):
                    return f"  {name} = invalid ({unit})"
                else:
                    return f"  {name} = {val:.4f} ({unit})"
            
            print(f"\nApproximate Continuous-Time Parameters:")
            print(safe_print(k_s, "k_s", "1/ps") + (f" → τ_s = {1/k_s:.2f} ps" if 0 < k_s < 1e6 else ""))
            print(safe_print(k_h, "k_h", "1/ps") + (f" → τ_h = {1/k_h:.2f} ps" if 0 < k_h < 1e6 else ""))
            print(safe_print(g_sh, "g_sh", "signal ← hot coupling"))
            print(safe_print(g_hs, "g_hs", "hot ← signal coupling"))
            print(safe_print(b_s, "b_s", "signal ← bath coupling"))
            print(safe_print(b_h, "b_h", "hot ← bath coupling"))
            
            if n_states == 3:
                k_kin = -A_c_approx[2, 2]  # Kinetic decay rate
                g_skin = A_c_approx[0, 2]  # Signal-kinetic cross-coupling
                g_hkin = A_c_approx[1, 2]  # Hot-kinetic cross-coupling
                g_kins = A_c_approx[2, 0]  # Kinetic-signal cross-coupling
                g_kinh = A_c_approx[2, 1]  # Kinetic-hot cross-coupling
                b_kin = B_c_approx[2, 0]   # Kinetic-bath coupling
                print(safe_print(k_kin, "k_kin", "1/ps") + (f" → τ_kin = {1/k_kin:.2f} ps" if 0 < k_kin < 1e6 else ""))
                print(safe_print(g_skin, "g_skin", "signal ← kinetic coupling"))
                print(safe_print(g_hkin, "g_hkin", "hot ← kinetic coupling"))
                print(safe_print(g_kins, "g_kins", "kinetic ← signal coupling"))
                print(safe_print(g_kinh, "g_kinh", "kinetic ← hot coupling"))
                print(safe_print(b_kin, "b_kin", "kinetic ← bath coupling"))
                
        except Exception as e:
            print(f"\n⚠ Warning: Continuous-time parameter extraction failed: {e}")
            print(f"  (This is non-critical - discrete-time parameters are used for control)")
        
        # Save parameters to file
        # Formulation string depends on whether kinetic tracking is enabled
        if self.track_kinetic_temp:
            formulation_str = "LQG: 3D state [x_s, x_h, x_kin], 1D control [u_bath]"
        else:
            formulation_str = "LQG: 2D state [x_s, x_h], 1D control [u_bath]"
        
        params = {
            "identification_info": {
                "timestamp": datetime.now().isoformat(),
                "initial_temperature": self.target_temperature,
                "quench_temperature": self.system_id_temp_K,
                "duration_ps": self.system_id_duration_ps,
                "sample_time_ps": self.dt,
                "num_data_points": N,
                "formulation": formulation_str
            },
            "system_matrices": {
                "A_discrete": self.A_d.tolist(),
                "B_discrete": self.B_d.flatten().tolist()
            },
            "continuous_parameters": {
                "k_s": float(k_s),
                "k_h": float(k_h),
                "g_sh": float(g_sh),
                "g_hs": float(g_hs),
                "b_s": float(b_s),
                "b_h": float(b_h)
            },
            "fit_quality": {
                "r_squared_raw": {
                    "signal": float(r_squared[0]),
                    "hot": float(r_squared[1])
                },
                "r_squared_smoothed": {
                    "signal": float(r_squared_smooth[0]),
                    "hot": float(r_squared_smooth[1])
                },
                "correlation": {
                    "signal": float(corr[0]),
                    "hot": float(corr[1])
                },
                "mae_K": {
                    "signal": float(mae[0]),
                    "hot": float(mae[1])
                }
            }
        }
        
        try:
            # Save with unique ID number
            base_filename = self.system_id_file.replace('.json', '')
            unique_filename = f"{base_filename}_id{self.system_id_count}.json"
            with open(unique_filename, 'w') as f:
                json.dump(params, f, indent=2)
            print(f"\nSystem parameters saved to: {unique_filename}")
            
            # Also save "latest" version that gets overwritten
            with open(self.system_id_file, 'w') as f:
                json.dump(params, f, indent=2)
            print(f"  (Latest version also saved to: {self.system_id_file})")
        except Exception as e:
            print(f"Warning: Could not save system parameters: {e}")
    
    def _plot_system_id_fit(self, X_actual, X_pred, u, time):
        """
        Plot actual vs predicted trajectories from system identification.
        
        LQG formulation: 2D or 3D state system
        
        Parameters
        ----------
        X_actual : np.ndarray
            Actual state trajectories (n_states, N): [signal, hot] or [signal, hot, kinetic]
        X_pred : np.ndarray
            Predicted state trajectories (n_states, N): [signal, hot] or [signal, hot, kinetic]
        u : np.ndarray
            Control input trajectory (N,): bath temperature deviation
        time : np.ndarray
            Time vector (N,)
        """
        try:
            import matplotlib
            matplotlib.use('Agg')  # Non-interactive backend
            import matplotlib.pyplot as plt
        except ImportError:
            print("Warning: matplotlib not available, skipping system ID fit plot")
            return
        
        # Smooth trajectories for visualization
        def smooth_1d(data, window=5):
            return np.convolve(data, np.ones(window)/window, mode='same')
        
        n_states = X_actual.shape[0]
        n_plots = n_states + 1  # States + control input
        fig, axes = plt.subplots(n_plots, 1, figsize=(12, 3*n_plots), sharex=True)
        
        # Title
        state_str = f"{n_states}D State " + ("[x_s, x_h]" if n_states == 2 else "[x_s, x_h, x_kin]")
        axes[0].set_title(f'LQG System ID: {state_str} + 1D Control [u_bath]', fontsize=14, fontweight='bold')
        
        # Signal temperature (state 1)
        axes[0].plot(time, X_actual[0, :] + self.target_temperature, 'b-', linewidth=1, label='Actual (raw)', alpha=0.3)
        axes[0].plot(time, X_pred[0, :] + self.target_temperature, 'r-', linewidth=1, label='Predicted (raw)', alpha=0.3)
        axes[0].plot(time, smooth_1d(X_actual[0, :], 10) + self.target_temperature, 'b-', linewidth=2.5, label='Actual (smoothed)', alpha=0.9)
        axes[0].plot(time, smooth_1d(X_pred[0, :], 10) + self.target_temperature, 'r--', linewidth=2.5, label='Predicted (smoothed)')
        axes[0].set_ylabel('Signal Temp (K)', fontsize=12, fontweight='bold')
        axes[0].legend(loc='best', fontsize=9, ncol=2)
        axes[0].grid(True, alpha=0.3)
        
        # Hot temperature (state 2)
        axes[1].plot(time, X_actual[1, :] + self.target_temperature, 'b-', linewidth=1, label='Actual (raw)', alpha=0.3)
        axes[1].plot(time, X_pred[1, :] + self.target_temperature, 'r-', linewidth=1, label='Predicted (raw)', alpha=0.3)
        axes[1].plot(time, smooth_1d(X_actual[1, :], 10) + self.target_temperature, 'b-', linewidth=2.5, label='Actual (smoothed)', alpha=0.9)
        axes[1].plot(time, smooth_1d(X_pred[1, :], 10) + self.target_temperature, 'r--', linewidth=2.5, label='Predicted (smoothed)')
        axes[1].set_ylabel('Hot Temp (K)', fontsize=12, fontweight='bold')
        axes[1].legend(loc='best', fontsize=9, ncol=2)
        axes[1].grid(True, alpha=0.3)
        
        # Kinetic temperature (state 3, if 3D)
        if n_states == 3:
            axes[2].plot(time, X_actual[2, :] + self.target_temperature, 'b-', linewidth=1, label='Actual (raw)', alpha=0.3)
            axes[2].plot(time, X_pred[2, :] + self.target_temperature, 'r-', linewidth=1, label='Predicted (raw)', alpha=0.3)
            axes[2].plot(time, smooth_1d(X_actual[2, :], 10) + self.target_temperature, 'b-', linewidth=2.5, label='Actual (smoothed)', alpha=0.9)
            axes[2].plot(time, smooth_1d(X_pred[2, :], 10) + self.target_temperature, 'r--', linewidth=2.5, label='Predicted (smoothed)')
            axes[2].set_ylabel('Kinetic Temp (K)', fontsize=12, fontweight='bold')
            axes[2].legend(loc='best', fontsize=9, ncol=2)
            axes[2].grid(True, alpha=0.3)
        
        # Control input (bath temperature deviation) - last subplot
        ctrl_ax = axes[-1]
        ctrl_ax.plot(time, u, 'g-', linewidth=2, label='Control Input u (bath deviation)', alpha=0.8)
        ctrl_ax.plot(time, u + self.target_temperature, 'c--', linewidth=1.5, label='Bath Temperature', alpha=0.6)
        ctrl_ax.set_ylabel('Temperature (K)', fontsize=12, fontweight='bold')
        ctrl_ax.set_xlabel('Time (ps)', fontsize=12, fontweight='bold')
        ctrl_ax.legend(loc='best', fontsize=10)
        ctrl_ax.grid(True, alpha=0.3)
        ctrl_ax.axhline(0, color='k', linestyle=':', linewidth=1, label='Target')
        
        plt.tight_layout()
        
        # Save plot with unique ID number
        base_filename = self.system_id_file.replace('.json', '')
        plot_filename = f"{base_filename}_fit_id{self.system_id_count}.png"
        plt.savefig(plot_filename, dpi=150, bbox_inches='tight')
        print(f"✓ System ID fit plot saved to: {plot_filename}")
        
        # Also save a "latest" version that gets overwritten
        latest_filename = f"{base_filename}_fit.png"
        plt.savefig(latest_filename, dpi=150, bbox_inches='tight')
        print(f"  (Latest version also saved to: {latest_filename})")
        
        plt.close()
    
    def _design_lqr_controller(self):
        """Design LQR controller with integral action, Kalman filter, and cross-coupling.
        
        LQG formulation:
        - Base system: x = [x_s, x_h] (2D) or [x_s, x_h, x_kin] (3D)
        - Augmented system adds integral: [x_s, x_h, ξ] (3D) or [x_s, x_h, x_kin, ξ] (4D)
        - LQR acts on augmented system: u = -K @ x_a
        - Kalman filter estimates base system only (separation principle)
        
        Cross-coupling (3D mode only):
        - Q matrix includes penalties on (T_i - T_j)² to enforce thermal equilibration
        - Examples: (T_signal - T_kinetic)², (T_signal - T_hot)², (T_hot - T_kinetic)²
        - Ensures Q remains PSD via rank-1 outer product construction
        """
        print(f"\n{'='*70}")
        print(f"Designing LQG Controller with Cross-Coupling")
        print(f"{'='*70}")
        
        n_states = self.state_dim  # 2 or 3
        n_aug = n_states + 1  # 3 or 4 (states + integral)
        
        # Augment state_dim system with integrator: x_a = [x_s; x_h; (x_kin); ξ]
        A_aug = np.zeros((n_aug, n_aug))
        A_aug[:n_states, :n_states] = self.A_d  # Copy state_dim x state_dim state matrix
        A_aug[n_states, 0] = 1.0  # ξ[k+1] = ξ[k] + x_s[k]
        A_aug[n_states, n_states] = 1.0
        
        B_aug = np.zeros((n_aug, 1))
        B_aug[:n_states, :] = self.B_d  # Copy state_dim x 1 input matrix
        
        # LQR weighting matrices with cross-coupling for thermal equilibration
        # Start with individual state penalties (diagonal)
        if n_states == 2:
            Q = np.diag([self.lqr_weight_signal,   # Penalize signal deviation
                         self.lqr_weight_hot,       # Penalize hot deviation
                         self.lqr_weight_integral])  # Penalize integral (tracking)
        else:  # n_states == 3
            # Individual penalties
            Q_individual = np.diag([self.lqr_weight_signal,   # Penalize signal deviation
                                   self.lqr_weight_hot,       # Penalize hot deviation
                                   self.weight_kinetic,       # Penalize kinetic deviation
                                   self.lqr_weight_integral])  # Penalize integral (tracking)
            
            # Cross-coupling penalties: penalize (T_i - T_j)²
            # This enforces thermal equilibration between subsystems
            def pairwise_penalty(i, j, n):
                """Construct PSD matrix that penalizes (x_i - x_j)²"""
                e_diff = np.zeros((n, 1))
                e_diff[i, 0] = 1.0
                e_diff[j, 0] = -1.0
                return e_diff @ e_diff.T  # Rank-1 PSD matrix
            
            # Add cross-coupling terms (only if weights > 0)
            Q = Q_individual.copy()
            if self.cross_coupling_signal_kinetic > 0:
                Q[:n_aug, :n_aug] += self.cross_coupling_signal_kinetic * pairwise_penalty(0, 2, n_aug)
            if self.cross_coupling_signal_hot > 0:
                Q[:n_aug, :n_aug] += self.cross_coupling_signal_hot * pairwise_penalty(0, 1, n_aug)
            if self.cross_coupling_hot_kinetic > 0:
                Q[:n_aug, :n_aug] += self.cross_coupling_hot_kinetic * pairwise_penalty(1, 2, n_aug)
        
        R = np.array([[self.lqr_control_effort]])  # Penalize control effort
        
        # Verify Q is positive semi-definite (PSD check)
        Q_eigenvalues = np.linalg.eigvalsh(Q)  # Hermitian eigenvalue solver (more stable)
        min_eigenvalue = np.min(Q_eigenvalues)
        if min_eigenvalue < -1e-10:  # Allow small numerical errors
            raise ValueError(f"Q matrix is not PSD! Minimum eigenvalue: {min_eigenvalue:.6e}. "
                           "Reduce cross-coupling weights or increase individual state weights.")
        
        print(f"LQR Weighting Matrices:")
        print(f"  Q (state weights diagonal): {np.diag(Q)}")
        if n_states == 3 and (self.cross_coupling_signal_kinetic > 0 or 
                              self.cross_coupling_signal_hot > 0 or 
                              self.cross_coupling_hot_kinetic > 0):
            print(f"  Q (cross-coupling): s-kin={self.cross_coupling_signal_kinetic:.2f}, "
                  f"s-hot={self.cross_coupling_signal_hot:.2f}, hot-kin={self.cross_coupling_hot_kinetic:.2f}")
            print(f"  Q eigenvalues (PSD check): min={min_eigenvalue:.6e}, max={np.max(Q_eigenvalues):.6e}")
        print(f"  R (control effort): {R[0, 0]}")
        
        # Solve discrete-time algebraic Riccati equation for augmented system
        try:
            P = solve_discrete_are(A_aug, B_aug, Q, R)
            self.K = np.linalg.inv(R + B_aug.T @ P @ B_aug) @ (B_aug.T @ P @ A_aug)
            
            # K is now 1 x n_aug: acts on augmented state
            print(f"\nComputed LQR Gain (1x{n_aug}):")
            print(f"  K = {self.K.flatten()}")
            if n_states == 2:
                print(f"  K breakdown: [K_signal={self.K[0,0]:.4f}, K_hot={self.K[0,1]:.4f}, K_integral={self.K[0,2]:.4f}]")
            else:  # n_states == 3
                print(f"  K breakdown: [K_signal={self.K[0,0]:.4f}, K_hot={self.K[0,1]:.4f}, K_kinetic={self.K[0,2]:.4f}, K_integral={self.K[0,3]:.4f}]")
            
            # Check closed-loop stability of augmented system
            A_cl = A_aug - B_aug @ self.K
            eigenvalues = np.linalg.eigvals(A_cl)
            print(f"\nClosed-Loop Eigenvalues ({n_aug}D augmented system):")
            for i, eig in enumerate(eigenvalues):
                print(f"  λ_{i+1} = {eig:.4f} (|λ| = {abs(eig):.4f})")
            
            if np.all(np.abs(eigenvalues) < 1.0):
                print("✓ Closed-loop system is stable")
            else:
                print("⚠ Warning: Closed-loop system may be unstable!")
        
        except Exception as e:
            print(f"Error solving LQR problem: {e}")
            raise
        
        # Design Kalman filter for state_dim base system estimation
        # Estimate [x_s, x_h], [x_s, x_h, x_kin], or [x_s, x_h, x_kin, x_bath]
        # (integral ξ is deterministic, no noise)
        # Separation principle: Kalman filter on base system, LQR on augmented system
        if n_states == 2:
            Q_process = np.diag([self.process_noise_signal**2, 
                                self.process_noise_hot**2])  # Process noise (2D)
            R_meas = np.diag([self.estimated_noise_signal**2,
                             self.estimated_noise_hot**2])  # Measurement noise (2D)
        else:  # n_states == 3
            Q_process = np.diag([self.process_noise_signal**2, 
                                self.process_noise_hot**2,
                                self.process_noise_kinetic**2])  # Process noise (3D)
            R_meas = np.diag([self.estimated_noise_signal**2,
                             self.estimated_noise_hot**2,
                             self.measurement_noise_kinetic**2])  # Measurement noise (3D)
        
        print(f"\nKalman Filter Design ({n_states}D base system):")
        if n_states == 2:
            print(f"  Process noise (σ_w): [{self.process_noise_signal:.2f}, {self.process_noise_hot:.2f}] K")
            print(f"  Measurement noise (σ_v): [{self.estimated_noise_signal:.2f}, {self.estimated_noise_hot:.2f}] K")
        else:  # n_states == 3
            print(f"  Process noise (σ_w): [{self.process_noise_signal:.2f}, {self.process_noise_hot:.2f}, {self.process_noise_kinetic:.2f}] K")
            print(f"  Measurement noise (σ_v): [{self.estimated_noise_signal:.2f}, {self.estimated_noise_hot:.2f}, {self.measurement_noise_kinetic:.2f}] K")
        if self.enable_adaptation:
            print(f"  (Adaptive noise estimation enabled)")
        
        # Solve Riccati for Kalman gain (state_dim system)
        try:
            P_kf = solve_discrete_are(self.A_d.T, self.C.T, Q_process, R_meas)
            self.L = P_kf @ self.C.T @ np.linalg.inv(self.C @ P_kf @ self.C.T + R_meas)
            
            print(f"\nComputed Kalman Filter Gain ({n_states}x{n_states}):")
            print(f"  L (estimator gain):")
            print(f"{self.L}")
            
            # Check estimator stability (state_dim system)
            A_est = self.A_d - self.L @ self.C @ self.A_d
            eigenvalues_est = np.linalg.eigvals(A_est)
            print(f"\nEstimator Eigenvalues ({n_states}D base system):")
            for i, eig in enumerate(eigenvalues_est):
                print(f"  λ_{i+1} = {eig:.4f} (|λ| = {abs(eig):.4f})")
            
            if np.all(np.abs(eigenvalues_est) < 1.0):
                print("✓ Estimator is stable")
            else:
                print("⚠ Warning: Estimator may be unstable!")
        
        except Exception as e:
            print(f"Error solving Kalman filter problem: {e}")
            raise
        
        # Initialize EKF for parameter estimation if enabled
        # NOTE: EKF currently only supports 2D states [x_s, x_h]
        # Skip EKF initialization in 3D mode (when track_kinetic_temp=True)
        if self.use_ekf_adaptation and not self.track_kinetic_temp:
            print(f"\nInitializing EKF for Online Parameter Estimation:")
            self.parameter_ekf = ParameterEKF(
                dt=self.dt,
                process_noise_state=self.process_noise_signal,
                process_noise_param=self.ekf_process_noise_param,
                measurement_noise_signal=self.estimated_noise_signal,
                measurement_noise_hot=self.estimated_noise_hot,
                initial_covariance_state=1.0,
                initial_covariance_param=self.ekf_initial_covariance_param
            )
            
            # Initialize with current system ID results
            self.parameter_ekf.initialize(
                x_initial=self.x_hat,
                A_initial=self.A_d,
                B_initial=self.B_d
            )
            
            # Store initial parameters for adaptive redesign tracking
            self.last_lqr_redesign_params = np.concatenate([
                self.A_d.flatten(),
                self.B_d.flatten()
            ])
        elif self.use_ekf_adaptation and self.track_kinetic_temp:
            print(f"\nℹ EKF parameter estimation disabled in 3D mode (kinetic temp tracking enabled)")
            print(f"  (EKF only supports 2D states; use 2D mode for adaptive parameter estimation)")
        
        print(f"{'='*70}\n")
    
    def _trigger_periodic_system_id(self, current_time_ps):
        """Trigger a periodic system identification to re-learn system dynamics.
        
        This allows the controller to adapt to changing system dynamics over time.
        """
        print(f"\n{'='*70}")
        print(f"PERIODIC SYSTEM RE-IDENTIFICATION #{self.system_id_count + 1}")
        print(f"{'='*70}")
        print(f"Time since last system ID: {current_time_ps - self.last_system_id_completion_time:.1f} ps")
        print(f"Re-identifying system parameters at t = {current_time_ps:.2f} ps...")
        
        # CRITICAL: Reset ALL system ID state variables
        self.phase = 'system_id'
        self.system_id_started = False
        self.system_id_current_step = 0  # Reset step counter to start from beginning
        self.system_id_sample_count = 0  # Reset sample counter
        
        # Clear old system ID steps to force regeneration
        if hasattr(self, 'system_id_steps'):
            del self.system_id_steps
        
        # Clear old data buffers
        self.system_id_data = {
            'time': [],
            'x_s': [],
            'x_h': [],
            'x_kin': [],  # Kinetic temperature (only used if track_kinetic_temp=True)
            'x_c': [],
            'u': []
        }
        
        # Clear median filter buffers to start fresh
        self.temp_buffer_signal = []
        self.temp_buffer_hot = []
        
        print(f"System will perform {self.system_id_mode} excitation for {self.system_id_duration_ps:.1f} ps")
        print(f"Control will resume automatically after re-identification completes.")
        print(f"{'='*70}\n")
    
    def _execute_lqr_control(self, timestep, current_time_ps):
        """Execute LQR optimal control."""
        # Only update at specified intervals
        if current_time_ps - self.last_update_time < self.update_interval_ps:
            return
        
        dt_ps = current_time_ps - self.last_update_time if self.last_update_time > 0 else self.update_interval_ps
        self.last_update_time = current_time_ps
        
        # Measure temperatures (raw)
        signal_temp_raw = self._measure_temperature(self.signal_temperature_method, current_time_ps)
        hot_temp_raw = self._measure_temperature(self.hot_temperature_method, current_time_ps)
        bath_temp = self._get_current_bath_temperature()
        
        if signal_temp_raw is None or hot_temp_raw is None or bath_temp is None:
            return
        
        # Apply median filter for outlier rejection
        signal_temp = self._apply_median_filter(signal_temp_raw, self.temp_buffer_signal, self.median_filter_window)
        hot_temp_median = self._apply_median_filter(hot_temp_raw, self.temp_buffer_hot, self.median_filter_window)
        
        # Apply low-pass filter to T_h to suppress high-frequency oscillations
        # This prevents the controller from chasing natural thermal fluctuations
        if self.T_h_filter_enabled:
            # Initialize filter on first measurement
            if self.T_h_filtered is None:
                self.T_h_filtered = hot_temp_median
            
            # Exponential moving average: x_filt[k] = α·x[k] + (1-α)·x_filt[k-1]
            # Time constant τ determines cutoff frequency: α = Δt/(Δt + τ)
            alpha_filter = self.dt / (self.dt + self.T_h_filter_time_constant)
            self.T_h_filtered = alpha_filter * hot_temp_median + (1 - alpha_filter) * self.T_h_filtered
            
            # Use filtered temperature for control (log unfiltered for diagnostics)
            hot_temp = self.T_h_filtered
            
            # Diagnostic: show filtering effect occasionally
            if self.num_updates <= 5 or self.num_updates % 200 == 0:
                print(f"[T_h Low-Pass Filter] Raw: {hot_temp_median:.2f}K → Filtered: {self.T_h_filtered:.2f}K "
                      f"(τ={self.T_h_filter_time_constant:.1f}ps, suppresses >~{1000/(2*np.pi*self.T_h_filter_time_constant):.1f}Hz)")
        else:
            hot_temp = hot_temp_median
        
        # Increment update counter
        self.num_updates += 1
        
        # GENTLE STARTUP: Ramp up control authority gradually
        # Prevents violent reactions to initial transients
        # Formula: linear ramp from min_authority to 1.0 over startup_steps
        if self.num_updates <= self.gentle_startup_steps:
            # Linear interpolation: min_authority at step 1, 1.0 at step N
            progress = self.num_updates / self.gentle_startup_steps
            startup_factor = self.gentle_startup_min_authority + (1.0 - self.gentle_startup_min_authority) * progress
            startup_factor = min(1.0, startup_factor)
        else:
            startup_factor = 1.0
        
        # Log filtering effect occasionally
        if self.num_updates <= 5 or self.num_updates % 100 == 0:
            print(f"[Median Filter] Signal: {signal_temp_raw:.2f}K → {signal_temp:.2f}K | "
                  f"Hot: {hot_temp_raw:.2f}K → {hot_temp_median:.2f}K")
            if self.num_updates <= self.gentle_startup_steps:
                print(f"[Gentle Startup] Control authority: {startup_factor*100:.1f}% (step {self.num_updates}/{self.gentle_startup_steps})")
        
        # Measure kinetic temperature if tracking is enabled
        if self.track_kinetic_temp:
            kinetic_temp_raw = self._measure_temperature('kinetic', current_time_ps)
            if kinetic_temp_raw is None:
                if self.num_updates <= 5:
                    print(f"[WARNING] Kinetic temperature measurement failed at t={current_time_ps:.1f}ps")
                return
            # Apply median filter to kinetic temperature
            kinetic_temp = self._apply_median_filter(kinetic_temp_raw, self.temp_buffer_signal, self.median_filter_window)
        else:
            kinetic_temp = None
        
        # Deviations from target
        x_s_meas = signal_temp - self.target_temperature
        x_h_meas = hot_temp - self.target_temperature
        x_bath_meas = bath_temp - self.target_temperature
        
        # Measurement vector (2D or 3D depending on configuration)
        if self.track_kinetic_temp:
            # 3D mode: [x_s, x_h, x_kin]
            x_kin_meas = kinetic_temp - self.target_temperature
            y = np.array([[x_s_meas], [x_h_meas], [x_kin_meas]])
        else:
            # 2D mode: [x_s, x_h]
            y = np.array([[x_s_meas], [x_h_meas]])
        
        # ADAPTIVE: Update noise estimates from actual signal variance
        if self.enable_adaptation:
            if self.track_kinetic_temp:
                self._update_noise_estimates(signal_temp, hot_temp, current_time_ps, kinetic_temp)
            else:
                self._update_noise_estimates(signal_temp, hot_temp, current_time_ps)
        
        # Save current state for EKF/RLS (state_dim state + 1D control)
        if self.track_kinetic_temp:
            x_current = np.array([[x_s_meas], [x_h_meas], [x_kin_meas]])
        else:
            x_current = np.array([[x_s_meas], [x_h_meas]])
        
        # ===== KALMAN FILTER: Estimate base system [x_s, x_h] or [x_s, x_h, x_kin] =====
        # Prediction: x_hat[k|k-1] = A_d @ x_hat[k-1|k-1] + B_d @ u[k-1]
        if not hasattr(self, 'u_prev'):
            self.u_prev = 0.0  # Initialize if missing
        
        x_hat_pred = self.A_d @ self.x_hat + self.B_d * self.u_prev
        
        # Innovation: r = y - C @ x_hat[k|k-1]
        y_pred = self.C @ x_hat_pred
        innovation = y - y_pred
        
        # Update: x_hat[k|k] = x_hat[k|k-1] + L @ innovation
        self.x_hat = x_hat_pred + self.L @ innovation
        
        # ===== INTEGRAL STATE: Deterministic integration (no Kalman) =====
        # ξ[k+1] = ξ[k] + x_s[k] * dt
        # CRITICAL: Skip during gentle startup to avoid windup during Kalman settling
        kalman_settling_steps = min(3, self.gentle_startup_steps // 3)  # Settle for first 1/3 of gentle startup
        if self.num_updates >= self.gentle_startup_steps:
            self.xi_hat[0, 0] += self.x_hat[0, 0] * dt_ps
        elif self.num_updates < kalman_settling_steps:
            print(f"[Gentle Startup] Skipping integral update (Kalman settling)")
        
        # Anti-windup: clamp integral
        if abs(self.xi_hat[0, 0]) > self.integral_max:
            self.xi_hat[0, 0] = np.sign(self.xi_hat[0, 0]) * self.integral_max
        
        # ===== GAIN SCHEDULING: Adapt control authority based on operating regime =====
        if self.enable_gain_scheduling:
            # Compute temperature spread (measure of proximity to equilibrium)
            # In 3D mode: consider ALL temperature spreads (not just T_h vs T_s!)
            if self.track_kinetic_temp:
                # Maximum spread across all three temperatures
                temps = [signal_temp, hot_temp, kinetic_temp]
                delta_T = max(temps) - min(temps)
            else:
                # 2D mode: only signal vs hot
                delta_T = abs(hot_temp - signal_temp)
            
            # Gain scheduling: reduce gain near equilibrium (high sensitivity region)
            if delta_T > self.gain_schedule_far_threshold:
                # Far from equilibrium: full gain
                alpha_scheduled = 1.0
            elif delta_T < self.gain_schedule_near_threshold:
                # Near equilibrium: reduced gain (prevent aggressive control in sensitive region)
                alpha_scheduled = 0.4
            else:
                # Transition region: linear interpolation
                alpha_scheduled = 0.4 + 0.6 * (delta_T - self.gain_schedule_near_threshold) / (
                    self.gain_schedule_far_threshold - self.gain_schedule_near_threshold
                )
            
            # Smooth transition (exponential smoothing)
            alpha_smoothing = 0.1
            self.gain_schedule_alpha = alpha_smoothing * alpha_scheduled + (1 - alpha_smoothing) * self.gain_schedule_alpha
            
            # Log gain scheduling occasionally
            if self.num_updates % 50 == 0 and self.num_updates > 10:
                if self.track_kinetic_temp:
                    print(f"[Gain Scheduling] ΔT_max={delta_T:.1f}K (across T_s, T_h, T_kin) → α={self.gain_schedule_alpha:.3f}")
                else:
                    print(f"[Gain Scheduling] ΔT={delta_T:.1f}K → α={self.gain_schedule_alpha:.3f}")
        else:
            self.gain_schedule_alpha = 1.0
        
        # ===== LQR CONTROL LAW: u = -α(regime) * K @ [x_s; x_h; (x_kin); ξ] =====
        # Build augmented state vector (state_dim + 1)
        x_aug = np.vstack([self.x_hat, self.xi_hat])  # (state_dim+1, 1): [x_s, x_h, (x_kin), ξ]
        u_deviation = -self.gain_schedule_alpha * self.K @ x_aug  # (1, 1) with gain scheduling
        
        # Apply gentle startup scaling
        u_deviation_scaled = u_deviation * startup_factor
        u_bath_temp_desired = self.target_temperature + u_deviation_scaled[0, 0]
        
        # Rate limiting: limit how fast bath temperature can change
        if self.u_prev_temp is not None:
            max_change = self.max_control_rate * dt_ps
            change = u_bath_temp_desired - self.u_prev_temp
            if abs(change) > max_change:
                # Limit the change
                u_bath_temp_limited = self.u_prev_temp + np.sign(change) * max_change
            else:
                u_bath_temp_limited = u_bath_temp_desired
        else:
            # First call: no rate limiting
            u_bath_temp_limited = u_bath_temp_desired
        
        # Apply saturation limits
        if self.T_max is not None:
            u_bath_temp_final = max(self.T_min, min(self.T_max, u_bath_temp_limited))
        else:
            u_bath_temp_final = max(self.T_min, u_bath_temp_limited)
        
        # Anti-windup back-calculation: if saturated, reduce integral
        # This prevents integral windup when control saturates
        if u_bath_temp_final != u_bath_temp_desired:
            saturation_error = (u_bath_temp_final - u_bath_temp_desired)
            # Back-calculate integral adjustment (simple proportional back-calc)
            k_back = 0.5  # Back-calculation gain
            self.xi_hat[0, 0] += saturation_error * k_back * dt_ps
            # Re-clamp
            if abs(self.xi_hat[0, 0]) > self.integral_max:
                self.xi_hat[0, 0] = np.sign(self.xi_hat[0, 0]) * self.integral_max
        
        # Update thermostats
        self._update_thermostats(u_bath_temp_final)
        
        # Store for next iteration
        self.u_prev_temp = u_bath_temp_final  # For rate limiting
        self.u_prev = u_bath_temp_final - self.target_temperature  # For Kalman prediction
        
        # ADAPTIVE: EKF-based parameter estimation and adaptive LQR redesign
        if self.enable_adaptation and self.use_ekf_adaptation and self.parameter_ekf is not None:
            # EKF predict-update cycle for parameter estimation
            self.adaptation_counter += 1
            
            if self.adaptation_counter >= self.ekf_update_interval:
                self.adaptation_counter = 0
                
                # EKF Prediction: propagate parameters forward using control input
                self.parameter_ekf.predict(self.u_prev)
                
                # EKF Update: correct with measurement
                self.parameter_ekf.update(y)
                
                # Get estimated parameters
                A_ekf, B_ekf = self.parameter_ekf.get_parameter_estimates()
                
                # Check if parameters changed significantly → trigger LQR redesign (WITH SAFEGUARDS)
                if self.last_lqr_redesign_params is not None:
                    current_params = np.concatenate([A_ekf.flatten(), B_ekf.flatten()])
                    param_change_norm = np.linalg.norm(current_params - self.last_lqr_redesign_params)
                    param_base_norm = np.linalg.norm(self.last_lqr_redesign_params)
                    relative_change = param_change_norm / (param_base_norm + 1e-10)
                    
                    # SAFEGUARD 1: Stability verification (eigenvalues < 0.99 for discrete-time)
                    A_eigs = np.linalg.eigvals(A_ekf)
                    max_eig = np.max(np.abs(A_eigs))
                    is_stable = max_eig < 0.99
                    
                    # SAFEGUARD 2: Confidence gating (check EKF parameter uncertainty)
                    P_theta = self.parameter_ekf.P[2:, 2:]  # Parameter covariance block
                    param_uncertainty = np.trace(P_theta) / P_theta.shape[0]  # Average variance
                    uncertainty_threshold = 0.05  # Only redesign if σ²_θ < 0.05
                    is_confident = param_uncertainty < uncertainty_threshold
                    
                    # SAFEGUARD 3: Minimum convergence time (wait at least 10 EKF updates)
                    min_updates_before_redesign = 10
                    has_converged = self.parameter_ekf.num_updates >= min_updates_before_redesign
                    
                    # SAFEGUARD 4: Parameter validation (eigenvalues haven't jumped too much)
                    A_old_eigs = np.linalg.eigvals(self.A_d)
                    max_eig_old = np.max(np.abs(A_old_eigs))
                    eig_change = abs(max_eig - max_eig_old)
                    max_eig_change = 0.15  # Don't allow eigenvalues to change by more than 0.15
                    eig_reasonable = eig_change < max_eig_change
                    
                    # Decision logic: ALL safeguards must pass
                    if relative_change > self.adaptive_lqr_threshold and is_stable and is_confident and has_converged and eig_reasonable:
                        # All safeguards passed → safe to redesign LQR
                        print(f"\n[Adaptive LQG] Parameters changed by {relative_change*100:.1f}% → SAFE redesign")
                        print(f"  ✓ Stability: max|λ|={max_eig:.4f} < 0.99")
                        print(f"  ✓ Confidence: σ²_θ={param_uncertainty:.4f} < {uncertainty_threshold}")
                        print(f"  ✓ Convergence: {self.parameter_ekf.num_updates} updates ≥ {min_updates_before_redesign}")
                        print(f"  ✓ Eigenvalue change: Δ|λ|={eig_change:.4f} < {max_eig_change}")
                        
                        # Update system matrices with EKF estimates
                        self.A_d = A_ekf.copy()
                        self.B_d = B_ekf.copy()
                        
                        # Re-design LQR controller with new parameters
                        self._design_lqr_controller()
                        
                        # Update tracking
                        self.last_lqr_redesign_params = current_params.copy()
                        self.num_lqr_redesigns += 1
                        
                        print(f"✓ LQR redesigned (total redesigns: {self.num_lqr_redesigns})")
                        print(f"  New A_d eigenvalues: {A_eigs}")
                    else:
                        # Safeguards failed → reject redesign
                        reasons = []
                        if relative_change <= self.adaptive_lqr_threshold:
                            reasons.append(f"insufficient change ({relative_change*100:.1f}% ≤ {self.adaptive_lqr_threshold*100:.1f}%)")
                        if not is_stable:
                            reasons.append(f"unstable (max|λ|={max_eig:.4f} ≥ 0.99)")
                        if not is_confident:
                            reasons.append(f"low confidence (σ²_θ={param_uncertainty:.4f} ≥ {uncertainty_threshold})")
                        if not has_converged:
                            reasons.append(f"insufficient convergence ({self.parameter_ekf.num_updates} < {min_updates_before_redesign})")
                        if not eig_reasonable:
                            reasons.append(f"eigenvalue jump too large (Δ|λ|={eig_change:.4f} ≥ {max_eig_change})")
                        
                        if self.num_updates % 100 == 0 and len(reasons) > 0:
                            print(f"[Adaptive] Redesign blocked: {', '.join(reasons)}")
                
                # Log EKF statistics occasionally
                if self.num_updates % 100 == 0:
                    A_std, B_std = self.parameter_ekf.get_parameter_uncertainty()
                    innov_stats = self.parameter_ekf.get_innovation_statistics()
                    print(f"[EKF] Parameter uncertainty: A_std={np.mean(A_std):.4f}, B_std={np.mean(B_std):.4f}")
                    print(f"[EKF] Innovation: mean={innov_stats['mean']:.3f}K, std={innov_stats['std']:.3f}K")
        
        elif self.enable_adaptation and not self.use_ekf_adaptation:
            # Fallback to RLS if EKF disabled
            # Next state measurement (2D state only, measured on next step)
            x_next = x_current  # Current measurement becomes "next" for this step
            
            # Collect data (use self.u_prev which was just updated)
            self._collect_rls_data(x_current, self.u_prev, x_next)
            
            # Periodic RLS update
            self.adaptation_counter += 1
            if self.adaptation_counter >= self.adaptation_interval:
                self.adaptation_counter = 0
                self._perform_rls_update()
        
        # Console output
        if current_time_ps - self.last_console_output_time >= self.console_output_period_ps:
            if self.enable_adaptation and self.use_ekf_adaptation and self.parameter_ekf is not None:
                adapt_info = f" [EKF:{self.parameter_ekf.num_updates}upd"
                if self.num_lqr_redesigns > 0:
                    adapt_info += f",{self.num_lqr_redesigns}redesign"
                adapt_info += "]"
            elif self.enable_adaptation:
                adapt_info = f" [R²={self.last_r_squared[0]:.2f},{self.last_r_squared[1]:.2f}]"
            else:
                adapt_info = ""
            
            gain_info = f" α={self.gain_schedule_alpha:.2f}" if self.enable_gain_scheduling else ""
            
            if self.track_kinetic_temp:
                # 3D mode: show signal, hot, kinetic
                print(f"LQR Control | t={current_time_ps:.1f}ps | "
                      f"T_s={signal_temp:.2f}K | T_h={hot_temp:.2f}K | T_kin={kinetic_temp:.2f}K | T_bath={u_bath_temp_final:.2f}K | "
                      f"e_s={x_s_meas:.2f}K | e_kin={x_kin_meas:.2f}K | ξ={self.xi_hat[0,0]:.2f}{gain_info}{adapt_info}")
            else:
                # 2D mode: show signal, hot only
                print(f"LQR Control | t={current_time_ps:.1f}ps | "
                      f"T_s={signal_temp:.2f}K | T_h={hot_temp:.2f}K | T_bath={u_bath_temp_final:.2f}K | "
                      f"e_s={x_s_meas:.2f}K | ξ={self.xi_hat[0,0]:.2f}{gain_info}{adapt_info}")
            self.last_console_output_time = current_time_ps
        
        # Log data
        try:
            with open(self.output_file, 'a', encoding='utf-8') as f:
                # Build log string based on state dimensionality
                if self.track_kinetic_temp:
                    # 3D: log signal, hot, kinetic temps and estimates
                    f.write(f"{current_time_ps:.6f},control,{signal_temp:.6f},{hot_temp:.6f},{kinetic_temp:.6f},"
                           f"{bath_temp:.6f},{self.x_hat[0,0]:.6f},{self.x_hat[1,0]:.6f},{self.x_hat[2,0]:.6f},"
                           f"{self.xi_hat[0,0]:.6f},{u_deviation_scaled[0,0]:.6f},1\n")
                else:
                    # 2D: log signal, hot temps and estimates (kinetic=0 as placeholder)
                    f.write(f"{current_time_ps:.6f},control,{signal_temp:.6f},{hot_temp:.6f},0.0,"
                           f"{bath_temp:.6f},{self.x_hat[0,0]:.6f},{self.x_hat[1,0]:.6f},0.0,"
                           f"{self.xi_hat[0,0]:.6f},{u_deviation_scaled[0,0]:.6f},1\n")
        except Exception as e:
            if self.num_updates <= 5:
                print(f"[WARNING] Failed to log data: {e}")
    
    def _apply_median_filter(self, value: float, buffer: list, window_size: int) -> float:
        """Apply median filter to reject outliers.
        
        Parameters
        ----------
        value : float
            New measurement
        buffer : list
            Circular buffer of recent measurements
        window_size : int
            Number of samples to use for median
            
        Returns
        -------
        float
            Median-filtered value
        """
        # Add new measurement to buffer
        buffer.append(value)
        
        # Keep buffer at window size
        if len(buffer) > window_size:
            buffer.pop(0)
        
        # Return median of buffer
        return np.median(buffer)
    
    def _measure_temperature(self, method: str, current_time_ps: float) -> Optional[float]:
        """Measure temperature using specified method."""
        # Try temperature tracker first (preferred - has all fictive temperatures ready)
        if self.temperature_tracker is not None:
            try:
                # Map method names to temperature tracker attributes
                method_map = {
                    'lj_coulombic': 'lj_coulombic_fictive_temperature',
                    'harmonic': 'harmonic_fictive_temperature',
                    'harmonic_equipartition': 'harmonic_equipartition_temperature',
                    'kinetic': 'kinetic_temperature'
                }
                
                attr_name = method_map.get(method)
                if attr_name and hasattr(self.temperature_tracker, attr_name):
                    temp = getattr(self.temperature_tracker, attr_name)
                    if temp is not None and temp > 0:
                        return temp
            except Exception as e:
                print(f"[LQR] WARNING: Failed to get {method} from temperature_tracker: {e}")
        
        # Fallback: Original methods using energy_tracker
        try:
            from .analysis import EmpiricalTemperatureData
            
            if method == 'kinetic':
                sim = self._get_hoomd_simulation()
                if sim is None:
                    return None
                with sim.state.cpu_local_snapshot as snap:
                    velocities = snap.particles.velocity
                    masses = snap.particles.mass
                kinetic_energy = 0.5 * np.sum(masses[:, np.newaxis] * (velocities ** 2))
                N_particles = len(masses)
                kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                temperature_K = (2.0 * kinetic_energy) / (3.0 * N_particles * kB_hartree)
                return max(temperature_K, 0.0)
            
            elif method == 'lj_coulombic':
                if self.empirical_data is None:
                    return None
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                lj_energy = energy_dict.get('lj', 0.0)
                coulomb_energy = energy_dict.get('coulombic', 0.0)
                total_energy = lj_energy + coulomb_energy
                temperature = self.empirical_data.calculate_systemic_temperature(total_energy)
                return max(temperature, 0.0)
            
            elif method == 'harmonic' or method == 'harmonic_equipartition':
                if self.energy_tracker is None:
                    return None
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                harmonic_energy = energy_dict.get('harmonic', 0.0)
                sim = self._get_hoomd_simulation()
                if sim is None:
                    return None
                N_particles = sim.state.N_particles
                kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                # Equipartition: <E_harmonic> = (3/2) * N * kB * T
                temperature_K = (2.0 * harmonic_energy) / (3.0 * N_particles * kB_hartree)
                return max(temperature_K, 0.0)
            
            else:
                print(f"Warning: Unknown temperature method '{method}'")
                return None
        
        except Exception as e:
            print(f"Warning: Could not measure temperature with method '{method}': {e}")
            return None
    
    def _get_current_bath_temperature(self) -> Optional[float]:
        """Get current bath temperature from thermostats."""
        try:
            from .utils import PhysicalConstants
            
            # Try molecular thermostat first
            if self.molecular_thermostat is not None:
                try:
                    if hasattr(self.molecular_thermostat, 'kT'):
                        kT = self.molecular_thermostat.kT
                        # Handle variant objects
                        if hasattr(kT, 'value'):
                            kT = kT.value
                        elif callable(kT):
                            kT = kT()
                        return kT / PhysicalConstants.KB_HARTREE_PER_K
                    elif hasattr(self.molecular_thermostat, 'thermostat'):
                        nested = self.molecular_thermostat.thermostat
                        if hasattr(nested, 'kT'):
                            kT = nested.kT
                            if hasattr(kT, 'value'):
                                kT = kT.value
                            elif callable(kT):
                                kT = kT()
                            return kT / PhysicalConstants.KB_HARTREE_PER_K
                except:
                    pass
            
            # Try cavity thermostat as backup
            if self.cavity_thermostat is not None:
                try:
                    if hasattr(self.cavity_thermostat, 'kT'):
                        kT = self.cavity_thermostat.kT
                        if hasattr(kT, 'value'):
                            kT = kT.value
                        elif callable(kT):
                            kT = kT()
                        return kT / PhysicalConstants.KB_HARTREE_PER_K
                    elif hasattr(self.cavity_thermostat, 'thermostat'):
                        nested = self.cavity_thermostat.thermostat
                        if hasattr(nested, 'kT'):
                            kT = nested.kT
                            if hasattr(kT, 'value'):
                                kT = kT.value
                            elif callable(kT):
                                kT = kT()
                            return kT / PhysicalConstants.KB_HARTREE_PER_K
                except:
                    pass
            
            return None
        except Exception as e:
            print(f"Warning: Could not get current bath temperature: {e}")
            return None
    
    def _update_thermostats(self, temperature_K: float):
        """Update thermostat temperatures using HOOMD variant objects."""
        try:
            import hoomd
            from .utils import PhysicalConstants
            kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
            kT_hartree = temperature_K * kB_hartree
            
            # CRITICAL: Must use hoomd.variant.Constant() for updates to take effect!
            updated = []
            if self.apply_to in ['molecular', 'both'] and self.molecular_thermostat is not None:
                if hasattr(self.molecular_thermostat, 'kT'):
                    self.molecular_thermostat.kT = hoomd.variant.Constant(kT_hartree)
                    updated.append('molecular')
                elif hasattr(self.molecular_thermostat, 'thermostat'):
                    # Handle nested thermostat structure
                    nested = self.molecular_thermostat.thermostat
                    if hasattr(nested, 'kT'):
                        nested.kT = hoomd.variant.Constant(kT_hartree)
                        updated.append('molecular(nested)')
            else:
                if self.molecular_thermostat is None and self.apply_to in ['molecular', 'both']:
                    print(f"⚠ LQR WARNING: molecular_thermostat is None, cannot update!")
            
            if self.apply_to in ['cavity', 'both'] and self.cavity_thermostat is not None:
                if hasattr(self.cavity_thermostat, 'kT'):
                    self.cavity_thermostat.kT = hoomd.variant.Constant(kT_hartree)
                    updated.append('cavity')
                elif hasattr(self.cavity_thermostat, 'thermostat'):
                    # Handle nested thermostat structure
                    nested = self.cavity_thermostat.thermostat
                    if hasattr(nested, 'kT'):
                        nested.kT = hoomd.variant.Constant(kT_hartree)
                        updated.append('cavity(nested)')
            else:
                if self.cavity_thermostat is None and self.apply_to in ['cavity', 'both']:
                    print(f"⚠ LQR WARNING: cavity_thermostat is None, cannot update!")
            
            if not updated:
                print(f"⚠ LQR WARNING: No thermostats were updated! apply_to={self.apply_to}")
        
        except Exception as e:
            print(f"Warning: Failed to update thermostats: {e}")
            import traceback
            traceback.print_exc()
    
    def _update_noise_estimates(self, signal_temp: float, hot_temp: float, current_time_ps: float, kinetic_temp: float = None):
        """Estimate measurement noise from signal variance (adaptive noise tuning).
        
        Parameters
        ----------
        signal_temp : float
            Signal temperature (K)
        hot_temp : float
            Hot temperature (K)
        current_time_ps : float
            Current simulation time (ps)
        kinetic_temp : float, optional
            Kinetic temperature (K), only used if track_kinetic_temp=True
        """
        # Add to noise buffers
        self.noise_buffer_signal.append(signal_temp)
        self.noise_buffer_hot.append(hot_temp)
        
        # Add kinetic temp to buffer if tracking
        if self.track_kinetic_temp and kinetic_temp is not None:
            if not hasattr(self, 'noise_buffer_kinetic'):
                self.noise_buffer_kinetic = []
            self.noise_buffer_kinetic.append(kinetic_temp)
        
        # Keep buffers at fixed size
        if len(self.noise_buffer_signal) > self.noise_estimation_window:
            self.noise_buffer_signal.pop(0)
            self.noise_buffer_hot.pop(0)
            if self.track_kinetic_temp and hasattr(self, 'noise_buffer_kinetic'):
                if len(self.noise_buffer_kinetic) > self.noise_estimation_window:
                    self.noise_buffer_kinetic.pop(0)
        
        # Update estimates periodically
        if current_time_ps - self.last_noise_update >= self.noise_update_interval:
            if len(self.noise_buffer_signal) >= 50:  # Need enough data
                # Compute standard deviation (noise estimate)
                signal_std = np.std(self.noise_buffer_signal)
                hot_std = np.std(self.noise_buffer_hot)
                
                # Update with exponential smoothing to avoid jumps
                alpha = 0.3  # Smoothing factor
                self.estimated_noise_signal = alpha * signal_std + (1 - alpha) * self.estimated_noise_signal
                self.estimated_noise_hot = alpha * hot_std + (1 - alpha) * self.estimated_noise_hot
                
                # Update kinetic temp noise estimate if tracking
                if self.track_kinetic_temp and hasattr(self, 'noise_buffer_kinetic'):
                    if len(self.noise_buffer_kinetic) >= 50:
                        kinetic_std = np.std(self.noise_buffer_kinetic)
                        self.measurement_noise_kinetic = alpha * kinetic_std + (1 - alpha) * self.measurement_noise_kinetic
                
                self.last_noise_update = current_time_ps
                
                # Log noise update
                if self.num_adaptations % 10 == 0:  # Log occasionally
                    if self.track_kinetic_temp:
                        print(f"[Adaptive] Noise estimates updated: signal={self.estimated_noise_signal:.3f}K, hot={self.estimated_noise_hot:.3f}K, kinetic={self.measurement_noise_kinetic:.3f}K")
                    else:
                        print(f"[Adaptive] Noise estimates updated: signal={self.estimated_noise_signal:.3f}K, hot={self.estimated_noise_hot:.3f}K")
    
    def _collect_rls_data(self, x_current: np.ndarray, u_current: float, x_next: np.ndarray):
        """Collect data for Recursive Least Squares (RLS) adaptation."""
        self.rls_data_buffer['x'].append(x_current.copy())
        self.rls_data_buffer['u'].append(u_current)
        self.rls_data_buffer['x_next'].append(x_next.copy())
        
        # Keep buffer at fixed window size
        if len(self.rls_data_buffer['x']) > self.adaptation_window:
            self.rls_data_buffer['x'].pop(0)
            self.rls_data_buffer['u'].pop(0)
            self.rls_data_buffer['x_next'].pop(0)
    
    def _perform_rls_update(self):
        """Perform Recursive Least Squares update of system parameters (2D state system)."""
        try:
            # Need enough data
            N = len(self.rls_data_buffer['x'])
            if N < 50:  # Minimum data for reliable estimation
                return False
            
            # Stack data (2D state system)
            X = np.array(self.rls_data_buffer['x'])  # (N, 2, 1)
            U = np.array(self.rls_data_buffer['u'])  # (N,)
            X_next = np.array(self.rls_data_buffer['x_next'])  # (N, 2, 1)
            
            # Reshape for least squares: x[k+1] = A_d @ x[k] + B_d @ u[k]
            X_curr = X.reshape(N, 2).T  # (2, N)
            X_next_mat = X_next.reshape(N, 2).T  # (2, N)
            U_vec = U.reshape(1, N)  # (1, N)
            
            # Build regressor matrix: Z = [X_curr; U_vec]
            Z = np.vstack([X_curr, U_vec])  # (3, N)
            
            # Regularized least squares: [A_d, B_d] = X_next @ Z^T @ (Z @ Z^T + λI)^{-1}
            ZZT = Z @ Z.T
            ZZT_reg = ZZT + self.rls_regularization * np.eye(3)
            AB_new = X_next_mat @ Z.T @ np.linalg.inv(ZZT_reg)
            
            A_d_new = AB_new[:, :2]  # (2, 2)
            B_d_new = AB_new[:, 2:3]  # (2, 1)
            
            # Compute fit quality (R²)
            X_pred = AB_new @ Z
            residual = X_next_mat - X_pred
            ss_res = np.sum(residual**2, axis=1)
            ss_tot = np.sum((X_next_mat - X_next_mat.mean(axis=1, keepdims=True))**2, axis=1)
            r_squared = 1 - ss_res / (ss_tot + 1e-10)
            
            # Check if model improved
            improved = (r_squared[0] > self.last_r_squared[0] or 
                       r_squared[1] > self.last_r_squared[1])
            
            # Check stability: all eigenvalues must have magnitude < 1
            is_stable = self._check_model_stability(A_d_new)
            
            if improved and is_stable:
                # Accept update
                self.A_d = A_d_new
                self.B_d = B_d_new
                self.last_r_squared = r_squared.tolist()
                self.num_adaptations += 1
                
                # Re-design controller with new model
                self._design_lqr_controller()
                
                print(f"\n[Adaptive LQG] 2D Model updated (#{self.num_adaptations}):")
                print(f"  R² = [{r_squared[0]:.4f}, {r_squared[1]:.4f}]")
                print(f"  Eigenvalues: {np.linalg.eigvals(A_d_new)}")
                
                return True
            else:
                if not is_stable:
                    print(f"[Adaptive] Rejected unstable model (max eigenvalue: {np.max(np.abs(np.linalg.eigvals(A_d_new))):.3f})")
                return False
                
        except Exception as e:
            print(f"Warning: RLS update failed: {e}")
            return False
    
    def _check_model_stability(self, A_d: np.ndarray) -> bool:
        """Check if discrete-time system is stable (all eigenvalues inside unit circle)."""
        try:
            eigs = np.linalg.eigvals(A_d)
            max_eig_mag = np.max(np.abs(eigs))
            return max_eig_mag < 0.99  # Conservative threshold
        except:
            return False
    
    def _load_system_parameters(self):
        """Load system parameters from file (2D state system)."""
        print(f"Loading system parameters from: {self.system_id_file}")
        try:
            with open(self.system_id_file, 'r') as f:
                params = json.load(f)
            
            self.A_d = np.array(params['system_matrices']['A_discrete'])
            B_d_flat = np.array(params['system_matrices']['B_discrete'])
            
            # Auto-detect dimensions and reshape appropriately
            if self.A_d.shape == (2, 2):
                # New format: 2D state system
                self.B_d = B_d_flat.reshape(2, 1)
                print(f"✓ Loaded LQG 2D state system")
            elif self.A_d.shape == (3, 3):
                # Old format: convert 3D to 2D by dropping bath state
                print(f"⚠ Warning: Loading old 3D state system file")
                print(f"  Converting to 2D LQG format (dropping bath state)")
                self.A_d = self.A_d[:2, :2]  # Keep only signal and hot
                self.B_d = B_d_flat[:2].reshape(2, 1) if len(B_d_flat) >= 2 else B_d_flat.reshape(2, 1)
            else:
                raise ValueError(f"Unexpected A_d shape: {self.A_d.shape}")
            
            print(f"Loaded A_d (2x2):")
            print(f"{self.A_d}")
            print(f"Loaded B_d (2x1):")
            print(f"{self.B_d.flatten()}")
        
        except Exception as e:
            raise RuntimeError(f"Could not load system parameters: {e}")
    
    def _set_manual_parameters(self):
        """Set system parameters manually from dict (2D state system)."""
        if 'A_discrete' not in self.system_params or 'B_discrete' not in self.system_params:
            raise ValueError("system_params must contain 'A_discrete' and 'B_discrete'")
        
        self.A_d = np.array(self.system_params['A_discrete'])
        
        # Ensure 2D state system
        if self.A_d.shape != (2, 2):
            raise ValueError(f"Expected A_discrete to be (2, 2), got {self.A_d.shape}. "
                           "Use 2D state system: [signal, hot]")
        self.B_d = np.array(self.system_params['B_discrete']).reshape(2, 1)
        
        print(f"Using manual system parameters (2D LQG):")
        print(f"A_d (2x2) = {self.A_d}")
        print(f"B_d (2x1) = {self.B_d.flatten()}")
    
    def _initialize_from_physics(self):
        """Initialize system matrices from physics-based guesses (2D LQG).
        
        For glassy systems with slow drift:
        - τ_signal ~ 100-1000 ps (slow glassy relaxation)
        - τ_hot ~ 10-50 ps (faster harmonic equilibration)
        - Weak coupling between subsystems
        - Bath (control input) has significant authority on both states
        
        LQG formulation: 2D state [signal, hot], 1D control [bath]
        This allows the controller to start with reasonable parameters
        and let RLS refine them online from control actions.
        """
        # Physics-based time constants (continuous time)
        tau_signal = 200.0  # ps - slow glassy drift
        tau_hot = 30.0      # ps - faster harmonic
        coupling = 0.01     # Weak signal-hot coupling coefficient
        
        # Convert to rates (1/time constant)
        k_s = 1.0 / tau_signal
        k_h = 1.0 / tau_hot
        
        # Continuous-time 2D state system:
        # dx_s/dt = -k_s * x_s + coupling * x_h + b_s * u_bath
        # dx_h/dt = coupling * x_s - k_h * x_h + b_h * u_bath
        #
        # where b_s, b_h are bath input gains (to be learned by RLS)
        
        # Initial guess: signal and hot ARE influenced by bath control
        b_s = k_s * 0.5  # Signal has moderate bath sensitivity
        b_h = k_h * 0.8  # Hot has stronger bath sensitivity (faster response)
        
        A_cont = np.array([
            [-k_s, coupling],
            [coupling, -k_h]
        ])  # (2, 2)
        
        B_cont = np.array([[b_s], [b_h]])  # (2, 1)
        
        # Discretize using zero-order hold (first-order approximation)
        dt = self.update_interval_ps
        
        from scipy.linalg import expm
        # Exact discretization using matrix exponential
        n = A_cont.shape[0]
        M = np.zeros((n + 1, n + 1))
        M[:n, :n] = A_cont * dt
        M[:n, n:n+1] = B_cont * dt
        M_exp = expm(M)
        
        self.A_d = M_exp[:n, :n]  # (2, 2)
        self.B_d = M_exp[:n, n:n+1]  # (2, 1)
        
        print(f"\nPhysics-Based Initialization (2D LQG):")
        print(f"  Time constants:")
        print(f"    τ_signal = {tau_signal:.1f} ps (slow glassy)")
        print(f"    τ_hot = {tau_hot:.1f} ps (harmonic)")
        print(f"  Coupling = {coupling:.4f}")
        print(f"  Bath input gains (initial guess):")
        print(f"    b_s = {b_s:.6f} (signal ← bath)")
        print(f"    b_h = {b_h:.6f} (hot ← bath)")
        print(f"\n  Discrete-time matrices (dt={dt:.2f}ps):")
        print(f"  A_d =")
        for row in self.A_d:
            print(f"    [{', '.join([f'{x:+.6f}' for x in row])}]")
        print(f"  B_d = [{', '.join([f'{x:+.6f}' for x in self.B_d.flatten()])}]")
        print(f"\n  ✓ B ≠ 0: Controller has control authority!")
        print(f"  ✓ RLS will refine these values online from control actions\n")



class AdaptiveLQIController(hoomd.custom.Action):
    """
    Adaptive LQI Controller with mode-based two-output regulation.
    
    Regulates BOTH T_signal and T_hot to target temperature using a single
    bath actuator through mode-based control and online parameter learning.
    
    Theory
    ------
    State vector (7D):
        e = T_signal - T_target      (signal error)
        z = T_hot - T_target          (hot error)
        u = T_bath - T_target         (bath deviation)
        η_c = ∫(x_c dt)               (common mode integral)
        η_d = ∫(x_d dt)               (difference mode integral)
        d_L = drift in signal         (random walk)
        d_H = drift in hot            (random walk)
    
    Mode coordinates:
        x_c = (e + z)/2               (common/average mode)
        x_d = (e - z)/2               (difference mode)
    
    Parameters (3D):
        θ = [a_L, a_H, k]             (rates and coupling)
        a_L = 1/τ_L                   (signal time constant)
        a_H = 1/τ_H                   (hot time constant)
        k = coupling coefficient
    
    Dynamics:
        ė = -(a_L + k)*e + k*z + a_L*u + d_L
        ż = k*e - (a_H + k)*z + a_H*u + d_H
        u̇ = -(1/τ_b)*u + (1/τ_b)*v
        η_c_dot = x_c
        η_d_dot = x_d
        d_L_dot = w_L  (process noise)
        d_H_dot = w_H  (process noise)
    
    Features
    --------
    - **Two-output regulation** with one actuator
    - **Mode-based control**: Common mode drives both temps, difference mode keeps them together
    - **Feedforward compensation**: u_ff from drift estimates
    - **Extended Kalman Filter**: 7-state estimation with drift tracking
    - **Constrained RLS**: Parameter estimation with physical bounds
    - **Adaptive relinearization**: Re-solve DARE when θ changes significantly
    - **Anti-windup**: Back-calculation on both integrals
    - **Rate limiting**: Prevent sudden bath changes
    
    References
    ----------
    [1] Åström & Murray - Feedback Systems (Ch. 12: Frequency Domain Design)
    [2] Anderson & Moore - Optimal Control (Ch. 3: LQR with Integral Action)
    """
    
    def __init__(self,
                 signal_temperature_method: str,
                 hot_temperature_method: str,
                 time_tracker,
                 energy_tracker,
                 simulation=None,
                 molecular_thermostat=None,
                 cavity_thermostat=None,
                 temperature_tracker=None,
                 target_temperature: float = 300.0,
                 dynamic_target: bool = True,
                 dynamic_target_method: Optional[str] = None,
                 # Initial parameter guesses
                 tau_L_initial: float = 200.0,  # ps - slow glassy
                 tau_H_initial: float = 30.0,   # ps - harmonic
                 k_initial: float = 0.01,        # Weak coupling
                 tau_b: float = 1.0,             # ps - bath actuator
                 # LQI weights (mode-based)
                 q_common: float = 1000.0,       # Common mode state
                 q_diff: float = 100.0,          # Difference mode state
                 q_bath: float = 1.0,            # Bath state
                 q_eta_common: float = 0.0,      # Common integral (DISABLED - testing PD control)
                 q_eta_diff: float = 0.0,        # Difference integral (DISABLED)
                 control_effort: float = 100.0,  # Control penalty
                 # Kalman filter noise
                 process_noise_signal: float = 1.0,
                 process_noise_hot: float = 2.0,
                 process_noise_drift: float = 0.01,   # Drift process noise
                 measurement_noise_signal: float = 2.0,
                 measurement_noise_hot: float = 5.0,
                 # RLS parameters
                 rls_forgetting_factor: float = 0.998,
                 rls_regularization: float = 1e-6,
                 rls_update_interval: int = 50,      # Update every N samples
                 # Control limits and safety
                 max_control_rate: float = 10.0,     # K/ps
                 integral_max_common: float = 1000.0,
                 integral_max_diff: float = 100.0,
                 T_min: float = 30.0,
                 T_max: float = 300.0,
                 # Relinearization
                 theta_change_threshold: float = 0.05,  # 5% change triggers re-solve
                 # Timing and output
                 turn_on_time_ps: float = 0.0,
                 turn_off_time_ps: Optional[float] = None,
                 update_interval_ps: float = 1.0,
                 apply_to: str = 'both',
                 output_file: str = 'adaptive_lqi_controller.csv',
                 params_file: str = 'adaptive_lqi_params.json',
                 empirical_data_file: Optional[str] = None,
                 console_output_period_ps: float = 1.0):
        
        super().__init__()
        
        if not HAS_SCIPY:
            raise RuntimeError("Adaptive LQI controller requires scipy for DARE solver")
        
        # Store configuration
        self.signal_temperature_method = signal_temperature_method
        self.hot_temperature_method = hot_temperature_method
        self.time_tracker = time_tracker
        self.energy_tracker = energy_tracker
        self.temperature_tracker = temperature_tracker
        self.simulation = simulation
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        
        # Target temperature
        self.target_temperature = target_temperature
        self.dynamic_target = dynamic_target
        self.dynamic_target_method = dynamic_target_method or signal_temperature_method
        self.target_set = not dynamic_target
        
        # Parameters θ = [a_L, a_H, k]
        self.theta = np.array([1.0/tau_L_initial, 1.0/tau_H_initial, k_initial])
        self.theta_prev = self.theta.copy()
        self.tau_b = tau_b
        self.theta_change_threshold = theta_change_threshold
        
        # LQI weights
        self.q_common = q_common
        self.q_diff = q_diff
        self.q_bath = q_bath
        self.q_eta_common = q_eta_common
        self.q_eta_diff = q_eta_diff
        self.control_effort = control_effort
        
        # Noise parameters
        self.process_noise_signal = process_noise_signal
        self.process_noise_hot = process_noise_hot
        self.process_noise_drift = process_noise_drift
        self.measurement_noise_signal = measurement_noise_signal
        self.measurement_noise_hot = measurement_noise_hot
        
        # State vector [e, z, u, η_c, η_d, d_L, d_H] (7x1)
        self.x_hat = np.zeros((7, 1))
        self.P = np.eye(7) * 1.0  # EKF covariance
        
        # Control gains (to be computed)
        self.K = None      # State feedback gain (1x7)
        self.u_ff = 0.0    # Feedforward compensation
        
        # RLS for parameter estimation
        self.rls_forgetting = rls_forgetting_factor
        self.rls_regularization = rls_regularization
        self.rls_update_interval = rls_update_interval
        self.rls_counter = 0
        self.rls_data_buffer = {'e': [], 'z': [], 'u': [], 'e_next': [], 'z_next': []}
        self.theta_cov = np.eye(3) * 100.0  # Large initial uncertainty
        
        # Control limits and safety
        self.max_control_rate = max_control_rate
        self.integral_max_common = integral_max_common
        self.integral_max_diff = integral_max_diff
        self.T_min = T_min
        self.T_max = T_max
        self.u_prev_temp = None
        self.is_saturated = False  # Track saturation for anti-windup
        
        # Timing
        self.turn_on_time_ps = turn_on_time_ps
        self.turn_off_time_ps = turn_off_time_ps
        self.update_interval_ps = update_interval_ps
        self.last_update_time = -float('inf')
        self.last_console_output = 0.0
        self.console_output_period_ps = console_output_period_ps
        
        # Apply to
        self.apply_to = apply_to
        
        # Output files
        self.output_file = output_file
        self.params_file = params_file
        
        # Empirical data
        self.empirical_data = None
        if empirical_data_file:
            from .analysis import EmpiricalTemperatureData
            self.empirical_data = EmpiricalTemperatureData(empirical_data_file)
        
        # Controller state
        self.is_active = False
        self.phase = 'control'  # Start immediately in control phase
        self.num_updates = 0
        self.num_adaptations = 0
        
        # Initialize controller
        print("\n" + "="*70)
        print("ADAPTIVE LQI CONTROLLER - Mode-Based Two-Output Regulation")
        print("="*70)
        self._initialize_controller()
        print("="*70 + "\n")
        
        # Initialize output file
        self._initialize_output_file()
    
    def _initialize_controller(self):
        """Initialize controller with physics-based parameters."""
        a_L, a_H, k = self.theta
        tau_L = 1.0 / a_L
        tau_H = 1.0 / a_H
        
        print(f"Initial Parameters:")
        print(f"  τ_L = {tau_L:.1f} ps (signal time constant)")
        print(f"  τ_H = {tau_H:.1f} ps (hot time constant)")
        print(f"  k = {k:.4f} (coupling coefficient)")
        print(f"  τ_b = {self.tau_b:.1f} ps (bath actuator)")
        
        print(f"\nMode-Based LQI Weights:")
        print(f"  q_common = {self.q_common:.0f} (average mode)")
        print(f"  q_diff = {self.q_diff:.0f} (difference mode)")
        print(f"  q_η_common = {self.q_eta_common:.0e} (STRONG integral on average)")
        print(f"  q_η_diff = {self.q_eta_diff:.0f} (weak integral on difference)")
        print(f"  R = {self.control_effort:.0f} (control penalty)")
        
        # Compute initial LQI gains
        self._compute_lqi_gains()
        
        print(f"\n✓ Controller initialized - ready for two-output regulation")
        print(f"✓ Both T_signal and T_hot will track T_target")
    
    def _compute_lqi_gains(self):
        """Compute LQI gains using current parameter estimates."""
        # Build continuous-time system matrices
        a_L, a_H, k = self.theta
        tau_b_inv = 1.0 / self.tau_b
        
        # State: [e, z, u, η_c, η_d, d_L, d_H]
        # But for control, we only care about first 5: [e, z, u, η_c, η_d]
        
        # Continuous A matrix (5x5 for control-relevant states)
        A_cont = np.array([
            [-(a_L + k), k, a_L, 0, 0],
            [k, -(a_H + k), a_H, 0, 0],
            [0, 0, -tau_b_inv, 0, 0],
            [0.5, 0.5, 0, 0, 0],  # η_c_dot = x_c = (e+z)/2
            [0.5, -0.5, 0, 0, 0]  # η_d_dot = x_d = (e-z)/2
        ])
        
        B_cont = np.array([[0], [0], [tau_b_inv], [0], [0]])
        
        # Discretize
        from scipy.linalg import expm
        dt = self.update_interval_ps
        M = np.block([[A_cont * dt, B_cont * dt],
                      [np.zeros((1, 6))]])
        M_exp = expm(M)
        A_d = M_exp[:5, :5]
        B_d = M_exp[:5, 5:6]
        
        # LQI cost matrices in original coordinates [e, z, u, η_c, η_d]
        Q = np.diag([
            self.q_common + self.q_diff,  # e weight (combines modes)
            self.q_common + self.q_diff,  # z weight (combines modes)
            self.q_bath,
            self.q_eta_common,
            self.q_eta_diff
        ])
        R = np.array([[self.control_effort]])
        
        # Solve DARE
        try:
            P = solve_discrete_are(A_d, B_d, Q, R)
            self.K = np.linalg.inv(R + B_d.T @ P @ B_d) @ (B_d.T @ P @ A_d)
            
            # Pad K to work with full 7-state vector (ignore d_L, d_H)
            self.K_full = np.hstack([self.K, np.zeros((1, 2))])
            
            # Check stability
            A_cl = A_d - B_d @ self.K
            eigenvalues = np.linalg.eigvals(A_cl)
            max_eig = np.max(np.abs(eigenvalues))
            
            if self.num_updates % 20 == 0:  # Print occasionally
                print(f"\n[Adaptive LQI] Gains recomputed (update #{self.num_updates})")
                print(f"  K = [{', '.join([f'{x:.4f}' for x in self.K.flatten()])}]")
                print(f"  Max |λ| = {max_eig:.4f} {'✓ stable' if max_eig < 1.0 else '⚠ UNSTABLE'}")
        
        except Exception as e:
            print(f"[Adaptive LQI] Error computing gains: {e}")
            if self.K is None:
                # Emergency fallback: conservative gains
                self.K = np.array([[0.1, 0.1, 0.05, 0.01, 0.005]])
                self.K_full = np.hstack([self.K, np.zeros((1, 2))])
                print(f"  Using fallback conservative gains")
    
    def _compute_feedforward(self):
        """Compute feedforward from drift estimates.
        
        Finds u_ff such that steady-state (e, z) = (0, 0) given drifts.
        From dynamics at steady state:
            0 = -(a_L + k)*0 + k*0 + a_L*u_ff + d_L  →  u_ff_L = -d_L/a_L
            0 = k*0 - (a_H + k)*0 + a_H*u_ff + d_H  →  u_ff_H = -d_H/a_H
        
        If drifts differ, use least-squares solution.
        """
        d_L = self.x_hat[5, 0]
        d_H = self.x_hat[6, 0]
        a_L, a_H, k = self.theta
        
        if abs(a_L) > 1e-6 and abs(a_H) > 1e-6:
            u_ff_L = -d_L / a_L
            u_ff_H = -d_H / a_H
            # Average (simple least-squares)
            self.u_ff = 0.5 * (u_ff_L + u_ff_H)
        else:
            self.u_ff = 0.0
    
    def act(self, timestep):
        """Execute adaptive LQI control at each timestep."""
        try:
            current_time_ps = self.time_tracker.elapsed_time
            
            # Check if controller should be active
            should_be_active = (current_time_ps >= self.turn_on_time_ps and
                              (self.turn_off_time_ps is None or current_time_ps < self.turn_off_time_ps))
            
            if not should_be_active:
                if self.is_active:
                    self.is_active = False
                    print(f"[Adaptive LQI] Turned OFF at t = {current_time_ps:.2f} ps")
                return
            
            # Activate controller
            if not self.is_active:
                self.is_active = True
                print(f"\n{'='*70}")
                print(f"[Adaptive LQI] Activated at t = {current_time_ps:.2f} ps")
                print(f"{'='*70}")
                
                # Set dynamic target
                if self.dynamic_target and not self.target_set:
                    measured_temp = self._measure_temperature(self.dynamic_target_method, current_time_ps)
                    if measured_temp is not None:
                        self.target_temperature = measured_temp
                        self.target_set = True
                        print(f"  Dynamic target: T* = {self.target_temperature:.2f} K")
            
            # Execute control
            self._execute_control(timestep, current_time_ps)
        
        except Exception as e:
            print(f"[Adaptive LQI] Error: {e}")
            import traceback
            traceback.print_exc()
    
    def _execute_control(self, timestep, current_time_ps):
        """Execute one control cycle."""
        # Only update at specified intervals
        if current_time_ps - self.last_update_time < self.update_interval_ps:
            return
        
        dt_ps = current_time_ps - self.last_update_time if self.last_update_time > 0 else self.update_interval_ps
        self.last_update_time = current_time_ps
        self.num_updates += 1
        
        # Measure temperatures
        T_signal = self._measure_temperature(self.signal_temperature_method, current_time_ps)
        T_hot = self._measure_temperature(self.hot_temperature_method, current_time_ps)
        T_bath = self._get_current_bath_temperature()
        
        # DIAGNOSTIC: Report measurement failures on first few attempts
        if T_signal is None or T_hot is None or T_bath is None:
            if self.num_updates <= 5:  # Only print first 5 failures to avoid spam
                print(f"[Adaptive LQI] WARNING at t={current_time_ps:.1f}ps: Temperature measurement failed")
                print(f"  T_signal ({self.signal_temperature_method}): {T_signal}")
                print(f"  T_hot ({self.hot_temperature_method}): {T_hot}")
                print(f"  T_bath: {T_bath}")
                if self.energy_tracker is None:
                    print(f"  ⚠ energy_tracker is None! Controller cannot measure fictive temperatures.")
                    print(f"  ⚠ Make sure --enable-energy-tracker is set in your run script!")
                if self.molecular_thermostat is None and self.cavity_thermostat is None:
                    print(f"  ⚠ Both thermostats are None! Controller cannot read bath temperature.")
            return
        
        # Error signals
        e_meas = T_signal - self.target_temperature
        z_meas = T_hot - self.target_temperature
        u_meas = T_bath - self.target_temperature
        
        # Mode coordinates
        x_c = 0.5 * (e_meas + z_meas)  # Common mode
        x_d = 0.5 * (e_meas - z_meas)  # Difference mode
        
        # Save for RLS
        e_prev = self.x_hat[0, 0]
        z_prev = self.x_hat[1, 0]
        u_prev = self.x_hat[2, 0]
        
        # EKF Prediction
        self._ekf_predict(dt_ps)
        
        # EKF Update (measurement)
        y_meas = np.array([[e_meas], [z_meas]])
        self._ekf_update(y_meas)
        
        # ANTI-WINDUP: Conditional integration (integrator clamping)
        # Only integrate if NOT currently saturated from previous step
        is_saturated = hasattr(self, 'is_saturated') and self.is_saturated
        
        if not is_saturated:
            # Update integrals normally
            self.x_hat[3, 0] += 0.5 * (self.x_hat[0, 0] + self.x_hat[1, 0]) * dt_ps  # η_c
            self.x_hat[4, 0] += 0.5 * (self.x_hat[0, 0] - self.x_hat[1, 0]) * dt_ps  # η_d
        else:
            # Saturated - FREEZE integrators to prevent windup
            if self.num_updates <= 10:
                print(f"[Anti-Windup] t={current_time_ps:.1f}ps: Integrators FROZEN (saturated)")
        
        # Clamp integrals (safety backup)
        if abs(self.x_hat[3, 0]) > self.integral_max_common:
            self.x_hat[3, 0] = np.sign(self.x_hat[3, 0]) * self.integral_max_common
        if abs(self.x_hat[4, 0]) > self.integral_max_diff:
            self.x_hat[4, 0] = np.sign(self.x_hat[4, 0]) * self.integral_max_diff
        
        # Compute feedforward
        self._compute_feedforward()
        
        # Control law: v = -K @ x_hat + u_ff
        v_feedback = -self.K_full @ self.x_hat
        v_desired = v_feedback[0, 0] + self.u_ff
        u_bath_desired = self.target_temperature + v_desired
        
        # DIAGNOSTIC: Print control calculation on first few steps
        if self.num_updates <= 10:
            print(f"[Adaptive LQI DEBUG] t={current_time_ps:.1f}ps")
            print(f"  x_hat = {self.x_hat.T}")
            print(f"  K_full = {self.K_full}")
            print(f"  v_feedback = {v_feedback[0, 0]:.2f}")
            print(f"  u_ff = {self.u_ff:.2f}")
            print(f"  v_desired = {v_desired:.2f}")
            print(f"  T_target = {self.target_temperature:.2f}")
            print(f"  u_bath_desired = {u_bath_desired:.2f}")
        
        # Rate limiting
        if self.u_prev_temp is not None:
            max_change = self.max_control_rate * dt_ps
            change = u_bath_desired - self.u_prev_temp
            if abs(change) > max_change:
                u_bath_limited = self.u_prev_temp + np.sign(change) * max_change
            else:
                u_bath_limited = u_bath_desired
        else:
            u_bath_limited = u_bath_desired
        
        # Saturation
        u_bath_final = np.clip(u_bath_limited, self.T_min, self.T_max)
        
        # Track saturation state for next iteration (for conditional integration)
        self.is_saturated = (u_bath_final != u_bath_limited)
        
        # Anti-windup back-calculation
        if u_bath_final != u_bath_desired:
            sat_error = u_bath_final - u_bath_desired
            k_back = 0.5
            # Back-calculate on both integrals
            self.x_hat[3, 0] += sat_error * k_back * dt_ps * 0.5  # Common mode
            self.x_hat[4, 0] += sat_error * k_back * dt_ps * 0.1  # Diff mode (weaker)
            # Re-clamp
            if abs(self.x_hat[3, 0]) > self.integral_max_common:
                self.x_hat[3, 0] = np.sign(self.x_hat[3, 0]) * self.integral_max_common
            if abs(self.x_hat[4, 0]) > self.integral_max_diff:
                self.x_hat[4, 0] = np.sign(self.x_hat[4, 0]) * self.integral_max_diff
        
        # Apply control
        self._update_thermostats(u_bath_final)
        self.u_prev_temp = u_bath_final
        
        # Update state estimate
        self.x_hat[2, 0] = u_bath_final - self.target_temperature
        
        # RLS adaptation
        self._collect_rls_data(e_prev, z_prev, u_prev, e_meas, z_meas)
        if self.rls_counter >= self.rls_update_interval:
            self._perform_rls_update()
            self.rls_counter = 0
            # Check if relinearization needed
            theta_change = np.linalg.norm(self.theta - self.theta_prev) / (np.linalg.norm(self.theta_prev) + 1e-6)
            if theta_change > self.theta_change_threshold:
                self._compute_lqi_gains()
                self.theta_prev = self.theta.copy()
        self.rls_counter += 1
        
        # Logging
        self._log_output(current_time_ps, T_signal, T_hot, T_bath, x_c, x_d)
    
    def _ekf_predict(self, dt):
        """EKF prediction step."""
        a_L, a_H, k = self.theta
        tau_b_inv = 1.0 / self.tau_b
        
        # Current state
        e, z, u, eta_c, eta_d, d_L, d_H = self.x_hat.flatten()
        
        # Prediction (Euler for simplicity - could use RK4)
        e_pred = e + dt * (-(a_L + k)*e + k*z + a_L*u + d_L)
        z_pred = z + dt * (k*e - (a_H + k)*z + a_H*u + d_H)
        u_pred = u + dt * (-tau_b_inv * u)  # Will be overwritten by control
        eta_c_pred = eta_c  # Updated separately
        eta_d_pred = eta_d  # Updated separately
        d_L_pred = d_L  # Random walk
        d_H_pred = d_H  # Random walk
        
        self.x_hat = np.array([[e_pred], [z_pred], [u_pred], [eta_c_pred], [eta_d_pred], [d_L_pred], [d_H_pred]])
        
        # Covariance prediction (simplified)
        Q_process = np.diag([
            self.process_noise_signal**2 * dt,
            self.process_noise_hot**2 * dt,
            1e-6 * dt,  # Bath (control input, low noise)
            1e-9 * dt,  # η_c (deterministic integral)
            1e-9 * dt,  # η_d (deterministic integral)
            self.process_noise_drift**2 * dt,  # d_L
            self.process_noise_drift**2 * dt   # d_H
        ])
        self.P = self.P + Q_process
    
    def _ekf_update(self, y_meas):
        """EKF measurement update."""
        # Measurement model: y = [e, z]
        C = np.array([[1, 0, 0, 0, 0, 0, 0],
                      [0, 1, 0, 0, 0, 0, 0]])
        
        R = np.diag([self.measurement_noise_signal**2,
                     self.measurement_noise_hot**2])
        
        # Innovation
        y_pred = C @ self.x_hat
        innovation = y_meas - y_pred
        S = C @ self.P @ C.T + R
        K_ekf = self.P @ C.T @ np.linalg.inv(S)
        
        # Update
        self.x_hat = self.x_hat + K_ekf @ innovation
        self.P = (np.eye(7) - K_ekf @ C) @ self.P
    
    def _collect_rls_data(self, e_prev, z_prev, u_prev, e_curr, z_curr):
        """Collect data for RLS parameter estimation."""
        self.rls_data_buffer['e'].append(e_prev)
        self.rls_data_buffer['z'].append(z_prev)
        self.rls_data_buffer['u'].append(u_prev)
        self.rls_data_buffer['e_next'].append(e_curr)
        self.rls_data_buffer['z_next'].append(z_curr)
        
        # Keep window
        max_window = 200
        if len(self.rls_data_buffer['e']) > max_window:
            for key in self.rls_data_buffer:
                self.rls_data_buffer[key].pop(0)
    
    def _perform_rls_update(self):
        """Perform RLS update for θ = [a_L, a_H, k]."""
        N = len(self.rls_data_buffer['e'])
        if N < 50:
            return
        
        try:
            # Extract data
            e = np.array(self.rls_data_buffer['e'])
            z = np.array(self.rls_data_buffer['z'])
            u = np.array(self.rls_data_buffer['u'])
            e_next = np.array(self.rls_data_buffer['e_next'])
            z_next = np.array(self.rls_data_buffer['z_next'])
            
            dt = self.update_interval_ps
            
            # Approximate derivatives
            e_dot = (e_next - e) / dt
            z_dot = (z_next - z) / dt
            
            # Regressor for e: ė ≈ -(a_L+k)*e + k*z + a_L*u
            # Regressor for z: ż ≈ k*e - (a_H+k)*z + a_H*u
            
            # Stack into linear system: [ė; ż] = Φ @ θ
            # where Φ rows are [(-e+u), 0, (-e+z)] for e equation
            #                  [0, (-z+u), (e-z)] for z equation
            
            Phi = np.zeros((2*N, 3))
            y_vec = np.zeros(2*N)
            
            for i in range(N):
                # e equation
                Phi[2*i, :] = [-e[i] + u[i], 0, -e[i] + z[i]]
                y_vec[2*i] = e_dot[i]
                # z equation  
                Phi[2*i+1, :] = [0, -z[i] + u[i], e[i] - z[i]]
                y_vec[2*i+1] = z_dot[i]
            
            # Regularized least squares
            theta_new = np.linalg.lstsq(
                Phi.T @ Phi + self.rls_regularization * np.eye(3),
                Phi.T @ y_vec,
                rcond=None
            )[0]
            
            # Project to physical bounds
            theta_new = np.clip(theta_new, [0.001, 0.001, 0.0], [1.0, 1.0, 0.5])
            
            # Check if improvement using robust metrics (less sensitive to noise)
            pred_old = Phi @ self.theta
            pred_new = Phi @ theta_new
            
            # 1. Smoothed predictions (reduce noise impact)
            def smooth_1d(data, window=5):
                if len(data) < window:
                    return data
                return np.convolve(data, np.ones(window)/window, mode='same')
            
            y_vec_smooth = smooth_1d(y_vec, window=5)
            pred_old_smooth = smooth_1d(pred_old, window=5)
            pred_new_smooth = smooth_1d(pred_new, window=5)
            
            # 2. Smoothed R²
            ss_res_old = np.sum((y_vec_smooth - pred_old_smooth)**2)
            ss_res_new = np.sum((y_vec_smooth - pred_new_smooth)**2)
            ss_tot = np.sum((y_vec_smooth - y_vec_smooth.mean())**2)
            r2_old = 1 - ss_res_old / (ss_tot + 1e-10)
            r2_new = 1 - ss_res_new / (ss_tot + 1e-10)
            
            # 3. Correlation coefficient
            corr_old = np.corrcoef(y_vec_smooth, pred_old_smooth)[0, 1] if len(y_vec_smooth) > 1 else 0
            corr_new = np.corrcoef(y_vec_smooth, pred_new_smooth)[0, 1] if len(y_vec_smooth) > 1 else 0
            
            # Accept if: (better smoothed R² OR better correlation) AND not much worse on raw metric
            raw_ratio = np.linalg.norm(pred_new - y_vec) / (np.linalg.norm(pred_old - y_vec) + 1e-10)
            improvement = (r2_new > r2_old * 1.02) or (corr_new > corr_old * 1.01)
            not_much_worse = raw_ratio < 1.5  # Don't accept if raw fit degrades too much
            
            if improvement and not_much_worse:
                self.theta = theta_new
                self.num_adaptations += 1
                
                if self.num_adaptations % 10 == 0:
                    tau_L, tau_H, k = 1.0/self.theta[0], 1.0/self.theta[1], self.theta[2]
                    print(f"[Adaptive LQI] Parameters updated (#{self.num_adaptations}): "
                          f"τ_L={tau_L:.1f}ps, τ_H={tau_H:.1f}ps, k={k:.4f} "
                          f"[R²_smooth={r2_new:.3f}, corr={corr_new:.3f}]")
        
        except Exception as e:
            if self.num_adaptations < 5:  # Only warn initially
                print(f"[Adaptive LQI] RLS update failed: {e}")
    
    def _measure_temperature(self, method: str, current_time_ps: float) -> Optional[float]:
        """Measure temperature using specified method."""
        # Try temperature tracker first (preferred - has all fictive temperatures ready)
        if self.temperature_tracker is not None:
            try:
                # Map method names to temperature tracker attributes
                method_map = {
                    'lj_coulombic': 'lj_coulombic_fictive_temperature',
                    'harmonic': 'harmonic_fictive_temperature',
                    'harmonic_equipartition': 'harmonic_equipartition_temperature',
                    'kinetic': 'kinetic_temperature'
                }
                
                attr_name = method_map.get(method)
                if attr_name and hasattr(self.temperature_tracker, attr_name):
                    temp = getattr(self.temperature_tracker, attr_name)
                    if temp is not None and temp > 0:
                        return temp
            except Exception as e:
                if self.num_updates <= 5:
                    print(f"[Adaptive LQI] WARNING: Failed to get {method} from temperature_tracker: {e}")
        
        # Fallback: Try energy tracker with manual conversion
        if self.energy_tracker is not None and self.empirical_data is not None:
            try:
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                if energy_dict is None:
                    return None
                
                # Extract the right energy component
                if method == 'lj_coulombic':
                    energy = energy_dict.get('lj_coulombic', 0.0)
                elif method in ['harmonic', 'harmonic_equipartition']:
                    energy = energy_dict.get('harmonic', 0.0)
                else:
                    return None
                
                if energy <= 0:
                    return None
                
                # Convert energy to temperature using empirical data
                temp = self.empirical_data.calculate_systemic_temperature(energy)
                return temp if temp > 0 else None
                
            except Exception as e:
                if self.num_updates <= 5:
                    print(f"[Adaptive LQI] WARNING: Failed to measure {method} via energy_tracker: {e}")
                return None
        
        # No valid measurement source
        if self.num_updates <= 5:
            print(f"[Adaptive LQI] WARNING: No valid source for {method} temperature")
            print(f"  temperature_tracker: {self.temperature_tracker}")
            print(f"  energy_tracker: {self.energy_tracker}")
            print(f"  empirical_data: {self.empirical_data}")
        return None
    
    def _get_current_bath_temperature(self) -> Optional[float]:
        """Get current bath temperature from thermostats."""
        try:
            from .utils import PhysicalConstants
            
            # Try molecular thermostat first
            if self.molecular_thermostat is not None:
                try:
                    if hasattr(self.molecular_thermostat, 'kT'):
                        kT = self.molecular_thermostat.kT
                        # Handle variant objects
                        if hasattr(kT, 'value'):
                            kT = kT.value
                        elif callable(kT):
                            kT = kT()
                        return kT / PhysicalConstants.KB_HARTREE_PER_K
                except Exception as e:
                    if self.num_updates <= 5:
                        print(f"[Adaptive LQI] Failed to read molecular thermostat: {e}")
            
            # Try cavity thermostat
            if self.cavity_thermostat is not None:
                try:
                    if hasattr(self.cavity_thermostat, 'kT'):
                        kT = self.cavity_thermostat.kT
                        # Handle variant objects
                        if hasattr(kT, 'value'):
                            kT = kT.value
                        elif callable(kT):
                            kT = kT()
                        return kT / PhysicalConstants.KB_HARTREE_PER_K
                except Exception as e:
                    if self.num_updates <= 5:
                        print(f"[Adaptive LQI] Failed to read cavity thermostat: {e}")
            
            if self.num_updates <= 5:
                print(f"[Adaptive LQI] WARNING: No valid thermostat found")
                print(f"  molecular_thermostat: {self.molecular_thermostat}")
                print(f"  cavity_thermostat: {self.cavity_thermostat}")
            return None
            
        except Exception as e:
            if self.num_updates <= 5:
                print(f"[Adaptive LQI] Unexpected error reading bath temperature: {e}")
            return None
    
    def _update_thermostats(self, temperature_K: float):
        """Update thermostat temperatures."""
        from .utils import PhysicalConstants
        kT = temperature_K * PhysicalConstants.KB_HARTREE_PER_K
        
        if self.apply_to in ['both', 'molecular'] and self.molecular_thermostat is not None:
            try:
                if hasattr(self.molecular_thermostat, 'kT'):
                    if hasattr(self.molecular_thermostat.kT, 'value'):
                        self.molecular_thermostat.kT.value = kT
                    else:
                        self.molecular_thermostat.kT = kT
            except Exception as e:
                if self.num_updates <= 5:
                    print(f"[Adaptive LQI] Failed to update molecular thermostat: {e}")
        
        if self.apply_to in ['both', 'cavity'] and self.cavity_thermostat is not None:
            try:
                if hasattr(self.cavity_thermostat, 'kT'):
                    if hasattr(self.cavity_thermostat.kT, 'value'):
                        self.cavity_thermostat.kT.value = kT
                    else:
                        self.cavity_thermostat.kT = kT
            except Exception as e:
                if self.num_updates <= 5:
                    print(f"[Adaptive LQI] Failed to update cavity thermostat: {e}")
    
    def _initialize_output_file(self):
        """Initialize CSV output file."""
        try:
            with open(self.output_file, 'w') as f:
                f.write("# Adaptive LQI Controller Output\n")
                f.write("# Mode-based two-output regulation\n")
                f.write("time_ps,T_signal_K,T_hot_K,T_bath_K,x_c,x_d,eta_c,eta_d,u_ff,control_v\n")
        except Exception as e:
            print(f"Warning: Could not initialize output file: {e}")
    
    def _log_output(self, time_ps, T_signal, T_hot, T_bath, x_c, x_d):
        """Log control data."""
        try:
            with open(self.output_file, 'a') as f:
                eta_c = self.x_hat[3, 0]
                eta_d = self.x_hat[4, 0]
                u_dev = self.x_hat[2, 0]
                f.write(f"{time_ps:.3f},{T_signal:.4f},{T_hot:.4f},{T_bath:.4f},"
                       f"{x_c:.4f},{x_d:.4f},{eta_c:.4f},{eta_d:.4f},{self.u_ff:.4f},{u_dev:.4f}\n")
            
            # Console output
            if time_ps - self.last_console_output >= self.console_output_period_ps:
                tau_L, tau_H = 1.0/self.theta[0], 1.0/self.theta[1]
                print(f"[Adaptive LQI] t={time_ps:.1f}ps | "
                      f"T_s={T_signal:.2f}K T_h={T_hot:.2f}K T_bath={T_bath:.2f}K | "
                      f"x_c={x_c:.2f}K x_d={x_d:.2f}K | "
                      f"τ_L={tau_L:.0f}ps τ_H={tau_H:.0f}ps")
                self.last_console_output = time_ps
        except:
            pass

