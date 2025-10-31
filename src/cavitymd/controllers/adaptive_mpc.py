"""
Adaptive Model Predictive Control (MPC) Controller

This module implements an adaptive MPC controller that simultaneously controls
coupling strength (λ) and bath temperature (T_bath) to regulate lj_coulombic
temperature towards a target temperature T_star.

The controller operates in three phases:
1. Pre-initialization: No control action
2. System identification: Apply randomized step changes to identify linear model
3. Control: Run MPC with online model adaptation via RLS

Author: AI Assistant
Date: 2025-10-30
"""

import numpy as np
import hoomd
from typing import Optional, List
import time as wall_time

try:
    import cvxpy as cp
    CVXPY_AVAILABLE = True
except ImportError:
    CVXPY_AVAILABLE = False
    print("WARNING: cvxpy not available. AdaptiveMPCController will not work.")
    print("Install with: pip install cvxpy")

from ..utils import PhysicalConstants


class AdaptiveMPCController(hoomd.custom.Action):
    """
    Adaptive Model Predictive Control controller for temperature regulation.
    
    Controls both coupling strength λ and bath temperature T_bath to track
    a target temperature T_star by regulating the lj_coulombic temperature.
    
    Uses a linear state-space model identified via randomized step changes
    with online adaptation via Recursive Least Squares (RLS).
    """
    
    def __init__(self,
                 simulation,
                 time_tracker,
                 temperature_tracker,
                 target_temperature: float = 100.0,
                 turn_on_time_ps: float = 0.0,
                 turn_off_time_ps: Optional[float] = None,
                 system_id_duration_ps: float = 50.0,
                 update_interval_ps: float = 0.1,
                 prediction_horizon: int = 10,
                 control_horizon: int = 5,
                 output_weight: float = 100.0,
                 control_effort_weight: List[float] = None,
                 rate_penalty_weight: List[float] = None,
                 lambda_min: float = 0.0,
                 lambda_max: float = 1e-2,
                 T_bath_min: float = 0.1,
                 T_bath_max: float = 500.0,
                 delta_lambda_max: float = 1e-4,
                 delta_T_bath_max: float = 10.0,
                 apply_to: str = 'both',
                 system_id_step_duration_ps: float = 5.0,
                 system_id_seed: int = 42,
                 rls_forgetting_factor: float = 0.995,
                 rls_initial_covariance: float = 100.0,
                 model_update_interval: int = 10,
                 empirical_data_file: Optional[str] = None,
                 output_file: str = 'adaptive_mpc_control.csv',
                 console_output_period_ps: float = 1.0,
                 regularization_param: float = 1e-3,
                 use_scaling: bool = True,
                 debug_mode: bool = False):
        """
        Initialize Adaptive MPC Controller.
        
        Parameters
        ----------
        simulation : CavityMDSimulation
            Reference to simulation object
        time_tracker : ElapsedTimeTracker
            Time tracking object
        temperature_tracker : TemperatureTracker
            Temperature tracking object
        target_temperature : float
            Target temperature T_star in K
        turn_on_time_ps : float
            Time to start system ID phase in ps
        turn_off_time_ps : float, optional
            Time to stop MPC (None = never stop)
        system_id_duration_ps : float
            Duration of system ID phase in ps
        update_interval_ps : float
            MPC update period in ps
        prediction_horizon : int
            Prediction horizon Np
        control_horizon : int
            Control horizon Nc (must be <= Np)
        output_weight : float
            Output tracking weight Q
        control_effort_weight : list of float
            Control effort weights [R_lambda, R_Tbath]
        rate_penalty_weight : list of float
            Rate-of-change penalties [S_lambda, S_Tbath]
        lambda_min, lambda_max : float
            Coupling strength bounds
        T_bath_min, T_bath_max : float
            Bath temperature bounds in K
        delta_lambda_max : float
            Maximum rate of λ change per update
        delta_T_bath_max : float
            Maximum rate of T_bath change in K per update
        apply_to : str
            Apply to 'molecular', 'cavity', or 'both'
        system_id_step_duration_ps : float
            Duration of each step change in ps
        system_id_seed : int
            Random seed for step change generation
        rls_forgetting_factor : float
            RLS forgetting factor (0.9-1.0)
        rls_initial_covariance : float
            Initial RLS covariance scaling
        model_update_interval : int
            Update model every N steps
        empirical_data_file : str, optional
            Empirical data file for T_lj calculation
        output_file : str
            CSV output file path
        console_output_period_ps : float
            Console logging period in ps
        regularization_param : float
            Tikhonov regularization parameter for model identification (default: 1e-3)
        use_scaling : bool
            Enable state/control scaling for better conditioning (default: True)
        debug_mode : bool
            Enable detailed debugging output (default: False)
        """
        if not CVXPY_AVAILABLE:
            raise ImportError("cvxpy is required for AdaptiveMPCController. Install with: pip install cvxpy")
        
        self.simulation = simulation
        self.time_tracker = time_tracker
        self.temperature_tracker = temperature_tracker
        
        # Target and timing
        self.target_temperature = target_temperature
        self.turn_on_time = turn_on_time_ps
        self.turn_off_time = turn_off_time_ps
        self.system_id_duration = system_id_duration_ps
        self.update_interval = update_interval_ps
        
        # MPC parameters
        self.prediction_horizon = prediction_horizon
        self.control_horizon = min(control_horizon, prediction_horizon)
        self.output_weight = output_weight
        self.control_effort_weight = np.array(control_effort_weight if control_effort_weight else [1.0, 0.1])
        self.rate_penalty_weight = np.array(rate_penalty_weight if rate_penalty_weight else [10.0, 1.0])
        
        # Constraints
        self.lambda_min = lambda_min
        self.lambda_max = lambda_max
        self.T_bath_min = T_bath_min
        self.T_bath_max = T_bath_max
        self.delta_lambda_max = delta_lambda_max
        self.delta_T_bath_max = delta_T_bath_max
        
        # Control application
        self.apply_to = apply_to
        
        # System ID parameters
        self.system_id_step_duration = system_id_step_duration_ps
        self.system_id_seed = system_id_seed
        
        # Online learning parameters
        self.rls_forgetting = rls_forgetting_factor
        self.rls_initial_cov = rls_initial_covariance
        self.model_update_interval = model_update_interval
        
        # Model identification parameters
        self.regularization_param = regularization_param
        self.use_scaling = use_scaling
        
        # Scaling factors (computed during system ID)
        self.state_scale = np.ones(3)
        self.state_offset = np.zeros(3)
        self.control_scale = np.ones(2)
        self.control_offset = np.zeros(2)
        
        # Debug mode
        self.debug_mode = debug_mode
        
        # Output
        self.output_file = output_file
        self.console_output_period = console_output_period_ps
        
        # Phase state
        self.phase = 'pre_init'  # 'pre_init', 'system_id', 'control'
        self.phase_start_time = 0.0
        
        # Model: x(k+1) = A*x(k) + B*u(k) + d, y = C*x
        self.model_A = None  # 3x3
        self.model_B = None  # 3x2
        self.model_d = None  # 3x1
        self.model_C = np.array([[1.0, 0.0, 0.0]])  # Fixed: observe T_lj only
        
        # RLS state (one per output dimension)
        self.rls_P = [None, None, None]  # Covariance matrices for each state
        self.theta_vec = [None, None, None]  # Parameter vectors for each state
        
        # System ID data collection
        self.state_history = []
        self.control_history = []
        self.time_history = []
        
        # System ID step sequence
        self.step_sequence_lambda = []
        self.step_sequence_T_bath = []
        self.current_step_index = 0
        
        # Control state
        self.u_prev = None  # Previous control [lambda, T_bath]
        self.x_prev = None  # Previous state [T_lj, T_harmonic, T_kinetic]
        
        # MPC solution cache
        self.mpc_last_solution = None
        
        # Timing
        self.last_update_time = -1.0
        self.last_console_time = -1.0
        self.step_counter = 0
        
        # Statistics
        self.mpc_solve_times = []
        self.mpc_costs = []
        
        # Initialize output file
        self._initialize_output_file()
        
        print(f"AdaptiveMPCController initialized:")
        print(f"  Target temperature: {self.target_temperature:.2f} K")
        print(f"  Turn-on time: {self.turn_on_time:.2f} ps")
        print(f"  System ID duration: {self.system_id_duration:.2f} ps")
        print(f"  Prediction horizon: {self.prediction_horizon}")
        print(f"  Control horizon: {self.control_horizon}")
        print(f"  Lambda range: [{self.lambda_min:.2e}, {self.lambda_max:.2e}]")
        print(f"  T_bath range: [{self.T_bath_min:.2f}, {self.T_bath_max:.2f}] K")
    
    def _initialize_output_file(self):
        """Initialize CSV output file with header."""
        try:
            with open(self.output_file, 'w') as f:
                f.write("time_ps,phase,T_lj,T_harmonic,T_kinetic,T_target,")
                f.write("lambda_cmd,T_bath_cmd,mpc_cost,A_trace,B_norm,solve_time_ms\n")
        except Exception as e:
            raise RuntimeError(f"FATAL: Could not initialize MPC output file '{self.output_file}': {e}") from e
    
    def act(self, timestep):
        """
        Main control loop called every timestep.
        
        Parameters
        ----------
        timestep : int
            Current simulation timestep
        """
        current_time = self.time_tracker.elapsed_time
        
        # Check if we should update (based on update_interval)
        if current_time - self.last_update_time < self.update_interval:
            return
        
        self.last_update_time = current_time
        self.step_counter += 1
        
        # Phase management
        if current_time < self.turn_on_time:
            # Pre-initialization phase: do nothing
            return
        
        elif current_time < self.turn_on_time + self.system_id_duration:
            # System ID phase
            if self.phase == 'pre_init':
                self._enter_system_id_phase(current_time)
            self._system_id_phase(current_time, timestep)
        
        else:
            # Control phase
            if self.phase == 'system_id':
                self._enter_control_phase(current_time)
            
            if self.turn_off_time is None or current_time < self.turn_off_time:
                self._control_phase(current_time, timestep)
    
    def _enter_system_id_phase(self, current_time):
        """Enter system ID phase and generate step sequence."""
        self.phase = 'system_id'
        self.phase_start_time = current_time
        
        # Generate random step sequence
        self._generate_step_sequence()
        
        # Initialize control to mid-range
        lambda_init = 0.5 * (self.lambda_min + self.lambda_max)
        T_bath_init = 0.5 * (self.T_bath_min + self.T_bath_max)
        self.u_prev = np.array([lambda_init, T_bath_init])
        self.x_prev = self._measure_state()
        
        print(f"\n[{current_time:.2f} ps] Entering SYSTEM ID phase")
        print(f"  Generated {len(self.step_sequence_lambda)} step changes")
        print(f"  Step duration: {self.system_id_step_duration:.2f} ps")
    
    def _generate_step_sequence(self):
        """
        Generate random step change sequence for system ID.
        
        Random step changes are limited to reasonable ranges around nominal values:
        - Lambda: ±0.05 around midpoint
        - Temperature: ±20 K around midpoint
        """
        np.random.seed(self.system_id_seed)
        
        # Calculate number of steps
        n_steps = int(self.system_id_duration / self.system_id_step_duration)
        
        # Calculate midpoints
        lambda_mid = 0.5 * (self.lambda_min + self.lambda_max)
        T_bath_mid = 0.5 * (self.T_bath_min + self.T_bath_max)
        
        # Define limited ranges for excitation (to avoid extreme changes)
        # Lambda: ±0.05 around midpoint
        lambda_excitation_range = 0.05
        lambda_low = max(self.lambda_min, lambda_mid - lambda_excitation_range)
        lambda_high = min(self.lambda_max, lambda_mid + lambda_excitation_range)
        
        # Temperature: ±20 K around midpoint
        T_bath_excitation_range = 20.0  # K
        T_bath_low = max(self.T_bath_min, T_bath_mid - T_bath_excitation_range)
        T_bath_high = min(self.T_bath_max, T_bath_mid + T_bath_excitation_range)
        
        # Generate random levels for lambda and T_bath within limited ranges
        self.step_sequence_lambda = np.random.uniform(
            lambda_low, lambda_high, size=n_steps
        )
        self.step_sequence_T_bath = np.random.uniform(
            T_bath_low, T_bath_high, size=n_steps
        )
        
        self.current_step_index = 0
        
        # Report ranges
        print(f"  System ID excitation ranges:")
        print(f"    Lambda: [{lambda_low:.2e}, {lambda_high:.2e}] (±{lambda_excitation_range:.2e})")
        print(f"    T_bath: [{T_bath_low:.1f}, {T_bath_high:.1f}] K (±{T_bath_excitation_range:.1f} K)")
    
    def _system_id_phase(self, current_time, timestep):
        """
        Apply randomized step changes and collect data.
        
        Parameters
        ----------
        current_time : float
            Current simulation time in ps
        timestep : int
            Current timestep
        """
        # Determine which step we're on
        time_in_phase = current_time - self.phase_start_time
        step_index = int(time_in_phase / self.system_id_step_duration)
        step_index = min(step_index, len(self.step_sequence_lambda) - 1)
        
        # Get control commands for this step
        lambda_cmd = self.step_sequence_lambda[step_index]
        T_bath_cmd = self.step_sequence_T_bath[step_index]
        
        # Apply controls
        self._apply_controls(lambda_cmd, T_bath_cmd)
        
        # Measure state
        state = self._measure_state()
        
        # Store data
        self.state_history.append(state.copy())
        self.control_history.append(np.array([lambda_cmd, T_bath_cmd]))
        self.time_history.append(current_time)
        
        # Update previous state
        self.x_prev = state
        self.u_prev = np.array([lambda_cmd, T_bath_cmd])
        
        # Console output
        if current_time - self.last_console_time >= self.console_output_period:
            self.last_console_time = current_time
            print(f"[{current_time:.2f} ps] System ID: step {step_index+1}/{len(self.step_sequence_lambda)}, "
                  f"λ={lambda_cmd:.2e}, T_bath={T_bath_cmd:.1f} K, T_lj={state[0]:.2f} K")
        
        # Log to CSV
        self._log_control_data(current_time, lambda_cmd, T_bath_cmd, 0.0, 0.0)
    
    def _enter_control_phase(self, current_time):
        """Enter control phase and fit initial model."""
        self.phase = 'control'
        self.phase_start_time = current_time
        
        print(f"\n[{current_time:.2f} ps] Entering CONTROL phase")
        print(f"  Collected {len(self.state_history)} data points")
        
        # Fit initial model from collected data
        self._fit_initial_model()
        
        # Initialize RLS
        self._initialize_rls()
        
        print(f"  Model identified:")
        print(f"    A eigenvalues: {np.linalg.eigvals(self.model_A)}")
        print(f"    A trace: {np.trace(self.model_A):.4f}")
        print(f"    B norm: {np.linalg.norm(self.model_B):.4f}")
    
    def _fit_initial_model(self):
        """Fit state-space model from collected data using regularized least squares."""
        if len(self.state_history) < 10:
            print(f"WARNING: Only {len(self.state_history)} data points collected. Using default model.")
            self._use_default_model()
            return
        
        # Convert to arrays
        X = np.array(self.state_history[:-1])  # (N-1) x 3
        U = np.array(self.control_history[:-1])  # (N-1) x 2
        X_next = np.array(self.state_history[1:])  # (N-1) x 3
        
        N = X.shape[0]
        
        # Compute scaling factors if enabled
        if self.use_scaling:
            self.state_offset = np.mean(X, axis=0)
            self.state_scale = np.std(X, axis=0) + 1e-6  # Avoid division by zero
            self.control_offset = np.mean(U, axis=0)
            self.control_scale = np.std(U, axis=0) + 1e-6
            
            # Normalize data
            X_norm = (X - self.state_offset) / self.state_scale
            U_norm = (U - self.control_offset) / self.control_scale
            X_next_norm = (X_next - self.state_offset) / self.state_scale
        else:
            X_norm = X
            U_norm = U
            X_next_norm = X_next
        
        # Augment regressor: Phi = [X_norm, U_norm, ones]
        Phi = np.hstack([X_norm, U_norm, np.ones((N, 1))])  # (N-1) x 6
        
        # Solve regularized least squares: X_next_norm = Phi * Theta^T
        # Add Tikhonov regularization: (Phi^T Phi + lambda*I) Theta = Phi^T X_next_norm
        try:
            PhiT_Phi = Phi.T @ Phi
            PhiT_X = Phi.T @ X_next_norm
            
            # Add regularization to diagonal (exclude bias term)
            reg_matrix = self.regularization_param * np.eye(6)
            reg_matrix[5, 5] = 0  # Don't regularize bias
            
            Theta = np.linalg.solve(PhiT_Phi + reg_matrix, PhiT_X)
            
            # Extract matrices (in normalized space)
            A_norm = Theta[:3, :].T  # 3 x 3
            B_norm = Theta[3:5, :].T  # 3 x 2
            d_norm = Theta[5, :]  # 3 x 1
            
            # Transform back to original space if scaling was used
            if self.use_scaling:
                # Original model: x_next = A*x + B*u + d
                # Normalized: (x_next - x0)/sx = A_norm*(x - x0)/sx + B_norm*(u - u0)/su + d_norm
                # Multiply by sx: x_next - x0 = A_norm*(x - x0) + B_norm*sx/su*(u - u0) + d_norm*sx
                # x_next = A_norm*x + B_norm*sx/su*u + (d_norm*sx - A_norm*x0 + B_norm*sx/su*u0 + x0)
                
                self.model_A = A_norm
                self.model_B = B_norm * (self.state_scale[:, None] / self.control_scale[None, :])
                self.model_d = (d_norm * self.state_scale 
                               - A_norm @ self.state_offset 
                               + self.model_B @ self.control_offset 
                               + self.state_offset)
            else:
                self.model_A = A_norm
                self.model_B = B_norm
                self.model_d = d_norm
            
            # Clip B matrix to prevent extremely large values
            B_norm_before = np.linalg.norm(self.model_B)
            max_B_element = 50.0  # Maximum allowed element in B
            if np.max(np.abs(self.model_B)) > max_B_element:
                self.model_B = np.clip(self.model_B, -max_B_element, max_B_element)
                print(f"  WARNING: B matrix clipped (norm before: {B_norm_before:.2f}, after: {np.linalg.norm(self.model_B):.2f})")
            
            # Check stability
            eigvals = np.linalg.eigvals(self.model_A)
            max_eigval = np.max(np.abs(eigvals))
            
            if max_eigval > 1.05:
                print(f"  WARNING: Identified model is unstable (max |eig(A)| = {max_eigval:.3f} > 1.0)")
                print(f"  Projecting A eigenvalues to unit circle")
                # Project to stable region
                D, V = np.linalg.eig(self.model_A)
                D_stable = np.where(np.abs(D) > 0.99, 0.99 * D / np.abs(D), D)
                self.model_A = (V @ np.diag(D_stable) @ np.linalg.inv(V)).real
            
            # Validate model with data
            predictions = []
            x = X[0]
            for i in range(min(20, N-1)):
                x = self.model_A @ x + self.model_B @ U[i] + self.model_d
                predictions.append(x[0])  # T_lj prediction
            
            actuals = X_next[:min(20, N-1), 0]
            prediction_error = np.mean(np.abs(np.array(predictions) - actuals))
            print(f"  Model validation: Mean prediction error = {prediction_error:.2f} K")
            
            if prediction_error > 50.0:
                print(f"  WARNING: Large prediction error. Model may be poor quality.")
            
        except np.linalg.LinAlgError as e:
            raise RuntimeError(f"FATAL: System identification failed - least squares solver error: {e}") from e
    
    def _use_default_model(self):
        """Use a simple default model if identification fails."""
        # Default: weak coupling, diagonal A
        self.model_A = np.diag([0.95, 0.95, 0.95])
        self.model_B = np.array([[0.01, 0.02],
                                  [0.005, 0.01],
                                  [0.005, 0.03]])
        self.model_d = np.zeros(3)
    
    def _initialize_rls(self):
        """Initialize RLS covariance matrices and parameter vectors."""
        # Flatten model parameters into theta vectors (one per state dimension)
        # theta[i] = [A[i,:], B[i,:], d[i]]^T (6 x 1 for each i=0,1,2)
        
        for i in range(3):
            self.theta_vec[i] = np.concatenate([
                self.model_A[i, :],  # A row i
                self.model_B[i, :],  # B row i
                [self.model_d[i]]    # d element i
            ])
            self.rls_P[i] = self.rls_initial_cov * np.eye(6)
    
    def _control_phase(self, current_time, timestep):
        """
        Run MPC optimization and apply control.
        
        Parameters
        ----------
        current_time : float
            Current simulation time in ps
        timestep : int
            Current timestep
        """
        # Update model online via RLS
        if self.step_counter % self.model_update_interval == 0:
            self._update_model_online()
        
        # Solve MPC optimization
        start_time = wall_time.time()
        u_opt, cost = self._solve_mpc(debug=self.debug_mode)
        solve_time_ms = (wall_time.time() - start_time) * 1000.0
        
        self.mpc_solve_times.append(solve_time_ms)
        self.mpc_costs.append(cost)
        
        # Apply first control move
        lambda_cmd = u_opt[0, 0]
        T_bath_cmd = u_opt[1, 0]
        self._apply_controls(lambda_cmd, T_bath_cmd)
        
        # Update state
        self.x_prev = self._measure_state()
        self.u_prev = np.array([lambda_cmd, T_bath_cmd])
        
        # Console output
        if current_time - self.last_console_time >= self.console_output_period:
            self.last_console_time = current_time
            T_lj = self.x_prev[0]
            error = T_lj - self.target_temperature
            print(f"[{current_time:.2f} ps] MPC: T_lj={T_lj:.2f} K (error={error:+.2f} K), "
                  f"λ={lambda_cmd:.2e}, T_bath={T_bath_cmd:.1f} K, cost={cost:.2e}, "
                  f"solve={solve_time_ms:.2f} ms")
        
        # Log to CSV
        self._log_control_data(current_time, lambda_cmd, T_bath_cmd, cost, solve_time_ms)
    
    def _update_model_online(self):
        """Update model parameters using RLS."""
        if self.x_prev is None or self.u_prev is None:
            return
        
        # Get current state
        x_current = self._measure_state()
        
        # Regressor: phi = [x_prev, u_prev, 1]
        phi = np.concatenate([self.x_prev, self.u_prev, [1.0]])  # 6 x 1
        
        # Update each state dimension independently
        for i in range(3):
            # Prediction error
            y_pred = phi @ self.theta_vec[i]
            y_actual = x_current[i]
            error = y_actual - y_pred
            
            # RLS gain
            P_phi = self.rls_P[i] @ phi
            denominator = self.rls_forgetting + phi @ P_phi
            K = P_phi / denominator
            
            # Parameter update
            self.theta_vec[i] += K * error
            
            # Covariance update
            self.rls_P[i] = (self.rls_P[i] - np.outer(K, P_phi)) / self.rls_forgetting
        
        # Unpack theta back to A, B, d
        self._unpack_theta_to_model()
    
    def _unpack_theta_to_model(self):
        """Unpack theta_vec back to model matrices A, B, d."""
        for i in range(3):
            self.model_A[i, :] = self.theta_vec[i][:3]
            self.model_B[i, :] = self.theta_vec[i][3:5]
            self.model_d[i] = self.theta_vec[i][5]
    
    def _solve_mpc(self, debug=False):
        """
        Solve MPC quadratic program using cvxpy.
        
        Parameters
        ----------
        debug : bool, optional
            If True, print detailed debug information
        
        Returns
        -------
        U : ndarray (2, Nc)
            Optimal control sequence
        cost : float
            Optimal cost value
        """
        # Get current state
        x0 = self._measure_state()
        
        if debug:
            print("\n" + "="*80)
            print("MPC DEBUG OUTPUT")
            print("="*80)
            print(f"Current state x0 = {x0}")
            print(f"Target temperature = {self.target_temperature:.2f} K")
            print(f"Previous control u_prev = {self.u_prev}")
            print(f"\nModel matrices:")
            print(f"  A =\n{self.model_A}")
            print(f"  B =\n{self.model_B}")
            print(f"  d = {self.model_d}")
            print(f"  C = {self.model_C}")
            print(f"\nWeights:")
            print(f"  output_weight = {self.output_weight}")
            print(f"  control_effort_weight = {self.control_effort_weight}")
            print(f"  rate_penalty_weight = {self.rate_penalty_weight}")
            print(f"\nControl bounds:")
            print(f"  lambda: [{self.lambda_min:.4e}, {self.lambda_max:.4e}]")
            print(f"  T_bath: [{self.T_bath_min:.2f}, {self.T_bath_max:.2f}] K")
            print(f"  delta_lambda_max: {self.delta_lambda_max:.4e}")
            print(f"  delta_T_bath_max: {self.delta_T_bath_max:.2f} K")
        
        # Decision variable: U = [u(0), u(1), ..., u(Nc-1)]
        U = cp.Variable((2, self.control_horizon))
        
        # Build cost function and constraints
        cost = 0
        constraints = []
        
        # Simulate forward
        x_pred = x0
        for i in range(self.prediction_horizon):
            # Control at time i
            if i < self.control_horizon:
                u_i = U[:, i]
            else:
                u_i = U[:, -1]  # Hold last control
            
            # State prediction: x(k+1) = A*x(k) + B*u(k) + d
            x_pred = self.model_A @ x_pred + self.model_B @ u_i + self.model_d
            
            # Output: y = C*x (T_lj only)
            y_pred = self.model_C @ x_pred
            
            # Tracking cost: ||y - T_star||^2 * Q
            cost += self.output_weight * cp.square(y_pred - self.target_temperature)
            
            # Control effort cost (only for control horizon)
            if i < self.control_horizon:
                cost += self.control_effort_weight[0] * cp.square(u_i[0])
                cost += self.control_effort_weight[1] * cp.square(u_i[1])
            
            # Rate-of-change cost
            if i < self.control_horizon - 1:
                delta_u = U[:, i+1] - U[:, i]
                cost += self.rate_penalty_weight[0] * cp.square(delta_u[0])
                cost += self.rate_penalty_weight[1] * cp.square(delta_u[1])
            elif i == 0 and self.u_prev is not None:
                delta_u = U[:, 0] - self.u_prev
                cost += self.rate_penalty_weight[0] * cp.square(delta_u[0])
                cost += self.rate_penalty_weight[1] * cp.square(delta_u[1])
        
        # Constraints on control inputs
        for i in range(self.control_horizon):
            # Box constraints
            constraints.append(U[0, i] >= self.lambda_min)
            constraints.append(U[0, i] <= self.lambda_max)
            constraints.append(U[1, i] >= self.T_bath_min)
            constraints.append(U[1, i] <= self.T_bath_max)
            
            # Rate constraints
            if i > 0:
                constraints.append(cp.abs(U[0, i] - U[0, i-1]) <= self.delta_lambda_max)
                constraints.append(cp.abs(U[1, i] - U[1, i-1]) <= self.delta_T_bath_max)
            elif self.u_prev is not None:
                constraints.append(cp.abs(U[0, 0] - self.u_prev[0]) <= self.delta_lambda_max)
                constraints.append(cp.abs(U[1, 0] - self.u_prev[1]) <= self.delta_T_bath_max)
        
        if debug:
            print(f"\nProblem setup:")
            print(f"  Number of constraints: {len(constraints)}")
            print(f"  Decision variable shape: U = (2, {self.control_horizon})")
            print(f"  Prediction horizon: {self.prediction_horizon}")
        
        # Formulate and solve problem
        problem = cp.Problem(cp.Minimize(cost), constraints)
        
        # Warm start if we have a previous solution
        if self.mpc_last_solution is not None:
            U.value = self.mpc_last_solution
            if debug:
                print(f"\nWarm starting with previous solution:")
                print(f"  U_warmstart[:, 0] = {self.mpc_last_solution[:, 0]}")
        
        try:
            if debug:
                print(f"\nSolving with OSQP...")
            
            solve_start = wall_time.time()
            
            # Solve with improved solver settings
            problem.solve(
                solver=cp.OSQP,
                warm_start=True,
                verbose=debug,  # Show solver output if debug=True
                max_iter=20000,           # Increase max iterations
                eps_abs=1e-4,             # Looser absolute tolerance
                eps_rel=1e-4,             # Looser relative tolerance
                polish=True,              # Refine solution
                adaptive_rho=True,        # Adaptive step size
                rho=0.1                   # Initial step size
            )
            
            solve_time = (wall_time.time() - solve_start) * 1000  # ms
            
            if debug:
                print(f"\nSolver results:")
                print(f"  Status: {problem.status}")
                print(f"  Optimal cost: {problem.value}")
                print(f"  Solve time: {solve_time:.2f} ms")
                if U.value is not None:
                    print(f"  Optimal U[:, 0] = {U.value[:, 0]}")
                    print(f"  Optimal U (all) =\n{U.value}")
                else:
                    print(f"  U.value is None!")
                
                # Additional solver info
                try:
                    print(f"\nSolver statistics:")
                    if hasattr(problem, 'solver_stats'):
                        stats = problem.solver_stats
                        print(f"  setup_time: {stats.setup_time}")
                        print(f"  solve_time: {stats.solve_time}")
                        print(f"  num_iters: {stats.num_iters}")
                except Exception as e:
                    print(f"  Could not get solver stats: {e}")
            
            if problem.status in ['optimal', 'optimal_inaccurate']:
                self.mpc_last_solution = U.value.copy()
                if debug:
                    print(f"\n✓ MPC solve successful")
                    print("="*80)
                return U.value, problem.value
            else:
                print(f"\n{'!'*80}")
                print(f"WARNING: MPC solve failed with status: {problem.status}")
                print(f"  Solve time: {solve_time:.2f} ms")
                print(f"  Current state: {x0}")
                print(f"  Target: {self.target_temperature:.2f} K")
                if self.u_prev is not None:
                    print(f"  Previous control: λ={self.u_prev[0]:.4e}, T_bath={self.u_prev[1]:.2f} K")
                print(f"{'!'*80}\n")
                
                # Fallback: hold previous control (small change allowed)
                if self.u_prev is not None:
                    # Allow small decay towards mid-range
                    u_mid = np.array([0.5 * (self.lambda_min + self.lambda_max),
                                      0.5 * (self.T_bath_min + self.T_bath_max)])
                    u_fallback = 0.95 * self.u_prev + 0.05 * u_mid
                    if debug:
                        print(f"Using fallback: u_fallback = {u_fallback.flatten()}")
                    return u_fallback.reshape(2, 1), 0.0
                else:
                    # Default: mid-range
                    u_default = np.array([[0.5 * (self.lambda_min + self.lambda_max)],
                                          [0.5 * (self.T_bath_min + self.T_bath_max)]])
                    if debug:
                        print(f"Using default: u_default = {u_default.flatten()}")
                    return u_default, 0.0
        
        except Exception as e:
            print(f"\n{'!'*80}")
            print(f"ERROR in MPC solve: {e}")
            import traceback
            traceback.print_exc()
            print(f"{'!'*80}\n")
            raise RuntimeError(f"FATAL: MPC optimization failed: {e}") from e
    
    def _measure_state(self):
        """
        Get current state [T_lj, T_harmonic, T_kinetic].
        
        Returns
        -------
        state : ndarray (3,)
            Current state vector
        """
        tt = self.temperature_tracker
        T_lj = tt.lj_coulombic_fictive_temperature if tt.lj_coulombic_fictive_temperature is not None else 0.0
        T_harmonic = tt.harmonic_equipartition_temperature if tt.harmonic_equipartition_temperature is not None else 0.0
        T_kinetic = tt.kinetic_temperature if tt.kinetic_temperature is not None else 0.0
        return np.array([T_lj, T_harmonic, T_kinetic])
    
    def _apply_controls(self, lambda_cmd, T_bath_cmd):
        """
        Apply control commands to simulation.
        
        Parameters
        ----------
        lambda_cmd : float
            Coupling strength command (dimensionless lambda)
        T_bath_cmd : float
            Bath temperature command in K
        """
        # Update cavity coupling strength
        try:
            if hasattr(self.simulation, 'cavity_force'):
                # Update the lambda_coupling_variant with new value
                from hoomd.variant import Constant
                new_coupling = Constant(float(lambda_cmd))
                self.simulation.cavity_force.lambda_coupling_variant = new_coupling
                
                # Update params in the C++ force implementation
                if hasattr(self.simulation.cavity_force, '_force_impl') and self.simulation.cavity_force._force_impl:
                    if hasattr(self.simulation.cavity_force._force_impl, 'setParams'):
                        self.simulation.cavity_force._force_impl.setParams(
                            self.simulation.cavity_force.omegac,
                            new_coupling,
                            self.simulation.cavity_force.phmass
                        )
        except Exception as e:
            raise RuntimeError(f"FATAL: Could not apply coupling strength control: {e}") from e
        
        # Update bath temperatures
        kT_au = T_bath_cmd * PhysicalConstants.KB_HARTREE_PER_K
        
        if self.apply_to in ['molecular', 'both']:
            if hasattr(self.simulation, 'molecular_thermostat_obj') and self.simulation.molecular_thermostat_obj:
                self.simulation.molecular_thermostat_obj.kT = kT_au
        
        if self.apply_to in ['cavity', 'both']:
            if hasattr(self.simulation, 'cavity_thermostat_obj') and self.simulation.cavity_thermostat_obj:
                self.simulation.cavity_thermostat_obj.kT = kT_au
    
    def _log_control_data(self, time, lambda_cmd, T_bath_cmd, cost, solve_time_ms):
        """
        Write control data to CSV file.
        
        Parameters
        ----------
        time : float
            Current time in ps
        lambda_cmd : float
            Applied coupling strength
        T_bath_cmd : float
            Applied bath temperature in K
        cost : float
            MPC cost value
        solve_time_ms : float
            Solve time in milliseconds
        """
        try:
            state = self._measure_state()
            A_trace = np.trace(self.model_A) if self.model_A is not None else 0.0
            B_norm = np.linalg.norm(self.model_B) if self.model_B is not None else 0.0
            
            with open(self.output_file, 'a') as f:
                f.write(f"{time:.6f},{self.phase},{state[0]:.6f},{state[1]:.6f},{state[2]:.6f},")
                f.write(f"{self.target_temperature:.6f},{lambda_cmd:.8e},{T_bath_cmd:.6f},")
                f.write(f"{cost:.8e},{A_trace:.6f},{B_norm:.6f},{solve_time_ms:.4f}\n")
        except Exception as e:
            raise RuntimeError(f"FATAL: Could not write MPC control data to '{self.output_file}': {e}") from e

