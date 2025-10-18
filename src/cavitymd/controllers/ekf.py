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

from ..utils import PhysicalConstants


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


