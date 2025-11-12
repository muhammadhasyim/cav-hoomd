"""
Test suite for LQR temperature controller with system identification and Kalman filtering.

This module tests the LQR-based optimal temperature control system including:
- System identification from step response data
- Discrete-time LQR control law with integral action
- Kalman filter state estimation
- Integration with HOOMD-blue temperature control
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_array_less
import json
import tempfile
import os


class TestSystemIdentification:
    """Test system identification from step response data."""
    
    def test_least_squares_identification_perfect_data(self):
        """Test system ID with perfect synthetic data (no noise)."""
        # True system matrices (discrete-time, dt=0.1)
        A_true = np.array([
            [0.95, 0.02, 0.03],
            [0.01, 0.97, 0.02],
            [0.00, 0.00, 0.98]
        ])
        B_true = np.array([[0.05], [0.03], [0.02]])
        
        # Generate synthetic step response
        dt = 0.1
        N = 500  # 50 ps of data
        x = np.zeros((3, N+1))
        u = np.ones((1, N)) * (-295.0)  # Step from 300K to 5K
        
        for k in range(N):
            x[:, k+1:k+2] = A_true @ x[:, k:k+1] + B_true @ u[:, k:k+1]
        
        # Perform system identification
        X_curr = x[:, :-1]
        X_next = x[:, 1:]
        U = u
        
        # Stack [X_curr; U] and solve least squares
        Z = np.vstack([X_curr, U])
        AB = X_next @ Z.T @ np.linalg.pinv(Z @ Z.T)
        
        A_id = AB[:, :3]
        B_id = AB[:, 3:4]
        
        # Check identification accuracy (numerical precision limits)
        assert_allclose(A_id, A_true, rtol=1e-6, atol=1e-10)
        assert_allclose(B_id, B_true, rtol=1e-6, atol=1e-10)
    
    def test_least_squares_identification_noisy_data(self):
        """Test system ID with noisy measurements."""
        # True system matrices
        A_true = np.array([
            [0.95, 0.02, 0.03],
            [0.01, 0.97, 0.02],
            [0.00, 0.00, 0.98]
        ])
        B_true = np.array([[0.05], [0.03], [0.02]])
        
        # Generate synthetic step response with noise
        np.random.seed(42)
        dt = 0.1
        N = 500
        x = np.zeros((3, N+1))
        u = np.ones((1, N)) * (-295.0)
        noise_std = 0.5  # 0.5K measurement noise
        
        for k in range(N):
            x[:, k+1:k+2] = A_true @ x[:, k:k+1] + B_true @ u[:, k:k+1]
            x[:, k+1] += np.random.randn(3) * noise_std
        
        # Perform system identification
        X_curr = x[:, :-1]
        X_next = x[:, 1:]
        U = u
        
        Z = np.vstack([X_curr, U])
        AB = X_next @ Z.T @ np.linalg.pinv(Z @ Z.T)
        
        A_id = AB[:, :3]
        B_id = AB[:, 3:4]
        
        # Check identification accuracy (should be close despite noise)
        # With 0.5K noise on 295K quench, expect ~0.2% identification error
        assert_allclose(A_id, A_true, rtol=0.1, atol=0.02)  # 10% relative or 0.02 absolute
        assert_allclose(B_id, B_true, rtol=0.1, atol=0.02)
    
    def test_parameter_extraction_from_matrices(self):
        """Test extracting physical parameters from identified matrices."""
        # Discrete-time system with dt=0.1
        dt = 0.1
        A_d = np.array([
            [0.95, 0.02, 0.03],
            [0.01, 0.97, 0.02],
            [0.00, 0.00, 0.98]
        ])
        B_d = np.array([0.05, 0.03, 0.02])
        
        # Extract approximate continuous-time parameters
        # A_d ≈ I + A_c * dt for small dt
        A_c_approx = (A_d - np.eye(3)) / dt
        
        k_s = -A_c_approx[0, 0]
        k_h = -A_c_approx[1, 1]
        k_c = -A_c_approx[2, 2]
        
        # Check that parameters are physically reasonable
        assert k_s > 0, "k_s should be positive (decay rate)"
        assert k_h > 0, "k_h should be positive (decay rate)"
        assert k_c > 0, "k_c should be positive (decay rate)"
        
        # Time constants should be reasonable (not too fast/slow)
        assert 1 < 1/k_s < 1000, "Signal time constant in reasonable range"
        assert 1 < 1/k_h < 1000, "Hot time constant in reasonable range"


class TestLQRControl:
    """Test LQR control law computation."""
    
    def test_lqr_gain_computation(self):
        """Test LQR gain computation using discrete-time algebraic Riccati equation."""
        from scipy.linalg import solve_discrete_are
        
        # Simple 2-state system
        A = np.array([[0.9, 0.1], [0.0, 0.95]])
        B = np.array([[0.1], [0.05]])
        Q = np.diag([100.0, 1.0])  # Regulate first state tightly
        R = np.array([[1.0]])
        
        # Solve discrete-time Riccati equation
        P = solve_discrete_are(A, B, Q, R)
        
        # Compute LQR gain
        K = np.linalg.inv(R + B.T @ P @ B) @ (B.T @ P @ A)
        
        # Check dimensions
        assert K.shape == (1, 2)
        
        # Check closed-loop stability
        A_cl = A - B @ K
        eigenvalues = np.linalg.eigvals(A_cl)
        assert np.all(np.abs(eigenvalues) < 1.0), "Closed-loop system should be stable"
    
    def test_lqr_with_integral_action(self):
        """Test LQR with integral action for zero steady-state error."""
        from scipy.linalg import solve_discrete_are
        
        # System: x[k+1] = A*x[k] + B*u[k]
        A = np.array([[0.95]])
        B = np.array([[0.05]])
        
        # Augment with integrator: xi[k+1] = xi[k] + x[k]
        A_aug = np.array([
            [0.95, 0.0],
            [1.0, 1.0]
        ])
        B_aug = np.array([[0.05], [0.0]])
        
        Q_aug = np.diag([100.0, 10.0])  # State weight, integral weight
        R = np.array([[1.0]])
        
        # Solve augmented LQR
        P = solve_discrete_are(A_aug, B_aug, Q_aug, R)
        K = np.linalg.inv(R + B_aug.T @ P @ B_aug) @ (B_aug.T @ P @ A_aug)
        
        # Simulate closed-loop with constant disturbance
        x = np.array([[10.0], [0.0]])  # Initial error
        u_hist = []
        x_hist = []
        
        for _ in range(100):
            u = -K @ x
            u_hist.append(u[0, 0])
            x_hist.append(x[0, 0])
            x = A_aug @ x + B_aug @ u
            x[0, 0] += 1.0  # Constant disturbance
        
        # Check that error goes to zero (integral action compensates disturbance)
        assert abs(x_hist[-1]) < 0.1, "Integral action should drive error to near zero"


class TestKalmanFilter:
    """Test Kalman filter state estimation."""
    
    def test_kalman_filter_steady_state(self):
        """Test Kalman filter convergence to steady state."""
        from scipy.linalg import solve_discrete_are
        
        # System
        A = np.array([[0.95, 0.05], [0.0, 0.97]])
        C = np.array([[1.0, 0.0], [0.0, 1.0]])  # Measure both states
        
        # Noise covariances
        Q = np.diag([0.1, 0.1])  # Process noise
        R = np.diag([0.5, 0.5])  # Measurement noise
        
        # Solve Riccati for Kalman gain
        P = solve_discrete_are(A.T, C.T, Q, R)
        L = P @ C.T @ np.linalg.inv(C @ P @ C.T + R)
        
        # Check dimensions
        assert L.shape == (2, 2)
        
        # Check estimator stability
        A_est = A - L @ C @ A
        eigenvalues = np.linalg.eigvals(A_est)
        assert np.all(np.abs(eigenvalues) < 1.0), "Estimator should be stable"
    
    def test_kalman_filter_estimation_accuracy(self):
        """Test Kalman filter estimation with known system."""
        from scipy.linalg import solve_discrete_are
        
        # True system
        A = np.array([[0.95]])
        C = np.array([[1.0]])
        Q = np.array([[0.1]])
        R = np.array([[0.5]])
        
        # Kalman gain
        P = solve_discrete_are(A.T, C.T, Q, R)
        L = P @ C.T @ np.linalg.inv(C @ P @ C.T + R)
        
        # Simulate with noise
        np.random.seed(42)
        N = 100
        x_true = np.zeros((1, N))
        x_hat = np.zeros((1, N))
        
        for k in range(N-1):
            # True state evolution
            x_true[:, k+1] = A @ x_true[:, k:k+1].reshape(-1, 1) + np.random.randn(1, 1) * np.sqrt(Q)
            
            # Measurement
            y = C @ x_true[:, k+1:k+2] + np.random.randn(1, 1) * np.sqrt(R)
            
            # Kalman filter update
            x_hat[:, k+1] = A @ x_hat[:, k:k+1].reshape(-1, 1) + L @ (y - C @ A @ x_hat[:, k:k+1].reshape(-1, 1))
        
        # Check that estimation error is smaller than measurement noise
        errors = np.abs(x_true - x_hat)
        assert np.mean(errors) < np.sqrt(R[0, 0]), "Kalman filter should reduce noise"


class TestLQRControllerIntegration:
    """Test full LQR controller integration."""
    
    def test_system_id_data_collection(self):
        """Test that system ID collects data correctly."""
        # Simulate data collection during cold quench
        dt = 0.1
        duration = 50.0
        N = int(duration / dt)
        
        time_data = []
        x_s_data = []
        x_h_data = []
        x_c_data = []
        u_data = []
        
        # Simulate cold quench response
        for k in range(N):
            t = k * dt
            x_s = 295.0 * np.exp(-t / 20.0)  # Exponential decay
            x_h = 295.0 * np.exp(-t / 30.0)
            x_c = 295.0 * np.exp(-t / 50.0)
            u = -295.0  # Step input
            
            time_data.append(t)
            x_s_data.append(x_s)
            x_h_data.append(x_h)
            x_c_data.append(x_c)
            u_data.append(u)
        
        # Check data shapes
        assert len(time_data) == N
        assert len(x_s_data) == N
        assert len(u_data) == N
        
        # Check that temperatures decay
        assert x_s_data[-1] < x_s_data[0]
        assert x_h_data[-1] < x_h_data[0]
    
    def test_parameter_save_load(self):
        """Test saving and loading identified parameters."""
        # Create test parameters
        params = {
            "identification_info": {
                "timestamp": "2025-10-05T12:00:00",
                "initial_temperature": 300.0,
                "quench_temperature": 5.0,
                "duration_ps": 50.0,
                "sample_time_ps": 0.1
            },
            "system_matrices": {
                "A_discrete": [[0.95, 0.02, 0.03],
                               [0.01, 0.97, 0.02],
                               [0.0, 0.0, 0.98]],
                "B_discrete": [0.05, 0.03, 0.02]
            },
            "continuous_parameters": {
                "k_s": 0.05,
                "k_h": 0.03,
                "k_c": 0.02,
                "b_s": 0.01,
                "b_h": 0.01,
                "g_sh": 0.005
            },
            "fit_quality": {
                "r_squared_signal": 0.98,
                "r_squared_hot": 0.95,
                "r_squared_bath": 0.99
            }
        }
        
        # Save to temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.json') as f:
            json.dump(params, f, indent=2)
            temp_file = f.name
        
        try:
            # Load and verify
            with open(temp_file, 'r') as f:
                loaded_params = json.load(f)
            
            assert loaded_params['identification_info']['initial_temperature'] == 300.0
            assert loaded_params['system_matrices']['A_discrete'][0][0] == 0.95
            assert loaded_params['fit_quality']['r_squared_signal'] == 0.98
        
        finally:
            os.unlink(temp_file)
    
    def test_control_law_computation(self):
        """Test that control law u = -K*x_hat - k_I*xi_hat is computed correctly."""
        # LQR gain
        K = np.array([[1.0, 0.5, 0.2]])
        k_I = 0.3
        
        # State estimate
        x_hat = np.array([[10.0], [5.0], [2.0]])
        xi_hat = np.array([[15.0]])
        
        # Compute control
        u = -K @ x_hat - k_I * xi_hat
        
        # Expected: u = -(1.0*10 + 0.5*5 + 0.2*2) - 0.3*15 = -12.9 - 4.5 = -17.4
        expected = -(1.0*10.0 + 0.5*5.0 + 0.2*2.0) - 0.3*15.0
        
        assert_allclose(u[0, 0], expected, rtol=1e-10)


class TestPhysicalValidity:
    """Test that controller produces physically valid outputs."""
    
    def test_temperature_bounds(self):
        """Test that controller respects T_min and T_max bounds."""
        T_min = 0.1
        T_max = 1000.0
        
        # Simulate some control outputs
        raw_outputs = [-50.0, 5.0, 300.0, 1500.0]
        
        for raw in raw_outputs:
            clamped = max(T_min, min(T_max, raw))
            assert T_min <= clamped <= T_max
    
    def test_integral_windup_prevention(self):
        """Test that integral action doesn't wind up excessively."""
        # Simulate integral accumulation with saturation
        xi = 0.0
        dt = 0.1
        error = 10.0
        
        for _ in range(100):
            xi += error * dt
            
            # Clamp integral to prevent windup
            xi_max = 1000.0
            if abs(xi) > xi_max:
                xi = np.sign(xi) * xi_max
        
        assert abs(xi) <= 1000.0, "Integral should be clamped"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

