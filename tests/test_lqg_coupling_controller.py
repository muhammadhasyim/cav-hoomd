#!/usr/bin/env python3
"""
Comprehensive test suite for LQG coupling controller implementation.

This module provides unit tests and integration tests for all components
of the LQG coupling controller system including:
- System identification logic
- Controllability and observability analysis  
- Kalman filter state estimation
- LQR control law computation
- Anti-windup and saturation handling
- Lambda to epsilon conversion consistency
- Numerical stability checks
"""

import pytest
import numpy as np
import tempfile
import json
import os
from unittest.mock import Mock, MagicMock, patch
from pathlib import Path
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for testing
import matplotlib.pyplot as plt

# Import test target - this will fail initially (TDD approach)
try:
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))
    from cavitymd.controllers.lqg_coupling import LQGCouplingController
    HAS_LQG_COUPLING = True
except ImportError:
    HAS_LQG_COUPLING = False


class TestSystemIdentification:
    """Test system identification functionality."""
    
    def test_step_excitation_design(self):
        """Test step excitation signal generation for system ID."""
        if not HAS_LQG_COUPLING:
            pytest.skip("LQGCouplingController not yet implemented")
            
        # Mock controller parameters
        controller = LQGCouplingController(
            temperature_methods=['lj_coulombic', 'harmonic_equipartition', 'kinetic'],
            bath_temperature_method='kinetic',
            equilibrium_duration_ps=5.0,
            step_duration_ps=5.0,
            n_steps=2,
            coupling_step_size=1e-4,
            target_temperature=100.0,
            update_interval_ps=0.1
        )
        
        # Test step excitation generation
        times = np.linspace(0, 30, 301)  # 0.1 ps timestep
        lambda_values = controller._generate_step_excitation(times)
        
        # Check phases
        assert np.all(lambda_values[:50] == 0.0)  # Equilibrium phase (0-5 ps)
        assert np.all(lambda_values[50:100] == 5e-5)  # Step 1 (5-10 ps)
        assert np.all(lambda_values[100:150] == 1e-4)  # Step 2 (10-15 ps) 
        assert np.all(lambda_values[150:200] == 5e-5)  # Step down 1 (15-20 ps)
        assert np.all(lambda_values[200:250] == 0.0)  # Step down 2 (20-25 ps)
        assert np.all(lambda_values[250:] == 0.0)  # Final equilibrium (25-30 ps)
    
    def test_system_identification_with_synthetic_data(self):
        """Test system ID with known synthetic data."""
        if not HAS_LQG_COUPLING:
            pytest.skip("LQGCouplingController not yet implemented")
            
        # Generate synthetic system response
        dt = 0.1
        n_samples = 300
        
        # True system matrices (4-state: T_s, T_harmonic, T_kinetic, T_bath)
        A_true = np.array([
            [0.95, 0.02, 0.01, 0.02],  # T_s dynamics
            [0.01, 0.97, 0.01, 0.01],  # T_harmonic dynamics
            [0.02, 0.01, 0.96, 0.01],  # T_kinetic dynamics
            [0.00, 0.00, 0.02, 0.98]   # T_bath dynamics
        ])
        
        B_true = np.array([
            [0.05],  # T_s responds to coupling rate
            [-0.03], # T_harmonic (negative - coupling heats vibrations)
            [0.02],  # T_kinetic
            [0.01]   # T_bath (weak response)
        ])
        
        C_true = np.eye(4)  # All states observable
        
        # Generate synthetic data
        x = np.zeros((4, n_samples + 1))
        u = np.random.randn(1, n_samples) * 0.1  # Random coupling rate inputs
        y = np.zeros((4, n_samples))
        
        for k in range(n_samples):
            x[:, k+1:k+2] = A_true @ x[:, k:k+1] + B_true @ u[:, k:k+1]
            y[:, k:k+1] = C_true @ x[:, k:k+1]
        
        # Test system identification
        controller = LQGCouplingController(
            temperature_methods=['lj_coulombic', 'harmonic_equipartition', 'kinetic'],
            bath_temperature_method='kinetic',
            target_temperature=100.0,
            update_interval_ps=0.1
        )
        
        A_id, B_id, C_id = controller._identify_system_matrices(
            y[:, :-1], y[:, 1:], u[:, :-1], dt
        )
        
        # Check identification accuracy
        np.testing.assert_allclose(A_id, A_true, rtol=1e-2)
        np.testing.assert_allclose(B_id, B_true, rtol=1e-2)
        np.testing.assert_allclose(C_id, C_true, rtol=1e-3)


class TestControllabilityObservability:
    """Test controllability and observability analysis."""
    
    def test_controllability_matrix_computation(self):
        """Test controllability matrix computation and rank checking."""
        if not HAS_LQG_COUPLING:
            pytest.skip("LQGCouplingController not yet implemented")
            
        # Test system matrices
        A = np.array([
            [0.95, 0.02, 0.01, 0.02],
            [0.01, 0.97, 0.01, 0.01],
            [0.02, 0.01, 0.96, 0.01],
            [0.00, 0.00, 0.02, 0.98]
        ])
        
        B = np.array([
            [0.05],
            [-0.03],
            [0.02],
            [0.01]
        ])
        
        controller = LQGCouplingController(
            temperature_methods=['lj_coulombic', 'harmonic_equipartition', 'kinetic'],
            bath_temperature_method='kinetic',
            target_temperature=100.0,
            update_interval_ps=0.1
        )
        
        # Test controllability
        is_controllable, ctrb_matrix = controller._check_controllability(A, B)
        
        # Should be controllable with this system
        assert is_controllable
        assert ctrb_matrix.shape == (4, 4)
        assert np.linalg.matrix_rank(ctrb_matrix) == 4
    
    def test_observability_matrix_computation(self):
        """Test observability matrix computation and rank checking.""" 
        if not HAS_LQG_COUPLING:
            pytest.skip("LQGCouplingController not yet implemented")
            
        # Test system matrices
        A = np.array([
            [0.95, 0.02, 0.01, 0.02],
            [0.01, 0.97, 0.01, 0.01],
            [0.02, 0.01, 0.96, 0.01],
            [0.00, 0.00, 0.02, 0.98]
        ])
        
        C = np.eye(4)  # All states observable
        
        controller = LQGCouplingController(
            temperature_methods=['lj_coulombic', 'harmonic_equipartition', 'kinetic'],
            bath_temperature_method='kinetic',
            target_temperature=100.0,
            update_interval_ps=0.1
        )
        
        # Test observability
        is_observable, obsv_matrix = controller._check_observability(A, C)
        
        # Should be observable with this system
        assert is_observable
        assert obsv_matrix.shape == (16, 4)  # 4*4 x 4
        assert np.linalg.matrix_rank(obsv_matrix) == 4


class TestKalmanFilter:
    """Test Kalman filter implementation."""
    
    def test_kalman_filter_initialization(self):
        """Test proper Kalman filter initialization."""
        if not HAS_LQG_COUPLING:
            pytest.skip("LQGCouplingController not yet implemented")
            
        controller = LQGCouplingController(
            temperature_methods=['lj_coulombic', 'harmonic_equipartition', 'kinetic'],
            bath_temperature_method='kinetic',
            target_temperature=100.0,
            update_interval_ps=0.1,
            process_noise_std=0.1,
            measurement_noise_std=0.5
        )
        
        # Test initialization
        A = np.eye(4)
        C = np.eye(4)
        controller._initialize_kalman_filter(A, C)
        
        # Check state and covariance initialization
        assert controller._kalman_state.shape == (4, 1)
        assert controller._kalman_covariance.shape == (4, 4)
        assert np.all(controller._kalman_state == 0)
        assert np.all(np.diag(controller._kalman_covariance) > 0)
    
    def test_kalman_filter_update(self):
        """Test Kalman filter prediction and update steps."""
        if not HAS_LQG_COUPLING:
            pytest.skip("LQGCouplingController not yet implemented")
            
        # System matrices
        A = np.array([
            [0.95, 0.02],
            [0.01, 0.97]
        ])
        B = np.array([[0.05], [0.03]])
        C = np.eye(2)
        
        controller = LQGCouplingController(
            temperature_methods=['lj_coulombic', 'harmonic_equipartition'],
            target_temperature=100.0,
            update_interval_ps=0.1,
            process_noise_std=0.1,
            measurement_noise_std=0.5
        )
        
        # Initialize with 2-state system for testing
        controller._n_states = 2
        controller._initialize_kalman_filter(A, C)
        
        # Test prediction step
        u = np.array([[0.1]])  # Control input
        controller._kalman_predict(A, B, u)
        
        # Test update step
        y = np.array([[100.5], [200.3]])  # Measurements
        controller._kalman_update(C, y)
        
        # State should be updated
        assert controller._kalman_state.shape == (2, 1)
        assert controller._kalman_covariance.shape == (2, 2)


class TestLQRController:
    """Test LQR control law computation."""
    
    def test_lqr_gain_computation(self):
        """Test LQR optimal gain computation."""
        if not HAS_LQG_COUPLING:
            pytest.skip("LQGCouplingController not yet implemented")
            
        # Test system
        A = np.array([
            [0.95, 0.02],
            [0.01, 0.97]
        ])
        B = np.array([[0.05], [0.03]])
        
        # LQR weights
        Q = np.diag([100.0, 1.0])  # Regulate first state tightly
        R = np.array([[1.0]])      # Control effort penalty
        
        controller = LQGCouplingController(
            temperature_methods=['lj_coulombic', 'harmonic_equipartition'],
            target_temperature=100.0,
            update_interval_ps=0.1
        )
        
        # Compute LQR gains
        K = controller._compute_lqr_gains(A, B, Q, R)
        
        # Check gain dimensions and properties
        assert K.shape == (1, 2)  # 1 control input, 2 states
        assert np.all(np.isfinite(K))  # No NaN/inf values
    
    def test_control_law_computation(self):
        """Test control law with anti-windup."""
        if not HAS_LQG_COUPLING:
            pytest.skip("LQGCouplingController not yet implemented")
            
        controller = LQGCouplingController(
            temperature_methods=['lj_coulombic', 'harmonic_equipartition'],
            target_temperature=100.0,
            update_interval_ps=0.1,
            lambda_min=0.0,
            lambda_max=1e-2
        )
        
        # Mock LQR gains and state
        controller._lqr_gains = np.array([[0.5, 0.3]])
        controller._integral_state = 0.1
        controller._integral_gain = 2.0
        
        # Test control computation
        state_error = np.array([[1.0], [0.5]])  # Temperature errors
        control = controller._compute_control_law(state_error)
        
        # Check control output
        assert isinstance(control, float)
        assert controller.lambda_min <= control <= controller.lambda_max


class TestAntiWindup:
    """Test anti-windup and saturation handling."""
    
    def test_saturation_handling(self):
        """Test smooth saturation function."""
        if not HAS_LQG_COUPLING:
            pytest.skip("LQGCouplingController not yet implemented")
            
        controller = LQGCouplingController(
            temperature_methods=['lj_coulombic'],
            target_temperature=100.0,
            update_interval_ps=0.1,
            lambda_min=0.0,
            lambda_max=1e-2
        )
        
        # Test saturation function
        test_values = [-0.1, -1e-3, 0.0, 5e-3, 1e-2, 2e-2]
        for val in test_values:
            saturated = controller._apply_saturation(val)
            assert controller.lambda_min <= saturated <= controller.lambda_max
    
    def test_anti_windup_logic(self):
        """Test anti-windup integrator clamping."""
        if not HAS_LQG_COUPLING:
            pytest.skip("LQGCouplingController not yet implemented")
            
        controller = LQGCouplingController(
            temperature_methods=['lj_coulombic'],
            target_temperature=100.0,
            update_interval_ps=0.1,
            lambda_min=0.0,
            lambda_max=1e-2
        )
        
        # Initialize integrator
        controller._integral_state = 0.0
        
        # Test normal operation (no saturation)
        error = 1.0
        dt = 0.1
        controller._update_integral_state(error, dt, saturated=False)
        assert controller._integral_state == error * dt
        
        # Test anti-windup (saturated)
        initial_integral = controller._integral_state
        controller._update_integral_state(error, dt, saturated=True)
        # Integral should not increase when saturated
        assert controller._integral_state == initial_integral


class TestNumericalStability:
    """Test numerical stability and conditioning."""
    
    def test_matrix_conditioning_check(self):
        """Test matrix conditioning checks."""
        if not HAS_LQG_COUPLING:
            pytest.skip("LQGCouplingController not yet implemented")
            
        controller = LQGCouplingController(
            temperature_methods=['lj_coulombic'],
            target_temperature=100.0,
            update_interval_ps=0.1
        )
        
        # Test well-conditioned matrix
        A_good = np.eye(3) * 0.95 + np.ones((3, 3)) * 0.01
        is_stable = controller._check_matrix_stability(A_good)
        assert is_stable
        
        # Test ill-conditioned matrix
        A_bad = np.array([[1.0, 1e-12], [1e12, 1.0]])
        is_stable = controller._check_matrix_stability(A_bad)
        assert not is_stable
    
    def test_joseph_form_covariance_update(self):
        """Test Joseph form for numerical stability."""
        if not HAS_LQG_COUPLING:
            pytest.skip("LQGCouplingController not yet implemented")
            
        controller = LQGCouplingController(
            temperature_methods=['lj_coulombic', 'harmonic_equipartition'],
            target_temperature=100.0,
            update_interval_ps=0.1
        )
        
        # Test covariance update
        P = np.eye(2) * 2.0  # Prior covariance
        K = np.array([[0.5], [0.3]])  # Kalman gain
        C = np.eye(2)  # Measurement matrix
        R = np.eye(2) * 0.1  # Measurement noise
        
        P_updated = controller._joseph_form_update(P, K, C, R)
        
        # Check positive definiteness
        assert P_updated.shape == (2, 2)
        eigenvals = np.linalg.eigvals(P_updated)
        assert np.all(eigenvals > 0)  # Positive definite


class TestLambdaToEpsilonConversion:
    """Test lambda to epsilon conversion consistency."""
    
    def test_lambda_integration(self):
        """Test lambda integration from coupling rate."""
        if not HAS_LQG_COUPLING:
            pytest.skip("LQGCouplingController not yet implemented")
            
        controller = LQGCouplingController(
            temperature_methods=['lj_coulombic'],
            target_temperature=100.0,
            update_interval_ps=0.1
        )
        
        # Initialize lambda value
        controller._current_lambda = 0.0
        
        # Test integration
        dlambda_dt = 0.001  # Coupling rate
        dt = 0.1
        
        new_lambda = controller._integrate_lambda(dlambda_dt, dt)
        expected_lambda = controller._current_lambda + dlambda_dt * dt
        
        assert abs(new_lambda - expected_lambda) < 1e-12
    
    def test_epsilon_scaling(self):
        """Test epsilon = lambda * omega_c scaling."""
        if not HAS_LQG_COUPLING:
            pytest.skip("LQGCouplingController not yet implemented")
            
        controller = LQGCouplingController(
            temperature_methods=['lj_coulombic'],
            target_temperature=100.0,
            update_interval_ps=0.1
        )
        
        # Test scaling
        lambda_val = 1e-4
        omega_c = 1560.0  # cm^-1 in atomic units
        
        epsilon = controller._lambda_to_epsilon(lambda_val, omega_c)
        expected_epsilon = lambda_val * omega_c
        
        assert abs(epsilon - expected_epsilon) < 1e-12


class TestSystemIDOutput:
    """Test system identification analysis and output."""
    
    def test_pdf_generation(self):
        """Test PDF figure generation for system ID results."""
        if not HAS_LQG_COUPLING:
            pytest.skip("LQGCouplingController not yet implemented")
            
        controller = LQGCouplingController(
            temperature_methods=['lj_coulombic', 'harmonic_equipartition'],
            target_temperature=100.0,
            update_interval_ps=0.1
        )
        
        # Mock system ID data
        times = np.linspace(0, 10, 100)
        temperatures = {
            'lj_coulombic': 100 + np.random.randn(100) * 0.5,
            'harmonic_equipartition': 200 + np.random.randn(100) * 1.0
        }
        lambda_excitation = np.where(times > 5, 1e-4, 0.0)
        
        with tempfile.NamedTemporaryFile(suffix='.pdf', delete=False) as f:
            pdf_path = f.name
        
        try:
            # Test PDF generation
            controller._generate_system_id_report(
                times, temperatures, lambda_excitation, 
                r2_scores={'raw': 0.95, 'smoothed': 0.98},
                output_path=pdf_path
            )
            
            # Check file was created
            assert os.path.exists(pdf_path)
            assert os.path.getsize(pdf_path) > 0
            
        finally:
            if os.path.exists(pdf_path):
                os.unlink(pdf_path)
    
    def test_r2_calculation(self):
        """Test R-squared goodness-of-fit calculation."""
        if not HAS_LQG_COUPLING:
            pytest.skip("LQGCouplingController not yet implemented")
            
        controller = LQGCouplingController(
            temperature_methods=['lj_coulombic'],
            target_temperature=100.0,
            update_interval_ps=0.1
        )
        
        # Perfect fit test
        y_true = np.array([1, 2, 3, 4, 5])
        y_pred = np.array([1, 2, 3, 4, 5])
        r2 = controller._calculate_r2(y_true, y_pred)
        assert abs(r2 - 1.0) < 1e-10
        
        # No fit test (predict mean)
        y_pred_mean = np.full_like(y_true, np.mean(y_true))
        r2 = controller._calculate_r2(y_true, y_pred_mean)
        assert abs(r2 - 0.0) < 1e-10


class TestPerformanceMonitoring:
    """Test performance monitoring and metrics."""
    
    def test_performance_metrics_calculation(self):
        """Test control performance indicators."""
        if not HAS_LQG_COUPLING:
            pytest.skip("LQGCouplingController not yet implemented")
            
        controller = LQGCouplingController(
            temperature_methods=['lj_coulombic'],
            target_temperature=100.0,
            update_interval_ps=0.1
        )
        
        # Mock temperature data showing step response
        times = np.linspace(0, 10, 100)
        target = 100.0
        # Simulated step response with overshoot and settling
        temperatures = target + 5 * np.exp(-times) * np.cos(2*times)
        
        # Calculate performance metrics
        metrics = controller._calculate_performance_metrics(times, temperatures, target)
        
        # Check metric types
        assert 'settling_time' in metrics
        assert 'overshoot_percent' in metrics
        assert 'steady_state_error' in metrics
        assert all(isinstance(v, (int, float)) for v in metrics.values())
    
    def test_kalman_filter_health_monitoring(self):
        """Test Kalman filter diagnostics."""
        if not HAS_LQG_COUPLING:
            pytest.skip("LQGCouplingController not yet implemented")
            
        controller = LQGCouplingController(
            temperature_methods=['lj_coulombic', 'harmonic_equipartition'],
            target_temperature=100.0,
            update_interval_ps=0.1
        )
        
        # Mock Kalman filter state
        controller._kalman_state = np.array([[100.0], [200.0]])
        controller._kalman_covariance = np.eye(2) * 0.5
        controller._innovation_history = [0.1, 0.2, 0.15, 0.08]
        
        # Test health metrics
        health = controller._assess_kalman_health()
        
        assert 'state_magnitude' in health
        assert 'covariance_trace' in health  
        assert 'innovation_variance' in health
        assert all(isinstance(v, (int, float)) for v in health.values())


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
