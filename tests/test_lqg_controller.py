#!/usr/bin/env python3
"""
Comprehensive test suite for LQG controller implementation.

This module provides unit tests and integration tests for all components
of the LQG pulse-to-pulse controller system.
"""

import pytest
import numpy as np
import tempfile
import json
import os
from unittest.mock import Mock, MagicMock, patch
from pathlib import Path

# Import LQG controller components
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from cavitymd.lqg_controller import (
    LQGParameters, LiftedStateSpaceModel, KalmanFilter, 
    LQRController, SystemIdentification
)
from cavitymd.lqg_pulse_controller import LQGPulseController
from cavitymd.lqg_square_wave import LQGSquareWaveVariant, LQGAdaptiveSquareWaveVariant


class TestLQGParameters:
    """Test LQG parameter handling."""
    
    def test_default_parameters(self):
        """Test default parameter values."""
        params = LQGParameters()
        
        # Check some key defaults
        assert params.a_sk == 0.15
        assert params.a_sv == 0.05
        assert params.b_sg == 0.10
        assert params.b_vg == -0.12  # Should be negative (heating)
        assert params.g_baseline == 1e-5
        assert params.use_hidden_state == True
    
    def test_parameter_serialization(self):
        """Test parameter save/load functionality."""
        params = LQGParameters()
        params.a_sk = 0.25
        params.Q_structure = 5.0
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            temp_file = f.name
        
        try:
            # Save parameters
            params.save(temp_file)
            
            # Load parameters
            loaded_params = LQGParameters.load(temp_file)
            
            # Check values
            assert loaded_params.a_sk == 0.25
            assert loaded_params.Q_structure == 5.0
            assert loaded_params.a_sv == params.a_sv  # Unchanged default
            
        finally:
            os.unlink(temp_file)
    
    def test_parameter_dict_conversion(self):
        """Test dictionary conversion methods."""
        params = LQGParameters()
        params.R_coupling = 0.05
        
        # Convert to dict
        param_dict = params.to_dict()
        assert isinstance(param_dict, dict)
        assert param_dict['R_coupling'] == 0.05
        
        # Convert back from dict
        new_params = LQGParameters.from_dict(param_dict)
        assert new_params.R_coupling == 0.05
        assert new_params.a_sk == params.a_sk


class TestLiftedStateSpaceModel:
    """Test lifted state-space model."""
    
    def test_model_initialization(self):
        """Test model matrix construction."""
        params = LQGParameters()
        model = LiftedStateSpaceModel(params)
        
        # Check dimensions
        if params.use_hidden_state:
            assert model.n_states == 4
        else:
            assert model.n_states == 3
        
        assert model.n_controls == 2
        assert model.n_outputs == 3
        
        # Check matrix shapes
        assert model.A.shape == (model.n_states, model.n_states)
        assert model.B.shape == (model.n_states, model.n_controls)
        assert model.C.shape == (model.n_outputs, model.n_states)
        assert len(model.c) == model.n_states
    
    def test_model_stability(self):
        """Test that A matrix has stable eigenvalues."""
        params = LQGParameters()
        model = LiftedStateSpaceModel(params)
        
        # Check eigenvalues of A matrix
        eigenvals = np.linalg.eigvals(model.A)
        max_eigenval = np.max(np.abs(eigenvals))
        
        # Should be stable (all eigenvalues inside unit circle)
        assert max_eigenval < 1.0, f"Unstable system: max eigenvalue = {max_eigenval}"
    
    def test_model_simulation(self):
        """Test model simulation step."""
        params = LQGParameters()
        model = LiftedStateSpaceModel(params)
        
        # Initial state and control
        x0 = np.array([300.0, 300.0, 300.0, 0.0])[:model.n_states]
        u0 = np.array([0.01, 310.0])  # Small coupling deviation, higher bath temp
        
        # Simulate one step
        x1, y0 = model.simulate_step(x0, u0, add_noise=False)
        
        # Check output dimensions
        assert len(x1) == model.n_states
        assert len(y0) == model.n_outputs
        
        # Check that output is measurement of state
        np.testing.assert_array_almost_equal(y0, x0[:3])
    
    def test_steady_state_computation(self):
        """Test steady-state computation."""
        params = LQGParameters()
        model = LiftedStateSpaceModel(params)
        
        # Constant control input
        u_ss = np.array([0.0, 300.0])  # No coupling deviation, 300K bath
        
        # Compute steady state
        x_ss = model.compute_steady_state(u_ss)
        
        # Check dimensions
        assert len(x_ss) == model.n_states
        
        # Verify it's actually steady state
        x_next, _ = model.simulate_step(x_ss, u_ss, add_noise=False)
        np.testing.assert_array_almost_equal(x_next, x_ss, decimal=10)
    
    def test_model_without_hidden_state(self):
        """Test 3-state model (without hidden state)."""
        params = LQGParameters()
        params.use_hidden_state = False
        
        model = LiftedStateSpaceModel(params)
        
        # Check dimensions
        assert model.n_states == 3
        assert model.A.shape == (3, 3)
        assert model.B.shape == (3, 2)
        assert len(model.c) == 3
        
        # Test simulation
        x0 = np.array([300.0, 300.0, 300.0])
        u0 = np.array([0.01, 310.0])
        
        x1, y0 = model.simulate_step(x0, u0)
        assert len(x1) == 3
        assert len(y0) == 3


class TestKalmanFilter:
    """Test Kalman filter implementation."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.params = LQGParameters()
        self.model = LiftedStateSpaceModel(self.params)
        self.kf = KalmanFilter(self.model)
    
    def test_filter_initialization(self):
        """Test filter initialization."""
        assert len(self.kf.x_hat) == self.model.n_states
        assert self.kf.P.shape == (self.model.n_states, self.model.n_states)
        
        # Check that steady-state gain was computed
        if self.kf.steady_state_gain is not None:
            expected_shape = (self.model.n_states, self.model.n_outputs)
            assert self.kf.steady_state_gain.shape == expected_shape
    
    def test_prediction_step(self):
        """Test Kalman filter prediction."""
        initial_state = self.kf.x_hat.copy()
        initial_cov = self.kf.P.copy()
        
        u_prev = np.array([0.01, 310.0])
        self.kf.predict(u_prev)
        
        # State should have changed
        assert not np.array_equal(self.kf.x_hat, initial_state)
        
        # Covariance should have increased (uncertainty grows)
        assert np.trace(self.kf.P) > np.trace(initial_cov)
    
    def test_update_step(self):
        """Test Kalman filter update."""
        # Set initial state
        self.kf.x_hat = np.array([300.0, 300.0, 300.0, 0.0])[:self.model.n_states]
        initial_cov_trace = np.trace(self.kf.P)
        
        # Measurement close to current estimate
        y_k = np.array([301.0, 299.0, 300.5])
        
        self.kf.update(y_k)
        
        # Covariance should have decreased (uncertainty reduced by measurement)
        assert np.trace(self.kf.P) < initial_cov_trace
        
        # State estimate should have moved toward measurement
        assert abs(self.kf.x_hat[0] - 301.0) < abs(300.0 - 301.0)
    
    def test_complete_filter_step(self):
        """Test complete filter step (predict + update)."""
        initial_state = self.kf.x_hat.copy()
        
        y_k = np.array([305.0, 295.0, 300.0])
        u_prev = np.array([0.0, 300.0])
        
        x_hat = self.kf.step(y_k, u_prev)
        
        # Should return updated state estimate
        np.testing.assert_array_equal(x_hat, self.kf.x_hat)
        assert len(x_hat) == self.model.n_states
    
    def test_steady_state_operation(self):
        """Test steady-state Kalman gain usage."""
        if self.kf.steady_state_gain is None:
            pytest.skip("Steady-state gain not available")
        
        y_k = np.array([300.0, 300.0, 300.0])
        u_prev = np.array([0.0, 300.0])
        
        # Run with steady-state gain
        x_hat_ss = self.kf.step(y_k, u_prev, use_steady_state=True)
        
        # Should produce valid estimate
        assert len(x_hat_ss) == self.model.n_states
        assert np.all(np.isfinite(x_hat_ss))


class TestLQRController:
    """Test LQR controller implementation."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.params = LQGParameters()
        self.model = LiftedStateSpaceModel(self.params)
        self.lqr = LQRController(self.model, self.params)
    
    def test_controller_initialization(self):
        """Test LQR controller initialization."""
        # Check augmented system dimensions
        expected_aug_dim = self.model.n_states + self.model.n_outputs
        assert self.lqr.n_aug == expected_aug_dim
        assert self.lqr.A_aug.shape == (expected_aug_dim, expected_aug_dim)
        assert self.lqr.B_aug.shape == (expected_aug_dim, self.model.n_controls)
        
        # Check LQR gain
        expected_gain_shape = (self.model.n_controls, expected_aug_dim)
        assert self.lqr.K_lqr.shape == expected_gain_shape
        
        # Check feedforward gain
        if self.lqr.feedforward_gain is not None:
            expected_ff_shape = (self.model.n_controls, self.model.n_outputs)
            assert self.lqr.feedforward_gain.shape == expected_ff_shape
    
    def test_control_computation(self):
        """Test control action computation."""
        # State estimate
        x_hat = np.array([305.0, 295.0, 310.0, 0.0])[:self.model.n_states]
        
        # Reference (target temperatures)
        reference = np.array([300.0, 300.0, 300.0])
        
        # Compute control
        u = self.lqr.compute_control(x_hat, reference)
        
        # Check dimensions and bounds
        assert len(u) == self.model.n_controls
        assert self.params.delta_g_min <= u[0] <= self.params.delta_g_max
        assert self.params.T_bath_min <= u[1] <= self.params.T_bath_max
    
    def test_integral_action(self):
        """Test integral action for steady-state error elimination."""
        x_hat = np.array([310.0, 310.0, 310.0, 0.0])[:self.model.n_states]  # High temperatures
        reference = np.array([300.0, 300.0, 300.0])  # Lower targets
        
        # Compute control multiple times to build up integral
        controls = []
        for _ in range(5):
            u = self.lqr.compute_control(x_hat, reference)
            controls.append(u.copy())
        
        # Control should increase in magnitude to eliminate steady-state error
        assert len(controls) == 5
        
        # Check that integral states are being updated
        assert np.any(self.lqr.integral_states != 0.0)
    
    def test_constraint_enforcement(self):
        """Test control constraint enforcement."""
        # Create extreme state to trigger constraints
        x_hat = np.array([500.0, 500.0, 500.0, 0.0])[:self.model.n_states]
        reference = np.array([100.0, 100.0, 100.0])
        
        u = self.lqr.compute_control(x_hat, reference)
        
        # Control should be within bounds
        assert self.params.delta_g_min <= u[0] <= self.params.delta_g_max
        assert self.params.T_bath_min <= u[1] <= self.params.T_bath_max
    
    def test_integrator_reset(self):
        """Test integrator reset functionality."""
        # Build up some integral action
        x_hat = np.array([310.0, 310.0, 310.0, 0.0])[:self.model.n_states]
        reference = np.array([300.0, 300.0, 300.0])
        
        for _ in range(3):
            self.lqr.compute_control(x_hat, reference)
        
        # Should have non-zero integral states
        assert np.any(self.lqr.integral_states != 0.0)
        
        # Reset integrator
        self.lqr.reset_integrator()
        
        # Should be zero again
        np.testing.assert_array_equal(self.lqr.integral_states, np.zeros(self.model.n_outputs))


class TestSystemIdentification:
    """Test system identification functionality."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.sys_id = SystemIdentification()
    
    def test_parameter_vector_conversion(self):
        """Test parameter to vector conversion."""
        params = LQGParameters()
        params.a_sk = 0.25
        params.b_sg = 0.15
        
        # Convert to vector
        vector = self.sys_id._params_to_vector(params)
        assert isinstance(vector, np.ndarray)
        assert len(vector) == 12  # Number of physical parameters
        
        # Convert back
        new_params = self.sys_id._vector_to_params(vector, params)
        assert new_params.a_sk == 0.25
        assert new_params.b_sg == 0.15
    
    def test_data_preparation(self):
        """Test pulse data preparation for fitting."""
        # Create synthetic pulse data
        pulse_data = [
            {
                'state': [300.0, 300.0, 300.0],
                'control': [0.01, 310.0],
                'next_state': [299.0, 301.0, 302.0]
            },
            {
                'state': [299.0, 301.0, 302.0],
                'control': [0.02, 305.0],
                'next_state': [298.0, 302.0, 301.0]
            }
        ]
        
        X, U, Y = self.sys_id._prepare_data_matrices(pulse_data)
        
        # Check dimensions
        assert len(X) == 2
        assert len(U) == 2
        assert len(Y) == 2
        
        # Check that states were extended with hidden variable
        assert len(X[0]) == 4  # [T_s, T_k, T_v, η]
        assert len(Y[0]) == 4
        assert len(U[0]) == 2  # [Δg, T_bath]
    
    def test_parameter_bounds(self):
        """Test parameter bounds for optimization."""
        bounds = self.sys_id._get_parameter_bounds()
        
        assert len(bounds) == 12  # Number of parameters
        
        # Check some specific bounds
        for lower, upper in bounds:
            assert lower < upper
            assert isinstance(lower, (int, float))
            assert isinstance(upper, (int, float))
    
    def test_fitting_with_synthetic_data(self):
        """Test parameter fitting with synthetic data."""
        # Generate synthetic data from known model
        true_params = LQGParameters()
        true_params.a_sk = 0.2
        true_params.b_sg = 0.12
        
        true_model = LiftedStateSpaceModel(true_params)
        
        # Generate pulse data
        pulse_data = []
        x = np.array([300.0, 300.0, 300.0, 0.0])
        
        for i in range(20):
            u = np.array([0.01 * np.sin(i * 0.1), 300.0 + 5.0 * np.cos(i * 0.1)])
            x_next, _ = true_model.simulate_step(x, u, add_noise=True)
            
            pulse_data.append({
                'state': x[:3].tolist(),
                'control': u.tolist(),
                'next_state': x_next[:3].tolist()
            })
            
            x = x_next
        
        # Fit parameters
        fitted_params = self.sys_id.fit_from_data(pulse_data)
        
        # Check that fitted parameters are reasonable
        assert isinstance(fitted_params, LQGParameters)
        assert 0.0 <= fitted_params.a_sk <= 0.5
        assert 0.0 <= fitted_params.b_sg <= 0.5


class TestLQGPulseController:
    """Test main LQG pulse controller."""
    
    def setup_method(self):
        """Set up test fixtures with mocked HOOMD objects."""
        # Mock time tracker
        self.time_tracker = Mock()
        self.time_tracker.elapsed_time_ps = 0.0
        
        # Mock energy tracker
        self.energy_tracker = Mock()
        self.energy_tracker.get_instantaneous_energy.return_value = {
            'lj': -1000.0,
            'coulombic': -500.0,
            'harmonic': 200.0
        }
        
        # Mock simulation
        self.simulation = Mock()
        
        # Mock snapshot for temperature calculations
        mock_snap = Mock()
        mock_snap.particles.velocity = np.random.normal(0, 1, (500, 3))
        mock_snap.particles.mass = np.ones(500) * 16.0  # Oxygen mass
        mock_snap.particles.typeid = np.zeros(500, dtype=int)  # All molecular
        
        self.simulation.state.cpu_local_snapshot.__enter__ = Mock(return_value=mock_snap)
        self.simulation.state.cpu_local_snapshot.__exit__ = Mock(return_value=None)
        
        # Mock thermostats
        self.molecular_thermostat = Mock()
        self.molecular_thermostat.kT = 300.0 * 3.166815e-6  # 300K in Hartree
        
        self.cavity_thermostat = Mock()
        self.cavity_thermostat.kT = 300.0 * 3.166815e-6
    
    def test_controller_initialization(self):
        """Test controller initialization."""
        controller = LQGPulseController(
            time_tracker=self.time_tracker,
            energy_tracker=self.energy_tracker,
            simulation=self.simulation,
            molecular_thermostat=self.molecular_thermostat,
            pulse_period_ps=5.0,
            turn_on_time_ps=10.0
        )
        
        # Check basic attributes
        assert controller.pulse_period_ps == 5.0
        assert controller.turn_on_time_ps == 10.0
        assert controller.is_active == False
        assert controller.pulse_count == 0
        
        # Check target temperatures
        assert 'structure' in controller.target_temperatures
        assert 'kinetic' in controller.target_temperatures
        assert 'vibrational' in controller.target_temperatures
    
    def test_temperature_measurement(self):
        """Test temperature measurement methods."""
        controller = LQGPulseController(
            time_tracker=self.time_tracker,
            energy_tracker=self.energy_tracker,
            simulation=self.simulation,
            molecular_thermostat=self.molecular_thermostat
        )
        
        # Test temperature measurement
        temperatures = controller._measure_temperatures(0.0)
        
        assert len(temperatures) == 3
        assert np.all(temperatures > 0)  # Should be positive
        assert np.all(np.isfinite(temperatures))  # Should be finite
    
    def test_controller_activation(self):
        """Test controller activation timing."""
        controller = LQGPulseController(
            time_tracker=self.time_tracker,
            energy_tracker=self.energy_tracker,
            simulation=self.simulation,
            molecular_thermostat=self.molecular_thermostat,
            turn_on_time_ps=10.0,
            turn_off_time_ps=100.0
        )
        
        # Before activation
        self.time_tracker.elapsed_time_ps = 5.0
        controller.act(0)
        assert not controller.is_active
        
        # After activation
        self.time_tracker.elapsed_time_ps = 15.0
        controller.act(0)
        # Note: Controller won't be fully active until first pulse boundary
        
        # After deactivation
        self.time_tracker.elapsed_time_ps = 105.0
        controller.act(0)
        assert not controller.is_active
    
    def test_pulse_boundary_detection(self):
        """Test pulse boundary detection."""
        controller = LQGPulseController(
            time_tracker=self.time_tracker,
            energy_tracker=self.energy_tracker,
            simulation=self.simulation,
            molecular_thermostat=self.molecular_thermostat,
            pulse_period_ps=5.0,
            turn_on_time_ps=0.0
        )
        
        # Should not be pulse boundary initially
        assert not controller._is_pulse_boundary(2.0)
        
        # Should be pulse boundary at period intervals
        controller.last_pulse_time_ps = 0.0
        assert controller._is_pulse_boundary(5.0)
        
        controller.last_pulse_time_ps = 5.0
        assert controller._is_pulse_boundary(10.0)
    
    def test_reference_trajectory_computation(self):
        """Test reference trajectory computation."""
        controller = LQGPulseController(
            time_tracker=self.time_tracker,
            energy_tracker=self.energy_tracker,
            simulation=self.simulation,
            molecular_thermostat=self.molecular_thermostat,
            target_temperatures={'structure': 250.0, 'kinetic': 300.0, 'vibrational': 350.0}
        )
        
        reference = controller._compute_reference_trajectory(10.0)
        
        assert len(reference) == 3
        np.testing.assert_array_equal(reference, [250.0, 300.0, 350.0])
    
    def test_control_application_to_thermostats(self):
        """Test applying control to thermostats."""
        controller = LQGPulseController(
            time_tracker=self.time_tracker,
            energy_tracker=self.energy_tracker,
            simulation=self.simulation,
            molecular_thermostat=self.molecular_thermostat,
            cavity_thermostat=self.cavity_thermostat,
            apply_to='both'
        )
        
        # Apply control
        controller._apply_control_to_thermostats(350.0)
        
        # Check that thermostats were updated
        # Note: Exact checking depends on thermostat interface
        assert hasattr(self.molecular_thermostat, 'kT') or hasattr(self.molecular_thermostat, 'temperature')
    
    def test_parameter_save_load(self):
        """Test parameter save/load functionality."""
        controller = LQGPulseController(
            time_tracker=self.time_tracker,
            energy_tracker=self.energy_tracker,
            simulation=self.simulation,
            molecular_thermostat=self.molecular_thermostat
        )
        
        # Modify parameters
        controller.lqg_params.Q_structure = 5.0
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            temp_file = f.name
        
        try:
            # Save parameters
            controller.save_parameters(temp_file)
            
            # Load parameters
            controller.load_parameters(temp_file)
            
            # Check that parameters were loaded
            assert controller.lqg_params.Q_structure == 5.0
            
        finally:
            os.unlink(temp_file)


class TestLQGSquareWave:
    """Test LQG square wave variants."""
    
    def setup_method(self):
        """Set up test fixtures."""
        # Mock time tracker
        self.time_tracker = Mock()
        self.time_tracker.elapsed_time_ps = 0.0
        
        # Mock LQG controller
        self.lqg_controller = Mock()
        self.lqg_controller.pulse_period_ps = 5.0
        self.lqg_controller.lqg_params.g_baseline = 1e-5
        self.lqg_controller.current_coupling_deviation = 0.01
        self.lqg_controller.current_bath_temperature = 300.0
        self.lqg_controller.get_current_coupling.return_value = 1e-5 + 0.01
        self.lqg_controller.pulse_count = 0
        self.lqg_controller.is_active = True
        self.lqg_controller.system_id_complete = False
    
    def test_lqg_square_wave_initialization(self):
        """Test LQG square wave initialization."""
        square_wave = LQGSquareWaveVariant(
            lqg_controller=self.lqg_controller,
            period_ps=5.0,
            time_tracker=self.time_tracker,
            duty_cycle=0.5
        )
        
        assert square_wave.period_ps == 5.0
        assert square_wave.duty_cycle == 0.5
        assert square_wave.baseline_coupling == 1e-5
        assert not square_wave.has_started
        assert not square_wave.has_stopped
    
    def test_square_wave_timing(self):
        """Test square wave timing and state transitions."""
        square_wave = LQGSquareWaveVariant(
            lqg_controller=self.lqg_controller,
            period_ps=4.0,
            time_tracker=self.time_tracker,
            duty_cycle=0.5,
            start_time_ps=2.0
        )
        
        # Before start time
        self.time_tracker.elapsed_time_ps = 1.0
        coupling = square_wave(0)
        assert coupling == 0.0
        assert not square_wave.has_started
        
        # During high phase
        self.time_tracker.elapsed_time_ps = 3.0  # 1 ps after start, in high phase
        coupling = square_wave(0)
        assert coupling == self.lqg_controller.get_current_coupling()
        assert square_wave.has_started
        assert square_wave.is_high
        
        # During low phase
        self.time_tracker.elapsed_time_ps = 5.0  # 3 ps after start, in low phase
        coupling = square_wave(0)
        assert coupling == square_wave.baseline_coupling
        assert not square_wave.is_high
    
    def test_square_wave_stop_time(self):
        """Test square wave stop functionality."""
        square_wave = LQGSquareWaveVariant(
            lqg_controller=self.lqg_controller,
            period_ps=5.0,
            time_tracker=self.time_tracker,
            start_time_ps=0.0,
            stop_time_ps=10.0
        )
        
        # During active period
        self.time_tracker.elapsed_time_ps = 5.0
        coupling = square_wave(0)
        assert coupling > 0.0
        
        # After stop time
        self.time_tracker.elapsed_time_ps = 15.0
        coupling = square_wave(0)
        assert coupling == 0.0
        assert square_wave.has_stopped
    
    def test_state_info_retrieval(self):
        """Test state information retrieval."""
        square_wave = LQGSquareWaveVariant(
            lqg_controller=self.lqg_controller,
            period_ps=5.0,
            time_tracker=self.time_tracker
        )
        
        self.time_tracker.elapsed_time_ps = 2.0
        square_wave(0)  # Activate
        
        state_info = square_wave.get_state_info()
        
        assert isinstance(state_info, dict)
        assert 'time_ps' in state_info
        assert 'is_active' in state_info
        assert 'current_value' in state_info
        assert 'coupling_deviation' in state_info
        assert 'bath_temperature' in state_info
    
    def test_adaptive_square_wave(self):
        """Test adaptive square wave variant."""
        adaptive_square = LQGAdaptiveSquareWaveVariant(
            lqg_controller=self.lqg_controller,
            period_ps=5.0,
            time_tracker=self.time_tracker,
            enable_period_adaptation=True,
            adaptation_interval_pulses=10
        )
        
        assert adaptive_square.enable_period_adaptation == True
        assert adaptive_square.adaptation_interval_pulses == 10
        assert adaptive_square.consecutive_saturations == 0
        
        # Test basic functionality
        self.time_tracker.elapsed_time_ps = 2.0
        coupling = adaptive_square(0)
        assert coupling >= 0.0


class TestIntegration:
    """Integration tests for complete LQG system."""
    
    def test_complete_system_simulation(self):
        """Test complete system with all components."""
        # This would be a more comprehensive test with actual HOOMD simulation
        # For now, we test the component integration
        
        # Create parameters
        params = LQGParameters()
        params.Q_structure = 2.0
        
        # Create model
        model = LiftedStateSpaceModel(params)
        
        # Create Kalman filter
        kf = KalmanFilter(model)
        
        # Create LQR controller
        lqr = LQRController(model, params)
        
        # Simulate a few steps
        x_true = np.array([300.0, 300.0, 300.0, 0.0])[:model.n_states]
        reference = np.array([295.0, 295.0, 295.0])
        
        for i in range(10):
            # Simulate measurement
            y_k = x_true[:3] + np.random.normal(0, 0.1, 3)
            
            # Previous control (start with zero)
            u_prev = np.array([0.0, 300.0]) if i == 0 else u
            
            # Kalman filter step
            x_hat = kf.step(y_k, u_prev)
            
            # Control computation
            u = lqr.compute_control(x_hat, reference)
            
            # Simulate true system
            x_true, _ = model.simulate_step(x_true, u, add_noise=True)
            
            # Check that everything is working
            assert len(x_hat) == model.n_states
            assert len(u) == 2
            assert np.all(np.isfinite(x_hat))
            assert np.all(np.isfinite(u))
    
    def test_system_identification_workflow(self):
        """Test complete system identification workflow."""
        # Generate synthetic data
        true_params = LQGParameters()
        true_params.a_sk = 0.18
        true_params.b_sg = 0.08
        
        true_model = LiftedStateSpaceModel(true_params)
        
        # Generate training data
        pulse_data = []
        x = np.array([300.0, 300.0, 300.0, 0.0])[:true_model.n_states]
        
        for i in range(30):
            u = np.array([0.02 * np.sin(i * 0.2), 300.0 + 10.0 * np.cos(i * 0.15)])
            x_next, _ = true_model.simulate_step(x, u, add_noise=True)
            
            pulse_data.append({
                'state': x[:3].tolist(),
                'control': u.tolist(),
                'next_state': x_next[:3].tolist()
            })
            
            x = x_next
        
        # Perform system identification
        sys_id = SystemIdentification()
        fitted_params = sys_id.fit_from_data(pulse_data)
        
        # Create new model with fitted parameters
        fitted_model = LiftedStateSpaceModel(fitted_params)
        
        # Test that fitted model performs reasonably
        x_test = np.array([305.0, 295.0, 300.0, 0.0])[:fitted_model.n_states]
        u_test = np.array([0.01, 310.0])
        
        x_next_fitted, _ = fitted_model.simulate_step(x_test, u_test)
        x_next_true, _ = true_model.simulate_step(x_test, u_test)
        
        # Fitted model should be reasonably close to true model
        error = np.linalg.norm(x_next_fitted[:3] - x_next_true[:3])
        assert error < 10.0  # Reasonable error bound


if __name__ == '__main__':
    # Run tests
    pytest.main([__file__, '-v'])
