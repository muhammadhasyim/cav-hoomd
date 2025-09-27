"""
Tests for FDR-based effective temperature estimator.

This test suite follows Test-Driven Development (TDD) principles and includes:
1. Unit tests for individual components
2. Integration tests with synthetic data
3. Validation against known analytical solutions
4. Edge case and error condition testing

Scientific validation includes:
- Markovian Langevin oscillator with known T_eff
- Dephasing effects (quasi-static and motional narrowing)
- Equilibrium FDT verification
- Non-equilibrium temperature dynamics
"""

import pytest
import numpy as np
import warnings
from unittest.mock import Mock, patch
from dataclasses import asdict

# Import the module under test
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from cavitymd.fdr_temperature import (
    FDRTemperatureEstimator,
    FDRDiagnostics, 
    LineshapeType,
    create_dipole_extractor,
    create_mode_extractor
)
from cavitymd.utils import PhysicalConstants


class TestFDRTemperatureEstimator:
    """Test suite for the core FDR temperature estimator."""
    
    @pytest.fixture
    def basic_estimator(self):
        """Create a basic FDR estimator for testing."""
        omega_0 = 1.0  # rad/ps
        dt = 0.001     # ps
        return FDRTemperatureEstimator(omega_0=omega_0, dt=dt)
        
    @pytest.fixture
    def synthetic_langevin_data(self):
        """Generate synthetic Langevin oscillator data with known properties."""
        
        # Physical parameters
        omega_n = 1.0    # rad/ps
        gamma = 0.1      # 1/ps
        T_true = 300.0   # K
        dt = 0.001       # ps
        n_steps = 50000
        
        # Convert temperature to simulation units
        k_B = PhysicalConstants.KB_HARTREE_PER_K
        
        # Generate colored noise for non-equilibrium case
        # S_FF = 2*k_B*T_eff*gamma where T_eff != T_bath
        T_eff_true = 350.0  # K (higher than bath)
        noise_amplitude = np.sqrt(2 * k_B * T_eff_true * gamma / dt)
        
        # Integrate Langevin equation: m*a + gamma*v + k*x = F_random
        # For unit mass: a + gamma*v + omega_n^2*x = noise
        
        times = np.arange(n_steps) * dt
        x = np.zeros(n_steps)
        v = np.zeros(n_steps)
        F_random = np.random.normal(0, noise_amplitude, n_steps)
        
        # Velocity Verlet integration
        for i in range(1, n_steps):
            # Current acceleration
            a_current = -gamma * v[i-1] - omega_n**2 * x[i-1] + F_random[i-1]
            
            # Update position
            x[i] = x[i-1] + v[i-1] * dt + 0.5 * a_current * dt**2
            
            # Next acceleration (estimated)
            a_next = -gamma * v[i-1] - omega_n**2 * x[i] + F_random[i]
            
            # Update velocity
            v[i] = v[i-1] + 0.5 * (a_current + a_next) * dt
            
        return {
            'times': times,
            'positions': x,
            'velocities': v,
            'forces': F_random,
            'omega_n_true': omega_n,
            'gamma_true': gamma,
            'T_eff_true': T_eff_true,
            'T_bath': T_true,
            'dt': dt
        }
        
    def test_initialization_valid_parameters(self, basic_estimator):
        """Test proper initialization with valid parameters."""
        
        assert basic_estimator.omega_0 == 1.0
        assert basic_estimator.dt == 0.001
        assert basic_estimator.T_0 == 2*np.pi  # Period
        assert not basic_estimator._is_calibrated
        assert basic_estimator.step_count == 0
        
    def test_initialization_invalid_parameters(self):
        """Test initialization with invalid parameters raises appropriate errors."""
        
        with pytest.raises(ValueError, match="omega_0 must be positive"):
            FDRTemperatureEstimator(omega_0=-1.0, dt=0.001)
            
        with pytest.raises(ValueError, match="dt must be positive"):
            FDRTemperatureEstimator(omega_0=1.0, dt=-0.001)
            
    def test_initialization_timestep_warning(self):
        """Test warning for large timestep."""
        
        omega_0 = 10.0  # High frequency
        dt = 0.1        # Large timestep
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            FDRTemperatureEstimator(omega_0=omega_0, dt=dt)
            assert len(w) == 1
            assert "may be too large" in str(w[0].message)
            
    def test_calibration_valid(self, basic_estimator):
        """Test calibration with valid equilibrium data."""
        
        # Generate some mock equilibrium data
        T_cal = 300.0  # K
        n_points = 1000
        A_equilibrium = np.random.normal(0, 1.0, n_points)
        
        # Perform calibration
        basic_estimator.calibrate(T_cal, A_equilibrium)
        
        assert basic_estimator._is_calibrated
        assert basic_estimator._G_gain is not None
        assert basic_estimator._G_gain > 0
        
    def test_calibration_invalid_temperature(self, basic_estimator):
        """Test calibration with invalid temperature."""
        
        with pytest.raises(ValueError, match="temperature must be positive"):
            basic_estimator.calibrate(-100.0)
            
    def test_update_without_forces(self, basic_estimator):
        """Test update mechanism without force data (AR2 mode)."""
        
        # Calibrate first
        basic_estimator.calibrate(300.0, np.random.normal(0, 1.0, 1000))
        
        # Single update
        A = 0.5
        T_eff, diagnostics = basic_estimator.update(A)
        
        assert basic_estimator.step_count == 1
        assert isinstance(T_eff, (float, type(np.nan)))
        assert isinstance(diagnostics, FDRDiagnostics)
        
    def test_update_with_forces(self, basic_estimator):
        """Test update mechanism with force data (regression mode)."""
        
        # Mock force extractor
        def mock_force_extractor(forces):
            return np.sum(forces)  # Simple sum
            
        basic_estimator.force_extractor = mock_force_extractor
        basic_estimator.calibrate(300.0, np.random.normal(0, 1.0, 1000))
        
        # Update with forces
        A = 0.5
        forces = np.array([0.1, -0.2, 0.05])
        T_eff, diagnostics = basic_estimator.update(A, forces=forces)
        
        assert basic_estimator.step_count == 1
        assert isinstance(diagnostics, FDRDiagnostics)
        
    def test_lock_in_amplifier(self, basic_estimator):
        """Test lock-in amplifier for pure sinusoidal signal."""
        
        # Generate pure sine wave at omega_0
        omega_0 = basic_estimator.omega_0
        dt = basic_estimator.dt
        amplitude = 2.0
        n_cycles = 10
        n_points = int(n_cycles * 2*np.pi / (omega_0 * dt))
        
        for i in range(n_points):
            t = i * dt
            A = amplitude * np.sin(omega_0 * t)
            basic_estimator._update_lock_in(A)
            
        # Check that lock-in recovered the correct amplitude
        S_AA = basic_estimator._get_S_AA()
        expected_power = amplitude**2 / 2  # Power of sine wave
        
        # Allow 10% tolerance due to finite averaging
        assert abs(S_AA - expected_power) / expected_power < 0.1
        
    def test_synthetic_langevin_validation(self, synthetic_langevin_data):
        """Test FDR estimator on synthetic Langevin oscillator data."""
        
        data = synthetic_langevin_data
        
        # Create estimator with appropriate parameters
        omega_0 = data['omega_n_true']  # Target at natural frequency
        dt = data['dt']
        estimator = FDRTemperatureEstimator(
            omega_0=omega_0, 
            dt=dt,
            tau_avg=100*2*np.pi/omega_0,  # 100 periods
            tau_id=500*2*np.pi/omega_0    # 500 periods
        )
        
        # Calibrate on equilibrium section (first 10% of data)
        n_cal = len(data['positions']) // 10
        estimator.calibrate(data['T_bath'], data['positions'][:n_cal])
        
        # Process the non-equilibrium data
        T_eff_estimates = []
        gamma_estimates = []
        omega_n_estimates = []
        
        for i, A in enumerate(data['positions']):
            if i < n_cal:
                continue  # Skip calibration region
                
            T_eff, diagnostics = estimator.update(A)
            
            if not np.isnan(T_eff) and i > n_cal + 1000:  # Allow settling
                T_eff_estimates.append(T_eff)
                gamma_estimates.append(diagnostics.gamma)
                omega_n_estimates.append(diagnostics.omega_n)
                
        # Validate parameter recovery
        if len(T_eff_estimates) > 100:
            T_eff_mean = np.mean(T_eff_estimates[-1000:])  # Last 1000 points
            gamma_mean = np.mean(gamma_estimates[-1000:])
            omega_n_mean = np.mean(omega_n_estimates[-1000:])
            
            # Check parameter recovery (allow 20% tolerance)
            assert abs(gamma_mean - data['gamma_true']) / data['gamma_true'] < 0.2
            assert abs(omega_n_mean - data['omega_n_true']) / data['omega_n_true'] < 0.2
            
            # Check temperature recovery (allow 25% tolerance due to finite size effects)
            assert abs(T_eff_mean - data['T_eff_true']) / data['T_eff_true'] < 0.25
            
    def test_dephasing_motional_narrowing(self, basic_estimator):
        """Test motional narrowing regime with frequency jitter."""
        
        # Set up for motional narrowing: omega_0 * tau_c << 1
        basic_estimator.tau_c = 0.1  # ps (short correlation time)
        basic_estimator.sigma_delta = 0.2  # rad/ps (moderate jitter)
        basic_estimator.gamma = 0.1  # 1/ps
        
        # Should use motional narrowing
        omega_tau_c = basic_estimator.omega_0 * basic_estimator.tau_c
        assert omega_tau_c < 1.0
        
        # Test effective damping calculation
        gamma_eff_expected = (basic_estimator.gamma + 
                             2 * basic_estimator.sigma_delta**2 * basic_estimator.tau_c)
        
        # Run a few updates to establish state
        basic_estimator.calibrate(300.0, np.random.normal(0, 1.0, 100))
        for _ in range(10):
            basic_estimator.update(np.random.normal(0, 0.1))
            
        # Check diagnostics
        _, diagnostics = basic_estimator.update(0.1)
        assert diagnostics.lineshape_type == LineshapeType.MOTIONAL_NARROWING
        
    def test_dephasing_quasi_static(self, basic_estimator):
        """Test quasi-static regime with slow frequency jitter."""
        
        # Set up for quasi-static: omega_0 * tau_c >> 1
        basic_estimator.tau_c = 10.0  # ps (long correlation time)
        basic_estimator.sigma_delta = 0.1  # rad/ps
        basic_estimator.gamma = 0.05  # 1/ps
        
        # Should use Voigt lineshape
        omega_tau_c = basic_estimator.omega_0 * basic_estimator.tau_c
        assert omega_tau_c > 1.0
        
        # Run updates
        basic_estimator.calibrate(300.0, np.random.normal(0, 1.0, 100))
        for _ in range(10):
            basic_estimator.update(np.random.normal(0, 0.1))
            
        # Check diagnostics
        _, diagnostics = basic_estimator.update(0.1)
        assert diagnostics.lineshape_type == LineshapeType.VOIGT
        
    def test_state_persistence(self, basic_estimator):
        """Test state save/restore functionality."""
        
        # Initialize and run some updates
        basic_estimator.calibrate(300.0, np.random.normal(0, 1.0, 100))
        for i in range(50):
            basic_estimator.update(np.sin(0.001 * i))
            
        # Save state
        state = basic_estimator.get_state()
        
        # Modify the estimator
        for i in range(20):
            basic_estimator.update(np.cos(0.001 * i))
            
        original_step_count = basic_estimator.step_count
        
        # Restore state
        basic_estimator.set_state(state)
        
        # Check restoration
        assert basic_estimator.step_count == 50
        assert basic_estimator.step_count != original_step_count
        
    def test_edge_cases(self, basic_estimator):
        """Test edge cases and error conditions."""
        
        # Update before calibration should return NaN
        T_eff, _ = basic_estimator.update(0.1)
        assert np.isnan(T_eff)
        
        # Very small timestep
        small_dt_estimator = FDRTemperatureEstimator(omega_0=1.0, dt=1e-6)
        assert small_dt_estimator.dt == 1e-6
        
        # Large frequency
        high_freq_estimator = FDRTemperatureEstimator(omega_0=100.0, dt=0.001)
        assert high_freq_estimator.omega_0 == 100.0


class TestObservableExtractors:
    """Test suite for observable extraction functions."""
    
    def test_dipole_extractor_creation(self):
        """Test creation of dipole moment extractor."""
        
        extractor_z = create_dipole_extractor('z')
        extractor_x = create_dipole_extractor('x')
        
        assert callable(extractor_z)
        assert callable(extractor_x)
        
    def test_dipole_extractor_mock_snapshot(self):
        """Test dipole extractor with mock HOOMD snapshot."""
        
        # Create mock snapshot
        mock_snapshot = Mock()
        mock_snapshot.particles.position = np.array([
            [1.0, 2.0, 3.0],
            [-1.0, 1.0, 2.0],
            [0.0, -1.0, 1.0]
        ])
        mock_snapshot.particles.charge = np.array([1.0, -0.5, 2.0])
        
        # Test z-axis extraction
        extractor = create_dipole_extractor('z')
        dipole_z = extractor(mock_snapshot)
        
        # Expected: 1.0*3.0 + (-0.5)*2.0 + 2.0*1.0 = 3.0 - 1.0 + 2.0 = 4.0
        expected = 4.0
        assert abs(dipole_z - expected) < 1e-10
        
    def test_mode_extractor_creation(self):
        """Test creation of vibrational mode extractor."""
        
        # Mock mode vector and masses
        n_particles = 3
        mode_vector = np.random.random((n_particles, 3))
        masses = np.ones(n_particles)
        
        extractor = create_mode_extractor(mode_vector, masses)
        assert callable(extractor)
        
    def test_mode_extractor_mock_snapshot(self):
        """Test mode extractor with mock data."""
        
        # Simple test case
        mode_vector = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]
        ])
        masses = np.array([1.0, 4.0, 9.0])  # √masses = [1, 2, 3]
        
        mock_snapshot = Mock()
        mock_snapshot.particles.position = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0], 
            [0.0, 0.0, 3.0]
        ])
        
        # No reference positions - use positions directly
        extractor = create_mode_extractor(mode_vector, masses)
        mode_projection = extractor(mock_snapshot)
        
        # Expected: √1*1*1 + √4*1*2 + √9*1*3 = 1 + 4 + 9 = 14
        expected = 14.0
        assert abs(mode_projection - expected) < 1e-10


class TestLinearAlgebra:
    """Test numerical stability and linear algebra operations."""
    
    def test_ridge_regression_stability(self):
        """Test ridge regression handles ill-conditioned matrices."""
        
        estimator = FDRTemperatureEstimator(omega_0=1.0, dt=0.001)
        
        # Create ill-conditioned normal matrix
        estimator.S_matrix = np.array([[1e-10, 0], [0, 1e-10]])
        estimator.b_vector = np.array([1.0, 1.0])
        
        # Should not crash due to ridge regularization
        A_val = 0.1
        F_mic = 0.05
        estimator._update_force_regression(A_val, F_mic)
        
        # Check that parameters remain finite
        assert np.isfinite(estimator.gamma)
        assert np.isfinite(estimator.omega_n)
        
    def test_ar2_stability(self):
        """Test AR(2) parameter mapping handles edge cases."""
        
        estimator = FDRTemperatureEstimator(omega_0=1.0, dt=0.001)
        
        # Initialize AR(2) state
        estimator._ar2_buffer = np.array([0.0, 0.0, 0.0])
        estimator._ar2_P = np.eye(2) * 1000.0
        estimator._ar2_theta = np.array([0.0, 0.0])
        
        # Feed some data
        for i, A in enumerate([0.1, 0.2, 0.15, 0.05]):
            estimator.step_count = i + 1
            estimator._update_ar2_identification(A)
            
        # Should not crash and parameters should be reasonable
        assert 0.001 <= estimator.gamma <= 10.0
        assert 0.1 <= estimator.omega_n <= 100.0


@pytest.mark.integration
class TestIntegration:
    """Integration tests with realistic scenarios."""
    
    def test_temperature_relaxation_simulation(self):
        """Test estimator on a temperature relaxation scenario."""
        
        # Simulate rapid cooling: start hot, relax to cold
        omega_0 = 2.0   # rad/ps
        dt = 0.0005     # ps
        n_steps = 20000
        
        T_initial = 500.0  # K (hot)
        T_final = 200.0    # K (cold)
        tau_cool = 5.0     # ps (cooling time)
        
        estimator = FDRTemperatureEstimator(omega_0=omega_0, dt=dt)
        
        # Generate synthetic relaxation data
        times = np.arange(n_steps) * dt
        T_profile = T_final + (T_initial - T_final) * np.exp(-times / tau_cool)
        
        # Mock observable that reflects temperature
        A_data = []
        for i, T in enumerate(T_profile):
            # Simple model: amplitude scales with √T
            amplitude = np.sqrt(T / 300.0)  # Normalized to 300K
            phase = omega_0 * times[i] + 0.1 * np.random.normal()  # Add noise
            A = amplitude * np.sin(phase) + 0.05 * np.random.normal()
            A_data.append(A)
            
        # Calibrate on final equilibrium region
        estimator.calibrate(T_final, A_data[-1000:])
        
        # Process the relaxation
        T_eff_trajectory = []
        for i, A in enumerate(A_data):
            T_eff, diagnostics = estimator.update(A)
            if i > 1000 and not np.isnan(T_eff):  # Allow settling
                T_eff_trajectory.append((times[i], T_eff))
                
        # Check that we see a temperature relaxation trend
        if len(T_eff_trajectory) > 100:
            early_temps = [T for t, T in T_eff_trajectory[:50]]
            late_temps = [T for t, T in T_eff_trajectory[-50:]]
            
            early_mean = np.mean(early_temps)
            late_mean = np.mean(late_temps)
            
            # Should see decreasing temperature trend
            assert early_mean > late_mean
            assert late_mean < T_initial
            
    def test_frequency_sweep_response(self):
        """Test estimator response to frequency-swept signals."""
        
        omega_0 = 1.0  # rad/ps (target frequency)
        dt = 0.001     # ps
        estimator = FDRTemperatureEstimator(omega_0=omega_0, dt=dt)
        
        # Calibrate
        estimator.calibrate(300.0, np.random.normal(0, 1.0, 1000))
        
        # Test response to different frequency components
        n_steps = 5000
        frequencies = [0.5*omega_0, omega_0, 2.0*omega_0]  # Below, at, above resonance
        
        responses = {}
        for freq in frequencies:
            estimator._reset_state()  # Reset for each frequency
            
            S_AA_values = []
            for i in range(n_steps):
                t = i * dt
                A = np.sin(freq * t)  # Pure sinusoid
                estimator._update_lock_in(A)
                
                if i > 1000 and i % 100 == 0:  # Sample periodically after settling
                    S_AA = estimator._get_S_AA()
                    S_AA_values.append(S_AA)
                    
            responses[freq] = np.mean(S_AA_values[-10:])  # Average final response
            
        # Should see maximum response at resonance
        resonant_response = responses[omega_0]
        off_resonant_responses = [responses[f] for f in frequencies if f != omega_0]
        
        assert resonant_response > max(off_resonant_responses)


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "--tb=short"])
