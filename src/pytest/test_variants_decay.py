"""Test the decay variant classes from the cavitymd plugin.

This module tests time-decay coupling protocols including:
- ExponentialDecayVariant: Exponential decay to baseline
- DecayingSquareWaveVariant: Square wave with exponential amplitude decay
"""

import pytest
import numpy as np
import hoomd
from unittest.mock import Mock, MagicMock

# Try to import from the plugin
try:
    from hoomd.cavitymd.variants import (
        ExponentialDecayVariant,
        DecayingSquareWaveVariant
    )
    from hoomd.cavitymd.analysis import ElapsedTimeTracker
    from hoomd.cavitymd.utils import PhysicalConstants
    PLUGIN_AVAILABLE = True
except ImportError:
    PLUGIN_AVAILABLE = False


@pytest.mark.skipif(not PLUGIN_AVAILABLE, reason='Plugin not available')
class TestExponentialDecayVariant:
    """Test suite for ExponentialDecayVariant class."""
    
    @pytest.fixture
    def mock_time_tracker(self):
        """Create a mock time tracker for testing."""
        tracker = Mock(spec=ElapsedTimeTracker)
        tracker.elapsed_time = 0.0
        return tracker
    
    def test_initialization_basic(self, mock_time_tracker):
        """Test basic initialization of ExponentialDecayVariant."""
        variant = ExponentialDecayVariant(
            initial_value=0.01,
            final_value=0.001,
            time_tracker=mock_time_tracker,
            decay_time_ps=100.0
        )
        
        assert variant.initial_value == 0.01
        assert variant.final_value == 0.001
        assert variant.decay_time_ps == 100.0
        assert variant.time_tracker is mock_time_tracker
    
    def test_initialization_with_start_time(self, mock_time_tracker):
        """Test initialization with delayed start time."""
        variant = ExponentialDecayVariant(
            initial_value=0.01,
            final_value=0.001,
            time_tracker=mock_time_tracker,
            decay_time_ps=100.0,
            start_time_ps=10.0
        )
        
        assert variant.start_time_ps == 10.0
    
    def test_value_at_t_equals_0(self, mock_time_tracker):
        """Test value at t=0 equals initial value."""
        variant = ExponentialDecayVariant(
            initial_value=0.01,
            final_value=0.001,
            time_tracker=mock_time_tracker,
            decay_time_ps=100.0,
            start_time_ps=0.0
        )
        
        mock_time_tracker.elapsed_time = 0.0
        value = variant(0)
        assert abs(value - 0.01) < 1e-10
    
    def test_value_after_long_time(self, mock_time_tracker):
        """Test value approaches final value after long time."""
        variant = ExponentialDecayVariant(
            initial_value=0.01,
            final_value=0.001,
            time_tracker=mock_time_tracker,
            decay_time_ps=100.0,
            start_time_ps=0.0
        )
        
        # After 5 time constants, should be very close to final value
        # exp(-5) ≈ 0.0067, so essentially at final value
        mock_time_tracker.elapsed_time = 500.0
        value = variant(0)
        assert abs(value - 0.001) < 1e-6
    
    def test_exponential_decay_behavior(self, mock_time_tracker):
        """Test exponential decay follows expected curve."""
        variant = ExponentialDecayVariant(
            initial_value=0.01,
            final_value=0.001,
            time_tracker=mock_time_tracker,
            decay_time_ps=100.0,
            start_time_ps=0.0
        )
        
        # Sample at multiple time points
        times = [0, 50, 100, 200, 500]
        values = []
        
        for t in times:
            mock_time_tracker.elapsed_time = t
            values.append(variant(0))
        
        values = np.array(values)
        
        # Check monotonic decrease
        assert all(values[i] >= values[i+1] for i in range(len(values)-1))
        
        # Check initial and final values
        assert abs(values[0] - 0.01) < 1e-10
        assert abs(values[-1] - 0.001) < 1e-6
        
        # At t = decay_time, should be at ~63% decay
        # value(t) = final + (initial - final) * exp(-t/tau)
        # At t = tau: value = final + (initial - final) * exp(-1)
        # For initial=0.01, final=0.001: value(100) = 0.001 + 0.009 * 0.368 ≈ 0.00433
        mock_time_tracker.elapsed_time = 100.0
        value_at_tau = variant(0)
        expected = 0.001 + (0.01 - 0.001) * np.exp(-1)
        assert abs(value_at_tau - expected) < 1e-6
    
    def test_before_start_time(self, mock_time_tracker):
        """Test behavior before start time."""
        variant = ExponentialDecayVariant(
            initial_value=0.01,
            final_value=0.001,
            time_tracker=mock_time_tracker,
            decay_time_ps=100.0,
            start_time_ps=50.0
        )
        
        # Before start, should return initial value
        mock_time_tracker.elapsed_time = 25.0
        value = variant(0)
        assert abs(value - 0.01) < 1e-10


@pytest.mark.skipif(not PLUGIN_AVAILABLE, reason='Plugin not available')
class TestDecayingSquareWaveVariant:
    """Test suite for DecayingSquareWaveVariant class."""
    
    @pytest.fixture
    def mock_time_tracker(self):
        """Create a mock time tracker for testing."""
        tracker = Mock(spec=ElapsedTimeTracker)
        tracker.elapsed_time = 0.0
        return tracker
    
    def test_initialization_basic(self, mock_time_tracker):
        """Test basic initialization of DecayingSquareWaveVariant."""
        variant = DecayingSquareWaveVariant(
            initial_high=0.01,
            final_high=0.001,
            g_low=0.0001,
            time_tracker=mock_time_tracker,
            period_ps=10.0,
            duty_cycle=0.5,
            decay_time_ps=100.0
        )
        
        assert variant.initial_high == 0.01
        assert variant.final_high == 0.001
        assert variant.g_low == 0.0001
        assert variant.period_ps == 10.0
        assert variant.duty_cycle == 0.5
        assert variant.decay_time_ps == 100.0
        assert variant.time_tracker is mock_time_tracker
    
    def test_high_value_decays_over_time(self, mock_time_tracker):
        """Test that HIGH value decays exponentially."""
        variant = DecayingSquareWaveVariant(
            initial_high=0.01,
            final_high=0.001,
            g_low=0.0001,
            time_tracker=mock_time_tracker,
            period_ps=10.0,
            duty_cycle=0.5,
            decay_time_ps=100.0,
            start_time_ps=0.0
        )
        
        # At t=2.5ps (HIGH state), should be close to initial
        mock_time_tracker.elapsed_time = 2.5
        value_early = variant(0)
        assert abs(value_early - 0.01) < 1e-6
        
        # At t=252.5ps (HIGH state after decay), should be closer to final
        mock_time_tracker.elapsed_time = 252.5
        value_late = variant(0)
        
        # Should have decayed significantly
        assert value_late < value_early
        assert value_late > 0.001
    
    def test_low_value_remains_constant(self, mock_time_tracker):
        """Test that LOW value remains constant over time."""
        variant = DecayingSquareWaveVariant(
            initial_high=0.01,
            final_high=0.001,
            g_low=0.0001,
            time_tracker=mock_time_tracker,
            period_ps=10.0,
            duty_cycle=0.5,
            decay_time_ps=100.0,
            start_time_ps=0.0
        )
        
        # At t=7.5ps (LOW state, early)
        mock_time_tracker.elapsed_time = 7.5
        value_early = variant(0)
        assert abs(value_early - 0.0001) < 1e-10
        
        # At t=257.5ps (LOW state, late)
        mock_time_tracker.elapsed_time = 257.5
        value_late = variant(0)
        assert abs(value_late - 0.0001) < 1e-10
        
        # LOW value should not change
        assert abs(value_early - value_late) < 1e-10
    
    def test_square_wave_continues_cycling(self, mock_time_tracker):
        """Test that square wave continues to cycle as amplitude decays."""
        variant = DecayingSquareWaveVariant(
            initial_high=0.01,
            final_high=0.001,
            g_low=0.0001,
            time_tracker=mock_time_tracker,
            period_ps=10.0,
            duty_cycle=0.5,
            decay_time_ps=100.0,
            start_time_ps=0.0
        )
        
        # First cycle HIGH
        mock_time_tracker.elapsed_time = 2.5
        value1_high = variant(0)
        
        # First cycle LOW
        mock_time_tracker.elapsed_time = 7.5
        value1_low = variant(0)
        
        # Second cycle HIGH (amplitude should have decayed slightly)
        mock_time_tracker.elapsed_time = 12.5
        value2_high = variant(0)
        
        # Second cycle LOW
        mock_time_tracker.elapsed_time = 17.5
        value2_low = variant(0)
        
        # HIGH values should both be above LOW
        assert value1_high > value1_low
        assert value2_high > value2_low
        
        # Second HIGH should be less than first due to decay
        assert value2_high < value1_high
        
        # LOW values should be equal
        assert abs(value1_low - value2_low) < 1e-10
    
    def test_asymmetric_duty_cycle_with_decay(self, mock_time_tracker):
        """Test decaying square wave with asymmetric duty cycle."""
        variant = DecayingSquareWaveVariant(
            initial_high=0.01,
            final_high=0.001,
            g_low=0.0001,
            time_tracker=mock_time_tracker,
            period_ps=10.0,
            duty_cycle=0.3,  # 30% HIGH, 70% LOW
            decay_time_ps=100.0,
            start_time_ps=0.0
        )
        
        # Within 30% of period, should be HIGH
        mock_time_tracker.elapsed_time = 1.5
        value_high = variant(0)
        assert value_high > 0.001
        
        # Beyond 30% of period, should be LOW
        mock_time_tracker.elapsed_time = 5.0
        value_low = variant(0)
        assert abs(value_low - 0.0001) < 1e-10


@pytest.mark.skipif(not PLUGIN_AVAILABLE, reason='Plugin not available')
class TestDecayVariantIntegration:
    """Integration tests for decay variants with HOOMD simulation."""
    
    def test_exponential_decay_in_simulation(self):
        """Test that ExponentialDecayVariant works in actual simulation."""
        device = hoomd.device.CPU()
        sim = hoomd.Simulation(device=device, seed=42)
        
        # Create simple system
        snapshot = hoomd.Snapshot(device.communicator)
        snapshot.configuration.box = [10, 10, 10, 0, 0, 0]
        snapshot.particles.N = 2
        snapshot.particles.types = ['A']
        snapshot.particles.position[0] = [0, 0, 0]
        snapshot.particles.position[1] = [2.0, 0, 0]
        snapshot.particles.velocity[0] = [0, 0, 0]
        snapshot.particles.velocity[1] = [0, 0, 0]
        snapshot.particles.mass[0] = 1.0
        snapshot.particles.mass[1] = 1.0
        
        sim.create_state_from_snapshot(snapshot)
        
        # Create time tracker
        time_tracker = ElapsedTimeTracker()
        
        # Create decay variant
        variant = ExponentialDecayVariant(
            initial_value=0.01,
            final_value=0.001,
            time_tracker=time_tracker,
            decay_time_ps=100.0
        )
        
        # Set up integrator
        integrator = hoomd.md.Integrator(dt=0.001)
        integrator.methods.append(
            hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All())
        )
        sim.operations.integrator = integrator
        
        # Attach time tracker
        sim.operations.writers.append(time_tracker)
        
        # Run simulation
        sim.run(100)
        
        # Verify variant can be called and returns reasonable values
        value = variant(0)
        assert isinstance(value, (int, float))
        assert 0.001 <= value <= 0.01


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

