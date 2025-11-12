"""Test the periodic variant classes from the cavitymd plugin.

This module tests periodic time-varying coupling protocols including:
- PeriodicVariant: Sinusoidal oscillations
- SquareWaveVariant: Square wave modulation
"""

import pytest
import numpy as np
import hoomd
from unittest.mock import Mock, MagicMock

# Try to import from the plugin
try:
    from hoomd.cavitymd.variants import PeriodicVariant, SquareWaveVariant
    from hoomd.cavitymd.analysis import ElapsedTimeTracker
    from hoomd.cavitymd.utils import PhysicalConstants
    PLUGIN_AVAILABLE = True
except ImportError:
    PLUGIN_AVAILABLE = False


@pytest.mark.skipif(not PLUGIN_AVAILABLE, reason='Plugin not available')
class TestPeriodicVariant:
    """Test suite for PeriodicVariant class."""
    
    @pytest.fixture
    def mock_time_tracker(self):
        """Create a mock time tracker for testing."""
        tracker = Mock(spec=ElapsedTimeTracker)
        tracker.elapsed_time = 0.0  # Start at t=0
        return tracker
    
    def test_initialization_basic(self, mock_time_tracker):
        """Test basic initialization of PeriodicVariant."""
        variant = PeriodicVariant(
            amplitude=0.001,
            time_tracker=mock_time_tracker,
            period_ps=1.0
        )
        
        assert variant.amplitude == 0.001
        assert variant.period_ps == 1.0
        assert variant.phase_offset == 0.0
        assert variant.start_time_ps == 0.0
        assert variant.time_tracker is mock_time_tracker
        
        # Check derived quantities
        assert abs(variant.frequency_hz - 1e12) < 1e9  # 1 THz ± 1 GHz
        assert abs(variant.angular_frequency - 2*np.pi) < 1e-10
        
        # Check initial state
        assert variant.current_value == 0.0
        assert not variant.oscillation_started
    
    def test_initialization_with_phase(self, mock_time_tracker):
        """Test initialization with phase offset."""
        variant = PeriodicVariant(
            amplitude=0.001,
            time_tracker=mock_time_tracker,
            period_ps=1.0,
            phase_offset=np.pi/2  # 90 degree phase shift
        )
        
        assert variant.phase_offset == np.pi/2
    
    def test_initialization_with_start_time(self, mock_time_tracker):
        """Test initialization with delayed start time."""
        variant = PeriodicVariant(
            amplitude=0.001,
            time_tracker=mock_time_tracker,
            period_ps=1.0,
            start_time_ps=10.0
        )
        
        assert variant.start_time_ps == 10.0
        assert not variant.oscillation_started
    
    def test_oscillation_value_at_t_equals_0(self, mock_time_tracker):
        """Test oscillation value at t=0."""
        variant = PeriodicVariant(
            amplitude=0.001,
            time_tracker=mock_time_tracker,
            period_ps=1.0,
            phase_offset=0.0
        )
        
        # At t=0, sin(0) = 0
        value = variant(0)
        assert abs(value) < 1e-10
    
    def test_oscillation_value_at_quarter_period(self, mock_time_tracker):
        """Test oscillation value at t=T/4."""
        variant = PeriodicVariant(
            amplitude=0.001,
            time_tracker=mock_time_tracker,
            period_ps=1.0,
            phase_offset=0.0
        )
        
        # Simulate time passing to T/4
        mock_time_tracker.elapsed_time = 0.25
        
        # At T/4, sin(π/2) = 1, so value should be amplitude
        value = variant(0)  # timestep arg not used, uses time_tracker
        assert abs(value - 0.001) < 1e-10
    
    def test_oscillation_value_at_half_period(self, mock_time_tracker):
        """Test oscillation value at t=T/2."""
        variant = PeriodicVariant(
            amplitude=0.001,
            time_tracker=mock_time_tracker,
            period_ps=1.0,
            phase_offset=0.0
        )
        
        # Simulate time passing to T/2
        mock_time_tracker.elapsed_time = 0.5
        
        # At T/2, sin(π) ≈ 0
        value = variant(0)
        assert abs(value) < 1e-10
    
    def test_oscillation_full_cycle(self, mock_time_tracker):
        """Test oscillation over a full period."""
        variant = PeriodicVariant(
            amplitude=0.001,
            time_tracker=mock_time_tracker,
            period_ps=1.0,
            phase_offset=0.0
        )
        
        # Sample at multiple points in the period
        times = np.linspace(0, 1.0, 20)
        values = []
        
        for t in times:
            mock_time_tracker.elapsed_time = t
            values.append(variant(0))
        
        values = np.array(values)
        
        # Check that max value is approximately amplitude
        assert abs(np.max(values) - 0.001) < 1e-6
        
        # Check that min value is approximately -amplitude
        assert abs(np.min(values) + 0.001) < 1e-6
        
        # Check oscillation completes a cycle (value at t=0 ≈ value at t=T)
        mock_time_tracker.elapsed_time = 1.0
        final_value = variant(0)
        assert abs(final_value - 0.0) < 1e-6


@pytest.mark.skipif(not PLUGIN_AVAILABLE, reason='Plugin not available')
class TestSquareWaveVariant:
    """Test suite for SquareWaveVariant class."""
    
    @pytest.fixture
    def mock_time_tracker(self):
        """Create a mock time tracker for testing."""
        tracker = Mock(spec=ElapsedTimeTracker)
        tracker.elapsed_time = 0.0
        return tracker
    
    def test_initialization_basic(self, mock_time_tracker):
        """Test basic initialization of SquareWaveVariant."""
        variant = SquareWaveVariant(
            g_high=0.001,
            g_low=0.0001,
            time_tracker=mock_time_tracker,
            period_ps=10.0,
            duty_cycle=0.5
        )
        
        assert variant.g_high == 0.001
        assert variant.g_low == 0.0001
        assert variant.period_ps == 10.0
        assert variant.duty_cycle == 0.5
        assert variant.time_tracker is mock_time_tracker
    
    def test_initialization_with_start_stop_times(self, mock_time_tracker):
        """Test initialization with start and stop times."""
        variant = SquareWaveVariant(
            g_high=0.001,
            g_low=0.0001,
            time_tracker=mock_time_tracker,
            period_ps=10.0,
            duty_cycle=0.5,
            start_time_ps=10.0,
            stop_time_ps=100.0
        )
        
        assert variant.start_time_ps == 10.0
        assert variant.stop_time_ps == 100.0
    
    def test_high_state_during_duty_cycle(self, mock_time_tracker):
        """Test that variant is HIGH during duty cycle."""
        variant = SquareWaveVariant(
            g_high=0.001,
            g_low=0.0001,
            time_tracker=mock_time_tracker,
            period_ps=10.0,
            duty_cycle=0.5,
            start_time_ps=0.0
        )
        
        # At t=2.5ps (within first half of 10ps period), should be HIGH
        mock_time_tracker.elapsed_time = 2.5
        value = variant(0)
        assert value == 0.001
    
    def test_low_state_after_duty_cycle(self, mock_time_tracker):
        """Test that variant is LOW after duty cycle."""
        variant = SquareWaveVariant(
            g_high=0.001,
            g_low=0.0001,
            time_tracker=mock_time_tracker,
            period_ps=10.0,
            duty_cycle=0.5,
            start_time_ps=0.0
        )
        
        # At t=7.5ps (second half of 10ps period), should be LOW
        mock_time_tracker.elapsed_time = 7.5
        value = variant(0)
        assert value == 0.0001
    
    def test_cycle_repeats(self, mock_time_tracker):
        """Test that square wave cycles repeat."""
        variant = SquareWaveVariant(
            g_high=0.001,
            g_low=0.0001,
            time_tracker=mock_time_tracker,
            period_ps=10.0,
            duty_cycle=0.5,
            start_time_ps=0.0
        )
        
        # First cycle: HIGH at t=2.5ps
        mock_time_tracker.elapsed_time = 2.5
        assert variant(0) == 0.001
        
        # First cycle: LOW at t=7.5ps
        mock_time_tracker.elapsed_time = 7.5
        assert variant(0) == 0.0001
        
        # Second cycle: HIGH at t=12.5ps (2.5 + 10)
        mock_time_tracker.elapsed_time = 12.5
        assert variant(0) == 0.001
        
        # Second cycle: LOW at t=17.5ps (7.5 + 10)
        mock_time_tracker.elapsed_time = 17.5
        assert variant(0) == 0.0001
    
    def test_before_start_time(self, mock_time_tracker):
        """Test behavior before start time."""
        variant = SquareWaveVariant(
            g_high=0.001,
            g_low=0.0001,
            time_tracker=mock_time_tracker,
            period_ps=10.0,
            duty_cycle=0.5,
            start_time_ps=10.0
        )
        
        # Before start, should return g_low
        mock_time_tracker.elapsed_time = 5.0
        value = variant(0)
        assert value == 0.0001
    
    def test_after_stop_time(self, mock_time_tracker):
        """Test behavior after stop time."""
        variant = SquareWaveVariant(
            g_high=0.001,
            g_low=0.0001,
            time_tracker=mock_time_tracker,
            period_ps=10.0,
            duty_cycle=0.5,
            start_time_ps=0.0,
            stop_time_ps=50.0
        )
        
        # After stop, should return g_low
        mock_time_tracker.elapsed_time = 60.0
        value = variant(0)
        assert value == 0.0001
    
    def test_asymmetric_duty_cycle(self, mock_time_tracker):
        """Test square wave with asymmetric duty cycle."""
        variant = SquareWaveVariant(
            g_high=0.001,
            g_low=0.0001,
            time_tracker=mock_time_tracker,
            period_ps=10.0,
            duty_cycle=0.3,  # 30% HIGH, 70% LOW
            start_time_ps=0.0
        )
        
        # At t=1.5ps (within 30% of period), should be HIGH
        mock_time_tracker.elapsed_time = 1.5
        assert variant(0) == 0.001
        
        # At t=5.0ps (beyond 30% of period), should be LOW
        mock_time_tracker.elapsed_time = 5.0
        assert variant(0) == 0.0001


@pytest.mark.skipif(not PLUGIN_AVAILABLE, reason='Plugin not available')
class TestPeriodicVariantIntegration:
    """Integration tests for periodic variants with HOOMD simulation."""
    
    def test_periodic_variant_in_simulation(self):
        """Test that PeriodicVariant works in actual simulation."""
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
        
        # Create periodic variant
        variant = PeriodicVariant(
            amplitude=0.001,
            time_tracker=time_tracker,
            period_ps=1.0
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
        
        # Verify variant can be called
        value = variant(0)
        assert isinstance(value, (int, float))


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

