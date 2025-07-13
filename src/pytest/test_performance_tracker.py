"""Test the PerformanceTracker class from the cavitymd plugin.

This module tests the performance tracking functionality including:
- Performance metrics calculation (ns/day, ETA)
- Integration with HOOMD simulation framework
- Logging capabilities
"""

from unittest.mock import patch

import hoomd

import pytest

# Try to import from the plugin
try:
    from hoomd.cavitymd import ElapsedTimeTracker, PerformanceTracker, PhysicalConstants

    PLUGIN_AVAILABLE = True
except ImportError:
    PLUGIN_AVAILABLE = False


@pytest.mark.skipif(not PLUGIN_AVAILABLE, reason='Plugin not available')
class TestPerformanceTracker:
    """Test suite for PerformanceTracker class."""

    def setup_method(self):
        """Set up test fixtures."""
        # Create a minimal simulation for testing
        self.device = hoomd.device.CPU()
        self.sim = hoomd.Simulation(device=self.device, seed=42)

        # Create a simple system (single particle for testing)
        self.snapshot = hoomd.Snapshot(self.device.communicator)
        self.snapshot.configuration.box = [10, 10, 10, 0, 0, 0]
        self.snapshot.particles.N = 1
        self.snapshot.particles.types = ['A']
        self.snapshot.particles.position[0] = [0, 0, 0]
        self.snapshot.particles.velocity[0] = [0, 0, 0]
        self.snapshot.particles.mass[0] = 1.0

        self.sim.create_state_from_snapshot(self.snapshot)

        # Set up a minimal integrator
        self.integrator = hoomd.md.Integrator(dt=0.001)
        self.integrator.methods.append(
            hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All())
        )
        self.sim.operations.integrator = self.integrator

        # Test parameters
        self.runtime_ps = 100.0

    def test_performance_tracker_initialization(self):
        """Test PerformanceTracker initialization."""
        tracker = PerformanceTracker(
            simulation=self.sim,
            runtime_ps=self.runtime_ps,
            time_tracker=None,
            output_prefix='test_performance',
            output_period_steps=1000,
        )

        assert tracker.sim == self.sim
        assert tracker.runtime_ps == self.runtime_ps
        assert tracker.current_ns_per_day == 0.0
        assert tracker.current_eta == ''
        assert tracker.output_prefix == 'test_performance'
        assert tracker.output_period_steps == 1000

    def test_performance_tracker_with_time_tracker(self):
        """Test PerformanceTracker with ElapsedTimeTracker."""
        time_tracker = ElapsedTimeTracker(simulation=self.sim, runtime=self.runtime_ps)

        tracker = PerformanceTracker(
            simulation=self.sim, runtime_ps=self.runtime_ps, time_tracker=time_tracker
        )

        assert tracker.time_tracker == time_tracker

    def test_performance_metrics_calculation(self):
        """Test performance metrics calculation."""
        tracker = PerformanceTracker(simulation=self.sim, runtime_ps=self.runtime_ps)

        # Mock time progression
        with patch('time.time') as mock_time:
            # Set initial time
            mock_time.return_value = 1000.0
            tracker.start_time = 1000.0

            # Simulate some time passing
            mock_time.return_value = 1010.0  # 10 seconds later

            # Mock simulation progress
            tracker.act(timestep=1000)

            # Check that performance metrics are calculated
            assert tracker.current_ns_per_day >= 0.0
            assert isinstance(tracker.current_eta, str)

    def test_logging_methods(self):
        """Test logging methods return correct types."""
        tracker = PerformanceTracker(simulation=self.sim, runtime_ps=self.runtime_ps)

        # Test that logging methods return strings (they are properties due to HOOMD's decorator)
        ns_per_day = tracker.ns_per_day
        eta = tracker.eta_remaining

        assert isinstance(ns_per_day, str)
        assert isinstance(eta, str)

    def test_performance_tracker_with_different_periods(self):
        """Test PerformanceTracker with different output periods."""
        # Test with step-based period
        tracker_steps = PerformanceTracker(
            simulation=self.sim, runtime_ps=self.runtime_ps, output_period_steps=500
        )
        assert tracker_steps.output_period_steps == 500
        assert not tracker_steps.use_time_based_output

        # Test with time-based period
        tracker_time = PerformanceTracker(
            simulation=self.sim, runtime_ps=self.runtime_ps, output_period_ps=1.0
        )
        assert tracker_time.output_period_ps == 1.0
        assert tracker_time.use_time_based_output

    def test_ns_per_day_formatting(self):
        """Test ns/day formatting for different values."""
        tracker = PerformanceTracker(simulation=self.sim, runtime_ps=self.runtime_ps)

        # Test small values (should have 4 decimal places)
        tracker.current_ns_per_day = 0.005
        result = tracker.ns_per_day  # Property access, not method call
        assert '.0050' in result

        # Test larger values (should have 2 decimal places)
        tracker.current_ns_per_day = 1.234
        result = tracker.ns_per_day  # Property access, not method call
        assert '1.23' in result

    def test_eta_calculation_with_mock_time(self):
        """Test ETA calculation with mocked time progression."""
        tracker = PerformanceTracker(simulation=self.sim, runtime_ps=self.runtime_ps)

        with patch('time.time') as mock_time:
            # Set initial time
            mock_time.return_value = 1000.0
            tracker.start_time = 1000.0

            # Simulate 10 seconds passing, 1 ps simulated
            mock_time.return_value = 1010.0

            # Mock dt and timestep for realistic calculation
            self.integrator.dt = PhysicalConstants.ps_to_atomic_units(0.001)

            tracker.act(timestep=1000)

            # Check that ETA is calculated (should be a valid time string)
            eta = tracker.eta_remaining  # Property access, not method call
            assert isinstance(eta, str)
            assert ':' in eta  # Time format should contain colons

    def test_performance_tracker_integration_with_hoomd_logging(self):
        """Test that PerformanceTracker integrates with HOOMD logging."""
        tracker = PerformanceTracker(simulation=self.sim, runtime_ps=self.runtime_ps)

        # Check that logging properties are accessible
        assert hasattr(tracker, 'ns_per_day')
        assert hasattr(tracker, 'eta_remaining')

        # Check that they return strings
        assert isinstance(tracker.ns_per_day, str)
        assert isinstance(tracker.eta_remaining, str)

    def test_performance_tracker_with_zero_runtime(self):
        """Test PerformanceTracker behavior with zero runtime."""
        tracker = PerformanceTracker(simulation=self.sim, runtime_ps=0.0)

        # Should not crash with zero runtime
        tracker.act(timestep=1)

        # Should return valid strings even with zero runtime
        assert isinstance(tracker.ns_per_day, str)  # Property access
        assert isinstance(tracker.eta_remaining, str)  # Property access


@pytest.mark.skipif(not PLUGIN_AVAILABLE, reason='Plugin not available')
class TestPerformanceTrackerEdgeCases:
    """Test edge cases and error conditions for PerformanceTracker."""

    def test_performance_tracker_with_invalid_simulation(self):
        """Test PerformanceTracker with invalid simulation object."""
        # This should work during initialization but fail during act()
        tracker = PerformanceTracker(simulation=None, runtime_ps=100.0)

        # Should raise an error when trying to act
        with pytest.raises(AttributeError):
            tracker.act(timestep=1)

    def test_performance_tracker_with_negative_runtime(self):
        """Test PerformanceTracker with negative runtime."""
        device = hoomd.device.CPU()
        sim = hoomd.Simulation(device=device, seed=42)

        # Set up minimal integrator to prevent None reference
        snapshot = hoomd.Snapshot(device.communicator)
        snapshot.configuration.box = [10, 10, 10, 0, 0, 0]
        snapshot.particles.N = 1
        snapshot.particles.types = ['A']
        snapshot.particles.position[0] = [0, 0, 0]
        snapshot.particles.velocity[0] = [0, 0, 0]
        snapshot.particles.mass[0] = 1.0

        sim.create_state_from_snapshot(snapshot)

        integrator = hoomd.md.Integrator(dt=0.001)
        integrator.methods.append(
            hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All())
        )
        sim.operations.integrator = integrator

        # Should not crash with negative runtime
        tracker = PerformanceTracker(simulation=sim, runtime_ps=-10.0)

        # Should handle negative runtime gracefully
        tracker.act(timestep=1)
        assert isinstance(tracker.ns_per_day, str)
        assert isinstance(tracker.eta_remaining, str)


if __name__ == '__main__':
    pytest.main([__file__])
