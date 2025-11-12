"""Test the AdaptiveTimestepUpdater class from the cavitymd plugin.

This module tests the adaptive timestep functionality including:
- Adaptive timestep adjustment based on error tolerance
- Shock dampening for sudden coupling changes
- Integration with HOOMD simulation framework
"""

import pytest
import numpy as np
import hoomd
from unittest.mock import Mock, MagicMock, patch

# Try to import from the plugin
try:
    from hoomd.cavitymd.simulation import AdaptiveTimestepUpdater
    from hoomd.cavitymd.analysis import ElapsedTimeTracker
    from hoomd.cavitymd.utils import PhysicalConstants
    PLUGIN_AVAILABLE = True
except ImportError:
    PLUGIN_AVAILABLE = False


@pytest.mark.skipif(not PLUGIN_AVAILABLE, reason='Plugin not available')
class TestAdaptiveTimestepUpdater:
    """Test suite for AdaptiveTimestepUpdater class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        # Create a minimal simulation for testing
        self.device = hoomd.device.CPU()
        self.sim = hoomd.Simulation(device=self.device, seed=42)
        
        # Create a simple system (two particles for testing)
        self.snapshot = hoomd.Snapshot(self.device.communicator)
        self.snapshot.configuration.box = [10, 10, 10, 0, 0, 0]
        self.snapshot.particles.N = 2
        self.snapshot.particles.types = ['A']
        self.snapshot.particles.position[0] = [0, 0, 0]
        self.snapshot.particles.position[1] = [1.5, 0, 0]
        self.snapshot.particles.velocity[0] = [1.0, 0, 0]
        self.snapshot.particles.velocity[1] = [-1.0, 0, 0]
        self.snapshot.particles.mass[0] = 1.0
        self.snapshot.particles.mass[1] = 1.0
        
        self.sim.create_state_from_snapshot(self.snapshot)
        
        # Set up a minimal integrator
        self.integrator = hoomd.md.Integrator(dt=0.001)
        self.integrator.methods.append(
            hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All())
        )
        
        # Add a simple LJ force
        lj = hoomd.md.pair.LJ(nlist=hoomd.md.nlist.Cell(buffer=0.4))
        lj.params[('A', 'A')] = {'epsilon': 1.0, 'sigma': 1.0}
        lj.r_cut[('A', 'A')] = 2.5
        self.integrator.forces.append(lj)
        
        self.sim.operations.integrator = self.integrator
        
        # Create time tracker
        self.time_tracker = ElapsedTimeTracker()
        
    def teardown_method(self):
        """Clean up after tests."""
        self.sim = None
        self.device = None
    
    def test_initialization_basic(self):
        """Test basic initialization of AdaptiveTimestepUpdater."""
        updater = AdaptiveTimestepUpdater(
            state=self.sim.state,
            integrator=self.integrator,
            error_tolerance=1e-6,
            time_constant_ps=50.0,
            time_tracker=self.time_tracker
        )
        
        assert updater.state is self.sim.state
        assert updater.integrator is self.integrator
        assert updater.target_error_tolerance == 1e-6
        assert updater.time_constant_ps == 50.0
        assert updater.time_tracker is self.time_tracker
    
    def test_initialization_with_shock_dampening(self):
        """Test initialization with shock dampening parameters."""
        updater = AdaptiveTimestepUpdater(
            state=self.sim.state,
            integrator=self.integrator,
            error_tolerance=1e-6,
            shock_dampening_enabled=True,
            shock_dampening_factor=1e-3,
            switch_time_ps=100.0,
            time_tracker=self.time_tracker
        )
        
        assert updater.shock_dampening_enabled == True
        assert updater.shock_dampening_factor == 1e-3
        assert updater.switch_time_ps == 100.0
        assert updater.shock_error_tolerance == 1e-6 * 1e-3
    
    def test_dynamic_coupling_detection(self):
        """Test dynamic coupling change detection."""
        updater = AdaptiveTimestepUpdater(
            state=self.sim.state,
            integrator=self.integrator,
            error_tolerance=1e-6,
            dynamic_coupling_detection=True,
            coupling_change_threshold=1e-5,
            time_tracker=self.time_tracker
        )
        
        assert updater.dynamic_coupling_detection == True
        assert updater.coupling_change_threshold == 1e-5
        assert updater.last_coupling_value is None
        assert updater.coupling_recovery_active == False
    
    def test_timestep_adjustment(self):
        """Test that timestep is adjusted during simulation."""
        updater = AdaptiveTimestepUpdater(
            state=self.sim.state,
            integrator=self.integrator,
            error_tolerance=1e-6,
            time_tracker=self.time_tracker
        )
        
        # Attach updater to simulation
        self.sim.operations.writers.append(updater)
        self.sim.operations.writers.append(self.time_tracker)
        
        # Get initial timestep
        initial_dt = self.integrator.dt
        
        # Run simulation
        self.sim.run(1000)
        
        # Check that timestep was adjusted (may be same or different)
        final_dt = self.integrator.dt
        
        # At minimum, check that dt is positive and reasonable
        assert final_dt > 0
        assert final_dt < 1.0  # Should be less than 1 a.u. (reasonable upper bound)
    
    def test_error_tolerance_ramping(self):
        """Test error tolerance ramping over time."""
        updater = AdaptiveTimestepUpdater(
            state=self.sim.state,
            integrator=self.integrator,
            error_tolerance=1e-6,
            initial_fraction=0.01,
            time_constant_ps=10.0,
            adaptiveerror=True,
            time_tracker=self.time_tracker
        )
        
        # Initial error tolerance should be reduced
        assert updater.current_error_tolerance == 1e-6 * 0.01
        
        # Attach to simulation
        self.sim.operations.writers.append(updater)
        self.sim.operations.writers.append(self.time_tracker)
        
        # Run simulation for some time
        self.sim.run(100)
        
        # Error tolerance should have increased (ramped up)
        assert updater.current_error_tolerance > 1e-6 * 0.01
        assert updater.current_error_tolerance <= 1e-6
    
    def test_conservative_timestep_parameters(self):
        """Test conservative timestep change parameters."""
        updater = AdaptiveTimestepUpdater(
            state=self.sim.state,
            integrator=self.integrator,
            error_tolerance=1e-6,
            timestep_change_threshold=0.1,
            max_timestep_change_factor=1.5,
            time_tracker=self.time_tracker
        )
        
        assert updater.timestep_change_threshold == 0.1
        assert updater.max_timestep_change_factor == 1.5
        assert updater.min_update_interval == 1000
    
    def test_switch_detection(self):
        """Test switch time detection."""
        updater = AdaptiveTimestepUpdater(
            state=self.sim.state,
            integrator=self.integrator,
            error_tolerance=1e-6,
            switch_time_ps=10.0,
            shock_dampening_enabled=True,
            shock_dampening_factor=1e-3,
            time_tracker=self.time_tracker
        )
        
        assert updater.switch_time_ps == 10.0
        assert updater.switch_detected == False
        
        # Attach to simulation
        self.sim.operations.writers.append(updater)
        self.sim.operations.writers.append(self.time_tracker)
        
        # Run until past switch time
        # Note: Exact behavior depends on timestep and update interval
        # Just verify no crashes occur
        self.sim.run(100)
    
    def test_elapsed_time_tracking(self):
        """Test elapsed time tracking property."""
        updater = AdaptiveTimestepUpdater(
            state=self.sim.state,
            integrator=self.integrator,
            error_tolerance=1e-6,
            time_tracker=self.time_tracker
        )
        
        # Attach to simulation
        self.sim.operations.writers.append(updater)
        self.sim.operations.writers.append(self.time_tracker)
        
        # Run simulation
        self.sim.run(100)
        
        # Check elapsed time is tracked
        elapsed = updater.elapsed_time_ps()
        assert elapsed >= 0
        assert isinstance(elapsed, (int, float))


@pytest.mark.skipif(not PLUGIN_AVAILABLE, reason='Plugin not available')
class TestAdaptiveTimestepIntegration:
    """Integration tests for AdaptiveTimestepUpdater with various force scenarios."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.device = hoomd.device.CPU()
        self.sim = hoomd.Simulation(device=self.device, seed=42)
        
        # Create a system with varying force magnitudes
        self.snapshot = hoomd.Snapshot(self.device.communicator)
        self.snapshot.configuration.box = [10, 10, 10, 0, 0, 0]
        self.snapshot.particles.N = 5
        self.snapshot.particles.types = ['A']
        
        # Position particles at varying separations
        for i in range(5):
            self.snapshot.particles.position[i] = [i*1.1, 0, 0]
            self.snapshot.particles.velocity[i] = [0, 0, 0]
            self.snapshot.particles.mass[i] = 1.0
        
        self.sim.create_state_from_snapshot(self.snapshot)
        
        # Set up integrator
        self.integrator = hoomd.md.Integrator(dt=0.001)
        self.integrator.methods.append(
            hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All())
        )
        
        # Add LJ force with strong interaction
        lj = hoomd.md.pair.LJ(nlist=hoomd.md.nlist.Cell(buffer=0.4))
        lj.params[('A', 'A')] = {'epsilon': 10.0, 'sigma': 1.0}
        lj.r_cut[('A', 'A')] = 2.5
        self.integrator.forces.append(lj)
        
        self.sim.operations.integrator = self.integrator
        
        # Create time tracker
        self.time_tracker = ElapsedTimeTracker()
        
    def teardown_method(self):
        """Clean up after tests."""
        self.sim = None
        self.device = None
    
    def test_adaptive_response_to_high_forces(self):
        """Test that timestep adapts to high force magnitudes."""
        updater = AdaptiveTimestepUpdater(
            state=self.sim.state,
            integrator=self.integrator,
            error_tolerance=1e-5,
            time_tracker=self.time_tracker
        )
        
        # Attach to simulation
        self.sim.operations.writers.append(updater)
        self.sim.operations.writers.append(self.time_tracker)
        
        # Initial timestep
        initial_dt = self.integrator.dt
        
        # Run simulation - should adapt to high forces
        self.sim.run(1000)
        
        # Final timestep should be adjusted
        final_dt = self.integrator.dt
        
        # With high forces, timestep should typically decrease
        # (though exact behavior depends on force magnitude)
        assert final_dt > 0
        assert final_dt < 1.0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

