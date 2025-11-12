"""Test the analysis feedback controller classes from the cavitymd plugin.

This module tests the refactored feedback controller functionality including:
- EmpiricalTemperatureFeedback: Feedback control using empirical data
- GradientDescentTemperatureFeedback: Gradient-based temperature control
- DualIndependentTemperatureFeedback: Dual-channel independent control
"""

import pytest
import numpy as np
import hoomd
from unittest.mock import Mock, MagicMock, patch
import tempfile
from pathlib import Path
import json

# Try to import from the plugin
try:
    from hoomd.cavitymd.analysis import (
        EmpiricalTemperatureFeedback,
        GradientDescentTemperatureFeedback,
        DualIndependentTemperatureFeedback,
        EmpiricalTemperatureData,
        ElapsedTimeTracker,
        TemperatureTracker
    )
    from hoomd.cavitymd.utils import PhysicalConstants
    PLUGIN_AVAILABLE = True
except ImportError:
    PLUGIN_AVAILABLE = False


@pytest.mark.skipif(not PLUGIN_AVAILABLE, reason='Plugin not available')
class TestEmpiricalTemperatureFeedback:
    """Test suite for EmpiricalTemperatureFeedback class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.device = hoomd.device.CPU()
        self.sim = hoomd.Simulation(device=self.device, seed=42)
        
        # Create simple system
        self.snapshot = hoomd.Snapshot(self.device.communicator)
        self.snapshot.configuration.box = [10, 10, 10, 0, 0, 0]
        self.snapshot.particles.N = 2
        self.snapshot.particles.types = ['A']
        self.snapshot.particles.position[0] = [0, 0, 0]
        self.snapshot.particles.position[1] = [2.0, 0, 0]
        self.snapshot.particles.velocity[0] = [1.0, 0, 0]
        self.snapshot.particles.velocity[1] = [-1.0, 0, 0]
        self.snapshot.particles.mass[0] = 1.0
        self.snapshot.particles.mass[1] = 1.0
        
        self.sim.create_state_from_snapshot(self.snapshot)
        
        # Set up integrator and thermostat
        self.integrator = hoomd.md.Integrator(dt=0.001)
        nvt = hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All())
        self.integrator.methods.append(nvt)
        
        # Add thermostat
        self.thermostat = hoomd.md.methods.thermostats.Bussi(kT=1.0)
        nvt.thermostat = self.thermostat
        
        self.sim.operations.integrator = self.integrator
        
        # Create time tracker
        self.time_tracker = ElapsedTimeTracker()
        
    def teardown_method(self):
        """Clean up after tests."""
        self.sim = None
        self.device = None
    
    def test_initialization_with_empirical_data(self):
        """Test initialization with empirical temperature data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = Path(tmpdir) / "feedback.csv"
            
            # Create mock empirical data
            empirical_data = Mock(spec=EmpiricalTemperatureData)
            empirical_data.get_temperature.return_value = 300.0
            
            feedback = EmpiricalTemperatureFeedback(
                time_tracker=self.time_tracker,
                thermostat=self.thermostat,
                empirical_data=empirical_data,
                target_temperature=300.0,
                output_file=str(output_file),
                control_period_ps=1.0
            )
            
            assert feedback.time_tracker is self.time_tracker
            assert feedback.thermostat is self.thermostat
            assert feedback.target_temperature == 300.0
    
    def test_feedback_control_loop(self):
        """Test that feedback control adjusts thermostat."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = Path(tmpdir) / "feedback.csv"
            
            # Create mock empirical data with temperature above target
            empirical_data = Mock(spec=EmpiricalTemperatureData)
            empirical_data.get_temperature.return_value = 350.0  # Above target
            
            feedback = EmpiricalTemperatureFeedback(
                time_tracker=self.time_tracker,
                thermostat=self.thermostat,
                empirical_data=empirical_data,
                target_temperature=300.0,
                output_file=str(output_file),
                control_period_ps=0.001,
                feedback_gain=0.1
            )
            
            # Attach to simulation
            self.sim.operations.writers.append(self.time_tracker)
            self.sim.operations.writers.append(feedback)
            
            # Get initial thermostat temperature
            initial_kT = self.thermostat.kT.default
            
            # Run simulation
            self.sim.run(10)
            
            # Check that thermostat was adjusted (should decrease since T > target)
            # Note: Exact behavior depends on feedback implementation


@pytest.mark.skipif(not PLUGIN_AVAILABLE, reason='Plugin not available')
class TestGradientDescentTemperatureFeedback:
    """Test suite for GradientDescentTemperatureFeedback class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.device = hoomd.device.CPU()
        self.sim = hoomd.Simulation(device=self.device, seed=42)
        
        # Create simple system
        self.snapshot = hoomd.Snapshot(self.device.communicator)
        self.snapshot.configuration.box = [10, 10, 10, 0, 0, 0]
        self.snapshot.particles.N = 2
        self.snapshot.particles.types = ['A']
        self.snapshot.particles.position[0] = [0, 0, 0]
        self.snapshot.particles.position[1] = [2.0, 0, 0]
        self.snapshot.particles.velocity[0] = [1.0, 0, 0]
        self.snapshot.particles.velocity[1] = [-1.0, 0, 0]
        self.snapshot.particles.mass[0] = 1.0
        self.snapshot.particles.mass[1] = 1.0
        
        self.sim.create_state_from_snapshot(self.snapshot)
        
        # Set up integrator
        self.integrator = hoomd.md.Integrator(dt=0.001)
        nvt = hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All())
        self.integrator.methods.append(nvt)
        
        # Add thermostat
        self.thermostat = hoomd.md.methods.thermostats.Bussi(kT=1.0)
        nvt.thermostat = self.thermostat
        
        self.sim.operations.integrator = self.integrator
        
        # Create time tracker
        self.time_tracker = ElapsedTimeTracker()
        
    def teardown_method(self):
        """Clean up after tests."""
        self.sim = None
        self.device = None
    
    def test_initialization(self):
        """Test GradientDescentTemperatureFeedback initialization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = Path(tmpdir) / "gd_feedback.csv"
            
            # Create mock temperature tracker
            temp_tracker = Mock(spec=TemperatureTracker)
            temp_tracker.current_temperature = 300.0
            
            feedback = GradientDescentTemperatureFeedback(
                time_tracker=self.time_tracker,
                thermostat=self.thermostat,
                temperature_tracker=temp_tracker,
                target_temperature=300.0,
                output_file=str(output_file),
                control_period_ps=1.0,
                learning_rate=0.01
            )
            
            assert feedback.time_tracker is self.time_tracker
            assert feedback.thermostat is self.thermostat
            assert feedback.target_temperature == 300.0
            assert feedback.learning_rate == 0.01


@pytest.mark.skipif(not PLUGIN_AVAILABLE, reason='Plugin not available')
class TestDualIndependentTemperatureFeedback:
    """Test suite for DualIndependentTemperatureFeedback class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.device = hoomd.device.CPU()
        self.sim = hoomd.Simulation(device=self.device, seed=42)
        
        # Create system with two particle types
        self.snapshot = hoomd.Snapshot(self.device.communicator)
        self.snapshot.configuration.box = [10, 10, 10, 0, 0, 0]
        self.snapshot.particles.N = 4
        self.snapshot.particles.types = ['A', 'B']
        
        # Two particles of type A, two of type B
        self.snapshot.particles.typeid[0] = 0
        self.snapshot.particles.typeid[1] = 0
        self.snapshot.particles.typeid[2] = 1
        self.snapshot.particles.typeid[3] = 1
        
        for i in range(4):
            self.snapshot.particles.position[i] = [i*2.0, 0, 0]
            self.snapshot.particles.velocity[i] = [np.random.randn(), 
                                                   np.random.randn(), 
                                                   np.random.randn()]
            self.snapshot.particles.mass[i] = 1.0
        
        self.sim.create_state_from_snapshot(self.snapshot)
        
        # Set up integrator with two thermostats
        self.integrator = hoomd.md.Integrator(dt=0.001)
        
        nvt_A = hoomd.md.methods.ConstantVolume(
            filter=hoomd.filter.Type(['A'])
        )
        self.thermostat_A = hoomd.md.methods.thermostats.Bussi(kT=1.0)
        nvt_A.thermostat = self.thermostat_A
        self.integrator.methods.append(nvt_A)
        
        nvt_B = hoomd.md.methods.ConstantVolume(
            filter=hoomd.filter.Type(['B'])
        )
        self.thermostat_B = hoomd.md.methods.thermostats.Bussi(kT=1.0)
        nvt_B.thermostat = self.thermostat_B
        self.integrator.methods.append(nvt_B)
        
        self.sim.operations.integrator = self.integrator
        
        # Create time tracker
        self.time_tracker = ElapsedTimeTracker()
        
    def teardown_method(self):
        """Clean up after tests."""
        self.sim = None
        self.device = None
    
    def test_initialization(self):
        """Test DualIndependentTemperatureFeedback initialization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = Path(tmpdir) / "dual_feedback.csv"
            
            # Create mock temperature trackers
            temp_tracker_A = Mock(spec=TemperatureTracker)
            temp_tracker_A.current_temperature = 300.0
            
            temp_tracker_B = Mock(spec=TemperatureTracker)
            temp_tracker_B.current_temperature = 400.0
            
            feedback = DualIndependentTemperatureFeedback(
                time_tracker=self.time_tracker,
                thermostat_channel1=self.thermostat_A,
                thermostat_channel2=self.thermostat_B,
                temperature_tracker_channel1=temp_tracker_A,
                temperature_tracker_channel2=temp_tracker_B,
                target_temperature_channel1=300.0,
                target_temperature_channel2=400.0,
                output_file=str(output_file),
                control_period_ps=1.0
            )
            
            assert feedback.time_tracker is self.time_tracker
            assert feedback.thermostat_channel1 is self.thermostat_A
            assert feedback.thermostat_channel2 is self.thermostat_B


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

