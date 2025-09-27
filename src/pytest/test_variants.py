"""
Tests for the variants module.

This module tests all custom variant classes for time-varying parameters
in cavity molecular dynamics simulations.
"""

import pytest
import numpy as np
import hoomd
from unittest.mock import Mock, MagicMock
import tempfile
from pathlib import Path

from hoomd.cavitymd.variants import StepVariant, ConstantVariant, PeriodicVariant
from hoomd.cavitymd.analysis import ElapsedTimeTracker


class TestPeriodicVariant:
    """Test PeriodicVariant class functionality."""
    
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
    
    def test_initialization_with_all_parameters(self, mock_time_tracker):
        """Test initialization with all optional parameters."""
        variant = PeriodicVariant(
            amplitude=0.002,
            time_tracker=mock_time_tracker,
            period_ps=2.5,
            phase_offset=np.pi/4,
            start_time_ps=1.0
        )
        
        assert variant.amplitude == 0.002
        assert variant.period_ps == 2.5
        assert variant.phase_offset == np.pi/4
        assert variant.start_time_ps == 1.0
        
        expected_frequency = 1.0 / (2.5 * 1e-12)  # Hz
        assert abs(variant.frequency_hz - expected_frequency) < expected_frequency * 1e-10
    
    def test_initialization_validation(self, mock_time_tracker):
        """Test parameter validation during initialization."""
        # Test negative period
        with pytest.raises(ValueError, match="Period must be positive"):
            PeriodicVariant(
                amplitude=0.001,
                time_tracker=mock_time_tracker,
                period_ps=-1.0
            )
        
        # Test zero period
        with pytest.raises(ValueError, match="Period must be positive"):
            PeriodicVariant(
                amplitude=0.001,
                time_tracker=mock_time_tracker,
                period_ps=0.0
            )
        
        # Test negative amplitude
        with pytest.raises(ValueError, match="Amplitude must be non-negative"):
            PeriodicVariant(
                amplitude=-0.001,
                time_tracker=mock_time_tracker,
                period_ps=1.0
            )
        
        # Test negative start time
        with pytest.raises(ValueError, match="Start time must be non-negative"):
            PeriodicVariant(
                amplitude=0.001,
                time_tracker=mock_time_tracker,
                period_ps=1.0,
                start_time_ps=-1.0
            )
    
    def test_call_before_start_time(self, mock_time_tracker):
        """Test variant evaluation before start time."""
        variant = PeriodicVariant(
            amplitude=0.001,
            time_tracker=mock_time_tracker,
            period_ps=1.0,
            start_time_ps=2.0
        )
        
        # Set time before start
        mock_time_tracker.elapsed_time = 1.5
        
        value = variant(1000)  # timestep argument is ignored
        assert value == 0.0
        assert variant.current_value == 0.0
        assert not variant.oscillation_started
    
    def test_call_basic_oscillation(self, mock_time_tracker):
        """Test basic sinusoidal oscillation."""
        variant = PeriodicVariant(
            amplitude=0.001,
            time_tracker=mock_time_tracker,
            period_ps=1.0,
            phase_offset=0.0,
            start_time_ps=0.0
        )
        
        # Test at t=0 (should be 0 due to sin(0))
        mock_time_tracker.elapsed_time = 0.0
        value = variant(0)
        assert abs(value) < 1e-15  # Should be exactly 0
        assert variant.oscillation_started
        
        # Test at t=T/4 (should be maximum)
        mock_time_tracker.elapsed_time = 0.25
        value = variant(1)
        expected = 0.001 * np.sin(2*np.pi * 0.25)  # sin(π/2) = 1
        assert abs(value - 0.001) < 1e-10
        assert abs(value - expected) < 1e-15
        
        # Test at t=T/2 (should be 0)
        mock_time_tracker.elapsed_time = 0.5
        value = variant(2)
        expected = 0.001 * np.sin(2*np.pi * 0.5)  # sin(π) = 0
        assert abs(value) < 1e-15
        
        # Test at t=3T/4 (should be minimum)
        mock_time_tracker.elapsed_time = 0.75
        value = variant(3)
        expected = 0.001 * np.sin(2*np.pi * 0.75)  # sin(3π/2) = -1
        assert abs(value - (-0.001)) < 1e-10
        
        # Test at t=T (should be 0 again)
        mock_time_tracker.elapsed_time = 1.0
        value = variant(4)
        expected = 0.001 * np.sin(2*np.pi * 1.0)  # sin(2π) = 0
        assert abs(value) < 1e-15
    
    def test_call_with_phase_offset(self, mock_time_tracker):
        """Test oscillation with phase offset."""
        variant = PeriodicVariant(
            amplitude=0.001,
            time_tracker=mock_time_tracker,
            period_ps=1.0,
            phase_offset=np.pi/2,  # Start at maximum
            start_time_ps=0.0
        )
        
        # Test at t=0 (should be maximum due to phase offset)
        mock_time_tracker.elapsed_time = 0.0
        value = variant(0)
        expected = 0.001 * np.sin(np.pi/2)  # sin(π/2) = 1
        assert abs(value - 0.001) < 1e-10
        
        # Test at t=T/4 (should be 0)
        mock_time_tracker.elapsed_time = 0.25
        value = variant(1)
        expected = 0.001 * np.sin(2*np.pi * 0.25 + np.pi/2)  # sin(π) = 0
        assert abs(value) < 1e-15
    
    def test_call_with_start_time(self, mock_time_tracker):
        """Test oscillation starting at non-zero time."""
        variant = PeriodicVariant(
            amplitude=0.001,
            time_tracker=mock_time_tracker,
            period_ps=2.0,
            phase_offset=0.0,
            start_time_ps=1.0
        )
        
        # Before start time
        mock_time_tracker.elapsed_time = 0.5
        value = variant(0)
        assert value == 0.0
        assert not variant.oscillation_started
        
        # At start time
        mock_time_tracker.elapsed_time = 1.0
        value = variant(1)
        assert abs(value) < 1e-15  # sin(0) = 0
        assert variant.oscillation_started
        
        # After start time (T/4 into oscillation)
        mock_time_tracker.elapsed_time = 1.5  # 0.5 ps since start
        value = variant(2)
        # Phase = 2π * 0.5 / 2.0 = π/2
        expected = 0.001 * np.sin(np.pi/2)  # sin(π/2) = 1
        assert abs(value - 0.001) < 1e-10
    
    def test_multiple_periods(self, mock_time_tracker):
        """Test oscillation across multiple periods."""
        variant = PeriodicVariant(
            amplitude=0.5,
            time_tracker=mock_time_tracker,
            period_ps=1.0,
            phase_offset=0.0,
            start_time_ps=0.0
        )
        
        # Test across several periods
        times = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25]
        expected_values = [
            0.0,    # t=0: sin(0) = 0
            0.5,    # t=T/4: sin(π/2) = 1
            0.0,    # t=T/2: sin(π) = 0  
            -0.5,   # t=3T/4: sin(3π/2) = -1
            0.0,    # t=T: sin(2π) = 0
            0.5,    # t=5T/4: sin(5π/2) = 1
            0.0,    # t=3T/2: sin(3π) = 0
            -0.5,   # t=7T/4: sin(7π/2) = -1
            0.0,    # t=2T: sin(4π) = 0
            0.5     # t=9T/4: sin(9π/2) = 1
        ]
        
        for i, (time_val, expected_val) in enumerate(zip(times, expected_values)):
            mock_time_tracker.elapsed_time = time_val
            value = variant(i)
            assert abs(value - expected_val) < 1e-10, f"Failed at t={time_val}"
    
    def test_min_max_values(self, mock_time_tracker):
        """Test _min and _max methods."""
        variant = PeriodicVariant(
            amplitude=0.003,
            time_tracker=mock_time_tracker,
            period_ps=1.0
        )
        
        assert variant._min() == -0.003
        assert variant._max() == 0.003
    
    def test_frequency_conversion_methods(self, mock_time_tracker):
        """Test frequency conversion utility methods."""
        # 1 ps period = 1 THz frequency
        variant = PeriodicVariant(
            amplitude=0.001,
            time_tracker=mock_time_tracker,
            period_ps=1.0
        )
        
        # Check THz conversion
        freq_thz = variant.get_frequency_thz()
        assert abs(freq_thz - 1.0) < 1e-10
        
        # Check cm⁻¹ conversion
        freq_cm = variant.get_frequency_cm_minus1()
        c_light = 2.99792458e10  # cm/s
        expected_cm = 1e12 / c_light  # Convert 1 THz to cm⁻¹
        assert abs(freq_cm - expected_cm) < expected_cm * 1e-10
        
        # Test different period
        variant2 = PeriodicVariant(
            amplitude=0.001,
            time_tracker=mock_time_tracker,
            period_ps=0.5  # 2 THz
        )
        
        freq_thz2 = variant2.get_frequency_thz()
        assert abs(freq_thz2 - 2.0) < 1e-10
    
    def test_zero_amplitude(self, mock_time_tracker):
        """Test oscillation with zero amplitude."""
        variant = PeriodicVariant(
            amplitude=0.0,
            time_tracker=mock_time_tracker,
            period_ps=1.0
        )
        
        # Should always return 0 regardless of time
        for time_val in [0.0, 0.25, 0.5, 0.75, 1.0]:
            mock_time_tracker.elapsed_time = time_val
            value = variant(0)
            assert value == 0.0
        
        assert variant._min() == 0.0
        assert variant._max() == 0.0
    
    def test_very_small_period(self, mock_time_tracker):
        """Test oscillation with very small period (high frequency)."""
        variant = PeriodicVariant(
            amplitude=0.001,
            time_tracker=mock_time_tracker,
            period_ps=0.01  # 100 THz
        )
        
        # Test at various fractions of the small period
        mock_time_tracker.elapsed_time = 0.0025  # T/4
        value = variant(0)
        expected = 0.001 * np.sin(np.pi/2)
        assert abs(value - 0.001) < 1e-10
        
        freq_thz = variant.get_frequency_thz()
        assert abs(freq_thz - 100.0) < 1e-8
    
    def test_consistency_across_calls(self, mock_time_tracker):
        """Test that repeated calls at same time give same result."""
        variant = PeriodicVariant(
            amplitude=0.001,
            time_tracker=mock_time_tracker,
            period_ps=1.0
        )
        
        mock_time_tracker.elapsed_time = 0.3
        
        # Call multiple times at same timestep
        value1 = variant(100)
        value2 = variant(100)
        value3 = variant(100)
        
        assert value1 == value2 == value3
        assert variant.current_value == value1


class TestStepVariant:
    """Test StepVariant class functionality."""
    
    @pytest.fixture
    def mock_time_tracker(self):
        """Create a mock time tracker for testing."""
        tracker = Mock(spec=ElapsedTimeTracker)
        tracker.elapsed_time = 0.0
        return tracker
    
    def test_step_basic_functionality(self, mock_time_tracker):
        """Test basic step function without decay."""
        variant = StepVariant(
            target_value=0.002,
            switch_time_ps=1.0,
            time_tracker=mock_time_tracker
        )
        
        # Before switch
        mock_time_tracker.elapsed_time = 0.5
        value = variant(100)
        assert value == 0.0
        assert not variant.has_switched
        
        # At switch time
        mock_time_tracker.elapsed_time = 1.0
        value = variant(200)
        assert value == 0.002
        assert variant.has_switched
        
        # After switch
        mock_time_tracker.elapsed_time = 2.0
        value = variant(300)
        assert value == 0.002
    
    def test_step_with_decay(self, mock_time_tracker):
        """Test step function with exponential decay."""
        variant = StepVariant(
            target_value=0.002,
            switch_time_ps=1.0,
            time_tracker=mock_time_tracker,
            decay_time_constant_ps=1.0
        )
        
        # Before switch
        mock_time_tracker.elapsed_time = 0.5
        value = variant(100)
        assert value == 0.0
        
        # At switch time
        mock_time_tracker.elapsed_time = 1.0
        value = variant(200)
        assert value == 0.002
        
        # After switch (one time constant)
        mock_time_tracker.elapsed_time = 2.0
        value = variant(300)
        expected = 0.002 * np.exp(-1.0)  # exp(-1) ≈ 0.368
        assert abs(value - expected) < 1e-10


class TestConstantVariant:
    """Test ConstantVariant class functionality."""
    
    def test_constant_functionality(self):
        """Test constant variant returns same value always."""
        variant = ConstantVariant(0.0015)
        
        # Should return same value regardless of timestep
        for timestep in [0, 100, 1000, 10000]:
            value = variant(timestep)
            assert value == 0.0015
        
        assert variant._min() == 0.0015
        assert variant._max() == 0.0015
    
    def test_constant_zero(self):
        """Test constant variant with zero value."""
        variant = ConstantVariant(0.0)
        
        value = variant(500)
        assert value == 0.0
        assert variant._min() == 0.0
        assert variant._max() == 0.0
    
    def test_constant_negative(self):
        """Test constant variant with negative value."""
        variant = ConstantVariant(-0.001)
        
        value = variant(500)
        assert value == -0.001
        assert variant._min() == -0.001
        assert variant._max() == -0.001


class TestVariantIntegration:
    """Test integration of variants with HOOMD system."""
    
    @pytest.fixture
    def mock_time_tracker(self):
        """Create a mock time tracker for testing."""
        tracker = Mock(spec=ElapsedTimeTracker)
        tracker.elapsed_time = 0.0
        return tracker
    
    def test_variant_inheritance(self, mock_time_tracker):
        """Test that variants properly inherit from hoomd.variant.Variant."""
        periodic = PeriodicVariant(0.001, mock_time_tracker)
        step = StepVariant(0.001, 1.0, mock_time_tracker)
        constant = ConstantVariant(0.001)
        
        # All should inherit from hoomd.variant.Variant
        assert isinstance(periodic, hoomd.variant.Variant)
        assert isinstance(step, hoomd.variant.Variant)
        assert isinstance(constant, hoomd.variant.Variant)
        
        # All should be callable
        assert callable(periodic)
        assert callable(step)
        assert callable(constant)
        
        # All should have required methods
        for variant in [periodic, step, constant]:
            assert hasattr(variant, '_min')
            assert hasattr(variant, '_max')
            assert callable(variant._min)
            assert callable(variant._max)
    
    def test_variant_serialization_compatibility(self, mock_time_tracker):
        """Test that variants can be created and called like HOOMD variants."""
        # This test ensures our variants work with HOOMD's internal systems
        periodic = PeriodicVariant(0.001, mock_time_tracker, period_ps=0.5)
        
        # Test that we can call with integer timesteps (as HOOMD does)
        timesteps = [0, 100, 1000, 10000]
        values = []
        
        for i, ts in enumerate(timesteps):
            mock_time_tracker.elapsed_time = i * 0.1  # Advance time
            value = periodic(ts)
            values.append(value)
            assert isinstance(value, float)
        
        # Values should change over time (not constant)
        assert not all(v == values[0] for v in values)
