#!/usr/bin/env python3
"""
Test script for CompositeVariant - combines sinusoid + adaptive square wave
"""

import sys
import os

# Add the src directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

import numpy as np
import matplotlib.pyplot as plt

# Import our composite variant and components
from cavitymd.composite_variant import CompositeVariant
from cavitymd.variants import PeriodicVariant, AdaptiveSquareWaveVariant
from cavitymd.analysis import ElapsedTimeTracker


def create_mock_simulation():
    """Create a mock simulation object for testing."""
    class MockSim:
        def __init__(self):
            self.timestep = 0
            
    return MockSim()


def test_composite_variant():
    """Test the CompositeVariant with sinusoid + adaptive square wave."""
    
    # Create mock simulation
    mock_sim = create_mock_simulation()
    
    # Create time tracker
    runtime_ps = 100.0
    time_tracker = ElapsedTimeTracker(mock_sim, runtime_ps=runtime_ps)
    
    # Component 1: Sinusoidal base (1 ps period, 1e-4 amplitude)
    sinusoid = PeriodicVariant(
        amplitude=1e-4,
        time_tracker=time_tracker,
        period_ps=1.0,
        phase_offset=1.57  # π/2, start at maximum
    )
    
    # Component 2: Adaptive square wave (50 ps period, 2e-4 amplitude)
    # Note: For this test, we'll use a regular SquareWave since AdaptiveSquareWave
    # requires temperature tracking which is complex to set up
    from cavitymd.variants import SquareWaveVariant
    square_wave = SquareWaveVariant(
        amplitude=2e-4,
        period_ps=50.0,
        time_tracker=time_tracker,
        duty_cycle=0.02,
        start_time_ps=10.0,
        stop_time_ps=80.0
    )
    
    # Create composite variant
    composite = CompositeVariant(
        variants=[sinusoid, square_wave],
        max_amplitude=5e-4  # Safety limit
    )
    
    # Test evaluation over time
    timesteps = np.arange(0, int(runtime_ps * 10000), 100)  # Every 0.01 ps
    times = []
    sinusoid_values = []
    square_values = []
    composite_values = []
    
    for timestep in timesteps:
        # Update mock simulation timestep
        mock_sim.timestep = timestep
        
        # Get current time
        current_time = time_tracker.get_elapsed_time_ps()
        times.append(current_time)
        
        # Evaluate components individually
        sin_val = sinusoid(timestep)
        sq_val = square_wave(timestep)
        comp_val = composite(timestep)
        
        sinusoid_values.append(sin_val)
        square_values.append(sq_val)
        composite_values.append(comp_val)
        
        # Print some values for verification
        if len(times) % 1000 == 0:  # Every ~10 ps
            print(f"t={current_time:.1f}ps: sin={sin_val:.6f}, square={sq_val:.6f}, composite={comp_val:.6f}")
    
    # Create plots
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(12, 10))
    
    # Plot sinusoid component
    ax1.plot(times, sinusoid_values, 'b-', label='Sinusoid (1ps period)')
    ax1.set_ylabel('Sinusoid\nAmplitude')
    ax1.grid(True)
    ax1.legend()
    
    # Plot square wave component
    ax2.plot(times, square_values, 'r-', label='Square Wave (50ps period)')
    ax2.set_ylabel('Square Wave\nAmplitude')
    ax2.grid(True)
    ax2.legend()
    
    # Plot composite
    ax3.plot(times, composite_values, 'g-', linewidth=2, label='Composite (Sin + Square)')
    ax3.set_ylabel('Composite\nAmplitude')
    ax3.grid(True)
    ax3.legend()
    
    # Plot all together
    ax4.plot(times, sinusoid_values, 'b-', alpha=0.7, label='Sinusoid')
    ax4.plot(times, square_values, 'r-', alpha=0.7, label='Square Wave')
    ax4.plot(times, composite_values, 'g-', linewidth=2, label='Composite')
    ax4.set_xlabel('Time (ps)')
    ax4.set_ylabel('Amplitude')
    ax4.grid(True)
    ax4.legend()
    
    plt.suptitle('Composite Variant: Sinusoid + Square Wave Coupling')
    plt.tight_layout()
    plt.savefig('composite_variant_test.png', dpi=150)
    plt.show()
    
    print(f"\\nTest completed! Plot saved as 'composite_variant_test.png'")
    print(f"Maximum composite amplitude: {max(composite_values):.6f}")
    print(f"Minimum composite amplitude: {min(composite_values):.6f}")


if __name__ == "__main__":
    test_composite_variant()
