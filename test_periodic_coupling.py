#!/usr/bin/env python3
"""
Validation script for periodic coupling constant functionality.

This script demonstrates and validates the new periodic coupling protocol
where the coupling constant oscillates sinusoidally with a specified period.

The coupling strength varies as:
g(t) = A * sin(2π*t/T + φ)

where:
- A is the amplitude (the couplstr parameter)
- T is the period in picoseconds  
- φ is the phase offset
- t is the elapsed time
"""

import numpy as np
import matplotlib.pyplot as plt
from unittest.mock import Mock
import sys
import os

def test_periodic_variant_basic():
    """Test basic PeriodicVariant functionality."""
    print("=" * 60)
    print("Testing PeriodicVariant Basic Functionality")
    print("=" * 60)
    
    try:
        from hoomd.cavitymd.variants import PeriodicVariant
        from hoomd.cavitymd.analysis import ElapsedTimeTracker
        
        # Create mock time tracker
        time_tracker = Mock()
        
        # Test 1: Basic oscillation with 1 ps period, 0.001 amplitude
        print("Test 1: Basic oscillation (period=1.0 ps, amplitude=0.001)")
        variant = PeriodicVariant(
            amplitude=0.001,
            time_tracker=time_tracker,
            period_ps=1.0,
            phase_offset=0.0,
            start_time_ps=0.0
        )
        
        print(f"  Amplitude: {variant.amplitude}")
        print(f"  Period: {variant.period_ps} ps")
        print(f"  Frequency: {variant.get_frequency_thz():.1f} THz")
        print(f"  Phase offset: {variant.phase_offset:.3f} rad")
        
        # Test values at key time points
        times = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
        expected = 0.001 * np.sin(2 * np.pi * times)
        
        for i, t in enumerate(times):
            time_tracker.elapsed_time = t
            value = variant(i)
            print(f"  t={t:.2f} ps: g(t)={value:.6f} (expected: {expected[i]:.6f})")
            assert abs(value - expected[i]) < 1e-12, f"Value mismatch at t={t}"
        
        print("  ✓ Basic oscillation test passed")
        
        # Test 2: Different period and phase
        print("\nTest 2: Different period and phase (period=2.0 ps, phase=π/2)")
        variant2 = PeriodicVariant(
            amplitude=0.002,
            time_tracker=time_tracker,
            period_ps=2.0,
            phase_offset=np.pi/2,
            start_time_ps=0.0
        )
        
        # At t=0, with phase π/2, should be at maximum
        time_tracker.elapsed_time = 0.0
        value = variant2(0)
        expected = 0.002 * np.sin(np.pi/2)  # = 0.002
        print(f"  t=0.0 ps with π/2 phase: g(t)={value:.6f} (expected: {expected:.6f})")
        assert abs(value - expected) < 1e-12
        
        print("  ✓ Phase offset test passed")
        
        # Test 3: Start time delay
        print("\nTest 3: Delayed start (start_time=1.0 ps)")
        variant3 = PeriodicVariant(
            amplitude=0.001,
            time_tracker=time_tracker,
            period_ps=1.0,
            phase_offset=0.0,
            start_time_ps=1.0
        )
        
        # Before start time - should be zero
        time_tracker.elapsed_time = 0.5
        value = variant3(0)
        print(f"  t=0.5 ps (before start): g(t)={value:.6f} (expected: 0.0)")
        assert value == 0.0
        
        # At start time - should begin oscillation
        time_tracker.elapsed_time = 1.25  # 0.25 ps after start
        value = variant3(1)
        expected = 0.001 * np.sin(2 * np.pi * 0.25)
        print(f"  t=1.25 ps (0.25 ps after start): g(t)={value:.6f} (expected: {expected:.6f})")
        assert abs(value - expected) < 1e-12
        
        print("  ✓ Delayed start test passed")
        
        print("\n✓ All PeriodicVariant tests passed!")
        return True
        
    except Exception as e:
        print(f"✗ PeriodicVariant test failed: {e}")
        return False


def test_simulation_integration():
    """Test integration with CavityMDSimulation class."""
    print("\n" + "=" * 60)
    print("Testing CavityMDSimulation Integration")
    print("=" * 60)
    
    try:
        from hoomd.cavitymd.simulation import CavityMDSimulation
        
        # Test creating a simulation with periodic coupling
        print("Creating CavityMDSimulation with periodic coupling...")
        
        sim = CavityMDSimulation(
            job_dir="test_periodic",
            replica=1,
            freq=2000.0,  # 2000 cm⁻¹
            couplstr=0.001,  # Amplitude
            incavity=True,
            runtime_ps=10.0,
            periodic_coupling=True,
            periodic_period_ps=1.0,
            periodic_phase_offset=0.0,
            periodic_start_time_ps=0.0,
            enable_energy_tracking=False,
            enable_fkt=False
        )
        
        print(f"  ✓ Simulation created successfully")
        print(f"  Periodic coupling: {sim.periodic_coupling}")
        print(f"  Period: {sim.periodic_period_ps} ps")
        print(f"  Phase offset: {sim.periodic_phase_offset} rad")
        print(f"  Start time: {sim.periodic_start_time_ps} ps")
        
        # Test parameter validation
        print("\nTesting parameter validation...")
        
        # Test negative period (should raise error when creating PeriodicVariant)
        try:
            bad_sim = CavityMDSimulation(
                job_dir="test_bad",
                replica=1,
                freq=2000.0,
                couplstr=0.001,
                incavity=True,
                runtime_ps=10.0,
                periodic_coupling=True,
                periodic_period_ps=-1.0,  # Invalid
                enable_energy_tracking=False,
                enable_fkt=False
            )
            print("  ✗ Should have failed with negative period")
            return False
        except:
            print("  ✓ Correctly rejected negative period")
        
        print("\n✓ CavityMDSimulation integration test passed!")
        return True
        
    except Exception as e:
        print(f"✗ Integration test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def plot_periodic_coupling_examples():
    """Create plots showing different periodic coupling patterns."""
    print("\n" + "=" * 60)
    print("Generating Periodic Coupling Example Plots")
    print("=" * 60)
    
    try:
        from hoomd.cavitymd.variants import PeriodicVariant
        from unittest.mock import Mock
        
        # Create time points
        time_ps = np.linspace(0, 5, 1000)
        
        # Create mock time tracker
        time_tracker = Mock()
        
        # Different periodic coupling examples
        variants = [
            {
                'name': 'Standard (1 ps period)',
                'amplitude': 0.001,
                'period_ps': 1.0,
                'phase_offset': 0.0,
                'start_time_ps': 0.0,
                'style': '-',
                'color': 'blue'
            },
            {
                'name': 'High frequency (0.5 ps period)',
                'amplitude': 0.001,
                'period_ps': 0.5,
                'phase_offset': 0.0,
                'start_time_ps': 0.0,
                'style': '--',
                'color': 'red'
            },
            {
                'name': 'Phase shifted (π/2)',
                'amplitude': 0.001,
                'period_ps': 1.0,
                'phase_offset': np.pi/2,
                'start_time_ps': 0.0,
                'style': '-.',
                'color': 'green'
            },
            {
                'name': 'Delayed start (2 ps)',
                'amplitude': 0.001,
                'period_ps': 1.0,
                'phase_offset': 0.0,
                'start_time_ps': 2.0,
                'style': ':',
                'color': 'orange'
            }
        ]
        
        plt.figure(figsize=(12, 8))
        
        for variant_params in variants:
            variant = PeriodicVariant(
                amplitude=variant_params['amplitude'],
                time_tracker=time_tracker,
                period_ps=variant_params['period_ps'],
                phase_offset=variant_params['phase_offset'],
                start_time_ps=variant_params['start_time_ps']
            )
            
            # Calculate coupling values
            coupling_values = []
            for t in time_ps:
                time_tracker.elapsed_time = t
                coupling_values.append(variant(0))
            
            plt.plot(time_ps, coupling_values, 
                    linestyle=variant_params['style'],
                    color=variant_params['color'],
                    linewidth=2,
                    label=variant_params['name'])
        
        plt.xlabel('Time (ps)', fontsize=12)
        plt.ylabel('Coupling Strength g(t) (a.u.)', fontsize=12)
        plt.title('Periodic Coupling Constant Examples', fontsize=14, fontweight='bold')
        plt.legend(fontsize=10)
        plt.grid(True, alpha=0.3)
        plt.xlim(0, 5)
        
        # Add annotations
        plt.text(0.5, 0.0008, 'g(t) = A sin(2πt/T + φ)', 
                fontsize=11, bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue"))
        
        plt.tight_layout()
        
        # Save plot
        output_file = 'periodic_coupling_examples.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"  ✓ Plot saved as: {output_file}")
        
        plt.show()
        
        return True
        
    except Exception as e:
        print(f"✗ Plotting failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all validation tests."""
    print("Periodic Coupling Validation Script")
    print("=" * 60)
    
    # Run tests
    success = True
    
    success &= test_periodic_variant_basic()
    success &= test_simulation_integration()
    
    # Only try plotting if matplotlib is available
    try:
        import matplotlib.pyplot as plt
        success &= plot_periodic_coupling_examples()
    except ImportError:
        print("\nSkipping plots (matplotlib not available)")
    
    # Final summary
    print("\n" + "=" * 60)
    if success:
        print("✓ ALL TESTS PASSED - Periodic coupling implementation is working correctly!")
        print("\nUsage example:")
        print("""
from hoomd.cavitymd.simulation import CavityMDSimulation

# Create simulation with periodic coupling
sim = CavityMDSimulation(
    job_dir="periodic_experiment",
    replica=1,
    freq=2000.0,             # Cavity frequency (cm⁻¹)
    couplstr=0.001,          # Amplitude of oscillation (a.u.)
    incavity=True,
    runtime_ps=100.0,
    periodic_coupling=True,
    periodic_period_ps=1.0,  # 1 ps period (1 THz)
    periodic_phase_offset=0.0,     # Start at zero
    periodic_start_time_ps=0.0,    # Start immediately
    enable_energy_tracking=True
)

# Run the simulation
exit_code = sim.run()
        """)
    else:
        print("✗ SOME TESTS FAILED - Please check the implementation")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
