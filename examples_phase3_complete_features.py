#!/usr/bin/env python3
"""
Comprehensive examples demonstrating all Phase 3 cavity MD enhancements.

This script provides complete working examples for:
1. All coupling variants (step, periodic, exponential, square wave)
2. Enhanced laser drive with timing control
3. PI feedback controller with different temperature methods
4. Combined multi-feature simulations
5. GPU vs CPU performance comparisons

Author: AI Assistant
Date: September 26, 2025
"""

import numpy as np
from pathlib import Path

def example_1_coupling_variants():
    """Example 1: All coupling variant types."""
    print("=" * 70)
    print("EXAMPLE 1: All Coupling Variants")
    print("=" * 70)
    
    examples = [
        {
            "name": "Constant Coupling",
            "params": {
                "coupling_variant_type": "constant",
                "couplstr": 0.001,
            },
            "description": "Standard constant coupling for equilibrium studies"
        },
        {
            "name": "Step Coupling",
            "params": {
                "coupling_variant_type": "step",
                "couplstr": 0.001,
                "switch_time_ps": 10.0,
                "decay_time_constant_ps": 20.0,  # Optional decay
            },
            "description": "Step coupling with exponential decay for non-equilibrium studies"
        },
        {
            "name": "Periodic Coupling",
            "params": {
                "coupling_variant_type": "periodic",
                "couplstr": 0.001,
                "periodic_period_ps": 5.0,
                "periodic_phase_offset": 0.0,
                "periodic_start_time_ps": 5.0,
                "periodic_stop_time_ps": 100.0,
            },
            "description": "Periodic coupling for dynamic cavity-molecule interactions"
        },
        {
            "name": "Exponential Decay",
            "params": {
                "coupling_variant_type": "exponential",
                "exponential_amplitude": 0.002,
                "exponential_decay_time_ps": 15.0,
                "exponential_turn_on_time_ps": 2.0,
                "exponential_turn_off_time_ps": 80.0,
            },
            "description": "Exponential decay coupling for studying decoherence effects"
        },
        {
            "name": "Square Wave",
            "params": {
                "coupling_variant_type": "square",
                "couplstr": 0.001,
                "square_period_ps": 8.0,
                "square_duty_cycle": 0.3,  # 30% on, 70% off
                "square_phase_offset": np.pi/4,
                "square_start_time_ps": 1.0,
                "square_stop_time_ps": 50.0,
            },
            "description": "Square wave coupling for digital control protocols"
        }
    ]
    
    print("Available coupling variants:\n")
    
    for i, example in enumerate(examples, 1):
        print(f"{i}. {example['name']}")
        print(f"   Description: {example['description']}")
        print(f"   Key parameters:")
        for key, value in example['params'].items():
            print(f"     {key}: {value}")
        print()
        
        # Show how to create the simulation
        print(f"   Usage example:")
        print(f"   ```python")
        print(f"   from cavitymd.simulation import CavityMDSimulation")
        print(f"   ")
        print(f"   sim = CavityMDSimulation(")
        print(f"       job_dir='{example['name'].lower().replace(' ', '_')}_run',")
        print(f"       replica=1,")
        print(f"       freq=2000.0,  # 2000 cm⁻¹")
        print(f"       incavity=True,")
        print(f"       runtime_ps=100.0,")
        for key, value in example['params'].items():
            if isinstance(value, str):
                print(f"       {key}='{value}',")
            else:
                print(f"       {key}={value},")
        print(f"   )")
        print(f"   exit_code = sim.run()")
        print(f"   ```")
        print("-" * 50)
    
    return examples

def example_2_laser_timing():
    """Example 2: Enhanced laser drive with timing control."""
    print("=" * 70)
    print("EXAMPLE 2: Enhanced Laser Drive with Timing")
    print("=" * 70)
    
    laser_examples = [
        {
            "name": "Basic Laser Pulse",
            "params": {
                "laser_enabled": True,
                "laser_amplitude": 1e-6,
                "laser_kvector": [1.0, 0.0, 0.0],
                "laser_start_time_ps": 10.0,
                "laser_stop_time_ps": 30.0,
            },
            "description": "Short laser pulse for perturbation studies"
        },
        {
            "name": "Continuous Laser",
            "params": {
                "laser_enabled": True,
                "laser_amplitude": 5e-7,
                "laser_kvector": [0.0, 1.0, 0.0],
                "laser_start_time_ps": 5.0,
                "laser_stop_time_ps": None,  # Runs until end
            },
            "description": "Continuous laser for steady-state studies"
        },
        {
            "name": "Multiple Pulse Sequence",
            "params": {
                "laser_enabled": True,
                "laser_amplitude": 2e-6,
                "laser_kvector": [1.0, 1.0, 0.0],
                "laser_start_time_ps": 20.0,
                "laser_stop_time_ps": 25.0,
            },
            "description": "Short pulse (can be repeated with multiple simulations)"
        }
    ]
    
    print("Enhanced laser drive examples:\n")
    
    for i, example in enumerate(laser_examples, 1):
        print(f"{i}. {example['name']}")
        print(f"   Description: {example['description']}")
        print(f"   Configuration:")
        for key, value in example['params'].items():
            print(f"     {key}: {value}")
        print()
        
        print(f"   Usage example:")
        print(f"   ```python")
        print(f"   sim = CavityMDSimulation(")
        print(f"       job_dir='{example['name'].lower().replace(' ', '_')}_laser',")
        print(f"       replica=1,")
        print(f"       freq=2000.0,")
        print(f"       couplstr=0.001,")
        print(f"       incavity=True,")
        print(f"       runtime_ps=50.0,")
        for key, value in example['params'].items():
            if value is None:
                print(f"       {key}=None,")
            elif isinstance(value, list):
                print(f"       {key}={value},")
            elif isinstance(value, str):
                print(f"       {key}='{value}',")
            else:
                print(f"       {key}={value},")
        print(f"   )")
        print(f"   exit_code = sim.run()")
        print(f"   ```")
        print("-" * 50)
    
    return laser_examples

def example_3_pi_feedback():
    """Example 3: PI feedback controller for all temperature methods."""
    print("=" * 70)
    print("EXAMPLE 3: PI Feedback Controller")
    print("=" * 70)
    
    feedback_examples = [
        {
            "name": "Kinetic Temperature Control",
            "params": {
                "enable_pi_feedback": True,
                "pi_target_temperature": 120.0,
                "pi_temperature_method": "kinetic",
                "pi_turn_on_time_ps": 10.0,
                "pi_turn_off_time_ps": 80.0,
                "pi_update_interval_ps": 2.0,
                "pi_Kc": None,  # Auto-tune
                "pi_Ti": None,  # Auto-tune
            },
            "description": "Control using kinetic temperature (from particle velocities)"
        },
        {
            "name": "LJ+Coulombic Fictive Temperature",
            "params": {
                "enable_pi_feedback": True,
                "pi_target_temperature": 100.0,
                "pi_temperature_method": "lj_coulombic",
                "pi_turn_on_time_ps": 5.0,
                "pi_update_interval_ps": 5.0,
                "pi_Kc": 1.5,  # Manual tuning
                "pi_Ti": 8.0,  # Manual tuning
                "pi_beta": 0.8,  # Setpoint weighting
            },
            "description": "Control using LJ+Coulombic potential energy with Rosenfeld scaling"
        },
        {
            "name": "Harmonic Fictive Temperature",
            "params": {
                "enable_pi_feedback": True,
                "pi_target_temperature": 110.0,
                "pi_temperature_method": "harmonic",
                "pi_turn_on_time_ps": 0.0,  # Start immediately
                "pi_update_interval_ps": 3.0,
                "pi_molecular_tau_ps": 6.0,  # For auto-tuning
                "pi_T_min": 80.0,
                "pi_T_max": 140.0,
            },
            "description": "Control using harmonic potential energy with T^(3/2) scaling"
        }
    ]
    
    print("PI feedback controller examples:\n")
    
    for i, example in enumerate(feedback_examples, 1):
        print(f"{i}. {example['name']}")
        print(f"   Description: {example['description']}")
        print(f"   Configuration:")
        for key, value in example['params'].items():
            print(f"     {key}: {value}")
        print()
        
        print(f"   Usage example:")
        print(f"   ```python")
        print(f"   sim = CavityMDSimulation(")
        print(f"       job_dir='{example['name'].lower().replace(' ', '_')}_feedback',")
        print(f"       replica=1,")
        print(f"       freq=2000.0,")
        print(f"       couplstr=0.001,")
        print(f"       incavity=True,")
        print(f"       runtime_ps=100.0,")
        print(f"       enable_energy_tracking=True,  # Required for fictive temp")
        for key, value in example['params'].items():
            if value is None:
                print(f"       {key}=None,")
            elif isinstance(value, str):
                print(f"       {key}='{value}',")
            else:
                print(f"       {key}={value},")
        print(f"   )")
        print(f"   exit_code = sim.run()")
        print(f"   ```")
        print("-" * 50)
    
    return feedback_examples

def example_4_combined_features():
    """Example 4: Combined multi-feature simulations."""
    print("=" * 70)
    print("EXAMPLE 4: Combined Multi-Feature Simulations")
    print("=" * 70)
    
    combined_examples = [
        {
            "name": "Periodic Coupling + Laser + PI Feedback",
            "description": "Dynamic cavity coupling with laser perturbation and temperature control",
            "features": ["Periodic coupling", "Laser timing", "PI feedback"],
            "params": {
                # Periodic coupling
                "coupling_variant_type": "periodic",
                "couplstr": 0.001,
                "periodic_period_ps": 10.0,
                "periodic_start_time_ps": 5.0,
                "periodic_stop_time_ps": 80.0,
                
                # Enhanced laser
                "laser_enabled": True,
                "laser_amplitude": 1e-6,
                "laser_start_time_ps": 20.0,
                "laser_stop_time_ps": 60.0,
                
                # PI feedback
                "enable_pi_feedback": True,
                "pi_target_temperature": 105.0,
                "pi_temperature_method": "kinetic",
                "pi_turn_on_time_ps": 0.0,
                
                # Required for complex simulations
                "enable_energy_tracking": True,
                "runtime_ps": 100.0,
            }
        },
        {
            "name": "Exponential Decay + Harmonic Temperature Control",
            "description": "Studying decoherence with harmonic temperature feedback",
            "features": ["Exponential decay coupling", "Harmonic fictive temperature"],
            "params": {
                # Exponential coupling
                "coupling_variant_type": "exponential",
                "exponential_amplitude": 0.0015,
                "exponential_decay_time_ps": 20.0,
                "exponential_turn_on_time_ps": 10.0,
                
                # PI feedback with harmonic method
                "enable_pi_feedback": True,
                "pi_target_temperature": 100.0,
                "pi_temperature_method": "harmonic",
                "pi_update_interval_ps": 2.0,
                
                # Enhanced tracking
                "enable_energy_tracking": True,
                "enable_fkt": True,
                "runtime_ps": 150.0,
            }
        },
        {
            "name": "Square Wave + Laser + Advanced Tracking",
            "description": "Digital coupling protocol with laser and comprehensive analysis",
            "features": ["Square wave coupling", "Laser pulses", "F(k,t) analysis"],
            "params": {
                # Square wave coupling
                "coupling_variant_type": "square",
                "couplstr": 0.001,
                "square_period_ps": 6.0,
                "square_duty_cycle": 0.25,  # 25% on
                
                # Laser perturbation
                "laser_enabled": True,
                "laser_amplitude": 5e-7,
                "laser_start_time_ps": 30.0,
                "laser_stop_time_ps": 40.0,
                
                # Advanced analysis
                "enable_energy_tracking": True,
                "enable_fkt": True,
                "fkt_log_time_spacing": True,
                "fkt_min_log_time_ps": 0.1,
                "fkt_max_log_time_ps": 50.0,
                "runtime_ps": 80.0,
            }
        }
    ]
    
    print("Combined feature examples:\n")
    
    for i, example in enumerate(combined_examples, 1):
        print(f"{i}. {example['name']}")
        print(f"   Description: {example['description']}")
        print(f"   Features: {', '.join(example['features'])}")
        print()
        
        print(f"   Complete simulation setup:")
        print(f"   ```python")
        print(f"   from cavitymd.simulation import CavityMDSimulation")
        print(f"   ")
        print(f"   # {example['description']}")
        print(f"   sim = CavityMDSimulation(")
        print(f"       job_dir='{example['name'].lower().replace(' ', '_')}_combined',")
        print(f"       replica=1,")
        print(f"       freq=2000.0,  # 2000 cm⁻¹ cavity frequency")
        print(f"       incavity=True,")
        print(f"       device='GPU',  # Use GPU for performance")
        
        for key, value in example['params'].items():
            if isinstance(value, str):
                print(f"       {key}='{value}',")
            elif isinstance(value, bool):
                print(f"       {key}={value},")
            else:
                print(f"       {key}={value},")
        
        print(f"   )")
        print(f"   ")
        print(f"   # Run the simulation")
        print(f"   print(' Starting {example['name']} simulation...')")
        print(f"   exit_code = sim.run()")
        print(f"   print(f' Simulation completed with exit code {{exit_code}}')")
        print(f"   ```")
        print("-" * 50)
    
    return combined_examples

def example_5_gpu_performance():
    """Example 5: GPU vs CPU performance comparison."""
    print("=" * 70)
    print("EXAMPLE 5: GPU vs CPU Performance Optimization")
    print("=" * 70)
    
    print("Performance optimization guidelines:\n")
    
    performance_tips = [
        {
            "category": "Device Selection",
            "tips": [
                "Use device='GPU' for large systems (>1000 particles)",
                "Use device='CPU' for small systems or debugging",
                "GPU shows best performance with complex force calculations",
                "All new features work identically on CPU and GPU"
            ]
        },
        {
            "category": "Coupling Variants",
            "tips": [
                "Constant coupling: Fastest (no time-dependent calculations)",
                "Step coupling: Good performance, one-time setup cost",
                "Periodic/Square: Moderate overhead from trigonometric functions", 
                "Exponential: Minimal overhead from exponential calculation"
            ]
        },
        {
            "category": "Feedback Controllers",
            "tips": [
                "PI feedback: Minimal performance impact (Python-side only)",
                "Kinetic temperature: Fastest calculation method",
                "Fictive temperatures: Require energy tracking (slight overhead)",
                "Update intervals >1 ps recommended for efficiency"
            ]
        },
        {
            "category": "Laser Drive",
            "tips": [
                "Laser force: GPU-optimized kernel, excellent performance",
                "Timing control: Zero performance impact (Python logic)",
                "Multiple laser forces: Linear scaling with number of forces"
            ]
        }
    ]
    
    for category in performance_tips:
        print(f"**{category['category']}:**")
        for tip in category['tips']:
            print(f"  • {tip}")
        print()
    
    print("Performance comparison example:")
    print("```python")
    print("import time")
    print("from cavitymd.simulation import CavityMDSimulation")
    print("")
    print("# Test configuration")
    print("base_params = {")
    print("    'replica': 1,")
    print("    'freq': 2000.0,")
    print("    'couplstr': 0.001,")
    print("    'incavity': True,")
    print("    'runtime_ps': 50.0,")
    print("    'coupling_variant_type': 'periodic',")
    print("    'periodic_period_ps': 5.0,")
    print("    'laser_enabled': True,")
    print("    'laser_amplitude': 1e-6,")
    print("    'enable_pi_feedback': True,")
    print("}")
    print("")
    print("# GPU performance test")
    print("start_time = time.time()")
    print("gpu_sim = CavityMDSimulation(")
    print("    job_dir='gpu_performance_test',")
    print("    device='GPU',")
    print("    **base_params")
    print(")")
    print("gpu_exit = gpu_sim.run()")
    print("gpu_time = time.time() - start_time")
    print("")
    print("# CPU performance test")
    print("start_time = time.time()")
    print("cpu_sim = CavityMDSimulation(")
    print("    job_dir='cpu_performance_test',")
    print("    device='CPU',")
    print("    **base_params")
    print(")")
    print("cpu_exit = cpu_sim.run()")
    print("cpu_time = time.time() - start_time")
    print("")
    print("# Performance comparison")
    print("speedup = cpu_time / gpu_time if gpu_time > 0 else float('inf')")
    print("print(f'GPU time: {gpu_time:.2f} seconds')")
    print("print(f'CPU time: {cpu_time:.2f} seconds')")
    print("print(f'GPU speedup: {speedup:.2f}x')")
    print("```")
    
    return performance_tips

def main():
    """Run all Phase 3 examples."""
    print(" COMPREHENSIVE PHASE 3 CAVITY MD EXAMPLES")
    print("=" * 80)
    print("This script demonstrates all enhanced features from Phase 3 development:")
    print("• Enhanced coupling variants (step, periodic, exponential, square wave)")
    print("• Laser drive with timing control (start/stop times)")
    print("• PI feedback controller (kinetic, LJ+Coulombic, harmonic temperature)")
    print("• Combined multi-feature simulations")
    print("• GPU vs CPU performance optimization")
    print("=" * 80)
    print()
    
    examples = []
    
    # Run all examples
    examples.append(("Coupling Variants", example_1_coupling_variants()))
    examples.append(("Laser Timing", example_2_laser_timing()))
    examples.append(("PI Feedback", example_3_pi_feedback()))
    examples.append(("Combined Features", example_4_combined_features()))
    examples.append(("GPU Performance", example_5_gpu_performance()))
    
    # Summary
    print("=" * 80)
    print(" SUMMARY OF ALL EXAMPLES")
    print("=" * 80)
    
    total_examples = sum(len(ex[1]) if isinstance(ex[1], list) else 1 for ex in examples)
    
    print(f"Total examples demonstrated: {total_examples}")
    print()
    print("Example categories:")
    for category, example_list in examples:
        if isinstance(example_list, list):
            print(f"  • {category}: {len(example_list)} examples")
        else:
            print(f"  • {category}: 1 comprehensive guide")
    
    print()
    print(" NEXT STEPS:")
    print("1. Choose an example that matches your research needs")
    print("2. Modify parameters for your specific system")
    print("3. Run simulations and analyze results")
    print("4. Combine features for advanced experiments")
    print("5. Use GPU acceleration for large-scale studies")
    
    print()
    print(" All features are production-ready and GPU-optimized!")
    print("   For questions, see documentation or example scripts.")

if __name__ == "__main__":
    main()
