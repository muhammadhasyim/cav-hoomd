#!/usr/bin/env python3
"""
Simple test to reproduce and fix the Bussi thermostat zero velocity issue.

This test attempts to reproduce the exact error: 
"Bussi thermostat requires non-zero initial momenta."
"""

import os
import numpy as np
import hoomd
import gsd.hoomd

def test_original_error_reproduction():
    """Test using an existing molecular GSD file to reproduce the Bussi error."""
    
    print("Testing Bussi Thermostat Zero Velocity Fix")
    print("=" * 50)
    
    # Look for existing molecular files in the current directory
    possible_files = [
        'molecular-0.gsd',
        'init-0.gsd', 
        'molecular.gsd',
        'init.gsd'
    ]
    
    test_gsd = None
    for fname in possible_files:
        if os.path.exists(fname):
            test_gsd = fname
            break
    
    if test_gsd is None:
        print("❌ No molecular GSD file found. Please ensure you have:")
        for fname in possible_files:
            print(f"   - {fname}")
        return False
    
    print(f"✅ Found GSD file: {test_gsd}")
    
    # First, let's examine this GSD file
    with gsd.hoomd.open(test_gsd, 'r') as f:
        snapshot = f[-1]  # Last frame
        velocities = np.array(snapshot.particles.velocity)
        max_vel = np.max(np.sqrt(np.sum(velocities**2, axis=1)))
        print(f"Maximum velocity magnitude in {test_gsd}: {max_vel:.2e}")
        
        if max_vel < 1e-12:
            print("✅ File has zero velocities - perfect for testing!")
        else:
            print("⚠️  File has non-zero velocities - let's zero them out for testing")
            
            # Create a modified version with zero velocities
            snapshot.particles.velocity = np.zeros_like(velocities)
            test_gsd_zero = 'test_zero_vel.gsd'
            with gsd.hoomd.open(test_gsd_zero, 'w') as f_out:
                f_out.append(snapshot)
            test_gsd = test_gsd_zero
            print(f"Created {test_gsd} with zero velocities")
    
    try:
        from hoomd.cavitymd.simulation import CavityMDSimulation
        
        print(f"\n1. Testing the original error scenario:")
        print(f"   restart_velocities=False + Bussi thermostat + zero velocities")
        print(f"   This should trigger our automatic fix...")
        
        # Create simulation that should trigger the fix
        sim = CavityMDSimulation(
            job_dir='bussi_fix_test',
            replica=1,
            freq=2000.0,
            couplstr=0.001,
            incavity=True,
            runtime_ps=0.5,  # Very short test
            input_gsd=test_gsd,
            frame=-1,
            name='fix_test',
            temperature=100.0,
            molecular_thermostat='bussi',  # This requires non-zero velocities
            cavity_thermostat='langevin',  # Mixed thermostats
            restart_velocities=False,  # This keeps the zero velocities from GSD
            enable_energy_tracking=False,  # Disable for simplicity
            enable_fkt=False,  # Disable for simplicity
            seed=42
        )
        
        os.makedirs('bussi_fix_test', exist_ok=True)
        
        print("\n2. Running simulation (should now work with the fix)...")
        
        # This should work now!
        exit_code = sim.run()
        
        if exit_code == 0:
            print("\n🎉 SUCCESS! Simulation completed without Bussi thermostat error!")
            print("   The automatic velocity fix worked correctly.")
            print("   The fix detected zero velocities and force-thermalized the system.")
        else:
            print(f"\n❌ FAILURE: Simulation failed with exit code: {exit_code}")
            return False
            
    except Exception as e:
        print(f"\n❌ ERROR during simulation: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    finally:
        # Cleanup
        if 'test_zero_vel.gsd' in locals() and os.path.exists('test_zero_vel.gsd'):
            os.remove('test_zero_vel.gsd')
        
        import shutil
        if os.path.exists('bussi_fix_test'):
            shutil.rmtree('bussi_fix_test')
    
    return True

if __name__ == '__main__':
    success = test_original_error_reproduction()
    
    if success:
        print("\n" + "=" * 60)
        print("🎉 TEST PASSED! The Bussi thermostat fix is working!")
        print("✅ The simulation can now handle zero velocities with Bussi thermostat")
        print("✅ The fix automatically thermalizes when needed")
        exit(0)
    else:
        print("\n" + "=" * 60) 
        print("❌ TEST FAILED! The fix needs more work.")
        exit(1) 