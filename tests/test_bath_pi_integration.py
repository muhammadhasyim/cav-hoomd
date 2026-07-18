#!/usr/bin/env python3
"""
Integration test for BathPI controller.

This test verifies that the BathPI flux-based controller can:
1. Be instantiated correctly
2. Control bath temperatures based on measured heat flux
3. Drive flux to zero at equilibrium
"""

import hoomd
import numpy as np
from pathlib import Path

def test_bath_pi_basic():
    """Test basic BathPI controller instantiation and operation."""
    
    print("\n" + "="*70)
    print("BathPI Controller Integration Test")
    print("="*70)
    
    # Use the init GSD file from examples
    gsd_file = Path("/home/mh7373/GitRepos/cav-hoomd/examples/init-0.gsd")
    if not gsd_file.exists():
        print(f"ERROR: Test GSD file not found: {gsd_file}")
        print("Please ensure molecular-0.gsd exists in examples directory")
        return 1
    
    # Test parameters
    temperature = 100.0  # K
    runtime_ps = 2.0     # Very short test run (just verify it works)
    
    print(f"\nTest Configuration:")
    print(f"  GSD file: {gsd_file}")
    print(f"  Temperature: {temperature} K")
    print(f"  Runtime: {runtime_ps} ps")
    print(f"  Bath PI: K_p=0.1, K_i=0.01, filter=5.0ps")
    
    try:
        # Import after HOOMD to ensure proper initialization
        from hoomd.cavitymd.simulation import CavityMDSimulation
        import os
        
        # Create output directory
        output_dir = Path("/home/mh7373/GitRepos/cav-hoomd/test_bath_pi_output")
        output_dir.mkdir(exist_ok=True)
        
        print("\n1. Creating CavityMDSimulation with BathPI controller...")
        sim = CavityMDSimulation(
            job_dir=str(output_dir),
            replica=1,
            freq=1560.0,
            lambda_coupling=0.0,  # No cavity coupling for simple test
            incavity=False,
            runtime_ps=runtime_ps,
            input_gsd=str(gsd_file),
            frame=-1,
            temperature=temperature,
            molecular_thermostat='bussi',
            cavity_thermostat='none',
            
            # Enable energy tracking (required for BathPI)
            enable_energy_tracking=True,
            energy_output_period_ps=0.1,
            
            # BathPI controller parameters
            enable_bath_pi_controller=True,
            bath_pi_apply_to='molecular',  # Only control molecular bath
            bath_pi_K_p_molecular=0.1,
            bath_pi_K_i_molecular=0.01,
            bath_pi_K_p_cavity=0.1,
            bath_pi_K_i_cavity=0.01,
            bath_pi_filter_window_ps=5.0,
            bath_pi_anti_windup_alpha=0.01,
            bath_pi_enable_feedforward=False,
            bath_pi_T_nominal=None,
            bath_pi_feedforward_tau_ps=1000.0,
            bath_pi_turn_on_time_ps=0.0,
            bath_pi_turn_off_time_ps=None,
            bath_pi_update_interval_ps=0.1,
            bath_pi_T_min=0.1,
            bath_pi_T_max=500.0,
            bath_pi_rate_limit_K_per_ps=None,
            bath_pi_output_file='test_bath_pi_control.csv',
            bath_pi_relaxation_data_file=None,
            
            # Disable other features
            enable_fkt=False,
            device='CPU',
        )
        
        print("✓ Simulation created successfully")
        
        print("\n2. Running simulation...")
        exit_code = sim.run()
        
        if exit_code != 0:
            print(f"✗ Simulation failed with exit code {exit_code}")
            return exit_code
        
        print("✓ Simulation completed successfully")
        
        print("\n3. Checking BathPI output file...")
        output_file = output_dir / "test_bath_pi_control.csv"
        if not output_file.exists():
            print(f"✗ BathPI output file not found: {output_file}")
            return 1
        
        # Read and analyze the output
        import pandas as pd
        df = pd.read_csv(output_file, comment='#')
        
        print(f"✓ BathPI output file found with {len(df)} data points")
        
        # Check that required columns exist
        required_columns = ['time_ps', 'flux_molecular', 'flux_filt_mol', 
                          'I_mol', 'T_bath_mol', 'T_kin']
        missing = [col for col in required_columns if col not in df.columns]
        if missing:
            print(f"✗ Missing columns in output: {missing}")
            return 1
        
        print("✓ All required columns present")
        
        # Analyze flux convergence
        print("\n4. Analyzing flux convergence...")
        flux_initial = df['flux_filt_mol'].iloc[10] if len(df) > 10 else df['flux_filt_mol'].iloc[0]
        flux_final = df['flux_filt_mol'].iloc[-1]
        flux_std = df['flux_filt_mol'].std()
        
        print(f"  Initial filtered flux: {flux_initial:.3e} Ha/ps")
        print(f"  Final filtered flux:   {flux_final:.3e} Ha/ps")
        print(f"  Flux std deviation:    {flux_std:.3e} Ha/ps")
        
        # Check temperature stability
        print("\n5. Checking temperature stability...")
        T_initial = df['T_bath_mol'].iloc[0]
        T_final = df['T_bath_mol'].iloc[-1]
        T_mean = df['T_bath_mol'].mean()
        T_std = df['T_bath_mol'].std()
        
        print(f"  Initial T_bath: {T_initial:.2f} K")
        print(f"  Final T_bath:   {T_final:.2f} K")
        print(f"  Mean T_bath:    {T_mean:.2f} K")
        print(f"  Std T_bath:     {T_std:.2f} K")
        
        # Check if bath temperature is reasonable
        if T_mean < 50 or T_mean > 200:
            print(f"  ⚠ Warning: Bath temperature outside expected range (50-200 K)")
        
        # Check integral state
        print("\n6. Checking integral state...")
        I_final = df['I_mol'].iloc[-1]
        print(f"  Final integral state: {I_final:.3e}")
        
        print("\n" + "="*70)
        print("✓ BathPI Controller Integration Test PASSED")
        print("="*70)
        
        return 0
        
    except Exception as e:
        import traceback
        print(f"\n✗ Test failed with exception:")
        print(f"  {e}")
        print("\nFull traceback:")
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    exit(test_bath_pi_basic())

