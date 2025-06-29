#!/usr/bin/env python3
"""
Test script to verify that the corrected energy tracking conserves energy properly.
This runs a short simulation and checks that universe total energy is conserved.
"""

import sys
import numpy as np
from run_cavity_experiments import run_single_experiment_enhanced

def test_energy_conservation():
    """Test that the corrected energy tracking conserves energy."""
    print("="*60)
    print("TESTING CORRECTED ENERGY CONSERVATION")
    print("="*60)
    
    # Run a short test simulation with corrected energy tracking
    success = run_single_experiment_enhanced(
        exp_name='bussi_langevin_finiteq',
        molecular_thermo='bussi',
        cavity_thermo='langevin', 
        finite_q=True,
        coupling=1e-3,
        temperature=100.0,
        frequency=1560.0,
        replica=999,  # Use special test replica number
        frame=0,
        runtime_ps=1.0,  # Very short simulation
        molecular_tau=5.0,
        cavity_tau=5.0,
        log_to_file=False,
        log_to_console=True,
        enable_fkt=False,  # Disable F(k,t) for speed
        fkt_kmag=1.0,
        fkt_wavevectors=50,
        fkt_ref_interval=1.0,
        fkt_max_refs=10,
        max_energy_output_time=None,
        device='CPU',
        gpu_id=0,
        incavity=True,
        fixed_timestep=False,
        timestep_fs=1.0,
        enable_energy_tracking=True,  # Enable energy tracking
        checkpoint_interval=None,
        energy_output_period_ps=0.01,  # Frequent energy output
        fkt_output_period_ps=1.0,
        gsd_output_period_ps=50.0,
        console_output_period_ps=0.1
    )
    
    if not success:
        print("❌ Test simulation failed!")
        return False
        
    # Check energy conservation by reading the energy file
    try:
        energy_file = 'bussi_langevin_finiteq_coupling_1eneg03/prod-999-energy.txt'
        
        # Read energy data
        with open(energy_file, 'r') as f:
            lines = f.readlines()
        
        # Find header and data
        header_line = None
        data_lines = []
        for line in lines:
            if line.startswith('timestep time(ps)'):
                header_line = line.strip()
            elif not line.startswith('#') and line.strip():
                data_lines.append(line.strip())
        
        if header_line is None:
            print("❌ Could not find header in energy file!")
            return False
            
        if len(data_lines) < 10:
            print("❌ Not enough energy data points!")
            return False
        
        # Parse header to find universe total energy column
        header_parts = header_line.split()
        try:
            universe_col = header_parts.index('universe_total_energy')
        except ValueError:
            print("❌ Could not find universe_total_energy column!")
            print(f"Available columns: {header_parts}")
            return False
        
        # Extract universe total energies
        universe_energies = []
        times = []
        for line in data_lines:
            parts = line.split()
            if len(parts) > universe_col:
                times.append(float(parts[1]))  # time in ps
                universe_energies.append(float(parts[universe_col]))
        
        if len(universe_energies) < 5:
            print("❌ Not enough universe energy values!")
            return False
        
        # Check energy conservation
        universe_energies = np.array(universe_energies)
        initial_energy = universe_energies[0]
        final_energy = universe_energies[-1]
        energy_drift = final_energy - initial_energy
        max_fluctuation = np.max(np.abs(universe_energies - initial_energy))
        std_fluctuation = np.std(universe_energies)
        
        print(f"\n🔍 ENERGY CONSERVATION ANALYSIS:")
        print(f"  Simulation time: {times[0]:.4f} - {times[-1]:.4f} ps")
        print(f"  Initial universe energy: {initial_energy:.6f} a.u.")
        print(f"  Final universe energy: {final_energy:.6f} a.u.")
        print(f"  Total energy drift: {energy_drift:.6f} a.u.")
        print(f"  Max energy fluctuation: {max_fluctuation:.6f} a.u.")
        print(f"  Energy standard deviation: {std_fluctuation:.6f} a.u.")
        
        # Define conservation criteria
        ACCEPTABLE_DRIFT = 0.01  # a.u.
        ACCEPTABLE_FLUCTUATION = 0.01  # a.u.
        
        drift_ok = abs(energy_drift) < ACCEPTABLE_DRIFT
        fluctuation_ok = max_fluctuation < ACCEPTABLE_FLUCTUATION
        
        print(f"\n📊 CONSERVATION CHECKS:")
        print(f"  Energy drift < {ACCEPTABLE_DRIFT:.3f} a.u.: {'✅' if drift_ok else '❌'} ({abs(energy_drift):.6f})")
        print(f"  Max fluctuation < {ACCEPTABLE_FLUCTUATION:.3f} a.u.: {'✅' if fluctuation_ok else '❌'} ({max_fluctuation:.6f})")
        
        if drift_ok and fluctuation_ok:
            print(f"\n🎉 ENERGY CONSERVATION TEST PASSED!")
            print(f"   Universe total energy is properly conserved.")
            return True
        else:
            print(f"\n⚠️  ENERGY CONSERVATION TEST FAILED!")
            print(f"   Universe total energy shows significant drift/fluctuations.")
            print(f"   This indicates the corrected energy tracking still has issues.")
            return False
            
    except Exception as e:
        print(f"❌ Error analyzing energy file: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == '__main__':
    success = test_energy_conservation()
    sys.exit(0 if success else 1) 