#!/usr/bin/env python3
"""
Fixed thermalization example for cavity MD simulations.

This script demonstrates the proper way to thermalize cavity MD systems
to prevent temperature drift from 100K to 120K.
"""

import sys
import os

# Add the src directory to the path to import project modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from cavitymd.simulation import CavityMDSimulation

def run_fixed_simulation():
    """
    Run cavity MD simulation with fixed thermalization.
    
    Key fixes:
    1. Use matched thermostats (both Bussi)
    2. Consistent thermalization method
    3. Proper initial conditions
    """
    
    print("FIXED CAVITY MD SIMULATION")
    print("="*50)
    print("Fixes applied:")
    print("1. Matched thermostats (both Bussi)")
    print("2. Consistent HOOMD thermalization")
    print("3. Short test to verify temperature stability")
    print("="*50)
    
    # Create output directory
    os.makedirs('fixed_test', exist_ok=True)
    
    # Fixed simulation parameters
    params = {
        'job_dir': 'fixed_test',
        'replica': 1,
        'freq': 2000.0,
        'couplstr': 0.001,
        'incavity': True,
        'temperature': 100.0,
        'runtime_ps': 10.0,  # Short test
        'input_gsd': 'test.gsd',
        'name': 'fixed',
        'device': 'CPU',
        'seed': 42,
        
        # FIX #1: Use matched thermostats
        'molecular_thermostat': 'bussi',
        'cavity_thermostat': 'bussi',  # SAME as molecular
        
        # FIX #2: Enable energy tracking to monitor temperature
        'enable_energy_tracking': True,
        'energy_output_period_ps': 0.1,
        
        # Use consistent thermalization (handled by modified simulation class)
        'restart_velocities': True
    }
    
    try:
        # Create and run simulation with fixes
        sim = FixedCavityMDSimulation(**params)
        exit_code = sim.run()
        
        if exit_code == 0:
            print("\n Fixed simulation completed successfully!")
            
            # Analyze results
            energy_file = f"fixed_test/fixed-1_energy_tracker.txt"
            if os.path.exists(energy_file):
                analyze_temperature_stability(energy_file)
            else:
                print("Energy file not found - check simulation output")
        else:
            print(" Fixed simulation failed")
            
    except Exception as e:
        print(f" Error: {e}")
        import traceback
        traceback.print_exc()

def analyze_temperature_stability(filepath):
    """Analyze temperature stability from energy output."""
    try:
        import numpy as np
        data = np.loadtxt(filepath, skiprows=1)
        
        if len(data.shape) == 1:
            data = data.reshape(1, -1)
        
        # Temperature is typically the last column
        temperatures = data[:, -1]
        times = data[:, 0]
        
        initial_temp = temperatures[0]
        final_temp = temperatures[-1]
        max_temp = np.max(temperatures)
        min_temp = np.min(temperatures)
        mean_temp = np.mean(temperatures)
        std_temp = np.std(temperatures)
        
        print(f"\n TEMPERATURE ANALYSIS:")
        print(f"   Initial temperature: {initial_temp:.1f} K")
        print(f"   Final temperature:   {final_temp:.1f} K")
        print(f"   Mean temperature:    {mean_temp:.1f} ± {std_temp:.1f} K")
        print(f"   Temperature range:   {min_temp:.1f} - {max_temp:.1f} K")
        print(f"   Total drift:         {final_temp - initial_temp:+.1f} K")
        
        # Assessment
        drift = abs(final_temp - initial_temp)
        if drift < 2.0:
            print(f"    EXCELLENT: Temperature drift < 2K")
        elif drift < 5.0:
            print(f"    GOOD: Temperature drift < 5K")
        elif drift < 10.0:
            print(f"     MODERATE: Temperature drift 5-10K")
        else:
            print(f"    POOR: Temperature drift > 10K")
        
        # Check if problem is solved
        if drift < 5.0 and mean_temp > 95.0 and mean_temp < 105.0:
            print(f"\n SUCCESS: Temperature drift issue appears to be FIXED!")
        else:
            print(f"\n Further investigation needed")
            
    except Exception as e:
        print(f"Error analyzing temperature: {e}")


class FixedCavityMDSimulation(CavityMDSimulation):
    """
    Fixed version of CavityMDSimulation with proper thermalization.
    
    Key fixes:
    1. Consistent thermalization for all particles
    2. Better cavity particle initialization
    """
    
    def thermalize_system(self):
        """
        FIXED: Use consistent thermalization for all particles.
        
        This fixes the main issue causing temperature drift by using
        HOOMD's built-in thermalization for ALL particles instead of
        mixing HOOMD and manual methods.
        """
        
        if not self.restart_velocities:
            self.log_info("Velocity restart disabled - keeping existing velocities from GSD file")
            return
        
        kT = self.kB * self.temperature
        
        # Set numpy seed for reproducible thermalization if seed is provided
        if self.seed is not None:
            import numpy as np
            np.random.seed(self.seed + 1)
            self.log_info(f"Using seed {self.seed + 1} for thermalization randomness")
        
        # FIX: Use HOOMD's thermalization for ALL particles (consistent method)
        import hoomd
        all_filter = hoomd.filter.All()
        self.sim.state.thermalize_particle_momenta(kT=kT, filter=all_filter)
        
        self.log_info(" FIXED: Thermalized ALL particles using consistent HOOMD State method")
        self.log_info(f"   Target temperature: {self.temperature:.1f} K")
        self.log_info(f"   Thermal energy kT: {kT:.6f} a.u.")
        self.log_info("   This should prevent temperature drift issues")
        
        # Additional verification
        self._verify_thermalization(kT)
    
    def _verify_thermalization(self, kT):
        """Verify that thermalization was successful."""
        try:
            with self.get_local_snapshot() as snap:
                # Get velocities and masses for all particles
                velocities = self.to_numpy(snap.particles.velocity)
                masses = self.to_numpy(snap.particles.mass)
                
                # Calculate kinetic energy
                ke_per_particle = 0.5 * masses[:, np.newaxis] * velocities**2
                total_ke = np.sum(ke_per_particle)
                
                # Calculate temperature
                N_dof = 3 * len(masses)
                calculated_temp = (2.0 * total_ke) / (N_dof * self.kB)
                
                self.log_info(f"   Verification: Calculated temperature = {calculated_temp:.1f} K")
                
                error = abs(calculated_temp - self.temperature) / self.temperature
                if error < 0.1:  # 10% tolerance
                    self.log_info(f"    Thermalization verification PASSED (error: {error:.1%})")
                else:
                    self.log_warning(f"     Thermalization verification WARNING (error: {error:.1%})")
                    
        except Exception as e:
            self.log_warning(f"   Could not verify thermalization: {e}")


def main():
    """Main function to run the fixed simulation."""
    
    print("HOOMD-blue Cavity MD: Fixed Thermalization Test")
    print("="*60)
    print("Testing the fix for temperature drift from 100K to 120K")
    print("="*60)
    
    # Check if test GSD file exists
    if not os.path.exists('test.gsd'):
        print(" ERROR: test.gsd file not found!")
        print("Please run: cp init-0.gsd test.gsd")
        return 1
    
    try:
        run_fixed_simulation()
        return 0
        
    except Exception as e:
        print(f" Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code) 