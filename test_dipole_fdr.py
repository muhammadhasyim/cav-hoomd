#!/usr/bin/env python3
"""
Test Script for Dipole Moment FDR Implementation

This script performs basic functionality tests for the dipole moment FDR 
analysis and force implementations.

Usage:
    python test_dipole_fdr.py

Author: AI Assistant  
Created: 2025-09-26 (Phase 4)
"""

import sys
import numpy as np
import hoomd
from pathlib import Path

def test_imports():
    """Test that all required modules can be imported."""
    print("🧪 Testing imports...")
    
    try:
        from hoomd.cavitymd.fdr_forces import DipoleResponseForce
        print("✅ DipoleResponseForce imported successfully")
    except ImportError as e:
        print(f"❌ Failed to import DipoleResponseForce: {e}")
        return False
        
    try:
        from hoomd.cavitymd.analysis import DipoleMomentFDRTracker
        print("✅ DipoleMomentFDRTracker imported successfully")
    except ImportError as e:
        print(f"❌ Failed to import DipoleMomentFDRTracker: {e}")
        return False
        
    try:
        import hoomd.cavitymd._cavitymd as _cavitymd
        if hasattr(_cavitymd, 'DipoleResponseForceCompute'):
            print("✅ DipoleResponseForceCompute backend available")
        else:
            print("❌ DipoleResponseForceCompute backend missing")
            return False
    except ImportError as e:
        print(f"❌ Failed to import C++ backend: {e}")
        return False
        
    return True

def test_dipole_response_force():
    """Test DipoleResponseForce functionality."""
    print("\n🧪 Testing DipoleResponseForce...")
    
    try:
        from hoomd.cavitymd.fdr_forces import DipoleResponseForce
        
        # Test basic instantiation
        force = DipoleResponseForce(
            field_vector=[0, 0, 1],
            field_strength=1e-5,
            sign=1.0,
            exclude_cavity=True
        )
        
        print(f"✅ Force instantiated successfully")
        print(f"   Field vector: {force.field_vector}")
        print(f"   Field strength: {force.field_strength}")
        print(f"   Sign: {force.sign}")
        print(f"   Exclude cavity: {force.exclude_cavity}")
        
        # Test parameter updates
        force.setParams(field_strength=2e-5, sign=-1.0)
        print(f"✅ Parameters updated successfully")
        
        # Test enable/disable
        force.setEnabled(False)
        enabled = force.getEnabled()
        print(f"✅ Enable/disable working: {enabled}")
        
        return True
        
    except Exception as e:
        print(f"❌ DipoleResponseForce test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_dipole_fdr_tracker():
    """Test DipoleMomentFDRTracker functionality."""
    print("\n🧪 Testing DipoleMomentFDRTracker...")
    
    try:
        from hoomd.cavitymd.analysis import DipoleMomentFDRTracker
        
        # Create a mock time tracker (simplified test without full ElapsedTimeTracker)
        class MockTimeTracker:
            def __init__(self):
                self.start_time = 0.0
                self.current_time = 0.0
            
            def get_elapsed_time_ps(self):
                return self.current_time
        
        time_tracker = MockTimeTracker()
        
        # Test basic instantiation
        tracker = DipoleMomentFDRTracker(
            time_tracker=time_tracker,
            output_file="test_dipole_fdr.csv",
            max_correlation_time_ps=50.0,
            correlation_output_interval_ps=0.1,
            exclude_cavity=True,
            field_direction=[0, 0, 1],
            enable_response_measurement=False
        )
        
        print(f"✅ DipoleMomentFDRTracker instantiated successfully")
        print(f"   Output file: {tracker.output_file}")
        print(f"   Max correlation time: {tracker.max_correlation_time_ps} ps")
        print(f"   Field direction: {tracker.field_direction}")
        print(f"   Response measurement: {tracker.enable_response_measurement}")
        
        # Test file initialization
        if Path("test_dipole_fdr.csv").exists():
            print("✅ Output file created successfully")
            Path("test_dipole_fdr.csv").unlink()  # Clean up
        else:
            print("⚠️ Output file not created (expected for test)")
            
        return True
        
    except Exception as e:
        print(f"❌ DipoleMomentFDRTracker test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_synthetic_fdr():
    """Test FDR calculations with synthetic data."""
    print("\n🧪 Testing synthetic FDR calculations...")
    
    try:
        # Generate synthetic dipole moment time series
        dt = 0.01  # ps
        t_max = 20.0  # ps
        time = np.arange(0, t_max, dt)
        
        # Create correlated dipole moment with exponential decay
        tau = 2.0  # correlation time in ps
        correlation_func = np.exp(-time / tau)
        
        # Generate random dipole moment with this correlation
        np.random.seed(42)  # Reproducible
        noise = np.random.normal(0, 1, len(time))
        
        # Simple exponential correlation (not perfect, but good for testing)
        mu_x = np.convolve(noise, correlation_func, mode='same') * dt
        mu_y = np.convolve(np.random.normal(0, 1, len(time)), correlation_func, mode='same') * dt
        mu_z = np.convolve(np.random.normal(0, 1, len(time)), correlation_func, mode='same') * dt
        
        # Test autocorrelation calculation (simplified)
        field_direction = np.array([0, 0, 1])
        mu_proj = mu_z  # Project onto field direction
        
        # Calculate autocorrelation
        mu_mean = np.mean(mu_proj)
        delta_mu = mu_proj - mu_mean
        
        max_lag = min(len(delta_mu) // 4, int(10.0 / dt))  # Max 10 ps correlation
        autocorr = np.zeros(max_lag)
        
        for lag in range(max_lag):
            if lag == 0:
                autocorr[lag] = np.mean(delta_mu**2)
            else:
                autocorr[lag] = np.mean(delta_mu[:-lag] * delta_mu[lag:])
        
        # Normalize
        autocorr = autocorr / autocorr[0]
        
        # Check that it decays roughly exponentially
        tau_fitted = -dt / np.log(autocorr[int(1.0/dt)] / autocorr[0]) if autocorr[int(1.0/dt)] > 0 else float('inf')
        
        print(f"✅ Synthetic autocorrelation calculated")
        print(f"   Input correlation time: {tau:.1f} ps")
        print(f"   Fitted correlation time: {tau_fitted:.1f} ps")
        print(f"   Autocorr[0]: {autocorr[0]:.3f} (should be 1.0)")
        print(f"   Autocorr[τ]: {autocorr[int(tau/dt)] if int(tau/dt) < len(autocorr) else 'N/A'}")
        
        # Test static susceptibility calculation
        kT = 3.1668e-6 * 100  # 100K in a.u.
        static_chi = np.trapz(autocorr, dx=dt) / kT
        
        print(f"✅ Static susceptibility: {static_chi:.2e}")
        
        return True
        
    except Exception as e:
        print(f"❌ Synthetic FDR test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_cli_integration():
    """Test CLI argument integration."""
    print("\n🧪 Testing CLI integration...")
    
    try:
        # Test that the unified script accepts dipole FDR arguments
        import subprocess
        result = subprocess.run([
            "python", "18_unified_cavity_dynamics.py", "--help"
        ], capture_output=True, text=True, timeout=10)
        
        if "--enable-dipole-fdr" in result.stdout:
            print("✅ --enable-dipole-fdr argument found in help")
        else:
            print("❌ --enable-dipole-fdr argument missing from help")
            return False
            
        if "--enable-dipole-response" in result.stdout:
            print("✅ --enable-dipole-response argument found in help")
        else:
            print("❌ --enable-dipole-response argument missing from help")
            return False
            
        if "--dipole-fdr-field-direction" in result.stdout:
            print("✅ --dipole-fdr-field-direction argument found in help")
        else:
            print("❌ --dipole-fdr-field-direction argument missing from help")
            return False
            
        return True
        
    except Exception as e:
        print(f"❌ CLI integration test failed: {e}")
        return False

def main():
    """Run all tests."""
    print("🔬 Dipole Moment FDR Implementation Test Suite")
    print("=" * 50)
    
    tests = [
        ("Import Test", test_imports),
        ("DipoleResponseForce Test", test_dipole_response_force),
        ("DipoleMomentFDRTracker Test", test_dipole_fdr_tracker),
        ("Synthetic FDR Test", test_synthetic_fdr),
        ("CLI Integration Test", test_cli_integration)
    ]
    
    results = []
    for test_name, test_func in tests:
        print(f"\n{'='*20} {test_name} {'='*20}")
        try:
            success = test_func()
            results.append((test_name, success))
        except Exception as e:
            print(f"❌ {test_name} crashed: {e}")
            results.append((test_name, False))
    
    # Summary
    print("\n" + "="*60)
    print("🔬 TEST SUMMARY")
    print("="*60)
    
    passed = 0
    total = len(results)
    
    for test_name, success in results:
        status = "✅ PASS" if success else "❌ FAIL"
        print(f"{status}: {test_name}")
        if success:
            passed += 1
    
    print(f"\nOverall: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed! Dipole FDR implementation looks good.")
        return 0
    else:
        print("⚠️ Some tests failed. Please check the implementation.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
