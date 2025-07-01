#!/usr/bin/env python3
"""
Test script to verify that the plugin imports work correctly.
This simulates the environment that Read the Docs will have.
"""

import os
import sys
from pathlib import Path
from unittest.mock import MagicMock
import types

def test_plugin_imports():
    """Test that the plugin modules can be imported correctly."""
    
    # Simulate Read the Docs environment
    os.environ['READTHEDOCS'] = 'True'
    
    # Get the project root
    project_root = Path(__file__).parent
    sys.path.insert(0, str(project_root / 'src'))
    
    print("Setting up plugin import test...")
    
    # Mock the C++ extensions
    sys.modules['hoomd.cavitymd._cavitymd'] = MagicMock()
    sys.modules['hoomd.bussi_reservoir._bussi_reservoir'] = MagicMock()
    print("✅ Mocked C++ extensions")
    
    # Verify plugin directories exist
    cavitymd_path = project_root / 'src' / 'cavitymd'
    bussi_path = project_root / 'src' / 'bussi_reservoir'
    
    print(f"Cavity MD path exists: {cavitymd_path.exists()}")
    print(f"Bussi Reservoir path exists: {bussi_path.exists()}")
    
    # Try to establish HOOMD namespace (will fail, so create minimal one)
    try:
        import hoomd
        print("✅ HOOMD base package available")
    except ImportError:
        print("Creating minimal HOOMD namespace...")
        hoomd = types.ModuleType('hoomd')
        sys.modules['hoomd'] = hoomd
    
    # Try to import the plugin modules
    try:
        import hoomd.cavitymd
        import hoomd.bussi_reservoir
        
        print("✅ Plugin modules imported successfully!")
        
        # Test specific classes
        if hasattr(hoomd.cavitymd, 'CavityForce'):
            print("  ✅ CavityForce found")
        else:
            print("  ❌ CavityForce not found")
            
        if hasattr(hoomd.bussi_reservoir, 'BussiReservoir'):
            print("  ✅ BussiReservoir found")
        else:
            print("  ❌ BussiReservoir not found")
            
        return True
        
    except Exception as e:
        print(f"❌ Failed to import plugin modules: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_plugin_imports()
    if success:
        print("\n🎉 Plugin import test PASSED!")
        print("The Read the Docs build should now succeed.")
    else:
        print("\n💥 Plugin import test FAILED!")
        print("Additional fixes may be needed.")
    
    sys.exit(0 if success else 1) 