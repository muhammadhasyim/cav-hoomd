# Module Structure

This repository has been reorganized into two separate HOOMD-blue plugins:

## hoomd.bussi_reservoir

**Purpose**: Bussi thermostat implementation with reservoir energy tracking

**Location**: `src/hoomd/bussi_reservoir/`

**Components**:
- `thermostats.py` - Python interface for BussiReservoir thermostat
- `module.cc` - C++ bindings for thermostat
- `BussiReservoir.cc/.h` - C++ implementation

**Usage**:
```python
from hoomd.bussi_reservoir import BussiReservoir

thermostat = BussiReservoir(kT=0.01, tau=1.0)
```

## hoomd.cavitymd

**Purpose**: High-performance cavity molecular dynamics force implementation

**Location**: `src/hoomd/cavitymd/`

**Components**:
- `forces.py` - Unified Python interface with automatic C++/Python fallback
- `cavity_force_python.py` - Pure Python fallback implementation
- `module.cc` - C++ bindings for cavity force
- `CavityForceCompute.cc/.h` - C++ CPU implementation
- `CavityForceComputeGPU.cc/.h/.cu/.cuh` - CUDA GPU implementation

**Usage**:
```python
from hoomd.cavitymd import CavityForce

cavity_force = CavityForce(
    kvector=np.array([0, 0, 1]),
    couplstr=0.001,
    omegac=0.01,
    phmass=1.0
)
```

## Build System

The CMakeLists.txt now builds two separate libraries:
- `_bussi_reservoir.so` - Installed to `hoomd/bussi_reservoir/`
- `_cavitymd.so` - Installed to `hoomd/cavitymd/`

Each module is self-contained and can be imported independently.

## Migration Guide

### Old imports:
```python
from hoomd.bussi_reservoir import BussiReservoir, CavityForce
```

### New imports:
```python
from hoomd.bussi_reservoir import BussiReservoir
from hoomd.cavitymd import CavityForce
```

This change allows for cleaner module organization and follows Python packaging best practices. 