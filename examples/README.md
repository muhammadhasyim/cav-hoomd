# CavityMD Examples

This directory contains 4 focused examples demonstrating the key capabilities of the CavityMD plugin for HOOMD-blue.

## Examples Overview

### 1. Basic Cavity Simulation (`01_basic_cavity_simulation.py`)
**Purpose:** Demonstrates the fundamental cavity MD simulation setup
**Duration:** 5 ps
**System:** 3 water molecules + 1 cavity particle (10 particles total)
**Features:**
- Basic cavity-molecule coupling mechanics
- Simple molecular system with realistic charges
- NVT thermostat for molecules, NVE for cavity
- No observable tracking - pure simulation
- Shows energy component evolution

**Usage:**
```bash
python 01_basic_cavity_simulation.py
```

### 2. Force Validation (`02_force_validation.py`)
**Purpose:** Validates consistency between Python, C++ CPU, and CUDA GPU implementations
**Systems:** Simple (3 particles) and complex (13 particles)
**Features:**
- Tests all three cavity force implementations
- Compares forces on each particle with high precision
- Compares energy components between implementations
- Creates detailed comparison plots
- Validates numerical accuracy within tolerance (1e-10)

**Usage:**
```bash
python 02_force_validation.py
```

**Output:**
- `force_validation_comparison.png` - Force and energy comparison plots
- Console output with detailed numerical comparisons

### 3. Energy Tracking (`03_energy_tracking_1ps.py`)
**Purpose:** Demonstrates comprehensive energy tracking over 1 ps
**Duration:** 1 ps (2000 steps at 0.5 fs)
**System:** 6 water molecules + 1 cavity particle (19 particles total)
**Features:**
- Uses `EnergyContributionTracker` for detailed energy logging
- Tracks cavity energy components (harmonic, coupling, dipole self-energy)
- Monitors molecular kinetic/potential energies
- Records thermostat reservoir energies
- Creates real-time energy analysis plots

**Usage:**
```bash
python 03_energy_tracking_1ps.py
```

**Output:**
- `energy_tracking_1ps.txt` - Detailed energy time series data
- `energy_tracking_1ps_analysis.png` - Energy evolution plots

### 4. Structure Tracking (`04_structure_tracking_10ps.py`)
**Purpose:** Advanced structure and correlation analysis over 10 ps
**Duration:** 10 ps (20000 steps at 0.5 fs)
**System:** 12 water molecules + 1 cavity particle (37 particles total)
**Features:**
- Uses `DensityCorrelationTracker` for F(k,t) calculation at multiple k-values
- Uses `DipoleAutocorrelation` for dipole dynamics
- Uses `CavityModeTracker` for cavity position/velocity tracking
- Long-time simulation to capture correlation timescales
- Comprehensive structure analysis plots

**Usage:**
```bash
python 04_structure_tracking_10ps.py
```

**Output:**
- `structure_k_*.txt` - F(k,t) data for different k-values (0.5, 1.0, 1.5, 2.0 Å⁻¹)
- `dipole_autocorr_10ps.txt` - Dipole autocorrelation function
- `cavity_mode_10ps.txt` - Cavity mode dynamics
- `structure_tracking_10ps_analysis.png` - Comprehensive analysis plots

## System Requirements

- HOOMD-blue with CavityMD plugin compiled
- Python packages: `numpy`, `matplotlib`, `gsd`
- For GPU validation: CUDA-capable GPU
- For C++ validation: Compiled C++ cavity force module

## Running Examples

Each example is self-contained and can be run independently:

```bash
cd examples
python 01_basic_cavity_simulation.py
python 02_force_validation.py
python 03_energy_tracking_1ps.py
python 04_structure_tracking_10ps.py
```

## Understanding the Physics

### Cavity Force Hamiltonian
The examples demonstrate the cavity-molecule interaction based on:
```
H = (1/2) * K * q² + g * q · d + (g²/2K) * d²
```
Where:
- `q` = cavity mode position
- `d` = molecular dipole moment  
- `g` = coupling strength
- `K = phmass * omegac²` = spring constant

### Key Parameters
- **`kvector`**: Cavity mode wave vector (currently [1,0,0])
- **`couplstr`** (g): Coupling strength between cavity and dipole
- **`omegac`**: Cavity frequency in atomic units (Hartree)
- **`phmass`**: Photon mass determining spring constant

### Energy Components
1. **Harmonic energy**: `(1/2) * K * q²` - cavity oscillator energy
2. **Coupling energy**: `g * (q · d)` - direct cavity-dipole interaction
3. **Dipole self-energy**: `(g²/2K) * d²` - self-interaction correction

## Analysis Capabilities

The examples showcase the refactored analysis framework:

- **`EnergyContributionTracker`**: Comprehensive energy logging
- **`DensityCorrelationTracker`**: Intermediate scattering function F(k,t)
- **`DipoleAutocorrelation`**: Dipole moment correlation functions
- **`CavityModeTracker`**: Cavity position and velocity dynamics

## Troubleshooting

**Force validation fails:**
- Check that C++ module is properly compiled
- Verify GPU drivers and CUDA installation
- Try reducing tolerance if small numerical differences occur

**Simulation crashes:**
- Reduce timestep (`dt`) for stability
- Check particle overlaps in initial configuration
- Verify molecular force parameters

**Analysis files not created:**
- Check write permissions in directory
- Verify tracker periods are set correctly
- Ensure simulation runs long enough for data collection

## Next Steps

After running these examples:
1. Modify cavity parameters (`couplstr`, `omegac`) to see different physics
2. Change system size and composition
3. Experiment with different thermostats and ensembles
4. Use the analysis data for your own scientific studies

For production simulations, see the `run_cavity_experiments.py` script for automated experimental workflows. 