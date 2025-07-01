# Cavity HOOMD - Molecular Dynamics with Optical Cavity Coupling

[![License](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](LICENSE)

A HOOMD-blue plugin for cavity-coupled molecular dynamics simulations. Study how molecular systems interact with optical cavity modes using high-performance MD simulations.

## Quick Start

### Installation

```bash
# Build and install
cmake -B build -S .
cmake --build build
cmake --install build
```

### Run Your First Simulation

```bash
# Basic cavity-coupled simulation
python examples/05_advanced_run.py --coupling 1e-3 --runtime 1000

# Control simulation (no cavity)
python examples/05_advanced_run.py --no-cavity --runtime 1000
```

That's it! The simulation will create output files with trajectory data, energy tracking, and cavity mode information.

## What It Does

The plugin adds cavity-molecule coupling to HOOMD-blue simulations using the interaction:

```
H = (1/2) * K * q² + g * q · d + (g²/2K) * d²
```

Where:
- `q` is the cavity mode position  
- `d` is the molecular dipole moment
- `g` is the coupling strength
- `K` is the cavity spring constant

## Example Usage

The main interface is the `examples/05_advanced_run.py` script:

```bash
# Run with specific parameters
python examples/05_advanced_run.py \
    --coupling 1e-3 \
    --temperature 100 \
    --frequency 2000 \
    --runtime 1000

# Parameter sweep over coupling strengths  
python examples/05_advanced_run.py \
    --coupling 1e-3,1e-4,1e-5 \
    --runtime 500

# Multiple independent replicas
python examples/05_advanced_run.py \
    --coupling 1e-3 \
    --replicas 1-5 \
    --runtime 1000
```

### Available Options

**Basic Parameters:**
- `--coupling` - Coupling strength (e.g., 1e-3)
- `--temperature` - Temperature in K (default: 100)
- `--frequency` - Cavity frequency in cm⁻¹ (default: 2000)  
- `--runtime` - Simulation time in ps (default: 500)
- `--no-cavity` - Run without cavity (control simulation)

**Advanced:**
- `--replicas` - Run multiple replicas (e.g., "1-5" or "1,3,5")
- `--device GPU` - Use GPU acceleration
- `--enable-energy-tracker` - Detailed energy tracking
- `--enable-fkt` - F(k,t) correlation analysis

## Output Files

Each simulation creates:
- `prod-X.gsd` - Trajectory file
- `prod-X-energy.txt` - Energy tracking data
- `prod-X-cavity_mode.txt` - Cavity mode properties

## Jupyter Notebook

See `examples/05_advanced_run.ipynb` for an interactive example that shows the same simulation setup and analysis.

## Python API

For programmatic access:

```python
import hoomd
from hoomd.cavitymd import CavityForce
from hoomd.bussi_reservoir import BussiReservoir

# Create simulation
sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=42)
sim.create_state_from_gsd('system.gsd')

# Add cavity force
cavity_force = CavityForce(
    kvector=[0, 0, 1],
    couplstr=0.001,
    omegac=0.01,
    phmass=1.0
)

# Setup integrator with thermostat
integrator = hoomd.md.Integrator(dt=0.001)
integrator.forces.append(cavity_force)

# Add thermostats for molecules and cavity
molecular_thermostat = BussiReservoir(kT=0.01, tau=1.0)
molecular_filter = hoomd.filter.Type(['O', 'N'])
integrator.methods.append(
    hoomd.md.methods.ConstantVolume(
        filter=molecular_filter, 
        thermostat=molecular_thermostat
    )
)

sim.operations.integrator = integrator
sim.run(10000)
```

## Installation Options

**With GPU support (default):**
```bash
cmake -B build -S . -DENABLE_GPU=ON
cmake --build build
cmake --install build
```

**CPU only:**
```bash
cmake -B build -S . -DENABLE_GPU=OFF  
cmake --build build
cmake --install build
```

## Requirements

- HOOMD-blue 4.0+
- Python 3.8+
- NumPy
- For GPU support: CUDA toolkit

## Citation

If you use this software in your research, please cite:

```bibtex
@software{cavity_hoomd_2025,
  title={Cavity HOOMD: Molecular Dynamics with Optical Cavity Coupling},
  author={Development Team},
  year={2025},
  url={https://github.com/yourusername/cavity-hoomd}
}
```

## License

BSD 3-Clause License - see [LICENSE](LICENSE) file.
