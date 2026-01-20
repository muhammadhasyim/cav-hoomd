# OpenMM CavityMD Plugin

Cavity-molecule coupling forces for OpenMM molecular dynamics simulations.

## Overview

This plugin implements cavity molecular dynamics (cavity-MD) for OpenMM, enabling
simulations of molecules coupled to optical cavity modes. The implementation includes:

- **CavityForce**: Single-mode cavity-molecule interaction force
- **Bussi Thermostat**: Stochastic velocity rescaling thermostat with reservoir energy tracking
- **RPMD Support**: Ring polymer molecular dynamics with centroid cavity coupling
- **Time-Varying Coupling**: Support for dynamic coupling protocols

## Physics

The cavity Hamiltonian is:

```
H = (1/2) K q² + g q·d + (g²/2K) d²
```

where:
- `q` = cavity mode coordinate
- `d` = molecular dipole moment (Σᵢ qᵢ rᵢ)
- `g` = coupling strength (λ × ωc)
- `K` = spring constant (m_ph × ω²)

## Requirements

- OpenMM 8.0 or later
- Python 3.8+
- CUDA toolkit (optional, for GPU acceleration)
- CMake 3.17+

## Building

```bash
# From the cav-hoomd root directory
./build_openmm.sh

# Or manually:
cd src/openmm-cavitymd
mkdir build && cd build
cmake .. -DOPENMM_DIR=/path/to/openmm
make -j$(nproc)
make install
make PythonInstall
```

## Usage

```python
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from openmm.cavitymd import CavityForce

# Create system
pdb = app.PDBFile('system.pdb')
forcefield = app.ForceField('amber14-all.xml')
system = forcefield.createSystem(pdb.topology)

# Add cavity particle
cavity_index = system.addParticle(1.0 * unit.dalton)

# Add cavity force
cavity_force = CavityForce(
    omegac=0.00913,  # 2000 cm^-1 in a.u.
    lambda_coupling=0.001,
    phmass=1.0
)
cavity_force.setCavityParticleIndex(cavity_index)
system.addForce(cavity_force)

# Create integrator and run
integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
simulation = app.Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.step(10000)
```

## Testing

```bash
cd tests
pytest -v
```

## License

BSD 3-Clause License

## Citation

If you use this software, please cite:

```bibtex
@software{openmm_cavitymd_2025,
  title={OpenMM CavityMD: Cavity Molecular Dynamics Plugin},
  author={Development Team},
  year={2025}
}
```
