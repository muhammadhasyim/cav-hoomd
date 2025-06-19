# HOOMD-blue Bussi Reservoir Energy Plugin

This plugin extends the Bussi stochastic velocity rescaling thermostat in HOOMD-blue to track the energy dumped into or taken from the thermal reservoir during molecular dynamics simulations.

## Features

- **Reservoir Energy Tracking**: Monitor cumulative and instantaneous energy exchange with the thermal bath
- **Separate Tracking**: Track translational and rotational contributions separately
- **HOOMD Logging Integration**: All quantities are available for logging and analysis
- **Reset Functionality**: Reset counters to measure energy over specific time periods
- **Drop-in Replacement**: Uses the same interface as the standard Bussi thermostat

## Installation

### Prerequisites

- HOOMD-blue 4.x installed from source
- CMake 3.16 or newer
- C++17 compatible compiler

### Building the Plugin

1. Clone this repository:
   ```bash
   git clone <this-repository>
   cd hoomd-bussi-reservoir
   ```

2. Configure the build:
   ```bash
   cmake -B build -S .
   ```

3. Build and install:
   ```bash
   cmake --build build
   cmake --install build
   ```

## Usage

### Basic Usage

```python
import hoomd
import hoomd.bussi_reservoir as bussi_res

# Create simulation
simulation = hoomd.Simulation(device=hoomd.device.CPU(), seed=42)

# Set up your system (particles, box, etc.)
# ... simulation setup code ...

# Create the Bussi reservoir thermostat
bussi = bussi_res.thermostats.BussiReservoir(
    kT=1.5,  # Temperature
    tau=simulation.operations.integrator.dt * 20  # Time constant
)

# Use with ConstantVolume method
nve = hoomd.md.methods.ConstantVolume(
    filter=hoomd.filter.All(),
    thermostat=bussi
)

simulation.operations.integrator = hoomd.md.Integrator(
    dt=0.001,
    methods=[nve]
)

# Run simulation
simulation.run(10000)

# Check reservoir energies
print(f"Total reservoir energy: {bussi.total_reservoir_energy}")
print(f"Translational component: {bussi.reservoir_energy_translational}")
print(f"Rotational component: {bussi.reservoir_energy_rotational}")
```

### Logging Reservoir Energies

```python
import hoomd.write

# Create logger
logger = hoomd.logging.Logger()

# Add reservoir energy quantities
logger.add(bussi, quantities=[
    'total_reservoir_energy',
    'reservoir_energy_translational',
    'reservoir_energy_rotational',
    'instantaneous_reservoir_total',
    'instantaneous_reservoir_translational',
    'instantaneous_reservoir_rotational'
])

# Write to file
gsd_writer = hoomd.write.GSD(
    filename='trajectory.gsd',
    trigger=hoomd.trigger.Periodic(1000),
    logger=logger
)
simulation.operations.writers.append(gsd_writer)

# Or write to table
table_writer = hoomd.write.Table(
    output=open('reservoir_energy.log', 'w'),
    trigger=hoomd.trigger.Periodic(100),
    logger=logger
)
simulation.operations.writers.append(table_writer)
```

### Measuring Energy Over Specific Periods

```python
# Reset counters before measurement period
bussi.reset_reservoir_energy()

# Run equilibration
simulation.run(5000)

# Reset again before production
bussi.reset_reservoir_energy()

# Run production
simulation.run(50000)

# Measure reservoir energy for production period only
production_reservoir_energy = bussi.total_reservoir_energy
print(f"Reservoir energy during production: {production_reservoir_energy}")
```

## Theory

The reservoir energy represents the energy exchange between the simulated system and the thermal bath during stochastic velocity rescaling. At each timestep, the kinetic energy is rescaled by a factor α:

```
KE_new = α² × KE_old
```

The energy dumped to the reservoir is:

```
ΔE_reservoir = KE_old - KE_new = KE_old × (1 - α²)
```

- **Positive values**: Energy dumped to reservoir (system cooling)
- **Negative values**: Energy taken from reservoir (system heating)
- **Cumulative tracking**: Total energy exchange over simulation time
- **Instantaneous tracking**: Energy exchange in the last timestep

## API Reference

### BussiReservoir Class

#### Constructor
```python
BussiReservoir(kT, tau=0.0)
```

**Parameters:**
- `kT`: Temperature set point (energy units)
- `tau`: Thermostat time constant (time units), default 0.0 for instant thermalization

#### Properties (Read-only)
- `reservoir_energy_translational`: Cumulative energy from translational DOF
- `reservoir_energy_rotational`: Cumulative energy from rotational DOF  
- `total_reservoir_energy`: Total cumulative energy
- `instantaneous_reservoir_translational`: Last timestep translational energy
- `instantaneous_reservoir_rotational`: Last timestep rotational energy
- `instantaneous_reservoir_total`: Last timestep total energy

#### Methods
- `reset_reservoir_energy()`: Reset all energy counters to zero

## Testing

Run the test suite:

```bash
cd build
python -m pytest ../src/pytest/
```

## Physical Interpretation

The reservoir energy provides insights into:

1. **Thermalization efficiency**: How much energy is exchanged to maintain temperature
2. **System behavior**: Whether the system tends to heat up or cool down naturally
3. **Thermostat performance**: Monitoring for proper statistical mechanics
4. **Energy conservation**: Tracking total energy flow in the system

## Comparison with Standard Bussi Thermostat

This plugin provides the exact same thermodynamic behavior as the standard Bussi thermostat, with the additional capability to track reservoir energies. The computational overhead is minimal (just a few floating-point operations per timestep).

## Citation

If you use this plugin in your research, please cite:

- The original Bussi et al. paper: [Canonical sampling through velocity rescaling](https://doi.org/10.1063/1.2408420)
- HOOMD-blue: [HOOMD-blue: A Python package for high-performance molecular dynamics and hard particle Monte Carlo simulations](https://doi.org/10.1016/j.cpc.2022.108363)

## License

This plugin is released under the same BSD 3-Clause License as HOOMD-blue.

## Contributing

Contributions are welcome! Please submit issues and pull requests on the GitHub repository.
