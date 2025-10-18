# Molecular Temperature Decomposition for Diatomic Molecules

## Overview

The `DiatomicMolecularTemperatures` class provides kinetic temperature estimates for different degrees of freedom in a system of diatomic molecules (A₂ and B₂ type).

## Temperature Components

For diatomic molecules, the total kinetic energy can be decomposed into three components:

### 1. Translational Temperature (T_trans)
- **Definition**: Temperature associated with center-of-mass motion
- **Theory**: Each molecule has 3 translational degrees of freedom
- **Relation**: `<KE_trans> = (3/2) * N * k_B * T_trans`
- **Formula**: `T_trans = (2 * KE_trans) / (3 * N * k_B)`

Where:
- `KE_trans = (1/2) * M * v_cm²`
- `M = m₁ + m₂` (total mass of dimer)
- `v_cm = (m₁*v₁ + m₂*v₂) / M` (center of mass velocity)

### 2. Rotational Temperature (T_rot)
- **Definition**: Temperature associated with rotation around center of mass
- **Theory**: Each molecule has 2 rotational degrees of freedom in 3D
- **Relation**: `<KE_rot> = N * k_B * T_rot`
- **Formula**: `T_rot = KE_rot / (N * k_B)`

Where:
- `KE_rot = (1/2) * μ * v_perp²`
- `μ = m₁*m₂/(m₁+m₂)` (reduced mass)
- `v_perp = v_rel - v_bond * r̂` (velocity perpendicular to bond)
- `v_rel = v₂ - v₁` (relative velocity)

### 3. Vibrational Temperature (T_vib)
- **Definition**: Temperature associated with motion along the bond direction
- **Theory**: Each molecule has 1 vibrational degree of freedom
- **Relation**: `<KE_vib> = (1/2) * N * k_B * T_vib`
- **Formula**: `T_vib = (2 * KE_vib) / (N * k_B)`

Where:
- `KE_vib = (1/2) * μ * v_bond²`
- `v_bond = v_rel · r̂` (relative velocity component along bond)
- `r̂ = (r₂ - r₁) / |r₂ - r₁|` (unit vector along bond)

## Usage

### Basic Integration

```python
from cavitymd.molecular_temperatures import DiatomicMolecularTemperatures
from cavitymd.analysis import ElapsedTimeTracker

# In your simulation setup:
# ... (after creating simulation and time_tracker)

# Create molecular temperature tracker
mol_temp_tracker = DiatomicMolecularTemperatures(
    simulation=sim,
    time_tracker=time_tracker,
    output_period_ps=1.0,  # Output every 1 ps
    output_file="molecular_temperatures.csv",
    debug=False  # Set to True for detailed output
)

# Add to simulation
sim.operations.writers.append(mol_temp_tracker)

# Run simulation
sim.run(N_steps)

# Access current temperatures
print(f"Translational: {mol_temp_tracker.translational_temp:.2f} K")
print(f"Rotational: {mol_temp_tracker.rotational_temp:.2f} K")
print(f"Vibrational: {mol_temp_tracker.vibrational_temp:.2f} K")
```

### Integration with Existing Simulation Scripts

To add this to your `18_unified_cavity_dynamics.py` or similar scripts:

```python
# Add import at the top
from cavitymd.molecular_temperatures import DiatomicMolecularTemperatures

# In run_single_experiment function, after creating time_tracker:
if args.track_molecular_temps:
    mol_temp_tracker = DiatomicMolecularTemperatures(
        simulation=sim,
        time_tracker=time_tracker,
        output_period_ps=args.temp_output_period,
        output_file=os.path.join(output_dir, "molecular_temperatures.csv"),
        debug=args.debug
    )
    sim.operations.writers.append(mol_temp_tracker)
```

## Output Format

The tracker outputs a CSV file with the following columns:

| Column | Description |
|--------|-------------|
| `time_ps` | Simulation time in picoseconds |
| `T_trans` | Overall translational temperature (K) |
| `T_rot` | Overall rotational temperature (K) |
| `T_vib` | Overall vibrational temperature (K) |
| `T_kinetic_total` | Total kinetic temperature from all velocities (K) |
| `T_trans_O2` | Translational temperature for O-O molecules only (K) |
| `T_trans_N2` | Translational temperature for N-N molecules only (K) |
| `T_rot_O2` | Rotational temperature for O-O molecules only (K) |
| `T_rot_N2` | Rotational temperature for N-N molecules only (K) |
| `T_vib_O2` | Vibrational temperature for O-O molecules only (K) |
| `T_vib_N2` | Vibrational temperature for N-N molecules only (K) |

## Physical Interpretation

### Equilibrium Expectations
At thermal equilibrium at temperature T, you should expect:
- `T_trans ≈ T`
- `T_rot ≈ T`
- `T_vib ≈ T` (classically; quantum effects may cause deviations)

### Non-Equilibrium Signatures

**Cooling/Heating Dynamics**:
- Different components may equilibrate at different rates
- Typically: `T_trans` equilibrates fastest, `T_vib` slowest

**Cavity Coupling Effects**:
- Selective coupling may preferentially affect certain DOF
- Monitor deviations from equipartition: `T_trans ≈ T_rot ≈ T_vib`

**Vibrational Heating**:
- If `T_vib > T_trans, T_rot`: Energy is preferentially deposited into bond vibrations
- Common in systems with strong intramolecular forces

## Energy Conservation Check

The total kinetic energy should equal the sum of components:
```
KE_total = KE_trans + KE_rot + KE_vib
```

This is verified internally and can be checked via:
```
3*N*T_total ≈ 3*N*T_trans + 2*N*T_rot + N*T_vib
```
Where N is the number of molecules.

## Implementation Details

### Molecular System Requirements
- System must contain diatomic molecules (A₂ and B₂ type)
- Bonds must be defined connecting pairs of atoms
- Currently supports:
  - O-O dimers (mass = 29150 each)
  - N-N dimers (mass = 25527 each)
  - Mixed O-N dimers

### Periodic Boundary Conditions
- Bond vectors are properly handled across periodic boundaries
- Uses minimum image convention

### Performance
- Computational cost: O(N_molecules) per timestep
- Minimal overhead compared to force calculations
- All calculations use NumPy for efficiency

## Example Analysis Script

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load data
data = pd.read_csv('molecular_temperatures.csv', comment='#')

# Plot temperature evolution
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(data['time_ps'], data['T_trans'], label='Translational')
ax.plot(data['time_ps'], data['T_rot'], label='Rotational')
ax.plot(data['time_ps'], data['T_vib'], label='Vibrational')
ax.plot(data['time_ps'], data['T_kinetic_total'], 
        label='Total Kinetic', ls='--', alpha=0.7)

ax.set_xlabel('Time (ps)')
ax.set_ylabel('Temperature (K)')
ax.legend()
ax.grid(True, alpha=0.3)
plt.savefig('temperature_evolution.pdf')
plt.show()

# Check equipartition
T_avg = data[['T_trans', 'T_rot', 'T_vib']].mean()
print("Average Temperatures:")
print(T_avg)
print("\nEquipartition Check:")
print(f"Std Dev / Mean: {T_avg.std() / T_avg.mean() * 100:.2f}%")
```

## References

- Frenkel, D., & Smit, B. (2002). *Understanding Molecular Simulation*. Academic Press.
- McQuarrie, D. A. (2000). *Statistical Mechanics*. University Science Books.

## Troubleshooting

**Issue**: Temperatures are all zero
- **Solution**: Check that velocities are initialized and simulation has run

**Issue**: Temperatures are negative or unphysical
- **Solution**: Check that bond topology is correct; enable debug mode to see detailed output

**Issue**: Large discrepancy between components
- **Solution**: This may be physical (non-equilibrium); verify with longer equilibration

**Issue**: `T_kinetic_total` doesn't match sum of components
- **Solution**: Check that `3*T_trans + 2*T_rot + T_vib ≈ 6*T_total` (accounting for DOF)

