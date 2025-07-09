# Time-Varying Coupling Strength in Cavity MD

The cavity molecular dynamics implementation now supports time-varying coupling strength using HOOMD-blue's `hoomd.variant` interface. This allows you to study how changing light-matter coupling affects molecular dynamics over the course of a simulation.

## Overview

Instead of using a constant coupling strength `g`, you can now provide any `hoomd.variant.variant_like` object that defines how the coupling strength changes over time. This is particularly useful for:

- **Adiabatic switching** of cavity coupling
- **Periodic modulation** of coupling strength
- **Ramping experiments** to study coupling threshold effects
- **Control pulse sequences** in cavity QED

## Basic Usage

### Constant Coupling (Traditional)

```python
from cavitymd.forces import CavityForce
import hoomd.variant

# Traditional constant coupling
cavity_force = CavityForce(
    kvector=[0, 0, 1],
    couplstr=1e-3,  # Constant value
    omegac=0.01,
    phmass=1.0
)

# Equivalent using variant (explicit)
cavity_force = CavityForce(
    kvector=[0, 0, 1],
    couplstr=hoomd.variant.Constant(1e-3),
    omegac=0.01,
    phmass=1.0
)
```

### Time-Varying Coupling

```python
import hoomd.variant

# Linear ramp from g=0 to g=1e-3 over 10,000 steps
ramp_coupling = hoomd.variant.Ramp(
    A=0.0,           # Initial coupling
    B=1e-3,          # Final coupling  
    t_start=0,       # Start time
    t_ramp=10000     # Ramp duration
)

cavity_force = CavityForce(
    kvector=[0, 0, 1],
    couplstr=ramp_coupling,  # Time-varying coupling
    omegac=0.01,
    phmass=1.0
)
```

## Available Variant Types

### 1. Constant
Maintains a fixed value throughout the simulation.

```python
constant = hoomd.variant.Constant(value=1e-3)
```

### 2. Ramp
Linear interpolation between two values over a specified time period.

```python
ramp = hoomd.variant.Ramp(
    A=0.0,           # Start value
    B=2e-3,          # End value
    t_start=1000,    # When to start ramping
    t_ramp=5000      # Duration of ramp
)
```

**Timeline:**
- Steps 0-999: g = 0.0
- Steps 1000-5999: g linearly increases from 0.0 to 2e-3
- Steps 6000+: g = 2e-3

### 3. Cycle
Oscillates between two values with configurable hold and transition times.

```python
cycle = hoomd.variant.Cycle(
    A=5e-4,          # First value
    B=1.5e-3,        # Second value
    t_start=0,       # When to start cycling
    t_A=1000,        # Hold time at A
    t_AB=500,        # Transition time A→B
    t_B=1000,        # Hold time at B
    t_BA=500         # Transition time B→A
)
```

**Timeline repeats:** A (1000 steps) → A→B (500 steps) → B (1000 steps) → B→A (500 steps) → ...

### 4. Power
Non-linear approach following a power law.

```python
power = hoomd.variant.Power(
    A=1e-4,          # Start value
    B=1e-3,          # End value
    power=0.5,       # Power (0.5 = square root)
    t_start=0,       # Start time
    t_ramp=8000      # Duration
)
```

The coupling follows: g(t) = A + (B-A) × (t/t_ramp)^power

## Physics Considerations

### Adiabatic vs. Sudden Changes

For physically meaningful results, consider the timescales in your system:

```python
# Adiabatic switching (slow compared to molecular motion)
# Typical molecular vibrations ~ few fs, so ramp over ps timescales
adiabatic_ramp = hoomd.variant.Ramp(
    A=0.0, B=1e-3,
    t_start=0, t_ramp=100000  # Slow ramp
)

# Sudden switching (fast compared to molecular motion)
sudden_switch = hoomd.variant.Ramp(
    A=0.0, B=1e-3, 
    t_start=10000, t_ramp=10  # Very fast ramp
)
```

### Energy Conservation

Time-varying coupling adds/removes energy from the system. Monitor the total energy:

```python
# Add energy logging
logger = hoomd.logging.Logger()
logger.add(cavity_force, quantities=[
    'harmonic_energy', 
    'coupling_energy', 
    'dipole_self_energy',
    'total_cavity_energy'
])

# Custom logger for coupling strength
class CouplingLogger:
    def __init__(self, variant):
        self.variant = variant
        
    @property 
    def coupling_strength(self):
        return self.variant(sim.timestep)

coupling_logger = CouplingLogger(ramp_coupling)
logger.add(coupling_logger, quantities=['coupling_strength'])
```

## Advanced Examples

### Control Pulse Sequence

```python
# Square pulse: 0 → g → 0
pulse = hoomd.variant.Cycle(
    A=0.0,           # Off
    B=1e-3,          # On
    t_start=5000,    # Start pulse at step 5000
    t_A=0,           # No hold at 0
    t_AB=100,        # Fast turn-on
    t_B=2000,        # Pulse duration
    t_BA=100         # Fast turn-off
)
# After one cycle, stays at A=0
```

### Exponential-like Approach

```python
# Approximate exponential using power function
exponential_like = hoomd.variant.Power(
    A=0.0,
    B=1e-3, 
    power=2.0,       # Quadratic approach
    t_start=0,
    t_ramp=10000
)
```

### Multi-Stage Protocol

```python
# Combine multiple variants using custom class
class MultiStageVariant(hoomd.variant.Variant):
    def __init__(self):
        super().__init__()
        
    def __call__(self, timestep):
        if timestep < 5000:
            # Stage 1: Linear ramp up
            return min(timestep / 5000.0 * 1e-3, 1e-3)
        elif timestep < 15000:
            # Stage 2: Constant
            return 1e-3
        else:
            # Stage 3: Exponential decay
            t_decay = timestep - 15000
            return 1e-3 * np.exp(-t_decay / 2000.0)
    
    def _min(self):
        return 0.0
        
    def _max(self): 
        return 1e-3

multi_stage = MultiStageVariant()
```

## Example Script

A complete example demonstrating time-varying coupling is provided in:
```
examples/time_varying_coupling_example.py
```

Run with different variants:
```bash
# Constant coupling
python time_varying_coupling_example.py --variant constant --coupling 1e-3

# Linear ramp
python time_varying_coupling_example.py --variant ramp --coupling 1e-3 --final_coupling 2e-3 --ramp_time 10000

# Oscillating coupling
python time_varying_coupling_example.py --variant cycle --coupling 5e-4 --final_coupling 1.5e-3 --cycle_time 3000

# Power law approach
python time_varying_coupling_example.py --variant power --coupling 1e-4 --final_coupling 1e-3 --power 0.5

# Generate plots
python time_varying_coupling_example.py --variant ramp --plot
```

## Performance Notes

- **Evaluation overhead**: Variants are evaluated at each timestep. For complex custom variants, consider caching if the function is expensive.
- **GPU compatibility**: All variant types work with both CPU and GPU implementations.
- **Numerical precision**: Ensure your variant doesn't have discontinuities that could cause numerical issues.

## Troubleshooting

**Issue**: Sudden energy jumps
**Solution**: Use smoother transitions (longer ramp times) or add appropriate thermostats.

**Issue**: Variant not changing coupling
**Solution**: Verify the variant is being evaluated correctly by logging the coupling strength.

**Issue**: Simulation becomes unstable
**Solution**: Reduce the rate of change or add damping to the cavity mode. 