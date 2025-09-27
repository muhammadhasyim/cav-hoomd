# FDR-Based Effective Temperature Measurement

## Overview

The FDR (Fluctuation-Dissipation Ratio) temperature estimator provides a **physics-based method** for measuring mode-specific effective temperatures in non-equilibrium molecular dynamics simulations. Unlike empirical approaches, this method directly measures violations of fluctuation-dissipation relations to determine `T_eff(ω₀,t)`.

### Key Formula

```
T_eff(ω₀,t) = (ω₀/2k_B) × S_AA(ω₀,t) / χ''(ω₀,t)
```

Where:
- `S_AA(ω₀,t)`: Power spectral density at frequency ω₀
- `χ''(ω₀,t)`: Imaginary susceptibility at frequency ω₀
- This ratio measures deviations from equilibrium FDT

## Quick Integration with Existing Simulations

### 1. Add FDR Monitor to Your Simulation

```python
from cavitymd.fdr_integration import create_cavity_fdr_monitor

# For your cavity coupling simulations
fdr_monitor = create_cavity_fdr_monitor(
    cavity_frequency_cm=2000.0,  # Use YOUR cavity frequency
    simulation_dt=0.001,         # Your timestep in ps
    observable_type='dipole',    # Monitor total dipole moment
    axis='z',                    # Project along z-axis
    output_file='fdr_temp_trajectory.dat'
)

# Attach to your existing simulation
fdr_monitor.attach_to_simulation(sim)
```

### 2. Calibrate Once at Equilibrium

```python
# Run equilibration and calibrate
T_bath = 300.0  # Your bath temperature
fdr_monitor.calibrate_equilibrium(T_cal=T_bath, n_steps=10000)
```

### 3. Monitor During Production

```python
# In your simulation loop
for step in range(production_steps):
    sim.run(10)  # Run 10 steps
    
    # Update FDR temperature estimate
    T_eff, diagnostics = fdr_monitor.update()
    
    # T_eff is the effective temperature at cavity frequency
    # diagnostics contains detailed information
```

## Avoiding Interference with Running Simulations

### Safe Frequency Choices
Choose a frequency **different** from your current simulations:

```python
# If you're running 1e-3 coupling, use different frequency for testing
test_frequencies = [1500.0, 2500.0, 3000.0]  # cm^-1, avoid your current setup

# For mode-specific analysis, extract vibrational modes first
fdr_monitor = create_cavity_fdr_monitor(
    cavity_frequency_cm=1500.0,  # Different from your current work
    simulation_dt=your_dt,
    observable_type='dipole',
    axis='z'
)
```

### Minimal Computational Overhead
- **O(1) per timestep** - constant memory usage
- **No large buffers** - exponential forgetting
- **Parallel safe** - doesn't modify simulation state

## Scientific Applications

### 1. Cavity QED Non-Equilibrium Dynamics
```python
# Monitor how cavity coupling affects molecular thermalization
cavity_frequencies = [1000, 2000, 3000]  # cm^-1
for freq in cavity_frequencies:
    monitor = create_cavity_fdr_monitor(cavity_frequency_cm=freq, ...)
    # Compare T_eff at different frequencies
```

### 2. Mode-Specific Temperature Measurement
```python
# Extract specific vibrational mode
mode_vector = get_normal_mode(molecule, mode_index=5)  # Your function
masses = get_particle_masses(molecule)

fdr_monitor = FDRTemperatureMonitor(
    omega_0=mode_frequency,
    dt=sim_dt,
    observable_type='mode',
    mode_vector=mode_vector,
    masses=masses
)
```

### 3. Validation of Thermalization Assumptions
```python
# Compare FDR temperature with equipartition
analyzer = CavityFDRAnalyzer(fdr_monitor)
comparison = analyzer.compare_with_equipartition(bath_temp)

if comparison['is_non_equilibrium']:
    print(f"Mode shows non-equilibrium: {comparison['relative_deviation']*100:.1f}% deviation")
```

## Output Files and Analysis

### FDR Trajectory File Format
```
# Columns: step, time_ps, T_eff_K, omega_n_rad_ps, gamma_inv_ps, 
#          sigma_delta_rad_ps, S_AA, chi_imag, snr, lineshape_type
     1000     1.000   305.23    1.000234    0.098123    0.012345    1.23e-04    2.34e-05    12.5  lorentzian
     1010     1.010   304.87    1.000156    0.098234    0.012234    1.22e-04    2.33e-05    12.8  lorentzian
```

### Key Diagnostics
- **T_eff**: Effective temperature at target frequency
- **omega_n**: Identified natural frequency
- **gamma**: Damping coefficient  
- **SNR**: Signal-to-noise ratio
- **lineshape_type**: Lorentzian vs Voigt (dephasing regime)

## Integration with Your Current Workflow

### With Temperature Trackers
```python
# Add FDR alongside your existing temperature tracking
from cavitymd.analysis import TemperatureTracker

# Your existing tracker
temp_tracker = TemperatureTracker(...)

# Add FDR monitor
fdr_monitor = create_cavity_fdr_monitor(...)

# Both run simultaneously
for step in simulation_loop:
    # Your current analysis
    temp_tracker.update(sim.state)
    
    # Add FDR analysis
    T_eff, diag = fdr_monitor.update()
```

### With GSD Trajectory Analysis
```python
# Post-process existing trajectories
import gsd.hoomd

trajectory = gsd.hoomd.open('trajectory.gsd', 'r')
fdr_monitor = create_cavity_fdr_monitor(...)

# Calibrate on equilibrium frames
equilibrium_frames = trajectory[-1000:]  # Last 1000 frames
# ... calibration process ...

# Analyze trajectory
for frame in trajectory:
    # Extract observable and update FDR
    # ... processing ...
```

## Performance Considerations

### Computational Cost
- **Per-step overhead**: ~10⁻⁶ seconds (negligible)
- **Memory usage**: ~1 KB (constant)
- **Frequency**: Update every 1-10 MD steps

### Recommended Settings
```python
# For production simulations
fdr_monitor = FDRTemperatureMonitor(
    omega_0=your_frequency,
    dt=your_dt,
    tau_avg=100*T_period,    # 100 vibrational periods
    tau_id=500*T_period,     # 500 periods for parameter ID
    log_interval=1000        # Log every 1000 steps
)
```

## Troubleshooting

### Common Issues
1. **NaN temperatures**: Check calibration and signal amplitude
2. **Unstable parameters**: Increase `tau_id` for more averaging
3. **Poor SNR**: Ensure observable couples to target frequency

### Validation Checks
```python
# Check diagnostics
_, diag = fdr_monitor.update()
if diag.snr < 3.0:
    print("Warning: Low signal-to-noise ratio")
if diag.effective_dof < 10:
    print("Warning: Insufficient statistics")
```

## Next Steps

1. **Test with synthetic data** using the validation tests
2. **Integrate gradually** with one frequency first
3. **Compare with existing methods** to validate
4. **Scale to multiple modes** for comprehensive analysis

The FDR estimator provides a rigorous, physics-based approach to temperature measurement that complements your existing cavity QED simulation infrastructure.
