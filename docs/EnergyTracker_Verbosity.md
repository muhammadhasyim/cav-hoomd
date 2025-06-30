# EnergyTracker Verbosity Control

The `EnergyTracker` class now supports verbosity control to manage the amount of debug output printed to the console during simulations.

## Verbosity Levels

### `'quiet'`
- **Minimal output**: Only essential messages and errors are printed
- **Use case**: Production runs where you want clean console output
- **Behavior**: 
  - No energy component details during simulation
  - Error messages are still printed
  - Energy data is still written to output files

### `'normal'` (default)
- **Moderate output**: Setup information and error messages are shown
- **Use case**: Regular simulations where you want some feedback
- **Behavior**:
  - Shows EnergyTracker setup information
  - Shows error messages and warnings
  - No detailed energy component output during simulation
  - Energy data is still written to output files

### `'verbose'`
- **Full debug output**: Shows all detailed energy information (original behavior)
- **Use case**: Debugging and detailed monitoring of energy components
- **Behavior**:
  - Shows detailed energy components for each timestep
  - Shows all intermediate calculations
  - Shows timing and status information
  - Identical to the original EnergyTracker behavior

## Usage Examples

### Quiet Mode (Recommended for production)
```python
energy_tracker = EnergyTracker(
    simulation=sim,
    components=['harmonic', 'lj', 'ewald_short', 'ewald_long', 'cavity'],
    force_objects=force_objects,
    thermostat_objects=thermostat_objects,
    verbose='quiet'  # Minimal console output
)
```

### Normal Mode (Default)
```python
energy_tracker = EnergyTracker(
    simulation=sim,
    components=['harmonic', 'lj', 'ewald_short', 'ewald_long', 'cavity'],
    force_objects=force_objects,
    thermostat_objects=thermostat_objects,
    verbose='normal'  # or omit verbose parameter (default)
)
```

### Verbose Mode (For debugging)
```python
energy_tracker = EnergyTracker(
    simulation=sim,
    components=['harmonic', 'lj', 'ewald_short', 'ewald_long', 'cavity'],
    force_objects=force_objects,
    thermostat_objects=thermostat_objects,
    verbose='verbose'  # Full debug output
)
```

## Migration Guide

### For existing code:
- **No changes needed**: Existing code will continue to work with 'normal' verbosity level
- **To reduce output**: Add `verbose='quiet'` to existing EnergyTracker instantiations
- **To maintain old behavior**: Add `verbose='verbose'` to existing EnergyTracker instantiations

### Output file behavior:
- **No changes**: All verbosity levels write the same detailed energy data to output files
- **Only console output is affected** by the verbosity setting
- File format and content remain identical regardless of verbosity level

## Typical Console Output Comparison

### Quiet Mode
```
[Minimal or no output during simulation]
```

### Normal Mode  
```
CORRECTED EnergyTracker (following working code pattern):
  Output file: prod-0-energy_energy_tracker.txt
  Components: ['harmonic', 'lj', 'ewald_short', 'ewald_long', 'cavity']
  ...
EnergyTracker: Successfully created output file prod-0-energy_energy_tracker.txt
```

### Verbose Mode
```
CORRECTED EnergyTracker (following working code pattern):
  Output file: prod-0-energy_energy_tracker.txt
  ...

=== ENERGY TRACKER DEBUG - Timestep 1000 ===
Current time from time_tracker: 0.073147 ps
=== POTENTIAL ENERGY COMPONENTS ===
Harmonic energy: 0.131957 Hartree
LJ energy: -0.719588 Hartree
Ewald short energy: -0.342615 Hartree
Ewald long energy: -0.403265 Hartree
...
```

This verbosity control helps maintain clean console output while preserving all the detailed energy tracking functionality in the output files. 