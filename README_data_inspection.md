# HDF5 Data Inspection Scripts

This directory contains two Python scripts for inspecting and analyzing HDF5 simulation data from CavityMD simulations.

## Scripts Overview

### 1. `quick_inspect.py` - Quick Data Overview
Simple script for immediate data viewing and basic statistics.

**Usage:**
```bash
python quick_inspect.py <hdf5_file>
```

**Example:**
```bash
python quick_inspect.py observables_replica_0.h5
```

**Features:**
- File attributes and metadata
- Time data range and sampling
- First 5 energy components with statistics
- All temperature components with statistics
- Quick data preview (first/last 5 values)

### 2. `inspect_hdf5_data.py` - Comprehensive Analysis
Full-featured script with plotting, export, and detailed analysis capabilities.

**Usage:**
```bash
python inspect_hdf5_data.py <hdf5_file> [OPTIONS]
```

**Options:**
- `--show-values`: Show actual data values in structure view
- `--max-values N`: Maximum values to show per dataset (default: 10)
- `--plot`: Generate energy and temperature evolution plots
- `--export-csv`: Export all data to CSV file
- `--print-arrays`: Print all data arrays with values
- `--output-dir DIR`: Output directory for plots (default: plots)

**Examples:**
```bash
# Basic file structure inspection
python inspect_hdf5_data.py observables_replica_0.h5

# Generate plots and export CSV
python inspect_hdf5_data.py observables_replica_0.h5 --plot --export-csv

# Print all data arrays
python inspect_hdf5_data.py observables_replica_0.h5 --print-arrays

# Full analysis with custom output directory
python inspect_hdf5_data.py observables_replica_0.h5 --plot --export-csv --print-arrays --output-dir my_plots
```

## Output Files

### Generated Plots
- `energy_evolution.png`: 4-panel plot showing:
  - Cavity energies (coupling, dipole_self, harmonic, kinetic, reservoir)
  - Molecular energies (LJ, harmonic, Ewald short/long)
  - Total energies (kinetic, potential, reservoir, system)
  - System temperature
- `temperature_evolution.png`: 2-panel plot showing:
  - All temperature components
  - Key temperature components (kinetic, cavity_bath, molecular_bath)

### CSV Export
- `simulation_data.csv`: All datasets exported as columns with time as first column
- Shape: (N_timepoints, N_datasets)
- Compatible with pandas, Excel, and other analysis tools

## Data Structure

### Energy Components (`energies/` group)
- `cavity_*`: Cavity-related energies (coupling, dipole_self, harmonic, kinetic, reservoir, total_potential)
- `molecular_*`: Molecular system energies (kinetic, reservoir)
- `ewald_*`: Electrostatic energies (short, long)
- `harmonic`, `lj`: Bonded and non-bonded molecular interactions
- `*_total`: System totals (kinetic, potential, reservoir, system, universe)
- `temperature`: Temperature from kinetic energy

### Temperature Components (`temperatures/` group)
- `kinetic`: Direct kinetic temperature
- `*_fictive`: Fictive temperatures from potential energies
- `*_bath`: Thermostat bath temperatures
- `*_equipartition`: Equipartition-based temperatures

### Time Data
- `time`: Simulation time in picoseconds
- Regular sampling intervals (typically 0.01 ps)

## Common Issues and Solutions

### HDF5 File Access Issues
If you get "unable to open file" errors:
```bash
# Clear HDF5 consistency flags
h5clear -s your_file.h5

# Then try accessing again
python quick_inspect.py your_file.h5
```

### Missing Dependencies
```bash
# Install required packages
pip install h5py numpy matplotlib pandas
```

### Large Files
For very large files, use:
- `--max-values 5` to limit printed values
- Skip `--print-arrays` to avoid memory issues
- Use `--export-csv` to work with data in external tools

## Example Workflow

1. **Quick overview:**
   ```bash
   python quick_inspect.py observables_replica_0.h5
   ```

2. **Generate visualizations:**
   ```bash
   python inspect_hdf5_data.py observables_replica_0.h5 --plot
   ```

3. **Export for further analysis:**
   ```bash
   python inspect_hdf5_data.py observables_replica_0.h5 --export-csv
   ```

4. **Full analysis:**
   ```bash
   python inspect_hdf5_data.py observables_replica_0.h5 --plot --export-csv --print-arrays
   ```

## Tips

- Use `quick_inspect.py` first to get an overview
- Generate plots to visualize trends and identify interesting features
- Export to CSV for detailed analysis in Python/R/Excel
- Check file attributes for simulation parameters
- Look for NaN values in temperature data (normal at t=0)
- Cavity energies should be zero before coupling activation
- System total energy should be conserved (with small fluctuations)
