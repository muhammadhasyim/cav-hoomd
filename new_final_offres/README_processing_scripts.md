# Cavity MD Processing Scripts

This directory contains two specialized processing scripts for analyzing cavity MD simulation data:

## 1. `process_fskt_only.py` - F(k,t) Processing and Analysis

Processes F(k,t) (structure factor autocorrelation) files from cavity MD simulations.

### File Pattern
- Input files: `prod*_ref*.txt` (where `*` can be replica indices and reference numbers)
- Example: `prod-0_ref1.txt`, `prod-1_ref1.txt`, etc.

### Key Features
- Processes F(k,t) correlation functions
- Creates master files by averaging across replicas for each reference group
- Generates analysis plots showing F(k,t) decay over time
- Configurable timestep resolution and time ranges

### Usage
```bash
python process_fskt_only.py --exp_dir ./my_experiment
python process_fskt_only.py --exp_dir ./cavity_coupling_1eneg03_switch_1.0ps --fskt_dt 0.5 --max_time 100.0
```

### Output Files
- `master_fskt_ref*.txt` - Master F(k,t) files for each reference group
- `fskt_overview_analysis.png` - Individual F(k,t) plots for each reference
- `fskt_combined_analysis.png` - Combined plot showing all references
- `fskt_analysis_summary.txt` - Summary report

---

## 2. `process_dipole_autocorr.py` - Dipole Autocorrelation and IR Spectrum Analysis

Processes dipole autocorrelation files and generates IR spectra using Discrete Cosine Transform (DCT).

### File Pattern
- Input files: `prod-i_dipole_autocorr_j.txt` where:
  - `i` = replica index
  - `j` = reference index
- Example: `prod-0_dipole_autocorr_1.txt`, `prod-1_dipole_autocorr_1.txt`, etc.

### Key Features
- Processes dipole moment autocorrelation functions
- Creates master autocorrelation files by averaging across replicas
- **Calculates IR spectra using DCT with quantum corrections**
- Applies proper temperature-dependent corrections for realistic IR intensities
- Configurable temperature for quantum corrections

### Usage
```bash
python process_dipole_autocorr.py --exp_dir ./my_experiment
python process_dipole_autocorr.py --exp_dir ./cavity_coupling_1eneg03_switch_1.0ps --temperature 300 --autocorr_dt 1.0
```

### Output Files
- `master_dipole_autocorr_ref*.txt` - Master autocorrelation files for each reference group
- `ir_spectrum_ref*.txt` - IR spectrum data files for each reference
- `dipole_autocorr_overview.png` - Individual autocorrelation plots
- `dipole_autocorr_and_ir_spectra.png` - **Combined plot showing both autocorrelations AND their corresponding IR spectra**
- `dipole_autocorr_analysis_summary.txt` - Summary report

---

## Common Features

Both scripts share the following features:

### Processing Options
- `--exp_dir` (required): Experiment directory to process
- `--job_name`: Job prefix for input files (default: 'prod')
- `--max_time`: Maximum time for analysis (auto-detect if not specified)
- `--skip_processing`: Skip processing, only create plots from existing master files
- `--running_sim`: Optimized for running simulations (conservative time detection)

### Data Handling
- Robust error handling for incomplete or malformed files
- Automatic time range detection for running simulations
- Configurable timestep resolution for interpolation
- Online averaging algorithm for memory efficiency
- Sample count tracking for data quality assessment

### File Organization
- Groups files by reference number automatically
- Creates comprehensive master files with detailed headers
- Generates sample count files for quality control
- Produces publication-ready plots and summary reports

---

## IR Spectrum Calculation Details

The `process_dipole_autocorr.py` script implements a complete IR spectrum calculation pipeline:

1. **Autocorrelation Processing**: Averages dipole autocorrelation functions across replicas
2. **DCT Transform**: Uses `scipy.fftpack.dct(type=1)` to convert time-domain to frequency-domain
3. **Quantum Corrections**: Applies temperature-dependent quantum mechanical corrections
4. **Field Description**: Includes proper classical-to-quantum field transformation
5. **Frequency Conversion**: Converts from Hz to cm⁻¹ wavenumbers

The resulting IR spectra are quantum-corrected and suitable for comparison with experimental infrared spectroscopy data.

---

## Example Workflow

1. **For F(k,t) analysis**:
   ```bash
   python process_fskt_only.py --exp_dir ./my_cavity_experiment
   ```

2. **For IR spectrum generation**:
   ```bash
   python process_dipole_autocorr.py --exp_dir ./my_cavity_experiment --temperature 300
   ```

Both scripts can be run independently and will automatically detect and process the appropriate file types in your experiment directory. 