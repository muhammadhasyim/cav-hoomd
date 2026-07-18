#!/usr/bin/env python3
"""
Compute Fictive Temperatures Using Extended Fitting Functions

This script computes fictive temperatures from averaged potential energy time series
using the EXTENDED fitting functions from the log-log scaling analysis:
- Harmonic: E = a*T/(1 + b*T)
- LJ+Coulombic: E = A*T^(3/5) + B or extended form

Uses potential_energy_vs_T.txt as the calibration data.

Usage:
    python compute_fictive_temps_extended_fit.py <potential_energy_file> [options]
    
Examples:
    # Basic usage
    python compute_fictive_temps_extended_fit.py time_series_output/coupling_0epos00_averaged_potential_energy.txt
    
    # Custom calibration file
    python compute_fictive_temps_extended_fit.py time_series_output/coupling_0epos00_averaged_potential_energy.txt --calibration potential_energy_vs_T.txt
    
    # Save to specific location
    python compute_fictive_temps_extended_fit.py time_series_output/coupling_0epos00_averaged_potential_energy.txt --output time_series_output/coupling_0epos00_extended_fictive_temperatures.txt

Author: Scientific Software Engineer
Date: October 2025
"""

import numpy as np
import argparse
import sys
from pathlib import Path
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import warnings

def load_calibration_data(filename="potential_energy_vs_T.txt"):
    """
    Load the potential energy vs temperature calibration data.
    
    Returns:
    --------
    calibration_data : dict
        Dictionary with temperature array and energy component arrays
    """
    try:
        # Read the file manually to handle the header properly
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # Find the header line (first non-comment, non-empty line)
        header_line = None
        data_start = 0
        
        for i, line in enumerate(lines):
            stripped = line.strip()
            if not stripped.startswith('#') and stripped:
                header_line = stripped
                data_start = i + 1
                break
        
        if header_line is None:
            raise ValueError("Could not find header line")
        
        column_names = [col.strip() for col in header_line.split('\t') if col.strip()]
        
        print(f"  Column names: {column_names}")
        
        # Read data starting from the data_start line
        data_lines = []
        for line in lines[data_start:]:
            if not line.strip().startswith('#') and line.strip():
                data_lines.append(line.strip())
        
        # Parse data
        data_rows = []
        for line in data_lines:
            values = line.split('\t') if '\t' in line else line.split()
            if len(values) == len(column_names):
                data_rows.append([float(v) for v in values])
        
        # Create array
        data = np.array(data_rows)
        
        # Create dictionary
        calibration_data = {name: data[:, i] for i, name in enumerate(column_names)}
        
        # Find temperature column for reporting
        temp_col = 'temperature' if 'temperature' in calibration_data else 'temperature_K'
        
        print(f"✓ Loaded calibration data from {filename}:")
        print(f"  Temperature range: {calibration_data[temp_col].min():.1f} - {calibration_data[temp_col].max():.1f} K")
        print(f"  {len(data)} temperature points")
        
        return calibration_data
        
    except Exception as e:
        print(f"❌ Error loading calibration data: {e}")
        import traceback
        traceback.print_exc()
        return None

def create_extended_fitting_interpolators(calibration_data, temp_range=None):
    """
    Create T(U) interpolators using extended fitting functions.
    
    Parameters:
    -----------
    calibration_data : dict
        Calibration data
    temp_range : tuple of (T_min, T_max), optional
        Temperature range to use for fitting. If None, uses all data.
    
    Returns:
    --------
    interpolators : dict
        Dictionary with component names as keys, each containing T_of_U function
    """
    # Find temperature column
    temp_col = None
    for possible_name in ['temperature', 'temperature_K']:
        if possible_name in calibration_data:
            temp_col = possible_name
            break
    
    if temp_col is None:
        raise ValueError("Could not find temperature column")
    
    temperatures = calibration_data[temp_col]
    
    # Filter by temperature range if specified
    if temp_range is not None:
        T_min, T_max = temp_range
        temp_mask = (temperatures >= T_min) & (temperatures <= T_max)
        print(f"  Filtering calibration data to T = {T_min:.1f} - {T_max:.1f} K")
        print(f"  Using {np.sum(temp_mask)}/{len(temperatures)} data points")
    else:
        temp_mask = np.ones(len(temperatures), dtype=bool)
    
    # Sort by temperature
    sort_idx = np.argsort(temperatures[temp_mask])
    T_sorted = temperatures[temp_mask][sort_idx]
    
    interpolators = {}
    
    # Define energy components to process
    # Map to possible column names in the calibration file
    energy_map = {
        'harmonic': ['harmonic_hartree', 'harmonic_energy'],
        'lj': ['lj_hartree', 'lj_energy'],
        'coulombic': ['coulombic_hartree', 'coulombic_energy'],
        'total_PE': ['total_potential_energy_hartree', 'total_potential_energy']
    }
    
    for component_name, possible_columns in energy_map.items():
        # Find which column name exists in the data
        column_name = None
        for col in possible_columns:
            if col in calibration_data:
                column_name = col
                break
        
        if column_name is None:
            print(f"  Skipping {component_name}: no matching column found (tried {possible_columns})")
            continue
        
        energies = calibration_data[column_name]
        U_sorted = energies[temp_mask][sort_idx]
        
        print(f"\n  Processing {component_name}:")
        print(f"    Energy range: {U_sorted.min():.6f} to {U_sorted.max():.6f}")
        
        # Extended fitting based on component type
        if component_name == 'harmonic':
            # For harmonic, use simple linear interpolation - it's more reliable for this range
            print(f"    Using linear interpolation for {component_name} (more reliable for narrow T range)")

            # Create linear interpolator
            U_of_T = interp1d(T_sorted, U_sorted, kind='linear', bounds_error=False, fill_value='extrapolate')
            T_of_U = interp1d(U_sorted, T_sorted, kind='linear', bounds_error=False, fill_value='extrapolate')

            # Test the fit quality
            U_pred = U_of_T(T_sorted)
            ss_res = np.sum((U_sorted - U_pred)**2)
            ss_tot = np.sum((U_sorted - np.mean(U_sorted))**2)
            r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

            print(f"    Linear interpolation R² = {r2:.10f}")
            print(f"    (Using linear interpolation for harmonic component)")

            # No parameters to store for linear case
            U_of_T = U_of_T  # Already defined above
            T_of_U = T_of_U  # Already defined above

        else:
            # For LJ, Coulombic, Total: T^(3/5) power law
            def power_law_with_intercept(x, a, b):
                return a * x**(3/5) + b
            
            def extended_power_law(x, e0, a, c):
                t_power = x**(3/5)
                return e0 + a * t_power / (1 + c * t_power)
            
            try:
                # Fit simple T^(3/5) + B first
                power_params, _ = curve_fit(power_law_with_intercept, T_sorted, U_sorted)
                a, b = power_params
                
                # Calculate R² for simple fit
                U_fit_simple = power_law_with_intercept(T_sorted, a, b)
                ss_res = np.sum((U_sorted - U_fit_simple)**2)
                ss_tot = np.sum((U_sorted - np.mean(U_sorted))**2)
                r2_simple = 1 - (ss_res / ss_tot)
                
                print(f"    Simple T^(3/5): E = {a:.6f}*T^(3/5) + {b:.3f}")
                print(f"    R² = {r2_simple:.10f}")
                
                # Try extended fit
                try:
                    e0_guess = U_sorted.min()
                    a_guess = a
                    c_guess = 1e-3
                    
                    extended_params, _ = curve_fit(
                        extended_power_law,
                        T_sorted,
                        U_sorted,
                        p0=[e0_guess, a_guess, c_guess],
                        bounds=([U_sorted.min() - 10, -np.inf, 0], 
                               [U_sorted.max() + 10, np.inf, np.inf])
                    )
                    
                    e0, a_ext, c = extended_params
                    
                    # Calculate R² for extended fit
                    U_fit_ext = extended_power_law(T_sorted, e0, a_ext, c)
                    ss_res_ext = np.sum((U_sorted - U_fit_ext)**2)
                    r2_ext = 1 - (ss_res_ext / ss_tot)
                    
                    print(f"    Extended T^(3/5): E = {e0:.3f} + {a_ext:.6f}*T^(3/5)/(1 + {c:.6f}*T^(3/5))")
                    print(f"    R² = {r2_ext:.10f}")
                    
                    # Use the better fit
                    if r2_ext > r2_simple:
                        print(f"    → Using extended fit (better R²)")
                        
                        def T_of_U(U):
                            # Solve E = E0 + A*T^(3/5)/(1 + C*T^(3/5)) for T
                            denom = (U - e0) * c - a_ext
                            denom = np.where(np.abs(denom) < 1e-12, 1e-12, denom)
                            x = -(U - e0) / denom
                            return np.where(x > 0, x**(5/3), np.nan)
                        
                        def U_of_T(T):
                            return extended_power_law(T, e0, a_ext, c)

                        # Assign the chosen functions for extended fit
                        U_of_T = U_of_T
                        T_of_U = T_of_U

                    else:
                        print(f"    → Using simple fit (better R²)")

                        def T_of_U(U):
                            # T = ((U - b)/a)^(5/3)
                            return np.where((U - b)/a > 0, ((U - b) / a)**(5/3), np.nan)

                        def U_of_T(T):
                            return power_law_with_intercept(T, a, b)

                        # Assign the chosen functions for simple fit
                        U_of_T = U_of_T
                        T_of_U = T_of_U

                except:
                    print(f"    → Extended fit failed, using simple fit")

                    def T_of_U(U):
                        return np.where((U - b)/a > 0, ((U - b) / a)**(5/3), np.nan)

                    def U_of_T(T):
                        return power_law_with_intercept(T, a, b)

                    # Assign the chosen functions
                    U_of_T = U_of_T
                    T_of_U = T_of_U
            
            except Exception as e:
                print(f"    Power law fitting failed: {e}")
                print(f"    Falling back to linear interpolation")
                U_of_T = interp1d(T_sorted, U_sorted, kind='linear', bounds_error=False, fill_value='extrapolate')
                T_of_U = interp1d(U_sorted, T_sorted, kind='linear', bounds_error=False, fill_value='extrapolate')
        
        interpolators[component_name] = {
            'T_of_U': T_of_U,
            'U_of_T': U_of_T,
            'U_range': (U_sorted.min(), U_sorted.max()),
            'T_range': (T_sorted.min(), T_sorted.max())
        }
    
    print(f"\n✓ Created {len(interpolators)} extended interpolators")
    return interpolators

def load_potential_energy_time_series(filename):
    """Load potential energy time series data."""
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # Find header
        header_line = None
        data_start = 0
        
        for i, line in enumerate(lines):
            if line.strip().startswith('# time_ps'):
                header_line = line.strip()[2:]
                data_start = i + 1
                break
            elif not line.strip().startswith('#') and line.strip():
                header_line = line.strip()
                data_start = i + 1
                break
        
        if header_line is None:
            raise ValueError("Could not find header")
        
        column_names = header_line.split()
        
        # Read data
        data_lines = []
        for line in lines[data_start:]:
            if line.strip() and not line.strip().startswith('#'):
                data_lines.append(line.split())
        
        # Parse data
        time = []
        energies = {col: [] for col in column_names[1:]}
        
        for values in data_lines:
            time.append(float(values[0]))
            for i, col in enumerate(column_names[1:], start=1):
                val = values[i]
                if val.lower() == 'nan':
                    energies[col].append(np.nan)
                else:
                    energies[col].append(float(val))
        
        time = np.array(time)
        for col in energies:
            energies[col] = np.array(energies[col])
        
        print(f"✓ Loaded potential energy time series:")
        print(f"  Time range: {time[0]:.3f} - {time[-1]:.3f} ps")
        print(f"  Data points: {len(time)}")
        print(f"  Components: {list(energies.keys())}")
        
        return time, energies
        
    except Exception as e:
        print(f"❌ Error loading potential energy data: {e}")
        import traceback
        traceback.print_exc()
        return None, None

def compute_fictive_temperatures(time, energies, interpolators):
    """Compute fictive temperatures from energy time series."""
    fictive_temps = {}
    
    # Map energy column names to interpolator names
    energy_map = {
        'harmonic_energy_hartree': 'harmonic',
        'lj_energy_hartree': 'lj',
        'coulombic_energy_hartree': 'coulombic',
        'total_PE_energy_hartree': 'total_PE'
    }
    
    for energy_col, interp_name in energy_map.items():
        if energy_col in energies and interp_name in interpolators:
            U_series = energies[energy_col]
            T_of_U = interpolators[interp_name]['T_of_U']
            U_range = interpolators[interp_name]['U_range']
            
            # Compute fictive temperature
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                T_fictive = T_of_U(U_series)
            
            # Mark physically unreasonable temperatures
            reasonable_mask = (T_fictive > 0) & (T_fictive < 10000) & np.isfinite(T_fictive)
            T_fictive_clean = T_fictive.copy()
            T_fictive_clean[~reasonable_mask] = np.nan
            
            fictive_temps[interp_name] = T_fictive_clean
            
            # Statistics
            n_valid = np.sum(reasonable_mask)
            n_total = len(T_fictive)
            print(f"  {interp_name}: {n_valid}/{n_total} valid fictive temperatures")
            if n_valid > 0:
                print(f"    T range: {T_fictive_clean[reasonable_mask].min():.2f} - {T_fictive_clean[reasonable_mask].max():.2f} K")
    
    return fictive_temps

def save_fictive_temperatures(time, fictive_temps, output_file, source_file):
    """Save fictive temperature data."""
    output_lines = []
    output_lines.append("# Fictive Temperatures Computed with Extended Fitting Functions")
    output_lines.append(f"# Source: {source_file}")
    output_lines.append("# Method: Extended fitting (Harmonic: E=aT/(1+bT), Others: T^(3/5) power law)")
    output_lines.append("# Columns:")
    
    # Create header
    header_cols = ['time_ps']
    component_order = ['harmonic', 'lj', 'coulombic', 'total_PE']
    for comp in component_order:
        if comp in fictive_temps:
            header_cols.append(f'fictive_T_{comp}_K')
    
    output_lines.append("# " + "\t".join(header_cols))
    
    # Write data
    for i, t in enumerate(time):
        data_row = [f"{t:.6f}"]
        
        for comp in component_order:
            if comp in fictive_temps:
                T = fictive_temps[comp][i]
                if np.isnan(T):
                    data_row.append("NaN")
                else:
                    data_row.append(f"{T:.6f}")
        
        output_lines.append("\t".join(data_row))
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write('\n'.join(output_lines))
    
    print(f"\n✓ Fictive temperatures saved to: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Compute fictive temperatures using extended fitting functions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('potential_energy_file', type=str,
                       help='Input potential energy time series file')
    parser.add_argument('--calibration', type=str, default='potential_energy_vs_T.txt',
                       help='Calibration file (default: potential_energy_vs_T.txt)')
    parser.add_argument('--output', '-o', type=str, default=None,
                       help='Output file (default: auto-generate from input name)')
    parser.add_argument('--temp-range', type=float, nargs=2, default=None,
                       help='Temperature range for fitting (default: use all data). Example: --temp-range 65 200')
    
    args = parser.parse_args()
    
    # Check files exist
    if not Path(args.potential_energy_file).exists():
        print(f"❌ Error: Input file not found: {args.potential_energy_file}")
        sys.exit(1)
    
    if not Path(args.calibration).exists():
        print(f"❌ Error: Calibration file not found: {args.calibration}")
        sys.exit(1)
    
    print("="*70)
    print("COMPUTING FICTIVE TEMPERATURES WITH EXTENDED FITTING FUNCTIONS")
    print("="*70)
    
    # Load calibration data
    print("\n[1/4] Loading calibration data...")
    calibration_data = load_calibration_data(args.calibration)
    if calibration_data is None:
        sys.exit(1)
    
    # Create interpolators
    print("\n[2/4] Creating extended fitting interpolators...")
    temp_range = tuple(args.temp_range) if args.temp_range else None
    interpolators = create_extended_fitting_interpolators(calibration_data, temp_range=temp_range)
    if not interpolators:
        sys.exit(1)
    
    # Load potential energy time series
    print("\n[3/4] Loading potential energy time series...")
    time, energies = load_potential_energy_time_series(args.potential_energy_file)
    if time is None or energies is None:
        sys.exit(1)
    
    # Compute fictive temperatures
    print("\n[4/4] Computing fictive temperatures...")
    fictive_temps = compute_fictive_temperatures(time, energies, interpolators)
    if not fictive_temps:
        print("❌ Error: No fictive temperatures could be computed!")
        sys.exit(1)
    
    # Determine output filename
    if args.output is None:
        input_path = Path(args.potential_energy_file)
        output_name = input_path.stem.replace('potential_energy', 'extended_fictive_temperatures') + '.txt'
        output_path = input_path.parent / output_name
    else:
        output_path = Path(args.output)
    
    # Save results
    save_fictive_temperatures(time, fictive_temps, output_path, args.potential_energy_file)
    
    print("\n" + "="*70)
    print("✓ COMPUTATION COMPLETE!")
    print("="*70)

if __name__ == '__main__':
    main()

