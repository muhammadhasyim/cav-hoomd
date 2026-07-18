#!/usr/bin/env python3
"""
Reconstruct Averaged Potential Energy Time Series from Fictive Temperature Data

This script reverses the fictive temperature calculation to recover the averaged
potential energy time series. It uses the OLD linear interpolation method 
(not the extended fitting functions) to match the original calculation.

Usage:
    python reconstruct_averaged_potential_energy.py <fictive_temp_file> [options]
    
Examples:
    # Basic usage
    python reconstruct_averaged_potential_energy.py time_series_output/coupling_0epos00_averaged_fictive_temperatures.txt
    
    # Specify output file
    python reconstruct_averaged_potential_energy.py time_series_output/coupling_0epos00_averaged_fictive_temperatures.txt --output coupling_0epos00_averaged_potential_energy.txt
    
    # Custom calibration file
    python reconstruct_averaged_potential_energy.py time_series_output/coupling_0epos00_averaged_fictive_temperatures.txt --temp-components my_calibration.txt

Author: Scientific Software Engineer
Date: October 2025
"""

import numpy as np
import argparse
import sys
from pathlib import Path
from scipy.interpolate import interp1d

def load_temperature_energy_calibration(filename="potential_energy_components_vs_temperature.txt"):
    """
    Load the temperature-energy calibration data.
    
    Returns:
    --------
    calibration_data : dict
        Dictionary with component names as keys, each containing:
        - 'T': temperature array
        - 'U': energy array  
        - 'U_of_T': interpolation function U(T)
    """
    try:
        # Read the file
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # Find the header line
        header_line = None
        data_start_line = 0
        
        for i, line in enumerate(lines):
            line = line.strip()
            if not line.startswith('#') and line:
                header_line = line
                data_start_line = i + 1
                break
        
        if header_line is None:
            raise ValueError("Could not find column headers")
        
        column_names = [col.strip() for col in header_line.split()]
        temp_columns = {name: i for i, name in enumerate(column_names)}
        
        # Read the data - handle tab or space separated
        data_lines = []
        for line in lines[data_start_line:]:
            if not line.strip().startswith('#') and line.strip():
                values = line.split('\t') if '\t' in line else line.split()
                if len(values) == len(column_names):
                    data_lines.append([float(v) for v in values])
        
        if not data_lines:
            raise ValueError("No valid data lines found")
        
        temp_data = np.array(data_lines)
        
        print(f"✓ Loaded temperature-energy calibration data:")
        print(f"  Temperature range: {temp_data[:, temp_columns['temperature']].min():.1f} - {temp_data[:, temp_columns['temperature']].max():.1f} K")
        print(f"  {len(temp_data)} temperature points")
        
        # Create interpolators for U(T) - OLD SCHOOL LINEAR INTERPOLATION
        temperatures = temp_data[:, temp_columns['temperature']]
        
        calibration_data = {}
        
        # Map display names to column names
        energy_components = {
            'total_PE': 'total_PE_hartree',
            'harmonic': 'harmonic_hartree',
            'lj': 'lj_hartree',
            'coulombic': 'coulombic_hartree'
        }
        
        for display_name, col_name in energy_components.items():
            if col_name in temp_columns:
                energies = temp_data[:, temp_columns[col_name]]
                
                # Sort by temperature
                sort_idx = np.argsort(temperatures)
                T_sorted = temperatures[sort_idx]
                U_sorted = energies[sort_idx]
                
                # Create OLD SCHOOL linear interpolator U(T) with extrapolation
                U_of_T = interp1d(T_sorted, U_sorted, kind='linear', 
                                 bounds_error=False, fill_value='extrapolate')
                
                calibration_data[display_name] = {
                    'T': T_sorted,
                    'U': U_sorted,
                    'U_of_T': U_of_T,
                    'T_range': (T_sorted.min(), T_sorted.max()),
                    'U_range': (U_sorted.min(), U_sorted.max())
                }
                
                print(f"  Created U(T) interpolator for {display_name}")
        
        return calibration_data
        
    except Exception as e:
        print(f"Error loading calibration data: {e}")
        import traceback
        traceback.print_exc()
        return None

def load_fictive_temperature_data(filename):
    """
    Load fictive temperature time series data.
    
    Returns:
    --------
    time : np.array
        Time array
    fictive_temps : dict
        Dictionary with component names as keys, values are temperature arrays
    """
    try:
        # Read the file
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # Find the header line
        header_line = None
        data_start = 0
        
        for i, line in enumerate(lines):
            if line.strip().startswith('# time_ps'):
                header_line = line.strip()[2:]  # Remove '# '
                data_start = i + 1
                break
            elif not line.strip().startswith('#') and line.strip():
                # Header without # prefix
                header_line = line.strip()
                data_start = i + 1
                break
        
        if header_line is None:
            raise ValueError("Could not find header line")
        
        column_names = header_line.split()
        print(f"✓ Found columns: {column_names}")
        
        # Read data
        data_lines = []
        for line in lines[data_start:]:
            if line.strip() and not line.strip().startswith('#'):
                values = line.split()
                data_lines.append(values)
        
        if not data_lines:
            raise ValueError("No data found")
        
        # Convert to arrays
        time = []
        fictive_temps = {col: [] for col in column_names[1:]}  # Skip time_ps
        
        for values in data_lines:
            time.append(float(values[0]))
            for i, col_name in enumerate(column_names[1:], start=1):
                val = values[i]
                if val.lower() == 'nan':
                    fictive_temps[col_name].append(np.nan)
                else:
                    fictive_temps[col_name].append(float(val))
        
        time = np.array(time)
        for col_name in fictive_temps:
            fictive_temps[col_name] = np.array(fictive_temps[col_name])
        
        print(f"✓ Loaded fictive temperature data:")
        print(f"  Time range: {time[0]:.3f} - {time[-1]:.3f} ps")
        print(f"  Data points: {len(time)}")
        print(f"  Components: {list(fictive_temps.keys())}")
        
        return time, fictive_temps
        
    except Exception as e:
        print(f"Error loading fictive temperature data: {e}")
        import traceback
        traceback.print_exc()
        return None, None

def reconstruct_potential_energy(time, fictive_temps, calibration_data):
    """
    Reconstruct averaged potential energy from fictive temperatures.
    
    Parameters:
    -----------
    time : np.array
        Time array
    fictive_temps : dict
        Dictionary of fictive temperature arrays
    calibration_data : dict
        Calibration data with U(T) interpolators
        
    Returns:
    --------
    reconstructed_energies : dict
        Dictionary with component names as keys, values are energy arrays
    """
    reconstructed_energies = {}
    
    # Map fictive temperature column names to calibration component names
    component_mapping = {
        'fictive_T_total_PE_K': 'total_PE',
        'fictive_T_harmonic_K': 'harmonic',
        'fictive_T_lj_K': 'lj',
        'fictive_T_coulombic_K': 'coulombic'
    }
    
    for fictive_col, calib_comp in component_mapping.items():
        if fictive_col in fictive_temps and calib_comp in calibration_data:
            T_fictive = fictive_temps[fictive_col]
            U_of_T = calibration_data[calib_comp]['U_of_T']
            
            # Reconstruct energy: U = U_of_T(T_fictive)
            U_reconstructed = np.zeros_like(T_fictive)
            
            for i, T in enumerate(T_fictive):
                if np.isnan(T):
                    U_reconstructed[i] = np.nan
                else:
                    U_reconstructed[i] = U_of_T(T)
            
            reconstructed_energies[calib_comp] = U_reconstructed
            
            # Statistics
            valid_mask = ~np.isnan(U_reconstructed)
            n_valid = np.sum(valid_mask)
            n_total = len(U_reconstructed)
            
            print(f"  Reconstructed {calib_comp}: {n_valid}/{n_total} valid points")
            if n_valid > 0:
                print(f"    Energy range: {U_reconstructed[valid_mask].min():.6f} to {U_reconstructed[valid_mask].max():.6f} Hartree")
    
    return reconstructed_energies

def save_reconstructed_energy(time, reconstructed_energies, output_file, source_file):
    """Save reconstructed potential energy data to file."""
    
    output_lines = []
    output_lines.append("# Reconstructed Averaged Potential Energy Time Series")
    output_lines.append(f"# Reconstructed from fictive temperature data: {source_file}")
    output_lines.append("# Method: OLD SCHOOL linear interpolation (not extended fitting)")
    output_lines.append("# Units: Time in ps, Energy in Hartree")
    output_lines.append("#")
    output_lines.append("# Columns:")
    
    # Create header
    header_cols = ['time_ps']
    component_order = ['harmonic', 'lj', 'coulombic', 'total_PE']
    for comp in component_order:
        if comp in reconstructed_energies:
            header_cols.append(f'{comp}_energy_hartree')
    
    output_lines.append("# " + "\t".join(header_cols))
    
    # Write data
    for i, t in enumerate(time):
        data_row = [f"{t:.6f}"]
        
        for comp in component_order:
            if comp in reconstructed_energies:
                U = reconstructed_energies[comp][i]
                if np.isnan(U):
                    data_row.append("NaN")
                else:
                    data_row.append(f"{U:.12e}")
        
        output_lines.append("\t".join(data_row))
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write('\n'.join(output_lines))
    
    print(f"\n✓ Reconstructed potential energy saved to: {output_file}")
    print(f"  Time range: {time[0]:.3f} - {time[-1]:.3f} ps")
    print(f"  Data points: {len(time)}")
    print(f"  Components: {list(reconstructed_energies.keys())}")

def main():
    parser = argparse.ArgumentParser(
        description='Reconstruct averaged potential energy from fictive temperature data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('fictive_temp_file', type=str,
                       help='Input fictive temperature file')
    parser.add_argument('--temp-components', type=str, 
                       default='potential_energy_components_vs_temperature.txt',
                       help='Temperature-energy calibration file (default: potential_energy_components_vs_temperature.txt)')
    parser.add_argument('--output', '-o', type=str, default=None,
                       help='Output file for reconstructed energy (default: auto-generate from input name)')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not Path(args.fictive_temp_file).exists():
        print(f"Error: Input file not found: {args.fictive_temp_file}")
        sys.exit(1)
    
    if not Path(args.temp_components).exists():
        print(f"Error: Calibration file not found: {args.temp_components}")
        sys.exit(1)
    
    print("="*70)
    print("RECONSTRUCTING AVERAGED POTENTIAL ENERGY FROM FICTIVE TEMPERATURE")
    print("="*70)
    
    # Load calibration data
    print("\n[1/4] Loading temperature-energy calibration data...")
    calibration_data = load_temperature_energy_calibration(args.temp_components)
    if calibration_data is None:
        sys.exit(1)
    
    # Load fictive temperature data
    print("\n[2/4] Loading fictive temperature data...")
    time, fictive_temps = load_fictive_temperature_data(args.fictive_temp_file)
    if time is None or fictive_temps is None:
        sys.exit(1)
    
    # Reconstruct potential energy
    print("\n[3/4] Reconstructing potential energy time series...")
    reconstructed_energies = reconstruct_potential_energy(time, fictive_temps, calibration_data)
    
    if not reconstructed_energies:
        print("Error: No energies could be reconstructed!")
        sys.exit(1)
    
    # Determine output filename
    if args.output is None:
        input_path = Path(args.fictive_temp_file)
        # Replace 'fictive_temperatures' with 'potential_energy' in filename
        output_name = input_path.stem.replace('fictive_temperatures', 'potential_energy') + '.txt'
        output_path = input_path.parent / output_name
    else:
        output_path = Path(args.output)
    
    # Save results
    print("\n[4/4] Saving reconstructed potential energy...")
    save_reconstructed_energy(time, reconstructed_energies, output_path, args.fictive_temp_file)
    
    print("\n" + "="*70)
    print("✓ RECONSTRUCTION COMPLETE!")
    print("="*70)

if __name__ == '__main__':
    main()

