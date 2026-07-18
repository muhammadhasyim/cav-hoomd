#!/usr/bin/env python3
"""
Recalculate Fictive Temperatures Using Empirical U(T) Relationships

This script recalculates fictive temperatures from energy tracker files using
proper empirical U(T) relationships instead of simple equipartition or other
approximations. This ensures that the fictive temperature plots use the same
sophisticated methods as plot_temperature_feedback.py.

The script:
1. Loads empirical energy-temperature calibration data from potential_energy_components_vs_temperature.txt
2. Creates proper U(T) -> T interpolators for each energy component
3. Loads raw energy data from coupling directories
4. Calculates fictive temperatures using the empirical relationships
5. Saves results as new averaged_fictive_temperatures_empirical.txt files

Usage:
    python recalculate_fictive_temperatures_empirical.py [--base_dir .] [--output_dir time_series_output]
"""

import numpy as np
import pandas as pd
from pathlib import Path
import argparse
from scipy import interpolate
import warnings

def load_temperature_energy_components(filename="potential_energy_components_vs_temperature.txt"):
    """Load empirical energy-temperature calibration data."""
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        header_line = None
        data_start_line = 0
        
        for i, line in enumerate(lines):
            line = line.strip()
            if not line.startswith('#') and line:
                header_line = line
                data_start_line = i + 1
                break
        
        if header_line is None:
            raise ValueError("Could not find column headers in temperature components file")
        
        column_names = header_line.split()
        temp_columns = {name: i for i, name in enumerate(column_names)}
        temp_data = np.loadtxt(filename, skiprows=data_start_line)
        
        return temp_data, temp_columns
        
    except Exception as e:
        print(f"Error: Could not load temperature components file {filename}: {e}")
        return None, None

def create_empirical_interpolators(temp_data, temp_columns):
    """Create T(U) interpolators from empirical calibration data."""
    if temp_data is None or temp_columns is None:
        return {}
    
    interpolators = {}
    temperatures = temp_data[:, temp_columns['temperature']]
    
    energy_components = {
        'harmonic': 'harmonic_hartree',
        'lj': 'lj_hartree', 
        'coulombic': 'coulombic_hartree',
        'total': 'total_PE_hartree'
    }
    
    for comp_name, col_name in energy_components.items():
        if col_name in temp_columns:
            energies = temp_data[:, temp_columns[col_name]]
            
            # Sort by energy for proper interpolation
            sort_idx = np.argsort(energies)
            U_sorted = energies[sort_idx]
            T_sorted = temperatures[sort_idx]
            
            # Create T(U) interpolator with extrapolation
            T_of_U = interpolate.interp1d(U_sorted, T_sorted, kind='linear', 
                                        bounds_error=False, fill_value='extrapolate')
            
            interpolators[comp_name] = {
                'T_of_U': T_of_U,
                'U_range': (U_sorted.min(), U_sorted.max()),
                'T_range': (T_sorted.min(), T_sorted.max())
            }
            
            print(f"  {comp_name}: U range = [{U_sorted.min():.6f}, {U_sorted.max():.6f}] Hartree")
            print(f"            T range = [{T_sorted.min():.1f}, {T_sorted.max():.1f}] K")
    
    return interpolators

def calculate_empirical_fictive_temperatures(energy_dict, interpolators):
    """Calculate fictive temperatures using empirical U(T) relationships."""
    fictive_temps = {}
    
    for comp_name, energies in energy_dict.items():
        if comp_name in interpolators:
            interp_data = interpolators[comp_name]
            T_of_U = interp_data['T_of_U']
            U_range = interp_data['U_range']
            
            # Calculate fictive temperatures
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fictive_T = T_of_U(energies)
            
            # Mask physically unreasonable temperatures
            reasonable_mask = (fictive_T > 0) & (fictive_T < 10000)
            fictive_T[~reasonable_mask] = np.nan
            
            # Count interpolated vs extrapolated points
            interpolated_mask = (energies >= U_range[0]) & (energies <= U_range[1])
            n_interpolated = np.sum(interpolated_mask & reasonable_mask)
            n_extrapolated = np.sum(~interpolated_mask & reasonable_mask)
            n_invalid = np.sum(~reasonable_mask)
            
            print(f"    {comp_name}: {n_interpolated} interpolated, {n_extrapolated} extrapolated, {n_invalid} invalid")
            
            fictive_temps[comp_name] = fictive_T
        else:
            print(f"    Warning: No interpolator for {comp_name}")
            fictive_temps[comp_name] = np.full_like(energies, np.nan)
    
    return fictive_temps

def load_and_average_energy_data(coupling_dir):
    """Load and average energy tracker data from a coupling directory."""
    coupling_path = Path(coupling_dir)
    energy_files = list(coupling_path.glob('prod-*_energy_tracker.txt'))
    
    if not energy_files:
        return None, None
    
    print(f"  Found {len(energy_files)} energy tracker files")
    
    all_energy_data = []
    for energy_file in energy_files:
        try:
            data = pd.read_csv(energy_file, sep=r'\s+', comment='#')
            all_energy_data.append(data)
        except Exception as e:
            print(f"    Warning: Could not read {energy_file}: {e}")
            continue
    
    if not all_energy_data:
        return None, None
    
    # Use shortest time range for consistency
    min_length = min(len(df) for df in all_energy_data)
    common_time = all_energy_data[0]['time(ps)'].iloc[:min_length].values
    
    # Average energy components across replicas
    energy_components = {}
    required_columns = ['harmonic_energy', 'lj_energy', 'ewald_short_energy', 'ewald_long_energy']
    
    for col in required_columns:
        if all(col in df.columns for df in all_energy_data):
            values = np.array([df[col].iloc[:min_length].values for df in all_energy_data])
            energy_components[col] = np.mean(values, axis=0)
    
    # Create combined components for fictive temperature calculation
    if 'harmonic_energy' in energy_components:
        energy_components['harmonic'] = energy_components['harmonic_energy']
    
    if 'lj_energy' in energy_components:
        energy_components['lj'] = energy_components['lj_energy']
    
    if 'ewald_short_energy' in energy_components and 'ewald_long_energy' in energy_components:
        energy_components['coulombic'] = energy_components['ewald_short_energy'] + energy_components['ewald_long_energy']
    
    if all(comp in energy_components for comp in ['harmonic_energy', 'lj_energy', 'ewald_short_energy', 'ewald_long_energy']):
        energy_components['total'] = (energy_components['harmonic_energy'] + 
                                     energy_components['lj_energy'] + 
                                     energy_components['ewald_short_energy'] + 
                                     energy_components['ewald_long_energy'])
    
    print(f"  Averaged {len(common_time)} time points from {len(all_energy_data)} replicas")
    
    return common_time, energy_components

def parse_coupling_strength(dir_name):
    """Parse coupling strength from directory name."""
    coupling_part = dir_name.replace('cavity_coupling_', '').split('_switch_')[0]
    
    if '0epos00' in coupling_part:
        return 0.0, '0.0'
    elif 'eneg' in coupling_part:
        parts = coupling_part.split('eneg')
        if len(parts) == 2:
            try:
                mantissa = int(parts[0])
                exponent = int(parts[1])
                value = mantissa * (10 ** (-exponent))
                return value, f'{mantissa}\\times10^{{-{exponent}}}'
            except ValueError:
                return float('nan'), coupling_part
    
    return float('nan'), coupling_part

def apply_sliding_window_filter(time, values, window_ps=5.0):
    """Apply sliding window filter to match the original processing."""
    if len(time) < 2:
        return values
    
    dt = np.mean(np.diff(time))
    window_points = max(1, int(window_ps / dt))
    
    if window_points <= 3:
        return values
    
    from scipy import ndimage
    
    values_array = np.array(values)
    nan_mask = np.isnan(values_array)
    
    if np.all(nan_mask):
        return values_array
    
    # Interpolate NaN values for filtering
    if np.any(nan_mask):
        valid_indices = np.where(~nan_mask)[0]
        if len(valid_indices) < 2:
            return values_array
        
        values_interp = np.copy(values_array)
        values_interp[nan_mask] = np.interp(
            np.where(nan_mask)[0], 
            valid_indices, 
            values_array[valid_indices]
        )
    else:
        values_interp = values_array
    
    # Apply uniform filter
    filtered = ndimage.uniform_filter1d(values_interp, size=window_points, mode='nearest')
    
    # Restore NaN values
    if np.any(nan_mask):
        filtered[nan_mask] = np.nan
    
    return filtered

def save_empirical_fictive_temperatures(time, fictive_temps, coupling_value, output_file):
    """Save recalculated fictive temperatures to file."""
    
    # Apply sliding window filter (5 ps window as in original files)
    filtered_fictive_temps = {}
    for comp_name, temps in fictive_temps.items():
        if not np.all(np.isnan(temps)):
            filtered_fictive_temps[comp_name] = apply_sliding_window_filter(time, temps, window_ps=5.0)
        else:
            filtered_fictive_temps[comp_name] = temps
    
    # Create header
    header_lines = [
        "# Empirically Recalculated Fictive Temperatures from Energy Components",
        f"# Coupling strength: {coupling_value:.6e} (atomic units)",
        "# Generated using empirical U(T) relationships from potential_energy_components_vs_temperature.txt",
        "# Filter window: 5.0 ps sliding window (matching original processing)",
        "# Method: Empirical interpolation (NOT equipartition approximation)",
        "#",
        "# Column 1: time_ps | Column 2: fictive_T_coulombic_K | Column 3: fictive_T_harmonic_K | Column 4: fictive_T_lj_K | Column 5: fictive_T_total_PE_K",
        "#",
        "time_ps\tfictive_T_coulombic_K\tfictive_T_harmonic_K\tfictive_T_lj_K\tfictive_T_total_PE_K"
    ]
    
    # Write data
    with open(output_file, 'w') as f:
        for line in header_lines:
            f.write(line + '\n')
        
        for i, t in enumerate(time):
            coulombic = filtered_fictive_temps.get('coulombic', np.full_like(time, np.nan))[i]
            harmonic = filtered_fictive_temps.get('harmonic', np.full_like(time, np.nan))[i]
            lj = filtered_fictive_temps.get('lj', np.full_like(time, np.nan))[i]
            total = filtered_fictive_temps.get('total', np.full_like(time, np.nan))[i]
            
            f.write(f"{t:.6f}\t{coulombic:.6f}\t{harmonic:.6f}\t{lj:.6f}\t{total:.6f}\n")
    
    print(f"  Saved: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Recalculate fictive temperatures using empirical U(T) relationships',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('--base_dir', default='.', 
                       help='Base directory containing coupling directories (default: .)')
    parser.add_argument('--output_dir', default='time_series_output', 
                       help='Output directory for empirical fictive temperature files (default: time_series_output)')
    parser.add_argument('--empirical_data_file', default='potential_energy_components_vs_temperature.txt',
                       help='Path to empirical energy-temperature calibration data')
    
    args = parser.parse_args()
    
    print("Empirical Fictive Temperature Recalculation Script")
    print("=" * 60)
    print(f"Base directory: {Path(args.base_dir).resolve()}")
    print(f"Output directory: {Path(args.output_dir).resolve()}")
    print(f"Empirical data file: {args.empirical_data_file}")
    
    # Load empirical calibration data
    print("\nLoading empirical energy-temperature relationships...")
    temp_data, temp_columns = load_temperature_energy_components(args.empirical_data_file)
    
    if temp_data is None:
        print("Error: Could not load empirical data. Exiting.")
        return 1
    
    # Create interpolators
    print("Creating empirical U(T) -> T interpolators...")
    interpolators = create_empirical_interpolators(temp_data, temp_columns)
    
    if not interpolators:
        print("Error: Could not create interpolators. Exiting.")
        return 1
    
    print(f"Successfully created interpolators for: {list(interpolators.keys())}")
    
    # Find coupling directories
    base_path = Path(args.base_dir)
    coupling_dirs = [d for d in base_path.iterdir() if d.is_dir() and d.name.startswith('cavity_coupling_')]
    
    if not coupling_dirs:
        print(f"Error: No coupling directories found in {args.base_dir}")
        return 1
    
    print(f"\nFound {len(coupling_dirs)} coupling directories")
    
    # Create output directory
    output_path = Path(args.output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Process each coupling directory
    processed_count = 0
    for coupling_dir in sorted(coupling_dirs):
        dir_name = coupling_dir.name
        coupling_value, coupling_label = parse_coupling_strength(dir_name)
        
        print(f"\nProcessing: {dir_name} (ε = {coupling_label})")
        
        # Load energy data
        time, energy_components = load_and_average_energy_data(coupling_dir)
        
        if time is None or not energy_components:
            print(f"  Warning: Could not load energy data from {coupling_dir}")
            continue
        
        # Calculate empirical fictive temperatures
        print(f"  Calculating empirical fictive temperatures...")
        fictive_temps = calculate_empirical_fictive_temperatures(energy_components, interpolators)
        
        if not fictive_temps:
            print(f"  Warning: Could not calculate fictive temperatures")
            continue
        
        # Save results
        output_filename = f"coupling_{dir_name.replace('cavity_coupling_', '').split('_switch_')[0]}_averaged_fictive_temperatures_empirical.txt"
        output_file = output_path / output_filename
        
        save_empirical_fictive_temperatures(time, fictive_temps, coupling_value, output_file)
        processed_count += 1
    
    print(f"\n" + "=" * 60)
    print(f"Recalculation complete! Processed {processed_count} coupling strengths")
    print(f"Empirical fictive temperature files saved to: {output_path}")
    print("\nThese files use proper empirical U(T) relationships and can be used")
    print("with plot_fictive_temperature_components.py for accurate analysis.")
    
    return 0

if __name__ == "__main__":
    exit(main())
