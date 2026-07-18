#!/usr/bin/env python3
"""
HDF5 Data Inspector for CavityMD Simulation Results
==================================================
This script provides comprehensive analysis and visualization of HDF5 simulation data.
"""

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for headless systems
import matplotlib.pyplot as plt
import argparse
import sys
from pathlib import Path

def print_dataset_summary(name, dataset, show_values=False, max_values=10):
    """Print summary information about a dataset."""
    print(f"\n📄 Dataset: {name}")
    print(f"   ├─ Shape: {dataset.shape}")
    print(f"   ├─ Data type: {dataset.dtype}")
    print(f"   ├─ Size: {dataset.nbytes / (1024**2):.3f} MB")
    
    if len(dataset.shape) == 1 and dataset.shape[0] > 0:
        data = dataset[:]
        print(f"   ├─ Min: {np.min(data):.6e}")
        print(f"   ├─ Max: {np.max(data):.6e}")
        print(f"   ├─ Mean: {np.mean(data):.6e}")
        print(f"   ├─ Std: {np.std(data):.6e}")
        
        if show_values:
            if len(data) <= max_values:
                print(f"   └─ All values: {data}")
            else:
                print(f"   ├─ First {max_values//2}: {data[:max_values//2]}")
                print(f"   └─ Last {max_values//2}: {data[-max_values//2:]}")
    else:
        print(f"   └─ Multi-dimensional data (not showing statistics)")

def print_group_structure(group, indent=0, show_values=False, max_values=10):
    """Recursively print the structure of HDF5 groups."""
    prefix = "  " * indent
    
    for key in group.keys():
        item = group[key]
        if isinstance(item, h5py.Group):
            print(f"{prefix}📁 Group: {key} ({len(item.keys())} items)")
            if len(item.keys()) > 0:
                print_group_structure(item, indent + 1, show_values, max_values)
            else:
                print(f"{prefix}   └─ (empty)")
        else:
            print_dataset_summary(f"{prefix[2:]}{key}", item, show_values, max_values)

def extract_time_series_data(filename):
    """Extract all time series data into a dictionary."""
    data_dict = {}
    
    with h5py.File(filename, 'r') as f:
        # Get time array
        if 'time' in f:
            data_dict['time'] = f['time'][:]
        
        # Recursively extract all datasets
        def extract_datasets(name, obj):
            if isinstance(obj, h5py.Dataset):
                # Clean up the name (remove leading slash)
                clean_name = name.lstrip('/')
                data_dict[clean_name] = obj[:]
        
        f.visititems(extract_datasets)
    
    return data_dict

def plot_energy_evolution(data_dict, output_dir="plots"):
    """Create plots of energy evolution."""
    Path(output_dir).mkdir(exist_ok=True)
    
    if 'time' not in data_dict:
        print("⚠️  No time data found, skipping plots")
        return
    
    time = data_dict['time']
    
    # Energy components plot
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('Energy Evolution', fontsize=16)
    
    # Cavity energies
    ax = axes[0, 0]
    cavity_keys = [k for k in data_dict.keys() if k.startswith('energies/cavity')]
    for key in cavity_keys:
        label = key.replace('energies/cavity_', '').replace('_', ' ').title()
        ax.plot(time, data_dict[key], label=label, linewidth=1.5)
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Energy (a.u.)')
    ax.set_title('Cavity Energies')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Molecular energies
    ax = axes[0, 1]
    molecular_keys = ['energies/lj', 'energies/harmonic', 'energies/ewald_short', 'energies/ewald_long']
    for key in molecular_keys:
        if key in data_dict:
            label = key.replace('energies/', '').replace('_', ' ').title()
            ax.plot(time, data_dict[key], label=label, linewidth=1.5)
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Energy (a.u.)')
    ax.set_title('Molecular Energies')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Total energies
    ax = axes[1, 0]
    total_keys = [k for k in data_dict.keys() if 'total' in k and k.startswith('energies/')]
    for key in total_keys:
        label = key.replace('energies/', '').replace('_', ' ').title()
        ax.plot(time, data_dict[key], label=label, linewidth=1.5)
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Energy (a.u.)')
    ax.set_title('Total Energies')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Temperature
    ax = axes[1, 1]
    if 'energies/temperature' in data_dict:
        ax.plot(time, data_dict['energies/temperature'], 'r-', linewidth=2, label='System Temperature')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Temperature (K)')
    ax.set_title('System Temperature')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/energy_evolution.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"✅ Saved energy evolution plot: {output_dir}/energy_evolution.png")

def plot_temperature_evolution(data_dict, output_dir="plots"):
    """Create plots of temperature evolution."""
    Path(output_dir).mkdir(exist_ok=True)
    
    if 'time' not in data_dict:
        return
    
    time = data_dict['time']
    
    # Temperature components
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    fig.suptitle('Temperature Evolution', fontsize=16)
    
    # All temperatures
    temp_keys = [k for k in data_dict.keys() if k.startswith('temperatures/') and not k.endswith('molecular')]
    colors = plt.cm.tab10(np.linspace(0, 1, len(temp_keys)))
    
    for i, key in enumerate(temp_keys):
        if key in data_dict:
            label = key.replace('temperatures/', '').replace('_', ' ').title()
            # Skip NaN values for plotting
            data = data_dict[key]
            valid_mask = ~np.isnan(data)
            if np.any(valid_mask):
                ax1.plot(time[valid_mask], data[valid_mask], 
                        color=colors[i], linewidth=2, label=label, marker='o', markersize=2)
    
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Temperature (K)')
    ax1.set_title('All Temperature Components')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Focus on key temperatures
    key_temps = ['temperatures/kinetic', 'temperatures/cavity_bath', 'temperatures/molecular_bath']
    for key in key_temps:
        if key in data_dict:
            label = key.replace('temperatures/', '').replace('_', ' ').title()
            data = data_dict[key]
            valid_mask = ~np.isnan(data)
            if np.any(valid_mask):
                ax2.plot(time[valid_mask], data[valid_mask], 
                        linewidth=2, label=label, marker='o', markersize=2)
    
    ax2.set_xlabel('Time (ps)')
    ax2.set_ylabel('Temperature (K)')
    ax2.set_title('Key Temperature Components')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/temperature_evolution.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"✅ Saved temperature evolution plot: {output_dir}/temperature_evolution.png")

def export_data_to_csv(data_dict, output_file="simulation_data.csv"):
    """Export all data to a CSV file."""
    import pandas as pd
    
    # Create DataFrame
    df = pd.DataFrame(data_dict)
    
    # Save to CSV
    df.to_csv(output_file, index=False)
    print(f"✅ Exported data to CSV: {output_file}")
    print(f"   └─ Shape: {df.shape}")
    print(f"   └─ Columns: {list(df.columns)}")

def print_data_arrays(data_dict, max_rows=20):
    """Print data arrays in a formatted way."""
    print("\n" + "="*80)
    print("📊 DATA ARRAYS CONTENT")
    print("="*80)
    
    for key, data in data_dict.items():
        print(f"\n🔍 {key}")
        print("-" * (len(key) + 4))
        
        if isinstance(data, np.ndarray):
            print(f"Shape: {data.shape}, Type: {data.dtype}")
            
            if len(data.shape) == 1:
                if len(data) <= max_rows:
                    for i, val in enumerate(data):
                        print(f"  [{i:4d}] {val:12.6e}")
                else:
                    print(f"  First {max_rows//2} values:")
                    for i in range(max_rows//2):
                        print(f"  [{i:4d}] {data[i]:12.6e}")
                    print("  ...")
                    print(f"  Last {max_rows//2} values:")
                    for i in range(-max_rows//2, 0):
                        print(f"  [{len(data)+i:4d}] {data[i]:12.6e}")
            else:
                print(f"  Multi-dimensional array (shape: {data.shape})")
                print(f"  First few elements: {data.flat[:min(10, data.size)]}")

def main():
    parser = argparse.ArgumentParser(description='Inspect HDF5 simulation data')
    parser.add_argument('filename', help='HDF5 file to inspect')
    parser.add_argument('--show-values', action='store_true', help='Show actual data values')
    parser.add_argument('--max-values', type=int, default=10, help='Maximum values to show per dataset')
    parser.add_argument('--plot', action='store_true', help='Generate plots')
    parser.add_argument('--export-csv', action='store_true', help='Export data to CSV')
    parser.add_argument('--print-arrays', action='store_true', help='Print all data arrays')
    parser.add_argument('--output-dir', default='plots', help='Output directory for plots')
    
    args = parser.parse_args()
    
    if not Path(args.filename).exists():
        print(f"❌ File not found: {args.filename}")
        sys.exit(1)
    
    print(f"🔍 Inspecting HDF5 file: {args.filename}")
    print("="*60)
    
    try:
        # Open and inspect file structure
        with h5py.File(args.filename, 'r') as f:
            print("📋 File Attributes:")
            for key, value in f.attrs.items():
                print(f"  {key}: {value}")
            
            print(f"\n📁 File Structure:")
            print_group_structure(f, show_values=args.show_values, max_values=args.max_values)
        
        # Extract data for further analysis
        print(f"\n📊 Extracting data...")
        data_dict = extract_time_series_data(args.filename)
        print(f"✅ Extracted {len(data_dict)} datasets")
        
        # Print arrays if requested
        if args.print_arrays:
            print_data_arrays(data_dict)
        
        # Generate plots if requested
        if args.plot:
            print(f"\n📈 Generating plots...")
            plot_energy_evolution(data_dict, args.output_dir)
            plot_temperature_evolution(data_dict, args.output_dir)
        
        # Export to CSV if requested
        if args.export_csv:
            print(f"\n💾 Exporting to CSV...")
            export_data_to_csv(data_dict)
        
        print(f"\n✅ Analysis complete!")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
