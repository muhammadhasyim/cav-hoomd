#!/usr/bin/env python3
"""
Quick HDF5 Data Inspector - Simple version for immediate data viewing
"""

import h5py
import numpy as np
import sys

def print_group_structure(group, indent=0):
    """Print the structure of HDF5 groups and datasets."""
    prefix = "  " * indent
    
    for key in group.keys():
        item = group[key]
        if isinstance(item, h5py.Group):
            print(f"{prefix}📁 Group: {key} ({len(item.keys())} items)")
            if len(item.keys()) > 0:
                print_group_structure(item, indent + 1)
            else:
                print(f"{prefix}   └─ (empty)")
        else:
            shape = item.shape
            dtype = item.dtype
            size_mb = item.nbytes / (1024 * 1024)
            print(f"{prefix}📄 Dataset: {key}")
            print(f"{prefix}   ├─ Shape: {shape}")
            print(f"{prefix}   ├─ Type: {dtype}")
            print(f"{prefix}   └─ Size: {size_mb:.3f} MB")

def quick_inspect(filename):
    """Quick inspection of HDF5 file with data arrays printed."""
    
    print(f"🔍 Quick inspection of: {filename}")
    print("="*60)
    
    with h5py.File(filename, 'r') as f:
        # File info
        print("📋 File Info:")
        for key, value in f.attrs.items():
            print(f"  {key}: {value}")
        
        # Data structure overview
        print(f"\n📁 Data Structure:")
        print_group_structure(f)
        
        # Summary statistics
        total_datasets = 0
        total_size_mb = 0
        max_length = 0
        
        def count_datasets(name, obj):
            nonlocal total_datasets, total_size_mb, max_length
            if isinstance(obj, h5py.Dataset):
                total_datasets += 1
                total_size_mb += obj.nbytes / (1024 * 1024)
                if len(obj.shape) > 0:
                    max_length = max(max_length, obj.shape[0])
        
        f.visititems(count_datasets)
        
        print(f"\n📊 Summary:")
        print(f"  Total datasets: {total_datasets}")
        print(f"  Total data size: {total_size_mb:.2f} MB")
        print(f"  Maximum data length: {max_length} points")
        
        # Time data first
        if 'time' in f:
            time_data = f['time'][:]
            print(f"\n⏰ Time Data (shape: {time_data.shape}):")
            print(f"  Range: {time_data[0]:.6f} to {time_data[-1]:.6f} ps")
            if len(time_data) > 1:
                print(f"  Step: {time_data[1] - time_data[0]:.6f} ps")
                total_time = time_data[-1] - time_data[0]
                print(f"  Total simulation time: {total_time:.3f} ps")
            print(f"  First 10: {time_data[:10]}")
            print(f"  Last 10:  {time_data[-10:]}")
        
        # Energy data
        print(f"\n🔋 Energy Data:")
        energy_keys = [k for k in f['energies'].keys()] if 'energies' in f else []
        for key in sorted(energy_keys)[:5]:  # Show first 5 energy components
            data = f[f'energies/{key}'][:]
            print(f"\n  📊 {key}:")
            print(f"    Range: {np.min(data):.6e} to {np.max(data):.6e}")
            print(f"    Mean: {np.mean(data):.6e}, Std: {np.std(data):.6e}")
            print(f"    First 5: {data[:5]}")
            print(f"    Last 5:  {data[-5:]}")
        
        if len(energy_keys) > 5:
            print(f"    ... and {len(energy_keys)-5} more energy components")
        
        # Temperature data
        print(f"\n🌡️  Temperature Data:")
        temp_keys = [k for k in f['temperatures'].keys() if isinstance(f[f'temperatures/{k}'], h5py.Dataset)] if 'temperatures' in f else []
        for key in sorted(temp_keys):
            data = f[f'temperatures/{key}'][:]
            valid_data = data[~np.isnan(data)]
            if len(valid_data) > 0:
                print(f"\n  📊 {key}:")
                print(f"    Range: {np.min(valid_data):.2f} to {np.max(valid_data):.2f} K")
                print(f"    Mean: {np.mean(valid_data):.2f} K, Std: {np.std(valid_data):.2f} K")
                print(f"    First 5: {data[:5]}")
                print(f"    Last 5:  {data[-5:]}")
            else:
                print(f"\n  📊 {key}: (all NaN values)")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python quick_inspect.py <hdf5_file>")
        sys.exit(1)
    
    filename = sys.argv[1]
    try:
        quick_inspect(filename)
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)
