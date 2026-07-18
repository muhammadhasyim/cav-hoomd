#!/usr/bin/env python3
"""
Debug script to analyze bond topology in GSD files.
This script will help identify the "same particle bonded to itself" issue.
"""

import gsd.hoomd
import numpy as np
import argparse
import sys

def analyze_bond_topology(gsd_filename, frame_index=-1):
    """
    Analyze bond topology in a GSD file to identify potential issues.
    
    Args:
        gsd_filename (str): Path to GSD file
        frame_index (int): Frame index to analyze (-1 for last frame)
    """
    print(f"Analyzing bond topology in: {gsd_filename}")
    print(f"Frame index: {frame_index}")
    print("="*60)
    
    try:
        with gsd.hoomd.open(gsd_filename, 'r') as f:
            n_frames = len(f)
            print(f"Total frames in file: {n_frames}")
            
            if frame_index < 0:
                frame_index = n_frames + frame_index
            
            if frame_index < 0 or frame_index >= n_frames:
                print(f"ERROR: Invalid frame index {frame_index}")
                return False
            
            frame = f[frame_index]
            print(f"Analyzing frame {frame_index}")
            print()
            
            # Analyze particles
            print("PARTICLE ANALYSIS:")
            print(f"Total particles: {frame.particles.N}")
            print(f"Particle types: {frame.particles.types}")
            print(f"Particle type IDs: {np.unique(frame.particles.typeid)}")
            
            # Count particles by type
            for i, ptype in enumerate(frame.particles.types):
                count = np.sum(frame.particles.typeid == i)
                print(f"  {ptype} particles (type {i}): {count}")
            print()
            
            # Analyze bonds
            print("BOND ANALYSIS:")
            if not hasattr(frame, 'bonds') or frame.bonds.N == 0:
                print("No bonds found in this frame")
                return True
            
            print(f"Total bonds: {frame.bonds.N}")
            if hasattr(frame.bonds, 'types'):
                print(f"Bond types: {frame.bonds.types}")
            print()
            
            # Check for self-bonding (particle bonded to itself)
            print("CHECKING FOR SELF-BONDING:")
            bond_groups = frame.bonds.group
            bond_types = frame.bonds.typeid if hasattr(frame.bonds, 'typeid') else np.zeros(len(bond_groups))
            
            self_bonds = []
            for i, (p1, p2) in enumerate(bond_groups):
                if p1 == p2:
                    bond_type = bond_types[i] if i < len(bond_types) else 'unknown'
                    self_bonds.append((i, p1, p2, bond_type))
                    print(f"  SELF-BOND FOUND: Bond {i} connects particle {p1} to itself (type {bond_type})")
            
            if self_bonds:
                print(f"\n*** PROBLEM IDENTIFIED: {len(self_bonds)} self-bonds found! ***")
                print("Self-bonds are invalid - a particle cannot be bonded to itself.")
                return False
            else:
                print("  No self-bonds found - bond topology looks valid")
            
            print()
            
            # Analyze bond statistics
            print("BOND STATISTICS:")
            unique_bond_types = np.unique(bond_types)
            for bond_type in unique_bond_types:
                count = np.sum(bond_types == bond_type)
                if hasattr(frame.bonds, 'types') and bond_type < len(frame.bonds.types):
                    type_name = frame.bonds.types[bond_type]
                    print(f"  Bond type {bond_type} ({type_name}): {count} bonds")
                else:
                    print(f"  Bond type {bond_type}: {count} bonds")
            
            # Check for out-of-range particle indices
            print("\nCHECKING FOR OUT-OF-RANGE PARTICLE INDICES:")
            max_particle_id = frame.particles.N - 1
            invalid_bonds = []
            
            for i, (p1, p2) in enumerate(bond_groups):
                if p1 < 0 or p1 > max_particle_id or p2 < 0 or p2 > max_particle_id:
                    bond_type = bond_types[i] if i < len(bond_types) else 'unknown'
                    invalid_bonds.append((i, p1, p2, bond_type))
                    print(f"  INVALID BOND: Bond {i} references non-existent particle(s) {p1}-{p2} (type {bond_type})")
            
            if invalid_bonds:
                print(f"\n*** PROBLEM IDENTIFIED: {len(invalid_bonds)} bonds with invalid particle indices! ***")
                return False
            else:
                print("  All bond particle indices are valid")
            
            print()
            print("SUMMARY:")
            if self_bonds or invalid_bonds:
                print("❌ BOND TOPOLOGY IS INVALID")
                print("   This will cause HOOMD to fail with bond topology errors")
                return False
            else:
                print("✅ BOND TOPOLOGY APPEARS VALID")
                return True
                
    except Exception as e:
        print(f"ERROR analyzing file: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Debug bond topology in GSD files")
    parser.add_argument('gsd_file', help='Path to GSD file to analyze')
    parser.add_argument('--frame', type=int, default=-1, help='Frame index to analyze (default: -1 for last frame)')
    
    args = parser.parse_args()
    
    success = analyze_bond_topology(args.gsd_file, args.frame)
    
    if not success:
        print("\nRECOMMENDATIONS:")
        print("1. Check if the GSD file was corrupted during transfer")
        print("2. Try using a different frame from the GSD file")
        print("3. Regenerate the GSD file from the original source")
        print("4. Check HOOMD version compatibility between systems")
        sys.exit(1)
    else:
        print("\nThe bond topology looks valid. The error might be caused by:")
        print("1. Different HOOMD versions between systems")
        print("2. Different compilation options or floating point precision")
        print("3. Race conditions if multiple processes access the same file")

if __name__ == "__main__":
    main() 