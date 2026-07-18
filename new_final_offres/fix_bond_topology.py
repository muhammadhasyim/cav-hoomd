#!/usr/bin/env python3
"""
Script to fix bond topology issues in GSD files.
This script removes invalid bonds (self-bonds and out-of-range particle indices).
"""

import gsd.hoomd
import numpy as np
import argparse
import sys
import shutil

def fix_bond_topology(input_gsd, output_gsd, frame_index=-1, verbose=True):
    """
    Fix bond topology in a GSD file by removing invalid bonds.
    
    Args:
        input_gsd (str): Path to input GSD file
        output_gsd (str): Path to output GSD file (can be same as input)
        frame_index (int): Frame index to fix (-1 for all frames)
        verbose (bool): Print detailed information
        
    Returns:
        bool: True if successful, False otherwise
    """
    if verbose:
        print(f"Fixing bond topology in: {input_gsd}")
        print(f"Output file: {output_gsd}")
        print("="*60)
    
    try:
        # Read the GSD file
        with gsd.hoomd.open(input_gsd, 'r') as f_in:
            n_frames = len(f_in)
            
            if frame_index >= 0:
                frames_to_fix = [frame_index]
            else:
                frames_to_fix = list(range(n_frames))
            
            if verbose:
                print(f"Total frames in file: {n_frames}")
                print(f"Frames to fix: {len(frames_to_fix)}")
            
            # Create output file
            with gsd.hoomd.open(output_gsd, 'w') as f_out:
                fixed_frames = 0
                total_bonds_removed = 0
                
                for frame_idx in range(n_frames):
                    frame = f_in[frame_idx]
                    
                    if frame_idx in frames_to_fix:
                        # Fix this frame
                        if verbose and len(frames_to_fix) <= 10:
                            print(f"\nFixing frame {frame_idx}:")
                        
                        # Check if frame has bonds
                        if hasattr(frame, 'bonds') and frame.bonds.N > 0:
                            bond_groups = frame.bonds.group
                            bond_types = frame.bonds.typeid if hasattr(frame.bonds, 'typeid') else np.zeros(len(bond_groups))
                            
                            # Find invalid bonds
                            max_particle_id = frame.particles.N - 1
                            valid_bonds = []
                            valid_bond_types = []
                            removed_count = 0
                            
                            for i, (p1, p2) in enumerate(bond_groups):
                                # Check for self-bonds
                                if p1 == p2:
                                    if verbose and len(frames_to_fix) <= 10:
                                        print(f"  Removing self-bond: {p1}-{p2}")
                                    removed_count += 1
                                    continue
                                
                                # Check for out-of-range particle indices
                                if p1 < 0 or p1 > max_particle_id or p2 < 0 or p2 > max_particle_id:
                                    if verbose and len(frames_to_fix) <= 10:
                                        print(f"  Removing out-of-range bond: {p1}-{p2}")
                                    removed_count += 1
                                    continue
                                
                                # Bond is valid
                                valid_bonds.append([p1, p2])
                                if i < len(bond_types):
                                    valid_bond_types.append(bond_types[i])
                                else:
                                    valid_bond_types.append(0)  # Default bond type
                            
                            # Update frame with valid bonds only
                            if len(valid_bonds) < len(bond_groups):
                                if verbose and len(frames_to_fix) <= 10:
                                    print(f"  Removed {removed_count} invalid bonds")
                                    print(f"  Kept {len(valid_bonds)} valid bonds")
                                
                                frame.bonds.N = len(valid_bonds)
                                if len(valid_bonds) > 0:
                                    frame.bonds.group = np.array(valid_bonds, dtype=np.uint32)
                                    frame.bonds.typeid = np.array(valid_bond_types, dtype=np.uint32)
                                else:
                                    # No bonds left - need to handle this carefully
                                    frame.bonds.group = np.empty((0, 2), dtype=np.uint32)
                                    frame.bonds.typeid = np.empty(0, dtype=np.uint32)
                                
                                fixed_frames += 1
                                total_bonds_removed += removed_count
                            elif verbose and len(frames_to_fix) <= 10:
                                print(f"  No invalid bonds found")
                        elif verbose and len(frames_to_fix) <= 10:
                            print(f"  No bonds in this frame")
                    
                    # Write frame to output
                    f_out.append(frame)
                
                if verbose:
                    print(f"\nSUMMARY:")
                    print(f"  Fixed {fixed_frames} frames")
                    print(f"  Removed {total_bonds_removed} invalid bonds total")
                    print(f"  Output written to: {output_gsd}")
                
                return True
                
    except Exception as e:
        print(f"ERROR fixing bond topology: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    parser = argparse.ArgumentParser(description="Fix bond topology issues in GSD files")
    parser.add_argument('input_gsd', help='Path to input GSD file')
    parser.add_argument('-o', '--output', help='Path to output GSD file (default: overwrite input)')
    parser.add_argument('--frame', type=int, default=-1, help='Frame index to fix (-1 for all frames)')
    parser.add_argument('--backup', action='store_true', help='Create backup of input file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # Set output file
    if args.output:
        output_gsd = args.output
    else:
        output_gsd = args.input_gsd
    
    # Create backup if requested
    if args.backup and args.input_gsd == output_gsd:
        backup_file = args.input_gsd + '.backup'
        if args.verbose:
            print(f"Creating backup: {backup_file}")
        shutil.copy2(args.input_gsd, backup_file)
    
    # Fix the bond topology
    success = fix_bond_topology(args.input_gsd, output_gsd, args.frame, args.verbose)
    
    if success:
        print("\n✅ Bond topology fixed successfully!")
        if args.verbose:
            print("\nTo verify the fix, run:")
            print(f"python debug_bond_topology.py {output_gsd}")
        sys.exit(0)
    else:
        print("\n❌ Failed to fix bond topology!")
        sys.exit(1)

if __name__ == "__main__":
    main() 