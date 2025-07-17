#!/usr/bin/env python3
"""
Script to analyze cavity particle displacement from GSD trajectory.

This script reads a GSD trajectory file and plots the unwrapped cavity particle
positions over time to visualize the finite-q displacement when coupling switches on.
"""

import numpy as np
import matplotlib.pyplot as plt
import gsd.hoomd
import sys
import os
import argparse
from pathlib import Path

def unwrap_positions(positions, images, box_lengths):
    """Unwrap particle positions across periodic boundaries."""
    pos = np.asarray(positions)
    img = np.asarray(images)
    box = np.asarray(box_lengths)
    return pos + img * box

def find_cavity_particle(snapshot):
    """Find the cavity particle index (type 'L', typeid=2)."""
    cavity_indices = np.where(snapshot.particles.typeid == 2)[0]
    if len(cavity_indices) == 0:
        raise ValueError("No cavity particle found in snapshot")
    return cavity_indices[0]

def analyze_cavity_trajectory(gsd_file, switch_time_ps=1.0):
    """Analyze cavity particle trajectory and detect displacement."""
    
    print(f"Analyzing cavity trajectory from: {gsd_file}")
    
    # Read trajectory
    with gsd.hoomd.open(gsd_file, 'r') as f:
        n_frames = len(f)
        print(f"Number of frames: {n_frames}")
        
        # Find cavity particle in first frame
        cavity_idx = find_cavity_particle(f[0])
        print(f"Cavity particle found at index: {cavity_idx}")
        
        # Extract data for all frames
        times = []
        positions = []
        box_lengths = []
        
        for i, frame in enumerate(f):
            # Get simulation time (assuming log data contains elapsed time)
            if hasattr(frame, 'log') and 'Time/elapsed_ps' in frame.log:
                time_ps = float(frame.log['Time/elapsed_ps'])
            else:
                # Fallback: estimate from frame number (assuming fixed timestep)
                # This is approximate and may not be accurate for adaptive timestep
                time_ps = i * 0.1  # Rough estimate
            
            times.append(time_ps)
            
            # Get cavity particle position and image
            pos = frame.particles.position[cavity_idx]
            img = frame.particles.image[cavity_idx]
            box = frame.configuration.box[:3]
            
            # Unwrap position
            unwrapped_pos = unwrap_positions([pos], [img], box)[0]
            
            positions.append(unwrapped_pos)
            box_lengths.append(box)
        
        times = np.array(times)
        positions = np.array(positions)
        box_lengths = np.array(box_lengths)
    
    # Detect displacement event
    displacement_magnitude = np.linalg.norm(positions - positions[0], axis=1)
    significant_displacement = displacement_magnitude > 1.0  # Threshold for significant displacement
    
    if np.any(significant_displacement):
        displacement_frame = np.where(significant_displacement)[0][0]
        displacement_time = float(times[displacement_frame])  # Ensure scalar
        displacement_pos = positions[displacement_frame]
        initial_pos = positions[0]
        
        print(f"\nDisplacement detected:")
        print(f"  Frame: {displacement_frame}")
        print(f"  Time: {displacement_time:.6f} ps (target: {switch_time_ps:.6f} ps)")
        print(f"  Initial position: {initial_pos}")
        print(f"  Final position: {displacement_pos}")
        print(f"  Displacement: {displacement_pos - initial_pos}")
        print(f"  Magnitude: {np.linalg.norm(displacement_pos - initial_pos):.6f}")
    else:
        print("\nNo significant displacement detected")
        displacement_frame = None
        displacement_time = None
    
    return times, positions, displacement_frame, displacement_time

def plot_cavity_positions(times, positions, displacement_frame=None, displacement_time=None, 
                         switch_time_ps=1.0, output_file=None):
    """Plot cavity particle positions over time."""
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot x, y, z positions
    axes[0, 0].plot(times, positions[:, 0], 'b-', linewidth=1.5, label='x position')
    axes[0, 0].set_xlabel('Time (ps)')
    axes[0, 0].set_ylabel('x position (σ)')
    axes[0, 0].set_title('Cavity Particle x Position')
    axes[0, 0].grid(True, alpha=0.3)
    
    axes[0, 1].plot(times, positions[:, 1], 'r-', linewidth=1.5, label='y position')
    axes[0, 1].set_xlabel('Time (ps)')
    axes[0, 1].set_ylabel('y position (σ)')
    axes[0, 1].set_title('Cavity Particle y Position')
    axes[0, 1].grid(True, alpha=0.3)
    
    axes[1, 0].plot(times, positions[:, 2], 'g-', linewidth=1.5, label='z position')
    axes[1, 0].set_xlabel('Time (ps)')
    axes[1, 0].set_ylabel('z position (σ)')
    axes[1, 0].set_title('Cavity Particle z Position')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Plot displacement magnitude
    displacement_magnitude = np.linalg.norm(positions - positions[0], axis=1)
    axes[1, 1].plot(times, displacement_magnitude, 'k-', linewidth=1.5, label='Displacement magnitude')
    axes[1, 1].set_xlabel('Time (ps)')
    axes[1, 1].set_ylabel('Displacement magnitude (σ)')
    axes[1, 1].set_title('Cavity Particle Displacement Magnitude')
    axes[1, 1].grid(True, alpha=0.3)
    
    # Mark switch time and displacement time
    for ax in axes.flat:
        # Mark expected switch time
        ax.axvline(x=switch_time_ps, color='orange', linestyle='--', alpha=0.7, 
                  label=f'Expected switch ({switch_time_ps:.1f} ps)')
        
        # Mark actual displacement time if detected
        if displacement_time is not None:
            ax.axvline(x=displacement_time, color='red', linestyle='-', alpha=0.8, 
                      label=f'Actual displacement ({displacement_time:.3f} ps)')
    
    # Add legends
    for ax in axes.flat:
        ax.legend()
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {output_file}")
    
    plt.show()

def main():
    """Main function to analyze cavity displacement."""
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Analyze cavity particle displacement from GSD trajectory')
    parser.add_argument('--dir', dest='directory', type=str, 
                       help='Specific cavity coupling directory to analyze')
    args = parser.parse_args()
    
    # Look for GSD files in the cavity coupling directory
    if args.directory:
        # Use specified directory
        cavity_dir = Path(args.directory)
        if not cavity_dir.exists():
            print(f"Specified directory not found: {cavity_dir}")
            sys.exit(1)
        print(f"Using specified directory: {cavity_dir}")
    else:
        # Auto-detect directories
        cavity_dirs = [d for d in Path('.').glob('cavity_coupling_*') if d.is_dir()]
        
        if not cavity_dirs:
            print("No cavity coupling directories found!")
            print("Make sure you're running this script from the examples directory")
            sys.exit(1)
        
        # Use the most recent directory
        cavity_dir = sorted(cavity_dirs)[-1]
        print(f"Auto-selected directory: {cavity_dir}")
    
    gsd_file = cavity_dir / 'prod-0.gsd'
    
    if not gsd_file.exists():
        print(f"GSD file not found: {gsd_file}")
        sys.exit(1)
    
    # Extract switch time from directory name if possible
    switch_time_ps = 1.0  # Default
    if 'switch' in cavity_dir.name:
        parts = cavity_dir.name.split('switch_')
        if len(parts) > 1:
            switch_part = parts[1].replace('ps', '')
            if switch_part.replace('.', '').isdigit():
                switch_time_ps = float(switch_part)
                print(f"Extracted switch time: {switch_time_ps} ps")
            else:
                print(f"Could not parse switch time from '{switch_part}', using default: {switch_time_ps} ps")
        else:
            print(f"Directory contains 'switch' but couldn't extract time, using default: {switch_time_ps} ps")
    
    # Analyze trajectory
    times, positions, displacement_frame, displacement_time = analyze_cavity_trajectory(
        gsd_file, switch_time_ps=switch_time_ps
    )
    
    # Plot results
    output_file = cavity_dir / 'cavity_displacement_analysis.png'
    plot_cavity_positions(times, positions, displacement_frame, displacement_time, 
                         switch_time_ps, output_file=output_file)
    
    # Print summary
    print(f"\nAnalysis Summary:")
    print(f"  Total simulation time: {times[-1]:.3f} ps")
    print(f"  Number of frames: {len(times)}")
    print(f"  Expected switch time: {switch_time_ps:.3f} ps")
    if displacement_time is not None:
        print(f"  Actual displacement time: {displacement_time:.3f} ps")
        print(f"  Timing error: {abs(displacement_time - switch_time_ps):.3f} ps")
        print(f"  Displacement successful: {'✓' if abs(displacement_time - switch_time_ps) < 0.1 else '✗'}")
    else:
        print(f"  No displacement detected - check if simulation completed")

if __name__ == '__main__':
    main() 