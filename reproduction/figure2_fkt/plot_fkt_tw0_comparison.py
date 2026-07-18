#!/usr/bin/env python3
"""
Plot F(k,t) at waiting times t_mathrm{w} = 0.0, 800, and 1400 ps for zero coupling vs 0.0007 a.u. coupling
Three panels stacked vertically
Matches the style of fkt_coupling_1x10neg3_filtered.pdf
Uses color scheme from normalized_relaxation_analysis_filtered.pdf
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import pandas as pd
from pathlib import Path
import subprocess
import sys
from pathlib import Path

REPRO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPRO_ROOT / "shared"))
from latex_config_adobe import setup_latex_fonts, latex_safe

def process_fkt_data(time, fkt, normalization_value=None):
    """
    Process F(k,t) data: normalize and remove values below 0.001.
    """
    time = np.array(time, dtype=np.float64)
    fkt = np.array(fkt, dtype=np.float64)

    # Remove NaN values and zero values
    valid_mask = ~(np.isnan(time) | np.isnan(fkt)) & (fkt != 0)
    time_clean = time[valid_mask]
    fkt_clean = fkt[valid_mask]
    
    if len(time_clean) < 2:
        return None, None, None
    
    # Sort by time to ensure monotonic behavior
    sort_indices = np.argsort(time_clean)
    time_sorted = time_clean[sort_indices]
    fkt_sorted = fkt_clean[sort_indices]
    
    # Normalize using provided normalization value
    if normalization_value is not None and normalization_value != 0:
        fkt_normalized = fkt_sorted / normalization_value
    else:
        # Fallback: normalize by first value
        if len(fkt_sorted) > 0:
            fkt_normalized = fkt_sorted / fkt_sorted[0]
        else:
            return None, None, None
    
    # Remove values below 0.001 after normalization
    above_threshold_mask = fkt_normalized >= 0.001
    if not np.any(above_threshold_mask):
        return None, None, None
    
    time_filtered = time_sorted[above_threshold_mask]
    fkt_filtered = fkt_normalized[above_threshold_mask]
    
    # Auto-detect maximum time for plotting
    max_time = time_filtered[-1] if len(time_filtered) > 0 else None
    
    return time_filtered, fkt_filtered, max_time

def read_fkt_data(file_path):
    """Read F(k,t) data from file."""
    try:
        data = pd.read_csv(file_path, sep=r'\s+', comment='#')
        time = data.iloc[:, 0].values
        fkt = data.iloc[:, 1].values
        return time, fkt
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None, None

def get_coupling_color(coupling_value, coupling_range=(0.0, 0.001)):
    """
    Get color for coupling strength from coolwarm colormap.
    Uses normalized position in the coupling range.
    """
    cmap = plt.colormaps.get_cmap('coolwarm')
    norm = Normalize(vmin=coupling_range[0], vmax=coupling_range[1])
    return cmap(norm(coupling_value))

def main():
    """Main plotting function."""
    # Set matplotlib style to classic and enable LaTeX
    plt.style.use('classic')
    try:
        # Test if LaTeX is available by trying to compile a simple expression
        subprocess.check_output(['latex', '--version'], stderr=subprocess.DEVNULL)
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Computer Modern Roman']
        plt.rcParams['font.size'] = 10
        # Add preamble to load amsmath for \boldsymbol and ensure proper font rendering
        plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb}'
        # Ensure PDF font embedding
        plt.rcParams['pdf.fonttype'] = 42  # TrueType fonts
        plt.rcParams['ps.fonttype'] = 42   # TrueType fonts
        use_latex = True
        print("Using LaTeX for mathematical typesetting with Computer Modern fonts")
    except (subprocess.CalledProcessError, FileNotFoundError, ImportError):
        print("Warning: LaTeX not available, using default matplotlib rendering")
        plt.rcParams['text.usetex'] = False
        plt.rcParams['font.family'] = 'DejaVu Sans'
        use_latex = False
    
    print("Creating F(k,t) comparison plot for t_mathrm{w} = 0.0, 800, and 1400 ps (ref1, ref4, ref7)")
    print("=" * 60)
    
    base_dir = Path('.')
    
    # File paths for the two coupling strengths and three waiting times
    # t_mathrm{w} = 0.0 ps (ref1)
    zero_coupling_ref1 = base_dir / 'cavity_coupling_0epos00_switch_200.0ps' / 'master_fskt_ref1.txt'
    coupling_7e4_ref1 = base_dir / 'cavity_coupling_7eneg04_switch_200.0ps' / 'master_fskt_ref1.txt'
    
    # ref4 (t_mathrm{w} = 4*200 ps = 800 ps)
    zero_coupling_ref4 = base_dir / 'cavity_coupling_0epos00_switch_200.0ps' / 'master_fskt_ref4.txt'
    coupling_7e4_ref4 = base_dir / 'cavity_coupling_7eneg04_switch_200.0ps' / 'master_fskt_ref4.txt'
    
    # ref7 (t_mathrm{w} = 7*200 ps = 1400 ps)
    zero_coupling_ref7 = base_dir / 'cavity_coupling_0epos00_switch_200.0ps' / 'master_fskt_ref7.txt'
    coupling_7e4_ref7 = base_dir / 'cavity_coupling_7eneg04_switch_200.0ps' / 'master_fskt_ref7.txt'
    
    # Check if files exist
    for file_path in [zero_coupling_ref1, coupling_7e4_ref1, zero_coupling_ref4, coupling_7e4_ref4, 
                      zero_coupling_ref7, coupling_7e4_ref7]:
        if not file_path.exists():
            print(f"ERROR: File not found: {file_path}")
            return 1
    
    # Read data for ref1 (t_mathrm{w} = 0.0 ps)
    print(f"\nReading ref1 (t_mathrm{{w}} = 0.0 ps) data...")
    time_zero_ref1, fkt_zero_ref1 = read_fkt_data(zero_coupling_ref1)
    time_7e4_ref1, fkt_7e4_ref1 = read_fkt_data(coupling_7e4_ref1)
    
    # Read data for ref4 (t_mathrm{w} = 800 ps)
    print(f"Reading ref4 (t_mathrm{{w}} = 800 ps) data...")
    time_zero_ref4, fkt_zero_ref4 = read_fkt_data(zero_coupling_ref4)
    time_7e4_ref4, fkt_7e4_ref4 = read_fkt_data(coupling_7e4_ref4)
    
    # Read data for ref7 (t_mathrm{w} = 1400 ps)
    print(f"Reading ref7 (t_mathrm{{w}} = 1400 ps) data...")
    time_zero_ref7, fkt_zero_ref7 = read_fkt_data(zero_coupling_ref7)
    time_7e4_ref7, fkt_7e4_ref7 = read_fkt_data(coupling_7e4_ref7)
    
    if any(x is None for x in [time_zero_ref1, fkt_zero_ref1, time_7e4_ref1, fkt_7e4_ref1,
                                 time_zero_ref4, fkt_zero_ref4, time_7e4_ref4, fkt_7e4_ref4,
                                 time_zero_ref7, fkt_zero_ref7, time_7e4_ref7, fkt_7e4_ref7]):
        print("ERROR: Could not read one or more F(k,t) data files")
        return 1
    
    # Process data for ref1
    print("\nProcessing ref1 data...")
    norm_zero_ref1 = fkt_zero_ref1[fkt_zero_ref1 != 0][0] if np.any(fkt_zero_ref1 != 0) else 1.0
    norm_7e4_ref1 = fkt_7e4_ref1[fkt_7e4_ref1 != 0][0] if np.any(fkt_7e4_ref1 != 0) else 1.0
    
    time_zero_ref1_proc, fkt_zero_ref1_proc, max_time_zero_ref1 = process_fkt_data(time_zero_ref1, fkt_zero_ref1, norm_zero_ref1)
    time_7e4_ref1_proc, fkt_7e4_ref1_proc, max_time_7e4_ref1 = process_fkt_data(time_7e4_ref1, fkt_7e4_ref1, norm_7e4_ref1)
    
    # Process data for ref4
    print("Processing ref4 data...")
    norm_zero_ref4 = fkt_zero_ref4[fkt_zero_ref4 != 0][0] if np.any(fkt_zero_ref4 != 0) else 1.0
    norm_7e4_ref4 = fkt_7e4_ref4[fkt_7e4_ref4 != 0][0] if np.any(fkt_7e4_ref4 != 0) else 1.0
    
    time_zero_ref4_proc, fkt_zero_ref4_proc, max_time_zero_ref4 = process_fkt_data(time_zero_ref4, fkt_zero_ref4, norm_zero_ref4)
    time_7e4_ref4_proc, fkt_7e4_ref4_proc, max_time_7e4_ref4 = process_fkt_data(time_7e4_ref4, fkt_7e4_ref4, norm_7e4_ref4)
    
    # Process data for ref7
    print("Processing ref7 data...")
    norm_zero_ref7 = fkt_zero_ref7[fkt_zero_ref7 != 0][0] if np.any(fkt_zero_ref7 != 0) else 1.0
    norm_7e4_ref7 = fkt_7e4_ref7[fkt_7e4_ref7 != 0][0] if np.any(fkt_7e4_ref7 != 0) else 1.0
    
    time_zero_ref7_proc, fkt_zero_ref7_proc, max_time_zero_ref7 = process_fkt_data(time_zero_ref7, fkt_zero_ref7, norm_zero_ref7)
    time_7e4_ref7_proc, fkt_7e4_ref7_proc, max_time_7e4_ref7 = process_fkt_data(time_7e4_ref7, fkt_7e4_ref7, norm_7e4_ref7)
    
    if any(x is None for x in [time_zero_ref1_proc, fkt_zero_ref1_proc, time_7e4_ref1_proc, fkt_7e4_ref1_proc,
                                time_zero_ref4_proc, fkt_zero_ref4_proc, time_7e4_ref4_proc, fkt_7e4_ref4_proc,
                                time_zero_ref7_proc, fkt_zero_ref7_proc, time_7e4_ref7_proc, fkt_7e4_ref7_proc]):
        print("ERROR: Processed F(k,t) data is empty after filtering")
        return 1
    
    print(f"  Ref1 - Zero coupling: {len(time_zero_ref1_proc)} points")
    print(f"  Ref1 - 0.0007 coupling: {len(time_7e4_ref1_proc)} points")
    print(f"  Ref4 - Zero coupling: {len(time_zero_ref4_proc)} points")
    print(f"  Ref4 - 0.0007 coupling: {len(time_7e4_ref4_proc)} points")
    print(f"  Ref7 - Zero coupling: {len(time_zero_ref7_proc)} points")
    print(f"  Ref7 - 0.0007 coupling: {len(time_7e4_ref7_proc)} points")
    
    # Get colors from coolwarm colormap based on coupling strength
    # Using the range from the normalized_relaxation_analysis_filtered.pdf
    color_zero = get_coupling_color(0.0, coupling_range=(0.0, 0.001))
    color_7e4 = get_coupling_color(0.0007, coupling_range=(0.0, 0.001))
    
    print(f"\nColor for zero coupling: {color_zero}")
    print(f"Color for 0.0007 coupling: {color_7e4}")
    
    # Create the plot with three panels stacked vertically
    print("\nCreating plot...")
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(2.2, 5))
    
    # Panel 1: t_mathrm{w} = 0.0 ps (ref1) - NO X-AXIS LABEL
    ax1.plot(time_zero_ref1_proc, fkt_zero_ref1_proc, 
            color=color_zero, linewidth=2, alpha=0.8, 
            label=r'\textrm{Out of cavity}')
    ax1.plot(time_7e4_ref1_proc, fkt_7e4_ref1_proc, 
            color=color_7e4, linewidth=2, alpha=0.8, 
            label=r'\textrm{In cavity}')
    ax1.set_ylabel(r'$\phi_{k}(t; t_{\mathrm{w}})$')
    ax1.set_title(r'$t_{\mathrm{w}} = 0.0$ \textrm{ps}')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8, loc='best')
    ax1.set_ylim(bottom=0.001)
    ax1.set_xlim(0, 1000)
    ax1.tick_params(labelbottom=False)  # Remove x-axis tick labels
    
    # Panel 2: ref4 (t_mathrm{w} = 800 ps) - NO X-AXIS LABEL
    ax2.plot(time_zero_ref4_proc, fkt_zero_ref4_proc, 
            color=color_zero, linewidth=2, alpha=0.8, 
            label=r'\textrm{Out of cavity}')
    ax2.plot(time_7e4_ref4_proc, fkt_7e4_ref4_proc, 
            color=color_7e4, linewidth=2, alpha=0.8, 
            label=r'\textrm{In cavity}')
    ax2.set_ylabel(r'$\phi_{k}(t; t_{\mathrm{w}})$')
    ax2.set_title(r'$t_{\mathrm{w}} = 800$ \textrm{ps}')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=8, loc='best')
    ax2.set_ylim(bottom=0.001)
    ax2.set_xlim(0, 1000)
    ax2.tick_params(labelbottom=False)  # Remove x-axis tick labels
    
    # Panel 3: ref7 (t_mathrm{w} = 1400 ps) - WITH X-AXIS LABEL
    ax3.plot(time_zero_ref7_proc, fkt_zero_ref7_proc, 
            color=color_zero, linewidth=2, alpha=0.8, 
            label=r'\textrm{Out of cavity}')
    ax3.plot(time_7e4_ref7_proc, fkt_7e4_ref7_proc, 
            color=color_7e4, linewidth=2, alpha=0.8, 
            label=r'\textrm{In cavity}')
    ax3.set_xlabel(r'$t-t_{\mathrm{w}}$ \textrm{(ps)}')
    ax3.set_ylabel(r'$\phi_{k}(t; t_{\mathrm{w}})$')
    ax3.set_title(r'$t_{\mathrm{w}} = 1400$ \textrm{ps}')
    ax3.grid(True, alpha=0.3)
    ax3.legend(fontsize=8, loc='best')
    ax3.set_ylim(bottom=0.001)
    ax3.set_xlim(0, 1000)
    
    # Tight layout
    plt.tight_layout()
    
    # Save plot
    output_file = 'fkt_tw_comparison.pdf'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved: {output_file}")
    
    # Also save as PNG
    output_png = 'fkt_tw_comparison.png'
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_png}")
    
    plt.close()
    
    print("\n" + "=" * 60)
    print("Plot generation complete!")
    
    return 0

if __name__ == "__main__":
    exit(main())

