#!/usr/bin/env python3
"""
Individual F(k,t) Coupling Plots Script

This script creates individual F(k,t) vs time plots for each coupling strength,
showing different waiting times (tw) on each plot with a colorbar. This is similar to the 
fkt_by_coupling_filtered.png functionality but creates separate files instead of
a multi-panel plot.

Each plot shows F(k,t) vs time for all waiting times at a specific coupling strength.

Output files: fkt_coupling_<coupling>_filtered.pdf (e.g., fkt_coupling_1x10neg3_filtered.pdf)

Usage: python plot_individual_fkt_coupling.py [--base_dir ./] [--output_dir ./]
"""

import argparse
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for remote servers
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path
import re
from matplotlib.colors import Normalize

REPRO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPRO_ROOT / "shared"))

# Import the functions we need from the main analysis script
from plot_fkt_analysis import (
    collect_data, 
    filter_coupling_strengths, 
    process_fkt_data,
    parse_coupling_strength
)

# Import improved LaTeX configuration for Adobe Illustrator compatibility
from latex_config_adobe import setup_latex_fonts, latex_safe, format_coupling, format_time, format_fkt

def create_coupling_strength_label(coupling_value):
    """Create a properly formatted coupling strength label."""
    if coupling_value == 0:
        return '0'
    elif coupling_value >= 1e-3:
        return f'{coupling_value:.0e}'
    else:
        # Format in scientific notation like 7.0 × 10⁻⁴
        exp_str = f'{coupling_value:.1e}'
        mantissa, exp_part = exp_str.split('e')
        exponent = int(exp_part)
        return f'{mantissa} × 10^{{{exponent}}}'

def plot_individual_fkt_coupling(data_dict, coupling_info, output_dir='.'):
    """
    Create individual F(k,t) vs time plots for each coupling strength.
    
    Each plot shows F(k,t) decay curves for all waiting times at a specific coupling strength.
    
    Parameters:
    -----------
    data_dict : dict
        Dictionary containing F(k,t) data for each coupling
    coupling_info : dict
        Information about coupling strengths
    output_dir : str
        Output directory for plots
    """
    
    print("Creating individual F(k,t) coupling plots...")
    
    # Set up LaTeX fonts for Adobe Illustrator compatibility
    use_latex = setup_latex_fonts()
    if use_latex:
        print("  Using LaTeX with Computer Modern fonts (Adobe Illustrator compatible)")
    else:
        print("  Using fallback fonts")
    
    # Filter to specific coupling strengths
    filtered_coupling_info = filter_coupling_strengths(coupling_info)
    
    # Get all available reference times
    all_refs = set()
    for coupling_name in filtered_coupling_info.keys():
        if coupling_name in data_dict:
            all_refs.update(data_dict[coupling_name].keys())
    
    all_refs = sorted(all_refs)
    
    if not all_refs:
        print("No reference data to plot")
        return []
    
    # Set up color mapping for waiting times
    ref_time_values = [ref * 200 for ref in all_refs]  # Convert to ps
    norm = Normalize(vmin=min(ref_time_values), vmax=max(ref_time_values))
    cmap = plt.colormaps.get_cmap('viridis')
    
    generated_files = []
    
    # Create individual plot for each coupling strength
    for coupling_name in filtered_coupling_info.keys():
        if coupling_name not in data_dict:
            continue
            
        coupling_value, coupling_label = filtered_coupling_info[coupling_name]
        print(f"  Creating plot for {coupling_label}...")
        
        # Create individual figure with smaller size and space for colorbar
        fig, ax = plt.subplots(1, 1, figsize=(4, 3))
        
        has_data = False
        max_time_for_plot = 0
        
        # Get normalization value from ref0 of this coupling
        normalization_value = None
        if 0 in data_dict[coupling_name]:
            time_ref0, fkt_ref0 = data_dict[coupling_name][0]
            # Find first non-zero value from ref0 for normalization
            nonzero_mask = fkt_ref0 != 0
            if np.any(nonzero_mask):
                first_nonzero_idx = np.where(nonzero_mask)[0][0]
                normalization_value = fkt_ref0[first_nonzero_idx]
        
        # Plot F(k,t) for each waiting time at this coupling strength
        for ref_num in all_refs:
            if ref_num in data_dict[coupling_name]:
                time, fkt = data_dict[coupling_name][ref_num]
                ref_time_ps = ref_num * 200
                
                color = cmap(norm(ref_time_ps))
                
                # Process F(k,t) data with filtering and normalization
                time_processed, fkt_processed, max_time = process_fkt_data(time, fkt, normalization_value)
                
                if time_processed is not None and fkt_processed is not None:
                    ax.plot(time_processed, fkt_processed, color=color, linewidth=2, alpha=0.8)
                    has_data = True
                    
                    if max_time is not None:
                        max_time_for_plot = max(max_time_for_plot, max_time)
        
        if not has_data:
            ax.text(0.5, 0.5, 'No data available', ha='center', va='center', 
                   transform=ax.transAxes, fontsize=16)
        
        # Format the plot with larger fonts
        ax.set_xlabel(latex_safe('$t-t_{\\mathrm{w}}$ (ps)', 't-t_w (ps)', use_latex), fontsize=16)
        ax.set_ylabel(latex_safe('$F_k(t+t_{\\mathrm{w}},t_{\\mathrm{w}})$', 'F_k(t+t_w,t_w)', use_latex), fontsize=16)
        
        # Title with coupling strength
        if coupling_value == 0:
            title_text = latex_safe('$\\varepsilon = 0$ a.u.', 'ε = 0 a.u.', use_latex)
        else:
            exp_str = f'{coupling_value:.0e}'
            mantissa, exp_part = exp_str.split('e')
            exponent = int(exp_part)
            if use_latex:
                title_text = f'$\\varepsilon = {mantissa} \\times 10^{{{exponent}}}$ a.u.'
            else:
                title_text = f'ε = {mantissa} × 10^{exponent} a.u.'
        
        ax.set_title(title_text, fontsize=18, fontweight='bold')
        
        # Add dashed grid lines
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.tick_params(axis='x', which='major', rotation=0, pad=5)  # Prevent overlap
        
        # Set axis limits
        ax.set_ylim(bottom=0.001)  # Start y-axis at threshold
        ax.set_xlim(0, 1600)  # Set x-axis maximum to 1600 ps
        
        # Set x-axis ticks every 400 ps
        ax.xaxis.set_major_locator(ticker.MultipleLocator(400))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(200))
        
        # Set y-axis ticks every 0.2 units
        ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
        
        # Add colorbar for waiting times
        if has_data:
            from matplotlib.colorbar import ColorbarBase
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            
            # Create smaller colorbar
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="3%", pad=0.1)
            
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            
            cbar = plt.colorbar(sm, cax=cax)
            cbar_label = latex_safe('$t_{\\mathrm{w}}$ (ps)', 't_w (ps)', use_latex)
            cbar.set_label(cbar_label, fontsize=14, rotation=90, labelpad=15)
            cbar.ax.tick_params(labelsize=12)
            
            # Set colorbar ticks to match x-axis frequency (every 400 ps)
            # Create tick positions at multiples of 400 within the colorbar range
            cb_min, cb_max = norm.vmin, norm.vmax
            tick_positions = []
            tick_val = 0
            while tick_val <= cb_max:
                if tick_val >= cb_min:
                    tick_positions.append(tick_val)
                tick_val += 400
            cbar.set_ticks(tick_positions)
        
        plt.tight_layout()
        
        # Create filename based on coupling strength
        if coupling_value == 0:
            filename_base = 'fkt_coupling_0_filtered'
        else:
            # Format for filename: 1x10neg3 for 1e-3
            exp_str = f'{coupling_value:.0e}'
            mantissa, exp_part = exp_str.split('e')
            exponent = int(exp_part)
            if exponent < 0:
                exp_suffix = f'{abs(exponent)}'
                filename_base = f'fkt_coupling_{mantissa}x10neg{exp_suffix}_filtered'
            else:
                filename_base = f'fkt_coupling_{mantissa}x10pos{exponent}_filtered'
        
        # Save both PNG and PDF versions
        output_file_png = Path(output_dir) / f'{filename_base}.png'
        plt.savefig(output_file_png, dpi=300, bbox_inches='tight')
        
        output_file_pdf = Path(output_dir) / f'{filename_base}.pdf'
        plt.savefig(output_file_pdf, bbox_inches='tight')
        
        print(f"    Saved: {output_file_png}")
        print(f"    Saved: {output_file_pdf} (Adobe Illustrator compatible)")
        
        generated_files.append(output_file_png.name)
        generated_files.append(output_file_pdf.name)
        
        plt.close()
    
    return generated_files

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Create individual F(k,t) plots for each coupling strength')
    parser.add_argument('--base_dir', default='./', 
                       help='Base directory containing coupling folders (default: ./)')
    parser.add_argument('--output_dir', default='./', 
                       help='Output directory for plots (default: ./)')
    
    args = parser.parse_args()
    
    print("Individual F(k,t) Coupling Plots Script")
    print("=" * 50)
    print(f"Base directory: {Path(args.base_dir).resolve()}")
    print(f"Output directory: {Path(args.output_dir).resolve()}")
    
    # Collect all data
    print("\nCollecting F(k,t) data...")
    data_dict, coupling_info = collect_data(args.base_dir)
    
    if not data_dict:
        print("No data found! Make sure you have cavity_coupling_* folders with master_fskt_ref*.txt files.")
        return 1
    
    print(f"\nSummary:")
    print(f"  Found {len(data_dict)} coupling strengths")
    
    total_refs = set()
    for folder_data in data_dict.values():
        total_refs.update(folder_data.keys())
    
    print(f"  Total unique references: {sorted(total_refs)}")
    
    # Create output directory if it doesn't exist
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    # Generate individual plots
    generated_files = plot_individual_fkt_coupling(data_dict, coupling_info, args.output_dir)
    
    print("\n" + "=" * 50)
    print("Analysis complete! Generated individual F(k,t) plots:")
    
    # Group by reference for cleaner output
    png_files = [f for f in generated_files if f.endswith('.png')]
    pdf_files = [f for f in generated_files if f.endswith('.pdf')]
    
    print(f"\nPNG files ({len(png_files)}):")
    for f in png_files:
        print(f"  📊 {f}")
    
    print(f"\nPDF files ({len(pdf_files)}):")
    for f in pdf_files:
        print(f"  📊 {f}")
    
    print(f"\nTotal files generated: {len(generated_files)}")
    
    return 0

if __name__ == "__main__":
    exit(main())
