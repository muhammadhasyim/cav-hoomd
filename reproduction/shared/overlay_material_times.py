#!/usr/bin/env python3
"""
Overlay rescaled fictive material time with F(k,t) material time.

This script creates a beautiful overlay plot showing:
1. ξ_fictive (raw) 
2. ξ_fictive / ln(10) (rescaled)
3. ξ_F(k,t) (reference)

For all coupling strengths to demonstrate the universal 1/ln(10) scaling.

Author: Assistant
Date: September 12, 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
import subprocess
from compare_material_time_scales import load_fictive_data, load_fkt_material_time, calculate_fictive_material_time

def create_overlay_plot():
    """Create overlay plot of all material times."""
    
    # Directories
    base_dir = '/media/extradrive/Trajectories/final_nodiss_cavitymd'
    time_series_dir = Path(base_dir) / 'time_series_output'
    material_time_file = os.path.join(base_dir, 'material_time_all_couplings.txt')
    output_dir = base_dir
    
    # All available coupling values
    coupling_values = [0.0, 3e-4, 5e-4, 7e-4, 1e-3]
    
    # Theoretical scale factor
    ln10_factor = 1.0 / np.log(10)
    
    print("=" * 70)
    print("CREATING OVERLAY PLOT: FICTIVE vs F(k,t) MATERIAL TIMES")
    print("=" * 70)
    print(f"Using theoretical scale factor: 1/ln(10) = {ln10_factor:.6f}")
    
    # Set up LaTeX rendering
    plt.style.use('classic')
    use_latex = True
    try:
        subprocess.check_output(['latex', '--version'], stderr=subprocess.DEVNULL)
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Computer Modern Roman']
        plt.rcParams['mathtext.fontset'] = 'cm'
        plt.rcParams['text.latex.preamble'] = r'\usepackage{lmodern}'
        print("  Using LaTeX with Computer Modern fonts")
    except (subprocess.CalledProcessError, FileNotFoundError, ImportError):
        print("  Warning: LaTeX not available, using default matplotlib rendering")
        plt.rcParams['text.usetex'] = False
        plt.rcParams['font.family'] = 'DejaVu Sans'
        use_latex = False
    
    # Helper function for LaTeX formatting
    def latex_format(text, fallback_text=None):
        if use_latex:
            return text
        else:
            return fallback_text if fallback_text else text.replace('$', '').replace(r'\xi', 'ξ').replace(r'\varepsilon', 'ε')
    
    # Create figure with two panels: individual plots and overlay
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Color scheme - same as F(k,t) analysis
    from matplotlib.colors import Normalize
    coupling_norm = Normalize(vmin=0, vmax=max(coupling_values))
    coupling_cmap = plt.colormaps.get_cmap('coolwarm')
    
    def get_coupling_color(coupling_value):
        """Get color for a given coupling value."""
        return coupling_cmap(coupling_norm(coupling_value))
    
    # Process each coupling strength
    all_data = {}
    
    for coupling_value in coupling_values:
        print(f"\nProcessing ε = {coupling_value}...")
        
        # Load data
        fictive_data = load_fictive_data(time_series_dir, coupling_value)
        fkt_data = load_fkt_material_time(material_time_file, coupling_value)
        
        if fictive_data is None or fkt_data is None:
            print(f"   Could not load data for ε = {coupling_value}")
            continue
        
        # Calculate fictive material time
        time_fictive = fictive_data['time']
        tau_fictive = fictive_data['fictive_relaxation_time']
        xi_fictive = calculate_fictive_material_time(time_fictive, tau_fictive)
        
        # Apply 1/ln(10) scaling
        xi_fictive_scaled = xi_fictive * ln10_factor
        
        # Store data for plotting
        all_data[coupling_value] = {
            'time_fictive': time_fictive,
            'xi_fictive_raw': xi_fictive,
            'xi_fictive_scaled': xi_fictive_scaled,
            'time_fkt': fkt_data['time'],
            'xi_fkt': fkt_data['material_time'],
            'tau_fictive': tau_fictive
        }
        
        # Calculate comparison statistics
        time_common = np.linspace(max(time_fictive[0], fkt_data['time'][0]), 
                                 min(time_fictive[-1], fkt_data['time'][-1]), 500)
        
        xi_fictive_interp = np.interp(time_common, time_fictive, xi_fictive_scaled)
        xi_fkt_interp = np.interp(time_common, fkt_data['time'], fkt_data['material_time'])
        
        relative_difference = (xi_fictive_interp - xi_fkt_interp) / xi_fkt_interp * 100
        mean_rel_diff = np.mean(np.abs(relative_difference))
        
        print(f"   Mean relative difference with 1/ln(10) scaling: {mean_rel_diff:.2f}%")
    
    # Panel 1: Individual comparison for ε = 0 (cleanest case)
    print(f"\nCreating detailed comparison for ε = 0...")
    
    if 0.0 in all_data:
        data_eps0 = all_data[0.0]
        color_eps0 = get_coupling_color(0.0)
        
        # Plot all three curves for ε = 0
        ax1.plot(data_eps0['time_fictive'], data_eps0['xi_fictive_raw'], 
                'b-', linewidth=2, alpha=0.7, label=latex_format(r'$\xi_{\mathrm{fictive}}$ (raw)', 'ξ_fictive (raw)'))
        
        ax1.plot(data_eps0['time_fictive'], data_eps0['xi_fictive_scaled'], 
                'r--', linewidth=2, alpha=0.8, 
                label=latex_format(r'$\xi_{\mathrm{fictive}} / \ln(10)$ (scaled)', 'ξ_fictive / ln(10) (scaled)'))
        
        ax1.plot(data_eps0['time_fkt'], data_eps0['xi_fkt'], 
                'g-', linewidth=2, alpha=0.8, 
                label=latex_format(r'$\xi_{\mathrm{F(k,t)}}$ (reference)', 'ξ_F(k,t) (reference)'))
        
        ax1.set_xlabel(latex_format(r'Time (ps)'))
        ax1.set_ylabel(latex_format(r'Material Time $\xi$', 'Material Time ξ'))
        ax1.set_title(latex_format(r'Material Time Comparison ($\varepsilon = 0$)', 'Material Time Comparison (ε = 0)'))
        ax1.legend(loc='upper left')
        ax1.grid(True, alpha=0.3)
        
        # Add text box with statistics
        time_common = np.linspace(max(data_eps0['time_fictive'][0], data_eps0['time_fkt'][0]), 
                                 min(data_eps0['time_fictive'][-1], data_eps0['time_fkt'][-1]), 500)
        
        xi_fictive_interp = np.interp(time_common, data_eps0['time_fictive'], data_eps0['xi_fictive_scaled'])
        xi_fkt_interp = np.interp(time_common, data_eps0['time_fkt'], data_eps0['xi_fkt'])
        
        mean_abs_diff = np.mean(np.abs(xi_fictive_interp - xi_fkt_interp))
        mean_rel_diff = np.mean(np.abs((xi_fictive_interp - xi_fkt_interp) / xi_fkt_interp * 100))
        
        stats_text = f'1/ln(10) = {ln10_factor:.4f}\nMean |diff| = {mean_abs_diff:.4f}\nMean rel diff = {mean_rel_diff:.2f}%'
        ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Panel 2: Overlay of all coupling strengths (scaled fictive vs F(k,t))
    print(f"\nCreating overlay plot for all coupling strengths...")
    
    for coupling_value in coupling_values:
        if coupling_value not in all_data:
            continue
        
        data = all_data[coupling_value]
        color = get_coupling_color(coupling_value)
        
        # Format coupling label
        if coupling_value == 0.0:
            label_fictive = latex_format(r'$\varepsilon = 0$ (fictive)', 'ε = 0 (fictive)')
            label_fkt = latex_format(r'$\varepsilon = 0$ (F(k,t))', 'ε = 0 (F(k,t))')
        else:
            if coupling_value >= 1e-3:
                exp_str = f'{coupling_value:.0e}'
                label_fictive = f'$\\varepsilon = {exp_str}$ (fictive)'
                label_fkt = f'$\\varepsilon = {exp_str}$ (F(k,t))'
            else:
                scaled_val = coupling_value * 10000
                label_fictive = f'$\\varepsilon = {scaled_val:.0f}\\times10^{{-4}}$ (fictive)'
                label_fkt = f'$\\varepsilon = {scaled_val:.0f}\\times10^{{-4}}$ (F(k,t))'
            
            label_fictive = latex_format(label_fictive)
            label_fkt = latex_format(label_fkt)
        
        # Plot scaled fictive material time (solid lines)
        ax2.plot(data['time_fictive'], data['xi_fictive_scaled'], 
                '-', color=color, linewidth=2, alpha=0.8, label=label_fictive)
        
        # Plot F(k,t) material time (dashed lines)
        ax2.plot(data['time_fkt'], data['xi_fkt'], 
                '--', color=color, linewidth=2, alpha=0.8, label=label_fkt)
    
    ax2.set_xlabel(latex_format(r'Time (ps)'))
    ax2.set_ylabel(latex_format(r'Material Time $\xi$', 'Material Time ξ'))
    ax2.set_title(latex_format(r'Universal Scaling: $\xi_{\mathrm{F(k,t)}} = \xi_{\mathrm{fictive}} / \ln(10)$', 
                              'Universal Scaling: ξ_F(k,t) = ξ_fictive / ln(10)'))
    
    # Create organized legend for panel 2
    handles, labels = ax2.get_legend_handles_labels()
    
    # Separate fictive and F(k,t) entries
    fictive_handles = []
    fictive_labels = []
    fkt_handles = []
    fkt_labels = []
    
    for handle, label in zip(handles, labels):
        if 'fictive' in label:
            fictive_handles.append(handle)
            fictive_labels.append(label.replace(' (fictive)', ''))
        else:
            fkt_handles.append(handle)
            fkt_labels.append(label.replace(' (F(k,t))', ''))
    
    # Create two-column legend
    legend1 = ax2.legend(fictive_handles, fictive_labels, loc='upper left', 
                        title=latex_format(r'$\xi_{\mathrm{fictive}} / \ln(10)$ (lines)', 'ξ_fictive / ln(10) (lines)'),
                        fontsize=8, title_fontsize=9)
    ax2.add_artist(legend1)
    
    legend2 = ax2.legend(fkt_handles, fkt_labels, loc='center left', 
                        title=latex_format(r'$\xi_{\mathrm{F(k,t)}}$ (dashed)', 'ξ_F(k,t) (dashed)'),
                        fontsize=8, title_fontsize=9)
    
    ax2.grid(True, alpha=0.3)
    
    # Add summary text box
    summary_text = latex_format(r'Theoretical: $\xi_{\mathrm{F(k,t)}} = \frac{\xi_{\mathrm{fictive}}}{\ln(10)}$', 
                               'Theoretical: ξ_F(k,t) = ξ_fictive / ln(10)')
    ax2.text(0.02, 0.98, summary_text, transform=ax2.transAxes, fontsize=11,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))
    
    plt.tight_layout()
    
    # Save plot
    output_file = Path(output_dir) / 'material_time_overlay_ln10_scaling.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSaved: {output_file}")
    
    # Also save as PDF
    output_file_pdf = Path(output_dir) / 'material_time_overlay_ln10_scaling.pdf'
    plt.savefig(output_file_pdf, bbox_inches='tight')
    print(f"Saved: {output_file_pdf}")
    
    plt.show()
    
    # Reset matplotlib settings
    plt.rcParams.update(plt.rcParamsDefault)
    
    # Print summary
    print(f"\n" + "=" * 70)
    print("OVERLAY ANALYSIS SUMMARY")
    print("=" * 70)
    print(f"Universal scaling relationship: ξ_F(k,t) = ξ_fictive / ln(10)")
    print(f"Scale factor: 1/ln(10) = {ln10_factor:.6f}")
    print(f"")
    print(f"Visual interpretation:")
    print(f"  - Solid lines: Rescaled fictive material time")
    print(f"  - Dashed lines: F(k,t) reference material time")
    print(f"  - Close overlay: Universal scaling works!")
    print(f"  - Panel 1: Detailed ε=0 comparison showing all three curves")
    print(f"  - Panel 2: All coupling strengths overlaid")
    
    return all_data

if __name__ == "__main__":
    create_overlay_plot()
