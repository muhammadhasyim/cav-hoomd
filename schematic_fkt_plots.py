#!/usr/bin/env python3
"""
Generate schematic F(k,t) vs t plots with exponential decay.
Creates three small square plots with different relaxation times.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple

# Configure matplotlib for LaTeX rendering with Computer Modern font
plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'serif',
    'font.serif': ['Computer Modern'],
    'font.size': 12,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14
})

def exponential_decay(t: np.ndarray, tau: float, amplitude: float = 1.0) -> np.ndarray:
    """
    Generate exponential decay function.
    
    Parameters:
    -----------
    t : np.ndarray
        Time array
    tau : float
        Relaxation time constant
    amplitude : float
        Initial amplitude (default: 1.0)
        
    Returns:
    --------
    np.ndarray
        F(k,t) values
    """
    return amplitude * np.exp(-t / tau)

def create_schematic_plot(ax: plt.Axes, tau: float, t_max: float = 1.0) -> None:
    """
    Create a single schematic F(k,t) vs t plot.
    
    Parameters:
    -----------
    ax : plt.Axes
        Matplotlib axes object
    tau : float
        Relaxation time constant
    t_max : float
        Maximum time value for plotting
    """
    # Generate time array
    t = np.linspace(0, t_max, 200)
    
    # Calculate F(k,t) using exponential decay
    fkt = exponential_decay(t, tau)
    
    # Plot the curve
    ax.plot(t, fkt, 'k-', linewidth=2)
    
    # Set labels with LaTeX formatting
    ax.set_xlabel(r'$t$', fontsize=12)
    ax.set_ylabel(r'$F(k,t)$', fontsize=12)
    
    # Set square aspect ratio with equal ranges
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 1.0)
    ax.set_aspect('equal', adjustable='box')
    
    # Remove tick labels for schematic appearance
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Add grid for better visibility
    ax.grid(True, alpha=0.3)
    
    # Set frame
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

def create_three_schematic_plots() -> Tuple[plt.Figure, np.ndarray]:
    """
    Create three schematic plots with different relaxation times.
    
    Returns:
    --------
    Tuple[plt.Figure, np.ndarray]
        Figure and array of axes
    """
    # Create figure with three subplots
    fig, axes = plt.subplots(1, 3, figsize=(9, 3))
    
    # Define three different relaxation times (progressively shorter)
    relaxation_times = [0.4, 0.25, 0.1]  # tau values
    
    # Create each plot (without titles)
    for i, (ax, tau) in enumerate(zip(axes, relaxation_times)):
        create_schematic_plot(ax, tau)
    
    # Adjust layout
    plt.tight_layout()
    
    return fig, axes

def main():
    """Main function to generate and save the schematic plots."""
    print("Generating schematic F(k,t) vs t plots...")
    
    # Create the plots
    fig, axes = create_three_schematic_plots()
    
    # Save the figure as PDF
    output_filename = 'schematic_fkt_decay.pdf'
    fig.savefig(output_filename, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    
    print(f"Saved schematic plots to: {output_filename}")
    
    # Also save individual plots for inset use
    for tau in [0.4, 0.25, 0.1]:
        # Create individual figure

        fig_single, ax_single = plt.subplots(figsize=(3, 3))
        create_schematic_plot(ax_single, tau)
        
        # Save individual plot as PDF
        individual_filename = f'schematic_fkt_tau_{tau:.1f}.pdf'
        fig_single.savefig(individual_filename, bbox_inches='tight',
                          facecolor='white', edgecolor='none')
        
        plt.close(fig_single)
        print(f"Saved individual plot: {individual_filename}")
    
    # Show the combined plot
   # plt.show()

if __name__ == "__main__":
    main()
