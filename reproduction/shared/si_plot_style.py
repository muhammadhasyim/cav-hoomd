"""
SI Plot Style Configuration for Material Time Reconstruction

Defines matplotlib styling for publication-quality figures in single-column format
with Computer Modern fonts matching LaTeX document typography.

Author: Material Time Analysis
Date: November 2, 2025
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
from contextlib import contextmanager

# Single column width for PRX/PRL style journals
SINGLE_COL_WIDTH = 3.3  # inches
DOUBLE_COL_WIDTH = 7.0  # inches (for reference)

# Font sizes following APS guidelines
SMALL_SIZE = 9
MEDIUM_SIZE = 10  
LARGE_SIZE = 11
HUGE_SIZE = 12

def configure_si_style():
    """
    Configure matplotlib with Computer Modern fonts and appropriate sizing
    for single-column supplementary information figures.
    """
    # Font configuration - Computer Modern to match LaTeX
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['CMU Serif', 'Computer Modern Roman', 'DejaVu Serif']
    plt.rcParams['font.sans-serif'] = ['CMU Sans Serif', 'Computer Modern Sans Serif', 'DejaVu Sans']
    plt.rcParams['mathtext.fontset'] = 'cm'
    plt.rcParams['mathtext.rm'] = 'serif'
    
    # Font sizes
    plt.rcParams['font.size'] = MEDIUM_SIZE
    plt.rcParams['axes.labelsize'] = MEDIUM_SIZE
    plt.rcParams['axes.titlesize'] = LARGE_SIZE
    plt.rcParams['xtick.labelsize'] = SMALL_SIZE
    plt.rcParams['ytick.labelsize'] = SMALL_SIZE
    plt.rcParams['legend.fontsize'] = SMALL_SIZE
    plt.rcParams['legend.title_fontsize'] = MEDIUM_SIZE
    plt.rcParams['figure.titlesize'] = LARGE_SIZE
    
    # Figure size - default to single column
    plt.rcParams['figure.figsize'] = (SINGLE_COL_WIDTH, SINGLE_COL_WIDTH * 0.75)
    plt.rcParams['figure.dpi'] = 100
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams['savefig.pad_inches'] = 0.05
    
    # Line widths and marker sizes appropriate for print
    plt.rcParams['lines.linewidth'] = 1.5
    plt.rcParams['lines.markersize'] = 4
    plt.rcParams['patch.linewidth'] = 1.0
    plt.rcParams['axes.linewidth'] = 0.8
    plt.rcParams['grid.linewidth'] = 0.5
    plt.rcParams['xtick.major.width'] = 0.8
    plt.rcParams['ytick.major.width'] = 0.8
    plt.rcParams['xtick.minor.width'] = 0.6
    plt.rcParams['ytick.minor.width'] = 0.6
    
    # Grid and layout
    plt.rcParams['axes.grid'] = False  # Turn off by default, enable per-plot
    plt.rcParams['grid.alpha'] = 0.3
    plt.rcParams['axes.axisbelow'] = True
    
    # Legend
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['legend.framealpha'] = 0.9
    plt.rcParams['legend.fancybox'] = False
    plt.rcParams['legend.edgecolor'] = 'black'
    
    # Use LaTeX for text rendering (optional, can be slow)
    # plt.rcParams['text.usetex'] = True
    # plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'


@contextmanager
def si_figure(ncols=1, nrows=1, width=None, height=None, **kwargs):
    """
    Context manager for creating SI figures with proper styling.
    
    Parameters:
    -----------
    ncols : int
        Number of subplot columns
    nrows : int
        Number of subplot rows
    width : float or None
        Figure width in inches (default: SINGLE_COL_WIDTH)
    height : float or None
        Figure height in inches (default: auto-calculated)
    **kwargs : dict
        Additional arguments passed to plt.subplots()
        
    Yields:
    -------
    fig, axes : matplotlib figure and axes
    
    Example:
    --------
    >>> with si_figure(ncols=2, nrows=1) as (fig, axes):
    ...     axes[0].plot(x, y)
    ...     axes[1].scatter(x, z)
    ...     fig.savefig('figure.pdf')
    """
    configure_si_style()
    
    # Set figure dimensions
    if width is None:
        width = SINGLE_COL_WIDTH
    if height is None:
        # Auto-calculate height maintaining reasonable aspect ratio
        height = width * 0.75 * nrows / ncols
    
    # Create figure
    fig, axes = plt.subplots(nrows, ncols, figsize=(width, height), **kwargs)
    
    try:
        yield fig, axes
    finally:
        plt.close(fig)


def save_si_figure(fig, basename, output_dir='.'):
    """
    Save figure in both PDF and PNG formats with SI naming convention.
    
    Parameters:
    -----------
    fig : matplotlib.figure.Figure
        Figure to save
    basename : str
        Base filename (e.g., 'fig_s1_lcurve')
    output_dir : str
        Output directory path
    """
    import os
    
    pdf_path = os.path.join(output_dir, f'{basename}.pdf')
    png_path = os.path.join(output_dir, f'{basename}.png')
    
    fig.savefig(pdf_path, format='pdf', bbox_inches='tight', pad_inches=0.05)
    fig.savefig(png_path, format='png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    
    print(f"Saved: {pdf_path}")
    print(f"Saved: {png_path}")


# Pre-configure style on import
configure_si_style()

