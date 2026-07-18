import numpy as np
import matplotlib.pyplot as plt
import argparse

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Plot a single spectrum file with LaTeX styling.')
    parser.add_argument('filename', help='Path to the spectrum data file')
    parser.add_argument('--xlim', nargs=2, type=float, default=[1200, 3000],
                       help='x-axis limits [min max] (default: [1200 3000])')
    parser.add_argument('--normalize', action='store_true',
                       help='Normalize spectrum using trapz')
    
    args = parser.parse_args()

    # Set up classic style with LaTeX fonts
    plt.style.use('classic')
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "figure.figsize": [5, 3],
        "font.size": 12,
        "axes.linewidth": 1.5,
        "lines.linewidth": 2
    })

    # Load data
    data = np.loadtxt(args.filename)
    freq = data[:, 0]
    spectrum = data[:, 1]
    if len(data.shape) > 1 and data.shape[1] > 2:
        std_spectrum = data[:, 2]
    else:
        std_spectrum = None

    # Normalize if requested
    if args.normalize:
        spectrum = spectrum / np.trapz(spectrum, freq)
        if std_spectrum is not None:
            std_spectrum = std_spectrum / np.max(spectrum)

    # Create figure
    fig, ax = plt.subplots()

    # Plot spectrum
    ax.plot(freq/8, spectrum, 'k-', linewidth=2, label='Spectrum')
    
    # Add reference lines
    ax.axvline(1560, color='gray', linestyle='--', alpha=0.5, label='A-A, 1560 cm$^{-1}$')
    ax.axvline(2335, color='gray', linestyle='--', alpha=0.5, label='B-B, 2335 cm$^{-1}$')

    # Customize plot
    ax.set_xlabel(r'Frequency (cm$^{-1}$)', fontsize=12)
    ax.set_ylabel(r'Intensity (arb. units)', fontsize=12)
    ax.set_xlim(args.xlim)
    ax.set_title('IR Spectrum', fontsize=14, pad=20)
    
    # Add grid with light gray color
    ax.grid(True, color='gray', linestyle='--', alpha=0.3)
    
    # Add legend
    ax.legend(loc='upper right', fontsize=12)
    
    # Set background color to white
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save and show
    output_filename = 'spectrum.pdf'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_filename}")
    plt.show()

if __name__ == "__main__":
    main() 
