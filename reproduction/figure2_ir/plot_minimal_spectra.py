import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from tqdm import tqdm

# Set classic style
plt.style.use('classic')
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
})

def main():
    # Color map for different directories (0-10 only)
    n_dirs = 11  # Directories 0 through 10
    colors = plt.cm.coolwarm(np.linspace(0, 1, n_dirs))
    
    # Process each directory
    for dir_num in tqdm(range(n_dirs)):
        # Calculate coupling strength (without sqrt(N) factor)
        coupling = dir_num * 0.0001
        print(f"Processing directory {dir_num}: coupling strength ε = {coupling:.4f}")
        
        # Create new figure for each spectrum
        fig = plt.figure(figsize=(2.5, 1.5))
        ax = fig.add_subplot(111)
        
        try:
            # Load averaged spectrum data
            data = np.loadtxt(f'average_spectrum_dir_{dir_num}.txt')
            data = data[np.isfinite(data[:,1])]
            freq = data[:,0]
            mean_spec = data[:,1]
            
            # Filter frequency range
            mask = (freq >= 1300) & (freq <= 2100)
            freq = freq[mask]
            mean_spec = mean_spec[mask]
            
            # Normalize spectrum
            mean_spec = mean_spec / np.trapz(mean_spec, freq)
            
            # Plot spectrum
            ax.plot(freq, mean_spec, 
                   color=colors[dir_num], 
                   alpha=0.8,
                   linewidth=2)  # Doubled line thickness
            
            # Add vertical line at 1560 cm-1
            ax.axvline(1560, color='gray', linestyle='--', alpha=0.6, linewidth=1.5, label='1560 cm$^{-1}$')

        except Exception as e:
            print(f"Error processing directory {dir_num}: {e}")
            continue
        
        # Set plot limits
        ax.set_xlim([1300, 2100])
        ax.set_ylim([-0.0005, 0.03])
        
        # Add axis labels with larger font sizes
        ax.set_xlabel(r'$\omega$ $(\mathrm{cm^{-1}})$', fontsize=14)
        ax.set_ylabel(r'$n(\omega)\alpha(\omega)$', fontsize=14)
        
        # Keep minimal ticks but make them visible with larger font
        ax.tick_params(axis='both', which='major', labelsize=12)
        
        # Add more ticks for better readability
        ax.locator_params(axis='x', nbins=6)  # More x-axis ticks
        ax.locator_params(axis='y', nbins=5)  # More y-axis ticks
        
        # Add legend for the vertical line
        ax.legend(loc='upper right', fontsize=10, framealpha=0.8)
        
        # Remove spines initially (will be added back later)
        for spine in ax.spines.values():
            spine.set_visible(False)
        
        # Set background color to white
        ax.set_facecolor('white')
        fig.patch.set_facecolor('white')
        # Add black box around figure
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color('black')
            spine.set_linewidth(1.0)
        # Tight layout and save
        plt.tight_layout(pad=0)
        plt.savefig(f'minimal_spectrum_{dir_num:02d}.pdf', dpi=300, bbox_inches='tight')
        plt.close()  # Close the figure to free memory

if __name__ == "__main__":
    main() 