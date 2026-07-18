import numpy as np
import matplotlib.pyplot as plt
import os
from glob import glob
import argparse
from tqdm import tqdm
from scipy.signal import find_peaks, savgol_filter

def load_and_process_directory(dir_path):
    """Load all spectrum files from a directory and compute average and std"""
    spectrum_files = glob(os.path.join(dir_path, "spectrum_hoomd_*.txt"))
    print(len(spectrum_files))
    if not spectrum_files:
        return None, None, None
    
    # Load first file to get frequency axis
    first_data = np.loadtxt(spectrum_files[0])
    frequencies = first_data[:, 0]
    
    # Initialize array to store all spectra
    all_spectra = np.zeros((len(spectrum_files), len(frequencies)))
    
    # Load all spectra
    for i, file in enumerate(spectrum_files[50:]):
        try:
            data = np.loadtxt(file)
            all_spectra[i] = data[:, 1]  # Second column is intensity
        except:
            print(f"Error loading {file}")
            continue
    
    # Compute mean and standard deviation
    mean_spectrum = np.mean(all_spectra, axis=0)
    std_spectrum = np.std(all_spectra, axis=0)
    
    return frequencies, mean_spectrum, std_spectrum

def detect_peaks_constrained(frequencies, spectrum, dir_num=0, prev_second_peak=None):
    """Detect peaks in the spectrum with constraint that second peak decreases with directory number."""
    # Normalize spectrum for peak detection
    norm_spectrum = spectrum / np.max(spectrum)
    
    # Apply light smoothing for non-zero directories to help reveal buried peaks
    if dir_num != 0 and len(norm_spectrum) > 50:
        # Use Savitzky-Golay filter with small window to preserve peak shape
        window_length = min(21, len(norm_spectrum)//10)
        if window_length % 2 == 0:  # Must be odd
            window_length += 1
        if window_length >= 5:  # Minimum window size
            norm_spectrum = savgol_filter(norm_spectrum, window_length, 3)
    
    # Adjust parameters based on directory
    if dir_num == 0:
        # Directory 0 may have different peak structure
        height_threshold = 0.05
        distance_threshold = 40
        prominence_threshold = 0.03
    else:
        # Other directories: expect 4 peaks max, more sensitive to buried peaks
        height_threshold = 0.02  # Even lower threshold for buried peaks
        distance_threshold = 20   # Closer peaks allowed
        prominence_threshold = 0.015  # Lower prominence for subtle peaks
    
    # Find initial peaks with specified criteria
    peaks, properties = find_peaks(norm_spectrum, 
                                 height=height_threshold,
                                 distance=distance_threshold,
                                 prominence=prominence_threshold)
    
    # Get peak frequencies and intensities
    peak_frequencies = frequencies[peaks]
    peak_intensities = norm_spectrum[peaks]
    
    # Apply constraint for second peak if we have previous information
    if dir_num > 0 and prev_second_peak is not None and len(peak_frequencies) >= 2:
        # Sort peaks by frequency
        sorted_indices = np.argsort(peak_frequencies)
        sorted_freqs = peak_frequencies[sorted_indices]
        sorted_intensities = peak_intensities[sorted_indices]
        
        # Find potential second peaks that satisfy BOTH constraints:
        # 1. Frequency ≤ previous second peak (with small tolerance)
        # 2. Proximity constraint: within reasonable distance from previous peak
        max_frequency_jump = 50  # Maximum allowed jump in frequency
        potential_second_peaks = []
        
        for i, freq in enumerate(sorted_freqs[1:], 1):  # Skip first peak
            freq_constraint = freq <= prev_second_peak + 5  # Small tolerance for numerical precision
            proximity_constraint = abs(freq - prev_second_peak) <= max_frequency_jump
            
            if freq_constraint and proximity_constraint:
                # Score based on both intensity and proximity to previous peak
                intensity_score = sorted_intensities[i]
                proximity_score = 1.0 / (1.0 + abs(freq - prev_second_peak))  # Closer = higher score
                combined_score = intensity_score * 0.6 + proximity_score * 0.4  # Weight intensity more
                potential_second_peaks.append((i, freq, sorted_intensities[i], combined_score))
        
        if potential_second_peaks:
            # Choose the best candidate based on combined score (intensity + proximity)
            best_idx, best_freq, best_intensity, best_score = max(potential_second_peaks, key=lambda x: x[3])
            
            # Reconstruct peak arrays with constrained second peak
            final_peaks = [sorted_indices[0]]  # Keep first peak
            final_peaks.append(sorted_indices[best_idx])  # Add constrained second peak
            
            # Add remaining peaks (3rd, 4th) if they exist
            remaining_indices = [i for j, i in enumerate(sorted_indices[1:], 1) 
                               if j != best_idx and j < 4]  # Limit to 4 total peaks
            final_peaks.extend(remaining_indices[:2])  # Add up to 2 more peaks
            
            # Sort final peaks by frequency
            final_peaks = sorted(final_peaks, key=lambda x: frequencies[peaks[x]])
            peak_frequencies = frequencies[peaks[final_peaks]]
            peak_intensities = norm_spectrum[peaks[final_peaks]]
        else:
            # If no valid constrained peak found, use more aggressive search
            print(f"    ⚠️  No valid second peak found within constraints, using fallback search")
            # Look for peaks with relaxed proximity constraint but still decreasing
            fallback_candidates = []
            for i, freq in enumerate(sorted_freqs[1:], 1):
                if freq <= prev_second_peak + 5:  # Still must be decreasing
                    fallback_candidates.append((i, freq, sorted_intensities[i]))
            
            if fallback_candidates:
                # Among decreasing candidates, pick the one closest to previous peak
                best_idx, best_freq, best_intensity = min(fallback_candidates, 
                                                         key=lambda x: abs(x[1] - prev_second_peak))
                final_peaks = [sorted_indices[0], sorted_indices[best_idx]]
                # Add remaining peaks
                remaining_indices = [i for j, i in enumerate(sorted_indices[1:], 1) 
                                   if j != best_idx and j < 4]
                final_peaks.extend(remaining_indices[:2])
                final_peaks = sorted(final_peaks, key=lambda x: frequencies[peaks[x]])
                peak_frequencies = frequencies[peaks[final_peaks]]
                peak_intensities = norm_spectrum[peaks[final_peaks]]
        
    # Limit to 4 peaks maximum for non-zero directories
    if dir_num != 0 and len(peak_frequencies) > 4:
        # Sort by intensity and keep top 4
        top_indices = np.argsort(peak_intensities)[-4:]
        sorted_by_freq = sorted(top_indices, key=lambda x: peak_frequencies[x])
        peak_frequencies = peak_frequencies[sorted_by_freq]
        peak_intensities = peak_intensities[sorted_by_freq]
    
    return peak_frequencies, peak_intensities, None

def print_peak_analysis(dir_num, peak_frequencies, peak_intensities):
    """Print detailed peak analysis for a directory."""
    epsilon = dir_num * 0.0001
    expected_peaks = "variable" if dir_num == 0 else "≤4"
    print(f"\n--- Directory {dir_num} (ε = {epsilon:.4f}) [Expected peaks: {expected_peaks}] ---")
    
    if len(peak_frequencies) == 0:
        print("No significant peaks detected.")
        return
    
    print(f"Detected {len(peak_frequencies)} peaks (with constraint):")
    
    # Sort peaks by frequency for consistent ordering
    sorted_indices = np.argsort(peak_frequencies)
    sorted_freqs = peak_frequencies[sorted_indices]
    sorted_intensities = peak_intensities[sorted_indices]
    
    # Reference peak positions - expand for better analysis
    ref_peaks = {'A-A': 1560, 'B-B': 2335}
    
    for i, (freq, intensity) in enumerate(zip(sorted_freqs, sorted_intensities)):
        peak_label = f"  Peak {i+1}"
        if i == 1 and dir_num > 0:  # Highlight second peak for tracking
            peak_label += " [TRACKED]"
        print(f"{peak_label}: {freq:.1f} cm⁻¹ (intensity: {intensity:.3f})")
        
        # Check proximity to reference peaks
        for ref_name, ref_freq in ref_peaks.items():
            diff = abs(freq - ref_freq)
            if diff < 100:  # Within 100 cm⁻¹ of reference
                print(f"    -> Close to {ref_name} reference ({ref_freq} cm⁻¹), diff: {diff:.1f} cm⁻¹")
    
    # Find the most intense peak
    max_idx = np.argmax(sorted_intensities)
    print(f"  Strongest peak: {sorted_freqs[max_idx]:.1f} cm⁻¹ (intensity: {sorted_intensities[max_idx]:.3f})")
    
    # Warn if we detect fewer than expected peaks for non-zero directories
    if dir_num != 0 and len(peak_frequencies) < 2:
        print(f"  ⚠️  Warning: Only {len(peak_frequencies)} peak(s) detected, some peaks may be buried")
    elif dir_num != 0 and len(peak_frequencies) == 4:
        print(f"  ✓ Found maximum expected peaks (4)")
    
    # Show second peak tracking info
    if len(sorted_freqs) >= 2 and dir_num > 0:
        print(f"  📍 Second peak for next constraint: {sorted_freqs[1]:.1f} cm⁻¹")

def main():
    # Set up classic style with LaTeX fonts
    plt.style.use('classic')
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "figure.figsize": [8, 8],  # Even larger figure size for poster
        "font.size": 24,  # Much larger base font size
        "axes.linewidth": 2.5,  # Thicker axes
        "lines.linewidth": 3.5,  # Thicker lines
        "xtick.labelsize": 20,
        "ytick.labelsize": 20
    })
    
    # Create figure with space for colorbar
    fig = plt.figure()
    gs = fig.add_gridspec(1, 2, width_ratios=[20, 1], wspace=0.05)
    ax = fig.add_subplot(gs[0])
    cax = fig.add_subplot(gs[1])
    
    # Color map for different directories
    colors = plt.cm.coolwarm(np.linspace(0, 1, 21))
    
    # Define offset between spectra
    offset = 0.004  # Adjust this value to change spacing between spectra
    
    # Create a list to store lines for colorbar mapping
    lines = []
    
    # Storage for peak analysis summary
    all_peaks_summary = []
    prev_second_peak = None  # Track second peak from previous directory
    
    print("=== CONSTRAINED PEAK DETECTION ANALYSIS ===")
    print("Note: Second peak frequency should decrease or stay constant across directories")
    print("Using ε (epsilon) = directory_number × 0.0001")
    
    # Process each directory
    for dir_num in tqdm(range(16)):  # 0 to 20
        dir_path = str(dir_num)
        coupling = dir_num * 0.0001  # Calculate actual coupling strength
        
        if not os.path.exists(dir_path):
            continue
            
        freq, mean_spec, std_spec = load_and_process_directory(dir_path)
        if freq is None:
            continue
        
        # Detect and print peaks before normalization using constrained search
        peak_freqs, peak_intensities, _ = detect_peaks_constrained(freq, mean_spec, dir_num, prev_second_peak)
        
        # Update tracking of second peak for next iteration
        if len(peak_freqs) >= 2:
            current_second_peak = sorted(peak_freqs)[1]  # Second peak by frequency
            if prev_second_peak is not None:
                constraint_satisfied = current_second_peak <= prev_second_peak + 5
                print(f"Second peak constraint: {current_second_peak:.1f} ≤ {prev_second_peak:.1f} → {'✓' if constraint_satisfied else '✗'}")
            prev_second_peak = current_second_peak
        
        print_peak_analysis(dir_num, peak_freqs, peak_intensities)
        
        # Store peak information for summary
        epsilon = dir_num * 0.0001
        all_peaks_summary.append({
            'dir': dir_num,
            'epsilon': epsilon,
            'peaks': list(zip(peak_freqs, peak_intensities))
        })
        
        # Save average and std for each directory
        output_data = np.column_stack((freq, mean_spec, std_spec))
        np.savetxt(f'average_spectrum_dir_{dir_num}.txt', 
                   output_data,
                   header='Frequency(cm^-1) Average_Intensity Std_Intensity',
                   delimiter='\t')
            
        # Normalize spectrum to [0,1] range
        mean_spec = mean_spec / np.trapz(mean_spec, freq)
        std_spec = std_spec / np.max(mean_spec)
        
        # Add vertical offset based on directory number
        vertical_offset = dir_num * offset
        
        # Plot mean spectrum with shaded error bars
        line = ax.plot(freq, mean_spec + vertical_offset, 
                      color=colors[dir_num], 
                      alpha=0.8)[0]
        # ax.fill_between(freq/8, 
        #                mean_spec - std_spec + vertical_offset, 
        #                mean_spec + std_spec + vertical_offset,
        #                color=colors[dir_num], 
        #                alpha=0.2)
        lines.append(line)
    
    # Print summary analysis
    print("\n=== CONSTRAINED PEAK ANALYSIS SUMMARY ===")
    print("Peak evolution across epsilon values (with second peak constraint):")
    
    # Track second peak evolution specifically
    print(f"\nSecond Peak Evolution (constraint: decreasing frequency):")
    second_peaks = []
    for summary in all_peaks_summary:
        dir_num, epsilon = summary['dir'], summary['epsilon']
        if len(summary['peaks']) >= 2:
            # Sort peaks by frequency and get second one
            sorted_peaks = sorted(summary['peaks'], key=lambda x: x[0])
            second_peak = sorted_peaks[1]
            second_peaks.append((dir_num, epsilon, second_peak[0], second_peak[1]))
            
    for i, (dir_num, epsilon, freq, intensity) in enumerate(second_peaks):
        trend = ""
        if i > 0:
            prev_freq = second_peaks[i-1][2]
            if freq <= prev_freq:
                trend = "✓ (decreasing/constant)"
            else:
                trend = f"✗ (increased by {freq-prev_freq:.1f})"
        print(f"  Dir {dir_num:2d} (ε={epsilon:.4f}): {freq:.1f} cm⁻¹ (intensity: {intensity:.3f}) {trend}")
    
    # Analyze peak trends for reference regions
    ref_peaks = {'A-A': 1560, 'B-B': 2335}
    for ref_name, ref_freq in ref_peaks.items():
        print(f"\n{ref_name} Peak Region ({ref_freq} cm⁻¹ ± 100):")
        for summary in all_peaks_summary:
            dir_num, epsilon = summary['dir'], summary['epsilon']
            close_peaks = [peak for peak in summary['peaks'] 
                         if abs(peak[0] - ref_freq) < 100]
            if close_peaks:
                closest_peak = min(close_peaks, key=lambda x: abs(x[0] - ref_freq))
                freq_shift = closest_peak[0] - ref_freq
                print(f"  Dir {dir_num:2d} (ε={epsilon:.4f}): {closest_peak[0]:.1f} cm⁻¹ "
                      f"(shift: {freq_shift:+.1f}, intensity: {closest_peak[1]:.3f})")
            else:
                print(f"  Dir {dir_num:2d} (ε={epsilon:.4f}): No peak detected")
    
    print("\n" + "="*60)
    
    # Save peak data to text file
    print("\nSaving peak data to 'peak_analysis.txt'...")
    with open('peak_analysis.txt', 'w') as f:
        f.write("# Peak Analysis Results\n")
        f.write("# Directory, Epsilon, Peak2_Freq, Peak2_Intensity, Peak3_Freq, Peak3_Intensity, Peak4_Freq, Peak4_Intensity\n")
        f.write("# Frequencies in cm⁻¹, intensities are normalized\n")
        f.write("# -1 indicates peak not detected\n")
        f.write("Dir\tEpsilon\tPeak2_Freq\tPeak2_Int\tPeak3_Freq\tPeak3_Int\tPeak4_Freq\tPeak4_Int\n")
        
        for summary in all_peaks_summary:
            dir_num = summary['dir']
            epsilon = summary['epsilon']
            peaks = summary['peaks']
            
            # Sort peaks by frequency
            sorted_peaks = sorted(peaks, key=lambda x: x[0])
            
            # Extract 2nd, 3rd, 4th peaks (pad with -1 if missing)
            peak_data = [-1, -1, -1, -1, -1, -1]  # [freq2, int2, freq3, int3, freq4, int4]
            
            for i in range(1, min(4, len(sorted_peaks))):  # Peaks 2, 3, 4 (indices 1, 2, 3)
                peak_data[(i-1)*2] = sorted_peaks[i][0]      # Frequency
                peak_data[(i-1)*2 + 1] = sorted_peaks[i][1]  # Intensity
            
            f.write(f"{dir_num:d}\t{epsilon:.6f}\t{peak_data[0]:.1f}\t{peak_data[1]:.6f}\t"
                   f"{peak_data[2]:.1f}\t{peak_data[3]:.6f}\t{peak_data[4]:.1f}\t{peak_data[5]:.6f}\n")
    
    print("Peak data saved successfully!")
    print(f"Average spectra saved as 'average_spectrum_dir_{{0-{len(all_peaks_summary)-1}}}.txt'")
    
    ax.axvline(1560, color='gray', linestyle='--', alpha=0.5,label='A-A, 1560 cm$^{-1}$')
    ax.axvline(2335, color='gray', linestyle='--', alpha=0.5,label='B-B, 2335 cm$^{-1}$')
    ax.legend(loc='upper left', fontsize=20)

    # Customize plot
    ax.set_xlabel(r'Frequency (cm$^{-1}$)', fontsize=24)
    ax.set_ylabel(r'Intensity (arb. units)', fontsize=24)
    #ax.set_ylim([0, 0.035])
    ax.set_xlim([1200, 3000])
    #ax.set_title(r'Stacked IR Spectra by Coupling Strength', fontsize=14, pad=20)
    
    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=plt.cm.coolwarm, 
                              norm=plt.Normalize(0, 19 * 0.0001))
    cbar = plt.colorbar(sm, cax=cax)
    cbar.set_label(r'$\epsilon$', fontsize=24)
    cbar.ax.tick_params(labelsize=20)
    
    # Add grid with light gray color
    ax.grid(True, color='gray', linestyle='--', alpha=0.3)
    
    # Set background color to white
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')
    
    # Remove y-axis ticks since they're not meaningful with the offset
    ax.set_yticks([])
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout(pad=2.0)  # Add more padding
    plt.subplots_adjust(right=0.85)  # Make room for colorbar label
    plt.savefig('allspectra_res_poster.pdf', dpi=300, bbox_inches='tight', pad_inches=0.5) 
    plt.show()
            
if __name__ == "__main__":
    main()
