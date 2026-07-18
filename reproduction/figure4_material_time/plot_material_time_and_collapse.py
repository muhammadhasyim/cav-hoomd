#!/usr/bin/env python3
"""
Three-Panel Material Time and Data Collapse Analysis Script

Creates a three-panel horizontal plot:
- Left Panel: Fictive material time ξ_S(t) evolution from fictive temperature analysis
- Middle Panel: Material time ξ(t) evolution for different coupling strengths  
- Right Panel: Collapsed F(k,t) data Φ(ξ) with stretched exponential fits

Uses consistent coolwarm colormap for coupling strength visualization.

Author: Cavity MD Analysis Suite
Date: 2025
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for remote servers
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
import matplotlib.style
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
import re
import sys
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

REPRO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPRO_ROOT / "shared"))
sys.path.insert(0, str(REPRO_ROOT / "figure4_material_time"))

# Set matplotlib to use classic style
mpl.style.use('classic')

# Import the MaterialTimeAnalyzer class
from material_time_analysis import MaterialTimeAnalyzer
from latex_config_adobe import setup_latex_fonts, latex_safe

# Import fictive material time calculation functions
from overlay_shifted_material_times import calculate_time_shifted_fictive_material_time

# Import corrected material time reconstruction
from material_time_correct import MaterialTimeReconstructor
from run_corrected_analysis import CorrectedMaterialTimeAnalyzer

class MaterialTimeCollapseVisualizer:
    """
    Visualizes material time evolution and collapsed F(k,t) data
    in a two-panel format with consistent styling.
    """
    
    def __init__(self, base_dir='.', use_latex=True, profile=None):
        """
        Initialize the visualizer.
        
        Parameters:
        -----------
        base_dir : str
            Base directory containing data files
        use_latex : bool
            Whether to use LaTeX rendering for text
        profile : DatasetProfile, optional
            Dataset profile for in-repo aging data
        """
        self.base_dir = Path(base_dir)
        self.use_latex = use_latex
        self.profile = profile
        
        # Set up matplotlib styling
        self._setup_matplotlib()
        
        if profile is not None:
            self.coupling_mapping = {
                entry.axis_value: (idx, f"{entry.axis_value:g}")
                for idx, entry in enumerate(profile.couplings)
            }
            self.axis = profile.axis
            self._tag_map = profile.tag_map()
        else:
            # Define coupling strength mapping
            self.coupling_mapping = {
                0.0: (0, '0.0'),
                3e-4: (1, '3 \\times 10^{-4}'),
                5e-4: (2, '5 \\times 10^{-4}'),
                7e-4: (3, '7 \\times 10^{-4}'),
                1e-3: (4, '1 \\times 10^{-3}')
            }
            self.axis = "epsilon"
            self._tag_map = None
        
        # Set up consistent colormap
        self.coupling_values = list(self.coupling_mapping.keys())
        self.norm = Normalize(vmin=0, vmax=max(self.coupling_values))
        self.cmap = plt.colormaps.get_cmap('coolwarm')
        
        # Initialize corrected material time analyzer
        self.corrected_analyzer = CorrectedMaterialTimeAnalyzer(
            criterion_value=0.1,
            alpha_smooth=0.01,  # Small regularization for smoothness
            n_grid_points=500,  # High resolution for smooth curves
            verbose=False  # Reduce verbosity for plotting
        )
        
    def _setup_matplotlib(self):
        """Set up matplotlib style and LaTeX rendering to match fictive temperature plots."""
        # Use classic style and LaTeX setup
        import matplotlib
        matplotlib.use('Agg')  # Use non-interactive backend
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        mpl.style.use('classic')
        
        if self.use_latex:
            try:
                # Use the same LaTeX setup as fictive temperature plots
                self.use_latex = setup_latex_fonts()
                print("Using LaTeX rendering with Adobe Illustrator compatibility")
                print("LaTeX successfully enabled with Computer Modern fonts")
            except Exception as e:
                print(f"LaTeX not available, using default rendering: {e}")
                self.use_latex = False
                self._setup_default_fonts()
        else:
            self._setup_default_fonts()
    
    def _setup_default_fonts(self):
        """Set up default font rendering."""
        plt.rcParams['text.usetex'] = False
        plt.rcParams['font.family'] = 'DejaVu Sans'
        plt.rcParams['mathtext.fontset'] = 'dejavusans'
        
    def _format_text(self, text, fallback=None):
        """Format text with LaTeX or fallback using latex_safe function."""
        return latex_safe(text, fallback if fallback else text, self.use_latex)
    
    def get_coupling_color(self, coupling_value):
        """Get color for a coupling value."""
        return self.cmap(self.norm(coupling_value))
    
    def get_coupling_color_viridis(self, coupling_value):
        """Get viridis color for a coupling value (for collapse panel)."""
        viridis_cmap = plt.colormaps.get_cmap('viridis')
        return viridis_cmap(self.norm(coupling_value))
    
    def convert_epsilon_to_lambda(self, epsilon):
        """Convert coupling constant ε to λ = ε/ωc."""
        OMEGA_C_CM_INV = 1560  # cm^-1
        OMEGA_C_AU = OMEGA_C_CM_INV * 4.556335e-6  # Convert to atomic units
        return epsilon / OMEGA_C_AU
    
    def load_material_time_data(self):
        """
        Compute material time data for all coupling strengths using CORRECTED method.
        
        Returns:
        --------
        dict
            Dictionary with coupling values as keys and (time, xi) tuples as values
        """
        print("Computing material time data using CORRECTED reconstruction method...")
        
        # Define coupling directories
        coupling_dirs = {
            0.0: self.base_dir / "cavity_coupling_0epos00_switch_200.0ps",
            3e-4: self.base_dir / "cavity_coupling_3eneg04_switch_200.0ps", 
            5e-4: self.base_dir / "cavity_coupling_5eneg04_switch_200.0ps",
            7e-4: self.base_dir / "cavity_coupling_7eneg04_switch_200.0ps",
            1e-3: self.base_dir / "cavity_coupling_1eneg03_switch_200.0ps"
        }
        
        material_time_data = {}
        
        for coupling_value, coupling_dir in coupling_dirs.items():
            if not coupling_dir.exists():
                print(f"  Warning: {coupling_dir} does not exist, skipping coupling ε = {coupling_value}")
                continue
                
            try:
                print(f"  Processing coupling ε = {coupling_value}...")
                
                # Load and normalize F(k,t) data using corrected analyzer
                fkt_files = self.corrected_analyzer.load_fkt_data(str(coupling_dir))
                waiting_times = self.corrected_analyzer.calculate_waiting_times(fkt_files)
                normalized_fkt = self.corrected_analyzer.normalize_fkt_data(fkt_files)
                
                # Calculate material time using CORRECTED method
                time_grid, xi_grid = self.corrected_analyzer.calculate_material_time(normalized_fkt)
                
                # Store results
                material_time_data[coupling_value] = (time_grid, xi_grid)
                
                print(f"    ✓ Computed {len(time_grid)} points, ξ range: [{xi_grid.min():.3f}, {xi_grid.max():.3f}]")
                
            except Exception as e:
                print(f"    ✗ Error processing coupling ε = {coupling_value}: {e}")
                # Create fallback linear material time for this coupling
                time_fallback = np.linspace(0, 2800, 1000)
                xi_fallback = time_fallback / np.max(time_fallback) * 2.0
                material_time_data[coupling_value] = (time_fallback, xi_fallback)
                print(f"    Using linear fallback for coupling ε = {coupling_value}")
                
        return material_time_data
    
    def load_collapsed_data(self):
        """
        Compute collapsed F(k,t) data using the CORRECTED material time method.
        
        Returns:
        --------
        dict
            Dictionary with coupling values as keys and collapsed data as values
        """
        print("Computing collapsed F(k,t) data using CORRECTED material time method...")
        
        # Define coupling directories
        coupling_dirs = {
            0.0: self.base_dir / "cavity_coupling_0epos00_switch_200.0ps",
            3e-4: self.base_dir / "cavity_coupling_3eneg04_switch_200.0ps", 
            5e-4: self.base_dir / "cavity_coupling_5eneg04_switch_200.0ps",
            7e-4: self.base_dir / "cavity_coupling_7eneg04_switch_200.0ps",
            1e-3: self.base_dir / "cavity_coupling_1eneg03_switch_200.0ps"
        }
        
        collapsed_data = {}
        
        for coupling_value, coupling_dir in coupling_dirs.items():
            if not coupling_dir.exists():
                print(f"  Warning: {coupling_dir} does not exist, skipping coupling ε = {coupling_value}")
                continue
                
            print(f"  Processing coupling ε = {coupling_value}...")
            
            try:
                # Load and normalize F(k,t) data using corrected analyzer
                fkt_files = self.corrected_analyzer.load_fkt_data(str(coupling_dir))
                waiting_times = self.corrected_analyzer.calculate_waiting_times(fkt_files)
                normalized_fkt = self.corrected_analyzer.normalize_fkt_data(fkt_files)
                
                # Calculate material time using CORRECTED method
                time_grid, xi_grid = self.corrected_analyzer.calculate_material_time(normalized_fkt)
                
                # Collapse data using corrected material time
                analyzer_collapsed_data = self.corrected_analyzer.collapse_data(
                    (time_grid, xi_grid), 
                    normalized_fkt, 
                    waiting_times
                )
                
                # Convert to the format expected by our plotting function
                coupling_collapsed_data = []
                for ref_num, data in analyzer_collapsed_data.items():
                    xi_coords = data['xi']
                    normalized_fkt_data = data['R_fkt']
                    coupling_collapsed_data.append((xi_coords, normalized_fkt_data, ref_num))
                
                collapsed_data[coupling_value] = coupling_collapsed_data
                print(f"    ✓ Computed {len(coupling_collapsed_data)} collapsed reference curves")
                
            except Exception as e:
                print(f"    ✗ Error processing coupling ε = {coupling_value}: {e}")
                continue
                
        return collapsed_data
    
    def _parse_coupling_from_dirname(self, dirname):
        """Parse coupling value from directory name."""
        if "0epos00" in dirname:
            return 0.0
        elif "3eneg04" in dirname:
            return 3e-4
        elif "5eneg04" in dirname:
            return 5e-4
        elif "7eneg04" in dirname:
            return 7e-4
        elif "1eneg03" in dirname:
            return 1e-3
        else:
            return None
    
    
    def load_stretched_exponential_fits(self):
        """
        Load stretched exponential fit parameters.
        
        Returns:
        --------
        dict
            Dictionary with coupling values as keys and (beta, r_squared) tuples as values
        """
        print("Loading stretched exponential fit parameters...")
        
        fit_file = self.base_dir / 'stretched_exponential_fit_results.txt'
        
        if not fit_file.exists():
            print("  Stretched exponential fit file not found")
            return {}
        
        try:
            # Read fit results - skip comments and header properly
            with open(fit_file, 'r') as f:
                lines = f.readlines()
            
            # Find data start (after header lines)
            data_start = 0
            for i, line in enumerate(lines):
                if line.strip().startswith('#') or line.strip() == '' or 'Coupling_au' in line:
                    data_start = i + 1
                else:
                    break
            
            # Load only the numeric columns (skip the Label column)
            data = []
            for line in lines[data_start:]:
                if line.strip():
                    parts = line.strip().split()
                    if len(parts) >= 7:
                        try:
                            # Extract only numeric columns: coupling_au, A, A_err, beta, beta_err, R_squared
                            coupling_au = float(parts[0])
                            beta = float(parts[4])
                            r_squared = float(parts[6])
                            data.append([coupling_au, beta, r_squared])
                        except ValueError:
                            continue
            
            data = np.array(data)
            
            fit_params = {}
            for row in data:
                coupling_au = row[0]
                beta = row[1]
                r_squared = row[2]
                
                fit_params[coupling_au] = (beta, r_squared)
                
            print(f"  Loaded fit parameters for {len(fit_params)} coupling strengths")
            return fit_params
            
        except Exception as e:
            print(f"  Error loading fit parameters: {e}")
            return {}
    
    def load_fictive_material_time_data(self):
        """
        Load material time data from fictive temperature analysis using proper methodology.
        
        Returns:
        --------
        material_time_data : dict
            Dictionary with coupling values as keys and (time, xi, aging_rate, time_lab, xi_scaled) tuples as values
        """
        print("Loading fictive material time data using proper methodology...")
        
        material_time_data = {}
        ln10_factor = 1.0 / np.log(10)
        
        for coupling_value in self.coupling_values:
            if self._tag_map is not None:
                coupling_str = self._tag_map[coupling_value]
            elif coupling_value == 0.0:
                coupling_str = '0epos00'
            else:
                exp_str = f'{coupling_value:.0e}'
                if 'e-04' in exp_str:
                    if '3e-04' in exp_str:
                        coupling_str = '3eneg04'
                    elif '5e-04' in exp_str:
                        coupling_str = '5eneg04'
                    elif '7e-04' in exp_str:
                        coupling_str = '7eneg04'
                    else:
                        coupling_str = exp_str.replace('e-04', 'eneg04')
                elif 'e-03' in exp_str:
                    if '1e-03' in exp_str:
                        coupling_str = '1eneg03'
                    else:
                        coupling_str = exp_str.replace('e-03', 'eneg03')
                else:
                    coupling_str = exp_str.replace('e', 'epos').replace('epos-', 'eneg')
            
            filename = f'coupling_{coupling_str}_fictive_time_series.txt'
            filepath = self.base_dir / 'time_series_output' / filename
            
            if not filepath.exists():
                print(f"  Warning: File {filename} not found, skipping coupling {coupling_value}")
                continue
            
            try:
                # Load the fictive time series data with proper header handling
                # Skip 12 rows to get past all the header comments and column names
                data = np.loadtxt(filepath, skiprows=12, usecols=(0, 2))  # Only load time and tau_fictive columns
                time_ps = data[:, 0]  # Time in ps
                fictive_relaxation_time = data[:, 1]  # τ_fictive in ps
                
                # Calculate time-shifted fictive material time
                time_shifted, xi_fictive_scaled = calculate_time_shifted_fictive_material_time(
                    time_ps, fictive_relaxation_time, t_shift=200.0
                )
                
                # Apply 1/ln(10) scaling to match proper methodology
                xi_fictive_scaled *= ln10_factor
                
                # Calculate aging rate γ(t) = dξ/dt
                aging_rate = np.gradient(xi_fictive_scaled, time_shifted)
                
                # Store all data for this coupling
                material_time_data[coupling_value] = (
                    time_shifted,      # time
                    xi_fictive_scaled, # xi (material time)
                    aging_rate,        # aging rate
                    time_ps,          # original time (time_lab)
                    xi_fictive_scaled  # xi_scaled (same as xi for fictive)
                )
                
                print(f"  Loaded ε = {coupling_value}: {len(time_shifted)} points, max ξ = {np.max(xi_fictive_scaled):.6f}")
                
            except Exception as e:
                print(f"  Error loading data for coupling {coupling_value}: {e}")
                continue
        
        return material_time_data
    
    def create_three_panel_plot(self, output_file='material_time_and_collapse.pdf'):
        """
        Create the main three-panel horizontal plot.
        
        Parameters:
        -----------
        output_file : str
            Output filename
        """
        print("\nCreating three-panel horizontal plot...")
        
        # Load all data
        fictive_material_time_data = self.load_fictive_material_time_data()
        material_time_data = self.load_material_time_data()
        collapsed_data = self.load_collapsed_data()
        fit_params = self.load_stretched_exponential_fits()
        
        # Create figure with three horizontal panels - even smaller size for larger fonts
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 3), sharex=False)
        
        # Left panel: Material time evolution
        self._plot_material_time_panel(ax1, material_time_data)
        
        # Middle panel: Collapsed data and fits
        self._plot_collapse_panel(ax2, collapsed_data, fit_params)
        
        # Right panel: Fictive material time evolution
        self._plot_fictive_material_time_panel(ax3, fictive_material_time_data)
        
        # Ensure both material time panels have the same y-limits
        y1_min, y1_max = ax1.get_ylim()
        y3_min, y3_max = ax3.get_ylim()
        
        # Use the same y-limits for both panels
        common_y_max = max(y1_max, y3_max)
        ax1.set_ylim(0, common_y_max)
        ax3.set_ylim(0, common_y_max)
        
        # Adjust layout and save
        plt.tight_layout()
        
        output_path = self.base_dir / output_file
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.savefig(output_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
        
        print(f"Saved: {output_path}")
        print(f"Saved: {output_path.with_suffix('.png')}")
        
        plt.show()
        
    def _plot_fictive_material_time_panel(self, ax, material_time_data):
        """Plot the right panel: Fictive material time evolution."""
        print("  Creating right panel: Fictive material time evolution...")
        
        for coupling_value in sorted(self.coupling_values):
            if coupling_value not in material_time_data:
                continue
            
            # Unpack the extended data structure
            time, xi, aging_rate, time_lab, xi_scaled = material_time_data[coupling_value]
            color = self.get_coupling_color(coupling_value)
            label = self.coupling_mapping[coupling_value][1]
            
            # Format label properly for LaTeX - use lambda instead of epsilon
            if self.use_latex:
                if coupling_value == 0.0:
                    lambda_val = self.convert_epsilon_to_lambda(coupling_value)
                    label_formatted = r'$\lambda = 0$'
                else:
                    # Convert to lambda and use simpler scientific notation
                    lambda_val = self.convert_epsilon_to_lambda(coupling_value)
                    label_formatted = f'$\\lambda = {lambda_val:.3f}$'
            else:
                lambda_val = self.convert_epsilon_to_lambda(coupling_value)
                label_formatted = f'λ = {lambda_val:.3f}'
            
            # Plot with dashed lines
            ax.plot(time, xi, color=color, linewidth=2, alpha=0.8, linestyle='--', label=label_formatted)
        
        # Formatting to match the material time panel (same axis limits)
        ax.set_xlabel(self._format_text(r'$t$ (ps)', 't (ps)'), fontsize=14)
        ax.set_ylabel(latex_safe(r'$\xi_{\mathrm{TN}}(t)$', 'ξ_TN(t)', self.use_latex), fontsize=16)
        # Remove title
        
        ax.legend(fontsize=10, loc='best')
        
        # Use same limits as material time panel (get y-limits from left panel)
        ax.set_xlim(0, 2600)
        ax.set_ylim(0, None)  # Will be set to match left panel after plotting
        
        # Add grid lines every 800 ps (reduced tick marks to prevent overlap)
        import matplotlib.ticker as ticker
        ax.xaxis.set_major_locator(ticker.MultipleLocator(800))
        ax.grid(True, alpha=0.3)
        
    def _plot_material_time_panel(self, ax, material_time_data):
        """Plot the left panel with material time evolution."""
        print("  Creating left panel: Material time evolution...")
        
        for coupling_value in sorted(self.coupling_values):
            if coupling_value not in material_time_data:
                continue
                
            time, xi = material_time_data[coupling_value]
            color = self.get_coupling_color(coupling_value)
            label = self.coupling_mapping[coupling_value][1]
            
            # Format label properly for LaTeX - use lambda instead of epsilon
            if self.use_latex:
                if coupling_value == 0.0:
                    label_formatted = r'$\lambda = 0$'
                else:
                    lambda_val = self.convert_epsilon_to_lambda(coupling_value)
                    label_formatted = f'$\\lambda = {lambda_val:.3f}$'
            else:
                lambda_val = self.convert_epsilon_to_lambda(coupling_value)
                label_formatted = f'λ = {lambda_val:.3f}'
            
            # Plot with solid lines only (no markers)
            ax.plot(time, xi, color=color, linewidth=3, alpha=0.8, label=label_formatted)
        
        # Formatting - use proper LaTeX labels, remove title
        ax.set_xlabel(self._format_text(r'$t$ (ps)', 't (ps)'), fontsize=14)
        ax.set_ylabel(self._format_text(r'$\xi(t)$', 'ξ(t)'), fontsize=14)
        # Remove title
        
        ax.legend(fontsize=10, loc='best')
        
        # Set reasonable limits with max time at 2600 ps
        ax.set_xlim(0, 2600)
        ax.set_ylim(0, None)
        
        # Add grid lines every 800 ps (reduced tick marks to prevent overlap)
        import matplotlib.ticker as ticker
        ax.xaxis.set_major_locator(ticker.MultipleLocator(800))
        ax.grid(True, alpha=0.3)
        
    def _plot_collapse_panel(self, ax, collapsed_data, fit_params):
        """Plot the middle panel: Data collapse and fits (overlapped style)"""
        print("  Creating middle panel: Data collapse and fits (overlapped style)...")
        
        # Reset fit plotting flag
        if hasattr(self, '_fit_plotted'):
            delattr(self, '_fit_plotted')
        
        # Determine max xi range
        all_xi_values = []
        for coupling_value in self.coupling_values:
            if coupling_value in collapsed_data:
                for xi_coords, normalized_fkt, ref_num in collapsed_data[coupling_value]:
                    all_xi_values.extend(xi_coords)
        
        # Set max xi to 3.0 as requested
        max_xi = 3.0
        
        legend_elements = []
        
        # Plot data for each coupling strength
        for coupling_value in sorted(self.coupling_values):
            if coupling_value not in collapsed_data:
                continue
                
            color = self.get_coupling_color(coupling_value)  # Use coolwarm colormap (back to original)
            label = self.coupling_mapping[coupling_value][1]
            coupling_curves = collapsed_data[coupling_value]
            
            # Plot raw data points (semi-transparent) - exactly like the original
            for xi_coords, normalized_fkt, ref_num in coupling_curves:
                # Filter to plotting range and remove invalid points
                mask = (xi_coords <= max_xi) & (xi_coords > 0) & (normalized_fkt > 0) & (normalized_fkt <= 1)
                if np.any(mask):
                    ax.plot(xi_coords[mask], normalized_fkt[mask], color=color, alpha=0.2, linewidth=1)
            
            # Plot only one stretched exponential fit curve (first coupling with fit data)
            if coupling_value in fit_params and not hasattr(self, '_fit_plotted'):
                beta, r_squared = fit_params[coupling_value]
                
                # Create fit line with ln(10) factor
                xi_fit = np.linspace(0.01, max_xi, 200)
                # Correct form: R(ξ) = exp(-ln(10) * ξ^β)
                R_fit = np.exp(-np.log(10) * xi_fit**beta)
                
                line = ax.plot(xi_fit, R_fit, color='black', linewidth=2, alpha=1.0, zorder=10)[0]  # Thicker line
                
                # Mark that we've plotted the fit to avoid duplicates
                self._fit_plotted = True
        
        # Add legend with stretched exponential formula only
        if fit_params:  # Only add legend if we have fits
            # Add a legend entry for the stretched exponential formula
            if self.use_latex:
                formula_label = r'$\Phi_k(\xi) = e^{-\xi^\beta}$'
            else:
                formula_label = 'Φ(ξ) = e^(-ξ^β)'
            
            # Create a dummy line for the formula legend entry
            formula_line = ax.plot([], [], color='black', linewidth=2, alpha=1.0)[0]
            
            # Only show the formula in legend (no β value)
            ax.legend([formula_line], [formula_label], fontsize=14, loc='upper right')
        
        # Format the plot with proper LaTeX labels, remove title
        ax.set_xlabel(self._format_text(r'$\xi$', 'ξ'), fontsize=14)
        ax.set_ylabel(self._format_text(r'$\Phi_k(\xi)$', 'Φ_k(ξ)'), fontsize=14)
        # Remove title
        
        # Add grid lines every 1 unit of material time
        import matplotlib.ticker as ticker
        ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, max_xi)
        ax.set_ylim(0, 1.05)

def main():
    """Main function."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Create three-panel horizontal material time and collapse plot')
    parser.add_argument('--profile', default='paper', help='Dataset profile (default: paper)')
    parser.add_argument('--staged-root', type=Path, default=None, help='Override staged data root')
    parser.add_argument('--base_dir', default=None,
                       help='Base directory containing data files')
    parser.add_argument('--output', default=None,
                       help='Output filename')
    parser.add_argument('--no-latex', action='store_true',
                       help='Disable LaTeX rendering for text')
    
    args = parser.parse_args()

    from repro_bootstrap import setup_profile

    profile = setup_profile(args, default="paper")
    if args.base_dir is not None:
        base_dir = args.base_dir
    elif profile.name != "paper":
        base_dir = str(profile.staged_root)
    else:
        base_dir = '/media/extradrive/Trajectories/final_nodiss_cavitymd'

    if args.output is None:
        if profile.figures_dir is not None:
            args.output = str(profile.figures_dir / "material_time_and_collapse.pdf")
        else:
            args.output = "material_time_and_collapse.pdf"
    
    # Create visualizer with LaTeX enabled by default
    visualizer = MaterialTimeCollapseVisualizer(
        base_dir=base_dir, 
        use_latex=not args.no_latex,
        profile=profile if profile.name != "paper" else None,
    )
    
    # Create the three-panel plot
    visualizer.create_three_panel_plot(output_file=args.output)
    
    print("\nThree-panel horizontal plot created successfully!")

if __name__ == "__main__":
    main()
