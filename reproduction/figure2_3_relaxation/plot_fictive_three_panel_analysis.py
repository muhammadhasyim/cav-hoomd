#!/usr/bin/env python3
"""
Fictive Temperature Analysis Figure Generator

Creates comprehensive fictive temperature analysis figures with two modes:

DEFAULT MODE (Separate Figures):
1. Standalone Material Time Figure: Material time evolution from fictive temperature analysis
2. Two-Panel Relaxation Figure: 
   - Panel 1: Normalization time vs coupling strengths (different waiting times overlaid)
   - Panel 2: Relaxation time prediction vs waiting times (different coupling strengths overlaid)

LEGACY MODE (--three-panel flag):
3-Panel Combined Figure: All three panels in a single figure (deprecated but maintained for compatibility)

This script combines data from fictive temperature analysis to show the complete
picture of how cavity coupling affects material time, normalization, and relaxation dynamics.

Author: Cavity MD Analysis Suite
Date: 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
from matplotlib.colors import Normalize
from matplotlib import cm
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
import os
import sys
from pathlib import Path

REPRO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPRO_ROOT / "shared"))

# Import functions from the overlay script for proper fictive calculations
try:
    from overlay_shifted_material_times import (
        calculate_time_shifted_fictive_material_time,
        generate_fkt_from_fictive,
        find_relaxation_time_01_criterion
    )
    from compare_material_time_t200 import calculate_fictive_material_time_t200
except ImportError as e:
    print(f"Warning: Could not import fictive calculation functions: {e}")
    print("Will use simplified fallback methods")

from direct_material_time import (
    direct_material_time_from_masters,
    staged_coupling_dir,
    tn_material_time_from_series,
)

# Import our LaTeX configuration
try:
    from latex_config_adobe import setup_latex_fonts, latex_safe
except ImportError:
    print("Warning: latex_config_adobe not found, using matplotlib defaults")
    def latex_safe(latex_text, fallback_text, use_latex_flag):
        return latex_text if use_latex_flag else fallback_text
    def setup_latex_fonts():
        pass

# Set up matplotlib style
mpl.style.use('classic')

class FictiveThreePanelAnalyzer:
    """
    Analyzer for creating three-panel fictive temperature analysis figures.
    """
    
    def __init__(self, base_dir='.', use_latex=True, profile=None):
        """
        Initialize the analyzer.
        
        Parameters:
        -----------
        base_dir : str
            Base directory containing data files
        use_latex : bool
            Whether to use LaTeX rendering
        profile : DatasetProfile, optional
            Dataset profile (overrides hardcoded paper coupling catalogue)
        """
        self.base_dir = Path(base_dir)
        self.use_latex = use_latex
        self.profile = profile

        if profile is not None:
            self.coupling_values = [entry.axis_value for entry in profile.couplings]
            self.axis = profile.axis
            self._tag_map = profile.tag_map()
        else:
            self.coupling_values = [0.0, 3e-4, 5e-4, 7e-4, 1e-3]
            self.axis = "epsilon"
            self._tag_map = {
                0.0: "0epos00",
                3e-4: "3eneg04",
                5e-4: "5eneg04",
                7e-4: "7eneg04",
                1e-3: "1eneg03",
            }

        if self.axis == "lambda":
            self.coupling_values_scaled = self.coupling_values
            self.norm = Normalize(vmin=0, vmax=max(self.coupling_values))
        else:
            self.coupling_values_scaled = [val * 1e4 for val in self.coupling_values]
            self.norm = Normalize(vmin=0, vmax=max(self.coupling_values))
        self.cmap = plt.colormaps.get_cmap('coolwarm')
        
        # Coupling strength mapping for labels
        self.coupling_mapping = {
            0.0: '0',
            3e-4: '3 \\times 10^{-4}',
            5e-4: '5 \\times 10^{-4}', 
            7e-4: '7 \\times 10^{-4}',
            1e-3: '1 \\times 10^{-3}'
        }
        
        # Constants for ε to λ conversion (matching plot_fkt_analysis.py)
        self.OMEGA_C_CM_INV = 1560  # cm^-1
        self.OMEGA_C_AU = self.OMEGA_C_CM_INV * 4.556335e-6  # Convert to atomic units

    def _coupling_str(self, coupling_value: float) -> str:
        return self._tag_map[coupling_value]

    def _axis_value_label(self, value: float) -> str:
        if self.axis == "lambda":
            return self.format_lambda_label(value, self.use_latex)
        return self.format_coupling_label(value)
    
    def convert_epsilon_to_lambda(self, epsilon):
        """Convert coupling constant ε to λ = ε/ωc.

        For profiles whose control axis is already λ (e.g. aging_weak_lambda),
        the stored coupling values ARE λ, so no rescaling is applied.
        """
        if self.axis == "lambda":
            return epsilon
        return epsilon / self.OMEGA_C_AU
    
    def format_lambda_label(self, lambda_val, use_latex=True):
        """Format λ label for legends."""
        if use_latex:
            return f'$\\lambda = {lambda_val:.3f}$'
        else:
            return f'λ = {lambda_val:.3f}'

    def format_lambda_axis_label(self, use_latex=True):
        """Format λ axis label."""
        if use_latex:
            return r'$\lambda$ (a.u.)'
        return 'λ (a.u.)'
    
    def format_coupling_axis_label(self, use_latex=True):
        if self.axis == "lambda":
            return self.format_lambda_axis_label(use_latex)
        return latex_safe(r'$\varepsilon$ ($\times 10^{-4}$ a.u.)', 'ε (×10⁻⁴ a.u.)', use_latex)

    def coupling_axis_xlim(self):
        if self.axis == "lambda":
            return (0.0, max(self.coupling_values) * 1.05)
        return (0, 10)
    
    def get_coupling_color(self, coupling_value):
        """Get color for a coupling value."""
        return self.cmap(self.norm(coupling_value))
    
    def format_coupling_label(self, coupling_value):
        """Format coupling value for legends."""
        if self.axis == "lambda":
            return self.format_lambda_label(coupling_value, self.use_latex)
        if coupling_value == 0:
            return latex_safe('$\\varepsilon = 0$', 'ε = 0', self.use_latex)
        label = self.coupling_mapping.get(coupling_value, f'{coupling_value:.1e}')
        return latex_safe(f'$\\varepsilon = {label}$', f'ε = {label}', self.use_latex)

    def _coupling_plot_storage(self, coupling_value: float) -> float:
        """Store coupling on the plot axis (λ directly; ε in ×10⁻⁴ units)."""
        if self.axis == "lambda":
            return coupling_value
        return coupling_value * 1e4

    def _coupling_from_storage(self, stored_value: float) -> float:
        """Recover the profile coupling value from stored plot data."""
        if self.axis == "lambda":
            return stored_value
        return stored_value / 1e4

    def _coupling_plot_x(self, coupling_value: float) -> float:
        """X-axis coordinate for coupling-strength panels."""
        return self.convert_epsilon_to_lambda(coupling_value)

    def load_measured_material_time_data(self):
        """
        Load F(k,t)-reconstructed material time from staged layout.

        Time in the file is already referenced to cavity turn-on (t = 0 at switch).

        Returns
        -------
        dict
            coupling_value -> (time_ps, xi) arrays
        """
        mt_path = self.base_dir / "material_time_all_couplings.txt"
        if not mt_path.exists():
            print(f"  Warning: measured material time file not found: {mt_path}")
            return {}

        header_cols = None
        rows: list[list[float]] = []
        with open(mt_path, encoding="utf-8") as fh:
            for line in fh:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                tokens = stripped.split()
                try:
                    float(tokens[0])
                except ValueError:
                    header_cols = tokens
                    continue
                rows.append([float(x) for x in tokens])

        if header_cols is None or not rows:
            print(f"  Warning: could not parse {mt_path}")
            return {}

        col_index = {name: idx for idx, name in enumerate(header_cols)}
        data = np.asarray(rows, dtype=float)
        measured: dict[float, tuple[np.ndarray, np.ndarray]] = {}

        for coupling_value in self.coupling_values:
            tag = f"{coupling_value:g}".replace(".", "p")
            t_col_name = f"t_{tag}_ps"
            xi_col_name = f"xi_{tag}"
            if t_col_name not in col_index or xi_col_name not in col_index:
                print(f"  Warning: columns for {self.axis}={coupling_value} missing in {mt_path.name}")
                continue
            measured[coupling_value] = (
                data[:, col_index[t_col_name]],
                data[:, col_index[xi_col_name]],
            )
            print(
                f"  Loaded measured ξ for {self.axis}={coupling_value}: "
                f"max={measured[coupling_value][1].max():.3f}"
            )

        return measured

    def load_direct_measured_material_time(self):
        """
        Load measured material time as direct integral dt/tau_0.1 from master F(k,t).

        Returns
        -------
        dict
            coupling_value -> (time_since_turn_on_ps, xi)
        """
        print("Loading direct measured material time from master F(k,t)...")
        switch_ps = self.profile.switch_time_ps if self.profile is not None else 200.0
        measured: dict[float, tuple[np.ndarray, np.ndarray]] = {}

        if self.profile is not None:
            coupling_dirs = [
                (entry.axis_value, self.base_dir / entry.staged_dir_name)
                for entry in self.profile.couplings
            ]
        else:
            coupling_dirs = [
                (
                    coupling_value,
                    staged_coupling_dir(self.base_dir, self._coupling_str(coupling_value), switch_ps),
                )
                for coupling_value in self.coupling_values
            ]

        for coupling_value, coupling_dir in coupling_dirs:
            time_grid, xi_grid = direct_material_time_from_masters(
                coupling_dir,
                switch_time_ps=switch_ps,
            )
            if time_grid.size == 0:
                print(f"  Warning: no direct measured material time for {self.axis}={coupling_value}")
                continue
            measured[coupling_value] = (time_grid, xi_grid)
            print(
                f"  Loaded direct measured xi for {self.axis}={coupling_value}: "
                f"max={xi_grid.max():.3f}, late slope={np.polyfit(time_grid[time_grid >= time_grid.max()-400], xi_grid[time_grid >= time_grid.max()-400], 1)[0]:.2e}"
            )

        return measured

    def load_tn_material_time(self):
        """
        Load TN material time from fictive temperature series (no 1/ln(10) scaling).

        Returns
        -------
        dict
            coupling_value -> (time_since_turn_on_ps, xi, aging_rate)
        """
        print("Loading TN material time from fictive temperature series...")
        switch_ps = self.profile.switch_time_ps if self.profile is not None else 200.0
        tn_data: dict[float, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

        for coupling_value in self.coupling_values:
            coupling_str = self._coupling_str(coupling_value)
            time_series_file = (
                self.base_dir / "time_series_output" / f"coupling_{coupling_str}_fictive_time_series.txt"
            )
            if not time_series_file.exists():
                print(f"  Warning: TN file not found: {time_series_file}")
                continue

            try:
                raw = np.loadtxt(time_series_file, skiprows=12)
                time_ps = raw[:, 0]
                tau_fictive = raw[:, 2]
                time_shifted, xi_shifted, aging_rate = tn_material_time_from_series(
                    time_ps,
                    tau_fictive,
                    switch_time_ps=switch_ps,
                )
                tn_data[coupling_value] = (time_shifted, xi_shifted, aging_rate)
                print(
                    f"  Loaded TN xi for {self.axis}={coupling_value}: "
                    f"max={xi_shifted.max():.3f}, late slope={np.polyfit(time_shifted[time_shifted >= time_shifted.max()-400], xi_shifted[time_shifted >= time_shifted.max()-400], 1)[0]:.2e}"
                )
            except Exception as exc:
                print(f"  Error loading TN material time for {self.axis}={coupling_value}: {exc}")

        return tn_data
    
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
            coupling_str = self._coupling_str(coupling_value)
            
            # Try to load time series data
            time_series_file = self.base_dir / 'time_series_output' / f'coupling_{coupling_str}_fictive_time_series.txt'
            
            if time_series_file.exists():
                try:
                    # Load data - format: time_ps, fictive_temperature_K, fictive_relaxation_time_ps, predicted_tau_alpha_ps
                    data = np.loadtxt(time_series_file, skiprows=12)  # Skip header lines
                    time = data[:, 0]  # Time in ps
                    tau_fictive = data[:, 2]  # Fictive relaxation time in ps
                    
                    # Calculate aging rate for all time points (1/τ_fictive exists everywhere)
                    aging_rate_full = 1.0 / tau_fictive
                    
                    # Use the proper time-shifted fictive material time calculation
                    try:
                        time_shifted, xi_shifted = calculate_time_shifted_fictive_material_time(
                            time, tau_fictive, t_shift=200.0)
                        
                        # Apply 1/ln(10) scaling for consistency with F(k,t) analysis
                        xi_scaled = xi_shifted * ln10_factor
                        
                        # Convert back to lab time for consistency with other panels
                        time_lab = time_shifted + 200.0
                        
                        # Pad with zeros before coupling turn-on for plotting
                        xi_full = np.zeros_like(time)
                        time_mask = time >= 200.0
                        if np.any(time_mask):
                            # Interpolate scaled material time back to full time grid
                            xi_interp = interp1d(time_lab, xi_scaled, kind='linear', 
                                               bounds_error=False, fill_value=0.0)
                            xi_full[time_mask] = xi_interp(time[time_mask])
                        
                        # Store both versions: for plotting (xi_full) and for F(k,t) generation (xi_scaled)
                        material_time_data[coupling_value] = (time, xi_full, aging_rate_full, time_lab, xi_scaled)
                        print(f"  Loaded {self.axis} = {coupling_value}: {len(time)} points, max ξ = {xi_scaled.max():.6f}")
                        
                    except Exception as e:
                        print(f"    Error in fictive calculation for ε = {coupling_value}: {e}")
                        print("    Falling back to simple integration...")
                        # Fallback to simple integration
                        dt = np.diff(time)
                        dt = np.append(dt, dt[-1])
                        time_mask = time >= 200.0
                        if np.any(time_mask):
                            integrand = 1.0 / tau_fictive[time_mask]
                            xi_integrated = np.cumsum(integrand * dt[time_mask])
                            xi_scaled = xi_integrated * ln10_factor
                            xi_full = np.zeros_like(time)
                            xi_full[time_mask] = xi_scaled
                            time_lab = time[time_mask]
                            material_time_data[coupling_value] = (time, xi_full, aging_rate_full, time_lab, xi_scaled)
                            print(f"  Fallback loaded ε = {coupling_value}: max ξ = {xi_scaled.max():.6f}")
                        
                except Exception as e:
                    print(f"    Error loading {time_series_file}: {e}")
                    continue
            else:
                print(f"  Warning: File not found: {time_series_file}")
        
        return material_time_data
    
    def calculate_fictive_relaxation_times_vs_coupling(self, material_time_data):
        """
        Calculate relaxation times from fictive material time for different coupling strengths and waiting times.
        
        Parameters:
        -----------
        material_time_data : dict
            Dictionary with fictive material time data
            
        Returns:
        --------
        normalization_data : dict
            Dictionary with relaxation time data vs coupling strength for different waiting times
        """
        print("Calculating relaxation times from fictive material time vs coupling strength...")
        
        normalization_data = {}
        beta = 0.55  # Stretched exponential exponent
        waiting_times = [0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]  # ps
        
        for tw in waiting_times:
            print(f"  Processing waiting time tw = {tw} ps...")
            
            coupling_strengths = []
            relaxation_times = []
            
            for coupling_value in self.coupling_values:
                if coupling_value not in material_time_data:
                    continue
                    
                try:
                    # Get material time data for this coupling
                    time, xi_full, aging_rate, time_lab, xi_scaled = material_time_data[coupling_value]
                    
                    # Generate F(k,t) from fictive material time
                    time_lag, fkt_fictive = generate_fkt_from_fictive(
                        time_lab, xi_scaled, tw=tw, beta=beta)
                    
                    # Find relaxation time where F(k,t) = 0.1
                    tau_relaxation = find_relaxation_time_01_criterion(time_lag, fkt_fictive)
                    
                    if not np.isnan(tau_relaxation) and tau_relaxation > 0:
                        coupling_strengths.append(self._coupling_plot_storage(coupling_value))
                        relaxation_times.append(tau_relaxation)
                    else:
                        print(f"    Invalid relaxation time for ε = {coupling_value} at tw = {tw}: τ = {tau_relaxation}")
                        
                except Exception as e:
                    print(f"    Error processing ε = {coupling_value} for tw = {tw}: {e}")
                    continue
            
            if len(relaxation_times) > 0:
                normalization_data[tw] = {
                    'coupling_strengths': np.array(coupling_strengths),
                    'normalization_times': np.array(relaxation_times)
                }
                print(f"    Calculated {len(relaxation_times)} relaxation times for tw = {tw} ps")
            else:
                print(f"    No valid relaxation times for tw = {tw} ps")
        
        print(f"  Calculated fictive relaxation data for {len(normalization_data)} waiting times")
        return normalization_data
    
    def _create_synthetic_normalization_data(self):
        """Create synthetic normalization data as fallback."""
        normalization_data = {}
        waiting_times = [0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]
        
        for tw in waiting_times:
            coupling_strengths = []
            norm_times = []
            
            for coupling_value in self.coupling_values:
                base_norm_time = 50.0 + tw * 0.1
                coupling_effect = coupling_value * 10000
                norm_time = base_norm_time * (1 + coupling_effect)
                
                coupling_strengths.append(self._coupling_plot_storage(coupling_value))
                norm_times.append(norm_time)
            
            normalization_data[tw] = {
                'coupling_strengths': np.array(coupling_strengths),
                'normalization_times': np.array(norm_times)
            }
        
        return normalization_data
    
    def calculate_fictive_relaxation_times_vs_waiting_time(self, material_time_data):
        """
        Calculate relaxation times from fictive material time vs waiting time for different coupling strengths.
        
        Parameters:
        -----------
        material_time_data : dict
            Dictionary with fictive material time data
            
        Returns:
        --------
        relaxation_data : dict
            Dictionary with relaxation time data vs waiting time for different coupling strengths
        """
        print("Calculating relaxation times from fictive material time vs waiting time...")
        
        relaxation_data = {}
        beta = 0.55  # Stretched exponential exponent
        waiting_times = [0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]  # ps
        
        for coupling_value in self.coupling_values:
            if coupling_value not in material_time_data:
                continue
                
            print(f"  Processing coupling ε = {coupling_value}...")
            
            waiting_times_valid = []
            relaxation_times_valid = []
            
            # Get material time data for this coupling
            time, xi_full, aging_rate, time_lab, xi_scaled = material_time_data[coupling_value]
            
            for tw in waiting_times:
                try:
                    # Generate F(k,t) from fictive material time
                    time_lag, fkt_fictive = generate_fkt_from_fictive(
                        time_lab, xi_scaled, tw=tw, beta=beta)
                    
                    # Find relaxation time where F(k,t) = 0.1
                    tau_relaxation = find_relaxation_time_01_criterion(time_lag, fkt_fictive)
                    
                    if not np.isnan(tau_relaxation):
                        waiting_times_valid.append(tw)
                        relaxation_times_valid.append(tau_relaxation)
                        
                except Exception as e:
                    print(f"    Error processing tw = {tw} for ε = {coupling_value}: {e}")
                    continue
            
            if len(relaxation_times_valid) > 0:
                relaxation_data[coupling_value] = {
                    'waiting_times': np.array(waiting_times_valid),
                    'relaxation_times': np.array(relaxation_times_valid)
                }
                print(f"    Calculated {len(relaxation_times_valid)} relaxation times for ε = {coupling_value}")
            else:
                print(f"    No valid relaxation times for ε = {coupling_value}")
        
        print(f"  Calculated fictive relaxation data for {len(relaxation_data)} coupling strengths")
        return relaxation_data
    
    def _create_synthetic_relaxation_data(self):
        """Create synthetic relaxation data as fallback."""
        relaxation_data = {}
        waiting_times = np.array([0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        
        for coupling_value in self.coupling_values:
            base_relaxation = 75.0
            aging_effect = waiting_times * 0.05
            coupling_effect = coupling_value * 5000
            
            relaxation_times = base_relaxation + aging_effect + coupling_effect * (1 + waiting_times / 1000)
            
            relaxation_data[coupling_value] = {
                'waiting_times': waiting_times,
                'relaxation_times': relaxation_times
            }
        
        return relaxation_data
    
    def create_material_time_figure(self, output_file='fictive_material_time_evolution.pdf'):
        """
        Create standalone material time evolution figure (Panel 1).
        
        Parameters:
        -----------
        output_file : str
            Output filename for the figure
        """
        print("\nCreating standalone material time evolution figure...")
        
        # Set matplotlib style and configure LaTeX fonts for Adobe Illustrator
        plt.style.use('classic')
        if self.use_latex:
            setup_latex_fonts()
            print("  Using LaTeX Computer Modern fonts for all text elements")
        
        # Load direct measured and TN material time data
        measured_material_time_data = self.load_direct_measured_material_time()
        tn_material_time_data = self.load_tn_material_time()
        
        # Create figure with single panel
        fig, ax = plt.subplots(figsize=(6, 3.5))
        
        # Material time: measured F(k,t) integral vs TN prediction
        self._plot_material_time_panel(ax, measured_material_time_data, tn_material_time_data)
        
        # Adjust layout manually instead of tight_layout to avoid warnings
        plt.subplots_adjust(top=0.90, bottom=0.15, left=0.12, right=0.95)
        
        output_path = self.base_dir / output_file
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.savefig(output_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
        
        print(f"Saved: {output_path}")
        print(f"Saved: {output_path.with_suffix('.png')}")
        
        plt.show()
        
        # Reset matplotlib settings
        plt.rcParams.update(plt.rcParamsDefault)
    
    def create_two_panel_relaxation_figure(self, output_file='fictive_two_panel_relaxation_analysis.pdf'):
        """
        Create two-panel relaxation analysis figure (Panels 2 & 3).
        Replicates exact styling from plot_fkt_analysis.py relaxation_analysis_filtered.pdf
        
        Parameters:
        -----------
        output_file : str
            Output filename for the figure
        """
        print("\nCreating two-panel relaxation analysis figure...")
        
        # Set matplotlib style and configure LaTeX fonts for Adobe Illustrator
        plt.style.use('classic')
        if self.use_latex:
            setup_latex_fonts()
            print("  Using LaTeX Computer Modern fonts for all text elements")
        
        # Load all data - now using fictive predictions for panels 2 and 3
        material_time_data = self.load_fictive_material_time_data()
        normalization_data = self.calculate_fictive_relaxation_times_vs_coupling(material_time_data)
        relaxation_data = self.calculate_fictive_relaxation_times_vs_waiting_time(material_time_data)
        
        # Create figure with exact same size and layout as relaxation_analysis_filtered.pdf
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        
        # Panel 1: Relaxation time vs coupling strength for different waiting times
        self._plot_normalization_panel_styled(ax1, normalization_data)
        
        # Panel 2: Relaxation time vs waiting times for different coupling strengths  
        self._plot_relaxation_panel_styled(ax2, relaxation_data)
        
        # Use exact same layout adjustment as original
        plt.tight_layout(pad=2.0)  # Add padding for colorbars above plots
        
        output_path = self.base_dir / output_file
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.savefig(output_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
        
        print(f"Saved: {output_path}")
        print(f"Saved: {output_path.with_suffix('.png')}")
        
        plt.show()
        
        # Reset matplotlib settings
        plt.rcParams.update(plt.rcParamsDefault)

    def create_three_panel_figure(self, output_file='fictive_three_panel_analysis.pdf'):
        """
        Create the three-panel figure combining all analyses.
        [DEPRECATED: Use create_material_time_figure() and create_two_panel_relaxation_figure() instead]
        
        Parameters:
        -----------
        output_file : str
            Output filename for the figure
        """
        print("\nCreating three-panel fictive temperature analysis figure...")
        
        # Set matplotlib style and configure LaTeX fonts for Adobe Illustrator
        plt.style.use('classic')
        if self.use_latex:
            setup_latex_fonts()
            print("  Using LaTeX Computer Modern fonts for all text elements")
        
        # Load all data - now using fictive predictions for panels 2 and 3
        material_time_data = self.load_fictive_material_time_data()
        measured_material_time_data = self.load_direct_measured_material_time()
        tn_material_time_data = self.load_tn_material_time()
        normalization_data = self.calculate_fictive_relaxation_times_vs_coupling(material_time_data)
        relaxation_data = self.calculate_fictive_relaxation_times_vs_waiting_time(material_time_data)
        
        # Create figure with three panels - wider figure with smaller fonts
        fig = plt.figure(figsize=(11.5, 3.5))
        gs = gridspec.GridSpec(1, 3, width_ratios=[1.4, 1.1, 1.1], wspace=0.3)
        
        # Main plot axes
        ax1 = fig.add_subplot(gs[0])  # Material time panel
        ax2 = fig.add_subplot(gs[1])  # Middle panel  
        ax3 = fig.add_subplot(gs[2])  # Right panel
        
        # Panel 1: Material time evolution (direct measured + TN)
        self._plot_material_time_panel(ax1, measured_material_time_data, tn_material_time_data)
        
        # Panel 2: Normalization time vs coupling strengths (different waiting times)
        self._plot_normalization_panel(ax2, normalization_data)
        
        # Panel 3: Relaxation time vs waiting times (different coupling strengths)
        self._plot_relaxation_panel(ax3, relaxation_data)
        
        # Adjust layout manually to accommodate colorbars
        plt.subplots_adjust(top=0.82, bottom=0.15, left=0.08, right=0.95, wspace=0.35)
        
        output_path = self.base_dir / output_file
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.savefig(output_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
        
        print(f"Saved: {output_path}")
        print(f"Saved: {output_path.with_suffix('.png')}")
        
        plt.show()
        
        # Reset matplotlib settings
        plt.rcParams.update(plt.rcParamsDefault)
    
    def _plot_normalization_panel_styled(self, ax, normalization_data):
        """Plot Panel 1: Relaxation time vs coupling strength (styled like relaxation_analysis_filtered.pdf)."""
        print("  Creating Panel 1: Relaxation time vs coupling strength...")
        
        # Extract waiting times from keys and filter out tw=0 as in the original
        all_waiting_times = sorted([tw for tw in normalization_data.keys() if tw > 0])
        
        if not all_waiting_times:
            print("    Warning: No valid waiting times found!")
            return
        
        # Create color mapping for waiting times (Panel 1) - use shifted times
        waiting_times_shifted = [(tw - 200) for tw in all_waiting_times]  # Shift by -200 ps
        waiting_time_norm = Normalize(vmin=0, vmax=max(waiting_times_shifted))
        waiting_time_cmap = plt.colormaps.get_cmap('viridis')
        
        # Store data for Panel 1 colorbar
        panel1_scatter_points = []
        
        # Find ε=0 coupling for normalization
        epsilon_zero_coupling = 0.0
        
        # Plot normalized data for Panel 1
        for tw in all_waiting_times:
            if tw not in normalization_data:
                continue
                
            data = normalization_data[tw]
            coupling_strengths = data['coupling_strengths']  # Already in 10^-4 units
            relaxation_times = data['normalization_times']
            
            coupling_values_orig = [
                self._coupling_from_storage(val) for val in coupling_strengths
            ]
            coupling_values_lambda = [
                self._coupling_plot_x(val) for val in coupling_values_orig
            ]
            
            # Find λ=0 / ε=0 reference value for normalization
            epsilon_zero_value = None
            for i, coupling_val in enumerate(coupling_values_orig):
                if abs(coupling_val - epsilon_zero_coupling) < 1e-8:
                    epsilon_zero_value = relaxation_times[i]
                    break
            
            if epsilon_zero_value is None or epsilon_zero_value <= 0:
                print(f"    Warning: No valid ε=0 reference for tw={tw} ps, skipping this waiting time")
                continue
            
            # Normalize all coupling values by ε=0 for this waiting time
            normalized_rel_times = relaxation_times / epsilon_zero_value
            
            if len(normalized_rel_times) > 0:  # Only plot if we have data
                tw_shifted = tw - 200  # Shift by -200 ps
                color = waiting_time_cmap(waiting_time_norm(tw_shifted))
                
                # Plot with both lines and markers in one call (matching original style)
                line = ax.plot(coupling_values_lambda, normalized_rel_times, '-D', 
                               color=color, linewidth=3,
                               markersize=8, markerfacecolor='white', markeredgecolor=color,
                               markeredgewidth=2.0,
                               label=latex_safe(f'$t_{{\\mathrm{{w}}}} = {tw_shifted}$ ps', f't_w = {tw_shifted} ps', self.use_latex))
                
                panel1_scatter_points.append(line[0])  # line returns a list of Line2D objects
                
                print(f"    tw={tw_shifted} ps: {len(normalized_rel_times)} valid coupling points")
        
        # Formatting exactly like original
        ax.set_xlabel(self.format_lambda_axis_label(self.use_latex), fontsize=16)
        ax.set_ylabel(latex_safe(r'$\tilde{\tau}_{\mathrm{s}}$', 'τ̃_s', self.use_latex), fontsize=16)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.tick_params(axis='both', which='major', labelsize=12)
        
        # Add horizontal dashed line at y=1 to show normalization reference
        ax.axhline(y=1, color='gray', linestyle='--', alpha=0.7, linewidth=1.5)
        
        # Set x-axis to start at 0; upper limit follows the profile's coupling range
        _lam_max = max(self.convert_epsilon_to_lambda(v) for v in self.coupling_values)
        ax.set_xlim(left=0, right=_lam_max * 1.1 if _lam_max > 0 else 0.15)
        
        # Set y-axis limits to match reference script (auto-scaling upper limit)
        ax.set_ylim(bottom=0.7)
        
        # Add colorbar for Panel 1 (waiting times) - above the plot
        if panel1_scatter_points:
            # Create a proper colorbar using the normalization and colormap
            import matplotlib.cm as cm
            sm1 = cm.ScalarMappable(norm=waiting_time_norm, cmap=waiting_time_cmap)
            sm1.set_array([])
            cbar1 = plt.colorbar(sm1, ax=ax, location='top', shrink=0.8, pad=0.1)
            cbar1.set_label(latex_safe(r'$t_{\mathrm{w}}$ (ps)', 't_w (ps)', self.use_latex), fontsize=12, labelpad=10)
            cbar1.ax.tick_params(labelsize=10)
    
    def _plot_relaxation_panel_styled(self, ax, relaxation_data):
        """Plot Panel 2: Relaxation time vs waiting times (styled like relaxation_analysis_filtered.pdf).""" 
        print("  Creating Panel 2: Relaxation time vs waiting times...")
        
        # Create color mapping for coupling strengths (Panel 2) - use λ values
        coupling_values_sorted = sorted(self.coupling_values)
        coupling_values_sorted_lambda = [self.convert_epsilon_to_lambda(val) for val in coupling_values_sorted]
        coupling_norm = Normalize(vmin=0, vmax=max(coupling_values_sorted_lambda))
        coupling_cmap = plt.colormaps.get_cmap('coolwarm')
        
        # Store data for Panel 2 colorbar
        panel2_scatter_points = []
        
        # Find ε=0 coupling for normalization
        epsilon_zero_coupling = None
        for coupling_val in self.coupling_values:
            if coupling_val == 0:
                epsilon_zero_coupling = coupling_val
                break
        
        if epsilon_zero_coupling is None:
            print("    Warning: No ε=0 coupling found! Cannot normalize data.")
            return
        
        # Plot normalized data for Panel 2
        for coupling_val in self.coupling_values:
            if coupling_val not in relaxation_data:
                continue
                
            data = relaxation_data[coupling_val]
            waiting_times = data['waiting_times']
            relaxation_times = data['relaxation_times']
            
            # Filter out tw=0 and shift waiting times by -200 ps
            valid_indices = [i for i, tw in enumerate(waiting_times) if tw > 0]
            if not valid_indices:
                continue
                
            filtered_waiting_times = [(waiting_times[i] - 200) for i in valid_indices]  # Shift by -200 ps
            filtered_relaxation_times = [relaxation_times[i] for i in valid_indices]
            
            # Get ε=0 reference data for normalization
            if epsilon_zero_coupling not in relaxation_data:
                print(f"    Warning: No ε=0 reference data found for normalization")
                continue
                
            ref_data = relaxation_data[epsilon_zero_coupling]
            ref_waiting_times = ref_data['waiting_times']
            ref_relaxation_times = ref_data['relaxation_times']
            
            # Normalize by ε=0 values at corresponding waiting times
            normalized_relaxation_times = []
            final_waiting_times = []
            
            for i, tw in enumerate(filtered_waiting_times):
                original_tw = tw + 200  # Convert back to original scale for lookup
                
                # Find corresponding ε=0 value
                ref_value = None
                for j, ref_tw in enumerate(ref_waiting_times):
                    if abs(ref_tw - original_tw) < 1e-6:  # Match within tolerance
                        ref_value = ref_relaxation_times[j]
                        break
                
                if ref_value is not None and ref_value > 0 and not np.isnan(filtered_relaxation_times[i]):
                    normalized_rel_time = filtered_relaxation_times[i] / ref_value
                    normalized_relaxation_times.append(normalized_rel_time)
                    final_waiting_times.append(tw)
            
            if normalized_relaxation_times:  # Only plot if we have data
                lambda_val = self.convert_epsilon_to_lambda(coupling_val)
                color = coupling_cmap(coupling_norm(lambda_val))
                
                # Plot with both lines and markers in one call (matching original style)
                line = ax.plot(final_waiting_times, normalized_relaxation_times, '-D',
                               color=color, linewidth=3,
                               markersize=8, markerfacecolor='white', markeredgecolor=color,
                               markeredgewidth=2.0,
                               label=self.format_lambda_label(lambda_val, self.use_latex))
                
                panel2_scatter_points.append(line[0])  # line returns a list of Line2D objects
                
                print(f"    Coupling λ={lambda_val:.3f}: {len(normalized_relaxation_times)} valid waiting time points")
        
        # Formatting exactly like original
        ax.set_xlabel(latex_safe(r'$t_{\mathrm{w}}$ (ps)', 't_w (ps)', self.use_latex), fontsize=16)
        ax.set_ylabel(latex_safe(r'$\tilde{\tau}_{\mathrm{s}}$', 'τ̃_s', self.use_latex), fontsize=16)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.legend(fontsize=10, loc='best')
        
        # Add horizontal dashed line at y=1 to show normalization reference
        ax.axhline(y=1, color='gray', linestyle='--', alpha=0.7, linewidth=1.5)
        
        # Set x-axis to start at 0 for Panel 2 (waiting times in ps, so different scale)
        ax.set_xlim(left=0, right=1800)  # Waiting times go up to ~1800 ps
        
        # Set y-axis limits to match reference script (auto-scaling upper limit)
        ax.set_ylim(bottom=0.7)
        
        
        # Colorbar for Panel 2 removed per user request (matching original)

    def _plot_material_time_panel(self, ax, measured_material_time_data, tn_material_time_data=None):
        """Plot material time: direct measured F(k,t) (solid) and TN prediction (dashed)."""
        print("  Creating material time panel (direct measured + TN prediction)...")
        if tn_material_time_data is None:
            tn_material_time_data = {}

        xi_values: list[float] = []
        t_max = 0.0

        for coupling_value in sorted(self.coupling_values):
            color = self.get_coupling_color(coupling_value)
            label = self.format_coupling_label(coupling_value)

            if coupling_value in measured_material_time_data:
                time_meas, xi_meas = measured_material_time_data[coupling_value]
                ax.plot(
                    time_meas,
                    xi_meas,
                    color=color,
                    linewidth=2.5,
                    alpha=0.9,
                    linestyle="-",
                    label=label,
                )
                xi_values.extend(xi_meas[np.isfinite(xi_meas)])
                t_max = max(t_max, float(np.nanmax(time_meas)))

            if coupling_value in tn_material_time_data:
                time_tn, xi_tn, _aging_rate = tn_material_time_data[coupling_value]
                if time_tn.size:
                    ax.plot(
                        time_tn,
                        xi_tn,
                        color=color,
                        linewidth=2.0,
                        alpha=0.85,
                        linestyle="--",
                    )
                    xi_values.extend(xi_tn[np.isfinite(xi_tn)])
                    t_max = max(t_max, float(np.nanmax(time_tn)))

        style_proxy = [
            plt.Line2D([0], [0], color="black", linewidth=2.5, linestyle="-", label=latex_safe(
                r"Measured $\xi_{\mathrm{S}}$ from $F(k,t)$",
                "Measured xi_S from F(k,t)",
                self.use_latex,
            )),
            plt.Line2D([0], [0], color="black", linewidth=2.0, linestyle="--", label=latex_safe(
                r"TN prediction $\xi_{\mathrm{TN}}$",
                "TN prediction xi_TN",
                self.use_latex,
            )),
        ]

        ax.set_xlabel(
            latex_safe(r"$t$ since cavity turn-on (ps)", "t since cavity turn-on (ps)", self.use_latex),
            fontsize=16,
        )
        ax.set_ylabel(latex_safe(r"Material time $\xi(t)$", "Material time xi(t)", self.use_latex), fontsize=16)
        ax.grid(True, alpha=0.3)
        ax.tick_params(axis="both", which="major", labelsize=14)
        ax.set_xlim(0, t_max * 1.02 if t_max > 0 else 2300)

        if xi_values:
            y_max = float(np.nanmax(xi_values))
            ax.set_ylim(0, y_max * 1.08)
        else:
            ax.set_ylim(bottom=0)

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=handles + style_proxy, fontsize=9, loc="upper left")

        inset_ax = inset_axes(ax, width="35%", height="35%", loc="upper right")
        for coupling_value in sorted(self.coupling_values):
            if coupling_value not in tn_material_time_data:
                continue
            time_tn, _xi_tn, aging_rate = tn_material_time_data[coupling_value]
            color = self.get_coupling_color(coupling_value)
            inset_ax.plot(
                time_tn,
                aging_rate,
                color=color,
                linewidth=1.8,
                alpha=0.85,
                linestyle="--",
            )

        inset_ax.set_xlabel(
            latex_safe(r"$t$ (ps)", "t (ps)", self.use_latex),
            fontsize=10,
        )
        inset_ax.set_ylabel(
            latex_safe(r"$\gamma(t)$ (ps$^{-1}$)", "gamma(t) (ps^-1)", self.use_latex),
            fontsize=10,
        )
        inset_ax.grid(True, alpha=0.3)
        inset_ax.tick_params(axis="both", which="major", labelsize=9)
        inset_ax.set_xlim(0, t_max * 1.02 if t_max > 0 else 2300)
        for spine in inset_ax.spines.values():
            spine.set_linewidth(1.0)
            spine.set_edgecolor("gray")
    
    def _plot_normalization_panel(self, ax, normalization_data):
        """Plot Panel 2: Relaxation time vs coupling strength (normalized by ε=0)."""
        print("  Creating Panel 2: Relaxation time vs coupling strength...")
        
        # Color scheme for waiting times (reference times)
        waiting_times = sorted(normalization_data.keys())
        tw_norm = Normalize(vmin=min(waiting_times), vmax=max(waiting_times))
        tw_cmap = plt.colormaps.get_cmap('viridis')
        
        # First pass: find ε=0 reference values for each waiting time
        reference_values = {}
        for tw in waiting_times:
            if tw in normalization_data:
                data = normalization_data[tw]
                coupling_strengths = data['coupling_strengths']
                relaxation_times = data['normalization_times']
                
                # Find ε=0 value (coupling_strength = 0)
                zero_mask = np.abs(coupling_strengths) < 1e-6
                if np.any(zero_mask):
                    reference_values[tw] = relaxation_times[zero_mask][0]
                else:
                    # Skip waiting times without ε=0 reference
                    print(f"    Warning: No ε=0 reference for tw={tw} ps, skipping this waiting time")
                    continue
        
        # Second pass: plot normalized data
        for tw in waiting_times:
            if tw not in normalization_data or tw not in reference_values:
                continue
                
            data = normalization_data[tw]
            ref_value = reference_values[tw]
            color = tw_cmap(tw_norm(tw))
            
            # Normalize by ε=0 value for this waiting time
            normalized_times = data['normalization_times'] / ref_value
            
            # Plot normalized relaxation time vs coupling strength with diamond markers
            ax.plot(data['coupling_strengths'], normalized_times,
                   '-D', color=color, linewidth=3, alpha=0.7,
                   markerfacecolor='white', markeredgecolor=color, 
                   markeredgewidth=2.0, markersize=10,
                   label=latex_safe(f'$t_{{\\mathrm{{w}}}} = {tw:.0f}$ ps', f't_w = {tw:.0f} ps', self.use_latex))
        
        # Formatting to match relaxation_analysis_filtered style
        ax.set_xlabel(self.format_coupling_axis_label(self.use_latex), fontsize=16)
        ax.set_ylabel(latex_safe(r'$\tilde{\tau}_{\mathrm{s}}$', 'τ̃_s', self.use_latex), fontsize=16)
        
        ax.grid(True, alpha=0.3, linestyle='--')
        # Remove legend from middle plot
        ax.tick_params(axis='both', which='major', labelsize=14)  # Match tick label size
        ax.set_xlim(*self.coupling_axis_xlim())
        
        # Add colorbar using make_axes_locatable for proper positioning
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("top", size="6%", pad=0.1)
        
        sm = cm.ScalarMappable(norm=tw_norm, cmap=tw_cmap)
        sm.set_array([])  # Fix gradient fill issue for PDF output
        cbar = plt.colorbar(sm, cax=cax, orientation='horizontal')
        cbar.set_label(latex_safe(r'$t_{\mathrm{w}}$ (ps)', 't_w (ps)', self.use_latex), fontsize=14)
        cbar.ax.tick_params(labelsize=12)
        cax.xaxis.set_label_position('top')
        cax.xaxis.tick_top()
    
    def _plot_relaxation_panel(self, ax, relaxation_data):
        """Plot Panel 3: Relaxation time vs waiting times (like Panel 2 of relaxation_time_comparison)."""
        print("  Creating Panel 3: Relaxation time vs waiting times...")
        
        # Normalize data by epsilon=0 values for each waiting time
        normalized_data = {}
        ref_coupling = 0.0  # Reference coupling for normalization
        
        if ref_coupling in relaxation_data:
            ref_data = relaxation_data[ref_coupling]
            ref_waiting_times = ref_data['waiting_times']
            ref_relaxation_times = ref_data['relaxation_times']
            
            # Create interpolation function for reference data
            ref_interp = interp1d(ref_waiting_times, ref_relaxation_times, 
                                bounds_error=False, fill_value='extrapolate')
            
            # Normalize all coupling data
            for coupling_value in sorted(relaxation_data.keys()):
                data = relaxation_data[coupling_value]
                waiting_times = data['waiting_times']
                relaxation_times = data['relaxation_times']
                
                # Get reference values at the same waiting times
                ref_values = ref_interp(waiting_times)
                
                # Normalize
                normalized_relaxation_times = relaxation_times / ref_values
                
                normalized_data[coupling_value] = {
                    'waiting_times': waiting_times,
                    'relaxation_times': normalized_relaxation_times
                }
        else:
            # If no reference coupling, use original data
            normalized_data = relaxation_data
        
        for coupling_value in sorted(normalized_data.keys()):
            data = normalized_data[coupling_value]
            color = self.get_coupling_color(coupling_value)
            
            # Format coupling value for legend with scientific notation
            if coupling_value == 0:
                label = latex_safe('$\\varepsilon = 0$', 'ε = 0', self.use_latex)
            else:
                exp_str = f'{coupling_value:.1e}'
                mantissa, exp_part = exp_str.split('e')
                exponent = int(exp_part)
                if self.use_latex:
                    label = f'$\\varepsilon = {mantissa} \\times 10^{{{exponent}}}$'
                else:
                    label = f'ε = {mantissa} × 10^{exponent}'
            
            # Set first value at tw=0 to 1
            relaxation_times = data['relaxation_times'].copy()
            relaxation_times[0] = 1.0
            
            ax.plot(data['waiting_times'], relaxation_times,
                   '-D', color=color, linewidth=3, alpha=0.7,
                   markerfacecolor='white', markeredgecolor=color,
                   markeredgewidth=2.0, markersize=10, label=label)
        
        # Add coupling turn-on line and pre-coupling shaded region (like original)
        ax.axvline(x=200, color='red', linestyle='--', alpha=0.7, linewidth=2,
                  label=latex_safe(r'$t_0 = 200$ ps', 't₀ = 200 ps', self.use_latex))
        ax.axvspan(0, 200, alpha=0.2, color='gray')
        
        # Formatting to match relaxation_analysis_filtered style
        ax.set_xlabel(latex_safe(r'$t_{\mathrm{w}}$ (ps)', 't_w (ps)', self.use_latex), fontsize=16)
        ax.set_ylabel(latex_safe(r'$\tilde{\tau}_{\mathrm{s}}$', 'τ̃_s', self.use_latex), fontsize=16)
        
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.legend(fontsize=10, loc='upper right')
        ax.tick_params(axis='both', which='major', labelsize=14)  # Match tick label size
        ax.set_xlim(0, 2100)
        
        # Add colorbar using make_axes_locatable for proper positioning
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("top", size="6%", pad=0.1)
        
        sm = cm.ScalarMappable(norm=self.norm, cmap=self.cmap)
        sm.set_array([])  # Fix gradient fill issue for PDF output
        cbar = plt.colorbar(sm, cax=cax, orientation='horizontal')
        cbar.set_label(self.format_coupling_axis_label(self.use_latex), fontsize=14)
        
        # Set colorbar ticks
        cbar.set_ticks(self.coupling_values)
        if self.axis == "lambda":
            cbar.set_ticklabels([f"{val:g}" if val != 0 else "0" for val in self.coupling_values_scaled])
        else:
            cbar.set_ticklabels([f"{val:.1f}" if val != 0 else "0" for val in self.coupling_values_scaled])
        cbar.ax.tick_params(labelsize=12)
        cax.xaxis.set_label_position('top')
        cax.xaxis.tick_top()


def main():
    """Main function."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Create fictive temperature analysis figures')
    parser.add_argument('--profile', default='paper', help='Dataset profile (default: paper)')
    parser.add_argument('--staged-root', type=Path, default=None, help='Override staged data root')
    parser.add_argument('--base_dir', default=None,
                       help='Base directory containing data files (default: profile staged_root)')
    parser.add_argument('--material-time-output', default='fictive_material_time_evolution.pdf',
                       help='Output filename for material time figure')
    parser.add_argument('--relaxation-output', default='fictive_two_panel_relaxation_analysis.pdf',
                       help='Output filename for two-panel relaxation figure')
    parser.add_argument('--three-panel', action='store_true',
                       help='Create legacy three-panel figure instead of separate figures')
    parser.add_argument('--three-panel-output', default='fictive_three_panel_analysis.pdf',
                       help='Output filename for three-panel figure (when --three-panel is used)')
    parser.add_argument('--no-latex', action='store_true',
                       help='Disable LaTeX rendering for text')
    
    args = parser.parse_args()

    from repro_bootstrap import setup_profile

    profile = setup_profile(args, default="paper")
    base_dir = args.base_dir
    if base_dir is None:
        base_dir = str(profile.staged_root if profile.name != "paper" else Path("."))

    # Create analyzer
    analyzer = FictiveThreePanelAnalyzer(
        base_dir=base_dir,
        use_latex=not args.no_latex,
        profile=profile if profile.name != "paper" else None,
    )
    
    if args.three_panel:
        # Create the legacy three-panel figure
        analyzer.create_three_panel_figure(output_file=args.three_panel_output)
        print("\nThree-panel analysis complete!")
        print("\nPanel descriptions:")
        print("  Panel 1: Material time ξ_fictive(t) evolution from fictive temperature analysis")
        print("  Panel 2: Normalization time vs coupling strength for different waiting times")  
        print("  Panel 3: Relaxation time prediction vs waiting time for different coupling strengths")
    else:
        # Create separate figures
        analyzer.create_material_time_figure(output_file=args.material_time_output)
        analyzer.create_two_panel_relaxation_figure(output_file=args.relaxation_output)
        
        print("\nSeparate figure analysis complete!")
        print("\nFigure descriptions:")
        print("  Figure 1 (Material Time): Material time ξ_fictive(t) evolution from fictive temperature analysis")
        print("  Figure 2 (Two-Panel Relaxation): Normalization time vs coupling strength + Relaxation time prediction vs waiting time")


if __name__ == "__main__":
    main()
