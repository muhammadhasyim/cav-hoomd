
# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Analysis and tracking tools for cavity molecular dynamics simulations."""

from typing import Optional, List, Dict, Union, Any, Tuple
import hoomd
import numpy as np
import logging
import sys
import os
import time
from pathlib import Path
import enum
import scipy.fft

# CuPy import with fallback for CPU/GPU agnostic code
try:
    import cupy as cp
    HAS_CUPY = True
except ImportError:
    cp = None
    HAS_CUPY = False

from ..utils import PhysicalConstants, unwrap_positions

class EmpiricalTemperatureData:
    """
    Handles empirical potential energy vs temperature data with extended function fitting.
    
    This class loads empirical energy-temperature relationships from equilibrium
    simulations and provides interpolation to determine systemic temperatures
    from instantaneous potential energies using extended functions:
    - Extended harmonic: E = aT/(1+bT) for harmonic energy
    - Extended T^(3/5): E = E₀ + AT^(3/5)/(1+CT^(3/5)) for LJ+Coulombic energy
    
    Parameters
    ----------
    data_file_path : str
        Path to file containing temperature and potential energy data
    energy_component : str, optional
        Which energy component to use ('lj_coulombic', 'harmonic', etc.). Default: 'lj_coulombic'
    use_direct_harmonic : bool, optional
        If True and energy_component='harmonic', use direct calculation T = 4*E/(N*kB). Default: False
    create_plots : bool, optional
        If True, create diagnostic plots showing fitted functions vs data. Default: False
        
    Attributes
    ----------
    temperatures : np.ndarray
        Temperature values from empirical data (K)
    energies : np.ndarray  
        Energy values from empirical data (hartree)
    has_extended_harmonic_fit : bool
        Whether extended harmonic fitting has been performed
    has_extended_t35_fit : bool
        Whether extended T^(3/5) fitting has been performed
    extended_harmonic_fit : dict
        Parameters from extended harmonic fitting: {'a': float, 'b': float, 'r2': float}
    extended_t35_fit : dict
        Parameters from extended T^(3/5) fitting: {'e0': float, 'a': float, 'c': float, 'r2': float}
    """
    
    def __init__(self, data_file_path: str, energy_component: str = 'lj_coulombic', use_direct_harmonic: bool = False, create_plots: bool = False):
        self.data_file_path = Path(data_file_path)
        self.energy_component = energy_component
        self.use_direct_harmonic = use_direct_harmonic
        self.create_plots = create_plots
        
        # Initialize extended fitting attributes (no more simple power law)
        self.has_extended_harmonic_fit = False
        self.has_extended_t35_fit = False
        self.extended_harmonic_fit = {}
        self.extended_t35_fit = {}
        
        self.load_empirical_data()
        
        # Use extended functions only (no more simple power law)
        if not self.use_direct_harmonic or self.energy_component != 'harmonic':
            if self.energy_component == 'harmonic':
                self.fit_extended_harmonic_function()
            else:
                self.fit_extended_t35_function()
            
            # Create diagnostic plots if requested
            if self.create_plots:
                self.plot_fits()
    
    def load_empirical_data(self):
        """Load empirical energy-temperature data from file."""
        if not self.data_file_path.exists():
            raise FileNotFoundError(f"Empirical data file not found: {self.data_file_path}")
        
        try:
            import pandas as pd
            data = pd.read_csv(self.data_file_path, sep=r'\s+', comment='#')
            
            # Extract temperature data
            if 'temperature' in data.columns:
                self.temperatures = data['temperature'].values
            else:
                raise ValueError("Temperature column not found in empirical data")
            
            # Extract energy data based on component
            if self.energy_component == 'lj_coulombic':
                if 'lj_hartree' in data.columns and 'coulombic_hartree' in data.columns:
                    self.energies = data['lj_hartree'].values + data['coulombic_hartree'].values
                else:
                    raise ValueError("LJ and Coulombic energy columns not found")
            elif self.energy_component == 'total_PE':
                if 'total_potential_energy_hartree' in data.columns:
                    self.energies = data['total_potential_energy_hartree'].values
                else:
                    raise ValueError("Total potential energy column not found")
            elif self.energy_component == 'harmonic':
                if 'harmonic_hartree' in data.columns:
                    self.energies = data['harmonic_hartree'].values
                else:
                    raise ValueError("Harmonic energy column not found")
            else:
                raise ValueError(f"Unknown energy component: {self.energy_component}")
                
            print(f"Loaded {len(self.temperatures)} empirical data points")
            print(f"Temperature range: {self.temperatures.min():.1f} - {self.temperatures.max():.1f} K")
            print(f"Energy range: {self.energies.min():.6f} - {self.energies.max():.6f} Ha")
            
        except Exception as e:
            raise RuntimeError(f"Failed to load empirical data: {e}")
    
    # Removed fit_rosenfeld_function - using only extended functions for better accuracy
    
    def fit_extended_harmonic_function(self):
        """Fit extended harmonic function: E = aT/(1+bT)"""
        try:
            from scipy.optimize import curve_fit
            
            def extended_harmonic_model(T, a, b):
                """Extended harmonic model: E(T) = aT/(1+bT)"""
                return a * T / (1 + b * T)
            
            # Initial parameter guesses
            a_guess = self.energies.max() / self.temperatures.max()  # Linear coefficient
            b_guess = 1e-3  # Small positive value
            
            # Perform fitting
            self.extended_harmonic_params, _ = curve_fit(
                extended_harmonic_model,
                self.temperatures,
                self.energies,
                p0=[a_guess, b_guess],
                bounds=([0, 0], [np.inf, np.inf])  # Both parameters must be positive
            )
            
            # Calculate R-squared
            predicted = extended_harmonic_model(self.temperatures, *self.extended_harmonic_params)
            residuals = self.energies - predicted
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((self.energies - np.mean(self.energies))**2)
            r2 = 1 - (ss_res / ss_tot)
            
            # Store fitted parameters
            self.extended_harmonic_fit = {
                'a': self.extended_harmonic_params[0],
                'b': self.extended_harmonic_params[1],
                'r2': r2
            }
            
            self.has_extended_harmonic_fit = True
            
            print(f"Extended harmonic fitting completed:")
            print(f"  a = {self.extended_harmonic_fit['a']:.6f} Ha·K^(-1)")
            print(f"  b = {self.extended_harmonic_fit['b']:.6f} K^(-1)")
            print(f"  R² = {self.extended_harmonic_fit['r2']:.12f}")
            
        except Exception as e:
            print(f"Warning: Extended harmonic fitting failed: {e}")
            self.has_extended_harmonic_fit = False
    
    def fit_extended_t35_function(self):
        """Fit extended T^(3/5) function: E = E0 + AT^(3/5)/(1+CT^(3/5))"""
        try:
            from scipy.optimize import curve_fit
            
            def extended_t35_model(T, e0, a, c):
                """Extended T^(3/5) model: E(T) = E₀ + AT^(3/5)/(1+CT^(3/5))"""
                t_power = T**(3/5)
                return e0 + a * t_power / (1 + c * t_power)
            
            # Initial parameter guesses
            e0_guess = self.energies.min()
            a_guess = (self.energies.max() - self.energies.min()) / (self.temperatures.max()**(3/5))
            c_guess = 1e-3
            
            # Perform fitting
            self.extended_t35_params, _ = curve_fit(
                extended_t35_model,
                self.temperatures,
                self.energies,
                p0=[e0_guess, a_guess, c_guess],
                bounds=([self.energies.min() - 10, -np.inf, 0], 
                       [self.energies.max() + 10, np.inf, np.inf])
            )
            
            # Calculate R-squared
            predicted = extended_t35_model(self.temperatures, *self.extended_t35_params)
            residuals = self.energies - predicted
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((self.energies - np.mean(self.energies))**2)
            r2 = 1 - (ss_res / ss_tot)
            
            # Store fitted parameters
            self.extended_t35_fit = {
                'e0': self.extended_t35_params[0],
                'a': self.extended_t35_params[1],
                'c': self.extended_t35_params[2],
                'r2': r2
            }
            
            self.has_extended_t35_fit = True
            
            print(f"Extended T^(3/5) fitting completed:")
            print(f"  E₀ = {self.extended_t35_fit['e0']:.6f} Ha")
            print(f"  A = {self.extended_t35_fit['a']:.6f} Ha·K^(-3/5)")
            print(f"  C = {self.extended_t35_fit['c']:.6f} K^(-3/5)")
            print(f"  R² = {self.extended_t35_fit['r2']:.12f}")
            
        except Exception as e:
            print(f"Warning: Extended T^(3/5) fitting failed: {e}")
            self.has_extended_t35_fit = False
    
    def plot_fits(self, output_file: str = None, show_plot: bool = False):
        """
        Plot empirical data and fitted functions for visual validation.
        
        Parameters
        ----------
        output_file : str, optional
            Output file path for saving plot. If None, uses component-based naming.
        show_plot : bool, optional
            Whether to display plot interactively. Default: False
        """
        try:
            import matplotlib.pyplot as plt
            import matplotlib.style
            
            # Set up scientific plotting style
            plt.style.use('default')
            plt.rcParams.update({
                'font.size': 12,
                'font.family': 'serif',
                'mathtext.fontset': 'cm',
                'axes.grid': True,
                'grid.alpha': 0.3,
                'figure.dpi': 150,
                'savefig.dpi': 300,
                'savefig.bbox': 'tight'
            })
            
            # Create figure with 2x2 subplots
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
            
            # Generate fine temperature array for smooth curves
            T_fine = np.linspace(self.temperatures.min(), self.temperatures.max(), 1000)
            
            # Plot 1: Linear scale
            ax1.scatter(self.temperatures, self.energies, color='blue', alpha=0.7, s=50, 
                       label='Empirical Data', zorder=3)
            
            # Plot 2: Log-log scale (same data)
            # Only plot positive values for log scale
            pos_mask = (self.temperatures > 0) & (self.energies > 0)
            if np.any(pos_mask):
                ax2.loglog(self.temperatures[pos_mask], self.energies[pos_mask], 'bo', 
                          alpha=0.7, markersize=6, label='Empirical Data')
            
            # Plot fitted functions on both linear and log-log scales
            if self.has_extended_harmonic_fit and self.energy_component == 'harmonic':
                a, b = self.extended_harmonic_fit['a'], self.extended_harmonic_fit['b']
                r2 = self.extended_harmonic_fit['r2']
                
                # Extended harmonic function
                E_fit = a * T_fine / (1 + b * T_fine)
                ax1.plot(T_fine, E_fit, 'r-', linewidth=2, alpha=0.8,
                        label=f'Extended Harmonic: $E = \\frac{{aT}}{{1+bT}}$\n'
                              f'a = {a:.6f} Ha/K, b = {b:.6f} K$^{{-1}}$\nR² = {r2:.6f}')
                
                # Linear approximation for comparison
                E_linear = a * T_fine
                ax1.plot(T_fine, E_linear, 'g--', linewidth=1, alpha=0.6,
                        label=f'Linear: $E = aT$ (a = {a:.6f})')
                
                # Log-log plots (only positive values)
                pos_fit_mask = E_fit > 0
                if np.any(pos_fit_mask):
                    ax2.loglog(T_fine[pos_fit_mask], E_fit[pos_fit_mask], 'r-', 
                              linewidth=2, alpha=0.8, label='Extended Harmonic')
                    ax2.loglog(T_fine, E_linear, 'g--', linewidth=1, alpha=0.6, 
                              label='Linear: $E \\propto T$')
                        
            elif self.has_extended_t35_fit and self.energy_component != 'harmonic':
                e0, a, c = self.extended_t35_fit['e0'], self.extended_t35_fit['a'], self.extended_t35_fit['c']
                r2 = self.extended_t35_fit['r2']
                
                # Extended T^(3/5) function
                t_power = T_fine**(3/5)
                E_fit = e0 + a * t_power / (1 + c * t_power)
                ax1.plot(T_fine, E_fit, 'r-', linewidth=2, alpha=0.8,
                        label=f'Extended T$^{{3/5}}$: $E = E_0 + \\frac{{AT^{{3/5}}}}{{1+CT^{{3/5}}}}$\n'
                              f'E₀ = {e0:.6f} Ha, A = {a:.6f} Ha·K$^{{-3/5}}$\n'
                              f'C = {c:.6f} K$^{{-3/5}}$, R² = {r2:.6f}')
                
                # Simple T^(3/5) for comparison if c is small
                if c < 1e-3:
                    E_simple = e0 + a * t_power
                    ax1.plot(T_fine, E_simple, 'g--', linewidth=1, alpha=0.6,
                            label=f'Simple T$^{{3/5}}$: $E = E_0 + AT^{{3/5}}$')
                
                # Log-log plots - need to handle negative energies
                # Shift data to make it positive for log plotting
                E_shift = E_fit - e0  # Remove baseline
                E_simple_shift = a * t_power if c < 1e-3 else None
                
                pos_shift_mask = E_shift > 0
                if np.any(pos_shift_mask):
                    ax2.loglog(T_fine[pos_shift_mask], E_shift[pos_shift_mask], 'r-', 
                              linewidth=2, alpha=0.8, 
                              label=f'Extended T$^{{3/5}}$ (shifted: $E - E_0$)')
                    
                    if E_simple_shift is not None:
                        ax2.loglog(T_fine, E_simple_shift, 'g--', linewidth=1, alpha=0.6,
                                  label='Pure T$^{{3/5}}$: $E \\propto T^{{3/5}}$')
            
            # Format linear plot
            ax1.set_xlabel('Temperature (K)')
            ax1.set_ylabel('Energy (Ha)')
            ax1.set_title(f'Linear Scale: {self.energy_component.replace("_", "+").title()}')
            ax1.legend(fontsize=9)
            ax1.grid(True, alpha=0.3)
            
            # Format log-log plot
            ax2.set_xlabel('Temperature (K)')
            ax2.set_ylabel('Energy (Ha)' if self.energy_component == 'harmonic' else 'Energy - E₀ (Ha)')
            ax2.set_title(f'Log-Log Scale: Power Law Analysis')
            ax2.legend(fontsize=9)
            ax2.grid(True, alpha=0.3)
            
            # Plot 3: Residuals
            if self.has_extended_harmonic_fit and self.energy_component == 'harmonic':
                a, b = self.extended_harmonic_fit['a'], self.extended_harmonic_fit['b']
                E_pred = a * self.temperatures / (1 + b * self.temperatures)
            elif self.has_extended_t35_fit and self.energy_component != 'harmonic':
                e0, a, c = self.extended_t35_fit['e0'], self.extended_t35_fit['a'], self.extended_t35_fit['c']
                t_power = self.temperatures**(3/5)
                E_pred = e0 + a * t_power / (1 + c * t_power)
            else:
                E_pred = np.interp(self.temperatures, self.temperatures, self.energies)  # No fit available
            
            residuals = self.energies - E_pred
            ax3.scatter(self.temperatures, residuals, color='red', alpha=0.7, s=50)
            ax3.axhline(y=0, color='black', linestyle='--', alpha=0.5)
            ax3.set_xlabel('Temperature (K)')
            ax3.set_ylabel('Residuals (Ha)')
            ax3.set_title('Fit Residuals')
            ax3.grid(True, alpha=0.3)
            
            # Add statistics text
            rmse = np.sqrt(np.mean(residuals**2))
            max_abs_residual = np.max(np.abs(residuals))
            ax3.text(0.05, 0.95, f'RMSE: {rmse:.2e} Ha\nMax |residual|: {max_abs_residual:.2e} Ha',
                    transform=ax3.transAxes, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            
            # Plot 4: Power law analysis
            if self.energy_component == 'harmonic':
                # For harmonic: analyze E vs T scaling
                if np.any(pos_mask):
                    # Fit power law: log(E) = log(A) + α*log(T)
                    log_T = np.log10(self.temperatures[pos_mask])
                    log_E = np.log10(self.energies[pos_mask])
                    
                    # Linear fit in log space
                    coeffs = np.polyfit(log_T, log_E, 1)
                    alpha_fit = coeffs[0]  # Power law exponent
                    log_A_fit = coeffs[1]  # log(prefactor)
                    A_fit = 10**log_A_fit
                    
                    # Calculate R²
                    log_E_pred = coeffs[0] * log_T + coeffs[1]
                    ss_res = np.sum((log_E - log_E_pred)**2)
                    ss_tot = np.sum((log_E - np.mean(log_E))**2)
                    r2_power = 1 - (ss_res / ss_tot)
                    
                    # Plot power law fit
                    T_power_fit = A_fit * T_fine**alpha_fit
                    ax4.loglog(T_fine, T_power_fit, 'purple', linewidth=2, alpha=0.8,
                              label=f'Power Law: $E = AT^{{\\alpha}}$\n'
                                    f'A = {A_fit:.2e} Ha, α = {alpha_fit:.3f}\n'
                                    f'R² = {r2_power:.6f}')
                    
                    # Add theoretical lines
                    E_theory_linear = (self.energies[pos_mask].mean() / self.temperatures[pos_mask].mean()) * T_fine
                    ax4.loglog(T_fine, E_theory_linear, 'k--', alpha=0.5, 
                              label='Linear: $E \\propto T^1$ (theory)')
                    
                    ax4.loglog(self.temperatures[pos_mask], self.energies[pos_mask], 'bo', 
                              alpha=0.7, markersize=6, label='Data')
                    
                    ax4.set_title(f'Power Law Fit: α = {alpha_fit:.3f} (expect 1.0)')
                    
            else:
                # For LJ+Coulombic: analyze (E - E₀) vs T^(3/5) scaling
                if self.has_extended_t35_fit:
                    e0, a, c = self.extended_t35_fit['e0'], self.extended_t35_fit['a'], self.extended_t35_fit['c']
                    E_shifted = self.energies - e0
                    pos_shifted_mask = E_shifted > 0
                    
                    if np.any(pos_shifted_mask):
                        T_shifted = self.temperatures[pos_shifted_mask]
                        E_shifted_pos = E_shifted[pos_shifted_mask]
                        
                        # Fit power law to shifted data
                        log_T = np.log10(T_shifted)
                        log_E_shift = np.log10(E_shifted_pos)
                        
                        coeffs = np.polyfit(log_T, log_E_shift, 1)
                        alpha_fit = coeffs[0]
                        log_A_fit = coeffs[1]
                        A_fit = 10**log_A_fit
                        
                        # Calculate R²
                        log_E_pred = coeffs[0] * log_T + coeffs[1]
                        ss_res = np.sum((log_E_shift - log_E_pred)**2)
                        ss_tot = np.sum((log_E_shift - np.mean(log_E_shift))**2)
                        r2_power = 1 - (ss_res / ss_tot)
                        
                        # Plot power law fit
                        T_power_fit = A_fit * T_fine**alpha_fit
                        ax4.loglog(T_fine, T_power_fit, 'purple', linewidth=2, alpha=0.8,
                                  label=f'Power Law: $(E-E_0) = AT^{{\\alpha}}$\n'
                                        f'A = {A_fit:.2e} Ha, α = {alpha_fit:.3f}\n'
                                        f'R² = {r2_power:.6f}')
                        
                        # Add theoretical T^(3/5) line
                        E_theory_t35 = a * T_fine**(3/5)
                        ax4.loglog(T_fine, E_theory_t35, 'k--', alpha=0.5,
                                  label='Theory: $E \\propto T^{3/5}$ (α = 0.6)')
                        
                        ax4.loglog(T_shifted, E_shifted_pos, 'bo', alpha=0.7, markersize=6, 
                                  label='Data (E - E₀)')
                        
                        ax4.set_title(f'Power Law Fit: α = {alpha_fit:.3f} (expect 0.6)')
            
            ax4.set_xlabel('Temperature (K)')
            ax4.set_ylabel('Energy (Ha)')
            ax4.legend(fontsize=9)
            ax4.grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            # Save plot
            if output_file is None:
                output_file = f'empirical_fit_{self.energy_component}.png'
            
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"✓ Saved empirical fit plot to {output_file}")
            
            if show_plot:
                plt.show()
            else:
                plt.close()
                
        except ImportError:
            print("Warning: matplotlib not available, cannot create fit plots")
        except Exception as e:
            print(f"Warning: Failed to create fit plot: {e}")
    
    def calculate_systemic_temperature(self, instantaneous_energy_hartree: float, num_particles: int = None) -> float:
        """
        Calculate systemic temperature from instantaneous potential energy.
        
        Uses extended fitting functions when available for higher accuracy:
        - Extended harmonic: E = aT/(1+bT) -> T = E/(a-bE) 
        - Extended T^(3/5): E = E₀ + AT^(3/5)/(1+CT^(3/5)) -> solve numerically
    
    Parameters
    ----------
        instantaneous_energy_hartree : float
            Measured potential energy in hartree
        num_particles : int, optional
            Number of particles (required for harmonic direct calculation)
            
        Returns
    -------
        float
            Systemic temperature in Kelvin
        """
        if self.use_direct_harmonic and self.energy_component == 'harmonic':
            # Direct harmonic calculation: T = 4*E/(N*kB)
            if num_particles is None:
                raise ValueError("num_particles required for direct harmonic calculation")
            
            # Convert energy from hartree to Kelvin: E_hartree * hartree_to_K
            from ..utils import PhysicalConstants
            energy_kelvin = instantaneous_energy_hartree * PhysicalConstants.hartree_to_kelvin()
            temperature = 4.0 * energy_kelvin / num_particles  # 4*E/(N*kB), kB=1 in these units
            return max(temperature, 0.0)  # Ensure positive temperature
        
        # Prioritize extended fits for better accuracy
        elif self.has_extended_harmonic_fit and self.energy_component == 'harmonic':
            # Extended harmonic: E = aT/(1+bT) -> T = E/(a-bE)
            a, b = self.extended_harmonic_fit['a'], self.extended_harmonic_fit['b']
            E = instantaneous_energy_hartree
            
            if E <= 0 or (a - b * E) <= 0:
                return 0.0
            
            try:
                temperature = E / (a - b * E)
                return max(temperature, 0.0)
            except (ValueError, ZeroDivisionError):
                return 0.0
        
        elif self.has_extended_t35_fit and self.energy_component != 'harmonic':
            # Extended T^(3/5): E = E₀ + AT^(3/5)/(1+CT^(3/5))
            # Solve numerically using scipy.optimize
            try:
                from scipy.optimize import fsolve
                
                e0, a, c = self.extended_t35_fit['e0'], self.extended_t35_fit['a'], self.extended_t35_fit['c']
                E = instantaneous_energy_hartree
                
                if E <= e0:
                    return 0.0
                
                def equation(T):
                    if T <= 0:
                        return float('inf')
                    t_power = T**(3/5)
                    return e0 + a * t_power / (1 + c * t_power) - E
                
                # Initial guess based on simple T^(3/5) scaling
                T_guess = max(((E - e0) / a) ** (5/3), 1.0) if a > 0 else 300.0
                
                solution = fsolve(equation, T_guess, full_output=True)
                temperature = solution[0][0]
                
                # Check if solution converged
                if solution[2] == 1 and temperature > 0:
                    return temperature
                else:
                    # Extended fit failed to converge, fall back to linear interpolation
                    pass
                    
            except Exception as e:
                print(f"Warning: Extended T^(3/5) temperature calculation failed: {e}")
                # Fall back to linear interpolation
                pass
        
        # Fallback to linear interpolation if no extended fits are available
        # This should rarely happen since we now always try to fit extended functions
        print(f"Warning: No extended fits available for {self.energy_component}, using linear interpolation")
        temperature = np.interp(instantaneous_energy_hartree, self.energies, self.temperatures)
        return max(temperature, 0.0)


