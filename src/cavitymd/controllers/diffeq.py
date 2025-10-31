"""
Advanced control systems for cavity molecular dynamics simulations.

This module provides optimal and robust control algorithms for temperature regulation
in coupled thermal systems, including:
- LQR (Linear-Quadratic Regulator) with integral action
- Kalman filtering for state estimation
- In-situ system identification
"""

from typing import Optional, Dict, Any, Tuple
import hoomd
import numpy as np
import json
from datetime import datetime
from pathlib import Path

try:
    from scipy.linalg import solve_discrete_are
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("Warning: scipy not available. LQR controller will not be functional.")

from ..utils import PhysicalConstants
from .empirical import EmpiricalTemperatureData

# Physical constants
k_B = 8.617333e-5  # Boltzmann constant in eV/K


class RelaxationTimeModel:
    """
    Temperature-dependent relaxation time model using dual-regime fitting.
    
    Implements DF theory of glassy dynamics with:
    - Arrhenius regime: ln(τ) = ln(τ₀) + Eₐ·(β-β₀) for T > T_onset  
    - Parabolic regime: ln(τ) = ln(τ₀) + Eₐ·(β-β₀) + J²·(β-β₀)² for T < T_onset
    
    Used for exact cancellation term β(t) = 1/T(T_s(t)) + 1/τ in PI control.
    """
    
    def __init__(self, data_file_path: str):
        """
        Initialize relaxation time model from data file.
        
        Parameters
        ----------
        data_file_path : str
            Path to relaxation_times_vs_temperature.txt file
        """
        self.data_file_path = data_file_path
        self.is_fitted = False
        self.T_onset = None
        self.fit_results = {}
        
        if data_file_path and Path(data_file_path).exists():
            self._load_and_fit_data()
        else:
            print(f"Warning: Relaxation data file not found: {data_file_path}")
            print("PI control will use fallback relaxation time calculation")
    
    def _load_and_fit_data(self):
        """Load relaxation time data and perform dual-regime fitting."""
        try:
            # Load data (skip header lines, only load temperature and relaxation time columns)
            data = np.loadtxt(self.data_file_path, skiprows=3, usecols=(0, 2))
            temperatures = data[:, 0]  # Temperature in K
            relaxation_times = data[:, 1]  # tau_relax in ps (now column 1 since we only loaded 2 columns)
            
            print(f"DEBUG RelaxationModel: Loaded {len(temperatures)} data points")
            print(f"DEBUG RelaxationModel: T range: {temperatures.min():.1f} - {temperatures.max():.1f} K")
            print(f"DEBUG RelaxationModel: τ range: {relaxation_times.min():.2f} - {relaxation_times.max():.2f} ps")
            
            # Perform dual-regime fitting
            self.T_onset, self.fit_results = self._find_onset_temperature(temperatures, relaxation_times)
            
            if not np.isnan(self.T_onset):
                self.is_fitted = True
                print(f"RelaxationTimeModel: Successfully fitted with T_onset = {self.T_onset:.1f} K")
                if self.fit_results:
                    print(f"DEBUG RelaxationModel: Arrhenius R² = {self.fit_results.get('arrhenius', {}).get('r2', 'N/A'):.3f}")
                    print(f"DEBUG RelaxationModel: Parabolic R² = {self.fit_results.get('parabolic', {}).get('r2', 'N/A'):.3f}")
                # Automatically save Arrhenius plot
                self.save_arrhenius_plot()
            else:
                print("RelaxationTimeModel: Could not fit dual-regime model")
                
        except Exception as e:
            print(f"RelaxationTimeModel: Error loading data: {e}")
            import traceback
            traceback.print_exc()
    
    def _fit_arrhenius_regime(self, beta: np.ndarray, log_tau: np.ndarray) -> Tuple[float, float, float]:
        """Fit Arrhenius regime: ln(τ) = ln(τ₀) + Eₐ·(β-β₀)"""
        try:
            from scipy.optimize import curve_fit
        except ImportError:
            print("Warning: scipy not available for relaxation time fitting, using fallback")
            return 1.0, np.mean(beta), np.mean(log_tau)
        
        def arrhenius_func(beta_arr, Ea, beta_0, ln_tau_0):
            return ln_tau_0 + Ea * (beta_arr - beta_0)
        
        # Initial guess: use linear fit as starting point
        coeffs = np.polyfit(beta, log_tau, 1)
        slope_guess, intercept_guess = coeffs
        beta_0_guess = np.mean(beta)
        ln_tau_0_guess = intercept_guess + slope_guess * beta_0_guess
        
        try:
            popt, _ = curve_fit(arrhenius_func, beta, log_tau, 
                               p0=[slope_guess, beta_0_guess, ln_tau_0_guess],
                               bounds=([-np.inf, np.min(beta), -np.inf], 
                                      [np.inf, np.max(beta), np.inf]))
            return popt[0], popt[1], popt[2]  # Ea, beta_0, ln_tau_0
        except:
            # Fallback to linear fit
            return slope_guess, beta_0_guess, ln_tau_0_guess

    def _fit_parabolic_regime(self, beta: np.ndarray, log_tau: np.ndarray, beta_0_fixed: float, 
                            arrhenius_params: Tuple[float, float]) -> Tuple[float, float, float]:
        """Fit parabolic regime: ln(τ) = ln(τ₀) + Eₐ·(β-β₀) + J²·(β-β₀)²"""
        try:
            from scipy.optimize import curve_fit
        except ImportError:
            return arrhenius_params[0], 1e-6, arrhenius_params[1]
        
        Ea_arr, ln_tau_0_arr = arrhenius_params
        
        def parabolic_func(beta_arr, Ea, J):
            # Force continuity: at β₀, parabolic = Arrhenius
            delta_beta = beta_arr - beta_0_fixed
            return ln_tau_0_arr + Ea * delta_beta + J**2 * delta_beta**2
        
        # Initial guess - start with Arrhenius Ea, small J
        delta_beta = beta - beta_0_fixed
        
        # For initial J guess, fit residuals after linear term
        linear_part = ln_tau_0_arr + Ea_arr * delta_beta
        residuals = log_tau - linear_part
        J_guess = np.mean(residuals / (delta_beta**2 + 1e-10))  # Avoid division by zero
        
        try:
            popt, _ = curve_fit(parabolic_func, beta, log_tau, 
                               p0=[Ea_arr, J_guess],  # Start with Arrhenius Ea
                               bounds=([-np.inf, 0],  # J must be positive 
                                      [np.inf, np.inf]))
            return popt[0], popt[1], ln_tau_0_arr  # Ea, J, ln_tau_0 (shared with Arrhenius)
        except:
            # Fallback: use Arrhenius Ea and estimate J
            return Ea_arr, max(J_guess, 1e-6), ln_tau_0_arr

    def _find_onset_temperature(self, temperatures: np.ndarray, relaxation_times: np.ndarray) -> Tuple[float, Dict]:
        """Find onset temperature where Arrhenius behavior transitions to parabolic."""
        # Filter valid data
        valid_mask = ~np.isnan(relaxation_times)
        T_valid = temperatures[valid_mask]
        tau_valid = relaxation_times[valid_mask]
        
        if len(T_valid) < 6:
            print("Warning: Not enough data points for onset temperature analysis")
            return np.nan, {}
        
        # Sort by temperature (ascending)
        sort_idx = np.argsort(T_valid)
        T_sorted = T_valid[sort_idx]
        tau_sorted = tau_valid[sort_idx]
        
        beta_sorted = 1.0 / (k_B * T_sorted)  # Inverse temperature β = 1/(k_B T) in eV⁻¹
        log_tau_sorted = np.log(tau_sorted)
        
        # Try different split points (minimum 5 points per regime for robust fitting)
        min_points = 5
        
        best_onset_idx = None
        best_total_r2 = -np.inf
        best_fits = {}
        for split_idx in range(min_points, len(T_sorted) - min_points + 1):
            # T_sorted is ascending (low to high T)
            # High T regime (Arrhenius): indices split_idx to end (high temperatures)
            # Low T regime (Parabolic): indices 0 to split_idx-1 (low temperatures)
            
            low_T_beta = beta_sorted[:split_idx]     # Low T = high β
            low_T_log_tau = log_tau_sorted[:split_idx]
            high_T_beta = beta_sorted[split_idx:]    # High T = low β  
            high_T_log_tau = log_tau_sorted[split_idx:]
            
            try:
                # Fit Arrhenius regime (high T, low β)
                Ea_arr, beta_0_arr, ln_tau_0_arr = self._fit_arrhenius_regime(high_T_beta, high_T_log_tau)
                
                # Calculate R² for Arrhenius fit
                high_T_pred = ln_tau_0_arr + Ea_arr * (high_T_beta - beta_0_arr)
                ss_res_high = np.sum((high_T_log_tau - high_T_pred)**2)
                ss_tot_high = np.sum((high_T_log_tau - np.mean(high_T_log_tau))**2)
                r2_high = 1 - (ss_res_high / (ss_tot_high + 1e-10))
                
                # Fit parabolic regime (low T, high β) with fixed β₀ and continuity constraint
                Ea_par, J_par, ln_tau_0_par = self._fit_parabolic_regime(low_T_beta, low_T_log_tau, beta_0_arr, 
                                                                       (Ea_arr, ln_tau_0_arr))
                
                # Calculate R² for parabolic fit
                delta_beta_low = low_T_beta - beta_0_arr
                low_T_pred = ln_tau_0_par + Ea_par * delta_beta_low + J_par**2 * delta_beta_low**2
                ss_res_low = np.sum((low_T_log_tau - low_T_pred)**2)
                ss_tot_low = np.sum((low_T_log_tau - np.mean(low_T_log_tau))**2)
                r2_low = 1 - (ss_res_low / (ss_tot_low + 1e-10))
                
                # Combined R² (weighted by number of points)
                n_high = len(high_T_beta)
                n_low = len(low_T_beta)
                total_r2 = (n_high * r2_high + n_low * r2_low) / (n_high + n_low)
                
                # Add penalty for very unbalanced splits (favor more balanced regimes)
                balance_ratio = min(n_high, n_low) / max(n_high, n_low)
                balance_penalty = 0.95 + 0.05 * balance_ratio  # Small penalty for imbalance
                adjusted_r2 = total_r2 * balance_penalty
                
                # Additional constraints: J > 0 and reasonable balance
                if J_par > 0 and adjusted_r2 > best_total_r2 and n_high >= min_points and n_low >= min_points:
                    best_total_r2 = adjusted_r2
                    best_onset_idx = split_idx
                    best_fits = {
                        'beta_0': beta_0_arr,
                        'arrhenius': {
                            'Ea': Ea_arr,
                            'ln_tau_0': ln_tau_0_arr,
                            'r2': r2_high,
                            'n_points': n_high,
                            'T_range': (T_sorted[split_idx], T_sorted[-1])
                        },
                        'parabolic': {
                            'Ea': Ea_par,
                            'J': J_par,
                            'ln_tau_0': ln_tau_0_par,
                            'r2': r2_low,
                            'n_points': n_low,
                            'T_range': (T_sorted[0], T_sorted[split_idx-1])
                        },
                        'total_r2': total_r2
                    }
                
            except Exception as e:
                continue
        
        if best_onset_idx is None:
            print("Warning: Could not find suitable onset temperature")
            return np.nan, {}
        
        # Calculate the true onset temperature from the mathematical reference point
        T_onset = 1.0 / (k_B * best_fits['beta_0'])  # Mathematical reference point for continuity
        
        return T_onset, best_fits
    
    def get_relaxation_time(self, temperature_K: float) -> float:
        """
        Get relaxation time T(T_s) at given temperature.
        
        Parameters
        ----------
        temperature_K : float
            Temperature in Kelvin
            
        Returns
        -------
        float
            Relaxation time in picoseconds
        """
        if not self.is_fitted:
            # Fallback: simple temperature-dependent relaxation time
            # Use rough estimate: τ ∝ exp(E/kT) with E ~ 0.1 eV
            tau = 100.0 * np.exp(0.1 / (k_B * temperature_K))
            return tau
        
        beta = 1.0 / (k_B * temperature_K)  # Inverse temperature in eV⁻¹
        beta_0 = self.fit_results['beta_0']
        
        if temperature_K > self.T_onset:
            # High temperature: Arrhenius regime
            arr_params = self.fit_results['arrhenius']
            ln_tau = arr_params['ln_tau_0'] + arr_params['Ea'] * (beta - beta_0)
            regime = "Arrhenius"
        else:
            # Low temperature: Parabolic regime
            par_params = self.fit_results['parabolic']
            delta_beta = beta - beta_0
            ln_tau = par_params['ln_tau_0'] + par_params['Ea'] * delta_beta + par_params['J']**2 * delta_beta**2
            regime = "Parabolic"
        
        tau = np.exp(ln_tau)  # Return τ in ps
        
        # Occasional debug output (every ~100 calls to avoid spam)
        if hasattr(self, '_call_count'):
            self._call_count += 1
        else:
            self._call_count = 1
            
        if self._call_count % 100 == 1:  # Print every 100 calls
            print(f"DEBUG RelaxationModel: T={temperature_K:.1f}K → τ={tau:.2f}ps ({regime}, T_onset={self.T_onset:.1f}K)")
        
        return tau

    def save_arrhenius_plot(self, output_path: str = "relaxation_arrhenius_fit.png"):
        """
        Generate and save Arrhenius plot showing the dual-regime fit.
        
        Parameters
        ----------
        output_path : str
            Path to save the PNG plot (default: "relaxation_arrhenius_fit.png")
        """
        if not self.is_fitted:
            print("RelaxationTimeModel: Cannot create plot - model not fitted")
            return
            
        try:
            import matplotlib.pyplot as plt
            
            # Load original data
            data = np.loadtxt(self.data_file_path, skiprows=3, usecols=(0, 2))
            temperatures = data[:, 0]  # Temperature in K
            relaxation_times = data[:, 1]  # tau_relax in ps
            
            # Create temperature range for smooth curves
            T_min, T_max = temperatures.min(), temperatures.max()
            T_plot = np.linspace(T_min, T_max, 200)
            beta_plot = 1.0 / (k_B * T_plot)
            
            # Calculate fitted curves
            tau_fit = np.array([self.get_relaxation_time(T) for T in T_plot])
            
            # Create Arrhenius plot: ln(τ) vs 1/T
            fig, ax = plt.subplots(figsize=(10, 6))
            
            # Plot experimental data
            beta_data = 1.0 / (k_B * temperatures)
            ax.semilogy(beta_data, relaxation_times, 'o', color='blue', markersize=6, 
                       alpha=0.7, label='Experimental Data')
            
            # Plot fitted curve
            ax.semilogy(beta_plot, tau_fit, '-', color='red', linewidth=2, 
                       label=f'Dual-Regime Fit (T_onset = {self.T_onset:.1f} K)')
            
            # Add onset temperature line
            beta_onset = 1.0 / (k_B * self.T_onset)
            tau_onset = self.get_relaxation_time(self.T_onset)
            ax.axvline(beta_onset, color='green', linestyle='--', alpha=0.7, 
                      label=f'T_onset = {self.T_onset:.1f} K')
            
            # Formatting
            ax.set_xlabel('1/T (eV⁻¹)', fontsize=12)
            ax.set_ylabel('Relaxation Time τ (ps)', fontsize=12)
            ax.set_title('Relaxation Time vs Temperature (Arrhenius Plot)', fontsize=14)
            ax.grid(True, alpha=0.3)
            ax.legend()
            
            # Add regime annotations
            if 'arrhenius' in self.fit_results and 'parabolic' in self.fit_results:
                arr_r2 = self.fit_results['arrhenius'].get('r2', 0)
                par_r2 = self.fit_results['parabolic'].get('r2', 0)
                textstr = f'Arrhenius R² = {arr_r2:.3f}\nParabolic R² = {par_r2:.3f}'
                props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
                ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
                       verticalalignment='top', bbox=props)
            
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"RelaxationTimeModel: Arrhenius plot saved to {output_path}")
            plt.close()
            
        except ImportError:
            print("RelaxationTimeModel: matplotlib not available - cannot create plot")
        except Exception as e:
            print(f"RelaxationTimeModel: Error creating plot: {e}")


class DiffEqController(hoomd.custom.Action):
    """
    Differential equation temperature controller for cavity MD simulations.
    
    This controller implements a first-order differential equation to control
    the bath temperature based on a signal temperature:
    
    dT_bath(t)/dt = -(T_bath(t) - T_signal(t)) / τ
    
    where:
    - T_bath(t) is the bath temperature
    - T_signal(t) is the measured signal temperature (kinetic, LJ+Coulombic, etc.)
    - τ is the time constant controlling response speed
    
    The discretized update rule is:
    T_bath(t+dt) = T_bath(t) + dt * (-(T_bath(t) - T_signal(t)) / τ)
    T_bath(t+dt) = T_bath(t) - (dt/τ) * (T_bath(t) - T_signal(t))
    
    Parameters
    ----------
    temperature_method : str
        Temperature calculation method for signal ('kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition', 'lj_coulombic_bath')
    time_constant_ps : float
        Time constant τ in picoseconds (controls response speed)
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing
    energy_tracker : EnergyTracker
        Energy tracker for temperature calculations
    molecular_thermostat : hoomd.md.methods.thermostats.Thermostat or None
        Molecular thermostat to control
    cavity_thermostat : hoomd.md.methods.thermostats.Thermostat or None
        Cavity thermostat to control
    apply_to : str, optional
        Which thermostats to control ('molecular', 'cavity', 'both'). Default: 'both'
    turn_on_time_ps : float, optional
        Time to start control in picoseconds. Default: 0.0
    turn_off_time_ps : float, optional
        Time to stop control in picoseconds. Default: None (never turn off)
    update_interval_ps : float, optional
        Interval between control updates in picoseconds. Default: 0.1
    T_min : float, optional
        Minimum allowed bath temperature in Kelvin. Default: 0.0
    T_max : float, optional
        Maximum allowed bath temperature in Kelvin. Default: None (no upper limit)
    rate_limit_K_per_ps : float, optional
        Maximum rate of temperature change in K/ps. Default: None (no rate limit)
    output_file : str, optional
        Output file for control data. Default: 'diffeq_feedback.csv'
    empirical_data_file : str, optional
        Path to empirical data file (required for certain methods)
    console_output_period_ps : float, optional
        Console output period in picoseconds. Default: 1.0
        
    Examples
    --------
    **Basic differential equation control following kinetic temperature:**
    
    >>> from hoomd.cavitymd.analysis import DiffEqController
    >>> 
    >>> diffeq_controller = DiffEqController(
    ...     temperature_method='kinetic',
    ...     time_constant_ps=5.0,  # 5 ps response time
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     apply_to='molecular'
    ... )
    
    **Fast response to LJ+Coulombic signal:**
    
    >>> diffeq_controller = DiffEqController(
    ...     temperature_method='lj_coulombic_bath',
    ...     time_constant_ps=1.0,  # Fast 1 ps response
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     cavity_thermostat=cavity_thermo,
    ...     empirical_data_file='lj_coulombic_calibration.txt',
    ...     apply_to='both'
    ... )
    
    Notes
    -----
    - Smaller time constants lead to faster response to signal changes
    - The bath temperature exponentially approaches the signal temperature
    - Unlike feedback controllers, this directly tracks the signal without a target
    - Compatible with all existing temperature calculation methods
    """
    
    def __init__(self,
                 temperature_method: str,
                 time_constant_ps: float,
                 time_tracker,
                 energy_tracker,
                 simulation=None,
                 time_constant_auto: bool = False,
                 molecular_thermostat=None,
                 cavity_thermostat=None,
                 apply_to: str = 'both',
                 turn_on_time_ps: float = 0.0,
                 turn_off_time_ps: Optional[float] = None,
                 update_interval_ps: float = 0.1,
                 T_min: float = 0.0,
                 T_max: Optional[float] = None,
                 rate_limit_K_per_ps: Optional[float] = None,
                 output_file: str = 'diffeq_feedback.csv',
                 empirical_data_file: Optional[str] = None,
                 console_output_period_ps: float = 1.0,
                 enable_bias_estimation: bool = True,
                 bias_process_noise: float = 1e-6,
                 bias_initial_covariance: float = 100.0,
                 # Exact-cancellation + PI control parameters
                 enable_pi_control: bool = False,
                 pi_rho: float = 1.0,
                 pi_epsilon: float = 0.5,
                 pi_zeta: float = 0.8,
                 relaxation_data_file: Optional[str] = None,
                 filter_window_ps: float = 0.0,
                 # Adaptive bias cancellation parameters (mutually exclusive with PI control)
                 enable_bias_cancellation: bool = False,
                 bias_tau_b_ps: float = 50.0,  # tau_b > tau_s for robustness  
                 bias_tau_b_auto: bool = False,  # compute tau_b = prefactor * tau_s from relaxation model
                 bias_kappa: float = 0.01,     # small for time-scale separation
                 bias_kappa_auto: bool = False,  # compute kappa = 1/(prefactor * tau_s) for proper hierarchy
                 bias_tau_b_prefactor: float = 5.0,  # prefactor for tau_b = prefactor * tau_s
                 bias_kappa_prefactor: float = 50.0,  # prefactor for kappa = 1/(prefactor * tau_s)
                 bias_calibration_time_ps: float = 10.0,  # time to estimate initial bias
                 enable_csv_output: bool = False):
        
        # Store core parameters
        self.temperature_method = temperature_method
        self.time_constant_ps = float(time_constant_ps)
        self.time_constant_auto = bool(time_constant_auto)
        self.time_tracker = time_tracker
        self.energy_tracker = energy_tracker
        self.simulation = simulation
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        
        # Control parameters
        self.apply_to = apply_to
        self.turn_on_time_ps = float(turn_on_time_ps)
        self.turn_off_time_ps = float(turn_off_time_ps) if turn_off_time_ps is not None else None
        self.update_interval_ps = float(update_interval_ps)
        self.T_min = float(T_min)
        self.T_max = float(T_max) if T_max is not None else None
        self.rate_limit_K_per_ps = float(rate_limit_K_per_ps) if rate_limit_K_per_ps is not None else None
        self.output_file = output_file
        self.console_output_period_ps = float(console_output_period_ps)
        self.enable_csv_output = enable_csv_output
        
        # Control history for HDF5 output
        self.control_history = []
        
        # Kalman filter parameters for bias estimation
        self.enable_bias_estimation = bool(enable_bias_estimation)
        self.bias_process_noise = float(bias_process_noise)  # Q_bias: process noise for bias (should be small)
        self.bias_initial_covariance = float(bias_initial_covariance)  # P_bias initial uncertainty
        
        # Kalman filter state
        self.bias_estimate = 0.0  # Estimated measurement bias (K)
        self.bias_covariance = self.bias_initial_covariance  # Uncertainty in bias estimate
        
        # Load empirical data if needed
        self.empirical_data = None
        if empirical_data_file is not None and temperature_method in ['lj_coulombic', 'harmonic', 'lj_coulombic_bath']:
            try:
                energy_component = 'lj_coulombic' if 'lj_coulombic' in temperature_method else 'harmonic'
                self.empirical_data = EmpiricalTemperatureData(
                    data_file_path=empirical_data_file,
                    energy_component=energy_component
                )
            except Exception as e:
                print(f"Warning: Could not load empirical data file {empirical_data_file}: {e}")
        
        # Exact-cancellation + PI control parameters
        self.enable_pi_control = bool(enable_pi_control)
        self.pi_rho = float(pi_rho)  # Expected disturbance rate (K/ps)
        self.pi_epsilon = float(pi_epsilon)  # Target tolerance (K)
        self.pi_zeta = float(pi_zeta)  # Damping coefficient (fixed at 0.8)
        
        # Calculate PI gains using design recipe
        self.k_i = self.pi_rho / self.pi_epsilon  # k_i ≥ ρ/ε
        self.k_p = 2.0 * self.pi_zeta * np.sqrt(self.k_i)  # k_p = 2ζ√k_i
        
        # PI control state
        self.integral_state = 0.0  # Integral state z
        self.last_error = 0.0  # For anti-windup
        
        # Adaptive bias cancellation parameters (must be set before RelaxationTimeModel)
        self.enable_bias_cancellation = enable_bias_cancellation
        self.bias_tau_b_ps = bias_tau_b_ps
        self.bias_tau_b_auto = bias_tau_b_auto
        self.bias_kappa = bias_kappa
        self.bias_kappa_auto = bias_kappa_auto
        self.bias_tau_b_prefactor = bias_tau_b_prefactor
        self.bias_kappa_prefactor = bias_kappa_prefactor
        self.bias_calibration_time_ps = bias_calibration_time_ps
        
        # Mutual exclusion check
        if self.enable_pi_control and self.enable_bias_cancellation:
            raise ValueError("PI control and bias cancellation cannot both be enabled")
        
        # Initialize relaxation time model
        self.relaxation_model = None
        if relaxation_data_file and (self.enable_pi_control or self.time_constant_auto or self.bias_tau_b_auto or self.bias_kappa_auto):
            self.relaxation_model = RelaxationTimeModel(relaxation_data_file)
            if self.enable_pi_control:
                print(f"PI gains calculated - k_i={self.k_i:.3f}, k_p={self.k_p:.3f} (ρ={self.pi_rho}, ε={self.pi_epsilon}, ζ={self.pi_zeta})")
            if self.bias_tau_b_auto:
                print(f"Auto tau_b mode enabled - will compute tau_b = 5 * tau_s from relaxation model")
            if self.bias_kappa_auto:
                print(f"Auto kappa mode enabled - will compute kappa = 1/(50 * tau_s) for proper time-scale hierarchy")
        elif self.enable_pi_control:
            print(f"PI control enabled but no relaxation data file provided - using fallback model")
        elif self.time_constant_auto:
            print(f"Adaptive time constant enabled but no relaxation data file provided - using fallback calculation")
        elif self.bias_tau_b_auto:
            print(f"Auto tau_b mode enabled but no relaxation data file provided - using fallback calculation")
        elif self.bias_kappa_auto:
            print(f"Auto kappa mode enabled but no relaxation data file provided - using fallback calculation")
        
        # Temperature filtering parameters
        self.filter_window_ps = filter_window_ps
        self.temperature_buffer = []  # Circular buffer for temperature measurements
        self.time_buffer = []        # Corresponding time stamps
        self.is_filter_ready = False  # Flag to indicate when we have enough measurements
        
        # Bias cancellation state (parameters already set above)
        self.bias_estimate = 0.0  # b_hat
        self.bias_calibration_samples = []
        self.is_calibration_complete = False
        
        # Control state
        self.is_active = False
        self.current_bath_temperature = None
        self.last_update_time = 0.0
        self.last_console_output_time = 0.0
        
        # Initialize output file (only if CSV enabled)
        if self.enable_csv_output:
            self._initialize_output_file()
        else:
            self.output_file = None
    
    def _get_effective_bias_tau_b(self, temperature_K: float) -> float:
        """
        Get effective tau_b for bias cancellation.
        
        If bias_tau_b_auto is enabled, computes tau_b = 5 * tau_s from relaxation model.
        Otherwise, uses the fixed bias_tau_b_ps value.
        
        Parameters
        ----------
        temperature_K : float
            Current bath temperature in Kelvin
            
        Returns
        -------
        float
            Effective tau_b in picoseconds
        """
        if not self.bias_tau_b_auto:
            return self.bias_tau_b_ps
            
        if self.relaxation_model and self.relaxation_model.is_fitted:
            tau_s = self.relaxation_model.get_relaxation_time(temperature_K)
            tau_b_auto = self.bias_tau_b_prefactor * tau_s  # Configurable prefactor
            return tau_b_auto
        else:
            # Fallback when no relaxation model available
            # Use a simple exponential estimate: tau_s ≈ 100 * exp(E_a / (k_B * T))
            # with E_a ≈ 0.1 eV for typical glass-forming liquids
            k_B_eV_per_K = 8.617333e-5  # Boltzmann constant in eV/K
            E_a_eV = 0.1  # Activation energy in eV
            tau_s_fallback = 100.0 * np.exp(E_a_eV / (k_B_eV_per_K * temperature_K))
            tau_b_fallback = self.bias_tau_b_prefactor * tau_s_fallback
            print(f"Warning: No relaxation model available for auto tau_b - using fallback: tau_b = {tau_b_fallback:.1f} ps (T = {temperature_K:.1f} K)")
            return tau_b_fallback

    def _get_effective_bias_kappa(self, temperature_K: float) -> float:
        """
        Get effective kappa for bias cancellation.
        
        If bias_kappa_auto is enabled, computes kappa = 1/(prefactor * tau_s) for proper time-scale hierarchy.
        Otherwise, uses the fixed bias_kappa value.
        
        Parameters
        ----------
        temperature_K : float
            Current bath temperature in Kelvin
            
        Returns
        -------
        float
            Effective kappa (dimensionless)
        """
        if not self.bias_kappa_auto:
            return self.bias_kappa
            
        if self.relaxation_model and self.relaxation_model.is_fitted:
            tau_s = self.relaxation_model.get_relaxation_time(temperature_K)
            kappa_auto = 1.0 / (self.bias_kappa_prefactor * tau_s)  # Configurable hierarchy: tau_s << tau_b << 1/kappa
            return kappa_auto
        else:
            # Fallback when no relaxation model available
            # Use a simple exponential estimate: tau_s ≈ 100 * exp(E_a / (k_B * T))
            # with E_a ≈ 0.1 eV for typical glass-forming liquids
            k_B_eV_per_K = 8.617333e-5  # Boltzmann constant in eV/K
            E_a_eV = 0.1  # Activation energy in eV
            tau_s_fallback = 100.0 * np.exp(E_a_eV / (k_B_eV_per_K * temperature_K))
            kappa_fallback = 1.0 / (self.bias_kappa_prefactor * tau_s_fallback)
            print(f"Warning: No relaxation model available for auto kappa - using fallback: kappa = {kappa_fallback:.6f} (T = {temperature_K:.1f} K)")
            return kappa_fallback

    def _get_effective_time_constant(self, temperature_K: float) -> float:
        """
        Get effective time constant for the differential equation controller.
        
        If time_constant_auto is enabled, uses T[T_bath] from relaxation model.
        Otherwise, uses the fixed time constant.
        
        Parameters
        ----------
        temperature_K : float
            Bath temperature in Kelvin
            
        Returns
        -------
        float
            Effective time constant in ps
        """
        if self.time_constant_auto and self.relaxation_model is not None:
            # Use adaptive time constant from relaxation model
            tau_adaptive = self.relaxation_model.get_relaxation_time(temperature_K)
            return tau_adaptive
        elif self.time_constant_auto:
            # Fallback when auto is enabled but no relaxation model
            # Use simple temperature-dependent estimate
            tau_fallback = 100.0 * np.exp(0.1 / (k_B * temperature_K))
            return tau_fallback
        else:
            # Use fixed time constant
            return self.time_constant_ps
    
    def _initialize_output_file(self):
        """Initialize CSV output file with headers."""
        try:
            with open(self.output_file, 'w', encoding='utf-8') as f:
                f.write("# Differential Equation Temperature Controller Output with Kalman Filter Bias Estimation\n")
                f.write("# time_ps: simulation time in picoseconds\n")
                f.write("# signal_temp_raw_K: raw measured signal temperature\n")
                f.write("# signal_temp_corrected_K: bias-corrected signal temperature\n")
                f.write("# bias_estimate_K: estimated measurement bias\n")
                f.write("# bias_covariance: uncertainty in bias estimate\n")
                f.write("# bath_temp_K: bath temperature setting\n")
                f.write("# dt_ps: timestep used for integration\n")
                f.write("# derivative_K_per_ps: dT_bath/dt calculated\n")
                f.write("# active: whether controller is active (1=active, 0=inactive)\n")
                f.write("time_ps,signal_temp_raw_K,signal_temp_corrected_K,bias_estimate_K,bias_covariance,"
                       f"bath_temp_K,dt_ps,derivative_K_per_ps,active\n")
        except Exception as e:
            print(f"Warning: Failed to initialize differential equation controller output file: {e}")
    
    def _calculate_system_temperature(self, current_time_ps: float) -> Optional[float]:
        """Calculate system temperature using the specified method."""
        # Use the same implementation as OffsetTemperatureController
        if self.temperature_method == 'kinetic':
            # Calculate kinetic temperature from particle velocities
            import numpy as np
            from ..utils import PhysicalConstants
            
            try:
                kinetic_energy = 0.0
                N_molecules = 0
                with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                    for i, ptype in enumerate(snap.particles.typeid):
                        # Exclude cavity particle (assume it's the last particle)
                        if i < len(snap.particles.mass) - 1:
                            mass = snap.particles.mass[i]
                            velocity = np.array(snap.particles.velocity[i])
                            ke_particle = 0.5 * mass * np.sum(velocity**2)
                            kinetic_energy += ke_particle
                            N_molecules += 1
                
                if N_molecules > 0:
                    kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                    temperature_K = (2.0/3.0) * kinetic_energy / (N_molecules * kB_hartree)
                    return max(temperature_K, 0.0)
                else:
                    return None
                    
            except Exception as e:
                print(f"Warning: Could not calculate kinetic temperature: {e}")
                return None
        
        elif self.temperature_method == 'harmonic_equipartition':
            # Calculate harmonic temperature using equipartition theorem
            try:
                energy_data = self.energy_tracker.get_instantaneous_energy()
                if energy_data is None:
                    return None
                
                harmonic_energy_hartree = energy_data.get('harmonic', 0.0)
                
                # Use equipartition theorem: E_harmonic = (3/2) * N * kB * T
                # Assume N=500 molecular particles (excluding cavity)
                N_molecules = 500
                from ..utils import PhysicalConstants
                kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                
                if harmonic_energy_hartree > 0:
                    temperature_K = (2.0/3.0) * harmonic_energy_hartree / (N_molecules * kB_hartree)
                    return max(temperature_K, 0.0)
                else:
                    return None
                    
            except Exception as e:
                print(f"Warning: Could not calculate harmonic equipartition temperature: {e}")
                return None
        
        elif self.temperature_method in ['lj_coulombic', 'lj_coulombic_bath']:
            # Calculate LJ+Coulombic fictive temperature using empirical data
            if self.empirical_data is None:
                print(f"Warning: No empirical data available for method '{self.temperature_method}'")
                return None
            
            try:
                energy_data = self.energy_tracker.get_instantaneous_energy()
                if energy_data is None:
                    return None
                
                lj_energy = energy_data.get('lj', 0.0)
                coulomb_energy = energy_data.get('coulombic', 0.0)
                total_lj_coulomb = lj_energy + coulomb_energy
                
                temperature_K = self.empirical_data.calculate_systemic_temperature(
                    instantaneous_energy_hartree=total_lj_coulomb,
                    num_particles=500
                )
                return max(temperature_K, 0.0)
                
            except Exception as e:
                print(f"Warning: Could not calculate LJ+Coulombic temperature: {e}")
                return None
        
        elif self.temperature_method == 'harmonic':
            # Calculate harmonic fictive temperature using empirical data
            if self.empirical_data is None:
                print(f"Warning: No empirical data available for method '{self.temperature_method}'")
                return None
            
            try:
                energy_data = self.energy_tracker.get_instantaneous_energy()
                if energy_data is None:
                    return None
                
                harmonic_energy = energy_data.get('harmonic', 0.0)
                temperature_K = self.empirical_data.calculate_systemic_temperature(
                    instantaneous_energy_hartree=harmonic_energy,
                    num_particles=500
                )
                return max(temperature_K, 0.0)
                
            except Exception as e:
                print(f"Warning: Could not calculate harmonic temperature: {e}")
                return None
        
        else:
            print(f"Warning: Unknown temperature method '{self.temperature_method}'")
            return None
    
    def _get_current_thermostat_temperature(self) -> float:
        """Get current thermostat temperature for initialization."""
        try:
            from ..utils import PhysicalConstants
            
            # Try molecular thermostat first
            if self.molecular_thermostat is not None:
                try:
                    if hasattr(self.molecular_thermostat, 'kT'):
                        kT = self.molecular_thermostat.kT
                        if hasattr(kT, 'value'):
                            kT = kT.value
                        elif callable(kT):
                            kT = kT()
                        return kT / PhysicalConstants.KB_HARTREE_PER_K
                    elif hasattr(self.molecular_thermostat, 'thermostat'):
                        nested = self.molecular_thermostat.thermostat
                        if hasattr(nested, 'kT'):
                            kT = nested.kT
                            if hasattr(kT, 'value'):
                                kT = kT.value
                            elif callable(kT):
                                kT = kT()
                            return kT / PhysicalConstants.KB_HARTREE_PER_K
                except:
                    pass
            
            # Try cavity thermostat as backup
            if self.cavity_thermostat is not None:
                try:
                    if hasattr(self.cavity_thermostat, 'kT'):
                        kT = self.cavity_thermostat.kT
                        if hasattr(kT, 'value'):
                            kT = kT.value
                        elif callable(kT):
                            kT = kT()
                        return kT / PhysicalConstants.KB_HARTREE_PER_K
                    elif hasattr(self.cavity_thermostat, 'thermostat'):
                        nested = self.cavity_thermostat.thermostat
                        if hasattr(nested, 'kT'):
                            kT = nested.kT
                            if hasattr(kT, 'value'):
                                kT = kT.value
                            elif callable(kT):
                                kT = kT()
                            return kT / PhysicalConstants.KB_HARTREE_PER_K
                except:
                    pass
            
            # Fallback to 300K
            return 300.0
            
        except Exception as e:
            print(f"Warning: Could not get current thermostat temperature: {e}")
            return 300.0
    

    def _apply_rate_limit(self, new_output: float, dt_ps: float) -> float:
        """Apply rate limiting to controller output."""
        if self.rate_limit_K_per_ps is None or self.current_bath_temperature is None:
            return new_output
        
        max_change = self.rate_limit_K_per_ps * dt_ps
        change = new_output - self.current_bath_temperature
        
        if abs(change) <= max_change:
            return new_output
        else:
            # Limit the change
            limited_change = max_change if change > 0 else -max_change
            return self.current_bath_temperature + limited_change
    
    def _update_thermostats(self, bath_temperature: float):
        """Update thermostat temperatures."""
        try:
            import hoomd
            from ..utils import PhysicalConstants
            kT_hartree = bath_temperature * PhysicalConstants.KB_HARTREE_PER_K
            
            if self.apply_to in ['molecular', 'both'] and self.molecular_thermostat is not None:
                if hasattr(self.molecular_thermostat, 'kT'):
                    self.molecular_thermostat.kT = hoomd.variant.Constant(kT_hartree)
                elif hasattr(self.molecular_thermostat, 'thermostat'):
                    nested = self.molecular_thermostat.thermostat
                    if hasattr(nested, 'kT'):
                        nested.kT = hoomd.variant.Constant(kT_hartree)
            
            if self.apply_to in ['cavity', 'both'] and self.cavity_thermostat is not None:
                if hasattr(self.cavity_thermostat, 'kT'):
                    self.cavity_thermostat.kT = hoomd.variant.Constant(kT_hartree)
                elif hasattr(self.cavity_thermostat, 'thermostat'):
                    nested = self.cavity_thermostat.thermostat
                    if hasattr(nested, 'kT'):
                        nested.kT = hoomd.variant.Constant(kT_hartree)
                        
        except Exception as e:
            print(f"Warning: Failed to update thermostats: {e}")
            import traceback
            traceback.print_exc()
    
    def _update_temperature_filter(self, temperature_K: float, current_time_ps: float):
        """
        Update the moving window filter with a new temperature measurement.
        
        Args:
            temperature_K: New temperature measurement in K
            current_time_ps: Current simulation time in ps
        """
        if self.filter_window_ps <= 0.0:
            # No filtering - use direct measurement
            self.filtered_temperature = temperature_K
            self.is_filter_ready = True
            return
        
        # Add new measurement to buffers
        self.temperature_buffer.append(temperature_K)
        self.time_buffer.append(current_time_ps)
        
        # Check if we have enough measurements to start filtering BEFORE removing old ones
        if self.time_buffer and (current_time_ps - self.time_buffer[0]) >= self.filter_window_ps:
            self.is_filter_ready = True
            # Calculate filtered temperature as simple moving average
            self.filtered_temperature = sum(self.temperature_buffer) / len(self.temperature_buffer)
        else:
            self.is_filter_ready = False
        
        # Remove old measurements outside the time window (after checking readiness)
        cutoff_time = current_time_ps - self.filter_window_ps
        while (self.time_buffer and self.time_buffer[0] < cutoff_time):
            self.temperature_buffer.pop(0)
            self.time_buffer.pop(0)
    
    def act(self, timestep):
        """Execute differential equation control at each timestep."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if we should be active
        should_be_active = (current_time_ps >= self.turn_on_time_ps and 
                          (self.turn_off_time_ps is None or current_time_ps < self.turn_off_time_ps))
        
        if should_be_active and not self.is_active:
            self.is_active = True
            # Reset PI control state when turning on
            self.integral_state = 0.0
            self.last_error = 0.0
            
            # Initialize current_bath_temperature from thermostats when first activating
            if self.current_bath_temperature is None:
                self.current_bath_temperature = self._get_current_thermostat_temperature()
                print(f"Differential equation controller turned ON at t = {current_time_ps:.2f} ps")
                print(f"Initial bath temperature: {self.current_bath_temperature:.2f} K")
                if self.time_constant_auto:
                    effective_tau = self._get_effective_time_constant(self.current_bath_temperature)
                    print(f"Time constant: ADAPTIVE τ = T[T_bath] = {effective_tau:.2f} ps (auto mode)")
                    if self.relaxation_model is None:
                        print(f"  Warning: No relaxation model - using fallback τ(T) calculation")
                else:
                    print(f"Time constant: FIXED τ = {self.time_constant_ps:.2f} ps")
                if self.enable_pi_control:
                    print(f"PI control ENABLED: k_p={self.k_p:.3f}, k_i={self.k_i:.3f}")
                    if self.relaxation_model and self.relaxation_model.is_fitted:
                        print(f"Relaxation model loaded with T_onset={self.relaxation_model.T_onset:.1f}K")
            else:
                print(f"Differential equation controller turned ON at t = {current_time_ps:.2f} ps")
        elif not should_be_active and self.is_active:
            self.is_active = False
            # Reset PI control state when turning off
            self.integral_state = 0.0
            self.last_error = 0.0
            print(f"Differential equation controller turned OFF at t = {current_time_ps:.2f} ps")
        
        # Skip if not active or not initialized
        if not self.is_active or self.current_bath_temperature is None:
            return
        
        # Only update at specified intervals
        if current_time_ps - self.last_update_time < self.update_interval_ps:
            return
        
        # Calculate dt_ps BEFORE updating last_update_time
        dt_ps = current_time_ps - self.last_update_time if self.last_update_time > 0 else self.update_interval_ps
        self.last_update_time = current_time_ps
        
        # Calculate raw signal temperature (may contain bias)
        signal_temperature_raw = self._calculate_system_temperature(current_time_ps)
        
        if signal_temperature_raw is None:
            return
        
        # Update temperature filter with raw measurement
        self._update_temperature_filter(signal_temperature_raw, current_time_ps)
        
        # If filtering is enabled but not ready, delay control action
        if self.filter_window_ps > 0.0 and not self.is_filter_ready:
            if current_time_ps - self.last_console_output_time >= self.console_output_period_ps:
                remaining_time = self.filter_window_ps - (current_time_ps - (self.time_buffer[0] if self.time_buffer else current_time_ps))
                print(f"DiffEq Controller | Collecting measurements for filter (window={self.filter_window_ps:.1f}ps, remaining={remaining_time:.1f}ps)")
                self.last_console_output_time = current_time_ps
            return
        
        # Use filtered temperature if available, otherwise use raw
        if hasattr(self, 'filtered_temperature'):
            signal_temperature_corrected = self.filtered_temperature
        else:
            signal_temperature_corrected = signal_temperature_raw
        
        # Calculate error signal e = T_signal - T_bath
        error = signal_temperature_corrected - self.current_bath_temperature
        
        if self.enable_pi_control:
            # Exact-cancellation + PI control
            # v = β(t)*e + k_p*e + k_i*z, where β(t) = 1/T(T_s) + 1/τ
            
            # Calculate exact cancellation term β(t)
            if self.relaxation_model is not None:
                relaxation_time_ps = self.relaxation_model.get_relaxation_time(signal_temperature_corrected)
                print(f"DEBUG PI: T_signal={signal_temperature_corrected:.2f}K, τ_relax={relaxation_time_ps:.2f}ps (from model)")
            else:
                # Fallback: use simple temperature-dependent estimate
                relaxation_time_ps = 100.0 * np.exp(0.1 / (k_B * signal_temperature_corrected))
                print(f"DEBUG PI: T_signal={signal_temperature_corrected:.2f}K, τ_relax={relaxation_time_ps:.2f}ps (fallback)")
            
            # Get effective time constant (fixed or adaptive)
            effective_time_constant = self._get_effective_time_constant(self.current_bath_temperature)
            beta_t = 1.0 / relaxation_time_ps + 1.0 / effective_time_constant
            
            # PI control terms
            proportional_term = self.k_p * error
            integral_term = self.k_i * self.integral_state
            exact_cancellation_term = beta_t * error
            
            # Control law: v = β(t)*e + k_p*e + k_i*z
            derivative_K_per_ps = exact_cancellation_term + proportional_term + integral_term
            
            print(f"DEBUG PI: e={error:+.2f}K, β(t)={beta_t:.4f}, β·e={exact_cancellation_term:+.3f}, k_p·e={proportional_term:+.3f}, k_i·z={integral_term:+.3f}")
            print(f"DEBUG PI: Total dT/dt={derivative_K_per_ps:+.3f}K/ps, T_bath={self.current_bath_temperature:.2f}K")
            if self.time_constant_auto:
                print(f"DEBUG PI: Adaptive τ_ctrl={effective_time_constant:.2f}ps (from T[T_bath={self.current_bath_temperature:.1f}K])")
            
            # Update integral state: ż = e (with anti-windup)
            # Only integrate if not saturated
            raw_output_test = self.current_bath_temperature + dt_ps * derivative_K_per_ps
            will_saturate = (self.T_max is not None and raw_output_test > self.T_max) or (raw_output_test < self.T_min)
            
            if not will_saturate:
                self.integral_state += dt_ps * error
                print(f"DEBUG PI: Integral updated: z={self.integral_state:.3f}")
            else:
                print(f"DEBUG PI: Anti-windup active - integral frozen at z={self.integral_state:.3f}")
            # Anti-windup: freeze integrator when saturated
            
            raw_output = self.current_bath_temperature + dt_ps * derivative_K_per_ps
            
            # Additional safety check
            if abs(derivative_K_per_ps) > 1000.0:  # Extremely large derivative
                print(f"WARNING PI: Excessive derivative {derivative_K_per_ps:.1f}K/ps - clamping to ±100K/ps")
                derivative_K_per_ps = np.sign(derivative_K_per_ps) * 100.0
                raw_output = self.current_bath_temperature + dt_ps * derivative_K_per_ps
            
        elif self.enable_bias_cancellation:
            # Adaptive bias cancellation: 
            # dT_b/dt = -((T_hat_s - b_hat) - T_b) / tau_b
            # db_hat/dt = kappa * (T_hat_s - T_b)
            
            # Check if we're in calibration phase
            if not self.is_calibration_complete:
                if current_time_ps < self.bias_calibration_time_ps:
                    # Collect calibration samples
                    bias_sample = signal_temperature_corrected - self.current_bath_temperature
                    self.bias_calibration_samples.append(bias_sample)
                    
                    # No control action during calibration - just track signal
                    derivative_K_per_ps = 0.0
                    raw_output = self.current_bath_temperature
                    
                    print(f"Bias Calibration | T_s={signal_temperature_corrected:.1f}K, T_bath={self.current_bath_temperature:.1f}K, bias_sample={bias_sample:+.2f}K, samples={len(self.bias_calibration_samples)}")
                else:
                    # End calibration phase - estimate initial bias
                    if self.bias_calibration_samples:
                        self.bias_estimate = np.mean(self.bias_calibration_samples)
                        print(f"Bias calibration complete: initial b_hat = {self.bias_estimate:.2f}K (from {len(self.bias_calibration_samples)} samples)")
                    else:
                        self.bias_estimate = 0.0
                        print(f"Bias calibration complete: no samples collected, b_hat = 0.0K")
                    self.is_calibration_complete = True
                    
                    # Start bias cancellation
                    derivative_K_per_ps = 0.0  # No change at transition
                    raw_output = self.current_bath_temperature
            else:
                # Bias cancellation active - implement coupled differential equations
                # dT_b/dt = -((T_s - b_hat) - T_b) / tau_b
                # db_hat/dt = kappa * (T_s - T_b)
                
                # Get effective parameters (auto mode or fixed)
                effective_tau_b = self._get_effective_bias_tau_b(self.current_bath_temperature)
                effective_kappa = self._get_effective_bias_kappa(self.current_bath_temperature)
                
                # Update bias estimate: db_hat/dt = kappa * (T_s - T_b)
                db_hat_dt = effective_kappa * (signal_temperature_corrected - self.current_bath_temperature)
                self.bias_estimate += dt_ps * db_hat_dt
                
                # CRITICAL FIX: Saturate bias estimate to prevent runaway heating
                # The bias should be reasonable - typically within ±50K of zero
                max_bias = 50.0  # Maximum reasonable bias magnitude (K)
                self.bias_estimate = max(-max_bias, min(max_bias, self.bias_estimate))
                
                # Control law: dT_b/dt = -((T_s - b_hat) - T_b) / tau_b
                unbiased_signal = signal_temperature_corrected - self.bias_estimate
                
                # CRITICAL FIX: Saturate unbiased signal to prevent extreme temperatures
                # Keep unbiased signal within reasonable bounds of initial temperature (300K)
                initial_temp = 300.0  # Your initial temperature
                max_temp_deviation = 100.0  # Maximum allowed deviation (K)
                min_allowed_temp = initial_temp - max_temp_deviation  # 200K
                max_allowed_temp = initial_temp + max_temp_deviation  # 400K
                unbiased_signal = max(min_allowed_temp, min(max_allowed_temp, unbiased_signal))
                
                derivative_K_per_ps = -(self.current_bath_temperature - unbiased_signal) / effective_tau_b
                raw_output = self.current_bath_temperature + dt_ps * derivative_K_per_ps
                
                # Print debug information with auto mode indicators
                auto_tau_b_str = f"Auto τ_b={effective_tau_b:.1f}ps" if self.bias_tau_b_auto else f"τ_b={effective_tau_b:.1f}ps"
                auto_kappa_str = f"Auto κ={effective_kappa:.6f}" if self.bias_kappa_auto else f"κ={effective_kappa:.6f}"
                
                # Check if saturation was applied
                bias_saturated = abs(self.bias_estimate) >= max_bias - 0.1  # Close to saturation
                temp_saturated = (unbiased_signal <= min_allowed_temp + 0.1) or (unbiased_signal >= max_allowed_temp - 0.1)
                saturation_warning = " [SATURATED]" if (bias_saturated or temp_saturated) else ""
                
                if self.bias_tau_b_auto or self.bias_kappa_auto:
                    print(f"Bias Cancellation ({auto_tau_b_str}, {auto_kappa_str}) | T_s={signal_temperature_corrected:.1f}K, b_hat={self.bias_estimate:+.2f}K, T_s_unbiased={unbiased_signal:.1f}K, T_bath={self.current_bath_temperature:.1f}K, dT_b/dt={derivative_K_per_ps:+.3f}K/ps, db_hat/dt={db_hat_dt:+.4f}K/ps{saturation_warning}")
                else:
                    print(f"Bias Cancellation | T_s={signal_temperature_corrected:.1f}K, b_hat={self.bias_estimate:+.2f}K, T_s_unbiased={unbiased_signal:.1f}K, T_bath={self.current_bath_temperature:.1f}K, dT_b/dt={derivative_K_per_ps:+.3f}K/ps, db_hat/dt={db_hat_dt:+.4f}K/ps{saturation_warning}")
        
        else:
            # Original differential equation: dT_bath/dt = -(T_bath - T_signal_corrected) / τ
            # Discretized: T_bath(t+dt) = T_bath(t) + dt * (-(T_bath(t) - T_signal(t)) / τ)
            effective_time_constant = self._get_effective_time_constant(self.current_bath_temperature)
            derivative_K_per_ps = -(self.current_bath_temperature - signal_temperature_corrected) / effective_time_constant
            raw_output = self.current_bath_temperature + dt_ps * derivative_K_per_ps
            
            if self.time_constant_auto:
                print(f"DEBUG DiffEq: Adaptive τ={effective_time_constant:.2f}ps (from T[T_bath={self.current_bath_temperature:.1f}K]), dT/dt={derivative_K_per_ps:+.3f}K/ps")
        
        # Apply saturation limits
        if self.T_max is not None:
            saturated_output = max(self.T_min, min(self.T_max, raw_output))
        else:
            saturated_output = max(self.T_min, raw_output)
        
        # Apply rate limiting
        rate_limited_output = self._apply_rate_limit(saturated_output, dt_ps)
        
        # Update current bath temperature and thermostats
        self.current_bath_temperature = rate_limited_output
        self._update_thermostats(self.current_bath_temperature)
        
        # Console output
        if current_time_ps - self.last_console_output_time >= self.console_output_period_ps:
            if self.filter_window_ps > 0.0:
                # Show both raw and filtered temperatures
                filter_info = f"[n={len(self.temperature_buffer)}, Δ={signal_temperature_raw - signal_temperature_corrected:+.1f}K]"
                if self.enable_pi_control:
                    print(f"DiffEq+PI Controller | {self.temperature_method}→{signal_temperature_raw:.1f}K→{signal_temperature_corrected:.1f}K {filter_info} | {self.current_bath_temperature:.1f}K | e={error:+.1f}K | z={self.integral_state:.2f} | dT/dt: {derivative_K_per_ps:+.2f}K/ps")
                else:
                    print(f"DiffEq Controller | {self.temperature_method}→{signal_temperature_raw:.1f}K→{signal_temperature_corrected:.1f}K {filter_info} | {self.current_bath_temperature:.1f}K | dT/dt: {derivative_K_per_ps:+.2f}K/ps")
            else:
                # No filtering - show original format
                if self.enable_pi_control:
                    print(f"DiffEq+PI Controller | {self.temperature_method}→{signal_temperature_raw:.1f}K | {signal_temperature_corrected:.1f}K | {self.current_bath_temperature:.1f}K | e={error:+.1f}K | z={self.integral_state:.2f} | dT/dt: {derivative_K_per_ps:+.2f}K/ps")
                else:
                    print(f"DiffEq Controller | {self.temperature_method}→{signal_temperature_raw:.1f}K | {signal_temperature_corrected:.1f}K | {self.current_bath_temperature:.1f}K | dT/dt: {derivative_K_per_ps:+.2f}K/ps")
            self.last_console_output_time = current_time_ps

