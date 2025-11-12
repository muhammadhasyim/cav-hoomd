"""
Integration of FDR Temperature Estimator with CavityMD Simulation Framework

This module provides seamless integration of the fluctuation-dissipation ratio (FDR)
based effective temperature estimator with the existing CavityMD simulation infrastructure.

Key Features:
- Automatic observable extraction during simulation
- Real-time temperature monitoring and logging
- Integration with existing analysis workflows
- Support for multiple observables (dipole, vibrational modes)
- CPU/GPU agnostic operation
"""

from typing import Optional, Dict, Any, List, Callable, Union
import numpy as np
import logging
import hoomd
from pathlib import Path
import time

# CuPy import with fallback for CPU/GPU agnostic code
try:
    import cupy as cp
    HAS_CUPY = True
except ImportError:
    cp = None
    HAS_CUPY = False

from .fdr_temperature import (
    FDRTemperatureEstimator, 
    FDRDiagnostics,
    create_dipole_extractor,
    create_mode_extractor
)
from .analysis import TemperatureTracker
from .utils import PhysicalConstants


class FDRTemperatureMonitor:
    """
    Real-time FDR temperature monitoring integrated with HOOMD simulations.
    
    This class acts as a bridge between the FDR temperature estimator and
    the HOOMD simulation, automatically extracting observables and computing
    effective temperatures during the simulation run.
    
    Examples
    --------
    Basic usage with dipole moment observable:
    
    >>> monitor = FDRTemperatureMonitor(
    ...     omega_0=2000.0,  # cm^-1
    ...     dt=0.001,        # ps
    ...     observable_type='dipole',
    ...     axis='z'
    ... )
    >>> monitor.attach_to_simulation(sim)
    >>> monitor.calibrate_equilibrium(T_cal=300.0, n_steps=10000)
    >>> # Simulation runs with automatic FDR monitoring
    
    With vibrational mode projection:
    
    >>> mode_vector = np.array([...])  # Normal mode eigenvector
    >>> masses = np.array([...])       # Particle masses
    >>> monitor = FDRTemperatureMonitor(
    ...     omega_0=mode_frequency,
    ...     dt=sim.dt,
    ...     observable_type='mode',
    ...     mode_vector=mode_vector,
    ...     masses=masses
    ... )
    """
    
    def __init__(
        self,
        omega_0: float,
        dt: float,
        observable_type: str = 'dipole',
        axis: str = 'z',
        mode_vector: Optional[np.ndarray] = None,
        masses: Optional[np.ndarray] = None,
        tau_avg: Optional[float] = None,
        tau_id: Optional[float] = None,
        output_file: Optional[str] = None,
        log_interval: int = 1000
    ):
        """
        Initialize FDR temperature monitor.
        
        Parameters
        ----------
        omega_0 : float
            Target frequency for temperature measurement. If > 100, assumed to be
            in cm^-1 and converted to rad/ps. Otherwise assumed to be rad/ps.
        dt : float
            Simulation timestep (ps)
        observable_type : str
            Type of observable ('dipole' or 'mode')
        axis : str
            Spatial axis for dipole projection ('x', 'y', 'z')
        mode_vector : array_like, optional
            Mass-weighted normal mode eigenvector (required if observable_type='mode')
        masses : array_like, optional  
            Particle masses (required if observable_type='mode')
        tau_avg : float, optional
            Lock-in averaging time constant (ps)
        tau_id : float, optional
            Parameter identification time constant (ps)
        output_file : str, optional
            Output file for FDR temperature trajectory
        log_interval : int
            Interval for logging diagnostics (steps)
        """
        
        # Convert frequency if needed
        if omega_0 > 100:  # Assume cm^-1
            self.omega_0_rad_ps = omega_0 / PhysicalConstants.HARTREE_TO_CM_MINUS1 * 2 * np.pi
            self.omega_0_cm = omega_0
        else:  # Assume rad/ps
            self.omega_0_rad_ps = omega_0
            self.omega_0_cm = omega_0 * PhysicalConstants.HARTREE_TO_CM_MINUS1 / (2 * np.pi)
            
        self.dt = dt
        self.observable_type = observable_type
        self.axis = axis
        self.log_interval = log_interval
        self.output_file = output_file
        
        # Create observable extractor
        if observable_type == 'dipole':
            self.observable_extractor = create_dipole_extractor(axis)
        elif observable_type == 'mode':
            if mode_vector is None or masses is None:
                raise ValueError("mode_vector and masses required for mode observable")
            self.observable_extractor = create_mode_extractor(mode_vector, masses)
        else:
            raise ValueError(f"Unknown observable type: {observable_type}")
            
        # Create FDR estimator
        self.fdr_estimator = FDRTemperatureEstimator(
            omega_0=self.omega_0_rad_ps,
            dt=dt,
            tau_avg=tau_avg,
            tau_id=tau_id,
            observable_extractor=self.observable_extractor
        )
        
        # State tracking
        self.simulation = None
        self.step_count = 0
        self.is_calibrated = False
        self.T_eff_history = []
        self.diagnostics_history = []
        
        # Logging
        self.logger = logging.getLogger(__name__)
        self._setup_output_file()
        
    def _setup_output_file(self):
        """Setup output file for FDR temperature trajectory."""
        
        if self.output_file is not None:
            self.output_path = Path(self.output_file)
            self.output_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Write header
            with open(self.output_path, 'w') as f:
                f.write("# FDR Effective Temperature Trajectory\n")
                f.write(f"# Target frequency: {self.omega_0_cm:.1f} cm^-1 ({self.omega_0_rad_ps:.6f} rad/ps)\n")
                f.write(f"# Observable: {self.observable_type}")
                if self.observable_type == 'dipole':
                    f.write(f" (axis: {self.axis})")
                f.write("\n")
                f.write("# Columns: step, time_ps, T_eff_K, omega_n_rad_ps, gamma_inv_ps, ")
                f.write("sigma_delta_rad_ps, S_AA, chi_imag, snr, lineshape_type\n")
                
    def attach_to_simulation(self, simulation: hoomd.Simulation) -> None:
        """
        Attach the FDR monitor to a HOOMD simulation.
        
        Parameters
        ----------
        simulation : hoomd.Simulation
            HOOMD simulation object
        """
        
        self.simulation = simulation
        
        # Validate simulation setup
        if not hasattr(simulation.state, 'thermodynamic_quantities'):
            self.logger.warning("Simulation state missing thermodynamic quantities")
            
        self.logger.info(f"FDR temperature monitor attached to simulation")
        self.logger.info(f"Target frequency: {self.omega_0_cm:.1f} cm^-1")
        self.logger.info(f"Observable: {self.observable_type}")
        
    def calibrate_equilibrium(
        self, 
        T_cal: float, 
        n_steps: int = 10000,
        force_recalibration: bool = False
    ) -> None:
        """
        Calibrate the FDR estimator using equilibrium simulation data.
        
        Parameters
        ----------
        T_cal : float
            Known equilibrium temperature (K)
        n_steps : int
            Number of equilibration steps for calibration
        force_recalibration : bool
            Force recalibration even if already calibrated
        """
        
        if self.is_calibrated and not force_recalibration:
            self.logger.info("FDR estimator already calibrated")
            return
            
        if self.simulation is None:
            raise RuntimeError("Must attach to simulation before calibration")
            
        self.logger.info(f"Calibrating FDR estimator at T = {T_cal:.1f} K")
        self.logger.info(f"Running {n_steps} equilibration steps...")
        
        # Collect equilibrium data
        equilibrium_data = []
        initial_timestep = self.simulation.timestep
        
        for i in range(n_steps):
            # Extract observable
            snapshot = self.simulation.state.get_snapshot()
            A_val = self.observable_extractor(snapshot)
            equilibrium_data.append(A_val)
            
            # Run one step
            self.simulation.run(1)
            
            # Progress reporting
            if (i + 1) % (n_steps // 10) == 0:
                progress = (i + 1) / n_steps * 100
                self.logger.info(f"Calibration progress: {progress:.1f}%")
                
        # Perform calibration
        self.fdr_estimator.calibrate(T_cal, np.array(equilibrium_data))
        self.is_calibrated = True
        
        final_timestep = self.simulation.timestep
        self.logger.info(f"Calibration complete (steps {initial_timestep} to {final_timestep})")
        
    def update(self) -> tuple[float, FDRDiagnostics]:
        """
        Update FDR temperature estimate for current simulation state.
        
        Returns
        -------
        T_eff : float
            Current effective temperature estimate (K)
        diagnostics : FDRDiagnostics
            Detailed diagnostic information
        """
        
        if self.simulation is None:
            raise RuntimeError("Must attach to simulation before updating")
            
        # Extract current observable
        snapshot = self.simulation.state.get_snapshot()
        A_val = self.observable_extractor(snapshot)
        
        # Update FDR estimator
        current_time = self.simulation.timestep * self.dt
        T_eff, diagnostics = self.fdr_estimator.update(A_val, t=current_time)
        
        self.step_count += 1
        
        # Store history
        self.T_eff_history.append(T_eff)
        self.diagnostics_history.append(diagnostics)
        
        # Log diagnostics periodically
        if self.step_count % self.log_interval == 0:
            self._log_diagnostics(T_eff, diagnostics)
            
        # Write to output file
        if self.output_file is not None:
            self._write_output_line(T_eff, diagnostics)
            
        return T_eff, diagnostics
        
    def _log_diagnostics(self, T_eff: float, diagnostics: FDRDiagnostics) -> None:
        """Log diagnostic information."""
        
        timestep = self.simulation.timestep
        time_ps = timestep * self.dt
        
        self.logger.info(f"Step {timestep:8d} | Time {time_ps:8.3f} ps | T_eff {T_eff:6.1f} K")
        self.logger.debug(f"  ωn = {diagnostics.omega_n:.6f} rad/ps, γ = {diagnostics.gamma:.6f} 1/ps")
        self.logger.debug(f"  S_AA = {diagnostics.S_AA:.6e}, χ'' = {diagnostics.chi_imag:.6e}")
        self.logger.debug(f"  SNR = {diagnostics.snr:.3f}, lineshape = {diagnostics.lineshape_type.value}")
        
    def _write_output_line(self, T_eff: float, diagnostics: FDRDiagnostics) -> None:
        """Write current data to output file."""
        
        if self.output_file is None:
            return
            
        timestep = self.simulation.timestep
        time_ps = timestep * self.dt
        
        with open(self.output_path, 'a') as f:
            f.write(f"{timestep:8d} {time_ps:12.6f} {T_eff:8.2f} ")
            f.write(f"{diagnostics.omega_n:12.6f} {diagnostics.gamma:12.6f} ")
            f.write(f"{diagnostics.sigma_delta:12.6f} {diagnostics.S_AA:12.6e} ")
            f.write(f"{diagnostics.chi_imag:12.6e} {diagnostics.snr:8.3f} ")
            f.write(f"{diagnostics.lineshape_type.value}\n")
            
    def get_recent_temperature(self, n_points: int = 100) -> float:
        """
        Get recent average effective temperature.
        
        Parameters
        ----------
        n_points : int
            Number of recent points to average
            
        Returns
        -------
        float
            Recent average effective temperature (K)
        """
        
        if len(self.T_eff_history) == 0:
            return np.nan
            
        recent_temps = self.T_eff_history[-n_points:]
        valid_temps = [T for T in recent_temps if not np.isnan(T)]
        
        if len(valid_temps) == 0:
            return np.nan
            
        return np.mean(valid_temps)
        
    def get_temperature_statistics(self) -> Dict[str, float]:
        """
        Get statistical summary of temperature trajectory.
        
        Returns
        -------
        dict
            Statistics including mean, std, min, max, trend
        """
        
        if len(self.T_eff_history) == 0:
            return {}
            
        valid_temps = [T for T in self.T_eff_history if not np.isnan(T)]
        
        if len(valid_temps) == 0:
            return {}
            
        temps = np.array(valid_temps)
        
        # Basic statistics
        stats = {
            'mean': np.mean(temps),
            'std': np.std(temps),
            'min': np.min(temps),
            'max': np.max(temps),
            'n_points': len(temps)
        }
        
        # Trend analysis (simple linear fit)
        if len(temps) > 10:
            x = np.arange(len(temps))
            slope, intercept = np.polyfit(x, temps, 1)
            stats['trend_slope'] = slope
            stats['trend_intercept'] = intercept
            
        return stats
        
    def reset(self) -> None:
        """Reset the FDR estimator state."""
        
        self.step_count = 0
        self.T_eff_history.clear()
        self.diagnostics_history.clear()
        self.fdr_estimator._reset_state()
        self.is_calibrated = False
        
        self.logger.info("FDR temperature monitor reset")


class CavityFDRAnalyzer:
    """
    Advanced FDR analysis for cavity-coupled molecular dynamics.
    
    This class provides comprehensive analysis tools for understanding
    effective temperature dynamics in cavity QED systems, including
    mode-specific thermalization and non-equilibrium dynamics.
    """
    
    def __init__(self, fdr_monitor: FDRTemperatureMonitor):
        """
        Initialize cavity FDR analyzer.
        
        Parameters
        ----------
        fdr_monitor : FDRTemperatureMonitor
            FDR temperature monitor instance
        """
        
        self.monitor = fdr_monitor
        self.logger = logging.getLogger(__name__)
        
    def analyze_thermalization_dynamics(
        self, 
        window_size: int = 1000
    ) -> Dict[str, Any]:
        """
        Analyze thermalization dynamics from FDR temperature trajectory.
        
        Parameters
        ----------
        window_size : int
            Window size for moving averages
            
        Returns
        -------
        dict
            Analysis results including relaxation times, equilibrium values
        """
        
        if len(self.monitor.T_eff_history) < window_size:
            self.logger.warning("Insufficient data for thermalization analysis")
            return {}
            
        temps = np.array([T for T in self.monitor.T_eff_history if not np.isnan(T)])
        times = np.arange(len(temps)) * self.monitor.dt
        
        # Moving average
        smoothed_temps = np.convolve(temps, np.ones(window_size)/window_size, mode='valid')
        smoothed_times = times[window_size//2:len(smoothed_temps)+window_size//2]
        
        # Fit exponential relaxation
        try:
            from scipy.optimize import curve_fit
            
            def exponential_relaxation(t, T_eq, T_initial, tau):
                return T_eq + (T_initial - T_eq) * np.exp(-t / tau)
                
            # Initial guess
            T_eq_guess = smoothed_temps[-window_size//4:].mean()
            T_initial_guess = smoothed_temps[:window_size//4].mean()
            tau_guess = smoothed_times[-1] / 3
            
            popt, pcov = curve_fit(
                exponential_relaxation,
                smoothed_times,
                smoothed_temps,
                p0=[T_eq_guess, T_initial_guess, tau_guess],
                maxfev=1000
            )
            
            T_eq, T_initial, tau = popt
            param_errors = np.sqrt(np.diag(pcov))
            
            results = {
                'equilibrium_temperature': T_eq,
                'initial_temperature': T_initial,
                'relaxation_time': tau,
                'temperature_errors': param_errors,
                'fit_quality': np.corrcoef(
                    smoothed_temps,
                    exponential_relaxation(smoothed_times, *popt)
                )[0, 1]**2
            }
            
        except Exception as e:
            self.logger.warning(f"Exponential fit failed: {e}")
            results = {}
            
        return results
        
    def compare_with_equipartition(self, bath_temperature: float) -> Dict[str, float]:
        """
        Compare FDR temperature with equipartition temperature.
        
        Parameters
        ----------
        bath_temperature : float
            Bath temperature for comparison (K)
            
        Returns
        -------
        dict
            Comparison metrics
        """
        
        recent_T_eff = self.monitor.get_recent_temperature()
        
        if np.isnan(recent_T_eff):
            return {}
            
        deviation = recent_T_eff - bath_temperature
        relative_deviation = deviation / bath_temperature
        
        return {
            'fdr_temperature': recent_T_eff,
            'bath_temperature': bath_temperature,
            'absolute_deviation': deviation,
            'relative_deviation': relative_deviation,
            'is_non_equilibrium': abs(relative_deviation) > 0.05  # 5% threshold
        }


def create_cavity_fdr_monitor(
    cavity_frequency_cm: float,
    simulation_dt: float,
    observable_type: str = 'dipole',
    **kwargs
) -> FDRTemperatureMonitor:
    """
    Convenience function to create FDR monitor for cavity simulations.
    
    Parameters
    ----------
    cavity_frequency_cm : float
        Cavity frequency in cm^-1
    simulation_dt : float
        Simulation timestep in ps
    observable_type : str
        Type of observable ('dipole' or 'mode')
    **kwargs
        Additional arguments for FDRTemperatureMonitor
        
    Returns
    -------
    FDRTemperatureMonitor
        Configured FDR monitor
    """
    
    return FDRTemperatureMonitor(
        omega_0=cavity_frequency_cm,
        dt=simulation_dt,
        observable_type=observable_type,
        **kwargs
    )
