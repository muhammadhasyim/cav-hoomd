"""
Fluctuation-Dissipation Ratio (FDR) Based Effective Temperature Estimator

This module implements a streaming algorithm for measuring single-frequency 
effective temperature T_eff(ω₀,t) from molecular dynamics trajectories using 
violations of fluctuation-dissipation relations.

The core principle is:
    T_eff(ω₀,t) = (ω₀/2k_B) × S_AA(ω₀,t) / χ''(ω₀,t)

Where:
- S_AA(ω₀,t): Fluctuation power spectrum at frequency ω₀
- χ''(ω₀,t): Imaginary susceptibility at frequency ω₀

Scientific References:
- Fluctuation-Dissipation Theorem and harmonic oscillator theory
- Non-equilibrium statistical mechanics
- Cavity quantum electrodynamics and polariton physics

Author: Computational Chemistry Team
Date: 2025
"""

from typing import Optional, Union, Tuple, Dict, Any, Callable
import numpy as np
import logging
from dataclasses import dataclass
from enum import Enum
import warnings

# CuPy import with fallback for CPU/GPU agnostic code
try:
    import cupy as cp
    HAS_CUPY = True
except ImportError:
    cp = None
    HAS_CUPY = False

from .utils import PhysicalConstants


class LineshapeType(Enum):
    """Lineshape models for susceptibility calculation."""
    LORENTZIAN = "lorentzian"
    VOIGT = "voigt"
    MOTIONAL_NARROWING = "motional_narrowing"


@dataclass
class FDRDiagnostics:
    """Diagnostic information from FDR temperature estimation."""
    
    # Mechanical parameters
    omega_n: float = 0.0        # Natural frequency (rad/ps)
    gamma: float = 0.0          # Damping coefficient (1/ps)
    gamma_eff: float = 0.0      # Effective damping including dephasing
    
    # Dephasing statistics
    sigma_delta: float = 0.0    # RMS frequency jitter (rad/ps)
    tau_c: float = 0.0          # Correlation time of frequency noise (ps)
    D_phi: float = 0.0          # Phase diffusion coefficient
    
    # Spectral quantities
    S_AA: float = 0.0           # Fluctuation power at ω₀
    chi_imag: float = 0.0       # Imaginary susceptibility at ω₀
    
    # Quality metrics
    snr: float = 0.0            # Signal-to-noise ratio
    residual_rms: float = 0.0   # RMS residual from operator identification
    lineshape_type: LineshapeType = LineshapeType.LORENTZIAN
    
    # Uncertainty estimates
    T_eff_uncertainty: float = 0.0
    effective_dof: float = 0.0  # Effective degrees of freedom for statistics


class FDRTemperatureEstimator:
    """
    Streaming FDR-based effective temperature estimator.
    
    This class implements the complete algorithm from the pedagogical notes,
    including three streaming estimators:
    1. Lock-in amplifier for S_AA(ω₀,t)
    2. Operator identification for (ωn, γ)
    3. Baseband phase tracking for dephasing statistics
    
    The estimator is designed to be O(1) per timestep with exponential forgetting,
    requiring no large buffers and maintaining constant memory usage.
    """
    
    def __init__(
        self,
        omega_0: float,
        dt: float,
        tau_avg: Optional[float] = None,
        tau_id: Optional[float] = None,
        observable_extractor: Optional[Callable] = None,
        force_extractor: Optional[Callable] = None,
        ridge_regularization: float = 1e-8,
        baseband_cutoff_factor: float = 10.0
    ):
        """
        Initialize FDR temperature estimator.
        
        Parameters
        ----------
        omega_0 : float
            Target frequency for temperature measurement (rad/ps)
        dt : float  
            MD timestep (ps)
        tau_avg : float, optional
            Averaging time constant for lock-in (ps). Default: 100/γ estimated
        tau_id : float, optional
            Parameter identification time constant (ps). Default: 500 periods
        observable_extractor : callable, optional
            Function to extract observable A(t) from simulation state
        force_extractor : callable, optional
            Function to extract microscopic force F_A^(mic)(t)
        ridge_regularization : float
            Ridge parameter for stable regression (default: 1e-8)
        baseband_cutoff_factor : float
            Baseband filter cutoff as factor × mechanical linewidth
        """
        
        # Validate inputs
        if omega_0 <= 0:
            raise ValueError("Target frequency omega_0 must be positive")
        if dt <= 0:
            raise ValueError("Timestep dt must be positive")
        if dt > 2*np.pi/omega_0/30:
            warnings.warn(f"Timestep dt={dt:.6f} ps may be too large for frequency "
                         f"omega_0={omega_0:.3f} rad/ps. Recommend dt ≤ {2*np.pi/omega_0/30:.6f} ps")
        
        # Core parameters
        self.omega_0 = omega_0
        self.dt = dt
        self.T_0 = 2*np.pi/omega_0  # Period
        
        # Time constants (will be finalized after initial parameter estimates)
        self._tau_avg_target = tau_avg
        self._tau_id_target = tau_id
        
        # Function extractors
        self.observable_extractor = observable_extractor
        self.force_extractor = force_extractor
        
        # Numerical parameters
        self.ridge_regularization = ridge_regularization
        self.baseband_cutoff_factor = baseband_cutoff_factor
        
        # State variables for streaming estimators
        self._reset_state()
        
        # Calibration
        self._G_gain = None  # Will be set during calibration
        self._is_calibrated = False
        
        # Logging
        self.logger = logging.getLogger(__name__)
        
    def _reset_state(self):
        """Reset all internal state variables."""
        
        # Lock-in amplifier state for S_AA estimation
        self.I = 0.0  # In-phase component
        self.Q = 0.0  # Quadrature component  
        self.W = 0.0  # Normalization weight
        
        # Parameter identification state
        self.S_matrix = np.zeros((2, 2))  # Normal matrix for [k_eff, gamma]
        self.b_vector = np.zeros(2)       # RHS vector
        self.param_history = []           # For adaptive time constants
        
        # Current parameter estimates
        self.omega_n = self.omega_0  # Start near target frequency
        self.gamma = 0.1             # Initial damping guess
        
        # Baseband phase tracking
        self.u_complex = 0.0 + 0.0j   # Complex envelope
        self.phi_prev = 0.0           # Previous phase
        self.phi_unwrapped = 0.0      # Unwrapped phase
        self.phase_diff_sq_avg = 0.0  # For phase diffusion
        self.sigma_delta = 0.0        # Frequency jitter RMS
        
        # Dephasing time constant (rough estimate)
        self.tau_c = 10 * self.T_0    # Start with 10 periods
        
        # Adaptive time constants
        self._update_time_constants()
        
        # Step counter
        self.step_count = 0
        
        # Savitzky-Golay differentiation buffer (causal, short window)
        self._sg_buffer = np.zeros(7)  # 7-point window
        self._sg_coeffs = np.array([-3, -2, -1, 0, 1, 2, 3]) / 28.0  # SG coeffs for order 1
        
    def _update_time_constants(self):
        """Update forgetting factors based on current parameter estimates."""
        
        # Lock-in averaging time constant
        if self._tau_avg_target is None:
            self.tau_avg = max(50 * self.T_0, 100.0 / max(self.gamma, 0.01))
        else:
            self.tau_avg = self._tau_avg_target
            
        # Parameter identification time constant  
        if self._tau_id_target is None:
            self.tau_id = max(200 * self.T_0, 500.0 / max(self.gamma, 0.01))
        else:
            self.tau_id = self._tau_id_target
            
        # Convert to forgetting factors
        self.lambda_avg = np.exp(-self.dt / self.tau_avg)
        self.lambda_id = np.exp(-self.dt / self.tau_id)
        
        # Ensure reasonable forgetting rates (prevent too slow averaging)
        self.lambda_avg = min(self.lambda_avg, 0.999)  # At least 0.1% forgetting per step
        self.lambda_id = min(self.lambda_id, 0.9999)   # Slower for parameter ID
        
        # Baseband filter parameters (single-pole IIR)
        gamma_est = max(self.gamma, 0.001)  # Avoid division by zero
        baseband_cutoff = gamma_est / self.baseband_cutoff_factor
        self.lambda_baseband = np.exp(-baseband_cutoff * self.dt)
        
    def calibrate(self, T_cal: float, A_data: Optional[np.ndarray] = None) -> None:
        """
        Perform one-time gain calibration using equilibrium data.
        
        Parameters
        ----------
        T_cal : float
            Known calibration temperature (K)
        A_data : array_like, optional
            Equilibrium observable data for calibration. If None, uses current
            window averages.
        """
        if T_cal <= 0:
            raise ValueError("Calibration temperature must be positive")
            
        if A_data is not None:
            # Process calibration data
            for A_val in A_data:
                self._update_lock_in(A_val)
                
        # Get current spectral estimates
        S_AA_cal = self._get_S_AA()
        chi_imag_cal = self._get_chi_imaginary_shape()  # Without gain factor
        
        if S_AA_cal <= 0 or chi_imag_cal <= 0:
            raise ValueError("Invalid spectral estimates for calibration")
            
        # Compute gain factor
        k_B = PhysicalConstants.KB_HARTREE_PER_K  # Boltzmann constant in Hartree/K
        self._G_gain = (self.omega_0 / (2 * k_B * T_cal)) * (S_AA_cal / chi_imag_cal)
        
        self._is_calibrated = True
        self.logger.info(f"FDR estimator calibrated with gain G = {self._G_gain:.6e}")
        
    def update(
        self, 
        A: float, 
        forces: Optional[np.ndarray] = None, 
        t: Optional[float] = None
    ) -> Tuple[float, FDRDiagnostics]:
        """
        Single-step update of the FDR temperature estimator.
        
        Parameters
        ----------
        A : float
            Current value of observable
        forces : array_like, optional
            Microscopic forces for force-based parameter identification
        t : float, optional
            Current time (for reference, not used in computation)
            
        Returns
        -------
        T_eff : float
            Current effective temperature estimate (K)
        diagnostics : FDRDiagnostics
            Detailed diagnostic information
        """
        
        self.step_count += 1
        current_time = self.step_count * self.dt
        
        # 1. Update lock-in amplifier for S_AA
        self._update_lock_in(A)
        
        # 2. Update parameter identification 
        if forces is not None and self.force_extractor is not None:
            F_A_mic = self.force_extractor(forces)
            self._update_force_regression(A, F_A_mic)
        else:
            self._update_ar2_identification(A)
            
        # 3. Update baseband phase tracking
        self._update_baseband_phase(A)
        
        # 4. Adaptive time constant updates (every 10 steps)
        if self.step_count % 10 == 0:
            self._update_time_constants()
            
        # 5. Compute susceptibility
        chi_imag = self._get_chi_imaginary()
        
        # 6. Compute effective temperature
        S_AA = self._get_S_AA()
        T_eff = self._compute_T_eff(S_AA, chi_imag)
        
        # 7. Prepare diagnostics
        diagnostics = self._get_diagnostics(S_AA, chi_imag, T_eff)
        
        return T_eff, diagnostics
        
    def _update_lock_in(self, A: float) -> None:
        """Update lock-in amplifier for S_AA estimation."""
        
        phase = self.omega_0 * self.step_count * self.dt
        
        # Update I/Q components with exponential forgetting
        cos_term = A * np.cos(phase)
        sin_term = A * np.sin(phase)
        
        self.I = self.lambda_avg * self.I + (1 - self.lambda_avg) * cos_term
        self.Q = self.lambda_avg * self.Q + (1 - self.lambda_avg) * sin_term
        self.W = self.lambda_avg * self.W + (1 - self.lambda_avg)
        
    def _update_force_regression(self, A: float, F_A_mic: float) -> None:
        """Update force-based parameter identification."""
        
        # Causal Savitzky-Golay differentiation
        self._sg_buffer[:-1] = self._sg_buffer[1:]  # Shift buffer
        self._sg_buffer[-1] = A
        
        if self.step_count < 7:
            A_dot = 0.0  # Not enough history
        else:
            A_dot = np.sum(self._sg_coeffs * self._sg_buffer) / self.dt
            
        # Regression: -F_A_mic = k_eff * A + gamma * A_dot
        phi = np.array([A, A_dot])
        y = -F_A_mic
        
        # Update normal equations with exponential forgetting
        self.S_matrix = self.lambda_id * self.S_matrix + np.outer(phi, phi)
        self.b_vector = self.lambda_id * self.b_vector + phi * y
        
        # Solve with ridge regularization
        S_reg = self.S_matrix + self.ridge_regularization * np.eye(2)
        try:
            theta = np.linalg.solve(S_reg, self.b_vector)
            k_eff, gamma = theta
            
            # Update parameters
            self.gamma = max(gamma, 0.001)  # Ensure positive damping
            self.omega_n = np.sqrt(max(k_eff, 0.0))
            
        except np.linalg.LinAlgError:
            # Keep previous estimates if solve fails
            pass
            
    def _update_ar2_identification(self, A: float) -> None:
        """Update AR(2) parameter identification (kinematic method)."""
        
        # Store recent A values for AR(2) fitting
        if not hasattr(self, '_ar2_buffer'):
            self._ar2_buffer = np.zeros(3)
            self._ar2_P = np.eye(2) * 1000.0  # Initial covariance
            self._ar2_theta = np.array([0.0, 0.0])  # [a1, a2]
            
        # Shift buffer
        self._ar2_buffer[:-1] = self._ar2_buffer[1:]
        self._ar2_buffer[-1] = A
        
        if self.step_count < 3:
            return  # Need at least 3 points
            
        # RLS update for A_n + a1*A_{n-1} + a2*A_{n-2} = epsilon_n
        x = np.array([self._ar2_buffer[1], self._ar2_buffer[0]])  # [A_{n-1}, A_{n-2}]
        y = self._ar2_buffer[2]  # A_n
        
        # RLS update
        k = self._ar2_P @ x / (self.lambda_id + x @ self._ar2_P @ x)
        self._ar2_theta = self._ar2_theta + k * (y + x @ self._ar2_theta)
        self._ar2_P = (self._ar2_P - np.outer(k, x @ self._ar2_P)) / self.lambda_id
        
        # Map AR(2) coefficients to physical parameters
        a1, a2 = self._ar2_theta
        
        # Characteristic equation: z^2 + a1*z + a2 = 0
        discriminant = a1**2 - 4*a2
        
        if discriminant < 0 and a2 > 0:  # Complex roots (oscillatory)
            r = np.sqrt(a2)
            if r < 1.0:  # Stable
                theta = np.arccos(-a1 / (2*r))
                
                # Map to physical parameters
                Omega_d = theta / self.dt
                gamma_mapped = (2.0 / self.dt) * np.log(1.0 / r)
                omega_n_mapped = np.sqrt(Omega_d**2 + (gamma_mapped/2)**2)
                
                # Update with validation
                if 0.001 < gamma_mapped < 10.0 and 0.1 < omega_n_mapped < 100.0:
                    self.gamma = gamma_mapped
                    self.omega_n = omega_n_mapped
                    
    def _update_baseband_phase(self, A: float) -> None:
        """Update baseband complex envelope and phase tracking."""
        
        # Demodulate to baseband
        phase = self.omega_0 * self.step_count * self.dt
        u_raw = A * np.exp(-1j * phase)
        
        # Low-pass filter (single-pole IIR)
        self.u_complex = self.lambda_baseband * self.u_complex + (1 - self.lambda_baseband) * u_raw
        
        # Extract phase and unwrap
        if abs(self.u_complex) > 1e-10:
            phi_current = np.angle(self.u_complex)
            
            # Unwrap phase
            phi_diff = phi_current - self.phi_prev
            phi_diff = ((phi_diff + np.pi) % (2*np.pi)) - np.pi  # Wrap to [-π, π]
            self.phi_unwrapped += phi_diff
            self.phi_prev = phi_current
            
            # Update phase diffusion statistics
            if self.step_count > 1:
                phase_diff_rate = phi_diff / self.dt
                self.phase_diff_sq_avg = (self.lambda_avg * self.phase_diff_sq_avg + 
                                        (1 - self.lambda_avg) * phase_diff_rate**2)
                
                # Phase diffusion coefficient
                D_phi = self.phase_diff_sq_avg / (2 * self.dt)
                
                # RMS frequency jitter
                self.sigma_delta = np.sqrt(self.phase_diff_sq_avg)
                
    def _get_S_AA(self) -> float:
        """Compute current S_AA estimate from lock-in."""
        if self.W > 0:
            return 2.0 * (self.I**2 + self.Q**2)  # Power spectral density
        return 0.0
        
    def _get_chi_imaginary_shape(self) -> float:
        """Compute susceptibility lineshape (without gain factor G)."""
        
        # Determine lineshape type
        omega_tau_c = self.omega_0 * self.tau_c
        
        if omega_tau_c < 1.0:  # Motional narrowing regime
            gamma_eff = self.gamma + 2 * self.sigma_delta**2 * self.tau_c
            return self._lorentzian_susceptibility(gamma_eff, gain=1.0)
        else:  # Quasi-static regime - use Voigt
            return self._voigt_susceptibility(gain=1.0)
            
    def _get_chi_imaginary(self) -> float:
        """Compute susceptibility with gain factor."""
        if not self._is_calibrated:
            warnings.warn("FDR estimator not calibrated - using unit gain")
            gain = 1.0
        else:
            gain = self._G_gain
            
        return gain * self._get_chi_imaginary_shape()
        
    def _lorentzian_susceptibility(self, gamma_eff: float, gain: float) -> float:
        """Compute Lorentzian susceptibility at omega_0."""
        
        numerator = gain * gamma_eff * self.omega_0
        delta_omega_sq = (self.omega_n**2 - self.omega_0**2)**2
        gamma_term_sq = (gamma_eff * self.omega_0)**2
        
        return numerator / (delta_omega_sq + gamma_term_sq)
        
    def _voigt_susceptibility(self, gain: float) -> float:
        """Compute Voigt susceptibility (approximation for single frequency)."""
        
        # Simplified Voigt approximation - full implementation would use Faddeeva function
        # For now, use Lorentzian with effective broadening
        sigma_omega = self.sigma_delta
        gamma_voigt = self.gamma + 0.5 * sigma_omega  # Rough approximation
        
        return self._lorentzian_susceptibility(gamma_voigt, gain)
        
    def _compute_T_eff(self, S_AA: float, chi_imag: float) -> float:
        """Compute effective temperature from S_AA and chi_imag."""
        
        if not self._is_calibrated:
            return np.nan
            
        if chi_imag <= 0 or S_AA <= 0:
            return np.nan
            
        k_B = PhysicalConstants.KB_HARTREE_PER_K
        T_eff = (self.omega_0 / (2 * k_B)) * (S_AA / chi_imag)
        
        return T_eff
        
    def _get_diagnostics(self, S_AA: float, chi_imag: float, T_eff: float) -> FDRDiagnostics:
        """Prepare comprehensive diagnostics."""
        
        # Determine effective damping and lineshape type
        omega_tau_c = self.omega_0 * self.tau_c
        if omega_tau_c < 1.0:
            gamma_eff = self.gamma + 2 * self.sigma_delta**2 * self.tau_c
            lineshape = LineshapeType.MOTIONAL_NARROWING
        else:
            gamma_eff = self.gamma + 0.5 * self.sigma_delta  # Voigt approximation
            lineshape = LineshapeType.VOIGT
            
        # SNR estimate
        snr = S_AA / max(1e-10, np.sqrt(S_AA / max(self.W, 1.0)))
        
        # Effective degrees of freedom for uncertainty
        bandwidth = 1.0 / (2 * np.pi * self.tau_avg)
        T_obs = self.step_count * self.dt
        eff_dof = 2 * bandwidth * T_obs
        
        # Temperature uncertainty (rough estimate)
        T_eff_uncertainty = T_eff / np.sqrt(max(eff_dof, 1.0))
        
        return FDRDiagnostics(
            omega_n=self.omega_n,
            gamma=self.gamma,
            gamma_eff=gamma_eff,
            sigma_delta=self.sigma_delta,
            tau_c=self.tau_c,
            D_phi=self.phase_diff_sq_avg / (2 * self.dt),
            S_AA=S_AA,
            chi_imag=chi_imag,
            snr=snr,
            residual_rms=0.0,  # TODO: implement residual tracking
            lineshape_type=lineshape,
            T_eff_uncertainty=T_eff_uncertainty,
            effective_dof=eff_dof
        )
        
    def get_state(self) -> Dict[str, Any]:
        """Get complete internal state for checkpointing."""
        return {
            'step_count': self.step_count,
            'I': self.I, 'Q': self.Q, 'W': self.W,
            'S_matrix': self.S_matrix, 'b_vector': self.b_vector,
            'omega_n': self.omega_n, 'gamma': self.gamma,
            'u_complex': self.u_complex, 'phi_prev': self.phi_prev,
            'phi_unwrapped': self.phi_unwrapped,
            'phase_diff_sq_avg': self.phase_diff_sq_avg,
            'sigma_delta': self.sigma_delta, 'tau_c': self.tau_c,
            'G_gain': self._G_gain, 'is_calibrated': self._is_calibrated
        }
        
    def set_state(self, state: Dict[str, Any]) -> None:
        """Restore complete internal state from checkpoint."""
        for key, value in state.items():
            if hasattr(self, key):
                setattr(self, key, value)
            elif key == 'G_gain':
                self._G_gain = value
            elif key == 'is_calibrated':
                self._is_calibrated = value


def create_dipole_extractor(axis: str = 'z') -> Callable:
    """
    Create observable extractor for total dipole moment projection.
    
    Parameters
    ----------
    axis : str
        Spatial axis for dipole projection ('x', 'y', or 'z')
        
    Returns
    -------
    extractor : callable
        Function that extracts A(t) = ê·M(t) from simulation state
    """
    
    axis_index = {'x': 0, 'y': 1, 'z': 2}[axis.lower()]
    
    def dipole_extractor(snapshot):
        """Extract dipole moment projection from HOOMD snapshot."""
        
        # Get array module (NumPy or CuPy)
        if hasattr(snapshot.particles.position, 'device'):
            xp = cp.get_array_module(snapshot.particles.position)
        else:
            xp = np
            
        positions = snapshot.particles.position
        charges = snapshot.particles.charge
        
        # Total dipole moment: M = Σᵢ qᵢ rᵢ
        dipole_moment = xp.sum(charges[:, None] * positions, axis=0)
        
        # Project onto specified axis
        A = dipole_moment[axis_index]
        
        # Convert to NumPy if needed
        if xp is cp:
            A = cp.asnumpy(A)
            
        return float(A)
        
    return dipole_extractor


def create_mode_extractor(mode_vector: np.ndarray, masses: np.ndarray) -> Callable:
    """
    Create observable extractor for vibrational mode projection.
    
    Parameters
    ----------
    mode_vector : array_like, shape (N, 3)
        Mass-weighted eigenvector for the target mode
    masses : array_like, shape (N,)
        Particle masses
        
    Returns
    -------
    extractor : callable
        Function that extracts A(t) = Σᵢ √mᵢ eᵢ,k·uᵢ(t)
    """
    
    mode_vector = np.asarray(mode_vector)
    masses = np.asarray(masses)
    
    def mode_extractor(snapshot, reference_positions=None):
        """Extract mode projection from HOOMD snapshot."""
        
        # Get array module
        if hasattr(snapshot.particles.position, 'device'):
            xp = cp.get_array_module(snapshot.particles.position)
        else:
            xp = np
            
        positions = snapshot.particles.position
        
        # Compute displacements from reference
        if reference_positions is not None:
            if xp is cp:
                ref_pos = cp.asarray(reference_positions)
            else:
                ref_pos = reference_positions
            displacements = positions - ref_pos
        else:
            displacements = positions
            
        # Mode projection: A = Σᵢ √mᵢ eᵢ,k·uᵢ
        if xp is cp:
            masses_gpu = cp.asarray(masses)
            mode_gpu = cp.asarray(mode_vector)
            sqrt_masses = cp.sqrt(masses_gpu)
            A = cp.sum(sqrt_masses[:, None] * mode_gpu * displacements)
            A = cp.asnumpy(A)
        else:
            sqrt_masses = np.sqrt(masses)
            A = np.sum(sqrt_masses[:, None] * mode_vector * displacements)
            
        return float(A)
        
    return mode_extractor
