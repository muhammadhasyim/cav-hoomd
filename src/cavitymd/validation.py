"""
Parameter validation framework for cavity molecular dynamics simulations.

This module provides comprehensive validation for all simulation parameters,
following scientific best practices and ensuring reproducible results.
"""

import logging
from dataclasses import dataclass, field
from typing import Optional, Union, List, Dict, Any
from pathlib import Path
import hoomd
import numpy as np
from .utils import PhysicalConstants


@dataclass
class CavitySimulationParams:
    """
    Validated parameters for cavity MD simulations.
    
    This class encapsulates all simulation parameters with comprehensive validation
    to ensure scientific accuracy and reproducibility.
    
    Parameters
    ----------
    temperature : float
        Temperature in Kelvin (must be positive)
    frequency : float
        Cavity mode frequency in cm^-1 (must be positive)
    coupling_strength : float
        Light-matter coupling strength in atomic units
    runtime_ps : float
        Simulation runtime in picoseconds (must be positive)
    molecular_thermostat : str
        Thermostat for molecular particles ('bussi', 'langevin', 'none')
    cavity_thermostat : str
        Thermostat for cavity particle ('bussi', 'langevin', 'none')
    molecular_thermostat_tau : float
        Molecular thermostat time constant in picoseconds (must be positive)
    cavity_thermostat_tau : float
        Cavity thermostat time constant in picoseconds (must be positive)
    
    Raises
    ------
    ValueError
        If any parameter fails validation
    """
    
    # Core simulation parameters
    temperature: float
    frequency: float
    coupling_strength: float
    runtime_ps: float
    
    # Thermostat configuration
    molecular_thermostat: str = 'bussi'
    cavity_thermostat: str = 'langevin'
    molecular_thermostat_tau: float = 5.0
    cavity_thermostat_tau: float = 5.0
    
    # Cavity configuration
    incavity: bool = True
    finite_q: bool = False
    add_cavity_particle: bool = True
    cavity_damping_factor: float = 1.0
    
    # Time-varying parameters
    switch_time_ps: Optional[float] = None
    dissipation: float = 0.0
    
    # Simulation control
    error_tolerance: float = 0.01
    dt_fs: Optional[float] = None
    device: str = 'CPU'
    gpu_id: int = 0
    seed: Optional[int] = None
    restart_velocities: bool = True
    
    # Input/output configuration
    job_dir: str = '.'
    replica: int = 1
    input_gsd: str = 'molecular-0.gsd'
    frame: int = -1
    name: str = 'prod'
    
    # Analysis and output periods (in ps)
    energy_output_period_ps: float = 0.1
    fkt_output_period_ps: float = 1.0
    gsd_output_period_ps: float = 50.0
    console_output_period_ps: float = 1.0
    
    # F(k,t) analysis parameters
    enable_fkt: bool = True
    fkt_kmag: float = 1.0
    fkt_num_wavevectors: int = 50
    fkt_reference_interval_ps: float = 1.0
    fkt_max_references: int = 10
    
    # Energy tracking
    enable_energy_tracking: bool = False
    max_energy_output_time_ps: Optional[float] = None
    
    # Logging configuration
    log_level: str = 'INFO'
    custom_log_file: Optional[str] = None
    enable_text_output: bool = False
    text_output_file: Optional[str] = None
    truncate_gsd: bool = False
    
    # Brownian dynamics (legacy)
    use_brownian_overdamped: bool = True
    
    def __post_init__(self):
        """Validate all parameters after initialization."""
        self.validate()
    
    def validate(self) -> None:
        """
        Comprehensive parameter validation.
        
        Raises
        ------
        ValueError
            If any parameter fails validation
        """
        # Core parameter validation
        self._validate_temperature()
        self._validate_frequency()
        self._validate_coupling_strength()
        self._validate_runtime()
        
        # Thermostat validation
        self._validate_thermostats()
        
        # Cavity configuration validation
        self._validate_cavity_config()
        
        # Time-varying parameters validation
        self._validate_time_varying_params()
        
        # Simulation control validation
        self._validate_simulation_control()
        
        # Output configuration validation
        self._validate_output_config()
        
        # Analysis parameters validation
        self._validate_analysis_params()
        
        # Input/output validation
        self._validate_io_config()
        
        # Log successful validation
        logging.info("Parameter validation completed successfully")
    
    def _validate_temperature(self) -> None:
        """Validate temperature parameter."""
        if not isinstance(self.temperature, (int, float)):
            raise ValueError(f"Temperature must be a number, got {type(self.temperature)}")
        if self.temperature <= 0:
            raise ValueError(f"Temperature must be positive, got {self.temperature} K")
        if self.temperature > 10000:
            logging.warning(f"Very high temperature: {self.temperature} K")
    
    def _validate_frequency(self) -> None:
        """Validate cavity frequency parameter."""
        if not isinstance(self.frequency, (int, float)):
            raise ValueError(f"Frequency must be a number, got {type(self.frequency)}")
        if self.frequency <= 0:
            raise ValueError(f"Frequency must be positive, got {self.frequency} cm^-1")
        if self.frequency > 10000:
            logging.warning(f"Very high frequency: {self.frequency} cm^-1")
    
    def _validate_coupling_strength(self) -> None:
        """Validate coupling strength parameter."""
        if not isinstance(self.coupling_strength, (int, float)):
            raise ValueError(f"Coupling strength must be a number, got {type(self.coupling_strength)}")
        if self.coupling_strength < 0:
            raise ValueError(f"Coupling strength must be non-negative, got {self.coupling_strength}")
        if self.coupling_strength > 1.0:
            logging.warning(f"Very strong coupling: {self.coupling_strength} a.u.")
    
    def _validate_runtime(self) -> None:
        """Validate simulation runtime."""
        if not isinstance(self.runtime_ps, (int, float)):
            raise ValueError(f"Runtime must be a number, got {type(self.runtime_ps)}")
        if self.runtime_ps <= 0:
            raise ValueError(f"Runtime must be positive, got {self.runtime_ps} ps")
        if self.runtime_ps > 1000000:
            logging.warning(f"Very long runtime: {self.runtime_ps} ps")
    
    def _validate_thermostats(self) -> None:
        """Validate thermostat configurations."""
        valid_thermostats = ['bussi', 'langevin', 'none']
        
        if self.molecular_thermostat not in valid_thermostats:
            raise ValueError(f"Invalid molecular thermostat: {self.molecular_thermostat}. "
                           f"Must be one of {valid_thermostats}")
        
        if self.cavity_thermostat not in valid_thermostats:
            raise ValueError(f"Invalid cavity thermostat: {self.cavity_thermostat}. "
                           f"Must be one of {valid_thermostats}")
        
        # Validate thermostat time constants
        if self.molecular_thermostat_tau <= 0:
            raise ValueError(f"Molecular thermostat tau must be positive, got {self.molecular_thermostat_tau} ps")
        
        if self.cavity_thermostat_tau <= 0:
            raise ValueError(f"Cavity thermostat tau must be positive, got {self.cavity_thermostat_tau} ps")
        
        # Special validation for Langevin thermostats
        if self.molecular_thermostat == 'langevin' and self.molecular_thermostat_tau <= 0:
            raise ValueError(
                f"Langevin molecular thermostat requires tau > 0, got {self.molecular_thermostat_tau} ps. "
                f"For overdamped dynamics, use Brownian dynamics instead."
            )
        
        if self.cavity_thermostat == 'langevin' and self.cavity_thermostat_tau <= 0:
            raise ValueError(
                f"Langevin cavity thermostat requires tau > 0, got {self.cavity_thermostat_tau} ps. "
                f"For overdamped dynamics, use Brownian dynamics instead."
            )
    
    def _validate_cavity_config(self) -> None:
        """Validate cavity configuration parameters."""
        if not isinstance(self.incavity, bool):
            raise ValueError(f"incavity must be boolean, got {type(self.incavity)}")
        
        if not isinstance(self.finite_q, bool):
            raise ValueError(f"finite_q must be boolean, got {type(self.finite_q)}")
        
        if not isinstance(self.add_cavity_particle, bool):
            raise ValueError(f"add_cavity_particle must be boolean, got {type(self.add_cavity_particle)}")
        
        if not isinstance(self.cavity_damping_factor, (int, float)):
            raise ValueError(f"Cavity damping factor must be a number, got {type(self.cavity_damping_factor)}")
        
        if self.cavity_damping_factor <= 0:
            raise ValueError(f"Cavity damping factor must be positive, got {self.cavity_damping_factor}")
    
    def _validate_time_varying_params(self) -> None:
        """Validate time-varying parameter configurations."""
        if self.switch_time_ps is not None:
            if not isinstance(self.switch_time_ps, (int, float)):
                raise ValueError(f"Switch time must be a number, got {type(self.switch_time_ps)}")
            if self.switch_time_ps < 0:
                raise ValueError(f"Switch time must be non-negative, got {self.switch_time_ps} ps")
            if self.switch_time_ps >= self.runtime_ps:
                raise ValueError(f"Switch time ({self.switch_time_ps} ps) must be less than runtime ({self.runtime_ps} ps)")
        
        if not isinstance(self.dissipation, (int, float)):
            raise ValueError(f"Dissipation must be a number, got {type(self.dissipation)}")
        if self.dissipation < 0:
            raise ValueError(f"Dissipation must be non-negative, got {self.dissipation}")
    
    def _validate_simulation_control(self) -> None:
        """Validate simulation control parameters."""
        if not isinstance(self.error_tolerance, (int, float)):
            raise ValueError(f"Error tolerance must be a number, got {type(self.error_tolerance)}")
        if self.error_tolerance < 0:
            raise ValueError(f"Error tolerance must be non-negative, got {self.error_tolerance}")
        
        if self.dt_fs is not None:
            if not isinstance(self.dt_fs, (int, float)):
                raise ValueError(f"Timestep must be a number, got {type(self.dt_fs)}")
            if self.dt_fs <= 0:
                raise ValueError(f"Timestep must be positive, got {self.dt_fs} fs")
            if self.dt_fs > 10.0:
                logging.warning(f"Large timestep: {self.dt_fs} fs")
        
        if self.device not in ['CPU', 'GPU']:
            raise ValueError(f"Device must be 'CPU' or 'GPU', got {self.device}")
        
        if not isinstance(self.gpu_id, int) or self.gpu_id < 0:
            raise ValueError(f"GPU ID must be a non-negative integer, got {self.gpu_id}")
        
        if self.seed is not None and (not isinstance(self.seed, int) or self.seed < 0):
            raise ValueError(f"Seed must be a non-negative integer, got {self.seed}")
    
    def _validate_output_config(self) -> None:
        """Validate output configuration parameters."""
        periods = [
            ('energy_output_period_ps', self.energy_output_period_ps),
            ('fkt_output_period_ps', self.fkt_output_period_ps),
            ('gsd_output_period_ps', self.gsd_output_period_ps),
            ('console_output_period_ps', self.console_output_period_ps)
        ]
        
        for name, period in periods:
            if not isinstance(period, (int, float)):
                raise ValueError(f"{name} must be a number, got {type(period)}")
            if period <= 0:
                raise ValueError(f"{name} must be positive, got {period} ps")
    
    def _validate_analysis_params(self) -> None:
        """Validate analysis parameter configurations."""
        if not isinstance(self.enable_fkt, bool):
            raise ValueError(f"enable_fkt must be boolean, got {type(self.enable_fkt)}")
        
        if not isinstance(self.fkt_kmag, (int, float)) or self.fkt_kmag <= 0:
            raise ValueError(f"F(k,t) k-magnitude must be positive, got {self.fkt_kmag}")
        
        if not isinstance(self.fkt_num_wavevectors, int) or self.fkt_num_wavevectors <= 0:
            raise ValueError(f"Number of wavevectors must be positive integer, got {self.fkt_num_wavevectors}")
        
        if not isinstance(self.fkt_reference_interval_ps, (int, float)) or self.fkt_reference_interval_ps <= 0:
            raise ValueError(f"F(k,t) reference interval must be positive, got {self.fkt_reference_interval_ps} ps")
        
        if not isinstance(self.fkt_max_references, int) or self.fkt_max_references <= 0:
            raise ValueError(f"Max references must be positive integer, got {self.fkt_max_references}")
        
        if not isinstance(self.enable_energy_tracking, bool):
            raise ValueError(f"enable_energy_tracking must be boolean, got {type(self.enable_energy_tracking)}")
        
        if self.max_energy_output_time_ps is not None:
            if not isinstance(self.max_energy_output_time_ps, (int, float)):
                raise ValueError(f"Max energy output time must be a number, got {type(self.max_energy_output_time_ps)}")
            if self.max_energy_output_time_ps <= 0:
                raise ValueError(f"Max energy output time must be positive, got {self.max_energy_output_time_ps} ps")
    
    def _validate_io_config(self) -> None:
        """Validate input/output configuration."""
        if not isinstance(self.job_dir, str):
            raise ValueError(f"Job directory must be a string, got {type(self.job_dir)}")
        
        if not isinstance(self.replica, int) or self.replica < 0:
            raise ValueError(f"Replica must be a non-negative integer, got {self.replica}")
        
        if not isinstance(self.input_gsd, str):
            raise ValueError(f"Input GSD must be a string, got {type(self.input_gsd)}")
        
        if not isinstance(self.frame, int):
            raise ValueError(f"Frame must be an integer, got {type(self.frame)}")
        
        if not isinstance(self.name, str):
            raise ValueError(f"Name must be a string, got {type(self.name)}")
        
        if self.log_level not in ['DEBUG', 'INFO', 'WARNING', 'ERROR']:
            raise ValueError(f"Log level must be one of ['DEBUG', 'INFO', 'WARNING', 'ERROR'], got {self.log_level}")
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert parameters to dictionary for serialization."""
        return {
            field.name: getattr(self, field.name)
            for field in self.__dataclass_fields__.values()
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'CavitySimulationParams':
        """Create parameters from dictionary."""
        return cls(**data)
    
    def get_physical_constants(self) -> Dict[str, float]:
        """Get relevant physical constants for the simulation."""
        return {
            'kB': PhysicalConstants.KB_HARTREE_PER_K,
            'kT': PhysicalConstants.KB_HARTREE_PER_K * self.temperature,
            'omegac_au': self.frequency / PhysicalConstants.HARTREE_TO_CM_MINUS1,
            'molecular_tau_au': PhysicalConstants.ps_to_atomic_units(self.molecular_thermostat_tau),
            'cavity_tau_au': PhysicalConstants.ps_to_atomic_units(self.cavity_thermostat_tau),
        }
    
    def get_summary(self) -> str:
        """Get a formatted summary of the simulation parameters."""
        summary = []
        summary.append("=== Cavity MD Simulation Parameters ===")
        summary.append(f"Temperature: {self.temperature} K")
        summary.append(f"Runtime: {self.runtime_ps} ps")
        summary.append(f"Cavity: {'Enabled' if self.incavity else 'Disabled'}")
        
        if self.incavity:
            summary.append(f"  Frequency: {self.frequency} cm^-1")
            summary.append(f"  Coupling: {self.coupling_strength} a.u.")
            summary.append(f"  Finite-q: {self.finite_q}")
        
        summary.append(f"Molecular thermostat: {self.molecular_thermostat} (tau={self.molecular_thermostat_tau} ps)")
        if self.incavity:
            summary.append(f"Cavity thermostat: {self.cavity_thermostat} (tau={self.cavity_thermostat_tau} ps)")
        
        summary.append(f"Device: {self.device}")
        summary.append(f"Seed: {self.seed if self.seed is not None else 'replica-based'}")
        
        return "\n".join(summary)


def validate_hoomd_device(device_type: str, gpu_id: int = 0) -> hoomd.device.Device:
    """
    Validate and create a HOOMD device.
    
    Parameters
    ----------
    device_type : str
        Device type ('CPU' or 'GPU')
    gpu_id : int
        GPU ID for GPU devices
        
    Returns
    -------
    hoomd.Device
        Validated HOOMD device
        
    Raises
    ------
    ValueError
        If device configuration is invalid
    """
    device_type = device_type.upper()
    
    if device_type not in ['CPU', 'GPU']:
        raise ValueError(f"Device type must be 'CPU' or 'GPU', got {device_type}")
    
    if device_type == 'CPU':
        return hoomd.device.CPU()
    else:
        if not isinstance(gpu_id, int) or gpu_id < 0:
            raise ValueError(f"GPU ID must be a non-negative integer, got {gpu_id}")
        
        try:
            # Try different GPU initialization methods for different HOOMD versions
            try:
                # First try with gpu_ids parameter (newer HOOMD versions)
                return hoomd.device.GPU(gpu_ids=[gpu_id])
            except TypeError:
                # Fall back to older parameter name
                return hoomd.device.GPU(gpu_id=gpu_id)
        except Exception as e:
            raise ValueError(f"Failed to initialize GPU device {gpu_id}: {e}")


def validate_input_file(filepath: str, frame: int = -1) -> None:
    """
    Validate input GSD file existence and frame access.
    
    Parameters
    ----------
    filepath : str
        Path to GSD file
    frame : int
        Frame index to validate
        
    Raises
    ------
    ValueError
        If file doesn't exist or frame is invalid
    FileNotFoundError
        If file doesn't exist
    """
    import gsd.hoomd
    
    file_path = Path(filepath)
    if not file_path.exists():
        raise FileNotFoundError(f"Input GSD file not found: {filepath}")
    
    try:
        with gsd.hoomd.open(filepath, 'r') as f:
            num_frames = len(f)
            if frame < 0:
                actual_frame = num_frames + frame
            else:
                actual_frame = frame
            
            if actual_frame < 0 or actual_frame >= num_frames:
                raise ValueError(f"Invalid frame index {frame} for file with {num_frames} frames")
                
    except Exception as e:
        raise ValueError(f"Error reading GSD file {filepath}: {e}")


def validate_directory_structure(job_dir: str, create_if_missing: bool = True) -> Path:
    """
    Validate and optionally create job directory structure.
    
    Parameters
    ----------
    job_dir : str
        Job directory path
    create_if_missing : bool
        Whether to create directory if it doesn't exist
        
    Returns
    -------
    Path
        Validated job directory path
        
    Raises
    ------
    ValueError
        If directory doesn't exist and create_if_missing is False
    """
    job_path = Path(job_dir)
    
    if not job_path.exists():
        if create_if_missing:
            job_path.mkdir(parents=True, exist_ok=True)
            logging.info(f"Created job directory: {job_path}")
        else:
            raise ValueError(f"Job directory does not exist: {job_dir}")
    
    return job_path 