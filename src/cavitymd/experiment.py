#!/usr/bin/env python3
"""
High-level experiment framework for cavity MD simulations.

This module provides the CavityMDExperiment class which encapsulates all
experiment logic and replaces the monolithic run script functionality.
"""

import argparse
import logging
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import replace

import hoomd
import numpy as np

from .simulation import CavityMDSimulation
from .validation import CavitySimulationParams
from .analysis import PhysicalConstants
from .utils import format_coupling_strength


class CavityMDExperiment:
    """
    High-level experiment runner that orchestrates cavity MD simulations.
    
    This class provides the same interface as the original run script but with
    proper separation of concerns and leveraging the plugin framework.
    
    Features:
    - SLURM array job support and local multi-replica execution
    - Comprehensive CLI interface with backward compatibility
    - Proper error handling and logging
    - Integration with the plugin-based simulation framework
    """
    
    def __init__(self, **kwargs):
        """
        Initialize the experiment with simulation parameters.
        
        Args:
            **kwargs: Simulation parameters matching the original run script interface
        """
        self.start_time = time.time()
        self.successful_experiments = 0
        self.failed_experiments = 0
        
        # Convert kwargs to validated parameters
        self.params = self._convert_kwargs_to_params(kwargs)
        
        # Determine replica list
        self.replica_list = self._determine_replica_list(
            kwargs.get('replicas'), 
            kwargs.get('slurm_task_id'), 
            kwargs.get('slurm_job_id')
        )
        
        # Setup logging
        self._setup_logging()
        
    def _convert_kwargs_to_params(self, kwargs: Dict) -> CavitySimulationParams:
        """
        Convert experiment kwargs to validated CavitySimulationParams.
        
        This method maps the original run script parameters to the new
        validation framework, ensuring backward compatibility.
        """
        # Map original parameter names to new parameter names
        param_mapping = {
            'input_gsd': 'input_gsd',
            'molecular_bath': 'molecular_thermostat',
            'cavity_bath': 'cavity_thermostat',
            'molecular_tau': 'molecular_thermostat_tau',
            'cavity_tau': 'cavity_thermostat_tau',
            'no_cavity': 'incavity',  # Will be inverted
            'finite_q': 'finite_q',
            'coupling': 'coupling_strength',
            'switch_time': 'switch_time_ps',
            'decay_time_constant': 'decay_time_constant_ps',
            'temperature': 'temperature',
            'frequency': 'frequency',
            'runtime': 'runtime_ps',
            'device': 'device',
            'gpu_id': 'gpu_id',
            'timestep': 'dt_fs',
            'enable_energy_tracker': 'enable_energy_tracking',
            'energy_output_period_ps': 'energy_output_period_ps',
            'fkt_output_period_ps': 'fkt_output_period_ps',
            'gsd_output_period_ps': 'gsd_output_period_ps',
            'console_output_period_ps': 'console_output_period_ps',
            'enable_fkt': 'enable_fkt',
            'fkt_kmag': 'fkt_kmag',
            'fkt_wavevectors': 'fkt_num_wavevectors',
            'fkt_ref_interval': 'fkt_reference_interval_ps',
            'fkt_max_refs': 'fkt_max_references',
            'max_energy_output_time': 'max_energy_output_time_ps',
            'truncate_gsd': 'truncate_gsd',
            'seed': 'seed',
            'no_restart_velocities': 'restart_velocities',  # Will be inverted
        }
        
        # Convert parameters
        converted_params = {}
        for orig_name, new_name in param_mapping.items():
            if orig_name in kwargs:
                value = kwargs[orig_name]
                # Handle inverted boolean parameters
                if orig_name == 'no_cavity':
                    converted_params['incavity'] = not value
                elif orig_name == 'no_restart_velocities':
                    converted_params['restart_velocities'] = not value
                else:
                    converted_params[new_name] = value
        
        # Handle fixed_timestep parameter separately
        if 'fixed_timestep' in kwargs and kwargs['fixed_timestep']:
            # For fixed timestep mode, use the specified timestep
            if 'timestep' in kwargs:
                converted_params['dt_fs'] = kwargs['timestep']
            else:
                converted_params['dt_fs'] = 1.0  # Default timestep
        else:
            # For adaptive timestep mode, dt_fs should be None
            converted_params['dt_fs'] = None
        
        # Calculate dissipation from damping ratio if provided
        damping_ratio = kwargs.get('damping_ratio', 0.0)
        if damping_ratio > 0:
            phmass = 1.0  # Photon mass is 1.0 in a.u.
            frequency = kwargs.get('frequency', 2000.0)
            omegac = frequency / PhysicalConstants.HARTREE_TO_CM_MINUS1
            dissipation = 2 * damping_ratio * phmass * omegac
            converted_params['dissipation'] = dissipation
        else:
            converted_params['dissipation'] = 0.0
        
        # Set defaults for parameters not in kwargs
        defaults = {
            'molecular_thermostat': 'bussi',
            'cavity_thermostat': 'langevin',
            'incavity': True,
            'finite_q': False,
            'coupling_strength': 1e-3,
            'temperature': 100.0,
            'frequency': 2000.0,
            'runtime_ps': 500.0,
            'molecular_thermostat_tau': 5.0,
            'cavity_thermostat_tau': 5.0,
            'device': 'CPU',
            'gpu_id': 0,
            'enable_energy_tracking': False,
            'energy_output_period_ps': 0.1,
            'fkt_output_period_ps': 1.0,
            'gsd_output_period_ps': 50.0,
            'console_output_period_ps': 1.0,
            'enable_fkt': False,
            'fkt_kmag': 1.0,
            'fkt_num_wavevectors': 50,
            'fkt_reference_interval_ps': 1.0,
            'fkt_max_references': 10,
            'truncate_gsd': False,
            'restart_velocities': True,
            'input_gsd': 'molecular-0.gsd',
            'frame': -1,
            'name': 'prod',
            'error_tolerance': 0.01,
            'log_level': 'INFO',
        }
        
        # Merge defaults with converted parameters
        for key, default_value in defaults.items():
            if key not in converted_params:
                converted_params[key] = default_value
        
        # Handle special parameter calculations
        if not converted_params['incavity']:
            converted_params['cavity_thermostat'] = 'none'
        
        # Set error tolerance based on timestep mode
        if converted_params['dt_fs'] is not None:
            # Fixed timestep mode
            converted_params['error_tolerance'] = 0.0
        else:
            # Adaptive timestep mode
            converted_params['error_tolerance'] = 1.0
        
        # Store damping ratio separately for access in experiment
        self.damping_ratio = damping_ratio
        
        return CavitySimulationParams(**converted_params)
    
    def _determine_replica_list(self, replicas_str: Optional[str], 
                               slurm_task_id: Optional[int], 
                               slurm_job_id: Optional[int]) -> List[int]:
        """
        Determine the list of replicas to run based on SLURM or local execution.
        
        Args:
            replicas_str: Replica specification string (e.g., "1,2,3" or "1-5")
            slurm_task_id: SLURM array task ID (if running under SLURM)
            slurm_job_id: SLURM job ID (if running under SLURM)
            
        Returns:
            List of replica IDs to execute
        """
        if slurm_task_id is not None:
            # Running under SLURM array job
            return [slurm_task_id]
        else:
            # Local execution - parse replicas
            return self._parse_replicas(replicas_str)
    
    def _parse_replicas(self, replicas_str: Optional[str]) -> List[int]:
        """
        Parse replica specification string into a list of replica IDs.
        
        Args:
            replicas_str: Replica specification (e.g., "1,2,3" or "1-5")
            
        Returns:
            List of replica IDs
        """
        if not replicas_str:
            return [1]  # Default to single replica
        
        replicas = []
        parts = replicas_str.split(',')
        
        for part in parts:
            part = part.strip()
            if '-' in part:
                # Handle range specification (e.g., "1-5")
                start, end = part.split('-')
                replicas.extend(range(int(start), int(end) + 1))
            else:
                # Handle single replica
                replicas.append(int(part))
        
        return sorted(list(set(replicas)))  # Remove duplicates and sort
    
    def _setup_logging(self):
        """Setup logging configuration for the experiment."""
        logging.basicConfig(
            level=getattr(logging, self.params.log_level.upper()),
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger('CavityMDExperiment')
    
    def _create_experiment_directory(self, replica: int) -> Path:
        """
        Create experiment directory with appropriate naming.
        
        Args:
            replica: Replica ID
            
        Returns:
            Path to the experiment directory
        """
        if self.params.incavity:
            # For cavity simulations, include coupling strength in directory name
            coupling_str = format_coupling_strength(self.params.coupling_strength)
            if self.params.switch_time_ps is not None:
                # Include switch time in directory name for time-varying simulations
                switch_str = f"_switch_{self.params.switch_time_ps}ps"
                exp_dir = Path(f"{coupling_str}{switch_str}")
            else:
                exp_dir = Path(f"{coupling_str}")
        else:
            # For non-cavity simulations
            exp_dir = Path("no_cavity")
        
        exp_dir.mkdir(exist_ok=True)
        return exp_dir
    
    def _run_single_experiment(self, replica: int) -> bool:
        """
        Run a single experiment for the given replica.
        
        Args:
            replica: Replica ID
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Create experiment directory
            exp_dir = self._create_experiment_directory(replica)
            
            self.logger.info(f"Starting replica {replica} in {exp_dir}")
            self.logger.info(f"Cavity coupling: {'Enabled' if self.params.incavity else 'Disabled'}")
            
            if self.params.incavity:
                self.logger.info(f"Coupling strength: {self.params.coupling_strength}")
                self.logger.info(f"Frequency: {self.params.frequency} cm⁻¹")
                self.logger.info(f"Finite-q mode: {self.params.finite_q}")
                self.logger.info(f"Cavity thermostat: {self.params.cavity_thermostat}")
                if self.params.switch_time_ps is not None:
                    self.logger.info(f"Switch time: {self.params.switch_time_ps} ps")
                    self.logger.info(f"Damping ratio: {self.damping_ratio}")
                    self.logger.info(f"Dissipation: {self.params.dissipation:.4e} a.u.")
            
            self.logger.info(f"Molecular thermostat: {self.params.molecular_thermostat}")
            self.logger.info(f"Temperature: {self.params.temperature} K")
            self.logger.info(f"Runtime: {self.params.runtime_ps} ps")
            self.logger.info(f"Device: {self.params.device}")
            
            # Update parameters for this replica
            replica_params = replace(
                self.params,
                job_dir=str(exp_dir),
                replica=replica,
                frame=self.params.frame,  # Use frame from CLI args instead of replica number
                input_gsd=self.params.input_gsd  # Use CLI-provided input GSD path directly
            )
            
            # Create and run simulation
            simulation = CavityMDSimulation(replica_params)
            simulation.setup()
            simulation.run()
            
            self.logger.info(f"SUCCESS: Replica {replica} completed successfully")
            return True
            
        except Exception as e:
            self.logger.error(f"ERROR: Replica {replica} failed: {e}")
            return False
    
    def run(self) -> int:
        """
        Run the complete experiment for all replicas.
        
        Returns:
            Exit code (0 for success, 1 for failure)
        """
        print("Advanced Cavity MD Experiment Runner")
        print("=" * 50)
        
        # Print experiment configuration
        print(f"\nSimulation Configuration:")
        print(f"  Cavity coupling: {'Enabled' if self.params.incavity else 'Disabled'}")
        if self.params.incavity:
            print(f"    Coupling strength: {self.params.coupling_strength}")
            print(f"    Frequency: {self.params.frequency} cm⁻¹")
            print(f"    Finite-q mode: {self.params.finite_q}")
            print(f"    Cavity thermostat: {self.params.cavity_thermostat}")
        print(f"  Molecular thermostat: {self.params.molecular_thermostat}")
        print(f"  Temperature: {self.params.temperature} K")
        print(f"  Runtime: {self.params.runtime_ps} ps")
        print(f"  Device: {self.params.device}")
        if self.params.device == 'GPU':
            print(f"    GPU ID: {self.params.gpu_id}")
        print(f"  Random seed: {self.params.seed if self.params.seed is not None else 'replica-based (deterministic)'}")
        print(f"  Velocity restart: {'Disabled - using GSD velocities' if not self.params.restart_velocities else 'Enabled - thermalizing velocities'}")
        
        print(f"\nStarting execution for {len(self.replica_list)} replica(s): {self.replica_list}")
        print("=" * 50)
        
        # Run replicas
        for replica in self.replica_list:
            print(f"\nRunning replica {replica}...")
            
            success = self._run_single_experiment(replica)
            
            if success:
                self.successful_experiments += 1
                print(f"SUCCESS: Replica {replica} completed successfully")
            else:
                self.failed_experiments += 1
                print(f"ERROR: Replica {replica} failed")
        
        # Final summary
        self._print_summary()
        
        return 0 if self.failed_experiments == 0 else 1
    
    def _print_summary(self):
        """Print final execution summary."""
        end_time = time.time()
        total_wall_time = end_time - self.start_time
        total_experiments = len(self.replica_list)
        
        print("\n" + "=" * 50)
        print("Execution Summary")
        print("=" * 50)
        print(f"Total replicas: {total_experiments}")
        print(f"Successful: {self.successful_experiments}")
        print(f"Failed: {self.failed_experiments}")
        print(f"Wall time: {total_wall_time:.2f} seconds")
        
        if self.failed_experiments > 0:
            print(f"\nWARNING: {self.failed_experiments} replicas failed - check individual logs for details")
        else:
            print("\nAll replicas completed successfully!")
    
    @classmethod
    def from_cli_args(cls, args: argparse.Namespace) -> 'CavityMDExperiment':
        """
        Create experiment from command-line arguments.
        
        Args:
            args: Parsed command-line arguments
            
        Returns:
            CavityMDExperiment instance
        """
        # Get SLURM information
        slurm_task_id, slurm_job_id = cls._get_slurm_info()
        
        # Convert argparse Namespace to kwargs
        kwargs = vars(args).copy()
        kwargs['slurm_task_id'] = slurm_task_id
        kwargs['slurm_job_id'] = slurm_job_id
        
        return cls(**kwargs)
    
    @staticmethod
    def _get_slurm_info() -> Tuple[Optional[int], Optional[int]]:
        """
        Get SLURM job information from environment variables.
        
        Returns:
            Tuple of (task_id, job_id) or (None, None) if not running under SLURM
        """
        try:
            task_id = int(os.environ.get('SLURM_ARRAY_TASK_ID', ''))
            job_id = int(os.environ.get('SLURM_ARRAY_JOB_ID', ''))
            return task_id, job_id
        except (ValueError, TypeError):
            return None, None 