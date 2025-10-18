#!/usr/bin/env python3
"""
LQG Controller Calibration and System Identification Tools

This script provides tools for calibrating the LQG controller parameters
from experimental or simulation data. It includes:

1. Data collection from existing simulations
2. System identification and parameter fitting
3. Controller validation and performance analysis
4. Parameter optimization and tuning

Usage:
    python lqg_calibration.py --mode collect --input simulation_data.csv
    python lqg_calibration.py --mode identify --data pulse_data.json
    python lqg_calibration.py --mode validate --params lqg_params.json
    python lqg_calibration.py --mode optimize --target performance_spec.json
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import csv
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import warnings
from dataclasses import asdict

# Add src to path for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from cavitymd.lqg_controller import (
    LQGParameters, LiftedStateSpaceModel, KalmanFilter, 
    LQRController, SystemIdentification
)


class DataCollector:
    """
    Collect and preprocess data for system identification.
    
    This class processes simulation output files to extract pulse-to-pulse
    data suitable for system identification.
    """
    
    def __init__(self, pulse_period_ps: float = 5.0, duty_cycle: float = 0.5):
        """
        Initialize data collector.
        
        Parameters
        ----------
        pulse_period_ps : float
            Expected pulse period in picoseconds
        duty_cycle : float
            Expected duty cycle
        """
        self.pulse_period_ps = pulse_period_ps
        self.duty_cycle = duty_cycle
    
    def collect_from_csv(self, csv_file: str, 
                        time_col: str = 'time_ps',
                        temp_cols: Optional[List[str]] = None,
                        control_cols: Optional[List[str]] = None) -> List[Dict[str, Any]]:
        """
        Collect pulse data from CSV simulation output.
        
        Parameters
        ----------
        csv_file : str
            Path to CSV file with simulation data
        time_col : str
            Name of time column
        temp_cols : List[str], optional
            Names of temperature columns [T_s, T_k, T_v]
        control_cols : List[str], optional
            Names of control columns [coupling, bath_temp]
            
        Returns
        -------
        pulse_data : List[Dict]
            List of pulse data points
        """
        if temp_cols is None:
            temp_cols = ['T_structure', 'T_kinetic', 'T_vibrational']
        
        if control_cols is None:
            control_cols = ['coupling_strength', 'bath_temperature']
        
        # Read CSV data
        df = pd.read_csv(csv_file)
        
        # Validate columns
        required_cols = [time_col] + temp_cols + control_cols
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing columns in CSV: {missing_cols}")
        
        # Detect pulse boundaries
        pulse_times = self._detect_pulse_boundaries(df[time_col].values)
        
        print(f"Detected {len(pulse_times)} pulse boundaries")
        
        # Extract pulse data
        pulse_data = []
        for i in range(len(pulse_times) - 1):
            current_time = pulse_times[i]
            next_time = pulse_times[i + 1]
            
            # Find data points closest to pulse boundaries
            current_idx = np.argmin(np.abs(df[time_col] - current_time))
            next_idx = np.argmin(np.abs(df[time_col] - next_time))
            
            # Extract state and control data
            current_state = [df[col].iloc[current_idx] for col in temp_cols]
            next_state = [df[col].iloc[next_idx] for col in temp_cols]
            
            # Average control over the pulse period
            period_mask = (df[time_col] >= current_time) & (df[time_col] < next_time)
            avg_control = [df[col][period_mask].mean() for col in control_cols]
            
            # Convert coupling to deviation (subtract baseline)
            if len(avg_control) >= 1:
                baseline_coupling = 1e-5  # Assumed baseline
                coupling_deviation = avg_control[0] - baseline_coupling
                avg_control[0] = coupling_deviation
            
            pulse_data.append({
                'time_ps': current_time,
                'state': current_state,
                'control': avg_control,
                'next_state': next_state,
                'period_ps': next_time - current_time
            })
        
        return pulse_data
    
    def _detect_pulse_boundaries(self, times: np.ndarray) -> np.ndarray:
        """
        Detect pulse boundaries from time series.
        
        Parameters
        ----------
        times : np.ndarray
            Array of time values
            
        Returns
        -------
        pulse_times : np.ndarray
            Array of pulse boundary times
        """
        # Simple method: assume regular spacing
        start_time = times[0]
        end_time = times[-1]
        
        # Generate expected pulse times
        pulse_times = []
        current_time = start_time
        
        while current_time <= end_time:
            pulse_times.append(current_time)
            current_time += self.pulse_period_ps
        
        return np.array(pulse_times)
    
    def collect_from_energy_tracker(self, energy_file: str, 
                                   temp_tracker_file: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Collect data from energy tracker and temperature tracker outputs.
        
        Parameters
        ----------
        energy_file : str
            Path to energy tracker CSV file
        temp_tracker_file : str, optional
            Path to temperature tracker CSV file
            
        Returns
        -------
        pulse_data : List[Dict]
            List of pulse data points
        """
        # Read energy data
        energy_df = pd.read_csv(energy_file)
        
        # Read temperature data if available
        if temp_tracker_file and Path(temp_tracker_file).exists():
            temp_df = pd.read_csv(temp_tracker_file)
            # Merge on time
            df = pd.merge_asof(energy_df.sort_values('time_ps'), 
                              temp_df.sort_values('time_ps'), 
                              on='time_ps', direction='nearest')
        else:
            df = energy_df
            # Calculate kinetic temperature from energy data
            df['T_kinetic'] = self._calculate_kinetic_temp_from_energy(df)
        
        # Use energy-based temperature calculations
        df['T_structure'] = self._calculate_structure_temp_from_energy(df)
        df['T_vibrational'] = self._calculate_vibrational_temp_from_energy(df)
        
        # Assume constant control (would need to be provided separately)
        df['coupling_strength'] = 1e-4  # Placeholder
        df['bath_temperature'] = 300.0  # Placeholder
        
        # Collect pulse data
        return self.collect_from_csv(
            csv_file=None,  # We already have the DataFrame
            time_col='time_ps',
            temp_cols=['T_structure', 'T_kinetic', 'T_vibrational'],
            control_cols=['coupling_strength', 'bath_temperature']
        )
    
    def _calculate_kinetic_temp_from_energy(self, df: pd.DataFrame) -> np.ndarray:
        """Calculate kinetic temperature from energy data."""
        # Placeholder - would use actual kinetic energy if available
        return np.full(len(df), 300.0)
    
    def _calculate_structure_temp_from_energy(self, df: pd.DataFrame) -> np.ndarray:
        """Calculate structure temperature from LJ+Coulombic energy."""
        if 'lj' in df.columns and 'coulombic' in df.columns:
            total_energy = df['lj'] + df['coulombic']
            # Simple linear relationship (would use empirical calibration)
            return 200.0 + 0.1 * total_energy
        else:
            return np.full(len(df), 300.0)
    
    def _calculate_vibrational_temp_from_energy(self, df: pd.DataFrame) -> np.ndarray:
        """Calculate vibrational temperature from harmonic energy."""
        if 'harmonic' in df.columns:
            # Simple relationship
            return 250.0 + 0.2 * df['harmonic']
        else:
            return np.full(len(df), 300.0)
    
    def save_pulse_data(self, pulse_data: List[Dict[str, Any]], output_file: str) -> None:
        """Save pulse data to JSON file."""
        with open(output_file, 'w') as f:
            json.dump(pulse_data, f, indent=2)
        
        print(f"Saved {len(pulse_data)} pulse data points to {output_file}")
    
    def load_pulse_data(self, input_file: str) -> List[Dict[str, Any]]:
        """Load pulse data from JSON file."""
        with open(input_file, 'r') as f:
            pulse_data = json.load(f)
        
        print(f"Loaded {len(pulse_data)} pulse data points from {input_file}")
        return pulse_data


class ParameterOptimizer:
    """
    Optimize LQG controller parameters for desired performance.
    
    This class provides tools for tuning controller parameters based on
    performance specifications and constraints.
    """
    
    def __init__(self):
        """Initialize parameter optimizer."""
        self.optimization_history = []
    
    def optimize_for_tracking(self, pulse_data: List[Dict[str, Any]],
                            target_temperatures: Dict[str, float],
                            performance_weights: Optional[Dict[str, float]] = None) -> LQGParameters:
        """
        Optimize parameters for tracking performance.
        
        Parameters
        ----------
        pulse_data : List[Dict]
            Training data for optimization
        target_temperatures : Dict[str, float]
            Target temperatures for each output
        performance_weights : Dict[str, float], optional
            Weights for different performance metrics
            
        Returns
        -------
        optimized_params : LQGParameters
            Optimized controller parameters
        """
        if performance_weights is None:
            performance_weights = {
                'tracking_error': 1.0,
                'control_effort': 0.1,
                'stability_margin': 0.5
            }
        
        # Initial parameter guess
        initial_params = LQGParameters()
        
        # Set up optimization problem
        def objective(param_vector):
            params = self._vector_to_params(param_vector, initial_params)
            return self._evaluate_performance(params, pulse_data, target_temperatures, performance_weights)
        
        # Parameter bounds
        bounds = self._get_optimization_bounds()
        
        # Initial parameter vector
        x0 = self._params_to_vector(initial_params)
        
        # Optimize using scipy
        from scipy.optimize import minimize
        
        result = minimize(
            objective, x0, method='L-BFGS-B', bounds=bounds,
            options={'maxiter': 100, 'disp': True}
        )
        
        if not result.success:
            warnings.warn(f"Optimization failed: {result.message}")
        
        # Convert back to parameters
        optimized_params = self._vector_to_params(result.x, initial_params)
        
        print(f"Optimization completed with final cost: {result.fun:.4f}")
        
        return optimized_params
    
    def _evaluate_performance(self, params: LQGParameters, 
                            pulse_data: List[Dict[str, Any]],
                            target_temperatures: Dict[str, float],
                            weights: Dict[str, float]) -> float:
        """Evaluate controller performance for given parameters."""
        try:
            # Create controller components
            model = LiftedStateSpaceModel(params)
            kf = KalmanFilter(model)
            lqr = LQRController(model, params)
            
            # Simulate controller on data
            total_tracking_error = 0.0
            total_control_effort = 0.0
            n_points = 0
            
            reference = np.array([
                target_temperatures.get('structure', 300.0),
                target_temperatures.get('kinetic', 300.0),
                target_temperatures.get('vibrational', 300.0)
            ])
            
            for i, data_point in enumerate(pulse_data[:50]):  # Limit for speed
                state = np.array(data_point['state'])
                if len(state) == 3:
                    state = np.append(state, 0.0)  # Add hidden state
                
                # Kalman filter step
                y_k = state[:3]
                u_prev = np.array([0.0, 300.0]) if i == 0 else u
                x_hat = kf.step(y_k, u_prev)
                
                # Control computation
                u = lqr.compute_control(x_hat, reference)
                
                # Compute metrics
                tracking_error = np.linalg.norm(y_k - reference)
                control_effort = np.linalg.norm(u)
                
                total_tracking_error += tracking_error
                total_control_effort += control_effort
                n_points += 1
            
            if n_points == 0:
                return 1e6  # Large penalty
            
            # Average metrics
            avg_tracking_error = total_tracking_error / n_points
            avg_control_effort = total_control_effort / n_points
            
            # Stability check
            eigenvals = np.linalg.eigvals(model.A)
            max_eigenval = np.max(np.abs(eigenvals))
            stability_penalty = max(0, max_eigenval - 0.95) * 1000  # Penalty if close to instability
            
            # Combined cost
            cost = (weights['tracking_error'] * avg_tracking_error +
                   weights['control_effort'] * avg_control_effort +
                   weights['stability_margin'] * stability_penalty)
            
            return cost
            
        except Exception as e:
            warnings.warn(f"Performance evaluation failed: {e}")
            return 1e6  # Large penalty for failed evaluation
    
    def _params_to_vector(self, params: LQGParameters) -> np.ndarray:
        """Convert parameters to optimization vector."""
        return np.array([
            params.Q_structure, params.Q_kinetic, params.Q_vibrational,
            params.R_coupling, params.R_bath,
            params.W_structure, params.W_kinetic, params.W_vibrational,
            params.V_structure, params.V_kinetic, params.V_vibrational
        ])
    
    def _vector_to_params(self, vector: np.ndarray, base_params: LQGParameters) -> LQGParameters:
        """Convert optimization vector to parameters."""
        new_params = LQGParameters(**base_params.to_dict())
        
        new_params.Q_structure = vector[0]
        new_params.Q_kinetic = vector[1]
        new_params.Q_vibrational = vector[2]
        new_params.R_coupling = vector[3]
        new_params.R_bath = vector[4]
        new_params.W_structure = vector[5]
        new_params.W_kinetic = vector[6]
        new_params.W_vibrational = vector[7]
        new_params.V_structure = vector[8]
        new_params.V_kinetic = vector[9]
        new_params.V_vibrational = vector[10]
        
        return new_params
    
    def _get_optimization_bounds(self) -> List[Tuple[float, float]]:
        """Get bounds for optimization variables."""
        return [
            # Q weights (tracking)
            (0.1, 100.0),   # Q_structure
            (0.01, 10.0),   # Q_kinetic
            (0.01, 10.0),   # Q_vibrational
            # R weights (control effort)
            (0.001, 1.0),   # R_coupling
            (0.001, 1.0),   # R_bath
            # Process noise
            (1e-6, 1e-2),   # W_structure
            (1e-6, 1e-2),   # W_kinetic
            (1e-6, 1e-2),   # W_vibrational
            # Measurement noise
            (1e-6, 1e-2),   # V_structure
            (1e-6, 1e-2),   # V_kinetic
            (1e-6, 1e-2),   # V_vibrational
        ]


class ControllerValidator:
    """
    Validate LQG controller performance and robustness.
    
    This class provides tools for testing controller performance under
    various conditions and validating against specifications.
    """
    
    def __init__(self):
        """Initialize controller validator."""
        self.validation_results = {}
    
    def validate_tracking_performance(self, params: LQGParameters,
                                    pulse_data: List[Dict[str, Any]],
                                    target_temperatures: Dict[str, float]) -> Dict[str, Any]:
        """
        Validate tracking performance on test data.
        
        Parameters
        ----------
        params : LQGParameters
            Controller parameters to validate
        pulse_data : List[Dict]
            Test data for validation
        target_temperatures : Dict[str, float]
            Target temperatures
            
        Returns
        -------
        results : Dict[str, Any]
            Validation results
        """
        # Create controller components
        model = LiftedStateSpaceModel(params)
        kf = KalmanFilter(model)
        lqr = LQRController(model, params)
        
        # Reference trajectory
        reference = np.array([
            target_temperatures.get('structure', 300.0),
            target_temperatures.get('kinetic', 300.0),
            target_temperatures.get('vibrational', 300.0)
        ])
        
        # Simulation results
        times = []
        states = []
        estimates = []
        controls = []
        tracking_errors = []
        
        for i, data_point in enumerate(pulse_data):
            time_ps = data_point.get('time_ps', i * 5.0)
            state = np.array(data_point['state'])
            
            if len(state) == 3:
                state = np.append(state, 0.0)  # Add hidden state
            
            # Kalman filter step
            y_k = state[:3]
            u_prev = np.array([0.0, 300.0]) if i == 0 else u
            x_hat = kf.step(y_k, u_prev)
            
            # Control computation
            u = lqr.compute_control(x_hat, reference)
            
            # Store results
            times.append(time_ps)
            states.append(state.copy())
            estimates.append(x_hat.copy())
            controls.append(u.copy())
            tracking_errors.append(np.linalg.norm(y_k - reference))
        
        # Compute performance metrics
        tracking_errors = np.array(tracking_errors)
        controls = np.array(controls)
        
        results = {
            'mean_tracking_error': np.mean(tracking_errors),
            'max_tracking_error': np.max(tracking_errors),
            'rms_tracking_error': np.sqrt(np.mean(tracking_errors**2)),
            'mean_control_effort': np.mean(np.linalg.norm(controls, axis=1)),
            'max_coupling_deviation': np.max(np.abs(controls[:, 0])),
            'bath_temp_range': [np.min(controls[:, 1]), np.max(controls[:, 1])],
            'settling_time_ps': self._estimate_settling_time(tracking_errors, times),
            'steady_state_error': np.mean(tracking_errors[-10:]) if len(tracking_errors) >= 10 else np.nan,
            'times': times,
            'states': states,
            'estimates': estimates,
            'controls': controls,
            'tracking_errors': tracking_errors.tolist()
        }
        
        return results
    
    def _estimate_settling_time(self, errors: np.ndarray, times: List[float], 
                               tolerance: float = 0.05) -> float:
        """Estimate settling time (time to reach within tolerance of steady state)."""
        if len(errors) < 10:
            return np.nan
        
        steady_state_error = np.mean(errors[-10:])
        threshold = steady_state_error + tolerance * np.max(errors)
        
        # Find last time error exceeded threshold
        for i in range(len(errors) - 1, -1, -1):
            if errors[i] > threshold:
                return times[i] if i < len(times) else np.nan
        
        return times[0] if times else 0.0
    
    def validate_robustness(self, params: LQGParameters,
                          disturbance_levels: List[float] = [0.1, 0.2, 0.5]) -> Dict[str, Any]:
        """
        Validate controller robustness to model uncertainties and disturbances.
        
        Parameters
        ----------
        params : LQGParameters
            Controller parameters
        disturbance_levels : List[float]
            Levels of model uncertainty to test
            
        Returns
        -------
        robustness_results : Dict[str, Any]
            Robustness analysis results
        """
        nominal_model = LiftedStateSpaceModel(params)
        
        robustness_results = {
            'nominal_stability': self._check_stability(nominal_model),
            'disturbance_tests': []
        }
        
        for disturbance_level in disturbance_levels:
            # Create perturbed model
            perturbed_params = self._perturb_parameters(params, disturbance_level)
            perturbed_model = LiftedStateSpaceModel(perturbed_params)
            
            # Test stability
            is_stable = self._check_stability(perturbed_model)
            
            # Test performance degradation
            performance_degradation = self._assess_performance_degradation(
                nominal_model, perturbed_model, params
            )
            
            robustness_results['disturbance_tests'].append({
                'disturbance_level': disturbance_level,
                'stable': is_stable,
                'performance_degradation': performance_degradation
            })
        
        return robustness_results
    
    def _check_stability(self, model: LiftedStateSpaceModel) -> bool:
        """Check if model is stable."""
        eigenvals = np.linalg.eigvals(model.A)
        max_eigenval = np.max(np.abs(eigenvals))
        return max_eigenval < 1.0
    
    def _perturb_parameters(self, params: LQGParameters, level: float) -> LQGParameters:
        """Create perturbed parameters for robustness testing."""
        perturbed = LQGParameters(**params.to_dict())
        
        # Add random perturbations to physical parameters
        perturbed.a_sk *= (1 + level * np.random.normal(0, 0.1))
        perturbed.a_sv *= (1 + level * np.random.normal(0, 0.1))
        perturbed.a_kv *= (1 + level * np.random.normal(0, 0.1))
        perturbed.b_sg *= (1 + level * np.random.normal(0, 0.1))
        perturbed.b_kg *= (1 + level * np.random.normal(0, 0.1))
        perturbed.b_vg *= (1 + level * np.random.normal(0, 0.1))
        
        # Ensure parameters stay within reasonable bounds
        perturbed.a_sk = np.clip(perturbed.a_sk, 0.01, 0.5)
        perturbed.a_sv = np.clip(perturbed.a_sv, 0.01, 0.3)
        perturbed.a_kv = np.clip(perturbed.a_kv, 0.01, 0.3)
        
        return perturbed
    
    def _assess_performance_degradation(self, nominal_model: LiftedStateSpaceModel,
                                      perturbed_model: LiftedStateSpaceModel,
                                      params: LQGParameters) -> float:
        """Assess performance degradation due to model mismatch."""
        # Simple metric: compare step responses
        x0 = np.array([310.0, 290.0, 305.0, 0.0])[:nominal_model.n_states]
        u = np.array([0.01, 300.0])
        
        # Simulate both models
        x_nom, _ = nominal_model.simulate_step(x0, u)
        x_pert, _ = perturbed_model.simulate_step(x0, u)
        
        # Compute relative error
        error = np.linalg.norm(x_nom[:3] - x_pert[:3]) / np.linalg.norm(x_nom[:3])
        
        return error
    
    def generate_validation_report(self, results: Dict[str, Any], output_file: str) -> None:
        """Generate validation report."""
        with open(output_file, 'w') as f:
            f.write("LQG Controller Validation Report\n")
            f.write("=" * 40 + "\n\n")
            
            if 'tracking_performance' in results:
                tp = results['tracking_performance']
                f.write("Tracking Performance:\n")
                f.write(f"  Mean tracking error: {tp['mean_tracking_error']:.3f} K\n")
                f.write(f"  RMS tracking error: {tp['rms_tracking_error']:.3f} K\n")
                f.write(f"  Max tracking error: {tp['max_tracking_error']:.3f} K\n")
                f.write(f"  Settling time: {tp['settling_time_ps']:.1f} ps\n")
                f.write(f"  Steady-state error: {tp['steady_state_error']:.3f} K\n\n")
            
            if 'robustness' in results:
                rob = results['robustness']
                f.write("Robustness Analysis:\n")
                f.write(f"  Nominal stability: {rob['nominal_stability']}\n")
                
                for test in rob['disturbance_tests']:
                    f.write(f"  Disturbance {test['disturbance_level']:.1%}: ")
                    f.write(f"Stable={test['stable']}, ")
                    f.write(f"Degradation={test['performance_degradation']:.3f}\n")
        
        print(f"Validation report saved to {output_file}")


def main():
    """Main function for command-line interface."""
    parser = argparse.ArgumentParser(description="LQG Controller Calibration Tools")
    parser.add_argument('--mode', choices=['collect', 'identify', 'optimize', 'validate'], 
                       required=True, help='Operation mode')
    parser.add_argument('--input', help='Input file path')
    parser.add_argument('--output', help='Output file path')
    parser.add_argument('--params', help='LQG parameters file')
    parser.add_argument('--data', help='Pulse data file')
    parser.add_argument('--period', type=float, default=5.0, help='Pulse period (ps)')
    parser.add_argument('--duty-cycle', type=float, default=0.5, help='Duty cycle')
    parser.add_argument('--target-temp', type=float, default=300.0, help='Target temperature (K)')
    
    args = parser.parse_args()
    
    if args.mode == 'collect':
        # Data collection mode
        if not args.input:
            print("Error: --input required for collect mode")
            return
        
        collector = DataCollector(args.period, args.duty_cycle)
        
        if args.input.endswith('.csv'):
            pulse_data = collector.collect_from_csv(args.input)
        else:
            print("Error: Only CSV input supported for collect mode")
            return
        
        output_file = args.output or 'pulse_data.json'
        collector.save_pulse_data(pulse_data, output_file)
    
    elif args.mode == 'identify':
        # System identification mode
        if not args.data:
            print("Error: --data required for identify mode")
            return
        
        collector = DataCollector()
        pulse_data = collector.load_pulse_data(args.data)
        
        sys_id = SystemIdentification()
        fitted_params = sys_id.fit_from_data(pulse_data)
        
        output_file = args.output or 'identified_params.json'
        fitted_params.save(output_file)
        print(f"Identified parameters saved to {output_file}")
    
    elif args.mode == 'optimize':
        # Parameter optimization mode
        if not args.data:
            print("Error: --data required for optimize mode")
            return
        
        collector = DataCollector()
        pulse_data = collector.load_pulse_data(args.data)
        
        target_temps = {
            'structure': args.target_temp,
            'kinetic': args.target_temp,
            'vibrational': args.target_temp
        }
        
        optimizer = ParameterOptimizer()
        optimized_params = optimizer.optimize_for_tracking(pulse_data, target_temps)
        
        output_file = args.output or 'optimized_params.json'
        optimized_params.save(output_file)
        print(f"Optimized parameters saved to {output_file}")
    
    elif args.mode == 'validate':
        # Validation mode
        if not args.params or not args.data:
            print("Error: --params and --data required for validate mode")
            return
        
        # Load parameters and data
        params = LQGParameters.load(args.params)
        collector = DataCollector()
        pulse_data = collector.load_pulse_data(args.data)
        
        target_temps = {
            'structure': args.target_temp,
            'kinetic': args.target_temp,
            'vibrational': args.target_temp
        }
        
        # Run validation
        validator = ControllerValidator()
        
        tracking_results = validator.validate_tracking_performance(
            params, pulse_data, target_temps
        )
        
        robustness_results = validator.validate_robustness(params)
        
        # Generate report
        results = {
            'tracking_performance': tracking_results,
            'robustness': robustness_results
        }
        
        output_file = args.output or 'validation_report.txt'
        validator.generate_validation_report(results, output_file)
        
        # Save detailed results
        json_output = output_file.replace('.txt', '_detailed.json')
        with open(json_output, 'w') as f:
            # Convert numpy arrays to lists for JSON serialization
            json_results = {}
            for key, value in results.items():
                if isinstance(value, dict):
                    json_results[key] = {}
                    for k, v in value.items():
                        if isinstance(v, np.ndarray):
                            json_results[key][k] = v.tolist()
                        elif isinstance(v, list) and len(v) > 0 and isinstance(v[0], np.ndarray):
                            json_results[key][k] = [arr.tolist() for arr in v]
                        else:
                            json_results[key][k] = v
                else:
                    json_results[key] = value
            
            json.dump(json_results, f, indent=2)
        
        print(f"Detailed results saved to {json_output}")


if __name__ == '__main__':
    main()
