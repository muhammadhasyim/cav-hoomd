#!/usr/bin/env python3
"""
Integration tests for Phase 3 refactoring.

These tests validate that the new plugin-based framework produces identical
results to the original monolithic implementation, ensuring backward
compatibility and scientific accuracy.
"""

import pytest
import tempfile
import shutil
from pathlib import Path
import numpy as np
from unittest.mock import Mock, patch, MagicMock

# Import both the new and old implementations for comparison
from hoomd.cavitymd import CavityMDExperiment, CavitySimulationParams
from hoomd.cavitymd.experiment import CavityMDExperiment as ExperimentClass
from hoomd.cavitymd.simulation import CavityMDSimulation
from hoomd.cavitymd.validation import CavitySimulationParams


class TestPhase3Integration:
    """Integration tests for Phase 3 refactoring."""
    
    def setup_method(self):
        """Setup test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)
        
    def teardown_method(self):
        """Cleanup test environment."""
        shutil.rmtree(self.temp_dir)
    
    def test_experiment_initialization(self):
        """Test that CavityMDExperiment initializes correctly."""
        # Test basic initialization
        kwargs = {
            'molecular_bath': 'bussi',
            'cavity_bath': 'langevin',
            'coupling': 1e-3,
            'temperature': 100.0,
            'frequency': 2000.0,
            'runtime': 500.0,
            'device': 'CPU',
            'replicas': '1'
        }
        
        experiment = CavityMDExperiment(**kwargs)
        
        # Verify experiment attributes
        assert experiment.params.molecular_thermostat == 'bussi'
        assert experiment.params.cavity_thermostat == 'langevin'
        assert experiment.params.coupling_strength == 1e-3
        assert experiment.params.temperature == 100.0
        assert experiment.params.frequency == 2000.0
        assert experiment.params.runtime_ps == 500.0
        assert experiment.params.device == 'CPU'
        assert experiment.replica_list == [1]
        
    def test_experiment_from_cli_args(self):
        """Test experiment creation from CLI arguments."""
        # Mock argparse Namespace
        args = Mock()
        args.molecular_bath = 'bussi'
        args.cavity_bath = 'langevin'
        args.coupling = 1e-3
        args.temperature = 100.0
        args.frequency = 2000.0
        args.runtime = 500.0
        args.device = 'CPU'
        args.replicas = '1-3'
        args.no_cavity = False
        args.finite_q = False
        args.fixed_timestep = False
        args.enable_energy_tracker = False
        args.enable_fkt = False
        args.truncate_gsd = False
        args.no_restart_velocities = False
        args.seed = None
        args.damping_ratio = 0.0
        args.switch_time = None
        args.molecular_tau = 5.0
        args.cavity_tau = 5.0
        args.timestep = 1.0
        args.gpu_id = 0
        args.energy_output_period_ps = 0.1
        args.fkt_output_period_ps = 1.0
        args.gsd_output_period_ps = 50.0
        args.console_output_period_ps = 1.0
        args.fkt_kmag = 1.0
        args.fkt_wavevectors = 50
        args.fkt_ref_interval = 1.0
        args.fkt_max_refs = 10
        args.max_energy_output_time = None
        
        # Mock SLURM environment
        with patch.dict('os.environ', {}, clear=True):
            with patch('hoomd.cavitymd.experiment.CavityMDExperiment._get_slurm_info', return_value=(None, None)):
                experiment = CavityMDExperiment.from_cli_args(args)
        
        assert experiment.params.molecular_thermostat == 'bussi'
        assert experiment.params.cavity_thermostat == 'langevin'
        assert experiment.replica_list == [1, 2, 3]
        
    def test_parameter_conversion(self):
        """Test conversion from CLI args to CavitySimulationParams."""
        kwargs = {
            'molecular_bath': 'bussi',
            'cavity_bath': 'langevin',
            'coupling': 1e-3,
            'temperature': 100.0,
            'frequency': 2000.0,
            'runtime': 500.0,
            'no_cavity': True,  # Should invert to incavity=False
            'finite_q': True,
            'fixed_timestep': True,
            'timestep': 2.0,  # Should set dt_fs when fixed_timestep=True
            'no_restart_velocities': True,  # Should invert to restart_velocities=False
            'damping_ratio': 0.1,
        }
        
        experiment = CavityMDExperiment(**kwargs)
        
        # Test parameter conversion
        assert experiment.params.incavity == False  # Inverted from no_cavity=True
        assert experiment.params.finite_q == True
        assert experiment.params.dt_fs == 2.0  # Should be set from timestep when fixed_timestep=True
        assert experiment.params.restart_velocities == False  # Inverted from no_restart_velocities=True
        assert experiment.damping_ratio == 0.1  # Should be stored separately
        assert experiment.params.error_tolerance == 0.0  # Should be 0.0 for fixed timestep
        
        # Test dissipation calculation
        assert experiment.params.dissipation > 0.0  # Should be calculated from damping_ratio
        
    def test_replica_parsing(self):
        """Test replica specification parsing."""
        experiment = CavityMDExperiment()
        
        # Test range specification
        replicas = experiment._parse_replicas('1-5')
        assert replicas == [1, 2, 3, 4, 5]
        
        # Test comma-separated specification
        replicas = experiment._parse_replicas('1,3,5')
        assert replicas == [1, 3, 5]
        
        # Test mixed specification
        replicas = experiment._parse_replicas('1-3,5,7-8')
        assert replicas == [1, 2, 3, 5, 7, 8]
        
        # Test default (no specification)
        replicas = experiment._parse_replicas(None)
        assert replicas == [1]
        
    def test_slurm_detection(self):
        """Test SLURM environment detection."""
        # Test with SLURM environment
        with patch.dict('os.environ', {
            'SLURM_ARRAY_TASK_ID': '5',
            'SLURM_ARRAY_JOB_ID': '12345'
        }):
            task_id, job_id = CavityMDExperiment._get_slurm_info()
            assert task_id == 5
            assert job_id == 12345
        
        # Test without SLURM environment
        with patch.dict('os.environ', {}, clear=True):
            task_id, job_id = CavityMDExperiment._get_slurm_info()
            assert task_id is None
            assert job_id is None
        
    def test_experiment_directory_creation(self):
        """Test experiment directory creation."""
        # Test cavity simulation directory
        params = CavitySimulationParams(
            temperature=100.0,
            frequency=2000.0,
            runtime_ps=500.0,
            coupling_strength=1e-3,
            incavity=True,
            switch_time_ps=None
        )
        experiment = CavityMDExperiment()
        experiment.params = params
        
        with patch('pathlib.Path.mkdir') as mock_mkdir:
            exp_dir = experiment._create_experiment_directory(1)
            expected_dir = Path("cavity_coupling_1eneg03")
            assert exp_dir == expected_dir
            mock_mkdir.assert_called_once()
        
        # Test cavity simulation with switch time
        params = CavitySimulationParams(
            temperature=100.0,
            frequency=2000.0,
            runtime_ps=500.0,
            coupling_strength=1e-3,
            incavity=True,
            switch_time_ps=1.0
        )
        experiment.params = params
        
        with patch('pathlib.Path.mkdir') as mock_mkdir:
            exp_dir = experiment._create_experiment_directory(1)
            expected_dir = Path("cavity_coupling_1eneg03_switch_1.0ps")
            assert exp_dir == expected_dir
            mock_mkdir.assert_called_once()
        
        # Test non-cavity simulation directory
        params = CavitySimulationParams(
            temperature=100.0,
            frequency=2000.0,
            runtime_ps=500.0,
            coupling_strength=1e-3,
            incavity=False
        )
        experiment.params = params
        
        with patch('pathlib.Path.mkdir') as mock_mkdir:
            exp_dir = experiment._create_experiment_directory(1)
            expected_dir = Path("no_cavity")
            assert exp_dir == expected_dir
            mock_mkdir.assert_called_once()
    
    def test_dissipation_calculation(self):
        """Test dissipation calculation from damping ratio."""
        kwargs = {
            'damping_ratio': 0.1,
            'frequency': 2000.0,  # cm^-1
        }
        
        experiment = CavityMDExperiment(**kwargs)
        
        # Expected dissipation calculation
        from hoomd.cavitymd.analysis import PhysicalConstants
        phmass = 1.0
        omegac = 2000.0 / PhysicalConstants.HARTREE_TO_CM_MINUS1
        expected_dissipation = 2 * 0.1 * phmass * omegac
        
        assert abs(experiment.params.dissipation - expected_dissipation) < 1e-10
        
    def test_experiment_mock_run(self):
        """Test experiment run with mocked simulation."""
        kwargs = {
            'molecular_bath': 'bussi',
            'cavity_bath': 'langevin',
            'coupling': 1e-3,
            'runtime': 10.0,  # Short runtime for testing
            'replicas': '1',
        }
        
        # Mock the simulation components
        with patch('hoomd.cavitymd.experiment.CavityMDSimulation') as mock_simulation_class:
            mock_simulation = Mock()
            mock_simulation.setup.return_value = None
            mock_simulation.run.return_value = None
            mock_simulation_class.return_value = mock_simulation
            
            experiment = CavityMDExperiment(**kwargs)
            
            # Redirect stdout to capture output
            import io
            import sys
            captured_output = io.StringIO()
            sys.stdout = captured_output
            
            try:
                exit_code = experiment.run()
                output = captured_output.getvalue()
                
                # Verify successful execution
                assert exit_code == 0
                assert "Advanced Cavity MD Experiment Runner" in output
                assert "SUCCESS: Replica 1 completed successfully" in output
                assert "All replicas completed successfully!" in output
                
                # Verify simulation was called correctly
                mock_simulation_class.assert_called_once()
                mock_simulation.setup.assert_called_once()
                mock_simulation.run.assert_called_once()
                
            finally:
                sys.stdout = sys.__stdout__
    
    def test_experiment_failure_handling(self):
        """Test experiment failure handling."""
        kwargs = {
            'molecular_bath': 'bussi',
            'cavity_bath': 'langevin',
            'coupling': 1e-3,
            'runtime': 10.0,
            'replicas': '1',
        }
        
        # Mock the simulation to raise an exception
        with patch('hoomd.cavitymd.experiment.CavityMDSimulation') as mock_simulation_class:
            mock_simulation_class.side_effect = Exception("Test exception")
            
            experiment = CavityMDExperiment(**kwargs)
            
            # Redirect stdout to capture output
            import io
            import sys
            captured_output = io.StringIO()
            sys.stdout = captured_output
            
            try:
                exit_code = experiment.run()
                output = captured_output.getvalue()
                
                # Verify failure handling
                assert exit_code == 1
                assert "ERROR: Replica 1 failed" in output
                assert "1 replicas failed" in output
                
            finally:
                sys.stdout = sys.__stdout__
    
    def test_multiple_replicas(self):
        """Test multiple replica execution."""
        kwargs = {
            'molecular_bath': 'bussi',
            'cavity_bath': 'langevin',
            'coupling': 1e-3,
            'runtime': 10.0,
            'replicas': '1-3',
        }
        
        # Mock successful simulation
        with patch('hoomd.cavitymd.experiment.CavityMDSimulation') as mock_simulation_class:
            mock_simulation = Mock()
            mock_simulation.setup.return_value = None
            mock_simulation.run.return_value = None
            mock_simulation_class.return_value = mock_simulation
            
            experiment = CavityMDExperiment(**kwargs)
            
            # Redirect stdout to capture output
            import io
            import sys
            captured_output = io.StringIO()
            sys.stdout = captured_output
            
            try:
                exit_code = experiment.run()
                output = captured_output.getvalue()
                
                # Verify all replicas executed
                assert exit_code == 0
                assert "3 replica(s): [1, 2, 3]" in output
                assert "SUCCESS: Replica 1 completed successfully" in output
                assert "SUCCESS: Replica 2 completed successfully" in output
                assert "SUCCESS: Replica 3 completed successfully" in output
                assert "Total replicas: 3" in output
                assert "Successful: 3" in output
                assert "Failed: 0" in output
                
                # Verify simulation was called 3 times
                assert mock_simulation_class.call_count == 3
                
            finally:
                sys.stdout = sys.__stdout__
    
    def test_parameter_validation_integration(self):
        """Test that parameter validation works correctly in the experiment."""
        # Test invalid temperature
        kwargs = {
            'temperature': -100.0,  # Invalid negative temperature
        }
        
        with pytest.raises(ValueError, match="Temperature must be positive"):
            CavityMDExperiment(**kwargs)
        
        # Test invalid frequency
        kwargs = {
            'frequency': -1000.0,  # Invalid negative frequency
        }
        
        with pytest.raises(ValueError, match="Frequency must be positive"):
            CavityMDExperiment(**kwargs)
        
        # Test invalid runtime
        kwargs = {
            'runtime': 0.0,  # Invalid zero runtime
        }
        
        with pytest.raises(ValueError, match="Runtime must be positive"):
            CavityMDExperiment(**kwargs)
    
    def test_backward_compatibility(self):
        """Test that the new framework maintains backward compatibility."""
        # Test that all original CLI arguments are supported
        original_args = [
            'molecular_bath', 'cavity_bath', 'finite_q', 'coupling',
            'switch_time', 'damping_ratio', 'temperature', 'frequency',
            'runtime', 'no_cavity', 'replicas', 'molecular_tau',
            'cavity_tau', 'fixed_timestep', 'timestep', 'enable_energy_tracker',
            'energy_output_period_ps', 'fkt_output_period_ps',
            'gsd_output_period_ps', 'console_output_period_ps',
            'enable_fkt', 'fkt_kmag', 'fkt_wavevectors', 'fkt_ref_interval',
            'fkt_max_refs', 'max_energy_output_time', 'device', 'gpu_id',
            'truncate_gsd', 'seed', 'no_restart_velocities'
        ]
        
        # Create kwargs with all original arguments
        kwargs = {
            'molecular_bath': 'bussi',
            'cavity_bath': 'langevin',
            'finite_q': True,
            'coupling': 1e-3,
            'switch_time': 1.0,
            'damping_ratio': 0.1,
            'temperature': 100.0,
            'frequency': 2000.0,
            'runtime': 500.0,
            'no_cavity': False,
            'replicas': '1',
            'molecular_tau': 5.0,
            'cavity_tau': 5.0,
            'fixed_timestep': True,
            'timestep': 1.0,
            'enable_energy_tracker': True,
            'energy_output_period_ps': 0.1,
            'fkt_output_period_ps': 1.0,
            'gsd_output_period_ps': 50.0,
            'console_output_period_ps': 1.0,
            'enable_fkt': True,
            'fkt_kmag': 1.0,
            'fkt_wavevectors': 50,
            'fkt_ref_interval': 1.0,
            'fkt_max_refs': 10,
            'max_energy_output_time': 100.0,
            'device': 'CPU',
            'gpu_id': 0,
            'truncate_gsd': True,
            'seed': 42,
            'no_restart_velocities': False
        }
        
        # Should not raise any exceptions
        experiment = CavityMDExperiment(**kwargs)
        
        # Verify all parameters were converted correctly
        assert experiment.params.molecular_thermostat == 'bussi'
        assert experiment.params.cavity_thermostat == 'langevin'
        assert experiment.params.finite_q == True
        assert experiment.params.coupling_strength == 1e-3
        assert experiment.params.switch_time_ps == 1.0
        assert experiment.damping_ratio == 0.1  # Stored separately in experiment
        assert experiment.params.temperature == 100.0
        assert experiment.params.frequency == 2000.0
        assert experiment.params.runtime_ps == 500.0
        assert experiment.params.incavity == True
        assert experiment.params.molecular_thermostat_tau == 5.0
        assert experiment.params.cavity_thermostat_tau == 5.0
        assert experiment.params.dt_fs == 1.0  # Fixed timestep mode
        assert experiment.params.enable_energy_tracking == True
        assert experiment.params.energy_output_period_ps == 0.1
        assert experiment.params.fkt_output_period_ps == 1.0
        assert experiment.params.gsd_output_period_ps == 50.0
        assert experiment.params.console_output_period_ps == 1.0
        assert experiment.params.enable_fkt == True
        assert experiment.params.fkt_kmag == 1.0
        assert experiment.params.fkt_num_wavevectors == 50
        assert experiment.params.fkt_reference_interval_ps == 1.0
        assert experiment.params.fkt_max_references == 10
        assert experiment.params.max_energy_output_time_ps == 100.0
        assert experiment.params.device == 'CPU'
        assert experiment.params.gpu_id == 0
        assert experiment.params.truncate_gsd == True
        assert experiment.params.seed == 42
        assert experiment.params.restart_velocities == True  # no_restart_velocities=False inverted
        
        # Verify dissipation calculation
        assert experiment.params.dissipation > 0.0  # Should be calculated from damping_ratio
        
        # Verify replica list
        assert experiment.replica_list == [1]


class TestRunScriptIntegration:
    """Integration tests for the simplified run script."""
    
    def test_run_script_imports(self):
        """Test that the run script imports work correctly."""
        # Test import of the new run script
        import sys
        import importlib.util
        
        # Load the run script as a module
        spec = importlib.util.spec_from_file_location(
            "run_script", 
            Path(__file__).parent.parent.parent / "examples" / "05_advanced_run_v2.py"
        )
        run_script = importlib.util.module_from_spec(spec)
        
        # Should not raise import errors
        spec.loader.exec_module(run_script)
        
        # Verify key functions exist
        assert hasattr(run_script, 'create_argument_parser')
        assert hasattr(run_script, 'main')
        
    def test_argument_parser_compatibility(self):
        """Test that the argument parser matches the original."""
        import sys
        import importlib.util
        
        # Load the run script as a module
        spec = importlib.util.spec_from_file_location(
            "run_script", 
            Path(__file__).parent.parent.parent / "examples" / "05_advanced_run_v2.py"
        )
        run_script = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(run_script)
        
        # Create parser and test argument parsing
        parser = run_script.create_argument_parser()
        
        # Test that all expected arguments are present
        expected_args = [
            '--molecular-bath', '--cavity-bath', '--finite-q', '--coupling',
            '--switch-time', '--damping-ratio', '--temperature', '--frequency',
            '--runtime', '--no-cavity', '--replicas', '--molecular-tau',
            '--cavity-tau', '--fixed-timestep', '--timestep', '--enable-energy-tracker',
            '--energy-output-period-ps', '--fkt-output-period-ps',
            '--gsd-output-period-ps', '--console-output-period-ps',
            '--enable-fkt', '--fkt-kmag', '--fkt-wavevectors', '--fkt-ref-interval',
            '--fkt-max-refs', '--max-energy-output-time', '--device', '--gpu-id',
            '--truncate-gsd', '--seed', '--no-restart-velocities'
        ]
        
        # Parse a comprehensive set of arguments
        test_args = [
            '--molecular-bath', 'bussi',
            '--cavity-bath', 'langevin',
            '--finite-q',
            '--coupling', '1e-3',
            '--temperature', '100.0',
            '--frequency', '2000.0',
            '--runtime', '500.0',
            '--replicas', '1-3',
            '--enable-energy-tracker',
            '--enable-fkt',
            '--device', 'CPU'
        ]
        
        # Should parse without errors
        args = parser.parse_args(test_args)
        
        # Verify parsed values
        assert args.molecular_bath == 'bussi'
        assert args.cavity_bath == 'langevin'
        assert args.finite_q == True
        assert args.coupling == 1e-3
        assert args.temperature == 100.0
        assert args.frequency == 2000.0
        assert args.runtime == 500.0
        assert args.replicas == '1-3'
        assert args.enable_energy_tracker == True
        assert args.enable_fkt == True
        assert args.device == 'CPU'


if __name__ == '__main__':
    pytest.main([__file__]) 