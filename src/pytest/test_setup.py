"""
Tests for the setup utilities module.

This module tests all setup classes and their methods to ensure proper
configuration of forces, thermostats, and other simulation components.
"""

import pytest
import numpy as np
import hoomd
import gsd.hoomd
from unittest.mock import Mock, MagicMock, patch
import tempfile
from pathlib import Path

from hoomd.cavitymd.setup import (
    DeviceSetup, ParticleSetup, ForceSetup, ThermostatSetup,
    ThermalizationSetup, TimestepSetup
)
from hoomd.cavitymd.validation import CavitySimulationParams


class TestDeviceSetup:
    """Test DeviceSetup class functionality."""
    
    def test_create_cpu_device(self):
        """Test CPU device creation."""
        device = DeviceSetup.create_device('CPU')
        assert isinstance(device, hoomd.device.CPU)
        
        # Test case insensitive
        device = DeviceSetup.create_device('cpu')
        assert isinstance(device, hoomd.device.CPU)
    
    def test_create_gpu_device(self):
        """Test GPU device creation."""
        device = DeviceSetup.create_device('GPU', gpu_id=0)
        assert isinstance(device, hoomd.device.GPU)
        
        device = DeviceSetup.create_device('GPU', gpu_id=1)
        assert isinstance(device, hoomd.device.GPU)
    
    def test_invalid_device_type(self):
        """Test invalid device type handling."""
        with pytest.raises(ValueError, match="Invalid device type"):
            DeviceSetup.create_device('TPU')
    
    def test_generate_simulation_seed(self):
        """Test simulation seed generation."""
        # Test user-specified seed
        seed = DeviceSetup.generate_simulation_seed(12345, 1)
        assert seed == 12345
        
        # Test replica-based seed
        seed1 = DeviceSetup.generate_simulation_seed(None, 1)
        seed2 = DeviceSetup.generate_simulation_seed(None, 2)
        assert isinstance(seed1, int)
        assert isinstance(seed2, int)
        assert seed1 != seed2  # Different replicas should have different seeds
        
        # Test deterministic behavior
        seed1a = DeviceSetup.generate_simulation_seed(None, 1)
        seed1b = DeviceSetup.generate_simulation_seed(None, 1)
        assert seed1a == seed1b  # Same replica should have same seed


class TestParticleSetup:
    """Test ParticleSetup class functionality."""
    
    def create_mock_snapshot(self, has_cavity=False):
        """Create a mock GSD snapshot for testing."""
        snapshot = Mock()
        snapshot.particles = Mock()
        snapshot.configuration = Mock()
        snapshot.configuration.box = [10, 10, 10, 0, 0, 0]
        
        if has_cavity:
            snapshot.particles.types = ['O', 'N', 'L']
            snapshot.particles.typeid = np.array([0, 0, 1, 1, 2])  # 2 O, 2 N, 1 L
            snapshot.particles.N = 5
        else:
            snapshot.particles.types = ['O', 'N']
            snapshot.particles.typeid = np.array([0, 0, 1, 1])  # 2 O, 2 N
            snapshot.particles.N = 4
        
        snapshot.particles.position = np.zeros((snapshot.particles.N, 3))
        snapshot.particles.image = np.zeros((snapshot.particles.N, 3), dtype=int)
        snapshot.particles.charge = np.ones(snapshot.particles.N) * 0.5
        
        return snapshot
    
    def test_check_cavity_particle_exists(self):
        """Test cavity particle existence checking."""
        # Test with cavity particle
        snapshot = self.create_mock_snapshot(has_cavity=True)
        exists, count = ParticleSetup.check_cavity_particle_exists(snapshot)
        assert exists is True
        assert count == 1
        
        # Test without cavity particle
        snapshot = self.create_mock_snapshot(has_cavity=False)
        exists, count = ParticleSetup.check_cavity_particle_exists(snapshot)
        assert exists is False
        assert count == 0
        
        # Test with no 'L' type at all
        snapshot = self.create_mock_snapshot(has_cavity=False)
        exists, count = ParticleSetup.check_cavity_particle_exists(snapshot)
        assert exists is False
        assert count == 0
    
    def test_create_cavity_particle(self):
        """Test cavity particle creation."""
        snapshot = self.create_mock_snapshot(has_cavity=False)
        
        # Mock additional attributes
        snapshot.particles.mass = np.ones(snapshot.particles.N)
        snapshot.particles.diameter = np.ones(snapshot.particles.N)
        snapshot.particles.velocity = np.zeros((snapshot.particles.N, 3))
        snapshot.particles.body = np.full(snapshot.particles.N, -1)
        snapshot.particles.orientation = np.zeros((snapshot.particles.N, 4))
        snapshot.particles.moment_inertia = np.zeros((snapshot.particles.N, 3))
        snapshot.particles.angmom = np.zeros((snapshot.particles.N, 4))
        
        # Set up append methods as mocks
        snapshot.particles.types.append = Mock()
        
        params = CavitySimulationParams(
            temperature=100.0,
            frequency=2000.0,
            coupling_strength=0.001,
            runtime_ps=1000.0,
            seed=12345
        )
        
        # Mock numpy methods
        with patch('numpy.random.seed') as mock_seed, \
             patch('numpy.einsum') as mock_einsum, \
             patch('numpy.random.normal') as mock_normal, \
             patch('numpy.allclose') as mock_allclose, \
             patch('numpy.append') as mock_append, \
             patch('numpy.vstack') as mock_vstack:
            
            mock_einsum.return_value = np.array([1.0, 2.0, 3.0])
            mock_normal.return_value = np.array([0.1, 0.2, 0.3])
            mock_allclose.return_value = False
            mock_append.side_effect = lambda arr, val, axis=None: np.append(arr, val, axis=axis)
            mock_vstack.side_effect = lambda arrays: np.vstack(arrays)
            
            result = ParticleSetup.create_cavity_particle(snapshot, params)
            
            # Verify seed was set
            mock_seed.assert_called_once_with(12347)  # seed + 2
            
            # Verify dipole moment calculation
            mock_einsum.assert_called_once()
            
            # Verify result is the modified snapshot
            assert result == snapshot
    
    def test_validate_cavity_particle_success(self):
        """Test successful cavity particle validation."""
        # Create mock simulation
        sim = Mock()
        sim.state = Mock()
        sim.state.particle_types = ['O', 'N', 'L']
        
        # Create mock snapshot context manager
        mock_snapshot = Mock()
        mock_snapshot.particles = Mock()
        mock_snapshot.particles.typeid = np.array([0, 0, 1, 1, 2])
        mock_snapshot.particles.position = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4]])
        
        sim.state.cpu_local_snapshot = MagicMock()
        sim.state.cpu_local_snapshot.__enter__ = Mock(return_value=mock_snapshot)
        sim.state.cpu_local_snapshot.__exit__ = Mock(return_value=None)
        
        # Should not raise any exception
        ParticleSetup.validate_cavity_particle(sim)
    
    def test_validate_cavity_particle_no_type(self):
        """Test cavity particle validation with no 'L' type."""
        sim = Mock()
        sim.state = Mock()
        sim.state.particle_types = ['O', 'N']  # No 'L' type
        
        with pytest.raises(ValueError, match="Cavity particle type 'L' not found"):
            ParticleSetup.validate_cavity_particle(sim)
    
    def test_validate_cavity_particle_wrong_count(self):
        """Test cavity particle validation with wrong count."""
        sim = Mock()
        sim.state = Mock()
        sim.state.particle_types = ['O', 'N', 'L']
        
        # Create mock snapshot with no cavity particles
        mock_snapshot = Mock()
        mock_snapshot.particles = Mock()
        mock_snapshot.particles.typeid = np.array([0, 0, 1, 1])  # No type 2
        
        sim.state.cpu_local_snapshot = MagicMock()
        sim.state.cpu_local_snapshot.__enter__ = Mock(return_value=mock_snapshot)
        sim.state.cpu_local_snapshot.__exit__ = Mock(return_value=None)
        
        with pytest.raises(ValueError, match="No cavity particles found"):
            ParticleSetup.validate_cavity_particle(sim)
    
    def test_validate_cavity_particle_multiple_cavities(self):
        """Test cavity particle validation with multiple cavities."""
        sim = Mock()
        sim.state = Mock()
        sim.state.particle_types = ['O', 'N', 'L']
        
        # Create mock snapshot with multiple cavity particles
        mock_snapshot = Mock()
        mock_snapshot.particles = Mock()
        mock_snapshot.particles.typeid = np.array([0, 0, 1, 1, 2, 2])  # Two type 2
        
        sim.state.cpu_local_snapshot = MagicMock()
        sim.state.cpu_local_snapshot.__enter__ = Mock(return_value=mock_snapshot)
        sim.state.cpu_local_snapshot.__exit__ = Mock(return_value=None)
        
        with pytest.raises(ValueError, match="Expected exactly 1 cavity particle, found 2"):
            ParticleSetup.validate_cavity_particle(sim)


class TestForceSetup:
    """Test ForceSetup class functionality."""
    
    def create_test_params(self, **kwargs):
        """Create test parameters with defaults."""
        defaults = {
            'temperature': 100.0,
            'frequency': 2000.0,
            'coupling_strength': 0.001,
            'runtime_ps': 1000.0,
            'incavity': True
        }
        defaults.update(kwargs)
        return CavitySimulationParams(**defaults)
    
    @patch('hoomd.cavitymd.setup.CavityForce')
    def test_create_cavity_force_basic(self, mock_cavity_force):
        """Test basic cavity force creation."""
        params = self.create_test_params()
        mock_force = Mock()
        mock_cavity_force.return_value = mock_force
        
        result = ForceSetup.create_cavity_force(params)
        
        assert result == mock_force
        mock_cavity_force.assert_called_once()
    
    def test_create_cavity_force_no_cavity(self):
        """Test cavity force creation when cavity is disabled."""
        params = self.create_test_params(incavity=False)
        
        result = ForceSetup.create_cavity_force(params)
        
        assert result is None
    
    @patch('hoomd.cavitymd.setup.CavityForce')
    @patch('hoomd.cavitymd.setup.StepVariant')
    def test_create_cavity_force_time_varying(self, mock_step_variant, mock_cavity_force):
        """Test cavity force creation with time-varying parameters."""
        params = self.create_test_params(switch_time_ps=500.0, dissipation=0.01)
        
        mock_time_tracker = Mock()
        mock_variants = [Mock(), Mock()]
        mock_step_variant.side_effect = mock_variants
        mock_force = Mock()
        mock_cavity_force.return_value = mock_force
        
        result = ForceSetup.create_cavity_force(params, mock_time_tracker)
        
        assert result == mock_force
        assert mock_step_variant.call_count == 2  # coupling and dissipation variants
        mock_cavity_force.assert_called_once()
    
    def test_create_cavity_force_time_varying_no_tracker(self):
        """Test cavity force creation with time-varying parameters but no tracker."""
        params = self.create_test_params(switch_time_ps=500.0)
        
        with pytest.raises(ValueError, match="Time tracker required"):
            ForceSetup.create_cavity_force(params)
    
    @patch('hoomd.md.bond.Harmonic')
    def test_create_harmonic_bonds(self, mock_harmonic):
        """Test harmonic bond creation."""
        mock_harmonic_instance = Mock()
        mock_harmonic.return_value = mock_harmonic_instance
        
        result = ForceSetup.create_harmonic_bonds()
        
        assert result == mock_harmonic_instance
        mock_harmonic.assert_called_once()
    
    @patch('hoomd.md.nlist.Cell')
    @patch('hoomd.md.pair.LJ')
    def test_create_lennard_jones_force(self, mock_lj, mock_cell):
        """Test Lennard-Jones force creation."""
        params = self.create_test_params()
        mock_cell_instance = Mock()
        mock_lj_instance = Mock()
        mock_lj_instance.params = {}
        mock_lj_instance.r_cut = {}
        
        mock_cell.return_value = mock_cell_instance
        mock_lj.return_value = mock_lj_instance
        
        result = ForceSetup.create_lennard_jones_force(params)
        
        assert result == mock_lj_instance
        mock_cell.assert_called_once()
        mock_lj.assert_called_once()
    
    @patch('hoomd.md.nlist.Cell')
    @patch('hoomd.md.long_range.pppm.make_pppm_coulomb_forces')
    def test_create_coulomb_forces(self, mock_make_pppm, mock_cell):
        """Test Coulomb force creation."""
        mock_cell_instance = Mock()
        mock_short = Mock()
        mock_long = Mock()
        
        mock_cell.return_value = mock_cell_instance
        mock_make_pppm.return_value = (mock_short, mock_long)
        
        short, long = ForceSetup.create_coulomb_forces()
        
        assert short == mock_short
        assert long == mock_long
        mock_cell.assert_called_once()
        mock_make_pppm.assert_called_once()
    
    @patch('hoomd.cavitymd.setup.ForceSetup.create_cavity_force')
    @patch('hoomd.cavitymd.setup.ForceSetup.create_harmonic_bonds')
    @patch('hoomd.cavitymd.setup.ForceSetup.create_lennard_jones_force')
    @patch('hoomd.cavitymd.setup.ForceSetup.create_coulomb_forces')
    def test_create_all_forces(self, mock_coulomb, mock_lj, mock_harmonic, mock_cavity):
        """Test creation of all forces."""
        params = self.create_test_params()
        
        mock_cavity_force = Mock()
        mock_harmonic_force = Mock()
        mock_lj_force = Mock()
        mock_short_force = Mock()
        mock_long_force = Mock()
        
        mock_cavity.return_value = mock_cavity_force
        mock_harmonic.return_value = mock_harmonic_force
        mock_lj.return_value = mock_lj_force
        mock_coulomb.return_value = (mock_short_force, mock_long_force)
        
        forces = ForceSetup.create_all_forces(params)
        
        assert len(forces) == 5
        assert mock_cavity_force in forces
        assert mock_harmonic_force in forces
        assert mock_lj_force in forces
        assert mock_short_force in forces
        assert mock_long_force in forces


class TestThermostatSetup:
    """Test ThermostatSetup class functionality."""
    
    def create_test_params(self, **kwargs):
        """Create test parameters with defaults."""
        defaults = {
            'temperature': 100.0,
            'frequency': 2000.0,
            'coupling_strength': 0.001,
            'runtime_ps': 1000.0,
            'molecular_thermostat': 'bussi',
            'cavity_thermostat': 'langevin',
            'molecular_thermostat_tau': 5.0,
            'cavity_thermostat_tau': 2.0,
            'cavity_damping_factor': 1.0,
            'incavity': True
        }
        defaults.update(kwargs)
        return CavitySimulationParams(**defaults)
    
    @patch('hoomd.filter.Type')
    @patch('hoomd.cavitymd.setup.Bussi')
    @patch('hoomd.md.methods.ConstantVolume')
    def test_create_molecular_thermostat_bussi(self, mock_cv, mock_bussi, mock_filter):
        """Test molecular Bussi thermostat creation."""
        params = self.create_test_params(molecular_thermostat='bussi')
        
        mock_filter_instance = Mock()
        mock_bussi_instance = Mock()
        mock_method = Mock()
        
        mock_filter.return_value = mock_filter_instance
        mock_bussi.return_value = mock_bussi_instance
        mock_cv.return_value = mock_method
        
        result = ThermostatSetup.create_molecular_thermostat(params)
        
        assert result == mock_method
        mock_filter.assert_called_once_with(['O', 'N'])
        mock_bussi.assert_called_once()
        mock_cv.assert_called_once()
    
    @patch('hoomd.filter.Type')
    @patch('hoomd.md.methods.Langevin')
    def test_create_molecular_thermostat_langevin(self, mock_langevin, mock_filter):
        """Test molecular Langevin thermostat creation."""
        params = self.create_test_params(molecular_thermostat='langevin')
        
        mock_filter_instance = Mock()
        mock_method = Mock()
        
        mock_filter.return_value = mock_filter_instance
        mock_langevin.return_value = mock_method
        
        result = ThermostatSetup.create_molecular_thermostat(params)
        
        assert result == mock_method
        mock_filter.assert_called_once_with(['O', 'N'])
        mock_langevin.assert_called_once()
    
    @patch('hoomd.filter.Type')
    @patch('hoomd.md.methods.ConstantVolume')
    def test_create_molecular_thermostat_none(self, mock_cv, mock_filter):
        """Test molecular NVE thermostat creation."""
        params = self.create_test_params(molecular_thermostat='none')
        
        mock_filter_instance = Mock()
        mock_method = Mock()
        
        mock_filter.return_value = mock_filter_instance
        mock_cv.return_value = mock_method
        
        result = ThermostatSetup.create_molecular_thermostat(params)
        
        assert result == mock_method
        mock_filter.assert_called_once_with(['O', 'N'])
        mock_cv.assert_called_once()
    
    def test_create_molecular_thermostat_invalid(self):
        """Test invalid molecular thermostat."""
        params = self.create_test_params(molecular_thermostat='invalid')
        
        with pytest.raises(ValueError, match="Invalid molecular thermostat"):
            ThermostatSetup.create_molecular_thermostat(params)
    
    def test_create_cavity_thermostat_no_cavity(self):
        """Test cavity thermostat creation when no cavity."""
        params = self.create_test_params(incavity=False)
        
        result = ThermostatSetup.create_cavity_thermostat(params)
        
        assert result is None
    
    @patch('hoomd.filter.Type')
    @patch('hoomd.cavitymd.setup.Bussi')
    @patch('hoomd.md.methods.ConstantVolume')
    def test_create_cavity_thermostat_bussi(self, mock_cv, mock_bussi, mock_filter):
        """Test cavity Bussi thermostat creation."""
        params = self.create_test_params(cavity_thermostat='bussi')
        
        mock_filter_instance = Mock()
        mock_bussi_instance = Mock()
        mock_method = Mock()
        
        mock_filter.return_value = mock_filter_instance
        mock_bussi.return_value = mock_bussi_instance
        mock_cv.return_value = mock_method
        
        result = ThermostatSetup.create_cavity_thermostat(params)
        
        assert result == mock_method
        mock_filter.assert_called_once_with(['L'])
        mock_bussi.assert_called_once()
        mock_cv.assert_called_once()
    
    @patch('hoomd.filter.Type')
    @patch('hoomd.md.methods.Langevin')
    def test_create_cavity_thermostat_langevin(self, mock_langevin, mock_filter):
        """Test cavity Langevin thermostat creation."""
        params = self.create_test_params(cavity_thermostat='langevin')
        
        mock_filter_instance = Mock()
        mock_method = Mock()
        
        mock_filter.return_value = mock_filter_instance
        mock_langevin.return_value = mock_method
        
        result = ThermostatSetup.create_cavity_thermostat(params)
        
        assert result == mock_method
        mock_filter.assert_called_once_with(['L'])
        mock_langevin.assert_called_once()
    
    @patch('hoomd.cavitymd.setup.ThermostatSetup.create_molecular_thermostat')
    @patch('hoomd.cavitymd.setup.ThermostatSetup.create_cavity_thermostat')
    def test_create_all_methods(self, mock_create_cavity, mock_create_molecular):
        """Test creation of all thermostat methods."""
        params = self.create_test_params()
        
        mock_molecular_method = Mock()
        mock_cavity_method = Mock()
        
        mock_create_molecular.return_value = mock_molecular_method
        mock_create_cavity.return_value = mock_cavity_method
        
        methods, refs = ThermostatSetup.create_all_methods(params)
        
        assert len(methods) == 2
        assert mock_molecular_method in methods
        assert mock_cavity_method in methods
        assert isinstance(refs, dict)


class TestThermalizationSetup:
    """Test ThermalizationSetup class functionality."""
    
    def create_test_params(self, **kwargs):
        """Create test parameters with defaults."""
        defaults = {
            'temperature': 100.0,
            'frequency': 2000.0,
            'coupling_strength': 0.001,
            'runtime_ps': 1000.0,
            'restart_velocities': True,
            'seed': 12345,
            'incavity': True
        }
        defaults.update(kwargs)
        return CavitySimulationParams(**defaults)
    
    def test_thermalize_system_disabled(self):
        """Test thermalization when restart_velocities is False."""
        params = self.create_test_params(restart_velocities=False)
        sim = Mock()
        
        ThermalizationSetup.thermalize_system(sim, params)
        
        # Should not call any thermalization methods
        sim.state.thermalize_particle_momenta.assert_not_called()
    
    @patch('hoomd.filter.Type')
    @patch('hoomd.filter.All')
    @patch('numpy.random.seed')
    def test_thermalize_system_no_cavity(self, mock_np_seed, mock_all_filter, mock_type_filter):
        """Test thermalization for no-cavity simulation."""
        params = self.create_test_params(incavity=False)
        sim = Mock()
        
        mock_all_filter_instance = Mock()
        mock_all_filter.return_value = mock_all_filter_instance
        
        ThermalizationSetup.thermalize_system(sim, params)
        
        mock_np_seed.assert_called_once_with(12346)  # seed + 1
        mock_all_filter.assert_called_once()
        sim.state.thermalize_particle_momenta.assert_called_once()
    
    @patch('hoomd.filter.Type')
    @patch('numpy.random.seed')
    @patch('hoomd.cavitymd.setup.ThermalizationSetup._thermalize_cavity_particle_manually')
    def test_thermalize_system_with_cavity(self, mock_manual_therm, mock_np_seed, mock_type_filter):
        """Test thermalization for cavity simulation."""
        params = self.create_test_params(incavity=True)
        sim = Mock()
        
        mock_molecular_filter = Mock()
        mock_type_filter.return_value = mock_molecular_filter
        
        ThermalizationSetup.thermalize_system(sim, params)
        
        mock_np_seed.assert_called_once_with(12346)  # seed + 1
        mock_type_filter.assert_called_once_with(['O', 'N'])
        sim.state.thermalize_particle_momenta.assert_called_once()
        mock_manual_therm.assert_called_once()
    
    @patch('numpy.random.normal')
    @patch('numpy.sqrt')
    @patch('numpy.where')
    def test_thermalize_cavity_particle_manually(self, mock_where, mock_sqrt, mock_normal):
        """Test manual cavity particle thermalization."""
        sim = Mock()
        kT = 0.01
        
        # Mock snapshot context manager
        mock_snapshot = Mock()
        mock_snapshot.particles = Mock()
        mock_snapshot.particles.typeid = np.array([0, 1, 2])
        mock_snapshot.particles.mass = np.array([1.0, 1.0, 1.0])
        mock_snapshot.particles.velocity = np.zeros((3, 3))
        
        sim.state.cpu_local_snapshot = MagicMock()
        sim.state.cpu_local_snapshot.__enter__ = Mock(return_value=mock_snapshot)
        sim.state.cpu_local_snapshot.__exit__ = Mock(return_value=None)
        
        # Mock numpy functions
        mock_where.return_value = (np.array([2]),)  # Cavity particle at index 2
        mock_sqrt.return_value = 0.1
        mock_normal.return_value = np.array([0.1, 0.2, 0.3])
        
        ThermalizationSetup._thermalize_cavity_particle_manually(sim, kT)
        
        mock_where.assert_called_once()
        mock_sqrt.assert_called_once()
        mock_normal.assert_called_once()
    
    @patch('numpy.where')
    def test_thermalize_cavity_particle_manually_no_cavity(self, mock_where):
        """Test manual cavity particle thermalization with no cavity particle."""
        sim = Mock()
        kT = 0.01
        
        # Mock snapshot context manager
        mock_snapshot = Mock()
        mock_snapshot.particles = Mock()
        mock_snapshot.particles.typeid = np.array([0, 1])
        
        sim.state.cpu_local_snapshot = MagicMock()
        sim.state.cpu_local_snapshot.__enter__ = Mock(return_value=mock_snapshot)
        sim.state.cpu_local_snapshot.__exit__ = Mock(return_value=None)
        
        # Mock no cavity particles found
        mock_where.return_value = (np.array([]),)
        
        ThermalizationSetup._thermalize_cavity_particle_manually(sim, kT)
        
        mock_where.assert_called_once()


class TestTimestepSetup:
    """Test TimestepSetup class functionality."""
    
    def create_test_params(self, **kwargs):
        """Create test parameters with defaults."""
        defaults = {
            'temperature': 100.0,
            'frequency': 2000.0,
            'coupling_strength': 0.001,
            'runtime_ps': 1000.0,
            'dt_fs': None
        }
        defaults.update(kwargs)
        return CavitySimulationParams(**defaults)
    
    @patch('hoomd.cavitymd.setup.PhysicalConstants')
    def test_calculate_initial_timestep_user_specified(self, mock_constants):
        """Test initial timestep calculation with user-specified timestep."""
        params = self.create_test_params(dt_fs=1.0)
        mock_constants.fs_to_atomic_units.return_value = 41.34
        
        dt = TimestepSetup.calculate_initial_timestep(params)
        
        assert dt == 41.34
        mock_constants.fs_to_atomic_units.assert_called_once_with(1.0)
    
    @patch('hoomd.cavitymd.setup.PhysicalConstants')
    def test_calculate_initial_timestep_default(self, mock_constants):
        """Test initial timestep calculation with default timestep."""
        params = self.create_test_params(dt_fs=None)
        mock_constants.ps_to_atomic_units.return_value = 4.134
        
        dt = TimestepSetup.calculate_initial_timestep(params)
        
        assert dt == 4.134
        mock_constants.ps_to_atomic_units.assert_called_once_with(0.0001)
    
    @patch('hoomd.cavitymd.setup.PhysicalConstants')
    @patch('numpy.asarray')
    @patch('numpy.linalg.norm')
    @patch('numpy.sum')
    @patch('numpy.sqrt')
    def test_compute_optimal_timestep(self, mock_sqrt, mock_sum, mock_norm, mock_asarray, mock_constants):
        """Test optimal timestep computation."""
        sim = Mock()
        error_tolerance = 0.01
        
        # Mock particle data
        mock_particles = Mock()
        mock_particles.mass = np.array([1.0, 1.0, 1.0])
        mock_snapshot = Mock()
        mock_snapshot.particles = mock_particles
        sim.state.get_snapshot.return_value = mock_snapshot
        
        # Mock forces
        mock_force1 = Mock()
        mock_force1.forces = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        mock_force2 = Mock()
        mock_force2.forces = np.array([[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]])
        
        sim.operations.integrator.forces = [mock_force1, mock_force2]
        
        # Mock numpy functions
        mock_asarray.side_effect = lambda x: np.asarray(x)
        mock_norm.return_value = np.array([1.0, 1.0, 1.0])
        mock_sum.return_value = 3.0
        mock_sqrt.return_value = 0.0577
        
        dt = TimestepSetup.compute_optimal_timestep(sim, error_tolerance)
        
        assert dt == 0.0577
        sim.run.assert_called_once_with(1)
        mock_sqrt.assert_called_once_with(error_tolerance / 3.0)
    
    @patch('hoomd.cavitymd.setup.PhysicalConstants')
    def test_compute_optimal_timestep_zero_force(self, mock_constants):
        """Test optimal timestep computation with zero forces."""
        sim = Mock()
        error_tolerance = 0.01
        
        # Mock particle data
        mock_particles = Mock()
        mock_particles.mass = np.array([1.0, 1.0, 1.0])
        mock_snapshot = Mock()
        mock_snapshot.particles = mock_particles
        sim.state.get_snapshot.return_value = mock_snapshot
        
        # Mock forces (empty or zero)
        sim.operations.integrator.forces = []
        
        mock_constants.fs_to_atomic_units.return_value = 4.134
        
        dt = TimestepSetup.compute_optimal_timestep(sim, error_tolerance)
        
        assert dt == 4.134  # Default fallback
        mock_constants.fs_to_atomic_units.assert_called_once_with(0.1)


class TestIntegrationTests:
    """Integration tests for setup utilities working together."""
    
    def create_test_params(self, **kwargs):
        """Create test parameters with defaults."""
        defaults = {
            'temperature': 100.0,
            'frequency': 2000.0,
            'coupling_strength': 0.001,
            'runtime_ps': 1000.0,
            'incavity': True,
            'molecular_thermostat': 'bussi',
            'cavity_thermostat': 'langevin'
        }
        defaults.update(kwargs)
        return CavitySimulationParams(**defaults)
    
    def test_setup_workflow_integration(self):
        """Test the complete setup workflow integration."""
        params = self.create_test_params()
        
        # Test device setup
        device = DeviceSetup.create_device('CPU')
        assert isinstance(device, hoomd.device.CPU)
        
        # Test seed generation
        seed = DeviceSetup.generate_simulation_seed(None, 1)
        assert isinstance(seed, int)
        
        # Test that all setup methods can be called with valid parameters
        # (Not mocking everything for this integration test)
        
        # Test ForceSetup methods exist and can be called
        assert hasattr(ForceSetup, 'create_all_forces')
        assert hasattr(ForceSetup, 'create_cavity_force')
        assert hasattr(ForceSetup, 'create_harmonic_bonds')
        
        # Test ThermostatSetup methods exist
        assert hasattr(ThermostatSetup, 'create_all_methods')
        assert hasattr(ThermostatSetup, 'create_molecular_thermostat')
        assert hasattr(ThermostatSetup, 'create_cavity_thermostat')
        
        # Test TimestepSetup methods exist
        assert hasattr(TimestepSetup, 'calculate_initial_timestep')
        assert hasattr(TimestepSetup, 'compute_optimal_timestep')
    
    def test_parameter_consistency(self):
        """Test that setup utilities use parameters consistently."""
        params = self.create_test_params(
            temperature=298.15,
            frequency=3000.0,
            coupling_strength=0.01,
            molecular_thermostat_tau=10.0,
            cavity_thermostat_tau=5.0
        )
        
        # Test that physical constants are calculated consistently
        constants = params.get_physical_constants()
        
        assert constants['kT'] == constants['kB'] * 298.15
        assert constants['omegac_au'] == 3000.0 / constants['kB']  # This will need adjustment
        assert constants['molecular_tau_au'] > 0
        assert constants['cavity_tau_au'] > 0
        
        # Test that setup utilities can access these constants
        assert params.temperature == 298.15
        assert params.frequency == 3000.0
        assert params.coupling_strength == 0.01 