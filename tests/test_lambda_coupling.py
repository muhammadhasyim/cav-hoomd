"""Tests for lambda-based coupling constant implementation."""

import pytest
import numpy as np
from hoomd.cavitymd.variants import LambdaScaledVariant
from hoomd.variant import Constant


class TestLambdaScaledVariant:
    """Test the LambdaScaledVariant wrapper class."""
    
    def test_constant_lambda_scaling(self):
        """Test that constant lambda is correctly scaled by omega_c."""
        lambda_val = 0.001
        omegac = 0.00913  # 2000 cm^-1 in a.u.
        
        variant = LambdaScaledVariant(lambda_val, omegac)
        
        # At any timestep, should return lambda * omega_c
        epsilon = variant(0)
        expected = lambda_val * omegac
        
        assert np.isclose(epsilon, expected), f"Expected {expected}, got {epsilon}"
    
    def test_variant_lambda_scaling(self):
        """Test that variant lambda is correctly scaled by omega_c."""
        lambda_variant = Constant(0.002)
        omegac = 0.00913
        
        variant = LambdaScaledVariant(lambda_variant, omegac)
        
        epsilon = variant(0)
        expected = 0.002 * omegac
        
        assert np.isclose(epsilon, expected), f"Expected {expected}, got {epsilon}"
    
    def test_float_to_constant_conversion(self):
        """Test that float lambda is converted to Constant variant."""
        lambda_val = 0.001
        omegac = 0.00913
        
        variant = LambdaScaledVariant(lambda_val, omegac)
        
        # Should work without error
        epsilon = variant(0)
        assert epsilon > 0
    
    def test_invalid_omegac(self):
        """Test that invalid omega_c raises error."""
        with pytest.raises(ValueError):
            LambdaScaledVariant(0.001, -0.00913)
        
        with pytest.raises(ValueError):
            LambdaScaledVariant(0.001, 0.0)
    
    def test_invalid_lambda_variant(self):
        """Test that invalid lambda variant raises error."""
        with pytest.raises(ValueError):
            LambdaScaledVariant("invalid", 0.00913)
    
    def test_multiple_timesteps(self):
        """Test that scaling is consistent across timesteps."""
        lambda_val = 0.001
        omegac = 0.00913
        
        variant = LambdaScaledVariant(lambda_val, omegac)
        
        # Should return same value for all timesteps (constant lambda)
        for timestep in [0, 100, 1000, 10000]:
            epsilon = variant(timestep)
            expected = lambda_val * omegac
            assert np.isclose(epsilon, expected)


class TestCavityMDSimulationLambdaCoupling:
    """Test CavityMDSimulation with lambda_coupling parameter."""
    
    def test_lambda_coupling_parameter_required(self):
        """Test that lambda_coupling is required when couplstr is not provided."""
        from hoomd.cavitymd.simulation import CavityMDSimulation

        # lambda_coupling is required
        with pytest.raises(ValueError, match="Must specify 'lambda_coupling'"):
            sim = CavityMDSimulation(
                job_dir="test",
                replica=1,
                freq=2000.0,
                incavity=True,
                # lambda_coupling not provided
            )
    
    def test_lambda_coupling_accepted(self):
        """Test that lambda_coupling parameter is accepted."""
        from hoomd.cavitymd.simulation import CavityMDSimulation
        
        # Should not raise error
        sim = CavityMDSimulation(
            job_dir="test",
            replica=1,
            freq=2000.0,
            incavity=True,
            lambda_coupling=0.001,
            runtime_ps=1.0
        )
        
        assert sim.lambda_coupling == 0.001
        assert sim.use_deprecated_couplstr is False
    
    def test_couplstr_deprecated_warning(self):
        """Test that using couplstr raises deprecation warning."""
        from hoomd.cavitymd.simulation import CavityMDSimulation
        import warnings
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            sim = CavityMDSimulation(
                job_dir="test",
                replica=1,
                freq=2000.0,
                incavity=True,
                couplstr=0.001,
                runtime_ps=1.0
            )
            
            # Check that deprecation warning was raised
            assert len(w) == 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert "deprecated" in str(w[0].message).lower()
    
    def test_cannot_specify_both_parameters(self):
        """Test that specifying both couplstr and lambda_coupling raises error."""
        from hoomd.cavitymd.simulation import CavityMDSimulation
        import warnings
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with pytest.raises(ValueError, match="Cannot specify both"):
                sim = CavityMDSimulation(
                    job_dir="test",
                    replica=1,
                    freq=2000.0,
                    incavity=True,
                    couplstr=0.001,
                    lambda_coupling=0.001,
                    runtime_ps=1.0
                )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

