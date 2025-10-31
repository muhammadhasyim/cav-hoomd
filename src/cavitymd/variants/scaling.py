# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Scaling variants for cavity molecular dynamics."""

from typing import Union, Optional
import hoomd.variant


class LambdaScaledVariant(hoomd.variant.Variant):
    r"""
    Wrapper variant that scales a base variant by omega_c to compute epsilon.
    
    This variant implements the relationship:
    
    .. math::
        
        \varepsilon(t) = \lambda(t) \cdot \omega_c
    
    where:
    - :math:`\lambda(t)` is the user-provided dimensionless coupling parameter (time-varying or constant)
    - :math:`\omega_c` is the cavity characteristic frequency in atomic units
    - :math:`\varepsilon(t)` is the effective coupling strength fed to the molecular dynamics
    
    **Physical Interpretation:**
    
    This scaling allows the user to work with a dimensionless coupling parameter :math:`\lambda`
    that is independent of the cavity frequency. The actual coupling strength :math:`\varepsilon`
    is automatically computed as the product of :math:`\lambda` and the cavity frequency.
    
    This is useful because:
    - :math:`\lambda` is a dimensionless parameter that characterizes the strength of coupling
    - Different cavity frequencies can be studied with the same :math:`\lambda` values
    - The relationship :math:`\varepsilon = \lambda \omega_c` is physically motivated
    
    Parameters
    ----------
    lambda_variant : float or hoomd.variant.Variant
        The dimensionless coupling parameter :math:`\lambda`. Can be:
        - A constant float value
        - A hoomd.variant.Variant object for time-varying coupling
    omegac : float
        The cavity characteristic frequency :math:`\omega_c` in atomic units (Hartree).
        This is typically the cavity mode frequency.
    
    Attributes
    ----------
    lambda_variant : hoomd.variant.Variant
        The wrapped lambda variant
    omegac : float
        The cavity frequency in atomic units
    
    Examples
    --------
    **Constant lambda with fixed cavity frequency:**
    
    >>> from hoomd.cavitymd.variants import LambdaScaledVariant
    >>> 
    >>> # Set lambda = 0.001, omega_c = 0.00913 a.u. (2000 cm⁻¹)
    >>> # This gives epsilon = 0.001 * 0.00913 = 9.13e-6 a.u.
    >>> scaled_coupling = LambdaScaledVariant(
    ...     lambda_variant=0.001,
    ...     omegac=0.00913
    ... )
    
    **Time-varying lambda with step function:**
    
    >>> from hoomd.cavitymd.variants import StepVariant, LambdaScaledVariant
    >>> from hoomd.cavitymd.analysis import ElapsedTimeTracker
    >>> 
    >>> # Create time tracker
    >>> time_tracker = ElapsedTimeTracker(sim, runtime_ps=10.0)
    >>> 
    >>> # Create step variant for lambda: switch at 1 ps to lambda = 0.001
    >>> lambda_variant = StepVariant(
    ...     target_value=0.001,
    ...     switch_time_ps=1.0,
    ...     time_tracker=time_tracker
    ... )
    >>> 
    >>> # Wrap with cavity frequency scaling
    >>> scaled_coupling = LambdaScaledVariant(
    ...     lambda_variant=lambda_variant,
    ...     omegac=0.00913  # 2000 cm⁻¹
    ... )
    >>> 
    >>> # Use in cavity force
    >>> cavity_force = CavityForce(
    ...     kvector=[0, 0, 1],
    ...     couplstr=scaled_coupling,
    ...     omegac=0.00913
    ... )
    
    **Periodic lambda with frequency scaling:**
    
    >>> from hoomd.cavitymd.variants import PeriodicVariant, LambdaScaledVariant
    >>> 
    >>> # Create periodic lambda: oscillate with amplitude 0.001, period 1 ps
    >>> lambda_variant = PeriodicVariant(
    ...     amplitude=0.001,
    ...     time_tracker=time_tracker,
    ...     period_ps=1.0
    ... )
    >>> 
    >>> # Scale by cavity frequency
    >>> scaled_coupling = LambdaScaledVariant(
    ...     lambda_variant=lambda_variant,
    ...     omegac=0.00913
    ... )
    
    Notes
    -----
    - The lambda parameter is dimensionless and independent of cavity frequency
    - The effective coupling epsilon = lambda * omega_c is what gets used in MD
    - This allows studying the same lambda values across different cavity frequencies
    - The wrapper is transparent to the rest of the simulation code
    - Works with any hoomd.variant.Variant subclass for lambda
    
    See Also
    --------
    hoomd.cavitymd.forces.CavityForce : Force class that uses the scaled coupling
    hoomd.cavitymd.variants.StepVariant : For time-varying lambda
    hoomd.cavitymd.variants.PeriodicVariant : For periodic lambda
    """
    
    def __init__(self, 
                 lambda_variant: Union[float, hoomd.variant.Variant],
                 omegac: float) -> None:
        """
        Initialize the lambda-scaled variant.
        
        Parameters
        ----------
        lambda_variant : float or hoomd.variant.Variant
            The dimensionless coupling parameter lambda
        omegac : float
            The cavity frequency in atomic units
        """
        super().__init__()
        
        # Convert constant float to Constant variant if needed
        if isinstance(lambda_variant, (int, float)):
            self.lambda_variant = hoomd.variant.Constant(float(lambda_variant))
        elif hasattr(lambda_variant, '__call__'):  # It's a variant
            self.lambda_variant = lambda_variant
        else:
            raise ValueError("lambda_variant must be a number or hoomd.variant.Variant")
        
        self.omegac = float(omegac)
        
        if self.omegac <= 0:
            raise ValueError(f"omegac must be positive, got {self.omegac}")
    
    def __call__(self, timestep: int) -> float:
        """
        Compute the scaled coupling strength at a given timestep.
        
        Parameters
        ----------
        timestep : int
            The current timestep
            
        Returns
        -------
        float
            The effective coupling strength epsilon = lambda(timestep) * omega_c
        """
        lambda_value = self.lambda_variant(timestep)
        return lambda_value * self.omegac

