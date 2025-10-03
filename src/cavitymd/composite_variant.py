"""
Composite variant for combining multiple coupling variants additively.
"""

import hoomd
import numpy as np
from typing import List, Optional
from .analysis import ElapsedTimeTracker


class CompositeVariant(hoomd.variant.Variant):
    r"""
    Composite variant that additively combines multiple coupling variants.
    
    This variant enables superimposition of different coupling protocols by
    adding their individual contributions:
    
    .. math::
        
        g(t) = \sum_{i=1}^{N} g_i(t)
    
    where each :math:`g_i(t)` is an individual variant (sinusoidal, square wave, etc.).
    
    **Physical Interpretation:**
    
    Composite coupling simulates experimental scenarios with multiple simultaneous
    driving protocols:
    
    - Sinusoidal base coupling + adaptive square wave modulation
    - Multiple frequency components (e.g., fundamental + harmonics)
    - Combined AC/DC driving protocols
    - Simultaneous cavity + mechanical modulation
    - Multi-mode cavity coupling with different temporal patterns
    
    **Energy Conservation:**
    
    Total energy input is the sum of contributions from all component variants.
    Care should be taken to ensure the combined amplitude doesn't become
    unphysically large.
    
    Parameters
    ----------
    variants : List[hoomd.variant.Variant]
        List of individual variants to combine additively
    max_amplitude : float, optional
        Maximum allowed total amplitude. If exceeded, all variants are
        scaled proportionally. Default: None (no limit)
    
    Attributes
    ----------
    variants : List[hoomd.variant.Variant]
        Component variants being combined
    max_amplitude : float or None
        Maximum amplitude limit
    current_value : float
        Current combined value
    current_components : List[float]
        Current values of individual component variants
        
    Examples
    --------
    **Sinusoid + Adaptive Square Wave:**
    
    >>> from hoomd.cavitymd.composite_variant import CompositeVariant
    >>> from hoomd.cavitymd.variants import PeriodicVariant, AdaptiveSquareWaveVariant
    >>> from hoomd.cavitymd.analysis import ElapsedTimeTracker
    >>> 
    >>> # Create time tracker
    >>> time_tracker = ElapsedTimeTracker(sim, runtime_ps=100.0)
    >>> 
    >>> # Create sinusoidal base (1 ps period, 1e-4 amplitude)
    >>> sinusoid = PeriodicVariant(
    ...     amplitude=1e-4,
    ...     time_tracker=time_tracker,
    ...     period_ps=1.0
    ... )
    >>> 
    >>> # Create adaptive square wave (50 ps period, 2e-4 amplitude)
    >>> square_wave = AdaptiveSquareWaveVariant(
    ...     base_amplitude=2e-4,
    ...     time_tracker=time_tracker,
    ...     period_ps=50.0,
    ...     duty_cycle=0.02,
    ...     start_time_ps=10.0
    ... )
    >>> 
    >>> # Combine them additively
    >>> composite_coupling = CompositeVariant(
    ...     variants=[sinusoid, square_wave],
    ...     max_amplitude=5e-4  # Safety limit
    ... )
    
    **Multiple Frequency Components:**
    
    >>> # Fundamental + harmonic
    >>> fundamental = PeriodicVariant(
    ...     amplitude=1e-3,
    ...     time_tracker=time_tracker,
    ...     period_ps=2.0
    ... )
    >>> 
    >>> harmonic = PeriodicVariant(
    ...     amplitude=2e-4,
    ...     time_tracker=time_tracker,
    ...     period_ps=1.0,  # 2x frequency
    ...     phase_offset=np.pi/4
    ... )
    >>> 
    >>> multi_freq = CompositeVariant(
    ...     variants=[fundamental, harmonic]
    ... )
    
    Notes
    -----
    - All component variants must be compatible with the same simulation
    - Component variants are evaluated independently at each timestep
    - Time tracking must be consistent across all components
    - Auto-stop signals from any component will affect the composite
    """
    
    def __init__(self, 
                 variants: List[hoomd.variant.Variant],
                 max_amplitude: Optional[float] = None) -> None:
        
        # Call parent class constructor
        super().__init__()
        
        if not variants:
            raise ValueError("At least one variant must be provided")
        
        self.variants = variants
        self.max_amplitude = max_amplitude
        self.current_value = 0.0
        self.current_components = [0.0] * len(variants)
        
        # Track initialization
        self._initialized = False
        
        print(f"CompositeVariant initialized:")
        print(f"   Number of components: {len(self.variants)}")
        print(f"   Max amplitude limit: {self.max_amplitude}")
        for i, variant in enumerate(self.variants):
            variant_type = type(variant).__name__
            print(f"   Component {i+1}: {variant_type}")
    
    def __call__(self, timestep: int) -> float:
        """
        Evaluate the composite variant by summing all component variants.
        
        Parameters
        ----------
        timestep : int
            Current simulation timestep
            
        Returns
        -------
        float
            Combined coupling strength
        """
        if not self._initialized:
            print(f"CompositeVariant: Starting evaluation at timestep {timestep}")
            self._initialized = True
        
        # Evaluate all component variants
        total_value = 0.0
        for i, variant in enumerate(self.variants):
            component_value = variant(timestep)
            self.current_components[i] = component_value
            total_value += component_value
        
        # Apply amplitude limit if specified
        if self.max_amplitude is not None and abs(total_value) > self.max_amplitude:
            # Scale all components proportionally
            if total_value != 0:
                scale_factor = self.max_amplitude / abs(total_value)
                total_value *= scale_factor
                for i in range(len(self.current_components)):
                    self.current_components[i] *= scale_factor
                
                if not hasattr(self, '_amplitude_warning_shown'):
                    print(f"CompositeVariant: Amplitude limited to {self.max_amplitude} (scale factor: {scale_factor:.3f})")
                    self._amplitude_warning_shown = True
        
        self.current_value = total_value
        return self.current_value
    
    def _min(self) -> float:
        """Return minimum possible value."""
        return sum(getattr(variant, '_min', lambda: 0.0)() for variant in self.variants)
    
    def _max(self) -> float:
        """Return maximum possible value."""
        if self.max_amplitude is not None:
            return self.max_amplitude
        return sum(getattr(variant, '_max', lambda: 1.0)() for variant in self.variants)
