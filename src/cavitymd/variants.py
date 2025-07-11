# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Custom variant classes for cavity molecular dynamics."""

import hoomd.variant
import numpy as np


class StepVariant(hoomd.variant.Variant):
    """
    Step function variant that switches at a specified time.
    
    This variant provides a completely discontinuous step function that switches
    from 0.0 to a target value at a specified switch time. Uses ElapsedTimeTracker
    for accurate time tracking in adaptive timestep simulations.
    
    Args:
        target_value (float): The value to switch to at switch_time
        switch_time_ps (float): The switch time in picoseconds
        time_tracker (ElapsedTimeTracker): Time tracker for accurate timing
    """
    
    def __init__(self, target_value, switch_time_ps, time_tracker):
        hoomd.variant.Variant.__init__(self)
        self.target_value = float(target_value)
        self.switch_time_ps = float(switch_time_ps)
        self.time_tracker = time_tracker
        
        # Store for debugging/logging
        self.current_value = 0.0
        self._has_switched = False
        
        print(f"StepVariant initialized:")
        print(f"  Target value: {self.target_value}")
        print(f"  Switch time: {self.switch_time_ps} ps")
        print(f"  Using ElapsedTimeTracker for accurate timing")
    
    def __call__(self, timestep):
        """Evaluate the step function at the given timestep."""
        # Get current elapsed time from the time tracker
        if self.time_tracker is not None:
            current_time_ps = self.time_tracker.elapsed_time
        else:
            # Fallback to timestep-based calculation (less accurate for adaptive timestep)
            # This should not be used in normal operation
            current_time_ps = 0.0
            
        # Step function: 0 before switch_time, target_value after
        if current_time_ps >= self.switch_time_ps:
            self.current_value = self.target_value
            if not self._has_switched:
                self._has_switched = True
                print(f"StepVariant: Switched ON at t = {current_time_ps:.6f} ps (target: {self.switch_time_ps:.6f} ps)")
        else:
            self.current_value = 0.0
            
        return self.current_value
    
    def _min(self):
        """Return the minimum value of this variant."""
        return 0.0
    
    def _max(self):
        """Return the maximum value of this variant."""
        return abs(self.target_value)


class ConstantVariant(hoomd.variant.Variant):
    """
    Simple constant variant for backward compatibility.
    
    Args:
        value (float): The constant value to return
    """
    
    def __init__(self, value):
        hoomd.variant.Variant.__init__(self)
        self.value = float(value)
    
    def __call__(self, timestep):
        """Return the constant value."""
        return self.value
    
    def _min(self):
        """Return the minimum value."""
        return self.value
    
    def _max(self):
        """Return the maximum value."""
        return self.value 