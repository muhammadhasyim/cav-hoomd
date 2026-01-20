
# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Analysis and tracking tools for cavity molecular dynamics simulations."""

from typing import Optional, List, Dict, Union, Any, Tuple
import hoomd
import numpy as np
import logging
import sys
import os
import time
from pathlib import Path
import enum
import scipy.fft

# CuPy import with fallback for CPU/GPU agnostic code
try:
    import cupy as cp
    HAS_CUPY = True
except ImportError:
    cp = None
    HAS_CUPY = False

from ..utils import PhysicalConstants, unwrap_positions


class Status:
    """Simulation status tracker that provides real-time status information."""
    
    def __init__(self, sim, runtime_ps, time_tracker):
        """Initialize status tracker.
        
        Parameters
        ----------
        sim : hoomd.Simulation
            HOOMD simulation object
        runtime_ps : float
            Total runtime in picoseconds
        time_tracker : ElapsedTimeTracker
            Time tracking object
        """
        self.sim = sim
        self.runtime_ps = runtime_ps
        self.time_tracker = time_tracker
        self._status = "initializing"
    
    @property 
    def status(self):
        """Current simulation status."""
        return self._status
    
    @status.setter
    def status(self, value):
        """Set simulation status."""
        self._status = value
        
    def update(self):
        """Update status based on current simulation state."""
        if hasattr(self.time_tracker, 'total_time'):
            elapsed_ps = PhysicalConstants.atomic_units_to_ps(self.time_tracker.total_time)
            if elapsed_ps >= self.runtime_ps:
                self._status = "completed"
            elif elapsed_ps > 0:
                self._status = "running"
        
    def get_progress(self):
        """Get current progress as a percentage."""
        if hasattr(self.time_tracker, 'total_time'):
            elapsed_ps = PhysicalConstants.atomic_units_to_ps(self.time_tracker.total_time)
            return min(100.0, (elapsed_ps / self.runtime_ps) * 100.0)
        return 0.0
    

class TimestepFormatter:
    """Utility class for formatting timestep-related output."""
    
    def __init__(self, integrator=None):
        """Initialize timestep formatter.
        
        Parameters
        ----------
        integrator : hoomd.md.Integrator, optional
            HOOMD integrator object for timestep information
        """
        self.integrator = integrator
    
    @staticmethod
    def format_timestep(timestep: int) -> str:
        """Format timestep with appropriate scale."""
        if timestep >= 1e6:
            return f"{timestep/1e6:.1f}M"
        elif timestep >= 1e3:
            return f"{timestep/1e3:.1f}k"
        else:
            return str(timestep)
    
    @staticmethod
    def format_time(time_ps: float) -> str:
        """Format time in appropriate units."""
        if time_ps >= 1000:
            return f"{time_ps/1000:.1f} ns"
        else:
            return f"{time_ps:.1f} ps"
            
    def get_current_timestep_info(self):
        """Get current timestep information from integrator."""
        if self.integrator is not None:
            return {
                'dt': self.integrator.dt,
                'dt_ps': PhysicalConstants.atomic_units_to_ps(self.integrator.dt)
            }
        return None
    
    # Note: dt_fs method temporarily removed due to HOOMD logging metaclass issues
    # TODO: Re-add when logging issue is resolved


class ElapsedTimeTracker(hoomd.custom.Action):
    """Tracks the total elapsed time in a simulation with variable timesteps."""
    
    def __init__(self, simulation, runtime):
        super().__init__()
        self.simulation = simulation
        self.total_time = 0.0
        self.runtime = runtime
        self.last_timestep = 0  # Start from 0, not simulation.timestep
        self.initial_timestep = None  # Track the starting timestep to handle inheritance

    def act(self, timestep):
        """Update the total elapsed time by accumulating time increments."""
        # Get current timestep size
        dt = self.simulation.operations.integrator.dt
        
        # For the first call, handle initialization
        if self.last_timestep == 0:
            # Initialize - record the starting timestep but don't add its time
            self.initial_timestep = timestep
            self.last_timestep = timestep
            self.total_time = 0.0  # Always start elapsed time from 0, regardless of inherited timestep
            if timestep > 0:
                print(f"NOTICE: Starting from inherited timestep {timestep}")
                print(f"  Elapsed time will start from 0, not from inherited simulation time")
            return

        # Calculate time increment since last update
        if timestep > self.last_timestep:
            timestep_increment = timestep - self.last_timestep
            time_increment = timestep_increment * dt
            self.total_time += time_increment
            
        # Update last timestep for next iteration
        self.last_timestep = timestep
        
        # Check if we've reached the runtime and signal to stop
        if PhysicalConstants.atomic_units_to_ps(self.total_time) >= self.runtime:
            print(f"Runtime {self.runtime} ps reached. Exiting simulation.")
            # Raise StopIteration to signal simulation should stop
            # This will be caught by the simulation.run() loop
            raise StopIteration("Runtime reached")

    @hoomd.logging.log(category="scalar")
    def elapsed_time(self):
        """Return elapsed time in picoseconds."""
        return PhysicalConstants.atomic_units_to_ps(self.total_time)


