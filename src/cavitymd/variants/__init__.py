# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Variants submodule for time-dependent coupling protocols."""

from .basic import StepVariant, ConstantVariant
from .periodic import PeriodicVariant, SquareWaveVariant
from .decay import ExponentialDecayVariant, DecayingSquareWaveVariant
from .adaptive import AdaptiveSquareWaveVariant, ExponentialWaveVariant
from .scaling import LambdaScaledVariant

__all__ = [
    # Basic variants
    'StepVariant',
    'ConstantVariant',

    # Periodic variants
    'PeriodicVariant',
    'SquareWaveVariant',

    # Decay variants
    'ExponentialDecayVariant',
    'DecayingSquareWaveVariant',

    # Adaptive variants
    'AdaptiveSquareWaveVariant',
    'ExponentialWaveVariant',

    # Scaling variants
    'LambdaScaledVariant',
]

