"""Pure-Python helpers for aging campaign checkpoint and resume I/O."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

try:
    from cavitymd.analysis.fkt_state import load_fkt_state, save_fkt_state
except ImportError:  # installed plugin lives under hoomd.cavitymd
    from hoomd.cavitymd.analysis.fkt_state import load_fkt_state, save_fkt_state

REFERENCE_TIME_PATTERN = re.compile(
    r"Reference time:\s*([0-9.eE+-]+)\s*ps"
)
FKT_REF_INDEX_PATTERN = re.compile(r"_fkt_ref_(\d+)\.txt$")


def parse_fkt_reference_time(path: Path) -> float | None:
    """Return the reference laboratory time encoded in an F(k,t) header."""
    try:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                match = REFERENCE_TIME_PATTERN.search(line)
                if match is not None:
                    return float(match.group(1))
    except OSError:
        return None
    return None


def parse_fkt_reference_index(path: Path) -> int | None:
    """Return the zero-based reference index from a production F(k,t) filename."""
    match = FKT_REF_INDEX_PATTERN.search(path.name)
    if match is None:
        return None
    return int(match.group(1))


def max_lag_for_reference(
    *,
    reference_time_ps: float,
    runtime_ps: float,
) -> float:
    """Return the maximum lag window available for one reference."""
    return max(0.0, runtime_ps - reference_time_ps)


def normalize_fkt_series(
    lag_times_ps: np.ndarray,
    values: np.ndarray,
    *,
    reference_time_ps: float,
    runtime_ps: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Normalize F(k,t) to F/F(0) and clip to the available lag window."""
    lag_times = np.asarray(lag_times_ps, dtype=float)
    series = np.asarray(values, dtype=float)
    if lag_times.size == 0:
        return lag_times, series

    order = np.argsort(lag_times)
    lag_times = lag_times[order]
    series = series[order]

    positive = series != 0.0
    if np.any(positive):
        normalization = series[np.argmax(positive)]
        if normalization != 0.0:
            series = series / normalization

    max_lag = max_lag_for_reference(
        reference_time_ps=reference_time_ps,
        runtime_ps=runtime_ps,
    )
    mask = (lag_times >= 0.0) & (lag_times <= max_lag + 1e-9)
    return lag_times[mask], series[mask]


def infer_hdf5_write_index(time_values: np.ndarray) -> int:
    """Infer the next append index for a monotonic HDF5 ``/time`` dataset."""
    times = np.asarray(time_values, dtype=float)
    if times.size == 0:
        return 0

    index = 0
    current_max = -np.inf
    for position, value in enumerate(times):
        if not np.isfinite(value):
            break
        if position > 0 and value <= 0.0 and current_max > 0.0:
            break
        if position > 0 and value + 1e-9 < current_max:
            break
        current_max = max(current_max, float(value))
        index = position + 1
    return index


__all__ = [
    "infer_hdf5_write_index",
    "load_fkt_state",
    "max_lag_for_reference",
    "normalize_fkt_series",
    "parse_fkt_reference_index",
    "parse_fkt_reference_time",
    "save_fkt_state",
]
