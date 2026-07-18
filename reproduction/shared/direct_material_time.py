"""Direct material-time construction from measured F(k,t) and TN fictive series.

Material time is defined as

    xi(t) = integral_{t_w=0}^{t} dt' / tau_0.1(t')

where tau_0.1 is the structural relaxation time at waiting time t_w obtained from
the F(k,t)/F(k,0) = 0.1 criterion on replica-averaged master curves.

Waiting time follows the paper staging convention: t_w = ref_num * switch_time_ps
with ref_000 at cavity turn-on (t_w = 0).

F(k,0) is the first nonzero F(k,t) value from master_fskt_ref_000, matching
``plot_fkt_analysis.find_relaxation_time`` and ``build_aging_repro_layout``.

The TN prediction uses the same integral with tau_fictive from the equilibrium
RelaxationTimeModel applied to the fictive structural temperature.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

DEFAULT_SWITCH_TIME_PS = 200.0
DEFAULT_FKT_THRESHOLD = 0.1


def staged_coupling_dir_name(tag: str, switch_time_ps: float = DEFAULT_SWITCH_TIME_PS) -> str:
    """Return the flat staged coupling directory name used by reproduction layouts."""
    return f"cavity_coupling_{tag}_switch_{switch_time_ps:.1f}ps"


def staged_coupling_dir(
    data_root: Path,
    tag: str,
    switch_time_ps: float = DEFAULT_SWITCH_TIME_PS,
) -> Path:
    """Return path to a staged coupling directory under ``data_root``."""
    return Path(data_root) / staged_coupling_dir_name(tag, switch_time_ps)


def load_master_fskt(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load lag time and F(k,t) from a master_fskt_ref_*.txt file."""
    lag_times: list[float] = []
    fkt_values: list[float] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or stripped.startswith("lag"):
                continue
            parts = stripped.split()
            lag_times.append(float(parts[0]))
            fkt_values.append(float(parts[1]))
    return np.asarray(lag_times, dtype=float), np.asarray(fkt_values, dtype=float)


def find_tau_f_over_f0_threshold(
    lag_time: np.ndarray,
    fkt: np.ndarray,
    threshold: float = DEFAULT_FKT_THRESHOLD,
    normalization_value: float | None = None,
) -> float:
    """Return lag time where F(k,t)/F(k,0) crosses ``threshold`` by linear interpolation."""
    if lag_time.size < 2:
        return float("nan")
    f0 = normalization_value if normalization_value not in (None, 0.0) else fkt[0]
    if f0 == 0:
        return float("nan")
    normalized = fkt / f0
    for idx in range(1, normalized.size):
        prev_val = normalized[idx - 1]
        curr_val = normalized[idx]
        if prev_val > threshold >= curr_val or prev_val >= threshold > curr_val:
            denom = curr_val - prev_val
            if abs(denom) < 1e-15:
                return float(lag_time[idx])
            fraction = (threshold - prev_val) / denom
            return float(lag_time[idx - 1] + fraction * (lag_time[idx] - lag_time[idx - 1]))
    return float("nan")


def _parse_ref_number(path: Path) -> int | None:
    """Extract reference index from a master_fskt_ref_* filename."""
    match = re.search(r"ref[_]?(\d+)", path.name)
    if not match:
        return None
    return int(match.group(1))


def ref0_normalization_value(coupling_dir: Path) -> float | None:
    """
    Return F(k,0) from ref_000: first nonzero F(k,t) in the pre-switch master curve.

    Matches ``build_aging_repro_layout`` and ``plot_fkt_analysis`` normalization.
    """
    coupling_dir = Path(coupling_dir)
    ref0_candidates = sorted(
        path
        for path in coupling_dir.glob("master_fskt_ref*.txt")
        if "_sample_counts" not in path.name and _parse_ref_number(path) == 0
    )
    if not ref0_candidates:
        return None

    _, fkt = load_master_fskt(ref0_candidates[0])
    nonzero = fkt[fkt != 0]
    if nonzero.size == 0:
        return None
    return float(nonzero[0])


def _find_relaxation_time_paper(
    lag_time: np.ndarray,
    fkt: np.ndarray,
    normalization_value: float | None,
    threshold: float,
) -> float:
    """Paper-style tau extraction via ``plot_fkt_analysis.find_relaxation_time``."""
    import sys

    shared_dir = Path(__file__).resolve().parent
    if str(shared_dir) not in sys.path:
        sys.path.insert(0, str(shared_dir))
    from plot_fkt_analysis import find_relaxation_time

    return float(
        find_relaxation_time(
            lag_time,
            fkt,
            target_value=threshold,
            normalization_value=normalization_value,
        )
    )


def collect_tau_vs_waiting_time(
    coupling_dir: Path,
    switch_time_ps: float = DEFAULT_SWITCH_TIME_PS,
    threshold: float = DEFAULT_FKT_THRESHOLD,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Collect tau_0.1 versus waiting time for all master references.

    Waiting time convention: t_w = ref_num * switch_time_ps with ref_000 at turn-on.
    """
    coupling_dir = Path(coupling_dir)
    norm_value = ref0_normalization_value(coupling_dir)
    waiting_times: list[float] = []
    relaxation_times: list[float] = []

    ref_paths = sorted(
        path
        for path in coupling_dir.glob("master_fskt_ref*.txt")
        if "_sample_counts" not in path.name
    )
    for ref_path in ref_paths:
        ref_num = _parse_ref_number(ref_path)
        if ref_num is None:
            continue

        lag_time, fkt = load_master_fskt(ref_path)
        tau = _find_relaxation_time_paper(lag_time, fkt, norm_value, threshold)
        if np.isnan(tau):
            continue

        waiting_times.append(ref_num * switch_time_ps)
        relaxation_times.append(tau)

    if not waiting_times:
        return np.array([]), np.array([])

    order = np.argsort(waiting_times)
    return (
        np.asarray(waiting_times, dtype=float)[order],
        np.asarray(relaxation_times, dtype=float)[order],
    )


def integrate_material_time(
    waiting_time: np.ndarray,
    relaxation_time: np.ndarray,
) -> np.ndarray:
    """Trapezoidal cumulative integral of 1/tau over waiting time."""
    if waiting_time.size == 0:
        return np.array([])
    if waiting_time.size == 1:
        return np.zeros(1, dtype=float)

    xi = np.zeros_like(waiting_time, dtype=float)
    for idx in range(1, waiting_time.size):
        dt = waiting_time[idx] - waiting_time[idx - 1]
        integrand = 0.5 * (1.0 / relaxation_time[idx - 1] + 1.0 / relaxation_time[idx])
        xi[idx] = xi[idx - 1] + integrand * dt
    return xi


def interpolate_material_time_grid(
    waiting_time: np.ndarray,
    xi: np.ndarray,
    n_points: int = 2500,
) -> tuple[np.ndarray, np.ndarray]:
    """Interpolate discrete (t_w, xi) onto a uniform grid for plotting."""
    if waiting_time.size == 0:
        return np.array([]), np.array([])
    if waiting_time.size == 1:
        return waiting_time.copy(), xi.copy()

    time_grid = np.linspace(0.0, waiting_time[-1], n_points)
    xi_grid = np.interp(time_grid, waiting_time, xi)
    return time_grid, xi_grid


def direct_material_time_from_masters(
    coupling_dir: Path,
    switch_time_ps: float = DEFAULT_SWITCH_TIME_PS,
    threshold: float = DEFAULT_FKT_THRESHOLD,
    n_interp: int = 2500,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build measured material time xi(t) from master F(k,t) files in ``coupling_dir``.

    Returns
    -------
    time_ps, xi
        Uniform grid in ps after cavity turn-on and corresponding material time.
    """
    coupling_dir = Path(coupling_dir)
    if not coupling_dir.is_dir():
        return np.array([]), np.array([])

    waiting_time, relaxation_time = collect_tau_vs_waiting_time(
        coupling_dir,
        switch_time_ps=switch_time_ps,
        threshold=threshold,
    )
    if waiting_time.size == 0:
        return np.array([]), np.array([])

    xi = integrate_material_time(waiting_time, relaxation_time)
    return interpolate_material_time_grid(waiting_time, xi, n_points=n_interp)


def tn_material_time_from_series(
    time_ps: np.ndarray,
    tau_fictive_ps: np.ndarray,
    switch_time_ps: float = DEFAULT_SWITCH_TIME_PS,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Integrate TN material time from fictive relaxation times without 1/ln(10) scaling.

    Returns
    -------
    time_shifted, xi_shifted, aging_rate
        Time since cavity turn-on (ps), material time referenced to turn-on, and 1/tau.
    """
    time_ps = np.asarray(time_ps, dtype=float)
    tau_fictive_ps = np.asarray(tau_fictive_ps, dtype=float)

    safe_tau = np.where(tau_fictive_ps > 0.0, tau_fictive_ps, np.nan)
    inv_tau = np.where(np.isfinite(safe_tau), 1.0 / safe_tau, 0.0)

    switch_idx = int(np.argmin(np.abs(time_ps - switch_time_ps)))
    switch_time = time_ps[switch_idx]

    xi = np.zeros_like(time_ps, dtype=float)
    for idx in range(switch_idx + 1, time_ps.size):
        dt = time_ps[idx] - time_ps[idx - 1]
        xi[idx] = xi[idx - 1] + 0.5 * (inv_tau[idx - 1] + inv_tau[idx]) * dt

    post_switch = time_ps >= switch_time
    time_shifted = time_ps[post_switch] - switch_time
    xi_shifted = xi[post_switch] - xi[switch_idx]
    aging_rate = inv_tau[post_switch]
    return time_shifted, xi_shifted, aging_rate


def load_direct_measured_mt_all(
    data_root: Path,
    coupling_values: list[float],
    coupling_tags: dict[float, str],
    switch_time_ps: float = DEFAULT_SWITCH_TIME_PS,
    n_interp: int = 2500,
) -> dict[float, tuple[np.ndarray, np.ndarray]]:
    """Load direct measured material time for all couplings under a staged layout."""
    data_root = Path(data_root)
    result: dict[float, tuple[np.ndarray, np.ndarray]] = {}

    for axis_value in coupling_values:
        tag = coupling_tags[axis_value]
        coupling_dir = staged_coupling_dir(data_root, tag, switch_time_ps)
        time_grid, xi_grid = direct_material_time_from_masters(
            coupling_dir,
            switch_time_ps=switch_time_ps,
            n_interp=n_interp,
        )
        if time_grid.size:
            result[axis_value] = (time_grid, xi_grid)

    return result
