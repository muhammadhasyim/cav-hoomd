"""
Stretched-exponential (KWW) fitting for structural relaxation times.

Fits ``phi(t) = A * exp(-(t / tau_K)^beta)`` to normalized F(k,t) curves,
optionally with a shared beta across waiting times of one coupling, and reports
the analytic F=0.1 crossing of the fitted curve as ``tau_s``.

Sample counts enter as heteroscedastic weights (``sigma ~ 1/sqrt(N)``) so
low-replica references do not dominate the estimate.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional, Sequence

import numpy as np
from scipy.optimize import least_squares

# Fit window: must extend below the reporting threshold (0.1) so the crossing
# is interpolated from the model, not extrapolated beyond the data.
PHI_FIT_LO = 0.03
PHI_FIT_HI = 0.90
MIN_FIT_POINTS = 8
DEFAULT_TARGET = 0.1


def kww(t: np.ndarray | float, A: float, tau_K: float, beta: float) -> np.ndarray | float:
    """Kohlrausch–Williams–Watts stretched exponential.

    Parameters
    ----------
    t : array_like
        Lag time (same units as ``tau_K``).
    A : float
        Amplitude (plateau height), dimensionless.
    tau_K : float
        KWW characteristic time.
    beta : float
        Stretching exponent (typically 0 < beta <= 1).

    Returns
    -------
    array_like
        ``A * exp(-(t / tau_K)^beta)``.
    """
    t_arr = np.asarray(t, dtype=np.float64)
    tau_safe = max(float(tau_K), 1e-30)
    # Clip the exponent argument to avoid overflow for large t/tau.
    arg = np.clip((t_arr / tau_safe) ** float(beta), 0.0, 700.0)
    return A * np.exp(-arg)


def relaxation_time_from_kww(
    A: float,
    tau_K: float,
    beta: float,
    target: float = DEFAULT_TARGET,
) -> float:
    """Analytic lag time where the KWW curve equals ``target``.

    Solves ``A * exp(-(t/tau_K)^beta) = target`` for ``t``:
    ``t = tau_K * (ln(A/target))^(1/beta)``.

    Returns NaN if ``A <= target``, or if parameters are non-positive.
    """
    if not (A > target > 0.0 and tau_K > 0.0 and beta > 0.0):
        return float("nan")
    return float(tau_K * (np.log(A / target)) ** (1.0 / beta))


def _prepare_fit_arrays(
    time: np.ndarray,
    phi: np.ndarray,
    weights: Optional[np.ndarray] = None,
    phi_lo: float = PHI_FIT_LO,
    phi_hi: float = PHI_FIT_HI,
) -> Optional[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Filter and sort data for KWW fitting; return (t, phi, sigma) or None."""
    time = np.asarray(time, dtype=np.float64)
    phi = np.asarray(phi, dtype=np.float64)
    if weights is None:
        weights = np.ones_like(time)
    else:
        weights = np.asarray(weights, dtype=np.float64)
        if weights.shape != time.shape:
            weights = np.ones_like(time)

    mask = (
        np.isfinite(time)
        & np.isfinite(phi)
        & (time > 0)
        & (phi >= phi_lo)
        & (phi <= phi_hi)
        & (weights > 0)
    )
    if np.count_nonzero(mask) < MIN_FIT_POINTS:
        return None

    t = time[mask]
    p = phi[mask]
    w = weights[mask]
    order = np.argsort(t)
    t, p, w = t[order], p[order], w[order]
    sigma = 1.0 / np.sqrt(w)
    return t, p, sigma


def _r_squared(y: np.ndarray, yhat: np.ndarray) -> float:
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    if ss_tot <= 0:
        return 0.0
    return 1.0 - ss_res / ss_tot


def _pack_result(
    A: float,
    tau_K: float,
    beta: float,
    r_squared: float,
    target: float = DEFAULT_TARGET,
) -> dict[str, float]:
    return {
        "A": float(A),
        "tau_K": float(tau_K),
        "beta": float(beta),
        "r_squared": float(r_squared),
        "tau_s": relaxation_time_from_kww(A, tau_K, beta, target=target),
    }


def fit_kww_single(
    time: np.ndarray,
    phi: np.ndarray,
    weights: Optional[np.ndarray] = None,
    frozen_beta: Optional[float] = None,
    target: float = DEFAULT_TARGET,
) -> Optional[dict[str, float]]:
    """Fit a single KWW curve with optional frozen beta.

    Parameters
    ----------
    time, phi :
        Lag time and normalized F(k,t).
    weights :
        Sample counts (or any positive relative weights). Converted to
        ``sigma = 1/sqrt(weight)`` for weighted least squares.
    frozen_beta :
        If provided, beta is held fixed at this value.
    target :
        Crossing level used to report ``tau_s`` (default 0.1).

    Returns
    -------
    dict or None
        Keys: ``A``, ``tau_K``, ``beta``, ``r_squared``, ``tau_s``.
        None if there are too few points or the fit fails.
    """
    prepared = _prepare_fit_arrays(time, phi, weights)
    if prepared is None:
        return None
    t, p, sigma = prepared

    # Initial guesses from the data.
    A0 = float(np.clip(np.max(p), 0.2, 1.2))
    # Rough tau from where phi ~ A0/e.
    idx_e = int(np.argmin(np.abs(p - A0 / np.e)))
    tau0 = float(max(t[idx_e], t[len(t) // 4], 1.0))
    beta0 = float(frozen_beta) if frozen_beta is not None else 0.65

    if frozen_beta is not None:
        def residuals(params: np.ndarray) -> np.ndarray:
            A, tau_K = params
            return (kww(t, A, tau_K, frozen_beta) - p) / sigma

        x0 = np.array([A0, tau0], dtype=np.float64)
        bounds = ([0.2, 1e-3], [1.5, 1e6])
    else:
        def residuals(params: np.ndarray) -> np.ndarray:
            A, tau_K, beta = params
            return (kww(t, A, tau_K, beta) - p) / sigma

        x0 = np.array([A0, tau0, beta0], dtype=np.float64)
        bounds = ([0.2, 1e-3, 0.15], [1.5, 1e6, 1.5])

    try:
        sol = least_squares(
            residuals, x0, bounds=bounds, method="trf",
            max_nfev=5000, ftol=1e-10, xtol=1e-10,
        )
    except Exception:
        return None

    if not sol.success and sol.cost > 1e3:
        return None

    if frozen_beta is not None:
        A_fit, tau_fit = sol.x
        beta_fit = float(frozen_beta)
    else:
        A_fit, tau_fit, beta_fit = sol.x

    yhat = kww(t, A_fit, tau_fit, beta_fit)
    return _pack_result(A_fit, tau_fit, beta_fit, _r_squared(p, yhat), target=target)


def _weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    cum = np.cumsum(w)
    cutoff = 0.5 * cum[-1]
    idx = int(np.searchsorted(cum, cutoff))
    return float(v[min(idx, len(v) - 1)])


def fit_kww_shared_beta(
    curves: Sequence[Mapping[str, Any]],
    target: float = DEFAULT_TARGET,
) -> Optional[list[dict[str, float]]]:
    """Joint KWW fit across multiple curves with one shared beta.

    Each curve is a mapping with keys ``time``, ``phi``, and optional
    ``weights``.  Parameter vector is ``[beta, A_0..A_R, tauK_0..tauK_R]``.

    Falls back to a two-stage procedure (free per-curve fits -> weighted-median
    beta -> refit with frozen beta) if the joint solve fails.

    Returns
    -------
    list of dict or None
        One result dict per input curve (same order), or None if all curves
        are unusable.
    """
    if not curves:
        return None

    prepared: list[Optional[tuple[np.ndarray, np.ndarray, np.ndarray]]] = []
    for curve in curves:
        prepared.append(
            _prepare_fit_arrays(
                curve["time"],
                curve["phi"],
                curve.get("weights"),
            )
        )

    usable_idx = [i for i, p in enumerate(prepared) if p is not None]
    if not usable_idx:
        return None

    n_u = len(usable_idx)
    # Initial guesses from independent fits.
    free_results: list[Optional[dict[str, float]]] = [None] * len(curves)
    beta_vals = []
    beta_wts = []
    for i in usable_idx:
        t, p, sigma = prepared[i]  # type: ignore[misc]
        w = 1.0 / (sigma ** 2)
        free = fit_kww_single(t, p, weights=w, target=target)
        free_results[i] = free
        if free is not None and np.isfinite(free["beta"]):
            beta_vals.append(free["beta"])
            beta_wts.append(float(np.sum(w)))

    if not beta_vals:
        return None

    beta_guess = _weighted_median(np.asarray(beta_vals), np.asarray(beta_wts))

    # Build joint residual: [beta, A_0.., tau_0..]
    A0 = []
    tau0 = []
    for i in usable_idx:
        fr = free_results[i]
        if fr is not None:
            A0.append(fr["A"])
            tau0.append(fr["tau_K"])
        else:
            t, p, _ = prepared[i]  # type: ignore[misc]
            A0.append(float(np.clip(np.max(p), 0.2, 1.2)))
            tau0.append(float(max(t[len(t) // 4], 1.0)))

    x0 = np.array([beta_guess] + A0 + tau0, dtype=np.float64)
    lb = np.array([0.15] + [0.2] * n_u + [1e-3] * n_u)
    ub = np.array([1.5] + [1.5] * n_u + [1e6] * n_u)

    def joint_residuals(params: np.ndarray) -> np.ndarray:
        beta = params[0]
        As = params[1 : 1 + n_u]
        taus = params[1 + n_u :]
        chunks = []
        for j, i in enumerate(usable_idx):
            t, p, sigma = prepared[i]  # type: ignore[misc]
            chunks.append((kww(t, As[j], taus[j], beta) - p) / sigma)
        return np.concatenate(chunks)

    joint_ok = False
    try:
        sol = least_squares(
            joint_residuals, x0, bounds=(lb, ub), method="trf",
            max_nfev=8000, ftol=1e-10, xtol=1e-10,
        )
        joint_ok = sol.success or sol.cost < joint_residuals(x0).dot(joint_residuals(x0))
    except Exception:
        joint_ok = False
        sol = None

    results: list[Optional[dict[str, float]]] = [None] * len(curves)

    if joint_ok and sol is not None:
        beta_fit = float(sol.x[0])
        As = sol.x[1 : 1 + n_u]
        taus = sol.x[1 + n_u :]
        for j, i in enumerate(usable_idx):
            t, p, _ = prepared[i]  # type: ignore[misc]
            yhat = kww(t, As[j], taus[j], beta_fit)
            results[i] = _pack_result(
                As[j], taus[j], beta_fit, _r_squared(p, yhat), target=target
            )
    else:
        # Two-stage fallback: freeze weighted-median beta, refit A and tau_K.
        for i in usable_idx:
            t, p, sigma = prepared[i]  # type: ignore[misc]
            w = 1.0 / (sigma ** 2)
            refit = fit_kww_single(t, p, weights=w, frozen_beta=beta_guess, target=target)
            results[i] = refit

    # Ensure every usable curve got a result with the same beta when possible.
    betas_final = [r["beta"] for r in results if r is not None]
    if betas_final and len(set(np.round(betas_final, decimals=10))) > 1:
        # Force shared beta via frozen refit for consistency.
        shared = float(np.median(betas_final))
        for i in usable_idx:
            t, p, sigma = prepared[i]  # type: ignore[misc]
            w = 1.0 / (sigma ** 2)
            results[i] = fit_kww_single(
                t, p, weights=w, frozen_beta=shared, target=target
            )

    if all(r is None for r in results):
        return None
    # Return list aligned with input; replace None with NaN-filled dicts so
    # callers can zip safely.
    out: list[dict[str, float]] = []
    for r in results:
        if r is None:
            out.append(
                {
                    "A": float("nan"),
                    "tau_K": float("nan"),
                    "beta": float("nan"),
                    "r_squared": float("nan"),
                    "tau_s": float("nan"),
                }
            )
        else:
            out.append(r)
    return out
