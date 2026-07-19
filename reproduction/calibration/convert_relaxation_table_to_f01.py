#!/usr/bin/env python3
"""
Convert equilibrium tau(T) from F/F0=1/e to F/F0=0.1 for TN consistency.

The aging material-time analysis defines
    xi(t) = integral dt / tau_0.1(t)
with tau_0.1 from the F(k,t)/F(k,0)=0.1 crossing.

The historical equilibrium table
    average_trajectory/relaxation_times_vs_temperature.txt
uses the 1/e criterion. Feeding that table into RelaxationTimeModel makes
TN material time ~5x steeper than measured.

This script converts each table entry under a stationary KWW assumption
    Phi(t) = exp(-(t/tau_K)^beta)
for which
    tau(F=0.1) / tau(F=1/e) = (-ln(0.1))^(1/beta) / (-ln(1/e))^(1/beta)
                             = (-ln(0.1))^(1/beta)

Optional anchoring rescales the converted table so tau(T_anchor) matches a
measured F=0.1 value (e.g. lambda=0 late-time structural tau), preserving the
relative T dependence.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

F01_THRESHOLD = 0.1
ONE_OVER_E = float(np.e) ** (-1.0)
DEFAULT_BETA = 0.55


def kww_1e_to_f01_factor(beta: float) -> float:
    """Return tau(F=0.1)/tau(F=1/e) for stationary KWW with exponent beta."""
    if beta <= 0.0:
        raise ValueError(f"beta must be positive, got {beta}")
    t_1e = (-np.log(ONE_OVER_E)) ** (1.0 / beta)  # = 1
    t_01 = (-np.log(F01_THRESHOLD)) ** (1.0 / beta)
    return float(t_01 / t_1e)


def convert_tau_1e_to_f01(tau_1e: np.ndarray, beta: float = DEFAULT_BETA) -> np.ndarray:
    """Convert an array of 1/e relaxation times to F=0.1 times."""
    tau_1e = np.asarray(tau_1e, dtype=float)
    return tau_1e * kww_1e_to_f01_factor(beta)


def _parse_data_rows(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if stripped.lower().startswith("temperature"):
                continue
            parts = stripped.split()
            if len(parts) < 3:
                continue
            try:
                float(parts[0])
                float(parts[2])
            except ValueError:
                continue
            rows.append(parts)
    return rows


def convert_table(
    source: Path,
    destination: Path,
    *,
    beta: float = DEFAULT_BETA,
    anchor_temperature_k: float | None = None,
    anchor_tau_ps: float | None = None,
) -> list[tuple[float, float, float]]:
    """
    Convert a relaxation_times_vs_temperature table to F=0.1.

    Returns list of (T_K, inv_T, tau_f01_ps).
    """
    raw_rows = _parse_data_rows(source)
    if not raw_rows:
        raise ValueError(f"No numeric rows found in {source}")

    temperatures = np.array([float(parts[0]) for parts in raw_rows], dtype=float)
    inv_t = np.array([float(parts[1]) for parts in raw_rows], dtype=float)
    tau_1e = np.array([float(parts[2]) for parts in raw_rows], dtype=float)
    tau_f01 = convert_tau_1e_to_f01(tau_1e, beta=beta)

    if anchor_temperature_k is not None or anchor_tau_ps is not None:
        if anchor_temperature_k is None or anchor_tau_ps is None:
            raise ValueError("anchor_temperature_k and anchor_tau_ps must be set together")
        if anchor_tau_ps <= 0.0:
            raise ValueError("anchor_tau_ps must be positive")
        idx = int(np.argmin(np.abs(temperatures - anchor_temperature_k)))
        scale = float(anchor_tau_ps) / float(tau_f01[idx])
        tau_f01 = tau_f01 * scale
        anchor_note = (
            f"# Anchored: tau({temperatures[idx]:.2f} K) set to {anchor_tau_ps:.6f} ps "
            f"(scale={scale:.6f} after KWW conversion)\n"
        )
    else:
        anchor_note = ""

    factor = kww_1e_to_f01_factor(beta)
    destination.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "# Relaxation time analysis from F(k,t) using F/F0 = 0.1 criterion\n"
        f"# Converted from 1/e table: {source.resolve()}\n"
        f"# Conversion: tau_0.1 = tau_1e * (-ln(0.1))^(1/beta) with beta={beta:g} "
        f"(factor={factor:.6f})\n"
        f"{anchor_note}"
        "# Temperature(K) 1/T(1/K) tau_relax(ps) F_initial F_final_norm decay_extent method success\n"
    )
    # RelaxationTimeModel uses skiprows=3, usecols=(0,2). Keep >=3 header lines
    # so the first temperature point is not dropped.
    while header.count("\n") < 3:
        header = "#\n" + header

    with open(destination, "w", encoding="utf-8") as fh:
        fh.write(header)
        for parts, t_k, inv, tau in zip(raw_rows, temperatures, inv_t, tau_f01):
            rest = parts[3:] if len(parts) > 3 else []
            # Rewrite criterion metadata if present
            if rest and rest[-1] in {"True", "False"}:
                # method field often second-to-last
                pass
            line_parts = [
                f"{t_k:.2f}",
                f"{inv:.6f}",
                f"{tau:.6f}",
            ]
            if rest:
                line_parts.extend(rest)
            else:
                line_parts.extend(["0", "0", "1", "converted_f01", "True"])
            fh.write(" ".join(line_parts) + "\n")

    return list(zip(temperatures.tolist(), inv_t.tolist(), tau_f01.tolist()))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        default=Path(__file__).resolve().parents[2]
        / "average_trajectory"
        / "relaxation_times_vs_temperature.txt",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[2]
        / "average_trajectory"
        / "relaxation_times_vs_temperature_f01.txt",
    )
    parser.add_argument("--beta", type=float, default=DEFAULT_BETA)
    parser.add_argument(
        "--anchor-temperature-k",
        type=float,
        default=None,
        help="Optional T at which to match measured F=0.1 tau",
    )
    parser.add_argument(
        "--anchor-tau-ps",
        type=float,
        default=None,
        help="Measured F=0.1 tau (ps) at the anchor temperature",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = convert_table(
        args.source,
        args.output,
        beta=args.beta,
        anchor_temperature_k=args.anchor_temperature_k,
        anchor_tau_ps=args.anchor_tau_ps,
    )
    near_100 = min(rows, key=lambda item: abs(item[0] - 100.0))
    print(f"Wrote {args.output}")
    print(f"  points={len(rows)}")
    print(f"  beta={args.beta:g}  KWW factor={kww_1e_to_f01_factor(args.beta):.4f}")
    print(f"  nearest T={near_100[0]:.2f} K -> tau_0.1={near_100[2]:.3f} ps")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
