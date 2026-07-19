#!/usr/bin/env python3
"""Background daemon for continuous IR DACF accumulation and dashboard updates."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_ir_from_dipole import CAMPAIGN_LAMBDAS, DEFAULT_CAVITY_FREQUENCY_CM1
from ir_dacf_accumulator import (
    ensemble_mesa_spectrum_from_state,
    refresh_accumulator,
    spectrum_from_mean_acf,
)
import json


def plot_latest_spectra(
    derived_dir: Path,
    *,
    min_freq_cm1: float = 1200.0,
    max_freq_cm1: float = 2200.0,
    method: str = "dct",
    temperature_k: float = 100.0,
) -> None:
    """Plot normalized ensemble IR spectra from accumulated data."""
    state_dir = derived_dir / "accumulator"
    available: list[tuple[float, str, np.ndarray, np.ndarray]] = []
    for lam, tag in CAMPAIGN_LAMBDAS:
        state_path = state_dir / f"lambda{tag}_state.json"
        if not state_path.is_file():
            continue
        state = json.loads(state_path.read_text(encoding="utf-8"))
        mean_acf = np.asarray(state.get("mean_acf", []), dtype=float)
        dt_fs = state.get("dt_fs")
        if mean_acf.size < 4 or dt_fs is None:
            continue
        if method == "mesa":
            freqs, spec = ensemble_mesa_spectrum_from_state(state)
        else:
            freqs, spec = spectrum_from_mean_acf(mean_acf, float(dt_fs))
        if freqs.size == 0:
            continue
        available.append((lam, tag, freqs, spec))

    output_pdf = derived_dir / "ir_spectra_latest.pdf"
    if not available:
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.text(0.5, 0.5, "No accumulated spectra yet", ha="center", va="center")
        fig.savefig(output_pdf, bbox_inches="tight")
        plt.close(fig)
        return

    ncols = len(available)
    fig, axes = plt.subplots(1, ncols, figsize=(3.2 * ncols, 3.2), sharey=True)
    if ncols == 1:
        axes = [axes]
    for ax, (lam, _tag, freqs, spec) in zip(axes, available):
        mask = (freqs >= min_freq_cm1) & (freqs <= max_freq_cm1)
        if not np.any(mask):
            continue
        peak = max(spec[mask].max(), 1e-30)
        ax.plot(freqs[mask], spec[mask] / peak, lw=1.5)
        ax.axvline(
            DEFAULT_CAVITY_FREQUENCY_CM1,
            color="gray",
            ls="--",
            lw=1.0,
        )
        ax.set_title(rf"$\lambda$={lam:g}")
        ax.set_xlim(min_freq_cm1, max_freq_cm1)
        ax.set_ylim(0.0, 1.05)
        ax.set_xlabel(r"frequency (cm$^{-1}$)")
    axes[0].set_ylabel(r"$n(\omega)\alpha(\omega)$ (norm.)")
    fig.suptitle(
        f"Continuous IR accumulation ({method}, ensemble-mean ACF)",
        y=1.02,
        fontsize=10,
    )
    fig.tight_layout()
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path("aging_weak_lambda_ir"),
    )
    parser.add_argument(
        "--derived-dir",
        type=Path,
        default=Path("aging_weak_lambda_ir/derived/continuous"),
    )
    parser.add_argument(
        "--registry",
        type=Path,
        default=Path("aging_weak_lambda_ir/logs/ir_continuous/active_runs.jsonl"),
    )
    parser.add_argument(
        "--poll-interval-s",
        type=float,
        default=120.0,
    )
    parser.add_argument(
        "--method",
        choices=("dct", "mesa"),
        default="dct",
    )
    parser.add_argument(
        "--temperature-k",
        type=float,
        default=100.0,
    )
    parser.add_argument(
        "--once",
        action="store_true",
        help="Run one refresh cycle and exit",
    )
    args = parser.parse_args()

    args.derived_dir.mkdir(parents=True, exist_ok=True)
    log_path = args.derived_dir / "daemon.log"
    while True:
        summary = refresh_accumulator(
            base_dir=args.base_dir,
            registry_path=args.registry,
            derived_dir=args.derived_dir,
            method=args.method,
            temperature_k=args.temperature_k,
        )
        plot_latest_spectra(
            args.derived_dir,
            method=args.method,
            temperature_k=args.temperature_k,
        )
        line = (
            f"{summary['updated_at']} active={len(summary['active_runs'])} "
            f"lambdas={len(summary['lambdas'])}\n"
        )
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write(line)
        print(line.strip())
        if args.once:
            break
        time.sleep(args.poll_interval_s)


if __name__ == "__main__":
    main()
