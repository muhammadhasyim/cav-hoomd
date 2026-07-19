#!/usr/bin/env python3
"""Aggregate smoke-test F(k,t) at |k|=1 a.u. and compare early decay vs legacy data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import sys

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "reproduction" / "shared"))
from relaxation_fit import fit_kww_single  # noqa: E402


def load_fkt_txt(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load two-column lag_time / fskt file, skipping comment and header lines."""
    rows: list[tuple[float, float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) != 2:
            continue
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except ValueError:
            continue
    if not rows:
        raise ValueError(f"no numeric rows in {path}")
    arr = np.asarray(rows, dtype=np.float64)
    return arr[:, 0], arr[:, 1]


def aggregate_ref000(run_dir: Path, n_replicas: int) -> tuple[np.ndarray, np.ndarray, int]:
    """Mean normalized F(k,t) for ref_000 across replicas."""
    curves: list[tuple[np.ndarray, np.ndarray]] = []
    for replica in range(n_replicas):
        path = run_dir / f"prod-{replica}_fkt_ref_000.txt"
        if not path.is_file():
            continue
        t, f = load_fkt_txt(path)
        if f[0] <= 0.0:
            continue
        curves.append((t, f / f[0]))
    if not curves:
        raise FileNotFoundError(f"no ref_000 FKT files under {run_dir}")
    t_grid = curves[0][0]
    stack = np.vstack([np.interp(t_grid, t, fn) for t, fn in curves])
    return t_grid, stack.mean(axis=0), stack.shape[0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--smoke-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "smoke_fkt_ic_init",
    )
    parser.add_argument("--n-replicas", type=int, default=50)
    parser.add_argument(
        "--legacy-master",
        type=Path,
        default=REPO
        / "aging_weak_lambda/lambda0/cavity_coupling_0epos00_switch_200.0ps/master_fskt_ref_000.txt",
        help="Legacy lambda=0 master F(k,t) (|k|=6.02 a.u. campaign reference)",
    )
    args = parser.parse_args()

    run_dir = args.smoke_dir / "lambda0/coupling_000000e+00_switch_200.0ps"
    t, fn, n_used = aggregate_ref000(run_dir, args.n_replicas)

    fit = fit_kww_single(t, fn)
    checkpoints = {int(x): float(np.interp(x, t, fn)) for x in (50, 100, 200, 300, 400)}

    legacy: dict[str, float] | None = None
    if args.legacy_master.is_file():
        tl, fl = load_fkt_txt(args.legacy_master)
        fln = fl / fl[0]
        legacy = {int(x): float(np.interp(x, tl, fln)) for x in (50, 100, 200, 300, 400)}

    summary = {
        "n_replicas_used": n_used,
        "k_au": 1.0,
        "lambda": 0.0,
        "checkpoints_normalized": checkpoints,
        "kww_fit": fit,
        "legacy_k6_master_checkpoints_normalized": legacy,
    }
    out_json = args.smoke_dir / "analysis" / "smoke_fkt_summary.json"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(t, fn, lw=2, label=f"IC-init smoke (k=1.0, n={n_used})")
    if legacy is not None:
        tl, fl = load_fkt_txt(args.legacy_master)
        ax.plot(tl, fl / fl[0], lw=1.5, ls="--", alpha=0.7, label="Legacy master (k=6.02)")
    ax.axvline(100.0, color="gray", ls=":", alpha=0.6)
    ax.axvline(200.0, color="gray", ls="-.", alpha=0.6)
    ax.text(102, 0.05, "100 ps", fontsize=8, color="gray")
    ax.text(202, 0.05, "200 ps", fontsize=8, color="gray")
    ax.set_xlabel("lag time (ps)")
    ax.set_ylabel(r"$F(k,t)/F(k,0)$")
    ax.set_xlim(0, min(600, t[-1]))
    ax.set_ylim(0, 1.05)
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title(r"Smoke $\lambda=0$, $|k|=1.0$ a.u., ref$_0$")
    fig.tight_layout()
    out_png = args.smoke_dir / "analysis" / "smoke_fkt_ref000.png"
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

    print(json.dumps(summary, indent=2))
    print(f"wrote {out_json}")
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
