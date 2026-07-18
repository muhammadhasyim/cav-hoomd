#!/usr/bin/env python3
"""Average F(k,t), restage, audit, and plot for the |k|=6.02 pilot campaign.

Parameters
----------
Requires complete enough ``prod-*_fkt_ref_*.txt`` files under
``aging_weak_lambda_k6/lambda*/cavity_coupling_*``.  Uses the same
``process_fskt_only.py`` averaging path as the legacy campaign.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
PY = Path("/scratch/mh7373/miniforge3/envs/hoomd/bin/python")
FKT_SCRIPT = REPO / "new_final_offres" / "process_fskt_only.py"
BASE = REPO / "aging_weak_lambda_k6"
ANALYSIS = BASE / "analysis"
STAGED = BASE / "derived" / "repro_layout"
FIGURES = BASE / "derived" / "repro_figures"
AUDIT_SCRIPT = REPO / "aging_weak_lambda" / "analysis" / "structural_relaxation_audit.py"
PAPER_KMAG = 6.02
PAPER_ANCHOR = 2.35
SUCCESS_MAX_TILDE = 2.6  # near/below paper lambda=0.042 peak


def _coupling_dirs() -> list[Path]:
    dirs = []
    for pattern in (
        "lambda0/cavity_coupling_*_switch_200.0ps",
        "lambda0p03/cavity_coupling_*_switch_200.0ps",
    ):
        matches = sorted(BASE.glob(pattern))
        if matches:
            dirs.append(matches[0])
    return dirs


def count_complete_h5(run_dir: Path, min_bytes: int = 1_380_000_000) -> int:
    """Count scientifically complete observables HDF5 files."""
    return sum(
        1
        for path in run_dir.glob("observables_replica_*.h5")
        if path.is_file() and path.stat().st_size >= min_bytes
    )


def average_fkt(run_dir: Path, max_time: float = 2500.0, dt: float = 1.0) -> None:
    """Rebuild master F(k,t) averages for one coupling directory."""
    for path in run_dir.glob("master_fskt_ref_*.txt"):
        path.unlink()
    subprocess.run(
        [
            str(PY),
            str(FKT_SCRIPT),
            "--exp_dir",
            str(run_dir),
            "--fskt_dt",
            str(dt),
            "--max_time",
            str(max_time),
        ],
        check=True,
    )
    counts_dir = run_dir / "fkt_sample_counts"
    counts_dir.mkdir(exist_ok=True)
    for path in run_dir.glob("master_fskt_ref_*_sample_counts.txt"):
        path.rename(counts_dir / path.name)


def restage_and_plot() -> None:
    """Build staged layout and regenerate the measured tau figure."""
    FIGURES.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            str(PY),
            str(REPO / "reproduction" / "adapters" / "build_aging_repro_layout.py"),
            "--profile",
            "aging_weak_lambda_k6",
            "--skip-material-time",
        ],
        check=True,
        cwd=str(REPO / "reproduction"),
    )
    subprocess.run(
        [
            str(PY),
            str(REPO / "reproduction" / "shared" / "plot_fkt_analysis.py"),
            "--profile",
            "aging_weak_lambda_k6",
            "--output_dir",
            str(FIGURES),
            "--fkt-kmag",
            f"{PAPER_KMAG:g}",
        ],
        check=True,
        cwd=str(REPO / "reproduction"),
    )


def run_audit() -> Path:
    """Write structural-relaxation audit for the k=6.02 staged masters."""
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            str(PY),
            str(AUDIT_SCRIPT),
            "--staged-root",
            str(STAGED),
            "--output-dir",
            str(ANALYSIS),
            "--fkt-kmag",
            f"{PAPER_KMAG:g}",
            "--coupling-cli",
            "0.03",
        ],
        check=True,
    )
    return ANALYSIS / "structural_relaxation_audit_report.txt"


def compare_to_paper(report_path: Path) -> dict[str, float]:
    """Parse peak tilde(tau) at t_w=0 for lambda=0.03 from the audit report."""
    text = report_path.read_text(encoding="utf-8")
    peaks: dict[str, float] = {}
    for line in text.splitlines():
        if "raw_per_ref_global_baseline" in line and "lambda=0.03" in line:
            # format: ... lambda=0.03: 3.770
            for part in line.split(","):
                if "lambda=0.03" in part:
                    peaks["raw_per_ref"] = float(part.split(":")[-1].strip())
        if "raw_ref0_global_baseline" in line and "lambda=0.03" in line:
            for part in line.split(","):
                if "lambda=0.03" in part:
                    peaks["raw_ref0"] = float(part.split(":")[-1].strip())
    return peaks


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--min-replicas",
        type=int,
        default=20,
        help="Minimum complete HDF5 replicas per coupling before averaging",
    )
    parser.add_argument(
        "--skip-average",
        action="store_true",
        help="Reuse existing master_fskt files",
    )
    args = parser.parse_args(argv)

    dirs = _coupling_dirs()
    if len(dirs) < 2:
        print("ERROR: expected lambda0 and lambda0p03 run directories", file=sys.stderr)
        return 2

    for run_dir in dirs:
        n = count_complete_h5(run_dir)
        print(f"{run_dir}: complete_h5={n}")
        if n < args.min_replicas:
            print(
                f"ERROR: need >= {args.min_replicas} complete replicas in {run_dir}",
                file=sys.stderr,
            )
            return 3
        if not args.skip_average:
            print(f"Averaging F(k,t) in {run_dir} ...")
            average_fkt(run_dir)

    print("Restaging and plotting ...")
    restage_and_plot()
    report = run_audit()
    peaks = compare_to_paper(report)
    summary = ANALYSIS / "k6_pilot_vs_paper.txt"
    lines = [
        "k=6.02 pilot comparison to paper Fig. 2b",
        f"paper_anchor_lambda=0.042 tilde~{PAPER_ANCHOR}",
        f"success_criterion: peak tilde(lambda=0.03, tw=0) <= {SUCCESS_MAX_TILDE}",
        f"audit_report={report}",
    ]
    for key, value in sorted(peaks.items()):
        ok = value <= SUCCESS_MAX_TILDE
        lines.append(f"{key}_peak_tilde_lambda0p03={value:.3f} pass={ok}")
    if not peaks:
        lines.append("WARNING: could not parse peak tilde from audit report")
    summary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(summary.read_text(encoding="utf-8"))

    # Copy report into figures for convenience.
    FIGURES.mkdir(parents=True, exist_ok=True)
    shutil.copy2(report, FIGURES / report.name)
    shutil.copy2(summary, FIGURES / summary.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
