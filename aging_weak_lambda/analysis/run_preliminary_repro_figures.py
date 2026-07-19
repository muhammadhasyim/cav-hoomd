#!/usr/bin/env python3
"""Generate preliminary reproduction figures from RF campaign replicas.

Mirrors ``reproduction/run_aging_weak_lambda.sh`` steps [2/7]–[6/7] with
``--profile aging_weak_lambda_preliminary`` so all outputs land under
``aging_weak_lambda/derived/preliminary/``.

RF-specific input prep (steps before plotting):
  - select scientifically complete replicas (e.g. ``--runtime-ps 1600`` /
    ``--max-references 8``)
  - stage F(k,t) masters and average energies from RF HDF5
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
REPRO_ROOT = REPO_ROOT / "reproduction"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from aging_weak_lambda.analysis.preliminary_repro_layout import (  # noqa: E402
    average_campaign_energies,
    build_preliminary_layout,
)

DEFAULT_PY = Path("/scratch/mh7373/miniforge3/envs/hoomd/bin/python")
PRELIM_ROOT = REPO_ROOT / "aging_weak_lambda" / "derived" / "preliminary"
STAGED_ROOT = PRELIM_ROOT / "layout"
ANALYSIS_DIR = PRELIM_ROOT / "analysis"
FIGURES_DIR = PRELIM_ROOT / "figures"
PROFILE = "aging_weak_lambda_preliminary"
FIGURE3_LAMBDAS = (0.0, 0.01, 0.016667, 0.023333, 0.03)
FIGURE3_DISPLAY_SMOOTH_WINDOW = 7
DEFAULT_FKT_KMAG = "1.0"
DEFAULT_MAX_TIME_PS = 2500.0
RF_CALIBRATION_PATH = (
    REPO_ROOT
    / "aging_weak_lambda"
    / "analysis"
    / "pe_vs_T_calib_rf"
    / "potential_energy_vs_T_rf.txt"
)
# HDF5 averaging: LJ only (RF Coulomb is noisy/inverted; keep Fig 3 LJ+Coulomb label).
DEFAULT_STRUCTURAL_ENERGY = "lj_only"
# EmpiricalTemperatureData channel: lj_hartree from RF PE-vs-T table.
STRUCTURAL_ENERGY_COMPONENT = "lj"


def _lambda_tag(lam: float) -> str:
    """Match ``run_aging_weak_lambda.sh``: ``echo $lam | tr '.' 'p'`` (lam=0 -> ``0``)."""
    if lam == 0.0:
        return "0"
    return str(lam).replace(".", "p")


def stage_time_series_output() -> None:
    """Write paper-format time series files from preliminary analysis CSVs."""
    sys.path.insert(0, str(REPRO_ROOT / "shared"))
    sys.path.insert(0, str(REPRO_ROOT / "adapters"))
    from build_aging_repro_layout import _write_time_series_output  # noqa: WPS433
    from dataset_profile import load_profile  # noqa: WPS433

    profile = load_profile(PROFILE)
    _write_time_series_output(profile)


def build_energy_pipeline_commands(
    *,
    python_executable: Path,
    analysis_dir: Path = ANALYSIS_DIR,
    empirical_data: Path = RF_CALIBRATION_PATH,
    structural_energy_component: str = STRUCTURAL_ENERGY_COMPONENT,
) -> list[list[str]]:
    """Return argv list(s) for energy -> fictive conversion (RF calibration)."""
    return [
        [
            str(python_executable),
            str(REPO_ROOT / "aging_weak_lambda" / "analysis" / "energies_to_fictive_temperatures.py"),
            "--analysis-dir",
            str(analysis_dir),
            "--empirical-data",
            str(empirical_data),
            "--structural-energy-component",
            structural_energy_component,
        ]
    ]


def run_energy_pipeline(
    python_executable: Path,
    output_base: Path,
    n_replicas: int,
    *,
    runtime_ps: float | None = None,
    max_references: int | None = None,
    use_all_available: bool = True,
    structural_energy: str = DEFAULT_STRUCTURAL_ENERGY,
    empirical_data: Path = RF_CALIBRATION_PATH,
    structural_energy_component: str = STRUCTURAL_ENERGY_COMPONENT,
) -> None:
    """Average energies, convert to fictive temperatures, and stage time series."""
    if not empirical_data.is_file():
        raise FileNotFoundError(
            f"RF calibration table missing: {empirical_data}. "
            "Run aging_weak_lambda/analysis/run_rf_pe_vs_T_calibration.sh first."
        )
    average_campaign_energies(
        output_base,
        ANALYSIS_DIR,
        n_replicas=n_replicas,
        runtime_ps=runtime_ps,
        max_references=max_references,
        use_all_available=use_all_available,
        structural_energy=structural_energy,
    )
    for argv in build_energy_pipeline_commands(
        python_executable=python_executable,
        analysis_dir=ANALYSIS_DIR,
        empirical_data=empirical_data,
        structural_energy_component=structural_energy_component,
    ):
        subprocess.run(argv, check=True)
    stage_time_series_output()


def build_figure_pipeline_steps(
    *,
    python_executable: Path,
    figures_dir: Path,
    profile: str = PROFILE,
    fkt_kmag: str = DEFAULT_FKT_KMAG,
    max_time_ps: float = DEFAULT_MAX_TIME_PS,
    figure3_lambdas: tuple[float, ...] = FIGURE3_LAMBDAS,
    figure3_smooth_window: int = FIGURE3_DISPLAY_SMOOTH_WINDOW,
    staged_coupling_dirs: list[Path] | None = None,
) -> list[tuple[str, list[str]]]:
    """Return (label, argv) steps matching ``run_aging_weak_lambda.sh`` [2/7]–[6/7].

    Parameters
    ----------
    python_executable
        Interpreter used to invoke reproduction scripts.
    figures_dir
        Output directory for PDFs/PNGs/CSVs (preliminary figures root).
    profile
        Dataset profile name (must be ``aging_weak_lambda_preliminary`` for RF).
    fkt_kmag
        Wavevector magnitude for F(k,t) plots (RF campaign uses ``1.0``).
    max_time_ps
        ``process_fskt_only.py --max_time`` (repro shell uses 2500).
    figure3_lambdas
        Couplings for Figure 3 panels.
    staged_coupling_dirs
        Staged coupling directories that already contain ``master_fskt_*.txt``.
        Required for RF preliminary data: ``process_fskt_only --profile``
        resolves ``data_root`` run dirs, which do not hold RF masters.
    """
    if staged_coupling_dirs is None:
        raise ValueError(
            "staged_coupling_dirs is required so step [2/7] uses staged masters "
            "(RF preliminary), not raw data_root run directories"
        )

    py = str(python_executable)
    figs = figures_dir
    steps: list[tuple[str, list[str]]] = []
    for staged_dir in staged_coupling_dirs:
        steps.append(
            (
                f"[2/7] Fig 2b F(k,t) process (staged {staged_dir.name})",
                [
                    py,
                    str(REPRO_ROOT / "figure2_fkt" / "process_fskt_only.py"),
                    "--exp_dir",
                    str(staged_dir),
                    "--skip_processing",
                    "--max_time",
                    str(max_time_ps),
                ],
            )
        )
    steps.extend(
        [
            (
                "[3/7] Fig 2b/c F(k,t) + measured relaxation",
                [
                    py,
                    str(REPRO_ROOT / "shared" / "plot_fkt_analysis.py"),
                    "--profile",
                    profile,
                    "--output_dir",
                    str(figs),
                    "--fkt-kmag",
                    fkt_kmag,
                ],
            ),
            (
                "[3b/7] Fig 2b/c fictive (TN) relaxation panels",
                [
                    py,
                    str(REPRO_ROOT / "figure2_3_relaxation" / "plot_fictive_three_panel_analysis.py"),
                    "--profile",
                    profile,
                    "--no-latex",
                    "--relaxation-output",
                    str(figs / "tn_model_relaxation_analysis.pdf"),
                    "--material-time-output",
                    str(figs / "fictive_material_time_evolution.pdf"),
                ],
            ),
        ]
    )
    for lam in figure3_lambdas:
        tag = _lambda_tag(lam)
        steps.append(
            (
                f"[4/7] Fig 3 energy/temperature lambda={lam:g}",
                [
                    py,
                    str(REPRO_ROOT / "figure3_energy" / "plot_figure3.py"),
                    "--profile",
                    profile,
                    "--coupling",
                    str(lam),
                    "--no-latex",
                    "--smooth-window",
                    str(figure3_smooth_window),
                    "--smooth-window-total",
                    str(figure3_smooth_window),
                    "--output",
                    str(figs / f"figure3_lambda_{tag}"),
                ],
            )
        )
    steps.extend(
        [
            (
                "[5/7] Fig 4 material-time corrected analysis",
                [
                    py,
                    str(REPRO_ROOT / "figure4_material_time" / "run_corrected_analysis.py"),
                    "--profile",
                    profile,
                ],
            ),
            (
                "[5/7] Fig 4 material_time_and_collapse",
                [
                    py,
                    str(REPRO_ROOT / "figure4_material_time" / "plot_material_time_and_collapse.py"),
                    "--profile",
                    profile,
                    "--no-latex",
                    "--output",
                    str(figs / "material_time_and_collapse.pdf"),
                ],
            ),
            (
                "[5b/7] Fig 4 four-panel composite",
                [
                    py,
                    str(REPRO_ROOT / "figure4_material_time" / "plot_figure4.py"),
                    "--profile",
                    profile,
                    "--output",
                    str(figs / "figure4.pdf"),
                ],
            ),
            (
                "[6/7] export figure3 CSV",
                [
                    py,
                    str(REPRO_ROOT / "csv_export" / "export_figure3_csv.py"),
                    "--profile",
                    profile,
                    "--coupling",
                    "0.03",
                    "--output",
                    str(figs / "figure3_lambda0p03.csv"),
                ],
            ),
            (
                "[6/7] export figure4 CSV",
                [
                    py,
                    str(REPRO_ROOT / "csv_export" / "export_figure4_csv.py"),
                    "--profile",
                    profile,
                    "--output",
                    str(figs / "figure4.csv"),
                ],
            ),
        ]
    )
    return steps


def _load_staged_coupling_dirs(profile_name: str = PROFILE) -> list[Path]:
    """Resolve staged coupling directories from the preliminary dataset profile."""
    sys.path.insert(0, str(REPRO_ROOT / "shared"))
    from dataset_profile import load_profile, staged_coupling_dirs  # noqa: WPS433

    return list(staged_coupling_dirs(load_profile(profile_name)))


def run_figure_pipeline(
    python_executable: Path,
    *,
    figures_dir: Path = FIGURES_DIR,
    fkt_kmag: str = DEFAULT_FKT_KMAG,
    max_time_ps: float = DEFAULT_MAX_TIME_PS,
) -> None:
    """Execute the full repro figure/CSV step list into ``figures_dir``."""
    figures_dir.mkdir(parents=True, exist_ok=True)
    staged_dirs = _load_staged_coupling_dirs(PROFILE)
    for label, argv in build_figure_pipeline_steps(
        python_executable=python_executable,
        figures_dir=figures_dir,
        profile=PROFILE,
        fkt_kmag=fkt_kmag,
        max_time_ps=max_time_ps,
        figure3_lambdas=FIGURE3_LAMBDAS,
        staged_coupling_dirs=staged_dirs,
    ):
        print(f"=== {label} ===", flush=True)
        subprocess.run(argv, check=True, cwd=str(REPRO_ROOT))


def write_summary(summary: dict[float, dict[str, int | str]]) -> None:
    """Write a small JSON/text summary of replica counts used."""
    PRELIM_ROOT.mkdir(parents=True, exist_ok=True)
    payload = {
        str(lam): values for lam, values in sorted(summary.items(), key=lambda item: item[0])
    }
    (PRELIM_ROOT / "replica_summary.json").write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    lines = [
        "Preliminary reproduction figure inputs",
        f"figures_dir={FIGURES_DIR}",
        f"staged_root={STAGED_ROOT}",
        f"analysis_dir={ANALYSIS_DIR}",
        f"profile={PROFILE}",
        "pipeline=run_aging_weak_lambda.sh steps [2/7]-[6/7]",
        "",
    ]
    total = 0
    for lam, values in sorted(summary.items(), key=lambda item: item[0]):
        count = int(values["complete_replicas"])
        total += count
        lines.append(
            f"lambda={lam:g}: complete_replicas={count} staged_dir={values['staged_dir']}"
        )
    lines.append(f"total_complete_replicas={total}")
    (PRELIM_ROOT / "replica_summary.txt").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )
    (PRELIM_ROOT / "structural_energy_proxy.txt").write_text(
        "lj_coul_ha column stores LJ only (RF Coulomb omitted: noisy and U(T)-inverted); "
        "Figure 3 still labels the structural channel as LJ+Coulomb. "
        f"Fictive mapping uses {RF_CALIBRATION_PATH.name} "
        f"via energy_component={STRUCTURAL_ENERGY_COMPONENT} (lj_hartree column). "
        "Kinetic temperature is temperatures/kinetic from HDF5 (Kelvin). "
        "Aging uses replica means at each lab time (no temporal block averaging). "
        f"Figure 3 applies display-only {FIGURE3_DISPLAY_SMOOTH_WINDOW}-sample moving "
        "average to LJ+Coulomb/Total curves; staged science CSVs stay unsmoothed.\n",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-base",
        type=Path,
        default=REPO_ROOT / "aging_weak_lambda",
    )
    parser.add_argument(
        "--python-executable",
        type=Path,
        default=DEFAULT_PY if DEFAULT_PY.is_file() else Path(sys.executable),
    )
    parser.add_argument("--n-replicas", type=int, default=500)
    parser.add_argument(
        "--runtime-ps",
        type=float,
        default=None,
        help=(
            "Completeness runtime for selecting replicas "
            "(default: campaign target, currently 2000 ps). "
            "Use 1600 while only primary-complete replicas exist."
        ),
    )
    parser.add_argument(
        "--max-references",
        type=int,
        default=None,
        help=(
            "Required F(k,t) reference count for completeness "
            "(default: campaign target, currently 10). Use 8 with --runtime-ps 1600. "
            "Ignored when --use-all-available is set (default)."
        ),
    )
    parser.add_argument(
        "--use-all-available",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Average the union of primary-complete (1600 ps / 8 refs) and "
            "near-target-complete (~1999 ps / 10 refs) replicas. Late times keep "
            "fewer contributors instead of dropping longer runs. Default: on."
        ),
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--skip-fkt-layout",
        action="store_true",
        help="Reuse existing staged F(k,t) masters",
    )
    parser.add_argument(
        "--skip-energy-pipeline",
        action="store_true",
        help="Skip energy averaging and fictive-temperature staging",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Only rebuild staged inputs (F(k,t) and/or energies)",
    )
    parser.add_argument(
        "--fkt-kmag",
        default=DEFAULT_FKT_KMAG,
        help="Wavevector magnitude for F(k,t) plots (RF campaign default: 1.0)",
    )
    parser.add_argument(
        "--max-time",
        type=float,
        default=DEFAULT_MAX_TIME_PS,
        help="process_fskt_only --max_time (default: 2500, same as repro shell)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.n_replicas <= 0:
        raise SystemExit("--n-replicas must be positive")

    if args.skip_fkt_layout:
        summary: dict[float, dict[str, int | str]] = {}
    else:
        print("=== [1/7] RF staged F(k,t) layout (complete replicas) ===", flush=True)
        summary = build_preliminary_layout(
            args.output_base,
            STAGED_ROOT,
            python_executable=args.python_executable,
            n_replicas=args.n_replicas,
            dry_run=args.dry_run,
            runtime_ps=args.runtime_ps,
            max_references=args.max_references,
            use_all_available=args.use_all_available,
        )
        write_summary(summary)
        print(json.dumps(summary, indent=2, sort_keys=True))

    if args.dry_run:
        return 0

    if not args.skip_energy_pipeline:
        print("=== [1b/7] RF energy average + fictive staging ===", flush=True)
        run_energy_pipeline(
            args.python_executable,
            args.output_base,
            args.n_replicas,
            runtime_ps=args.runtime_ps,
            max_references=args.max_references,
            use_all_available=args.use_all_available,
        )

    if args.skip_plots:
        return 0

    run_figure_pipeline(
        args.python_executable,
        figures_dir=FIGURES_DIR,
        fkt_kmag=args.fkt_kmag,
        max_time_ps=args.max_time,
    )
    print(f"=== [7/7] Done ===\nPreliminary figures written under: {FIGURES_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
