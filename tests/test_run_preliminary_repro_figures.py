"""Tests that preliminary figure runner mirrors run_aging_weak_lambda.sh."""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "aging_weak_lambda" / "analysis"))

from run_preliminary_repro_figures import (  # noqa: E402
    DEFAULT_STRUCTURAL_ENERGY,
    FIGURE3_DISPLAY_SMOOTH_WINDOW,
    FIGURE3_LAMBDAS,
    PROFILE,
    RF_CALIBRATION_PATH,
    STRUCTURAL_ENERGY_COMPONENT,
    build_energy_pipeline_commands,
    build_figure_pipeline_steps,
)


def _script_names(steps: list[tuple[str, list[str]]]) -> list[str]:
    names: list[str] = []
    for _label, argv in steps:
        # argv[0]=python, argv[1]=script path
        names.append(Path(argv[1]).name)
    return names


def test_figure_pipeline_matches_repro_shell_step_scripts() -> None:
    """Step scripts must match reproduction/run_aging_weak_lambda.sh [2/7]–[6/7]."""
    staged = [
        Path("/tmp/prelim_layout/cavity_coupling_0epos00_switch_200.0ps"),
        Path("/tmp/prelim_layout/cavity_coupling_1eneg02_switch_200.0ps"),
    ]
    steps = build_figure_pipeline_steps(
        python_executable=Path("/usr/bin/python"),
        figures_dir=Path("/tmp/prelim_figures"),
        profile=PROFILE,
        fkt_kmag="1.0",
        max_time_ps=2500.0,
        figure3_lambdas=FIGURE3_LAMBDAS,
        staged_coupling_dirs=staged,
    )
    scripts = _script_names(steps)
    # RF masters live under staged_root, so step [2/7] is one process_fskt_only
    # per staged coupling dir (not --profile -> data_root run dirs).
    assert scripts[:2] == ["process_fskt_only.py", "process_fskt_only.py"]
    assert scripts[2:] == [
        "plot_fkt_analysis.py",
        "plot_fictive_three_panel_analysis.py",
        "plot_figure3.py",
        "plot_figure3.py",
        "plot_figure3.py",
        "plot_figure3.py",
        "plot_figure3.py",
        "run_corrected_analysis.py",
        "plot_material_time_and_collapse.py",
        "plot_figure4.py",
        "export_figure3_csv.py",
        "export_figure4_csv.py",
    ]
    for _label, argv in steps:
        if Path(argv[1]).name == "plot_figure3.py":
            assert "--smooth-window" in argv
            sw_idx = argv.index("--smooth-window")
            assert argv[sw_idx + 1] == str(FIGURE3_DISPLAY_SMOOTH_WINDOW)
    for _label, argv in steps[:2]:
        assert "--exp_dir" in argv
        assert "--profile" not in argv
        assert any(str(path) in argv for path in staged)


def test_figure_pipeline_uses_preliminary_profile_and_output_dir() -> None:
    figures = Path("/tmp/prelim_figures")
    steps = build_figure_pipeline_steps(
        python_executable=Path("/usr/bin/python"),
        figures_dir=figures,
        profile=PROFILE,
        fkt_kmag="1.0",
        max_time_ps=2500.0,
        figure3_lambdas=FIGURE3_LAMBDAS,
        staged_coupling_dirs=[Path("/tmp/prelim_layout/cavity_coupling_0epos00_switch_200.0ps")],
    )
    assert PROFILE == "aging_weak_lambda_preliminary"
    for _label, argv in steps:
        joined = " ".join(argv)
        # process_fskt_only uses staged --exp_dir; later steps use --profile.
        if Path(argv[1]).name == "process_fskt_only.py":
            assert "--exp_dir" in argv
            assert "/tmp/prelim_layout/" in joined
            continue
        assert f"--profile {PROFILE}" in joined or PROFILE in argv
        if "--output" in argv or "--output_dir" in argv or "--material-time-output" in argv:
            assert str(figures) in joined


def test_preliminary_uses_rf_calibration_lj_channel() -> None:
    """Map structural fictive T with LJ column of the RF PE-vs-T table."""
    assert RF_CALIBRATION_PATH.name == "potential_energy_vs_T_rf.txt"
    assert "pe_vs_T_calib_rf" in str(RF_CALIBRATION_PATH)
    assert DEFAULT_STRUCTURAL_ENERGY == "lj_only"
    assert STRUCTURAL_ENERGY_COMPONENT == "lj"

    commands = build_energy_pipeline_commands(
        python_executable=Path("/usr/bin/python"),
        analysis_dir=Path("/tmp/prelim_analysis"),
        empirical_data=RF_CALIBRATION_PATH,
        structural_energy_component=STRUCTURAL_ENERGY_COMPONENT,
    )
    assert len(commands) == 1
    argv = commands[0]
    assert "energies_to_fictive_temperatures.py" in argv[1]
    assert str(RF_CALIBRATION_PATH) in argv
    assert "--structural-energy-component" in argv
    comp_idx = argv.index("--structural-energy-component")
    assert argv[comp_idx + 1] == "lj"
    # Must not keep the old PPPM/Ewald table or inverted RF lj+coulombic channel.
    assert "reproduction/calibration/potential_energy_vs_T.txt" not in " ".join(argv)
    assert argv[comp_idx + 1] != "lj_coulombic"


def test_figure_pipeline_includes_plot_figure4_and_csv_exports() -> None:
    steps = build_figure_pipeline_steps(
        python_executable=Path("/usr/bin/python"),
        figures_dir=Path("/tmp/prelim_figures"),
        profile=PROFILE,
        fkt_kmag="1.0",
        max_time_ps=2500.0,
        figure3_lambdas=FIGURE3_LAMBDAS,
        staged_coupling_dirs=[Path("/tmp/prelim_layout/cavity_coupling_0epos00_switch_200.0ps")],
    )
    labels = [label for label, _ in steps]
    assert any("5b" in label or "figure4" in label.lower() for label in labels)
    assert any("export" in label.lower() or "csv" in label.lower() for label in labels)
    scripts = set(_script_names(steps))
    assert "plot_figure4.py" in scripts
    assert "export_figure4_csv.py" in scripts
    assert "export_figure3_csv.py" in scripts
