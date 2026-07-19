#!/usr/bin/env python3
"""
Assemble Figure 4: four-panel aging / material-time analysis.

Uses the same data path and styling as the smooth TN figures:
  - Panel 1: measured xi(t) + TN xi_TN(t) from staged masters / fictive series
  - Panel 2: ISF collapse Phi_k(h) with Phi = exp(-h^beta), beta=0.55
  - Panel 3/4: measured tau_tilde_s from relaxation_vs_reference_time_data.txt

No subplot letter labels. Computer Modern fonts via setup_latex_fonts().
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable

REPRO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPRO_ROOT / "shared"))
sys.path.insert(0, str(REPRO_ROOT / "figure2_3_relaxation"))
sys.path.insert(0, str(REPRO_ROOT / "figure4_material_time"))

from latex_config_adobe import latex_safe  # noqa: E402
from plot_fictive_three_panel_analysis import FictiveThreePanelAnalyzer  # noqa: E402
from run_corrected_analysis import CorrectedMaterialTimeAnalyzer  # noqa: E402

# Paper thermal-aging value is 0.55; weak-lambda collapse overlays near ~0.72
# when using the F=0.1 / ln(10) material-time convention (see stretched_exponential_guide).
DEFAULT_BETA = 0.72
PAPER_BETA = 0.55
DEFAULT_WAITING_TIMES_PS = [0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]


def _register_matplotlib_cm_fonts() -> str:
    """Register matplotlib-bundled Computer Modern TTFs; return family name."""
    from matplotlib import font_manager

    ttf_dir = Path(matplotlib.__file__).resolve().parent / "mpl-data" / "fonts" / "ttf"
    roman = ttf_dir / "cmr10.ttf"
    if not roman.is_file():
        raise FileNotFoundError(f"matplotlib cmr10.ttf not found at {roman}")
    for fname in ("cmr10.ttf", "cmmi10.ttf", "cmsy10.ttf", "cmtt10.ttf"):
        path = ttf_dir / fname
        if path.is_file():
            font_manager.fontManager.addfont(str(path))
    return "cmr10"


def setup_figure4_fonts() -> None:
    """
    Computer Modern via matplotlib cmr10 + mathtext CM (same path as figure3).

    Do not call setup_latex_fonts()/usetex here: latex.fmt is missing on this
    cluster and the fallback silently becomes DejaVu sans.
    """
    mpl.style.use("classic")
    try:
        cm_family = _register_matplotlib_cm_fonts()
    except (OSError, FileNotFoundError):
        cm_family = "DejaVu Serif"
    plt.rcParams.update(
        {
            "text.usetex": False,
            "font.family": cm_family,
            "font.serif": [cm_family, "Computer Modern Roman", "DejaVu Serif"],
            "mathtext.fontset": "cm",
            "axes.formatter.use_mathtext": True,
            "font.size": 14,
            "axes.linewidth": 1.2,
            "axes.unicode_minus": False,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "legend.framealpha": 0.95,
            "axes.axisbelow": True,
            "axes.grid": True,
            "grid.alpha": 0.3,
        }
    )


def stretched_exponential_guide(h: np.ndarray, beta: float) -> np.ndarray:
    """
    Guide curve consistent with F/F0 = 0.1 material-time definition.

    Material time h is built from the F=0.1 crossing, so Phi(h=1) = 0.1.
    That requires Phi(h) = exp(-ln(10) * h^beta), not exp(-h^beta).
    The paper writes e^{-h^beta}; the ln(10) is absorbed into the operational
    definition of h used in this codebase (same as plot_material_time_and_collapse).
    """
    h = np.asarray(h, dtype=float)
    return np.exp(-np.log(10.0) * np.power(np.maximum(h, 0.0), beta))


def load_panel_b_csv(path: Path) -> list[dict[str, np.ndarray]]:
    """Load collapse trajectories from figure4_panel_b.csv if present."""
    if not path.is_file():
        return []
    trajectories: list[dict[str, object]] = []
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(row for row in fh if not row.startswith("#"))
        for row in reader:
            if row["series"] != "collapse_data":
                continue
            trajectories.append(
                {
                    "lambda_au": float(row["lambda_au"]),
                    "t_w_ps": float(row["t_w_ps"]),
                    "ref_num": int(float(row["ref_num"])),
                    "h": float(row["h"]),
                    "phi": float(row["Phi_k"]),
                }
            )
    grouped: dict[tuple[float, float, int], dict[str, list[float]]] = defaultdict(
        lambda: {"h": [], "phi": []}
    )
    meta: dict[tuple[float, float, int], float] = {}
    for item in trajectories:
        key = (item["lambda_au"], item["t_w_ps"], item["ref_num"])
        grouped[key]["h"].append(item["h"])
        grouped[key]["phi"].append(item["phi"])
        meta[key] = item["lambda_au"]
    merged: list[dict[str, np.ndarray]] = []
    for key, values in grouped.items():
        order = np.argsort(values["h"])
        merged.append(
            {
                "lambda_au": meta[key],
                "h": np.asarray(values["h"], dtype=float)[order],
                "phi": np.asarray(values["phi"], dtype=float)[order],
            }
        )
    return merged


def compute_collapse_from_staged(
    staged_root: Path,
    coupling_tags: dict[float, str],
    switch_time_ps: float = 200.0,
) -> list[dict[str, np.ndarray]]:
    """Compute ISF collapse curves from master F(k,t) using corrected material time."""
    analyzer = CorrectedMaterialTimeAnalyzer(
        criterion_value=0.1,
        alpha_smooth=0.01,
        n_grid_points=500,
        verbose=False,
    )
    collapsed: list[dict[str, np.ndarray]] = []
    for lam, tag in sorted(coupling_tags.items()):
        coupling_dir = staged_root / f"cavity_coupling_{tag}_switch_{switch_time_ps}ps"
        if not coupling_dir.is_dir():
            continue
        try:
            fkt_files = analyzer.load_fkt_data(str(coupling_dir))
            waiting_times = analyzer.calculate_waiting_times(fkt_files)
            normalized_fkt = analyzer.normalize_fkt_data(fkt_files)
            time_grid, xi_grid = analyzer.calculate_material_time(normalized_fkt)
            curves = analyzer.collapse_data(
                (time_grid, xi_grid), normalized_fkt, waiting_times
            )
            for ref_num, data in curves.items():
                collapsed.append(
                    {
                        "lambda_au": lam,
                        "h": np.asarray(data["xi"], dtype=float),
                        "phi": np.asarray(data["R_fkt"], dtype=float),
                    }
                )
        except Exception as exc:  # noqa: BLE001
            print(f"  Warning: collapse failed for lambda={lam}: {exc}")
    return collapsed


def build_measured_tau_tilde_tables(
    relaxation_rows: list[tuple[float, int, float, float]],
) -> tuple[dict[float, dict[str, np.ndarray]], dict[float, dict[str, np.ndarray]]]:
    """Build measured tau_tilde tables for panels (c) and (d).

    Uses the same F(k,t) relaxation table as ``export_figure4_csv.build_panels_cd``.
    """
    tau_dict: dict[float, dict[int, float]] = {}
    ref_times: dict[int, float] = {}
    for coupling_value, ref_num, ref_time_ps, relaxation_time in relaxation_rows:
        tau_dict.setdefault(coupling_value, {})[ref_num] = relaxation_time
        ref_times[ref_num] = ref_time_ps

    tau0 = tau_dict.get(0.0, {})
    panel_c: dict[float, dict[str, np.ndarray]] = defaultdict(
        lambda: {"lambda": [], "tau_tilde": []}
    )
    panel_d: dict[float, dict[str, np.ndarray]] = defaultdict(
        lambda: {"t_w": [], "tau_tilde": []}
    )

    for coupling_value, ref_num, ref_time_ps, relaxation_time in relaxation_rows:
        if ref_num == 0:
            continue
        tau_ref = tau0.get(ref_num)
        if tau_ref is None or tau_ref <= 0:
            continue
        tau_tilde = relaxation_time / tau_ref
        t_w = ref_time_ps - 200.0
        panel_c[t_w]["lambda"].append(coupling_value)
        panel_c[t_w]["tau_tilde"].append(tau_tilde)
        panel_d[coupling_value]["t_w"].append(t_w)
        panel_d[coupling_value]["tau_tilde"].append(tau_tilde)

    panel_c_out: dict[float, dict[str, np.ndarray]] = {}
    for tw, data in panel_c.items():
        order = np.argsort(data["lambda"])
        panel_c_out[float(tw)] = {
            "lambda": np.asarray(data["lambda"], dtype=float)[order],
            "tau_tilde": np.asarray(data["tau_tilde"], dtype=float)[order],
        }

    panel_d_out: dict[float, dict[str, np.ndarray]] = {}
    for lam, data in panel_d.items():
        order = np.argsort(data["t_w"])
        panel_d_out[float(lam)] = {
            "t_w": np.asarray(data["t_w"], dtype=float)[order],
            "tau_tilde": np.asarray(data["tau_tilde"], dtype=float)[order],
        }
    return panel_c_out, panel_d_out


def load_relaxation_table(data_root: Path, alt_root: Path | None = None) -> list[tuple[float, int, float, float]]:
    """Load relaxation_vs_reference_time_data.txt rows."""
    fpath = data_root / "relaxation_vs_reference_time_data.txt"
    if not fpath.is_file() and alt_root is not None:
        candidate = alt_root / "relaxation_vs_reference_time_data.txt"
        if candidate.is_file():
            fpath = candidate
    if not fpath.is_file():
        raise FileNotFoundError(
            f"relaxation_vs_reference_time_data.txt not found under {data_root}"
            + (f" or {alt_root}" if alt_root else "")
        )

    rows: list[tuple[float, int, float, float]] = []
    with open(fpath, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "coupling_name" in line or "ref_num" in line:
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                if parts[0].startswith("cavity_coupling_"):
                    coupling_value = float(parts[1])
                    ref_num = int(parts[3])
                    ref_time_ps = float(parts[4])
                    relaxation_time = float(parts[5])
                else:
                    ref_num = int(parts[0])
                    ref_time_ps = float(parts[1])
                    coupling_value = float(parts[4])
                    relaxation_time = float(parts[6])
            except ValueError:
                continue
            rows.append((coupling_value, ref_num, ref_time_ps, relaxation_time))
    return rows


def build_tn_tau_tilde_tables(
    analyzer: FictiveThreePanelAnalyzer,
) -> tuple[dict[float, dict[str, np.ndarray]], dict[float, dict[str, np.ndarray]]]:
    """
    Build TN-predicted tau_tilde tables for panels (c) and (d).

    Returns
    -------
    panel_c : dict[t_w_shifted] -> {lambda, tau_tilde}
    panel_d : dict[lambda] -> {t_w_shifted, tau_tilde}
    """
    material_time_data = analyzer.load_fictive_material_time_data()
    vs_coupling = analyzer.calculate_fictive_relaxation_times_vs_coupling(
        material_time_data
    )
    vs_waiting = analyzer.calculate_fictive_relaxation_times_vs_waiting_time(
        material_time_data
    )

    panel_c: dict[float, dict[str, np.ndarray]] = {}
    for tw_lab, data in vs_coupling.items():
        if tw_lab <= 0:
            continue
        couplings = np.asarray(
            [analyzer._coupling_from_storage(v) for v in data["coupling_strengths"]],
            dtype=float,
        )
        taus = np.asarray(data["normalization_times"], dtype=float)
        zero_mask = np.abs(couplings) < 1e-12
        if not np.any(zero_mask):
            continue
        tau0 = float(taus[zero_mask][0])
        if tau0 <= 0 or not np.isfinite(tau0):
            continue
        order = np.argsort(couplings)
        panel_c[float(tw_lab - 200.0)] = {
            "lambda": couplings[order],
            "tau_tilde": taus[order] / tau0,
        }

    panel_d: dict[float, dict[str, np.ndarray]] = {}
    ref = vs_waiting.get(0.0)
    if ref is None:
        return panel_c, panel_d
    ref_tw = np.asarray(ref["waiting_times"], dtype=float)
    ref_tau = np.asarray(ref["relaxation_times"], dtype=float)
    ref_map = {float(tw): float(tau) for tw, tau in zip(ref_tw, ref_tau) if tau > 0}

    for lam, data in vs_waiting.items():
        tws = np.asarray(data["waiting_times"], dtype=float)
        taus = np.asarray(data["relaxation_times"], dtype=float)
        xs: list[float] = []
        ys: list[float] = []
        for tw, tau in zip(tws, taus):
            if tw <= 0 or not np.isfinite(tau) or tau <= 0:
                continue
            tau_ref = ref_map.get(float(tw))
            if tau_ref is None:
                continue
            xs.append(float(tw - 200.0))
            ys.append(float(tau / tau_ref))
        if xs:
            order = np.argsort(xs)
            panel_d[float(lam)] = {
                "t_w": np.asarray(xs, dtype=float)[order],
                "tau_tilde": np.asarray(ys, dtype=float)[order],
            }
    return panel_c, panel_d


def plot_figure4(
    *,
    analyzer: FictiveThreePanelAnalyzer,
    output: Path,
    beta: float = DEFAULT_BETA,
    collapse_csv: Path | None = None,
    use_latex: bool = True,
    relaxation_rows: list[tuple[float, int, float, float]] | None = None,
) -> plt.Figure:
    """Render the four-panel Figure 4."""
    setup_figure4_fonts()
    # Always keep mathtext labels (Computer Modern math), even without usetex.
    use_latex = True

    measured = analyzer.load_direct_measured_material_time()
    tn = analyzer.load_tn_material_time()
    if relaxation_rows is None:
        figures_dir = Path(analyzer.base_dir).parent / "figures"
        relaxation_rows = load_relaxation_table(
            Path(analyzer.base_dir),
            alt_root=figures_dir,
        )
    panel_c, panel_d = build_measured_tau_tilde_tables(relaxation_rows)

    collapse_data: list[dict[str, np.ndarray]] = []
    if collapse_csv is not None and collapse_csv.is_file():
        collapse_data = load_panel_b_csv(collapse_csv)
    if not collapse_data:
        collapse_data = compute_collapse_from_staged(
            Path(analyzer.base_dir),
            analyzer._tag_map,
            switch_time_ps=analyzer.profile.switch_time_ps if analyzer.profile else 200.0,
        )

    lambdas = sorted(analyzer.coupling_values)
    lam_norm = Normalize(vmin=0.0, vmax=max(lambdas) if max(lambdas) > 0 else 1.0)
    cmap_lam = plt.colormaps.get_cmap("coolwarm")

    fig = plt.figure(figsize=(14.0, 3.6))
    gs = gridspec.GridSpec(1, 4, figure=fig, wspace=0.32)

    # Panel 1: material time
    ax_a = fig.add_subplot(gs[0, 0])
    for lam in lambdas:
        color = cmap_lam(lam_norm(analyzer.convert_epsilon_to_lambda(lam)))
        if lam in measured:
            t_m, xi_m = measured[lam]
            ax_a.plot(t_m, xi_m, color=color, lw=2.2, ls="-")
        if lam in tn:
            t_tn, xi_tn, _ = tn[lam]
            ax_a.plot(t_tn, xi_tn, color=color, lw=1.8, ls="--", alpha=0.9)
    ax_a.set_xlabel(latex_safe(r"$t$ (ps)", "t (ps)", use_latex), fontsize=14)
    ax_a.set_ylabel(
        latex_safe(
            r"$h_\lambda(t),\ h_{\lambda,\mathrm{TN}}(t)$",
            "h_lambda(t), h_lambda,TN(t)",
            use_latex,
        ),
        fontsize=14,
    )
    ax_a.set_xlim(left=0.0)
    ax_a.set_ylim(bottom=0.0)
    ax_a.grid(True, alpha=0.3, linestyle="--")
    ax_a.tick_params(labelsize=11)
    ax_a.legend(
        handles=[
            plt.Line2D([0], [0], color="black", lw=2.2, ls="-", label=latex_safe("Measured", "Measured", use_latex)),
            plt.Line2D([0], [0], color="black", lw=1.8, ls="--", label=latex_safe("TN model", "TN model", use_latex)),
        ],
        loc="upper left",
        fontsize=9,
        frameon=True,
    )
    divider_a = make_axes_locatable(ax_a)
    cax_a = divider_a.append_axes("top", size="6%", pad=0.12)
    sm_a = cm.ScalarMappable(norm=lam_norm, cmap=cmap_lam)
    sm_a.set_array([])
    cbar_a = fig.colorbar(sm_a, cax=cax_a, orientation="horizontal")
    cbar_a.set_label(latex_safe(r"$\lambda$ (a.u.)", "lambda (a.u.)", use_latex), fontsize=11)
    cax_a.xaxis.set_ticks_position("top")
    cax_a.xaxis.set_label_position("top")
    cbar_a.ax.tick_params(labelsize=9)

    # Panel 2: collapse
    ax_b = fig.add_subplot(gs[0, 1])
    for entry in collapse_data:
        color = cmap_lam(lam_norm(float(entry["lambda_au"])))
        mask = (entry["h"] >= 0.0) & (entry["h"] <= 2.0) & np.isfinite(entry["phi"])
        if np.any(mask):
            ax_b.plot(entry["h"][mask], entry["phi"][mask], color=color, alpha=0.18, lw=1.0)
    h_fit = np.linspace(0.0, 2.0, 300)
    phi_fit = stretched_exponential_guide(h_fit, beta)
    ax_b.plot(
        h_fit,
        phi_fit,
        color="black",
        lw=2.5,
        zorder=10,
        label=latex_safe(
            rf"$\Phi_k(h)=e^{{-h^{{\beta}}}}$, $\beta={beta:.2f}$",
            f"Phi_k(h)=exp(-h^beta), beta={beta:.2f}",
            use_latex,
        ),
    )
    ax_b.set_xlabel(latex_safe(r"$h$", "h", use_latex), fontsize=14)
    ax_b.set_ylabel(latex_safe(r"$\Phi_k(h)$", "Phi_k(h)", use_latex), fontsize=14)
    ax_b.set_xlim(0.0, 2.0)
    ax_b.set_ylim(0.0, 1.05)
    ax_b.grid(True, alpha=0.3, linestyle="--")
    ax_b.tick_params(labelsize=11)
    ax_b.legend(loc="upper right", fontsize=9, frameon=True)

    # Panel 3: measured tau_tilde vs lambda
    ax_c = fig.add_subplot(gs[0, 2])
    tw_keys = sorted(panel_c.keys())
    tw_norm = Normalize(vmin=min(tw_keys) if tw_keys else 0.0, vmax=max(tw_keys) if tw_keys else 1.0)
    cmap_tw = plt.colormaps.get_cmap("viridis")
    for tw in tw_keys:
        data = panel_c[tw]
        color = cmap_tw(tw_norm(tw))
        ax_c.plot(
            data["lambda"],
            data["tau_tilde"],
            "-D",
            color=color,
            lw=2.2,
            ms=6,
            markerfacecolor="white",
            markeredgecolor=color,
            markeredgewidth=1.6,
        )
    ax_c.axhline(1.0, color="gray", ls="--", lw=1.2, alpha=0.7)
    ax_c.set_xlabel(latex_safe(r"$\lambda$ (a.u.)", "lambda (a.u.)", use_latex), fontsize=14)
    ax_c.set_ylabel(latex_safe(r"$\tilde{\tau}_{s}$", "tau_tilde_s", use_latex), fontsize=14)
    ax_c.set_xlim(left=0.0)
    ax_c.set_ylim(bottom=0.7)
    ax_c.grid(True, alpha=0.3, linestyle="--")
    ax_c.tick_params(labelsize=11)
    divider_c = make_axes_locatable(ax_c)
    cax_c = divider_c.append_axes("top", size="6%", pad=0.12)
    sm_c = cm.ScalarMappable(norm=tw_norm, cmap=cmap_tw)
    sm_c.set_array([])
    cbar_c = fig.colorbar(sm_c, cax=cax_c, orientation="horizontal")
    cbar_c.set_label(latex_safe(r"$t_{w}$ (ps)", "t_w (ps)", use_latex), fontsize=11)
    cax_c.xaxis.set_ticks_position("top")
    cax_c.xaxis.set_label_position("top")
    cbar_c.ax.tick_params(labelsize=9)

    # Panel 4: measured tau_tilde vs waiting time
    ax_d = fig.add_subplot(gs[0, 3])
    for lam in sorted(panel_d.keys()):
        data = panel_d[lam]
        color = cmap_lam(lam_norm(analyzer.convert_epsilon_to_lambda(lam)))
        ax_d.plot(
            data["t_w"],
            data["tau_tilde"],
            "-D",
            color=color,
            lw=2.2,
            ms=6,
            markerfacecolor="white",
            markeredgecolor=color,
            markeredgewidth=1.6,
            label=latex_safe(
                rf"$\lambda={analyzer.convert_epsilon_to_lambda(lam):g}$",
                f"lambda={analyzer.convert_epsilon_to_lambda(lam):g}",
                use_latex,
            ),
        )
    ax_d.axhline(1.0, color="gray", ls="--", lw=1.2, alpha=0.7)
    ax_d.set_xlabel(latex_safe(r"$t_{w}$ (ps)", "t_w (ps)", use_latex), fontsize=14)
    ax_d.set_ylabel(latex_safe(r"$\tilde{\tau}_{s}$", "tau_tilde_s", use_latex), fontsize=14)
    ax_d.set_xlim(left=0.0, right=max(1800.0, max((max(data["t_w"]) for data in panel_d.values()), default=1600.0) + 50.0))
    ax_d.set_ylim(bottom=0.7)
    ax_d.grid(True, alpha=0.3, linestyle="--")
    ax_d.tick_params(labelsize=11)
    ax_d.legend(loc="upper right", fontsize=8, frameon=True)

    fig.subplots_adjust(top=0.86, bottom=0.16, left=0.06, right=0.98, wspace=0.35)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    return fig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", default="aging_weak_lambda")
    parser.add_argument("--staged-root", type=Path, default=None)
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output PDF path (default: profile figures_dir/figure4.pdf)",
    )
    parser.add_argument(
        "--collapse-csv",
        type=Path,
        default=None,
        help="Optional figure4_panel_b.csv; otherwise collapse is recomputed",
    )
    parser.add_argument("--beta", type=float, default=DEFAULT_BETA)
    parser.add_argument("--no-latex", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    from repro_bootstrap import setup_profile

    profile = setup_profile(args, default="aging_weak_lambda")
    staged_root = args.staged_root or profile.staged_root
    figures_dir = profile.figures_dir or staged_root
    output = args.output or (figures_dir / "figure4.pdf")
    collapse_csv = args.collapse_csv
    if collapse_csv is None:
        candidate = figures_dir / "figure4_panel_b.csv"
        if candidate.is_file():
            collapse_csv = candidate

    analyzer = FictiveThreePanelAnalyzer(
        base_dir=str(staged_root),
        use_latex=not args.no_latex,
        profile=profile if profile.name != "paper" else None,
    )
    plot_figure4(
        analyzer=analyzer,
        output=output,
        beta=args.beta,
        collapse_csv=collapse_csv,
        use_latex=not args.no_latex,
    )
    print(f"Wrote {output}")
    print(f"Wrote {output.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
