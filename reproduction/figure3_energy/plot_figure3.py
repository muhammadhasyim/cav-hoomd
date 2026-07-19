#!/usr/bin/env python3
"""
Plot Figure 3 panels (b) and (c): energy deviations and fictive/kinetic temperatures.

Publication styling:
  - Computer Modern LaTeX (via latex_config_adobe)
  - Panel (b): ΔV harmonic / LJ+Coulomb / total vs t
  - Panel (c): T_v (harmonic), T_s (LJ+Coulomb), T_k (kinetic bath) vs t

Data pipeline matches export_figure3_csv.py.
"""

from __future__ import annotations

import os
os.environ.pop("LD_LIBRARY_PATH", None)
os.environ.pop("LD_PRELOAD", None)

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

SCRIPT_DIR = Path(__file__).resolve().parent
REPRO_ROOT = SCRIPT_DIR.parents[0]
sys.path.insert(0, str(REPRO_ROOT / "shared"))
from latex_config_adobe import latex_safe, setup_illustrator_safe_fonts

ARCHIVE = Path("/media/extradrive/Trajectories/cavity_supercooled_archive")
DEFAULT_BASE = ARCHIVE / "final_production_run/data/derived/2026-01-29/time_series_output"

COUPLING_TAG = {
    0.0: "0epos00",
    3e-4: "3eneg04",
    5e-4: "5eneg04",
    7e-4: "7eneg04",
    1e-3: "1eneg03",
}


def coupling_tag_for(coupling: float, profile=None) -> str:
    if profile is not None:
        return profile.entry_for_axis_value(coupling).epsilon_tag
    if coupling not in COUPLING_TAG:
        raise ValueError(f"Unsupported coupling {coupling}; choose from {list(COUPLING_TAG)}")
    return COUPLING_TAG[coupling]


def parse_coupling(value: str) -> float:
    if "eneg" in value:
        parts = value.replace("coupling_", "").split("eneg")
        return int(parts[0]) * (10 ** (-int(parts[1])))
    return float(value)


def load_numeric(path: Path) -> np.ndarray:
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                float(line.split()[0])
            except ValueError:
                continue
            rows.append([float(x) for x in line.split()])
    return np.array(rows, dtype=np.float64)


def load_named_columns(path: Path) -> dict[str, np.ndarray]:
    """Load a whitespace table keyed by its bare (non-comment) header row.

    Robust to column-count drift: the potential-energy files may be written with
    a 5-column layout (time, harmonic, lj, coul, total) or a 6-column layout that
    also carries an explicit combined ``lj_coul_ha`` column. Reading by header
    name avoids silently mapping ``lj_coul`` onto ``total`` (or vice versa).
    """
    header: list[str] | None = None
    rows: list[list[float]] = []
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            tokens = stripped.split()
            try:
                float(tokens[0])
            except ValueError:
                header = tokens
                continue
            rows.append([float(x) for x in tokens])
    if header is None:
        raise ValueError(f"No column header row found in {path}")
    arr = np.array(rows, dtype=np.float64)
    if arr.shape[1] != len(header):
        raise ValueError(
            f"Header/data column mismatch in {path}: "
            f"{len(header)} names vs {arr.shape[1]} data columns"
        )
    return {name: arr[:, i] for i, name in enumerate(header)}


def align_to_master(master_t: np.ndarray, arr: np.ndarray, t_col: int = 0) -> np.ndarray:
    t_src = np.round(arr[:, t_col], 6)
    t_tgt = np.round(master_t, 6)
    idx = np.searchsorted(t_src, t_tgt)
    idx = np.clip(idx, 0, len(t_src) - 1)
    return arr[idx]


def interp_nans(col: np.ndarray) -> np.ndarray:
    col = col.copy()
    nans = np.isnan(col)
    if nans.any():
        idx = np.arange(len(col))
        col[nans] = np.interp(idx[nans], idx[~nans], col[~nans])
    return col


def moving_average(x: np.ndarray, window: int) -> np.ndarray:
    """Centered moving average with edge padding, preserving array length."""
    window = int(window)
    if window <= 1:
        return x
    pad = window // 2
    x_padded = np.pad(x, (pad, pad), mode="edge")
    kernel = np.ones(window) / window
    smoothed = np.convolve(x_padded, kernel, mode="valid")
    return smoothed[: len(x)]


def apply_display_smoothing(
    data: dict,
    *,
    smooth_window: int = 1,
    smooth_window_total: int = 1,
    smooth_harmonic: bool = False,
) -> dict:
    """Return a copy with display-only moving averages applied to noisy channels."""
    display = dict(data)
    display["dV_harm"] = (
        moving_average(data["dV_harm"], smooth_window)
        if smooth_harmonic
        else data["dV_harm"]
    )
    display["dV_ljc"] = moving_average(data["dV_ljc"], smooth_window)
    display["dV_tot"] = moving_average(data["dV_tot"], smooth_window_total)
    if data.get("dV_lj") is not None:
        display["dV_lj"] = moving_average(data["dV_lj"], smooth_window)
    if data.get("dV_coul") is not None:
        display["dV_coul"] = moving_average(data["dV_coul"], smooth_window)
    return display


def _load_analysis_energy_uncertainty(
    path: Path,
    master_t: np.ndarray,
    *,
    t_switch: float,
) -> dict[str, np.ndarray] | None:
    """Load replica SEM columns from avg_energies CSV when present."""
    import csv

    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []
        if "coulombic_std_ha" not in fieldnames:
            return None
        rows = list(reader)

    if not rows:
        return None

    time_ps = np.asarray([float(row["time_ps"]) for row in rows], dtype=float)
    harmonic_std = np.asarray([float(row["harmonic_std_ha"]) for row in rows], dtype=float)
    lj_std = np.asarray([float(row["lj_std_ha"]) for row in rows], dtype=float)
    coul_std = np.asarray([float(row["coulombic_std_ha"]) for row in rows], dtype=float)
    lj_coul_std = np.asarray([float(row["lj_coul_std_ha"]) for row in rows], dtype=float)

    aligned = np.column_stack(
        [
            time_ps,
            harmonic_std,
            lj_std,
            coul_std,
            lj_coul_std,
        ]
    )
    aligned = align_to_master(master_t, aligned)
    total_std = np.sqrt(
        aligned[:, 1] ** 2 + aligned[:, 2] ** 2 + aligned[:, 3] ** 2
    )
    pre = master_t < t_switch
    if pre.any():
        # Baseline subtraction is constant; replica std of ΔV equals component std.
        pass
    return {
        "dV_harm_sem": aligned[:, 1],
        "dV_lj_sem": aligned[:, 2],
        "dV_coul_sem": aligned[:, 3],
        "dV_ljc_sem": aligned[:, 4],
        "dV_tot_sem": total_std,
    }


def load_figure3_data(
    base_dir: Path,
    coupling: float,
    t_switch: float = 200.0,
    profile=None,
    analysis_energy_csv: Path | None = None,
):
    tag = coupling_tag_for(coupling, profile=profile)
    pe_file = base_dir / f"coupling_{tag}_averaged_potential_energy.txt"
    fict_file = base_dir / f"coupling_{tag}_averaged_fictive_temperatures.txt"
    fts_file = base_dir / f"coupling_{tag}_fictive_time_series.txt"

    for p in (pe_file, fict_file, fts_file):
        if not p.exists():
            raise FileNotFoundError(f"Missing required file: {p}")

    pe_cols = load_named_columns(pe_file)
    fict = load_numeric(fict_file)
    fts = load_numeric(fts_file)

    t = pe_cols["time_ps"]
    fict_a = align_to_master(t, fict)
    fts_a = align_to_master(t, fts)

    pre = t < t_switch
    if pre.sum() == 0:
        raise ValueError(f"No points before t={t_switch} ps")

    U_harm = interp_nans(pe_cols["harmonic_ha"])
    U_tot = interp_nans(pe_cols["total_ha"])

    eq_harm = np.mean(U_harm[pre])
    eq_tot = np.mean(U_tot[pre])

    dV_harm = U_harm - eq_harm
    dV_tot = U_tot - eq_tot

    # Combined LJ+Coulomb deviation: prefer an explicit combined column, else
    # sum the separate LJ and Coulomb columns.
    if "lj_coul_ha" in pe_cols:
        U_ljc = interp_nans(pe_cols["lj_coul_ha"])
        dV_ljc = U_ljc - np.mean(U_ljc[pre])
    else:
        U_lj = interp_nans(pe_cols["lj_ha"])
        U_coul = interp_nans(pe_cols.get("coul_ha", np.zeros_like(U_lj)))
        dV_ljc = (U_lj - np.mean(U_lj[pre])) + (U_coul - np.mean(U_coul[pre]))

    dV_lj = None
    dV_coul = None
    if "lj_ha" in pe_cols and "coul_ha" in pe_cols:
        U_lj_only = interp_nans(pe_cols["lj_ha"])
        U_coul_only = interp_nans(pe_cols["coul_ha"])
        # When Coulomb is omitted (LJ-only structural proxy), keep a single
        # LJ+Coulomb curve rather than plotting a flat zero Coulomb channel.
        if float(np.nanmax(np.abs(pe_cols["coul_ha"]))) > 1e-12:
            dV_lj = U_lj_only - np.mean(U_lj_only[pre])
            dV_coul = U_coul_only - np.mean(U_coul_only[pre])
            dV_ljc = dV_lj + dV_coul

    T_v = fict_a[:, 2]
    T_s = fts_a[:, 1]

    kt_file = base_dir / f"coupling_{tag}_kinetic_temperature.txt"
    if kt_file.is_file():
        kt_cols = load_named_columns(kt_file)
        T_k = kt_cols["kinetic_temp_K"]
        if len(T_k) != len(t):
            kt_aligned = align_to_master(t, np.column_stack([kt_cols["time_ps"], T_k]))
            T_k = kt_aligned[:, 1]
    else:
        T_k = np.full(len(t), 100.0)

    result = {
        "t": t,
        "dV_harm": dV_harm,
        "dV_ljc": dV_ljc,
        "dV_lj": dV_lj,
        "dV_coul": dV_coul,
        "dV_tot": dV_tot,
        "T_v": T_v,
        "T_s": T_s,
        "T_k": T_k,
    }
    if analysis_energy_csv is not None and analysis_energy_csv.is_file():
        unc = _load_analysis_energy_uncertainty(
            analysis_energy_csv,
            t,
            t_switch=t_switch,
        )
        if unc is not None:
            result.update(unc)
    return result


# Reference paper panels are individual near-square figures (~1024x844 px,
# aspect ratio 1024/844 = 1.213). We now render (b) and (c) as two axes inside
# ONE figure with a shared GridSpec, which is the only way to *guarantee*
# identical axes box dimensions, fonts, and legend styling between them --
# saving them as two separate bbox_inches="tight" files let each panel's
# content (e.g. the inset in (c)) crop differently and drift out of sync.
PANEL_FIGSIZE = (4.1, 3.38)  # inches; 4.1 / 3.38 = 1.213 (standalone single-panel use)
COMBINED_FIGSIZE = (9.6, 3.55)
XLABEL_FONTSIZE = 16
YLABEL_FONTSIZE = 18
LABEL_FONTSIZE = YLABEL_FONTSIZE  # legacy alias
TICK_FONTSIZE = 12
YTICK_FONTSIZE = 13.5
LEGEND_FONTSIZE = 11
PANEL_LETTER_FONTSIZE = 22
LINEWIDTH = 1.8


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


def setup_figure3_fonts(*, prefer_illustrator: bool = True) -> bool:
    """Configure Computer Modern fonts for Figure 3.

    Tries Illustrator-safe CMU OpenType first, then matplotlib's bundled
    ``cmr10`` + mathtext CM (needed when ``latex.fmt`` / CMU fonts are missing).
    Always keeps mathtext ``$...$`` labels.
    """
    plt.style.use("classic")
    if prefer_illustrator and setup_illustrator_safe_fonts():
        return True
    try:
        cm_family = _register_matplotlib_cm_fonts()
    except OSError:
        cm_family = "DejaVu Serif"
    plt.rcParams.update(
        {
            "text.usetex": False,
            "font.family": cm_family,
            "font.serif": [cm_family, "Computer Modern Roman", "DejaVu Serif"],
            "mathtext.fontset": "cm",
            "axes.formatter.use_mathtext": True,
            "font.size": 14,
            "axes.labelsize": YLABEL_FONTSIZE,
            "xtick.labelsize": TICK_FONTSIZE,
            "ytick.labelsize": YTICK_FONTSIZE,
            "legend.fontsize": LEGEND_FONTSIZE,
            "axes.linewidth": 1.2,
            "axes.unicode_minus": False,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "legend.framealpha": 0.95,
            "axes.axisbelow": True,
        }
    )
    return True


def _style_figure3_axes(ax) -> None:
    """Inward ticks and grid matching the paper screenshot."""
    ax.grid(True, alpha=0.3, linestyle="-")
    ax.tick_params(
        axis="x", which="major", direction="in", top=True, labelsize=TICK_FONTSIZE,
    )
    ax.tick_params(
        axis="y", which="major", direction="in", right=True, labelsize=YTICK_FONTSIZE,
    )
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)


def _with_ext(output_path: Path, ext: str) -> Path:
    """Append an extension without clobbering dotted stems (e.g. '0.01')."""
    if output_path.suffix.lower() in {".pdf", ".png", ".svg"}:
        return output_path.with_suffix(ext)
    return output_path.parent / (output_path.name + ext)


def _save_fig(fig, output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pdf_path = _with_ext(output_path, ".pdf")
    png_path = _with_ext(output_path, ".png")
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {pdf_path}")
    print(f"Saved {png_path}")


def _panel_label(ax, label: str):
    # Plain (non-bold) serif letter outside the axes box, top-left, matching
    # the reference figures. Positioned in axes-fraction coords so it lines
    # up identically regardless of where the axes sits in the figure.
    ax.text(-0.30, 1.08, label, transform=ax.transAxes, fontsize=PANEL_LETTER_FONTSIZE,
            ha="left", va="top")


STACKED_FIGSIZE = (10.0, 14.0)
STACKED_LABEL_FONTSIZE = 14
STACKED_TICK_FONTSIZE = 11


def _autoscale_ylim(ax, values: np.ndarray, ref_line: float | None = None, pad_frac: float = 0.08) -> None:
    """Set y-limits from finite data with modest padding around an optional reference."""
    vals = values[np.isfinite(values)]
    if vals.size == 0:
        return
    v_min = float(np.min(vals))
    v_max = float(np.max(vals))
    if ref_line is not None:
        v_min = min(v_min, ref_line)
        v_max = max(v_max, ref_line)
    span = max(v_max - v_min, 1e-12)
    ax.set_ylim(v_min - pad_frac * span, v_max + pad_frac * span)


def draw_stacked_figure3(
    data: dict,
    output_path: Path,
    t_switch: float,
    t_max: float,
    use_latex: bool,
    *,
    coupling_label: str | None = None,
    smooth_window: int = 1,
    smooth_window_total: int = 1,
) -> None:
    """Render five vertically stacked panels mirroring the aging analysis plots.

    Energy panels show deviations from the pre-switch equilibrium (ΔV in Hartree).
    Temperature panels show baseline-shifted fictive temperatures (100 K reference).
    """
    t = data["t"]
    mask = t <= t_max

    display = apply_display_smoothing(
        data,
        smooth_window=smooth_window,
        smooth_window_total=smooth_window_total,
    )
    dV_harm = display["dV_harm"]
    dV_ljc = display["dV_ljc"]
    dV_tot = display["dV_tot"]
    dV_lj = display.get("dV_lj")
    dV_coul = display.get("dV_coul")
    t_v = data["T_v"]
    t_s = data["T_s"]

    panels = [
        {
            "ylabel": latex_safe(r"$\Delta V_{\mathrm{harm}}$ (Ha)", "ΔV harmonic (Ha)", use_latex),
            "values": dV_harm[mask],
            "color": "blue",
            "ref": 0.0,
            "ref_style": {"color": "0.5", "ls": ":", "lw": 1, "alpha": 0.8},
        },
        {
            "ylabel": latex_safe(r"$\Delta V_{\mathrm{LJ+Coul}}$ (Ha)", "ΔV LJ+Coul (Ha)", use_latex),
            "values": dV_ljc[mask],
            "color": "red",
            "ref": 0.0,
            "ref_style": {"color": "0.5", "ls": ":", "lw": 1, "alpha": 0.8},
            "overlay": (
                [
                    (dV_lj[mask], "darkred", latex_safe(r"$\Delta V_{\mathrm{LJ}}$", "ΔV LJ", use_latex)),
                    (dV_coul[mask], "orange", latex_safe(r"$\Delta V_{\mathrm{Coul}}$", "ΔV Coul", use_latex)),
                ]
                if dV_lj is not None and dV_coul is not None
                else None
            ),
        },
        {
            "ylabel": latex_safe(r"$\Delta V_{\mathrm{tot}}$ (Ha)", "ΔV total (Ha)", use_latex),
            "values": dV_tot[mask],
            "color": "green",
            "ref": 0.0,
            "ref_style": {"color": "0.5", "ls": ":", "lw": 1, "alpha": 0.8},
        },
        {
            "ylabel": latex_safe(r"$T_f^{\mathrm{harm}}$ (K)", "Harmonic Tf (K)", use_latex),
            "values": t_v[mask],
            "color": "blue",
            "ref": 100.0,
            "ref_style": {"color": "0.5", "ls": ":", "lw": 1, "alpha": 0.8},
        },
        {
            "ylabel": latex_safe(r"$T_f^{\mathrm{LJ+Coul}}$ (K)", "LJ + Coul Tf (K)", use_latex),
            "values": t_s[mask],
            "color": "red",
            "ref": 100.0,
            "ref_style": {"color": "0.5", "ls": ":", "lw": 1, "alpha": 0.8},
        },
    ]

    fig, axes = plt.subplots(
        len(panels), 1, figsize=STACKED_FIGSIZE, sharex=True,
        layout="constrained",
    )
    if len(panels) == 1:
        axes = [axes]

    t_plot = t[mask]
    for ax, panel in zip(axes, panels):
        ax.plot(t_plot, panel["values"], color=panel["color"], lw=1.8)
        overlay = panel.get("overlay")
        if overlay:
            for values, color, label in overlay:
                ax.plot(t_plot, values, color=color, lw=1.4, ls="--", alpha=0.9, label=label)
            ax.legend(loc="best", fontsize=8, frameon=True)
        ax.axhline(panel["ref"], **panel["ref_style"])
        ax.axvline(t_switch, color="k", ls="--", lw=1, alpha=0.5)
        ax.set_ylabel(panel["ylabel"], fontsize=STACKED_LABEL_FONTSIZE)
        ax.set_xlim(0, t_max)
        _autoscale_ylim(ax, panel["values"], ref_line=panel["ref"])
        ax.tick_params(labelsize=STACKED_TICK_FONTSIZE)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel(latex_safe(r"$t$ (ps)", "Time (ps)", use_latex), fontsize=STACKED_LABEL_FONTSIZE)

    title = "Energy deviations and fictive temperatures"
    if coupling_label is not None:
        title += f" ({coupling_label})"
    axes[0].set_title(title, fontsize=STACKED_LABEL_FONTSIZE + 1)

    _save_fig(fig, output_path)


def draw_energy_panel(ax, data: dict, t_switch: float, t_max: float, use_latex: bool,
                       smooth_window: int = 1, smooth_window_total: int = 1,
                       auto_ylim: bool = False):
    t = data["t"]
    mask = t <= t_max

    display = apply_display_smoothing(
        data,
        smooth_window=smooth_window,
        smooth_window_total=smooth_window_total,
    )
    dV_harm = display["dV_harm"]
    dV_ljc = display["dV_ljc"]
    dV_tot = display["dV_tot"]
    dV_lj = display.get("dV_lj")
    dV_coul = display.get("dV_coul")

    ax.plot(t[mask], dV_harm[mask], color="blue", lw=LINEWIDTH, label="Harmonic")
    ax.plot(t[mask], dV_ljc[mask], color="red", lw=LINEWIDTH, label="LJ+Coulomb")
    if data.get("dV_ljc_sem") is not None:
        sem = data["dV_ljc_sem"]
        ax.fill_between(
            t[mask],
            dV_ljc[mask] - sem[mask],
            dV_ljc[mask] + sem[mask],
            color="red",
            alpha=0.15,
            linewidth=0,
            label="LJ+Coulomb ± SEM",
        )
    if dV_lj is not None and dV_coul is not None:
        ax.plot(t[mask], dV_lj[mask], color="darkred", lw=LINEWIDTH, ls="--", label="LJ")
        ax.plot(t[mask], dV_coul[mask], color="orange", lw=LINEWIDTH, ls="--", label="Coulomb")
        if data.get("dV_coul_sem") is not None:
            sem_coul = data["dV_coul_sem"]
            ax.fill_between(
                t[mask],
                dV_coul[mask] - sem_coul[mask],
                dV_coul[mask] + sem_coul[mask],
                color="orange",
                alpha=0.12,
                linewidth=0,
            )
    ax.plot(t[mask], dV_tot[mask], color="green", lw=LINEWIDTH, label="Total")
    if data.get("dV_tot_sem") is not None:
        sem_tot = data["dV_tot_sem"]
        ax.fill_between(
            t[mask],
            dV_tot[mask] - sem_tot[mask],
            dV_tot[mask] + sem_tot[mask],
            color="green",
            alpha=0.12,
            linewidth=0,
        )
    ax.axhline(0, color="0.4", ls="--", lw=1, label="Equilibrium")
    ax.axvline(t_switch, color="red", ls=":", lw=1.5, alpha=0.8)
    ax.set_xlabel(latex_safe(r"$t$ (ps)", "t (ps)", use_latex), fontsize=XLABEL_FONTSIZE)
    # Paper uses arbitrary units for ΔV; values are still Hartree internally.
    ax.set_ylabel(latex_safe(r"$\Delta V$ (a.u.)", "ΔV (a.u.)", use_latex), fontsize=YLABEL_FONTSIZE)
    ax.set_xlim(0, t_max)

    # Legend sits upper-right (empty region at large t), so the peak -- which
    # occurs at low t, upper-left -- only needs modest headroom above it,
    # letting the y-axis max come down close to the actual peak height.
    # For aging profiles, scale from harmonic + LJ only so the total spike
    # does not compress the component curves (paper 2-panel layout).
    if auto_ylim:
        v_all = np.concatenate([dV_harm[mask], dV_ljc[mask]])
    else:
        v_all = np.concatenate([dV_harm[mask], dV_ljc[mask], dV_tot[mask]])
    v_max, v_min = np.nanmax(v_all), np.nanmin(v_all)
    span = v_max - v_min
    ax.set_ylim(v_min - 0.10 * span, v_max + 0.18 * span)

    _style_figure3_axes(ax)
    ax.legend(fontsize=LEGEND_FONTSIZE, loc="upper right", framealpha=0.95,
              handlelength=1.3, labelspacing=0.3, borderpad=0.35)
    if smooth_window > 1 or smooth_window_total > 1:
        ax.text(
            0.02,
            0.02,
            f"display smooth: {smooth_window} ps"
            + (f" (total {smooth_window_total} ps)" if smooth_window_total != smooth_window else ""),
            transform=ax.transAxes,
            fontsize=8,
            color="0.35",
            ha="left",
            va="bottom",
        )


def draw_temperature_panel(ax, data: dict, t_switch: float, t_max: float, use_latex: bool,
                            auto_ylim: bool = False):
    t = data["t"]
    mask = t <= t_max

    ax.plot(t[mask], data["T_v"][mask], color="blue", lw=LINEWIDTH,
            label=latex_safe(r"$T_v$ Harmonic", "T_v Harmonic", use_latex))
    ax.plot(t[mask], data["T_s"][mask], color="red", lw=LINEWIDTH,
            label=latex_safe(r"$T_s$ LJ+Coulomb", "T_s LJ+Coulomb", use_latex))
    ax.plot(t[mask], data["T_k"][mask], color="orange", lw=LINEWIDTH,
            label=latex_safe(r"$T_k$ Kinetic", "T_k Kinetic", use_latex))
    ax.axhline(100, color="0.4", ls="--", lw=1,
               label=latex_safe(r"$100$ K", "100 K", use_latex))
    ax.axvline(t_switch, color="red", ls=":", lw=1.5, alpha=0.8)
    ax.set_xlabel(latex_safe(r"$t$ (ps)", "t (ps)", use_latex), fontsize=XLABEL_FONTSIZE)
    ax.set_ylabel(latex_safe(r"$T$ (K)", "T (K)", use_latex), fontsize=YLABEL_FONTSIZE)
    ax.set_xlim(0, t_max)
    if auto_ylim:
        vals = np.concatenate([data["T_v"][mask], data["T_s"][mask], data["T_k"][mask]])
        vals = vals[np.isfinite(vals)]
        v_max = max(float(np.max(vals)), 100.0)
        span = max(v_max - 100.0, 10.0)
        # Paper framing: baseline at bottom of the range (ymin = 0).
        ax.set_ylim(0.0, v_max + 0.15 * span)
    else:
        ax.set_ylim(0, 800.0)
    _style_figure3_axes(ax)
    ax.legend(fontsize=LEGEND_FONTSIZE, loc="upper left", framealpha=0.95,
              handlelength=1.3, labelspacing=0.3, borderpad=0.35)

    # Inset zoom on structural/kinetic recovery (paper: 250-750 ps, 70-100 K)
    inset = inset_axes(ax, width="36%", height="36%", loc="upper right",
                        bbox_to_anchor=(0.02, -0.04, 0.98, 1.0), bbox_transform=ax.transAxes)
    inset.plot(t[mask], data["T_s"][mask], color="red", lw=1.5)
    inset.plot(t[mask], data["T_k"][mask], color="orange", lw=1.5)
    inset.axhline(100, color="0.4", ls="--", lw=1)
    inset.axvline(t_switch, color="red", ls=":", lw=1, alpha=0.8)
    inset.set_xlim(250, 750)
    inset.set_ylim(70, 100)
    inset.set_xticks([250, 500, 750])
    inset.set_yticks([70, 80, 90, 100])
    inset.grid(True, alpha=0.3)
    inset.tick_params(labelsize=9, direction="in")
    for spine in inset.spines.values():
        spine.set_linewidth(1.0)
    mark_inset(ax, inset, loc1=2, loc2=4, fc="none", ec="0.5", ls="--", alpha=0.7)


def plot_panel_b(data: dict, output_path: Path, t_switch: float, t_max: float, use_latex: bool,
                  smooth_window: int = 1, smooth_window_total: int = 1, auto_ylim: bool = False):
    fig, ax = plt.subplots(figsize=PANEL_FIGSIZE)
    draw_energy_panel(ax, data, t_switch, t_max, use_latex, smooth_window, smooth_window_total,
                      auto_ylim=auto_ylim)
    _panel_label(ax, "(b)")
    fig.tight_layout(rect=(0.08, 0, 1, 1))
    _save_fig(fig, output_path)


def plot_panel_c(data: dict, output_path: Path, t_switch: float, t_max: float, use_latex: bool,
                  auto_ylim: bool = False):
    fig, ax = plt.subplots(figsize=PANEL_FIGSIZE)
    draw_temperature_panel(ax, data, t_switch, t_max, use_latex, auto_ylim=auto_ylim)
    _panel_label(ax, "(c)")
    fig.tight_layout(rect=(0.08, 0, 1, 1))
    _save_fig(fig, output_path)


def plot_combined_figure3(data: dict, output_path: Path, t_switch: float,
                           t_max_b: float, t_max_c: float, use_latex: bool,
                           smooth_window: int = 1, smooth_window_total: int = 1,
                           auto_ylim: bool = False, auto_ylim_energy: bool = False):
    fig, (ax_b, ax_c) = plt.subplots(
        1, 2, figsize=COMBINED_FIGSIZE, gridspec_kw={"wspace": 0.28},
    )

    draw_energy_panel(ax_b, data, t_switch, t_max_b, use_latex, smooth_window, smooth_window_total,
                      auto_ylim=auto_ylim_energy)
    draw_temperature_panel(ax_c, data, t_switch, t_max_c, use_latex, auto_ylim=auto_ylim)

    _panel_label(ax_b, "(b)")
    _panel_label(ax_c, "(c)")

    # Avoid tight_layout: inset_axes is incompatible and can shift panel boxes.
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.16, top=0.90, wspace=0.28)
    _save_fig(fig, output_path)


def plot_figure3(data: dict, output_path: Path, t_switch: float = 200.0,
                 t_max: float = 2000.0, t_max_b: float | None = None, use_latex: bool = True,
                 smooth_window: int = 1, smooth_window_total: int = 1, combined: bool = True,
                 auto_ylim: bool = False, auto_ylim_energy: bool = False,
                 layout: str = "paper", coupling_label: str | None = None):
    # Prefer Illustrator-safe CMU OpenType; otherwise matplotlib cmr10 + mathtext CM.
    # Always keep math-mode label strings (``latex_safe`` with use_latex=True).
    prefer_illustrator = bool(use_latex)
    setup_figure3_fonts(prefer_illustrator=prefer_illustrator)
    use_latex = True

    if layout == "stacked":
        draw_stacked_figure3(
            data, output_path, t_switch, t_max, use_latex,
            coupling_label=coupling_label,
            smooth_window=smooth_window,
            smooth_window_total=smooth_window_total,
        )
        return

    t_max_b = t_max if t_max_b is None else t_max_b

    if combined:
        plot_combined_figure3(data, output_path, t_switch, t_max_b, t_max, use_latex,
                               smooth_window, smooth_window_total,
                               auto_ylim=auto_ylim, auto_ylim_energy=auto_ylim_energy)
    else:
        stem = output_path
        plot_panel_b(data, stem.with_name(stem.name + "_panel_b"), t_switch, t_max_b, use_latex,
                     smooth_window, smooth_window_total, auto_ylim=auto_ylim_energy)
        plot_panel_c(data, stem.with_name(stem.name + "_panel_c"), t_switch, t_max, use_latex,
                     auto_ylim=auto_ylim)


def main():
    parser = argparse.ArgumentParser(description="Plot Figure 3 energy and temperature panels")
    parser.add_argument("--profile", default="paper", help="Dataset profile (default: paper)")
    parser.add_argument("--staged-root", type=Path, default=None, help="Override staged data root")
    parser.add_argument("--base-dir", type=Path, default=None, help="Time-series directory override")
    parser.add_argument("--coupling", type=str, default="1e-3",
                        help="Coupling strength (paper ε or profile axis value, e.g. 0.01)")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--t-max", type=float, default=2000.0,
                        help="x-axis upper limit (ps) for panel (c)")
    parser.add_argument("--t-max-b", type=float, default=2000.0,
                        help="x-axis upper limit (ps) for panel (b)")
    parser.add_argument("--t-switch", type=float, default=200.0)
    parser.add_argument("--smooth-window", type=int, default=1,
                        help="Display-only moving-average window (samples) for LJ+Coulomb (1 = none)")
    parser.add_argument("--smooth-window-total", type=int, default=1,
                        help="Display-only moving-average window (samples) for Total (1 = none)")
    parser.add_argument(
        "--analysis-energy-csv",
        type=Path,
        default=None,
        help="Optional avg_energies CSV for replica SEM bands (default: profile analysis_dir)",
    )
    parser.add_argument("--no-latex", action="store_true")
    parser.add_argument("--layout", choices=("paper", "stacked"), default="paper",
                        help="Figure layout: paper (2-panel b/c) or stacked (5-panel aging style)")
    parser.add_argument("--separate-panels", action="store_true",
                        help="Save (b) and (c) as two separate figure files instead of one combined figure")
    args = parser.parse_args()

    from repro_bootstrap import setup_profile

    profile = setup_profile(args, default="paper")
    if args.base_dir is not None:
        base_dir = args.base_dir
    elif profile.name != "paper":
        base_dir = profile.time_series_dir
    else:
        base_dir = DEFAULT_BASE

    coupling = parse_coupling(args.coupling)
    if profile.name == "paper" and coupling not in COUPLING_TAG:
        raise ValueError(f"Unsupported coupling {coupling}; choose from {list(COUPLING_TAG)}")
    profile.entry_for_axis_value(coupling)

    if args.output is None:
        if profile.axis == "lambda":
            lam_str = f"{coupling:g}".replace(".", "p")
            out_root = profile.figures_dir or profile.staged_root
            args.output = out_root / f"figure3_lambda_{lam_str}"
        else:
            lam_str = f"{coupling:.0e}".replace("+0", "").replace("-0", "-")
            out_dir = ARCHIVE / "final_production_run/figures/2026-01-29"
            args.output = out_dir / f"figure3_coupling_{lam_str}"

    analysis_csv = args.analysis_energy_csv
    if analysis_csv is None and profile.analysis_dir is not None:
        entry = profile.entry_for_axis_value(coupling)
        tag = getattr(entry, "analysis_tag", None)
        if tag:
            candidate = Path(profile.analysis_dir) / f"avg_energies_{tag}.csv"
            if candidate.is_file():
                analysis_csv = candidate

    data = load_figure3_data(
        base_dir,
        coupling,
        args.t_switch,
        profile=profile,
        analysis_energy_csv=analysis_csv,
    )

    t_max = args.t_max
    if args.layout == "stacked" and t_max == 2000.0:
        # Default stacked view uses the full replica-averaged time range.
        t_max = float(np.max(data["t"]))

    coupling_label = None
    if profile.axis == "lambda":
        coupling_label = f"λ={coupling:g}"

    plot_figure3(data, args.output, args.t_switch, t_max, args.t_max_b,
                 use_latex=not args.no_latex, smooth_window=args.smooth_window,
                 smooth_window_total=args.smooth_window_total,
                 combined=not args.separate_panels,
                 auto_ylim=(profile.axis == "lambda"),
                 auto_ylim_energy=(profile.axis == "lambda"),
                 layout=args.layout,
                 coupling_label=coupling_label)


if __name__ == "__main__":
    main()
