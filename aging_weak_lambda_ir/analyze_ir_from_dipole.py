#!/usr/bin/env python3
"""IR spectra from equilibrium dipole trajectories stored in observables HDF5.

In-plane convention (cavity mode in x/y):
  C(t) = <mu_x(t) mu_x(0)> + <mu_y(t) mu_y(0)>   (FFT or direct time average)

Primary IR pipeline (Efren Braun / ``compute_irspectrum.py``):
  DCT type-I on the in-plane ACF with omega^2 field factor.

Fallback IR pipeline (memspectrum MESA + FPE):
  Separate MESA solve on mu_x(t) and mu_y(t), sum spectra with omega^2.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import fftpack

BOLTZ = 1.38064852e-23
LIGHTSPEED = 299792458.0
REDUCED_PLANCK = 1.05457180013e-34
DEFAULT_TEMPERATURE_K = 100.0
DEFAULT_CAVITY_FREQUENCY_CM1 = 1560.0
DEFAULT_ACF_FRACTION = 1.0
IRMethod = str  # "auto" | "dct" | "mesa"

CAMPAIGN_LAMBDAS: list[tuple[float, str]] = [
    (0.0, "0"),
    (0.01, "0p01"),
    (0.016667, "0p016667"),
    (0.023333, "0p023333"),
    (0.03, "0p03"),
]


def _acf_component_fft(component: np.ndarray) -> np.ndarray:
    """FFT autocorrelation of one mean-subtracted scalar series."""
    signal = np.asarray(component, dtype=float)
    n_time = signal.size
    if n_time < 2:
        return np.array([])
    centered = signal - np.mean(signal)
    padded = np.zeros(2 * n_time, dtype=float)
    padded[:n_time] = centered
    fft_signal = np.fft.fft(padded)
    correlation = np.fft.ifft(fft_signal * np.conj(fft_signal)).real[:n_time]
    normalization = np.arange(n_time, 0, -1, dtype=float)
    return correlation / normalization


def _acf_component_direct(component: np.ndarray) -> np.ndarray:
    """Direct time-average autocorrelation for validation."""
    signal = np.asarray(component, dtype=float)
    n_time = signal.size
    if n_time < 2:
        return np.array([])
    centered = signal - np.mean(signal)
    acf = np.empty(n_time, dtype=float)
    for lag in range(n_time):
        acf[lag] = np.sum(centered[: n_time - lag] * centered[lag:]) / (n_time - lag)
    return acf


def dipole_acf_xy_fft(dipole: np.ndarray) -> np.ndarray:
    """In-plane dipole ACF: C(t) = C_x(t) + C_y(t) via FFT."""
    dipole = np.asarray(dipole, dtype=float)
    if dipole.ndim != 2 or dipole.shape[0] < 2 or dipole.shape[1] < 2:
        return np.array([])
    acf = np.zeros(dipole.shape[0], dtype=float)
    for dim in (0, 1):
        acf += _acf_component_fft(dipole[:, dim])
    return acf


def dipole_acf_xy_direct(dipole: np.ndarray) -> np.ndarray:
    """In-plane dipole ACF via explicit time averaging."""
    dipole = np.asarray(dipole, dtype=float)
    if dipole.ndim != 2 or dipole.shape[0] < 2 or dipole.shape[1] < 2:
        return np.array([])
    acf = np.zeros(dipole.shape[0], dtype=float)
    for dim in (0, 1):
        acf += _acf_component_direct(dipole[:, dim])
    return acf


def dipole_autocorrelation_fftconvolve(dipole: np.ndarray) -> np.ndarray:
    """Backward-compatible alias for the in-plane FFT ACF."""
    return dipole_acf_xy_fft(dipole)


def ir_spectrum_dct_from_acf(
    acf: np.ndarray,
    dt_fs: float,
) -> tuple[np.ndarray, np.ndarray]:
    """DCT IR spectrum from in-plane ACF (Efren Braun ``fft3`` convention)."""
    acf = np.asarray(acf, dtype=float)
    if acf.size < 4 or dt_fs <= 0.0:
        return np.array([]), np.array([])
    lineshape = fftpack.dct(acf, type=1)
    freq_au = np.linspace(0.0, 0.5 / dt_fs * 1.0e15, acf.size)
    freqs_cm1 = freq_au / (100.0 * LIGHTSPEED)
    spectrum = lineshape * freq_au**2
    return freqs_cm1, spectrum


def _mean_timestep_fs(times_ps: np.ndarray) -> float | None:
    """Return mean sampling interval in femtoseconds."""
    times_fs = np.asarray(times_ps, dtype=float) * 1000.0
    if times_fs.size < 2:
        return None
    dt_fs = float(np.mean(np.diff(times_fs)))
    return dt_fs if dt_fs > 0.0 else None


def _spectrum_is_valid(freqs: np.ndarray, spec: np.ndarray) -> bool:
    """Return True when a spectrum array is usable for plotting and I/O."""
    if freqs.size == 0 or spec.size == 0 or freqs.shape != spec.shape:
        return False
    if not np.all(np.isfinite(freqs)) or not np.all(np.isfinite(spec)):
        return False
    return bool(np.nanmax(spec) > 0.0)


def ir_spectrum_dct_from_dipole(
    dipole: np.ndarray,
    times_ps: np.ndarray,
    *,
    temperature_k: float = DEFAULT_TEMPERATURE_K,
    acf_fraction: float = DEFAULT_ACF_FRACTION,
) -> tuple[np.ndarray, np.ndarray]:
    """DCT IR spectrum from in-plane dipole ACF."""
    del temperature_k
    if dipole.shape[0] < 4:
        return np.array([]), np.array([])
    dt_fs = _mean_timestep_fs(times_ps)
    if dt_fs is None:
        return np.array([]), np.array([])

    autocorr_full = dipole_acf_xy_fft(dipole)
    n_acf = max(4, int(autocorr_full.size * acf_fraction))
    n_acf = min(n_acf, autocorr_full.size)
    return ir_spectrum_dct_from_acf(autocorr_full[:n_acf], dt_fs)


def _mesa_component_spectrum(
    component: np.ndarray,
    dt_s: float,
) -> tuple[np.ndarray, np.ndarray]:
    from memspectrum import MESA

    mesa = MESA()
    mesa.solve(
        np.asarray(component, dtype=float),
        method="Standard",
        optimisation_method="FPE",
    )
    freqs_hz, lineshape = mesa.spectrum(dt=dt_s, onesided=True)
    omega = 2.0 * np.pi * np.asarray(freqs_hz, dtype=float)
    freqs_cm1 = np.asarray(freqs_hz, dtype=float) / (100.0 * LIGHTSPEED)
    spectrum = np.asarray(lineshape, dtype=float) * omega**2
    return freqs_cm1, spectrum


def ir_spectrum_mesa_fpe_from_dipole_xy(
    dipole: np.ndarray,
    times_ps: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """MaxEnt IR spectrum from mu_x and mu_y separately, then summed."""
    try:
        from memspectrum import MESA  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "memspectrum is required for MESA IR analysis "
            "(pip install memspectrum)"
        ) from exc

    if dipole.shape[0] < 4 or dipole.shape[1] < 2:
        return np.array([]), np.array([])
    dt_fs = _mean_timestep_fs(times_ps)
    if dt_fs is None:
        return np.array([]), np.array([])

    dt_s = dt_fs * 1.0e-15
    spectrum_total: np.ndarray | None = None
    freqs_ref: np.ndarray | None = None
    for dim in (0, 1):
        freqs_cm1, component_spec = _mesa_component_spectrum(dipole[:, dim], dt_s)
        if spectrum_total is None:
            spectrum_total = component_spec
            freqs_ref = freqs_cm1
        elif freqs_ref.shape != freqs_cm1.shape:
            return np.array([]), np.array([])
        else:
            spectrum_total += component_spec

    assert freqs_ref is not None and spectrum_total is not None
    return freqs_ref, spectrum_total


def ir_spectrum_mesa_fpe_from_dipole(
    dipole: np.ndarray,
    times_ps: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Backward-compatible alias."""
    return ir_spectrum_mesa_fpe_from_dipole_xy(dipole, times_ps)


def ir_spectrum_from_dipole(
    dipole: np.ndarray,
    times_ps: np.ndarray,
    *,
    temperature_k: float = DEFAULT_TEMPERATURE_K,
    acf_fraction: float = DEFAULT_ACF_FRACTION,
    method: IRMethod = "auto",
) -> tuple[np.ndarray, np.ndarray, str]:
    """Compute IR spectrum with DCT, MESA/FPE, or automatic fallback."""
    if method not in {"auto", "dct", "mesa"}:
        raise ValueError(f"unsupported IR method: {method}")

    if method == "mesa":
        freqs, spec = ir_spectrum_mesa_fpe_from_dipole(dipole, times_ps)
        if not _spectrum_is_valid(freqs, spec):
            return np.array([]), np.array([]), "mesa"
        return freqs, spec, "mesa"

    freqs, spec = ir_spectrum_dct_from_dipole(
        dipole,
        times_ps,
        temperature_k=temperature_k,
        acf_fraction=acf_fraction,
    )
    if method == "dct" or _spectrum_is_valid(freqs, spec):
        return freqs, spec, "dct"

    freqs, spec = ir_spectrum_mesa_fpe_from_dipole(dipole, times_ps)
    if _spectrum_is_valid(freqs, spec):
        return freqs, spec, "mesa"
    return np.array([]), np.array([]), "mesa"


def trim_monotonic_trajectory(
    times_ps: np.ndarray,
    dipole: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Keep the first contiguous segment with strictly increasing laboratory time."""
    times_ps = np.asarray(times_ps, dtype=float)
    dipole = np.asarray(dipole, dtype=float)
    if times_ps.size < 2:
        return times_ps, dipole
    bad = np.where(np.diff(times_ps) <= 0.0)[0]
    end = int(bad[0] + 1) if bad.size else times_ps.size
    return times_ps[:end], dipole[:end]


def load_dipole_trajectory_hdf5(
    h5_path: Path,
) -> tuple[np.ndarray, np.ndarray]:
    """Load laboratory time (ps) and dipole components (e*Ang) from observables HDF5."""
    with h5py.File(h5_path, "r") as handle:
        times = np.asarray(handle["time"][:], dtype=float)
        dipole = np.asarray(
            handle["order_parameters/dipole/components"][:],
            dtype=float,
        )
    times, dipole = trim_monotonic_trajectory(times, dipole)
    if dipole.ndim != 2 or dipole.shape[1] < 3:
        raise ValueError(f"unexpected dipole shape in {h5_path}: {dipole.shape}")
    return times, dipole


def find_observables_h5(run_dir: Path, replica: int) -> Path | None:
    """Return observables HDF5 for a replica if present."""
    candidate = run_dir / f"observables_replica_{replica}.h5"
    if candidate.is_file():
        return candidate
    matches = sorted(run_dir.glob(f"**/observables_replica_{replica}.h5"))
    return matches[0] if matches else None


def resolve_run_dir(base_dir: Path, lam_tag: str, lam: float) -> Path | None:
    """Locate the coupling run directory under ``lambda{tag}/``."""
    lam_root = base_dir / f"lambda{lam_tag}"
    if not lam_root.is_dir():
        return None
    coupling_dirs = sorted(
        p for p in lam_root.iterdir() if p.is_dir() and p.name.startswith("coupling_")
    )
    if coupling_dirs:
        return coupling_dirs[-1]
    return lam_root


def ensemble_ir_spectrum(
    base_dir: Path,
    lam: float,
    lam_tag: str,
    replicas: list[int],
    *,
    temperature_k: float = DEFAULT_TEMPERATURE_K,
    acf_fraction: float = DEFAULT_ACF_FRACTION,
    method: IRMethod = "auto",
) -> tuple[np.ndarray, np.ndarray, str] | None:
    """Average IR spectra over replicas for one coupling."""
    run_dir = resolve_run_dir(base_dir, lam_tag, lam)
    if run_dir is None:
        return None
    specs: list[np.ndarray] = []
    freqs_ref: np.ndarray | None = None
    methods: list[str] = []
    for replica in replicas:
        h5_path = find_observables_h5(run_dir, replica)
        if h5_path is None:
            continue
        times, dipole = load_dipole_trajectory_hdf5(h5_path)
        freqs, spec, used_method = ir_spectrum_from_dipole(
            dipole,
            times,
            temperature_k=temperature_k,
            acf_fraction=acf_fraction,
            method=method,
        )
        if not _spectrum_is_valid(freqs, spec):
            continue
        if freqs_ref is None:
            freqs_ref = freqs
        elif freqs_ref.shape != freqs.shape:
            continue
        specs.append(spec)
        methods.append(used_method)
    if freqs_ref is None or not specs:
        return None
    resolved = "mesa" if "mesa" in methods else "dct"
    return freqs_ref, np.mean(np.stack(specs, axis=0), axis=0), resolved


def write_spectrum(path: Path, freqs: np.ndarray, spec: np.ndarray) -> None:
    """Write frequency / spectrum columns to text."""
    path.parent.mkdir(parents=True, exist_ok=True)
    data = np.column_stack([freqs, spec])
    np.savetxt(
        path,
        data,
        header="frequency_cm-1 spectrum_qm_au",
        comments="# ",
    )


def plot_average_spectra(
    base_dir: Path,
    output_pdf: Path,
    *,
    replicas: list[int],
    lambdas: list[tuple[float, str]] | None = None,
    min_freq_cm1: float = 1200.0,
    max_freq_cm1: float = 2200.0,
    temperature_k: float = DEFAULT_TEMPERATURE_K,
    acf_fraction: float = DEFAULT_ACF_FRACTION,
    method: IRMethod = "auto",
) -> None:
    """Plot normalized ensemble IR spectra for each lambda."""
    plot_lams = lambdas or CAMPAIGN_LAMBDAS
    available = [
        (lam, tag)
        for lam, tag in plot_lams
        if ensemble_ir_spectrum(
            base_dir,
            lam,
            tag,
            replicas,
            temperature_k=temperature_k,
            acf_fraction=acf_fraction,
            method=method,
        )
        is not None
    ]
    if not available:
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.text(0.5, 0.5, "No spectra available", ha="center", va="center")
        fig.savefig(output_pdf, bbox_inches="tight")
        plt.close(fig)
        return

    ncols = len(available)
    fig, axes = plt.subplots(1, ncols, figsize=(3.2 * ncols, 3.2), sharey=True)
    if ncols == 1:
        axes = [axes]
    for ax, (lam, tag) in zip(axes, available):
        result = ensemble_ir_spectrum(
            base_dir,
            lam,
            tag,
            replicas,
            temperature_k=temperature_k,
            acf_fraction=acf_fraction,
            method=method,
        )
        assert result is not None
        freqs, spec, _resolved = result
        ax.set_title(rf"$\lambda$={lam:g}")
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
        ax.set_xlim(min_freq_cm1, max_freq_cm1)
        ax.set_ylim(0.0, 1.05)
        ax.set_xlabel(r"frequency (cm$^{-1}$)")
        if lam == available[0][0]:
            ax.set_ylabel(r"$n(\omega)\alpha(\omega)$ (norm.)")
    fig.suptitle(
        f"Equilibrium IR (method={method}, ACF fraction={acf_fraction:g}, "
        f"replicas {replicas}, T={temperature_k:g} K)",
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
        help="Campaign root containing lambda* directories",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("aging_weak_lambda_ir/derived"),
        help="Directory for spectra text files and figures",
    )
    parser.add_argument(
        "--replicas",
        type=int,
        nargs="+",
        default=[0, 1],
        help="Replica indices to average (default: 0 1)",
    )
    parser.add_argument(
        "--temperature-k",
        type=float,
        default=DEFAULT_TEMPERATURE_K,
        help="Temperature for quantum correction (default: 100 K)",
    )
    parser.add_argument(
        "--acf-fraction",
        type=float,
        default=DEFAULT_ACF_FRACTION,
        help="Fraction of ACF to DCT (default: 1.0 = full trajectory)",
    )
    parser.add_argument(
        "--min-freq-cm1",
        type=float,
        default=1200.0,
    )
    parser.add_argument(
        "--max-freq-cm1",
        type=float,
        default=2200.0,
    )
    parser.add_argument(
        "--method",
        choices=("auto", "dct", "mesa"),
        default="auto",
        help=(
            "Spectrum estimator: DCT on ACF, MESA/FPE on dipole x/y, or auto "
            "(DCT with MESA fallback)"
        ),
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for lam, tag in CAMPAIGN_LAMBDAS:
        if args.method == "auto":
            method_results: list[tuple[str, np.ndarray, np.ndarray]] = []
            for spectrum_method in ("dct", "mesa"):
                result = ensemble_ir_spectrum(
                    args.base_dir,
                    lam,
                    tag,
                    args.replicas,
                    temperature_k=args.temperature_k,
                    acf_fraction=args.acf_fraction,
                    method=spectrum_method,
                )
                if result is None:
                    print(f"skip lambda={lam:g} ({spectrum_method}): no spectra")
                    continue
                freqs, spec, resolved_method = result
                method_results.append((resolved_method, freqs, spec))
                tagged_txt = args.output_dir / (
                    f"ir_spectrum_lambda{tag}_{resolved_method}.txt"
                )
                write_spectrum(tagged_txt, freqs, spec)
                print(f"wrote {tagged_txt} (method={resolved_method})")

            if not method_results:
                continue
            preferred = next(
                (item for item in method_results if item[0] == "mesa"),
                method_results[0],
            )
            resolved_method, freqs, spec = preferred
            out_txt = args.output_dir / f"ir_spectrum_lambda{tag}.txt"
            write_spectrum(out_txt, freqs, spec)
            print(f"wrote {out_txt} (primary method={resolved_method})")
            continue

        result = ensemble_ir_spectrum(
            args.base_dir,
            lam,
            tag,
            args.replicas,
            temperature_k=args.temperature_k,
            acf_fraction=args.acf_fraction,
            method=args.method,
        )
        if result is None:
            print(f"skip lambda={lam:g}: no spectra")
            continue
        freqs, spec, resolved_method = result
        out_txt = args.output_dir / f"ir_spectrum_lambda{tag}.txt"
        write_spectrum(out_txt, freqs, spec)
        print(f"wrote {out_txt} (method={resolved_method})")

    plot_average_spectra(
        args.base_dir,
        args.output_dir / "ir_spectra_smoke.pdf",
        replicas=args.replicas,
        min_freq_cm1=args.min_freq_cm1,
        max_freq_cm1=args.max_freq_cm1,
        temperature_k=args.temperature_k,
        acf_fraction=args.acf_fraction,
        method="mesa" if args.method == "auto" else args.method,
    )
    print(
        f"wrote {args.output_dir / 'ir_spectra_smoke.pdf'} "
        f"(method={'mesa' if args.method == 'auto' else args.method})"
    )
    if args.method == "auto":
        plot_average_spectra(
            args.base_dir,
            args.output_dir / "ir_spectra_smoke_dct.pdf",
            replicas=args.replicas,
            min_freq_cm1=args.min_freq_cm1,
            max_freq_cm1=args.max_freq_cm1,
            temperature_k=args.temperature_k,
            acf_fraction=args.acf_fraction,
            method="dct",
        )
        print(f"wrote {args.output_dir / 'ir_spectra_smoke_dct.pdf'}")


if __name__ == "__main__":
    main()
