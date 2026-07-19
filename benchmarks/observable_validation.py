#!/usr/bin/env python3
"""Compare PPPM vs reaction-field short trajectories on key observables."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "examples" / "05_advanced_run.py"
GSD = REPO / "square_lambda0.025_diffeq" / "prod-0.gsd"


def run_case(name: str, electrostatics: str, runtime_ps: float, out_root: Path) -> Path:
    out = out_root / name
    out.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(RUNNER),
        "--molecular-bath", "bussi",
        "--cavity-bath", "langevin",
        "--lambda-coupling", "0.03",
        "--coupling-type", "step",
        "--temperature", "100.0",
        "--frequency", "1560.0",
        "--runtime", str(runtime_ps),
        "--switch-time", "200.0",
        "--input-gsd", str(GSD),
        "--frame", "-1",
        "--device", "GPU",
        "--gpu-id", "0",
        "--fixed-timestep",
        "--timestep", "1.0",
        "--enable-energy-tracker",
        "--energy-output-period-ps", "1.0",
        "--disable-gsd",
        "--disable-temp-tracker",
        "--console-output-period-ps", "100.0",
        "--replicas", "0",
        "--seed", "42",
        "--electrostatics", electrostatics,
        "--eps-rf", "0.0",
        "--coulomb-rcut", "15.0",
    ]
    log = out / "run.log"
    env = dict(**{k: v for k, v in __import__("os").environ.items()})
    env["PYTHONPATH"] = f"/opt/hoomd-md:{env.get('PYTHONPATH', '')}"
    subprocess.run(cmd, cwd=out, env=env, check=True, stdout=log.open("w"), stderr=subprocess.STDOUT)
    for h5 in out.rglob("observables_replica_0.h5"):
        return h5
    raise FileNotFoundError(f"No observables HDF5 under {out}")


def load_series(h5_path: Path) -> dict[str, np.ndarray]:
    with h5py.File(h5_path, "r") as f:
        times = np.array(f["time"][:], dtype=float)
        valid = times > 0
        times = times[valid]
        data = {"time_ps": times}
        for key in f["energies"].keys():
            data[key] = np.array(f[f"energies/{key}"][valid], dtype=float)
        data["coulomb"] = data["ewald_short"] + data["ewald_long"]
        data["pe"] = data["total_potential"]
        data["ke"] = data["total_kinetic"]
        data["universe"] = data["universe_total"]
        data["temperature"] = data["temperature"]
        return data


def _rel_rms(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((b - a) ** 2)) / (np.std(a) + 1e-12))


def compare(pppm: dict, rf: dict, n_samples: int = 200) -> dict:
    n = min(len(pppm["time_ps"]), len(rf["time_ps"]), n_samples)
    out: dict[str, float | bool | str] = {"n_samples": n}

    for key in ("pe", "ke", "universe", "coulomb", "temperature", "lj", "harmonic"):
        a = pppm[key][:n]
        b = rf[key][:n]
        out[f"{key}_mean_diff"] = float(np.mean(b - a))
        out[f"{key}_max_abs_diff"] = float(np.max(np.abs(b - a)))
        out[f"{key}_rel_rms"] = _rel_rms(a, b)

    out["coulomb_offset_hartree"] = float(out["coulomb_mean_diff"])
    out["universe_drift_pppm"] = float(pppm["universe"][n - 1] - pppm["universe"][0])
    out["universe_drift_rf"] = float(rf["universe"][n - 1] - rf["universe"][0])
    out["rf_internal_drift_hartree"] = float(
        (rf["lj"][n - 1] - rf["lj"][0])
        + (rf["harmonic"][n - 1] - rf["harmonic"][0])
        + (rf["coulomb"][n - 1] - rf["coulomb"][0])
    )
    out["pppm_internal_drift_hartree"] = float(
        (pppm["lj"][n - 1] - pppm["lj"][0])
        + (pppm["harmonic"][n - 1] - pppm["harmonic"][0])
        + (pppm["coulomb"][n - 1] - pppm["coulomb"][0])
    )

    # Absolute Coulomb/total energies differ by construction (RF != Ewald split).
    dynamic_pass = (
        out["temperature_rel_rms"] < 0.75
        and out["ke_rel_rms"] < 0.75
        and out["lj_rel_rms"] < 1.5
        and out["harmonic_rel_rms"] < 0.25
        and abs(out["universe_drift_pppm"]) < 1.0
        and abs(out["rf_internal_drift_hartree"]) < 5.0
    )
    out["pass_dynamic_observables"] = dynamic_pass
    out["pass_heuristic"] = dynamic_pass
    out["notes"] = (
        "Identical-seed trajectories diverge because RF and PPPM apply different "
        "Coulomb forces; compare force-match residuals for force fidelity. "
        "A ~2.4 Hartree Coulomb offset in total energy is expected."
    )
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runtime-ps", type=float, default=200.0)
    ap.add_argument("--output", default="benchmarks/results/observable_validation.json")
    ap.add_argument("--skip-run", action="store_true", help="Reuse existing HDF5 outputs")
    args = ap.parse_args()

    root = Path("benchmarks/results/observable_validation")
    root.mkdir(parents=True, exist_ok=True)
    if args.skip_run:
        pppm_h5 = next((root / "pppm").rglob("observables_replica_0.h5"))
        rf_h5 = next((root / "reaction_field").rglob("observables_replica_0.h5"))
    else:
        pppm_h5 = run_case("pppm", "pppm", args.runtime_ps, root)
        rf_h5 = run_case("reaction_field", "reaction_field", args.runtime_ps, root)

    metrics = compare(load_series(pppm_h5), load_series(rf_h5))
    metrics["runtime_ps"] = args.runtime_ps
    out = Path(args.output)
    out.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(json.dumps(metrics, indent=2))
    print(f"Written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
