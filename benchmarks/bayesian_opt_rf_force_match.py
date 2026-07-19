#!/usr/bin/env python3
"""Bayesian optimization of reaction-field params to match PPPM Coulomb forces."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import gsd.hoomd
import hoomd
import numpy as np
from skopt import gp_minimize
from skopt.space import Real, Categorical

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "benchmarks"))

from electrostatic_force_match import (  # noqa: E402
    GSD_DEFAULT,
    MOL_TYPES,
    _compare,
    _coulomb_forces_pppm,
    _coulomb_forces_rf,
)

DT = 41.341379


def _mol_mask(type_names, typeid):
    names = list(type_names)
    return np.array([names[t] in MOL_TYPES for t in typeid])


def cache_pppm_forces(frames_path: str, frame_indices: list[int],
                      rcut: float = 15.0, grid: int = 32, order: int = 6):
    """Evaluate and cache PPPM reference forces for selected frames."""
    cached = []
    for idx in frame_indices:
        dev = hoomd.device.GPU()
        sim = hoomd.Simulation(device=dev, seed=1)
        sim.create_state_from_gsd(filename=frames_path, frame=idx)
        cell = hoomd.md.nlist.Cell(buffer=1.0)
        forces = _coulomb_forces_pppm(cell, rcut=rcut, grid=grid, order=order)
        integrator = hoomd.md.Integrator(dt=DT, forces=forces)
        sim.operations.integrator = integrator
        sim.run(0)
        with sim.state.cpu_local_snapshot as snap:
            forces_arr = np.array(snap.particles.net_force, dtype=float)
            typeid = np.array(snap.particles.typeid)
            mask = _mol_mask(sim.state.particle_types, typeid)
        cached.append(forces_arr[mask])
    return cached


def eval_rf_cosine(cached_refs: list[np.ndarray], frames_path: str,
                   frame_indices: list[int], eps_rf: float, rcut: float) -> dict:
    """Return force-match metrics for one RF parameter set."""
    rms_list, cos_list, max_list = [], [], []
    for idx, ref in zip(frame_indices, cached_refs):
        dev = hoomd.device.GPU()
        sim = hoomd.Simulation(device=dev, seed=1)
        sim.create_state_from_gsd(filename=frames_path, frame=idx)
        cell = hoomd.md.nlist.Cell(buffer=1.0)
        forces = _coulomb_forces_rf(cell, eps_rf=eps_rf, rcut=rcut)
        integrator = hoomd.md.Integrator(dt=DT, forces=forces)
        sim.operations.integrator = integrator
        sim.run(0)
        with sim.state.cpu_local_snapshot as snap:
            cand = np.array(snap.particles.net_force, dtype=float)
            typeid = np.array(snap.particles.typeid)
            mask = _mol_mask(sim.state.particle_types, typeid)
        rms, mx, cos = _compare(ref, cand[mask])
        rms_list.append(rms)
        max_list.append(mx)
        cos_list.append(cos)
    return {
        "mean_cosine": float(np.mean(cos_list)),
        "min_cosine": float(np.min(cos_list)),
        "mean_rms": float(np.mean(rms_list)),
        "max_rms": float(np.max(max_list)),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--frames", default=GSD_DEFAULT)
    ap.add_argument("--max-frames", type=int, default=30)
    ap.add_argument("--n-calls", type=int, default=40,
                    help="Total BO evaluations (including initial random)")
    ap.add_argument("--rcut-min", type=float, default=14.0)
    ap.add_argument("--rcut-max", type=float, default=18.0)
    ap.add_argument("--output", default="benchmarks/results/rf_bayes_opt.json")
    args = ap.parse_args()

    frames_path = str(Path(args.frames))
    with gsd.hoomd.open(frames_path, "r") as f:
        n_avail = len(f)
    frame_indices = list(range(min(args.max_frames, n_avail)))

    print(f"Caching PPPM references for {len(frame_indices)} frames...")
    cached = cache_pppm_forces(frames_path, frame_indices)

    # Baseline configs
    baselines = [
        ("conducting (eps_rf=0)", 0.0, 15.0),
        ("eps_rf=80", 80.0, 15.0),
    ]
    history = []
    for label, eps, rc in baselines:
        m = eval_rf_cosine(cached, frames_path, frame_indices, eps, rc)
        m.update({"label": label, "eps_rf": eps, "rcut": rc})
        history.append(m)
        print(f"  baseline {label}: cos={m['mean_cosine']:.6f} rms={m['mean_rms']:.4e}")

    # eps_rf: 0 = HOOMD conducting limit; otherwise finite dielectric
    space = [
        Categorical([0.0, 1], name="conducting"),  # 1 => use log_eps_rf
        Real(0.8, 2.5, name="log10_eps_rf"),  # 6.3..316 when conducting=1
        Real(args.rcut_min, args.rcut_max, name="rcut"),
    ]

    def objective(x):
        conducting, log_eps, rcut = x
        eps_rf = 0.0 if conducting == 0 else float(10 ** log_eps)
        metrics = eval_rf_cosine(cached, frames_path, frame_indices, eps_rf, rcut)
        entry = {
            "eps_rf": eps_rf,
            "rcut": float(rcut),
            "mean_cosine": metrics["mean_cosine"],
            "min_cosine": metrics["min_cosine"],
            "mean_rms": metrics["mean_rms"],
            "max_rms": metrics["max_rms"],
        }
        history.append(entry)
        print(f"  BO: eps_rf={eps_rf:.4g} rcut={rcut:.2f} "
              f"cos={metrics['mean_cosine']:.6f} rms={metrics['mean_rms']:.4e}")
        return -metrics["mean_cosine"]

    print(f"\nRunning Bayesian optimization ({args.n_calls} calls)...")
    result = gp_minimize(
        objective,
        space,
        n_calls=args.n_calls,
        n_random_starts=10,
        random_state=42,
        acq_func="EI",
    )

    best_idx = int(np.argmin(result.func_vals))
    conducting, log_eps, rcut = result.x_iters[best_idx]
    best_eps = 0.0 if conducting == 0 else float(10 ** log_eps)
    best_cos = -result.fun

    # Full evaluation on all available frames for the best BO point
    print(f"\nRe-evaluating best params on all {n_avail} frames...")
    with gsd.hoomd.open(frames_path, "r") as _:
        pass
    all_indices = list(range(n_avail))
    all_cached = cache_pppm_forces(frames_path, all_indices)
    full_metrics = eval_rf_cosine(
        all_cached, frames_path, all_indices, best_eps, float(rcut))

    payload = {
        "target_cosine": 0.9999,
        "best_bayes": {
            "eps_rf": best_eps,
            "rcut": float(rcut),
            "mean_cosine_train": best_cos,
            **{f"full_{k}": v for k, v in full_metrics.items()},
        },
        "baselines": history[:2],
        "all_evaluations": history[2:],
        "feasible_0.9999": full_metrics["mean_cosine"] >= 0.9999,
    }

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print("\n=== Results ===")
    print(f"Best BO params: eps_rf={best_eps:.6g}, rcut={rcut:.3f}")
    print(f"  train ({len(frame_indices)} frames): cos={best_cos:.6f}")
    print(f"  full  ({n_avail} frames): cos={full_metrics['mean_cosine']:.6f}, "
          f"rms={full_metrics['mean_rms']:.4e}")
    print(f"Target 0.9999 reachable: {payload['feasible_0.9999']}")
    print(f"Written to {out}")
    return 0 if payload["feasible_0.9999"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
