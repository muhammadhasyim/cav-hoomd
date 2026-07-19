#!/usr/bin/env python3
"""Compare reaction-field vs PPPM electrostatic forces over many configurations."""
from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass, asdict
from pathlib import Path

import gsd.hoomd
import hoomd
import numpy as np

GSD_DEFAULT = "/workspace/refs/cav-hoomd/benchmarks/results/force_match_frames.gsd"
DT = 41.341379
MOL_TYPES = {'O', 'N'}


@dataclass
class ForceMatchResult:
    eps_rf: float
    coulomb_rcut: float
    mean_rms: float
    max_rms: float
    mean_cosine: float
    min_cosine: float
    steps_per_s: float | None = None


def _mol_mask(type_names, typeid):
    names = list(type_names)
    return np.array([names[t] in MOL_TYPES for t in typeid])


def _coulomb_forces_pppm(nlist, rcut=15.0, grid=32, order=6):
    short, lng = hoomd.md.long_range.pppm.make_pppm_coulomb_forces(
        nlist=nlist,
        resolution=[grid, grid, grid],
        order=order,
        r_cut=rcut,
        alpha=0.0,
    )
    return [short, lng]


def _coulomb_forces_rf(nlist, eps_rf=0.0, rcut=15.0):
    rf = hoomd.md.pair.ReactionField(nlist=nlist, default_r_cut=rcut)
    for type_i in ['O', 'N', 'L']:
        for type_j in ['O', 'N', 'L']:
            if type_i == 'L' or type_j == 'L':
                rf.params[(type_i, type_j)] = dict(
                    epsilon=1.0, eps_rf=eps_rf, use_charge=True)
                rf.r_cut[(type_i, type_j)] = 0.0
            else:
                rf.params[(type_i, type_j)] = dict(
                    epsilon=1.0, eps_rf=eps_rf, use_charge=True)
                rf.r_cut[(type_i, type_j)] = rcut
    return [rf]


def _eval_coulomb_forces(frame_path: str, frame_idx: int, method: str,
                           eps_rf=0.0, rcut=15.0, grid=32, order=6):
    dev = hoomd.device.GPU()
    sim = hoomd.Simulation(device=dev, seed=1)
    sim.create_state_from_gsd(filename=frame_path, frame=frame_idx)
    cell = hoomd.md.nlist.Cell(buffer=1.0)
    if method == 'pppm':
        forces = _coulomb_forces_pppm(cell, rcut=rcut, grid=grid, order=order)
    else:
        forces = _coulomb_forces_rf(cell, eps_rf=eps_rf, rcut=rcut)

    integrator = hoomd.md.Integrator(dt=DT, forces=forces)
    sim.operations.integrator = integrator
    sim.run(0)
    with sim.state.cpu_local_snapshot as snap:
        forces_arr = np.array(snap.particles.net_force, dtype=float)
        typeid = np.array(snap.particles.typeid)
        mask = _mol_mask(sim.state.particle_types, typeid)
    return forces_arr[mask]


def _compare(ref: np.ndarray, cand: np.ndarray):
    diff = cand - ref
    rms = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))
    max_err = float(np.max(np.linalg.norm(diff, axis=1)))
    ref_norm = np.linalg.norm(ref, axis=1)
    cand_norm = np.linalg.norm(cand, axis=1)
    denom = ref_norm * cand_norm
    valid = denom > 1e-12
    if not np.any(valid):
        cos = 1.0
    else:
        cos = float(np.mean(np.sum(ref[valid] * cand[valid], axis=1) / denom[valid]))
    return rms, max_err, cos


def _bench_rf(rcut=15.0, eps_rf=0.0, steps=4000):
    dev = hoomd.device.GPU()
    sim = hoomd.Simulation(device=dev, seed=1)
    sim.create_state_from_gsd(
        filename="/workspace/refs/cav-hoomd/square_lambda0.025_diffeq/prod-0.gsd",
        frame=-1,
    )
    cell = hoomd.md.nlist.Cell(buffer=1.0, exclusions=('bond',))
    harmonic = hoomd.md.bond.Harmonic()
    harmonic.params['O-O'] = dict(k=2 * 0.36602, r0=2.281655158)
    harmonic.params['N-N'] = dict(k=2 * 0.71625, r0=2.0743522177)
    lj = hoomd.md.pair.LJ(nlist=cell, mode='shift')
    lj_rcut = 15.0
    lj.params[('O', 'O')] = dict(epsilon=0.00016685201, sigma=6.230426584)
    lj.r_cut[('O', 'O')] = lj_rcut
    lj.params[('N', 'N')] = dict(epsilon=0.000083426, sigma=5.48277488)
    lj.r_cut[('N', 'N')] = lj_rcut
    lj.params[('N', 'O')] = dict(epsilon=0.00025027802, sigma=4.9832074319)
    lj.r_cut[('N', 'O')] = lj_rcut
    for pair in [('L', 'N'), ('N', 'L'), ('O', 'L'), ('L', 'O'), ('L', 'L')]:
        lj.params[pair] = dict(epsilon=0.0, sigma=1.0)
        lj.r_cut[pair] = 0.0
    rf = _coulomb_forces_rf(cell, eps_rf=eps_rf, rcut=rcut)[0]
    integrator = hoomd.md.Integrator(dt=DT, forces=[harmonic, lj, rf])
    mol = hoomd.filter.Type(['O', 'N'])
    bussi = hoomd.md.methods.thermostats.Bussi(kT=KT, tau=DT * 100)
    integrator.methods.append(
        hoomd.md.methods.ConstantVolume(filter=mol, thermostat=bussi))
    sim.operations.integrator = integrator
    sim.run(200)
    t0 = time.perf_counter()
    sim.run(steps)
    return steps / (time.perf_counter() - t0)


KT = 100.0 * 3.166811563e-6


def run_force_match(frames_path: str, max_frames: int | None, bench: bool):
    with gsd.hoomd.open(frames_path, 'r') as f:
        n_frames = len(f) if max_frames is None else min(len(f), max_frames)

    cell = hoomd.md.nlist.Cell(buffer=1.0)
    _ = cell  # nlist built per-frame inside _eval_coulomb_forces

    candidates = [
        (0.0, 12.0),
        (0.0, 15.0),
        (80.0, 15.0),
        (1.0, 15.0),
    ]

    results: list[ForceMatchResult] = []
    for eps_rf, rcut in candidates:
        rms_list = []
        max_list = []
        cos_list = []
        for idx in range(n_frames):
            try:
                ref = _eval_coulomb_forces(frames_path, idx, 'pppm', rcut=15.0)
                cand = _eval_coulomb_forces(
                    frames_path, idx, 'reaction_field', eps_rf=eps_rf, rcut=rcut)
            except RuntimeError:
                continue
            rms, mx, cos = _compare(ref, cand)
            rms_list.append(rms)
            max_list.append(mx)
            cos_list.append(cos)
        if not rms_list:
            continue
        tps = _bench_rf(rcut=rcut, eps_rf=eps_rf) if bench else None
        results.append(ForceMatchResult(
            eps_rf=eps_rf,
            coulomb_rcut=rcut,
            mean_rms=float(np.mean(rms_list)),
            max_rms=float(np.max(max_list)),
            mean_cosine=float(np.mean(cos_list)),
            min_cosine=float(np.min(cos_list)),
            steps_per_s=tps,
        ))

    results.sort(key=lambda r: r.mean_rms)
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--frames', default=GSD_DEFAULT)
    ap.add_argument('--max-frames', type=int, default=None)
    ap.add_argument('--bench', action='store_true', help='Include throughput for each RF config')
    ap.add_argument('--output', default='benchmarks/results/force_match_results.json')
    args = ap.parse_args()

    frames = Path(args.frames)
    if not frames.exists():
        raise SystemExit(
            f"Missing {frames}. Run benchmarks/generate_force_match_frames.py first."
        )

    results = run_force_match(str(frames), args.max_frames, args.bench)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = [asdict(r) for r in results]
    out.write_text(json.dumps(payload, indent=2), encoding='utf-8')

    print(f"{'eps_rf':>8} {'rcut':>6} {'mean_rms':>12} {'max_rms':>12} "
          f"{'mean_cos':>10} {'min_cos':>10} {'tps':>10}")
    for r in results:
        tps = f"{r.steps_per_s:.1f}" if r.steps_per_s else 'n/a'
        print(f"{r.eps_rf:8.1f} {r.coulomb_rcut:6.1f} {r.mean_rms:12.6e} "
              f"{r.max_rms:12.6e} {r.mean_cosine:10.6f} {r.min_cosine:10.6f} {tps:>10}")
    print(f"\nBest match: eps_rf={results[0].eps_rf}, rcut={results[0].coulomb_rcut}")
    print(f"Results written to {out}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
