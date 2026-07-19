#!/usr/bin/env python3
"""Generate decorrelated PPPM trajectory frames for electrostatic force matching."""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import gsd.hoomd
import hoomd
import numpy as np

GSD_DEFAULT = "/workspace/refs/cav-hoomd/square_lambda0.025_diffeq/prod-0.gsd"
KT = 100.0 * 3.166811563e-6
DT = 41.341379


def build_physics_sim(input_gsd: str, frame: int = -1):
    dev = hoomd.device.GPU()
    sim = hoomd.Simulation(device=dev, seed=1)
    sim.create_state_from_gsd(filename=input_gsd, frame=frame)

    harmonic = hoomd.md.bond.Harmonic()
    harmonic.params['O-O'] = dict(k=2 * 0.36602, r0=2.281655158)
    harmonic.params['N-N'] = dict(k=2 * 0.71625, r0=2.0743522177)
    cell = hoomd.md.nlist.Cell(buffer=1.0, exclusions=('bond',))
    lj = hoomd.md.pair.LJ(nlist=cell, mode='shift')
    rcut = 15.0
    lj.params[('O', 'O')] = dict(epsilon=0.00016685201, sigma=6.230426584)
    lj.r_cut[('O', 'O')] = rcut
    lj.params[('N', 'N')] = dict(epsilon=0.000083426, sigma=5.48277488)
    lj.r_cut[('N', 'N')] = rcut
    lj.params[('N', 'O')] = dict(epsilon=0.00025027802, sigma=4.9832074319)
    lj.r_cut[('N', 'O')] = rcut
    for pair in [('L', 'N'), ('N', 'L'), ('O', 'L'), ('L', 'O'), ('L', 'L')]:
        lj.params[pair] = dict(epsilon=0.0, sigma=1.0)
        lj.r_cut[pair] = 0.0

    short, lng = hoomd.md.long_range.pppm.make_pppm_coulomb_forces(
        nlist=cell, resolution=[32, 32, 32], order=6, r_cut=rcut, alpha=0.0)

    integrator = hoomd.md.Integrator(dt=DT, forces=[harmonic, lj, short, lng])
    mol = hoomd.filter.Type(['O', 'N'])
    bussi = hoomd.md.methods.thermostats.Bussi(kT=KT, tau=DT * 100)
    integrator.methods.append(
        hoomd.md.methods.ConstantVolume(filter=mol, thermostat=bussi))
    sim.operations.integrator = integrator
    return sim


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input-gsd', default=GSD_DEFAULT)
    ap.add_argument('--output', default='benchmarks/results/force_match_frames.gsd')
    ap.add_argument('--n-frames', type=int, default=100)
    ap.add_argument('--steps-between-frames', type=int, default=50)
    args = ap.parse_args()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        out.unlink()

    sim = build_physics_sim(args.input_gsd)
    writer = hoomd.write.GSD(
        trigger=hoomd.trigger.Periodic(args.steps_between_frames),
        filename=str(out),
        mode='wb',
    )
    sim.operations.writers.append(writer)
    total_steps = args.n_frames * args.steps_between_frames
    print(f"Generating {args.n_frames} frames ({total_steps} steps) -> {out}")
    t0 = time.perf_counter()
    sim.run(total_steps)
    elapsed = time.perf_counter() - t0
    with gsd.hoomd.open(out, 'r') as f:
        print(f"Wrote {len(f)} frames in {elapsed:.1f}s")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
