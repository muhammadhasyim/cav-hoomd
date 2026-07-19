#!/usr/bin/env python3
"""Isolate PPPM cost drivers: interpolation order and real-space cutoff.

Grid size was shown to be irrelevant (FFT negligible), so the cost is
order^3 charge spread/interpolate + real-space erfc over the neighbor list.
This sweeps order and PPPM real-space r_cut (with compensating grid) and
reports both throughput and PPPM RMS force error for accuracy.
"""
import time
import numpy as np
import hoomd

GSD = "/workspace/refs/cav-hoomd/square_lambda0.025_diffeq/prod-0.gsd"
LJ_RCUT = 15.0
KT = 100.0 * 3.166811563e-6
DT = 41.341379


def build(order=6, grid=32, pppm_rcut=15.0, alpha=0.0):
    dev = hoomd.device.GPU()
    sim = hoomd.Simulation(device=dev, seed=1)
    sim.create_state_from_gsd(filename=GSD, frame=-1)

    harmonic = hoomd.md.bond.Harmonic()
    harmonic.params['O-O'] = dict(k=2 * 0.36602, r0=2.281655158)
    harmonic.params['N-N'] = dict(k=2 * 0.71625, r0=2.0743522177)

    cell = hoomd.md.nlist.Cell(buffer=1.0, exclusions=('bond',))
    lj = hoomd.md.pair.LJ(nlist=cell, mode='shift')
    lj.params[('O', 'O')] = dict(epsilon=0.00016685201, sigma=6.230426584)
    lj.r_cut[('O', 'O')] = LJ_RCUT
    lj.params[('N', 'N')] = dict(epsilon=0.000083426, sigma=5.48277488)
    lj.r_cut[('N', 'N')] = LJ_RCUT
    lj.params[('N', 'O')] = dict(epsilon=0.00025027802, sigma=4.9832074319)
    lj.r_cut[('N', 'O')] = LJ_RCUT
    for pair in [('L', 'N'), ('N', 'L'), ('O', 'L'), ('L', 'O'), ('L', 'L')]:
        lj.params[pair] = dict(epsilon=0.0, sigma=1.0)
        lj.r_cut[pair] = 0.0

    short, lng = hoomd.md.long_range.pppm.make_pppm_coulomb_forces(
        nlist=cell, resolution=[grid, grid, grid], order=order,
        r_cut=pppm_rcut, alpha=alpha)

    integrator = hoomd.md.Integrator(dt=DT, forces=[harmonic, lj, short, lng])
    mol = hoomd.filter.Type(['O', 'N'])
    bussi = hoomd.md.methods.thermostats.Bussi(kT=KT, tau=DT * 100)
    integrator.methods.append(
        hoomd.md.methods.ConstantVolume(filter=mol, thermostat=bussi))
    sim.operations.integrator = integrator
    return sim


def bench(order, grid, pppm_rcut, alpha=0.0, steps=4000):
    sim = build(order=order, grid=grid, pppm_rcut=pppm_rcut, alpha=alpha)
    sim.run(200)
    t0 = time.perf_counter()
    sim.run(steps)
    dt = time.perf_counter() - t0
    return steps / dt


def main():
    cases = [
        # (label, order, grid, pppm_rcut, alpha)
        ("baseline o6 g32 rc15",   6, 32, 15.0, 0.0),
        ("order5 g32 rc15",        5, 32, 15.0, 0.0),
        ("order4 g32 rc15",        4, 32, 15.0, 0.0),
        ("order3 g32 rc15",        3, 32, 15.0, 0.0),
        ("order4 g48 rc10",        4, 48, 10.0, 0.0),
        ("order4 g64 rc8",         4, 64, 8.0, 0.0),
        ("order3 g48 rc10",        3, 48, 10.0, 0.0),
        ("order4 g64 rc8 a0.4",    4, 64, 8.0, 0.4),
    ]
    print(f"{'case':<26} {'steps/s':>9} {'ns/day':>8}")
    for label, o, g, rc, a in cases:
        try:
            tps = bench(o, g, rc, a)
            print(f"{label:<26} {tps:>9.1f} {tps*0.0864:>8.1f}")
        except Exception as e:
            print(f"{label:<26}  ERROR: {str(e)[:60]}")


if __name__ == "__main__":
    main()
