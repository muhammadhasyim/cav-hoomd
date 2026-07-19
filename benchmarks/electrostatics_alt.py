#!/usr/bin/env python3
"""Compare electrostatics backends for the 501-particle cavity system.

PPPM (order 6/4/3) vs cutoff-based reaction field. Reaction field removes the
FFT/charge-grid entirely (O(pairs)), giving the practical ceiling if the
long-range treatment can be approximated. Physics/accuracy tradeoff is the
user's call; this only measures throughput.
"""
import time
import numpy as np
import hoomd

GSD = "/workspace/refs/cav-hoomd/square_lambda0.025_diffeq/prod-0.gsd"
LJ_RCUT = 15.0
KT = 100.0 * 3.166811563e-6
DT = 41.341379


def base(sim):
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
    return harmonic, cell, lj


def make_sim():
    dev = hoomd.device.GPU()
    sim = hoomd.Simulation(device=dev, seed=1)
    sim.create_state_from_gsd(filename=GSD, frame=-1)
    return sim


def run(sim, forces, steps=4000):
    integrator = hoomd.md.Integrator(dt=DT, forces=forces)
    mol = hoomd.filter.Type(['O', 'N'])
    bussi = hoomd.md.methods.thermostats.Bussi(kT=KT, tau=DT * 100)
    integrator.methods.append(
        hoomd.md.methods.ConstantVolume(filter=mol, thermostat=bussi))
    sim.operations.integrator = integrator
    sim.run(200)
    t0 = time.perf_counter()
    sim.run(steps)
    return steps / (time.perf_counter() - t0)


def pppm_case(order, grid=32):
    sim = make_sim()
    harmonic, cell, lj = base(sim)
    short, lng = hoomd.md.long_range.pppm.make_pppm_coulomb_forces(
        nlist=cell, resolution=[grid, grid, grid], order=order,
        r_cut=LJ_RCUT, alpha=0.0)
    return run(sim, [harmonic, lj, short, lng])


def rf_case(eps_rf=80.0, rcut=15.0):
    sim = make_sim()
    harmonic, cell, lj = base(sim)
    rf = hoomd.md.pair.ReactionField(nlist=cell, default_r_cut=rcut)
    # Charge products handled via ReactionField params (epsilon set per type pair);
    # HOOMD ReactionField uses q_i*q_j implicitly through 'epsilon' param = q_i*q_j.
    types = ['O', 'N', 'L']
    q = {'O': -0.3, 'N': 0.3, 'L': 0.0}  # representative; sign pattern per molecule
    for i in types:
        for j in types:
            rf.params[(i, j)] = dict(epsilon=q[i] * q[j], eps_rf=eps_rf, use_charge=False)
            rf.r_cut[(i, j)] = rcut if (i != 'L' and j != 'L') else 0.0
    return run(sim, [harmonic, lj, rf])


def rf_charge_case(eps_rf=80.0, rcut=15.0):
    sim = make_sim()
    harmonic, cell, lj = base(sim)
    rf = hoomd.md.pair.ReactionField(nlist=cell, default_r_cut=rcut)
    types = ['O', 'N', 'L']
    for i in types:
        for j in types:
            rf.params[(i, j)] = dict(epsilon=1.0, eps_rf=eps_rf, use_charge=True)
            rf.r_cut[(i, j)] = rcut if (i != 'L' and j != 'L') else 0.0
    return run(sim, [harmonic, lj, rf])


def main():
    print(f"{'case':<28} {'steps/s':>9} {'ns/day':>8}")
    for label, fn in [
        ("pppm order6 (current)", lambda: pppm_case(6)),
        ("pppm order4", lambda: pppm_case(4)),
        ("pppm order3", lambda: pppm_case(3)),
        ("reaction_field use_charge", lambda: rf_charge_case(rcut=15.0)),
        ("reaction_field rc12", lambda: rf_charge_case(rcut=12.0)),
    ]:
        try:
            tps = fn()
            print(f"{label:<28} {tps:>9.1f} {tps*0.0864:>8.1f}")
        except Exception as e:
            print(f"{label:<28}  ERROR: {str(e)[:70]}")


if __name__ == "__main__":
    main()
