#!/usr/bin/env python3
"""Standalone kernel-level throughput profiler for the cavity MD physics stack.

Bypasses CavityMDSimulation entirely: builds forces/integrator directly and
times raw sim.run() with NO Python custom updaters, isolating pure GPU kernel
cost. Used to find the physics-only ceiling and A/B individual components.
"""
import time
import argparse
import numpy as np
import hoomd

GSD = "/workspace/refs/cav-hoomd/square_lambda0.025_diffeq/prod-0.gsd"
RCUT = 15.0
KT = 100.0 * 3.166811563e-6  # 100 K in Hartree
DT = 41.341379  # 1 fs in atomic units (matches production)


def build_sim(nlist_type="cell", grid=32, order=6, use_pppm=True,
              use_cavity=True, use_lj=True, use_bonds=True,
              cavity_variant="constant", buffer=1.0, rcut=RCUT,
              integrate_photon=False):
    dev = hoomd.device.GPU()
    sim = hoomd.Simulation(device=dev, seed=1)
    sim.create_state_from_gsd(filename=GSD, frame=-1)

    forces = []

    if use_cavity:
        from hoomd.cavitymd.forces import CavityForce
        omegac = 1560.0 / 219474.6313702  # cm^-1 -> Hartree
        if cavity_variant == "step":
            # Python variant: measures GIL-crossing cost on force hot path
            from hoomd.cavitymd.variants.basic import StepVariant

            class _TT:
                elapsed_time = 0.0
            coupling = StepVariant(target_value=0.03, switch_time_ps=200.0,
                                   time_tracker=_TT())
        else:
            coupling = hoomd.variant.Constant(0.03)
        cav = CavityForce(kvector=np.array([0, 0, 1]), lambda_coupling=coupling,
                          omegac=omegac)
        forces.append(cav)

    if use_bonds:
        harmonic = hoomd.md.bond.Harmonic()
        harmonic.params['O-O'] = dict(k=2 * 0.36602, r0=2.281655158)
        harmonic.params['N-N'] = dict(k=2 * 0.71625, r0=2.0743522177)
        forces.append(harmonic)

    if nlist_type == "cell":
        cell = hoomd.md.nlist.Cell(buffer=buffer, exclusions=('bond',))
    elif nlist_type == "tree":
        cell = hoomd.md.nlist.Tree(buffer=buffer, exclusions=('bond',))
    elif nlist_type == "stencil":
        cell = hoomd.md.nlist.Stencil(buffer=buffer, cell_width=rcut / 2.0,
                                      exclusions=('bond',))
    else:
        raise ValueError(nlist_type)

    if use_lj:
        lj = hoomd.md.pair.LJ(nlist=cell, mode='shift')
        lj.params[('O', 'O')] = dict(epsilon=0.00016685201, sigma=6.230426584)
        lj.r_cut[('O', 'O')] = rcut
        lj.params[('N', 'N')] = dict(epsilon=0.000083426, sigma=5.48277488)
        lj.r_cut[('N', 'N')] = rcut
        lj.params[('N', 'O')] = dict(epsilon=0.00025027802, sigma=4.9832074319)
        lj.r_cut[('N', 'O')] = rcut
        for pair in [('L', 'N'), ('N', 'L'), ('O', 'L'), ('L', 'O'), ('L', 'L')]:
            lj.params[pair] = dict(epsilon=0.0, sigma=1.0)
            lj.r_cut[pair] = 0.0
        forces.append(lj)

    if use_pppm:
        short, lng = hoomd.md.long_range.pppm.make_pppm_coulomb_forces(
            nlist=cell, resolution=[grid, grid, grid], order=order,
            r_cut=rcut, alpha=0.0)
        forces.append(short)
        forces.append(lng)

    integrator = hoomd.md.Integrator(dt=DT, forces=forces)
    mol = hoomd.filter.Type(['O', 'N'])
    bussi = hoomd.md.methods.thermostats.Bussi(kT=KT, tau=DT * 100)
    integrator.methods.append(
        hoomd.md.methods.ConstantVolume(filter=mol, thermostat=bussi))
    if integrate_photon and 'L' in sim.state.particle_types:
        cavf = hoomd.filter.Type(['L'])
        integrator.methods.append(
            hoomd.md.methods.Langevin(filter=cavf, kT=KT, default_gamma=1.0))
    sim.operations.integrator = integrator
    return sim


def bench(name, warmup=200, steps=3000, **kw):
    sim = build_sim(**kw)
    sim.run(warmup)  # JIT autotuners, nlist build
    t0 = time.perf_counter()
    sim.run(steps)
    dt = time.perf_counter() - t0
    tps = steps / dt
    ns_day = tps * DT_TO_PS * 86400 * 1e-3  # steps/s * ps/step -> ps/s ... see below
    return tps


DT_TO_PS = 0.001  # 1 fs = 0.001 ps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", type=int, default=3000)
    args = ap.parse_args()

    cases = [
        ("full_cell_g32_const",   dict(nlist_type="cell", grid=32, cavity_variant="constant")),
        ("full_cell_g32_stepvar", dict(nlist_type="cell", grid=32, cavity_variant="step")),
        ("no_cavity_cell_g32",    dict(nlist_type="cell", grid=32, use_cavity=False)),
        ("tree_g32",              dict(nlist_type="tree", grid=32)),
        ("stencil_g32",           dict(nlist_type="stencil", grid=32)),
        ("cell_g24",              dict(nlist_type="cell", grid=24)),
        ("cell_g16",              dict(nlist_type="cell", grid=16)),
        ("no_pppm_lj_bonds",      dict(nlist_type="cell", use_pppm=False)),
        ("lj_only_no_bonds_nopppm", dict(nlist_type="cell", use_pppm=False, use_bonds=False, use_cavity=False)),
        ("no_forces_integ_only",  dict(nlist_type="cell", use_pppm=False, use_bonds=False, use_cavity=False, use_lj=False)),
    ]
    print(f"{'case':<28} {'steps/s':>10} {'ns/day':>10}")
    results = {}
    for name, kw in cases:
        try:
            tps = bench(name, steps=args.steps, **kw)
            ns = tps * DT_TO_PS * 86.4  # steps/s * 0.001 ps/step * 86400 s/day / 1000 = *0.0864
            print(f"{name:<28} {tps:>10.1f} {ns:>10.1f}")
            results[name] = tps
        except Exception as e:
            print(f"{name:<28}  ERROR: {e}")
    return results


if __name__ == "__main__":
    main()
