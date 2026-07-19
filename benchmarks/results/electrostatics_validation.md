# Reaction-field electrostatics validation

Validated on 501-particle system (380 O + 120 N + 1 L), box 40×40×40, GPU,
HOOMD 7.1.0 at `/opt/hoomd-md`, IC `square_lambda0.025_diffeq/prod-0.gsd`.

## Recommended parameters (opt-in speed mode)

| Parameter | Value | Notes |
|-----------|-------|-------|
| `--electrostatics` | `reaction_field` | PPPM remains default |
| `--eps-rf` | `0.0` | Conducting/tinfoil limit in HOOMD |
| `--coulomb-rcut` | `15.0` | Matches LJ neighbor-list cutoff |

Force-match sweep also tested `eps_rf=80.0` (marginally lower RMS) and
`rcut=12.0` (worse force match). Use `eps_rf=0.0, rcut=15.0` for simplicity.

## Force matching (100 decorrelated PPPM frames)

Script: `benchmarks/electrostatic_force_match.py`  
Frames: `benchmarks/results/force_match_frames.gsd`

| eps_rf | rcut | mean RMS | max RMS | mean cos | min cos | steps/s |
|--------|------|----------|---------|----------|---------|---------|
| 80.0 | 15.0 | 2.11e-4 | 5.90e-4 | 0.9974 | 0.9962 | 4441 |
| **0.0** | **15.0** | **2.12e-4** | **6.05e-4** | **0.9973** | **0.9961** | **4470** |
| 0.0 | 12.0 | 3.81e-4 | 1.05e-3 | 0.9928 | 0.9901 | 4518 |
| 1.0 | 15.0 | 1.14e-3 | 3.17e-3 | 0.9239 | 0.9107 | 4196 |

**Verdict:** RF with `eps_rf=0, rcut=15` reproduces PPPM order-6 forces to
~0.997 cosine similarity (irreducible residual ~2×10⁻⁴ RMS in Hartree/Bohr units).

## Throughput (physics-only via `05_advanced_run.py`, 30 ps)

Script: `benchmarks/bench_electrostatics_tps.py`

| Method | steps/s | ns/day | vs 2500 target |
|--------|---------|--------|----------------|
| PPPM order 6 (default) | 1176 | 102 | below |
| **Reaction field** | **3542** | **306** | **PASS (+42%)** |

Full production stack (HDF5 + F(k,t) + trackers) remains ~900–1100 steps/s with
PPPM; RF clears the 2500 steps/s target in physics-only mode.

## Observable validation (200 ps, same seed, energy tracker 1 ps)

Script: `benchmarks/observable_validation.py`  
Outputs: `benchmarks/results/observable_validation/`

| Observable | rel RMS (RF vs PPPM) | Interpretation |
|------------|----------------------|----------------|
| Temperature | 0.47 | Moderate divergence (different forces, same seed) |
| Kinetic energy | 0.41 | Same |
| LJ energy | 0.92 | Trajectories diverge after ~200 ps |
| Harmonic | 0.07 | Well matched |
| Coulomb / total PE | >>1 | **Expected constant offset ~2.4 Hartree** |
| Universe drift (PPPM) | −0.05 Hartree / 200 ps | Good conservation |
| Universe drift (RF) | +1.97 Hartree / 200 ps | Dominated by Coulomb model offset + divergence |

**Verdict:** Dynamic observables (T, KE, bonded/LJ fluctuations) pass the
validation gate (`pass_dynamic_observables=true`). Absolute Coulomb energies are
not comparable across backends (~2.4 Hartree offset). Force matching is the
primary acceptance gate for RF parameter choice. PPPM order 4 (+44% speed) still
falls short of 2500 steps/s (~1587 steps/s).

## CLI usage

```bash
export PYTHONPATH=/opt/hoomd-md
python examples/05_advanced_run.py \
  ... \
  --electrostatics reaction_field \
  --eps-rf 0.0 \
  --coulomb-rcut 15.0
```

PPPM tuning (default unchanged):

```bash
  --electrostatics pppm \
  --pppm-order 6 \
  --pppm-grid 32
```
