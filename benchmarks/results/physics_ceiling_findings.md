# Physics-only throughput investigation (501 particles, single GPU)

System: 40x40x40 box, 380 O + 120 N + 1 photon (L), r_cut=15 (2.4 sigma for O-O),
PPPM 32^3 order 6, LJ + harmonic bonds + Bussi(mol) + Langevin(cavity).

## Where the time goes (tracker-free, raw sim.run)

| stack                      | steps/s | ns/day |
|----------------------------|---------|--------|
| integrator only (no forces)| 26,305  | 2,273  |
| LJ + bonds (no PPPM)       |  6,649  |   574  |
| full physics (PPPM o6)     |  1,081  |    93  |

=> PPPM electrostatics is ~100% of the bottleneck (~775 us/step).
   nlist type (Cell/Tree/Stencil), grid size (16/24/32), and the cavity
   Python variant make ZERO measurable difference. FFT is negligible;
   cost is charge spread/interpolate (order^3 pts/particle) + real-space erfc.

## PPPM tuning

| config                | steps/s | ns/day | RMS force err |
|-----------------------|---------|--------|---------------|
| order 6 (current)     |  1,099  |    95  | 2.3e-8        |
| order 5               |  1,361  |   118  | 8.9e-8        |
| order 4               |  1,587  |   137  | 4.3e-7        |
| order 3               |  1,756  |   152  | 2.7e-6        |
| order4 grid48 rc10    |  1,572  |   136  | 5.2e-7        |
| order4 grid64 rc8     |  1,558  |   135  | 4.7e-7        |

- Reducing PPPM real-space r_cut gives NOTHING: the neighbor list is pinned at
  15 by the LJ cutoff (shared list), and box=40 => cell list is ~2x2x2 (all-pairs).
- Interpolation ORDER is the only PPPM knob that matters.

## Alternative electrostatics (cutoff-based, no FFT/grid)

| method                     | steps/s | ns/day |
|----------------------------|---------|--------|
| PPPM order 6               |  1,099  |    95  |
| reaction field (rc15)      |  4,471  |   386  |
| reaction field (rc12)      |  4,546  |   393  |

Reaction field EXCEEDS the ~2500 steps/s target (4x speedup) but replaces the
Ewald long-range treatment with a dielectric-continuum approximation beyond the
cutoff. This changes electrostatics/dielectric response and must be validated
against the science (IR spectra, aging dynamics, cavity dipole coupling) before
use in production.

## Recommendations (in order of safety)

1. PPPM order 6 -> 4: +44% (95 -> 137 ns/day), RMS 4e-7. Drop-in, same method.
2. PPPM order 6 -> 3: +59% (95 -> 152 ns/day), RMS 2.7e-6. Still accurate.
3. Reaction field: 4x (95 -> ~390 ns/day, >2500 steps/s). PHYSICS CHANGE - validate.

To exceed 2500 steps/s with proper Ewald long-range is not feasible for N=501 on
this GPU: even order-3 PPPM leaves ~420 us/step of irreducible solver overhead.
The 2500 target is reachable only by dropping the FFT-based long-range solver.
