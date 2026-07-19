# FOR_AGENTS.md — `aging_weak_lambda` cavity MD workflow (cav-hoomd)

This document is for agents (and humans) continuing work on the **weak-coupling aging campaign** in this repository. It summarizes the scientific objective, directory layout, how to launch and monitor runs, known pitfalls, and how to benchmark `CavityForce` GPU performance.

---

## Scientific objective

Reproduce and extend the **weak-coupling aging** study: supercooled molecular liquid (500 O/N particles) coupled to a single cavity photon (`L`) at **100 K**, with a **step turn-on** of dimensionless coupling **λ** at **200 ps** into a **2500 ps** production run.

Campaign λ values:

| λ | Directory tag |
|---|----------------|
| 0.0 | `lambda0` |
| 0.01 | `lambda0p01` |
| 0.016667 | `lambda0p016667` |
| 0.023333 | `lambda0p023333` |
| 0.03 | `lambda0p03` |

Fixed physical settings (production protocol):

- **Temperature:** 100 K
- **Cavity frequency:** 1560 cm⁻¹ (`--frequency 1560.0`)
- **Switch time:** 200 ps (`--switch-time 200.0`, step coupling)
- **Runtime:** 2500 ps
- **Timestep:** 1 fs fixed (`--fixed-timestep --timestep 1.0`)
- **Molecular bath:** Bussi NVT
- **Cavity bath:** Langevin on the `L` particle
- **F(k,t):** enabled in production (`|k| = 6.02` a.u. for paper Fig. 2; see `run_packed_aging_task.py`)
- **GSD trajectories:** disabled in production (`--disable-gsd`); observables go to HDF5
- **Electrostatics:** PPPM order 6 (default). Opt-in reaction field for ~3× speed when force fidelity is acceptable (see below).

Reference profile for analysis staging: `reproduction/profiles/aging_weak_lambda.yaml`.

---

## Repository map (what matters for runs)

```
/scratch/mh7373/projects/cav-hoomd/
├── build_install.sh              # Rebuild native plugin after C++/CUDA changes
├── examples/
│   ├── 05_advanced_run.py        # Main CLI entry point for one replica
│   ├── init-0.gsd                # Default production IC (molecular + cavity L)
│   └── slurm/
│       ├── submit_aging_campaign.sh      # SLURM packed campaign (fixed dt)
│       ├── submit_aging_adaptive.sh      # Adaptive-dt variant (separate output tree)
│       ├── run_packed_aging_task.py     # Runs 1–2 replicas per manifest row
│       └── aging_campaign_status.py      # Completion checks, cleanup, resume helpers
├── aging_weak_lambda/            # Production campaign outputs
│   └── lambda{tag}/cavity_coupling_{lam}_switch_200.0ps/
│       ├── observables_replica_{N}.h5    # Primary production artifact
│       ├── prod-{N}.gsd                  # Only if GSD enabled
│       └── master_fskt_ref_*.txt         # Derived F(k,t) analysis
├── aging_weak_lambda_ic_library/   # Molecular-only equilibration / IC generation
│   ├── init_molecular_only.gsd     # 500-particle molecular seed (no L)
│   ├── ic_checkpoints/             # Per-replica continue ICs (single-frame GSDs)
│   ├── replicas/replica_{N}/no_cavity/prod-{N}.gsd
│   └── run_local_ic_gsd_library.sh
├── src/                            # CavityForce CPU/GPU native code
├── src/pytest/test_cavity_force_gpu_optimization.py
├── reproduction/                   # Figure re-plotting from staged data
└── third_party_refs/openmm-lm/     # OpenMM reference implementation (comparison)
```

---

## Environment and plugin install

**Python:** `/scratch/mh7373/miniforge3/envs/hoomd/bin/python` (HOOMD 5.4.0 in current use).

After any change to `src/CavityForce*.cc|cu|h` or Python plugin code:

```bash
cd /scratch/mh7373/projects/cav-hoomd
./build_install.sh          # GPU build (default)
# ./build_install.sh --no-gpu
```

**Do not** `import` from source tree for tests; always reinstall then run pytest.

Regression tests for GPU `CavityForce`:

```bash
/scratch/mh7373/miniforge3/envs/hoomd/bin/python -m pytest \
  src/pytest/test_cavity_force_gpu_optimization.py -q
```

Collect errors from unrelated tests with `--ignore=tests/test_lqg_controller.py` if needed.

---

## Workflow A — Production cavity aging runs

### Launch (SLURM)

```bash
cd /scratch/mh7373/projects/cav-hoomd/examples/slurm
./submit_aging_campaign.sh [--test] [--dry-run] [--lam 0.03]
```

- **Output base:** `/scratch/mh7373/projects/cav-hoomd/aging_weak_lambda`
- **Input GSD:** `examples/init-0.gsd`
- **Packed runner:** `run_packed_aging_task.py` (up to 2 replicas per array task)
- **Target:** 500 valid replicas per λ (configurable via `--target-valid`)

Canonical per-replica command is built in `build_simulation_command()` inside `run_packed_aging_task.py`. Use that as ground truth when debugging CLI flags.

### Output layout

For λ = 0.03:

```
aging_weak_lambda/lambda0p03/cavity_coupling_3.000000e-02_switch_200.0ps/
```

Note: directory name uses **λ formatted in scientific notation**, not ε.

### Completion criterion

`examples/slurm/aging_campaign_status.py` treats a replica as **complete** when:

```text
observables_replica_{N}.h5  >=  ~1.38 GB  (COMPLETE_MIN_BYTES)
```

Check campaign status:

```bash
/scratch/mh7373/miniforge3/envs/hoomd/bin/python \
  examples/slurm/aging_campaign_status.py scan \
  --output-base /scratch/mh7373/projects/cav-hoomd/aging_weak_lambda
```

### Local / single-GPU smoke test

```bash
cd /scratch/mh7373/projects/cav-hoomd/examples
/scratch/mh7373/miniforge3/envs/hoomd/bin/python 05_advanced_run.py \
  --molecular-bath bussi --cavity-bath langevin \
  --lambda-coupling 0.03 --coupling-type step \
  --temperature 100.0 --frequency 1560.0 \
  --runtime 2500.0 --switch-time 200.0 \
  --input-gsd init-0.gsd --frame -1 \
  --device GPU --gpu-id 0 \
  --fixed-timestep --timestep 1.0 \
  --enable-energy-tracker --energy-output-period-ps 1.0 \
  --enable-fkt --fkt-kmag 6.02 --fkt-ref-interval 200.0 --fkt-max-refs 13 \
  --disable-gsd --console-output-period-ps 100.0 \
  --replicas 0 --seed 1
```

### Opt-in reaction-field electrostatics (speed mode)

PPPM is the production default. For physics-only throughput targets (~2500 steps/s),
use cutoff-based reaction field (validated 2026-07-18 on 501-particle GPU stack):

```bash
python examples/05_advanced_run.py \
  ... \
  --electrostatics reaction_field \
  --eps-rf 0.0 \
  --coulomb-rcut 15.0
```

| Flag | Default | Purpose |
|------|---------|---------|
| `--electrostatics {pppm,reaction_field}` | `pppm` | Coulomb backend |
| `--pppm-order` | `6` | PPPM interpolation order |
| `--pppm-grid` | `32` | PPPM mesh resolution |
| `--eps-rf` | `0.0` | RF dielectric (`0` = conducting/tinfoil limit in HOOMD) |
| `--coulomb-rcut` | `15.0` | Real-space Coulomb cutoff (Å) |

**Benchmarks:** `benchmarks/electrostatic_force_match.py`, `benchmarks/bench_electrostatics_tps.py`, `benchmarks/observable_validation.py`. Results in `benchmarks/results/electrostatics_validation.md`.

**Validated RF performance:** ~3540 steps/s physics-only vs ~1176 for PPPM order 6; force-match cosine ~0.997 vs PPPM reference. Absolute Coulomb energies differ by ~2.4 Hartree (expected); use force matching + dynamic observables, not raw Ewald energy equality, when judging RF.

---

## Workflow B — Molecular-only IC library (no cavity)

Purpose: generate **100 K molecular equilibrated structures** for future cavity runs or analysis. Eight replicas, **100 ns**, GSD every **200 ps** → **500 frames/replica** when append mode is correct.

### Launcher

```bash
cd /scratch/mh7373/projects/cav-hoomd/aging_weak_lambda_ic_library
bash run_local_ic_gsd_library.sh
```

Environment overrides: `N_REPLICAS`, `PER_GPU`, `RUNTIME_PS`, `GSD_PERIOD_PS`, `TEMP_K`.

### IC chain (current state as of 2026-07-18)

1. **Seed:** `init_molecular_only.gsd` (500 particles, types O/N only)
2. **First truncated runs:** accidentally used `--truncate-gsd` → only latest frame kept
3. **Checkpoints saved:** `ic_checkpoints/replica_{0..7}.gsd` (2.4–4.8 ns equilibrated at 100 K)
4. **Restart:** script now uses checkpoints as `--input-gsd` and **does not truncate**

Backups of truncated outputs: `replicas/replica_*/no_cavity/prod-*.truncated_backup.gsd`.

### Monitor frame accumulation

```bash
/scratch/mh7373/miniforge3/envs/hoomd/bin/python - <<'PY'
import gsd.hoomd
from pathlib import Path
root = Path("/scratch/mh7373/projects/cav-hoomd/aging_weak_lambda_ic_library/replicas")
for r in range(8):
    p = root / f"replica_{r}" / "no_cavity" / f"prod-{r}.gsd"
    if p.exists():
        with gsd.hoomd.open(p, "r") as f:
            print(f"replica_{r}: frames={len(f)} size={p.stat().st_size}")
PY
```

Expect **file size to grow** beyond ~32 KB and **frame count to increase** every 200 ps. Target: **500 frames** per replica (4000 total).

---

## Critical GSD pitfall — never truncate trajectory libraries

HOOMD `hoomd.write.GSD(truncate=True)` **truncates and writes frame 0 on every trigger**, not just at job start.

**Wrong for libraries:**

```bash
--truncate-gsd   # DO NOT use when accumulating trajectories
```

**Correct:** omit `--truncate-gsd` (default `truncate=False`, append mode).

Log line to verify:

```text
GSD truncate mode: False (append to existing file)
```

This is **not** a flushing issue: truncated files stay ~32 KB with exactly one readable frame.

---

## CavityForce GPU performance (agent benchmark guide)

### Objective

Production throughput should be **practically similar** with vs without cavity coupling. The bottleneck was the old GPU `CavityForce` path (host sync, per-step photon search, bad reduction). That has been rewritten; validate on each new GPU node.

### What NOT to benchmark

A minimal 501-particle HOOMD script with **only** `CavityForce` and no LJ/Ewald/thermostats is **not** representative. It can show ~1000+ ns/day but does not reflect production.

### What TO benchmark

1. **Full physics stack** via `CavityMDSimulation` / `05_advanced_run.py` on the real system (`init-0.gsd`, 500 molecules + L, PPPM, bonds, Bussi, etc.)
2. Compare **median TPS** for `--no-cavity` vs cavity-on with matched settings
3. Run on an **idle GPU** or use **paired synchronized processes** on the same GPU to cancel contention noise
4. For fair force-only comparison, detach heavy writers/updaters after setup **or** disable HDF5/temp trackers in code (they are **on by default** in `CavityMDSimulation` even when `--enable-energy-tracker` is false)

### ns/day conversion (1 fs timestep)

```text
ns/day = TPS * 0.0864
```

where TPS = integration steps per second.

### Expected implementation

- GPU class: `CavityForceComputeGPU` (`src/CavityForceComputeGPU.cc`, `.cu`)
- Python wrapper: `hoomd.cavitymd.forces.CavityForce` → should log `CUDA CavityForceComputeGPU initialized successfully`
- Hot path must **not** contain `hipDeviceSynchronize`, host `ArrayHandle` reads, or per-step `hipMemcpy` in `computeForces`
- Energy getters lazily copy a small device buffer only when queried

### Regression tests before trusting benchmarks

```bash
/scratch/mh7373/miniforge3/envs/hoomd/bin/python -m pytest \
  src/pytest/test_cavity_force_gpu_optimization.py -q
```

---

## Analysis and figure reproduction

After production data exist under `aging_weak_lambda/`:

```bash
cd /scratch/mh7373/projects/cav-hoomd/reproduction
/scratch/mh7373/miniforge3/envs/hoomd/bin/python adapters/build_aging_repro_layout.py --profile aging_weak_lambda
bash run_aging_weak_lambda.sh
```

Figures land in `aging_weak_lambda/derived/repro_figures/`.

See `reproduction/README.md` for panel-to-script mapping.

---

## Agent checklist before changing run behavior

1. **Rebuild plugin** if native code changed (`./build_install.sh`).
2. **Run GPU regression tests** (`test_cavity_force_gpu_optimization.py`).
3. **Confirm GPU implementation** in logs (`CavityForceComputeGPU`, not silent CPU fallback).
4. **Never add `--truncate-gsd`** to trajectory-library jobs.
5. **Distinguish** production (`aging_weak_lambda/`, cavity on, HDF5 observables) from IC library (`aging_weak_lambda_ic_library/`, `--no-cavity`, GSD output).
6. **Check GPU contention** (`nvidia-smi`) before interpreting ns/day numbers.
7. **Use `aging_campaign_status.py`** for replica completion, not guesswork from partial files.

---

## Known performance caveats (production driver)

Even with optimized `CavityForce`, the default `CavityMDSimulation` driver enables:

- HDF5 observable writer (default **0.01 ps** period)
- Comprehensive temperature tracker (default **0.1 ps**)

These dominate wall time compared to the cavity force kernel itself. Production campaign accepts this for science output; **throughput benchmarks** should either disable these in the driver or clear `sim.operations.updaters` after attach when measuring raw integrator speed.

---

## OpenMM reference

OpenMM cavity MD reference code lives at:

```text
third_party_refs/openmm-lm/
```

Useful for comparing force definitions and GPU scheduling patterns, not for running this HOOMD campaign directly.

---

## Quick reference commands

| Task | Command |
|------|---------|
| Rebuild plugin | `./build_install.sh` |
| GPU cavity tests | `python -m pytest src/pytest/test_cavity_force_gpu_optimization.py -q` |
| Submit production campaign | `examples/slurm/submit_aging_campaign.sh` |
| Campaign status | `python examples/slurm/aging_campaign_status.py scan --output-base aging_weak_lambda` |
| IC library (local 8×GPU) | `aging_weak_lambda_ic_library/run_local_ic_gsd_library.sh` |
| Reproduce figures | `reproduction/run_aging_weak_lambda.sh` |

---

## Contact context for future agents

If continuing **CavityForce performance** work: the user expects cavity-on production runs to be within roughly **~25% overhead** of no-cavity on the same hardware when measuring the integrator fairly, and far better than the pre-optimization ~30× slowdown caused by sync/host transfers. Always report **full-stack** numbers alongside any microbenchmark.

If continuing **IC library** work: verify GSD frame counts increase over time; 1 frame at multi-ns elapsed time means truncate or overwrite bug, not flush delay.
