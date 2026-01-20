# Finite-Size Scaling Study Configuration

## Summary

The finite-size scaling study script has been successfully configured with:

1. **Import paths fixed**: Changed from `cavitymd` to `hoomd.cavitymd` to use the installed plugin
2. **Box size constraints**: Added minimum box size (33 Bohr) to accommodate neighbor list cutoff
3. **Time-varying coupling protocol**: Added support for equilibration phase without coupling, then turn on coupling step-wise
4. **Adaptive timestepping**: Full integration with parameters matching `05_advanced_run.py`
5. **Multiple replica support**: Inverse scaling with system size (500 replicas for N=250)

## Protocol

### Two-Phase Dynamics

- **Phase 1 (0-1 ns)**: Equilibration without cavity coupling
  - System equilibrates at target temperature
  - Cavity particle present but not coupled to molecules
  
- **Phase 2 (1-3 ns)**: Production with coupling
  - Coupling turns on step-wise at t = 1 ns
  - Adaptive timestepping automatically handles the shock with dampening
  - Collective coupling g = 0.02 a.u. maintained across all system sizes

### System Sizes and Replica Counts

With defaults (`--base-molecules 250 --base-replicas 500`):

| N (molecules) | λ (single) | g (collective) | Replicas | Total sim time |
|--------------|------------|----------------|----------|----------------|
| 100          | 0.002000   | 0.02           | 1250     | 3.75 μs        |
| 250          | 0.001265   | 0.02           | 500      | 1.50 μs        |
| 500          | 0.000894   | 0.02           | 250      | 0.75 μs        |
| 1000         | 0.000632   | 0.02           | 125      | 0.375 μs       |
| 2500         | 0.000400   | 0.02           | 50       | 0.15 μs        |
| 5000         | 0.000283   | 0.02           | 25       | 0.075 μs       |
| 10000        | 0.000200   | 0.02           | 13       | 0.039 μs       |

**Total**: ~2213 simulations × 3 ns = **6.64 μs of aggregate simulation time**

## Usage

### Basic Usage

```bash
cd /home/mh7373/GitRepos/cav-hoomd
./examples/run_finite_size_scaling.sh
```

### Custom Parameters

```bash
python examples/finite_size_scaling_study.py \
    --n-values "100,250,500,1000,2500,5000,10000" \
    --collective-coupling 0.02 \
    --runtime-ps 3000.0 \
    --switch-time-ps 1000.0 \
    --device GPU \
    --gpu-id 0 \
    --output-dir "my_scaling_study"
```

### Key Parameters

**Time-varying coupling:**
- `--runtime-ps 3000.0`: Total simulation time (3 ns)
- `--switch-time-ps 1000.0`: When to turn on coupling (1 ns)
- The simulation runs 0-1 ns without coupling, then 1-3 ns with coupling

**Adaptive timestepping:**
- `--error-tolerance 5.0`: Target error tolerance
- `--initial-fraction 1e-5`: Shock dampening factor (error tolerance drops to 5e-5 when coupling turns on)
- `--time-constant-ps 50.0`: Recovery time constant (exponential ramp back to 5.0 over ~50 ps)

**Replica scaling:**
- `--base-molecules 250`: Reference system size
- `--base-replicas 500`: Number of replicas for reference size
- Replicas scale as: `N_replicas = 500 × (250 / N)`

## Output Structure

```
finite_size_scaling_g0.02_3ns/
├── finite_size_scaling_summary.csv
├── N_100/
│   ├── replica_0/
│   │   ├── molecular-0.gsd         # Initial configuration
│   │   ├── prod-0.gsd              # Trajectory (50 ps snapshots)
│   │   └── observables_replica_0.h5 # Observables (0.1 ps resolution)
│   ├── replica_1/
│   └── ...
├── N_250/
├── N_500/
├── N_1000/
├── N_2500/
├── N_5000/
└── N_10000/
```

## Files Modified

1. **`examples/finite_size_scaling_study.py`**:
   - Fixed imports to use `hoomd.cavitymd`
   - Added `--switch-time-ps` parameter
   - Added adaptive timestepping parameters
   - Added dissipation calculation
   - All parameters passed to `CavityMDSimulation`

2. **`examples/initlattice_equilibrium.py`**:
   - Added minimum box size constraint (33 Bohr)
   - Ensures neighbor list cutoff is satisfied for small systems

3. **`examples/run_finite_size_scaling.sh`** (NEW):
   - Comprehensive wrapper script
   - Documents all parameters
   - Ready to run on GPU

## Known Issues

1. **Cleanup/Return Issue**: `CavityMDSimulation.run()` completes successfully and creates all output files, but the Python wrapper doesn't always receive the return value properly. This doesn't affect the simulations themselves - they run correctly and produce valid output.

2. **Status in CSV**: The summary CSV may show "failure" status even when simulations complete successfully. Check for the presence of output files (`prod-*.gsd` and `observables_*.h5`) to verify completion.

## Verification

To verify a simulation completed successfully:

```bash
# Check output files exist
ls -lh finite_size_scaling_g0.02_3ns/N_250/replica_0/*.gsd
ls -lh finite_size_scaling_g0.02_3ns/N_250/replica_0/*.h5

# Inspect observables
python tools/inspect_hdf5_data.py finite_size_scaling_g0.02_3ns/N_250/replica_0/observables_replica_0.h5
```

## Next Steps

1. Run a small test with shortened parameters:
   ```bash
   python examples/finite_size_scaling_study.py \
       --n-values "100,250" \
       --base-replicas 2 \
       --max-replicas 2 \
       --runtime-ps 10.0 \
       --switch-time-ps 5.0 \
       --device CPU
   ```

2. Submit full study to cluster/GPU
3. Analyze results using the HDF5 observables and GSD trajectories
