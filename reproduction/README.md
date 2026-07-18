# Paper Figure Reproduction

Self-contained copy of all plotting and post-processing scripts needed to reproduce paper figures and CSV exports. **Scripts and calibration files only** — derived simulation data remain in their original archive locations.

Archive root: `/media/extradrive/Trajectories/cavity_supercooled_archive`

## Figure → script map

| Paper panel | Script(s) in `reproduction/` | Typical output |
|-------------|------------------------------|----------------|
| **Fig 2a** — IR spectrum | `figure2_ir/compute_irspectrum.py`, `plot_average_spectra.py`, `plot_minimal_spectra.py` | `average_spectrum_dir_*.txt`, IR PDFs |
| **Fig 2b** — F(k,t) vs t | `figure2_fkt/process_fskt_only.py`, `shared/plot_fkt_analysis.py`, `plot_individual_fkt_coupling.py` | `fkt_coupling_*_filtered.pdf` |
| **Fig 2b/c** — τ̃_s vs λ / t_w | `figure2_3_relaxation/plot_fictive_three_panel_analysis.py` | `fictive_two_panel_relaxation_analysis.pdf` |
| **Fig 3b/c** — ΔV and fictive/kinetic T | `figure3_energy/plot_figure3.py`, `plot_fictive_temperature_components.py` | `figure3_coupling_*.pdf` |
| **Fig 4a–d** — material time, collapse, τ̃_s | `figure4_material_time/run_corrected_analysis.py`, `plot_material_time_and_collapse.py` | `unified_material_time_collapse.png`, collapse fits |
| **CSV exports** | `csv_export/export_fig2_minimal_csvs.py`, `export_figure3_csv.py`, `export_figure4_csv.py` | `figure2a_irspectrum.csv`, etc. |

## Recommended run order

1. **FSKT masters** (if regenerating from replicas):  
   `python figure2_fkt/process_fskt_only.py --base_dir ../final_production_run/data/derived/2026-01-29`

2. **Relaxation table** (1/e criterion):  
   `python shared/plot_fkt_analysis.py` with `--base_dir` pointing at coupling directories

3. **Material time + collapse**:  
   `python figure4_material_time/run_corrected_analysis.py`  
   `python figure4_material_time/plot_material_time_and_collapse.py`

4. **Figure 3 energy/temperature**:  
   `python figure3_energy/plot_figure3.py --coupling 3e-4`

5. **Relaxation panels**:  
   `python figure2_3_relaxation/plot_fictive_three_panel_analysis.py`

6. **IR** (second_rigorous_check data):  
   `figure2_ir/run_ir_analysis.sh` or `compute_irspectrum.py` → `plot_average_spectra.py`

7. **CSV export**:  
   `python csv_export/export_fig2_minimal_csvs.py --out_dir ./csv_output`

## Data locations (not copied)

| Data | Path |
|------|------|
| Production FSKT, material time, relaxation | `final_production_run/data/derived/2026-01-29/` |
| Time series (energy, fictive T) | `.../time_series_output/coupling_*_*.txt` |
| IR aggregated spectra | `second_rigorous_check/data/derived/2025-10-03/average_spectrum_dir_*.txt` |
| Calibration (also in `reproduction/calibration/`) | `potential_energy_components_vs_temperature.txt`, `potential_energy_vs_T.txt` |

Coupling strengths (ε → λ): 0, 3×10⁻⁴→0.042, 5×10⁻⁴→0.070, 7×10⁻⁴→0.098, 1×10⁻³→0.141 a.u.

## Shared modules

Scripts that import `latex_config_adobe`, `material_time_*`, or `plot_fkt_analysis` bootstrap `reproduction/shared/` via:

```python
REPRO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPRO_ROOT / "shared"))
```

See `shared/repro_bootstrap.py` for a reusable helper.

## Running on aging_weak_lambda (in-repo)

The `aging_weak_lambda/` campaign uses a nested layout (`lambda*/cavity_coupling_*`) and sweeps
**λ = 0, 0.01, 0.016667, 0.023333, 0.03** (not the paper ε sweep). Use profile
`aging_weak_lambda` to stage a paper-compatible layout and run the same scripts.

| Item | Path |
|------|------|
| Profile | `profiles/aging_weak_lambda.yaml` |
| Staging adapter | `adapters/build_aging_repro_layout.py` |
| Staged layout (generated) | `aging_weak_lambda/derived/repro_layout/` |
| Figure outputs (generated) | `aging_weak_lambda/derived/repro_figures/` |

**Scope:** Fig 2b/c, Fig 3, Fig 4 (non-IR). Fig 2a IR is unavailable until IR data are added.

```bash
cd /scratch/mh7373/projects/cav-hoomd/reproduction

# 1. Build symlinks + time_series_output + relaxation table
python adapters/build_aging_repro_layout.py --profile aging_weak_lambda

# 2. Full orchestrator (staging + figures + CSV)
bash run_aging_weak_lambda.sh

# 3. Individual figures with profile flag
python figure3_energy/plot_figure3.py --profile aging_weak_lambda --coupling 0.01
python figure4_material_time/run_corrected_analysis.py --profile aging_weak_lambda
python csv_export/export_figure4_csv.py --profile aging_weak_lambda --output /tmp/figure4.csv
```

Paper-archive behaviour is unchanged when `--profile paper` (default) or when `--profile` is omitted.

## Example commands

```bash
cd /media/extradrive/Trajectories/cavity_supercooled_archive/reproduction

# Export all Fig 2 CSVs
python csv_export/export_fig2_minimal_csvs.py --out_dir /tmp/repro_csv

# Figure 3 (g = 3×10⁻⁴)
python figure3_energy/plot_figure3.py --coupling 3e-4

# Two-panel relaxation figure
python figure2_3_relaxation/plot_fictive_three_panel_analysis.py

# Figure 4 CSV bundle
python csv_export/export_figure4_csv.py --out_dir /tmp/repro_csv
```

## Known gaps

- **`analyze_minimum_fictive_temperatures.py` is missing** from the archive. Precomputed `coupling_*_fictive_time_series.txt` files exist in `time_series_output/` but cannot be regenerated from this archive alone.
- **Raw per-replica energy trackers** (`prod-*_energy_tracker.txt`) for the production coupling sweep are not retained; only derived/averaged products remain.
- **`coupling_*_averaged_potential_energy.txt`** in `time_series_output/` was reconstructed from fictive temperature (Oct 2025), not direct MD averaging.

## File inventory

Every copied file is listed in [`MANIFEST.txt`](MANIFEST.txt) with its original source path.
