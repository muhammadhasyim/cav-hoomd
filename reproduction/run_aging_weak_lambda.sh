#!/usr/bin/env bash
# End-to-end reproduction pipeline for aging_weak_lambda in-repo data.
set -euo pipefail

REPO=/scratch/mh7373/projects/cav-hoomd
REPRO=$REPO/reproduction
PY=${PY:-/scratch/mh7373/miniforge3/envs/hoomd/bin/python}
PROFILE=aging_weak_lambda
STAGED=$REPO/aging_weak_lambda/derived/repro_layout
FIGURES=$REPO/aging_weak_lambda/derived/repro_figures

mkdir -p "$FIGURES"

cd "$REPRO"

echo "=== [1/7] Build staged paper-format layout ==="
"$PY" adapters/build_aging_repro_layout.py --profile "$PROFILE"

echo "=== [2/7] Fig 2b F(k,t) plots (existing masters) ==="
# Masters should be built with --lag-extent-mode max (now the process_fskt_only
# default). To rebuild truncated masters, drop --skip_processing.
"$PY" figure2_fkt/process_fskt_only.py \
  --profile "$PROFILE" \
  --skip_processing \
  --max_time 2500

echo "=== [3/7] Fig 2b F(k,t) multi-panel + Fig 2c measured relaxation-time analysis ==="
"$PY" shared/plot_fkt_analysis.py \
  --profile "$PROFILE" \
  --output_dir "$FIGURES"

echo "=== [3b/7] Fig 2b/c fictive (predicted) relaxation panels ==="
"$PY" figure2_3_relaxation/plot_fictive_three_panel_analysis.py \
  --profile "$PROFILE" \
  --no-latex \
  --relaxation-output "$FIGURES/tn_model_relaxation_analysis.pdf" \
  --material-time-output "$FIGURES/fictive_material_time_evolution.pdf"

echo "=== [4/7] Fig 3 energy/temperature panels (all lambda) ==="
# Prefer --no-latex: plot_figure3 falls back to matplotlib cmr10 + mathtext CM
# (CMU OpenType / latex.fmt are unavailable on this cluster).
for lam in 0 0.01 0.016667 0.023333 0.03; do
  tag=$(echo "$lam" | tr '.' 'p')
  "$PY" figure3_energy/plot_figure3.py \
    --profile "$PROFILE" \
    --coupling "$lam" \
    --no-latex \
    --output "$FIGURES/figure3_lambda_${tag}"
done

echo "=== [5/7] Fig 4 material time + collapse ==="
"$PY" figure4_material_time/run_corrected_analysis.py --profile "$PROFILE"
"$PY" figure4_material_time/plot_material_time_and_collapse.py \
  --profile "$PROFILE" \
  --output "$FIGURES/material_time_and_collapse.pdf"

echo "=== [5b/7] Fig 4 four-panel composite ==="
"$PY" figure4_material_time/plot_figure4.py \
  --profile "$PROFILE" \
  --output "$FIGURES/figure4.pdf"

echo "=== [6/7] CSV exports ==="
"$PY" csv_export/export_figure3_csv.py \
  --profile "$PROFILE" \
  --coupling 0.03 \
  --output "$FIGURES/figure3_lambda0p03.csv"
"$PY" csv_export/export_figure4_csv.py \
  --profile "$PROFILE" \
  --output "$FIGURES/figure4.csv"

echo "=== [7/7] Done ==="
echo "Figures and CSVs written under: $FIGURES"
echo "Staged data layout: $STAGED"
