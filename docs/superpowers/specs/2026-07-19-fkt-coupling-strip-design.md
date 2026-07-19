# F(k,t) Coupling Strip Restyle Design

## Goal

Replace the existing multi-panel `fkt_by_coupling_filtered` figure with a
publication-style 1×N strip that matches the reference screenshot: shared y-axis,
viridis \(t_w\) coloring, single right colorbar, and Computer Modern typography.

## Scope

- Restyle `plot_fkt_by_coupling` in `reproduction/shared/plot_fkt_analysis.py`.
- Keep the output stem `fkt_by_coupling_filtered` (PNG) and also write PDF.
- Add focused unit tests for pure helpers (waiting-time mapping, colorbar ticks).
- Do not change F(k,t) data loading, normalization, or other plot types
  (`plot_fkt_by_ref`, diagnostic, relaxation panels).

## Hard requirements

### Typography (non-negotiable)

Use the same Computer Modern path already used by measured relaxation panels:

1. Prefer true LaTeX Computer Modern via `text.usetex=True` when a working
   LaTeX install is available (reuse the existing probe / `_apply_measured_relaxation_rcparams`).
2. If LaTeX is broken on the cluster, fall back to matplotlib Computer Modern
   mathtext + registered `cmr10` TTF fonts (`mathtext.fontset='cm'`).
3. All math labels must be TeX-style strings, e.g.
   - y: `$\phi_k(t; t_w)$` (or `$t_{\mathrm{w}}$` form consistent with sibling figures)
   - x: `$t - t_w$ (ps)` / `$t - t_{\mathrm{w}}$ (ps)`
   - colorbar: `$t_w$ (ps)`
   - panel titles: `$\lambda = \ldots$`

Never ship this figure with the default sans-serif matplotlib font stack.

### Layout

| Element | Spec |
|---|---|
| Grid | `1 × N`, `sharey=True`, wide landscape |
| Y | `[0, 1]`; label only on leftmost panel |
| X | fixed `[0, 1600]` ps; ticks `0, 400, 800, 1200, 1600` |
| Color | `viridis` by \(t_w = \mathrm{ref\_num} \times 200\) ps |
| Colorbar | one vertical bar on the right; ticks every 400 ps; range covering data (typically 0–2400) |
| Style | inward ticks on all sides; light dashed grid; no per-panel legends; no large bold suptitle |
| Panel ID | small `$\lambda = \ldots$` title above each panel |
| Outputs | `fkt_by_coupling_filtered.png` and `.pdf` |

### Data path

Unchanged: filtered couplings → existing `process_fkt_data` → plot. Only layout
and styling change.

## Non-goals

- Restyling `plot_fkt_by_ref` or diagnostic plots.
- Changing relaxation-time algorithms or master-file staging.
- Renaming pipeline entrypoints or profile wiring.

## Testing

- Unit-test pure helpers for \(t_w\) mapping and colorbar tick generation.
- Regeneration check: run `plot_fkt_analysis.py` against staged repro/preliminary
  data and confirm PNG/PDF exist with the new layout (manual visual check).

## Approval notes

- User chose: replace existing multi-panel (not add a second figure).
- User chose: keep small \(\lambda\) titles.
- User required: Computer Modern / LaTeX fonts.
