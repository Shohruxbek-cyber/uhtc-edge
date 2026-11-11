# Changelog

## 2.0.0 (11 Nov 2026)

**Code and data only.** The repository no longer holds the manuscript: its LaTeX source, PDFs, figure copies and build scripts, the
circulated copy, the highlights, and the number audit that linked each number of the manuscript to its source. `reproduce.py`
draws the figures into `figures/out/`, and its `numbers` target is gone.

**Added since 1.2.0.**

- Validation uncertainty of V6–V9 following ASME V&V 20 (`validation/VVU/`, `reproduce.py vvu`, `tests/test_vvu.py`). The plan was
  committed before any run. The verdicts of V6–V9 are unchanged.
- Version 3, three optional modules that are off by default. They are glass flow (`src/uhtc/glassflow.py`), a thermoelastic stress
  screen (`src/uhtc/thermoelastic.py`) and an integral flank-heating model. Their checks V10, V2b, V9b, V11 and V12 are in
  `validation/` with their plans, their runs are in `runs/v3/`, `runs/thermoelastic/` and new stages of `runs/envelope_2d/`, and
  `figures/make_figures_v3.py` draws their figures.
- `runs/journal/`:
  - a comparison of Sutton–Graves with Fay–Riddell (equilibrium and frozen edge) at the 30-minute boundary of the 1 mm tip
    (`heating_check.py`);
  - a local wall-pressure option for every surface node of the 2D model (`runs/common.py`, stage `lpgrid`), with its partial
    1 mm grid and the comparison with the baseline (`lpgrid_partial.py`).

**Unchanged.** With every optional module off, the 79 saved cases of `reproduce.py regression` are reproduced bit-exactly, so the
maps and results of 1.2.0 stand.
