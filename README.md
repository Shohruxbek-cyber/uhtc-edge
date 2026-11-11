# Oxidation and temperature limits of sharp ZrB₂–SiC nose tips in hypersonic flight

Code, data and saved results of the study of that name by Shohruxbek Bozorov
([ORCID 0009-0005-5165-093X](https://orcid.org/0009-0005-5165-093X)). The manuscript is not part of this repository.

Repository: https://github.com/Shohruxbek-cyber/uhtc-edge · Archive: https://doi.org/10.5281/zenodo.24099760

## What this is

A low-order model of the stagnation region of a sharp ZrB₂–SiC nose tip in hypersonic flight. It couples stagnation-point heating with
catalysis bounds, a surface energy balance with the growing oxide scale, conduction along the stagnation line (1D) or in the nose tip
(quasi-2D), and an oxide-scale model with three failure criteria: loss of the protective glass, SiC active oxidation and a temperature
limit. It maps the time to failure over Mach 5–12 and 20–40 km for nose radii of 1–10 mm.

**Read the results as a screening tool, not a prediction.** The furnace check (V6), both arc-jet checks (V7, V8) and the plasmatron check
(V9) fail against criteria set before they were run. With the fully catalytic bound the model's surface is too hot under a ground-test
heat flux. In V9 the model also keeps its glass longer than the test does.

## Quick start

Tested on Windows 11 with Python 3.14.6. It is pure Python and should run elsewhere, but no other platform has been tested.

```
python -m venv .venv
.venv\Scripts\activate            (Linux or macOS: source .venv/bin/activate)
pip install -r requirements.txt
python reproduce.py check
```

`check` downloads the NIST-JANAF tables, runs the unit tests, redraws every figure into `figures/out/`, rebuilds the validation and
gain tables and rechecks the reported numbers against the saved outputs. It takes about 5 minutes.

On Windows set `PYTHONUTF8=1` if you run scripts directly; `reproduce.py` sets it for you. Runs use 4 worker processes; set `UHTC_WORKERS`
to change this.

## One-command reproduction

`python reproduce.py <target>`, or `make <target>` if you have `make`. `python reproduce.py list` shows them all.

| Target | What it does | Time |
|---|---|---|
| `setup` | Downloads the NIST-JANAF tables and verifies their checksums | under 1 min |
| `test` | Unit tests (needs `setup` first: it downloads the NIST-JANAF tables) | under 1 min |
| `fig1` … `fig8`, `figS1` … `figS4`, `figures` | Redraws figures from saved outputs into `figures/out/` (`figures` also draws S5–S7) | seconds to 1 min |
| `table1`, `table2`, `table3` | Prints the inputs varied, the validation outcomes and the 1D to 2D gain | seconds |
| `trace` | Rechecks the reported numbers against saved outputs | under 1 min |
| `sample` | Reruns a fixed random sample of 79 stored cases and compares them | about 25 min |
| `v1` … `v9`, `validation` | Reruns the validation checks | about 13 min |
| `histories` | Temperature and layer histories, X-43A path | about 1 min |
| `sensitivity` | One-at-a-time and Latin-hypercube runs | about 4 min |
| `q5` | 1D against 2D at fixed conditions | about 20 min |
| `maps1d` | 5,040 1D map cases | about 35 min |
| `maps2d` | 2D maps, uncertainty band, variants | about 19 h |
| `estimates` | Glass flow, thermal stress, edge oxygen, flank ratio | about 10 min |
| `vvu` | Validation uncertainty of V6–V9 (ASME V&V 20, `validation/VVU/`) | about 2.5 h |
| `regression` | Version 3: reruns the 79 saved cases with every optional module off and checks bit-exact equality | about 13 min |
| `v3checks` | Version 3: checks V10, V2b, V9b, V11, V12 of the optional modules | about 30 min |
| `v3runs` | Version 3: glass flow on the band columns and the full maps, integral flank and geometry columns, stress, frozen shock | about 9 h |
| `v3analysis` | Version 3: boundary comparisons, glass-flow numbers, Figures S5–S7 | about 1 min |
| `all` | Everything from scratch | about 20 h |

Version 3 adds three optional modules: glass flow (`glass_flow="on"`, `src/uhtc/glassflow.py`), a thermoelastic stress screen
(`src/uhtc/thermoelastic.py`) and an integral flank-heating model (`heat_dist="integral"`). They are off in every reported map, and
`regression` checks that those results are reproduced bit-exactly with them off.

Times are for 4 workers on an 8-thread laptop. Saved outputs are included, and run scripts skip cases already stored in their `.jsonl` file.
To recompute a stage, move that file away first.

## Layout

| Path | Content |
|---|---|
| `src/uhtc/` | The model package |
| `tests/` | Unit tests |
| `validation/` | V1–V9, and V2b, V9b, V10–V12 of the optional modules, one folder and report each, with the plans; `validation/VVU/` the validation uncertainty |
| `runs/` | Run scripts and saved outputs; `runs/analysis/` holds the checks and summaries |
| `figures/make_figures.py`, `figures/make_figures_v3.py` | Draw every figure from saved outputs |
| `audit/` | The glass-flow, thermal-stress, edge-oxygen and flank-ratio estimates |
| `literature/digitized/`, `literature/notes/*.csv` | Digitized data points, each with its source figure |
| `INPUTS.md` | Every input with its source at page, table or figure level |
| `DECISIONS.md` | Every modelling decision with its reason |
| `CHANGELOG.md` | Changes of this version |

## What is not included

- **The manuscript.** Neither its source nor its PDF is part of this repository.
- **Source papers.** No third-party PDFs or text. `INPUTS.md` says where every value comes from.
- **NIST-JANAF tables.** NIST reserves copyright on its Standard Reference Data, so `setup` downloads them from
  [janaf.nist.gov](https://janaf.nist.gov) (doi 10.18434/T42S31) and checks them against the checksums of the tables used here.
- **Reading notes.** They quote their sources, so they stay private. `INPUTS.md` still names some of them.

## Licence

Code: MIT (`LICENSE`). Data and saved outputs created for this study: CC BY 4.0 (`LICENSE-DATA`). Digitized points remain attributed to the
sources named in each file.

## How to cite

See `CITATION.cff`.
