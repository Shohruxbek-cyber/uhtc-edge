# V7 — Arc-jet test: ZrB2–15SiC hemisphere (SPES, Savino et al. 2010)

Inputs, targets and tolerances: `run_v7.py` header. The table behind them, with the source location of every value
(`literature/V7_TARGET.md`), is a reading note that quotes its sources, so it is not part of the public archive.

Plateau temperature at the end of each step (K):

| Case | Step 1 | Step 2 | Step 3 | Step 4 | Step 5 | Steps within tolerance | Oxide at stagnation (µm) | Depleted (µm) |
|---|---|---|---|---|---|---|---|---|
| Measured (Fig. 8 inset) | 1573 | 1723 | 1843 | 1923 | 2053 | — | 150 (max) | none |
| Tolerance | ±80 | ±80 | ±60 | ±60 | ±60 | | ×1.5 → 100–225 | |
| 1d-fc | 2050 | 2277 | 2583 | 2780 | 2901 | 0/5 | 913 (out) | 59 |
| 1d-nc | 1824 | 2044 | 2187 | 2358 | 2588 | 0/5 | 384 (out) | 57 |
| 1d-switch | 1824 | 2044 | 2187 | 2358 | 2904 | 0/5 | 863 (out) | 57 |
| 1d-nc-eps0.6 | 1885 | 2206 | 2386 | 2610 | 2733 | 0/5 | 601 (out) | 57 |
| 2d-nc | 1685 | 1751 | 1859 | 2024 | 2164 | 2/5 | 100 (out) | 49 |
| 2d-fc | 1827 | 1982 | 2194 | 2431 | 2701 | 0/5 | 306 (out) | 57 |

Result: **FAIL** for the predictive cases (no parameter adjusted to this test). Cases meeting all temperature and oxide tolerances: none. Best temperature agreement: 2d-nc (2/5 steps).

Interpretation. The 1D slab over-predicts by 250-540 K (non-catalytic) because the specimen is a small blunt body heated mainly near the stagnation point (cosine distribution) that conducts heat to, and radiates from, its cooler flanks; the 2D axisymmetric model represents this and, without any adjusted parameter, comes within +112, +28, +16, +101, +111 K of the five plateaus. The remaining over-prediction at steps 4-5 (about +100 K) is consistent with what the model leaves out: the arc-jet stream is far from equilibrium (step-5 CFD: N2 vibrational temperature 8230 K, atomic O mass fraction 0.061) so part of H0 does not reach a weakly catalytic wall (the authors estimate that only about 1/3 of H0 is 'operative', Di Maso 2009 p. 46), and the stem/holder heat leak is not modelled. The fully catalytic bound is far too hot, as the authors also found.

Diagnostics (not predictions):
- D1: in 1D, a constant factor of 0.53 on the non-catalytic heating matches the five plateaus with an rms error of 104 K (plateaus 1395, 1730, 1876, 2033, 2149 K); the measured plateaus rise more gently with enthalpy than any constant factor can reproduce.
- D2: driven by the measured temperature history, the oxidation model gives 83 µm of oxide (ZrO2 52 µm + glass 31 µm) at the stagnation point vs 150 µm measured (maximum; tolerance 100–225 µm) → **outside**; SiC-depleted layer 47 µm vs none observed on the hemisphere.

