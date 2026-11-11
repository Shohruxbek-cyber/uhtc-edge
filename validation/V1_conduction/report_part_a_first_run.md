<!-- First run of V1 part a (2 Nov 2026), archived from commit a9912a1 in the final review (R5p2-03): its order test failed and was redesigned (Supplementary Section S3). The current report is report_part_a.md. -->
# V1 — Conduction solvers vs analytic solutions (part a: constant properties)

Method and exact solutions: `run_v1_analytic.py` header.

## Case A: semi-infinite slab, constant flux 5 MW/m², t = 2 s

| dx0 (µm) | dt (s) | max error / surface rise | energy residual |
|---|---|---|---|
| 80 | 2e-4 | 4.70e-04 | -1.9e-16 |
| 40 | 2e-4 | 4.14e-04 | 0.0e+00 |
| 20 | 2e-4 | 3.86e-04 | -3.7e-16 |
| 10 | 2e-4 | 3.72e-04 | -3.7e-16 |
| 10 | 0.04 | 2.11e-03 | 1.9e-16 |
| 10 | 0.02 | 8.65e-04 | 1.9e-16 |
| 10 | 0.01 | 2.41e-04 | -1.9e-16 |
| 10 | 0.005 | 7.18e-05 | -1.9e-16 |

Observed order, space: 0.18, 0.10, 0.05; time: 1.29, 1.84, 1.75. Finest errors: 3.7e-04 (space study), 7.2e-05 (time study) → **FAIL**.

## Case B: semi-infinite slab, surface step 300 → 2300 K, t = 2 s

| dx0 (µm) | dt (s) | max error / surface rise | energy residual |
|---|---|---|---|
| 80 | 2e-4 | 6.00e-05 | 2.1e-15 |
| 40 | 2e-4 | 7.07e-05 | -8.0e-15 |
| 20 | 2e-4 | 7.56e-05 | -2.3e-15 |
| 10 | 2e-4 | 7.81e-05 | 1.3e-15 |
| 10 | 0.04 | 2.53e-03 | -8.8e-15 |
| 10 | 0.02 | 1.22e-03 | 1.2e-14 |
| 10 | 0.01 | 5.65e-04 | 8.8e-15 |
| 10 | 0.005 | 2.39e-04 | 5.3e-15 |

Observed order, space: -0.24, -0.10, -0.05; time: 1.05, 1.11, 1.24. Finest errors: 7.8e-05 (space study), 2.4e-04 (time study) → **FAIL**.

## Case C: sphere and cylinder, radius 3 mm, uniform flux 2 MW/m², t = 0.12 s (Fo = 0.222)

Exact surface temperature: sphere 403.87 K, cylinder 382.71 K (series checks: T(a,0) − T0 = 6.1e-03 K and 6.1e-03 K).

| Model | mesh | dt (s) | surface T (K) | error / rise | spread over surface (K) | mean T − energy-balance mean (K) | energy residual |
|---|---|---|---|---|---|---|---|
| 2D axisymmetric (sphere) | 20×10 | 1e-3 | 403.78 | 7.84e-04 | 1.2e+00 | 2.5e-01 | -7.3e-14 |
| 2D axisymmetric (sphere) | 40×20 | 1e-3 | 403.85 | 1.33e-04 | 3.3e-01 | 6.2e-02 | -7.1e-14 |
| 2D axisymmetric (sphere) | 80×40 | 1e-3 | 403.86 | 8.98e-05 | 9.4e-02 | 1.5e-02 | 3.9e-14 |
| 2D axisymmetric (sphere): observed order | | | | 2.56, 0.57 | | | |
| 2D planar (cylinder) | 20×10 | 1e-3 | 382.38 | 4.02e-03 | 7.3e-01 | 1.6e-01 | 4.3e-14 |
| 2D planar (cylinder) | 40×20 | 1e-3 | 382.64 | 8.63e-04 | 1.8e-01 | 4.1e-02 | 1.2e-14 |
| 2D planar (cylinder) | 80×40 | 1e-3 | 382.67 | 4.00e-04 | 4.4e-02 | 1.0e-02 | -8.4e-14 |
| 2D planar (cylinder): observed order | | | | 2.22, 1.11 | | | |
| 1D spherical shell | dx0 40 µm | 1e-3 | 403.88 | 1.57e-04 | — | — | 0.0e+00 |
| 1D spherical shell | dx0 20 µm | 1e-3 | 403.87 | 5.03e-05 | — | — | 0.0e+00 |
| 1D spherical shell | dx0 10 µm | 1e-3 | 403.87 | 8.83e-06 | — | — | 0.0e+00 |

Case C → **PASS** (time step 1 ms; its first-order error is included in all rows, so the observed spatial order levels off once the time error dominates).

Overall V1 part a: **FAIL**.
