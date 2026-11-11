# V9b: the plasmatron check V9 with glass flow at the stagnation point

Plan: `validation/V9b_plan.md` (committed before this run). V9 model and inputs, plus the stagnation-point film-flow sink (LA viscosity, beta from the Fay-Riddell inversion of the measured flux). V9 values in brackets.

| Bound | Sample | Low plateau (K) | T at reference time (K) | Glass loss after switch (s) | Largest sink rate (um/s) | Film removed by flow (um) | Outcome | Timing | Temperature |
|---|---|---|---|---|---|---|---|---|---|
| fc | 1.8 | 2066.7 [2066.7] | 2387.4 [2387.4] | none by 480 s [none by 480 s] | 6.95e-04 | 5.86e-02 | yes | yes | no |
| fc | 1.5 | 2066.7 [2066.7] | 2444.4 [2444.4] | none by 480 s [none by 480 s] | 1.09e-03 | 9.84e-02 | no | no | no |
| fc | 1.6 | 2066.6 [2066.6] | 2475.0 [2475.0] | 70.8 [70.8] | 1.82e-03 | 5.88e-02 | yes | no | no |
| nc | 1.8 | 1733.0 [1733.0] | 1871.8 [1871.8] | none by 480 s [none by 480 s] | 3.07e-10 | 2.33e-08 | yes | yes | no |
| nc | 1.5 | 1733.0 [1733.0] | 1890.3 [1890.3] | none by 480 s [none by 480 s] | 6.00e-10 | 4.42e-08 | no | no | no |
| nc | 1.6 | 1733.0 [1733.0] | 1851.7 [1851.7] | none by 480 s [none by 480 s] | 6.47e-08 | 2.53e-06 | no | no | no |

**V9b FAILS** (fully catalytic misses, non-catalytic misses the criterion in all three samples).

Sensitivity, fully catalytic, not scored (glass loss after the switch in s; film removed by flow in um):

- beta_x0.5_1.8: none by 480 s; removed 2.86e-02
- beta_x0.5_1.5: none by 480 s; removed 4.82e-02
- beta_x0.5_1.6: 70.8; removed 2.90e-02
- beta_x2_1.8: none by 480 s; removed 1.23e-01
- beta_x2_1.5: none by 480 s; removed 2.05e-01
- beta_x2_1.6: 70.7; removed 1.21e-01
- F320_1.8: none by 480 s; removed 9.55e-02
- F320_1.5: none by 480 s; removed 1.45e-01
- F320_1.6: 70.8; removed 7.70e-02
