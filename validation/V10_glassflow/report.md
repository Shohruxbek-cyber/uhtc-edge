# V10: verification of the glass-flow module

Plan: `validation/V10_plan.md` (committed before this run). Criterion: every value within 1 % of the exact solution.

| Check | Point | Computed (um) | Exact (um) | Relative error | Within 1 % |
|---|---|---|---|---|---|
| 1 shear | t = 2 t_head | 10.0321 | 10.0000 | 3.21e-03 | yes |
| 1 shear | t = 4 t_head | 5.0285 | 5.0000 | 5.70e-03 | yes |
| 1 shear | t = 8 t_head | 2.5171 | 2.5000 | 6.84e-03 | yes |
| 1b pressure | t = 2 t_head | 14.1348 | 14.1421 | 5.16e-04 | yes |
| 1b pressure | t = 4 t_head | 10.0079 | 10.0000 | 7.88e-04 | yes |
| 1b pressure | t = 8 t_head | 7.0815 | 7.0711 | 1.48e-03 | yes |
| 2 steady, explicit | s = 0.25 L | 22.3194 | 22.3607 | 1.85e-03 | yes |
| 2 steady, explicit | s = 0.5 L | 31.5782 | 31.6228 | 1.41e-03 | yes |
| 2 steady, explicit | s = 1 L | 44.6689 | 44.7214 | 1.17e-03 | yes |
| 2 steady, implicit | s = 0.25 L | 22.3220 | 22.3607 | 1.73e-03 | yes |
| 2 steady, implicit | s = 0.5 L | 31.5808 | 31.6228 | 1.33e-03 | yes |
| 2 steady, implicit | s = 1 L | 44.6714 | 44.7214 | 1.12e-03 | yes |

**V10 PASSES: 12 of 12 values within 1 %.** Largest relative error 6.84e-03. Largest mass-audit error of the runs 1.7e-14.
