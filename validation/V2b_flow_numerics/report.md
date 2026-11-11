# V2b: 2D convergence with glass flow on

Plan: `validation/V2b_plan.md` (committed before this run). R_n = 1 mm, 30 km, fully catalytic, Lees distribution, flow on (LA, shear x1, pressure on). Production mesh and 2 s cap against the refined mesh and 0.5 s cap. Criterion: peak temperature within 0.5 %, same first failure, failure time within 2 %.

| Mach | Mode (prod / fine) | Peak T prod / fine (K) | Difference | Failure time prod / fine (s) | Difference | Film audit | Meets |
|---|---|---|---|---|---|---|---|
| 8.5 | none / none | 1827.9 / 1827.3 | 0.032 % | none |  | 8.1e-15 | yes |
| 9 | none / none | 2007.3 / 2003.2 | 0.203 % | none |  | 4.1e-15 | yes |
| 9.5 | glass / glass | 2264.1 / 2263.9 | 0.009 % | 1388.6 / 1399.6 | 0.79 % | 2.7e-15 | yes |
| 10 | glass / glass | 2369.6 / 2369.0 | 0.025 % | 200.6 / 202.7 | 1.04 % | 1.2e-15 | yes |

**V2b PASSES.**
