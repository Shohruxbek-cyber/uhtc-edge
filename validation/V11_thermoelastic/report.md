# V11: verification of the thermoelastic module

Plan: `validation/V11_plan.md` (committed before this run).

**Check 0, patch test:** largest element stress 1.5e-14 E alpha dT, displacement error 6.6e-15: meets 1e-9.

**Check 1, thick-walled sphere:** exact hoop stress -1993.4 MPa (inner) and 797.3 MPa (outer).

| Mesh | Inner hoop (MPa) | Error | Outer hoop (MPa) | Error |
|---|---|---|---|---|
| 40x8 | -2100.5 | 5.38 % | 1099.3 | 37.87 % |
| 80x16 | -2060.4 | 3.37 % | 797.8 | 0.06 % |
| 160x32 | -2030.6 | 1.87 % | 797.1 | 0.02 % |
| 320x64 | -2012.9 | 0.98 % | 797.1 | 0.03 % |

Finest mesh meets 2 %. Observed order of the inner-surface error over the last three meshes 0.89 (pairs 0.85, 0.93): meets 0.8.

**Check 2, free thin disc (h/b = 0.02, plane-stress solution):**

| r/b | Hoop (MPa) | Exact (MPa) | Error |
|---|---|---|---|
| 0.25 | 487.5 | 487.5 | 0.00 % |
| 0.5 | 150.0 | 150.0 | 0.00 % |
| 0.9 | -858.0 | -858.0 | 0.00 % |

**V11 PASSES.**
