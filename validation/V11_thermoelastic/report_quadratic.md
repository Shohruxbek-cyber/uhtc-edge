# V11 with quadratic triangles (amendment 1)

Check 0: stress 1.8e-13, displacement 5.6e-14 (meets 1e-9).

| Mesh | Inner hoop (MPa) | Error | Outer hoop (MPa) | Error |
|---|---|---|---|---|
| 40x8 | -1872.6 | 6.06 % | 787.7 | 1.21 % |
| 80x16 | -1928.9 | 3.23 % | 791.8 | 0.69 % |
| 160x32 | -1959.9 | 1.68 % | 794.4 | 0.37 % |
| 320x64 | -1976.4 | 0.85 % | 795.8 | 0.19 % |

Finest meets 2 %; order 0.96 (meets 0.8).

| r/b | Hoop (MPa) | Exact (MPa) | Error |
|---|---|---|---|
| 0.25 | 487.5 | 487.5 | 0.00 % |
| 0.5 | 150.0 | 150.0 | 0.00 % |
| 0.9 | -858.0 | -858.0 | 0.00 % |

**V11 with quadratic triangles PASSES.**
