# V3 — Radiative-equilibrium stagnation temperature vs Kolodziej (1997) Table 1

Case: 2D unswept leading edge, R_n = 0.762 mm, emissivity 0.80, limit 2149 °C, fully catalytic hot wall,
radiative equilibrium (no conduction). Velocity at which T_w reaches the limit:

| Altitude | Kolodziej | Hill | Model (S–G, 2D) | Model/Kolodziej | Fay–Riddell 2D | S–G with q ×0.85 / ×1.15 |
|---|---|---|---|---|---|---|
| 0 km | 2.30 km/s | 2.25 km/s | 2.296 km/s | 0.998 | 2.298 km/s | 2.301 / 2.292 km/s |
| 30 km | 2.53 km/s | 2.56 km/s | 2.523 km/s | 0.997 | 2.537 km/s | 2.556 / 2.498 km/s |

Digitized points of Kolodziej's Fig. 1 (higher altitudes, heating not blocked; ±0.05 km/s reading
uncertainty; `literature/digitized/Kolodziej1997_fig1.csv`):

| Wall | Altitude | Kolodziej Fig. 1 | Model | Model/Kolodziej |
|---|---|---|---|---|
| FCHW | 50 km | 3.04 km/s | 3.018 km/s | 0.993 |
| FCHW | 60 km | 3.43 km/s | 3.401 km/s | 0.992 |
| FCHW | 65 km | 3.68 km/s | 3.657 km/s | 0.994 |
| FCHW | 70 km | 4.01 km/s | 3.973 km/s | 0.991 |
| NCHW | 50 km | 3.23 km/s | 3.246 km/s | 1.005 |
| NCHW | 55 km | 3.55 km/s | 3.564 km/s | 1.004 |
| NCHW | 60 km | 3.82 km/s | 3.901 km/s | 1.021 |
| NCHW | 65 km | 4.12 km/s | 4.232 km/s | 1.027 |
| NCHW | 70 km | 4.56 km/s | 4.790 km/s | 1.051 |

Result (primary check, tabulated Table 1 values): **PASS** — largest |V_model/V_Kolodziej − 1| = 0.003 (criterion ≤ 0.05, set before the run).
Digitized Fig. 1 points: largest deviation 0.051; 1 of 9 outside ±0.05 (flagged; see note below).

The last column shows the velocity at the limit when the heating is changed by ±15 % (the spread of
stagnation-heating methods quoted by Tauber 1989), to show the size of the tolerance in velocity terms.

**Flagged point.** The non-catalytic (NCHW) limit at 70 km lies 5 % above Kolodziej's curve, the model
giving the higher velocity (lower non-catalytic heating). The deviation grows with altitude and enthalpy, where
h_D/H_s is largest, so it points to the treatment of the dissociation enthalpy in the non-catalytic bound
(Kolodziej's curve fits for the edge composition vs our JANAF equilibrium; frozen vs equilibrium wall enthalpy).
The point is outside the conditions of this study (z ≤ 40 km, V ≤ 3.8 km/s); at the digitized conditions
nearest to the study's (50–55 km, 3.2–3.6 km/s) the bound agrees within 0.5 %.
UNSURE: which of the two treatments is closer to reality at 70 km.

**Limitation of this check.** At both tabulated points the limit temperature (2422 K) is close to the
stagnation temperature of the flow, so the heating is nearly 'blocked' (H_s − h_w small) and the limiting
velocity depends mainly on the stagnation thermodynamics (H_s, equilibrium h_w) rather than on the heat-transfer
coefficient: a ±15 % change of the heating moves it by ≤ 1.3 %. The heat-transfer coefficient itself is checked
in V4 (Sutton–Graves vs Fay–Riddell) and, with the material response, in V7 (arc jet).
