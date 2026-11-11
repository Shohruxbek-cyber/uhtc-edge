# V6 — Oxidation sub-model vs furnace scale thicknesses

Result: **FAIL** — 24 of 37 independent points within a factor of 3 (criterion, set before the run: all). Geometric-mean model/measured ratio: independent 2.10, calibration set 1.22, glass thickness 1.41.

- Oxide scale, ZrO2 layer, recession and scale + depleted zone (independent): 15 of 18 within a factor of 3, geometric-mean ratio 1.09. The outliers are batch differences in opposite directions: Opila & Halbig's commercial ZrB2-20SiC at 1327 °C (10 and 100 min) recedes 4-8 times faster than the model (calibrated on Shugart's material), and the NASA-made batch of Opila et al. 2003 grew about 5 times less oxide at 1627 °C than the model (and about 3 times less than the commercial batch at the same condition).
- SiC-depleted layer at 1627-1700 °C (Shugart, digitized; independent of the 1800 °C calibration): 9 of 19 within a factor of 3 (zero values excluded from the ratio). The model has no incubation period, while Shugart reports that depletion needs 10-20 min to start at 1650 °C; it therefore over-predicts the depleted depth at short times. Consequence for the flight results: the SiO2 that the model adds to the glass from the depleted layer arrives too early, which makes glass survival somewhat optimistic; the sensitivity case f_SiO_lost = 1 (no SiO2 from the depleted layer) bounds this (Q4).
- Glass thickness depends on the assumed furnace mass-transfer conductance and is reported, not judged.

**With the incubation option** (D-OX-10: depletion starts after an incubation time taken from the onset brackets of Shugart's Fig. 40, bracket midpoints; not fitted to the thicknesses below, but from the same experiments as the depleted-layer points, so not independent of them): 29 of 37 independent points within a factor of 3 (oxide-type 15 of 18, depleted layer 14 of 19); geometric-mean ratio, depleted layer: 3.00 with the option, 3.93 without. Against the criterion set before the run (all within a factor of 3) V6 still fails with the option.

Method, data mapping and the reason for the tolerance: `run_v6.py` header. Furnace mass-transfer conductance g = 0.01 kg/m²/s (ASSUMPTION); it changes only glass-related quantities.

| Set | Source | SiC vol% | T (°C) | t (min) | Quantity | Measured (µm) | ± | Model (µm) | Model/measured |
|---|---|---|---|---|---|---|---|---|---|
| calibration | Shugart2014 | 30 | 1300 | 100 | ZrO2(+C) layer | 25.1 | 8.7 | 21.5 | 0.86 |
| calibration | Shugart2014 | 30 | 1400 | 100 | ZrO2(+C) layer | 15.3 | 13.3 | 23.1 | 1.51 |
| calibration | Shugart2014 | 30 | 1500 | 100 | ZrO2(+C) layer | 21.9 | 15.6 | 24.6 | 1.12 |
| calibration | Shugart2014 | 30 | 1550 | 100 | ZrO2(+C) layer | 23.2 | 10.7 | 25.3 | 1.09 |
| calibration | Shugart2014 | 30 | 1800 | 5 | SiC-depleted layer | 43 | 10 | 48.6 | 1.13 |
| calibration | Shugart2014 | 30 | 1800 | 10 | SiC-depleted layer | 66 | 2 (reading; error bar ~13) | 68.8 | 1.04 |
| calibration | Shugart2014 | 30 | 1800 | 10 | SiC-depleted layer | 23 | 2 (reading; error bar ~6) | 68.8 | 2.99 |
| calibration | Shugart2014 | 30 | 1800 | 15 | SiC-depleted layer | 74 | 6 | 84.3 | 1.14 |
| calibration | Shugart2014 | 30 | 1800 | 15 | SiC-depleted layer | 66 | 2 (reading; error bar ~13) | 84.3 | 1.28 |
| calibration | Shugart2014 | 30 | 1800 | 20 | SiC-depleted layer | 104 | 8 | 97.3 | 0.94 |
| glass (depends on furnace g) | Shugart2014 | 30 | 1300 | 100 | borosilicate layer | 7.7 | 6.0 | 37.3 | 4.85 |
| glass (depends on furnace g) | Shugart2014 | 30 | 1400 | 100 | borosilicate layer | 22.1 | 15.3 | 26.0 | 1.18 |
| glass (depends on furnace g) | Shugart2014 | 30 | 1500 | 20 | borosilicate layer | 5.5 | 0.8 | 9.1 | 1.65 |
| glass (depends on furnace g) | Shugart2014 | 30 | 1500 | 20 | borosilicate layer | 7.1 | 0.8 | 9.1 | 1.28 |
| glass (depends on furnace g) | Shugart2014 | 30 | 1500 | 100 | borosilicate layer | 27.6 | 20.5 | 16.8 | 0.61 |
| glass (depends on furnace g) | Shugart2014 | 30 | 1550 | 100 | borosilicate layer | 14.9 | 6.9 | 16.2 | 1.09 |
| independent | Levine2003 | 20 | 1627 | 100 | scale + depleted zones (face 1; 13 meas.) | 258 | 85 | 506.5 | 1.96 |
| independent | Levine2003 | 20 | 1627 | 100 | scale + depleted zone (face 2; 12 meas.) | 515 | 253 | 506.5 | 0.98 |
| independent | Levine2003 | 20 | 1627 | 100 | scale + depleted zones (both faces; 25 meas.) | 381 | 224 | 506.5 | 1.33 |
| independent | Opila2001 | 20 | 1327 | 10 | recession (incl. depleted zone; per face assu | 60 | 20 (reading) | 6.9 | 0.12 **out** |
| independent | Opila2001 | 20 | 1327 | 100 | oxide scale (ZrO2 with embedded unoxidized Si | 30 | approx. (no SD given) | 58.5 | 1.95 |
| independent | Opila2001 | 20 | 1327 | 100 | recession (incl. depleted zone; per face assu | 80 | 20 (reading) | 21.9 | 0.27 **out** |
| independent | Opila2001 | 20 | 1627 | 10 | recession (incl. depleted zone; per face assu | 100 | 20 (reading) | 106.8 | 1.07 |
| independent | Opila2001 | 20 | 1627 | 100 | oxide scale (SiO2 + ZrO2/glass) | 150 | approx. (no SD given) | 288.9 | 1.93 |
| independent | Opila2001 | 20 | 1627 | 100 | SiC-depleted layer | 100 | approx. (no SD given) | 217.6 | 2.18 |
| independent | Opila2001 | 20 | 1627 | 100 | recession (incl. depleted zone; per face assu | 260 | 20 (reading) | 337.8 | 1.30 |
| independent | Opila2001 | 20 | 1927 | 10 | recession (per face assumed) | 190 | 20 (reading) | 223.9 | 1.18 |
| independent | Opila2001 | 20 | 1927 | 100 | recession (per face assumed) | 960 | 20 (reading) | 708.2 | 0.74 |
| independent | Opila2003 | 20 | 1627 | 100 | total oxide (SiO2 + ZrO2) at one location | 55 | 3 (reading; single location) | 288.9 | 5.25 **out** |
| independent | Shugart2014 | 30 | 1300 | 100 | total oxide | 32.8 | 6.3 | 58.8 | 1.79 |
| independent | Shugart2014 | 30 | 1400 | 100 | total oxide | 37.5 | 17.3 | 49.1 | 1.31 |
| independent | Shugart2014 | 30 | 1500 | 20 | ZrO2(+C) layer | 13.1 | 1.6 | 11.0 | 0.84 |
| independent | Shugart2014 | 30 | 1500 | 20 | ZrO2(+C) layer | 11 | 2.1 | 11.0 | 1.00 |
| independent | Shugart2014 | 30 | 1500 | 100 | total oxide | 39.6 | 26.4 | 41.4 | 1.04 |
| independent | Shugart2014 | 30 | 1550 | 100 | total oxide | 38 | 9.9 | 41.5 | 1.09 |
| independent | Shugart2014 | 30 | 1627 | 60 | SiC-depleted layer | 4 | 2 | 168.5 | 42.13 **out** |
| independent | Shugart2014 | 30 | 1627 | 70 | SiC-depleted layer | 113 | 2 (reading; error bar ~12) | 182.0 | 1.61 |
| independent | Shugart2014 | 30 | 1650 | 35 | SiC-depleted layer | 17 | 2 (reading; error bar ~9) | 128.7 | 7.57 **out** |
| independent | Shugart2014 | 30 | 1650 | 45 | SiC-depleted layer | 12 | 2 | 145.9 | 12.16 **out** |
| independent | Shugart2014 | 30 | 1650 | 50 | SiC-depleted layer | 45 | 2 (reading; error bar ~6) | 153.8 | 3.42 **out** |
| independent | Shugart2014 | 30 | 1650 | 50 | SiC-depleted layer | 0 | 2 | 153.8 | nan **out** |
| independent | Shugart2014 | 30 | 1650 | 55 | SiC-depleted layer | 71 | 2 (reading; error bar ~11) | 161.4 | 2.27 |
| independent | Shugart2014 | 30 | 1650 | 60 | SiC-depleted layer | 75 | 2 (reading; error bar ~5) | 168.5 | 2.25 |
| independent | Shugart2014 | 30 | 1650 | 70 | SiC-depleted layer | 7 | 2 (reading; error bar ~6) | 182.0 | 26.00 **out** |
| independent | Shugart2014 | 30 | 1700 | 10 | SiC-depleted layer | 0 | 2 | 68.8 | nan **out** |
| independent | Shugart2014 | 30 | 1700 | 15 | SiC-depleted layer | 22 | 2 (reading; error bar ~10) | 84.3 | 3.83 **out** |
| independent | Shugart2014 | 30 | 1700 | 20 | SiC-depleted layer | 31 | 2 (reading; error bar ~5) | 97.3 | 3.14 **out** |
| independent | Shugart2014 | 30 | 1700 | 20 | SiC-depleted layer | 8 | 2 | 97.3 | 12.16 **out** |
| independent | Shugart2014 | 30 | 1700 | 25 | SiC-depleted layer | 30 | 2 (reading; error bar ~9) | 108.8 | 3.63 **out** |
| independent | Shugart2014 | 30 | 1700 | 30 | SiC-depleted layer | 70 | 2 (reading; error bar ~18) | 119.2 | 1.70 |
| independent | Shugart2014 | 30 | 1700 | 30 | SiC-depleted layer | 33 | 2 (reading; error bar ~9) | 119.2 | 3.61 **out** |
| independent | Shugart2014 | 30 | 1700 | 35 | SiC-depleted layer | 86 | 2 (reading; error bar ~7) | 128.7 | 1.50 |
| independent | Shugart2014 | 30 | 1700 | 35 | SiC-depleted layer | 70 | 2 (reading; error bar ~18) | 128.7 | 1.84 |
| independent | Shugart2014 | 30 | 1700 | 40 | SiC-depleted layer | 87 | 2 (reading; error bar ~16) | 137.6 | 1.58 |
| independent | Shugart2014 | 30 | 1700 | 40 | SiC-depleted layer | 62 | 2 (reading; error bar ~19) | 137.6 | 2.22 |

Not used (value known only from an abstract; author decision, 10 Nov 2026): Karlsdottir2007 oxide scale 102 µm at 1700 °C, 15 min (abstract (10.1111/j.1551-2916.2007.01861.x)); Karlsdottir2007 SiC-depleted zone 120 µm at 1700 °C, 15 min (abstract (10.1111/j.1551-2916.2007.01861.x)).
