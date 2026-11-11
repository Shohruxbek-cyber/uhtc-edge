# V8 — arc-jet / plasma tests with measured heat flux (ZrB2-SiC)

Inputs, targets and acceptance bands: `run_v8.py` header and `literature/V8_TARGET.md`. 1D model, no adjusted parameter.

## V8-1 NASA Ames AHF, ZrB2-20SiC flat-face coupon (US Patent 5,750,450)

| H0 (MJ/kg) | bound | passive peak T (°C) [1630 ± 80] | oxide 300 s (µm) [110, ×1.5] | depleted (µm) [40, ×2] | active peak T (°C) [2200] | active: T at 50 / 180 s (°C) | pass (passive) |
|---|---|---|---|---|---|---|---|
| 15 | fc | 2541 | 626 | 47 | 2739 | 2727 / 2638 | no |
| 15 | nc | 1868 | 121 | 45 | 2035 | 1961 / 2023 | no |
| 30 | fc | 2667 | 816 | 48 | 2907 | 2889 / 2761 | no |
| 30 | nc | 1673 | 58 | 41 | 1829 | 1629 / 1826 | no |
| 45 | fc | 2710 | 891 | 48 | 2971 | 2952 / 2808 | no |
| 45 | nc | 1801 | 97 | 44 | 1961 | 1827 / 1954 | no |

V8-1 (passive, quantitative bands): **FAIL** (at least one enthalpy and bound inside all bands). The active case is judged qualitatively in the text.

## V8-2 HIT plasma wind tunnel, ZrB2-SiC cylinders (Hu et al. 2013)

| bound | ZS30 T cond.1 (°C) [1960 ± 80] | ZS30 T cond.2 (°C) [1900 ± 80] | ZS30 oxide at 201/402/600 s (µm) [65/99/125] | oxide+depleted at 201/402/600 s (µm) | cond.2 oxide 600 s (µm) [90] | T ZS10 / ZS15 / ZS30, cond.1 (°C) | pass |
|---|---|---|---|---|---|---|---|
| fc | 2869 | 2672 | 661/1349/1848 | 698/1404/1915 | 1226 | 2882 / 2871 / 2869 | no |
| nc | 1967 | 1926 | 35/156/241 | 53/199/299 | 209 | 1915 / 1968 / 1967 | no |

V8-2: **FAIL** (at least one catalysis bound inside all bands).

## Diagnostic (not a prediction): oxidation driven by the measured surface temperature

Isothermal at the measured (peak or steady) temperature, no heat-up; separates the oxidation model from the heating error.

- V8-1 passive (peak T held), 1630 °C for 300 s: 201 s: oxide 50, depleted 40, 300 s: oxide 60, depleted 49 µm (measured: 110 µm oxide, 40 µm depleted).
- V8-2 ZS30 cond. 1 (steady T), 1960 °C for 600 s: 201 s: oxide 189, depleted 40, 300 s: oxide 229, depleted 49, 402 s: oxide 261, depleted 56, 600 s: oxide 316, depleted 69 µm (measured: 65/99/125 µm at 201/402/600 s).
- V8-2 ZS30 cond. 2 (steady T), 1900 °C for 600 s: 201 s: oxide 155, depleted 40, 300 s: oxide 189, depleted 49, 402 s: oxide 216, depleted 56, 600 s: oxide 263, depleted 69 µm (measured: ≈90 µm at 600 s).

