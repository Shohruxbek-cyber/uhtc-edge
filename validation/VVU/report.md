# VVU: validation uncertainty of V7–V9 (post hoc; the preset verdicts of Table 2 are unchanged)

Plan: `VVU_plan.md` (committed before any run, amendment 1 before any sample). Method: ASME V&V 20-2009 (R2021). E = S − D in K for temperatures and E = ln(S/D) for thicknesses and times; U_val = 2 u_val. V7 heating ratio r = 1.997. Failed runs: 4.

| Check | Quantity | Bound | S | D | E | u_num | u_BC | u_input | u_D | U_val | Label | Share of U_val² (num/BC/input/D) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V7 | plateau temperature, step 1 | fc | 1827 | 1573 | 254 | 5.273 | 160 | 27 | 20 | 327 | not distinguishable at the test uncertainty | 0/96/3/2 |
| V7 | plateau temperature, step 1 | nc | 1685 | 1573 | 112 | 4.863 | 147 | 36 | 20 | 306 | not distinguishable at the test uncertainty | 0/93/5/2 |
| V7 | plateau temperature, step 2 | fc | 1982 | 1723 | 259 | 5.720 | 193 | 39 | 22 | 397 | not distinguishable at the test uncertainty | 0/95/4/1 |
| V7 | plateau temperature, step 2 | nc | 1751 | 1723 | 28 | 5.055 | 164 | 36 | 22 | 340 | not distinguishable at the test uncertainty | 0/94/5/2 |
| V7 | plateau temperature, step 3 | fc | 2194 | 1843 | 351 | 6.335 | 250 | 50 | 24 | 513 | not distinguishable at the test uncertainty | 0/95/4/1 |
| V7 | plateau temperature, step 3 | nc | 1859 | 1843 | 16 | 5.367 | 191 | 34 | 24 | 391 | not distinguishable at the test uncertainty | 0/95/3/2 |
| V7 | plateau temperature, step 4 | fc | 2431 | 1923 | 508 | 7.019 | 295 | 83 | 25 | 614 | not distinguishable at the test uncertainty | 0/92/7/1 |
| V7 | plateau temperature, step 4 | nc | 2024 | 1923 | 101 | 5.842 | 229 | 46 | 25 | 470 | not distinguishable at the test uncertainty | 0/95/4/1 |
| V7 | plateau temperature, step 5 | fc | 2701 | 2053 | 648 | 7.797 | 308 | 97 | 27 | 649 | not distinguishable at the test uncertainty | 0/90/9/1 |
| V7 | plateau temperature, step 5 | nc | 2164 | 2053 | 111 | 6.246 | 297 | 56 | 27 | 607 | not distinguishable at the test uncertainty | 0/96/3/1 |
| V7 | oxide at the stagnation point | fc | 306 | 150 | 0.712 | 0.012 | 0.694 | 0.559 | 0.000* | 1.783 | not distinguishable at the test uncertainty | 0/61/39/0 |
| V7 | oxide at the stagnation point | nc | 100 | 150 | -0.407 | 0.012 | 1.192 | 0.498 | 0.000* | 2.583 | not distinguishable at the test uncertainty | 0/85/15/0 |
| V7 | SiC-depleted layer (none observed) | fc | no | measured | | | | | | | outcome not reproduced within input uncertainty (n = 97) | |
| V7 | SiC-depleted layer (none observed) | nc | no | measured | | | | | | | outcome reproduced in 11 % of the joint sample (n = 100) | |
| V8_1 | peak temperature (30 MJ/kg nominal) | fc | 2940 | 1903 | 1037 | 1.358 | 127 | 160 | 0.000* | 408 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/39/61/0 |
| V8_1 | peak temperature (30 MJ/kg nominal) | nc | 1946 | 1903 | 43 | 0.899 | 130 | 44 | 0.000* | 274 | not distinguishable at the test uncertainty | 0/90/10/0 |
| V8_1 | oxide at 300 s | fc | 816 | 110 | 2.003 | 0.000 | 0.242 | 0.621 | 0.000* | 1.332 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/13/87/0 |
| V8_1 | oxide at 300 s | nc | 58 | 110 | -0.633 | 0.000 | 0.362 | 0.422 | 0.000* | 1.111 | not distinguishable at the test uncertainty | 0/42/58/0 |
| V8_1 | SiC-depleted layer at 300 s | fc | 48 | 40 | 0.173 | 0.000 | 0.003 | 0.203 | 0.000* | 0.406 | not distinguishable at the test uncertainty | 0/0/100/0 |
| V8_1 | SiC-depleted layer at 300 s | nc | 41 | 40 | 0.024 | 0.000 | 0.042 | 0.194 | 0.000* | 0.397 | not distinguishable at the test uncertainty | 0/4/96/0 |
| V8_2 | steady temperature, condition 1 | fc | 3143 | 2233 | 909 | 1.452 | 53 | 116 | 0.000* | 255 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/18/82/0 |
| V8_2 | steady temperature, condition 1 | nc | 2241 | 2233 | 7.427 | 1.035 | 43 | 50 | 0.000* | 131 | not distinguishable at the test uncertainty | 0/42/57/0 |
| V8_2 | steady temperature, condition 2 | fc | 2945 | 2173 | 772 | 1.360 | 47 | 105 | 0.000* | 231 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/17/83/0 |
| V8_2 | steady temperature, condition 2 | nc | 2199 | 2173 | 26 | 1.016 | 40 | 48 | 0.000* | 126 | not distinguishable at the test uncertainty | 0/41/59/0 |
| V8_2 | reacted layer (oxide + depleted), condition 1, 201 s | fc | 698 | 65 | 2.374 | 0.000 | 0.101 | 0.413 | 0.013* | 0.850 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/6/94/0 |
| V8_2 | reacted layer (oxide + depleted), condition 1, 201 s | nc | 53 | 65 | -0.210 | 0.000 | 0.231 | 0.243 | 0.013* | 0.671 | not distinguishable at the test uncertainty | 0/47/52/0 |
| V8_2 | reacted layer (oxide only), condition 1, 201 s | fc | 661 | 65 | 2.319 | 0.000 | 0.105 | 0.446 | 0.013* | 0.916 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/5/95/0 |
| V8_2 | reacted layer (oxide only), condition 1, 201 s | nc | 35 | 65 | -0.614 | 0.000 | 0.274 | 0.376 | 0.013* | 0.931 | not distinguishable at the test uncertainty | 0/35/65/0 |
| V8_2 | reacted layer (oxide + depleted), condition 1, 402 s | fc | 1404 | 99 | 2.652 | 0.000 | 0.091 | 0.445 | 0.009* | 0.909 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/4/96/0 |
| V8_2 | reacted layer (oxide + depleted), condition 1, 402 s | nc | 199 | 99 | 0.699 | 0.000 | 0.122 | 0.332 | 0.009* | 0.707 | not distinguishable at the test uncertainty | 0/12/88/0 |
| V8_2 | reacted layer (oxide only), condition 1, 402 s | fc | 1349 | 99 | 2.612 | 0.000 | 0.094 | 0.471 | 0.009* | 0.961 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/4/96/0 |
| V8_2 | reacted layer (oxide only), condition 1, 402 s | nc | 156 | 99 | 0.452 | 0.000 | 0.144 | 0.458 | 0.009* | 0.961 | not distinguishable at the test uncertainty | 0/9/91/0 |
| V8_2 | reacted layer (oxide + depleted), condition 1, 600 s | fc | 1915 | 125 | 2.729 | 0.000 | 0.091 | 0.462 | 0.007* | 0.943 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/4/96/0 |
| V8_2 | reacted layer (oxide + depleted), condition 1, 600 s | nc | 299 | 125 | 0.873 | 0.000 | 0.088 | 0.362 | 0.007* | 0.745 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/6/94/0 |
| V8_2 | reacted layer (oxide only), condition 1, 600 s | fc | 1848 | 125 | 2.693 | 0.000 | 0.093 | 0.487 | 0.007* | 0.991 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/4/96/0 |
| V8_2 | reacted layer (oxide only), condition 1, 600 s | nc | 241 | 125 | 0.655 | 0.000 | 0.104 | 0.493 | 0.007* | 1.008 | not distinguishable at the test uncertainty | 0/4/96/0 |
| V8_2 | oxide, condition 2, 600 s | fc | 1226 | 90 | 2.612 | 0.000 | 0.102 | 0.505 | 0.000* | 1.031 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/4/96/0 |
| V8_2 | oxide, condition 2, 600 s | nc | 209 | 90 | 0.841 | 0.000 | 0.122 | 0.471 | 0.000* | 0.973 | not distinguishable at the test uncertainty | 0/6/94/0 |
| V8_2 | trend T(ZS10) > T(ZS15) > T(ZS30) | fc | yes | measured | | | | | | | outcome reproduced in 97 % of the joint sample (n = 100) | |
| V8_2 | trend T(ZS10) > T(ZS15) > T(ZS30) | nc | no | measured | | | | | | | outcome not reproduced within input uncertainty (n = 100) | |
| V9 | sample 1.8: low-stage plateau | fc | 2067 | 1895 | 172 | 0.955 | 27 | 45 | 10 | 106 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/26/71/4 |
| V9 | sample 1.8: low-stage plateau | nc | 1733 | 1895 | -162 | 0.800 | 23 | 37 | 10 | 90 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/26/69/5 |
| V9 | sample 1.8: high-stage temperature at the reference time | fc | 2387 | 2110 | 277 | 1.103 | 70 | 67 | 10 | 195 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/52/47/1 |
| V9 | sample 1.8: high-stage temperature at the reference time | nc | 1872 | 2110 | -238 | 0.865 | 37 | 43 | 10 | 115 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/41/55/3 |
| V9 | sample 1.8: jump outcome (no glass loss in the high stage) | fc | yes | measured | | | | | | | outcome reproduced in 58 % of the joint sample (n = 100) | |
| V9 | sample 1.8: jump outcome (no glass loss in the high stage) | nc | yes | measured | | | | | | | outcome reproduced in 100 % of the joint sample (n = 100) | |
| V9 | sample 1.5: low-stage plateau | fc | 2067 | 1825 | 242 | 0.955 | 27 | 45 | 10 | 106 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/26/71/4 |
| V9 | sample 1.5: low-stage plateau | nc | 1733 | 1825 | -92 | 0.800 | 23 | 37 | 10 | 90 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/26/69/5 |
| V9 | sample 1.5: high-stage temperature at the reference time | fc | 2444 | 2210 | 234 | 1.129 | 77 | 67 | 10 | 205 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/56/43/1 |
| V9 | sample 1.5: high-stage temperature at the reference time | nc | 1890 | 2210 | -320 | 0.873 | 41 | 43 | 10 | 121 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/46/51/3 |
| V9 | sample 1.5: glass loss after the switch | fc | — | 68 | > 0.593 | 0.003 | — | — | 0.026 | — | censored at the baseline (E > 0.59) and in 61 of 100 group samples; the rule of plan Section 3 cannot be applied | — |
| V9 | sample 1.5: glass loss after the switch | nc | — | 68 | > 0.593 | 0.003 | — | — | 0.026 | — | no glass loss by the end of the run in any sample (censored, E > 0.59); see the outcome row | — |
| V9 | sample 1.5: jump outcome (glass loss in the high stage) | fc | no | measured | | | | | | | outcome reproduced in 78 % of the joint sample (n = 100) | |
| V9 | sample 1.5: jump outcome (glass loss in the high stage) | nc | no | measured | | | | | | | outcome not reproduced within input uncertainty (n = 100) | |
| V9 | sample 1.6: low-stage plateau | fc | 2067 | 1860 | 207 | 0.955 | 27 | 45 | 10 | 106 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/26/71/4 |
| V9 | sample 1.6: low-stage plateau | nc | 1733 | 1860 | -127 | 0.800 | 23 | 37 | 10 | 90 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/26/69/5 |
| V9 | sample 1.6: high-stage temperature at the reference time | fc | 2475 | 2250 | 225 | 1.143 | 42 | 44 | 10 | 123 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/47/51/3 |
| V9 | sample 1.6: high-stage temperature at the reference time | nc | 1852 | 2250 | -398 | 0.855 | 40 | 40 | 10 | 115 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/49/48/3 |
| V9 | sample 1.6: glass loss after the switch | fc | 71 | 21 | 1.215 | 0.003 | 0.344 | 0.445 | 0.084 | 1.138 | model error indicated at the test uncertainty (stated uncertainties only: U_val is a lower bound) | 0/37/61/2 |
| V9 | sample 1.6: glass loss after the switch | nc | — | 21 | > 1.751 | 0.003 | — | — | 0.084 | — | no glass loss by the end of the run in any sample (censored, E > 1.75); see the outcome row | — |
| V9 | sample 1.6: jump outcome (glass loss in the high stage) | fc | yes | measured | | | | | | | outcome reproduced in 99 % of the joint sample (n = 100) | |
| V9 | sample 1.6: jump outcome (glass loss in the high stage) | nc | no | measured | | | | | | | outcome not reproduced within input uncertainty (n = 100) | |

\* u_D holds only what the source states. For every check u_BC holds only what the sources state (plan 5.2), so U_val is a lower bound throughout.

## Calorimeter catalysis (plan 5.2): S with a copper efficiency of 0.01 and 0.067, other inputs at baseline

| Check | Quantity | Bound | S (fully catalytic calorimeter) | γ_Cu 0.01 | γ_Cu 0.067 |
|---|---|---|---|---|---|
| V8_1 | peak temperature (30 MJ/kg nominal) | fc | 2940 | 3400 | 3041 |
| V8_1 | peak temperature (30 MJ/kg nominal) | nc | 1946 | 2235 | 2010 |
| V8_1 | oxide at 300 s | fc | 816 | 1907 | 1007 |
| V8_1 | oxide at 300 s | nc | 58 | 150 | 76 |
| V8_1 | SiC-depleted layer at 300 s | fc | 48 | 48 | 48 |
| V8_1 | SiC-depleted layer at 300 s | nc | 41 | 46 | 43 |
| V8_2 | steady temperature, condition 1 | fc | 3143 | 3342 | 3180 |
| V8_2 | steady temperature, condition 1 | nc | 2241 | 2400 | 2270 |
| V8_2 | steady temperature, condition 2 | fc | 2945 | 3119 | 2978 |
| V8_2 | steady temperature, condition 2 | nc | 2199 | 2348 | 2228 |
| V8_2 | reacted layer (oxide + depleted), condition 1, 201 s | fc | 698 | 1008 | 754 |
| V8_2 | reacted layer (oxide + depleted), condition 1, 201 s | nc | 53 | 122 | 64 |
| V8_2 | reacted layer (oxide only), condition 1, 201 s | fc | 661 | 970 | 717 |
| V8_2 | reacted layer (oxide only), condition 1, 201 s | nc | 35 | 95 | 44 |
| V8_2 | reacted layer (oxide + depleted), condition 1, 402 s | fc | 1404 | 1955 | 1501 |
| V8_2 | reacted layer (oxide + depleted), condition 1, 402 s | nc | 199 | 302 | 220 |
| V8_2 | reacted layer (oxide only), condition 1, 402 s | fc | 1349 | 1900 | 1446 |
| V8_2 | reacted layer (oxide only), condition 1, 402 s | nc | 156 | 254 | 175 |
| V8_2 | reacted layer (oxide + depleted), condition 1, 600 s | fc | 1915 | 2670 | 2047 |
| V8_2 | reacted layer (oxide + depleted), condition 1, 600 s | nc | 299 | 388 | 325 |
| V8_2 | reacted layer (oxide only), condition 1, 600 s | fc | 1848 | 2602 | 1980 |
| V8_2 | reacted layer (oxide only), condition 1, 600 s | nc | 241 | 326 | 265 |
| V8_2 | oxide, condition 2, 600 s | fc | 1226 | 1776 | 1323 |
| V8_2 | oxide, condition 2, 600 s | nc | 209 | 316 | 231 |
| V9 | sample 1.8: low-stage plateau | fc | 2067 | 2167 | 2086 |
| V9 | sample 1.8: low-stage plateau | nc | 1733 | 1818 | 1749 |
| V9 | sample 1.8: high-stage temperature at the reference time | fc | 2387 | 2665 | 2417 |
| V9 | sample 1.8: high-stage temperature at the reference time | nc | 1872 | 2007 | 1897 |
| V9 | sample 1.5: low-stage plateau | fc | 2067 | 2167 | 2086 |
| V9 | sample 1.5: low-stage plateau | nc | 1733 | 1818 | 1749 |
| V9 | sample 1.5: high-stage temperature at the reference time | fc | 2444 | 2742 | 2476 |
| V9 | sample 1.5: high-stage temperature at the reference time | nc | 1890 | 2040 | 1924 |
| V9 | sample 1.5: glass loss after the switch | fc | — | 39 | — |
| V9 | sample 1.5: glass loss after the switch | nc | — | — | — |
| V9 | sample 1.6: low-stage plateau | fc | 2067 | 2167 | 2086 |
| V9 | sample 1.6: low-stage plateau | nc | 1733 | 1818 | 1749 |
| V9 | sample 1.6: high-stage temperature at the reference time | fc | 2475 | 2637 | 2511 |
| V9 | sample 1.6: high-stage temperature at the reference time | nc | 1852 | 2003 | 1879 |
| V9 | sample 1.6: glass loss after the switch | fc | 71 | 21 | 51 |
| V9 | sample 1.6: glass loss after the switch | nc | — | — | — |
