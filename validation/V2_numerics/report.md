# V2: convergence of the coupled 1D model on production cases

Settings and criterion: `run_v2.py` header. Each case runs 600 s (no stop at failure).

## Mach 8.5, 30 km, R_n = 3 mm, fc

| Settings | T_s(10 s) | T_s(60 s) | T_s(300 s) | first failure (s, mode) | oxide at 600 s (µm) | steps | CPU (s) |
|---|---|---|---|---|---|---|---|
| baseline | 1380.6 | 2262.0 | 2368.5 | nan (none) | 324.80 | 675 | 0.6 |
| fine mesh | 1380.4 | 2262.1 | 2368.5 | nan (none) | 324.81 | 678 | 0.5 |
| fine time | 1381.2 | 2263.6 | 2368.5 | nan (none) | 324.83 | 1962 | 2.1 |
| finest | 1381.1 | 2263.6 | 2368.5 | nan (none) | 324.83 | 1962 | 3.3 |
| deeper slab | 1316.8 | 1938.0 | 2335.3 | nan (none) | 284.65 | 650 | 1.1 |

Baseline vs finest: temperatures 0.07 %, failure time 0.00 %, oxide 0.01 % → **PASS**.

## Mach 10.0, 30 km, R_n = 3 mm, nc

| Settings | T_s(10 s) | T_s(60 s) | T_s(300 s) | first failure (s, mode) | oxide at 600 s (µm) | steps | CPU (s) |
|---|---|---|---|---|---|---|---|
| baseline | 1683.9 | 2769.3 | 2875.9 | 30.33 (temperature) | 1174.21 | 743 | 1.1 |
| fine mesh | 1683.7 | 2769.3 | 2875.9 | 30.33 (temperature) | 1174.23 | 743 | 1.2 |
| fine time | 1684.5 | 2771.6 | 2875.9 | 30.22 (temperature) | 1174.51 | 2231 | 3.2 |
| finest | 1684.4 | 2771.6 | 2875.9 | 30.22 (temperature) | 1174.53 | 2235 | 3.2 |
| deeper slab | 1609.8 | 2376.8 | 2838.7 | 67.78 (temperature) | 878.06 | 720 | 1.1 |

Baseline vs finest: temperatures 0.08 %, failure time 0.36 %, oxide 0.03 % → **PASS**.

## Mach 9.0, 25 km, R_n = 1 mm, fc

| Settings | T_s(10 s) | T_s(60 s) | T_s(300 s) | first failure (s, mode) | oxide at 600 s (µm) | steps | CPU (s) |
|---|---|---|---|---|---|---|---|
| baseline | 2232.8 | 2739.1 | 2756.7 | 15.84 (temperature) | 897.61 | 749 | 1.0 |
| fine mesh | 2232.7 | 2739.2 | 2756.7 | 15.84 (temperature) | 897.63 | 749 | 1.1 |
| fine time | 2234.1 | 2740.0 | 2756.7 | 15.77 (temperature) | 897.81 | 2242 | 3.2 |
| finest | 2234.0 | 2740.0 | 2756.7 | 15.77 (temperature) | 897.82 | 2238 | 3.0 |
| deeper slab | 2162.6 | 2545.4 | 2742.9 | 28.43 (temperature) | 731.95 | 733 | 1.2 |

Baseline vs finest: temperatures 0.05 %, failure time 0.44 %, oxide 0.02 % → **PASS**.

Overall V2: **PASS**.
