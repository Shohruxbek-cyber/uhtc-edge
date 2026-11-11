# Final-review analyses (runs/final_review/analyses.py)

## Rarefaction (Kolodziej Eq. 7) and Knudsen number at the highest passing 2D conditions
- 1.0 mm 20-22 km: Re_t2 [np.float64(11939.622093020751), np.float64(15874.7373433586)], q_rare/q_cont [1.0478594833537864, 1.0596845063548466], Kn [0.0008035872386256937, 0.0011108355700331085]
- 1.0 mm 30-30 km: Re_t2 [np.float64(3823.209832239455), np.float64(3823.209832239455)], q_rare/q_cont [1.119794293623673, 1.119794293623673], Kn [0.003938707968186385, 0.003938707968186385]
- 1.0 mm 35-40 km: Re_t2 [np.float64(1098.6583267300919), np.float64(2006.262561071123)], q_rare/q_cont [1.1552777244099608, 1.181754813161593], Kn [0.008688198965865954, 0.01873222032676369]
- 3.0 mm 20-22 km: Re_t2 [np.float64(36894.8782015236), np.float64(50332.594776861006)], q_rare/q_cont [1.0202116303877942, 1.0243429912789006], Kn [0.00026786241287523124, 0.0003702785233443695]
- 3.0 mm 30-30 km: Re_t2 [np.float64(12318.249291646016), np.float64(12318.249291646016)], q_rare/q_cont [1.0583035402839587, 1.0583035402839587], Kn [0.0013129026560621281, 0.0013129026560621281]
- 3.0 mm 35-40 km: Re_t2 [np.float64(3521.3131049852987), np.float64(6470.067905993342)], q_rare/q_cont [1.0902577758254686, 1.1244797401758335], Kn [0.002896066321955318, 0.006244073442254563]
- 10.0 mm 20-22 km: Re_t2 [np.float64(130761.49248486139), np.float64(172713.10886356796)], q_rare/q_cont [1.0, 1.0], Kn [8.035872386256936e-05, 0.00011108355700331085]
- 10.0 mm 30-30 km: Re_t2 [np.float64(43997.70431444805), np.float64(43997.70431444805)], q_rare/q_cont [1.0217119411177538, 1.0217119411177538], Kn [0.00039387079681863846, 0.00039387079681863846]
- 10.0 mm 35-40 km: Re_t2 [np.float64(12734.922323972502), np.float64(23054.098388270693)], q_rare/q_cont [1.035236541111272, 1.0568542134837653], Kn [0.0008688198965865954, 0.001873222032676369]

## Laminar flow: free-stream Re over the wetted length (cap and 30 mm flank), 1, 3 and 10 mm tips
- 1 mm, 20 km, Mach 8.5 (2D boundary): Re_x = 4.94e+05
- 1 mm, 20 km, Mach 10.0: Re_x = 5.81e+05
- 1 mm, 20 km, Mach 12.0: Re_x = 6.97e+05
- 1 mm, 30 km, Mach 8.5: Re_x = 1.01e+05
- 1 mm, 30 km, Mach 9.25 (2D boundary): Re_x = 1.1e+05
- 1 mm, 30 km, Mach 10.0: Re_x = 1.19e+05
- 1 mm, 30 km, Mach 12.0: Re_x = 1.42e+05
- 1 mm, 40 km, Mach 8.5: Re_x = 2.12e+04
- 1 mm, 40 km, Mach 10.0: Re_x = 2.49e+04
- 1 mm, 40 km, Mach 10.25 (2D boundary): Re_x = 2.55e+04
- 1 mm, 40 km, Mach 12.0: Re_x = 2.99e+04
- 3 mm, 20 km, Mach 8.5: Re_x = 5.4e+05
- 3 mm, 20 km, Mach 9.0 (2D boundary): Re_x = 5.72e+05
- 3 mm, 20 km, Mach 10.0: Re_x = 6.36e+05
- 3 mm, 20 km, Mach 12.0: Re_x = 7.63e+05
- 3 mm, 30 km, Mach 8.5: Re_x = 1.1e+05
- 3 mm, 30 km, Mach 9.75 (2D boundary): Re_x = 1.26e+05
- 3 mm, 30 km, Mach 10.0: Re_x = 1.3e+05
- 3 mm, 30 km, Mach 12.0: Re_x = 1.56e+05
- 3 mm, 40 km, Mach 8.5: Re_x = 2.32e+04
- 3 mm, 40 km, Mach 10.0: Re_x = 2.73e+04
- 3 mm, 40 km, Mach 10.75 (2D boundary): Re_x = 2.93e+04
- 3 mm, 40 km, Mach 12.0: Re_x = 3.27e+04
- 10 mm, 20 km, Mach 8.5: Re_x = 7.03e+05
- 10 mm, 20 km, Mach 9.25 (2D boundary): Re_x = 7.65e+05
- 10 mm, 20 km, Mach 10.0: Re_x = 8.27e+05
- 10 mm, 20 km, Mach 12.0: Re_x = 9.93e+05
- 10 mm, 30 km, Mach 8.5: Re_x = 1.43e+05
- 10 mm, 30 km, Mach 10.0: Re_x = 1.69e+05
- 10 mm, 30 km, Mach 10.25 (2D boundary): Re_x = 1.73e+05
- 10 mm, 30 km, Mach 12.0: Re_x = 2.03e+05
- 10 mm, 40 km, Mach 8.5: Re_x = 3.02e+04
- 10 mm, 40 km, Mach 10.0: Re_x = 3.55e+04
- 10 mm, 40 km, Mach 11.5 (2D boundary): Re_x = 4.08e+04
- 10 mm, 40 km, Mach 12.0: Re_x = 4.26e+04
- largest at a 2D boundary: 7.65e+05; largest up to Mach 12: 9.93e+05

## Eq. 2 and its cap
- 1D settled state, Mach 9, 30 km, T_w 2700 K: Eq. 2/fc 1.198, frozen/fc 0.607, cap active True
- 2D boundary, Mach 8.5, 20 km, T_w 2300 K: Eq. 2/fc 1.015, frozen/fc 0.832, cap active True
- B1, Mach 8, 35 km, T_w 2320 K: Eq. 2/fc 1.051, frozen/fc 0.732, cap active True
- 2D boundary, Mach 10.25, 40 km, T_w 2300 K: Eq. 2/fc 0.685, frozen/fc 0.480, cap active False

## Interfacial gas pressure (kPa) by temperature: {'2073.0': 88.64210435936914, '2200.0': 279.34287969079486, '2300.0': 627.2999942962166, '2400.0': 1311.363126400134}
- front.jsonl, fully catalytic: {'n_cases': 126, 'n_exceed': np.int64(45), 'n_passing': 63, 'n_passing_exceed': np.int64(18), 'ratio_span': [np.float64(0.10995910293498308), np.float64(4.942564937066284)]}

## Inside the time tolerance: baseline columns 3 ([{'R': '1.0', 'cat': 'fc', 'alt': '25.0', 't': 1794.7534982949066}, {'R': '3.0', 'cat': 'nc', 'alt': '32.0', 't': 1786.7694404412814}, {'R': '10.0', 'cat': 'fc', 'alt': '26.0', 't': 1786.1215553505226}]); variant entries 14

## Converted plus depleted depth at the highest passing 2D conditions (um): {'1.0': [54.74813058680202, 212.9490394492014], '3.0': [203.26655700225498, 266.2624746243673], '10.0': [271.343627231489, 338.7924957174673]}

## Five-minute boundary minus 30-minute boundary (fully catalytic, grid step 0.5): {'10.0': {'gain': [0.0, 0.5], 'modes': ['none', 'temperature']}, '1.0': {'gain': [0.25, 0.75], 'modes': ['glass', 'temperature']}, '3.0': {'gain': [0.0, 0.75], 'modes': ['glass', 'temperature']}}

## PRCC with bootstrap 95 % intervals
- B1: emissivity, glass 0.69 [0.6, 0.77]; heating correlation (x) -0.68 [-0.76, -0.59]; mass-transfer coefficient (x) -0.62 [-0.75, -0.5]; evaporation coefficient -0.57 [-0.69, -0.43]; oxide growth rate kp (x) 0.47 [0.29, 0.62]; temperature limit (K) 0.46 [0.31, 0.57]; O2 pressure exponent n_p -0.34 [-0.48, -0.2]; SiO2 from depleted layer lost -0.23 [-0.4, -0.07]; passive-to-active offset (decades) 0.14 [-0.02, 0.29]; heat capacity (x) 0.13 [-0.01, 0.25]; depletion onset T_dep (K) -0.12 [-0.28, 0.04]; SiC depletion rate (x) 0.08 [-0.08, 0.25]; ZrO2 scale conductivity (W/m/K) -0.05 [-0.2, 0.11]; substrate conductivity (x) 0.03 [-0.13, 0.19]; emissivity, oxide (after glass loss) 0.01 [-0.13, 0.15]; catalysis bound -0.01 [-0.17, 0.14]; glass conductivity (W/m/K) 0.0 [-0.15, 0.17]
- B2: temperature limit (K) 0.98 [0.97, 0.98]; heating correlation (x) -0.87 [-0.9, -0.79]; emissivity, glass 0.7 [0.57, 0.76]; catalysis bound -0.32 [-0.44, -0.15]; substrate conductivity (x) -0.23 [-0.36, -0.06]; heat capacity (x) 0.19 [0.03, 0.32]; passive-to-active offset (decades) -0.18 [-0.3, -0.02]; mass-transfer coefficient (x) -0.18 [-0.32, -0.02]; depletion onset T_dep (K) 0.17 [0.02, 0.3]; evaporation coefficient -0.11 [-0.26, 0.06]; O2 pressure exponent n_p 0.1 [-0.07, 0.26]; oxide growth rate kp (x) -0.08 [-0.22, 0.08]; ZrO2 scale conductivity (W/m/K) 0.05 [-0.1, 0.21]; SiO2 from depleted layer lost -0.05 [-0.2, 0.12]; emissivity, oxide (after glass loss) -0.02 [-0.17, 0.15]; SiC depletion rate (x) 0.01 [-0.15, 0.16]; glass conductivity (W/m/K) 0.0 [-0.14, 0.16]
