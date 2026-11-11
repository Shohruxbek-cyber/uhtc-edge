# INPUTS

Every input value used by the code, with its source. Tags: **SOURCE** (value stated in the cited document, at the
location given), **DIGITIZED** (read off a figure; method and figure given), **ASSUMPTION** (no source; reason and
planned sensitivity case given), **DERIVED** (computed from sourced values; method given).
PDFs and extracted text of every open source are in `literature/pdf/` and `literature/text/`; notes in
`literature/notes/`. Paywalled sources are listed in `literature/NOT_OBTAINED.md`.

## 1. Physical constants

| Quantity | Value | Tag | Source |
|---|---|---|---|
| Molar gas constant R | 8.314462618 J mol⁻¹ K⁻¹ | SOURCE | CODATA 2018 (exact in the 2019 SI) |
| Stefan–Boltzmann constant σ | 5.670374419×10⁻⁸ W m⁻² K⁻⁴ | SOURCE | CODATA 2018 (exact in the 2019 SI) |
| Standard atomic weights N, O, Si, C, B, Zr, Ar | 14.007, 15.999, 28.085, 12.011, 10.81, 91.222, 39.95 g/mol | SOURCE | CIAAW abridged standard atomic weights, ciaaw.org/abridged-atomic-weights.htm (read 2 Nov 2026); molar masses computed from these in `thermo.py` |
| ft → m | 0.3048 (exact) | SOURCE | international foot definition |

## 2. Atmosphere (US Standard Atmosphere 1976, `literature/pdf/USSA1976.pdf`, NTRS 19770009539)

| Quantity | Value | Tag | Location |
|---|---|---|---|
| R* | 8.31432×10³ N m kmol⁻¹ K⁻¹ | SOURCE | doc p. 3 (Table 2 misprint corrected by the errata bound at PDF p. 242) |
| g0 | 9.80665 m s⁻² | SOURCE | Table 2, doc p. 2 |
| r0 | 6356.766 km | SOURCE | Table 2 (corrected by the errata, PDF p. 242) |
| P0, T0 | 101325 Pa, 288.15 K | SOURCE | Table 2 |
| γ | 1.40 | SOURCE | Table 2 |
| Sutherland β, S | 1.458×10⁻⁶ kg s⁻¹ m⁻¹ K⁻¹ᐟ², 110.4 K | SOURCE | Table 2; S corrected by the errata; viscosity law eq. (51), doc p. 19 ("fails for ... very high and very low temperatures" — used above 1000 K as an extrapolation, see D-HEAT-1) |
| Sea-level composition N2, O2, Ar, CO2, ... | Table 3 (F_i, M_i) | SOURCE | doc p. 3 |
| M0 | 28.9644 kg kmol⁻¹ | SOURCE | doc p. 9, eq. (21) |
| Layer base heights H_b and gradients L_M,b | 0, 11, 20, 32, 47, 51, 71, 84.852 km′; −6.5, 0, 1.0, 2.8, 0, −2.8, −2.0 K km′⁻¹ | SOURCE | Table 4, doc p. 3 |
| Check values Z = 20, 30, 40 km: T, P, ρ | 216.650 K, 5529.3 Pa, 8.8910×10⁻² kg m⁻³; 226.509 K, 1197.0 Pa, 1.8410×10⁻²; 250.350 K, 287.14 Pa, 3.9957×10⁻³ | SOURCE | Table I, doc pp. 61, 63, 65 (PDF pp. 76, 78, 80) |

## 3. Heating

| Quantity | Value | Tag | Location |
|---|---|---|---|
| Sutton–Graves form q = K (p_s/R_n)^½ (h_s − h_w) | eq. (33) | SOURCE | Sutton & Graves 1971, NASA TR R-376, doc p. 13 (PDF p. 15) |
| K for air (0.2320 O2 – 0.7680 N2 by mass) | 0.1113 kg s⁻¹ m⁻³ᐟ² atm⁻½; average error 3.3 %, max 9.8 % | SOURCE | Table II, doc p. 39 (PDF p. 41) |
| Range of the Sutton–Graves study | H 2.3–116.2 MJ/kg, p 0.001–100 atm, T_w 300 and 1111 K | SOURCE | abstract (NTRS record and report summary) |
| Wall Prandtl number range | 0.67–0.71 for planetary-entry gases | SOURCE | Sutton & Graves 1971, doc p. 15 |
| Fay–Riddell (equilibrium) q = 0.76 Pr^−0.6 (ρ_wμ_w)^0.1 (ρ_eμ_e)^0.4 [1 + (Le^0.52 − 1) h_D/H_s](H_s − h_w)(du_e/dx)^½ | eq. (6) | SOURCE | Hamilton 1982, NASA TM-84580, doc p. 4 (PDF p. 5); original Fay & Riddell 1958 paywalled (NOT_OBTAINED) |
| Newtonian velocity gradient (du_e/dx) = (1/r_n)[2(p_s − p_∞)/ρ_s]^½ | eq. (39a) | SOURCE | Tauber 1989, NASA TP-2914, doc p. 6 (PDF p. 10) |
| p_s = ρ_∞V²(1 − ε/2) + p_∞, ε = ρ_∞/ρ_s | — | SOURCE | Tauber 1989 doc p. 6 |
| Sphere (simplified) q = 1.83×10⁻⁴ (ρ_∞/r_n)^½ V³ (1 − h_w/H_s) | eq. (40) | SOURCE | Tauber 1989 doc p. 6 |
| Unswept cylinder q = 1.29×10⁻⁴ (ρ_∞/r_cyl)^½ V³ (1 − h_w/h_aw) | eq. (41), Λ = 0 | SOURCE | Tauber 1989 doc p. 7 (PDF p. 11) |
| Sharp cone q = 4.03×10⁻⁵ (ρ_∞ cos δ_c/x)^½ V^3.2 sin δ_c (1 − h_w/h_aw), h_aw = h_∞ + 0.40 V² | eq. (46) | SOURCE | Tauber 1989 doc p. 9 (PDF p. 13) |
| Wedge / flat plate q = 2.42×10⁻⁵ (ρ_∞ cos δ/x)^½ V^3.2 sin δ (1 − h_w/h_aw) | eq. (48) | SOURCE | Tauber 1989 doc p. 9 |
| Heating distribution on hemisphere / cylinder q/q_s = cos θ (to ≈70°) | — | SOURCE | Tauber 1989 doc p. 8, Figs. 5–6 (PDF p. 12) |
| Prandtl number Pr | 0.71 | ASSUMPTION | upper end of Sutton & Graves' wall range; Fay–Riddell enters only V4 (sensitivity Pr 0.67) |
| Lewis number Le | 1.4 | ASSUMPTION | commonly used for dissociated air; enters the non-catalytic factor and V4 (sensitivity Le 1.0) |
| **Heating distribution, baseline (revision, D-2D-7)** | cos φ on the arc to 70°; Lees / Kemp–Rose–Detra local similarity from the tangency on (one integral over arc and flank, modified-Newtonian p/p_s = (1 − ε)cos²ψ + ε, isentropic surface flow γ = 1.4, factor 1 + 0.1853 β^½); linear bridge 70°–tangency | SOURCE (method) / DERIVED (bridge) | Tauber 1989 p. 8; Cleary 1969 NASA TN D-5450 eqs. 2, 4, 5 (doc pp. 7–8); DeJarnette & Tai 1969 NASA CR-101678 eqs. 10–11 (doc p. 11), 50–51 (pp. 26–27); Crabtree et al. 1965 ARC R&M 3637 eqs. 3.4–3.5 (doc pp. 16–17); notes `HeatingDistribution_*.md` |
| Check data for the distribution | Crabtree Fig. 12 (5/10/15° hemisphere-cones, M 5–14, theory); Cleary Table III (15° sphere-cone, M 10.6, measured); Kolodziej et al. 1997 Fig. 2 (planar 5°, M 15, BLIMPK/PERFORM) | DIGITIZED / transcribed | `literature/digitized/Crabtree1965_fig12.csv`, `Cleary1969_TablesI-III_alpha0.csv`, `Kolodziej1997conf_fig2.csv` |
| Heating-distribution band | beyond 70°: ×0.66 (low) to ×1.5 (high; at least the Kemp–Rose–Detra shoulder curve of Tauber Fig. 5) | DERIVED (D-2D-7) | measured/theory 0.66–1.05 (Cleary Table III vs Crabtree Fig. 12c); mirrored for the unknown sign at 5°; `literature/digitized/Tauber1989_fig5.csv` (DIGITIZED) |

## 4. Thermochemistry (NIST-JANAF, NIST SRD 13, doi 10.18434/T42S31; fetched by `literature/janaf/fetch_janaf.py`, not redistributed)

All enthalpies, Gibbs energies of formation and equilibrium constants are interpolated (linear in T, 100 K
tables) from the JANAF text tables listed in `literature/janaf/README.md`: N2, O2, N, O, NO, CO, CO2, SiO, SiO2(g),
B2O3(g), B2O3(l), B2O3(cr,l), SiO2(l), SiO2(cr,l), SiC (β), ZrB2(cr,l), ZrO2(cr,l). Argon: ideal monatomic gas,
c_p = 5/2 R (exact for an ideal gas). DERIVED quantities: equilibrium air composition (N2, O2, NO, N, O, Ar),
vapour pressures p_B2O3(T), p_SiO(T, p_O2), p_SiO2(g)(T), reaction enthalpies.

## 5. Trajectory

| Quantity | Value | Tag | Location |
|---|---|---|---|
| X-43A flight 3 altitude and Mach number vs time, boost | curves | DIGITIZED | Karlgaard et al. 2005, AIAA 2005-5920 (NTRS 20050232845), Fig. 9(a,b), PDF p. 12; extracted from the PDF vector paths (solid line = BET), calibration residual ≤ 0.004 s, ≤ 0.001 axis units |
| Same, descent | curves | DIGITIZED | Fig. 11(a,b), PDF p. 15; calibration residual ≤ 0.05 s, ≤ 0.001 |
| Separation: t, altitude, Mach (BET) | 88.35 s, 109,653 ft, 9.702 | SOURCE | Table 3, PDF p. 13 |
| 88.8–111.8 s (engine test, between figures) | linear interpolation | ASSUMPTION | D-TRAJ-1 |

## 6. ZrB2–SiC composite properties (`src/uhtc/data.py`, `src/uhtc/materials.py`)

| Quantity | Value | Tag | Source / location | Sensitivity |
|---|---|---|---|---|
| k(T), ZrB2–20 vol% SiC | 98.7 W/m/K (373 K) … 69.0 (2273 K), table of 20 values; at 1273 K mean of 78.2 (cut-bar) and 82.0 (laser flash); constant outside 373–2273 K | SOURCE | NASA Ames TPSX database, material 28 (ManLabs data, AFML-TR-68-190), ±5 % to 2073 K, ±8 % above; `literature/notes/TPSX2006_ZrB2-20SiC_ManLabs.md` | ×0.70 / ×1.05 (INDEX_properties 2.1) |
| k composition factor | 1.17 (10 vol%), 1.00 (20 vol%), 0.92 (30 vol%) | DERIVED (revision; D-TH-4) | Van Laningham (2012) Fig. 27 (PDF p. 81), within-study ratios at 1200–2000 °C: 10.7/21.9 vol% SiC 1.15–1.20; 30 vol% linearly interpolated between 21.9 and 48.7 vol%: 0.94–0.91 (`literature/notes/CHECK_3_conductivity_vs_SiC_fraction.md`). Kim et al. 2012 (SPS, ≤1000 °C, Fig. 8): no drop with SiC (30/20 = 0.96–1.18). The draft's 0.75 (Zimmermann 2008 abstract / TPSX: different materials) is kept as a low bound | 0.75 (low bound); k multiplier ×0.70–1.05 |
| cp(T) | mass-weighted rule of mixtures of JANAF ZrB2(cr) and β-SiC; e.g. 20 vol%: 456 J/kg/K (298 K), 819 (2200 K) | DERIVED | JANAF tables; agrees with TPSX drop calorimetry (717–806 J/kg/K, 1273–2273 K, ±2 %) within ≈3 % | ±3 % |
| ρ (room temperature) | rule of mixtures of ZrB2 6090 and SiC 3210 kg/m³ (20 vol%: 5514) | DERIVED | ZrB2 theoretical density 6.09 g/cm³: Shugart 2014 thesis, PDF p. 16; SiC 3.21 g/cm³: TPSX K&K, PubChem 3.16–3.23 | — |
| ρ(T) | scaled with the measured expansion of ZrB2–20SiC: 5560 → 5460 → 5310 kg/m³ at 293, 1273, 2273 K | SOURCE | TPSX material 28 (density table) | — |

## 7. Oxide scale and oxidation kinetics (`src/uhtc/data.py`, `src/uhtc/oxidation.py`)

| Quantity | Value | Tag | Source / location | Sensitivity |
|---|---|---|---|---|
| ZrO2-layer parabolic constant, 1300/1400/1500 °C | 295, 280, 390 µm²/h (ZrB2–30SiC, stagnant air) | SOURCE | Shugart 2014, Table 9, PDF p. 75 | kp ×1/3, ×3 |
| Same, 1650/1700/1800 °C | 8 500, 28 000, 42 000 µm²/h (flowing O2, 900 sccm; R² 0.21–0.52) | SOURCE | Shugart 2014, Table 9, PDF p. 75 | kp ×1/3, ×3 |
| Arrhenius fits of the two sets | low: A = 8.5×10⁻¹³ m²/s, E = 31.5 kJ/mol; high: A = 2.2×10⁻³ m²/s, E = 326 kJ/mol | DERIVED (least squares in `data.py`) | — | — |
| Recession constants kp″ (ZrB2–20SiC, stagnant air) | 7.9×10⁻⁵, 1.3×10⁻³, 2.2×10⁻² mm²/min at 1327, 1627, 1927 °C | SOURCE | Opila & Halbig 2001, Table 1, PDF p. 6 (checked on the page image) | used for V6 and for the slope below 1573 K |
| Low-T continuation below 1573 K | E = 236 kJ/mol (Opila & Halbig, 1327–1627 °C points) | DERIVED | — | — |
| SiC-depletion onset T_dep | 1627 °C (1900 K) ± 25 K | SOURCE | Shugart 2014 pp. 79–87; Opila & Halbig 2001 (depletion at 1627 °C) | ±25 K |
| SiC-depleted layer growth | 43 ± 10, 74 ± 6, 104 ± 8 µm at 5, 15, 20 min, 1800 °C → k_d = 473 µm²/min (fit through origin, no incubation, no T dependence as Shugart concludes) | SOURCE + DERIVED | Shugart 2014, Fig. 44 caption, PDF p. 85 | k_d ×0.5, ×2 |
| Extrapolation band of k_p above 2073 K (D-OX-11) | E = 250 or 354 kJ/mol anchored at the fit value at 2073 K; k_p ×0.5 or ×3 | DERIVED (E) / SOURCE (level) | Oguri et al. 2011a, J. Japan Inst. Metals 75:207, Fig. 3 p. 210 (thicknesses; E DERIVED); Hu et al. 2013 Materials 6:1730 (ZS30 at 2233 K, 125 µm at 600 s vs 177 µm extrapolated); furnace points ≥ 1900 K 1.5–5× above (`digitized/HighT_oxidation_points.csv`) | band itself |
| Incubation of SiC depletion (option, off in the baseline; D-OX-10) | t_inc = 20, 15, 12.5, 2.5 min at 1627, 1650, 1700, 1800 °C (bracket midpoints), clock Σ dt/t_inc = 1 | DIGITIZED + DERIVED | Shugart 2014 Fig. 40 (PDF p. 81) and pp. 79–80: onset brackets ≤40, 10–20, 10–15, ≤5 min; `digitized/Shugart2014_fig40_SiC_depletion_map.csv` | with / without (V6, base cases) |
| Pressure exponent n_p of kp | 0 (calibration at 1 atm O2 applied unchanged) | ASSUMPTION (D-OX-6) | no usable evidence: the only point once quoted, Karlsdottir 2007, is known only from its abstract and is no longer used (author, 10 Nov 2026) | n_p = 0.5 |
| ZrB2 → ZrO2 volume expansion | 16.9 % | SOURCE | Campbell et al. 2021 (AIAA 2021-3168), p. 7 | — |
| Glass densities | SiO2 2200, B2O3 1800 kg/m³ (amorphous) | SOURCE | PubChem/HSDB (Merck Index) | — |
| Glass activity model | ideal solution (a = mole fraction) | ASSUMPTION (D-OX-3) | — | — |
| Glass kept in the ZrO2 skeleton | pore volume per m³ converted = 1 − (1 − f_SiC)(1.169) (rest goes to the outer layer) | DERIVED (D-OX-2) | — | — |
| SiC of the depleted layer → glass | all re-oxidised (f_SiO_lost = 0) | ASSUMPTION (D-OX-7) | Fahrenholtz 2007 mechanism (abstract) | f_SiO_lost = 1 |
| k of glass | 1.37 W/m/K (298 K value; no open high-T value) | SOURCE (RT), ASSUMPTION (high T) | PubChem/HSDB citing Ullmann | 0.7–2.8 |
| k of porous ZrO2 scale | 1.5 W/m/K | SOURCE (range) | Zhu & Miller 2002, NASA TM-2002-211481: 1.0–1.9 at 1316 °C | 1.0–2.5 |
| k of SiC-depleted layer | 30 W/m/K | ASSUMPTION | porous ZrB2 after SiC removal; minor effect | 10–60 (run: `run_extra.py`; < 0.1 %) |
| Evaporation coefficient α | 1 (diffusion-limited in all flight cases) | ASSUMPTION | — | 0.1 |
| Mass-transfer conductance | = heat-transfer conductance K_SG (p_s/R_n)^½ | ASSUMPTION (D-HEAT-3) | — | ×0.5, ×2 |

## 8. Surface: emissivity, catalysis, failure criteria

| Quantity | Value | Tag | Source / location | Sensitivity |
|---|---|---|---|---|
| Emissivity, glass-covered | 0.88 | ASSUMPTION within the sourced range 0.70–0.95 | Pellegrini 2022 (abstract, ≈0.90), Stewart 1997 RCG glass 0.89, Padovano 2018 0.877, Van Laningham 2012 0.73–0.92 (INDEX_properties 2.7) | 0.70–0.95 |
| Emissivity, after glass loss (ZrO2-rich) | 0.75 | ASSUMPTION within 0.50–0.85 | Liu 2021 (abstract) 0.85 at 1768 °C; Meng 2011 0.66; TPSX/Thomas 2002 0.50 at 2100 K | 0.50–0.85 |
| Catalysis | bounds: fully catalytic and non-catalytic (D-HEAT-2); variant: non-catalytic while glass exists, fully catalytic after loss (D-HEAT-4) | — | literature γO ≈ 10⁻³–0.1 on oxidized ZrB2–SiC (Scatteia 2008 abstract 0.08 at 1800 K; Stewart 1997 RCG fits ≤ 0.022) | bounds |
| SiC active-to-passive line (molecular O2) | log10 p_O2[bar] = 5.708 − 15936/T (45 points, 1473–1924 K, rms 0.28 decades; corrected in revision 1, see CHANGELOG.md) | DERIVED (our fit to DIGITIZED points) | Jacobson, Myers & Harder 2011 poster compilation (NTRS); points in `literature/notes/Jacobson2011_MoO8_poster_transition_points_DIGITIZED.csv`; spot-checked against the poster | — |
| Passive-to-active offset | 0.65 decade below the active-to-passive line | ASSUMPTION within the sourced range 0.3–1 | Jacobson et al. 2012, Fig. 9 (DIGITIZED) | 0.3, 1.0; atomic-O case +1.65 (Rosner & Allendorf points lie 1.6–1.7 decades below the fit; run: `run_extra.py`, no change) |
| Alternative passive-to-active line (sensitivity, D-OX-12) | p_O2 [Pa] = 10¹⁵ exp(−55850 K/T) | SOURCE (fitted law, C/SiC, air plasma; extrapolated in pressure) | Balat-Pichelin et al. 2015, 8th Eur. Symp. Aerothermodyn. Space Vehicles, HAL hal-02313822, PDF p. 4–5 (exponent checked on a rendered page image) | used instead of Eq. pa in 1D set D and the 2D "pa_balat" stage |
| Temperature limit T_limit | 2422 K (2149 °C) | SOURCE (design value) | Kolodziej 1997, NASA TM-112204, doc p. 2: reusable leading-edge temperature used for sharp-edge constraints | 2289 K (TPSX multiple-use limit for ZrB2/SiC, Thomas 2002); 2543 K (ZrB2–SiC eutectic, lower of 2543 ± 30 K Ordan'yan 1989 and 2570 K Tu 2008, `literature/notes/ZrB2-SiC_eutectic.md`). The earlier 3033 K (2760 °C single use, Kolodziej 1997 p. 5) was dropped: it is above the eutectic (D-FAIL-5) |
| Glass considered gone (draft rule, now a sensitivity case) | d_g < 0.1 µm after having exceeded it | ASSUMPTION | numerical threshold | 0.01–1 µm (run: `runs/sensitivity/run_extra.py`; 1 µm gives early failures, D-FAIL-3); 0.1/0.5/1 µm in the 2D boundary runs |
| Criterion (a), revision (D-FAIL-4) | d_c = K(T, x_B) p_O2,w / J_gas, J_gas = g ln(1 + Y_O2)/M_O2, Y_O2 = 0.2315 | DERIVED | Henry's law flux K p/d (Kajihara et al. 2004 p. 561); mass-transfer conductance g (D-HEAT-3) | permeability low/high below |
| O2 permeability of silica, base | K = 2.15×10⁻¹⁴ exp(−14169 K/T) mol m⁻¹ s⁻¹ Pa⁻¹ (= 0.75 N₁ B/p from CVD-SiC B = 285 exp(−117800/RT) µm²/h, 1473–1773 K) | DERIVED | Jacobson et al. 2005, NASA NTRS 20050060620, Table II p. 7 (citing Ogbuji & Opila 1995); `notes/GlassOxygenTransport_SiC-kp_protective-thickness.md` | low: Kajihara et al. 2004 eqs. 11–12 (D·S): 2.83×10⁻¹⁶ exp(−8471/T); high: Deal–Grove Si dry O2 (Jacobson Table II): 6.50×10⁻¹⁴ exp(−14914/T) |
| O2 permeability of B2O3 (liquid), base | K = 7.7×10⁻⁷ exp(−16600 K/T) mol m⁻¹ s⁻¹ Pa⁻¹ (O self-diffusion 1.44×10⁻⁹ cm²/s at 1010 K, E = 138 kJ/mol, × O2-equivalent density) | DIGITIZED + DERIVED, UNSURE | Schlichting 1982, Glastech. Ber. 55, Bild 3 (p. 169; Tokuda data), `digitized/Schlichting1982_Bild3_...csv`; `notes/GlassOxygenTransport_borosilicate.md` | high: Luthra as printed in Campbell et al. 2021 eq. 7 (sign-typo reading, UNSURE): 0.15/101325 exp(−16000/T) |
| Mixed B2O3–SiO2 glass | K linear in the B2O3 mole fraction | SOURCE (rule) | Campbell et al. 2021 eqs. 10–11; reproduces Schlichting Bild 4 within ×1.5 at 750 °C | — |

## 8a. Geometry and case inputs (added in revision 1, task 2.2)

Every input varied in Q4 or Q5 now has a row here or above and in Table 1 of the paper.

| Quantity | Baseline | Status | Source / justification | Sensitivity range |
|---|---|---|---|---|
| Nose radius R_n | 1, 3, 10 mm | study cases | chosen to span sharp edges (D-2D-4) | Q4: ×2/3, ×1.5 of the base case |
| SiC fraction | 20 vol% | study case | composition of the TPSX record; Shugart 2014 tests 30 vol% | 10, 30 vol% |
| Slab depth (1D) | 20 mm | ASSUMPTION (D-TH-2) | heat-soaked after about 20 s, the hotter limit | 50 mm (V2) |
| Flank length (2D) | 30 mm | ASSUMPTION (D-2D-4) | slender nose, adiabatic back plane at the end of the flank | 15, 60 mm |
| Cone half-angle (2D) | 5° | ASSUMPTION (D-2D-4) | slender nose | 10° |
| Flank heating multiplier (2D) | 1 | ASSUMPTION (D-2D-5) | scale on g(s) over the flank | ×0.5, ×2 |
| Evaporation coefficient, measured value for comparison | not used | SOURCE | SiO over Ta + Y2O3 + Y2O3·SiO2, 0.004–0.007 at 1560–1670 K, Jacobson 2015, NTRS 20150007506, Fig. 9 (DIGITIZED by eye): a different system, so the model keeps α = 1 | — |

## 8b. Glass flow (version 3, Task 1; `src/uhtc/glassflow.py`; off in every run of the paper)

| Quantity | Baseline | Status | Source / justification | Sensitivity range |
|---|---|---|---|---|
| Silica-melt viscosity, natural silica LA | log10 mu [Pa s] = -8.60 + 587.4 kJ/mol / (ln 10 R T) | SOURCE | Mauro et al. 2022, J. Am. Ceram. Soc. 105, 61-66, Table II, fibre-draw fit, valid 1973-2573 K (checked on arXiv:2108.09487, review/R6_citations.md); extrapolated above 2573 K (D-FLOW-5) | F320 law below |
| Silica-melt viscosity, fluorine-doped F320 | log10 mu [Pa s] = -6.80 + 495.8 kJ/mol / (ln 10 R T) | SOURCE | same table; the softest of the seven silicas at 2000-2200 K | low-viscosity case |
| Temperature of the viscosity | surface temperature of the upstream node | ASSUMPTION (D-FLOW-5) | the film is the outer layer, so T_s is its hottest point | bound in D-FLOW-5 |
| Prandtl number of the Reynolds analogy | 0.71 | ASSUMPTION | air; the value of the Supplement S9 estimate | none |
| Reynolds analogy c_f/2 = St Pr^(2/3) | laminar | SOURCE | Crabtree et al. 1965, ARC R&M 3637, Section 1.3 (D-FLOW-3) | wall shear x0.5-2 |
| Wall-shear multiplier | 1 | ASSUMPTION (D-FLOW-3) | same band as the flank heating, since the same g(s) carries both | x0.5, x2 |
| Edge expansion exponent gamma | 1.4 | ASSUMPTION | isentropic expansion from the stagnation point, as the S9 estimate and the Lees distribution | none |
| Pressure gradient | modified Newtonian, p/p_s = (1 - eps) cos^2(phi) + eps | SOURCE (model) | the ratio of the Lees distribution (D-FLOW-4) | term off |

## 8c. Thermoelastic module (version 3, Task 2; `src/uhtc/thermoelastic.py`; not coupled to the thermal solve)

| Quantity | Baseline | Status | Source / justification | Sensitivity range |
|---|---|---|---|---|
| Young's modulus E(T), ZrB2-20SiC | 530, 523, 446, 106 GPa at 296, 1073, 1673, 1873 K, linear between, held beyond | SOURCE | NASA TPSX record 28 (`literature/text/TPSX_ZrB2-20SiC_id28.txt`), 4-point bending in argon; as the Section 6 index (D-MECH-2) | Neuman set below |
| Young's modulus E(T), ZrB2-30SiC in air | 510, 210, 110 GPa at 296, 1773, 1873 K | SOURCE (abstract only) | Neuman et al. 2013, abstract values, as the Section 6 index | property set 2 |
| Thermal expansion alpha(T) | 5.00, 5.94, 6.33, 6.88e-6 /K at 673, 1073, 1473, 1873 K, held beyond | SOURCE | TPSX record 28; thermal strain alpha(T)(T - 300 K) (D-MECH-2) | used with both E sets |
| Poisson's ratio | 0.14 | ASSUMPTION (in the TPSX record) | TPSX record 28 | none |
| Bend strength sigma_f(T), mean and standard deviation | in-plane 430 +- 28.5, 548 +- 18.7, 463 +- 18.0, 342 +- 6.8, 239 +- 16.7 MPa at 296, 673, 1073, 1473, 1873 K | SOURCE | TPSX record 28 (4-point bend, 2.54 x 5.08 x 44.5 mm bars) | Neuman strengths 680, 750, 360, 360 MPa (abstract) for the strength ratio |
| Weibull modulus m | 6 and 20 (bracket) | DERIVED (D-MECH-3) | from the TPSX coefficients of variation by the exact two-parameter relation: 5.9-64, median about 25, specimen counts unknown | bracket itself |
| Effective volume of the strength test | Quinn's 4-point formula, outer span 38.1 mm, quarter- and third-point loading | ASSUMPTION (span not given) | D-MECH-3 | bracket itself |
| Restraint | free tip and restrained back plane | ASSUMPTION (D-MECH-1) | bracketing cases | both always reported |

## 8d. Integral flank model (version 3, Task 3; `conduction2d.integral_distribution`; not used in the baseline)

| Quantity | Baseline | Status | Source / justification | Sensitivity range |
|---|---|---|---|---|
| Laminar heating and momentum thickness | q = 0.22 Re_theta^-1 (rho*/rho_e)(mu*/mu_e) rho_e u_e (h_aw - h_w) Pr_w^-0.6; theta = 0.664 (int rho* mu* u_e r^2 ds)^1/2 / (rho_e u_e r) | SOURCE | Riley 1992, NASA TM-107630 (NTRS 19920016036), Eqs. 3.3.1-3.3.2, PDF pp. 45-46 (Zoby et al.) | model option itself |
| Intermediate (reference) enthalpy | i* = i_e + 0.5 (i_w - i_e) + 0.22 (i_r - i_e) | SOURCE | Crabtree et al. 1965, ARC R&M 3637, Eq. 1.3 (Eckert) | none |
| Laminar recovery factor | Pr^1/2, Pr = 0.71 | SOURCE | Crabtree et al. 1965, Eq. 1.4 | none |
| Edge state | isentropic expansion (gamma 1.4) from the stagnation point to the modified-Newtonian pressure | ASSUMPTION (D-HEAT-7, option (a)) | as the Lees distribution | sharp-cone model as the swallowed limit |

## 9. Validation targets

- V1: analytic solutions (no external data).
- V3: Kolodziej 1997, Table 1 (doc p. 3) and DIGITIZED Fig. 1 points (`literature/digitized/Kolodziej1997_fig1.csv`).
- V6: furnace data in `literature/notes/oxidation_furnace_validation_points.csv` (Shugart 2014; Opila & Halbig 2001;
  Levine & Opila 2003; Opila et al. 2003 preprint; the two Karlsdottir 2007 abstract values are no longer used, author 10 Nov
  2026); see `validation/V6_oxidation/report.md`.
- V7: Savino et al. 2010 (Open Aerosp. Eng. J. 3:20) and Di Maso 2009 thesis; full input/target table in
  `literature/V7_TARGET.md`; DIGITIZED temperature histories in `literature/digitized/Savino2010_TOAEJ_ArcJetUHTC_fig8.csv`.
- V4: no external data (correlation cross-check).
