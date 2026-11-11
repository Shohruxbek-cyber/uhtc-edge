# Thermal-stress screening during heat-up (task 2.6)

Indicative only, not a failure criterion. Script: `audit/tools/thermal_stress_screen.py`; numbers:
`audit/thermal_stress_screen.csv`.

## Result

For a cold start (300 K) at constant flight conditions on the fully catalytic 2D boundary, the fully restrained stress
index σ_th = E α ΔT/(1 − ν) exceeds the measured flexural strength at every condition, within 0.1 s of the start.

| Condition | Peak ΔT within 2 R_n | Peak index / strength, 2 R_n (TPSX; Neuman) | 5 R_n (TPSX; Neuman) | Peak reached at |
|---|---|---|---|---|
| Mach 8.5, 20 km, 1 mm | 363 K | 2.8; 1.4 | 4.5; 2.3 | 3 s |
| Mach 9.25, 30 km, 1 mm | 260 K | 2.2; 1.0 | 3.7; 1.7 | 10 s |
| Mach 10.25, 40 km, 1 mm | 214 K | 1.9; 0.8 | 3.4; 1.5 | 10–30 s |
| Mach 9.75, 30 km, 3 mm | 476 K | 3.8; 1.7 | 5.6; 3.0 | 10–30 s |
| Mach 10.25, 30 km, 10 mm | 825 K | 6.1; 2.9 | 6.9; 4.0 | 30 s |

Ratios use properties at the coldest point of the region. With TPSX properties at the region's mean temperature they reach
7.2 (2 R_n) and 7.9 (5 R_n). All are upper bounds for the elastic stress at a given ΔT, for three reasons:

- A real tip is only partly restrained, and σ_th is the fully restrained stress.
- A step entry is the fastest heat-up the model can give. A flown trajectory heats more slowly.
- The hot surface is in compression; tension appears in the cooler interior, and the comparison with bend strength ignores
  the actual stress state.

The screen therefore does not show that tips crack. It shows that thermal stress cannot be dismissed: the index would have
to overstate the real stress by a factor of up to about 7 for the tip to stay below its strength, so a proper thermoelastic
calculation is needed before a survival claim can include fracture. Larger radii score worse, the opposite of the oxidation
trend, because they hold larger temperature differences.

The simple critical temperature difference R = σ_f (1 − ν)/(E α) is 85–300 K from the TPSX data (lowest near 1673 K, highest
above 1873 K, where E collapses), against 160–1030 K of ΔT in the tip.

## Inputs

ZrB₂–20 vol% SiC, NASA TPSX record 28 (`literature/text/TPSX_ZrB2-20SiC_id28.txt`), interpolated, held constant past the
end points: E 530 GPa (296 K), 523 (1073 K), 446 (1673 K), 106 (1873 K); α 5.00–6.88×10⁻⁶ K⁻¹ (673–1873 K); ν 0.14 (marked
assumed in TPSX); in-plane bend strength 430, 548, 463, 342, 239 MPa (296, 673, 1073, 1473, 1873 K), 252 MPa through
thickness at 2073 K.

Cross-check, the second number in each table cell: Neuman et al. 2013 (ZrB₂–30SiC, air, abstract values, interpolated
linearly) give E of 510 GPa at room temperature, 210 at 1773 K and 110 at 1873 K, and strength of 680 MPa at room
temperature, 750 at 1073 K and about 360 at 1773–1873 K. That material is stronger and loses stiffness sooner than the TPSX
one, which halves the ratios. It has more SiC than the model's composite, so TPSX stays the primary set.

ΔT is the spread of nodal temperatures within a distance L of the stagnation point (L = 2 R_n, 5 R_n), from the production
2D mesh with the baseline heating distribution, at 0.1, 0.3, 1, 3, 10 and 30 s.
