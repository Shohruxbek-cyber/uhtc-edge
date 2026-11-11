# V9 plan: glass-loss clock against a plasmatron test

Written and committed before any V9 run. Nothing below may change after the run. No model parameter is adjusted.

## Why this check

V1–V8 test conduction, heating, oxide growth and surface temperature. None tests the glass-loss clock, which sets the 30-minute boundary in most columns of the 2D maps. Marschall et al. report spontaneous temperature jumps on ZrB₂–30SiC, which they associate with the loss of the protective silica glass, with enough detail to simulate.

## Source

Marschall, Pejaković, Fahrenholtz, Hilmas, Panerai, Chazot, "Temperature Jump Phenomenon during Plasmatron Testing of ZrB₂-SiC Ultra-High Temperature Ceramics". Read as the edited manuscript, Appendix 1 of DTIC ADA564348 (`literature/pdf/Marschall2011_AFOSR_final_report_ADA564348.pdf`). Published as J. Thermophys. Heat Transfer 26 (2012) 559, which was not seen. Page numbers are those of the appendix.

## Test case

ZrB₂–30SiC, VKI Plasmatron, static pressure 10 kPa, 16 g/s of air. Two power stages in one run: low power for about 6 min, then high power. Three samples are used, because Fig. 4 shows these with their exact traces:

| Sample | Low stage q_cw (W/cm²), H_e (MJ/kg), P_dyn (Pa) | High stage q_cw, H_e, P_dyn | Measured outcome |
|---|---|---|---|
| 1.8 | 109, 12.0, 30 | 185, 17.7, 60 | No jump |
| 1.5 | 109, 12.0, 30 | 202 (marked as estimated), 18.6, 75 | Jump |
| 1.6 | 109, 12.0, 30 | 236, 21.7, 90 | Jump |

Source: Table 1, A1-11. Sample 1.7 also jumps, but it is not in Fig. 4, so its times can only be read coarsely from Fig. 3; it is left out.

## Measured targets, read by eye from Fig. 4 (A1-13), about ±3 s and ±15 K

| Sample | Power switch (test time, s) | Low-stage plateau (K) | High-stage plateau before any jump (K) | Jump onset (test time, s) | Onset after the switch (s) |
|---|---|---|---|---|---|
| 1.8 | 358 | 1895 | 2110 | none up to 480 s | none in 122 s |
| 1.5 | 357 | 1825 | 2210 (text: near 2215, sample 1.5 T_L = 2206, A1-17) | 425 | 68 |
| 1.6 | 359 | 1860 | 2250 at onset | 380 | 21 |

## Inputs to the model

- `coupled.run_1d`, production settings, as in V8.
- `PrescribedConditions` with the measured cold-wall flux: Φ = q_cw/(H − h(300 K)) and mass-transfer conductance equal to Φ.
- Stagnation pressure p_s = 10 kPa + P_dyn of each stage.
- f_SiC = 0.3, `materials.composite(0.3)`, baseline `OxParams`.
- Slab depth 7 mm with an adiabatic back face, the target cap thickness (A1-5). Real caps are up to 2 mm thinner. The stem is ignored.
- Start at 300 K at t = 0, stage change at the switch time above, end at 480 s. Both catalysis bounds.
- O₂ counted as in the baseline: 0.2095 p_s. The source computes atomic oxygen at the surface; the model does not.

## Pass criterion

Following V8, a case passes if at least one catalysis bound meets every band, and V9 passes if all three samples pass with the same bound.

1. **Outcome.** Glass loss (criterion a) occurs in the high-power stage for samples 1.5 and 1.6, and no glass loss occurs by 480 s for sample 1.8.
2. **Timing.** For samples 1.5 and 1.6, the predicted time of glass loss after the switch is within a factor of 3 of the measured onset (the V6 band).
3. **Temperature.** The low-stage plateau, and the high-stage temperature at the measured onset time (1.5, 1.6) or at 480 s (1.8), are within ±80 K of the measured values (the V8 band).

## Reporting

The result goes into Table 2 and Section 3 whether it passes or fails, and into every place that summarizes validation: Abstract, Limitations and Conclusions. If V9 fails, it fails: the criterion stays as written here.
