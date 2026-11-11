# VVU plan: validation uncertainty of checks V6–V9

Written and committed before any run or calculation of this analysis (10 Nov 2026). The analysis is **post hoc**: it was
added after the first results of V6–V9 were known. The preset criteria and verdicts of V6–V9 (Table 2, `validation/V9_plan.md`,
the `report.md` and result files of each check) are unchanged, and the words "pass" and "fail" keep their meaning. Nothing
below may change after this commit. If a step proves impossible, a dated `## Amendment` section at the end says what changed
and why; the text above it stays as committed. No model code (`src/uhtc/`) changes and no parameter is adjusted.

## 1. Purpose

V6–V9 fail their criteria. This analysis asks, quantity by quantity, whether the model error is larger than the combined
uncertainty of the numerics, the test conditions, the model inputs and the measurement, or whether the test is too
uncertain to tell. It also reports each measured quantity on its own, in addition to the single verdict per check.

## 2. Method

ASME V&V 20-2009 (R2021), *Standard for Verification and Validation in Computational Fluid Dynamics and Heat Transfer*,
ASME, New York (reaffirmed 2021, as the ASME product page shows).

For each scalar quantity, with S the simulation and D the measurement:
- comparison error E = S − D;
- validation standard uncertainty u_val = √(u_num² + u_BC² + u_input² + u_D²), where u_BC and u_input split the input
  uncertainty of V&V 20 into the test boundary conditions and the model inputs;
- expanded uncertainty U_val = 2 u_val.

Temperatures are compared in kelvin. Times and thicknesses are positive and spread over factors, so for them
E = ln(S/D) and every uncertainty is a standard deviation of a logarithm.

S is the baseline simulation of the existing check (sample 0, Section 7). Each catalysis bound is a separate case, because
a bound is not a random input: E and U_val are reported for the fully catalytic and the non-catalytic bound side by side.

## 3. Decision rule (fixed now)

- |E| > U_val: "model error indicated at the test uncertainty". The interval E ± U_val brackets the model-form error.
- |E| ≤ U_val: "not distinguishable at the test uncertainty". The check cannot separate model error from uncertainty.
- Outcome quantities (Section 4): the share of the joint samples (Section 6) that reproduce the measured outcome. Below 5 %:
  "outcome not reproduced within input uncertainty".
- Where u_D or u_BC holds only part of the real uncertainty because the source states none (Section 5), U_val is a lower
  bound, and a label "model error indicated" then holds only at the stated uncertainty. The tables mark these quantities.

These labels are in addition to the preset pass/fail. They never replace it.

## 4. Quantities

| Check | Case (model) | Quantity | D | Kind |
|---|---|---|---|---|
| V7 SPES, ZrB2–15SiC hemisphere [Savino 2010; Di Maso 2009] | 2D nose-tip model of the hemisphere, as in `run_v7.py` (`2d-nc`, `2d-fc`) | plateau temperatures at the end of steps 1–5 | 1573, 1723, 1843, 1923, 2053 K | scalar, K |
| | | oxide (ZrO2 + glass) at the stagnation point at power-off | 150 µm (maximum) | scalar, log |
| | | SiC-depleted layer at the stagnation point | none observed | outcome: δ_d < 1 µm |
| V8-1 Ames, ZrB2–20SiC flat face [Bull 1998] | 1D, passive condition, 300 s, as `run_v8.py` | peak surface temperature | 1630 °C | scalar, K |
| | | oxide (ZrO2 + glass) at 300 s | 110 µm | scalar, log |
| | | SiC-depleted layer at 300 s | 40 µm | scalar, log |
| V8-2 HIT, ZrB2–30SiC [Hu 2013] | 1D, as `run_v8.py` | steady temperature (end of run), condition 1 | 1960 °C | scalar, K |
| | | steady temperature (end of run), condition 2 | 1900 °C | scalar, K |
| | | reacted layer, condition 1, at 201, 402 and 600 s | 65, 99, 125 µm | scalar, log; reported for oxide + depleted and for oxide alone, because the source does not define the layer |
| | | oxide, condition 2, at 600 s | 90 µm | scalar, log |
| | | trend with SiC content, condition 1 | T(ZS10) > T(ZS15) > T(ZS30) | outcome |
| V9 Plasmatron, ZrB2–30SiC [Marschall et al.] | 1D, as `run_v9.py`, samples 1.8, 1.5, 1.6 | low-stage plateau | 1895, 1825, 1860 K | scalar, K |
| | | high-stage temperature at the measured onset (1.5, 1.6) or at 480 s (1.8) | 2110, 2210, 2250 K | scalar, K |
| | | glass-loss time after the switch (1.5, 1.6) | 68, 21 s | scalar, log; a sample without glass loss by 480 s is censored at the end of the run and counted as such |
| | | jump outcome | glass loss in the high stage for 1.5 and 1.6, none by 480 s for 1.8 | outcome |
| V6 furnace | no runs | per-quantity statistics of the existing points (`per_quantity.py`, Section 10) | | |

V8-1: the enthalpy is not reported. S is given at the three bracketing enthalpies of the check (15, 30 and 45 MJ/kg), and E
and U_val use 30 MJ/kg as the nominal case, with the enthalpy sampled as a test condition (Section 5.2).

## 5. Uncertainty sources

Four groups, reported separately, because the split between them is the main result.

### 5.1 Numerical, u_num

The V2 results, each bound taken as the half-width of a uniform distribution, so u = bound/√3:
- 1D (V8, V9): the largest baseline-to-finest differences of V2: temperature 0.08 %, failure time 0.44 %, oxide thickness 0.03 %.
  The glass-loss time takes the failure-time value and every thickness the oxide value.
- 2D (V7): the V2 criteria, temperature 0.5 % and failure time 2 %. V2 did not refine the 2D oxide, so the 2D thicknesses
  take the 2 % failure-time value.

### 5.2 Test boundary conditions, u_BC

Only uncertainties that the sources state, or that follow from a rule written here. Where a source states none, the item
gets no spread and the table marks u_BC as a lower bound. Source locations are those of the source reports of this
analysis (`validation/VVU/sources.md`, written with this plan).

- **V7.** Total enthalpy and Pitot pressure: no uncertainty stated (Savino et al. Table 3; Di Maso p. 52). The heat flux
  was computed by the authors, with no stated accuracy: only for step 5, a CFD value of 5 MW/m² at the hemisphere
  stagnation point for a non-catalytic cold wall at 300 K (Di Maso p. 56, Fig. 51; Savino et al. Fig. 20). Rule: r is the
  ratio of that CFD flux to the model's non-catalytic cold-wall flux at the same total enthalpy and Pitot pressure (step 5,
  R = 5 mm). The model's heating is multiplied by a factor drawn uniformly from 1 − |1 − r| to 1 + |1 − r|, for both bounds and
  all steps. The non-catalytic flux is used because the CFD wall was non-catalytic. The Table 1 band of the heating constant
  (×0.85–1.15) is then left out of u_input for V7, so that the heating is counted once.
- **V8-1.** Enthalpy: not reported; uniform over the existing bracket, 15–45 MJ/kg. Calorimeter flux and pressure: no
  uncertainty stated ("approximately", "nominal"; patent col. 7). Calorimeter catalysis: see below.
- **V8-2.** Flux, enthalpy and pressure: no uncertainty stated (Hu et al., Section 2 and Table 1). Calorimeter catalysis:
  see below.
- **V9.** Flux, edge enthalpy, static and dynamic pressure: no uncertainty stated (Marschall et al., A1-6, A1-9, A1-11); the
  dynamic pressures are rounded to 5 Pa, which is not counted. Calorimeter catalysis: see below.
- **Calorimeter catalysis (V8, V9).** All three tests used a water-cooled copper calorimeter. The checks convert its cold-wall
  flux as if it were fully catalytic at 300 K (`coupled.PrescribedConditions`). Marschall et al. rebuilt their edge enthalpy
  with a copper recombination efficiency of 0.01 at 350 K (A1-9) and cite measured values of 0.014–0.067 for cold copper
  (Krassilchikoff et al., A1-26). The other two sources state no efficiency. Rule: the recombination efficiency γ_Cu of the
  calorimeter is drawn log-uniformly from 0.01 to 1, which spans the value of the V9 rebuild, the cited measurements and the
  fully catalytic assumption of the checks. A partly catalytic calorimeter reads less than a fully catalytic one, so the
  measured flux is divided by the ratio q(γ_Cu)/q_fc of the calorimeter, computed for a frozen boundary layer with first-order
  wall recombination (Goulard form):
  q(γ)/q_fc = f_nc + φ(γ)(1 − f_nc), φ = 1/(1 + Φ/(ρ_w k_w)), k_w = γ √(R T_cal/(2π M_O)),
  with f_nc the model's non-catalytic factor of the step, Φ the step's heat-transfer conductance from the measured flux,
  ρ_w the air density at the stagnation pressure and T_cal, M_O = 0.016 kg/mol and T_cal = 300 K (350 K for V9, the value
  of the source). The correction scales the heat- and mass-transfer conductances of each step. For V9 the edge enthalpy is
  held at its rebuilt value: the source gives the direction of its dependence on γ_Cu (higher γ_Cu, lower H_e) but not its size.
  A separate deterministic case at γ_Cu = 0.01 and 0.067 reports the sign and size of the correction on its own.
- **Pyrometry.** None of the four sources states an emissivity assumed for a two-colour reading, so no emissivity correction
  enters u_D. The Ames patent names "an optical pyrometer" with no type.

### 5.3 Model inputs, u_input

The ranges of Table 1 that act in these tests, with the distributions of Section 4.4 (`runs/sensitivity/run_sensitivity.py`,
PARAMS): uniform, or log-uniform where marked.

| Input | Range | Distribution |
|---|---|---|
| glass emissivity | 0.70–0.95 | uniform |
| oxide emissivity | 0.50–0.85 | uniform |
| substrate conductivity multiplier | 0.70–1.05 | uniform |
| heat capacity multiplier | 0.97–1.03 | uniform |
| k_p multiplier | 1/3–3 | log-uniform |
| k_p multiplier above 2073 K (extrapolation band) | 0.5–3 | log-uniform |
| k_d multiplier | 0.5–2 | log-uniform |
| depletion onset T_dep | 1875–1925 K | uniform |
| passive-to-active offset | 0.3–1.0 decade | uniform |
| mass-transfer conductance multiplier | 0.5–2 | log-uniform |
| glass conductivity | 0.7–2.8 W/m/K | uniform |
| ZrO2-layer conductivity | 1.0–2.5 W/m/K | uniform |
| SiC-depleted-layer conductivity | 10–60 W/m/K | uniform |
| O2-pressure exponent n_p | 0–0.5 | uniform |
| SiO2 of the depleted layer lost as gas | 0–1 | uniform |
| evaporation coefficient α | 0.1–1 | log-uniform |
| O2 permeability law of the glass | low, base, high | discrete, equal shares |

Left out, because they are fixed by the test or do not act on the measured quantities: the catalysis bound (a case),
the SiC fraction and the geometry (set by each test), the temperature limit (a failure clock only) and, for V8 and V9, the
heating constant (the flux is measured). For V7 the heating constant is counted in u_BC (Section 5.2). No input outside
Table 1 is added.

### 5.4 Measurement, u_D

Stated uncertainties, each taken as the half-width of a uniform distribution unless the source gives a standard
uncertainty, combined with the reading uncertainty of values read from figures:
- V7 temperatures: two-colour pyrometer, ±1 % of reading and a repeatability of about 2 % (Di Maso p. 23). Step 3 is given
  as 1843 K by Savino et al. and 1823 K by Di Maso (Table 12): ±10 K is added for that step. Oxide: no uncertainty stated;
  the value is the maximum of 150–50 µm over the hemisphere (Di Maso p. 55). u_D = 0, a lower bound.
- V8-1: no uncertainty stated for the temperature or the layers. u_D = 0, a lower bound.
- V8-2: temperatures are tabulated (Hu et al. Table 2) with no uncertainty: u_D = 0, a lower bound. Layers read from
  Fig. 10: the recorded reading uncertainty of ±1.5 µm (`literature/digitized/Hu2013_ZrB2SiC_HFPWT_fig10.csv`), a lower bound.
- V9: two-colour radiometer ±10 °C (A1-7), combined with the reading uncertainty of ±15 K (`V9_plan.md`). Times: reading
  ±3 s (`V9_plan.md`) and the 1 Hz sampling (A1-7), ±0.5 s.

## 6. Sampling

- Latin hypercube over groups (b) and (c) together, seed 20261110, N = 100 per case and bound (the joint sample).
- Each group alone, the other at baseline: N = 50 each, seeds 20261111 (b only) and 20261112 (c only), to split the variance.
- u_BC and u_input are the sample standard deviations of S over the b-only and c-only samples (of ln S for log quantities).
  The joint sample gives the outcome shares and a check that the two groups add.
- One draw of the inputs serves every run of a case: the five V7 steps, the three V9 samples, the four V8-2 runs.
- A run that fails numerically is reported and left out; it is not rerun with other settings.

## 7. Reproducibility gate

Sample 0 of every case has every input at its baseline and no calorimeter correction. It must reproduce the stored values
of `v7_results.json` (`2d-nc`, `2d-fc`: Tend, oxide_um, dep_um), `v8_results.json` (V8_1_H15/H30/H45 passive and V8_2: T_max,
T_end, oxide_um, dep_um) and `v9_results.json` (T_low, T_ref, t_glass) to a relative difference of 1e-9. If it does not, the
analysis stops and is reported; it does not continue.

## 8. Cost rule (fixed now)

Three samples of each case are timed first. If the full plan would take more than 8 hours on this machine, the joint
sample of the V7 (2D) cases is reduced to N = 50; the group-only samples stay at 50. Applying this rule is not a change of plan.

## 9. Outputs and what will be written, whatever the outcome

- `validation/VVU/run_vvu.py` (imports the run functions and constants of `run_v7.py`, `run_v8.py` and `run_v9.py`;
  a thin wrapper adds the sampled inputs), `vvu_samples.csv` (every sample, inputs and outputs), `vvu_results.json`, and
  `report.md` with, per check: quantity, bound, S, D, E, u_num, u_BC, u_input, u_D, U_val, the label of Section 3 and the
  share of U_val² from each group.
- Paper: Supplementary Section S10 "Validation uncertainty (added after the first results)", one paragraph in Section 3 of
  the main text, and one sentence each in Section 6 and in the Conclusions. These are written whichever way the result
  falls, and every place says that the analysis was added after the first results and that the preset verdicts stand.

## 10. Per-quantity statistics (no runs)

`validation/VVU/per_quantity.py` reads the existing result files and splits the points of each check by quantity class
(ZrO2 layer, total oxide, recession, SiC-depleted layer and, for V6 and outside its criterion, glass), keeping the
calibration set, the independent points and the two zero measurements apart. For each class: n, the number within the
check's own band, the geometric-mean ratio and the geometric standard deviation, with and without the incubation option
where it applies. It must reproduce the values in the text (V6: 15 of 18 oxide and recession points, geometric mean 1.09;
9 of 19 depleted-layer points); if it does not, the analysis stops and is reported.

## Amendment 1 (10 Nov 2026, after the reproducibility gate, before any sample was run)

The gate of Section 7 passed: sample 0 reproduced every stored value of V7, V8 and V9 within a relative 1e-9
(`vvu_gate.json`). The same step gave the V7 heating ratio of Section 5.2: the authors' CFD flux of 5 MW/m² is
r = 2.00 times the model's non-catalytic cold-wall flux of 2.50 MW/m² at step 5. The committed rule draws the heating
multiplier uniformly from 1 − |1 − r| to 1 + |1 − r|, which for r = 2.00 runs from 0.003 to 2.0, that is, down to no heating
at all. A ratio is symmetric in its logarithm, and a multiplier near zero has no physical meaning. The rule is therefore
changed to: the multiplier is drawn log-uniformly from 1/r to r (0.50 to 2.00). Nothing else changes. At the time of this
amendment no sample of any case had been run and no comparison error had been computed, other than the gate's baseline values,
which are those of the existing checks.

## Amendment 2 (10 Nov 2026, after the sample runs, before the final analysis)

Four of the 1612 runs did not terminate. Their workers ran at full CPU for 56 minutes without returning, while the slowest
completed run of the same case took 182 s (V7) and 5 s (V8-1). They are: V7, fully catalytic, joint sample, draws 4, 38 and
83 (heating multipliers 1.93, 1.70 and 1.88), and V8-1, fully catalytic, input-only sample, draw 30 (k_p multiplier 2.66).
Section 6 says that a run that fails numerically is reported and left out, and is not rerun with other settings, but it does not
say how long to wait. The four runs were stopped and are recorded in `vvu_runs.jsonl` as numerical failures with their inputs.
Nothing else changes. Effect: the V7 fully catalytic joint sample has 97 runs, which enter only the outcome share of the depleted
layer, and the V8-1 fully catalytic input-only sample has 49. The three V7 runs lie at the hot end of the heating range, so the
joint spread of V7 with the fully catalytic bound is somewhat narrower than the full sample would give. Before this amendment the analysis script had been run once on the partial output, as a
test of the script; no range, sample size or rule was changed after that run.

## Amendment 3 (10 Nov 2026, after the results and after the independent review)

This amendment records changes made after results existed, so that none of them is hidden. None changes a range, a sample size
of the plan or the decision rule of Section 3.
- **Correction to Amendment 2.** Its last sentence says that no rule was changed after the test run of the analysis script. That is
  not accurate. After that run the labelling code was changed so that U_val is marked as a lower bound for every check, which is
  what Section 5.2 states (the code had marked only quantities with an incomplete u_D). Amendment 2 also lowered two sample sizes:
  the V7 fully catalytic joint sample has 97 runs and the V8-1 fully catalytic input-only sample 49.
- **Censored glass-loss times.** Section 4 says that a V9 sample without glass loss by 480 s is censored and counted as such, but
  gives no rule for E and U_val. Where the baseline keeps its glass, E is reported as a lower bound, E > ln(time left/D), and the
  rule of Section 3 is not applied. U_val is not reported there, because the group samples that lose their glass are a biased subset.
- **Time limit for a reproduction.** `run_vvu.py` stops a run after 30 min of wall time and records it as Amendment 2 did. This
  affects only a rerun; the committed runs are unchanged.
- **Text.** Besides the places of Section 9, the paper has one sentence in the V6 part of Section 3, which the author's prompt
  allows. Outcome quantities carry the label of Section 3 in the tables.
