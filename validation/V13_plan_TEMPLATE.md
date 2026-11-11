# V13 plan (TEMPLATE): a pre-registered check of the glass-loss clock for version 2

Renamed from `V10_plan_TEMPLATE.md` on 10 Nov 2026: the author's version-3 prompt assigns V10 to the verification of the
glass-flow module, V11 to the thermoelastic module and V12 to the flank-heating comparison, so this data check is V13.

**Fill in and commit before any run. Nothing may change after the run.** If a step proves impossible, add a dated
`## Amendment` section at the end and leave the text above it as committed. No model parameter is adjusted to this test.

Why a new check: once the oxide model changes because of V9 (for example by the glass-flow term of
`docs/v2_glass_flow_design.md`), V9 is no longer an independent check of that model. Version 2 needs a check that played no
part in its design. This template follows the structure of `validation/V9_plan.md`.

## Why this check

_Which part of the model it tests and why no earlier check tests it independently._

## Source

_Full reference, the file in `literature/pdf/` and `literature/text/`, the pages, tables and figures used, and whether the
source was read in full, as an edited manuscript, or only as an abstract (abstract-only sources are not admissible)._

## Test case

_Material, facility, gas, static and stagnation pressure, cold-wall heat flux and how it was measured (calorimeter type and
material, catalytic efficiency used by the source), enthalpy and how it was obtained, sample geometry and mounting, timeline
of the power stages._

| Sample | Stage | q_cw (W/cm²) | H (MJ/kg) | p (Pa) | Duration (s) | Measured outcome |
|---|---|---|---|---|---|---|
| | | | | | | |

## Measured targets, with reading uncertainty

_Every target, how it was read (table value or figure reading), the stated instrument uncertainty with its location, and the
reading uncertainty. Where the source states none, write "not stated"._

| Sample | Quantity | Value | Stated uncertainty (location) | Reading uncertainty |
|---|---|---|---|---|
| | | | | |

## Inputs to the model

_Model version (commit), conditions object, conversion of the measured flux (state the calorimeter efficiency used and why),
pressures, material, slab depth or geometry, start state, end time, catalysis bounds, oxygen bookkeeping._

## Pass criterion

_For each quantity the band and the rule that combines them (for example: a case passes if at least one catalysis bound meets
every band, and the check passes if all samples pass with the same bound). State also the validation-uncertainty analysis
that will be run (ASME V&V 20, as in `validation/VVU/VVU_plan.md`), with its uncertainty sources fixed now._

## Reporting

The result goes into Table 2, Section 3 and every place that summarises validation (Abstract, Limitations, Conclusions)
whether it passes or fails. If V13 fails, it fails: the criterion stays as written here.

## Candidate test data (sources opened for this template, 10 Nov 2026)

| Source | What it gives | Status | Note |
|---|---|---|---|
| Marschall et al., Appendix 1 of DTIC ADA564348, sample 1.7 (ZrB₂–30SiC, VKI Plasmatron) | Table 1: low stage 108 W/cm², 12.2 MJ/kg, P_dyn 25 Pa; high stage 219 W/cm², 20.3 MJ/kg, P_dyn 80 Pa; static 10 kPa; jump reported (A1-12) | read (edited manuscript) | Jump timing only on Fig. 3, coarsely; same campaign and calorimeter rebuild as V9, so not independent of V9's test conditions |
| Mungiguerra et al. 2019, IAC paper on SiC content in ZrB₂ (SPES arc jet), `literature/text/Mungiguerra2019_IAC_SiC_concentration_ZrB2.txt` | ZrB₂ with 5–18 vol% SiC, enthalpies up to 20 MJ/kg, two-colour pyrometers; spontaneous jumps of about 400 K on the lower-SiC compositions, surface above 2800 K after the jump | read (text file) | Heat flux not measured (SPES, as V7); jump onset times to be read from its figures |
| Mungiguerra et al. 2019, arc-jet tests of UHTC matrix composites, `literature/text/Mungiguerra2019_arcjet_UHTCMC.txt` | Fibre-reinforced ZrB₂–SiC composites, stepwise enthalpy to about 20 MJ/kg, two-colour pyrometer at 100 Hz | read (text file) | Different material (composite); usable only with a model of that material |
| Hu et al. 2013, HIT plasma wind tunnel (`literature/text/Hu2013_ZrB2SiC_HFPWT.txt`) | Jumps of 300–500 K on ZrB₂–10SiC at 4.78 MW/m², 27.9 MJ/kg, 18 kPa | read (Europe PMC text) | Same source as V8-2, so not independent of V8 |
