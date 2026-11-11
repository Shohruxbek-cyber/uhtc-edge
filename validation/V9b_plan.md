# V9b plan: the plasmatron check V9 with glass flow at the stagnation point (version 3, Task 1)

Written and committed on 10 Nov 2026 before any V9b run. Nothing below may change after this commit; a dated `## Amendment`
section at the end records any change. V9b is a new check, not a replacement: the V9 result (`validation/V9_plan.md`, Table 2)
stands in the record as it is.

## Decision: is a flat-face flow analogue defensible?

The author's prompt asks whether a radial-flow analogue on the flat faces of the Marschall et al. samples is defensible
(stagnation-region shear on a disc), and to run V9b only if it is. The author asked that this decision be made and documented here
rather than wait for a reading.

**Decision: yes, for the stagnation point, and only there.** Reasons:

1. *The quantity V9 scores is local.* V9 scores the glass-loss clock and the temperature at the stagnation point, which the
   two-colour radiometer views on the flat face away from the rounded edges (source, A1-7). Near the stagnation point of any blunt
   body, flat or spherical, the edge velocity grows linearly with the radius r, u_e = beta r, so the wall shear is tau = tau' r
   with tau' = Phi Pr^(2/3) beta (Reynolds analogy, D-FLOW-3) and the pressure gradient is G = -dp/dr = rho_e beta^2 r. With a
   film thickness that is smooth at r = 0, the film equation at the stagnation point reduces to a local sink,

       w_flow = tau' delta^2 / mu + 2 rho_e beta^2 delta^3 / (3 mu),

   which depends on the geometry only through beta. The 2D sphere-cone module computes the same sink at its stagnation node with the
   sphere's beta. A flat face changes beta, not the law.
2. *beta can be set from the test's own data.* Marschall et al. give the edge velocity V_e for each stage (Table 1, A1-11) but not
   beta. Two estimates agree within a factor of 1.6 (`validation/V9b_glassflow/beta_estimate.py`, run before this plan was written,
   `beta_estimate.json`): inverting the Fay-Riddell relation on the measured cold-wall flux at the rebuilt edge state gives
   beta = 1.1 x 10^4 s^-1 at the low stage and 1.6-1.8 x 10^4 s^-1 at the high stage; V_e divided by the 15 mm radius of the 30 mm
   probe gives 0.7 and 1.1-1.4 x 10^4 s^-1. V9b uses the Fay-Riddell value, because V9 converts the same measured flux.
3. *What it cannot do.* The analogue says nothing about the flow over the rest of the face, the rounded edge or the holder, about
   glass that might reach the centre from elsewhere (on a flat face none can), or about a temperature that varies over the face.
   It is a check of the stagnation-point clock with flow, as V9 is without it.

**Expected size, stated before the run.** With these beta values the wall shear 1 mm from the centre is 0.8-1.6 Pa, against about
1700 Pa 0.2 mm from the stagnation point of a 1 mm tip at the flight boundary (20 km, Mach 8.75). The flow sink at 2250 K with the LA
viscosity is 2-4 x 10^-5 um/s for a 50 um film, about four orders of magnitude below the thinning that V9 needs (tens of um in about
a minute). V9b is therefore expected to reproduce the V9 results to within the step tolerance. If it does, the omitted glass flow
cannot explain why the model keeps its glass too long in the plasmatron. That conclusion is the purpose of V9b, and it is reported
whichever way the run falls.

## Model

The 1D stagnation-line model and the inputs of V9 exactly (`validation/V9_plasmatron/run_v9.py`: samples 1.8, 1.5 and 1.6,
measured cold-wall flux converted as if the calorimeter were fully catalytic at 300 K, both catalysis bounds, no adjusted parameter),
plus the stagnation-point flow sink above after every oxide step: tau' from the step's conductance Phi and beta, G' = rho_e beta^2
with the equilibrium edge density at (H_e, p_s), viscosity LA at the surface temperature (D-FLOW-5), composition kept (D-FLOW-6),
criterion (a) rechecked after the sink. beta per stage from the Fay-Riddell inversion of the stage's measured flux with a
calorimeter wall at 350 K. The sink is applied by a wrapper of the oxide model in `validation/V9b_glassflow/run_v9b.py`; the model
code in `src/uhtc/` is not changed for V9b.

Sensitivity, reported but not scored: beta x0.5 and x2, the F320 viscosity.

## Pass criterion (the V9 criterion, unchanged)

Following V8, a case passes if at least one catalysis bound meets every band, and V9b passes if all three samples pass with the
same bound.
1. **Outcome.** Glass loss (criterion a) occurs in the high-power stage for samples 1.5 and 1.6, and no glass loss occurs by 480 s
   for sample 1.8.
2. **Timing.** For samples 1.5 and 1.6, the predicted time of glass loss after the switch is within a factor of 3 of the measured
   onset (68 s and 21 s).
3. **Temperature.** The low-stage plateau, and the high-stage temperature at the measured onset time (1.5, 1.6) or at 480 s (1.8),
   are within +-80 K of the measured values.

## Reporting

`validation/V9b_glassflow/report.md` and `v9b_results.json`, with the V9 values beside the V9b values and the sink rate at the
glass-loss time. If V9b fails, it fails. Its result goes into the paper wherever the glass-flow results are reported, with V9 left
as it stands.
