# V12 plan: comparison protocol of the three flank-heating models (version 3, Task 3)

Written and committed on 10 Nov 2026 before any V12 run and before the integral model is used in any map. Nothing below may change
after this commit; a dated `## Amendment` section at the end records any change.

## The three models

- **baseline**: the Lees / Kemp-Rose-Detra local-similarity distribution on the flank (`conduction2d.lees_distribution`, D-2D-7),
  edge state from the normal-shock stagnation point by isentropic expansion to the modified-Newtonian pressure;
- **sharpcone**: the laminar sharp-cone correlation of Tauber 1989 (eq. 46) on the flank with the running length from the
  stagnation point (`heat_dist = "tauber"`, D-2D-12), whose edge state is that of a sharp cone: the limit in which the boundary layer
  has swallowed the whole entropy layer;
- **integral** (new): the laminar integral method of Zoby et al. as printed by Riley 1992 (NASA TM-107630, Eqs. 3.3.1-3.3.2),
  momentum thickness marched along s with Eckert's reference enthalpy (Crabtree et al. 1965 Eqs. 1.3-1.4), edge state option (a)
  of the prompt (the same isentropic expansion as the baseline) (`conduction2d.integral_distribution`, D-HEAT-7).

Option (b) of the prompt, a swallowed entropy-layer edge state by mass balance, needs a shock shape. No shock-shape correlation was
found in a source read for this project by the time of this plan (Billig's correlation is known only through secondary quotations),
so option (b) is not run. The sharp-cone model is the fully swallowed limit and the two others the unswallowed one, so together they
bracket option (b).

## Deliverables and protocol (fixed now)

1. **Overlay.** g(s) = q(s)/q_s of the three models at the map corners: Mach 6, 9 and 12 at 20, 30 and 40 km, for R_n = 1 and 10 mm,
   5 deg sphere-cone, 30 mm flank, at the tangency point and at 1/4, 1/2 and the end of the flank. Table and figure.
2. **Measured anchor.** The 15 deg cone of V5 (Cleary 1969, NASA TN D-5450, M = 10.6, the digitised points that V5 uses): ratio of
   each model to the measured q/q_s at every measured point on the cone. Reported as the range of ratios per model.
3. **Columns.** The 15 band columns (fully catalytic, 20/25/30/35/40 km, R_n 1/3/10 mm) re-searched with the integral model
   (2D runner stage `flankint`), reported against the baseline boundary as the existing flank cases are (`flank` stage).
4. **Recommendation on the x0.5-2 flank band (rule fixed now).** Let L and H be the smallest and largest ratio, over the 18 corner
   points of item 1 on the flank, of each alternative model to the baseline, and of the measured anchor of item 2 to the baseline.
   If L >= 0.67 and H <= 1.5 for both alternatives and for the anchor, the band can be narrowed to [L, H], each end rounded outward to
   0.05; the Limitations text then states the narrower band and that it rests on two engineering models and one cone at one Mach
   number. Otherwise the band cannot be narrowed, and the text says why (which model or point falls outside).

## What V12 supports and what it does not

It compares engineering models with each other and with one measured cone. It does not test the flank heating of a 5 deg cone in
flight, transition, real-gas effects on the flank or the entropy-swallowed edge state itself, and it does not make the flank heating
validated. The band of the maps changes only through item 4, and only for the Limitations text and a sensitivity case; the baseline
maps keep the x0.5-2 band as computed.

## Outputs

`validation/V12_flank/run_v12.py`, `v12_results.json`, `report.md`, figure `figures/v3_flank_models.png`.
