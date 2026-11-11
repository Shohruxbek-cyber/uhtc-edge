# DECISIONS

Every modelling or process decision, with its reason. IDs are cited in the code and the paper.
Newest entries at the end of each section.

## Index (updated 11 Nov 2026, version 3)

| Section | IDs |
|---|---|
| General | D-GEN-1 to D-GEN-9 |
| Atmosphere and flow | D-ATM-1, D-FLOW-1, D-FLOW-2 |
| Heating | D-HEAT-1 to D-HEAT-4 |
| Thermal model | D-TH-1 to D-TH-3, D-NUM-1 to D-NUM-3 |
| Quasi-2D model | D-2D-1 to D-2D-12 |
| Oxidation | D-OX-1 to D-OX-9 |
| Failure, runs and trajectory | D-FAIL-1 to D-FAIL-5, D-OX-10 to D-OX-14, D-TH-4, D-TRAJ-1 |
| Glass flow (version 3, Task 1) | D-FLOW-3 to D-FLOW-7 |
| Thermoelastic module (version 3, Task 2) | D-MECH-1 to D-MECH-4 |
| Flank heating, third model (version 3, Task 3) | D-HEAT-7 |

D-HEAT-5 and D-HEAT-6 are not used. Rule 0.2 of the version-3 prompt asks for the next free ID of a family, which was D-HEAT-5; the
integral flank model was given D-HEAT-7 before the family was checked, and it keeps that ID because the committed plan of V12
(`validation/V12_plan.md`), which may not change, cites it. D-FLOW-3 to D-FLOW-6 are the prompt's D-FLOW-1 to D-FLOW-4, renumbered
because D-FLOW-1 and D-FLOW-2 already existed (section Glass flow).

## General

- **D-GEN-1 Environment.** Pure Python on the existing interpreter (Python 3.14.6, numpy, scipy, matplotlib,
  numba, pytest, pymupdf already installed). `pandas` is not installed; rather than change the author's Python
  environment, the code uses the `csv` module and numpy instead (the task allowed pandas, it did not require it).
- **D-GEN-2 PDF build without a TeX engine.** No LaTeX engine is installed (checked 2 Nov 2026: no pdflatex,
  xelatex, latexmk or tectonic). Installing one would mean downloading and running an external binary, which
  the standing rules avoid. The paper's master source is LaTeX (`paper/main.tex`, compiles with pdflatex or
  Overleaf); the PDF in this repository is rendered from the same `.tex` by a small Python converter (a LaTeX
  subset → HTML → PDF with PyMuPDF; equations drawn with matplotlib mathtext), adapted from the author's
  `ife_wall/scripts/report_builder.py`. The typeset pdflatex version will look better; the content is identical.
- **D-GEN-3 Literature search in parallel.** Phase 1 was split among four AI research helpers (oxidation
  kinetics; properties/emissivity/catalysis/thermochemistry; arc-jet validation data; prior work). Each saved the
  PDFs and their extracted text in `literature/`, and every number with page/table location and a verbatim
  snippet, so that each value used in `INPUTS.md` was checked against the saved text before use.
- **D-GEN-4 Constants.** Physical constants are CODATA 2018 exact values (R, sigma). Molar masses are computed in
  `thermo.py` from the CIAAW abridged standard atomic weights (INPUTS.md section 1).

- **D-GEN-5 Literature files not in git.** The open-access PDFs (≈240 MB) and their text extracts stay in
  `literature/pdf/` and `literature/text/` on disk but are not tracked: they are large, and several are free to read
  but not licensed for redistribution, which matters if the repository is ever published. The notes, digitized data,
  JANAF tables (public domain) and `literature/MANIFEST_SHA256.txt` (file, checksum) are tracked.
- **D-GEN-6 UTF-8 console.** Scripts print non-ASCII symbols; on Windows run with `PYTHONUTF8=1` (as in the fusion project).
- **D-GEN-7 TeX engine (revision, 2 Nov 2026; supersedes the "no TeX" part of D-GEN-2).** The revision task asked for
  real typesetting with a user-level install. Installed Tectonic 0.17.0 (single binary, no admin rights, no system
  settings changed): `tectonic-0.17.0-x86_64-pc-windows-msvc.zip` from the official GitHub release
  (github.com/tectonic-typesetting/tectonic), SHA-256 f61ce51f0b0ade1015b7de7ef368541c5424e9756ecbd0d7af97d6d48030845f,
  identical to the digest GitHub publishes for the asset; unpacked to a user-level tools folder. On first use it
  downloads TeX packages into its user cache a user cache folder. To remove: delete
  that folder and a user cache folder (nothing else was touched; it is not on PATH).
  The Python renderer `paper/build_pdf.py` stays as a fallback only.
- **D-GEN-8 Page footer (3 Nov 2026).** elsarticle's preprint footer reads "Preprint submitted to Elsevier", which
  would be false on the arXiv version: no journal has been chosen and nothing has been submitted. The footer is
  redefined to "Preprint" and the date. When the author submits to a journal, `\journal{...}` can be restored.
- **D-GEN-9 Overfull lines (3 Nov 2026).** Two aligned equations (heating bound, surface balance) had their two
  rows joined by a lost line break and ran 3--29 pt into the margin; the break was restored. The two glass
  permeability laws moved from running text into a short display, unchanged in value.

## Atmosphere and flow

- **D-ATM-1** US Standard Atmosphere 1976 coded from the published document (Tables 2–4, eqs. 18, 21, 51)
  with the errata bound in the NTRS PDF; verified against Table I at 20, 30, 40 km (5 significant figures).
  No third-party atmosphere package, so the provenance of every constant is visible.
- **D-FLOW-1 Stagnation pressure.** Computed with an equilibrium-air normal shock (Rankine–Hugoniot with
  JANAF-based equilibrium air: N2, O2, NO, N, O, Ar), then p_s = p_inf + rho_inf V^2 (1 − eps/2) (Tauber 1989).
  The ideal-gas (gamma = 1.4) pitot value is kept for comparison; they differ by 0.3–4 % over Mach 5–12.
- **D-FLOW-2 Total enthalpy and wall enthalpy.** H_s = h(T_inf) + V^2/2 (exact energy conservation), from the
  JANAF enthalpy on the formation basis. Wall enthalpy: for the fully catalytic bound, equilibrium air at
  (T_w, p_s) (the gas at a catalytic wall is in equilibrium); for the non-catalytic bound, undissociated
  (frozen) air at T_w. Checked 2 Nov 2026: at T_w = 2500 K and p = 10 kPa the equilibrium value exceeds the
  frozen one by 0.25 MJ/kg (about 13 % of H_s − h_w at Mach 10, 30 km), so the distinction matters and is kept.
  Both are tabulated once per case at p_s (10 K grid) and interpolated in the time loop.

## Heating

- **D-HEAT-1 Correlation.** Baseline stagnation heating: Sutton & Graves (1971) eq. (33) with the air
  coefficient K = 0.1113 kg s^-1 m^-3/2 atm^-1/2 (Table II), which already contains the hot-wall factor
  (H_s − h_w). Fay–Riddell (Hamilton 1982 eq. 6) is the cross-check (V4). Sutton–Graves was derived for wall
  temperatures of 300 and 1111 K and H_s ≥ 2.3 MJ/kg; using it at T_w up to ~2800 K and H_s down to 1 MJ/kg is
  an extrapolation, quantified by the V4 comparison with Fay–Riddell evaluated with hot-wall properties.
- **D-HEAT-2 Catalysis bounds.** Fully catalytic = Sutton–Graves (equilibrium boundary layer) with the wall
  enthalpy of equilibrium air at (T_w, p_s). Non-catalytic bound = the Fay–Riddell bracket evaluated with Le = 0,
  as Kolodziej (1997, NASA TM-112204, doc p. 3) does for his "noncatalytic hot wall" limit:
  q_nc = q_fc (1 − h_D/H_s)/[1 + (Le^0.52 − 1) h_D/H_s], with the frozen-air wall enthalpy. Physically: atoms
  reach the wall but do not recombine there. Le = 1.4 (ASSUMPTION, see INPUTS.md). Real ZrB2–SiC oxide surfaces
  lie between the bounds; the literature recombination coefficients are listed in INPUTS.md.
  Cap (added 2 Nov 2026 after a production check): this multiplicative form comes from cold-wall theory; near
  convective blocking (wall temperature approaching the edge temperature, e.g. Mach 8, 30 km, T_w ≈ 2330 K vs
  T_e = 2575 K) it gave a non-catalytic flux above the catalytic one, which is impossible (the equilibrium wall
  enthalpy already contains the dissociation the hot wall cannot recombine). The non-catalytic flux is therefore
  capped at the fully catalytic value; near blocking the two bounds coincide. The cap is inactive at the V3
  comparison points (checked: V3 results unchanged).
- **D-HEAT-4 Catalysis tied to the oxidation state (variant).** Besides the two bounds, a variant is run in which the
  wall is non-catalytic while an outer glass exists and fully catalytic after it is lost. Reason: the novelty review
  found this coupling untreated in lifetime studies, while measured recombination coefficients are low on glass
  (≤ 0.02, RCG fits) and higher on oxidized ZrB2–SiC / zirconia (≈0.08 at 1800 K). It is a bounding variant, not a
  calibrated model. Implementation: the wall is catalytic when a glass has formed once and is now thinner than
  0.1 µm (d_fail); criterion (a) additionally needs 1 µm of converted composite. In the hottest map cells (failure
  within 6 s) the first thin glass briefly falls below 0.1 µm before that, so the switched wall turns catalytic early
  and the time to first failure differs from the non-catalytic bound by up to 1.5 % (found by the Phase 5 re-trace;
  reported in the paper, not changed).
- **D-HEAT-3 Mass-transfer coefficient.** The boundary-layer mass-transfer conductance for evaporation is taken
  equal to the heat-transfer conductance rho_e u_e C_H = q/(H_s − h_w) = K sqrt(p_s/R_n) (Reynolds/Chilton–Colburn
  analogy with Le = 1; ASSUMPTION with sensitivity ×0.5/×2). Revision 1, task 2.3: Eq. 2 uses Le = 1.4, so the two are
  not consistent. The correction factor is Le^0.52 = 1.19 (the Fay–Riddell exponent of Eq. 2) to Le^(2/3) = 1.25
  (Chilton–Colburn), inside the ×0.5–×2 band. Baseline kept pending the author's decision; run as
  `runs/sensitivity/run_lewis.py` (B1 711 → 548 and 510 s, B2 unchanged).

## Thermal model

- **D-TH-1 1D solver.** Adapted from the author's fusion-wall enthalpy solver (backward Euler, Newton,
  Lagrangian mass coordinate), see `src/uhtc/conduction1d.py` header for what was kept and changed.
- **D-TH-2 1D back boundary.** Planar slab 20 mm deep with an adiabatic back face. Reason: for exposures longer than
  about (20 mm)²/α ≈ 20 s the slab is heat-soaked and the surface tends to radiative equilibrium, which is the
  conservative (hotter) limit compared with a semi-infinite solid that keeps absorbing heat; the real sink for a sharp
  nose — conduction aft through a narrowing section — is what the quasi-2D model of Q5 resolves. V2 shows that
  the depth changes transient failure times by up to a factor of 2 (20 vs 50 mm) and steady states not at all, so the
  depth is a modelling choice reported with the results, not a numerical parameter.
- **D-TH-3 Oxide scale as a quasi-steady resistance.** The scale (glass + ZrO2 + depleted layer, at most a few
  hundred µm) is represented by its thermal resistance between the outer surface and the first solid node;
  its heat capacity is neglected. Reason: the volumetric heat capacities of ZrO2 and of the substrate are
  similar, and after 1–10 s the heated depth of the substrate (sqrt(alpha t), alpha ≈ 2e-5 m²/s) is 4–14 mm, so
  a 100 µm scale stores about 1–3 % of the heat; its resistance (k of order 1–2 W/m/K vs tens of W/m/K for the
  substrate) changes the surface temperature by tens to hundreds of K at MW/m² fluxes and must not be dropped.

- **D-NUM-1 Coupling and step control.** Each step solves conduction + surface balance implicitly with the oxide state
  frozen, then advances the oxide state with the new temperatures (operator splitting). Steps are limited to a surface
  temperature change of 10 K, 5 % change of the converted depth (evaluated 50 K above the current front temperature)
  and 5 % of the glass inventory (but at least 0.05 µm of glass). Checked by halving the limits (V2).
- **D-NUM-2 Speed-up for the 2D production runs (revision).** Profiling a 2D run showed 70 % of the time in the
  surface-chemistry term (evaluated by forward differences inside every surface Newton iteration of every node) and in
  `np.interp` calls on the uniform wall-enthalpy grid. The chemistry term is now compiled with numba (same arithmetic;
  the Python version is kept as `LayerOxidation.chem_reference`) and the wall-enthalpy lookup uses index arithmetic on
  its uniform 10 K grid. Unit tests compare both with the reference paths; a 1D base case (B1, 710.862… s) and a 2D case
  (Mach 8.5, 30 km, T_max 1845.2039… K) reproduce the saved results bit for bit. A 2D run is about 2.4× faster. A banded
  linear solver would save a further ~20 % but changes round-off and was not adopted.
- **D-NUM-3 2D production settings checked (revision).** `runs/envelope_2d/check_mesh.py`: production mesh (30 arc,
  40 flank, 30 normal cells, first cell 0.003 R-fraction) with a 2 s step cap against a mesh twice as fine in each
  direction (first cell halved) with a 0.5 s cap. Criterion set before the run: peak stagnation temperature within
  0.5 % and failure time within 2 %. Survivors (600 s, three radii): temperatures within 0.01 %. Failing cells (1800 s):
  temperatures within 0.02 %; glass-loss times within 0.8 %; one temperature-limit time differs by 3.5 %
  (48.1 vs 49.8 s) — the plateau there lies only 0.4 K above the limit, so the crossing time is ill-conditioned
  (a 0.4 K change of the plateau moves it by seconds), not inaccurate. Recorded as a FAIL of the time criterion for
  that case and reported; such near-tangent crossings are flagged in the maps by their peak temperature. A 5 s cap
  was as accurate but saves only 30 % on a 30-minute run (36 vs 51 s); the verified 2 s cap is kept. The full 2D grid
  (21 altitudes) is therefore affordable (about 6 h on 4 processes) and is run instead of the reduced 5-altitude grid.
  The checks above used the draft distribution; repeated with the sourced baseline (D-2D-7; `check_mesh_lees*.jsonl`)
  they pass: survivors within 0.27 % in peak temperature, failing cells within 0.03 % in temperature and 1.02 % in
  failure time (all three failing cells are glass-loss cases).

## Quasi-2D model

- **D-2D-1 Heating distribution.** cos(phi) on the nose arc (Tauber 1989 p. 8, Figs. 5–6: follows Lees' and
  Kemp et al.'s solutions to about 70°); on the flank the laminar sharp-cone (eq. 46) or wedge (eq. 48)
  correlation with the running length taken as the wetted distance from the stagnation point (ASSUMPTION; at
  the tangency point it matches the cosine law within about 15–20 % for theta_c = 5–10°; sensitivity ×0.5/×2
  on the flank heating).
- **D-2D-2 Back boundary.** Adiabatic plane normal to the flank at the end of the modelled flank length
  (conservative: no heat sink behind the tip region).

- **D-2D-3 Surface model in 2D.** Every surface node carries its own oxide state and solves the same surface
  balance as the 1D model; local heating and mass transfer follow g(s); the O2 partial pressure of the stagnation
  point is used on the whole surface (no pressure distribution is modelled; it matters only away from the tip, while
  Q5 compares the stagnation point). The failure clocks are those of the stagnation node. The 2D Newton tolerance is
  1e-4 K (1e-6 K in 1D) because the per-node chemistry derivative is a finite difference.
- **D-2D-4 Q5 geometry.** Half-angle 5° (slender nose / wing edge; 10° as a sensitivity case for R_n = 1 mm), flank
  length 30 mm, adiabatic back plane; nose tip = axisymmetric sphere-cone heated with the sphere stagnation value,
  wing edge = planar cylinder-wedge heated with the unswept-cylinder value (Tauber eq. 41 / eq. 40 ratio 0.705). The 1D
  reference for each uses the same stagnation heating.
- **D-2D-5 Flank sensitivity for Q5 (added in Phase 4).** The first Q5 runs showed the 2D stagnation temperature to be
  almost the same for R_n = 1 and 3 mm (about 1530 K at Mach 7, 25 km), i.e. the slender tip acts as one radiating fin
  and the flank assumptions of D-2D-1/D-2D-4 can matter. The planned x0.5/x2 flank-heating case (D-2D-1) was therefore
  run, plus flank lengths of 15 and 60 mm, for the 1 mm nose tip at Mach 8.5 and 9, 30 km (`runs/q5_2d/run_q5_sens.py`;
  `run_2d(..., flank_mult=)` scales the flank heating only, default 1 leaves all other runs unchanged).
- **D-2D-6 Coarse 2D survival map (added in Phase 4).** Q5 at fixed conditions showed 350–710 K lower stagnation
  temperatures in 2D and almost no failures within 600 s, so the answer to "how much does 2D extend survival" is better
  given as a shift of the 30-minute boundary. A full 2D map is too expensive (each 1800 s 2D run takes minutes), so a
  coarse one is run: nose tip, R_n = 1 and 3 mm, fully catalytic, Mach 9–12 in steps of 1, at 20, 30 and 40 km
  (`runs/q5_2d/run_q5_envelope.py`). Its Mach resolution (1) is coarser than the 1D maps (0.5); the boundary is
  reported as a bracket. Superseded in the revision by the full 2D maps (D-2D-9).
- **D-2D-7 Sourced heating distribution (revision; replaces the draft D-2D-1 as baseline).** Sources (open):
  `literature/notes/HeatingDistribution_formulas.md`, `_data.md`, `_baseline.md`. The draft flank model (sharp-cone
  correlation with running length = wetted distance) is not supported for blunted cones: Cleary 1969 (NASA TN D-5450,
  15 deg sphere-cone, M 10.6) measured 0.57–0.65 of it at s/R 1.9–5.6, and against the Crabtree et al. 1965 charts
  for 5 deg it was read as 1.1–3.8 times higher at M 8–14 (chart comparison, script not kept; recomputed in revision 1
  against the baseline distribution as 0.8–4.1 over the flank, `audit/tools/flank_ratio.py`) (Vallerani, via Tauber p. 9: blunting reduces cone heating for
  5–25 deg). Baseline now (`uhtc.conduction2d.baseline_distribution`): cos(phi) on the arc up to 70 deg (Tauber 1989
  p. 8: supported to ~70 deg); from the tangency on, the Lees / Kemp-Rose-Detra local-similarity result in one integral
  over arc and flank (Cleary eqs. 2, 4, 5; DeJarnette & Tai 1969 eqs. 10–11), modified-Newtonian pressure, isentropic
  surface flow, pressure-gradient factor of Crabtree eq. 3.5 (`lees_distribution`); between 70 deg and the tangency
  (sources disagree, 0.25–2 x cos) a linear bridge in phi. Checks: the implementation reproduces the Crabtree Fig. 12
  hemisphere-cone charts (5/10/15 deg, M 5–14) within 10 % (unit test), and Cleary's measured 15 deg data lie at
  0.72–1.15 of the baseline (V5). Near the stagnation point the Lees-Newtonian curve falls below cos(phi) by up to
  22 % at 50–70 deg; cos is kept there because Tauber's Fig. 5 data follow it.
  **Band (Part 1.3):** beyond 70 deg (shoulder and flank) x0.66 (low) to x1.5 (high), the high case also at least the
  Kemp-Rose-Detra shoulder curve of Tauber Fig. 5 (DIGITIZED, 0.25–0.16 at 75–90 deg). Reason: measured/Crabtree theory
  for the 15 deg cone spans 0.66–1.05 (Cleary Table III); at 5 deg the sign of the Newtonian-pressure error is unknown
  (the helper's analysis), so the low bound is mirrored (1/0.66 = 1.5). This replaces the arbitrary x0.5/x2 of D-2D-5.
  The draft distribution ("tauber") is kept as a comparison case.
- **D-2D-8 Wing edge (planar) not a 2D headline result (revision).** For the planar cylinder-wedge the same open
  method with Newtonian pressure gives 0.09–0.36 of the flank heating that Kolodziej et al. 1997 (Fig. 2) computed for a
  5 deg wedge with the BLIMPK boundary-layer code, and 0.14–0.46 of their PERFORM result (V5 part 2): a blunt planar body
  keeps a much higher pressure downstream than Newtonian theory gives. The draft "tauber" wedge flank is also low
  (1.5–7x, helper's analysis). No open planar distribution consistent with BLIMPK for our Mach range was found, so 2D
  wing-edge results would be non-conservative; the revision reports 2D results for the axisymmetric nose tip only and
  keeps the 1D cylinder stagnation line as the (conservative) wing-edge estimate.
- **D-2D-9 2D production maps (revision).** Nose tip, full 1D grid (Mach 5–12 step 0.5, 20–40 km step 1 km, R_n 1/3/10
  mm, both catalysis bounds, 20 vol% SiC), baseline distribution, production numerics (D-NUM-3), stopped at the first
  failure; half-step refinement at the 30-minute boundary; band runs (D-2D-7 low/high, emissivity low/high) and glass-
  criterion variants re-locate the boundary at 20/25/30/35/40 km for the fully catalytic bound
  (`runs/envelope_2d/run_envelope_2d.py`).
- **D-2D-10 Definition of the uncertainty band of the 30-minute boundary (revision, Part 1.3).** "Low/high heating"
  has two parts for a nose tip: the level at the stagnation point (Sutton–Graves constant ×0.85–1.15, INPUTS.md; the
  strongest single factor in the 1D sensitivity) and the distribution around the nose (D-2D-7 low/high). Both are run
  one at a time, as are the low/high emissivities (glass 0.70/0.95 with oxide 0.50/0.85). The reported band is the
  range spanned by the two combined corners — "worst" (high level, high distribution, low emissivities) and "best"
  (low level, low distribution, high emissivities) — so that it is not narrower than any single factor; the one-at-a-
  time cases show which factor dominates. Fully catalytic bound only (the conservative bound; the non-catalytic
  baseline map is reported alongside). A corner combines independent uncertainties at their extremes, so the band is
  wide by construction; it is a bound, not a confidence interval.
- **D-2D-11 Q5 flank sensitivity rerun with the baseline distribution (revision).** The flank sensitivity of D-2D-5
  (flank heating x0.5/x2, flank 15/60 mm) had been run with the draft distribution (`run_2d` default "tauber"), so it
  no longer described the baseline model. `runs/q5_2d/run_q5_sens.py` now passes the baseline (D-2D-7) explicitly and is
  rerun; the coarse 2D map of D-2D-6 and its table are dropped from the paper (the full maps of D-2D-9 replace them).
- **D-2D-12 Flank-heating dependence of the 1D to 2D gain (revision 1, task 2.7).** The largest effect of the paper rests on
  flank heating, tested until now only at Mach 8.5–9 and 30 km. Stage `flank` of `runs/envelope_2d/run_envelope_2d.py` relocates the
  boundary for ×0.5, ×2 and the sharp-cone running-length correlation at five altitudes, 1 and 3 mm, fully catalytic. These
  variants are wider than the ×0.66–1.5 band of D-2D-7, which stays the uncertainty band; they test the headline, not the band.

## Oxidation

- **D-OX-1 Oxygen at the wall.** p_O2,w = 0.2095 p_s: all oxygen counted as O2, i.e. atoms reaching the wall
  are assumed to recombine or to oxidise like O2. Atomic oxygen is known to shift the SiC active/passive
  transition (Balat 1996); this is a stated limitation with a sensitivity case.
- **D-OX-2 Layer bookkeeping.** The ZrO2 skeleton occupies the converted depth x_r (Shugart's measured "ZrO2 layer"
  is compared with x_r); its pore volume, 1 − (1 − f_SiC)(1 + 0.169), is filled with glass; the remaining glass from
  B2O3 (ZrB2) and SiO2 (SiC) joins the outer layer. Glass flow away from the stagnation point is neglected (zero shear
  at the stagnation point; stated limitation).
- **D-OX-3 Glass thermodynamics.** Ideal B2O3–SiO2 solution (activity = mole fraction) for the vapour pressures.
- **D-OX-4 Temperatures used.** Oxide growth and SiC depletion use the reaction-front temperature (first solid node
  below the scale); evaporation and the active-oxidation check use the outer-surface temperature. No change of the
  growth law after glass loss is modelled (no open data for unprotected ZrB2–SiC at these temperatures), so histories
  after the first failure are indicative only and are labelled so.
- **D-OX-5 Emissivity.** 0.88 while at least 1 µm of glass covers the surface, falling linearly to 0.75 as the glass
  disappears (literature ranges in INPUTS.md section 8); this reproduces in kind the emissivity drop reported with the
  ground-test temperature jumps.
- **D-OX-6 Oxygen-pressure dependence of the growth law.** n_p = 0: the high-temperature calibration (1 atm O2) is used
  unchanged at other O2 pressures. No usable evidence (10 Nov 2026: values known only from an abstract are not used). Earlier
  note, kept for the record: at 1700 °C Karlsdottir et al. (abstract) measured 102 µm in 15 min,
  more than the 65 µm the O2 constant gives. UNSURE: the abstract does not state the furnace atmosphere (assumed air
  when this was first written; corrected 2 Nov 2026) and the material is ZrB2–15SiC, so this is not evidence that the
  growth rate is pressure-independent. Sensitivity n_p = 0.5 (Q4: it roughly halves the glass-governed time to failure
  at 35 km, so it is a reported uncertainty, not a negligible one).
- **D-OX-7 SiC of the depleted layer.** Assumed to re-oxidise to SiO2 inside the scale (SiO(g) out-diffusing and
  oxidising; mechanism stated by Opila & Halbig 2001, p. 7; the Fahrenholtz 2007 abstract is no longer cited). Sensitivity: lost as gas (f_SiO_lost = 1). The depleted layer
  grows without incubation period: conservative for the depleted depth, but it adds SiO2 to the glass earlier than real
  material would (V6 shows the depleted depth over-predicted at short times at 1627–1700 °C), which makes glass survival
  somewhat optimistic; the f_SiO_lost = 1 case bounds this.
- **D-OX-8 Kinetic data set.** Baseline = Shugart 2014 (ZrB2–30SiC, layer-resolved ZrO2 growth in two regimes split at
  the SiC-depletion onset 1627 °C; depletion depth at 1800 °C) because it is the only open data set that separates the
  layers; below 1573 K the low-T law is continued with the Opila & Halbig activation energy (the Shugart low-T constants
  are nearly flat in T and must not be extrapolated down). The independent data (Opila & Halbig 2001, Levine & Opila
  2003, Opila et al. 2003) are kept for V6; the Karlsdottir et al. 2007 values, known only from its abstract, were dropped
  from V6 on 10 Nov 2026 (author decision: no validation data from an abstract alone). Above 2073 K the high-T Arrhenius law is an
  extrapolation (stated). The same constants are used for 10/20/30 vol% SiC (no open composition-resolved kinetics);
  the SiC fraction enters through glass production and properties. Sensitivity: kp ×1/3 and ×3 (batch-to-batch spread
  of about 3 at 1627 °C between Opila & Halbig 2001 and Opila et al. 2003).
- **D-OX-9 SiC active-oxidation criterion.** For a pre-oxidised surface that heats up, the relevant boundary is the
  passive-to-active one (scale breakdown), which lies 0.3–1 decade below the active-to-passive line in p_O2 (Jacobson
  et al. 2012); baseline 0.65 decade. The line is our fit to 45 literature points (1473–1924 K, rms 0.28 decade; revision 1 corrected one misread point and
  added two hidden ones, CHANGELOG.md) and is
  extrapolated to the 2300–2900 K reached here (stated). Sensitivity: offsets 0.3 and 1.0; the active-to-passive line
  itself (offset 0); atomic oxygen (+1.65 decades, from the Rosner & Allendorf points).

## Failure, runs and trajectory

- **D-FAIL-1 Failure criteria.** (a) no protective glass: > 1 µm of composite converted but < 0.1 µm of outer glass
  (lost after forming, or evaporating as fast as it forms); (b) SiC active oxidation: p_O2,w below the passive-to-active
  line at T_s; (c) T_s ≥ 2422 K (2149 °C, the reusable sharp-leading-edge limit used by Kolodziej 1997). Time to failure
  = the first of the three. The limit (c) is a design value, not a property; sensitivities 2289 K (TPSX multiple-use
  limit for ZrB2/SiC) and 2543 K (the eutectic, D-FAIL-5; 3033 K until revision 1). All three times are stored for every run.
- **D-FAIL-5 Upper end of the temperature-limit range (revision 1, task 2.1).** 3033 K (2760 °C single use) is above the
  ZrB2–SiC pseudo-binary eutectic, measured at 2543 ± 30 K (Ordan'yan et al. 1989, as reported by Tu et al. 2008) and 2570 K
  (Tu et al. 2008). First liquid forms there for any composition, so a higher limit is not physical for ZrB2–20SiC. The range
  now ends at 2543 K, the lower measurement, with no margin subtracted because it is a ceiling and not a design value. A
  2207 °C value quoted by one secondary paper for the same source is unresolved (`literature/notes/ZrB2-SiC_eutectic.md`).
- **D-FAIL-2 Exposure cap.** Constant-condition runs stop at the first failure or at 1800 s (30 min, long compared with
  hypersonic cruise segments of minutes); "≥ 1800 s" is reported as survival.
- **D-FAIL-3 Glass-loss threshold kept at 0.1 µm; its sensitivity reported (added in the final review).** The planned
  0.01–1 µm sensitivity (INPUTS section 8) had not been run; `runs/sensitivity/run_extra.py` ran it with the depleted-
  layer conductivity (10/60 W/m/K) and the atomic-oxygen offset (+1.65 decades). 0.01 µm, k_dep and atomic O change
  nothing (< 2 %). 1 µm turns both base cases into early glass-loss failures (B1 24 s, B2 13 s): below 1900 K the oxide
  grows slowly and the glass is only 0.1–0.5 µm thick for about 10–28 s of the heat-up. UNSURE: whether a film this
  thin protects the surface is not known; the baseline treats it as protective (it is a continuous film in the model).
  Reported in the paper (Q4, limitations) with a [CHECK] marker. The Lewis-number case (1.4 → 1.0) was not run: it moves
  the non-catalytic bound towards the catalytic one, so its effect is bounded by the catalysis-bound row.
  Superseded as baseline by D-FAIL-4; the fixed thresholds remain sensitivity cases.
- **D-FAIL-4 Physics-based criterion (a) (revision; primary).** The glass protects as long as it limits the O2 supply
  to the reaction front. The largest O2 flux the boundary layer can deliver is J_gas = g ln(1 + Y_O2)/M_O2 (wall O2
  mass fraction taken as 0; g = the mass-transfer conductance of D-HEAT-3; Y_O2 = 0.2315, all O counted as O2,
  D-OX-1); the film carries at most K(T, x_B) p_O2,w / d (Henry's law). Protection is lost when d < d_c = K p_O2,w /
  J_gas (and more than 1 µm of composite has been converted, as before). K: molecular O2 permeability, linear in the
  B2O3 mole fraction between silica and B2O3 (data.O2_PERM; sources in INPUTS section 8 and
  `literature/notes/GlassOxygenTransport_*.md`): silica base = implied by the CVD-SiC parabolic constant (Jacobson et al.
  2005 Table II, citing Ogbuji & Opila 1995; 1473–1773 K), low = Kajihara et al. 2004 (D·S, 1073–1373 K), high =
  Deal–Grove Si (1073–1473 K); B2O3 base = O self-diffusion in the melt (Tokuda via Schlichting 1982, 942–1010 K) times
  the O2-equivalent density, high = Luthra's value as printed by Campbell et al. 2021 (UNSURE reading). All are
  extrapolated to 1300–2450 K. Result: d_c is sub-nanometre for silica-rich glass and at most ~4 nm (base; ~12 nm high)
  for the B2O3-rich glass of the heat-up (1470–1550 K), when the glass is 0.25–0.3 µm thick; so the criterion fires only
  when the glass has essentially gone, and never during the heat-up (B1: 711.3 s vs 710.9 s with 0.1 µm; B2 unchanged).
  Fixed thresholds 0.1/0.5/1 µm are kept as sensitivity cases. Limitation: below ~30 nm a silica film need not follow
  the parabolic (permeation) law (Deal–Grove thin-oxide anomaly, Hollauer 2007 §2.6.6); continuity and bubbling of a
  real film are not modelled. The catalysis switch (D-HEAT-4) now uses the same test, so switch and criterion (a) agree.
- **D-OX-10 Incubation of SiC depletion (revision; model option, off in the baseline).** Shugart 2014 Fig. 40 (the only
  open source with onset times): onset brackets 1627 °C ≤ 40 min, 1650 °C 10–20 min, 1700 °C 10–15 min, 1800 °C ≤ 5 min
  (ZrB2-30SiC, flowing O2). Option: depletion starts when the incubation clock sum(dt/t_inc(T_f)) reaches 1 (additivity
  rule, our assumption), t_inc = bracket midpoints 20/15/12.5/2.5 min (our construction), held constant outside. Not
  fitted to any thickness; but it comes from the same experiments as V6's depleted-layer points, so V6 with the option
  is not an independent test of those points. Kept off in the baseline because the onset data exist only for flowing
  1-atm O2 and one material; reported as a sensitivity.
- **D-OX-11 Extrapolation band of the high-temperature growth law (revision, Part 3.3).** The fit's own confidence
  interval is useless (three points: the 95 % interval of E spans negative values), so the band comes from the open
  alternatives (`literature/notes/HighT_kinetics_alternatives.md`, `digitized/HighT_oxidation_points.csv`): activation
  energies on a k_p basis of ~250 kJ/mol (Oguri et al. 2011a Fig. 3, 1700–1900 °C, DERIVED) to 354 kJ/mol (same data,
  1800–1900 °C) — the baseline 326 lies inside; rate level: furnace data in air at ≥ 1900 K lie 1.5–5× above the law,
  the one controlled plasma series (Hu et al. 2013, ZS30 at 2233 K) grows about half as fast (125 vs 177 µm at 600 s).
  Band cases above the fit end (2073 K): E = 250 and 354 kJ/mol anchored at the fit's value at 2073 K, and k_p ×0.5 and
  ×3 (`OxParams.kp_E_hi`, `kp_mult_hi`); run on the glass-loss boundary cells. The ×3 upper multiplier follows the
  helper's proposal (within the 1.5–5× furnace spread); our choice.
- **D-OX-12 Alternative passive-to-active line (revision, Part 4, sensitivity only).** The [CHECK] search found a
  measured transition law for SiC-based composites in air plasma (Balat-Pichelin et al. 2015, HAL hal-02313822, PDF
  p. 4–5, fitted on solar-furnace and Plasmatron points at 200–5000 Pa total air): p_O2 [Pa] = 10¹⁵ exp(−55850 K/T)
  (`data.PA_BALAT_A`, `PA_BALAT_B`; `OxParams(pa_law="balat2015")`). At flight O2 pressures it puts the transition
  280–870 K below our line (Eq. pa: 2360–3380 K for 2–220 kPa O2; this law: ≈2070–2510 K, extrapolated 0.3–2.3
  decades in pressure). It is for C/SiC coatings, not ZrB2–SiC, and plasma tests include atomic oxygen, so it is not
  adopted as the baseline (the baseline keeps the molecular-O2 compilation, D-OX-9); it is run as a sensitivity case on
  the 1D maps (set D) and on the 2D grid ("pa_balat" stage), because it could make active oxidation, not glass loss,
  the governing mode near the boundary (example: at 10 kPa O2 the law gives 2205 K, below the 2320 K surface of the
  draft's case B1).
- **D-OX-14 Effective oxidant pressure with atomic oxygen (revision 1, task 2.4; sensitivity only).** p_eff =
  max(0.2095 p_s, (x_O2 + x_O) p_s) with the equilibrium edge composition, every O atom counted as one O2 molecule; 2D stage
  `oxidant`, the 15 fully catalytic band columns. The baseline stays 0.2095 p_s: a fully catalytic wall recombines the atoms.
- **D-OX-13 Where the growth law is extrapolated in 2D (revision, Part 3.3 follow-up).** The growth law is evaluated
  at the reaction-front temperature T_f (D-OX-4), not at the surface. In 1D with an adiabatic back the steady conducted
  flux is near zero, so T_f is close to T_s; in the 2D tip heat flows continuously through the scale into the flank,
  so T_f lies below T_s by the drop across the scale. The 2D kinetics band moved the boundary by at most 0.25 Mach in
  3 of 15 columns although every boundary point has T_s > 2073 K, which prompted a check: the grid records only T_s,
  so `runs/envelope_2d/run_front.py` reruns the highest surviving and lowest failing case of every bracketed column
  with unchanged inputs (it must reproduce the stored T_s_max and failure time) and records the largest T_f and the
  time with T_f > 2073 K. The map line "surface above 2073 K" is kept (it is what the grid records) and the paper
  states where the front, which is what the growth law sees, exceeds 2073 K.
  Outcome (4 Nov 2026): all 226 reruns reproduce the stored results exactly. The hottest front lies 160–310 K below
  the hottest surface at the boundary. 1 mm: front ≤ 2060 K at the boundary (failing points: > 2073 K for ≤ 103 s
  in the heat-up). 3 mm: > 2073 K in 19 of 38 columns, for ≤ 415 s. 10 mm: > 2073 K in every column for 630–1770 s.
  So the law is extrapolated in earnest at 10 mm, and the small boundary shift there (≤ 0.25 Mach, the search
  resolution) is a result, not an artefact of a cool front. Fig. 7 marks the columns whose highest surviving
  condition has T_f > 2073 K with a dot; the paper does not claim that the front stays cool in general.
- **D-TH-4 Conductivity factor for 30 vol% SiC (revision).** The draft's 0.75 divided a 30 vol% value from one study
  (Zimmermann 2008 abstract) by the ManLabs 20 vol% table, i.e. it mixed material and composition effects. Replaced by
  0.92, the within-study ratio of Van Laningham 2012 (Fig. 27, interpolated between 21.9 and 48.7 vol%: 0.94–0.91 at
  1200–2000 °C), the same source and method as the 1.17 used for 10 vol%; 0.75 kept as a low bound. Kim et al. 2012
  (SPS, ≤ 1000 °C) show no drop with SiC fraction.
- **D-TRAJ-1 X-43A trajectory.** The 88.8–111.8 s gap between the boost and descent figures (separation at 88.35 s,
  Table 3; engine test near Mach 9.7) is filled by linear interpolation; only the Mach ≥ 3 part (33.5–637.5 s) is
  simulated because the correlations are for hypersonic flow; the edge starts at 300 K at Mach 3 (the subsonic carry
  and early boost add little heat). The X-43A nose and leading edges were carbon-composite based (Peters et al. 2024,
  Fig. 1 caption: carbon composites, tungsten alloy, Ir coating); here the trajectory is used as
  a published, reproducible flight path for a hypothetical ZrB2–SiC edge.

## Glass flow (version 3, Task 1)

The author's version-3 prompt names these decisions D-FLOW-1 to D-FLOW-4. D-FLOW-1 and D-FLOW-2 already name the stagnation
pressure and enthalpy (section "Atmosphere and flow"), so by the prompt's own rule 0.2 (next free ID in the family) they are
D-FLOW-3 to D-FLOW-6 here, and D-FLOW-7 is added for the numerics. Code: `src/uhtc/glassflow.py`, `run_2d(glass_flow=...)`.
Off in every run of the paper (regression: `runs/analysis/regression_v3.py`).

- **D-FLOW-3 Wall shear (prompt D-FLOW-1).** tau_w(s) = Phi_s g(s) Pr^(2/3) u_e(s): the laminar Reynolds analogy
  c_f/2 = St Pr^(2/3) (Crabtree et al. 1965, Section 1.3), the form of the Supplement S9 estimate extended to every surface node,
  with Pr = 0.71, the heating conductance Phi_s of Eq. 1 (heating multiplier included), the heating distribution g(s) of the map
  (cosine on the arc, flank correlation on the flank) and u_e the isentropic expansion (gamma = 1.4) from the stagnation point to
  the modified-Newtonian wall pressure (D-FLOW-4). Face values are the mean of the two nodes. Band x0.5-2 on tau_w (case key
  flow_tau_mult), the same band as the flank heating, because the same distribution g(s) carries both. 10 Nov 2026.
- **D-FLOW-4 Pressure gradient (prompt D-FLOW-2).** G = -dp/ds from the modified-Newtonian ratio
  p/p_s = (1 - eps) cos^2(phi) + eps, eps = p_inf/p_s, the ratio of the Lees distribution and of the local_pressure check,
  differenced between neighbouring nodes. It is zero on the flank, where the pressure is constant. Included in the baseline;
  a sensitivity case switches it off (case key flow_pressure = False) to show whether shear dominates. 10 Nov 2026.
- **D-FLOW-5 Viscosity (prompt D-FLOW-3).** Silica melt, Mauro et al. 2022 Table II: natural silica LA (baseline) and
  fluorine-doped F320 (the softest of their seven at 2000-2200 K, low-viscosity case), the pair the Supplement S9 estimate uses.
  Evaluated at the surface temperature T_s of the node upstream of each face, not at a film-average temperature. The film is the
  outer layer, so T_s is its hottest point and the viscosity there its lowest: flow is overestimated, not underestimated. Bound:
  the drop across the film itself, film thickness times conducted flux over k_glass, is 4-70 K for 10-50 um, 0.5-2 MW/m2 and
  k_glass = 1.4 W/m/K, which raises the LA viscosity by less than a factor of 3 at 2300 K. The whole drop to the reaction front,
  160-310 K in 2D (Section 4), would raise it by factors of 8-190 (LA) and 6-85 (F320) at 2200-2400 K; that is the outer bound of
  the error, not its size. Above 2573 K the law is extrapolated. 10 Nov 2026.
- **D-FLOW-6 Composition not advected (prompt D-FLOW-4).** The flow moves glass volume only. A node that loses glass keeps its
  B2O3 fraction, and glass arriving on a bare node takes the composition of the node it came from. Reason given by the prompt:
  B2O3 leaves the film within seconds above 1800 K (Section 4.2), so advecting it would change only the first minute. A one-off
  diagnostic run at one low-altitude boundary point quantifies the neglected effect. 10 Nov 2026.
- **D-FLOW-7 Numerics.** Operator splitting inside each oxide step: LayerOxidation.advance applies the silica source and the
  evaporation node by node (one source term), then glassflow.film_step moves the film over the same step. Flux form on the lumped
  surface weights of the conduction solver (Solver2D.w_surf), so the glass volume is conserved exactly; first-order upwind in the
  flow direction (outward: tau_w >= 0 and G >= 0 everywhere, checked); q = 0 at the stagnation point; glass leaving the last flank
  node is shed and tallied. Explicit sub-steps with CFL 0.9 on dq/d(delta); above 10^4 sub-steps in one oxide step, a linearised
  implicit upwind sweep (positive and conservative for any step). Criterion (a) is rechecked after the flow and its time
  interpolated inside the step as in advance(). The film audit (film + evaporated + shed against the cumulative source) is
  evaluated at every output time. Verified by V10 (`validation/V10_plan.md`). 10 Nov 2026.

## Thermoelastic module (version 3, Task 2)

Code: `src/uhtc/thermoelastic.py`; science runs `runs/thermoelastic/run_stress.py`. The module replaces the fully restrained index of
Section 6 with computed stress fields. It is a screen upgraded from the index, not a fracture prediction.

- **D-MECH-1 Boundary conditions.** Two bracketing cases on every run. (a) Free: u_r = 0 on the axis, every other surface
  traction-free, including the back plane, and u_z fixed at one node only to remove the rigid axial translation: the
  free-standing tip, a lower bound of the restraint. (b) Restrained back plane: u_z = 0 on the back plane of the 2D mesh (the plane
  normal to the flank at the end of the 30 mm flank, D-2D-2): the mounted tip, an upper bound. Quasi-static, linear, uncoupled;
  stress-free at 300 K, the cold start of every run. 10 Nov 2026.
- **D-MECH-2 Properties.** The two property sets of the Section 6 index, unchanged: TPSX record 28 for ZrB2-20SiC (E 530, 523,
  446 and 106 GPa at 296, 1073, 1673 and 1873 K; alpha 5.00, 5.94, 6.33 and 6.88e-6 /K at 673, 1073, 1473 and 1873 K; nu = 0.14,
  assumed in the record), and Neuman et al. 2013 for ZrB2-30SiC in air (E 510, 210 and 110 GPa at 296, 1773 and 1873 K, abstract
  values) with the TPSX alpha and nu. TPSX gives sparse points: linear interpolation in T, held constant beyond the end points
  (as `audit/tools/thermal_stress_screen.py`). Thermal strain alpha(T)(T - 300 K) with alpha the tabulated mean coefficient, as
  the prompt sets. E and alpha at the element temperature. 10 Nov 2026.
- **D-MECH-3 Weibull parameters (derived).** No open Weibull modulus for ZrB2-SiC was found in the sources read for this
  project; TPSX gives the mean and standard deviation of the 4-point bend strength at five or six temperatures (in-plane and
  through-thickness). The modulus is derived from the coefficient of variation with the exact two-parameter Weibull relation
  CoV = sqrt(Gamma(1 + 2/m) - Gamma(1 + 1/m)^2) / Gamma(1 + 1/m), and the characteristic strength as mean / Gamma(1 + 1/m). The
  derived moduli range from 5.9 to 64 (median about 25) with unknown specimen counts, so the screen uses m = 6 and m = 20 as a
  bracket, sigma_0(T) from the in-plane means (linear in T, held beyond the end points), and the effective volume of the TPSX
  bend bar (width 5.08 mm, depth 2.54 mm; outer span 38.1 mm assumed, the record gives none) by the formula of Quinn,
  V_eff = b d L_o (m L_i/L_o + 1) / (2 (m + 1)^2), with quarter-point (L_i = L_o/2) and third-point (L_i = L_o/3) loading as a
  further bracket. Failure probability of the tip P_f = 1 - exp(-sum over tensile elements of (sigma_1/sigma_0(T))^m V/V_eff).
  The survival estimate is reported as a range, never as a point value. 10 Nov 2026.

- **D-MECH-4 Snapshot times to the end of the exposure (after the independent review).** The prompt's snapshots stop at 120 s plus the
  failure time, but the five boundary cases survive the 1800 s cap and have no failure time. The independent review found the
  restrained 10 mm case still rising at 120 s (whole-body stress at the restrained back plane), so the snapshots were extended to
  200, 300, 600, 1200 and 1800 s for all five cases (`run_stress.py late`). Only the 10 mm case changes: its failure probability
  peaks at 0.109 at 200 s with the back plane restrained (0.0989 at 120 s) and stays near 0.1 to 1800 s; the free 10 mm tip peaks
  at 0.059 within 10 s. The near-tip stresses (within 2 R_n) and the other four cases peak before 120 s. A mesh check of the
  10 mm case at 3 and 200 s (`run_stress.py meshcheck10 200`) changes its stresses by less than 1 % and its failure probability by
  less than 0.2 %. 11 Nov 2026.

## Flank heating, third model (version 3, Task 3)

- **D-HEAT-7 Integral flank model.** heat_dist = "integral" (`conduction2d.integral_distribution`), a third option beside the
  baseline (Lees local similarity, D-2D-7) and the sharp-cone correlation (D-2D-12); the baseline is unchanged. Laminar integral
  method of Zoby et al. as printed by Riley 1992 (PhD thesis, NC State, NASA TM-107630, NTRS 19920016036, Eqs. 3.3.1-3.3.2; the
  original Zoby, Moss and Sutton 1981 paper was not read): q = (0.22/0.664) rho* mu* u_e r (h_aw - h_w) Pr^-0.6 /
  (int rho* mu* u_e r^2 ds)^1/2, marched along s; Eckert's intermediate enthalpy and the laminar recovery factor Pr^1/2 as given by
  Crabtree et al. 1965 (Eqs. 1.3-1.4); rho* and mu* of equilibrium air (Sutherland viscosity); cold-wall ratio q/q_s normalised by
  the stagnation limit of the same equations, used on the flank and on the bridge from 70 deg to the tangency as the baseline does.
  Edge state option (a) of the prompt (isentropic expansion from the normal-shock stagnation point to the modified-Newtonian
  pressure). Option (b), an entropy-swallowed edge state, needs a shock shape for which no source was read; the sharp-cone model,
  the fully swallowed limit, and the two unswallowed models bracket it. Compared in V12 (`validation/V12_plan.md`). 10 Nov 2026.
