# Design note: glass flow in the oxide model (version 2)

Status: design only. Nothing in `src/uhtc/` changes in version 1, and no parameter is adjusted.

## Why

V9 fails on the glass-loss clock: in the plasmatron the model keeps its glass longer than the test does (Section 3). The
model omits glass flow, and the order-of-magnitude estimate of Supplementary Section S9 shows that near the 2D boundary flow
can thin the film as fast as evaporation once the film is 11–53 µm thick at 1 mm. Flow is therefore the first candidate for
a model change. Whether it explains V9 is untested: the post hoc uncertainty analysis (`validation/VVU/`) can only say whether
the V9 error exceeds the validation uncertainty.

## Physics, from the thin-film estimate of Supplement S9

- Pressure-driven film velocity: v_p = (δ_g²/μ)(p_s/R_n). Shear-driven: v_τ = τ_w δ_g/μ, with τ_w from the laminar Reynolds
  analogy (550–4800 Pa near the boundary).
- The film leaves a region of size R_n, so the thinning rate is w_flow = δ_g (v_p + v_τ)/R_n. The estimate is good to about a
  factor of 2.
- Viscosity of silica: log₁₀ μ = log₁₀ μ∞ + E/(2.303 R T), with the fibre-draw constants of Mauro et al. (natural silica LA:
  E = 587.4 kJ/mol, log₁₀ μ∞ = −8.60 Pa s; fluorine-doped F320: 495.8 kJ/mol, −6.80), valid over 1973–2573 K. Silica is the
  right end member above about 1800 K, where the glass holds less than 0.3 mol% B₂O₃.

## Which module changes and how

- `src/uhtc/oxidation.py`, class `LayerOxidation`: the outer-glass balance in `advance` gains a sink term,
  dδ_g/dt = (supply − evaporation)/ρ_g − w_flow(δ_g, T_s, p_s, τ_w, R_n). The flow term scales as δ_g³ (pressure) and δ_g²
  (shear), so it limits a thick film and vanishes for a thin one. The time step limit on the converted depth stays, and a
  limit on the relative change of δ_g per step is added for stiffness.
- `OxParams` gains `glass_flow: bool` (off in version 1 results), `visc_case: 'LA' | 'F320'` and `flow_mult` (sensitivity,
  ×0.5–2 for the factor-of-2 accuracy of the estimate).
- The conditions objects (`coupled.Conditions` and subclasses) pass τ_w and R_n to the oxide model; the 2D model passes the
  local values of each surface node.

## New inputs

| Input | Baseline | Range | Source |
|---|---|---|---|
| Silica viscosity law | LA (natural silica) | LA, F320 | Mauro et al. (fibre draw) |
| Flow multiplier | 1 | 0.5–2 | accuracy of the thin-film estimate (S9) |
| Wall shear stress | laminar Reynolds analogy | ×0.5–2 | Crabtree et al. |

## How it would be tested

1. Unit tests: the flow term against the analytic thinning of an isothermal film of fixed supply; energy and mass bookkeeping
   unchanged with `glass_flow = False` (all version 1 results reproduce exactly).
2. V1–V8 rerun with flow on; the verdicts are reported as they fall.
3. **V9 is no longer an independent check once the model changes because of V9.** Its result with flow is reported as a
   diagnostic only.
4. A new, pre-registered check V13 (`validation/V13_plan_TEMPLATE.md`; renamed from V10 on 10 Nov 2026, because version 3 uses V10 for the verification of the flow module), filled in and committed before any run, decides
   whether the changed clock is supported.
5. The validation-uncertainty analysis of `validation/VVU/` is repeated for V13 with its sources fixed in its plan.
