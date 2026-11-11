# Glass flow away from the stagnation point: order-of-magnitude check (task 2.5)

Not added to the model. Script: `audit/tools/glass_flow_estimate.py`; numbers: `audit/glass_flow_table.csv`.

## Result

**Flow can rival or exceed evaporation near the boundary, so this is flagged as a major limitation for the author.**

At the fully catalytic 2D boundary (20, 30 and 40 km; the last surviving and the first failing point of each column,
surface 2060–2420 K), the model's evaporation thins the glass at 0.004–0.16 µm/s. Flow thins it at the same rate once the
film has reached the thickness below:

| R_n | Film thickness at which flow equals evaporation, natural silica (LA) | Same, fluorine-doped silica (F320) | Model's film at the surviving point, end of run |
|---|---|---|---|
| 1 mm | 14–53 µm | 11–39 µm | 15–38 µm |
| 3 mm | 23–92 µm | 19–68 µm | 11–34 µm |
| 10 mm | 40–167 µm | 34–136 µm | 10–78 µm |

The ratio of flow to evaporation (LA silica) is 0.02–0.45 at 10 µm and 0.8–35 at 50 µm for 1 mm, 0.2–8 at 50 µm for 3 mm
and 0.05–1.8 at 50 µm for 10 mm. So:

- At 10 µm evaporation dominates everywhere: flow is at most 0.45 of it (LA), 0.83 with the softer F320 curve.
- At 50 µm flow is comparable to or faster than evaporation at 1 mm, and at the low altitudes for 3 and 10 mm.
- The model's own films reach this size. Near the 1 mm boundary they are at or above the matching thickness at 20 km
  (38 against 14 µm), about equal at 30 km (26 against 28 µm) and below it at 40 km (15 against 50 µm). 1D histories
  run thicker, at 40–140 µm.

## What this means for the paper

- **Direction: the omission is not conservative.** Flow takes glass away from the stagnation point and nothing replaces
  it from upstream, so the modelled glass-loss times are too long if flow matters. This is the opposite sign to most of
  the model's other simplifications.
- **Flow is self-limiting.** Its thinning rate scales as the square to the cube of the thickness, so it trims a thick film
  quickly and a thin one hardly at all. It cannot remove the film alone; it caps the film at roughly the thickness above,
  and evaporation then has less glass to consume.
- **Rough size of the effect on the boundary.** Evaporation typically grows by a factor of 1.5–2.5 over one 0.25 Mach step
  at these points (15 at the 1 mm, 40 km column, where the surface jumps 160 K across the step). An extra removal of the
  same order as evaporation is therefore worth about one search step. This is an inference from the table, not a run.
- Task 2.5 forbids adding flow to the model, so no boundary has been recomputed with it.

## Method

Velocities are the forms set in the brief, with the film leaving a region of size R_n:

- pressure driven: v = (δ²/μ)(p_s/R_n)
- shear driven: v = τ_w δ/μ
- thinning rate by continuity: w = δ v / R_n, summed over the two

The exact lubrication means are v/3 and v/2, and axisymmetric divergence adds a factor of 2, so w is good to about a
factor of 2 either way.

**Viscosity.** Measured fibre-draw constants of Mauro et al., J. Am. Ceram. Soc. 105 (2022) 61–66, Table II, valid over
1973–2573 K, so nothing is extrapolated: log₁₀ μ[Pa s] = log₁₀ μ∞ + E/(2.303 R T).

| Sample | E (kJ/mol) | log₁₀ μ∞ | μ at 2000 K | 2200 K | 2300 K | 2450 K |
|---|---|---|---|---|---|---|
| LA, natural silica | 587.4 | −8.60 | 5.5×10⁶ | 2.2×10⁵ | 5.5×10⁴ | 8.4×10³ Pa s |
| F320, fluorine doped | 495.8 | −6.80 | 1.4×10⁶ | 9.4×10⁴ | 2.9×10⁴ | 5.9×10³ Pa s |

Silica is the right end member here because the model's glass holds less than 0.3 mol% B₂O₃ at these temperatures
(`runs/histories/histories.csv`). The Urbain et al. 1982 data could not be read (paywalled); Shao et al. (arXiv:1305.2236)
quote their activation energy as about 551 kJ/mol over 1200–2400 °C, between the two samples used.

**Wall shear.** Reynolds analogy for a laminar boundary layer, c_f/2 = St Pr^(2/3) (Crabtree et al. 1965, Section 1.3),
with g = cos φ, p/p_s = cos²φ and isentropic edge velocity, Pr = 0.71, taking the largest value over φ ≤ 70°. This gives
τ_w = 550–4800 Pa.

## Caveats

- Real glass is not pure silica. Dissolved ZrO₂ would stiffen it and residual B₂O₃ would soften it; no data for such a melt
  at these temperatures were found. Across the seven silicas of Mauro et al., with impurities at the ppm to 0.15 % level, viscosity
  spreads by a factor of 4.7 at 2000 K and 2.4 at 2450 K.
- Part of the glass sits in the pores of the ZrO₂ skeleton and is held there; the estimate treats the outer film as free.
- Surface temperature falls away from the stagnation point, so viscosity rises downstream and the film would thicken
  there. The estimate uses the stagnation temperature throughout, which overstates flow.
- Pressure gradient p_s/R_n and shear are order-of-magnitude values, not local solutions.
