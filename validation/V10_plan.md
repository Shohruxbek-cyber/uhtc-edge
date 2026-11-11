# V10 plan: verification of the glass-flow module (version 3, Task 1)

Written and committed on 10 Nov 2026 before the check script exists and before any scoring run. Nothing below may change after
this commit; a dated `## Amendment` section at the end records any change, and the text above it stays as committed. V10 is a
verification of the numerics of `src/uhtc/glassflow.py` against exact solutions of the film equation. It checks no physics, and it
cannot support any statement about the glass-loss times of the maps.

## Why this check

The flow module (`uhtc.glassflow.film_step`) integrates

    d(delta)/dt + (1/r) d(r q)/ds = 0,   q = tau delta^2 / (2 mu) + G delta^3 / (3 mu),   G = -dp/ds,

on the surface nodes, first-order upwind, explicit under a CFL limit or with a linearised implicit sweep. Before the module is
used for any map, its discretisation must reproduce exact solutions of the same equation. The geometry is a cylinder of constant
radius R (flow along the axis), on which the axisymmetric operator reduces to the planar one, so the exact solutions below hold.

## Test cases and exact solutions

Common setup: cylinder radius R = 1 mm, length L = 10 mm, N = 2000 uniform cells (node spacing ds = L/N, lumped weights
w_i = R ds and R ds/2 at the two ends, face radius R), constant viscosity mu = 1e4 Pa s, CFL number 0.9.

**Check 1 (shear only).** tau = 1000 Pa, G = 0, no source, no evaporation; initial film delta(s, 0) = delta0 = 20 um on all nodes;
no glass enters at s = 0. The characteristic speed is dq/d(delta) = tau delta / mu, so the film thins by an expansion fan from s = 0:

    delta(s, t) = mu s / (tau t)          for s <= tau delta0 t / mu,
    delta(s, t) = delta0                  beyond the head of the fan.

Compared at s* = L/4 = 2.5 mm, inside the fan, at t = 2, 4 and 8 times t_head = mu s* / (tau delta0) (the time at which the head of
the fan passes s*): delta(s*, t) = delta0 t_head / t.

**Check 1b (pressure gradient only).** tau = 0, G = 1e6 Pa/m, otherwise as check 1. Characteristic speed G delta^2 / mu, fan

    delta(s, t) = sqrt(mu s / (G t))      for s <= G delta0^2 t / mu,

compared at s* = L/4 at t = 2, 4 and 8 times t_head = mu s* / (G delta0^2): delta(s*, t) = delta0 sqrt(t_head / t).

**Check 2 (steady state, source balancing outflow).** tau = 1000 Pa, G = 0, uniform source S = 1e-8 m/s (volume per area per time),
no evaporation, starting from a bare surface. The source is applied by operator splitting as in `run_2d` (the film gains S dt, then
moves). At steady state q(s) = S s, so

    delta(s) = sqrt(2 mu S s / tau),

compared at s = L/4, L/2 and L after 20 transit times of the film over L. Run twice: with the explicit sub-steps and with the
implicit sweep forced (max_explicit = 0, dt = 10 s).

## Criterion (fixed now)

Each compared value within 1 % of the exact solution: three times in check 1, three times in check 1b, three positions in check 2
for each of the two integrators. V10 passes only if all 15 values meet the criterion. The mass audit of every run (film + shed
against the initial film and the cumulative source) must close to 1e-10 relative; this is a unit test, not part of the criterion.

## What a pass supports and what it does not

A pass supports the discretisation of the film equation, the upwind direction, the source splitting and both integrators. It does
not support the wall shear, the pressure gradient or the viscosity used in the maps (D-FLOW-3 to D-FLOW-5), the neglect of the
B2O3 advection (D-FLOW-6), the lubrication approximation itself, or any glass-loss time. A failed first run is recorded as it falls;
a redesign is logged the way the V1 order-test redesign was.

## Outputs

`validation/V10_glassflow/run_v10.py`, `v10_results.json` and `report.md` with every compared value, its exact value and the
relative error.

## Amendment 1 (10 Nov 2026, after the first run)

The criterion paragraph says that V10 passes only if "all 15 values" meet the criterion. That count is an arithmetic slip: the
values the plan lists are three times in check 1, three in check 1b and three positions in check 2 for each of two integrators,
12 in all. The criterion, every listed value within 1 %, is unchanged, and the first run is reported as it fell (12 of 12).
