# V2b plan: mesh and step convergence of the 2D model with glass flow on (version 3, Task 1)

Written and committed on 10 Nov 2026 before any V2b run. Nothing below may change after this commit; a dated `## Amendment`
section at the end records any change. V2b repeats the 2D part of V2 (`runs/envelope_2d/check_mesh.py`) with the glass-flow
module of `src/uhtc/glassflow.py` switched on.

## Why this check

The flow module moves glass between surface nodes, so the glass-loss time at the stagnation node now depends on the surface
mesh and on the step length in a way that V2 did not test. The flow maps of Task 1 may be used only if the production settings
are converged with flow on to the same tolerances as without it.

## Case

One band column: R_n = 1 mm, 30 km, fully catalytic bound, Lees heating distribution (production baseline), 5 deg sphere-cone
with a 30 mm flank, 20 vol% SiC, glass flow on with its baseline settings (LA viscosity, wall-shear multiplier 1, pressure term
on), 1800 s cap, stopped at the first failure. Mach 8.5, 9.0, 9.5 and 10.0: the grid points around the flow-off boundary of this
column (highest surviving Mach 9.25, lowest failing 9.5), extended by one grid step on the failing side because flow can only
lower the boundary.

Settings compared, as in `check_mesh.py`:
- production: mesh `common.MESH_2D` (n_arc 30, n_flank 40, n_eta 30, first cell 0.003), step cap 2 s;
- refined: n_arc 60, n_flank 80, n_eta 45, first cell 0.0015, step cap 0.5 s.

## Criterion (fixed now)

For every Mach number: the peak stagnation-surface temperature of the production run within 0.5 % of the refined run, and the
same first-failure mode. Where both fail, the failure time within 2 % of the refined run. V2b passes only if all four Mach numbers
meet every applicable part. These are the V2 tolerances for 2D.

## What a pass supports and what it does not

A pass supports using the production mesh and step cap for the flow maps of Task 1. It does not support the flow physics, the
other band columns individually or any column with a different radius, which V2b does not run; a failure of V2b means the flow maps
must be run with the refined settings or reported with the measured numerical uncertainty.

## Outputs

`validation/V2b_flow_numerics/run_v2b.py`, `v2b.jsonl` and `report.md`.
