# V11 plan: verification of the thermoelastic module (version 3, Task 2)

Written and committed on 10 Nov 2026 before the module's check script exists and before any scoring run. Nothing below may change
after this commit; a dated `## Amendment` section at the end records any change. V11 verifies the numerics of the axisymmetric
thermoelastic solver (`src/uhtc/thermoelastic.py`, quasi-static, linear, uncoupled, linear triangles on the body-fitted meshes of the
2D model) against exact solutions. It says nothing about strength, fracture or the properties of ZrB2-SiC.

## Exact solutions (constant E, nu, alpha; stress-free at the reference temperature)

`validation/V11_thermoelastic/exact.py` evaluates the general solutions of the displacement equations with the constants fixed by the
boundary conditions (Timoshenko and Goodier, Theory of Elasticity, 3rd ed., Sections 150 and 152), and checks them before use:
equilibrium and boundary conditions hold to 3e-7 relative (finite-difference check), and the sphere's hoop stress equals the closed
form for steady conduction to machine precision.

Common properties: E = 400 GPa, nu = 0.14, alpha = 6e-6 /K.

**Check 0 (patch test).** A small irregular axisymmetric mesh, uniform temperature rise of 1000 K, only the rigid-body modes fixed
(u_r = 0 on the axis, u_z = 0 at one node). Exact solution: uniform strain alpha dT in every direction, zero stress.
Criterion: every element stress below 1e-9 E alpha dT, nodal displacements within 1e-9 relative of alpha dT x position.

**Check 1 (thick-walled sphere).** Hollow sphere a = 0.5 mm, b = 1 mm (full meridian, axisymmetric mesh between the two spheres),
steady radial temperature T(r) = 1000 K (1/r - 1/b)/(1/a - 1/b), traction-free surfaces. Hoop stress at the inner and outer surfaces
(-1993 and +797 MPa exactly), taken from the element nearest each surface on the equator and extrapolated linearly to the surface
from the two nearest element centroids along the radius. Criterion: within 2 % at both surfaces on the finest of the meshes below.

**Check 2 (free thin disc).** Disc of radius b = 1 mm and thickness 0.02 mm (h/b = 0.02), T(r) = 1000 K (r/b)^2, traction-free,
plane-stress solution: hoop stress +600 MPa at the centre and -1200 MPa at the rim. Compared at r = 0.25 b, 0.5 b and 0.9 b on the
mid-plane. Criterion: within 2 % at all three radii.

**Convergence.** Check 1 on four meshes, each doubling the cells along the meridian and through the wall. Observed order of the
hoop-stress error at the inner surface from the last three meshes. Criterion: order >= 0.8 (the V1 style).

V11 passes only if checks 0, 1 and 2 and the convergence criterion are all met.

## What a pass supports and what it does not

A pass supports the element formulation, the axisymmetric hoop term, the thermal load, the boundary conditions and the stress
recovery. It does not support the temperature-dependent properties of ZrB2-SiC (D-MECH-2), the restraint of a real tip (D-MECH-1),
the Weibull parameters (D-MECH-3) or any statement about fracture. A failed first run is recorded as it falls; a redesign is logged.

## Mesh check of the science runs (rule of the prompt, fixed now)

For the science runs the stress at the probe points (the largest first principal stress within 2 R_n of the stagnation point, and
its value at the stagnation point) is computed on the production mesh and on a mesh with twice the cells along the surface and
through the thickness, for one case. If they differ by more than 5 %, the elasticity solve moves to quadratic triangles on the same
geometry.

## Outputs

`validation/V11_thermoelastic/run_v11.py`, `v11_results.json`, `report.md`.

## Amendment 1 (10 Nov 2026, after the first science runs)

The mesh check of the science runs (1 mm, 20 km, Mach 8.5, 10 s) missed the 5 % rule with the linear triangles: the largest first
principal stress within 2 R_n was 16.4 MPa on the production mesh and 6.0 MPa on the doubled mesh, the peak lying at the centre of the
nose arc, where the arc stations of the mesh meet in a fan of thin triangles. As the plan requires, the elasticity solve moved to
6-node (quadratic) straight-sided triangles on the same geometry (`thermoelastic.solve_quadratic`). The checks of this plan were
repeated with them and pass (`run_v11_quadratic.py`, `report_quadratic.md`: patch test 2e-13, sphere 0.85 % inner and 0.19 % outer,
order 0.96, disc exact). With quadratic triangles the mesh check gives 4.9 against 4.6 MPa (TPSX properties, 6.5 %) and 4.91 against
4.70 MPa (Neuman, 4.5 %). The plan gives no further step, so the stresses of the science runs are reported with a mesh uncertainty
of +-7 %, and a second mesh check at 3 s, near the time of the peak index, is added and reported. The linear results are kept in
`runs/thermoelastic/stress_linear_superseded.jsonl` and not used.
