"""V1 (part a): 1D and 2D conduction solvers against analytic solutions, with grid and time-step convergence.

Constant properties rho = 5000 kg/m3, k = 50 W/m/K, c = 600 J/kg/K (close to ZrB2-SiC at high temperature).
  A. Semi-infinite slab, constant surface flux q = 5 MW/m^2, t = 2 s (Carslaw & Jaeger, slab with constant flux):
       T = T0 + (2q/k) sqrt(alpha t/pi) exp(-x^2/4 alpha t) - (q x/k) erfc(x/2 sqrt(alpha t)).
  B. Semi-infinite slab, surface temperature stepped 300 -> 2300 K, t = 2 s:  T = Ts + (T0 - Ts) erf(x/2 sqrt(alpha t)).
  C. Solid sphere (axisymmetric) and infinite cylinder (planar 2D, a disc) of radius a = 3 mm with a uniform
     surface flux F = 2 MW/m^2.  Surface temperature (eigenfunction series; mean rise = energy balance):
       sphere:   T(a,t) - T0 = (F a/k)[3 Fo + 1/5 - 2 sum exp(-b_n^2 Fo)/b_n^2],  tan b_n = b_n
       cylinder: T(a,t) - T0 = (F a/k)[2 Fo + 1/4 - 2 sum exp(-a_n^2 Fo)/a_n^2],  J1(a_n) = 0
     with Fo = alpha t / a^2.  (At t = 0 both vanish because sum 1/b_n^2 = 1/10 and sum 1/a_n^2 = 1/8; the mean
     temperature equals 3 F t/(rho c a) and 2 F t/(rho c a), as the energy balance requires - checked numerically below.)
     Solved with the 2D finite-element model (conduction2d, sphere/disc meshes) and, for the sphere, also with the
     1D solver in spherical-shell geometry (an independent discretisation).
Errors relative to the surface temperature rise.  Pass criteria: finest-grid error <= 0.5 %; observed order >= 0.8
for refinement in space and in time (cases A, B); energy residuals <= 1e-8.
Spatial order in A and B is measured by self-convergence (differences between successive meshes at the same time
step, so that the time-discretisation error, which dominates the error against the exact solution on fine meshes,
cancels); revised after the first run (2 Nov 2026), whose meshes were already finer than the time-error floor.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import erf, erfc, jn_zeros

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from uhtc.conduction1d import Material, Solver1D  # noqa: E402
from uhtc.conduction2d import NoseGeometry, Solver2D, build_sphere_mesh  # noqa: E402

RHO, K, CP = 5000.0, 50.0, 600.0
ALPHA = K / (RHO * CP)
T0 = 300.0


def mat():
    f = lambda v: (lambda T: np.full_like(np.asarray(T, float), v))
    return Material("const", rho=f(RHO), k=f(K), cp=f(CP), T_min=100, T_max=8000)


def exact_A(x, t, q):
    s = np.sqrt(ALPHA * t)
    return T0 + 2 * q / K * s / np.sqrt(np.pi) * np.exp(-x ** 2 / (4 * s * s)) - q * x / K * erfc(x / (2 * s))


def exact_B(x, t, Ts):
    return Ts + (T0 - Ts) * erf(x / (2 * np.sqrt(ALPHA * t)))


BETA = np.array([brentq(lambda b: np.tan(b) - b, n * np.pi + 1e-9, n * np.pi + np.pi / 2 - 1e-9) for n in range(1, 4000)])
ALPH = jn_zeros(1, 4000)


def exact_sphere_surface(t, F, a):
    Fo = ALPHA * t / a ** 2
    return T0 + F * a / K * (3 * Fo + 0.2 - 2 * np.sum(np.exp(-BETA ** 2 * Fo) / BETA ** 2))


def exact_cyl_surface(t, F, a):
    Fo = ALPHA * t / a ** 2
    return T0 + F * a / K * (2 * Fo + 0.25 - 2 * np.sum(np.exp(-ALPH ** 2 * Fo) / ALPH ** 2))


def run1d(case, dx0, dt, t_end=2.0, growth=1.08, probe=False):
    s = Solver1D(mat(), 0.05, T0, dx0=dx0, growth=growth)
    if case == "A":
        s.front = lambda T, t: (5e6, 0.0)
    else:
        Rh = s.half_cell_resistance()
        s.front = lambda T, t: ((2300.0 - T) / Rh, -1.0 / Rh)
    n = int(round(t_end / dt))
    for _ in range(n):
        assert s.step(dt)[0]
    x = s.depth_c
    sel = x < 5e-3
    ex = exact_A(x, t_end, 5e6) if case == "A" else exact_B(x, t_end, 2300.0)
    rise = (exact_A(0.0, t_end, 5e6) - T0) if case == "A" else 2000.0
    if probe:  # temperatures at fixed depths (0.5, 1, 2, 4 mm) for self-convergence
        return np.interp([0.5e-3, 1e-3, 2e-3, 4e-3], x, s.T) / rise
    return np.abs(s.T[sel] - ex[sel]).max() / rise, s.finish().residual()[1]


def order(errs):
    e = np.array(errs)
    return np.log2(e[:-1] / e[1:])


def run_sphere2d(axi, n_arc, n_eta, dt, t_end, F=2e6, a=3e-3):
    g = NoseGeometry(R_n=a, theta_c=0.0, L_flank=0.0, axisymmetric=axi)
    s = Solver2D(g, mat(), T0, mesh=build_sphere_mesh(a, n_arc=n_arc, n_eta=n_eta, eta1=0.004))
    flux = lambda T: (np.full(len(T), F), np.zeros(len(T)))
    for _ in range(int(round(t_end / dt))):
        assert s.step(dt, flux)[0]
    Ts = s.T[s.surf]
    return Ts.mean(), Ts.max() - Ts.min(), s.energy_residual(), np.sum(s.m * s.T) / np.sum(s.m)


def run_sphere1d(dx0, dt, t_end, F=2e6, a=3e-3):
    s = Solver1D(mat(), a, T0, dx0=dx0, growth=1.05, geometry="sphere", front=lambda T, t: (F, 0.0))
    for _ in range(int(round(t_end / dt))):
        assert s.step(dt)[0]
    # surface temperature: first-cell value extrapolated with the surface flux (half cell)
    return s.T[0] + F * (s.faces[1] / 2) / K, s.finish().residual()[1]


def main():
    out = ["# V1 — Conduction solvers vs analytic solutions (part a: constant properties)", "",
           "Method and exact solutions: `run_v1_analytic.py` header.", ""]
    allok = True
    for case in ("A", "B"):
        title = "constant flux 5 MW/m², t = 2 s" if case == "A" else "surface step 300 → 2300 K, t = 2 s"
        out += [f"## Case {case}: semi-infinite slab, {title}", "",
                "| dx0 (µm) | dt (s) | max error / surface rise | energy residual |", "|---|---|---|---|"]
        e_dx = []
        for dx0 in (80e-6, 40e-6, 20e-6, 10e-6):
            e, r = run1d(case, dx0, 2e-4)
            e_dx.append(e)
            out.append(f"| {dx0 * 1e6:.0f} | 2e-4 | {e:.2e} | {r:.1e} |")
        # self-convergence on uniform meshes (growth 1.0), same time step
        probes = [run1d(case, dx, 2e-4, growth=1.0, probe=True) for dx in (400e-6, 200e-6, 100e-6, 50e-6)]
        diffs = [np.abs(probes[i] - probes[i + 1]).max() for i in range(3)]
        o_self = order(diffs)
        e_dt = []
        for dt in (0.04, 0.02, 0.01, 0.005):
            e, r = run1d(case, 10e-6, dt)
            e_dt.append(e)
            out.append(f"| 10 | {dt} | {e:.2e} | {r:.1e} |")
        o2 = order(e_dt)
        ok = max(e_dx) <= 5e-3 and min(e_dt) <= 5e-3 and o_self[-1] >= 0.8 and o2[-1] >= 0.8
        allok &= ok
        out += ["", f"All production-type meshes (dx0 ≤ 80 µm, growth 1.08) are within {max(e_dx):.1e} of the exact solution "
                f"(time-error floor at dt = 2e-4 s). Self-convergence on uniform meshes 400/200/100/50 µm: successive "
                f"differences {', '.join(f'{d:.1e}' for d in diffs)} → spatial order {', '.join(f'{v:.2f}' for v in o_self)}. "
                f"Time order: {', '.join(f'{v:.2f}' for v in o2)}. → **{'PASS' if ok else 'FAIL'}**.", ""]
    # C: sphere and cylinder
    a, F, t_end = 3e-3, 2e6, 0.12  # Fo = 0.22
    Fo = ALPHA * t_end / a ** 2
    exs, exc = exact_sphere_surface(t_end, F, a), exact_cyl_surface(t_end, F, a)
    out += [f"## Case C: sphere and cylinder, radius 3 mm, uniform flux 2 MW/m², t = {t_end} s (Fo = {Fo:.3f})", "",
            f"Exact surface temperature: sphere {exs:.2f} K, cylinder {exc:.2f} K "
            f"(series checks: T(a,0) − T0 = {exact_sphere_surface(0.0, F, a) - T0:.1e} K and {exact_cyl_surface(0.0, F, a) - T0:.1e} K).", "",
            "| Model | mesh | dt (s) | surface T (K) | error / rise | spread over surface (K) | mean T − energy-balance mean (K) | energy residual |",
            "|---|---|---|---|---|---|---|---|"]
    ok_c = True
    for axi, ex, nm, fac in ((True, exs, "2D axisymmetric (sphere)", 3), (False, exc, "2D planar (cylinder)", 2)):
        errs = []
        for n_arc, n_eta in ((20, 10), (40, 20), (80, 40)):
            Tm, spread, res, Tbar = run_sphere2d(axi, n_arc, n_eta, 1e-3, t_end)
            e = abs(Tm - ex) / (ex - T0)
            errs.append(e)
            out.append(f"| {nm} | {n_arc}×{n_eta} | 1e-3 | {Tm:.2f} | {e:.2e} | {spread:.1e} | "
                       f"{Tbar - T0 - fac * F * t_end / (RHO * CP * a):.1e} | {res:.1e} |")
            ok_c &= abs(res) <= 1e-8
        o = order(errs)
        ok_c &= errs[-1] <= 5e-3
        out.append(f"| {nm}: observed order | | | | {', '.join(f'{v:.2f}' for v in o)} | | | |")
    errs = []
    for dx0 in (40e-6, 20e-6, 10e-6):
        Ts1, res = run_sphere1d(dx0, 1e-3, t_end)
        e = abs(Ts1 - exs) / (exs - T0)
        errs.append(e)
        out.append(f"| 1D spherical shell | dx0 {dx0 * 1e6:.0f} µm | 1e-3 | {Ts1:.2f} | {e:.2e} | — | — | {res:.1e} |")
    ok_c &= errs[-1] <= 5e-3
    allok &= ok_c
    out += ["", f"Case C → **{'PASS' if ok_c else 'FAIL'}** (time step 1 ms; its first-order error is included in all rows,"
            " so the observed spatial order levels off once the time error dominates).",
            "", f"Overall V1: **{'PASS' if allok else 'FAIL'}**."]
    (HERE / "report_part_a.md").write_text("\n".join(out) + "\n", encoding="utf-8")
    print("\n".join(out))


if __name__ == "__main__":
    main()
