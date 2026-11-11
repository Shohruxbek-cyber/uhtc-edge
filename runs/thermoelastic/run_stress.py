"""Version 3, Task 2: computed thermal-stress fields of the 2D nose tip (replaces the fully restrained index of Section 6).

Cases: the five fully catalytic boundary conditions of the Section 6 index (audit/tools/thermal_stress_screen.py), cold start
from 300 K at constant flight conditions, 2D model as the index ran it (Lees distribution, production mesh, step cap 0.5 s);
and the X-43A Mach 10 path for the 1 mm tip (runs/trajectory, 2D). Temperature fields at snapshot times (each from its own cold
start: run_2d to t), stress by uhtc.thermoelastic for the free and the restrained-back-plane tip (D-MECH-1) and the TPSX and
Neuman property sets (D-MECH-2), Weibull failure probability bracket (D-MECH-3).

Usage: python runs/thermoelastic/run_stress.py [cases|late|x43a|meshcheck|all]   -> runs/thermoelastic/stress.jsonl
       python runs/thermoelastic/run_stress.py meshcheck10 <t_s>   (mesh check of the 10 mm case at 3 s and at t_s)
"""
import json
import math
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "runs"))
sys.path.insert(0, str(ROOT / "src"))
import common  # noqa: E402
from uhtc import coupled as cp, materials as mt  # noqa: E402
from uhtc.conduction2d import NoseGeometry  # noqa: E402
from uhtc.oxidation import LayerOxidation, OxParams  # noqa: E402
from uhtc.thermoelastic import Props, back_plane_nodes, principal_max, solve_quadratic  # noqa: E402

# D-MECH-2: property sets of the Section 6 index (audit/tools/thermal_stress_screen.py)
E_T, E_V = [296.15, 1073.15, 1673.15, 1873.15], [530e9, 523e9, 446e9, 106e9]
A_T, A_V = [673.15, 1073.15, 1473.15, 1873.15], [5.00e-6, 5.94e-6, 6.33e-6, 6.88e-6]
EN_T, EN_V = [296.15, 1773.15, 1873.15], [510e9, 210e9, 110e9]
S_T, S_V = [296.15, 673.15, 1073.15, 1473.15, 1873.15], [430e6, 548e6, 463e6, 342e6, 239e6]  # TPSX in-plane bend means
SN_T, SN_V = [296.15, 1073.15, 1773.15, 1873.15], [680e6, 750e6, 360e6, 360e6]
NU = 0.14
PROPS = {"TPSX": Props(lambda T: np.interp(T, E_T, E_V), lambda T: np.interp(T, A_T, A_V), NU),
         "Neuman": Props(lambda T: np.interp(T, EN_T, EN_V), lambda T: np.interp(T, A_T, A_V), NU)}
STRENGTH = {"TPSX": (S_T, S_V), "Neuman": (SN_T, SN_V)}
CASES = [(8.5, 20e3, 1e-3), (9.25, 30e3, 1e-3), (10.25, 40e3, 1e-3), (9.75, 30e3, 3e-3), (10.25, 30e3, 10e-3)]
TIMES = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 60.0, 120.0]
# D-MECH-4 (11 Nov 2026, after the independent review): the five cases survive the 1800 s cap, so they have no failure time; the
# snapshots are extended to the end of the exposure, because the restrained 10 mm case was still rising at 120 s
TIMES_LATE = [200.0, 300.0, 600.0, 1200.0, 1800.0]
X43A_TIMES = [10.0, 30.0, 60.0, 100.0, 150.0, 200.0, 300.0, 400.0, 500.0, 600.0]
# D-MECH-3: Weibull bracket; TPSX bend bar 5.08 x 2.54 mm, outer span 38.1 mm (assumed), quarter- and third-point loading
BAR_B, BAR_D, BAR_LO = 5.08e-3, 2.54e-3, 38.1e-3
WEIBULL_M = (6.0, 20.0)


def v_eff(m, Li_over_Lo):
    return BAR_B * BAR_D * BAR_LO * (m * Li_over_Lo + 1.0) / (2.0 * (m + 1.0) ** 2)


def field(M, z, R, t_end, mesh=None, traj=False):
    p = OxParams(f_SiC=0.2)
    ox = LayerOxidation(p)
    mat = mt.composite(0.2)
    if traj:
        sys.path.insert(0, str(ROOT / "runs" / "trajectory"))
        import run_x43a
        t, zz, MM = run_x43a.load_traj()
        cond = cp.TrajectoryConditions(t - t[0], zz, MM, R, dt_eval=2.0, state_fn=common.stag_state, catalysis="fc")
    else:
        cond = cp.ConstantConditions(M, z, R, state=common.stag_state(M, z), catalysis="fc")
    g = NoseGeometry(R_n=R, theta_c=math.radians(5.0), L_flank=0.03, axisymmetric=True)
    r = cp.run_2d(cond, g, mat, ox, t_end, dt_max=0.5, mesh_kw=mesh or common.MESH_2D, heat_dist="lees")
    return r.extra["solver"], r


def analyse(sol, T, R):
    """Quadratic (6-node) triangles: the linear solve missed the 5 % mesh check (V11 plan, amendment 1). Stresses at the Gauss
    points; a point counts as 'at the surface' if its element has a surface node."""
    nodes, tri = sol.nodes, sol.tri
    back = back_plane_nodes(sol)
    tip = nodes[sol.surf[0]]
    surf_set = set(int(i) for i in sol.surf)
    n_el = len(tri)
    touches = np.array([any(int(v) in surf_set for v in t) for t in tri])
    out = {}
    for bc in ("free", "restrained"):
        for pname, props in PROPS.items():
            _, sig, cent, info = solve_quadratic(nodes, tri, T, props, bc=bc, back_nodes=back)
            s1 = principal_max(sig)
            d = np.hypot(*(cent - tip).T)
            Te = info["T"]
            near = d <= 2 * R
            k = int(np.argmax(np.where(near, s1, -np.inf)))
            kall = int(np.argmax(s1))
            on_surface = bool(touches[k % n_el])
            St, Sv = STRENGTH[pname]
            ratio = s1 / np.interp(Te, St, Sv)
            # Weibull (D-MECH-3): characteristic strength from the TPSX in-plane means for both sets (TPSX is the only scatter source)
            vol = 2 * np.pi * info["vol"]
            pf = {}
            for m in WEIBULL_M:
                s0 = np.interp(Te, S_T, S_V) / math.gamma(1.0 + 1.0 / m)
                for tag, li in (("quarter", 0.5), ("third", 1.0 / 3.0)):
                    risk = float(np.sum(np.where(s1 > 0, (np.maximum(s1, 0) / s0) ** m, 0.0) * vol) / v_eff(m, li))
                    pf[f"m{m:g}_{tag}"] = 1.0 - math.exp(-risk)
            out[f"{bc}_{pname}"] = dict(s1_max_2R_MPa=float(s1[k]) / 1e6, at_surface=bool(on_surface), dist_mm=float(d[k] * 1e3),
                                        T_at_K=float(Te[k]), ratio_max_2R=float(np.max(np.where(near, ratio, -np.inf))),
                                        s1_max_body_MPa=float(s1[kall]) / 1e6, dist_body_mm=float(d[kall] * 1e3),
                                        ratio_max_body=float(np.max(ratio)), pf_min=min(pf.values()), pf_max=max(pf.values()), pf=pf)
    return out


def one(args):
    kind, M, z, R, t = args
    if kind == "x43a":
        sol, r = field(None, None, R, t, traj=True)
    else:
        sol, r = field(M, z, R, t, mesh=None if kind != "fine" else dict(n_arc=60, n_flank=80, n_eta=60, eta1=0.0015))
    T = np.asarray(sol.T)
    return dict(kind=kind, mach=M, alt_km=None if z is None else z / 1e3, R_mm=R * 1e3, t_s=t, T_stag=float(r.T_s[-1]),
                T_max=float(T.max()), T_min_2R=float(T[np.hypot(*(sol.nodes - sol.nodes[sol.surf[0]]).T) <= 2 * R].min()),
                n_nodes=len(T), **analyse(sol, T, R))


def failure_time(M, z, R):
    rows = {(r["R_n"], r["alt"], r["mach"]): r for f in ("grid", "refine") for r in common.load(ROOT / f"runs/envelope_2d/{f}.jsonl")
            if r["catalysis"] == "fc"}
    r = rows.get((R, z, M))
    return None if r is None or r["t_fail"] != r["t_fail"] else float(r["t_fail"])


def main():
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    jobs = []
    if what in ("cases", "all"):
        for M, z, R in CASES:
            ts = list(TIMES)
            tf = failure_time(M, z, R)
            if tf is not None and tf not in ts:
                ts.append(math.floor(tf * 10) / 10)
            jobs += [("case", M, z, R, t) for t in ts]
    if what in ("late", "all"):
        jobs += [("case", M, z, R, t) for M, z, R in CASES for t in TIMES_LATE]
    if what == "meshcheck10":
        M, z, R = CASES[4]
        for t in (3.0, float(sys.argv[2])):
            jobs += [("case", M, z, R, t), ("fine", M, z, R, t)]
    if what in ("x43a", "all"):
        jobs += [("x43a", None, None, 1e-3, t) for t in X43A_TIMES]
    if what in ("meshcheck", "all"):
        M, z, R = CASES[0]
        jobs += [("case", M, z, R, 10.0), ("fine", M, z, R, 10.0), ("case", M, z, R, 3.0), ("fine", M, z, R, 3.0)]
    out = HERE / "stress.jsonl"
    done = set()
    if out.exists():
        done = {(r["kind"], r["mach"], r["alt_km"], r["R_mm"], r["t_s"]) for r in map(json.loads, out.read_text().splitlines())}
    todo = [j for j in jobs if (j[0], j[1], None if j[2] is None else j[2] / 1e3, j[3] * 1e3, j[4]) not in done]
    print(f"{len(jobs)} jobs, {len(todo)} to run", flush=True)
    with Pool(4) as pool, open(out, "a", encoding="utf-8") as f:
        for row in pool.imap_unordered(one, todo):
            f.write(json.dumps(row, default=float) + "\n")
            f.flush()
            print(row["kind"], row["mach"], row["alt_km"], row["R_mm"], row["t_s"], round(row["free_TPSX"]["s1_max_2R_MPa"]),
                  round(row["restrained_TPSX"]["s1_max_2R_MPa"]), flush=True)


if __name__ == "__main__":
    main()
