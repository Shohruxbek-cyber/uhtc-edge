"""V11 checks repeated with the 6-node (quadratic) triangles (V11 plan, amendment 1): the mesh check of the science runs moved the
elasticity solve to quadratic triangles, so the same exact solutions and criteria verify that solver.

Usage: python validation/V11_thermoelastic/run_v11_quadratic.py   -> v11_quadratic_results.json, report_quadratic.md
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_v11 as v  # noqa: E402
from uhtc.thermoelastic import solve_quadratic  # noqa: E402


def check0():
    nodes, tri = v.grid(4, 3, 1e-3, 0.7e-3, jitter=0.3)
    dT = 1000.0
    k = int(np.argmin(nodes[:, 0] + nodes[:, 1]))
    U, sig, _, info = solve_quadratic(nodes, tri, np.full(len(nodes), 300.0 + dT), v.P, bc="free", z_fix_node=k)
    n6 = info["nodes6"]
    s_err = float(np.abs(sig).max() / (v.E * v.AL * dT))
    ue = np.c_[v.AL * dT * (n6[:, 0] - nodes[k, 0]), v.AL * dT * n6[:, 1]]
    u_err = float(np.abs(U - ue).max() / (v.AL * dT * 1e-3))
    return dict(stress_rel=s_err, disp_rel=u_err, ok=s_err < 1e-9 and u_err < 1e-9)


def check1():
    a, b = 0.5e-3, 1e-3
    T = lambda r: 300.0 + 1000.0 * (1 / r - 1 / b) / (1 / a - 1 / b)  # noqa: E731
    ex = v.exact.sphere(T, a, b, v.E, v.NU, v.AL, T_ref=300.0)
    ea, eb = ex(a)[1], ex(b)[1]
    rows = []
    for n_arc, n_eta in ((40, 8), (80, 16), (160, 32), (320, 64)):
        nodes, tri = v.sphere_mesh(n_arc, n_eta, a, b)
        # quadratic mid-side nodes of a straight-sided mesh do not lie on the circles: use the exact T at every node
        _, sig, pos, _ = solve_quadratic(nodes, tri, T(np.hypot(nodes[:, 0], nodes[:, 1])), v.P, bc="free",
                                         z_fix_node=int(np.argmax(nodes[:, 0])))
        sa, sb = v.surface_hoop(pos, sig, a, b, "a"), v.surface_hoop(pos, sig, a, b, "b")
        rows.append(dict(mesh=f"{n_arc}x{n_eta}", h=b / n_eta, hoop_a=sa, hoop_b=sb, err_a=abs(sa / ea - 1), err_b=abs(sb / eb - 1)))
    h = np.array([r["h"] for r in rows[-3:]])
    e = np.array([r["err_a"] for r in rows[-3:]])
    order = float(np.polyfit(np.log(h), np.log(e), 1)[0])
    fin = rows[-1]
    return dict(exact_a=ea, exact_b=eb, meshes=rows, order=order, ok=fin["err_a"] <= 0.02 and fin["err_b"] <= 0.02, ok_order=order >= 0.8)


def check2():
    b, h = 1e-3, 0.02e-3
    T = lambda r: 300.0 + 1000.0 * (r / b) ** 2  # noqa: E731
    ex = v.exact.disc(T, b, v.E, v.NU, v.AL, T_ref=300.0)
    nodes, tri = v.grid(400, 4, b, h)
    k = int(np.argmin(np.hypot(nodes[:, 0] - h / 2, nodes[:, 1])))
    _, sig, pos, _ = solve_quadratic(nodes, tri, T(nodes[:, 1]), v.P, bc="free", z_fix_node=k)
    rows = []
    for frac in (0.25, 0.5, 0.9):
        r0 = frac * b
        sel = np.abs(pos[:, 1] - r0) <= b / 800
        s = float(np.mean(sig[sel, 2]))
        rc = float(np.mean(pos[sel, 1]))
        e = ex(rc)[1]
        rows.append(dict(r_over_b=frac, r=rc, hoop=s, exact=e, err=abs(s / e - 1)))
    return dict(points=rows, ok=all(r["err"] <= 0.02 for r in rows))


def main():
    c0, c1, c2 = check0(), check1(), check2()
    passed = c0["ok"] and c1["ok"] and c1["ok_order"] and c2["ok"]
    (HERE / "v11_quadratic_results.json").write_text(json.dumps(dict(passed=passed, check0=c0, check1=c1, check2=c2), indent=1, default=float))
    L = ["# V11 with quadratic triangles (amendment 1)", "", f"Check 0: stress {c0['stress_rel']:.1e}, displacement {c0['disp_rel']:.1e} "
         f"({'meets' if c0['ok'] else 'MISSES'} 1e-9).", "", "| Mesh | Inner hoop (MPa) | Error | Outer hoop (MPa) | Error |", "|---|---|---|---|---|"]
    L += [f"| {r['mesh']} | {r['hoop_a'] / 1e6:.1f} | {100 * r['err_a']:.2f} % | {r['hoop_b'] / 1e6:.1f} | {100 * r['err_b']:.2f} % |" for r in c1["meshes"]]
    L += ["", f"Finest {'meets' if c1['ok'] else 'MISSES'} 2 %; order {c1['order']:.2f} ({'meets' if c1['ok_order'] else 'MISSES'} 0.8).", "",
          "| r/b | Hoop (MPa) | Exact (MPa) | Error |", "|---|---|---|---|"]
    L += [f"| {r['r_over_b']:g} | {r['hoop'] / 1e6:.1f} | {r['exact'] / 1e6:.1f} | {100 * r['err']:.2f} % |" for r in c2["points"]]
    L += ["", f"**V11 with quadratic triangles {'PASSES' if passed else 'FAILS'}.**"]
    (HERE / "report_quadratic.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
