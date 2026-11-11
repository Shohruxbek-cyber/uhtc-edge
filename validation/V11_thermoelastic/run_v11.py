"""V11: verification of the thermoelastic module against exact solutions (plan: validation/V11_plan.md, committed before this run).

Usage: python validation/V11_thermoelastic/run_v11.py   -> v11_results.json, report.md
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))
from uhtc.conduction2d import _assemble  # noqa: E402
from uhtc.thermoelastic import Props, solve  # noqa: E402
import exact  # noqa: E402

E, NU, AL = 400e9, 0.14, 6e-6
P = Props(E, AL, NU)


def grid(nr, nz, R, H, jitter=0.0, seed=0):
    """Structured rectangle r in [0, R], z in [0, H] split into triangles; interior nodes optionally jittered."""
    r, z = np.meshgrid(np.linspace(0, R, nr + 1), np.linspace(0, H, nz + 1), indexing="ij")
    if jitter:
        rng = np.random.default_rng(seed)
        inner = (r > 0) & (r < R) & (z > 0) & (z < H)
        r = r + inner * rng.uniform(-jitter, jitter, r.shape) * R / nr
        z = z + inner * rng.uniform(-jitter, jitter, z.shape) * H / nz
    nodes = np.c_[z.ravel(), r.ravel()]  # [x axial, y radius]
    idx = np.arange((nr + 1) * (nz + 1)).reshape(nr + 1, nz + 1)
    a, b, c, d = idx[:-1, :-1].ravel(), idx[1:, :-1].ravel(), idx[1:, 1:].ravel(), idx[:-1, 1:].ravel()
    return nodes, np.r_[np.c_[a, b, c], np.c_[a, c, d]]


def check0():
    nodes, tri = grid(4, 3, 1e-3, 0.7e-3, jitter=0.3)
    dT = 1000.0
    k = int(np.argmin(nodes[:, 0] + nodes[:, 1]))
    U, sig, _ = solve(nodes, tri, np.full(len(nodes), 300.0 + dT), P, bc="free", z_fix_node=k)
    s_err = float(np.abs(sig).max() / (E * AL * dT))
    ue = np.c_[AL * dT * (nodes[:, 0] - nodes[k, 0]), AL * dT * nodes[:, 1]]
    u_err = float(np.abs(U - ue).max() / (AL * dT * 1e-3))
    return dict(stress_rel=s_err, disp_rel=u_err, ok=s_err < 1e-9 and u_err < 1e-9)


def sphere_mesh(n_arc, n_eta, a, b):
    phi = np.linspace(0.0, np.pi, n_arc + 1)
    xs, ys = b * np.cos(phi), b * np.sin(phi)
    xe, ye = a * np.cos(phi), a * np.sin(phi)
    ys[[0, -1]] = 0.0
    ye[[0, -1]] = 0.0
    nodes, tri, _ = _assemble(xs, ys, xe, ye, np.linspace(0.0, 1.0, n_eta + 1), b)
    return nodes, tri


def surface_hoop(cent, sig, a, b, side):
    """Hoop stress at the inner (side 'a') or outer ('b') surface on the equator: linear extrapolation in radius from the two
    element centroids nearest the surface within the equatorial band."""
    rr = np.hypot(cent[:, 0], cent[:, 1])
    band = np.abs(cent[:, 0]) <= 0.02 * b
    i = np.where(band)[0]
    order = i[np.argsort(rr[i])] if side == "a" else i[np.argsort(-rr[i])]
    # two nearest distinct radii
    r1, s1 = rr[order[0]], sig[order[0], 2]
    j = next(k for k in order[1:] if abs(rr[k] - r1) > 1e-6 * b)
    r2, s2 = rr[j], sig[j, 2]
    r0 = a if side == "a" else b
    return s1 + (s2 - s1) * (r0 - r1) / (r2 - r1)


def check1():
    a, b = 0.5e-3, 1e-3
    T = lambda r: 300.0 + 1000.0 * (1 / r - 1 / b) / (1 / a - 1 / b)  # noqa: E731
    ex = exact.sphere(T, a, b, E, NU, AL, T_ref=300.0)
    ea, eb = ex(a)[1], ex(b)[1]
    rows = []
    for n_arc, n_eta in ((40, 8), (80, 16), (160, 32), (320, 64)):
        nodes, tri = sphere_mesh(n_arc, n_eta, a, b)
        Tn = T(np.hypot(nodes[:, 0], nodes[:, 1]))
        k = int(np.argmax(nodes[:, 0]))
        _, sig, cent = solve(nodes, tri, Tn, P, bc="free", z_fix_node=k)
        sa, sb = surface_hoop(cent, sig, a, b, "a"), surface_hoop(cent, sig, a, b, "b")
        rows.append(dict(mesh=f"{n_arc}x{n_eta}", h=b / n_eta, hoop_a=sa, hoop_b=sb, err_a=abs(sa / ea - 1), err_b=abs(sb / eb - 1)))
    h = np.array([r["h"] for r in rows[-3:]])
    e = np.array([r["err_a"] for r in rows[-3:]])
    order = float(np.polyfit(np.log(h), np.log(e), 1)[0])
    pair = [float(math.log(rows[i]["err_a"] / rows[i + 1]["err_a"]) / math.log(2)) for i in range(1, 3)]
    fin = rows[-1]
    return dict(exact_a=ea, exact_b=eb, meshes=rows, order=order, order_pairs=pair,
                ok=fin["err_a"] <= 0.02 and fin["err_b"] <= 0.02, ok_order=order >= 0.8)


def check2():
    b, h = 1e-3, 0.02e-3
    T = lambda r: 300.0 + 1000.0 * (r / b) ** 2  # noqa: E731
    ex = exact.disc(T, b, E, NU, AL, T_ref=300.0)
    nodes, tri = grid(400, 4, b, h)
    k = int(np.argmin(np.hypot(nodes[:, 0] - h / 2, nodes[:, 1])))
    _, sig, cent = solve(nodes, tri, T(nodes[:, 1]), P, bc="free", z_fix_node=k)
    rows = []
    for frac in (0.25, 0.5, 0.9):
        r0 = frac * b
        sel = np.argsort(np.abs(cent[:, 1] - r0) + np.abs(cent[:, 0] - h / 2) * 0.0)[:16]
        # average the hoop stress of the elements in the cell column at r0 (mid-plane average over the thin thickness)
        col = sel[np.abs(cent[sel, 1] - r0) <= 1.01 * np.abs(cent[sel, 1] - r0).min() + b / 400]
        s = float(np.mean(sig[col, 2]))
        rc = float(np.mean(cent[col, 1]))
        e = ex(rc)[1]
        rows.append(dict(r_over_b=frac, r=rc, hoop=s, exact=e, err=abs(s / e - 1)))
    return dict(points=rows, ok=all(r["err"] <= 0.02 for r in rows))


def main():
    c0, c1, c2 = check0(), check1(), check2()
    passed = c0["ok"] and c1["ok"] and c1["ok_order"] and c2["ok"]
    out = dict(passed=passed, check0=c0, check1=c1, check2=c2)
    (HERE / "v11_results.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    L = ["# V11: verification of the thermoelastic module", "", "Plan: `validation/V11_plan.md` (committed before this run).", "",
         f"**Check 0, patch test:** largest element stress {c0['stress_rel']:.1e} E alpha dT, displacement error {c0['disp_rel']:.1e}: "
         f"{'meets' if c0['ok'] else 'MISSES'} 1e-9.", "",
         f"**Check 1, thick-walled sphere:** exact hoop stress {c1['exact_a'] / 1e6:.1f} MPa (inner) and {c1['exact_b'] / 1e6:.1f} MPa (outer).", "",
         "| Mesh | Inner hoop (MPa) | Error | Outer hoop (MPa) | Error |", "|---|---|---|---|---|"]
    L += [f"| {r['mesh']} | {r['hoop_a'] / 1e6:.1f} | {100 * r['err_a']:.2f} % | {r['hoop_b'] / 1e6:.1f} | {100 * r['err_b']:.2f} % |" for r in c1["meshes"]]
    L += ["", f"Finest mesh {'meets' if c1['ok'] else 'MISSES'} 2 %. Observed order of the inner-surface error over the last three meshes "
          f"{c1['order']:.2f} (pairs {c1['order_pairs'][0]:.2f}, {c1['order_pairs'][1]:.2f}): {'meets' if c1['ok_order'] else 'MISSES'} 0.8.", "",
          "**Check 2, free thin disc (h/b = 0.02, plane-stress solution):**", "", "| r/b | Hoop (MPa) | Exact (MPa) | Error |", "|---|---|---|---|"]
    L += [f"| {r['r_over_b']:g} | {r['hoop'] / 1e6:.1f} | {r['exact'] / 1e6:.1f} | {100 * r['err']:.2f} % |" for r in c2["points"]]
    L += ["", f"**V11 {'PASSES' if passed else 'FAILS'}.**"]
    (HERE / "report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
