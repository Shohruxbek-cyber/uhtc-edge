"""V12: comparison of the three flank-heating models (plan: validation/V12_plan.md, committed before this run).

Usage: python validation/V12_flank/run_v12.py   -> v12_results.json, report.md, figures/v3_flank_models.png
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "runs"))
import common  # noqa: E402
from uhtc import heating as ht, materials as mt  # noqa: E402
from uhtc.conduction2d import NoseGeometry, Solver2D, integral_distribution  # noqa: E402
from uhtc.coupled import FlowState  # noqa: E402

MODELS = ("lees", "tauber", "integral")
NAMES = {"lees": "baseline", "tauber": "sharpcone", "integral": "integral"}


def overlay():
    rows = []
    for R in (1e-3, 10e-3):
        g = NoseGeometry(R_n=R, theta_c=math.radians(5.0), L_flank=0.03, axisymmetric=True)
        sol = Solver2D(g, mt.composite(0.2), 300.0, mesh_kw=common.MESH_2D)
        s_t = R * (math.pi / 2 - g.theta_c)
        pts = {"tangency": s_t, "flank 1/4": s_t + 0.0075, "flank 1/2": s_t + 0.015, "flank end": s_t + 0.03}
        for M in (6.0, 9.0, 12.0):
            for z in (20e3, 30e3, 40e3):
                st = common.stag_state(M, z)
                flow = FlowState(st.H_s, st.p_s, st.p_O2_w, st.h_D, 1.0, 1.0, 1.0, M, z)
                q_s_cold = ht.tauber_sphere(st.rho_inf, st.V, R, 0.0, 1.0)
                gs = {m: sol.heating_shape(st.V, st.rho_inf, q_s_cold, dist=m, flow=flow) for m in MODELS}
                for name, sp in pts.items():
                    rows.append(dict(R_mm=R * 1e3, mach=M, alt_km=z / 1e3, point=name,
                                     **{NAMES[m]: float(np.interp(sp, sol.s, gs[m])) for m in MODELS}))
    return rows


def anchor():
    """Cleary 1969 15 deg sphere-cone, M 10.6, R = 1.1 in, R/r_b = 0.183 (the V5 points): measured q/q_s against the models."""
    v5 = [r for r in csv.DictReader(open(ROOT / "validation/V5_heating_distribution/v5_results.csv", encoding="utf-8"))
          if r["case"].startswith("Cleary")]
    M, T0, G, cp, Rg = 10.6, 2000.0 / 1.8, 1.4, 1004.5, 287.05
    T_inf = T0 / (1 + 0.5 * (G - 1) * M ** 2)
    V = M * math.sqrt(G * Rg * T_inf)
    p0 = 1200 * 6894.757  # 1200 psia (the middle Reynolds number); the shape does not depend on the pressure level
    p_inf = p0 / (1 + 0.5 * (G - 1) * M ** 2) ** (G / (G - 1))
    p_s = p_inf * ((G + 1) ** 2 * M ** 2 / (4 * G * M ** 2 - 2 * (G - 1))) ** (G / (G - 1)) * (1 - G + 2 * G * M ** 2) / (G + 1)
    rho_inf = p_inf / (Rg * T_inf)
    R = 1.1 * 0.0254
    g = NoseGeometry(R_n=R, theta_c=math.radians(15.0), L_flank=0.5, axisymmetric=True)
    s = np.array([float(r["s_over_R"]) for r in v5]) * R
    qi = integral_distribution(g, s, M, V, cp * T0, p_s, p_inf, n_fine=800, T_wall=0.3 * T0, perfect_gas=True)
    s_t = R * (math.pi / 2 - g.theta_c)
    qc = ht.tauber_cone(rho_inf, V, s, g.theta_c, 0.0, 1.0) / ht.tauber_sphere(rho_inf, V, R, 0.0, 1.0)
    out = []
    for r, a, c in zip(v5, qi, qc):
        meas, base = float(r["measured"]), float(r["model"])
        out.append(dict(s_over_R=float(r["s_over_R"]), Re_per_ft=r["Re_per_ft"], measured=meas, baseline=base, integral=float(a),
                        sharpcone=float(c), on_flank=bool(float(r["s_over_R"]) * R >= s_t)))
    return out


def main():
    ov, an = overlay(), anchor()
    fl = [r for r in ov if r["point"] != "tangency"]
    rat = {m: [r[m] / r["baseline"] for r in fl] for m in ("sharpcone", "integral")}
    meas = [a["measured"] / a["baseline"] for a in an if a["on_flank"]]
    L = min(min(v) for v in rat.values()) if rat else math.nan
    H = max(max(v) for v in rat.values())
    La, Ha = min(min(v) for v in [rat["sharpcone"], rat["integral"], meas]), max(max(v) for v in [rat["sharpcone"], rat["integral"], meas])
    narrow = La >= 0.67 and Ha <= 1.5
    rec = (dict(narrow=True, band=[math.floor(La * 20) / 20, math.ceil(Ha * 20) / 20]) if narrow
           else dict(narrow=False, band=None, reason=[k for k, v in (("sharpcone", rat["sharpcone"]), ("integral", rat["integral"]),
                                                                      ("measured anchor", meas)) if min(v) < 0.67 or max(v) > 1.5]))
    out = dict(overlay=ov, anchor=an, ratios_to_baseline={k: [min(v), max(v)] for k, v in {**rat, "measured": meas}.items()},
               recommendation=rec)
    (HERE / "v12_results.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    Lr = ["# V12: the three flank-heating models", "", "Plan: `validation/V12_plan.md` (committed before this run).", "",
          "## Overlay at the map corners, q/q_s", "", "| R (mm) | Mach | Altitude (km) | Point | Baseline (Lees) | Sharp cone | Integral |",
          "|---|---|---|---|---|---|---|"]
    Lr += [f"| {r['R_mm']:g} | {r['mach']:g} | {r['alt_km']:g} | {r['point']} | {r['baseline']:.4f} | {r['sharpcone']:.4f} | {r['integral']:.4f} |"
           for r in ov]
    Lr += ["", "## Measured anchor: Cleary 1969, 15 deg sphere-cone, M 10.6 (the V5 points)", "",
           "| s/R | Re/ft | Measured | Baseline | Sharp cone | Integral |", "|---|---|---|---|---|---|"]
    Lr += [f"| {a['s_over_R']:.3f} | {a['Re_per_ft']} | {a['measured']:.4f} | {a['baseline']:.4f} | {a['sharpcone']:.4f} | {a['integral']:.4f} |"
           for a in an]
    Lr += ["", "## Ratios to the baseline on the flank (corner points) and of the measured anchor", ""]
    Lr += [f"- {k}: {v[0]:.3f} to {v[1]:.3f}" for k, v in out["ratios_to_baseline"].items()]
    Lr += ["", "## Recommendation (rule of the plan, item 4)", "",
           (f"The flank band can be narrowed to x{rec['band'][0]:.2f}-{rec['band'][1]:.2f}." if rec["narrow"] else
            f"The flank band cannot be narrowed: outside 0.67-1.5: {', '.join(rec['reason'])}.")]
    (HERE / "report.md").write_text("\n".join(Lr) + "\n", encoding="utf-8")
    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    for k, R in enumerate((1.0, 10.0)):
        sel = [r for r in ov if r["R_mm"] == R and r["point"] != "tangency"]
        for m, mk in (("sharpcone", "o"), ("integral", "s")):
            ax[k].scatter([r["mach"] + (0.15 if m == "integral" else -0.15) for r in sel], [r[m] / r["baseline"] for r in sel], marker=mk,
                          facecolors="none", edgecolors="k" if m == "integral" else "0.5", label=f"{m} / baseline")
        ax[k].axhline(1.0, color="k", lw=0.8)
        for y in (0.5, 2.0):
            ax[k].axhline(y, color="k", lw=0.8, ls=":")
        ax[k].set_yscale("log")
        ax[k].set_xlabel("Mach number")
        ax[k].set_title(f"R_n = {R:g} mm, flank points, 20-40 km")
    ax[0].set_ylabel("q/q_s relative to the baseline")
    ax[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(ROOT / "figures/v3_flank_models.png", dpi=200)
    print("\n".join(Lr[-8:]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
