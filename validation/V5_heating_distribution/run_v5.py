"""V5: heating distribution around the nose (revision, D-2D-7) against independent open results.

(1) Axisymmetric: Cleary 1969 (NASA TN D-5450) measured laminar heating on a 15 deg sphere-cone, M 10.6, cold wall,
    R/r_b = 0.183 (Table III, transcribed: literature/digitized/Cleary1969_TablesI-III_alpha0.csv), normalised by the
    model's own stagnation value.  Model: baseline distribution (uhtc.conduction2d.baseline_distribution) for the same
    geometry, perfect gas (gamma 1.4) wind-tunnel stream.  NOTE: the band of D-2D-7 (x0.66-x1.5) was derived from these
    data, so (1) is a consistency check of the band, not an independent validation.
(2) Planar: Kolodziej et al. 1997 (NTRS 19970040120) Fig. 2, steady radiative-equilibrium temperatures of a 5 deg
    cylinder-wedge leading edge (R = 2.54 mm, Mach 15, 30.48 km) from two independent methods (BLIMPK boundary-layer code;
    PERFORM engineering method), DIGITIZED (literature/digitized/Kolodziej1997conf_fig2.csv).  Compared on a hot-wall
    basis: q(s)/q(0) = (T(s)/T0)^4 for a constant emissivity, against the model's g(s) (H_s - h_w(T(s)))/(H_s - h_w(T0))
    with the project's equilibrium stagnation state and wall enthalpy (fully catalytic).  Independent of the band.
Criterion (set before the run): (1) every measured point inside the band x0.66-x1.5 around the baseline;
(2) reported, not judged (single computed case outside the study's Mach range).
Output: report.md, v5_results.csv
"""
import csv
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "runs"))
from uhtc import atmosphere as atm  # noqa: E402
from uhtc.conduction2d import PHI_COS, NoseGeometry, baseline_distribution, lees_distribution  # noqa: E402
from uhtc.coupled import WallEnthalpy  # noqa: E402
import common  # noqa: E402

DIG = ROOT / "literature" / "digitized"
G = 1.4


def pitot(M):
    return ((G + 1) ** 2 * M ** 2 / (4 * G * M ** 2 - 2 * (G - 1))) ** (G / (G - 1)) * (1 - G + 2 * G * M ** 2) / (G + 1)


def rows(f):
    return [r for r in csv.DictReader(l for l in open(f, encoding="utf-8") if not l.startswith("#"))]


def main():
    out, lines = [], ["# V5 — heating distribution around the nose (D-2D-7)", ""]
    # (1) Cleary 15 deg sphere-cone, M 10.6 (all points lie on the flank: s/R >= 1.31 = tangency)
    M = 10.6
    geo = NoseGeometry(R_n=1.0, theta_c=math.radians(15.0), L_flank=30.0, axisymmetric=True)
    pts = [(float(r["s_over_R"]), float(r["q_over_qstag_own"]), r["Re_per_ft"]) for r in rows(DIG / "Cleary1969_TablesI-III_alpha0.csv")
           if r["R_over_rb"] == "0.1833" and r["q_over_qstag_own"]]
    s = np.array([p[0] for p in pts])
    phi = np.minimum(s, math.pi / 2 - geo.theta_c)
    base = baseline_distribution(geo, s, phi, s > (math.pi / 2 - geo.theta_c), M, M, 0.5 * M ** 2 * (1 + 2 / ((G - 1) * M ** 2)),
                                 pitot(M), 1.0)
    inside = 0
    for (sr, q, re), b in zip(pts, base):
        r = q / b
        ok = 0.66 <= r <= 1.5
        inside += ok
        out.append(dict(case="Cleary1969 15deg M10.6", Re_per_ft=re, s_over_R=sr, measured=q, model=round(float(b), 4),
                        ratio=round(r, 3), in_band=ok))
    # continuous baseline curve for the figure (same geometry and stream); band factors apply beyond PHI_COS
    sc = np.linspace(0.0, 17.0, 341)
    phic = np.minimum(sc, math.pi / 2 - geo.theta_c)
    bc = baseline_distribution(geo, sc, phic, sc > (math.pi / 2 - geo.theta_c), M, M,
                               0.5 * M ** 2 * (1 + 2 / ((G - 1) * M ** 2)), pitot(M), 1.0)
    with open(HERE / "v5_curve.csv", "w", newline="", encoding="utf-8") as f:
        f.write("# baseline distribution q/q0 for Cleary's 15 deg sphere-cone at M 10.6; band x0.66-x1.5 beyond "
                f"phi = {math.degrees(PHI_COS):.0f} deg (s/R = {PHI_COS:.4f})\n")
        f.write("s_over_R,q_over_q0\n")
        for a, b in zip(sc, bc):
            f.write(f"{a:.4f},{float(b):.5f}\n")
    rat = [o["ratio"] for o in out]
    lines += [f"(1) Cleary 1969, 15 deg sphere-cone, M 10.6 (measured, {len(rat)} points, s/R {min(s):.2f}-{max(s):.2f}): "
              f"measured/baseline = {min(rat):.2f}-{max(rat):.2f}; {inside} of {len(rat)} inside the band x0.66-x1.5 "
              f"({'consistent' if inside == len(rat) else 'NOT consistent'}; the band was derived from these data).", ""]
    # (2) Kolodziej 1997 planar 5 deg, M 15, 30.48 km, R = 2.54 mm
    Mk, zk, Rk = 15.0, 30480.0, 2.54e-3
    st = common.stag_state(Mk, zk)
    T_inf, p_inf, rho_inf = atm.state(zk)
    V = Mk * math.sqrt(atm.GAMMA * atm.R_AIR * float(T_inf))
    geo2 = NoseGeometry(R_n=Rk, theta_c=math.radians(5.0), L_flank=0.05, axisymmetric=False)
    K = rows(DIG / "Kolodziej1997conf_fig2.csv")
    hw = WallEnthalpy(st.p_s, True)
    T0 = float(K[0]["T_BLIMPK_K"])
    sK = np.array([float(r["s_over_R"]) for r in K]) * Rk
    phi2 = np.minimum(sK / Rk, math.pi / 2 - geo2.theta_c)
    g2 = baseline_distribution(geo2, sK, phi2, sK / Rk > (math.pi / 2 - geo2.theta_c), Mk, V, st.H_s, st.p_s, float(p_inf))
    lines.append(f"(2) Kolodziej et al. 1997 Fig. 2, planar 5 deg, M 15, 30.48 km, R 2.54 mm (computed by BLIMPK and "
                 f"PERFORM, radiative equilibrium; hot-wall ratios, constant emissivity):")
    lines.append("")
    lines.append("| s/R | BLIMPK q/q0 | PERFORM q/q0 | model q/q0 at BLIMPK T(s) | model / BLIMPK | model / PERFORM |")
    lines.append("|---|---|---|---|---|---|")
    rb, rp = [], []
    for r, gg in zip(K[1:], g2[1:]):
        Tb, Tp = float(r["T_BLIMPK_K"]), float(r["T_PERFORM_K"])
        qb, qp = (Tb / T0) ** 4, (Tp / T0) ** 4
        mb = gg * (st.H_s - hw(Tb)[0]) / (st.H_s - hw(T0)[0])
        mp = gg * (st.H_s - hw(Tp)[0]) / (st.H_s - hw(T0)[0])
        rb.append(mb / qb)
        rp.append(mp / qp)
        lines.append(f"| {float(r['s_over_R']):.2f} | {qb:.3f} | {qp:.3f} | {mb:.3f} | {mb / qb:.2f} | {mp / qp:.2f} |")
        out.append(dict(case="Kolodziej1997 planar 5deg M15", Re_per_ft="", s_over_R=float(r["s_over_R"]), measured=round(qb, 4),
                        model=round(float(mb), 4), ratio=round(mb / qb, 3), in_band=""))
    lines += ["", f"Model / BLIMPK = {min(rb):.2f}-{max(rb):.2f}; model / PERFORM = {min(rp):.2f}-{max(rp):.2f} over s/R "
              f"{float(K[1]['s_over_R']):.1f}-{float(K[-1]['s_over_R']):.1f}. Below 1 means the model's planar flank heating is "
              "lower than these methods (the wing-edge results are then not conservative)."]
    with open(HERE / "v5_results.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    (HERE / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
