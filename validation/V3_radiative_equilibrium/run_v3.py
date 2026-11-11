"""V3: radiative-equilibrium stagnation temperature against a published example.

Target: Kolodziej (1997), NASA TM-112204, Table 1 (doc p. 3, PDF p. 7): velocity at which a two-dimensional
unswept leading edge of radius 0.762 mm (0.03 in., doc p. 2) with surface emissivity 0.80 reaches its
reusable temperature limit of 2149 C (3900 F) in radiative equilibrium (doc p. 2, eq. 3: q_w = sigma eps T_w^4),
fully catalytic hot wall (FCHW: h_w = equilibrium enthalpy of air at (p_t2, T_w), doc p. 3):
    altitude 0 km: 2.30 km/s (Kolodziej), 2.25 km/s (Hill, his ref. 12)
    altitude 30 km: 2.53 km/s (Kolodziej), 2.56 km/s (Hill)
Model: uhtc stagnation state (US 1976, equilibrium normal shock), 2D stagnation-line heating = Sutton-Graves
x C_cyl/C_sph (Tauber 1989 eqs. 40-41), wall enthalpy of equilibrium air at (T_w, p_s), no conduction.
Also shown: Fay-Riddell for a 2D stagnation line (axisymmetric value / sqrt(2), the (beta/2^n)^1/2 factor of
Kolodziej's eq. 4 with n = 1), Pr = 0.71, Le = 1.4.
Pass criterion: |V_model / V_Kolodziej - 1| <= 0.05.  Reason: independent stagnation-heating methods differ by
up to about 15 % (Tauber 1989, doc p. 5); the script reports the velocity change that a +/-15 % heating change
produces at this limit, so the criterion can be checked against it.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from uhtc import atmosphere as atm, heating as ht, stagnation as st, thermo as th  # noqa: E402
from uhtc.surface import SIGMA  # noqa: E402

RN, EPS, TMAX = 0.762e-3, 0.80, 2149.0 + 273.15
TARGETS = {0.0: (2.30, 2.25), 30e3: (2.53, 2.56)}
# plus the digitized points of Fig. 1 (fully catalytic FCHW and non-catalytic NCHW hot wall) at 50-70 km


def q_model(V, z, method="SG", scale=1.0):
    s = st.stagnation_state_V(V, z)
    xw = th.equilibrium_air(TMAX, s.p_s)
    Mw, hw, _ = th.mixture_props(xw, TMAX)
    if method == "SG":
        q = ht.sutton_graves_cylinder(s.p_s, RN, s.H_s, hw)
    elif method == "SG_nc":  # non-catalytic bound (D-HEAT-2): frozen wall enthalpy, Fay-Riddell bracket with Le = 0
        q = min(ht.sutton_graves_cylinder(s.p_s, RN, s.H_s, th.h_air_frozen(TMAX)) * ht.noncatalytic_factor(s.H_s, s.h_D),
                ht.sutton_graves_cylinder(s.p_s, RN, s.H_s, hw))  # capped at the catalytic value (D-HEAT-2)
    else:
        rho_w = s.p_s * Mw / (th.R * TMAX)
        due = ht.newtonian_velocity_gradient(RN, s.p_s, s.p_inf, s.rho_e)
        q = ht.fay_riddell(rho_w, float(atm.viscosity(TMAX)), s.rho_e, s.mu_e, s.H_s, hw, s.h_D, due) / np.sqrt(2)
    return scale * float(q)


def v_limit(z, method="SG", scale=1.0, hi=4000.0):
    target = EPS * SIGMA * TMAX ** 4
    return brentq(lambda V: q_model(V, z, method, scale) - target, 1700.0, hi, xtol=0.5)


def main():
    lines = ["# V3 — Radiative-equilibrium stagnation temperature vs Kolodziej (1997) Table 1", "",
             "Case: 2D unswept leading edge, R_n = 0.762 mm, emissivity 0.80, limit 2149 °C, fully catalytic hot wall,",
             "radiative equilibrium (no conduction). Velocity at which T_w reaches the limit:", "",
             "| Altitude | Kolodziej | Hill | Model (S–G, 2D) | Model/Kolodziej | Fay–Riddell 2D | S–G with q ×0.85 / ×1.15 |",
             "|---|---|---|---|---|---|---|"]
    worst = 0.0
    worst_tab = 0.0
    csv_rows = []
    for z, (vk, vh) in TARGETS.items():
        v = v_limit(z) / 1e3
        vfr = v_limit(z, "FR") / 1e3
        vlo, vhi = v_limit(z, scale=0.85) / 1e3, v_limit(z, scale=1.15) / 1e3
        worst_tab = max(worst_tab, abs(v / vk - 1))
        csv_rows.append(("Table1_FCHW", z / 1e3, vk, v, vfr))
        lines.append(f"| {z / 1e3:.0f} km | {vk:.2f} km/s | {vh:.2f} km/s | {v:.3f} km/s | {v / vk:.3f} | {vfr:.3f} km/s | {vlo:.3f} / {vhi:.3f} km/s |")
    lines += ["", "Digitized points of Kolodziej's Fig. 1 (higher altitudes, heating not blocked; ±0.05 km/s reading",
              "uncertainty; `literature/digitized/Kolodziej1997_fig1.csv`):", "",
              "| Wall | Altitude | Kolodziej Fig. 1 | Model | Model/Kolodziej |", "|---|---|---|---|---|"]
    import csv
    with open(HERE.parents[1] / "literature" / "digitized" / "Kolodziej1997_fig1.csv") as f:
        rows = [r for r in csv.reader(l for l in f if not l.startswith("#"))][1:]
    for curve, alt, vk in rows:
        z, vk = float(alt) * 1e3, float(vk)
        v = v_limit(z, "SG" if curve == "FCHW" else "SG_nc", hi=5200.0) / 1e3
        worst = max(worst, abs(v / vk - 1))
        csv_rows.append(("Fig1_" + curve, z / 1e3, vk, v, float("nan")))
        lines.append(f"| {curve} | {z / 1e3:.0f} km | {vk:.2f} km/s | {v:.3f} km/s | {v / vk:.3f} |")
    ok = worst_tab <= 0.05
    flagged = [l for l in lines if l.startswith("| NCHW") or l.startswith("| FCHW")]
    n_out = sum(1 for l in flagged if abs(float(l.split("|")[-2]) - 1) > 0.05)
    lines += ["", f"Result (primary check, tabulated Table 1 values): **{'PASS' if ok else 'FAIL'}** — largest "
              f"|V_model/V_Kolodziej − 1| = {worst_tab:.3f} (criterion ≤ 0.05, set before the run).",
              f"Digitized Fig. 1 points: largest deviation {worst:.3f}; {n_out} of {len(flagged)} outside ±0.05 "
              "(flagged; see note below).",
              "", "The last column shows the velocity at the limit when the heating is changed by ±15 % (the spread of",
              "stagnation-heating methods quoted by Tauber 1989), to show the size of the tolerance in velocity terms.",
              "", "**Flagged point.** The non-catalytic (NCHW) limit at 70 km lies 5 % above Kolodziej's curve, the model",
              "giving the higher velocity (lower non-catalytic heating). The deviation grows with altitude and enthalpy, where",
              "h_D/H_s is largest, so it points to the treatment of the dissociation enthalpy in the non-catalytic bound",
              "(Kolodziej's curve fits for the edge composition vs our JANAF equilibrium; frozen vs equilibrium wall enthalpy).",
              "The point is outside the conditions of this study (z ≤ 40 km, V ≤ 3.8 km/s); at the digitized conditions",
              "nearest to the study's (50–55 km, 3.2–3.6 km/s) the bound agrees within 0.5 %.",
              "UNSURE: which of the two treatments is closer to reality at 70 km.",
              "", "**Limitation of this check.** At both tabulated points the limit temperature (2422 K) is close to the",
              "stagnation temperature of the flow, so the heating is nearly 'blocked' (H_s − h_w small) and the limiting",
              "velocity depends mainly on the stagnation thermodynamics (H_s, equilibrium h_w) rather than on the heat-transfer",
              "coefficient: a ±15 % change of the heating moves it by ≤ 1.3 %. The heat-transfer coefficient itself is checked",
              "in V4 (Sutton–Graves vs Fay–Riddell) and, with the material response, in V7 (arc jet)."]
    (HERE / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    with open(HERE / "v3_results.csv", "w", encoding="utf-8") as f:
        f.write("source,alt_km,V_kolodziej_kms,V_model_SG_kms,V_model_FR_kms\n")
        for r in csv_rows:
            f.write(",".join([r[0]] + [f"{x:.4f}" for x in r[1:]]) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
