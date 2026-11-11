"""V4: stagnation-point heating, Sutton-Graves (baseline) vs Fay-Riddell vs Tauber's simplified sphere formula.

Fay-Riddell (Hamilton 1982 eq. 6) is evaluated with: edge = equilibrium air at (p_s, H_s); wall = equilibrium air
at (T_w, p_s) (fully catalytic wall); viscosity from Sutherland's law (US 1976 eq. 51); Pr = 0.71, Le = 1.4;
Newtonian velocity gradient with rho_s = rho_e (Tauber 1989 eq. 39a).  Both correlations scale as R_n^-1/2, so
the ratio does not depend on the nose radius (R_n = 3 mm is used for the absolute values).
Pass criterion: |q_SG / q_FR - 1| <= 0.15, the spread among independent stagnation-heating calculations quoted
by Tauber (1989, doc p. 5: "maximum difference ... is about 15%").
Output: validation/V4_heating/v4_results.csv and report.md (written by this script).
"""
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

from uhtc import atmosphere as atm, heating as ht, stagnation as st, thermo as th  # noqa: E402

RN = 0.003
CASES = [(m, 30e3) for m in (5, 6, 8, 10, 12)] + [(10, 20e3), (10, 40e3), (8, 25e3), (6, 35e3)]
TW = (300.0, 1111.0, 2000.0)


def main():
    rows = []
    for mach, z in CASES:
        s = st.stagnation_state(mach, z)
        due = ht.newtonian_velocity_gradient(RN, s.p_s, s.p_inf, s.rho_e)
        for Tw in TW:
            if Tw > 0.95 * s.T_e:
                continue
            xw = th.equilibrium_air(Tw, s.p_s)
            Mw, hw, _ = th.mixture_props(xw, Tw)
            rho_w = s.p_s * Mw / (th.R * Tw)
            mu_w = float(atm.viscosity(Tw))
            q_sg = float(ht.sutton_graves(s.p_s, RN, s.H_s, hw))
            q_fr = float(ht.fay_riddell(rho_w, mu_w, s.rho_e, s.mu_e, s.H_s, hw, s.h_D, due))
            q_ta = float(ht.tauber_sphere(s.rho_inf, s.V, RN, hw, s.H_s))
            rows.append(dict(mach=mach, alt_km=z / 1e3, T_w=Tw, H_s_MJkg=s.H_s / 1e6, p_s_kPa=s.p_s / 1e3,
                             T_e=s.T_e, hD_over_Hs=s.h_D / s.H_s, q_SG_MWm2=q_sg / 1e6, q_FR_MWm2=q_fr / 1e6,
                             q_Tauber_MWm2=q_ta / 1e6, SG_over_FR=q_sg / q_fr, Tauber_over_FR=q_ta / q_fr))
    with open(HERE / "v4_results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()})
    ratios = [r["SG_over_FR"] for r in rows]
    worst = max(abs(x - 1) for x in ratios)
    ok = worst <= 0.15
    lines = ["# V4 — Stagnation-point heating: Sutton–Graves vs Fay–Riddell", "",
             f"Result: **{'PASS' if ok else 'FAIL'}** — largest |q_SG/q_FR − 1| = {worst:.3f} (criterion ≤ 0.15).", "",
             "Method: see `run_v4.py` header. R_n = 3 mm (ratios are independent of R_n).", "",
             "| Mach | z km | T_w K | H_s MJ/kg | p_s kPa | T_e K | h_D/H_s | q_SG MW/m² | q_FR MW/m² | q_Tauber MW/m² | SG/FR | Tauber/FR |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['mach']} | {r['alt_km']:.0f} | {r['T_w']:.0f} | {r['H_s_MJkg']:.2f} | {r['p_s_kPa']:.1f} | "
                     f"{r['T_e']:.0f} | {r['hD_over_Hs']:.3f} | {r['q_SG_MWm2']:.3f} | {r['q_FR_MWm2']:.3f} | "
                     f"{r['q_Tauber_MWm2']:.3f} | {r['SG_over_FR']:.3f} | {r['Tauber_over_FR']:.3f} |")
    lines += ["", f"SG/FR range: {min(ratios):.3f}–{max(ratios):.3f}.",
              "Cases with T_w > 0.95 T_e are skipped (heating tends to zero; ratios meaningless)."]
    (HERE / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
