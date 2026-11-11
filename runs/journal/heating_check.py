"""Journal revision, ASK 4 (author, 11 Nov 2026): Sutton-Graves against Fay-Riddell at the 30-minute boundary of the 1 mm tip.

At the five altitudes of the uncertainty band (20, 25, 30, 35, 40 km), fully catalytic 2D maps, for the highest passing and the
lowest failing Mach number (runs/analysis/boundaries_2d.csv), with a cold wall (300 K) and the peak surface temperature at the
highest passing condition. Methods as in check V4 (validation/V4_heating/run_v4.py): wall gas is equilibrium air at (T_w, p_s);
Fay-Riddell (Hamilton 1982 eq. 6, Pr 0.71, Le 1.4) with a Newtonian velocity gradient, edge (a) equilibrium air at (p_s, H_s), the
model's stagnation state, and (b) frozen, undissociated air at H_s (h_D = 0), as in runs/v3/frozen_shock.py. Sutton-Graves is
Eq. (1) of the paper with the same wall enthalpy. No map is rerun.
Output: runs/journal/heating_check.csv and a printed summary.
"""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "runs"))
import common  # noqa: E402
from uhtc import atmosphere as atm, heating as ht, thermo as th  # noqa: E402

M_AIR = 0.028965
R_N = 1e-3


def main():
    b = [r for r in csv.DictReader(open(ROOT / "runs/analysis/boundaries_2d.csv", encoding="utf-8"))
         if r["R_n_mm"] == "1.0" and r["catalysis"] == "fc" and float(r["alt_km"]) in (20.0, 25.0, 30.0, 35.0, 40.0)]
    rows = []
    for x in b:
        z, T_hot = float(x["alt_km"]) * 1e3, float(x["T_max_s_2d"])
        for point, M in (("highest passing", float(x["M_s_2d"])), ("lowest failing", float(x["M_f_2d"]))):
            s = common.stag_state(M, z)
            for Tw in (300.0, T_hot):
                xw = th.equilibrium_air(Tw, s.p_s)
                Mw, hw, _ = th.mixture_props(xw, Tw)
                rho_w, mu_w = s.p_s * Mw / (th.R * Tw), float(atm.viscosity(Tw))
                q_sg = float(ht.sutton_graves(s.p_s, R_N, s.H_s, hw))
                due_a = ht.newtonian_velocity_gradient(R_N, s.p_s, s.p_inf, s.rho_e)
                q_a = float(ht.fay_riddell(rho_w, mu_w, s.rho_e, s.mu_e, s.H_s, hw, s.h_D, due_a))
                Tf = float(th.T_air_frozen(s.H_s))
                rho_f, mu_f = s.p_s * M_AIR / (th.R * Tf), float(atm.viscosity(Tf))
                due_b = ht.newtonian_velocity_gradient(R_N, s.p_s, s.p_inf, rho_f)
                q_b = float(ht.fay_riddell(rho_w, mu_w, rho_f, mu_f, s.H_s, hw, 0.0, due_b))
                rows.append(dict(alt_km=z / 1e3, point=point, mach=M, T_wall=round(Tw, 1), H_s_MJkg=round(s.H_s / 1e6, 3),
                                 h_D_frac=round(s.h_D / s.H_s, 4), q_SG_MW=q_sg / 1e6, q_FR_eq_MW=q_a / 1e6,
                                 q_FR_frozen_MW=q_b / 1e6, SG_over_FR_eq=q_sg / q_a, SG_over_FR_frozen=q_sg / q_b,
                                 FR_frozen_over_eq=q_b / q_a))
    out = ROOT / "runs/journal/heating_check.csv"
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
    for key in ("SG_over_FR_eq", "SG_over_FR_frozen", "FR_frozen_over_eq"):
        for wall in ("cold", "hot"):
            v = [r[key] for r in rows if (r["T_wall"] == 300.0) == (wall == "cold")]
            print(f"{key:20} {wall} wall: {min(v):.4f} to {max(v):.4f}")
    print("hot-wall temperatures:", sorted({r["T_wall"] for r in rows if r["T_wall"] != 300.0}))


if __name__ == "__main__":
    main()
