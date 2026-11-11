"""Version 3, Task 4.2: frozen-shock-layer bound of the stagnation heating at the five 2D boundary points of the 1 mm tip.

Fully catalytic Fay-Riddell stagnation heating (uhtc.heating.fay_riddell, Le = 1.4) with the edge gas (a) in equilibrium at
(p_s, H_s), as the model's stagnation state, and (b) frozen at the free-stream composition (undissociated air at the stagnation
enthalpy, h_D = 0): edge temperature from the frozen-air enthalpy, density from the ideal gas of free-stream composition, Sutherland
viscosity; Newtonian velocity gradient with each edge density; cold wall (300 K) and the boundary wall temperature. The ratio (b)/(a)
bounds the effect of a shock layer too fast for chemistry to keep up. No maps are rerun.
Output: runs/v3/frozen_shock.csv and a printed table.
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


def main():
    b = [r for r in csv.DictReader(open(ROOT / "runs/analysis/boundaries_2d.csv", encoding="utf-8"))
         if r["R_n_mm"] == "1.0" and r["catalysis"] == "fc" and float(r["alt_km"]) in (20.0, 25.0, 30.0, 35.0, 40.0)]
    rows = []
    R = 1e-3
    for x in b:
        z, M = float(x["alt_km"]) * 1e3, float(x["M_f_2d"])
        s = common.stag_state(M, z)
        for Tw in (300.0, float(x["T_max_s_2d"])):
            p = s.p_s
            rho_w, mu_w = p * M_AIR / (th.R * Tw), float(atm.viscosity(Tw))
            h_w = float(th.h_air_frozen(Tw))
            # (a) equilibrium edge (the stagnation state of the model)
            beta_a = ht.newtonian_velocity_gradient(R, s.p_s, s.p_inf, s.rho_e)
            qa = ht.fay_riddell(rho_w, mu_w, s.rho_e, s.mu_e, s.H_s, h_w, s.h_D, beta_a)
            # (b) frozen edge: undissociated air at H_s
            Tf = float(th.T_air_frozen(s.H_s))
            rho_f, mu_f = p * M_AIR / (th.R * Tf), float(atm.viscosity(Tf))
            beta_b = ht.newtonian_velocity_gradient(R, s.p_s, s.p_inf, rho_f)
            qb = ht.fay_riddell(rho_w, mu_w, rho_f, mu_f, s.H_s, h_w, 0.0, beta_b)
            rows.append(dict(alt_km=z / 1e3, mach=M, T_wall=round(Tw, 1), T_e_eq=round(s.T_e, 1), T_e_frozen=round(Tf, 1),
                             h_D_frac=round(s.h_D / s.H_s, 4), q_eq_MW=qa / 1e6, q_frozen_MW=qb / 1e6, ratio=qb / qa))
    with open(ROOT / "runs/v3/frozen_shock.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})


if __name__ == "__main__":
    main()
