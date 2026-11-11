"""V9b pre-run estimate (written before the V9b plan): stagnation-point velocity gradient beta of the plasmatron samples and the
size of the film-flow sink there, compared with the flight boundary.

beta from the measured cold-wall flux by inverting Fay-Riddell (uhtc.heating.fay_riddell, fully catalytic, Le = 1.4) at the
rebuilt edge state (H_e, p = P_stat + P_dyn; equilibrium edge temperature of the model) and a calorimeter wall at 350 K, the value
of the source's rebuild; cross-check beta ~ V_e / R_probe with the 15 mm radius of the 30 mm probe (Table 1 of the source, A1-11).

Stagnation-point flow sink (axisymmetric, delta smooth at r = 0): w_flow = tau' delta^2 / mu + 2 G' delta^3 / (3 mu), with
tau = tau' r, tau' = Phi Pr^(2/3) beta (Reynolds analogy, D-FLOW-3) and G = -dp/dr = G' r, G' = rho_e beta^2.
Usage: python validation/V9b_glassflow/beta_estimate.py   -> beta_estimate.json and a printed table
"""
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "validation" / "V9_plasmatron"))
from uhtc import atmosphere as atm, coupled as cp, heating as ht, thermo as th  # noqa: E402
from uhtc.glassflow import PR, viscosity  # noqa: E402
import run_v9  # noqa: E402

V_E = {"1.8": (100.0, 159.0), "1.5": (100.0, 180.0), "1.6": (100.0, 207.0)}  # Table 1, A1-11, low and high stage, m/s
R_PROBE = 0.015


def edge(H, p):
    with th.extended_tables(cp.T_EDGE_MAX):
        Te = th.T_air_equilibrium(H, p, T_hi=cp.T_EDGE_MAX - 1.0)
        x = th.equilibrium_air(Te, p)
        M, _, hD = th.mixture_props(x, Te)
    return Te, p * M / (th.R * Te), float(atm.viscosity(Te)), hD


def main():
    out = []
    for name, (ts, lo, hi, *_rest) in run_v9.CASES.items():
        for stage, (q, H, pd), Ve in (("low", lo, V_E[name][0]), ("high", hi, V_E[name][1])):
            p = run_v9.P_STAT + pd
            Te, rho_e, mu_e, hD = edge(H, p)
            Tw = 350.0
            rho_w = p * 0.028965 / (th.R * Tw)
            mu_w = float(atm.viscosity(Tw))
            h_w = float(th.h_air_frozen(Tw))
            unit = ht.fay_riddell(rho_w, mu_w, rho_e, mu_e, H, h_w, hD, 1.0)
            beta = (q / unit) ** 2
            Phi = q / (H - float(th.h_air_frozen(300.0)))  # as PrescribedConditions
            tau1 = Phi * PR ** (2 / 3) * beta  # tau / r, Pa/m
            G1 = rho_e * beta ** 2  # G / r, Pa/m^2
            row = dict(sample=name, stage=stage, q_cw=q, H=H, p=p, T_e=Te, rho_e=rho_e, beta_FR=beta, beta_Ve=Ve / R_PROBE,
                       tau_at_1mm=tau1 * 1e-3, tau1=tau1, G1=G1)
            for T in (2100.0, 2250.0):
                mu = float(viscosity(T, "LA"))
                for d in (10e-6, 50e-6):
                    row[f"w_flow_um_s_T{T:.0f}_d{d * 1e6:.0f}"] = (tau1 * d * d / mu + 2 * G1 * d ** 3 / (3 * mu)) * 1e6
            out.append(row)
    # flight comparison: 1 mm sphere, 20 km, Mach 8.75 (first failing point of the fc 2D boundary)
    from uhtc import stagnation as st
    s = st.stagnation_state(8.75, 20e3)
    beta_f = ht.newtonian_velocity_gradient(1e-3, s.p_s, s.p_inf, s.rho_e)
    Phi_f = ht.K_SG_AIR * math.sqrt(s.p_s / ht.ATM / 1e-3)
    flight = dict(beta=beta_f, tau1=Phi_f * PR ** (2 / 3) * beta_f, G1=s.rho_e * beta_f ** 2, tau_at_0p2mm=Phi_f * PR ** (2 / 3) * beta_f * 2e-4)
    mu = float(viscosity(2380.0, "LA"))
    flight["w_flow_um_s_T2380_d10"] = (flight["tau1"] * 1e-10 / mu + 2 * flight["G1"] * 1e-15 / (3 * mu)) * 1e6
    (HERE / "beta_estimate.json").write_text(json.dumps(dict(plasmatron=out, flight_1mm_20km_M8p75=flight), indent=1), encoding="utf-8")
    for r in out:
        print(f"{r['sample']} {r['stage']:4s} beta_FR {r['beta_FR']:9.0f} beta_Ve {r['beta_Ve']:7.0f} tau(1mm) {r['tau_at_1mm']:6.3f} Pa "
              f"w_flow(2250 K, 50 um) {r['w_flow_um_s_T2250_d50']:.2e} um/s")
    print("flight", {k: f"{v:.3g}" for k, v in flight.items()})


if __name__ == "__main__":
    main()
