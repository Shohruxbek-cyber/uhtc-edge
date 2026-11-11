"""Task 2.5: order-of-magnitude check of glass flow away from the stagnation point.  NOT part of the model.

Compares, at the conditions of the fully catalytic 2D boundary, the thinning rate of the glass by flow with the
thinning rate by evaporation that the model computes (uhtc.oxidation.LayerOxidation.evap).

Film velocities (the forms set in the revision brief):
    pressure driven   v_p = (d^2 / mu) (dp/ds),  dp/ds ~ p_s / R_n
    shear driven      v_t = tau_w d / mu
Thinning of the film over the cap by continuity, film flux d v leaving a region of size R_n:
    w = d v / R_n
The exact lubrication means are v/3 and v/2 and the axisymmetric divergence adds a factor 2, so w is good to about a
factor of 2 either way; the estimate does not claim better.

Viscosity: measured fibre-draw constants of Mauro et al. 2022, Table II (1973-2573 K, so no extrapolation is needed).
    log10 mu [Pa s] = log10 mu_inf + E / (2.303 R T)
    LA, Type II natural silica, E = 587.4 kJ/mol, log10 mu_inf = -8.60
    F320, fluorine doped,        E = 495.8 kJ/mol, log10 mu_inf = -6.80 (the softest of the seven at 2000-2200 K)
Silica is the right end member: the model's own glass holds x_B2O3 < 0.003 at these temperatures (histories.csv).

Wall shear: Reynolds analogy for a laminar boundary layer, c_f/2 = St Pr^(2/3) (Crabtree et al. 1965, Section 1.3),
    tau_w(phi) = Phi_s g(phi) Pr^(2/3) u_e(phi),  g = cos(phi) on the cap (Section 2.2 of the paper),
with Phi_s the stagnation conductance of Eq. 1, p/p_s = cos^2(phi) and u_e the isentropic expansion from the stagnation
point (gamma = 1.4), as in the model.  The maximum over phi <= 70 deg is used.

Usage: python audit/tools/glass_flow_estimate.py   -> audit/glass_flow_table.csv and a printed summary
"""
import csv
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "runs"))
import common  # noqa: E402
from uhtc import coupled as cp, data as dt, heating as ht  # noqa: E402
from uhtc.oxidation import LayerOxidation, OxParams  # noqa: E402

R_GAS = 8.314462618
VISC = {"LA": (587.4e3, -8.60), "F320": (495.8e3, -6.80)}
PR = 0.71
DELTAS = (10e-6, 50e-6, 100e-6)  # the brief's 10-50 um, plus 100 um: the model's own films reach 90-140 um


def mu(T, sample):
    E, l0 = VISC[sample]
    return 10.0 ** (l0 + E / (math.log(10.0) * R_GAS * T))


def tau_w(state, Phi_s):
    """Largest wall shear on the cap, Pa (Reynolds analogy)."""
    phi = np.radians(np.linspace(1.0, 70.0, 139))
    pr = np.cos(phi) ** 2
    gam = 1.4
    ue = np.sqrt(2 * gam / (gam - 1) * state.p_s / state.rho_e * (1 - pr ** ((gam - 1) / gam)))
    return float(np.max(Phi_s * np.cos(phi) * PR ** (2 / 3) * ue))


def rows_near_boundary():
    b = list(csv.DictReader(open(ROOT / "runs/analysis/boundaries_2d.csv", encoding="utf-8")))
    runs = {}
    for f in ("grid", "refine"):
        for r in common.load(ROOT / f"runs/envelope_2d/{f}.jsonl"):
            runs[(r["R_n"], r["catalysis"], r["alt"], r["mach"])] = r
    out = []
    for x in b:
        R, z = float(x["R_n_mm"]) * 1e-3, float(x["alt_km"]) * 1e3
        if x["catalysis"] != "fc" or float(x["alt_km"]) not in (20.0, 30.0, 40.0):
            continue
        for tag, M, key in (("survives", float(x["M_s_2d"]), "T_s_end"), ("fails", float(x["M_f_2d"]), "T_s_max")):
            r = runs[(R, "fc", z, M)]
            out.append(dict(tag=tag, M=M, z=z, R=R, T=r[key], film_um=r["glass_um_end"], mode=r["mode"]))
    return out


def main():
    ox = LayerOxidation(OxParams())
    table = []
    for c in rows_near_boundary():
        s = common.stag_state(c["M"], c["z"])
        Phi = ht.K_SG_AIR * math.sqrt(s.p_s / ht.ATM / c["R"])
        fl = cp.FlowState(s.H_s, s.p_s, s.p_O2_w, s.h_D, Phi, 1.0, Phi, s.mach, s.altitude)
        _, jS, _ = ox.evap(dict(n_B=1e-3, n_S=1.0), c["T"], fl)
        w_ev = jS / dt.RHO_SiO2_GLASS
        tw = tau_w(s, Phi)
        row = dict(point=c["tag"], mach=c["M"], alt_km=c["z"] / 1e3, R_mm=c["R"] * 1e3, T_s=round(c["T"], 1),
                   p_s_kPa=round(s.p_s / 1e3, 1), tau_w_Pa=round(tw), film_um=round(c["film_um"], 1),
                   w_evap_um_s=w_ev * 1e6)
        for sample in VISC:
            m = mu(c["T"], sample)
            row[f"mu_{sample}_Pa_s"] = m
            for d in DELTAS:
                wp = d ** 3 * (s.p_s / c["R"]) / (m * c["R"])
                wt = d ** 2 * tw / (m * c["R"])
                row[f"w_{sample}_{d * 1e6:.0f}um_um_s"] = (wp + wt) * 1e6
                row[f"ratio_{sample}_{d * 1e6:.0f}um"] = (wp + wt) / w_ev if w_ev > 0 else math.inf
            # film thickness at which flow equals evaporation
            f = lambda d: d ** 3 * (s.p_s / c["R"]) / (m * c["R"]) + d ** 2 * tw / (m * c["R"]) - w_ev  # noqa: E731
            lo, hi = 1e-8, 1e-2
            for _ in range(200):
                mid = math.sqrt(lo * hi)
                lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
            row[f"d_equal_{sample}_um"] = mid * 1e6
        table.append(row)
    with open(ROOT / "audit/glass_flow_table.csv", "w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=list(table[0]))
        w.writeheader()
        for r in table:
            w.writerow({k: (f"{v:.3g}" if isinstance(v, float) and k not in ("mach", "alt_km", "R_mm", "T_s", "p_s_kPa", "film_um") else v) for k, v in r.items()})
    print(f"{'pt':8s}{'M':>6s}{'z':>4s}{'R':>4s}{'T':>7s}{'p_s':>7s}{'tau':>7s}{'film':>6s}{'w_ev':>9s}{'muLA':>9s}{'r10':>9s}{'r50':>9s}{'r100':>9s}{'dEq':>7s}{'muF':>9s}{'r50F':>9s}{'r100F':>9s}{'dEqF':>7s}")
    for r in table:
        print(f"{r['point']:8s}{r['mach']:6.2f}{r['alt_km']:4.0f}{r['R_mm']:4.0f}{r['T_s']:7.0f}{r['p_s_kPa']:7.1f}{r['tau_w_Pa']:7d}{r['film_um']:6.0f}{r['w_evap_um_s']:9.3f}"
              f"{r['mu_LA_Pa_s']:9.2e}{r['ratio_LA_10um']:9.1e}{r['ratio_LA_50um']:9.1e}{r['ratio_LA_100um']:9.1e}{r['d_equal_LA_um']:7.0f}"
              f"{r['mu_F320_Pa_s']:9.2e}{r['ratio_F320_50um']:9.1e}{r['ratio_F320_100um']:9.1e}{r['d_equal_F320_um']:7.0f}")
    for T in (2000, 2100, 2200, 2300, 2400, 2450):
        print(T, {k: f"{mu(T, k):.2e}" for k in VISC})


if __name__ == "__main__":
    main()
