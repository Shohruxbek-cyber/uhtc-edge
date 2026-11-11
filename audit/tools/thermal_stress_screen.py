"""Task 2.6: thermal-stress screening index during heat-up.  An indicative screen, not a failure criterion.

Sigma_th = E alpha dT / (1 - nu) is the stress of a fully restrained body; a real tip is partly free, so the index is an
upper bound on the elastic stress for a given dT.  The 2D model is started cold (300 K) at constant flight conditions, a
step entry, which is the severest heat-up the model can give.  The temperature field is taken at 0.1-30 s; dT is the
spread of nodal temperatures within a distance L of the stagnation point, for L = 2 R_n and 5 R_n.

Properties of ZrB2-20 vol% SiC, NASA TPSX record 28 (literature/text/TPSX_ZrB2-20SiC_id28.txt), interpolated linearly
in T and held constant beyond the end points:
    E      530 GPa 296 K, 523 1073 K, 446 1673 K, 106 1873 K
    alpha  5.00e-6 /K 673 K, 5.94e-6 1073 K, 6.33e-6 1473 K, 6.88e-6 1873 K
    nu     0.14 (assumed in the TPSX record)
    sigma_f (in-plane bend) 430 MPa 296 K, 548 673 K, 463 1073 K, 342 1473 K, 239 1873 K; 252 at 2073 K (through thickness)
Open cross-check (Neuman et al. 2013, ZrB2-30SiC in air, abstract): E 510 GPa at room temperature, 210 at 1773 K, 110 at
1873 K; sigma_f 680 MPa at room temperature, 750 at 1073 K, about 360 at 1773 and 1873 K.

Indices are evaluated with the properties at the coldest temperature in the region (the tensile side and the larger E)
and at the region's mean temperature.

Usage: python audit/tools/thermal_stress_screen.py   -> audit/thermal_stress_screen.csv
       python audit/tools/thermal_stress_screen.py --extend   -> adds 60-300 s for cases still rising at 30 s
"""
import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "runs"))
import common  # noqa: E402
from uhtc import coupled as cp, materials as mt  # noqa: E402
from uhtc.conduction2d import NoseGeometry  # noqa: E402
from uhtc.oxidation import LayerOxidation, OxParams  # noqa: E402

E_T, E_V = [296.15, 1073.15, 1673.15, 1873.15], [530e9, 523e9, 446e9, 106e9]
A_T, A_V = [673.15, 1073.15, 1473.15, 1873.15], [5.00e-6, 5.94e-6, 6.33e-6, 6.88e-6]
S_T, S_V = [296.15, 673.15, 1073.15, 1473.15, 1873.15, 2073.15], [430e6, 548e6, 463e6, 342e6, 239e6, 252e6]
NU = 0.14
EN_T, EN_V = [296.15, 1773.15, 1873.15], [510e9, 210e9, 110e9]  # Neuman et al. 2013, ZrB2-30SiC in air (abstract values)
SN_T, SN_V = [296.15, 1073.15, 1773.15, 1873.15], [680e6, 750e6, 360e6, 360e6]
CASES = [(8.5, 20e3, 1e-3), (9.25, 30e3, 1e-3), (10.25, 40e3, 1e-3), (9.75, 30e3, 3e-3), (10.25, 30e3, 10e-3)]  # 2D boundary, fc
TIMES = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0]
TIMES_EXT = [60.0, 100.0, 200.0, 300.0]  # --extend: only for cases whose index is largest at the last time of TIMES


def index(dT, T):
    E, a, s = np.interp(T, E_T, E_V), np.interp(T, A_T, A_V), np.interp(T, S_T, S_V)
    sig = E * a * dT / (1 - NU)
    return sig, s


def add_neuman(rows):
    """Same index with the E and strength of Neuman et al. 2013 (alpha from TPSX), at the coldest point of the region."""
    for r in rows:
        for name in ("2R", "5R"):
            T, dT = float(r[f"Tcold_{name}"]), float(r[f"dT_{name}"])
            sig = np.interp(T, EN_T, EN_V) * np.interp(T, A_T, A_V) * dT / (1 - NU)
            r[f"ratio_cold_neuman_{name}"] = round(float(sig / np.interp(T, SN_T, SN_V)), 3)


def write(rows):
    with open(ROOT / "audit/thermal_stress_screen.csv", "w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    if "--post" in sys.argv:  # recompute the Neuman columns from the saved temperatures
        rows = list(csv.DictReader(open(ROOT / "audit/thermal_stress_screen.csv", encoding="utf-8")))
        add_neuman(rows)
        write(rows)
        return
    if "--extend" in sys.argv:  # final review R4p2-11: cases whose index was still rising at the last time, sampled later
        rows = list(csv.DictReader(open(ROOT / "audit/thermal_stress_screen.csv", encoding="utf-8")))
        for M, z, R in CASES:
            mine = [r for r in rows if float(r["mach"]) == M and float(r["alt_km"]) == z / 1e3 and float(r["R_mm"]) == R * 1e3]
            last = max(mine, key=lambda r: float(r["t_s"]))
            if float(last["ratio_cold_2R"]) < max(float(r["ratio_cold_2R"]) for r in mine):
                continue  # peak already inside the sampled times
            for t_end in TIMES_EXT:
                if any(float(r["t_s"]) == t_end for r in mine):
                    continue
                row = run_row(M, z, R, t_end)
                rows.append({k: row.get(k, "") for k in rows[0]} | row)
                print(row, flush=True)
        add_neuman(rows)
        write(rows)
        return
    out = []
    for M, z, R in CASES:
        for t_end in TIMES:
            row = run_row(M, z, R, t_end)
            out.append(row)
            print(row, flush=True)
    add_neuman(out)
    write(out)


def run_row(M, z, R, t_end):
    state = common.stag_state(M, z)
    p = OxParams(f_SiC=0.2)
    ox = LayerOxidation(p)
    mat = mt.composite(0.2)
    cond = cp.ConstantConditions(M, z, R, state=state, catalysis="fc")
    g = NoseGeometry(R_n=R, theta_c=np.radians(5.0), L_flank=0.03, axisymmetric=True)
    r = cp.run_2d(cond, g, mat, ox, t_end, dt_max=0.5, mesh_kw=common.MESH_2D, heat_dist="lees")
    sol = r.extra["solver"]
    T = np.asarray(sol.T)
    d = np.hypot(*(sol.nodes - sol.nodes[sol.surf[0]]).T)
    row = dict(mach=M, alt_km=z / 1e3, R_mm=R * 1e3, t_s=t_end, T_stag=round(float(r.T_s[-1]), 1))
    for name, L in (("2R", 2 * R), ("5R", 5 * R)):
        m = d <= L
        Tm = float(np.sum(T[m] * sol.vol[m]) / np.sum(sol.vol[m]))
        Tc, dT = float(T[m].min()), float(T[m].max() - T[m].min())
        sc, fc = index(dT, Tc)
        sm, fm = index(dT, Tm)
        row.update({f"dT_{name}": round(dT, 1), f"Tcold_{name}": round(Tc, 1), f"Tmean_{name}": round(Tm, 1),
                    f"sig_cold_MPa_{name}": round(float(sc) / 1e6, 1), f"ratio_cold_{name}": round(float(sc / fc), 3),
                    f"sig_mean_MPa_{name}": round(float(sm) / 1e6, 1), f"ratio_mean_{name}": round(float(sm / fm), 3)})
    return row


if __name__ == "__main__":
    main()
