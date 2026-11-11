"""Task 2.4: equilibrium composition at the boundary-layer edge along the fully catalytic 2D boundary.

Prints, for each of the 15 band columns, the edge temperature, h_D/H_s and mole fractions, and the effective oxidant
pressure of the sensitivity case relative to the baseline 0.2095 p_s:
    p_eff = max(0.2095 p_s, (x_O2 + x_O) p_s)          every edge O atom counted as one O2 molecule
Writes audit/edge_oxygen.csv.
"""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "runs"))
import common  # noqa: E402

rows = list(csv.DictReader(open(ROOT / "runs/analysis/boundaries_2d.csv", encoding="utf-8")))
out = []
for r in rows:
    if r["catalysis"] != "fc" or float(r["alt_km"]) not in (20.0, 25.0, 30.0, 35.0, 40.0):
        continue
    for tag, M in (("survives", float(r["M_s_2d"])), ("fails", float(r["M_f_2d"]))):
        s = common.stag_state(M, float(r["alt_km"]) * 1e3)
        x = s.x_e
        out.append(dict(R_n_mm=r["R_n_mm"], alt_km=r["alt_km"], point=tag, mach=M, T_e=round(s.T_e), p_s_kPa=round(s.p_s / 1e3, 1),
                        hD_over_Hs=round(s.h_D / s.H_s, 3), x_O=round(x["O"], 4), x_O2=round(x["O2"], 4), x_NO=round(x["NO"], 4),
                        x_N=round(x["N"], 5), O_atom_share=round(x["O"] / (x["O"] + 2 * x["O2"] + x["NO"]), 3),
                        p_eff_over_base=round(max(1.0, (x["O"] + x["O2"]) / 0.2095), 3)))
with open(ROOT / "audit/edge_oxygen.csv", "w", newline="", encoding="utf-8") as h:
    w = csv.DictWriter(h, fieldnames=list(out[0]))
    w.writeheader()
    w.writerows(out)
for o in out:
    print(o)
print("x_O", min(o["x_O"] for o in out), max(o["x_O"] for o in out), "share", min(o["O_atom_share"] for o in out),
      max(o["O_atom_share"] for o in out), "p_eff", min(o["p_eff_over_base"] for o in out), max(o["p_eff_over_base"] for o in out),
      "hD", min(o["hD_over_Hs"] for o in out), max(o["hD_over_Hs"] for o in out))
