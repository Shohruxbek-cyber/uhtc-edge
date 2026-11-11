"""Audit item: ratio of the sharp-cone running-length correlation to the baseline distribution on a 5 deg flank.

The paper quoted 1.1-3.8 from a comparison whose script was not kept.  This recomputes it from the code at Mach 8-14
(30 km) for R_n = 1 and 3 mm, both over s/R_n = 1.5-6 and over the whole 30 mm flank.  Writes audit/flank_ratio.csv.
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
from uhtc import atmosphere as atm, coupled as cp, heating as ht, materials as mt  # noqa: E402
from uhtc.conduction2d import NoseGeometry, Solver2D  # noqa: E402

out = []
for R in (1e-3, 3e-3):
    g = NoseGeometry(R_n=R, theta_c=math.radians(5.0), L_flank=0.03, axisymmetric=True)
    sol = Solver2D(g, mt.composite(0.2), 300.0, mesh_kw=common.MESH_2D)
    for M in (8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0):
        for z in (20e3, 30e3, 40e3):
            st = common.stag_state(M, z)
            flow = cp.ConstantConditions(M, z, R, state=st, catalysis="fc").at(0.0)
            T, _, rho = atm.state(z)
            V = M * math.sqrt(atm.GAMMA * atm.R_AIR * T)
            q = ht.tauber_sphere(rho, V, R, 0.0, 1.0)
            gt = sol.heating_shape(V, rho, q, dist="tauber", flow=flow)
            gl = sol.heating_shape(V, rho, q, dist="lees", flow=flow)
            m = sol.on_flank
            r = gt[m] / gl[m]
            x = sol.s[m] / R
            k = (x >= 1.5) & (x <= 6.0)
            out.append(dict(R_mm=R * 1e3, mach=M, alt_km=z / 1e3, lo_1p5_6=round(float(r[k].min()), 3), hi_1p5_6=round(float(r[k].max()), 3),
                            lo_flank=round(float(r.min()), 3), hi_flank=round(float(r.max()), 3)))
with open(ROOT / "audit/flank_ratio.csv", "w", newline="", encoding="utf-8") as h:
    w = csv.DictWriter(h, fieldnames=list(out[0]))
    w.writeheader()
    w.writerows(out)
print("s/R 1.5-6:", min(o["lo_1p5_6"] for o in out), max(o["hi_1p5_6"] for o in out))
print("whole flank:", min(o["lo_flank"] for o in out), max(o["hi_flank"] for o in out))
for o in out[:6]:
    print(o)
_ = np
