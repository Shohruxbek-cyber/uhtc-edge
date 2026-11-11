"""Phase 4, Q5 sensitivity: how much do the 2D results depend on the assumed flank?

D-2D-1 plans a x0.5 / x2 sensitivity on the flank heating; D-2D-4 fixes the flank length at 30 mm with an adiabatic
back plane.  Both matter if the whole slender tip behaves as one radiating body.  Cases: axisymmetric nose tip,
R_n = 1 mm, half-angle 5 deg, fully catalytic, 20 vol% SiC, 600 s, at Mach 8.5 and 9 (30 km); variants: baseline,
flank heating x0.5 and x2, flank length 15 mm and 60 mm (flank cells scaled with the length).
Output: runs/q5_2d/q5_sensitivity.csv
"""
import csv
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common  # noqa: E402
from uhtc import coupled as cp, materials as mt  # noqa: E402
from uhtc.conduction2d import NoseGeometry  # noqa: E402
from uhtc.oxidation import LayerOxidation, OxParams  # noqa: E402

HERE = Path(__file__).resolve().parent
T_END = 600.0
VARIANTS = {"baseline": (0.03, 1.0), "flank heating x0.5": (0.03, 0.5), "flank heating x2": (0.03, 2.0),
            "flank 15 mm": (0.015, 1.0), "flank 60 mm": (0.06, 1.0)}


def run(args):
    M, z, name = args
    L, fm = VARIANTS[name]
    mesh = dict(n_arc=30, n_flank=max(10, int(round(40 * L / 0.03))), n_eta=30, eta1=0.003)
    ox = LayerOxidation(OxParams(f_SiC=0.2))
    cond = cp.ConstantConditions(M, z, 0.001, state=common.stag_state(M, z), catalysis="fc")
    g = NoseGeometry(R_n=0.001, theta_c=np.radians(5.0), L_flank=L, axisymmetric=True)
    r = cp.run_2d(cond, g, mt.composite(0.2), ox, T_END, record_every=1.0, dt_max=1.0, mesh_kw=mesh, flank_mult=fm,
                  heat_dist="lees")  # baseline distribution (D-2D-7, D-2D-11); the run_2d default is the draft law
    tf, mode = cp.first_failure(r.events)
    row = dict(mach=M, alt_km=z / 1e3, variant=name, L_flank_mm=L * 1e3, flank_mult=fm, T_max=float(np.max(r.T_s)),
               T_end=float(r.T_s[-1]), t_fail=tf, mode=mode, t_glass=r.events["t_glass_loss"],
               t_limit=r.events["t_limit"], E_res=r.energy_residual)
    print({k: (round(v, 1) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    return row


if __name__ == "__main__":
    jobs = [(M, 30e3, v) for M in (8.5, 9.0) for v in VARIANTS]
    with Pool(int(os.environ.get("UHTC_WORKERS", 4))) as pool:
        rows = pool.map(run, jobs, chunksize=1)
    with open(HERE / "q5_sensitivity.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
