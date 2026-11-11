"""Phase 4, Q5 extension: coarse 2D survival map for the nose tip, to express the gain from lateral conduction as a
shift of the 30-minute boundary (the 1D maps are runs/envelope/).

Axisymmetric sphere-cone (half-angle 5 deg, flank 30 mm, adiabatic back plane; D-2D-4), R_n = 1 and 3 mm, fully
catalytic, 20 vol% SiC, Mach 9-12 (step 1) at 20, 30 and 40 km, up to 1800 s, stopped at the first failure.
Output: runs/q5_2d/q5_envelope.csv; with the argument "refine", a second stage at half-Mach resolution around the
boundary (q5_envelope_refine.csv).
"""
import csv
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common  # noqa: E402
from uhtc import coupled as cp, materials as mt  # noqa: E402
from uhtc.conduction2d import NoseGeometry  # noqa: E402
from uhtc.oxidation import LayerOxidation, OxParams  # noqa: E402

HERE = Path(__file__).resolve().parent
MESH = dict(n_arc=30, n_flank=40, n_eta=30, eta1=0.003)


def run(args):
    M, z, R = args
    t0 = time.time()
    ox = LayerOxidation(OxParams(f_SiC=0.2))
    cond = cp.ConstantConditions(M, z, R, state=common.stag_state(M, z), catalysis="fc")
    g = NoseGeometry(R_n=R, theta_c=np.radians(5.0), L_flank=0.03, axisymmetric=True)
    r = cp.run_2d(cond, g, mt.composite(0.2), ox, common.T_CAP, record_every=5.0, dt_max=2.0, mesh_kw=MESH,
                  stop_on=common.stop_first)
    tf, mode = cp.first_failure(r.events)
    row = dict(mach=M, alt_km=z / 1e3, R_n_mm=R * 1e3, t_fail=tf, mode=mode, T_max=float(np.max(r.T_s)),
               t_end=float(r.t[-1]), t_glass=r.events["t_glass_loss"], t_active=r.events["t_active"],
               t_limit=r.events["t_limit"], E_res=r.energy_residual, cpu_s=round(time.time() - t0, 1))
    print({k: (round(v, 1) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    return row


def refine_jobs():
    """Second stage ('refine'): from q5_envelope.csv, bracket the 30-minute boundary at half-Mach resolution for each
    radius and altitude: the midpoint between the highest surviving and lowest failing Mach number, or, where even
    Mach 9 fails, Mach 8.0 and 8.5."""
    with open(HERE / "q5_envelope.csv", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    jobs = []
    for R in (1.0, 3.0):
        for z in (20.0, 30.0, 40.0):
            cell = [r for r in rows if float(r["R_n_mm"]) == R and float(r["alt_km"]) == z]
            ok = [float(r["mach"]) for r in cell if r["mode"] == "none"]
            bad = [float(r["mach"]) for r in cell if r["mode"] != "none"]
            if not ok:
                jobs += [(M, z * 1e3, R * 1e-3) for M in (8.0, 8.5)]
            elif bad and min(bad) - max(ok) > 0.5:
                jobs.append((max(ok) + 0.5, z * 1e3, R * 1e-3))
    return jobs


if __name__ == "__main__":
    refine = "refine" in sys.argv[1:]
    if refine:
        jobs = refine_jobs()
    else:
        jobs = [(M, z, R) for R in (0.001, 0.003) for z in (20e3, 30e3, 40e3) for M in (9.0, 10.0, 11.0, 12.0)]
    with Pool(int(os.environ.get("UHTC_WORKERS", 4))) as pool:
        rows = pool.map(run, jobs, chunksize=1)
    with open(HERE / ("q5_envelope_refine.csv" if refine else "q5_envelope.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
