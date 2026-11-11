"""Q5: how much does lateral conduction change the stagnation temperature?  1D stagnation line vs quasi-2D nose region.

Revision scope (D-2D-7, D-2D-8): axisymmetric sphere-cone nose tip only, half-angle 5 deg (sensitivity 10 deg for
R_n = 1 mm), flank 30 mm, adiabatic back plane, baseline heating distribution "lees" (and the draft "tauber" for
comparison at the constant conditions); R_n = 1 and 3 mm at five constant conditions (600 s, fully catalytic), and the
X-43A trajectory for R_n = 1, 3 and 10 mm and both catalysis bounds (Part 1.5).  The 1D reference uses the same
stagnation heating and the production 1D settings (20 mm slab, adiabatic back); 20 vol% SiC.
Each 2D run records the stagnation-node temperature and the three failure times there.
Output: runs/q5_2d/q5_results.csv, q5_histories.csv.
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
CONDS = [(7.0, 25e3), (8.0, 30e3), (8.5, 30e3), (9.0, 30e3), (10.0, 35e3)]
MESH = common.MESH_2D


def conditions(M, z, R, cat, traj=None):
    if traj is None:
        return cp.ConstantConditions(M, z, R, state=common.stag_state(M, z), catalysis=cat)
    t, zz, MM = traj
    return cp.TrajectoryConditions(t - t[0], zz, MM, R, dt_eval=2.0, state_fn=common.stag_state, catalysis=cat)


def run(args):
    kind, M, z, R, theta, cat, dist = args
    traj = None
    t_end = T_END
    if kind == "x43a":
        a = np.loadtxt(common.ROOT / "runs" / "trajectory" / "x43a_mach10.csv", delimiter=",")
        traj = (a[:, 0], a[:, 1], a[:, 2])
        t_end = a[-1, 0] - a[0, 0]
    out = {}
    hist = []
    for model in ("1d", "2d"):
        ox = LayerOxidation(OxParams(f_SiC=0.2))
        cond = conditions(M, z, R, cat, traj)
        if model == "1d":
            r = cp.run_1d(cond, mt.composite(0.2), ox, t_end, depth=0.02, record_every=1.0, dt_max=1.0)
        else:
            g = NoseGeometry(R_n=R, theta_c=np.radians(theta), L_flank=0.03, axisymmetric=True)
            r = cp.run_2d(cond, g, mt.composite(0.2), ox, t_end, record_every=1.0, dt_max=1.0, mesh_kw=MESH,
                          heat_dist=dist)
        tf, mode = cp.first_failure(r.events)
        out[model] = dict(T_max=float(np.max(r.T_s)), T_end=float(r.T_s[-1]), t_fail=tf, mode=mode,
                          t_glass=r.events["t_glass_loss"], t_active=r.events["t_active"], t_limit=r.events["t_limit"],
                          E_res=r.energy_residual)
        for ti, Ts in zip(r.t, r.T_s):
            hist.append(dict(kind=kind, mach=M, alt_km=z / 1e3, R_n_mm=R * 1e3, theta_deg=theta, catalysis=cat,
                             heat_dist=dist, model=model, t=ti, T_s=Ts))
    row = dict(kind=kind, mach=M, alt_km=z / 1e3, R_n_mm=R * 1e3, theta_deg=theta, catalysis=cat, heat_dist=dist)
    for m in ("1d", "2d"):
        for k, v in out[m].items():
            row[f"{m}_{k}"] = v
    print({k: (round(v, 1) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    return row, hist


if __name__ == "__main__":
    jobs = [("const", M, z, R, 5.0, "fc", d) for d in ("lees", "tauber") for (M, z) in CONDS for R in (0.001, 0.003)]
    jobs += [("const", M, z, 0.001, 10.0, "fc", "lees") for (M, z) in CONDS[1:4]]
    jobs += [("x43a", 0.0, 0.0, R, 5.0, cat, "lees") for R in (0.001, 0.003, 0.010) for cat in ("fc", "nc")]
    with Pool(int(os.environ.get("UHTC_WORKERS", 4))) as pool:
        res = pool.map(run, jobs, chunksize=1)
    rows = [r for r, _ in res]
    hist = [h for _, hs in res for h in hs]
    for name, data in (("q5_results.csv", rows), ("q5_histories.csv", hist)):
        with open(HERE / name, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0]))
            w.writeheader()
            w.writerows(data)
