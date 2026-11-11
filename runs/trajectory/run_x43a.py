"""Phase 4, Q1-Q2 along a published flight path: the X-43A Mach 10 flight (runs/trajectory/x43a_mach10.csv).

1D runs: R_n = 1, 3, 10 mm x catalysis fc / nc / switch x SiC 10 / 20 / 30 vol% (27 runs).  The whole Mach >= 3
part (33.5-637.5 s) is simulated without stopping at the first failure (histories after the first failure are
indicative only, DECISIONS D-OX-4).  Histories are stored at 1 s.
Output: runs/trajectory/x43a_histories.csv (long format), x43a_summary.csv.
"""
import csv
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common  # noqa: E402
from uhtc import coupled as cp, materials as mt  # noqa: E402
from uhtc.oxidation import LayerOxidation, OxParams  # noqa: E402

HERE = Path(__file__).resolve().parent


def load_traj():
    a = np.loadtxt(HERE / "x43a_mach10.csv", delimiter=",")
    return a[:, 0], a[:, 1], a[:, 2]


def run(args):
    R, cat, f = args
    t, z, M = load_traj()
    cond = cp.TrajectoryConditions(t - t[0], z, M, R, dt_eval=2.0, state_fn=common.stag_state, catalysis=cat)
    ox = LayerOxidation(OxParams(f_SiC=f))
    r = cp.run_1d(cond, mt.composite(f), ox, t[-1] - t[0], depth=0.02, record_every=1.0, dt_max=1.0)
    rows = []
    for ti, Ts, s in zip(r.t, r.T_s, r.states):
        L = ox.layers(s)
        fl = cond.at(ti)
        rows.append(dict(R_n_mm=R * 1e3, catalysis=cat, f_SiC=f, t_flight=ti + t[0], mach=fl.mach,
                         alt_km=fl.altitude / 1e3, T_s=Ts, glass_um=L["glass"] * 1e6, ZrO2_um=L["ZrO2"] * 1e6,
                         dep_um=L["depleted"] * 1e6, x_B2O3=L["x_B2O3"],
                         T_tr=ox.T_transition(fl.p_O2), p_O2_kPa=fl.p_O2 / 1e3))
    tf, mode = cp.first_failure(r.events)
    summ = dict(R_n_mm=R * 1e3, catalysis=cat, f_SiC=f, T_s_max=float(np.max(r.T_s)),
                t_at_Tmax=float(r.t[np.argmax(r.T_s)] + t[0]),
                t_glass=r.events["t_glass_loss"] + t[0], t_active=r.events["t_active"] + t[0],
                t_limit=r.events["t_limit"] + t[0], t_fail=tf + t[0], mode=mode,
                ZrO2_um_end=rows[-1]["ZrO2_um"], glass_um_end=rows[-1]["glass_um"], dep_um_end=rows[-1]["dep_um"],
                energy_res=r.energy_residual)
    return rows, summ


if __name__ == "__main__":
    combos = [(R, cat, f) for R in (0.001, 0.003, 0.010) for cat in ("fc", "nc", "switch") for f in (0.1, 0.2, 0.3)]
    t, z, M = load_traj()
    # warm the stagnation-state cache once (serially through the pool) for the 2 s evaluation grid
    te = np.arange(0.0, t[-1] - t[0] + 1e-9, 2.0)
    pts = sorted({(round(float(np.interp(x, t - t[0], M)), 4), round(float(np.interp(x, t - t[0], z)), 1)) for x in te})
    with Pool(4) as pool:
        pool.starmap(common.stag_state, pts, chunksize=8)
        results = pool.map(run, combos, chunksize=1)
    hist = [r for rows, _ in results for r in rows]
    summ = [s for _, s in results]
    for name, data in (("x43a_histories.csv", hist), ("x43a_summary.csv", summ)):
        with open(HERE / name, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0]))
            w.writeheader()
            w.writerows(data)
    for s in summ:
        print({k: (round(v, 1) if isinstance(v, float) else v) for k, v in s.items()})
