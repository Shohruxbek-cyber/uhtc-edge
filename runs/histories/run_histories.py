"""Phase 4, Q1-Q2: full time histories (surface temperature, layer thicknesses) for selected constant conditions.

Conditions (Mach, altitude): (7, 25 km), (8, 30 km), (9, 30 km), (10, 35 km) x R_n = 1, 3, 10 mm x fc / nc,
20 vol% SiC, 1800 s or until T_s > 3400 K, recorded every 2 s; the run is not stopped at the first failure (later
history indicative only, D-OX-4).
Output: runs/histories/histories.csv (long format) and events.csv.
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
CONDS = [(7.0, 25e3), (8.0, 30e3), (9.0, 30e3), (10.0, 35e3)]


def run(args):
    M, z, R, cat = args
    ox = LayerOxidation(OxParams(f_SiC=0.2))
    cond = cp.ConstantConditions(M, z, R, state=common.stag_state(M, z), catalysis=cat)
    r = cp.run_1d(cond, mt.composite(0.2), ox, common.T_CAP, depth=0.02, record_every=2.0,
                  stop_on=lambda e, T, t: "T > abort" if T > common.T_ABORT else None)
    rows = []
    for t, Ts, Tsub, q, s in zip(r.t, r.T_s, r.T_sub, r.q_conv, r.states):
        L = ox.layers(s)
        rows.append(dict(mach=M, alt_km=z / 1e3, R_n_mm=R * 1e3, catalysis=cat, t=t, T_s=Ts, T_front=Tsub,
                         q_conv_MW=q / 1e6, glass_um=L["glass"] * 1e6, ZrO2_um=L["ZrO2"] * 1e6,
                         dep_um=L["depleted"] * 1e6, x_B2O3=L["x_B2O3"], dcrit_nm=s.get("_dcrit", np.nan) * 1e9))
    tf, mode = cp.first_failure(r.events)
    ev = dict(mach=M, alt_km=z / 1e3, R_n_mm=R * 1e3, catalysis=cat, t_glass=r.events["t_glass_loss"],
              t_active=r.events["t_active"], t_limit=r.events["t_limit"], t_fail=tf, mode=mode,
              T_transition=ox.T_transition(cond.flow.p_O2), stopped=r.stopped)
    return rows, ev


if __name__ == "__main__":
    jobs = [(M, z, R, cat) for (M, z) in CONDS for R in (0.001, 0.003, 0.010) for cat in ("fc", "nc")]
    import os
    with Pool(int(os.environ.get("UHTC_WORKERS", "4"))) as pool:
        res = pool.map(run, jobs, chunksize=1)
    for name, data in (("histories.csv", [r for rows, _ in res for r in rows]), ("events.csv", [e for _, e in res])):
        with open(HERE / name, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0]))
            w.writeheader()
            w.writerows(data)
    for _, e in res:
        print(e)
