"""Temperature jump after glass loss: histories for the glass-loss-governed cases with the three catalysis options
(fully catalytic, non-catalytic, switched at glass loss, D-HEAT-4).  The switch cannot change the time to first
failure (it happens at failure (a)), only what follows; this run shows that part.
Cases: Mach 8 at 30 and 35 km with R_n = 1 mm (edge gas weakly dissociated, h_D/H_s = 0.06-0.08) and Mach 10 at
35 km and Mach 10.5 at 38 km with R_n = 10 mm (h_D/H_s = 0.23-0.30; glass-loss cells of the non-catalytic map);
20 vol% SiC, 1800 s, recorded every 2 s.
Output: runs/histories/jump_histories.csv
"""
import csv
import os
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common  # noqa: E402
from uhtc import coupled as cp, materials as mt  # noqa: E402
from uhtc.oxidation import LayerOxidation, OxParams  # noqa: E402

HERE = Path(__file__).resolve().parent


def run(args):
    M, z, R, cat = args
    ox = LayerOxidation(OxParams(f_SiC=0.2))
    cond = cp.ConstantConditions(M, z, R, state=common.stag_state(M, z), catalysis=cat)
    r = cp.run_1d(cond, mt.composite(0.2), ox, common.T_CAP, depth=0.02, record_every=2.0,
                  stop_on=lambda e, T, t: "T > abort" if T > common.T_ABORT else None)
    tf, mode = cp.first_failure(r.events)
    return [dict(mach=M, alt_km=z / 1e3, R_n_mm=R * 1e3, catalysis=cat, t=t, T_s=Ts,
                 glass_um=ox.layers(s)["glass"] * 1e6, t_fail=tf, mode=mode) for t, Ts, s in zip(r.t, r.T_s, r.states)]


if __name__ == "__main__":
    cells = [(8.0, 30e3, 0.001), (8.0, 35e3, 0.001), (10.0, 35e3, 0.010), (10.5, 38e3, 0.010)]
    jobs = [(M, z, R, c) for (M, z, R) in cells for c in ("fc", "nc", "switch")]
    with Pool(int(os.environ.get("UHTC_WORKERS", 4))) as pool:
        res = pool.map(run, jobs, chunksize=1)
    rows = [r for rr in res for r in rr]
    with open(HERE / "jump_histories.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for rr in res:
        print(rr[0]["mach"], rr[0]["alt_km"], rr[0]["R_n_mm"], rr[0]["catalysis"], "t_fail", round(rr[0]["t_fail"], 1), rr[0]["mode"],
              "T_s end", round(rr[-1]["T_s"]), "max", round(max(r["T_s"] for r in rr)))
