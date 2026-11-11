"""Version 3, Task 1, D-FLOW-6 diagnostic: how much does keeping the local B2O3 fraction (not advecting it with the flow) matter?

One low-altitude boundary point (1 mm, 20 km, Mach 8.75: the lowest failing Mach of the fully catalytic 2D column), production
2D settings, Lees distribution, glass flow on and off, run to 5, 10, 20, 30 and 60 s from a cold start. At each time and surface
node: the B2O3 fraction x_B of the film and the relative change of the film thickness caused by the flow, |delta_on - delta_off| /
delta_off. Their product bounds the B2O3 fraction carried by the flow that the composition rule misplaces. Output:
runs/v3/b2o3_advection.csv and a printed table.
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
from uhtc import coupled as cp, materials as mt  # noqa: E402
from uhtc.conduction2d import NoseGeometry  # noqa: E402
from uhtc.glassflow import FlowParams  # noqa: E402
from uhtc.oxidation import LayerOxidation, OxParams  # noqa: E402

M, Z, R = 8.75, 20e3, 1e-3


def states_at(t_end, flow):
    ox = LayerOxidation(OxParams(f_SiC=0.2))
    cond = cp.ConstantConditions(M, Z, R, state=common.stag_state(M, Z), catalysis="fc")
    g = NoseGeometry(R_n=R, theta_c=math.radians(5.0), L_flank=0.03, axisymmetric=True)
    r = cp.run_2d(cond, g, mt.composite(0.2), ox, t_end, dt_max=2.0, mesh_kw=common.MESH_2D, heat_dist="lees",
                  glass_flow=FlowParams() if flow else None)
    st = r.extra["states"]
    d = np.array([ox.glass_thickness(s) for s in st])
    xb = np.array([ox.x_B(s) for s in st])
    return d, xb


def main():
    rows = []
    for t in (5.0, 10.0, 20.0, 30.0, 60.0):
        d_on, xb_on = states_at(t, True)
        d_off, _ = states_at(t, False)
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(d_off > 1e-9, np.abs(d_on - d_off) / d_off, 0.0)
        prod = xb_on * rel
        k = int(np.argmax(prod))
        rows.append(dict(t_s=t, xB_max=float(xb_on.max()), rel_flow_change_max=float(rel.max()), bound_max=float(prod.max()),
                         node=k, xB_at=float(xb_on[k]), rel_at=float(rel[k]), film_um_stag=float(d_on[0] * 1e6)))
        print(rows[-1], flush=True)
    with open(ROOT / "runs/v3/b2o3_advection.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
