"""Reaction-front temperature at the 30-minute boundary of the 2D maps (diagnostic, D-OX-13).

The growth law k_p is evaluated at the reaction-front temperature T_f (the solid node under the oxide scale), not at
the surface temperature T_s (uhtc.oxidation, LayerOxidation.advance); the two differ by the drop across the scale.
The grid records only T_s, so this script reruns, for every bracketed (R_n, bound, altitude) column of grid + refine,
the highest surviving and the lowest failing case with the baseline inputs unchanged, and records the largest T_f at
the stagnation node and the time it spends above 2073 K (the end of the high-temperature fit of k_p).
The rerun must reproduce the stored T_s_max and failure time (checked by runs/analysis/analyze_2d.py).
Output: runs/envelope_2d/front.jsonl (+ .csv).
"""
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common  # noqa: E402
import run_envelope_2d as R2  # noqa: E402
from uhtc import coupled as cp, data as dt, materials as mt  # noqa: E402
from uhtc.conduction2d import NoseGeometry  # noqa: E402
from uhtc.oxidation import LayerOxidation, OxParams  # noqa: E402

HERE = Path(__file__).resolve().parent


class FrontOx(LayerOxidation):
    """LayerOxidation that also keeps the largest front temperature and the time spent above the fit range."""

    def advance(self, s, dt_, T_s, T_f, flow):
        s2 = super().advance(s, dt_, T_s, T_f, flow)
        if s.get("_Tf_max", -1.0) < T_f:
            s2["_Tf_max"], s2["_Ts_at_Tf_max"], s2["_t_Tf_max"] = T_f, T_s, s.get("_t", 0.0) + dt_
        s2["_t"] = s.get("_t", 0.0) + dt_
        s2["_t_Tf_hi"] = s.get("_t_Tf_hi", 0.0) + (dt_ if T_f > dt.KP_FIT_TMAX else 0.0)
        s2["_t_Ts_hi"] = s.get("_t_Ts_hi", 0.0) + (dt_ if T_s > dt.KP_FIT_TMAX else 0.0)
        return s2


def run_front_case(case):
    """As common.run_case_2d (same inputs, mesh, step cap and stop rule), with FrontOx."""
    t0 = time.time()
    ox = FrontOx(OxParams(f_SiC=case.get("f_SiC", 0.2), **case.get("ox", {})))
    mat = mt.composite(case.get("f_SiC", 0.2))
    state = common.stag_state(case["mach"], case["alt"])
    cond = cp.ConstantConditions(case["mach"], case["alt"], case["R_n"], state=state, catalysis=case["catalysis"])
    g = NoseGeometry(R_n=case["R_n"], theta_c=math.radians(5.0), L_flank=0.03, axisymmetric=True)
    try:
        r = cp.run_2d(cond, g, mat, ox, common.T_CAP, record_every=5.0, dt_max=2.0, mesh_kw=common.MESH_2D,
                      stop_on=common.stop_first, heat_dist=case["heat_dist"])
    except Exception as e:  # recorded, never silently dropped
        return dict(case_id=case["id"], mach=case["mach"], alt=case["alt"], R_n=case["R_n"],
                    catalysis=case["catalysis"], error=repr(e)[:200], cpu_s=round(time.time() - t0, 2))
    tf, mode = cp.first_failure(r.events)
    s = r.states[-1]
    return dict(case_id=case["id"], mach=case["mach"], alt=case["alt"], R_n=case["R_n"], catalysis=case["catalysis"],
                t_fail=tf, mode=mode, T_s_max=float(np.max(r.T_s)), t_end=float(r.t[-1]),
                T_f_max=s.get("_Tf_max", math.nan), T_s_at_T_f_max=s.get("_Ts_at_Tf_max", math.nan),
                t_T_f_max=s.get("_t_Tf_max", math.nan), t_T_f_above_fit=s.get("_t_Tf_hi", 0.0),
                t_T_s_above_fit=s.get("_t_Ts_hi", 0.0), error="", cpu_s=round(time.time() - t0, 2))


def boundary_cases():
    rows = common.load(HERE / "grid.jsonl") + common.load(HERE / "refine.jsonl")
    cols = R2.columns(rows)
    cases = []
    for (R, cat, z), (lo, hi) in sorted(cols.items()):
        if lo is None or hi is None:
            continue
        cases.append(R2.case(lo, z, R, cat, "front_s"))
        cases.append(R2.case(hi, z, R, cat, "front_f"))
    return cases


if __name__ == "__main__":
    common.run_all(boundary_cases(), HERE / "front.jsonl", runner=run_front_case)
