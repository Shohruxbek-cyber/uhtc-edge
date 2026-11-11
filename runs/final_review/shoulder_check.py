"""Shoulder check of the 2D boundary (external review, item A4 / run R1).

The 2D model gives every surface node its own oxide state, but uses the stagnation wall pressure and O2 pressure on the
whole surface (D-2D-3), and the failure clocks are read at the stagnation node only.  This script evaluates criteria
(a), (b) and (c) at every surface node, for
  * the 15 band columns of the 2D maps (fully catalytic, R_n = 1, 3, 10 mm, 20-40 km in 5 km steps), at the highest
    passing (M_s) and the lowest failing (M_f) Mach number of the baseline search (grid + refine, 0.25 Mach), and
  * the six 2D X-43A runs (R_n = 1, 3, 10 mm, fully catalytic and non-catalytic),
each in two modes:
  'stag'  the production model (stagnation pressure at every node), run to the end without stopping at the first
          stagnation failure;
  'local' the oxide model of each node sees the local wall pressure and O2 pressure, p/p_s from the modified-Newtonian
          law of the Lees distribution (uhtc.coupled.run_2d, local_pressure=True); heating unchanged.
Each node's failure times come from its own oxide state, with its own surface temperature and glass.  Where a node
other than the stagnation node fails within 1800 s at M_s, the Mach number is lowered in 0.25 steps until no node
fails, which gives the boundary with every node monitored.

Output: shoulder_check.jsonl (one row per run, with the per-node failure times), shoulder_check.json/.md (summary).
Run: python runs/final_review/shoulder_check.py [run|summary]   (UHTC_WORKERS processes, at most 4)
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "runs"))
sys.path.insert(0, str(ROOT / "runs" / "envelope_2d"))
import common  # noqa: E402
from uhtc import atmosphere as atm, coupled as cp, materials as mt  # noqa: E402
from uhtc.conduction2d import NoseGeometry  # noqa: E402
from uhtc.oxidation import LayerOxidation, OxParams  # noqa: E402

OUT = HERE / "shoulder_check.jsonl"
MODES = ("local", "stag")
STEP = 0.25


def stop_hot(events, T_s, t):
    """Run to the end; stop only when the stagnation surface leaves the property tables (as common.stop_first)."""
    return "T > abort limit" if T_s > common.T_ABORT else None


def job_id(j):
    if j["kind"] == "x43a":
        return f"{j['mode']}_x43a_R{j['R_n'] * 1e3:g}_{j['catalysis']}"
    return f"{j['mode']}_R{j['R_n'] * 1e3:g}_{j['catalysis']}_M{j['mach']:g}_z{j['alt'] / 1e3:g}"


def run_job(j):
    t0 = time.time()
    R, cat = j["R_n"], j["catalysis"]
    ox = LayerOxidation(OxParams(f_SiC=0.2))
    g = NoseGeometry(R_n=R, theta_c=math.radians(5.0), L_flank=0.03, axisymmetric=True)
    if j["kind"] == "const":
        state = common.stag_state(j["mach"], j["alt"])
        cond = cp.ConstantConditions(j["mach"], j["alt"], R, state=state, catalysis=cat)
        t_end, dt_max, rec = common.T_CAP, 2.0, 5.0  # production settings of common.run_case_2d
        eps = min(float(atm.state(j["alt"])[1]) / state.p_s, 1.0)
    else:
        a = np.loadtxt(ROOT / "runs" / "trajectory" / "x43a_mach10.csv", delimiter=",")
        cond = cp.TrajectoryConditions(a[:, 0] - a[0, 0], a[:, 1], a[:, 2], R, dt_eval=2.0,
                                       state_fn=common.stag_state, catalysis=cat)
        t_end, dt_max, rec = float(a[-1, 0] - a[0, 0]), 1.0, 1.0  # settings of runs/q5_2d/run_q5.py
        eps = math.nan
    try:
        r = cp.run_2d(cond, g, mt.composite(0.2), ox, t_end, record_every=rec, dt_max=dt_max, mesh_kw=common.MESH_2D,
                      stop_on=stop_hot, heat_dist="lees", local_pressure=(j["mode"] == "local"))
    except Exception as e:  # recorded, never dropped
        return dict(case_id=job_id(j), **j, error=repr(e)[:200], cpu_s=round(time.time() - t0, 1))
    sol, states = r.extra["solver"], r.extra["states"]
    nodes = []
    for i, s in enumerate(states):
        tf, mode = cp.first_failure(ox.events(s))
        nodes.append(dict(i=i, s_mm=float(sol.s[i] * 1e3), phi_deg=float(np.degrees(sol.phi[i])),
                          flank=bool(sol.on_flank[i]),
                          pr=float((1 - eps) * math.cos(sol.phi[i]) ** 2 + eps) if eps == eps else math.nan,
                          t_glass=s["t_glass_loss"], t_active=s["t_active"], t_limit=s["t_limit"], t_fail=tf,
                          mode=mode, Ts_max=float(s["Ts_max"]), glass_um_end=ox.glass_thickness(s) * 1e6,
                          t_end_node=float(s["t"])))
    return dict(case_id=job_id(j), **j, t_end=float(r.t[-1]) if len(r.t) else math.nan, t_end_target=t_end,
                stopped=r.stopped, eps=eps, nodes=nodes, error="", cpu_s=round(time.time() - t0, 1))


def load():
    if not OUT.exists():
        return {}
    return {r["case_id"]: r for r in (json.loads(l) for l in OUT.read_text(encoding="utf-8").splitlines() if l.strip())}


def run_jobs(jobs):
    done = load()
    todo = [j for j in jobs if job_id(j) not in done]
    print(f"{OUT.name}: {len(jobs)} jobs, {len(jobs) - len(todo)} done, {len(todo)} to run", flush=True)
    if not todo:
        return
    with Pool(min(int(os.environ.get("UHTC_WORKERS", "4")), 4)) as pool:
        t0 = time.time()
        with open(OUT, "a", encoding="utf-8") as f:
            for k, res in enumerate(pool.imap_unordered(run_job, todo, chunksize=1)):
                f.write(json.dumps(res, default=float) + "\n")
                f.flush()
                print(f"  {k + 1}/{len(todo)} {res['case_id']} {time.time() - t0:.0f} s", flush=True)


def baseline_columns():
    import run_envelope_2d as env
    rows = common.load(env.HERE / "grid.jsonl") + common.load(env.HERE / "refine.jsonl")
    cols = env.columns(rows)
    return {(R, z): cols[(R, "fc", z)] for R in env.RADII for z in env.BAND_ALTS}


def first_any(row, exclude_stag=False, t_max=common.T_CAP):
    """Earliest failure (t, mode, node) over the nodes of a run, within t_max; node 0 is the stagnation node."""
    best = (math.nan, "none", None)
    for n in row["nodes"]:
        if exclude_stag and n["i"] == 0:
            continue
        t = n["t_fail"]
        if t == t and t <= t_max and (best[0] != best[0] or t < best[0]):
            best = (t, n["mode"], n)
    return best


def const_job(mode, R, z, M):
    return dict(kind="const", mode=mode, R_n=R, catalysis="fc", mach=M, alt=z)


def main_run():
    cols = baseline_columns()
    jobs = []
    for mode in MODES:
        jobs += [const_job(mode, R, z, Ms) for (R, z), (Ms, Mf) in sorted(cols.items())]
    for mode in MODES:
        jobs += [const_job(mode, R, z, Mf) for (R, z), (Ms, Mf) in sorted(cols.items())]
    jobs += [dict(kind="x43a", mode=mode, R_n=R, catalysis=cat, mach=0.0, alt=0.0)
             for mode in MODES for R in (0.001, 0.003, 0.010) for cat in ("fc", "nc")]
    run_jobs(jobs)
    # every-node boundary: lower M in 0.25 steps where any node fails within 1800 s at the current highest pass.
    # Local mode only: in the 'stag' mode the flank nodes meet criterion (a) within the first minute because the
    # stagnation O2 pressure, applied where the boundary layer delivers far less O2, inflates their critical thickness
    # (Eq. 7); that is an artefact of the production assumption, so no boundary is searched for it (see summary).
    for _ in range(12):
        done = load()
        more = []
        for mode in ("local",):
            for (R, z), (Ms, Mf) in cols.items():
                M = Ms
                while True:
                    row = done.get(job_id(const_job(mode, R, z, M)))
                    if row is None:
                        more.append(const_job(mode, R, z, M))
                        break
                    if first_any(row)[0] != first_any(row)[0] or M - STEP < 5.0:
                        break  # no node fails: every-node boundary found
                    M -= STEP
        if not more:
            break
        run_jobs(more)


PA_LAWS = {"baseline": dict(), "refit": dict(pa_law="jacobson_fe"), "refit_off0.3": dict(pa_law="jacobson_fe", pa_offset=0.3),
           "balat": dict(pa_law="balat2015")}


def pa_margins(row):
    """Criterion (b) a posteriori at every node of a constant-condition run: margin T_transition(p_O2 at the node) minus
    the node's highest surface temperature, for the baseline line and the alternative lines (criterion (b) does not feed
    back into the solution, and p_O2 is constant in time, so the node turns active iff the margin is negative).  In the
    'stag' mode every node sees the stagnation p_O2; in 'local' mode p_O2 x p/p_s."""
    st = common.stag_state(row["mach"], row["alt"])
    out = {}
    for name, kw in PA_LAWS.items():
        ox = LayerOxidation(OxParams(f_SiC=0.2, **kw))
        m = []
        for n in row["nodes"]:
            p = st.p_O2_w * (n["pr"] if row["mode"] == "local" else 1.0)
            m.append((ox.T_transition(p) - n["Ts_max"], n["phi_deg"], n["i"]))
        off = min(m[1:])
        pr_min = min(n["pr"] for n in row["nodes"]) if row["mode"] == "local" else math.nan
        drop = ox.T_transition(st.p_O2_w) - ox.T_transition(st.p_O2_w * pr_min) if pr_min == pr_min else math.nan
        out[name] = dict(stag=m[0][0], off_min=off[0], off_phi=off[1], any_active=min(x[0] for x in m) < 0,
                         tr_drop_flank=drop)
    return out


def summary():
    cols = baseline_columns()
    done = load()
    out = dict(columns=[], x43a=[])
    for mode in MODES:
        for (R, z), (Ms, Mf) in sorted(cols.items()):
            rs, rf = done.get(job_id(const_job(mode, R, z, Ms))), done.get(job_id(const_job(mode, R, z, Mf)))
            if rs is None or rf is None:
                continue
            ts, ms, ns_ = first_any(rs)
            stag_f = rf["nodes"][0]
            tf_off, mf_off, nf_off = first_any(rf, exclude_stag=True)
            M = Ms
            while True:  # every-node boundary
                row = done.get(job_id(const_job(mode, R, z, M)))
                if row is None or first_any(row)[0] != first_any(row)[0]:
                    break
                M -= STEP
            out["columns"].append(dict(
                mode=mode, R_mm=R * 1e3, z_km=z / 1e3, M_s=Ms, M_f=Mf,
                Ms_any_node_fails=ts == ts, Ms_first_t=ts, Ms_first_mode=ms,
                Ms_first_phi=None if ns_ is None else ns_["phi_deg"], Ms_first_s_mm=None if ns_ is None else ns_["s_mm"],
                Mf_stag_t=stag_f["t_fail"], Mf_stag_mode=stag_f["mode"], Mf_off_t=tf_off, Mf_off_mode=mf_off,
                Mf_off_phi=None if nf_off is None else nf_off["phi_deg"],
                Mf_off_before_stag=bool(tf_off == tf_off and tf_off < stag_f["t_fail"]),
                M_s_every_node=M if (row is not None) else None, shift=(M - Ms) if row is not None else None,
                Ms_Tmax_off=max(n["Ts_max"] for n in rs["nodes"][1:]), Ms_Tmax_stag=rs["nodes"][0]["Ts_max"],
                Ms_min_pr=min(n["pr"] for n in rs["nodes"]), Ms_pa_margins=pa_margins(rs),
                errors=[rs["error"], rf["error"]]))
        for R in (0.001, 0.003, 0.010):
            for cat in ("fc", "nc"):
                r = done.get(job_id(dict(kind="x43a", mode=mode, R_n=R, catalysis=cat)))
                if r is None:
                    continue
                t, m, n = first_any(r, t_max=math.inf)
                out["x43a"].append(dict(mode=mode, R_mm=R * 1e3, catalysis=cat, any_node_fails=t == t, first_t=t,
                                        first_mode=m, first_phi=None if n is None else n["phi_deg"],
                                        stag_t=r["nodes"][0]["t_fail"], stag_mode=r["nodes"][0]["mode"],
                                        Tmax_stag=r["nodes"][0]["Ts_max"], error=r["error"]))
    (HERE / "shoulder_check.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    lines = ["# Shoulder check (external review A4, run R1)", "",
             "Generated by `runs/final_review/shoulder_check.py summary` from `shoulder_check.jsonl`. "
             "Mode `stag`: production model (stagnation pressure at every node); `local`: local wall and O2 pressure "
             "(modified Newtonian) in each node's oxide model. Times in s, angles in degrees from the axis.", "",
             "| mode | R (mm) | z (km) | M_s | any node fails at M_s | first (t, mode, phi) | M_f stagnation (t, mode) "
             "| M_f first off-stagnation (t, mode, phi) | every-node boundary | shift |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    f = lambda v: "-" if v is None or v != v else (f"{v:.0f}" if isinstance(v, float) else str(v))
    for c in out["columns"]:
        lines.append(f"| {c['mode']} | {c['R_mm']:g} | {c['z_km']:g} | {c['M_s']:g} | {c['Ms_any_node_fails']} | "
                     f"{f(c['Ms_first_t'])}, {c['Ms_first_mode']}, {f(c['Ms_first_phi'])} | "
                     f"{f(c['Mf_stag_t'])}, {c['Mf_stag_mode']} | {f(c['Mf_off_t'])}, {c['Mf_off_mode']}, "
                     f"{f(c['Mf_off_phi'])} | {'-' if c['M_s_every_node'] is None else format(c['M_s_every_node'], 'g')} | "
                     f"{'-' if c['shift'] is None else format(c['shift'], 'g')} |")
    lines += ["", "Criterion (b) a posteriori at M_s: transition temperature minus the highest surface temperature (K), at the "
              "stagnation node / lowest over the other nodes (body angle).", "",
              "| mode | R (mm) | z (km) | " + " | ".join(PA_LAWS) + " |", "|---|---|---|" + "---|" * len(PA_LAWS)]
    for c in out["columns"]:
        pm = c["Ms_pa_margins"]
        lines.append(f"| {c['mode']} | {c['R_mm']:g} | {c['z_km']:g} | " + " | ".join(
            f"{pm[k]['stag']:.0f} / {pm[k]['off_min']:.0f} ({pm[k]['off_phi']:.0f})" for k in PA_LAWS) + " |")
    lines += ["", "| mode | R (mm) | catalysis | any node fails | first (t, mode, phi) | stagnation (t, mode) | "
              "T_max stagnation (K) |", "|---|---|---|---|---|---|---|"]
    for c in out["x43a"]:
        lines.append(f"| {c['mode']} | {c['R_mm']:g} | {c['catalysis']} | {c['any_node_fails']} | {f(c['first_t'])}, "
                     f"{c['first_mode']}, {f(c['first_phi'])} | {f(c['stag_t'])}, {c['stag_mode']} | "
                     f"{f(c['Tmax_stag'])} |")
    (HERE / "shoulder_check.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "run"
    if what == "run":
        main_run()
        summary()
    elif what == "summary":
        summary()
    else:
        raise SystemExit(__doc__)
