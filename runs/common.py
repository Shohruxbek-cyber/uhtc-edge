"""Shared case runner for the production runs (runs/*.py).  Every run records the three failure times
(glass loss, SiC active oxidation, temperature limit), the first failure and its mode, peak temperatures and the
layer thicknesses at the first failure and at the end.  Parallel execution uses at most 4 worker processes
(standing rule 6).  All outputs are small CSV/JSON files in runs/<study>/.
"""
from __future__ import annotations

import csv
import dataclasses
import json
import math
import os
import pickle
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from uhtc import coupled as cp, materials as mt, stagnation as st  # noqa: E402
from uhtc.oxidation import LayerOxidation, OxParams  # noqa: E402

CACHE = ROOT / "runs" / "_cache"
T_CAP = 1800.0  # s, exposure cap (DECISIONS D-FAIL-2)
T_ABORT = 3400.0  # K, stop: all criteria long exceeded, beyond the property tables


def stag_state(mach, alt):
    """Stagnation state with a disk cache (pickle per (Mach, altitude))."""
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"stag_M{mach:.4f}_z{alt:.1f}.pkl"
    if f.exists():
        with open(f, "rb") as h:
            return pickle.load(h)
    s = st.stagnation_state(mach, alt)
    tmp = f.with_suffix(f".{os.getpid()}.tmp")
    with open(tmp, "wb") as h:
        pickle.dump(s, h)
    os.replace(tmp, f)
    return s


def stop_all_or_hot(events, T_s, t):
    """Stop when all three failure times are known, or the surface is hotter than T_ABORT."""
    if all(v == v for v in events.values()):
        return "all criteria met"
    if T_s > T_ABORT:
        return "T > abort limit"
    return None


def stop_first(events, T_s, t):
    if any(v == v for v in events.values()):
        return "first failure"
    if T_s > T_ABORT:
        return "T > abort limit"
    return None


def run_case(case: dict):
    """case keys: mach, alt, R_n, catalysis, f_SiC, and optional overrides:
    ox (dict of OxParams fields), k_mult, cp_mult, heat_mult, mtf, t_cap, stop ('all'|'first'), depth, cylinder."""
    t0 = time.time()
    p = OxParams(f_SiC=case.get("f_SiC", 0.2), **case.get("ox", {}))
    ox = LayerOxidation(p)
    mat = mt.composite(case.get("f_SiC", 0.2), k_mult=case.get("k_mult", 1.0), cp_mult=case.get("cp_mult", 1.0))
    state = stag_state(case["mach"], case["alt"])
    cond = cp.ConstantConditions(case["mach"], case["alt"], case["R_n"], state=state, catalysis=case["catalysis"],
                                 heat_mult=case.get("heat_mult", 1.0), mass_transfer_factor=case.get("mtf", 1.0),
                                 cylinder=case.get("cylinder", False))
    stop = stop_first if case.get("stop", "all") == "first" else stop_all_or_hot
    try:
        r = cp.run_1d(cond, mat, ox, case.get("t_cap", T_CAP), depth=case.get("depth", 0.02), stop_on=stop)
        err = ""
    except Exception as e:  # recorded, never silently dropped
        return dict(case_id=case.get("id", ""), **_flat(case), error=repr(e)[:200], cpu_s=round(time.time() - t0, 2))
    tf, mode = cp.first_failure(r.events)
    i_f = int(np.searchsorted(r.t, tf)) if tf == tf else len(r.t) - 1
    i_f = min(i_f, len(r.t) - 1)
    Lf = ox.layers(r.states[i_f])
    Le = ox.layers(r.states[-1])
    out = dict(case_id=case.get("id", ""), **_flat(case),
               p_s_kPa=state.p_s / 1e3, H_s_MJkg=state.H_s / 1e6, hD_frac=state.h_D / state.H_s, T_e=state.T_e,
               p_O2_kPa=state.p_O2_w / 1e3, T_transition=ox.T_transition(state.p_O2_w),
               t_glass=r.events["t_glass_loss"], t_active=r.events["t_active"], t_limit=r.events["t_limit"],
               t_fail=tf, mode=mode, T_s_max=float(np.max(r.T_s)), T_s_end=float(r.T_s[-1]), t_end=float(r.t[-1]),
               T_s_at_fail=float(r.T_s[i_f]) if tf == tf else math.nan,
               glass_um_fail=Lf["glass"] * 1e6, ZrO2_um_fail=Lf["ZrO2"] * 1e6, dep_um_fail=Lf["depleted"] * 1e6,
               glass_um_end=Le["glass"] * 1e6, ZrO2_um_end=Le["ZrO2"] * 1e6, dep_um_end=Le["depleted"] * 1e6,
               stopped=r.stopped, energy_res=r.energy_residual, error=err, cpu_s=round(time.time() - t0, 2))
    return out


MESH_2D = dict(n_arc=30, n_flank=40, n_eta=30, eta1=0.003)  # production 2D mesh (convergence: runs/envelope_2d)


def run_case_2d(case: dict):
    """2D nose region (uhtc.conduction2d) for one condition.  case keys as run_case plus: geometry ('nose' axisymmetric
    sphere-cone | 'wing' planar cylinder-wedge), theta_deg (5), L_flank (0.03 m), heat_dist (heating-distribution name,
    uhtc.conduction2d.HEAT_DISTRIBUTIONS), flank_mult (1), mesh (dict, default MESH_2D), dt_max (2 s).  Failure clocks
    and layers are those of the stagnation node; T_body_max is the hottest surface node."""
    from uhtc.conduction2d import NoseGeometry
    t0 = time.time()
    p = OxParams(f_SiC=case.get("f_SiC", 0.2), **case.get("ox", {}))
    ox = LayerOxidation(p)
    mat = mt.composite(case.get("f_SiC", 0.2), k_mult=case.get("k_mult", 1.0), cp_mult=case.get("cp_mult", 1.0))
    state = stag_state(case["mach"], case["alt"])
    if case.get("oxidant") == "O_as_O2":  # D-OX-14: sensitivity case, every edge O atom counted as one O2 molecule
        state = dataclasses.replace(state, p_O2_w=max(state.p_O2_w, (state.x_e["O2"] + state.x_e["O"]) * state.p_s))
    axi = case.get("geometry", "nose") == "nose"
    cond = cp.ConstantConditions(case["mach"], case["alt"], case["R_n"], state=state, catalysis=case["catalysis"],
                                 heat_mult=case.get("heat_mult", 1.0), mass_transfer_factor=case.get("mtf", 1.0),
                                 cylinder=not axi)
    g = NoseGeometry(R_n=case["R_n"], theta_c=math.radians(case.get("theta_deg", 5.0)),
                     L_flank=case.get("L_flank", 0.03), axisymmetric=axi)
    stop = stop_first if case.get("stop", "first") == "first" else stop_all_or_hot
    flow_kw = {}
    if case.get("glass_flow", "off") == "on":  # version 3, Task 1 (uhtc.glassflow); absent or "off" in every run of the paper
        from uhtc.glassflow import FlowParams
        flow_kw["glass_flow"] = FlowParams(law=case.get("flow_visc", "LA"), tau_mult=case.get("flow_tau_mult", 1.0),
                                           pressure=case.get("flow_pressure", True))
    lp_kw = {"local_pressure": True} if case.get("local_pressure", False) else {}  # journal revision, ASK 5 check only
    try:
        r = cp.run_2d(cond, g, mat, ox, case.get("t_cap", T_CAP), record_every=case.get("record_every", 5.0),
                      dt_max=case.get("dt_max", 2.0), mesh_kw=case.get("mesh", MESH_2D), stop_on=stop,
                      flank_mult=case.get("flank_mult", 1.0), heat_dist=case.get("heat_dist", "tauber"), **flow_kw, **lp_kw)
        err = ""
    except Exception as e:  # recorded, never silently dropped
        return dict(case_id=case.get("id", ""), **_flat(case), error=repr(e)[:200], cpu_s=round(time.time() - t0, 2))
    tf, mode = cp.first_failure(r.events)
    Le = ox.layers(r.states[-1])
    extra = {}
    if flow_kw:
        gfx = r.extra["glass_flow"]
        led = gfx["ledger"]
        extra = dict(flow_audit_max=gfx["audit_max"], flow_shed=led.shed, flow_source=led.source, flow_evap=led.evap,
                     flow_steps_explicit=led.n_explicit, flow_steps_implicit=led.n_implicit, flow_max_substeps=led.max_sub)
    return dict(case_id=case.get("id", ""), **_flat(case), **extra,
                p_s_kPa=state.p_s / 1e3, H_s_MJkg=state.H_s / 1e6, hD_frac=state.h_D / state.H_s, T_e=state.T_e,
                p_O2_kPa=state.p_O2_w / 1e3, T_transition=ox.T_transition(state.p_O2_w),
                t_glass=r.events["t_glass_loss"], t_active=r.events["t_active"], t_limit=r.events["t_limit"],
                t_fail=tf, mode=mode, T_s_max=float(np.max(r.T_s)), T_s_end=float(r.T_s[-1]), t_end=float(r.t[-1]),
                T_body_max=float(np.max(r.T_sub)) if len(r.T_sub) else math.nan,
                glass_um_end=Le["glass"] * 1e6, ZrO2_um_end=Le["ZrO2"] * 1e6, dep_um_end=Le["depleted"] * 1e6,
                stopped=r.stopped, energy_res=r.energy_residual, error=err, cpu_s=round(time.time() - t0, 2))


def _flat(case):
    d = {k: v for k, v in case.items() if k not in ("ox", "id", "mesh")}
    for k, v in case.get("ox", {}).items():
        d["ox_" + k] = v
    return d


def run_all(cases, out_jsonl: Path, workers=int(os.environ.get("UHTC_WORKERS", "4")), runner=None):
    """Run cases in parallel (<= 4 processes); results appended to a JSON-lines file (restartable: finished case ids
    are skipped), then exported to CSV next to it."""
    out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if out_jsonl.exists():
        done = {json.loads(l)["case_id"] for l in out_jsonl.read_text(encoding="utf-8").splitlines() if l.strip()}
    todo = [c for c in cases if c["id"] not in done]
    print(f"{out_jsonl.name}: {len(cases)} cases, {len(done)} done, {len(todo)} to run", flush=True)
    if todo:
        from multiprocessing import Pool
        uniq = sorted({(c["mach"], c["alt"]) for c in todo})
        with Pool(min(workers, 4)) as pool:
            pool.starmap(stag_state, uniq, chunksize=4)
            t0 = time.time()
            with open(out_jsonl, "a", encoding="utf-8") as f:
                fn = runner or run_case
                for i, res in enumerate(pool.imap_unordered(fn, todo, chunksize=2 if fn is run_case else 1)):
                    f.write(json.dumps(res, default=float) + "\n")
                    f.flush()
                    if (i + 1) % 50 == 0:
                        print(f"  {i + 1}/{len(todo)}  {time.time() - t0:.0f} s", flush=True)
    export_csv(out_jsonl)


def load(out_jsonl: Path):
    return [json.loads(l) for l in out_jsonl.read_text(encoding="utf-8").splitlines() if l.strip()]


def export_csv(out_jsonl: Path):
    rows = load(out_jsonl)
    fields = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with open(out_jsonl.with_suffix(".csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in sorted(rows, key=lambda r: r["case_id"]):
            w.writerow(r)


def write_meta(path: Path, **meta):
    meta = dict(meta)
    meta["oxparams_baseline"] = {k: v for k, v in asdict(OxParams()).items() if k != "sources"}
    meta["T_cap_s"] = T_CAP
    meta["T_abort_K"] = T_ABORT
    path.write_text(json.dumps(meta, indent=1, default=str), encoding="utf-8")
