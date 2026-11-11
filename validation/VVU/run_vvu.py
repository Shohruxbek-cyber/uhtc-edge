"""VVU: validation uncertainty of V7-V9 (ASME V&V 20 style). Plan: validation/VVU/VVU_plan.md, committed before any run.

Post hoc: added after the first results of V6-V9. The preset verdicts are unchanged and nothing in src/uhtc changes.
The run functions, constants and conditions of run_v7.py, run_v8.py and run_v9.py are imported, not copied; the thin
wrappers below only add the sampled inputs (OxParams fields, material multipliers, the mass-transfer multiplier, the
heating multiplier of V7 and the calorimeter correction of V8/V9). Sample 0 of every case is the baseline and must
reproduce the stored results of the checks (plan Section 7).

Usage: python validation/VVU/run_vvu.py gate        reproducibility gate and the V7 heating ratio r
       python validation/VVU/run_vvu.py time        time three samples of each case (plan Section 8)
       python validation/VVU/run_vvu.py run [n7]    all samples (n7 = joint N of V7, 100 or 50 by the cost rule)
Output: vvu_runs.jsonl (one line per run, resumable), then analyse_vvu.py writes vvu_samples.csv, vvu_results.json, report.md.
"""
import json
import math
import os
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
for d in ("V7_arcjet", "V8_measured_flux", "V9_plasmatron"):
    sys.path.insert(0, str(ROOT / "validation" / d))
import run_v7  # noqa: E402
import run_v8  # noqa: E402
import run_v9  # noqa: E402
from uhtc import coupled as cp, materials as mt, thermo as th  # noqa: E402
from uhtc.conduction2d import NoseGeometry  # noqa: E402
from uhtc.oxidation import LayerOxidation, OxParams  # noqa: E402

OUT = HERE / "vvu_runs.jsonl"
SEEDS = {"joint": 20261110, "bc": 20261111, "input": 20261112}

# group (c), plan Section 5.3: name -> (low, base, high, scale); base = the OxParams / material default
INPUTS = {
    "eps_glass": (0.70, 0.88, 0.95, "lin"),
    "eps_oxide": (0.50, 0.75, 0.85, "lin"),
    "k_mult": (0.70, 1.0, 1.05, "lin"),
    "cp_mult": (0.97, 1.0, 1.03, "lin"),
    "kp_mult": (1 / 3, 1.0, 3.0, "log"),
    "kp_mult_hi": (0.5, 1.0, 3.0, "log"),
    "kd_mult": (0.5, 1.0, 2.0, "log"),
    "T_dep": (1875.15, 1627.0 + 273.15, 1925.15, "lin"),  # base computed as in data.T_DEP_ONSET
    "pa_offset": (0.3, 0.65, 1.0, "lin"),
    "mtf": (0.5, 1.0, 2.0, "log"),
    "k_glass": (0.7, 1.37, 2.8, "lin"),
    "k_ZrO2": (1.0, 1.5, 2.5, "lin"),
    "k_dep": (10.0, 30.0, 60.0, "lin"),
    "n_p": (0.0, 0.0, 0.5, "lin"),
    "f_SiO_lost": (0.0, 0.0, 1.0, "lin"),
    "alpha_evap": (0.1, 1.0, 1.0, "log"),
    "perm_case": (None, "base", None, "discrete"),
}
OX_KEYS = ("eps_glass", "eps_oxide", "kp_mult", "kp_mult_hi", "kd_mult", "T_dep", "pa_offset", "k_glass", "k_ZrO2", "k_dep",
           "n_p", "f_SiO_lost", "alpha_evap", "perm_case")
PERM = ("low", "base", "high")
CASES = ("V7", "V8_1", "V8_2", "V9")


def bc_inputs(case, r7=None):
    """Group (b), plan Section 5.2: name -> (low, base, high, scale)."""
    if case == "V7":

        lo, hi = sorted((r7, 1.0 / r7))  # Amendment 1: log-uniform from 1/r to r
        return {"heat_bc": (lo, 1.0, hi, "log")}
    g = {"gamma_cu": (0.01, None, 1.0, "log")}  # base None = no calorimeter correction (the checks as run)
    if case == "V8_1":
        g["H_MJ"] = (15.0, 30.0, 45.0, "lin")
    return g


def lhs(spec, n, seed):
    """Latin hypercube: one stratified uniform per input, transformed by its scale."""
    rng = np.random.default_rng(seed)
    names = list(spec)
    u = (rng.permuted(np.tile(np.arange(n), (len(names), 1)), axis=1).T + rng.random((n, len(names)))) / n
    rows = [dict() for _ in range(n)]
    for j, k in enumerate(names):
        lo, _, hi, sc = spec[k]
        for i in range(n):
            if sc == "log":
                v = lo * (hi / lo) ** u[i, j]
            elif sc == "discrete":
                v = PERM[min(2, int(3 * u[i, j]))]
            else:
                v = lo + (hi - lo) * u[i, j]
            rows[i][k] = float(v) if sc != "discrete" else v
    return rows


def designs(case, n_joint, r7=None):
    """{set name: [input dicts]}; every dict holds the full set of inputs (baseline where not sampled)."""
    b, c = bc_inputs(case, r7), INPUTS
    base = {k: v[1] for k, v in {**c, **b}.items()}
    out = {"base": [dict(base)]}
    for name, spec, n in (("joint", {**b, **c}, n_joint), ("bc", b, 50), ("input", c, 50)):
        rows = lhs(spec, n, SEEDS[name] + CASES.index(case))
        out[name] = [{**base, **r} for r in rows]
    return out


def ox_params(f_SiC, x):
    return OxParams(f_SiC=f_SiC, **{k: x[k] for k in OX_KEYS})


def material(f_SiC, x):
    return mt.composite(f_SiC, k_mult=x["k_mult"], cp_mult=x["cp_mult"])


def calorimeter_ratio(H, p, q_cw, gamma, T_cal):
    """q(gamma)/q_fc of a copper calorimeter: frozen boundary layer, first-order wall recombination (plan Section 5.2)."""
    f = cp.PrescribedConditions([(0.0, H, p, q_cw)], R_n=0.01).flows[0]
    rho_w = p * 0.028965 / (8.314462618 * T_cal)
    k_w = gamma * math.sqrt(8.314462618 * T_cal / (2.0 * math.pi * 0.016))
    phi = 1.0 / (1.0 + f.Phi / (rho_w * k_w))
    return f.f_nc + phi * (1.0 - f.f_nc), f.f_nc, phi


def corrected(H, p, q_cw, x, T_cal=300.0):
    if x.get("gamma_cu") is None:
        return q_cw
    return q_cw / calorimeter_ratio(H, p, q_cw, x["gamma_cu"], T_cal)[0]


# ---------------------------------------------------------------- the four cases, one draw x each
def run_v7_case(cat, x):
    steps = [(t - run_v7.STARTS[0], H, p) for t, H, p in zip(run_v7.STARTS, run_v7.H0, run_v7.P02)]
    cond = cp.PrescribedConditions(steps, run_v7.RN, x_O2=run_v7.X_O2, catalysis=cat, heat_mult=x["heat_bc"],
                                   mass_transfer_factor=x["mtf"])
    ox = LayerOxidation(ox_params(0.15, x))
    g = NoseGeometry(R_n=run_v7.RN, theta_c=0.0, L_flank=5e-4, axisymmetric=True)
    r = cp.run_2d(cond, g, material(0.15, x), ox, run_v7.T_OFF - run_v7.STARTS[0], dt_max=1.0, record_every=1.0,
                  mesh_kw=dict(n_arc=30, n_flank=3, n_eta=30, eta1=0.004))
    t = r.t + run_v7.STARTS[0]
    ends = [run_v7.STARTS[i + 1] - 1.0 for i in range(4)] + [run_v7.T_OFF - 1.0]
    L = ox.layers(r.states[-1])
    return dict(Tend=[float(np.interp(te, t, r.T_s)) for te in ends], oxide_um=float((L["ZrO2"] + L["glass"]) * 1e6),
                dep_um=float(L["depleted"] * 1e6))


def run_1d(q_cw, H, p, t_end, f_SiC, depth, cat, x):
    """run_v8.run with the sampled inputs."""
    ox = LayerOxidation(ox_params(f_SiC, x))
    cond = cp.PrescribedConditions([(0.0, H, p, q_cw)], R_n=0.01, catalysis=cat, mass_transfer_factor=x["mtf"])
    r = cp.run_1d(cond, material(f_SiC, x), ox, t_end, depth=depth, record_every=1.0)
    L = [ox.layers(s) for s in r.states]
    return dict(t=r.t, T=r.T_s, oxide=np.array([(l["ZrO2"] + l["glass"]) * 1e6 for l in L]),
                dep=np.array([l["depleted"] * 1e6 for l in L]))


def run_v8_1_case(cat, x):
    H, p = x["H_MJ"] * 1e6, 0.042 * 101325.0
    a = run_1d(corrected(H, p, 3.0e6, x), H, p, 300.0, 0.2, 0.00635, cat, x)
    return dict(T_max=float(np.max(a["T"])), oxide_um=float(a["oxide"][-1]), dep_um=float(a["dep"][-1]))


def run_v8_2_case(cat, x):
    def one(q, H, t_end, f):
        return run_1d(corrected(H, 18e3, q, x), H, 18e3, t_end, f, 0.03, cat, x)
    c1, c2 = one(4.78e6, 27.9e6, 600.0, 0.3), one(3.80e6, 20.8e6, 600.0, 0.3)
    z10, z15 = one(4.78e6, 27.9e6, 300.0, 0.1), one(4.78e6, 27.9e6, 600.0, 0.15)
    at = lambda r, k, tt: float(np.interp(tt, r["t"], r[k]))  # noqa: E731
    return dict(T1=float(c1["T"][-1]), T2=float(c2["T"][-1]),
                oxide_c1=[at(c1, "oxide", tt) for tt in (201, 402, 600)],
                reacted_c1=[at(c1, "oxide", tt) + at(c1, "dep", tt) for tt in (201, 402, 600)],
                oxide_c2=float(c2["oxide"][-1]), T_zs10=float(z10["T"][-1]), T_zs15=float(z15["T"][-1]),
                trend=bool(z10["T"][-1] > z15["T"][-1] > c1["T"][-1]))


def run_v9_case(cat, x):
    res = {}
    for name, (ts, lo, hi, _, _, t_ref, _) in run_v9.CASES.items():
        steps = []
        for t0, (q, H, pd) in ((0.0, lo), (ts, hi)):
            p = run_v9.P_STAT + pd
            steps.append((t0, H, p, corrected(H, p, q, x, T_cal=350.0)))
        ox = LayerOxidation(ox_params(run_v9.F_SIC, x))
        cond = cp.PrescribedConditions(steps, R_n=0.01, catalysis=cat, mass_transfer_factor=x["mtf"])
        r = cp.run_1d(cond, material(run_v9.F_SIC, x), ox, run_v9.T_END, depth=run_v9.DEPTH, record_every=1.0)
        t, T = np.array(r.t), np.array(r.T_s)
        tg = float(r.events["t_glass_loss"])
        res[name] = dict(T_low=float(np.interp(ts - 1.0, t, T)), T_ref=float(np.interp(t_ref, t, T)), t_glass=tg,
                         after=(tg - ts) if not math.isnan(tg) else None, lost_in_high=bool(not math.isnan(tg) and tg > ts))
    return res


RUN = {"V7": run_v7_case, "V8_1": run_v8_1_case, "V8_2": run_v8_2_case, "V9": run_v9_case}


def task(args):
    case, cat, sset, i, x = args
    t0 = time.time()
    try:
        out, err = RUN[case](cat, x), None
    except Exception as e:  # a numerical failure is reported and left out (plan Section 6)
        out, err = None, f"{type(e).__name__}: {e}"
    return dict(case=case, bound=cat, set=sset, i=i, x=x, out=out, error=err, cpu_s=round(time.time() - t0, 2))


def v7_ratio():
    """r = CFD flux of step 5 (5 MW/m2, non-catalytic cold wall, 300 K) / model non-catalytic cold-wall flux (plan 5.2)."""
    H, p = run_v7.H0[-1], run_v7.P02[-1]
    f = cp.PrescribedConditions([(0.0, H, p)], run_v7.RN, x_O2=run_v7.X_O2).flows[0]
    q_nc = f.f_nc * f.Phi * (H - float(th.h_air_frozen(300.0)))
    return 5.0e6 / q_nc, q_nc, f.f_nc


def close(a, b):
    return abs(a - b) <= 1e-9 * max(1.0, abs(b))


def gate():
    """Plan Section 7: sample 0 reproduces the stored results."""
    v7 = json.loads((ROOT / "validation/V7_arcjet/v7_results.json").read_text())
    v8 = json.loads((ROOT / "validation/V8_measured_flux/v8_results.json").read_text())
    v9 = json.loads((ROOT / "validation/V9_plasmatron/v9_results.json").read_text())
    base = designs("V7", 1, r7=1.0)["base"][0]
    bad = []
    jobs = []
    for cat in ("fc", "nc"):
        jobs.append(("V7", cat, "gate", 0, dict(base)))
        for H in (15.0, 30.0, 45.0):
            jobs.append(("V8_1", cat, "gate", int(H), {**base, "H_MJ": H, "gamma_cu": None}))
        jobs.append(("V8_2", cat, "gate", 0, {**base, "gamma_cu": None}))
        jobs.append(("V9", cat, "gate", 0, {**base, "gamma_cu": None}))
    with ProcessPoolExecutor(max_workers=7) as ex:
        res = list(ex.map(task, jobs))
    for r in res:
        o, cat = r["out"], r["bound"]
        if r["case"] == "V7":
            s = v7[f"2d-{cat}"]
            pairs = list(zip(o["Tend"], s["Tend"])) + [(o["oxide_um"], s["oxide_um"]), (o["dep_um"], s["dep_um"])]
        elif r["case"] == "V8_1":
            s = v8[f"V8_1_H{r['i']}_{cat}"]["passive"]
            pairs = [(o["T_max"], s["T_max"]), (o["oxide_um"], s["oxide_um"][-1]), (o["dep_um"], s["dep_um"][-1])]
        elif r["case"] == "V8_2":
            s = v8[f"V8_2_{cat}"]
            pairs = [(o["T1"], s["c1"]["T_end"]), (o["T2"], s["c2"]["T_end"]), (o["oxide_c2"], s["c2"]["oxide_um"][-1]),
                     (o["T_zs10"], s["zs10"]["T_end"]), (o["T_zs15"], s["zs15"]["T_end"])]
            for k, tt in enumerate((201, 402, 600)):
                pairs.append((o["oxide_c1"][k], float(np.interp(tt, s["c1"]["t"], s["c1"]["oxide_um"]))))
        else:
            pairs = []
            for name in run_v9.CASES:
                s = v9[f"{cat}_{name}"]
                pairs += [(o[name]["T_low"], s["T_low"]), (o[name]["T_ref"], s["T_ref"])]
                if not math.isnan(s["t_glass"]):
                    pairs.append((o[name]["t_glass"], s["t_glass"]))
                elif not math.isnan(o[name]["t_glass"]):
                    bad.append(f"V9 {cat} {name}: glass lost in the gate run but not in the check")
        for a, b in pairs:
            if not close(a, b):
                bad.append(f"{r['case']} {cat} {r['i']}: {a!r} vs stored {b!r}")
        print(f"gate {r['case']:5s} {cat} {r['i']:2d}: {len(pairs)} values, cpu {r['cpu_s']} s", flush=True)
    r7, q_nc, f_nc = v7_ratio()
    print(f"V7 heating ratio r = CFD 5.0 MW/m2 / model non-catalytic cold-wall {q_nc / 1e6:.3f} MW/m2 = {r7:.4f} (f_nc {f_nc:.4f})")
    (HERE / "vvu_gate.json").write_text(json.dumps(dict(passed=not bad, problems=bad, r7=r7, q_nc_model=q_nc, f_nc=f_nc,
                                                         cpu=[(r["case"], r["bound"], r["i"], r["cpu_s"]) for r in res]), indent=1))
    print("GATE PASSED" if not bad else "GATE FAILED:\n  " + "\n  ".join(bad[:20]))
    return 0 if not bad else 1


def all_jobs(n7, r7):
    jobs = []
    for case in CASES:
        d = designs(case, n7 if case == "V7" else 100, r7=r7)
        for cat in ("fc", "nc"):
            for sset in ("joint", "bc", "input"):
                jobs += [(case, cat, sset, i, x) for i, x in enumerate(d[sset])]
            if case != "V7":  # deterministic calorimeter cases (plan Section 5.2): gamma_Cu = 0.01 and 0.067, rest baseline
                for g in (0.01, 0.067):
                    jobs.append((case, cat, "cal", int(round(g * 1000)), {**d["base"][0], "gamma_cu": g}))
    return jobs


def main():
    what = sys.argv[1] if len(sys.argv) > 1 else "gate"
    if what == "gate":
        return gate()
    g = json.loads((HERE / "vvu_gate.json").read_text())
    if not g["passed"]:
        raise SystemExit("the reproducibility gate has not passed: stop (plan Section 7)")
    r7 = g["r7"]
    if what == "time":
        jobs = []
        for case in CASES:
            d = designs(case, 100, r7=r7)
            jobs += [(case, "fc", "joint", i, d["joint"][i]) for i in range(3)]
        with ProcessPoolExecutor(max_workers=7) as ex:
            res = list(ex.map(task, jobs))
        cost = {}
        for r in res:
            cost.setdefault(r["case"], []).append(r["cpu_s"])
            if r["error"]:
                print("error:", r["case"], r["error"])
        per = {k: float(np.mean(v)) for k, v in cost.items()}
        n_runs = {"V7": 2 * 200, "V8_1": 2 * 204, "V8_2": 2 * 204, "V9": 2 * 204}
        total = sum(per[k] * n_runs[k] for k in per) / 7 / 3600
        print("mean cpu per sample (s):", {k: round(v, 1) for k, v in per.items()})
        print(f"estimated wall time with 7 workers and N7 = 100: {total:.1f} h ->",
              "reduce V7 joint N to 50 (plan Section 8)" if total > 8 else "keep N = 100")
        (HERE / "vvu_timing.json").write_text(json.dumps(dict(per_sample_cpu_s=per, est_hours_n100=total), indent=1))
        return 0
    n7 = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    jobs = all_jobs(n7, r7)
    done = set()
    if OUT.exists():
        for line in OUT.read_text(encoding="utf-8").splitlines():
            r = json.loads(line)
            done.add((r["case"], r["bound"], r["set"], r["i"]))
    todo = [j for j in jobs if (j[0], j[1], j[2], j[3]) not in done]
    # long 2D runs first so that the pool stays busy at the end
    todo.sort(key=lambda j: 0 if j[0] == "V7" else 1)
    print(f"{len(jobs)} jobs, {len(done)} done, {len(todo)} to run", flush=True)
    run_all(task, todo, OUT)
    return 0


# Plan Section 6 and Amendment 2: a run that does not return is a numerical failure. In the committed run four runs were
# stopped by hand after 56 min; for a reproduction the runner stops any run after RUN_LIMIT_S of wall time (more than nine
# times the slowest completed run, 182 s) and records it the same way, with its inputs and no outputs.
RUN_LIMIT_S = 1800.0


def run_all(fn, todo, out, limit_s=RUN_LIMIT_S, workers=7, poll_s=10.0):
    """Run fn over todo in a process pool, appending one JSON line per run to out; stop runs that exceed limit_s."""
    t0, n_done, n_all = time.time(), 0, len(todo)
    while todo:
        ex = ProcessPoolExecutor(max_workers=workers)
        futs = {ex.submit(fn, j): j for j in todo}
        pending, started, stopped = set(futs), {}, False
        with open(out, "a", encoding="utf-8") as f:
            while pending:
                ready, pending = wait(pending, timeout=poll_s, return_when=FIRST_COMPLETED)
                for fu in ready:
                    r = fu.result()
                    f.write(json.dumps(r, default=float) + "\n")
                    f.flush()
                    n_done += 1
                    if n_done % 50 == 0 or n_done == n_all:
                        print(f"{n_done}/{n_all} runs, {time.time() - t0:.0f} s", flush=True)
                now = time.time()
                for fu in pending:
                    if fu.running():
                        started.setdefault(fu, now)
                late = [fu for fu in pending if fu in started and now - started[fu] > limit_s]
                if late:
                    for fu in late:
                        case, cat, sset, i, x = futs[fu]
                        f.write(json.dumps(dict(case=case, bound=cat, set=sset, i=i, x=x, out=None, cpu_s=None,
                                                error=f"no result after {limit_s:.0f} s of wall time; stopped and left out "
                                                      "(plan Section 6, Amendment 2)"), default=float) + "\n")
                        f.flush()
                        n_done += 1
                        print(f"stopped after {limit_s:.0f} s: {case} {cat} {sset} {i}", flush=True)
                    for p in list(getattr(ex, "_processes", {}).values()):
                        p.kill()
                    ex.shutdown(wait=False, cancel_futures=True)
                    stopped = True
                    break
        if not stopped:
            ex.shutdown()
        written = {(r["case"], r["bound"], r["set"], r["i"])
                   for r in map(json.loads, Path(out).read_text(encoding="utf-8").splitlines())}
        todo = [j for j in todo if (j[0], j[1], j[2], j[3]) not in written]


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    sys.exit(main())
