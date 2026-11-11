"""VVU analysis (plan validation/VVU/VVU_plan.md, Sections 2-5): comparison error, uncertainty groups and labels per quantity.

Reads vvu_runs.jsonl (run_vvu.py) and the stored results of the checks (the baseline S, identical to sample 0 by the gate).
Writes vvu_samples.csv (every run: inputs and outputs), vvu_results.json and report.md.
"""
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
C2K = 273.15
R3 = math.sqrt(3.0)
# plan 5.1: bound (fraction), uniform half-width
NUM = {"1d": dict(T=0.0008, time=0.0044, thick=0.0003), "2d": dict(T=0.005, time=0.02, thick=0.02)}
# plan 5.2: every check has test conditions whose source states no uncertainty (V7 enthalpy and Pitot pressure; V8-1 flux
# and pressure; V8-2 flux, enthalpy and pressure; V9 flux, enthalpy and pressures), so u_BC is a lower bound for every check
BC_LOWER = {"V7": True, "V8_1": True, "V8_2": True, "V9": True}
# V9 glass-loss time censored at the end of the run (run_v9.T_END = 480 s; switches at 357 and 359 s): time left after the switch
CENS_AT = {"sample 1.5: glass loss after the switch": 480.0 - 357.0, "sample 1.6: glass loss after the switch": 480.0 - 359.0}


def runs():
    out = defaultdict(list)
    for line in (HERE / "vvu_runs.jsonl").read_text(encoding="utf-8").splitlines():
        r = json.loads(line)
        out[(r["case"], r["bound"], r["set"])].append(r)
    return out


def qty_defs():
    """(case, name, kind, D, u_D (K or log), u_D lower bound?, extract(out) -> value or None for censored)."""
    q = []
    T7, tol = [1573.0, 1723.0, 1843.0, 1923.0, 2053.0], None
    for k, D in enumerate(T7):
        uD = math.sqrt((0.01 * D / R3) ** 2 + (0.02 * D / R3) ** 2 + ((10.0 / R3) ** 2 if k == 2 else 0.0))
        q.append(("V7", f"plateau temperature, step {k + 1}", "T", D, uD, False, lambda o, k=k: o["Tend"][k]))
    q.append(("V7", "oxide at the stagnation point", "thick", 150.0, 0.0, True, lambda o: o["oxide_um"]))
    q.append(("V7", "SiC-depleted layer (none observed)", "outcome", None, None, None, lambda o: o["dep_um"] < 1.0))
    q.append(("V8_1", "peak temperature (30 MJ/kg nominal)", "T", 1630 + C2K, 0.0, True, lambda o: o["T_max"]))
    q.append(("V8_1", "oxide at 300 s", "thick", 110.0, 0.0, True, lambda o: o["oxide_um"]))
    q.append(("V8_1", "SiC-depleted layer at 300 s", "thick", 40.0, 0.0, True, lambda o: o["dep_um"]))
    q.append(("V8_2", "steady temperature, condition 1", "T", 1960 + C2K, 0.0, True, lambda o: o["T1"]))
    q.append(("V8_2", "steady temperature, condition 2", "T", 1900 + C2K, 0.0, True, lambda o: o["T2"]))
    for k, (t, D) in enumerate(((201, 65.0), (402, 99.0), (600, 125.0))):
        q.append(("V8_2", f"reacted layer (oxide + depleted), condition 1, {t} s", "thick", D, 1.5 / R3 / D, True,
                  lambda o, k=k: o["reacted_c1"][k]))
        q.append(("V8_2", f"reacted layer (oxide only), condition 1, {t} s", "thick", D, 1.5 / R3 / D, True,
                  lambda o, k=k: o["oxide_c1"][k]))
    q.append(("V8_2", "oxide, condition 2, 600 s", "thick", 90.0, 0.0, True, lambda o: o["oxide_c2"]))
    q.append(("V8_2", "trend T(ZS10) > T(ZS15) > T(ZS30)", "outcome", None, None, None, lambda o: o["trend"]))
    meas = {"1.8": (1895.0, 2110.0, None), "1.5": (1825.0, 2210.0, 68.0), "1.6": (1860.0, 2250.0, 21.0)}
    uT9 = math.sqrt((10.0 / R3) ** 2 + (15.0 / R3) ** 2)
    for s, (tl, th, on) in meas.items():
        q.append(("V9", f"sample {s}: low-stage plateau", "T", tl, uT9, False, lambda o, s=s: o[s]["T_low"]))
        q.append(("V9", f"sample {s}: high-stage temperature at the reference time", "T", th, uT9, False,
                  lambda o, s=s: o[s]["T_ref"]))
        if on is not None:
            q.append(("V9", f"sample {s}: glass loss after the switch", "time", on, math.sqrt((3 / R3) ** 2 + (0.5 / R3) ** 2) / on,
                      False, lambda o, s=s: o[s]["after"]))
        want = on is not None
        q.append(("V9", f"sample {s}: jump outcome ({'glass loss' if want else 'no glass loss'} in the high stage)", "outcome",
                  None, None, None, lambda o, s=s, want=want: o[s]["lost_in_high"] == want))
    return q


def baseline_outputs():
    """Stored results of the checks, in the output format of run_vvu.py (identical to sample 0 by the gate)."""
    v7 = json.loads((ROOT / "validation/V7_arcjet/v7_results.json").read_text())
    v8 = json.loads((ROOT / "validation/V8_measured_flux/v8_results.json").read_text())
    v9 = json.loads((ROOT / "validation/V9_plasmatron/v9_results.json").read_text())
    out = {}
    for cat in ("fc", "nc"):
        r = v7[f"2d-{cat}"]
        out[("V7", cat)] = dict(Tend=r["Tend"], oxide_um=r["oxide_um"], dep_um=r["dep_um"])
        for H in (15, 30, 45):
            p = v8[f"V8_1_H{H}_{cat}"]["passive"]
            out[(f"V8_1@{H}", cat)] = dict(T_max=p["T_max"], oxide_um=p["oxide_um"][-1], dep_um=p["dep_um"][-1])
        out[("V8_1", cat)] = out[("V8_1@30", cat)]
        h = v8[f"V8_2_{cat}"]
        it = lambda r, k, t: float(np.interp(t, r["t"], r[k]))  # noqa: E731
        out[("V8_2", cat)] = dict(T1=h["c1"]["T_end"], T2=h["c2"]["T_end"], oxide_c2=h["c2"]["oxide_um"][-1],
                                  oxide_c1=[it(h["c1"], "oxide_um", t) for t in (201, 402, 600)],
                                  reacted_c1=[it(h["c1"], "oxide_um", t) + it(h["c1"], "dep_um", t) for t in (201, 402, 600)],
                                  trend=bool(h["zs10"]["T_end"] > h["zs15"]["T_end"] > h["c1"]["T_end"]))
        o9 = {}
        for s, ts in (("1.8", 358.0), ("1.5", 357.0), ("1.6", 359.0)):
            r = v9[f"{cat}_{s}"]
            tg = r["t_glass"]
            o9[s] = dict(T_low=r["T_low"], T_ref=r["T_ref"], after=None if tg != tg else tg - ts, lost_in_high=bool(tg == tg and tg > ts))
        out[("V9", cat)] = o9
    return out


def sd(vals, log):
    v = [x for x in vals if x is not None and not (isinstance(x, float) and math.isnan(x))]
    if len(v) < 2:
        return float("nan"), len(vals) - len(v)
    a = np.log(v) if log else np.asarray(v, float)
    return float(np.std(a, ddof=1)), len(vals) - len(v)


def main():
    R, base, Q = runs(), baseline_outputs(), qty_defs()
    results, lines = [], []
    for case, name, kind, D, uD, lower, ex in Q:
        dim = "2d" if case == "V7" else "1d"
        for cat in ("fc", "nc"):
            S = ex(base[(case, cat)])
            joint = [ex(r["out"]) for r in R[(case, cat, "joint")] if r["out"] is not None]
            rec = dict(check=case, quantity=name, bound=cat, kind=kind)
            if kind == "outcome":
                share = float(np.mean([bool(v) for v in joint])) if joint else float("nan")
                rec.update(S=bool(S), share_reproducing=share, n=len(joint),
                           label="outcome not reproduced within input uncertainty" if share < 0.05 else
                           f"outcome reproduced in {100 * share:.0f} % of the joint sample")
                results.append(rec)
                continue
            log = kind in ("thick", "time")
            u_num = NUM[dim][kind] / R3
            bc = [ex(r["out"]) for r in R[(case, cat, "bc")] if r["out"] is not None]
            inp = [ex(r["out"]) for r in R[(case, cat, "input")] if r["out"] is not None]
            u_bc, cens_bc = sd(bc, log)
            u_in, cens_in = sd(inp, log)
            u_jt, cens_jt = sd(joint, log)
            if kind == "T":
                u_num_abs = u_num * S
                E = S - D
            else:
                u_num_abs = u_num
                E = math.log(S / D) if S is not None and S > 0 else float("nan")
            comps = dict(num=u_num_abs, bc=u_bc, input=u_in, D=uD)
            u_val = math.sqrt(sum((0.0 if v != v else v) ** 2 for v in comps.values()))
            U = 2.0 * u_val
            E_lower = None
            if E != E:
                # baseline censored (plan Section 4): no glass loss by the end of the run, so E > ln(time left / D)
                E_lower = math.log(CENS_AT[name] / D)
                n_c, n_g = cens_bc + cens_in, len(bc) + len(inp)
                if n_c == n_g:
                    label = (f"no glass loss by the end of the run in any sample (censored, E > {E_lower:.2f}); "
                             "see the outcome row")
                else:
                    label = (f"censored at the baseline (E > {E_lower:.2f}) and in {n_c} of {n_g} group samples; "
                             "the rule of plan Section 3 cannot be applied")
            else:
                label = ("model error indicated at the test uncertainty" if abs(E) > U else "not distinguishable at the test uncertainty")
                if (lower or BC_LOWER[case]) and abs(E) > U:
                    label += " (stated uncertainties only: U_val is a lower bound)"
            share = {k: (((0.0 if v != v else v) ** 2) / u_val ** 2 if u_val > 0 else float("nan")) for k, v in comps.items()}
            if E_lower is not None:
                # plan Amendment 3: the group samples that lose their glass are a biased subset, so no spread is reported
                u_bc = u_in = u_val = U = float("nan")
                share = {k: float("nan") for k in share}
            rec.update(S=S, D=D, E=E, E_lower=E_lower, unit="K" if kind == "T" else "ln", u_num=u_num_abs, u_BC=u_bc, u_input=u_in, u_D=uD,
                       u_D_lower_bound=lower, u_BC_lower_bound=BC_LOWER[case], u_val=u_val, U_val=U, label=label, share=share, u_joint=u_jt,
                       censored=dict(bc=cens_bc, input=cens_in, joint=cens_jt), n=dict(bc=len(bc), input=len(inp), joint=len(joint)))
            results.append(rec)
    # calorimeter cases (plan 5.2): change of S with gamma_Cu = 0.01 and 0.067
    cal = []
    for case, name, kind, D, uD, lower, ex in Q:
        if case == "V7" or kind == "outcome":
            continue
        for cat in ("fc", "nc"):
            S = ex(base[(case, cat)])
            row = dict(check=case, quantity=name, bound=cat, S=S)
            for r in R[(case, cat, "cal")]:
                v = ex(r["out"]) if r["out"] is not None else None
                row[f"gamma_{r['i'] / 1000:g}"] = v
            cal.append(row)
    # V8-1 at the three bracketing enthalpies (stored)
    v81 = {f"{H} MJ/kg, {cat}": base[(f"V8_1@{H}", cat)] for H in (15, 30, 45) for cat in ("fc", "nc")}
    gate = json.loads((HERE / "vvu_gate.json").read_text())
    errors = sum(1 for k, v in R.items() for r in v if r["error"])
    (HERE / "vvu_results.json").write_text(json.dumps(dict(results=results, calorimeter=cal, v8_1_bracket=v81, r7=gate["r7"],
                                                           runs={"/".join(k): len(v) for k, v in R.items()}, failed_runs=errors),
                                                      indent=1, default=float))
    # samples table
    keys = sorted({k for v in R.values() for r in v for k in r["x"]})
    with open(HERE / "vvu_samples.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["case", "bound", "set", "i", "error", "cpu_s"] + keys + ["outputs"])
        for (case, cat, sset), v in sorted(R.items()):
            for r in sorted(v, key=lambda r: r["i"]):
                w.writerow([case, cat, sset, r["i"], r["error"] or "", r["cpu_s"]] + [r["x"].get(k, "") for k in keys]
                           + [json.dumps(r["out"], default=float)])
    # report
    L = ["# VVU: validation uncertainty of V7–V9 (post hoc; the preset verdicts of Table 2 are unchanged)", "",
         "Plan: `VVU_plan.md` (committed before any run, amendment 1 before any sample). Method: ASME V&V 20-2009 (R2021). "
         "E = S − D in K for temperatures and E = ln(S/D) for thicknesses and times; U_val = 2 u_val. "
         f"V7 heating ratio r = {gate['r7']:.3f}. Failed runs: {errors}.", "",
         "| Check | Quantity | Bound | S | D | E | u_num | u_BC | u_input | u_D | U_val | Label | Share of U_val² (num/BC/input/D) |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    f2 = lambda x: "—" if x is None or x != x else (f"{x:.3f}" if abs(x) < 10 else f"{x:.0f}")  # noqa: E731
    for r in results:
        if r["kind"] == "outcome":
            L.append(f"| {r['check']} | {r['quantity']} | {r['bound']} | {'yes' if r['S'] else 'no'} | measured | | | | | | | "
                     f"{r['label']} (n = {r['n']}) | |")
            continue
        sh = "—" if r["E_lower"] is not None else "/".join(f"{100 * r['share'][k]:.0f}" for k in ("num", "bc", "input", "D"))
        L.append(f"| {r['check']} | {r['quantity']} | {r['bound']} | {f2(r['S'])} | {f2(r['D'])} | {f2(r['E']) if r['E_lower'] is None else '> ' + f2(r['E_lower'])} | {f2(r['u_num'])} | "
                 f"{f2(r['u_BC'])} | {f2(r['u_input'])} | {f2(r['u_D'])}{'*' if r['u_D_lower_bound'] else ''} | {f2(r['U_val'])} | "
                 f"{r['label']} | {sh} |")
    L += ["", "\\* u_D holds only what the source states. For every check u_BC holds only what the sources state (plan 5.2), so U_val is a lower bound throughout.", "",
          "## Calorimeter catalysis (plan 5.2): S with a copper efficiency of 0.01 and 0.067, other inputs at baseline", "",
          "| Check | Quantity | Bound | S (fully catalytic calorimeter) | γ_Cu 0.01 | γ_Cu 0.067 |", "|---|---|---|---|---|---|"]
    for c in cal:
        L.append(f"| {c['check']} | {c['quantity']} | {c['bound']} | {f2(c['S'])} | {f2(c.get('gamma_0.01'))} | {f2(c.get('gamma_0.067'))} |")
    (HERE / "report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    sys.exit(main())
