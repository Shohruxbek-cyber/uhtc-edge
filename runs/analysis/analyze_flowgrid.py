"""Version 3, Task 1: fully catalytic 2D maps with glass flow on against the baseline maps (flow off), all 21 altitudes.

Inputs : runs/envelope_2d/{grid,refine}.jsonl (baseline, fully catalytic columns)
         runs/envelope_2d/flowgrid_R{1,3,10}.jsonl + flowrefine_R{1,3,10}.jsonl (glass flow on, otherwise identical cases)
Boundary of a column as in analyze_2d.py (same function): highest surviving Mach below the lowest failing Mach.
Point-wise: for every grid point run in both, the change of the outcome (survive / first-failure mode) and of the
failure time.
Outputs: runs/analysis/flowgrid_boundaries.csv   one row per (R_n, altitude)
         runs/analysis/flowgrid_summary.json     counts used in the paper
         runs/analysis/flowgrid_summary.md       readable tables
Usage: python runs/analysis/analyze_flowgrid.py
"""
import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
E2 = ROOT / "runs" / "envelope_2d"
OUT = ROOT / "runs" / "analysis"
sys.path.insert(0, str(OUT))
import analyze_2d as A  # noqa: E402

RADII_MM = (1, 3, 10)


def key(r):
    return (round(r["R_n"] * 1e3, 6), round(r["alt"] / 1e3, 6))


def main():
    base = [r for r in A.load(E2 / "grid.jsonl") + A.load(E2 / "refine.jsonl") if r["catalysis"] == "fc"]
    base_grid = {(key(r), r["mach"]): r for r in A.load(E2 / "grid.jsonl") if r["catalysis"] == "fc" and A.ok(r)}
    bcols, fcols, fgrid, errors = {}, {}, {}, []
    for r in base:
        bcols.setdefault(key(r), []).append(r)
    for R in RADII_MM:
        g = A.load(E2 / f"flowgrid_R{R}.jsonl")
        f = A.load(E2 / f"flowrefine_R{R}.jsonl")
        errors += [r["case_id"] for r in g + f if not A.ok(r)]
        for r in g + f:
            fcols.setdefault(key(r), []).append(r)
        for r in g:
            if A.ok(r):
                fgrid[(key(r), r["mach"])] = r
    rows = []
    for k in sorted(fcols):
        b, v = A.boundary(bcols[k]), A.boundary(fcols[k])
        shift = None
        if b["M_s"] is not None and v["M_s"] is not None:
            shift = round((v["M_s"] - b["M_s"]) / 0.25, 2)
        at_old = next((r for r in fcols[k] if A.ok(r) and b["M_f"] is not None and r["mach"] == b["M_f"]), None)
        rows.append(dict(R_n_mm=k[0], alt_km=k[1], M_s_base=b["M_s"], M_f_base=b["M_f"], mode_base=b["mode_f"],
                         t_fail_base=b["t_fail_f"], M_s_flow=v["M_s"], M_f_flow=v["M_f"], mode_flow=v["mode_f"],
                         t_fail_flow=v["t_fail_f"], shift_steps=shift,
                         t_fail_flow_at_old=at_old["t_fail"] if at_old else None,
                         mode_flow_at_old=at_old["mode"] if at_old else None,
                         dt_at_old_pct=(100 * (at_old["t_fail"] / b["t_fail_f"] - 1)
                                        if at_old and b["t_fail_f"] and at_old["mode"] != "none" else None),
                         nonmono_base=b["nonmono"], nonmono_flow=v["nonmono"]))
    with open(OUT / "flowgrid_boundaries.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    # point-wise comparison on the plain grid
    pts = []
    for (k, M), rf in fgrid.items():
        rb = base_grid.get((k, M))
        if rb is None:
            continue
        pts.append(dict(R=k[0], alt=k[1], M=M, mode_b=rb["mode"], mode_f=rf["mode"], t_b=rb["t_fail"], t_f=rf["t_fail"],
                        Tmax_b=rb["T_s_max"], Tmax_f=rf["T_s_max"], glass_b=rb.get("glass_um_end"), glass_f=rf.get("glass_um_end"),
                        shed=rf.get("flow_shed")))
    summ = dict(n_flow_points=sum(len(v) for v in fcols.values()), n_paired_grid=len(pts), errors=errors, radii={})
    for R in RADII_MM:
        rr = [r for r in rows if r["R_n_mm"] == R]
        pp = [p for p in pts if p["R"] == R]
        both_fail = [p for p in pp if p["mode_b"] != "none" and p["mode_f"] != "none"]
        dts = [100 * (p["t_f"] / p["t_b"] - 1) for p in both_fail if p["t_b"]]
        glass_both = [p for p in both_fail if p["mode_b"] == "glass" and p["mode_f"] == "glass"]
        dtg = [100 * (p["t_f"] / p["t_b"] - 1) for p in glass_both if p["t_b"]]
        summ["radii"][R] = dict(
            columns=len(rr), shifts=dict(Counter(str(r["shift_steps"]) for r in rr)),
            mode_change=dict(Counter(f"{r['mode_base']}->{r['mode_flow']}" for r in rr)),
            shifted=[dict(alt=r["alt_km"], base=r["M_s_base"], flow=r["M_s_flow"], mode=f"{r['mode_base']}->{r['mode_flow']}")
                     for r in rr if r["shift_steps"] not in (0, 0.0, None)],
            dt_at_old_pct_range=[min((r["dt_at_old_pct"] for r in rr if r["dt_at_old_pct"] is not None), default=None),
                                 max((r["dt_at_old_pct"] for r in rr if r["dt_at_old_pct"] is not None), default=None)],
            paired_points=len(pp),
            outcome_changes=dict(Counter(f"{p['mode_b']}->{p['mode_f']}" for p in pp if p["mode_b"] != p["mode_f"])),
            dt_fail_pct_range_both_fail=[min(dts, default=None), max(dts, default=None)],
            dt_fail_pct_median_both_fail=(sorted(dts)[len(dts) // 2] if dts else None),
            dt_glass_pct_range=[min(dtg, default=None), max(dtg, default=None)], n_glass_both=len(glass_both),
            nonmono_flow=sum(r["nonmono_flow"] for r in rr), nonmono_base=sum(r["nonmono_base"] for r in rr))
        gl = [r["dt_at_old_pct"] for r in rr if r["mode_base"] == "glass" and r["mode_flow_at_old"] == "glass"]
        te = [r["dt_at_old_pct"] for r in rr if r["mode_base"] == "temperature" and r["mode_flow_at_old"] == "temperature"]
        summ["radii"][R].update(dt_at_old_glass_pct_range=[min(gl, default=None), max(gl, default=None)], n_glass_cols=len(gl),
                                dt_at_old_temperature_pct_range=[min(te, default=None), max(te, default=None)], n_temperature_cols=len(te))
    # the flow band stage (runs/envelope_2d/flow_flow.jsonl) ran some of the same points: they must be identical
    same, shared = 0, 0
    fmap = {}
    for R in RADII_MM:
        for r in A.load(E2 / f"flowgrid_R{R}.jsonl") + A.load(E2 / f"flowrefine_R{R}.jsonl"):
            fmap[(key(r), r["mach"])] = r
    eq = lambda a, b: a == b or (isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b))  # noqa: E731
    for r in A.load(E2 / "flow_flow.jsonl"):
        g = fmap.get((key(r), r["mach"]))
        if g is not None:
            shared += 1
            same += all(eq(r[x], g.get(x)) for x in r if x not in ("case_id", "cpu_s"))
    summ["band_stage_crosscheck"] = dict(shared=shared, identical=same)
    (OUT / "flowgrid_summary.json").write_text(json.dumps(summ, indent=1), encoding="utf-8")
    fm = lambda x: "none" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:g}" if isinstance(x, (int, float)) else str(x)  # noqa: E731
    L = ["# Fully catalytic 2D maps with glass flow on against the baseline (flow off)", "",
         f"Flow points run: {summ['n_flow_points']}; paired grid points: {summ['n_paired_grid']}; errors: {len(errors)}.", "",
         "| R (mm) | Alt (km) | Boundary base -> flow | Shift (0.25 steps) | First failure base -> flow | t_fail base -> flow at old M_f (s) |",
         "|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {fm(r['R_n_mm'])} | {fm(r['alt_km'])} | {fm(r['M_s_base'])} -> {fm(r['M_s_flow'])} | {fm(r['shift_steps'])} | "
                 f"{r['mode_base']} -> {r['mode_flow']} | {fm(r['t_fail_base'] and round(r['t_fail_base']))} -> "
                 f"{fm(r['t_fail_flow_at_old'] and round(r['t_fail_flow_at_old']))} ({r['mode_flow_at_old']}) |")
    L += ["", "## Per radius", "", "```", json.dumps(summ["radii"], indent=1), "```"]
    (OUT / "flowgrid_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
