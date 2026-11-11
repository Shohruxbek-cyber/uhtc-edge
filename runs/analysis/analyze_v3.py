"""Version 3: boundary of every variant column against the baseline 2D boundary (runs/analysis/boundaries_2d.csv, fully catalytic).

For each variant file runs/envelope_2d/<stage>_<variant>.jsonl and each (R_n, altitude) column: highest surviving and lowest
failing Mach number, the shift of the boundary (highest surviving Mach) in grid steps of 0.25 against the baseline, the first
failure mode at the lowest failing Mach against the baseline's, and the failure (glass-loss) time at the baseline's lowest failing
Mach where the variant ran that point. Reference: baseline columns, or for flow_alpha_0.005 the flow-off transport_alpha_0.005 stage.
Output: runs/analysis/v3_boundaries.csv and runs/analysis/v3_boundaries.md.
Usage: python runs/analysis/analyze_v3.py
"""
import csv
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
E2D = ROOT / "runs/envelope_2d"
STAGES = {"flow": ["flow", "flow_alpha_0.005", "flow_F320", "flow_tau_x0.5", "flow_tau_x2", "flow_no_pressure"],
          "flankint": ["flank_integral"], "geometry": ["theta_10", "flank_15mm", "flank_60mm"]}
REFS = {"flow_alpha_0.005": "transport_alpha_0.005"}


def load(f):
    return [json.loads(l) for l in open(f, encoding="utf-8") if l.strip()] if f.exists() else []


def columns(rows):
    out = {}
    for r in rows:
        if r.get("error"):
            continue
        k = (round(r["R_n"] * 1e3, 6), round(r["alt"] / 1e3, 6))
        c = out.setdefault(k, dict(lo=None, hi=None, rows={}))
        c["rows"][r["mach"]] = r
        if r["mode"] == "none":
            c["lo"] = r["mach"] if c["lo"] is None else max(c["lo"], r["mach"])
        else:
            c["hi"] = r["mach"] if c["hi"] is None else min(c["hi"], r["mach"])
    return out


def baseline():
    b = {}
    for r in csv.DictReader(open(ROOT / "runs/analysis/boundaries_2d.csv", encoding="utf-8")):
        if r["catalysis"] == "fc":
            b[(float(r["R_n_mm"]), float(r["alt_km"]))] = dict(lo=float(r["M_s_2d"]), hi=float(r["M_f_2d"]), mode=r["mode_2d"],
                                                               t_fail=float(r["t_fail_2d"]))
    return b


def main():
    base = baseline()
    out = []
    for stage, variants in STAGES.items():
        for v in variants:
            rows = load(E2D / f"{stage}_{v}.jsonl")
            if not rows:
                continue
            cols = columns(rows)
            ref = None
            if v in REFS:
                ref = {k: dict(lo=c["lo"], hi=c["hi"], mode=c["rows"][c["hi"]]["mode"] if c["hi"] in c["rows"] else "",
                               t_fail=c["rows"][c["hi"]]["t_fail"] if c["hi"] in c["rows"] else math.nan)
                       for k, c in columns(load(E2D / f"transport_{REFS[v][len('transport_'):]}.jsonl")).items()}
            for k, c in sorted(cols.items()):
                b = (ref or base).get(k)
                if b is None:
                    continue
                at_old = c["rows"].get(b["hi"])
                mode_v = c["rows"][c["hi"]]["mode"] if c["hi"] in c["rows"] else ""
                out.append(dict(stage=stage, variant=v, R_mm=k[0], alt_km=k[1], M_s_ref=b["lo"], M_f_ref=b["hi"], M_s=c["lo"], M_f=c["hi"],
                                shift_steps=(None if c["lo"] is None or b["lo"] is None else round((c["lo"] - b["lo"]) / 0.25, 2)),
                                mode_ref=b["mode"], mode=mode_v,
                                t_fail_ref_at_old=b["t_fail"], t_fail_at_old=(at_old["t_fail"] if at_old else None),
                                mode_at_old=(at_old["mode"] if at_old else None), ref="flow off, alpha 0.005" if ref else "baseline"))
    with open(ROOT / "runs/analysis/v3_boundaries.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    L = ["# Version-3 boundary columns against the reference (fully catalytic)", "",
         "| Stage | Variant | R (mm) | Alt (km) | Boundary ref -> variant (highest surviving Mach) | Shift (0.25 steps) | First failure ref -> variant | "
         "Failure time at the old lowest failing Mach, ref -> variant (s) |", "|---|---|---|---|---|---|---|---|"]
    for r in out:
        tf = lambda x: "none" if x is None or x != x else f"{x:.0f}"  # noqa: E731
        L.append(f"| {r['stage']} | {r['variant']} | {r['R_mm']:g} | {r['alt_km']:g} | {r['M_s_ref']} -> {r['M_s']} | {r['shift_steps']} | "
                 f"{r['mode_ref']} -> {r['mode']} | {tf(r['t_fail_ref_at_old'])} -> {tf(r['t_fail_at_old'])} |")
    (ROOT / "runs/analysis/v3_boundaries.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    from collections import Counter
    for (st, v), grp in sorted({(r["stage"], r["variant"]): None for r in out}.items()):
        rs = [r for r in out if r["stage"] == st and r["variant"] == v]
        print(st, v, len(rs), "columns; shifts", Counter(r["shift_steps"] for r in rs), "; modes", Counter(f"{r['mode_ref']}->{r['mode']}" for r in rs))


if __name__ == "__main__":
    main()
