"""Final review (8-9 Nov 2026): failure time of each 2D variant over the baseline failure time one grid step beyond the
baseline boundary (lowest failing baseline Mach number, fully catalytic band columns), per radius, and highest passing
baseline conditions that fail with the variant.  Usage: python runs/final_review/ratios_2d.py > runs/final_review/ratios_2d.md"""
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "runs" / "envelope_2d"
A = ROOT / "runs" / "analysis"


def load(p):
    return [json.loads(l) for l in open(p, encoding="utf-8") if l.strip()]


base = {}
for f in ("grid.jsonl", "refine.jsonl"):
    for r in load(D / f):
        if r.get("catalysis") == "fc" and not r.get("error"):
            base[(round(r["R_n"] * 1e3, 3), round(r["alt"] / 1e3, 3), r["mach"])] = r
import csv
bnd = {(float(r["R_n_mm"]), float(r["alt_km"])): r for r in csv.DictReader(open(A / "boundaries_2d.csv", encoding="utf-8"))
       if r["catalysis"] == "fc"}
VARIANTS = ["kinetics_kp_x0.5", "kinetics_kp_x3", "kinetics_kp_E250", "kinetics_kp_E354", "kinetics_incubation",
            "transport_mtf_x1.25", "transport_mtf_x2", "transport_alpha_0.1", "transport_alpha_0.005",
            "transport_kZrO2_1.0", "transport_kZrO2_2.5", "active_pa_balat", "active_pa_fe"]
for var in sys.argv[1:] or VARIANTS:
    rows = load(D / f"{var}.jsonl")
    ratios = defaultdict(list)
    lost = []
    for r in rows:
        key = (round(r["R_n"] * 1e3, 3), round(r["alt"] / 1e3, 3))
        b = bnd.get(key)
        if not b or b["M_f_2d"] in ("", "None"):
            continue
        Mf, Ms = float(b["M_f_2d"]), float(b["M_s_2d"])
        if abs(r["mach"] - Mf) < 1e-9:
            t0 = base[key + (Mf,)]["t_fail"]
            t1 = r["t_fail"] if r["t_fail"] == r["t_fail"] else math.inf
            ratios[key[0]].append((t1 / t0, key[1], r["mode"]))
        if abs(r["mach"] - Ms) < 1e-9 and r["t_fail"] == r["t_fail"]:
            lost.append((key, Ms, round(r["t_fail"]), r["mode"]))
    print(var)
    for R, v in sorted(ratios.items()):
        fin = [x for x, _, _ in v if x != math.inf]
        rng = f"{min(fin):.2f}..{max(fin):.2f}" if fin else "-"
        print(f"   R {R:g} mm: ratio {rng} (n={len(v)}, survive {sum(1 for x, _, _ in v if x == math.inf)}), "
              f"detail {[(round(x, 2), z, m) for x, z, m in v]}")
    print("   highest passing baseline condition now fails:", lost)
