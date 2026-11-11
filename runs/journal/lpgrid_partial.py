"""Journal revision, ASK 5: the local-pressure grid was stopped by the author on 11 Nov 2026 at 12:51 after 365 of the 630 runs
at 1 mm. Point-wise comparison of those runs with the baseline grid at the same points. Output: runs/journal/lpgrid_partial_R1.txt"""
import math
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "runs" / "analysis"))
import analyze_2d as A  # noqa: E402

E2 = ROOT / "runs" / "envelope_2d"


def key(r):
    return (round(r["R_n"] * 1e3, 6), r["catalysis"], round(r["alt"] / 1e3, 6), round(r["mach"], 6))


base = {key(r): r for r in A.load(E2 / "grid.jsonl") if A.ok(r)}
lp = [r for r in A.load(E2 / "lpgrid_R1.jsonl")]
bad = [r for r in lp if not A.ok(r)]
lp = [r for r in lp if A.ok(r)]
fail = lambda r: r["mode"] != "none"  # noqa: E731
n = Counter()
dt, modes, cols = [], Counter(), set()
for r in lp:
    b = base.get(key(r))
    if b is None:
        n["no baseline point"] += 1
        continue
    cols.add(key(r)[:3])
    n[(r["catalysis"], "same outcome" if fail(r) == fail(b) else f"changed {('fail' if fail(b) else 'pass')}->{('fail' if fail(r) else 'pass')}")] += 1
    if fail(r) and fail(b):
        modes[(b["mode"], r["mode"])] += 1
        if b["t_fail"] and not math.isnan(b["t_fail"]):
            dt.append((r["t_fail"] - b["t_fail"]) / b["t_fail"] * 100)
print("runs read:", len(lp) + len(bad), "errors:", len(bad))
print("columns touched:", len(cols), sorted({(c[1], c[2]) for c in cols})[:3], "...")
for k, v in sorted(n.items(), key=str):
    print(" ", k, v)
print("first-failure mode baseline -> local pressure:", dict(modes))
if dt:
    dt.sort()
    print(f"failure time change, failing in both: n={len(dt)}, min {dt[0]:+.2f} %, median {dt[len(dt)//2]:+.2f} %, max {dt[-1]:+.2f} %")
# columns complete in the partial grid: boundary at grid resolution (highest passing Mach below the lowest failing)
by = {}
for r in lp:
    by.setdefault(key(r)[:3], []).append(r)
same = diff = 0
for c, rs in sorted(by.items()):
    if len(rs) < 15:
        continue
    bs = [base[key(r)] for r in rs if key(r) in base]
    def bnd(rows):
        lo = min((x["mach"] for x in rows if fail(x)), default=99)
        return max((x["mach"] for x in rows if not fail(x) and x["mach"] < lo), default=None), lo
    if bnd(rs) == bnd(bs):
        same += 1
    else:
        diff += 1
        print("  column differs:", c, "baseline", bnd(bs), "local", bnd(rs))
print(f"complete columns (all 15 Mach numbers): {same + diff}, boundary at 0.5-Mach grid resolution unchanged in {same}")
# lowest failing grid point of each complete column (the paper's 1-4 % refers to the refined lowest failing conditions of the
# 15 band columns, a different set)
lf = []
for c, rs in by.items():
    f = [r for r in rs if fail(r)]
    if len(rs) < 15 or not f:
        continue
    m = min(f, key=lambda r: r["mach"])
    b = base[key(m)]
    lf.append((m["t_fail"] - b["t_fail"]) / b["t_fail"] * 100)
print(f"lowest failing grid point of each complete column: n={len(lf)}, failure time {min(lf):+.1f} to {max(lf):+.1f} %")
