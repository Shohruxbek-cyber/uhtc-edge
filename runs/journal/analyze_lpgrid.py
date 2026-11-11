"""Journal revision, ASK 5: the full 2D grid with the local wall pressure at every surface node against the baseline maps.

Inputs : runs/envelope_2d/{grid,refine}.jsonl (baseline, both catalysis bounds)
         runs/envelope_2d/lpgrid_R{1,3,10}.jsonl + lprefine_R{1,3,10}.jsonl (local pressure, otherwise identical cases)
Boundary of a column as in runs/analysis/analyze_2d.py (same function): highest passing Mach below the lowest failing Mach.
Point-wise: for every grid point run in both, the change of the outcome and of the failure time.
Outputs: runs/journal/lpgrid_boundaries.csv   one row per (R_n, catalysis bound, altitude)
         runs/journal/lpgrid_summary.json     counts used in the paper
Usage: python runs/journal/analyze_lpgrid.py   (works on the radii finished so far)
"""
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
E2 = ROOT / "runs" / "envelope_2d"
OUT = ROOT / "runs" / "journal"
sys.path.insert(0, str(ROOT / "runs" / "analysis"))
import analyze_2d as A  # noqa: E402

RADII_MM = (1, 3, 10)


def key(r):
    return (round(r["R_n"] * 1e3, 6), r["catalysis"], round(r["alt"] / 1e3, 6))


def main():
    bcols, lcols, bgrid, lgrid, errors, done = {}, {}, {}, {}, [], []
    for r in A.load(E2 / "grid.jsonl") + A.load(E2 / "refine.jsonl"):
        bcols.setdefault(key(r), []).append(r)
    for r in A.load(E2 / "grid.jsonl"):
        if A.ok(r):
            bgrid[(key(r), r["mach"])] = r
    for R in RADII_MM:
        gp, fp = E2 / f"lpgrid_R{R}.jsonl", E2 / f"lprefine_R{R}.jsonl"
        if not fp.exists():
            continue
        done.append(R)
        g, f = A.load(gp), A.load(fp)
        errors += [r["case_id"] for r in g + f if not A.ok(r)]
        for r in g + f:
            lcols.setdefault(key(r), []).append(r)
        for r in g:
            if A.ok(r):
                lgrid[(key(r), r["mach"])] = r
    rows = []
    for k in sorted(lcols):
        b, v = A.boundary(bcols[k]), A.boundary(lcols[k])
        shift = None
        if b["M_s"] is not None and v["M_s"] is not None:
            shift = round((v["M_s"] - b["M_s"]) / 0.25, 2)
        rows.append(dict(R_n_mm=k[0], catalysis=k[1], alt_km=k[2], M_s_base=b["M_s"], M_f_base=b["M_f"], mode_base=b["mode_f"],
                         t_fail_base=b["t_fail_f"], M_s_local=v["M_s"], M_f_local=v["M_f"], mode_local=v["mode_f"],
                         t_fail_local=v["t_fail_f"], shift_steps=shift))
    if not rows:
        print("no finished radius yet")
        return
    with open(OUT / "lpgrid_boundaries.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    pts = []
    for (k, M), rl in lgrid.items():
        rb = bgrid.get((k, M))
        if rb is not None:
            pts.append((k, M, rb, rl))
    changed_outcome = [(k, M) for k, M, rb, rl in pts if (rb["mode"] == "none") != (rl["mode"] == "none")]
    dts = [100 * (rl["t_fail"] / rb["t_fail"] - 1) for k, M, rb, rl in pts
           if rb["mode"] != "none" and rl["mode"] != "none" and rb["t_fail"] > 0]
    mode_changes = sum(1 for k, M, rb, rl in pts if rb["mode"] != rl["mode"])
    summ = dict(radii_done=done, n_columns=len(rows), errors=errors,
                columns_moved=sum(1 for r in rows if r["shift_steps"] not in (None, 0)),
                columns_down=sum(1 for r in rows if r["shift_steps"] is not None and r["shift_steps"] < 0),
                columns_up=sum(1 for r in rows if r["shift_steps"] is not None and r["shift_steps"] > 0),
                moved=[r for r in rows if r["shift_steps"] not in (None, 0)],
                n_paired_grid=len(pts), grid_outcome_changes=len(changed_outcome), grid_mode_changes=mode_changes,
                dt_pct_min=min(dts) if dts else None, dt_pct_max=max(dts) if dts else None)
    (OUT / "lpgrid_summary.json").write_text(json.dumps(summ, indent=1, default=str), encoding="utf-8")
    print(json.dumps({k: v for k, v in summ.items() if k != "moved"}, default=str))
    for r in summ["moved"]:
        print("  moved:", r)


if __name__ == "__main__":
    main()
