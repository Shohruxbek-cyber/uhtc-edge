"""VVU plan Section 10: per-quantity statistics of the existing check results (no new runs).

Reads validation/V6_oxidation/v6_results.csv and the stored V7, V8 and V9 results. For each check the points are split by
quantity class; for each class: n, the number within the check's own band, and for thicknesses and times the geometric-mean
model/measured ratio and the geometric standard deviation (for temperatures the mean and spread of model minus measured).
V6 keeps the calibration set, the independent points, the two zero measurements and the glass apart. The text values must be
reproduced (15 of 18 oxide-type points within a factor of 3, geometric mean 1.09; 9 of 19 depleted-layer points), or the
script stops (plan Section 10).
Output: per_quantity.json, per_quantity.md
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
C2K = 273.15


def gstats(ratios):
    lr = np.log(np.asarray(ratios, float))
    return float(np.exp(lr.mean())), (float(np.exp(lr.std(ddof=1))) if len(lr) > 1 else float("nan"))


def v6_class(q):
    q = q.lower()
    if q.startswith("sic-depleted"):
        return "SiC-depleted layer"
    if q.startswith("borosilicate"):
        return "glass (outside the criterion)"
    if q.startswith("zro2"):
        return "ZrO2 layer"
    if q.startswith("recession"):
        return "recession"
    if q.startswith("scale + depleted"):
        return "scale + depleted zone"
    return "total oxide"


def v6():
    rows = list(csv.DictReader(open(ROOT / "validation/V6_oxidation/v6_results.csv", encoding="utf-8")))
    out = []
    for setname in ("independent", "calibration", "glass (depends on furnace g)"):
        sub = [r for r in rows if r["set"] == setname]
        for cls in sorted({v6_class(r["quantity"]) for r in sub}):
            pts = [r for r in sub if v6_class(r["quantity"]) == cls]
            zero = [r for r in pts if float(r["measured_um"]) == 0.0]
            pos = [r for r in pts if float(r["measured_um"]) > 0.0]
            rat = [float(r["ratio"]) for r in pos]
            rat_i = [float(r["ratio_inc"]) for r in pos]
            gm, gsd = gstats(rat)
            gmi, gsdi = gstats(rat_i)
            out.append(dict(check="V6", set=setname, cls=cls, n=len(pos), n_zero=len(zero),
                            within=sum(1 / 3 <= x <= 3 for x in rat), gm=gm, gsd=gsd,
                            within_inc=sum(1 / 3 <= x <= 3 for x in rat_i), gm_inc=gmi, gsd_inc=gsdi, band="x3"))
    ind = [r for r in rows if r["set"] == "independent" and float(r["measured_um"]) > 0]
    oxq = [r for r in ind if v6_class(r["quantity"]) != "SiC-depleted layer"]
    dep = [r for r in ind if v6_class(r["quantity"]) == "SiC-depleted layer"]
    n_ox, n_dep = sum(1 / 3 <= float(r["ratio"]) <= 3 for r in oxq), sum(1 / 3 <= float(r["ratio"]) <= 3 for r in dep)
    gm_ox = gstats([float(r["ratio"]) for r in oxq])[0]
    summary = dict(oxide_type=f"{n_ox} of {len(oxq)}", gm_oxide_type=round(gm_ox, 2), depleted=f"{n_dep} of {len(dep)}")
    if (n_ox, len(oxq), round(gm_ox, 2), n_dep, len(dep)) != (15, 18, 1.09, 9, 19):
        raise SystemExit(f"V6 statistics do not reproduce the text (15 of 18, 1.09, 9 of 19): {summary}; stop (plan Section 10)")
    return out, summary


def temps(check, cls, case, pairs, band):
    d = [s - m for s, m in pairs]
    return dict(check=check, set=case, cls=cls, n=len(d), within=sum(abs(x) <= b for x, b in zip(d, band)),
                mean_diff_K=float(np.mean(d)), min_diff_K=float(min(d)), max_diff_K=float(max(d)), band="+-" + "/".join(
                    f"{b:.0f}" for b in sorted(set(band))) + " K")


def thick(check, cls, case, pairs, f):
    rat = [s / m for s, m in pairs]
    gm, gsd = gstats(rat)
    return dict(check=check, set=case, cls=cls, n=len(rat), within=sum(1 / f <= x <= f for x in rat), gm=gm, gsd=gsd,
                band=f"x{f:g}")


def v7_v9():
    out = []
    v7 = json.loads((ROOT / "validation/V7_arcjet/v7_results.json").read_text())
    T7, TOL7 = [1573.0, 1723.0, 1843.0, 1923.0, 2053.0], [80.0, 80.0, 60.0, 60.0, 60.0]
    for cat in ("nc", "fc"):
        r = v7[f"2d-{cat}"]
        out.append(temps("V7", "plateau temperature", f"2D, {cat}", list(zip(r["Tend"], T7)), TOL7))
        out.append(thick("V7", "oxide at the stagnation point", f"2D, {cat}", [(r["oxide_um"], 150.0)], 1.5))
        out.append(dict(check="V7", set=f"2D, {cat}", cls="SiC-depleted layer (none observed)", n=1,
                        within=int(r["dep_um"] < 1.0), band="< 1 um", model_um=r["dep_um"]))
    v8 = json.loads((ROOT / "validation/V8_measured_flux/v8_results.json").read_text())
    for cat in ("nc", "fc"):
        pas = [v8[f"V8_1_H{h}_{cat}"]["passive"] for h in (15, 30, 45)]
        out.append(temps("V8-1", "peak temperature", f"{cat}, 15/30/45 MJ/kg", [(p["T_max"], 1630 + C2K) for p in pas], [80.0] * 3))
        out.append(thick("V8-1", "oxide at 300 s", f"{cat}, 15/30/45 MJ/kg", [(p["oxide_um"][-1], 110.0) for p in pas], 1.5))
        out.append(thick("V8-1", "SiC-depleted layer at 300 s", f"{cat}, 15/30/45 MJ/kg", [(p["dep_um"][-1], 40.0) for p in pas], 2.0))
        h = v8[f"V8_2_{cat}"]
        out.append(temps("V8-2", "steady temperature", cat, [(h["c1"]["T_end"], 1960 + C2K), (h["c2"]["T_end"], 1900 + C2K)], [80.0] * 2))
        lay = [(float(np.interp(t, h["c1"]["t"], h["c1"]["oxide_um"])), m) for t, m in ((201, 65), (402, 99), (600, 125))]
        out.append(thick("V8-2", "oxide, conditions 1 (201/402/600 s) and 2", cat, lay + [(h["c2"]["oxide_um"][-1], 90.0)], 1.5))
    v9 = json.loads((ROOT / "validation/V9_plasmatron/v9_results.json").read_text())
    meas = {"1.8": (1895.0, 2110.0, None), "1.5": (1825.0, 2210.0, 68.0), "1.6": (1860.0, 2250.0, 21.0)}
    for cat in ("nc", "fc"):
        pairs = []
        for s, (tl, th, on) in meas.items():
            r = v9[f"{cat}_{s}"]
            pairs += [(r["T_low"], tl), (r["T_ref"], th)]
        out.append(temps("V9", "low- and high-stage temperature", cat, pairs, [80.0] * 6))
        times = [(v9[f"{cat}_{s}"]["after_switch"], on) for s, (_, _, on) in meas.items() if on is not None]
        lost = [(a, m) for a, m in times if a == a]
        rec = thick("V9", "glass-loss time after the switch", cat, lost, 3.0) if lost else dict(
            check="V9", set=cat, cls="glass-loss time after the switch", n=0, within=0, band="x3")
        rec["no_loss_by_480_s"] = sum(1 for a, _ in times if a != a)
        out.append(rec)
    return out


def main():
    v6rows, summary = v6()
    rows = v6rows + v7_v9()
    (HERE / "per_quantity.json").write_text(json.dumps(dict(v6_summary=summary, rows=rows), indent=1, default=float))
    L = ["# Per-quantity statistics of the existing checks (VVU plan Section 10; no new runs)", "",
         f"V6 reproduces the text: oxide-type independent points {summary['oxide_type']} within a factor of 3, geometric mean "
         f"{summary['gm_oxide_type']}; SiC-depleted layer {summary['depleted']}.", "",
         "| Check | Set / case | Quantity class | n | Within band | Band | Geometric mean ratio (GSD) or mean model − measured (range) | With incubation |",
         "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        if "gm" in r:
            stat = f"{r['gm']:.2f} ({r['gsd']:.2f})" if r["n"] > 1 else (f"{r['gm']:.2f}" if r["n"] else "—")
        elif "mean_diff_K" in r:
            stat = f"{r['mean_diff_K']:+.0f} K ({r['min_diff_K']:+.0f} to {r['max_diff_K']:+.0f})"
        else:
            stat = f"model {r['model_um']:.0f} µm" if "model_um" in r else "—"
        inc = f"{r['within_inc']} of {r['n']}, {r['gm_inc']:.2f}" if "within_inc" in r and r["n"] else ""
        extra = f"; no glass loss by 480 s in {r['no_loss_by_480_s']}" if "no_loss_by_480_s" in r else ""
        zero = f" (+{r['n_zero']} zero)" if r.get("n_zero") else ""
        L.append(f"| {r['check']} | {r['set']} | {r['cls']} | {r['n']}{zero} | {r['within']} | {r['band']} | {stat}{extra} | {inc} |")
    (HERE / "per_quantity.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    sys.exit(main())
