"""Summaries of the production runs used in the paper (every number quoted in the text comes from here or from the
run files themselves).  Output: runs/analysis/summary.json and envelope_boundaries.csv.

Envelope (sets A-C): for each radius, catalysis bound and SiC fraction
  * fraction of grid points surviving the 1800 s cap, and the split of failures by first mode;
  * per altitude: the highest Mach number that survives 1800 s and the highest that survives 60 s
    (largest grid Mach with t_fail >= threshold such that every lower Mach also does; grid step 0.5);
  * median time to failure among failed points.
Trajectory: summary per radius/catalysis/SiC fraction (copied from x43a_summary.csv).
Sensitivity: OAT effects and LHS rank statistics (also written to runs/sensitivity/lhs_stats.json for Fig. 8).
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs" / "analysis"
OUT.mkdir(parents=True, exist_ok=True)


def load_jsonl(p):
    return [json.loads(l) for l in open(p, encoding="utf-8") if l.strip()]


def tfail(r):
    t = r["t_fail"]
    return t if t == t else math.inf


def envelope():
    rows = []
    for s in ("A", "B", "C"):
        f = ROOT / "runs" / "envelope" / f"set_{s}.jsonl"
        if f.exists():
            rows += load_jsonl(f)
    errors = [r for r in rows if r.get("error")]
    groups = {}
    for r in rows:
        if r.get("error"):
            continue
        groups.setdefault((r["R_n"], r["catalysis"], r["f_SiC"]), []).append(r)
    summ, bnd = {}, []
    for (R, cat, f), g in sorted(groups.items()):
        key = f"R{R * 1e3:g}mm_{cat}_SiC{f:g}"
        n = len(g)
        surv = [r for r in g if tfail(r) >= 1800]
        failed = [r for r in g if tfail(r) < 1800]
        modes = {m: sum(1 for r in failed if r["mode"] == m) for m in ("glass", "active", "temperature")}
        alts = sorted({r["alt"] for r in g})
        per_alt = {}
        for z in alts:
            rz = sorted((r for r in g if r["alt"] == z), key=lambda r: r["mach"])
            out = {}
            for thr in (1800.0, 60.0):
                mmax = None
                for r in rz:
                    if tfail(r) >= thr:
                        mmax = r["mach"]
                    else:
                        break
                out[f"Mmax_{int(thr)}s"] = mmax
            per_alt[f"{z / 1e3:g}"] = out
            bnd.append(dict(R_n_mm=R * 1e3, catalysis=cat, f_SiC=f, alt_km=z / 1e3, **out))
        summ[key] = dict(n=n, n_survive=len(surv), frac_survive=len(surv) / n, failed_by_mode=modes,
                         median_tfail_failed=float(np.median([tfail(r) for r in failed])) if failed else None,
                         T_s_max_range=[min(r["T_s_max"] for r in g), max(r["T_s_max"] for r in g)],
                         per_alt=per_alt)
    with open(OUT / "envelope_boundaries.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(bnd[0]))
        w.writeheader()
        w.writerows(bnd)
    return summ, len(rows), len(errors), [r["case_id"] for r in errors][:20]


def sensitivity():
    import importlib.util
    spec = importlib.util.spec_from_file_location("rs", ROOT / "runs" / "sensitivity" / "run_sensitivity.py")
    rs = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(ROOT / "runs"))
    spec.loader.exec_module(rs)
    bases = json.loads((ROOT / "runs" / "sensitivity" / "base_cases.json").read_text(encoding="utf-8"))
    out, stats = {}, {}
    for bn in bases:
        f = ROOT / "runs" / "sensitivity" / f"oat_{bn}.jsonl"
        if not f.exists():
            continue
        rows = load_jsonl(f)
        base = [r for r in rows if r["case_id"] == f"{bn}_base"][0]
        tb = min(tfail(base), 1800.0)
        eff = {}
        for r in rows:
            if r["case_id"] == f"{bn}_base":
                continue
            k, tag = r["case_id"][len(bn) + 1:].rsplit("_", 1)
            eff.setdefault(k, {})[tag] = dict(t_fail=min(tfail(r), 1800.0), mode=r["mode"],
                                              dlog10=math.log10(min(tfail(r), 1800.0) / tb))
        out[bn] = dict(base_t_fail=tb, base_mode=base["mode"], effects=eff)
        fl = ROOT / "runs" / "sensitivity" / f"lhs_{bn}.jsonl"
        fd = ROOT / "runs" / "sensitivity" / f"lhs_design_{bn}.csv"
        if fl.exists() and fd.exists():
            res = {r["case_id"]: r for r in load_jsonl(fl)}
            design = list(csv.DictReader(open(fd, encoding="utf-8")))
            names = [k for k in design[0] if k != "case_id"]
            X, y = [], []
            for d in design:
                r = res.get(d["case_id"])
                if r is None or r.get("error"):
                    continue
                X.append([1.0 if d[k] == "fc" else (0.0 if d[k] == "nc" else float(d[k])) for k in names])
                y.append(math.log10(min(tfail(r), 1800.0)))
            X, y = np.array(X), np.array(y)
            from scipy.stats import rankdata, spearmanr
            rho = {k: float(spearmanr(X[:, j], y)[0]) for j, k in enumerate(names)}
            # partial rank correlation: correlate residuals of ranked y and ranked x_j after regressing on the others
            Xr = np.column_stack([rankdata(X[:, j]) for j in range(X.shape[1])])
            yr = rankdata(y)
            prcc = {}
            for j, k in enumerate(names):
                others = np.column_stack([np.ones(len(yr)), np.delete(Xr, j, axis=1)])
                bx = np.linalg.lstsq(others, Xr[:, j], rcond=None)[0]
                by = np.linalg.lstsq(others, yr, rcond=None)[0]
                rx, ry = Xr[:, j] - others @ bx, yr - others @ by
                prcc[k] = float(np.corrcoef(rx, ry)[0, 1])
            stats[bn] = dict(n=len(y), spearman=rho, prcc=prcc,
                             frac_survive=float(np.mean(np.array(y) >= math.log10(1799.0))),
                             tfail_quantiles_s=[float(10 ** q) for q in np.quantile(y, [0.05, 0.25, 0.5, 0.75, 0.95])])
    (ROOT / "runs" / "sensitivity" / "lhs_stats.json").write_text(json.dumps(stats, indent=1), encoding="utf-8")
    return out, stats


def main():
    env, n_rows, n_err, err_ids = envelope()
    summary = dict(envelope=env, envelope_rows=n_rows, envelope_errors=n_err, envelope_error_ids=err_ids)
    tf = ROOT / "runs" / "trajectory" / "x43a_summary.csv"
    if tf.exists():
        summary["x43a"] = list(csv.DictReader(open(tf, encoding="utf-8")))
    qf = ROOT / "runs" / "q5_2d" / "q5_results.csv"
    if qf.exists():
        summary["q5"] = list(csv.DictReader(open(qf, encoding="utf-8")))
    for name in ("q5_sensitivity", "q5_envelope", "q5_envelope_refine"):
        fq = ROOT / "runs" / "q5_2d" / f"{name}.csv"
        if fq.exists():
            summary[name] = list(csv.DictReader(open(fq, encoding="utf-8")))
    fx = ROOT / "runs" / "sensitivity" / "oat_extra.jsonl"
    if fx.exists():
        summary["oat_extra"] = load_jsonl(fx)
    if (ROOT / "runs" / "sensitivity" / "base_cases.json").exists():
        summary["oat"], summary["lhs"] = sensitivity()
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=str), encoding="utf-8")
    for k, v in env.items():
        print(k, f"survive {v['n_survive']}/{v['n']}", v["failed_by_mode"], "median t_fail", v["median_tfail_failed"])
    print("rows", n_rows, "errors", n_err, err_ids)


if __name__ == "__main__":
    main()
