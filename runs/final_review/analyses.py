"""Final review (8 Nov 2026): numbers asked for by the reviewers that need no new simulation.

Reads saved outputs (runs/analysis/boundaries_2d.csv, runs/envelope_2d/*.csv/jsonl, runs/sensitivity/*) and evaluates the
model's own thermochemistry and atmosphere. Writes runs/final_review/analyses.json and analyses.md.

  1 rarefaction   Kolodziej (1997, NASA TM-112204) Eq. 7: q_rare/q_cont as a polynomial in X = log10(Re_t2),
                  Re_t2 = rho_t2 R_n sqrt(H_s) / mu_t2 (nose tip, n = 0), at the highest passing 2D conditions (R2-04)
  2 laminar       free-stream Reynolds number over the wetted length to the end of the 30 mm flank (R2-10)
  3 eq2_cap       Eq. 2 against the fully catalytic potential and against the frozen sensible potential (R2-03)
  4 interfacial   gas pressure of SiC + 2 SiO2(l) = 3 SiO + CO at the reaction front against p_s (R1-04)
  5 tolerance     columns whose pass or fail is decided inside the 1.02 % time tolerance (R3-17)
  6 converted     converted plus SiC-depleted depth at the highest passing 2D conditions (R3-15)
  7 five_minutes  boundary for a 5-minute instead of a 30-minute exposure, from the saved 2D runs (R4-08)
  8 prcc_ci       bootstrap 95 % intervals of the partial rank correlations (R3-13)
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "runs"))
from uhtc import atmosphere as atm  # noqa: E402
from uhtc import heating as ht  # noqa: E402
from uhtc import thermo as th  # noqa: E402
import common  # noqa: E402

OUT = Path(__file__).resolve().parent
res = {}


def rows_csv(p):
    return list(csv.DictReader(open(p, encoding="utf-8")))


bnd = rows_csv(ROOT / "runs" / "analysis" / "boundaries_2d.csv")
fc = [r for r in bnd if r["catalysis"] == "fc"]

# 1 rarefaction
A = [-1.38151, 2.26375, -0.54005, -0.02164, 0.02016, -0.00170]
rare = []
for r in fc:
    R, z, M = float(r["R_n_mm"]) * 1e-3, float(r["alt_km"]) * 1e3, float(r["M_s_2d"])
    s = common.stag_state(M, z)
    Re = s.rho_e * R * math.sqrt(s.H_s) / s.mu_e
    X = math.log10(Re)
    ratio = sum(a * X ** i for i, a in enumerate(A)) if Re < 1e5 else 1.0
    T, p, rho = atm.state(z)
    mu = float(atm.viscosity(T))
    lam = mu / rho * math.sqrt(math.pi / (2 * 287.05 * T))
    rare.append(dict(R_mm=r["R_n_mm"], alt_km=r["alt_km"], mach=M, Re_t2=Re, q_ratio=ratio, Kn=lam / R))
res["rarefaction"] = rare


def span(vals):
    return [min(vals), max(vals)] if vals else None


summ = {}
for R in ("1.0", "3.0", "10.0"):
    for zl, zh in ((20, 22), (30, 30), (35, 40)):
        v = [x for x in rare if x["R_mm"] == R and zl <= float(x["alt_km"]) <= zh]
        summ[f"{R} mm {zl}-{zh} km"] = dict(Re_t2=span([x["Re_t2"] for x in v]), q_ratio=span([x["q_ratio"] for x in v]),
                                            Kn=span([x["Kn"] for x in v]))
res["rarefaction_summary"] = summ

# 2 laminar: wetted length from the stagnation point to the end of the flank, 5 deg cone, 30 mm flank; all three radii
# (third review pass R5p3-02: the first version used the 1 mm tip only), plus the value at each column's highest passing Mach number
lam = []
for R in (0.001, 0.003, 0.01):
    for z in (20e3, 30e3, 40e3):
        T, p, rho = atm.state(z)
        a = float(atm.speed_of_sound(z))
        mu = float(atm.viscosity(T))
        L = R * (math.pi / 2 - math.radians(5)) + 0.03
        Ms = [float(r["M_s_2d"]) for r in bnd if r["catalysis"] == "fc" and float(r["R_n_mm"]) == R * 1e3
              and float(r["alt_km"]) == z / 1e3 and r["M_s_2d"] not in ("", "None")]
        for M in sorted({8.5, 10.0, 12.0, *Ms}):
            lam.append(dict(R_mm=R * 1e3, alt_km=z / 1e3, mach=M, at_boundary=bool(Ms and M == Ms[0]),
                            Re_x=rho * M * a * L / mu, L_mm=L * 1e3))
res["laminar"] = lam

# 3 Eq. 2 and the cap
cap = []
for M, z, Tw, label in ((9.0, 30e3, 2700.0, "1D settled state, Mach 9, 30 km"),
                        (8.5, 20e3, 2300.0, "2D boundary, Mach 8.5, 20 km"),
                        (8.0, 35e3, 2320.0, "B1, Mach 8, 35 km"),
                        (10.25, 40e3, 2300.0, "2D boundary, Mach 10.25, 40 km")):
    s = common.stag_state(M, z)
    fcp = s.H_s - th.h_air_equilibrium(Tw, s.p_s)
    eq2 = ht.noncatalytic_factor(s.H_s, s.h_D) * (s.H_s - th.h_air_frozen(Tw))
    frozen = s.H_s - s.h_D - th.h_air_frozen(Tw)
    cap.append(dict(case=label, T_w=Tw, fc_MJkg=fcp / 1e6, eq2_uncapped_MJkg=eq2 / 1e6, frozen_MJkg=frozen / 1e6,
                    eq2_over_fc=eq2 / fcp, frozen_over_fc=frozen / fcp, cap_active=bool(eq2 > fcp)))
res["eq2_cap"] = cap

# 4 interfacial gas pressure, stoichiometric SiO:CO = 3:1, activities of SiC and SiO2(l) = 1
def p_interface(T):
    K = math.exp(-float(th.dG_reaction(T, {"SiC": 1, "SiO2(l)": 2}, {"SiO": 3, "CO": 1})) / (th.R * T))
    x = (K / 27.0) ** 0.25  # p_CO in bar; p_SiO = 3 p_CO
    return 4.0 * x * 1e5  # Pa


front = [json.loads(l) for l in open(ROOT / "runs" / "envelope_2d" / "front.jsonl", encoding="utf-8") if l.strip()]
fr = []
for f in front:
    if f.get("error") or f.get("catalysis") != "fc":
        continue
    s = common.stag_state(f["mach"], f["alt"])
    pi = p_interface(f["T_f_max"])
    fr.append(dict(case=f["case_id"], T_f_max=f["T_f_max"], p_int_kPa=pi / 1e3, p_s_kPa=s.p_s / 1e3, exceeds=pi > s.p_s))
res["interfacial_by_T"] = {str(T): p_interface(T) / 1e3 for T in (2073.0, 2200.0, 2300.0, 2400.0)}
passing = [x for x in fr if "front_s_" in x["case"] or "_s_" in x["case"]]
res["interfacial"] = dict(n_cases=len(fr), n_exceed=sum(x["exceeds"] for x in fr),
                          n_passing=len(passing), n_passing_exceed=sum(x["exceeds"] for x in passing),
                          ratio_span=span([x["p_int_kPa"] / x["p_s_kPa"] for x in fr]))

# 5 tolerance: pass or fail decided inside the 1.02 % time tolerance (failing within [1800/1.0102, 1800])
tol = 1800.0 / 1.0102
near = [dict(R=r["R_n_mm"], cat=r["catalysis"], alt=r["alt_km"], t=float(r["t_fail_2d"]))
        for r in bnd if r["t_fail_2d"] and float(r["t_fail_2d"]) > tol]
res["tolerance_baseline"] = near
var = rows_csv(ROOT / "runs" / "analysis" / "variants_2d.csv")
tcol = [k for k in var[0] if k.startswith("t_fail")]
res["tolerance_variants"] = sum(1 for r in var for k in tcol if r[k] and float(r[k]) > tol)

# 6 converted depth at the highest passing 2D conditions (grid and refine)
runs2d = rows_csv(ROOT / "runs" / "envelope_2d" / "grid.csv") + rows_csv(ROOT / "runs" / "envelope_2d" / "refine.csv")
idx = {(round(float(r["R_n"]) * 1e3, 3), r["catalysis"], round(float(r["alt"]) / 1e3, 3), float(r["mach"])): r for r in runs2d}
conv = {}
for r in fc:
    k = (float(r["R_n_mm"]), "fc", float(r["alt_km"]), float(r["M_s_2d"]))
    g = idx.get(k)
    if g:
        conv.setdefault(r["R_n_mm"], []).append(float(g["ZrO2_um_end"]) + float(g["dep_um_end"]))
res["converted_um"] = {R: span(v) for R, v in conv.items()}

# 7 five-minute boundary: highest grid Mach whose run lasts >= 300 s, against the 30-minute boundary
cols = {}
for r in runs2d:
    if r["catalysis"] != "fc" or r.get("error"):
        continue
    k = (round(float(r["R_n"]) * 1e3, 3), round(float(r["alt"]) / 1e3, 3))
    tf = float(r["t_fail"]) if r["t_fail"] not in ("", "nan") else math.inf
    cols.setdefault(k, []).append((float(r["mach"]), tf if not math.isnan(tf) else math.inf, r["mode"]))
five = {}
for (R, z), lst in cols.items():
    lst.sort()
    s30 = max([M for M, t, m in lst if t >= 1800.0 - 1e-6] or [None], key=lambda x: -1 if x is None else x)
    s5 = max([M for M, t, m in lst if t >= 300.0] or [None], key=lambda x: -1 if x is None else x)
    f5 = [(M, m) for M, t, m in lst if t < 300.0]
    mode5 = min(f5)[1] if f5 else "none"
    five.setdefault(R, []).append(dict(alt=z, M30=s30, M5=s5, gain=(s5 - s30) if (s5 is not None and s30 is not None) else None,
                                       mode_first_fail_5min=mode5))
res["five_minutes"] = {str(R): dict(gain=span([x["gain"] for x in v if x["gain"] is not None]),
                                    modes=sorted({x["mode_first_fail_5min"] for x in v})) for R, v in five.items()}

# 8 bootstrap of the partial rank correlations (same design and outputs as runs/analysis/summarize.py)
from scipy.stats import rankdata  # noqa: E402


def prcc(X, y):
    Xr = np.column_stack([rankdata(X[:, j]) for j in range(X.shape[1])])
    yr = rankdata(y)
    out = []
    for j in range(X.shape[1]):
        o = np.column_stack([np.ones(len(yr)), np.delete(Xr, j, axis=1)])
        rx = Xr[:, j] - o @ np.linalg.lstsq(o, Xr[:, j], rcond=None)[0]
        ry = yr - o @ np.linalg.lstsq(o, yr, rcond=None)[0]
        out.append(float(np.corrcoef(rx, ry)[0, 1]))
    return np.array(out)


ci = {}
for bn in ("B1", "B2"):
    resr = {json.loads(l)["case_id"]: json.loads(l) for l in open(ROOT / "runs" / "sensitivity" / f"lhs_{bn}.jsonl", encoding="utf-8") if l.strip()}
    design = rows_csv(ROOT / "runs" / "sensitivity" / f"lhs_design_{bn}.csv")
    names = [k for k in design[0] if k != "case_id"]
    X, y = [], []
    for d in design:
        r = resr.get(d["case_id"])
        if r is None or r.get("error"):
            continue
        t = r["t_fail"] if r["t_fail"] is not None and not (isinstance(r["t_fail"], float) and math.isnan(r["t_fail"])) else 1800.0
        X.append([1.0 if d[k] == "fc" else (0.0 if d[k] == "nc" else float(d[k])) for k in names])
        y.append(math.log10(min(t, 1800.0)))
    X, y = np.array(X), np.array(y)
    base = prcc(X, y)
    rng = np.random.default_rng(20261108)
    boots = np.array([prcc(X[i], y[i]) for i in (rng.integers(0, len(y), len(y)) for _ in range(1000))])
    lo, hi = np.percentile(boots, [2.5, 97.5], axis=0)
    ci[bn] = {k: dict(prcc=round(float(base[j]), 2), lo=round(float(lo[j]), 2), hi=round(float(hi[j]), 2),
                      excludes_zero=bool(lo[j] > 0 or hi[j] < 0)) for j, k in enumerate(names)}
res["prcc_ci"] = ci

(OUT / "analyses.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
md = ["# Final-review analyses (runs/final_review/analyses.py)", ""]
md.append("## Rarefaction (Kolodziej Eq. 7) and Knudsen number at the highest passing 2D conditions")
for k, v in summ.items():
    md.append(f"- {k}: Re_t2 {v['Re_t2']}, q_rare/q_cont {v['q_ratio']}, Kn {v['Kn']}")
md.append("\n## Laminar flow: free-stream Re over the wetted length (cap and 30 mm flank), 1, 3 and 10 mm tips")
for x in lam:
    md.append(f"- {x['R_mm']:g} mm, {x['alt_km']:.0f} km, Mach {x['mach']}{' (2D boundary)' if x['at_boundary'] else ''}: "
              f"Re_x = {x['Re_x']:.3g}")
mb = max(x["Re_x"] for x in lam if x["at_boundary"])
md.append(f"- largest at a 2D boundary: {mb:.3g}; largest up to Mach 12: {max(x['Re_x'] for x in lam):.3g}")
md.append("\n## Eq. 2 and its cap")
for x in cap:
    md.append(f"- {x['case']}, T_w {x['T_w']:.0f} K: Eq. 2/fc {x['eq2_over_fc']:.3f}, frozen/fc {x['frozen_over_fc']:.3f}, cap active {x['cap_active']}")
md.append(f"\n## Interfacial gas pressure (kPa) by temperature: {res['interfacial_by_T']}")
md.append(f"- front.jsonl, fully catalytic: {res['interfacial']}")
md.append(f"\n## Inside the time tolerance: baseline columns {len(near)} ({near}); variant entries {res['tolerance_variants']}")
md.append(f"\n## Converted plus depleted depth at the highest passing 2D conditions (um): {res['converted_um']}")
md.append(f"\n## Five-minute boundary minus 30-minute boundary (fully catalytic, grid step 0.5): {res['five_minutes']}")
md.append("\n## PRCC with bootstrap 95 % intervals")
for bn, d in ci.items():
    md.append(f"- {bn}: " + "; ".join(f"{k} {v['prcc']} [{v['lo']}, {v['hi']}]" for k, v in sorted(d.items(), key=lambda kv: -abs(kv[1]['prcc']))))
(OUT / "analyses.md").write_text("\n".join(md) + "\n", encoding="utf-8")
print("\n".join(md))
