"""Re-trace the numbers quoted in paper/main.tex to the saved run outputs (part of REVIEW.md).

Every quoted number is recomputed here from the raw outputs in runs/ and validation/, independently of the code that
wrote the paper text, and compared with the value in the text.  Conventions of the paper: failure times rounded down,
temperatures to the nearest 10 K.  A check passes when the quoted value is what that rounding gives (or, for ranges and
bounds, when the recomputed values lie inside the quoted range).
Each section runs on its own: an exception inside a section is recorded as a failed check and the others still run.
Output: runs/analysis/trace_numbers.md (table) and the exit code (number of failed checks).
Revision (Nov 2026): rewritten for the revised text (2D headline, criterion (a), validation table, Appendix A).
"""
import csv
import json
import math
import re
import sys
import traceback
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "runs"))
sys.path.insert(0, str(ROOT / "src"))
import common  # noqa: E402
from uhtc import data as uhtc_data  # noqa: E402

RESULTS = []
SECTIONS = []


def section(fn):
    SECTIONS.append(fn)
    return fn


def rows(path):
    with open(ROOT / path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def jl(path):
    with open(ROOT / path, encoding="utf-8") as fh:
        return [json.loads(x) for x in fh if x.strip()]


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def check(sec, claim, quoted, value, ok):
    RESULTS.append((sec, claim, str(quoted), value if isinstance(value, str) else f"{value}", bool(ok)))


def down(x, q=1.0):
    return math.floor(x / q + 1e-9) * q


def down2(x):
    """failure times >= 200 s are quoted rounded down to 10 s, shorter ones to 1 s"""
    return down(x, 10) if x >= 200 else down(x)


def near10(x):
    return round(x / 10.0) * 10.0


def tf(r):
    v = f(r["t_fail"])
    return 1800.0 if v != v else v


def isnan(x):
    return f(x) != f(x)


# ------------------------------------------------------------------------------------------------ Validation
@section
def validation():
    V6 = rows("validation/V6_oxidation/v6_results.csv")
    ind = [r for r in V6 if r["set"] == "independent" and f(r["measured_um"]) > 0]  # zero readings have no ratio
    w3 = lambda x: x == x and 1 / 3 <= x <= 3
    n = sum(w3(f(r["ratio"])) for r in ind)
    check("Val", "V6 independent points within x3", "24 of 37", f"{n} of {len(ind)}", n == 24 and len(ind) == 37)
    cal = [f(r["ratio"]) for r in V6 if r["set"] == "calibration"]
    gm = math.exp(np.mean(np.log(cal)))
    check("Val", "V6 calibration geometric-mean ratio", 1.22, f"{gm:.3f}", round(gm, 2) == 1.22)
    dep = [r for r in ind if r["quantity"].lower().startswith("sic-depleted")]
    ox = [r for r in ind if r not in dep]
    gmo = math.exp(np.mean(np.log([f(r["ratio"]) for r in ox])))
    no = sum(w3(f(r["ratio"])) for r in ox)
    check("Val", "V6 oxide-type points within x3, geometric mean", "15 of 18, 1.09", f"{no} of {len(ox)}, {gmo:.3f}",
          no == 15 and len(ox) == 18 and round(gmo, 2) == 1.09)
    nd = sum(w3(f(r["ratio"])) for r in dep)
    check("Val", "V6 depleted-layer points within x3", "9 of 19", f"{nd} of {len(dep)}", nd == 9 and len(dep) == 19)
    gmi = math.exp(np.mean(np.log([f(r["ratio"]) for r in ind])))
    check("Val", "V6 geometric-mean ratio of all independent points (Supplement S3)", "2.10", f"{gmi:.3f}", round(gmi, 2) == 2.10)
    srcs = sorted({r["source"] for r in V6})
    check("Val", "V6: four studies; no value known only from an abstract (Karlsdottir 2007 dropped, author 10 Nov 2026)",
          "4; none", srcs, len(srcs) == 4 and "Karlsdottir2007" not in srcs)
    sh_ox = [r for r in ox if r["source"] == "Shugart2014"]
    check("Val", "Supplement S3: Shugart points among the independent oxide-type points", "6 of 18", f"{len(sh_ox)} of {len(ox)}",
          len(sh_ox) == 6 and len(ox) == 18)
    ndi = sum(w3(f(r["ratio_inc"])) for r in dep)
    nai = sum(w3(f(r["ratio_inc"])) for r in ind)
    g0 = math.exp(np.mean(np.log([f(r["ratio"]) for r in dep])))
    g1 = math.exp(np.mean(np.log([f(r["ratio_inc"]) for r in dep])))
    check("Val", "V6 incubation: depleted within x3, GM with / without, all within x3",
          "14 of 19, 3.00 / 3.93, 29 of 37", f"{ndi}, {g1:.3f} / {g0:.3f}, {nai}",
          ndi == 14 and round(g1, 2) == 3.00 and round(g0, 2) == 3.93 and nai == 29)
    # V7 (SPES), computed flux
    V7 = json.loads((ROOT / "validation" / "V7_arcjet" / "v7_results.json").read_text(encoding="utf-8"))
    meas = [1573, 1723, 1843, 1923, 2053]  # SPES plateaus, Savino et al. 2010 Fig. 8 inset (validation/V7_arcjet/report.md)
    tol = [80, 80, 60, 60, 60]
    d1 = [a - m for a, m in zip(V7["1d-nc"]["Tend"], meas)]
    d2 = [a - m for a, m in zip(V7["2d-nc"]["Tend"], meas)]
    d3 = [a - m for a, m in zip(V7["2d-fc"]["Tend"], meas)]
    check("Val", "V7 1D non-catalytic over-prediction (K)", "250-540", f"{min(d1):.1f}-{max(d1):.1f}",
          near10(min(d1)) == 250 and near10(max(d1)) == 540)
    check("Val", "V7 2D non-catalytic over-prediction (K)", "16-112", f"{min(d2):.1f}-{max(d2):.1f}",
          round(min(d2)) == 16 and round(max(d2)) == 112)
    check("Val", "V7 2D fully catalytic over-prediction (K)", "250-650", f"{min(d3):.1f}-{max(d3):.1f}",
          near10(min(d3)) == 250 and near10(max(d3)) == 650)
    nok = sum(abs(x) <= t for x, t in zip(d2, tol))
    check("Val", "V7 2D non-catalytic plateaus within tolerance", "2 of 5", nok, nok == 2)
    check("Val", "V7 oxide 2D-nc / measured-T driven (um)", "100 / 83",
          f"{V7['2d-nc']['oxide_um']:.1f} / {V7['_D2']['oxide_um']:.1f}",
          round(V7["2d-nc"]["oxide_um"]) == 100 and round(V7["_D2"]["oxide_um"]) == 83)
    # V8, measured flux
    B = json.loads((ROOT / "validation" / "V8_measured_flux" / "v8_results.json").read_text(encoding="utf-8"))
    C2K = 273.15
    ames_nc = [B[f"V8_1_H{h}_nc"]["passive"]["T_max"] - C2K - 1630 for h in (15, 30, 45)]
    hit_nc = [B["V8_2_nc"]["c1"]["T_end"] - C2K - 1960, B["V8_2_nc"]["c2"]["T_end"] - C2K - 1900]
    fc = [B[f"V8_1_H{h}_fc"]["passive"]["T_max"] - C2K - 1630 for h in (15, 30, 45)] + \
         [B["V8_2_fc"]["c1"]["T_end"] - C2K - 1960, B["V8_2_fc"]["c2"]["T_end"] - C2K - 1900]
    # final review (R2-01): edge states beyond the 6000 K table end are no longer clipped (runs/final_review/edge_states.py)
    check("Val", "V8 non-catalytic over-prediction, Ames 30 MJ/kg (K)", "40", f"{ames_nc[1]:.1f}", near10(ames_nc[1]) == 40)
    check("Val", "V8 non-catalytic over-prediction, Ames 15 and 45 MJ/kg (K)", "170-240",
          f"{min(ames_nc[0], ames_nc[2]):.1f}-{max(ames_nc[0], ames_nc[2]):.1f}",
          near10(min(ames_nc[0], ames_nc[2])) == 170 and near10(max(ames_nc[0], ames_nc[2])) == 240)
    check("Val", "V8 non-catalytic over-prediction, HIT (K)", "7 and 26",
          f"{hit_nc[0]:.1f} and {hit_nc[1]:.1f}", round(hit_nc[0]) == 7 and round(hit_nc[1]) == 26)
    check("Val", "V8 fully catalytic over-prediction (K)", "770-1080", f"{min(fc):.1f}-{max(fc):.1f}",
          near10(min(fc)) == 770 and near10(max(fc)) == 1080)
    oxv = [B[f"V8_1_H{h}_nc"]["passive"]["oxide_um"][-1] for h in (15, 30, 45)]
    dpv = [B[f"V8_1_H{h}_nc"]["passive"]["dep_um"][-1] for h in (15, 30, 45)]
    check("Val", "V8 Ames nc oxide at 30 MJ/kg (um); inside 110 x1.5 at 15 and 45", "58; inside",
          f"{oxv[1]:.1f}; {oxv[0]:.1f}, {oxv[2]:.1f}",
          round(oxv[1]) == 58 and all(110 / 1.5 <= x <= 165 for x in (oxv[0], oxv[2])) and not 110 / 1.5 <= oxv[1])
    check("Val", "V8 Ames nc depleted layer, three enthalpies (um)", "41-45", f"{min(dpv):.1f}-{max(dpv):.1f}",
          round(min(dpv)) == 41 and round(max(dpv)) == 45)
    c1n = B["V8_2_nc"]["c1"]
    ox600 = float(np.interp(600.0, c1n["t"], c1n["oxide_um"]))
    zs = [B["V8_2_nc"][k]["T_end"] for k in ("zs10", "zs15", "c1")]
    check("Val", "V8 HIT nc oxide at 600 s (um); trend T(ZS10) > T(ZS15) > T(ZS30) not reproduced", "241; not reproduced",
          f"{ox600:.1f}; {[round(x - C2K) for x in zs]}", round(ox600) == 241 and not (zs[0] > zs[1] > zs[2]))
    rep = (ROOT / "validation" / "V8_measured_flux" / "report.md").read_text(encoding="utf-8")
    d_ames = float(re.search(r"V8-1 passive.*?300 s: oxide (\d+)", rep).group(1)) / 110.0
    c1 = [float(x) for x in re.search(r"cond\. 1 \(steady T\).*?201 s: oxide (\d+).*?402 s: oxide (\d+).*?600 s: oxide (\d+)",
                                      rep).groups()]
    c2 = float(re.search(r"cond\. 2 \(steady T\).*?600 s: oxide (\d+)", rep).group(1))
    rat = [c1[0] / 65, c1[1] / 99, c1[2] / 125, c2 / 90]
    check("Val", "V8 diagnostic (measured T): oxide ratio Ames / HIT", "0.55 / 2.5-2.9",
          f"{d_ames:.3f} / {min(rat):.2f}-{max(rat):.2f}",
          round(d_ames, 2) == 0.55 and round(min(rat), 1) == 2.5 and round(max(rat), 1) == 2.9)
    # V9, plasmatron (targets of validation/V9_plan.md)
    P = json.loads((ROOT / "validation" / "V9_plasmatron" / "v9_results.json").read_text(encoding="utf-8"))
    tg = {"1.8": (1895, 2110), "1.5": (1825, 2210), "1.6": (1860, 2250)}
    dfc = [P[f"fc_{k}"][q] - v for k, tv in tg.items() for q, v in zip(("T_low", "T_ref"), tv)]
    dnc = [P[f"nc_{k}"][q] - v for k, tv in tg.items() for q, v in zip(("T_low", "T_ref"), tv)]
    check("Val", "V9 fully catalytic over-prediction (K)", "170-280", f"{min(dfc):.1f}-{max(dfc):.1f}", near10(min(dfc)) == 170 and near10(max(dfc)) == 280)
    check("Val", "V9 non-catalytic under-prediction (K)", "90-400", f"{max(dnc):.1f} to {min(dnc):.1f}", near10(-max(dnc)) == 90 and near10(-min(dnc)) == 400)
    check("Val", "V9 glass loss after the switch, fc, sample 1.6 (s)", 71, f"{P['fc_1.6']['after_switch']:.1f}", round(P["fc_1.6"]["after_switch"]) == 71)
    lost = [k for k in P if P[k]["t_glass"] == P[k]["t_glass"]]
    check("Val", "V9 only fc sample 1.6 loses its glass", "fc_1.6", str(lost), lost == ["fc_1.6"])
    # V5
    V4 = rows("validation/V5_heating_distribution/v5_results.csv")
    cl = [f(r["ratio"]) for r in V4 if r["case"].startswith("Cleary")]
    check("Val", "V5 Cleary measured/baseline (27 points)", "0.72-1.15", f"{len(cl)}: {min(1 / x for x in cl):.3f}-{max(1 / x for x in cl):.3f}"
          if False else f"{len(cl)}: {min(cl):.3f}-{max(cl):.3f}", len(cl) == 27 and round(min(cl), 2) == 0.72
          and round(max(cl), 2) == 1.15)
    # V2 2D: production vs refined settings
    dT, dt = [], []
    for name in ("check_mesh_lees", "check_mesh_lees_failing"):
        rr = {r["case_id"]: r for r in jl(f"runs/envelope_2d/{name}.jsonl")}
        for cid, r in rr.items():
            if cid.endswith("_prod"):
                b = rr[cid[:-5] + "_fine"]
                dT.append(abs(r["T_s_max"] / b["T_s_max"] - 1))
                if r["t_fail"] == r["t_fail"]:
                    dt.append(abs(r["t_fail"] / b["t_fail"] - 1))
    check("Val", "V2 2D: peak T (6 cases), failure time (3 failing)", "<= 0.27 %, <= 1.02 %",
          f"{len(dT)}: {100 * max(dT):.3f} %, {len(dt)}: {100 * max(dt):.3f} %",
          len(dT) == 6 and round(100 * max(dT), 2) <= 0.27 and len(dt) == 3 and round(100 * max(dt), 2) <= 1.02)
    cr = [l.split("|") for l in (ROOT / "runs" / "analysis" / "check_results.md").read_text(encoding="utf-8").splitlines()
          if l.startswith("| ") and l.count("|") == 6 and "%" in l]
    pairs = [(float(x[2]), float(x[3])) for x in cr]
    pct = [float(x[5].strip().rstrip("%")) for x in cr]  # the report's own full-precision change per case
    rel = max(pct)
    check("Num", "refined reruns (20 failing map cases): max change, direction", "0.39 %, none longer",
          f"{len(pairs)} cases, {rel:.2f} %, shorter or equal: {all(b <= a for a, b in pairs)}",
          len(pairs) == 20 and round(rel, 2) == 0.39 and all(b <= a for a, b in pairs))


# ------------------------------------------------------------------------------------------------ Q1 (1D histories)
def histories():
    H = rows("runs/histories/histories.csv")
    E = {(r["mach"], r["alt_km"], r["R_n_mm"], r["catalysis"]): r for r in rows("runs/histories/events.csv")}
    hist = defaultdict(list)
    for r in H:
        hist[(r["mach"], r["alt_km"], r["R_n_mm"], r["catalysis"])].append(r)
    for k in hist:
        hist[k].sort(key=lambda r: f(r["t"]))
    return hist, E


@section
def q1():
    hist, E = histories()
    for R, qf, qn in (("1.0", 24, 27), ("3.0", 55, 65)):
        a, b = f(E[("9.0", "30.0", R, "fc")]["t_limit"]), f(E[("9.0", "30.0", R, "nc")]["t_limit"])
        check("Q1", f"M9/30 km, R {R} mm: t_limit fc / nc (s)", f"{qf} / {qn}", f"{a:.2f} / {b:.2f}",
              down(a) == qf and down(b) == qn)
    for R, cat, q in (("1.0", "fc", 2700), ("3.0", "fc", 2570), ("10.0", "fc", 2330), ("10.0", "nc", 2300),
                      ("1.0", "nc", 2700), ("3.0", "nc", 2570)):
        v = f(hist[("9.0", "30.0", R, cat)][-1]["T_s"])
        check("Q1", f"M9/30 km, R {R} mm {cat}: T at the end (K)", q, f"{v:.1f}", near10(v) == q)
    st = common.stag_state(9.0, 30e3)
    check("Q1", "M9/30 km edge temperature (K)", 2950, f"{st.T_e:.2f}", near10(st.T_e) == 2950)
    gap = st.T_e - f(hist[("9.0", "30.0", "1.0", "fc")][-1]["T_s"])
    check("Q1", "M9/30 km: 1 mm surface below edge T (K)", 260, f"{gap:.1f}", near10(gap) == 260)
    gaps = [st.T_e - f(hist[("9.0", "30.0", R, "fc")][-1]["T_s"]) for R in ("1.0", "3.0", "10.0")]
    check("Disc", "M9/30 km: surface below edge T, 1-10 mm (K)", "260-620", f"{min(gaps):.1f}-{max(gaps):.1f}",
          near10(min(gaps)) == 260 and near10(max(gaps)) == 620)
    for (M, z), q in ((("9.0", "30.0"), 120), (("10.0", "35.0"), 320)):
        t = np.linspace(0.05, 1800, 200001)
        dm = 0.0
        for R in ("1.0", "3.0", "10.0"):
            a = np.array([[f(r["t"]), f(r["T_s"])] for r in hist[(M, z, R, "fc")]])
            b = np.array([[f(r["t"]), f(r["T_s"])] for r in hist[(M, z, R, "nc")]])
            hi = min(a[-1, 0], b[-1, 0])
            tt = t[t <= hi]
            dm = max(dm, float(np.max(np.interp(tt, *a.T) - np.interp(tt, *b.T))))
        check("Q1", f"M{M}/{z} km: max fc-nc during heat-up (K)", f"up to {q}", f"{dm:.1f}", near10(dm) == q)
    X = rows("runs/trajectory/x43a_summary.csv")
    xs = {(r["R_n_mm"], r["catalysis"], r["f_SiC"]): r for r in X}
    traj = np.loadtxt(ROOT / "runs" / "trajectory" / "x43a_mach10.csv", delimiter=",")
    a, b = f(xs[("1.0", "fc", "0.2")]["t_limit"]), f(xs[("1.0", "nc", "0.2")]["t_limit"])
    check("Q1", "X-43A 1 mm 1D: t_limit fc / nc (s)", "75 / 77", f"{a:.2f} / {b:.2f}", down(a) == 75 and down(b) == 77)
    Mx, zx = np.interp(a, traj[:, 0], traj[:, 2]), np.interp(a, traj[:, 0], traj[:, 1]) / 1e3
    check("Q1", "X-43A 1 mm: Mach and altitude at t_limit", "9.9, 31.5 km", f"{Mx:.2f}, {zx:.2f}",
          round(Mx, 1) == 9.9 and round(zx, 1) == 31.5)
    check("Q1", "X-43A 1 mm fc peak (K)", 2820, xs[("1.0", "fc", "0.2")]["T_s_max"][:7],
          near10(f(xs[("1.0", "fc", "0.2")]["T_s_max"])) == 2820)
    a, b = f(xs[("3.0", "fc", "0.2")]["t_limit"]), f(xs[("3.0", "nc", "0.2")]["t_limit"])
    check("Q1", "X-43A 3 mm 1D: t_limit fc / nc (s)", "94 / 108", f"{a:.2f} / {b:.2f}", down(a) == 94 and down(b) == 108)
    for R, q in (("3.0", (2610, 2500)), ("10.0", (2280, 2200))):
        a, b = f(xs[(R, "fc", "0.2")]["T_s_max"]), f(xs[(R, "nc", "0.2")]["T_s_max"])
        check("Q1", f"X-43A {R} mm peak fc / nc (K)", f"{q[0]} / {q[1]}", f"{a:.1f} / {b:.1f}",
              near10(a) == q[0] and near10(b) == q[1])
    check("Q1", "X-43A 10 mm 1D: no failure (fc, nc, switch)", "survives",
          ",".join(xs[("10.0", c, "0.2")]["mode"] for c in ("fc", "nc", "switch")),
          all(xs[("10.0", c, "0.2")]["mode"] == "none" for c in ("fc", "nc", "switch")))
    dT = [f(xs[(R, c, "0.1")]["T_s_max"]) - f(xs[(R, c, "0.3")]["T_s_max"]) for R in ("1.0", "3.0", "10.0")
          for c in ("fc", "nc")]
    check("Q1", "X-43A: peak change 10 -> 30 vol% SiC (K)", "10-60", f"{min(dT):.1f}-{max(dT):.1f}",
          near10(min(dT)) == 10 and near10(max(dT)) == 60)
    dt = [abs(f(xs[(R, c, "0.1")]["t_limit"]) - f(xs[(R, c, "0.3")]["t_limit"])) for R in ("1.0", "3.0") for c in ("fc", "nc")]
    check("Q1", "X-43A: t_limit change 10 -> 30 vol% SiC (s)", "at most 1", f"{max(dt):.2f}", max(dt) <= 1.0)


# ------------------------------------------------------------------------------------------------ Q2
@section
def q2():
    hist, E = histories()
    h = hist[("8.0", "30.0", "1.0", "fc")]

    th = np.array([f(r["t"]) for r in h])

    def at(t0, key):  # interpolated at the quoted time
        return float(np.interp(t0, th, [f(r[key]) for r in h]))
    x8, T8, x16, T16 = at(8.0, "x_B2O3"), at(8.0, "T_s"), at(16.0, "x_B2O3"), at(16.0, "T_s")
    check("Q2", "M8/30/1 mm: x_B2O3 and T at 8 s", "0.69, 1480 K", f"{x8:.3f}, {T8:.0f}", round(x8, 2) == 0.69 and near10(T8) == 1480)
    check("Q2", "M8/30/1 mm: x_B2O3 and T at 16 s", "0.009, 1790 K", f"{x16:.4f}, {T16:.0f}",
          round(x16, 3) == 0.009 and near10(T16) == 1790)
    t_lt = min(f(r["t"]) for r in h if f(r["t"]) > 5 and f(r["x_B2O3"]) < 1e-3)
    check("Q2", "M8/30/1 mm: x_B2O3 < 1e-3 from (s)", 47, f"{t_lt:.2f}", math.ceil(t_lt) == 47)
    g = max(h, key=lambda r: f(r["glass_um"]))
    check("Q2", "M8/30/1 mm: glass maximum (um) and its time", "47 um, about 5 min",
          f"{f(g['glass_um']):.2f} at {f(g['t']):.0f} s", round(f(g["glass_um"])) == 47 and 270 <= f(g["t"]) <= 330)
    A = rows("runs/envelope/set_A.csv")
    env = {(r["mach"], r["alt"], r["R_n"], r["catalysis"], r["f_SiC"]): r for r in A}
    c = env[("8.0", "30000.0", "0.001", "fc", "0.2")]
    check("Q2", "M8/30/1 mm: glass gone (s), T then (K)", "1140, 2340", f"{tf(c):.1f}, {f(c['T_s_at_fail']):.1f}",
          down(tf(c), 10) == 1140 and near10(f(c["T_s_at_fail"])) == 2340)
    check("Q2", "M8/30/1 mm at glass loss: ZrO2 (mm), depleted (um)", "0.33, 94",
          f"{f(c['ZrO2_um_fail']) / 1e3:.3f}, {f(c['dep_um_fail']):.1f}",
          round(f(c["ZrO2_um_fail"]) / 1e3, 2) == 0.33 and round(f(c["dep_um_fail"])) == 94)
    h3 = hist[("9.0", "30.0", "3.0", "fc")]
    gm = max(f(r["glass_um"]) for r in h3)
    tg = f(E[("9.0", "30.0", "3.0", "fc")]["t_glass"])
    check("Q2", "M9/30/3 mm: glass max (um), gone (s)", "28, 132", f"{gm:.2f}, {tg:.2f}", round(gm) == 28 and down(tg) == 132)
    # criterion (a) during the heat-up: critical thickness in the 24 recorded histories (12 conditions x 2 bounds)
    pk, ratio, late, t1, t01, fall = [], [], [], [], [], []
    for k, hh in hist.items():
        ev = E[k]
        tfirst = min([f(ev[c]) for c in ("t_glass", "t_active", "t_limit") if not isnan(ev[c])] or [1800.0])
        dc = np.array([f(r["dcrit_nm"]) for r in hh])
        gl = np.array([f(r["glass_um"]) for r in hh]) * 1e3
        zr = np.array([f(r["ZrO2_um"]) for r in hh])
        t = np.array([f(r["t"]) for r in hh])
        T = np.array([f(r["T_s"]) for r in hh])
        xb = np.array([f(r["x_B2O3"]) for r in hh])
        ok = np.isfinite(dc) & (dc > 0)
        i = int(np.nanargmax(np.where(ok, dc, -1)))
        pk.append((dc[i], T[i], xb[i], gl[i] / 1e3, gl[i] / dc[i]))
        j = int(np.nanargmax(np.where(t <= tfirst, gl, -1)))
        act = ok & (zr >= 1.0) & (t <= tfirst)  # criterion (a) applies once 1 um of composite is converted
        ratio.append(np.min(gl[act] / dc[act]))
        post = act & (t > t[j])
        if post.any():
            late.append(dc[post].max())
        t1.append(t[(gl < 1000.0) & (t < t[j])].max(initial=0.0))
        t01.append(t[(gl < 100.0) & (t < t[j])].max(initial=0.0))
        lost = np.where((t > t[j]) & (gl <= 1.0))[0]  # glass practically gone (<= 1 nm)
        if lost.size:
            above = np.where((t > t[j]) & (t < t[lost[0]]) & (gl >= 1000.0))[0]
            fall.append(t[lost[0]] - t[above[-1]])
    p = np.array(pk)
    check("Q2", "criterion (a): recorded histories", 24, len(pk), len(pk) == 24)
    check("Q2", "peak critical thickness (nm), at T (K), x_B2O3, glass (um), glass/critical",
          "4-54 nm, 1530-1680 K, 0.13-0.68, 0.14-1.2 um, 12-39",
          f"{p[:, 0].min():.1f}-{p[:, 0].max():.1f}, {p[:, 1].min():.0f}-{p[:, 1].max():.0f}, {p[:, 2].min():.2f}-"
          f"{p[:, 2].max():.2f}, {p[:, 3].min():.3f}-{p[:, 3].max():.3f}, {p[:, 4].min():.1f}-{p[:, 4].max():.1f}",
          round(p[:, 0].min()) == 4 and round(p[:, 0].max()) == 54 and near10(p[:, 1].min()) == 1530
          and near10(p[:, 1].max()) == 1680 and round(p[:, 2].min(), 2) == 0.13 and round(p[:, 2].max(), 2) == 0.68
          and round(p[:, 3].min(), 2) == 0.14 and round(p[:, 3].max(), 1) == 1.2 and round(p[:, 4].min()) == 12
          and round(p[:, 4].max()) == 39)
    check("Q2", "glass / critical thickness while criterion (a) applies, minimum", ">= 22", f"{min(ratio):.1f}", min(ratio) >= 22)
    check("Q2", "critical thickness after the glass maximum (nm)", "< 1", f"{len(late)} cases, max {max(late):.2f}",
          max(late) < 1.0)
    check("Q2", "glass falls from 1 um to practically zero within (s)", "2-12", f"{len(fall)}: {min(fall):.1f}-{max(fall):.1f}",
          round(min(fall)) >= 2 and round(max(fall)) <= 12)
    check("Q2", "glass < 1 um for the first (s); < 0.1 um for the first (s)", "6-114; 2-34",
          f"{min(t1):.1f}-{max(t1):.1f}; {min(t01):.1f}-{max(t01):.1f}",
          down(min(t1)) == 6 and down(max(t1)) == 114 and down(min(t01)) == 2 and down(max(t01)) == 34)
    # active oxidation
    B, C = rows("runs/envelope/set_B.csv"), rows("runs/envelope/set_C.csv")
    allr = A + B + C
    check("Q2", "active oxidation first, baseline line (sets A-C)", "never (of 4095)",
          f"{sum(r['mode'] == 'active' for r in allr)} of {len(allr)}",
          sum(r["mode"] == "active" for r in allr) == 0 and len(allr) == 4095)
    po = [f(r["p_O2_kPa"]) for r in allr]
    tt = [f(r["T_transition"]) for r in allr]
    check("Q2", "wall O2 pressure range (kPa)", "2-220", f"{min(po):.2f}-{max(po):.1f}",
          round(min(po)) == 2 and near10(max(po)) == 220)
    check("Q2", "passive-to-active transition range (K)", "2360-3380", f"{min(tt):.1f}-{max(tt):.1f}",
          near10(min(tt)) == 2360 and near10(max(tt)) == 3380)
    low = [r for r in allr if f(r["T_transition"]) < 2422.15]
    cells = sorted({(f(r["mach"]), f(r["alt"]) / 1e3) for r in low})
    margin = min(f(r["T_transition"]) - f(r["T_s_max"]) for r in low)
    check("Q2", "transition below the limit only at", "Mach 5-6, 37-40 km, >= 780 K margin",
          f"Mach {min(m for m, _ in cells)}-{max(m for m, _ in cells)}, {min(z for _, z in cells)}-{max(z for _, z in cells)} km, {margin:.0f} K",
          min(m for m, _ in cells) >= 5 and max(m for m, _ in cells) <= 6 and min(z for _, z in cells) >= 37 and margin >= 780 and down(margin, 10) == 780)
    gap = [-uhtc_data.TR_B / (uhtc_data.TR_A - 0.65 - math.log10(p / 100.0)) - 55850.0 / math.log(1e15 / (p * 1e3)) for p in po]
    check("Model", "Balat-Pichelin line below Eq. (pa) at 2-220 kPa (K)", "280-870", f"{min(gap):.1f}-{max(gap):.1f}",
          near10(min(gap)) == 280 and near10(max(gap)) == 870)
    xs = {(r["R_n_mm"], r["catalysis"], r["f_SiC"]): r for r in rows("runs/trajectory/x43a_summary.csv")}
    check("Q2", "X-43A 1 mm fc: active oxidation at (s)", 94, xs[("1.0", "fc", "0.2")]["t_active"][:7],
          down(f(xs[("1.0", "fc", "0.2")]["t_active"])) == 94)
    check("Q2", "M10/35/1 mm fc: active oxidation after (s)", 25, E[("10.0", "35.0", "1.0", "fc")]["t_active"][:7],
          down(f(E[("10.0", "35.0", "1.0", "fc")]["t_active"])) == 25)


# ------------------------------------------------------------------------------------------------ Q4 (1D OAT, LHS)
def cols2d(rows_):
    """{(R_mm, cat, z_km): {Mach: row}} of 2D (or 1D) map rows"""
    c = defaultdict(dict)
    for r in rows_:
        c[(round(r["R_n"] * 1e3), r["catalysis"], round(r["alt"] / 1e3))][r["mach"]] = r
    return c


def bound2d(d):
    """(highest surviving Mach below the lowest failing one, lowest failing Mach) of one altitude column"""
    fails = sorted(M for M in d if d[M]["mode"] != "none")
    Mf = fails[0] if fails else None
    surv = [M for M in d if d[M]["mode"] == "none" and (Mf is None or M < Mf)]
    return (max(surv) if surv else None), Mf


@section
def q3():
    G, Rf = jl("runs/envelope_2d/grid.jsonl"), jl("runs/envelope_2d/refine.jsonl")
    A1 = [r for r in jl("runs/envelope/set_A.jsonl") if abs(r["f_SiC"] - 0.2) < 1e-9]
    check("Q3", "2D runs: grid, refinement, numerical failures", "1890, 113, none",
          f"{len(G)}, {len(Rf)}, {sum(1 for r in G + Rf if r['error'])}",
          len(G) == 1890 and len(Rf) == 113 and not any(r["error"] for r in G + Rf))
    c2, c1 = cols2d(G + Rf), cols2d(A1)
    Z = range(20, 41)
    B = {(R, c, z): bound2d(c2[(R, c, z)]) for R in (1, 3, 10) for c in ("fc", "nc") for z in Z}
    nonmono = sum(1 for k, (Ms, Mf) in B.items() if Mf is not None
                  and any(c2[k][M]["mode"] == "none" and M > Mf for M in c2[k]))
    res = {(k[0], k[1], k[2]): Mf - Ms for k, (Ms, Mf) in B.items() if Mf is not None}
    check("Q3", "bracketed columns, resolution 0.25, monotone", "113, 0.25",
          f"{len(res)}, {set(res.values())}, non-monotone {nonmono}", len(res) == 113 and set(res.values()) == {0.25}
          and nonmono == 0)
    for R, q in ((1, (8.5, 10.25)), (3, (9.0, 10.75)), (10, (9.25, 11.5))):
        v = [B[(R, "fc", z)][0] for z in Z]
        check("Q3", f"fc boundary {R} mm rises from 20 to 40 km", q, (v[0], v[-1], min(v), max(v)),
              (v[0], v[-1]) == q and min(v) == v[0] and max(v) == v[-1] and all(B[(R, "fc", z)][1] for z in Z))
    v = [(B[(R, "nc", 20)][0] == B[(R, "fc", 20)][0]) for R in (1, 3, 10)]
    check("Q3", "nc = fc boundary at 20 km", "same", v, all(v))
    v = tuple(B[(R, c, 30)][0] for c in ("nc", "fc") for R in (1, 3, 10))
    check("Q3", "30 km: nc 9.75/10.5/11.25 vs fc 9.25/9.75/10.25", (9.75, 10.5, 11.25, 9.25, 9.75, 10.25), v,
          v == (9.75, 10.5, 11.25, 9.25, 9.75, 10.25))
    op = {R: [z for z in Z if B[(R, "nc", z)][1] is None] for R in (1, 3, 10)}
    check("Q3", "nc: nothing fails up to Mach 12", "1 mm 40; 3 mm 37-40; 10 mm 33-40 km", op,
          op == {1: [40], 3: list(range(37, 41)), 10: list(range(33, 41))})
    nc_higher = all(B[(R, "nc", z)][1] is None or B[(R, "nc", z)][0] >= B[(R, "fc", z)][0] for R in (1, 3, 10) for z in Z)
    check("Q3", "nc boundary never below fc", "higher above 20 km", nc_higher, nc_higher)
    for R, q in ((1, ("63", "71")), (3, ("69", "79")), (10, ("76", "86"))):
        sh = tuple(f"{100 * sum(r['mode'] == 'none' for r in G if round(r['R_n'] * 1e3) == R and r['catalysis'] == c) / 315:.0f}"
                   for c in ("fc", "nc"))
        check("Q3", f"survival share {R} mm fc/nc (%)", q, sh, sh == q)
    modes = {R: {z: c2[(R, "fc", z)][B[(R, "fc", z)][1]]["mode"] for z in Z} for R in (1, 3, 10)}
    g = {R: sorted(z for z, m in modes[R].items() if m == "glass") for R in (1, 3, 10)}
    check("Q3", "glass first beyond fc boundary", "1 mm all 21; 3 mm all but 20-21; 10 mm 11 of 13 from 28 km",
          {R: len(g[R]) for R in g}, g[1] == list(Z) and g[3] == list(range(22, 41)) and len(g[10]) == 11
          and min(g[10]) == 28 and all(modes[10][z] == "temperature" for z in Z if z not in g[10]))
    na = sum(1 for r in G + Rf if r["mode"] == "active")
    check("Q3", "SiC active oxidation first in 2D maps", "never", na, na == 0)
    for R, q in ((1, (620, 1790)), (3, (720, 1760)), (10, (780, 1640))):
        t = [c2[(R, "fc", z)][B[(R, "fc", z)][1]]["t_fail"] for z in g[R]]
        check("Q3", f"glass lasts just beyond, {R} mm fc (s)", q, f"{min(t):.1f}-{max(t):.1f}",
              (down2(min(t)), down2(max(t))) == q)
    Tl = [c2[(R, "fc", z)][B[(R, "fc", z)][1]]["T_s_end"] for R in (1, 3, 10) for z in g[R]]
    check("Q3", "T at glass loss just beyond (K)", "2210-2420", f"{min(Tl):.1f}-{max(Tl):.1f}",
          (near10(min(Tl)), near10(max(Tl))) == (2210, 2420))
    for R, q in ((1, 17), (3, 47), (10, 410)):  # upper bounds, rounded up (Supplementary Table S3)
        m12 = [r for r in G if round(r["R_n"] * 1e3) == R and r["catalysis"] == "fc" and r["mach"] == 12.0]
        t = max(tf(r) for r in m12)
        up = math.ceil(t - 1e-9) if t < 200 else math.ceil(t / 10 - 1e-9) * 10
        check("Q3", f"Mach 12 fc, {R} mm: all fail within (s)", q, f"{t:.1f}, n={len(m12)}",
              all(r["mode"] != "none" for r in m12) and up == q and len(m12) == 21)
    for R, q in ((1, (0.5, 2.75)), (3, (1.0, 2.75)), (10, (0.75, 2.0))):
        d = [B[(R, "fc", z)][0] - bound2d(c1[(R, "fc", z)])[0] for z in Z]
        check("Q3", f"2D minus 1D boundary, {R} mm fc (Mach)", q, (min(d), max(d)), (min(d), max(d)) == q)
    g2 = {(round(r["R_n"] * 1e3), r["catalysis"], r["mach"], round(r["alt"] / 1e3)): r for r in G}
    bad = sum(1 for r in A1 if r["mode"] == "none"
              and g2.get((round(r["R_n"] * 1e3), r["catalysis"], r["mach"], round(r["alt"] / 1e3)), {}).get("mode", "none") != "none")
    check("Q3", "conditions surviving in 1D but failing in 2D", "none", bad, bad == 0)
    Tb = {k: c2[k][Ms]["T_s_max"] for k, (Ms, Mf) in B.items() if Mf is not None}
    check("Q3", "every bracketed column: boundary point above 2073 K", "113 of 113",
          f"{sum(v > 2073.15 for v in Tb.values())} of {len(Tb)}", all(v > 2073.15 for v in Tb.values()) and len(Tb) == 113)
    for R, q in ((1, 2310), (3, 2380), (10, 2410)):
        m = max(v for k, v in Tb.items() if k[0] == R and k[1] == "fc")
        check("Q3", f"peak surface at the fc boundary, {R} mm (K)", f"up to {q}", f"{m:.1f}", near10(m) == q)


@section
def appendix_switch():
    J = rows("runs/histories/jump_histories.csv")
    h = defaultdict(list)
    for r in J:
        h[(f(r["mach"]), f(r["alt_km"]), f(r["R_n_mm"]), r["catalysis"])].append(r)
    A1 = jl("runs/envelope/set_A.jsonl")

    def hd(M, z, R):
        return [r["hD_frac"] for r in A1 if r["mach"] == M and r["alt"] == z * 1e3 and abs(r["R_n"] * 1e3 - R) < 1e-9][0]

    def end(k):
        return f(h[k][-1]["T_s"])
    weak = [(8.0, 30.0, 1.0), (8.0, 35.0, 1.0)]
    strong = [(10.0, 35.0, 10.0), (10.5, 38.0, 10.0)]
    d = [abs(end(c + ("switch",)) - end(c + ("nc",))) for c in weak]
    check("AppA", "switch, 1 mm Mach 8: end temperature unchanged (K)", "unchanged", f"{max(d):.2f}", max(d) < 0.5)
    v = [hd(*c) for c in weak] + [hd(*c) for c in strong]
    check("AppA", "h_D/H_s weak / strong cells", "0.06-0.08 / 0.23-0.30", [round(x, 3) for x in v],
          round(min(v[:2]), 2) == 0.06 and round(max(v[:2]), 2) == 0.08 and round(min(v[2:]), 2) == 0.23
          and round(max(v[2:]), 2) == 0.30)
    t = [f(h[c + ("nc",)][0]["t_fail"]) for c in strong]
    mo = [h[c + ("nc",)][0]["mode"] for c in strong]
    check("AppA", "10 mm nc: temperature limit after (s)", "1560, 920", f"{t}, {mo}",
          [down2(x) for x in t] == [1560, 920] and mo == ["temperature", "temperature"])
    dd = [end(c + ("switch",)) - end(c + ("nc",)) for c in strong]
    check("AppA", "switch adds at 1800 s (K)", "about 100, 160", [round(x, 1) for x in dd],
          [near10(x) for x in dd] == [100, 160])
    e = [near10(end(c + (k,))) for k in ("switch", "nc") for c in strong]
    check("AppA", "end temperatures switched / nc (K)", "2560, 2630 / 2460, 2470", e, e == [2560, 2630, 2460, 2470])


@section
def q4():
    o = {r["case_id"]: r for r in jl("runs/sensitivity/oat_B1.jsonl")}
    X = {r["case_id"]: r for r in jl("runs/sensitivity/oat_extra.jsonl")}
    for bn, q, mode in (("B1", 711, "glass"), ("B2", 48, "temperature")):
        b = {r["case_id"]: r for r in jl(f"runs/sensitivity/oat_{bn}.jsonl")}[f"{bn}_base"]
        check("Q4", f"{bn} base time (s), mode", f"{q}, {mode}", f"{b['t_fail']:.2f}, {b['mode']}",
              down(b["t_fail"]) == q and b["mode"] == mode)
    check("Q4", "B1 surface temperature (K)", 2320, f"{o['B1_base']['T_s_max']:.1f}", near10(o["B1_base"]["T_s_max"]) == 2320)
    for cid, q in (("heat_mult_low", 1750), ("heat_mult_high", 370), ("f_SiC_low", 240), ("f_SiC_high", 1390),
                   ("kp_mult_low", 390), ("kp_mult_high", 1540), ("eps_glass_low", 270), ("eps_glass_high", 1060),
                   ("R_n_low", 230), ("R_n_high", None), ("mtf_high", 280), ("mtf_low", None), ("alpha_evap_low", None),
                   ("n_p_high", 380), ("f_SiO_lost_high", 440), ("kd_mult_low", 620), ("kd_mult_high", 840),
                   ("T_limit_low", 105)):
        v = o[f"B1_{cid}"]["t_fail"]
        ok = (v != v) if q is None else down2(v) == q
        check("Q4", f"B1 {cid} (s)", "no failure" if q is None else q, f"{v:.2f}", ok)
    check("Q4", "B1 T_limit_low mode", "temperature", o["B1_T_limit_low"]["mode"], o["B1_T_limit_low"]["mode"] == "temperature")
    base = o["B1_base"]["t_fail"]
    small = ["k_mult_low", "k_mult_high", "cp_mult_low", "cp_mult_high", "k_glass_low", "k_glass_high", "k_ZrO2_low",
             "k_ZrO2_high", "eps_oxide_low", "eps_oxide_high", "T_dep_low", "T_dep_high", "pa_offset_low", "pa_offset_high"]
    mx = max(abs(o[f"B1_{c}"]["t_fail"] / base - 1) for c in small)
    check("Q4", "B1 small effects (k, cp, k_glass, k_ZrO2, eps_oxide, T_dep, offset)", "< 2 %", f"{100 * mx:.2f} %", mx < 0.02)
    cv = abs(o["B1_catalysis_low"]["t_fail"] / base - 1)
    check("Q4", "B1 catalysis bound", "< 1 %", f"{100 * cv:.2f} %", cv < 0.01)
    mx = max(abs(X[f"B1_{c}"]["t_fail"] / base - 1) for c in ("k_dep_low", "k_dep_high", "pa_offset_atomicO"))
    check("Q4", "B1 k_dep 10/60, atomic-O offset", "< 2 %", f"{100 * mx:.3f} %", mx < 0.02)
    mx = max(abs(X[f"B1_{c}"]["t_fail"] / base - 1) for c in ("perm_low", "perm_high"))
    check("Q4", "B1 O2 permeability low / high", "< 0.01 %", f"{100 * mx:.4f} %", mx < 1e-4)
    inc = X["B1_incub_high"]["t_fail"]
    check("Q4", "B1 incubation: time (s), change", "670, -6 %", f"{inc:.1f}, {100 * (inc / base - 1):.1f} %",
          down2(inc) == 670 and round(100 * (inc / base - 1)) == -6)
    v = [X[f"B1_{c}"]["t_fail"] for c in ("kpE_low", "kpE_high", "kphi_low", "kphi_high")]
    check("Q4", "B1 growth law above 2073 K: E 250 / 354, x0.5 / x3 (s)", "550 / 780, 480 / 1530",
          " / ".join(f"{x:.1f}" for x in v), [down2(x) for x in v] == [550, 780, 480, 1530])
    whole = o["B1_kp_mult_high"]["t_fail"]
    check("Q4", "x3 above 2073 K vs x3 everywhere", "within 0.2 %", f"{100 * abs(v[3] / whole - 1):.3f} %",
          abs(v[3] / whole - 1) <= 0.002)
    check("Q4", "B1 runs 250 K beyond the fit end (K)", 250, f"{o['B1_base']['T_s_max'] - 2073.15:.1f}",
          near10(o["B1_base"]["T_s_max"] - 2073.15) == 250)
    check("Q4", "kinetics band factor up (down: 1.7 with the hot-end slope, [FR])", "2", f"{v[3] / base:.2f}",
          round(v[3] / base) == 2)
    b = X["B1_paBalat_high"]
    check("Q4", "B1 Balat line: mode, time (s), T (K)", "active, 59, 2210", f"{b['mode']}, {b['t_fail']:.2f}, {b['T_s_max']:.1f}",
          b["mode"] == "active" and down(b["t_fail"]) == 59 and near10(b["T_s_max"]) == 2210)
    a = [X[f"B1_dfix{d}_high"] for d in ("0.1", "0.5", "1.0")]
    check("Q4", "B1 fixed 0.5 / 1 um: time (s); 0.1 um: change", "24, 24; 0.1 %",
          f"{a[1]['t_fail']:.2f}, {a[2]['t_fail']:.2f}; {100 * abs(a[0]['t_fail'] / base - 1):.2f} %",
          down(a[1]["t_fail"]) == 24 and down(a[2]["t_fail"]) == 24 and round(100 * abs(a[0]["t_fail"] / base - 1), 1) == 0.1)
    o2 = {r["case_id"]: r for r in jl("runs/sensitivity/oat_B2.jsonl")}
    base2 = o2["B2_base"]["t_fail"]
    check("Q4", "B2 T_limit low / high", "29 s / glass 91 s",
          f"{o2['B2_T_limit_low']['t_fail']:.2f} / {o2['B2_T_limit_high']['mode']} {o2['B2_T_limit_high']['t_fail']:.2f}",
          down(o2["B2_T_limit_low"]["t_fail"]) == 29 and o2["B2_T_limit_high"]["mode"] == "glass"
          and down(o2["B2_T_limit_high"]["t_fail"]) == 91)
    for cid, q in ((("R_n_low", "R_n_high"), (35, 73)), (("heat_mult_low", "heat_mult_high"), (66, 39)),
                   (("eps_glass_low", "eps_glass_high"), (41, 52))):
        v = [o2[f"B2_{c}"]["t_fail"] for c in cid]
        check("Q4", f"B2 {cid[0][:-4]} (s)", q, f"{v[0]:.2f}, {v[1]:.2f}", (down(v[0]), down(v[1])) == q)
    check("Q4", "B2 k x0.7 (s)", 53, f"{o2['B2_k_mult_low']['t_fail']:.2f}", down(o2["B2_k_mult_low"]["t_fail"]) == 53)
    b2 = X["B2_paBalat_high"]
    check("Q4", "B2 Balat line: mode, time (s)", "active, 29", f"{b2['mode']}, {b2['t_fail']:.2f}",
          b2["mode"] == "active" and down(b2["t_fail"]) == 29)
    a = [X[f"B2_dfix{d}_high"] for d in ("0.1", "0.5", "1.0")]
    check("Q4", "B2 fixed 1 um: glass after (s); 0.1 / 0.5 um unchanged", "13; yes",
          f"{a[2]['mode']} {a[2]['t_fail']:.2f}; {a[0]['t_fail']:.2f}, {a[1]['t_fail']:.2f}",
          a[2]["mode"] == "glass" and down(a[2]["t_fail"]) == 13
          and all(abs(x["t_fail"] / base2 - 1) < 1e-6 for x in a[:2]))
    kin = [X[f"B2_{c}"]["t_fail"] / base2 - 1 for c in ("kpE_low", "kpE_high", "kphi_low", "kphi_high")]
    check("Q4", "B2 growth law above 2073 K (%)", "-6 to +3", f"{100 * min(kin):.1f} to {100 * max(kin):.1f}",
          round(100 * min(kin)) == -6 and round(100 * max(kin)) == 3)
    big = {"B2_T_limit_low", "B2_T_limit_high", "B2_R_n_low", "B2_R_n_high", "B2_heat_mult_low", "B2_heat_mult_high",
           "B2_eps_glass_low", "B2_eps_glass_high", "B2_k_mult_low", "B2_base"}
    rest = [abs(r["t_fail"] / base2 - 1) for k, r in o2.items() if k not in big]
    rest += [abs(r["t_fail"] / base2 - 1) for k, r in X.items() if k.startswith("B2_")
             and not k.startswith(("B2_paBalat", "B2_dfix1.0"))]
    check("Q4", "B2 every other input (incl. revision inputs)", "< 10 %", f"{100 * max(rest):.1f} %", max(rest) < 0.10)
    cv = abs(o2["B2_catalysis_low"]["t_fail"] / base2 - 1)
    check("Disc", "B2 catalysis bound", "< 4 %", f"{100 * cv:.2f} %", cv < 0.04)
    Lw = {r["case_id"]: r for r in jl("runs/sensitivity/oat_lewis.jsonl")}
    check("Q4", "Lewis factors", "1.19, 1.25", f"{Lw['B1_Le052']['mtf']:.3f}, {Lw['B1_Le23']['mtf']:.3f}",
          round(Lw["B1_Le052"]["mtf"], 2) == 1.19 and round(Lw["B1_Le23"]["mtf"], 2) == 1.25)
    check("Q4", "B1 with the Lewis factors (s)", "540, 510", f"{Lw['B1_Le052']['t_fail']:.1f}, {Lw['B1_Le23']['t_fail']:.1f}",
          down2(Lw["B1_Le052"]["t_fail"]) == 540 and down2(Lw["B1_Le23"]["t_fail"]) == 510)
    lw2 = max(abs(Lw[k]["t_fail"] / base2 - 1) for k in ("B2_Le052", "B2_Le23"))
    check("Q4", "B2 with the Lewis factors", "< 1 %", f"{100 * lw2:.2f} %", lw2 < 0.01)
    S = json.loads((ROOT / "runs" / "sensitivity" / "lhs_stats.json").read_text(encoding="utf-8"))
    L1, L2 = jl("runs/sensitivity/lhs_B1.jsonl"), jl("runs/sensitivity/lhs_B2.jsonl")
    q1 = S["B1"]["tfail_quantiles_s"]
    check("Q4", "B1 LHS median, 5 %, 95 % (s)", "410, 99, >1800", f"{q1[2]:.1f}, {q1[0]:.1f}, {q1[4]:.1f}",
          down2(q1[2]) == 410 and down(q1[0]) == 99 and q1[4] >= 1799)
    m1 = Counter(r["mode"] for r in L1)
    check("Q4", "B1 LHS: survive / glass / temperature", "23 / 142 / 35 of 200", dict(m1),
          (m1["none"], m1["glass"], m1["temperature"]) == (23, 142, 35) and len(L1) == 200)
    pr = S["B1"]["prcc"]
    want = {"emissivity, glass": 0.69, "heating correlation (x)": -0.68, "mass-transfer coefficient (x)": -0.62,
            "evaporation coefficient": -0.57, "oxide growth rate kp (x)": 0.47, "temperature limit (K)": 0.46}
    check("Q4", "B1 PRCC top six", want, {k: round(pr[k], 2) for k in want}, all(round(pr[k], 2) == v for k, v in want.items()))
    q2_ = S["B2"]["tfail_quantiles_s"]
    check("Q4", "B2 LHS median, 5-95 % (s)", "47, 31-112", f"{q2_[2]:.1f}, {q2_[0]:.1f}-{q2_[4]:.1f}",
          down(q2_[2]) == 47 and down(q2_[0]) == 31 and down(q2_[4]) == 112)
    check("Q4", "B2 LHS glass-loss samples", "26 of 200", f"{sum(1 for r in L2 if r['mode'] == 'glass')}",
          sum(1 for r in L2 if r["mode"] == "glass") == 26)
    pr = S["B2"]["prcc"]
    want = {"temperature limit (K)": 0.98, "heating correlation (x)": -0.87, "emissivity, glass": 0.7,
            "catalysis bound": -0.32}
    check("Q4", "B2 PRCC top four", want, {k: round(pr[k], 2) for k in want}, all(round(pr[k], 2) == v for k, v in want.items()))
    check("Q4", "LHS inputs", "16 continuous + catalysis", len(S["B1"]["prcc"]), len(S["B1"]["prcc"]) == 17)
    hdv = [s.h_D / s.H_s for s in (common.stag_state(8.0, 35e3), common.stag_state(8.5, 30e3))]
    check("Disc", "base cases: h_D/H_s (%)", "8-9", f"{100 * hdv[0]:.1f}, {100 * hdv[1]:.1f}",
          round(100 * min(hdv)) >= 8 and round(100 * max(hdv)) <= 9)


# ------------------------------------------------------------------------------------------------ Appendix A (1D maps)
def maps1d(rows_):
    m = defaultdict(list)
    for r in rows_:
        m[(round(r["R_n"] * 1e3), r["catalysis"], r["f_SiC"])].append(r)
    return m


def bound1d(rows_):
    """per altitude (km): highest surviving Mach below the lowest failing one, lowest failing Mach and its mode"""
    by = defaultdict(dict)
    for r in rows_:
        by[r["alt"] / 1e3][r["mach"]] = r
    out = {}
    for z, d in sorted(by.items()):
        fail = [M for M in sorted(d) if d[M]["mode"] != "none"]
        Mf = min(fail) if fail else None
        surv = [M for M in d if d[M]["mode"] == "none" and (Mf is None or M < Mf)]
        out[z] = (max(surv) if surv else None, Mf, d[Mf]["mode"] if Mf else None)
    return out


@section
def appendix_1d():
    A, B, C, D = (jl(f"runs/envelope/set_{s}.jsonl") for s in "ABCD")
    allr = A + B + C + D
    check("AppA", "1D runs, numerical failures", "5040, none", f"{len(allr)}, {sum(1 for r in allr if r['error'])}",
          len(allr) == 5040 and not any(r["error"] for r in allr))
    mA = maps1d(A)
    bd = {k: bound1d(v) for k, v in mA.items()}
    for R, q in ((1, (7.5, 8.0)), (3, (8.0, 8.5)), (10, (8.5, 9.5))):
        b = bd[(R, "fc", 0.2)]
        v = (min(x[0] for x in b.values()), max(x[0] for x in b.values()))
        check("AppA", f"fc boundary {R} mm (Mach)", q, v, v == q and len(b) == 21)
    for R, q in ((1, 0.0), (3, 0.5), (10, 1.0)):
        d = max(bd[(R, "nc", 0.2)][z][0] - bd[(R, "fc", 0.2)][z][0] for z in bd[(R, "fc", 0.2)])
        check("AppA", f"nc shift {R} mm (Mach)", f"up to {q}", d, d == q)
    b10 = bd[(10, "nc", 0.2)]
    v = (min(x[0] for x in b10.values()), max(x[0] for x in b10.values()))
    check("AppA", "10 mm nc boundary (Mach)", (8.5, 10.5), v, v == (8.5, 10.5))
    for R, q in ((1, "41"), (3, "47-50"), (10, "60-63")):
        fr = [100 * sum(r["mode"] == "none" for r in mA[(R, c, 0.2)]) / 315 for c in ("fc", "nc")]
        lo, hi = round(min(fr)), round(max(fr))
        check("AppA", f"survival share {R} mm (%)", q, f"{min(fr):.1f}-{max(fr):.1f}", q == (f"{lo}" if lo == hi else f"{lo}-{hi}"))
    tl = [r["t_fail"] for r in mA[(1, "fc", 0.2)] if r["mode"] == "temperature"]
    check("AppA", "1 mm fc temperature-limit times (s)", "0.4-97", f"{min(tl):.2f}-{max(tl):.2f}",
          round(min(tl), 1) == 0.4 and down(max(tl)) == 97)
    for R, q in ((1, 22), (3, 47), (10, 440)):  # upper bounds, rounded up (Supplementary Table S3)
        m = max(tf(r) for r in mA[(R, "fc", 0.2)] if r["mach"] >= 10.0)
        up = math.ceil(m - 1e-9) if m < 200 else math.ceil(m / 10 - 1e-9) * 10
        check("AppA", f"fc, Mach >= 10, {R} mm: longest time (s)", f"within {q}", f"{m:.2f}", up == q)
    for R, M, zz, n in ((1, (8.0, 8.0), (23, 40), 18), (3, (8.5, 9.0), (28, 40), 12), (10, (9.5, 9.5), (35, 35), 1)):
        g = [r for r in mA[(R, "fc", 0.2)] if r["mode"] == "glass"]
        ms, zs = [r["mach"] for r in g], [r["alt"] / 1e3 for r in g]
        check("AppA", f"glass band {R} mm fc", f"Mach {M}, {zz} km, n={n}",
              f"Mach {min(ms)}-{max(ms)}, {min(zs):.0f}-{max(zs):.0f} km, n={len(g)}",
              (min(ms), max(ms)) == M and (round(min(zs)), round(max(zs))) == zz and len(g) == n)
    g = [r for R in (1, 3, 10) for r in mA[(R, "fc", 0.2)] if r["mode"] == "glass"]
    gT, gt = [r["T_s_at_fail"] for r in g], [r["t_fail"] for r in g]
    check("AppA", "glass band: T at loss (K), glass lasts (s)", "2300-2420, 200-1790",
          f"{min(gT):.1f}-{max(gT):.1f}, {min(gt):.1f}-{max(gt):.1f}",
          near10(min(gT)) == 2300 and near10(max(gT)) == 2420 and down2(min(gt)) == 200 and down2(max(gt)) == 1790)
    for R, q in ((1, 18), (3, 12), (10, 1)):
        n = sum(x[2] == "glass" for x in bd[(R, "fc", 0.2)].values())
        check("AppA", f"altitudes where glass loss sets the fc boundary, {R} mm", q, n, n == q)
    mB = maps1d(B)
    sv = [sum(r["mode"] == "none" for r in (mB[(3, "fc", 0.1)], mA[(3, "fc", 0.2)], mB[(3, "fc", 0.3)])[i]) for i in range(3)]
    check("AppA", "3 mm fc survivors at 10/20/30 vol%", "147, 149, 161", sv, sv == [147, 149, 161])
    gc = [sum(r["mode"] == "glass" for r in (mB[(3, "fc", 0.1)], mA[(3, "fc", 0.2)], mB[(3, "fc", 0.3)])[i]) for i in range(3)]
    check("AppA", "3 mm fc glass-loss failures at 10/20/30 vol%", "15, 12, 1", gc, gc == [15, 12, 1])
    bb = [bound1d(x) for x in (mB[(3, "fc", 0.1)], mA[(3, "fc", 0.2)], mB[(3, "fc", 0.3)])]
    d = max(max(b[z][0] for b in bb) - min(b[z][0] for b in bb) for z in bb[0])
    check("AppA", "SiC fraction moves the 3 mm boundary by (Mach)", "at most 0.5", d, d <= 0.5)
    nc = {(round(r["R_n"] * 1e3), r["mach"], r["alt"]): r for r in A if r["catalysis"] == "nc"}
    same = all((r["t_fail"] == nc[(round(r["R_n"] * 1e3), r["mach"], r["alt"])]["t_fail"])
               or (isnan(r["t_fail"]) and isnan(nc[(round(r["R_n"] * 1e3), r["mach"], r["alt"])]["t_fail"])) for r in C)
    check("AppA", "switch = non-catalytic time to first failure, all conditions", "945, identical", f"{len(C)}, {same}",
          len(C) == 945 and same)
    mD = maps1d(D)
    for R, q, alts, nsurv in ((1, (157, 186), {0.0: 21}, 129), (3, (139, 168), {0.0: 19, -0.5: 2}, 147),
                              (10, (119, 147), {0.0: 6, -0.5: 10, -1.0: 5}, 168)):
        rr = mD[(R, "fc", 0.2)]
        na = sum(r["mode"] == "active" for r in rr)
        nf = sum(r["mode"] != "none" for r in rr)
        bD, bA = bound1d(rr), bd[(R, "fc", 0.2)]
        sh = dict(Counter(bD[z][0] - bA[z][0] for z in bD))
        ns = sum(r["mode"] == "none" for r in rr)
        check("AppA", f"Balat line {R} mm: active first / failing, boundary shifts, survivors",
              f"{q[0]} of {q[1]}, {alts}, {nsurv}", f"{na} of {nf}, {sh}, {ns}", (na, nf) == q and sh == alts and ns == nsurv)
    # the onset is the transition temperature; the end-of-step surface temperature can lie past the limit (line check, 10 Nov 2026)
    act = [f(r["T_transition"]) for r in D if r["mode"] == "active"]
    check("AppA", "Balat line: T at active onset (K)", "2150-2420", f"{min(act):.1f}-{max(act):.1f}",
          near10(min(act)) == 2150 and near10(max(act)) == 2420)


# ------------------------------------------------------------------------------------------------ Q1, 2D (X-43A)
def _q5rows():
    R = rows("runs/q5_2d/q5_results.csv")
    H = rows("runs/q5_2d/q5_histories.csv")
    return R, H


def _hist(H, **kw):
    out = [r for r in H if all(r[k] == v for k, v in kw.items())]
    return np.array([f(r["t"]) for r in out]), np.array([f(r["T_s"]) for r in out])


def _base2d():
    G, Rf = jl("runs/envelope_2d/grid.jsonl"), jl("runs/envelope_2d/refine.jsonl")
    c2 = cols2d(G + Rf)
    return c2, {k: bound2d(v) for k, v in c2.items()}


@section
def q1_2d():
    R, H = _q5rows()
    traj = np.loadtxt(ROOT / "runs" / "trajectory" / "x43a_mach10.csv", delimiter=",")
    t0 = traj[0, 0]
    X = {(r["R_n_mm"], r["catalysis"]): r for r in R if r["kind"] == "x43a"}
    ok = all(X[(Rn, c)]["2d_mode"] == "none" and isnan(X[(Rn, c)]["2d_t_glass"]) for Rn in ("1.0", "3.0", "10.0")
             for c in ("fc", "nc"))
    check("Q1", "X-43A 2D: every radius and bound survives, glass never lost", "survives", ok, ok and len(X) == 6)
    for Rn, q in (("1.0", (2120, 1920)), ("3.0", (1990, 1820)), ("10.0", (1870, 1750))):
        a, b = f(X[(Rn, "fc")]["2d_T_max"]), f(X[(Rn, "nc")]["2d_T_max"])
        check("Q1", f"X-43A 2D {Rn} mm peak fc (nc) (K)", f"{q[0]} ({q[1]})", f"{a:.1f} / {b:.1f}",
              near10(a) == q[0] and near10(b) == q[1])
    pk = [f(X[k]["2d_T_max"]) for k in X]
    lo, hi = 2422 - max(pk), 2422 - min(pk)
    check("Q1", "X-43A 2D: peaks below the limit (K)", "310-670", f"{lo:.1f}-{hi:.1f}",
          near10(lo) == 310 and near10(hi) == 670)
    t, T = _hist(H, kind="x43a", R_n_mm="1.0", catalysis="fc", model="2d")
    i = int(np.argmax(T))
    tp = t[i] + t0
    Mp, zp = np.interp(tp, traj[:, 0], traj[:, 2]), np.interp(tp, traj[:, 0], traj[:, 1]) / 1e3
    check("Q1", "X-43A 2D 1 mm fc: peak at flight time, Mach, altitude", "81 s, Mach 10, 32.5 km",
          f"{tp:.1f} s, {Mp:.2f}, {zp:.2f}", round(tp) == 81 and round(Mp) == 10 and round(zp, 1) == 32.5)
    hot = {}
    for k in X:
        t, T = _hist(H, kind="x43a", R_n_mm=k[0], catalysis=k[1], model="2d")
        hot[k] = float(np.sum(np.diff(t)[T[1:] > 2073.15]))
    others = all(v == 0 for k, v in hot.items() if k != ("1.0", "fc"))
    check("Q1", "X-43A 2D: only 1 mm fc above 2073 K, for (s)", 10, f"{hot[('1.0', 'fc')]:.1f}; others 0: {others}",
          round(hot[("1.0", "fc")]) == 10 and others)
    M, z = traj[:, 2], traj[:, 1]
    zz = z[M > 8] / 1e3
    check("Q1", "X-43A path altitude while faster than Mach 8 (km)", "29-33.5", f"{zz.min():.2f}-{zz.max():.2f}",
          math.floor(zz.min()) == 29 and round(zz.max(), 1) == 33.5)
    c2, B = _base2d()
    for Rn, q in ((1, 12), (3, 3)):
        ms = {zk: B[(Rn, "fc", zk)][0] for zk in range(20, 41)}
        mb = np.array([ms[min(max(int(round(zi / 1e3)), 20), 40)] if 19.5e3 <= zi <= 40.5e3 else np.inf for zi in z])
        over = M > mb
        dur = float(np.sum(np.diff(traj[:, 0])[over[1:]]))
        check("Q1", f"X-43A path beyond the 2D fc boundary, {Rn} mm (s)", q, f"{dur:.1f}", round(dur) == q)
    for Rn in (10,):
        ms = {zk: B[(Rn, "fc", zk)][0] for zk in range(20, 41)}
        mb = np.array([ms[min(max(int(round(zi / 1e3)), 20), 40)] if 19.5e3 <= zi <= 40.5e3 else np.inf for zi in z])
        check("Q1", "X-43A path never beyond the 10 mm fc boundary", "(implied)", bool(np.any(M > mb)), not np.any(M > mb))
    for Rn, q in (("1.0", 2000), ("3.0", 1910)):
        t, T = _hist(H, kind="const", mach="9.0", alt_km="30.0", R_n_mm=Rn, heat_dist="lees", theta_deg="5.0", model="2d")
        check("Q1", f"M9/30 km 2D {Rn} mm after 600 s (K)", q, f"{T[-1]:.1f}", near10(T[-1]) == q and t[-1] >= 599)


# ------------------------------------------------------------------------------------------------ Q3 band and variants
def _variant(stage, name):
    return cols2d(jl(f"runs/envelope_2d/{stage}_{name}.jsonl"))


@section
def q3_band():
    c2, B = _base2d()
    cols = [(R, z) for R in (1, 3, 10) for z in (20, 25, 30, 35, 40)]

    def shifts(stage, name):
        V = _variant(stage, name)
        out = {}
        for R, z in cols:
            Ms, Mf = bound2d(V[(R, "fc", z)])
            out[(R, z)] = (Ms - B[(R, "fc", z)][0], Mf is None, Ms, Mf, V[(R, "fc", z)])
        return out
    S = {n: shifts(st, n) for st, n in (("band", "worst"), ("band", "best"), ("band", "q0_lo"), ("band", "q0_hi"),
                                         ("band", "eps_lo"), ("band", "eps_hi"), ("band", "heat_lo"),
                                         ("band", "heat_hi"), ("glass", "dfail_0.1"), ("glass", "dfail_0.5"),
                                         ("glass", "dfail_1.0"), ("active", "pa_balat"), ("kinetics", "kp_E250"),
                                         ("kinetics", "kp_E354"), ("kinetics", "kp_x0.5"), ("kinetics", "kp_x3"),
                                         ("kinetics", "incubation"), ("kinetics", "perm_low"),
                                         ("kinetics", "perm_high"))}
    rng = lambda n: (min(v[0] for v in S[n].values()), max(v[0] for v in S[n].values()))
    check("Q3b", "worst corner shift (Mach)", "-1.0 to -0.25", rng("worst"), rng("worst") == (-1.0, -0.25))
    check("Q3b", "best corner shift (Mach; one open column)", "+0.25 to +1.0", rng("best"),
          rng("best") == (0.25, 1.0) and sum(v[1] for v in S["best"].values()) == 1)
    check("Q3b", "stagnation heating -15% / +15% shift", "+0.25..+0.5 / -0.5..0", (rng("q0_lo"), rng("q0_hi")),
          rng("q0_lo") == (0.25, 0.5) and rng("q0_hi") == (-0.5, 0.0))
    check("Q3b", "emissivity low / high shift", "-0.5..0 / 0..+0.25", (rng("eps_lo"), rng("eps_hi")),
          rng("eps_lo") == (-0.5, 0.0) and rng("eps_hi") == (0.0, 0.25))
    check("Q3b", "heating beyond 70 deg low / high shift", "at most 0.25", (rng("heat_lo"), rng("heat_hi")),
          max(abs(x) for x in rng("heat_lo") + rng("heat_hi")) == 0.25)
    for R, q20, q40 in ((1, (8.25, 9.0), (9.25, 11.25)), (3, (8.5, 9.25), (10.0, 11.75)), (10, (8.75, 9.5), (10.5, None))):
        w20, b20 = S["worst"][(R, 20)][2], S["best"][(R, 20)][2]
        w40, b40 = S["worst"][(R, 40)][2], S["best"][(R, 40)][2]
        b40v = None if S["best"][(R, 40)][1] else b40
        check("Q3b", f"band {R} mm at 20 and 40 km (worst-best)", f"{q20}, {q40}", f"{(w20, b20)}, {(w40, b40v)}",
              (w20, b20) == q20 and (w40, b40v) == q40)
    # failure-time ratios a quarter Mach beyond the baseline boundary

    def ratios(n, need_glass=True):
        V = S[n]
        out = []
        for R, z in cols:
            Mf = B[(R, "fc", z)][1]
            rb = c2[(R, "fc", z)][Mf]
            rv = V[(R, z)][4].get(Mf)
            if rv is None:
                continue
            if need_glass and not (rb["mode"] == "glass" and rv["mode"] == "glass"):
                continue
            out.append(tf(rv) / tf(rb))
        return out
    r1, r2 = ratios("q0_hi"), ratios("eps_lo")
    check("Q3b", "glass life a quarter beyond: +15% heating / low emissivity", "0.26-0.45 / 0.31-0.50",
          f"{min(r1):.3f}-{max(r1):.3f} / {min(r2):.3f}-{max(r2):.3f}",
          (round(min(r1), 2), round(max(r1), 2), round(min(r2), 2), round(max(r2), 2)) == (0.26, 0.45, 0.31, 0.50))
    nf = [S["q0_lo"][(R, z)][4][B[(R, "fc", z)][1]]["mode"] for R, z in cols if B[(R, "fc", z)][1] in S["q0_lo"][(R, z)][4]]
    check("Q3b", "-15% heating: none of these conditions fails", "none", f"{Counter(nf)}", set(nf) == {"none"} and len(nf) >= 14)
    check("Q3b", "fixed 0.1 um: same boundary in all 15 columns", "same", rng("dfail_0.1"), rng("dfail_0.1") == (0.0, 0.0))
    check("Q3b", "fixed 0.5 / 1 um: largest drop (Mach)", "1.75 / 3.25", (rng("dfail_0.5")[0], rng("dfail_1.0")[0]),
          (rng("dfail_0.5")[0], rng("dfail_1.0")[0]) == (-1.75, -3.25))
    for n, q in (("dfail_0.5", (22, 43)), ("dfail_1.0", (52, 260))):
        t = [tf(v[4][v[3]]) for v in S[n].values() if v[0] < 0]
        check("Q3b", f"{n}: added conditions fail after (s)", q, f"{min(t):.2f}-{max(t):.2f}",
              (down2(min(t)), down2(max(t))) == q)
    bal = S["pa_balat"]
    act = [(R, z) for (R, z), v in bal.items() if v[4][v[3]]["mode"] == "active"]
    check("Q3b", "Balat line: active first beyond the boundary (columns 1/3/10 mm)", "8 (1, 3, 4)",
          f"{len(act)}: {Counter(R for R, _ in act)}", len(act) == 8 and Counter(R for R, _ in act) == {1: 1, 3: 3, 10: 4})
    Ta = [r["T_s_max"] for st, n in (("active", "pa_balat"),) for r in jl(f"runs/envelope_2d/{st}_{n}.jsonl")
          if r["mode"] == "active"]
    check("Q3b", "Balat line: surface temperature at active failure (K)", "2190-2390", f"{min(Ta):.1f}-{max(Ta):.1f}",
          near10(min(Ta)) == 2190 and near10(max(Ta)) == 2390)
    moved = sorted((R, z, v[0]) for (R, z), v in bal.items() if v[0] != 0)
    check("Q3b", "Balat line: boundary lower in five columns by 0.25-0.75", "3 mm 35/40; 10 mm 30/35/40", moved,
          [(R, z) for R, z, _ in moved] == [(3, 35), (3, 40), (10, 30), (10, 35), (10, 40)]
          and min(-s for *_, s in moved) == 0.25 and max(-s for *_, s in moved) == 0.75)
    ta, tb = [], []
    for R, z in cols:
        Mf = B[(R, "fc", z)][1]
        rv = bal[(R, z)][4].get(Mf)
        if rv is not None and rv["mode"] == "active":
            ta.append(tf(rv))
            tb.append(tf(c2[(R, "fc", z)][Mf]))
    check("Q3b", "Balat: failure a quarter beyond, active vs baseline (s)", "76-320 vs 160-1640",
          f"{min(ta):.1f}-{max(ta):.1f} vs {min(tb):.1f}-{max(tb):.1f}",
          (down2(min(ta)), down2(max(ta)), down2(min(tb)), down2(max(tb))) == (76, 320, 160, 1640))
    # growth-law band and other kinetics variants
    for n in ("kp_E250", "kp_E354"):
        check("Q3b", f"{n}: boundary unchanged", 0, rng(n), rng(n) == (0.0, 0.0))
    mv = sorted((n, R, z, v[0]) for n in ("kp_x0.5", "kp_x3") for (R, z), v in S[n].items() if v[0] != 0)
    check("Q3b", "kp x0.5 / x3: boundary moves by 0.25 in 3 columns", "x3: 3 mm 20, 10 mm 40; x0.5: 10 mm 35", mv,
          mv == [("kp_x0.5", 10, 35, 0.25), ("kp_x3", 3, 20, -0.25), ("kp_x3", 10, 40, -0.25)])
    for n in ("incubation", "perm_low", "perm_high"):
        check("Q3b", f"{n}: boundary unchanged", 0, rng(n), rng(n) == (0.0, 0.0))
    for R, q in ((1, (0.90, 1.05)), (3, (0.67, 2.55)), (10, (0.24, 2.97))):
        rr = []
        for n in ("kp_x0.5", "kp_x3"):
            for z in (20, 25, 30, 35, 40):
                Mf = B[(R, "fc", z)][1]
                rv = S[n][(R, z)][4].get(Mf)
                if rv is not None and rv["mode"] != "none":
                    rr.append(tf(rv) / tf(c2[(R, "fc", z)][Mf]))
        check("Q3b", f"kp x0.5/x3: failure-time factor a quarter beyond, {R} mm", q, f"{min(rr):.3f}-{max(rr):.3f}",
              (round(min(rr), 2), round(max(rr), 2)) == q)
    ri, rp = [], []
    for R, z in cols:
        Mf = B[(R, "fc", z)][1]
        tb_ = tf(c2[(R, "fc", z)][Mf])
        rv = S["incubation"][(R, z)][4].get(Mf)
        if rv is not None:
            ri.append(tf(rv) / tb_ - 1)
        for n in ("perm_low", "perm_high"):
            rv = S[n][(R, z)][4].get(Mf)
            if rv is not None:
                rp.append(abs(tf(rv) / tb_ - 1))
    check("Q3b", "incubation: failure times change (%)", "-8 to +39", f"{100 * min(ri):.1f} to {100 * max(ri):.1f}",
          round(100 * min(ri)) == -8 and round(100 * max(ri)) == 39)
    check("Q3b", "permeability band: failure times change (%)", "< 0.1", f"{100 * max(rp):.4f}", max(rp) < 1e-3)


# ------------------------------------------------------------------------------------------------ Q5
@section
def q5():
    R, H = _q5rows()
    base = [r for r in R if r["kind"] == "const" and r["heat_dist"] == "lees" and r["theta_deg"] == "5.0"]
    d = [f(r["1d_T_max"]) - f(r["2d_T_max"]) for r in base]
    check("Q5", "peak drop 1D -> 2D, constant conditions (K)", "390-730", f"{min(d):.1f}-{max(d):.1f}",
          near10(min(d)) == 390 and near10(max(d)) == 730 and len(base) == 10)
    x = [f(r["1d_T_max"]) - f(r["2d_T_max"]) for r in R if r["kind"] == "x43a"]
    check("Q5", "peak drop 1D -> 2D, X-43A, all runs of the comparison (K)", "(not quoted)", f"{min(x):.1f}-{max(x):.1f}", True)
    xx = [f(r["1d_T_max"]) - f(r["2d_T_max"]) for r in R
          if r["kind"] == "x43a" and r["catalysis"] == "fc" and r["R_n_mm"] in ("1.0", "3.0")]
    check("Q5", "peak drop 1D -> 2D along the X-43A path, fc, 1 and 3 mm (K)", "620-710", f"{min(xx):.1f}-{max(xx):.1f}",
          len(xx) == 2 and near10(min(xx)) == 620 and near10(max(xx)) == 710)
    d2, d1 = [], []
    for M, z in (("7.0", "25.0"), ("8.0", "30.0"), ("8.5", "30.0"), ("9.0", "30.0"), ("10.0", "35.0")):
        a = {r["R_n_mm"]: r for r in base if r["mach"] == M and r["alt_km"] == z}
        d2.append(abs(f(a["1.0"]["2d_T_max"]) - f(a["3.0"]["2d_T_max"])))
        d1.append(abs(f(a["1.0"]["1d_T_max"]) - f(a["3.0"]["1d_T_max"])))
    check("Q5", "1 vs 3 mm: 2D < 90 K; 1D 60-170 K", "<90; 60-170", f"{max(d2):.1f}; {min(d1):.1f}-{max(d1):.1f}",
          max(d2) < 90 and near10(min(d1)) == 60 and near10(max(d1)) == 170)
    f1 = [r for r in base if r["1d_mode"] != "none"]
    t1 = [f(r["1d_t_fail"]) for r in f1]
    check("Q5", "1D: 5 of 10 reach the limit within 600 s, after (s); 2D none", "5, 14-55, none",
          f"{len(f1)}, {min(t1):.1f}-{max(t1):.1f}, {sum(r['2d_mode'] != 'none' for r in base)}",
          len(f1) == 5 and all(r["1d_mode"] == "temperature" for r in f1) and (down(min(t1)), down(max(t1))) == (14, 55)
          and all(r["2d_mode"] == "none" for r in base))
    hot, rate = [], []
    for r in base:
        kw = dict(kind="const", mach=r["mach"], alt_km=r["alt_km"], R_n_mm=r["R_n_mm"], heat_dist="lees", theta_deg="5.0")
        t2, T2 = _hist(H, model="2d", **kw)
        t1_, T1 = _hist(H, model="1d", **kw)
        hh = t2[T2 > np.interp(t2, t1_, T1)]
        hot.append(float(hh.max()) if len(hh) else 0.0)
        m = t2 >= 540
        rate.append(float(np.polyfit(t2[m], T2[m], 1)[0] * 60))
    check("Q5", "2D tip hotter only for the first (s)", "1-9", f"{min(hot):.1f}-{max(hot):.1f}",
          round(min(hot)) == 1 and round(max(hot)) == 9)
    check("Q5", "2D still rising at 600 s (K/min)", "0.3-6", f"{min(rate):.2f}-{max(rate):.2f}",
          round(min(rate), 1) == 0.3 and round(max(rate)) == 6)
    dr = {(r["mach"], r["alt_km"], r["R_n_mm"]): r for r in R if r["kind"] == "const" and r["heat_dist"] == "tauber"}
    dd = [f(dr[(r["mach"], r["alt_km"], r["R_n_mm"])]["2d_T_max"]) - f(r["2d_T_max"]) for r in base]
    check("Q5", "sharp-cone correlation: 2D peak higher by up to (K)", 180, f"{max(dd):.1f}", near10(max(dd)) == 180)
    g = dr[("10.0", "35.0", "1.0")]
    check("Q5", "sharp-cone correlation: 1 mm, Mach 10, 35 km loses glass after (s)", 580, g["2d_t_fail"][:7],
          g["2d_mode"] == "glass" and down2(f(g["2d_t_fail"])) == 580)
    t10 = {r["mach"]: r for r in R if r["kind"] == "const" and r["theta_deg"] == "10.0"}
    d10 = [f(t10[M]["2d_T_max"]) - f({(b["mach"]): b for b in base if b["R_n_mm"] == "1.0"}[M]["2d_T_max"])
           for M in ("8.0", "8.5", "9.0")]
    check("Q5", "10 deg half-angle: 2D peak higher (K)", "50-190", f"{min(d10):.1f}-{max(d10):.1f}",
          near10(min(d10)) == 50 and near10(max(d10)) == 190)
    A1 = [r for r in jl("runs/envelope/set_A.jsonl") if abs(r["f_SiC"] - 0.2) < 1e-9]
    c1 = cols2d(A1)
    c2, B = _base2d()
    for Rn, q in ((1, (7.5, 8.0)), (3, (8.0, 8.5)), (10, (8.5, 9.5))):
        v = [bound2d(c1[(Rn, "fc", z)])[0] for z in range(20, 41)]
        check("Q5", f"1D fc boundary range, {Rn} mm", q, (min(v), max(v)), (min(v), max(v)) == q)
    for zk, q in ((20, (0.5, 1.0)), (40, (2.0, 2.75))):
        v = [B[(Rn, "fc", zk)][0] - bound2d(c1[(Rn, "fc", zk)])[0] for Rn in (1, 3, 10)]
        check("Q5", f"gain 1D -> 2D at {zk} km (Mach)", q, (min(v), max(v)), (min(v), max(v)) == q)


# ------------------------------------------------------------------------------------------------ Discussion, Limitations
@section
def disc_limits():
    c2, B = _base2d()
    mx = {}
    for Rn in (1, 3, 10):
        mx[Rn] = max(B[(Rn, "nc", z)][0] - B[(Rn, "fc", z)][0] for z in range(21, 41) if B[(Rn, "nc", z)][1] is not None)
    check("Disc", "nc boundary above fc by up to (Mach), 1/3/10 mm, where nc fails below Mach 12", "1.5/1.25/1.25", mx,
          mx == {1: 1.5, 3: 1.25, 10: 1.25})
    lb = {Rn: [12.0 - B[(Rn, "fc", z)][0] for z in range(21, 41) if B[(Rn, "nc", z)][1] is None] for Rn in (1, 3, 10)}
    opr = {Rn: (min(v), max(v)) if v else None for Rn, v in lb.items()}  # the lower bounds of the gap where nc passes Mach 12
    eq = [(Rn, z) for Rn in (1, 3, 10) for z in range(21, 41)
          if B[(Rn, "nc", z)][1] is not None and B[(Rn, "nc", z)][0] == B[(Rn, "fc", z)][0]]
    lower = [(Rn, z) for Rn in (1, 3, 10) for z in range(20, 41) if B[(Rn, "nc", z)][0] < B[(Rn, "fc", z)][0]]
    check("Disc", "nc passes Mach 12: lower bound of the gap at 1 mm; at 3 mm (Mach); nc never lower", "1.75; 1.25-1.5; never lower",
          f"{opr}; lower {lower}; equal above 20 km {eq}", opr[1] == (1.75, 1.75) and opr[3] == (1.25, 1.5) and not lower)
    L = json.loads((ROOT / "runs" / "final_review" / "analyses.json").read_text(encoding="utf-8"))["laminar"]
    rb = max(x["Re_x"] for x in L if x["at_boundary"])
    r12 = max(x["Re_x"] for x in L)
    check("Lim", "free-stream Re over the wetted length: at the boundary below; up to Mach 12 below", "8e5; 1e6",
          f"{rb:.3g}; {r12:.3g}", rb < 8e5 and r12 < 1e6)
    G = jl("runs/envelope_2d/grid.jsonl")
    A1 = [r for r in jl("runs/envelope/set_A.jsonl") if abs(r["f_SiC"] - 0.2) < 1e-9]
    a1 = {(r["R_n"], r["catalysis"], r["alt"], r["mach"]): r for r in A1}
    ok = all(tf(r) <= tf(a1[(r["R_n"], "fc", r["alt"], r["mach"])]) for r in G
             if r["catalysis"] == "fc" and r["alt"] == 20e3 and r["mach"] >= 10.5)
    check("Disc", "20 km, Mach >= 10.5: 2D reaches the limit as fast as 1D or faster", "yes", ok, ok)
    H = rows("runs/histories/histories.csv")
    v = []
    for Rn in ("1.0", "3.0"):
        for c in ("fc", "nc"):
            hh = [r for r in H if r["mach"] == "10.0" and r["alt_km"] == "35.0" and r["R_n_mm"] == Rn and r["catalysis"] == c]
            r300 = min(hh, key=lambda r: abs(f(r["t"]) - 300))
            v.append(f(r300["ZrO2_um"]) / 1e3)
    check("Disc", "1D histories, Mach 10/35 km: converted depth at 300 s, 1-3 mm (mm)", "0.55-0.95",
          f"{min(v):.3f}-{max(v):.3f}", round(min(v), 2) == 0.55 and round(max(v), 2) == 0.95)


@section
def q3_front():
    F = jl("runs/envelope_2d/front.jsonl")
    base = {(r["R_n"], r["catalysis"], r["alt"], r["mach"]): r for r in jl("runs/envelope_2d/grid.jsonl")
            + jl("runs/envelope_2d/refine.jsonl")}
    same = all(abs(r["T_s_max"] - base[(r["R_n"], r["catalysis"], r["alt"], r["mach"])]["T_s_max"]) < 1e-6
               and r["mode"] == base[(r["R_n"], r["catalysis"], r["alt"], r["mach"])]["mode"]
               and tf(r) == tf(base[(r["R_n"], r["catalysis"], r["alt"], r["mach"])]) for r in F)
    check("Q3f", "front reruns reproduce the stored results (226 = 2 x 113)", "exactly", f"{len(F)}, {same}",
          same and len(F) == 226 and not any(r["error"] for r in F))
    S = [r for r in F if r["case_id"].startswith("front_s")]
    Fl = [r for r in F if r["case_id"].startswith("front_f")]
    d = [r["T_s_max"] - r["T_f_max"] for r in S]
    check("Q3f", "hottest front below hottest surface at the boundary (K)", "160-310", f"{min(d):.1f}-{max(d):.1f}",
          near10(min(d)) == 160 and near10(max(d)) == 310)
    s1 = [r for r in S if round(r["R_n"] * 1e3) == 1]
    f1 = [r for r in Fl if round(r["R_n"] * 1e3) == 1]
    m = max(r["T_f_max"] for r in s1)
    t = max(r["t_T_f_above_fit"] for r in f1)
    check("Q3f", "1 mm: front at the boundary at most (K); failing: above 2073 K for up to (s)", "2060; 103",
          f"{m:.1f}; {t:.1f}", m < 2073.15 and near10(m) == 2060 and round(t) == 103)
    s3 = [r for r in S if round(r["R_n"] * 1e3) == 3]
    h3 = [r for r in s3 if r["T_f_max"] > 2073.15]
    check("Q3f", "3 mm: columns with front > 2073 K at the boundary; longest (s); hottest (K)", "19 of 38; 415; 2160",
          f"{len(h3)} of {len(s3)}; {max(r['t_T_f_above_fit'] for r in h3):.1f}; {max(r['T_f_max'] for r in h3):.1f}",
          len(h3) == 19 and len(s3) == 38 and round(max(r["t_T_f_above_fit"] for r in h3)) == 415
          and near10(max(r["T_f_max"] for r in h3)) == 2160)
    s10 = [r for r in S if round(r["R_n"] * 1e3) == 10]
    tt = [r["t_T_f_above_fit"] for r in s10]
    check("Q3f", "10 mm: front > 2073 K in every column, for (s), up to (K)", "all; 630-1770; 2240",
          f"{sum(r['T_f_max'] > 2073.15 for r in s10)} of {len(s10)}; {min(tt):.1f}-{max(tt):.1f}; "
          f"{max(r['T_f_max'] for r in s10):.1f}",
          all(r["T_f_max"] > 2073.15 for r in s10) and down(min(tt), 10) == 630 and down(max(tt), 10) == 1770
          and near10(max(r["T_f_max"] for r in s10)) == 2240)


@section
def q5_flank_abs():
    S = {(r["mach"], r["variant"]): f(r["T_max"]) for r in rows("runs/q5_2d/q5_sensitivity.csv")}
    R, _ = _q5rows()
    one_d = {r["mach"]: f(r["1d_T_max"]) for r in R if r["kind"] == "const" and r["alt_km"] == "30.0"
             and r["R_n_mm"] == "1.0" and r["heat_dist"] == "lees" and r["theta_deg"] == "5.0"}
    base2 = {r["mach"]: f(r["2d_T_max"]) for r in R if r["kind"] == "const" and r["alt_km"] == "30.0"
             and r["R_n_mm"] == "1.0" and r["heat_dist"] == "lees" and r["theta_deg"] == "5.0"}
    q = {("8.5", "flank heating x0.5"): -100, ("8.5", "flank heating x2"): 150, ("9.0", "flank heating x0.5"): -150,
         ("9.0", "flank heating x2"): 240, ("8.5", "flank 15 mm"): 20, ("8.5", "flank 60 mm"): -10,
         ("9.0", "flank 15 mm"): 130, ("9.0", "flank 60 mm"): -60}
    got = {k: S[k] - S[(k[0], "baseline")] for k in q}
    same_base = all(abs(S[(M, "baseline")] - base2[M]) < 1e-6 for M in ("8.5", "9.0"))
    check("Q5", "flank heating x0.5/x2 and length 15/60 mm (K), baseline = Q5 run", q,
          {k: round(v, 1) for k, v in got.items()}, all(near10(got[k]) == q[k] for k in q) and same_base)
    gap = min(one_d[k[0]] - S[k] for k in S if k[0] in one_d)
    check("Q5", "flank variants: 2D at least below 1D (K)", 470, f"{gap:.1f}", near10(gap) == 470 and gap >= 465)
    check("Lim", "flank heating / length ranges (K)", "-150 to +240; -60 to +130",
          f"{min(got[k] for k in got if 'heating' in k[1]):.1f}..{max(got[k] for k in got if 'heating' in k[1]):.1f}; "
          f"{min(got[k] for k in got if 'mm' in k[1]):.1f}..{max(got[k] for k in got if 'mm' in k[1]):.1f}",
          near10(min(got[k] for k in got if "heating" in k[1])) == -150
          and near10(max(got[k] for k in got if "heating" in k[1])) == 240
          and near10(min(got[k] for k in got if "mm" in k[1])) == -60
          and near10(max(got[k] for k in got if "mm" in k[1])) == 130)
    # Abstract / Conclusions: derived statements
    c2, B = _base2d()
    t = [c2[(Rn, "fc", z)][B[(Rn, "fc", z)][1]]["t_fail"] for Rn in (1, 3, 10) for z in range(20, 41)
         if c2[(Rn, "fc", z)][B[(Rn, "fc", z)][1]]["mode"] == "glass"]
    check("Abs", "glass lost just beyond the boundary after (min)", "10-30", f"{min(t) / 60:.1f}-{max(t) / 60:.1f}",
          min(t) >= 600 and max(t) <= 1800)
    V = cols2d(jl("runs/envelope_2d/band_q0_hi.jsonl"))
    r = [tf(V[(Rn, "fc", z)][B[(Rn, "fc", z)][1]]) / tf(c2[(Rn, "fc", z)][B[(Rn, "fc", z)][1]])
         for Rn in (1, 3, 10) for z in (20, 25, 30, 35, 40)
         if B[(Rn, "fc", z)][1] in V[(Rn, "fc", z)] and c2[(Rn, "fc", z)][B[(Rn, "fc", z)][1]]["mode"] == "glass"
         and V[(Rn, "fc", z)][B[(Rn, "fc", z)][1]]["mode"] == "glass"]
    check("Abs", "+15 % heating, glass loss still first: glass life ratio (8 columns)", "0.26-0.45", f"{len(r)}: {min(r):.3f}-{max(r):.3f}",
          len(r) == 8 and round(min(r), 2) == 0.26 and round(max(r), 2) == 0.45)
    r2 = [tf(V[(Rn, "fc", z)][B[(Rn, "fc", z)][1]]) / tf(c2[(Rn, "fc", z)][B[(Rn, "fc", z)][1]])
          for Rn in (1, 3, 10) for z in (20, 25, 30, 35, 40)
          if B[(Rn, "fc", z)][1] in V[(Rn, "fc", z)] and c2[(Rn, "fc", z)][B[(Rn, "fc", z)][1]]["mode"] == "glass"
          and V[(Rn, "fc", z)][B[(Rn, "fc", z)][1]]["mode"] == "temperature"]
    check("Abs", "+15 % heating, temperature limit takes over: ratio (3 columns)", "0.05-0.28", f"{len(r2)}: {min(r2):.3f}-{max(r2):.3f}",
          len(r2) == 3 and round(min(r2), 2) == 0.05 and round(max(r2), 2) == 0.28)


@section
def q4_front1d():
    H = rows("runs/histories/histories.csv")
    d = [f(r["T_s"]) - f(r["T_front"]) for r in H if 300 <= f(r["t"]) <= 1800]
    check("Q4", "1D histories: front below surface after the first 5 min (K)", "< 35", f"max {max(d):.1f}", max(d) < 35)


# ------------------------------------------------------------------------------------------------ run and write
@section
def rev1_glass_flow():
    """Task 2.5: audit/glass_flow_table.csv against the Limitations text."""
    T = rows("audit/glass_flow_table.csv")
    w = [f(r["w_evap_um_s"]) for r in T]
    check("Lim", "evaporation thinning (um/s)", "0.004-0.16", f"{min(w):.4f}-{max(w):.3f}", round(min(w), 3) == 0.004 and round(max(w), 2) == 0.16)
    Ts = [f(r["T_s"]) for r in T]
    check("Lim", "surface temperature of the estimate (K)", "2060-2420", f"{min(Ts):.1f}-{max(Ts):.1f}", near10(min(Ts)) == 2060 and near10(max(Ts)) == 2420)
    mus = [f(r["mu_F320_Pa_s"]) for r in T] + [f(r["mu_LA_Pa_s"]) for r in T]
    check("Lim", "silica viscosity (Pa s)", "8e3-2e6", f"{min(mus):.3g}-{max(mus):.3g}", 7.5e3 <= min(mus) <= 8.5e3 and 1.9e6 <= max(mus) <= 2.1e6)
    for R_, q in (("1.0", (14, 53)), ("3.0", (23, 92)), ("10.0", (40, 167))):
        d = [f(r["d_equal_LA_um"]) for r in T if r["R_mm"] == R_]
        check("Lim", f"film thickness where flow = evaporation, R = {R_} mm (um)", q, f"{min(d):.1f}-{max(d):.1f}",
              round(min(d)) == q[0] and round(max(d)) == q[1])
    step = [f(b["w_evap_um_s"]) / f(a["w_evap_um_s"]) for a, b in zip(T[::2], T[1::2])]
    step.sort()
    check("Lim", "evaporation growth over one 0.25 Mach step, median", "about 2", f"{step[len(step) // 2]:.2f}", 1.5 <= step[len(step) // 2] <= 2.5)


@section
def rev1_thermal_stress():
    """Task 2.6: audit/thermal_stress_screen.csv against the Limitations text."""
    T = rows("audit/thermal_stress_screen.csv")
    groups = defaultdict(list)
    for r in T:
        groups[(r["mach"], r["alt_km"], r["R_mm"])].append(r)
    check("Lim", "thermal screen: boundary conditions", 5, len(groups), len(groups) == 5)
    dT = [max(f(r["dT_2R"]) for r in g) for g in groups.values()]
    check("Lim", "peak dT within 2 R_n (K)", "210-820", f"{min(dT):.1f}-{max(dT):.1f}", near10(min(dT)) == 210 and near10(max(dT)) == 820)
    a = [max(f(r["ratio_cold_2R"]) for r in g) for g in groups.values()]
    check("Lim", "peak index / strength, TPSX", "1.9-6.1", f"{min(a):.2f}-{max(a):.2f}", round(min(a), 1) == 1.9 and round(max(a), 1) == 6.1)
    b = [max(f(r["ratio_cold_neuman_2R"]) for r in g) for g in groups.values()]
    check("Lim", "peak index / strength, Neuman", "0.8-2.9", f"{min(b):.2f}-{max(b):.2f}", round(min(b), 1) == 0.8 and round(max(b), 1) == 2.9)
    tp = [f(max(g, key=lambda r: f(r["ratio_cold_2R"]))["t_s"]) for g in groups.values()]
    last = {k: max(g, key=lambda r: f(r["t_s"])) for k, g in groups.items()}
    falling = all(f(last[k]["ratio_cold_2R"]) < max(f(r["ratio_cold_2R"]) for r in g) for k, g in groups.items())
    check("Lim", "peak reached within the first minute; index falls after it in every case (s)", "3-30, falling",
          f"{min(tp)}-{max(tp)}, falling: {falling}", min(tp) == 3 and max(tp) <= 60 and falling)


@section
def rev1_audit():
    """Statements corrected after the Phase 1 audit."""
    T = rows("audit/flank_ratio.csv")
    lo, hi = min(f(r["lo_flank"]) for r in T), max(f(r["hi_flank"]) for r in T)
    X = json.loads((ROOT / "audit" / "tools" / "_extra.json").read_text(encoding="utf-8"))["draft_over_lees"]
    v = [x for k in ("5", "8", "10", "12") for x in X[k]]
    check("Heating", "sharp-cone correlation / baseline on the flank, Mach 5-12", "0.5-3.6", f"{min(v):.2f}-{max(v):.2f}",
          round(min(v), 1) == 0.5 and round(max(v), 1) == 3.6)
    tm = max(f(r["T_s"]) for r in rows("runs/histories/histories.csv"))
    check("Lim", "highest wall temperature in the 1D histories (K)", 2960, f"{tm:.1f}", near10(tm) == 2960)


@section
def rev1_oxidant():
    """Task 2.4: edge composition and the oxidant sensitivity case against the Limitations text."""
    T = rows("audit/edge_oxygen.csv")
    hd, xo, sh = ([f(r[k]) for r in T] for k in ("hD_over_Hs", "x_O", "O_atom_share"))
    check("Lim", "h_D/H_s near the boundary", "0.06-0.38", f"{min(hd):.3f}-{max(hd):.3f}", round(min(hd), 2) == 0.06 and round(max(hd), 2) == 0.38)
    check("Lim", "O mole fraction at the edge", "0.01-0.25", f"{min(xo):.4f}-{max(xo):.4f}", round(min(xo), 2) == 0.01 and round(max(xo), 2) == 0.25)
    check("Lim", "share of oxygen atoms free (%)", "2-68", f"{100 * min(sh):.1f}-{100 * max(sh):.1f}", round(100 * min(sh)) == 2 and round(100 * max(sh)) == 68)
    b = {(r["R_n"], r["alt"], r["mach"]): r for s in ("grid", "refine") for r in jl(f"runs/envelope_2d/{s}.jsonl") if r["catalysis"] == "fc"}
    v = [r for r in jl("runs/envelope_2d/oxidant_o_as_o2.jsonl") if (r["R_n"], r["alt"], r["mach"]) in b]
    pr = [r["p_O2_kPa"] / b[(r["R_n"], r["alt"], r["mach"])]["p_O2_kPa"] for r in v]
    check("Lim", "oxidant pressure raised by up to (%)", 40, f"{100 * (max(pr) - 1):.1f}", round(100 * (max(pr) - 1), -1) == 40)
    tr = [r["t_fail"] / b[(r["R_n"], r["alt"], r["mach"])]["t_fail"] for r in v if r["t_fail"] == r["t_fail"] and b[(r["R_n"], r["alt"], r["mach"])]["t_fail"] == b[(r["R_n"], r["alt"], r["mach"])]["t_fail"]]
    check("Lim", "failure later by (%)", "0-9", f"{100 * (min(tr) - 1):.2f}-{100 * (max(tr) - 1):.2f}", min(tr) >= 1.0 - 1e-9 and round(100 * (max(tr) - 1)) == 9)
    md = sum(1 for r in v if r["mode"] != b[(r["R_n"], r["alt"], r["mach"])]["mode"])
    check("Lim", "failure modes changed by the oxidant case", 0, md, md == 0)
    vr = [r for r in rows("runs/analysis/variants_2d.csv") if r["variant"] == "o_as_o2"]
    check("Lim", "boundary columns moved by the oxidant case (of 15)", 0, sum(f(r["shift"]) != 0 for r in vr), len(vr) == 15 and all(f(r["shift"]) == 0 for r in vr))


@section
def rev1_flank():
    """Task 2.7: gain from 1D to 2D for the flank-heating variants, Table tab:gain and its text."""
    B = {(r["R_n_mm"], r["alt_km"]): r for r in rows("runs/analysis/boundaries_2d.csv") if r["catalysis"] == "fc"}
    V = rows("runs/analysis/variants_2d.csv")
    table = {"flank_x0.5": {("1.0", "20.0"): .75, ("1.0", "25.0"): 1.5, ("1.0", "30.0"): 2.0, ("1.0", "35.0"): 2.5, ("1.0", "40.0"): 3.25,
                            ("3.0", "20.0"): 1.0, ("3.0", "25.0"): 1.5, ("3.0", "30.0"): 2.0, ("3.0", "35.0"): 2.5, ("3.0", "40.0"): 2.5},
             "flank_x2": {("1.0", "20.0"): .25, ("1.0", "25.0"): 1.0, ("1.0", "30.0"): 1.5, ("1.0", "35.0"): 2.0, ("1.0", "40.0"): 2.25,
                          ("3.0", "20.0"): .75, ("3.0", "25.0"): 1.0, ("3.0", "30.0"): 1.5, ("3.0", "35.0"): 2.0, ("3.0", "40.0"): 2.0},
             "sharp_cone": {("1.0", "20.0"): .5, ("1.0", "25.0"): 1.25, ("1.0", "30.0"): 1.75, ("1.0", "35.0"): 2.0, ("1.0", "40.0"): 2.5,
                            ("3.0", "20.0"): .75, ("3.0", "25.0"): 1.25, ("3.0", "30.0"): 1.5, ("3.0", "35.0"): 2.0, ("3.0", "40.0"): 2.0}}
    base = {("1.0", "20.0"): .5, ("1.0", "25.0"): 1.25, ("1.0", "30.0"): 1.75, ("1.0", "35.0"): 2.25, ("1.0", "40.0"): 2.75,
            ("3.0", "20.0"): 1.0, ("3.0", "25.0"): 1.25, ("3.0", "30.0"): 1.75, ("3.0", "35.0"): 2.25, ("3.0", "40.0"): 2.25}
    got0 = {k: f(B[k]["M_s_2d"]) - f(B[k]["M_s_1d"]) for k in base}
    check("Q5", "gain table, baseline columns", base, got0, got0 == base)
    allg, sh = [], {}
    for name, q in table.items():
        rv = {(r["R_n_mm"], r["alt_km"]): r for r in V if r["variant"] == name}
        got = {k: f(rv[k]["M_s"]) - f(B[k]["M_s_1d"]) for k in q}
        check("Q5", f"gain table, {name}", q, got, got == q)
        allg += list(got.values())
        sh[name] = [f(rv[k]["shift"]) for k in q]
    check("Q5", "gain range over the 30 columns", "0.25-3.25", f"{min(allg)}-{max(allg)}", min(allg) == .25 and max(allg) == 3.25 and len(allg) == 30)
    check("Q5", "shifts: x0.5 up to 0.5; x2 0.25-0.5 down; sharp cone at most 0.25 down", "see text", sh,
          max(sh["flank_x0.5"]) == .5 and min(sh["flank_x0.5"]) >= 0 and max(sh["flank_x2"]) == -.25 and min(sh["flank_x2"]) == -.5
          and min(sh["sharp_cone"]) == -.25 and max(sh["sharp_cone"]) == 0)


# ------------------------------------------------------------------------------------------------ final review (8-9 Nov 2026)
@section
def final_review():
    """Numbers added in the final review: 1D checks, edge states, 2D transport and active variants, the refitted Eq. 8."""
    K = json.loads((ROOT / "runs" / "final_review" / "checks_1d.json").read_text(encoding="utf-8"))
    b1, b2 = K["B1_base"]["t_fail"], K["B2_base"]["t_fail"]
    check("FR", "B1 / B2 base (s)", "710 / 48", f"{b1:.1f} / {b2:.1f}", down2(b1) == 710 and down2(b2) == 48)
    hs = K["hot_slope"]["t_fail"]
    check("FR", "B1 with E = 138 kJ/mol above 2073 K (s); factor downward", "410; about 1.7", f"{hs:.1f}; {b1 / hs:.2f}",
          down2(hs) == 410 and round(b1 / hs, 1) == 1.7)
    check("FR", "B1 with silica-filled pores (s)", "600", f"{K['pore_silica']['t_fail']:.1f}", down2(K["pore_silica"]["t_fail"]) == 600)
    s1, s2 = K["slab_50mm_B1"]["t_fail"], K["slab_50mm_B2"]["t_fail"]
    check("FR", "50 mm slab: B1 / B2 (s); factors", "970 / 160; 1.4 and 3.3", f"{s1:.1f} / {s2:.1f}; {s1 / b1:.2f}, {s2 / b2:.2f}",
          down2(s1) == 970 and down2(s2) == 160 and round(s1 / b1, 1) == 1.4 and round(s2 / b2, 1) == 3.3)
    check("FR", "B1 at the V2 finest settings: change of failure time", "0.12 %",
          f"{100 * abs(K['numerics_B1']['rel_diff_t_fail']):.3f} %", round(100 * abs(K["numerics_B1"]["rel_diff_t_fail"]), 2) == 0.12)
    ob = K["o2_budget_B1"]
    check("FR", "B1 O2 budget: median demand / permeation; gas supply / demand at least", "1400; 150",
          f"{ob['demand_over_permeation'][1]:.0f}; {ob['gas_supply_over_demand'][0]:.1f}",
          round(ob["demand_over_permeation"][1], -2) == 1400 and down(ob["gas_supply_over_demand"][0], 10) == 150)
    E = json.loads((ROOT / "runs" / "final_review" / "edge_states.json").read_text(encoding="utf-8"))
    bt = E["backtest"].values()
    worst = max(max(abs(v["h_rel"]), abs(v["hD_rel"])) for v in bt)
    fnc = max(abs(v["f_nc_diff"]) for v in bt)
    check("FR", "table extension back-test: h, h_D within; f_nc within", "0.7 %; 0.001", f"{100 * worst:.2f} %; {fnc:.4f}",
          worst <= 0.007 and fnc <= 0.001)
    sa = E["saha"]
    i8, i9 = sa["4.26 kPa, 8000 K"]["ionised_atom_fraction"], sa["4.26 kPa, 8999 K"]["ionised_atom_fraction"]
    check("FR", "Saha, 4.3 kPa: ionised atoms at 8000 / 9000 K", "1 % / 4 %", f"{100 * i8:.2f} % / {100 * i9:.2f} %",
          round(100 * i8) == 1 and round(100 * i9) == 4)
    ed = E["edges"]
    ames30, h1, h2 = ed["V8-1 30 MJ/kg"], ed["V8-2 condition 1"], ed["V8-2 condition 2"]
    fn = [ames30["f_nc"], h1["f_nc"], h2["f_nc"]]
    check("FR", "edge T: Ames 30 MJ/kg; HIT; f_nc of the three", "6250; 6100-6560; 0.25-0.34",
          f"{ames30['T_e']:.0f}; {h2['T_e']:.0f}-{h1['T_e']:.0f}; {min(fn):.3f}-{max(fn):.3f}",
          near10(ames30["T_e"]) == 6250 and near10(h2["T_e"]) == 6100 and near10(h1["T_e"]) == 6560
          and round(min(fn), 2) == 0.25 and round(max(fn), 2) == 0.34)
    held = sorted(k for k, v in ed.items() if v["held_at_T_max"])
    above = sorted(k for k, v in ed.items() if v["T_e"] > 6000.0)
    check("FR", "edge above 6000 K: all V8 but the lowest Ames enthalpy, and V9 1.6; held at 9000 K: Ames 45", "see text",
          f"above {above}; held {held}",
          set(above) == {k for k in ed if k.startswith("V8") and "15 MJ" not in k} | {"V9 1.6 high"}
          and held == ["V8-1 45 MJ/kg", "V8-1 active 45 MJ/kg"])
    # refitted Eq. 8 (study intercepts)
    from uhtc.oxidation import LayerOxidation, OxParams
    a, b = LayerOxidation(OxParams()), LayerOxidation(OxParams(pa_law="jacobson_fe"))
    d = [a.T_transition(p) - b.T_transition(p) for p in (2e3, 220e3)]
    check("FR", "refit: slope (K), scatter (decade), lower by (K) at 2-220 kPa", "21985; 0.23; 210-700",
          f"{-uhtc_data.TR_FE_B:.0f}; {uhtc_data.TR_FE_RMS:.3f}; {d[0]:.0f}-{d[1]:.0f}",
          round(-uhtc_data.TR_FE_B) == 21985 and round(uhtc_data.TR_FE_RMS, 2) == 0.23 and near10(d[0]) == 210
          and round(d[1], -2) == 700)
    # 2D variants (runs/analysis/variants_2d.csv)
    V = rows("runs/analysis/variants_2d.csv")
    def sh(n):
        return [f(r["shift"]) for r in V if r["variant"] == n]
    def cnt(n, sign):
        return sum(1 for x in sh(n) if x == x and (x > 0 if sign > 0 else x < 0))
    got = {n: (cnt(n, -1), cnt(n, 1), min(sh(n)), max(sh(n))) for n in
           ("mtf_x1.25", "mtf_x2", "alpha_0.1", "alpha_0.005", "kZrO2_1.0", "kZrO2_2.5", "pa_fe")}
    check("FR", "2D variants: columns down / up, range", "x1.25: 2 down; x2: 3 down; alpha 0.1: 7 up; alpha 0.005: 9 up, 0-0.5; "
          "k 1.0: 3 down; k 2.5: 11 up; refit: none", got,
          got["mtf_x1.25"][:2] == (2, 0) and got["mtf_x2"][:2] == (3, 0) and got["alpha_0.1"][:2] == (0, 7)
          and got["alpha_0.005"] == (0, 9, 0.0, 0.5) and got["kZrO2_1.0"][:2] == (3, 0) and got["kZrO2_2.5"][:2] == (0, 11)
          and got["pa_fe"] == (0, 0, 0.0, 0.0) and all(abs(v[2]) <= 0.25 and v[3] <= 0.25 for k, v in got.items() if k != "alpha_0.005"))
    kpa = {n: (cnt(n, -1), cnt(n, 1), min(sh(n)), max(sh(n)),
               sorted({r["R_n_mm"] for r in V if r["variant"] == n and f(r["shift"]) != 0}))
           for n in ("kp_all_x0.33", "kp_all_x3")}
    check("FR", "whole law /3: up one step in 8 of 15 at all three radii; x3: down one step in 4", "8 up; 4 down", kpa,
          kpa["kp_all_x0.33"][:4] == (0, 8, 0.0, 0.25) and len(kpa["kp_all_x0.33"][4]) == 3
          and kpa["kp_all_x3"][:4] == (4, 0, -0.25, 0.0))
    off = [r for r in V if r["variant"] == "pa_fe_off0.3"]
    act = sum(1 for r in off if r["mode"] == "active")
    fell = sorted((r["R_n_mm"], r["alt_km"], f(r["shift"])) for r in off if f(r["shift"]) < 0)
    check("FR", "refit at offset 0.3: SiC active first in; falls 0.25-0.5 in three columns, all at 3 or 10 mm; 1 mm unchanged",
          "6; 3 columns", f"{act}; {fell}", act == 6 and len(fell) == 3 and {d[0] for d in fell} <= {"3.0", "10.0"}
          and min(d[2] for d in fell) == -0.5 and max(d[2] for d in fell) == -0.25)
    modes = Counter(r["mode"] for r in V if r["variant"] == "alpha_0.005")
    fe_a = [(r["R_n_mm"], r["alt_km"]) for r in V if r["variant"] == "pa_fe" and r["mode"] == "active"]
    check("FR", "alpha 0.005: first failure beyond the boundary; refit: SiC active first", "temperature in all 15; 10 mm, 40 km",
          f"{dict(modes)}; {fe_a}", modes == Counter({"temperature": 15}) and fe_a == [("10.0", "40.0")])
    # V1 case C at 0.25 ms (runs/final_review/v1_case_c_dt.py)
    O = json.loads((ROOT / "runs" / "final_review" / "v1_case_c_dt.json").read_text(encoding="utf-8"))
    o_s, o_c = O["sphere_dt0.00025"]["order"], O["cylinder_dt0.00025"]["order"]
    check("FR", "V1 2D order at 0.25 ms: sphere; cylinder", "3.8 and 2.8; 2.6 and 2.1",
          f"{o_s[0]:.2f}, {o_s[1]:.2f}; {o_c[0]:.2f}, {o_c[1]:.2f}",
          [round(x, 1) for x in o_s] == [3.8, 2.8] and [round(x, 1) for x in o_c] == [2.6, 2.1])


# ------------------------------------------------------------------------------------------------ external review (9 Nov 2026)
@section
def external_review():
    """Numbers added in response to the external review (runs/final_review/kd_arrhenius.py, external_1d.py,
    external_calcs.py, shoulder_check.py; runs/envelope_2d stage 'external')."""
    FR = ROOT / "runs" / "final_review"
    K = json.loads((FR / "kd_arrhenius.json").read_text(encoding="utf-8"))
    check("EXT", "Arrhenius k_d: E (kJ/mol); k at 2073 K (um2/min)", "343; 348", f"{K['E_kJ_mol']:.1f}; {K['k_ref_um2_min']:.1f}",
          round(K["E_kJ_mol"]) == 343 and round(K["k_ref_um2_min"]) == 348)
    rb = K["ratio_base_to_fit"]
    check("EXT", "baseline / Arrhenius fit at 1650 and 1700 C inside the paper's 2-8", "2-8", f"{rb['1650']:.2f}, {rb['1700']:.2f}",
          all(2.0 <= v <= 8.0 for v in rb.values()))
    E1 = json.loads((FR / "external_1d.json").read_text(encoding="utf-8"))
    b1, ka, kc = E1["B1_base"], E1["B1_kd_arrh"], E1["B1_kd_arrh_cap"]
    keff = lambda r: r["dep_um_end"] ** 2 / (r["t_fail"] / 60.0)  # um^2/min, effective depletion rate up to the failure
    b2r = [E1[f"B2_{v}"]["t_fail"] / E1["B2_base"]["t_fail"] for v in ("kd_arrh", "kd_arrh_cap")]
    check("EXT", "B1 with the Arrhenius k_d: extrapolated; held at 2073 K; base (s); B2 change; effective rate x base",
          "1150; 660; 710; hardly; faster",
          f"{ka['t_fail']:.1f}; {kc['t_fail']:.1f}; {b1['t_fail']:.1f}; {min(b2r):.3f}-{max(b2r):.3f}; {keff(ka) / keff(b1):.1f}",
          down2(ka["t_fail"]) == 1150 and down2(kc["t_fail"]) == 660 and down2(b1["t_fail"]) == 710
          and all(abs(x - 1) < 0.01 for x in b2r) and keff(ka) > 2 * keff(b1) and ka["mode"] == kc["mode"] == "glass")
    gates = [E1[f"{b}_{g}"]["t_fail"] / E1[f"{b}_base"]["t_fail"] for b in ("B1", "B2") for g in ("gate_1nm", "gate_0.5", "gate_2")]
    g0 = [E1[f"{b}_gate_0"] for b in ("B1", "B2")]
    check("EXT", "gate 1 nm to 2 um: same failure times in B1 and B2; no gate: glass loss within 1 s below 550 K", "same; <1 s, <550 K",
          f"ratios {min(gates):.6f}-{max(gates):.6f}; t {[round(x['t_fail'], 2) for x in g0]}, T {[round(x['T_s_max']) for x in g0]}",
          max(abs(x - 1.0) for x in gates) < 1e-9 and all(x["t_fail"] < 1.0 and x["mode"] == "glass" and x["T_s_max"] < 550 for x in g0))
    # glass conductivity 0.7-2.8 W/m/K, one at a time (runs/sensitivity/oat_B1/B2.csv)
    dev = []
    for b in ("B1", "B2"):
        R = {r["case_id"]: r for r in rows(f"runs/sensitivity/oat_{b}.csv")}
        base = E1[f"{b}_base"]["t_fail"]
        dev += [abs(tf(R[f"{b}_k_glass_{s}"]) / base - 1.0) for s in ("low", "high")]
    check("EXT", "glass conductivity 0.7-2.8: failure times of B1 and B2 change by less than", "1 %", f"{100 * max(dev):.2f} %",
          max(dev) < 0.01)
    # V6: zero points counted as failures (the abstract-only Karlsdottir points are no longer in V6, author 10 Nov 2026)
    V6 = [r for r in rows("validation/V6_oxidation/v6_results.csv") if r["set"] == "independent"]
    sc = [r for r in V6 if not isnan(r["ratio"])]
    ok3 = lambda r: 1 / 3 <= f(r["ratio"]) <= 3
    zeros = [r for r in V6 if f(r["measured_um"]) == 0.0]
    check("EXT", "V6 with the zero points counted as outside; the 37 left out the zeros", "24 of 39; 37",
          f"{sum(map(ok3, sc))} of {len(sc) + len(zeros)}; {len(sc)}",
          sum(map(ok3, sc)) == 24 and len(sc) + len(zeros) == 39 and len(sc) == 37 and len(zeros) == 2)
    C = json.loads((FR / "external_calcs.json").read_text(encoding="utf-8"))
    kn = {(k["R_mm"], k["z_km"]): k for k in C["knudsen"]}
    k1 = kn[(1.0, 40.0)]
    check("EXT", "1 mm, 40 km boundary: Kn free stream; Kn at the edge", "0.019; 0.003", f"{k1['Kn_inf']:.4f}; {k1['Kn_e']:.4f}",
          round(k1["Kn_inf"], 3) == 0.019 and round(k1["Kn_e"], 3) == 0.003)
    q1, q3 = kn[(1.0, 40.0)]["q_rare_over_cont"], kn[(3.0, 40.0)]["q_rare_over_cont"]
    check("EXT", "Table 4 footnote: rarefaction heating at 40 km, 1 and 3 mm", "18 %; 12 %", f"{100 * (q1 - 1):.1f} %; {100 * (q3 - 1):.1f} %",
          round(100 * (q1 - 1)) == 18 and round(100 * (q3 - 1)) == 12)
    others = [k["q_rare_over_cont"] for (R, z), k in kn.items() if z == 40.0 and R == 10.0]
    check("EXT", "40 km, 10 mm stays under the 10 % footnote threshold", "< 10 %", f"{100 * (others[0] - 1):.1f} %", others[0] < 1.10)
    qd = C["dynamic_pressure"]
    q20 = [q["q_inf_kPa"] for q in qd if q["z_km"] == 20.0]
    q40 = [q["q_inf_kPa"] for q in qd if q["z_km"] == 40.0]
    xa = C["x43a_above_mach8"]
    check("EXT", "dynamic pressure, boundaries to Mach 12: 20 km; 40 km; X-43A above Mach 8 (kPa)", "280-560; 21-29; 46-65",
          f"{min(q20):.1f}-{max(q20):.1f}; {min(q40):.1f}-{max(q40):.1f}; {xa['q_min_kPa']:.1f}-{xa['q_max_kPa']:.1f}",
          round(min(q20), -1) == 280 and round(max(q20), -1) == 560 and round(min(q40)) == 21 and round(max(q40)) == 29
          and round(xa["q_min_kPa"]) == 46 and round(xa["q_max_kPa"]) == 65)
    # time tolerance (runs/final_review/analyses.json) and 10 mm 2D at Mach 9, 30 km
    A = json.loads((FR / "analyses.json").read_text(encoding="utf-8"))
    tb = [x["t"] for x in A["tolerance_baseline"]]
    check("EXT", "3 columns inside the time tolerance fail after (s)", "1786-1794", f"{len(tb)}: {down(min(tb))}-{down(max(tb))}",
          len(tb) == 3 and down(min(tb)) == 1786 and down(max(tb)) == 1794)
    g = [r for r in rows("runs/envelope_2d/grid.csv") if f(r["R_n"]) == 0.01 and r["catalysis"] == "fc" and f(r["alt"]) == 30e3
         and f(r["mach"]) == 9.0]
    check("EXT", "Table 3 note: 10 mm 2D map passes Mach 9 at 30 km", "passes", g[0]["mode"] if g else "missing",
          len(g) == 1 and g[0]["mode"] == "none")
    # shoulder check (runs/final_review/shoulder_check.json, recomputed from the jsonl by its summary)
    S = json.loads((FR / "shoulder_check.json").read_text(encoding="utf-8"))
    loc = [c for c in S["columns"] if c["mode"] == "local"]
    stg = [c for c in S["columns"] if c["mode"] == "stag"]
    xl = [c for c in S["x43a"] if c["mode"] == "local"]
    check("EXT", "local pressure: no node fails at any highest passing condition (15 columns) nor on the X-43A path (6 runs)",
          "none", f"{sum(c['Ms_any_node_fails'] for c in loc)} of {len(loc)}; {sum(c['any_node_fails'] for c in xl)} of {len(xl)}",
          len(loc) == 15 and not any(c["Ms_any_node_fails"] for c in loc) and len(xl) == 6
          and not any(c["any_node_fails"] for c in xl) and all(c["shift"] == 0 for c in loc))
    sm = {(c["R_mm"], c["z_km"]): c["Mf_stag_t"] for c in stg}
    rat = [c["Mf_stag_t"] / sm[(c["R_mm"], c["z_km"])] for c in loc]
    lag = [c["Mf_off_t"] - c["Mf_stag_t"] for c in loc]
    check("EXT", "lowest failing: stagnation first, sooner than production by (%); neighbours within (s)", "1-4 %; 5 s",
          f"{100 * (1 - max(rat)):.1f}-{100 * (1 - min(rat)):.1f} %; {min(lag):.2f}-{max(lag):.2f} s",
          not any(c["Mf_off_before_stag"] for c in loc) and round(100 * (1 - max(rat))) == 1 and round(100 * (1 - min(rat))) == 4
          and 0 <= min(lag) and max(lag) <= 5.0)
    fl = [c for c in stg if c["Ms_any_node_fails"]]
    check("EXT", "production pressure: flank flags (a) at M_s: columns; time (s); angle (deg)", "13 of 15; 13-48; 79-85",
          f"{len(fl)} of {len(stg)}; {min(c['Ms_first_t'] for c in fl):.1f}-{max(c['Ms_first_t'] for c in fl):.1f}; "
          f"{min(c['Ms_first_phi'] for c in fl):.1f}-{max(c['Ms_first_phi'] for c in fl):.1f}; modes {set(c['Ms_first_mode'] for c in fl)}",
          len(fl) == 13 and len(stg) == 15 and down(min(c["Ms_first_t"] for c in fl)) == 13 and down(max(c["Ms_first_t"] for c in fl)) == 48
          and round(min(c["Ms_first_phi"] for c in fl)) == 79 and round(max(c["Ms_first_phi"] for c in fl)) == 85
          and {c["Ms_first_mode"] for c in fl} == {"glass"})
    SJ = {r["case_id"]: r for r in jl("runs/final_review/shoulder_check.jsonl")}
    inv = []
    for c in fl:  # first flagged node of each column at M_s (production pressure): critical thickness scales with p_O2
        r = SJ[f"stag_R{c['R_mm']:g}_fc_M{c['M_s']:g}_z{c['z_km']:g}"]
        node = min((n for n in r["nodes"][1:] if n["t_fail"] == n["t_fail"] and n["t_fail"] <= 1800), key=lambda n: n["t_fail"])
        inv.append(1.0 / node["pr"])
    check("EXT", "critical thickness at the flagged nodes, production / local = 1/(p/p_s); flank p/p_s under", "21-76; 2 %",
          f"{min(inv):.1f}-{max(inv):.1f}; {100 * max(c['Ms_min_pr'] for c in loc):.2f} %",
          round(min(inv)) == 21 and round(max(inv)) == 76 and max(c["Ms_min_pr"] for c in loc) < 0.02)
    xs = [c for c in S["x43a"] if c["mode"] == "stag" and c["any_node_fails"]]
    check("EXT", "production pressure on the X-43A path: runs with flank flags; time (s)", "3 of 6; 47-52",
          f"{len(xs)} of {sum(c['mode'] == 'stag' for c in S['x43a'])}; {min(c['first_t'] for c in xs):.1f}-{max(c['first_t'] for c in xs):.1f}",
          len(xs) == 3 and down(min(c["first_t"] for c in xs)) == 47 and down(max(c["first_t"] for c in xs)) == 52)
    b = [c["Ms_pa_margins"]["baseline"] for c in loc]
    check("EXT", "criterion (b), Eq. 8: transition lower on the flank by (K); flank margin at least (K)", "600-860; 150",
          f"{min(x['tr_drop_flank'] for x in b):.0f}-{max(x['tr_drop_flank'] for x in b):.0f}; {min(x['off_min'] for x in b):.0f}",
          near10(min(x["tr_drop_flank"] for x in b)) == 600 and near10(max(x["tr_drop_flank"] for x in b)) == 860
          and min(x["off_min"] for x in b) >= 150 and not any(x["any_active"] for x in b))
    alt_ok = all((not c["Ms_pa_margins"][k]["any_active"]) or c["Ms_pa_margins"][k]["stag"] < 0
                 for c in loc for k in ("refit", "refit_off0.3", "balat"))
    check("EXT", "alternative lines: a node turns active only where the stagnation node does", "yes", alt_ok, alt_ok)
    # 2D variants of the external review on the 15 band columns (runs/envelope_2d/external_*.jsonl)
    bnd = {(f(r["R_n_mm"]), f(r["alt_km"])): r for r in rows("runs/analysis/boundaries_2d.csv") if r["catalysis"] == "fc"}
    base = {}
    for p in ("runs/envelope_2d/grid.jsonl", "runs/envelope_2d/refine.jsonl"):
        for r in jl(p):
            if r.get("catalysis") == "fc" and not r.get("error"):
                base[(round(r["R_n"] * 1e3, 3), round(r["alt"] / 1e3, 3), r["mach"])] = r

    def ext_variant(name):
        """(boundary shifts of the 15 columns, {R_mm: failure-time ratios one step beyond the baseline boundary})"""
        V = jl(f"runs/envelope_2d/external_{name}.jsonl")
        hi = defaultdict(lambda: -math.inf)
        for r in V:
            if r["mode"] == "none" and not r.get("error"):
                k = (round(r["R_n"] * 1e3, 3), round(r["alt"] / 1e3, 3))
                hi[k] = max(hi[k], r["mach"])
        shifts, rat = [], defaultdict(list)
        for k, b in bnd.items():
            if k[1] not in (20.0, 25.0, 30.0, 35.0, 40.0):
                continue
            shifts.append(hi[k] - f(b["M_s_2d"]))
            Mf = f(b["M_f_2d"])
            rv = [r for r in V if (round(r["R_n"] * 1e3, 3), round(r["alt"] / 1e3, 3)) == k and abs(r["mach"] - Mf) < 1e-9]
            if rv:
                rat[k[0]].append(tf(rv[0]) / tf(base[k + (Mf,)]))
        return shifts, rat
    sa, ra = ext_variant("kd_arrh")
    sc, rc = ext_variant("kd_arrh_cap")
    allr = [x for d in (ra, rc) for v in d.values() for x in v]
    check("EXT", "Arrhenius k_d (extrapolated, held): boundary moves in no column; failure-time factors one step beyond",
          "0 of 15; 0.80-1.05", f"{sum(x != 0 for x in sa)}, {sum(x != 0 for x in sc)} of {len(sa)}; {min(allr):.3f}-{max(allr):.3f}",
          len(sa) == len(sc) == 15 and not any(sa) and not any(sc) and round(min(allr), 2) == 0.80 and round(max(allr), 2) == 1.05)
    gate_dev = []
    for n in ("gate_0.5", "gate_2"):
        for r in jl(f"runs/envelope_2d/external_{n}.jsonl"):
            b = base.get((round(r["R_n"] * 1e3, 3), round(r["alt"] / 1e3, 3), r["mach"]))
            if b is not None:
                gate_dev.append(abs(tf(r) / tf(b) - 1.0) if r["mode"] == b["mode"] else math.inf)
    sg = [ext_variant(n)[0] for n in ("gate_0.5", "gate_2")]
    check("EXT", "2D gates 0.5 and 2 um: no boundary moves; every failure time as in the baseline run at the same condition",
          "none; same", f"shifts {[sum(x != 0 for x in s) for s in sg]}; {len(gate_dev)} runs, max |ratio - 1| {max(gate_dev):.2e}",
          all(len(s) == 15 and not any(s) for s in sg) and len(gate_dev) == 90 and max(gate_dev) < 1e-9)
    check("EXT", "1 mm: both versions give practically the same factors", "0.80-0.90",
          f"{[round(x, 4) for x in sorted(ra[1.0])]} / {[round(x, 4) for x in sorted(rc[1.0])]}",
          [round(x, 3) for x in sorted(ra[1.0])] == [round(x, 3) for x in sorted(rc[1.0])]
          and round(min(ra[1.0]), 2) == 0.80 and round(max(ra[1.0]), 2) == 0.90)
    # stress screen with the abstract values of Neuman et al.: time of the peak within 2 R_n (audit/thermal_stress_screen.csv)
    grp = defaultdict(list)
    for r in rows("audit/thermal_stress_screen.csv"):
        grp[(r["mach"], r["alt_km"], r["R_mm"])].append((f(r["ratio_cold_neuman_2R"]), f(r["t_s"])))
    tp = [max(v)[1] for v in grp.values()]
    check("EXT", "stress index with the Neuman values peaks after (s)", "0.1-3", f"{min(tp)}-{max(tp)}", min(tp) == 0.1 and max(tp) == 3.0)


if __name__ == "__main__":
    only = sys.argv[1:]
    for fn in SECTIONS:
        if only and fn.__name__ not in only:
            continue
        try:
            fn()
        except Exception as e:  # noqa: BLE001 - a crash is a failed check, the other sections still run
            check(fn.__name__, "section ran", "yes", f"{type(e).__name__}: {e}", False)
            traceback.print_exc()
    out = ["# Re-trace of the numbers in paper/main.tex", "",
           "Generated by `runs/analysis/trace_numbers.py` from the saved outputs. Rounding conventions of the paper: failure",
           "times rounded down, temperatures to the nearest 10 K.", "",
           "| Section | Claim | Quoted | Recomputed | OK |", "|---|---|---|---|---|"]
    for s, c, q, v, ok in RESULTS:
        out.append(f"| {s} | {c} | {q} | {v} | {'yes' if ok else '**NO**'} |")
    nf = sum(1 for r in RESULTS if not r[-1])
    out += ["", f"{len(RESULTS)} checks, {nf} failed."]
    if not only:
        (ROOT / "runs" / "analysis" / "trace_numbers.md").write_text("\n".join(out) + "\n", encoding="utf-8")
    for s, c, q, v, ok in RESULTS:
        if not ok:
            print("FAIL", s, c, "| quoted", q, "| got", v)
    print(f"{len(RESULTS)} checks, {nf} failed")
    sys.exit(nf)
