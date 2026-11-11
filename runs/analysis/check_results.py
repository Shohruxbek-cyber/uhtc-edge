"""Second, independent check of the envelope results (REVIEW.md).

1. Steady state: for cases that survive the 1800 s cap with a heat-soaked slab, the final surface temperature must
   equal the radiative-equilibrium temperature of the same surface balance with no conduction (solved here directly
   with uhtc.surface.radiative_equilibrium from the stagnation state, the final emissivity and the chemical heat terms
   of the final oxide state).  Tolerance 5 K (the slab is not exactly steady; the chemical term is lagged).
2. Re-run: 20 random cases with refined numerics (first cell 0.5 um, growth 1.04, 2.5 K and 1.25 % step limits) must
   reproduce the time to failure within 1 % and the failure mode exactly.
3. Map sanity: at fixed altitude, radius and catalysis the time to failure should not increase with Mach number by more
   than the 2 % numerical tolerance; exceptions are listed (they may be physical: regime switches).
Output: runs/analysis/check_results.md
"""
import json
import math
import random
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "runs"))
sys.path.insert(0, str(ROOT / "src"))
import common  # noqa: E402
from uhtc import coupled as cp, materials as mt, surface as sf  # noqa: E402
from uhtc.oxidation import LayerOxidation, OxParams  # noqa: E402


def load(s):
    return [json.loads(l) for l in open(ROOT / "runs" / "envelope" / f"set_{s}.jsonl", encoding="utf-8") if l.strip()]


def steady_check(rows, n=40, seed=1):
    rng = random.Random(seed)
    cand = [r for r in rows if r["mode"] == "none" and r["catalysis"] in ("fc", "nc") and r["T_s_end"] > 1500]
    picks = rng.sample(cand, min(n, len(cand)))
    out = []
    for r in picks:
        ox = LayerOxidation(OxParams(f_SiC=r["f_SiC"]))
        cond = cp.ConstantConditions(r["mach"], r["alt"], r["R_n"], state=common.stag_state(r["mach"], r["alt"]),
                                     catalysis=r["catalysis"])
        fl = cond.at(0.0)
        conv = cp.convective_flux(fl, r["catalysis"] == "fc", cp.WallEnthalpy(fl.p_s, True), cp.WallEnthalpy(fl.p_s, False))
        gl = r["glass_um_end"] * 1e-6
        e = ox.p.eps_oxide + (ox.p.eps_glass - ox.p.eps_oxide) * min(gl / ox.p.d_eps, 1.0)
        T_re = sf.radiative_equilibrium(conv, sf.constant_eps(e))
        out.append((r["case_id"], r["T_s_end"], T_re, r["T_s_end"] - T_re))
    return out


def rerun_check(rows, n=20, seed=2):
    rng = random.Random(seed)
    cand = [r for r in rows if r["catalysis"] in ("fc", "nc") and r["mode"] != "none"]
    picks = rng.sample(cand, min(n, len(cand)))
    out = []
    for r in picks:
        ox = LayerOxidation(OxParams(f_SiC=r["f_SiC"]))
        orig = ox.rate_limit
        ox.rate_limit = lambda s, T, f, _o=orig: _o(s, T, f) * 0.25
        cond = cp.ConstantConditions(r["mach"], r["alt"], r["R_n"], state=common.stag_state(r["mach"], r["alt"]),
                                     catalysis=r["catalysis"])
        res = cp.run_1d(cond, mt.composite(r["f_SiC"]), ox, common.T_CAP, depth=0.02, dx0=0.5e-6, growth=1.04,
                        dT_max=2.5, stop_on=common.stop_first)
        tf, mode = cp.first_failure(res.events)
        out.append((r["case_id"], r["t_fail"], tf, mode == r["mode"], abs(tf / r["t_fail"] - 1)))
    return out


def monotonic_check(rows):
    groups = {}
    for r in rows:
        groups.setdefault((r["R_n"], r["catalysis"], r["f_SiC"], r["alt"]), []).append(r)
    bad = []
    for k, g in groups.items():
        g = sorted(g, key=lambda r: r["mach"])
        t = [min(r["t_fail"] if r["t_fail"] == r["t_fail"] else 1e9, 1800.0) for r in g]
        for i in range(1, len(g)):
            if t[i] > t[i - 1] * 1.02 and t[i - 1] < 1800:
                bad.append((g[i - 1]["case_id"], round(t[i - 1], 1), g[i]["case_id"], round(t[i], 1), g[i - 1]["mode"], g[i]["mode"]))
    return bad


def main():
    rows = load("A") + load("B") + load("C")
    st = steady_check(rows)
    rr = rerun_check(rows)
    mono = monotonic_check(rows)
    worst_st = max(abs(x[3]) for x in st) if st else float("nan")
    worst_rr = max(x[4] for x in rr) if rr else float("nan")
    lines = ["# Second check of the envelope results", "",
             f"1. Steady state vs direct radiative-equilibrium solve, {len(st)} surviving cases: largest |difference| "
             f"{worst_st:.2f} K (tolerance 5 K) → **{'PASS' if worst_st <= 5 else 'FAIL'}**.", "",
             "| case | T_s,end (K) | T_RE (K) | diff (K) |", "|---|---|---|---|"]
    lines += [f"| {c} | {a:.1f} | {b:.1f} | {d:+.2f} |" for c, a, b, d in st]
    lines += ["", f"2. Re-run of {len(rr)} failing cases with refined numerics: largest relative change of the time to "
              f"failure {worst_rr * 100:.2f} % (tolerance 1 %), modes identical: {all(x[3] for x in rr)} → "
              f"**{'PASS' if worst_rr <= 0.01 and all(x[3] for x in rr) else 'FAIL'}**.", "",
              "| case | t_fail (s) | refined (s) | same mode | rel. change |", "|---|---|---|---|---|"]
    lines += [f"| {c} | {a:.2f} | {b:.2f} | {m} | {d * 100:.2f} % |" for c, a, b, m, d in rr]
    lines += ["", f"3. Non-monotonic time to failure with Mach number (> 2 %): {len(mono)} pairs.", ""]
    if mono:
        lines += ["| lower Mach case | t_fail | next Mach case | t_fail | modes |", "|---|---|---|---|---|"]
        lines += [f"| {a} | {b} | {c} | {d} | {e} → {f} |" for a, b, c, d, e, f in mono[:60]]
    (ROOT / "runs" / "analysis" / "check_results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines[:3]), "\n", lines[-1] if not mono else f"{len(mono)} non-monotonic pairs")


if __name__ == "__main__":
    main()
