"""V9: glass-loss clock against the VKI Plasmatron tests of Marschall et al. Plan, inputs and criterion: validation/V9_plan.md,
committed before the first run. No model parameter is adjusted.

Model: 1D, production settings, measured cold-wall flux in two stages (as V8). Slab 7 mm, adiabatic back, ZrB2-30SiC, both catalysis bounds.
Output: report.md, v9_results.json
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
from uhtc import coupled as cp, materials as mt  # noqa: E402
from uhtc.oxidation import LayerOxidation, OxParams  # noqa: E402

F_SIC, DEPTH, T_END, P_STAT = 0.3, 0.007, 480.0, 1.0e4
# sample: t_switch, (q_cw W/m2, H J/kg, P_dyn Pa) low, high; low plateau K; high temperature K at t_ref; t_ref s; measured jump onset s
CASES = {
    "1.8": (358.0, (109e4, 12.0e6, 30.0), (185e4, 17.7e6, 60.0), 1895.0, 2110.0, 480.0, None),
    "1.5": (357.0, (109e4, 12.0e6, 30.0), (202e4, 18.6e6, 75.0), 1825.0, 2210.0, 425.0, 425.0),
    "1.6": (359.0, (109e4, 12.0e6, 30.0), (236e4, 21.7e6, 90.0), 1860.0, 2250.0, 380.0, 380.0),
}


def run(case, cat):
    ts, lo, hi = case[:3]
    steps = [(0.0, lo[1], P_STAT + lo[2], lo[0]), (ts, hi[1], P_STAT + hi[2], hi[0])]
    ox = LayerOxidation(OxParams(f_SiC=F_SIC))
    cond = cp.PrescribedConditions(steps, R_n=0.01, catalysis=cat)
    r = cp.run_1d(cond, mt.composite(F_SIC), ox, T_END, depth=DEPTH, record_every=1.0)
    L = [ox.layers(s) for s in r.states]
    return r, [l["glass"] * 1e6 for l in L]


def main():
    out = {}
    lines = ["# V9: glass-loss clock against the Plasmatron tests of Marschall et al.", "",
             "Plan, inputs and criterion: `validation/V9_plan.md`, committed before this run. 1D model, no adjusted parameter.", "",
             "| Bound | Sample | Low plateau (K) [measured] | Peak glass (µm) | Glass at switch (µm) | Glass loss after switch (s) [measured onset] | "
             "T at reference time (K) [measured] | Outcome | Timing | Temperature |", "|---|---|---|---|---|---|---|---|---|---|"]
    verdict = {}
    for cat in ("fc", "nc"):
        ok_all = True
        for name, c in CASES.items():
            ts, _, _, T_lo, T_hi, t_ref, t_jump = c
            r, g = run(c, cat)
            t = np.array(r.t)
            T = np.array(r.T_s)
            T_low = float(np.interp(ts - 1.0, t, T))
            T_ref = float(np.interp(t_ref, t, T))
            tg = float(r.events["t_glass_loss"])
            lost = not math.isnan(tg)
            after = tg - ts if lost else math.nan
            if t_jump is None:
                out_ok = not lost
                tim_ok = True
                meas = "none by 480 s"
            else:
                out_ok = lost and tg > ts
                on = t_jump - ts
                tim_ok = out_ok and on / 3.0 <= after <= on * 3.0
                meas = f"{on:.0f}"
            tem_ok = abs(T_low - T_lo) <= 80.0 and abs(T_ref - T_hi) <= 80.0
            ok = out_ok and tim_ok and tem_ok
            ok_all &= ok
            out[f"{cat}_{name}"] = dict(T_low=T_low, T_ref=T_ref, t_glass=tg, after_switch=after, t_limit=float(r.events["t_limit"]), t_active=float(r.events["t_active"]),
                                        glass_peak_um=float(max(g)), glass_switch_um=float(np.interp(ts, t, g)), T_max=float(T.max()), outcome=out_ok, timing=tim_ok,
                                        temperature=tem_ok, t=t.tolist(), T=T.tolist(), glass_um=list(g))
            yn = lambda b: "pass" if b else "fail"  # noqa: E731
            aft = f"{after:.0f}" if lost else "no loss by 480 s"
            lines.append(f"| {cat} | {name} | {T_low:.0f} [{T_lo:.0f}] | {max(g):.1f} | {np.interp(ts, t, g):.1f} | {aft} [{meas}] | {T_ref:.0f} [{T_hi:.0f}] | "
                         f"{yn(out_ok)} | {yn(tim_ok)} | {yn(tem_ok)} |")
        verdict[cat] = ok_all
    res = any(verdict.values())
    lines += ["", f"V9: **{'PASS' if res else 'FAIL'}** (all three samples must pass with the same catalysis bound; fully catalytic {verdict['fc']}, non-catalytic {verdict['nc']}).", ""]
    (HERE / "v9_results.json").write_text(json.dumps(out), encoding="utf-8")
    (HERE / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
