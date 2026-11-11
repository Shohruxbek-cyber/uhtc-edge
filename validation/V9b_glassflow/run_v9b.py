"""V9b: the plasmatron check V9 with glass flow at the stagnation point (plan: validation/V9b_plan.md, committed before this run).

The V9 model and inputs exactly (validation/V9_plasmatron/run_v9.py), plus the stagnation-point flow sink
w_flow = tau' delta^2 / mu + 2 rho_e beta^2 delta^3 / (3 mu), tau' = Phi Pr^(2/3) beta, applied after every oxide step by a
wrapper of the oxide model (src/uhtc is not changed). beta from the Fay-Riddell inversion of the stage's measured cold-wall flux.
Usage: python validation/V9b_glassflow/run_v9b.py   -> v9b_results.json, report.md
"""
import json
import math
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "validation" / "V9_plasmatron"))
from uhtc import coupled as cp, materials as mt, thermo as th, heating as ht, atmosphere as atm  # noqa: E402
from uhtc.glassflow import PR, viscosity  # noqa: E402
from uhtc.oxidation import LayerOxidation, OxParams  # noqa: E402
import run_v9  # noqa: E402

V9 = json.loads((ROOT / "validation/V9_plasmatron/v9_results.json").read_text())
H300 = float(th.h_air_frozen(300.0))


@lru_cache(maxsize=64)
def edge_beta(H, p, q):
    """(rho_e, beta) of a stage: equilibrium edge at (H, p), beta from the Fay-Riddell inversion of q (calorimeter wall 350 K)."""
    with th.extended_tables(cp.T_EDGE_MAX):
        Te = th.T_air_equilibrium(H, p, T_hi=cp.T_EDGE_MAX - 1.0)
        M, _, hD = th.mixture_props(th.equilibrium_air(Te, p), Te)
    rho_e, mu_e = p * M / (th.R * Te), float(atm.viscosity(Te))
    Tw = 350.0
    unit = ht.fay_riddell(p * 0.028965 / (th.R * Tw), float(atm.viscosity(Tw)), rho_e, mu_e, H, float(th.h_air_frozen(Tw)), hD, 1.0)
    return rho_e, (q / unit) ** 2


class FlowSinkOxidation(LayerOxidation):
    """LayerOxidation with the stagnation-point film-flow sink applied after each step (V9b only)."""

    def __init__(self, p, law="LA", beta_mult=1.0):
        super().__init__(p)
        self.law, self.beta_mult = law, beta_mult
        self.max_rate = 0.0  # largest sink rate seen, m/s
        self.removed = 0.0  # film removed by the sink, m

    def advance(self, s, dt_, T_s, T_f, flow):
        b = s
        s = super().advance(s, dt_, T_s, T_f, flow)
        q = flow.Phi * (flow.H_s - H300)  # the measured cold-wall flux of the stage (PrescribedConditions)
        rho_e, beta = edge_beta(round(flow.H_s, 3), round(flow.p_s, 6), round(q, 6))
        beta *= self.beta_mult
        mu = float(viscosity(T_s, self.law))
        A = flow.Phi * PR ** (2.0 / 3.0) * beta / mu
        B = 2.0 * rho_e * beta ** 2 / (3.0 * mu)
        d0 = self.glass_thickness(s)
        if d0 <= 0.0:
            return s
        d, n = d0, 20  # d(delta)/dt = -(A delta^2 + B delta^3), RK4 over the step (the sink is many orders below the step limit)
        h = dt_ / n
        f = lambda x: -(A * x * x + B * x ** 3)  # noqa: E731
        for _ in range(n):
            k1 = f(d)
            k2 = f(max(d + 0.5 * h * k1, 0.0))
            k3 = f(max(d + 0.5 * h * k2, 0.0))
            k4 = f(max(d + h * k3, 0.0))
            d = max(d + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0, 0.0)
        self.max_rate = max(self.max_rate, -f(d0))
        self.removed += d0 - d
        fac = d / d0
        s["n_B"] *= fac
        s["n_S"] *= fac
        dg, dc = self.glass_thickness(s), s["_dcrit"]
        if dg <= dc and s["x_r"] - self.p.x_r0 > self.p.x_ox_min and math.isnan(s["t_glass_loss"]):
            s["t_glass_loss"] = s["t"]  # criterion (a) after the sink (end of step; the sink is negligible, see report)
        s["_dg_prev"] = dg
        return s


def run(case, cat, law="LA", beta_mult=1.0):
    ts, lo, hi = case[:3]
    steps = [(0.0, lo[1], run_v9.P_STAT + lo[2], lo[0]), (ts, hi[1], run_v9.P_STAT + hi[2], hi[0])]
    ox = FlowSinkOxidation(OxParams(f_SiC=run_v9.F_SIC), law, beta_mult)
    cond = cp.PrescribedConditions(steps, R_n=0.01, catalysis=cat)
    r = cp.run_1d(cond, mt.composite(run_v9.F_SIC), ox, run_v9.T_END, depth=run_v9.DEPTH, record_every=1.0)
    return r, ox


def score(name, case, r):
    ts, _, _, T_lo, T_hi, t_ref, t_jump = case
    t, T = np.array(r.t), np.array(r.T_s)
    T_low, T_ref = float(np.interp(ts - 1.0, t, T)), float(np.interp(t_ref, t, T))
    tg = float(r.events["t_glass_loss"])
    lost = not math.isnan(tg)
    after = tg - ts if lost else math.nan
    if t_jump is None:
        out_ok, tim_ok = not lost, True
    else:
        out_ok = lost and tg > ts
        meas = t_jump - ts
        tim_ok = out_ok and meas / 3.0 <= after <= 3.0 * meas
    tmp_ok = abs(T_low - T_lo) <= 80.0 and abs(T_ref - T_hi) <= 80.0
    return dict(T_low=T_low, T_ref=T_ref, t_glass=tg, after=after, outcome=out_ok, timing=tim_ok, temperature=tmp_ok,
                ok=out_ok and tim_ok and tmp_ok)


def main():
    res, lines = {}, ["# V9b: the plasmatron check V9 with glass flow at the stagnation point", "",
                      "Plan: `validation/V9b_plan.md` (committed before this run). V9 model and inputs, plus the stagnation-point "
                      "film-flow sink (LA viscosity, beta from the Fay-Riddell inversion of the measured flux). V9 values in brackets.", "",
                      "| Bound | Sample | Low plateau (K) | T at reference time (K) | Glass loss after switch (s) | Largest sink rate (um/s) | "
                      "Film removed by flow (um) | Outcome | Timing | Temperature |", "|---|---|---|---|---|---|---|---|---|---|"]
    verdict = {}
    for cat in ("fc", "nc"):
        ok_all = True
        for name, c in run_v9.CASES.items():
            r, ox = run(c, cat)
            sc = score(name, c, r)
            ok_all &= sc["ok"]
            v9 = V9[f"{cat}_{name}"]
            tg9 = v9["t_glass"]
            after9 = tg9 - c[0] if tg9 == tg9 else math.nan
            sc.update(sink_max_um_s=ox.max_rate * 1e6, removed_um=ox.removed * 1e6, v9_T_low=v9["T_low"], v9_T_ref=v9["T_ref"],
                      v9_after=after9)
            res[f"{cat}_{name}"] = sc
            fa = lambda x: "none by 480 s" if x != x else f"{x:.1f}"  # noqa: E731
            lines.append(f"| {cat} | {name} | {sc['T_low']:.1f} [{v9['T_low']:.1f}] | {sc['T_ref']:.1f} [{v9['T_ref']:.1f}] | "
                         f"{fa(sc['after'])} [{fa(after9)}] | {sc['sink_max_um_s']:.2e} | {sc['removed_um']:.2e} | "
                         f"{'yes' if sc['outcome'] else 'no'} | {'yes' if sc['timing'] else 'no'} | {'yes' if sc['temperature'] else 'no'} |")
        verdict[cat] = ok_all
    passed = any(verdict.values())
    # sensitivity (not scored): beta x0.5, x2 and the F320 viscosity, fully catalytic bound
    sens = {}
    for tag, kw in (("beta_x0.5", dict(beta_mult=0.5)), ("beta_x2", dict(beta_mult=2.0)), ("F320", dict(law="F320"))):
        for name, c in run_v9.CASES.items():
            r, ox = run(c, "fc", **kw)
            sc = score(name, c, r)
            sens[f"{tag}_{name}"] = dict(after=sc["after"], T_ref=sc["T_ref"], removed_um=ox.removed * 1e6)
    out = dict(passed=passed, by_bound=verdict, results=res, sensitivity=sens)
    (HERE / "v9b_results.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    lines += ["", f"**V9b {'PASSES' if passed else 'FAILS'}** (fully catalytic {'meets' if verdict['fc'] else 'misses'}, "
              f"non-catalytic {'meets' if verdict['nc'] else 'misses'} the criterion in all three samples).", "",
              "Sensitivity, fully catalytic, not scored (glass loss after the switch in s; film removed by flow in um):", ""]
    for k, v in sens.items():
        a = "none by 480 s" if v["after"] != v["after"] else f"{v['after']:.1f}"
        lines.append(f"- {k}: {a}; removed {v['removed_um']:.2e}")
    (HERE / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
