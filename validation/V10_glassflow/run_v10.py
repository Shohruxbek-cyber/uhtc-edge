"""V10: verification of the glass-flow module against exact solutions (plan: validation/V10_plan.md, committed before this run).

Usage: python validation/V10_glassflow/run_v10.py   -> v10_results.json, report.md
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
from uhtc.glassflow import FlowParams, Surface, film_step  # noqa: E402

R, L, N, MU, TOL = 1e-3, 10e-3, 2000, 1e4, 0.01


def cylinder():
    s = np.linspace(0.0, L, N + 1)
    ds = L / N
    w = np.full(N + 1, R * ds)
    w[0] = w[-1] = R * ds / 2
    return Surface(s, np.full(N + 1, R), w, np.zeros(N + 1), np.full(N + 1, R), R, True)


def run_to(surf, d, tau, G, t_end, params, source=0.0, dt=None):
    """Advance the film to t_end; returns delta and the mass-audit error. dt: oxide-step length (default: whole interval)."""
    mu = np.full(len(d), MU)
    tau_f, G_f = np.full(len(d), tau), np.full(len(d), G)
    G_f[-1] = 0.0
    film0 = float(np.sum(surf.w * d))
    t, shed, src = 0.0, 0.0, 0.0
    step = dt or t_end
    while t < t_end * (1 - 1e-12):
        h = min(step, t_end - t)
        if source:
            d = d + source * h
            src += source * h * float(np.sum(surf.w))
        d, sh, _ = film_step(d, surf, tau_f, G_f, mu, h, params)
        shed += sh
        t += h
    audit = abs(float(np.sum(surf.w * d)) + shed - film0 - src) / max(film0 + src, 1e-300)
    return d, audit


def main():
    surf = cylinder()
    s_star = L / 4
    i_star = int(round(s_star / (L / N)))
    out, rows = {}, []
    expl = FlowParams()
    d0 = 20e-6
    # check 1: shear only
    tau, G = 1000.0, 0.0
    t_head = MU * s_star / (tau * d0)
    d, t_prev, audits = np.full(N + 1, d0), 0.0, []
    for k in (2, 4, 8):
        d, a = run_to(surf, d, tau, G, k * t_head - t_prev, expl)
        t_prev = k * t_head
        audits.append(a)
        exact = d0 / k
        rows.append(("1 shear", f"t = {k} t_head", d[i_star], exact))
    # check 1b: pressure gradient only
    tau, G = 0.0, 1e6
    t_head = MU * s_star / (G * d0 ** 2)
    d, t_prev = np.full(N + 1, d0), 0.0
    for k in (2, 4, 8):
        d, a = run_to(surf, d, tau, G, k * t_head - t_prev, expl)
        t_prev = k * t_head
        audits.append(a)
        exact = d0 * math.sqrt(1.0 / k)
        rows.append(("1b pressure", f"t = {k} t_head", d[i_star], exact))
    # check 2: steady state with a uniform source, explicit and forced implicit
    tau, S = 1000.0, 1e-8
    d_end = math.sqrt(2 * MU * S * L / tau)
    transit = L / (tau * d_end / (2 * MU))
    for name, params, dt in (("2 steady, explicit", expl, 10.0), ("2 steady, implicit", FlowParams(max_explicit=0), 10.0)):
        d, a = run_to(surf, np.zeros(N + 1), tau, 0.0, 20 * transit, params, source=S, dt=dt)
        audits.append(a)
        for frac in (0.25, 0.5, 1.0):
            i = int(round(frac * N))
            rows.append((name, f"s = {frac:g} L", d[i], math.sqrt(2 * MU * S * surf.s[i] / tau)))
    res = [dict(check=c, point=p, value=float(v), exact=float(e), rel_err=float(abs(v - e) / e), ok=bool(abs(v - e) <= TOL * e))
           for c, p, v, e in rows]
    passed = all(r["ok"] for r in res)
    out = dict(passed=passed, n=len(res), n_ok=sum(r["ok"] for r in res), max_rel_err=max(r["rel_err"] for r in res),
               audit_max=max(audits), results=res, grid=dict(R=R, L=L, N=N, mu=MU), tolerance=TOL)
    (HERE / "v10_results.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    L_ = ["# V10: verification of the glass-flow module", "",
          "Plan: `validation/V10_plan.md` (committed before this run). Criterion: every value within 1 % of the exact solution.", "",
          "| Check | Point | Computed (um) | Exact (um) | Relative error | Within 1 % |", "|---|---|---|---|---|---|"]
    L_ += [f"| {r['check']} | {r['point']} | {r['value'] * 1e6:.4f} | {r['exact'] * 1e6:.4f} | {r['rel_err']:.2e} | {'yes' if r['ok'] else 'NO'} |"
           for r in res]
    L_ += ["", f"**V10 {'PASSES' if passed else 'FAILS'}: {out['n_ok']} of {out['n']} values within 1 %.** Largest relative error "
           f"{out['max_rel_err']:.2e}. Largest mass-audit error of the runs {out['audit_max']:.1e}."]
    (HERE / "report.md").write_text("\n".join(L_) + "\n", encoding="utf-8")
    print("\n".join(L_))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
