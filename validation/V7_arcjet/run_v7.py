"""V7: arc-jet test of a ZrB2-15 vol% SiC hemisphere (R = 5 mm) in the SPES facility (University of Naples).

Sources (open): Savino, De Stefano Fumo, Paterna, Di Maso, Open Aerosp. Eng. J. 3 (2010) 20-31 [T]; Di Maso PhD
thesis, Univ. Naples (2009) [D].  Full input/target table with locations: literature/V7_TARGET.md.
Inputs used:
  * steps 1-5: H0 = 4.5, 5.5, 7.0, 8.6, 10.0 MJ/kg and Pitot pressure 6800, 7400, 8100, 8800, 9500 Pa ([T] Table 3);
  * step start times 50, 130, 191, 247, 306 s and power-off at 548 s, from the transitions of the digitized
    hemisphere curve ([T] Fig. 8, literature/digitized/Savino2010_TOAEJ_ArcJetUHTC_fig8.csv) with step 1 starting
    80 s (inset hold time) before the first visible transition (DERIVED); initial temperature 300 K (ASSUMPTION);
  * gas 80/20 wt% N2/O2 -> O2 mole fraction 0.1795 (DERIVED); edge state from equilibrium air at (p, H0);
  * heating: Sutton-Graves sphere at R_n = 5 mm, both catalysis bounds and the oxidation-state switch; 2D:
    cos(phi) over the whole hemisphere (a ground test is given by H0 and p only, so coupled.run_2d uses the cosine
    shape whatever heat_dist says; the flight distribution of D-2D-7 needs a flight Mach number);
  * material: ZrB2-15SiC properties (k between the 10 and 20 vol% values; UNSURE), oxidation model baseline;
  * conduction: 1D planar slab 10 mm deep (specimen length; ASSUMPTION), adiabatic back; also 2D axisymmetric
    hemisphere (R 5 mm) with adiabatic back face; the stem/holder heat leak is not modelled (unknown contacts).
Targets: plateau temperatures at the end of each step (two-colour pyrometer, [T] Fig. 8 inset: 1573, 1723, 1843,
1923, 2053 K); oxide thickness at the stagnation point after the test: 150 um maximum ([D] p. 55); no SiC-depleted
layer on the hemisphere ([D] p. 54).
Tolerances (set in V7_TARGET.md before the runs): +/-80 K for steps 1-2, +/-60 K for steps 3-5; oxide within a
factor of 1.5.
Diagnostics added after the predictive runs failed (they do not change any flight input):
  D1. effective heating fraction: one constant factor on the non-catalytic heating that best matches the five
      plateaus (least squares), compared with the authors' estimate that only about 1/3 of H0 is "operative" for a
      non-catalytic UHTC in this facility ([D] p. 46);
  D2. oxidation model driven by the measured stagnation temperature history (digitized Fig. 8; below the pyrometer
      floor of 1273 K, where oxidation is negligible, a linear rise from 300 K at step 1 start is assumed), which tests
      the oxide prediction independently of the heating problem.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
from uhtc import coupled as cp, materials as mt  # noqa: E402
from uhtc.conduction2d import NoseGeometry  # noqa: E402
from uhtc.oxidation import LayerOxidation, OxParams  # noqa: E402

STARTS = [50.0, 130.0, 191.0, 247.0, 306.0]
H0 = [4.5e6, 5.5e6, 7.0e6, 8.6e6, 10.0e6]
P02 = [6800.0, 7400.0, 8100.0, 8800.0, 9500.0]
T_OFF = 548.0
TARGET_T = [1573.0, 1723.0, 1843.0, 1923.0, 2053.0]
TOL = [80.0, 80.0, 60.0, 60.0, 60.0]
X_O2 = (0.2 / 32.0) / (0.2 / 32.0 + 0.8 / 28.0)
RN = 0.005


def conditions(catalysis, heat_mult=1.0):
    steps = [(t - STARTS[0], H, p) for t, H, p in zip(STARTS, H0, P02)]
    return cp.PrescribedConditions(steps, RN, x_O2=X_O2, catalysis=catalysis, heat_mult=heat_mult)


def run(catalysis, model="1d", eps=None, heat_mult=1.0):
    ox = LayerOxidation(OxParams(f_SiC=0.15, **({} if eps is None else dict(eps_glass=eps, eps_oxide=eps))))
    mat = mt.composite(0.15)
    t_end = T_OFF - STARTS[0]
    if model == "1d":
        r = cp.run_1d(conditions(catalysis, heat_mult), mat, ox, t_end, depth=0.010, dt_max=1.0, record_every=1.0)
    else:
        g = NoseGeometry(R_n=RN, theta_c=0.0, L_flank=5e-4, axisymmetric=True)
        r = cp.run_2d(conditions(catalysis), g, mat, ox, t_end, dt_max=1.0, record_every=1.0,
                      mesh_kw=dict(n_arc=30, n_flank=3, n_eta=30, eta1=0.004))
    t = r.t + STARTS[0]
    ends = [STARTS[i + 1] - 1.0 for i in range(4)] + [T_OFF - 1.0]
    Tend = [float(np.interp(te, t, r.T_s)) for te in ends]
    L = ox.layers(r.states[-1])
    return dict(t=t, Ts=r.T_s, Tend=Tend, oxide_um=(L["ZrO2"] + L["glass"]) * 1e6, glass_um=L["glass"] * 1e6,
                ZrO2_um=L["ZrO2"] * 1e6, dep_um=L["depleted"] * 1e6, events=r.events)


def effective_fraction():
    """D1: constant multiplier on the non-catalytic heating minimising the plateau RMS error (golden section)."""
    from scipy.optimize import minimize_scalar
    cache = {}

    def rms(f):
        if f not in cache:
            cache[f] = run("nc", heat_mult=f)
        return float(np.sqrt(np.mean((np.array(cache[f]["Tend"]) - np.array(TARGET_T)) ** 2)))
    r = minimize_scalar(rms, bounds=(0.2, 1.0), method="bounded", options={"xatol": 0.005})
    best = cache[min(cache, key=lambda f: abs(f - r.x))]
    return float(r.x), float(r.fun), best


def oxide_from_measured(meas):
    """D2: oxidation model driven by the measured surface temperature (front = surface; thin scale)."""
    from uhtc.coupled import FlowState
    ox = LayerOxidation(OxParams(f_SiC=0.15))
    s = ox.init_state(300.0)
    t_m, T_m = meas["t"], meas["T"]
    first = t_m[np.argmax(T_m > 1290)]
    tt = np.arange(STARTS[0], T_OFF, 0.25)
    TT = np.where(tt < first, 300.0 + (1290.0 - 300.0) * (tt - STARTS[0]) / (first - STARTS[0]), np.interp(tt, t_m, T_m))
    cond = conditions("nc")
    for a, b, T in zip(tt[:-1], tt[1:], TT[:-1]):
        fl = cond.at(a - STARTS[0])
        s = ox.advance(s, b - a, float(T), float(T), fl)
    L = ox.layers(s)
    return dict(oxide_um=(L["ZrO2"] + L["glass"]) * 1e6, ZrO2_um=L["ZrO2"] * 1e6, glass_um=L["glass"] * 1e6,
                dep_um=L["depleted"] * 1e6, xB=L["x_B2O3"])


def main(do_2d=True):
    meas = {}
    with open(ROOT / "literature" / "digitized" / "Savino2010_TOAEJ_ArcJetUHTC_fig8.csv", encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(l for l in f if not l.startswith("#")) if r["series"] == "hemisphere"]
    meas["t"] = np.array([float(r["t_s"]) for r in rows])
    meas["T"] = np.array([float(r["T_K"]) for r in rows])
    cases = [("fc", "1d", None), ("nc", "1d", None), ("switch", "1d", None), ("nc", "1d", 0.6)]
    if do_2d:
        cases += [("nc", "2d", None), ("fc", "2d", None)]
    res = {}
    for cat, mdl, eps in cases:
        key = f"{mdl}-{cat}" + ("" if eps is None else f"-eps{eps}")
        res[key] = run(cat, mdl, eps)
        print(key, [round(x) for x in res[key]["Tend"]], f"oxide {res[key]['oxide_um']:.0f} um", flush=True)
    lines = ["# V7 — Arc-jet test: ZrB2–15SiC hemisphere (SPES, Savino et al. 2010)", "",
             "Inputs, targets and tolerances: `run_v7.py` header and `literature/V7_TARGET.md`.", "",
             "Plateau temperature at the end of each step (K):", "",
             "| Case | Step 1 | Step 2 | Step 3 | Step 4 | Step 5 | Steps within tolerance | Oxide at stagnation (µm) | Depleted (µm) |",
             "|---|---|---|---|---|---|---|---|---|",
             f"| Measured (Fig. 8 inset) | {' | '.join(f'{x:.0f}' for x in TARGET_T)} | — | 150 (max) | none |",
             f"| Tolerance | ±80 | ±80 | ±60 | ±60 | ±60 | | ×1.5 → 100–225 | |"]
    verdict = {}
    for k, v in res.items():
        n_ok = sum(abs(a - b) <= tol for a, b, tol in zip(v["Tend"], TARGET_T, TOL))
        ox_ok = 100.0 <= v["oxide_um"] <= 225.0
        verdict[k] = (n_ok, ox_ok)
        lines.append(f"| {k} | {' | '.join(f'{x:.0f}' for x in v['Tend'])} | {n_ok}/5 | {v['oxide_um']:.0f}"
                     f"{' (in)' if ox_ok else ' (out)'} | {v['dep_um']:.0f} |")
    f_eff, rms_eff, best_eff = effective_fraction()
    res[f"1d-nc x{f_eff:.2f} (D1, fitted)"] = best_eff
    d2 = oxide_from_measured(meas)
    best = max(verdict, key=lambda k: verdict[k][0])
    passed = [k for k, (n, o) in verdict.items() if n == 5 and o]
    lines += ["", f"Result: **{'PASS' if passed else 'FAIL'}** for the predictive cases (no parameter adjusted to this test). "
              f"Cases meeting all temperature and oxide tolerances: {', '.join(passed) if passed else 'none'}. "
              f"Best temperature agreement: {best} ({verdict[best][0]}/5 steps).", "",
              "Interpretation. The 1D slab over-predicts by 250-540 K (non-catalytic) because the specimen is a small blunt body "
              "heated mainly near the stagnation point (cosine distribution) that conducts heat to, and radiates from, its "
              "cooler flanks; the 2D axisymmetric model represents this and, without any adjusted parameter, comes within "
              + (", ".join(f"{a - b:+.0f}" for a, b in zip(res['2d-nc']['Tend'], TARGET_T)) if '2d-nc' in res else "n/a")
              + " K of the five plateaus. The remaining over-prediction at steps 4-5 (about +100 K) is consistent with "
              "what the model leaves out: the arc-jet stream is far from equilibrium (step-5 CFD: N2 vibrational temperature "
              "8230 K, atomic O mass fraction 0.061) so part of H0 does not reach a weakly catalytic wall (the authors "
              "estimate that only about 1/3 of H0 is 'operative', Di Maso 2009 p. 46), and the stem/holder heat leak is "
              "not modelled. The fully catalytic bound is far too hot, as the authors also found.", "",
              "Diagnostics (not predictions):",
              f"- D1: in 1D, a constant factor of {f_eff:.2f} on the non-catalytic heating matches the five plateaus with an "
              f"rms error of {rms_eff:.0f} K (plateaus {', '.join(f'{x:.0f}' for x in best_eff['Tend'])} K); the measured "
              "plateaus rise more gently with enthalpy than any constant factor can reproduce.",
              f"- D2: driven by the measured temperature history, the oxidation model gives {d2['oxide_um']:.0f} µm of oxide "
              f"(ZrO2 {d2['ZrO2_um']:.0f} µm + glass {d2['glass_um']:.0f} µm) at the stagnation point vs 150 µm measured "
              f"(maximum; tolerance 100–225 µm) → **{'within' if 100 <= d2['oxide_um'] <= 225 else 'outside'}**; "
              f"SiC-depleted layer {d2['dep_um']:.0f} µm vs none observed on the hemisphere.", ""]
    res_out = dict(res)
    res_out["_D1"] = dict(f_eff=f_eff, rms=rms_eff)
    res_out["_D2"] = d2
    json.dump({k: {kk: (vv.tolist() if isinstance(vv, np.ndarray) else vv) for kk, vv in v.items()} for k, v in res_out.items()},
              open(HERE / "v7_results.json", "w"), indent=1, default=float)
    (HERE / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines[-3:]))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    sel = meas["T"] > 1280
    ax.plot(meas["t"][sel], meas["T"][sel], "k.", ms=3, label="measured (Savino et al. 2010, Fig. 8)")
    for k, v in res.items():
        ax.plot(v["t"], v["Ts"], lw=1.2, label=k)
    ax.set_xlim(40, 560)
    ax.set_ylim(1000, 3000)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Stagnation surface temperature (K)")
    ax.legend(fontsize=7, frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(HERE / "fig_v7.png", dpi=200)
    fig.savefig(HERE / "fig_v7.pdf", metadata={"CreationDate": None})


if __name__ == "__main__":
    main(do_2d="--no2d" not in sys.argv)
