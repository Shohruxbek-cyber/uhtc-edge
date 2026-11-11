"""V8: arc-jet / plasma tests of ZrB2-SiC with a MEASURED cold-wall heat flux (revision, Part 3.2).

Targets and inputs: literature/V8_TARGET.md (open sources only; no open case has both a measured heat flux and a
surface-temperature history, so two complementary cases are used).
  V8-1  NASA Ames AHF arc jet, ZrB2-20 vol% SiC flat-face coupon (US Patent 5,750,450, Bull, White & Kaufman 1998):
         calorimeter on a Cu model of identical flat-face geometry, 3.0 MW/m^2 at 0.042 atm for 5 min ("passive") and
         4.0 MW/m^2 at 0.058 atm for 3 min ("active"); coupon 6.35 mm thick, Fiberform behind (adiabatic back face).
         Total enthalpy NOT given: bracket 15 / 30 / 45 MJ/kg (inside the Ames capability 9.3-46.5 MJ/kg, Ellerby 2003
         slide 8, DERIVED conversion); test gas air (ASSUMPTION).
  V8-2  HIT high-frequency plasma wind tunnel (Hu et al. 2013, Materials 6:1730): calorimeter 4.78 / 3.80 MW/m^2,
         27.9 / 20.8 MJ/kg, "18 kPa" (taken as the stagnation pressure), 600 s; flat face of a 30 mm cylinder (1D slab
         30 mm, adiabatic back); ZrB2-30SiC (targets) and ZS10/ZS15 (trend).
Model: 1D (a flat-faced coupon is outside the 2D nose geometry), production settings, measured q_cw on a fully
catalytic cold wall -> Phi = q_cw / (H - h(300 K)), mass-transfer conductance = Phi (D-HEAT-3); fully catalytic and
non-catalytic bounds; no parameter adjusted.
Acceptance bands, set before the run (helper's proposal in V8_TARGET.md, adopted unchanged):
  V8-1 passive peak T 1630 C +/- 80 K; oxide 110 um within x1.5; depletion present, 40 um within x2;
        active: loss of steady behaviour at 4.0 but not at 3.0 MW/m^2 (qualitative).
  V8-2 ZS30 steady T 1960 C (cond. 1) and 1900 C (cond. 2) +/- 80 K; ZS30 reacted layer 65/99/125 um at 201/402/600 s
        within x1.5 each; cond. 2 oxide 90 um within x1.5; trend T(ZS10) > T(ZS15) > T(ZS30).
A case passes if at least one catalysis bound meets every band (the test's catalytic state is unknown).
Output: report.md, v8_results.json
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
from uhtc import coupled as cp, materials as mt  # noqa: E402
from uhtc.oxidation import LayerOxidation, OxParams  # noqa: E402

C2K = 273.15


def run(q_cw, H, p, t_end, f_SiC, depth, cat, rec=1.0):
    ox = LayerOxidation(OxParams(f_SiC=f_SiC))
    cond = cp.PrescribedConditions([(0.0, H, p, q_cw)], R_n=0.01, catalysis=cat)
    r = cp.run_1d(cond, mt.composite(f_SiC), ox, t_end, depth=depth, record_every=rec)
    L = [ox.layers(s) for s in r.states]
    return dict(t=r.t.tolist(), T=r.T_s.tolist(), oxide_um=[(l["ZrO2"] + l["glass"]) * 1e6 for l in L],
                ZrO2_um=[l["ZrO2"] * 1e6 for l in L], dep_um=[l["depleted"] * 1e6 for l in L],
                events={k: float(v) for k, v in r.events.items()}, T_max=float(np.max(r.T_s)), T_end=float(r.T_s[-1]))


def at(res, key, t):
    return float(np.interp(t, res["t"], res[key]))


def main():
    out, lines = {}, ["# V8 — arc-jet / plasma tests with measured heat flux (ZrB2-SiC)", "",
                      "Inputs, targets and acceptance bands: `run_v8.py` header and `literature/V8_TARGET.md`. "
                      "1D model, no adjusted parameter.", ""]
    # ---------------- V8-1
    lines += ["## V8-1 NASA Ames AHF, ZrB2-20SiC flat-face coupon (US Patent 5,750,450)", "",
              "| H0 (MJ/kg) | bound | passive peak T (°C) [1630 ± 80] | oxide 300 s (µm) [110, ×1.5] | depleted (µm) [40, ×2] | "
              "active peak T (°C) [2200] | active: T at 50 / 180 s (°C) | pass (passive) |", "|---|---|---|---|---|---|---|---|"]
    ok1 = False
    for H in (15e6, 30e6, 45e6):
        for cat in ("fc", "nc"):
            a = run(3.0e6, H, 0.042 * 101325.0, 300.0, 0.2, 0.00635, cat)
            b = run(4.0e6, H, 0.058 * 101325.0, 180.0, 0.2, 0.00635, cat)
            out[f"V8_1_H{H / 1e6:.0f}_{cat}"] = dict(passive=a, active=b)
            Tp = a["T_max"] - C2K
            oxv, dep = a["oxide_um"][-1], a["dep_um"][-1]
            pas = abs(Tp - 1630) <= 80 and 110 / 1.5 <= oxv <= 110 * 1.5 and 20 <= dep <= 80
            ok1 |= pas
            lines.append(f"| {H / 1e6:.0f} | {cat} | {Tp:.0f} | {oxv:.0f} | {dep:.0f} | {b['T_max'] - C2K:.0f} | "
                         f"{at(b, 'T', 50) - C2K:.0f} / {b['T_end'] - C2K:.0f} | {'yes' if pas else 'no'} |")
    lines += ["", f"V8-1 (passive, quantitative bands): **{'PASS' if ok1 else 'FAIL'}** (at least one enthalpy and bound "
              "inside all bands). The active case is judged qualitatively in the text.", ""]
    # ---------------- V8-2
    lines += ["## V8-2 HIT plasma wind tunnel, ZrB2-SiC cylinders (Hu et al. 2013)", "",
              "| bound | ZS30 T cond.1 (°C) [1960 ± 80] | ZS30 T cond.2 (°C) [1900 ± 80] | ZS30 oxide at 201/402/600 s (µm) "
              "[65/99/125] | oxide+depleted at 201/402/600 s (µm) | cond.2 oxide 600 s (µm) [90] | T ZS10 / ZS15 / ZS30, cond.1 (°C) | pass |",
              "|---|---|---|---|---|---|---|---|"]
    ok2 = False
    for cat in ("fc", "nc"):
        c1 = run(4.78e6, 27.9e6, 18e3, 600.0, 0.3, 0.03, cat)
        c2 = run(3.80e6, 20.8e6, 18e3, 600.0, 0.3, 0.03, cat)
        z10 = run(4.78e6, 27.9e6, 18e3, 300.0, 0.1, 0.03, cat)
        z15 = run(4.78e6, 27.9e6, 18e3, 600.0, 0.15, 0.03, cat)
        out[f"V8_2_{cat}"] = dict(c1=c1, c2=c2, zs10=z10, zs15=z15)
        T1, T2 = c1["T_end"] - C2K, c2["T_end"] - C2K
        lay = [at(c1, "oxide_um", t) for t in (201, 402, 600)]
        layd = [at(c1, "oxide_um", t) + at(c1, "dep_um", t) for t in (201, 402, 600)]
        meas = [65, 99, 125]
        lay_ok = all(m / 1.5 <= v <= m * 1.5 for v, m in zip(lay, meas)) or all(m / 1.5 <= v <= m * 1.5 for v, m in zip(layd, meas))
        o2 = c2["oxide_um"][-1]
        trend = z10["T_end"] > z15["T_end"] > c1["T_end"]
        pas = abs(T1 - 1960) <= 80 and abs(T2 - 1900) <= 80 and lay_ok and 60 <= o2 <= 135 and trend
        ok2 |= pas
        lines.append(f"| {cat} | {T1:.0f} | {T2:.0f} | {'/'.join(f'{v:.0f}' for v in lay)} | {'/'.join(f'{v:.0f}' for v in layd)} | "
                     f"{o2:.0f} | {z10['T_end'] - C2K:.0f} / {z15['T_end'] - C2K:.0f} / {T1:.0f} | {'yes' if pas else 'no'} |")
    lines += ["", f"V8-2: **{'PASS' if ok2 else 'FAIL'}** (at least one catalysis bound inside all bands).", ""]
    # ---------------- diagnostic (not a prediction): oxidation model driven by the measured surface temperature
    from uhtc.coupled import FlowState
    diag = []
    for name, f_SiC, T_C, t_s, q_cw, H, p, meas in (
            ("V8-1 passive (peak T held)", 0.2, 1630.0, 300.0, 3.0e6, 30e6, 0.042 * 101325.0, "110 µm oxide, 40 µm depleted"),
            ("V8-2 ZS30 cond. 1 (steady T)", 0.3, 1960.0, 600.0, 4.78e6, 27.9e6, 18e3, "65/99/125 µm at 201/402/600 s"),
            ("V8-2 ZS30 cond. 2 (steady T)", 0.3, 1900.0, 600.0, 3.80e6, 20.8e6, 18e3, "≈90 µm at 600 s")):
        ox = LayerOxidation(OxParams(f_SiC=f_SiC))
        T = T_C + C2K
        s = ox.init_state(T)
        Phi = q_cw / H
        fl = FlowState(H_s=H, p_s=p, p_O2=0.2095 * p, h_D=0.0, Phi=Phi, f_nc=1.0, g_mass=Phi)
        edges = np.r_[0.0, np.geomspace(1e-3, t_s, 400)]
        vals = {}
        for t0_, t1_ in zip(edges[:-1], edges[1:]):
            s = ox.advance(s, t1_ - t0_, T, T, fl)
            for tt in (201.0, 402.0, 600.0, 300.0):
                if t0_ < tt <= t1_:
                    L = ox.layers(s)
                    vals[tt] = ((L["ZrO2"] + L["glass"]) * 1e6, L["depleted"] * 1e6)
        L = ox.layers(s)
        txt = ", ".join(f"{int(k)} s: oxide {v[0]:.0f}, depleted {v[1]:.0f}" for k, v in sorted(vals.items()) if k <= t_s)
        diag.append(f"- {name}, {T_C:.0f} °C for {t_s:.0f} s: {txt} µm (measured: {meas}).")
    lines += ["## Diagnostic (not a prediction): oxidation driven by the measured surface temperature", "",
              "Isothermal at the measured (peak or steady) temperature, no heat-up; separates the oxidation model from the "
              "heating error.", ""] + diag + [""]
    (HERE / "v8_results.json").write_text(json.dumps(out), encoding="utf-8")
    (HERE / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
