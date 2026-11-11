"""V6: oxidation sub-model against furnace scale thicknesses (isothermal exposures).

Data: literature/notes/oxidation_furnace_validation_points.csv (compiled from Shugart 2014, Opila & Halbig 2001,
Levine & Opila 2003, Opila et al. 2003 preprint; locations in the CSV).  Rows whose value comes only from an abstract
(value_type 'ABSTRACT TEXT': the two Karlsdottir et al. 2007 points, whose full text could not be obtained) are not used
(author decision, 10 Nov 2026); they stay in the CSV for the record and are listed at the end of report.md.
Model: uhtc.oxidation.LayerOxidation with the baseline parameters, run isothermally (surface = front = furnace T).
  furnace gas: p = 101.325 kPa; p_O2 = 21.2 kPa (air) or 101.325 kPa (flowing O2);
  mass-transfer conductance for evaporation in a furnace: g = 0.01 kg/m^2/s (ASSUMPTION, natural convection,
  h/cp ~ 10 W/m^2/K / 1000 J/kg/K; sensitivity 0.003 and 0.03) - affects only the glass thickness;
  cyclic exposures (10 x 10 min) treated as continuous (ASSUMPTION; each cycle includes heat-up, Opila 2010).
Mapping of measured quantities to model quantities:
  ZrO2 layer -> x_r;  borosilicate / glass -> d_g;  total oxide / oxide scale -> x_r + d_g;
  SiC-depleted layer -> d_dep;  recession -> x_r + d_dep;  scale + depleted -> x_r + d_g + d_dep.
Calibration set (fit check only, not validation): Shugart ZrO2 layers (Table 9 source runs: 1300-1550 C box
furnace, and the 1800 C depletion depths used for k_d).  Independent set: everything else.
Pass criterion (set before the run): every independent point within a factor of 3 of the measurement.  Reason: the
same nominal ZrB2-20SiC oxidised at 1627 C for 100 min gave 150 um (Opila & Halbig 2001) and about 55 um (Opila et al.
2003) of oxide, a factor of about 3 between batches; within-specimen scatter is 30-80 % (Shugart Table 3).
Outputs: v6_results.csv, report.md, fig_v6.png.
"""
import csv
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
from uhtc.coupled import FlowState  # noqa: E402
from uhtc.oxidation import LayerOxidation, OxParams  # noqa: E402

DATA = ROOT / "literature" / "notes" / "oxidation_furnace_validation_points.csv"
G_FURNACE = 0.01


def simulate(f_SiC, T, t_min, pO2, g=G_FURNACE, n_steps=400, incubation=False):
    p = OxParams(f_SiC=f_SiC, incubation=incubation)
    ox = LayerOxidation(p)
    s = ox.init_state(T)
    fl = FlowState(H_s=0.0, p_s=101325.0, p_O2=pO2, h_D=0.0, Phi=0.0, f_nc=1.0, g_mass=g)
    t_end = t_min * 60.0
    # geometric steps (parabolic transients)
    edges = np.r_[0.0, np.geomspace(1e-3, t_end, n_steps)]
    for dt in np.diff(edges):
        s = ox.advance(s, dt, T, T, fl)
    L = ox.layers(s)
    return dict(ZrO2=L["ZrO2"], glass=L["glass"], dep=L["depleted"], xB=L["x_B2O3"])


def model_value(q, m):
    q = q.lower()
    if "depleted" in q and "scale" in q:
        return m["ZrO2"] + m["glass"] + m["dep"]
    if q.startswith("recession"):
        return m["ZrO2"] + m["dep"]
    if "sic-depleted" in q or q.startswith("sic depleted"):
        return m["dep"]
    if q.startswith("zro2"):
        return m["ZrO2"]
    if "borosilicate" in q or q.startswith("glass"):
        return m["glass"]
    if "oxide" in q:
        return m["ZrO2"] + m["glass"]
    raise ValueError(q)


def is_calibration(r):
    if r["source"] != "Shugart2014":
        return False
    q = r["quantity"].lower()
    if q.startswith("zro2") and "box furnace" in r["set_up"] and float(r["t_min"]) == 100:
        return True
    if "depleted" in q and float(r["T_C"]) == 1800:
        return True
    return False


def main():
    allrows = [r for r in csv.DictReader(l for l in open(DATA, encoding="utf-8") if not l.startswith("#"))]
    excluded = [r for r in allrows if r["value_type"].upper().startswith("ABSTRACT")]  # value from an abstract only
    rows = [r for r in allrows if r not in excluded]
    out = []
    for r in rows:
        try:
            val = float(r["value_um"])
        except ValueError:
            continue
        T = float(r["T_C"]) + 273.15
        pO2 = 101325.0 if "O2" in r["atmosphere"] and "air" not in r["atmosphere"] else 0.2095 * 101325.0
        m = simulate(float(r["vol_pct_SiC"]) / 100.0, T, float(r["t_min"]), pO2)
        mv = model_value(r["quantity"], m) * 1e6
        glass_q = "borosilicate" in r["quantity"].lower()
        lo = hi = None
        if not glass_q:
            mlo = model_value(r["quantity"], simulate(float(r["vol_pct_SiC"]) / 100.0, T, float(r["t_min"]), pO2, g=0.03)) * 1e6
            mhi = model_value(r["quantity"], simulate(float(r["vol_pct_SiC"]) / 100.0, T, float(r["t_min"]), pO2, g=0.003)) * 1e6
            lo, hi = min(mlo, mhi), max(mlo, mhi)
        mi = model_value(r["quantity"], simulate(float(r["vol_pct_SiC"]) / 100.0, T, float(r["t_min"]), pO2,
                                                 incubation=True)) * 1e6  # option D-OX-10
        out.append(dict(source=r["source"], SiC=r["vol_pct_SiC"], T_C=r["T_C"], t_min=r["t_min"],
                        atmosphere=r["atmosphere"][:40], quantity=r["quantity"][:45], measured_um=val,
                        sd_um=r["sd_or_uncert_um"], model_um=round(mv, 2), ratio=round(mv / val, 3) if val > 0 else np.nan,
                        model_inc_um=round(mi, 2), ratio_inc=round(mi / val, 3) if val > 0 else np.nan,
                        set="calibration" if is_calibration(r) else ("glass (depends on furnace g)" if glass_q else "independent"),
                        location=r["location"]))
    with open(HERE / "v6_results.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    ind = [o for o in out if o["set"] == "independent" and o["measured_um"] > 0]
    cal = [o for o in out if o["set"] == "calibration"]
    gl = [o for o in out if o["set"].startswith("glass")]
    within = [o for o in ind if 1 / 3 <= o["ratio"] <= 3]
    ok = len(within) == len(ind)
    gm = lambda L: math.exp(np.mean([math.log(o["ratio"]) for o in L])) if L else float("nan")
    dep = [o for o in ind if o["quantity"].lower().startswith("sic-depleted")]
    oxq = [o for o in ind if o not in dep]
    n_in = lambda L: sum(1 for o in L if 1 / 3 <= o["ratio"] <= 3)
    n_in_i = lambda L: sum(1 for o in L if 1 / 3 <= o["ratio_inc"] <= 3)
    gm_i = lambda L: (math.exp(np.mean([math.log(o["ratio_inc"]) for o in L])) if L and all(o["ratio_inc"] > 0 for o in L)
                      else float("nan"))
    lines = ["# V6 — Oxidation sub-model vs furnace scale thicknesses", "",
             f"Result: **{'PASS' if ok else 'FAIL'}** — {len(within)} of {len(ind)} independent points within a factor of 3 "
             f"(criterion, set before the run: all). Geometric-mean model/measured ratio: independent {gm(ind):.2f}, "
             f"calibration set {gm(cal):.2f}, glass thickness {gm(gl):.2f}.", "",
             f"- Oxide scale, ZrO2 layer, recession and scale + depleted zone (independent): {n_in(oxq)} of {len(oxq)} within "
             f"a factor of 3, geometric-mean ratio {gm(oxq):.2f}. The outliers are batch differences in opposite "
             "directions: Opila & Halbig's commercial ZrB2-20SiC at 1327 °C (10 and 100 min) recedes 4-8 times faster than the model "
             "(calibrated on Shugart's material), and the NASA-made batch of Opila et al. 2003 grew about 5 times less "
             "oxide at 1627 °C than the model (and about 3 times less than the commercial batch at the same condition).",
             f"- SiC-depleted layer at 1627-1700 °C (Shugart, digitized; independent of the 1800 °C calibration): "
             f"{n_in(dep)} of {len(dep)} within a factor of 3 (zero values excluded from the ratio). The model has no "
             "incubation period, while Shugart reports that depletion needs 10-20 min to start at 1650 °C; it therefore "
             "over-predicts the depleted depth at short times. Consequence for the flight results: the SiO2 that the model "
             "adds to the glass from the depleted layer arrives too early, which makes glass survival somewhat optimistic; "
             "the sensitivity case f_SiO_lost = 1 (no SiO2 from the depleted layer) bounds this (Q4).",
             "- Glass thickness depends on the assumed furnace mass-transfer conductance and is reported, not judged.", "",
             "**With the incubation option** (D-OX-10: depletion starts after an incubation time taken from the onset "
             "brackets of Shugart's Fig. 40, bracket midpoints; not fitted to the thicknesses below, but from the same "
             f"experiments as the depleted-layer points, so not independent of them): {n_in_i(ind)} of {len(ind)} "
             f"independent points within a factor of 3 (oxide-type {n_in_i(oxq)} of {len(oxq)}, depleted layer "
             f"{n_in_i(dep)} of {len(dep)}); geometric-mean ratio, depleted layer: {gm_i(dep):.2f} with the option, "
             f"{gm(dep):.2f} without. Against the criterion set before the run (all within a factor of 3) V6 "
             f"{'passes' if n_in_i(ind) == len(ind) else 'still fails'} with the option.", "",
             "Method, data mapping and the reason for the tolerance: `run_v6.py` header. Furnace mass-transfer conductance "
             f"g = {G_FURNACE} kg/m²/s (ASSUMPTION); it changes only glass-related quantities.", "",
             "| Set | Source | SiC vol% | T (°C) | t (min) | Quantity | Measured (µm) | ± | Model (µm) | Model/measured |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for o in sorted(out, key=lambda o: (o["set"], o["source"], float(o["T_C"]), float(o["t_min"]))):
        flag = "" if o["set"] != "independent" or (1 / 3 <= o["ratio"] <= 3) else " **out**"
        lines.append(f"| {o['set']} | {o['source']} | {o['SiC']} | {o['T_C']} | {o['t_min']} | {o['quantity']} | "
                     f"{o['measured_um']:g} | {o['sd_um']} | {o['model_um']:.1f} | {o['ratio']:.2f}{flag} |")
    if excluded:
        lines += ["", "Not used (value known only from an abstract; author decision, 10 Nov 2026): " + "; ".join(
            f"{r['source']} {r['quantity']} {r['value_um']} µm at {r['T_C']} °C, {r['t_min']} min ({r['location']})" for r in excluded) + "."]
    (HERE / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines[:3]))
    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    for setname, mk, col in (("calibration", "s", "0.55"), ("independent", "o", "C0"), ("glass (depends on furnace g)", "^", "C1")):
        L = [o for o in out if o["set"] == setname and o["measured_um"] > 0]
        ax.loglog([o["measured_um"] for o in L], [o["model_um"] for o in L], mk, color=col, mfc="none" if setname != "independent" else col,
                  label=setname.split(" (")[0], ms=5)
    x = np.array([1, 3000])
    ax.plot(x, x, "k-", lw=0.8)
    ax.plot(x, 3 * x, "k:", lw=0.8)
    ax.plot(x, x / 3, "k:", lw=0.8)
    ax.set_xlim(2, 3000)
    ax.set_ylim(2, 3000)
    ax.set_xlabel("Measured thickness (µm)")
    ax.set_ylabel("Model thickness (µm)")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(HERE / "fig_v6.png", dpi=200)
    fig.savefig(HERE / "fig_v6.pdf", metadata={"CreationDate": None})


if __name__ == "__main__":
    main()
