"""Temperature-dependent SiC-depletion rate (external review, item A6 / run R2).

The baseline k_d = 473 um^2/min is fitted to the three printed depths at 1800 C (Shugart 2014, Fig. 44 caption) and does
not depend on temperature.  This script fits an Arrhenius law k_d(T) = k_ref exp(-E/R (1/T - 1/T_ref)) through all of
Shugart's depleted-layer depths at the four temperatures he reports (1627, 1650, 1700 and 1800 C; Fig. 44 caption and
the Fig. 45 readings in literature/notes/oxidation_furnace_validation_points.csv), in two steps:
  1. at each temperature, the parabolic constant through the origin, k_T = sum(x^2 t) / sum(t^2), the estimator of the
     baseline (uhtc.data.K_DEP), with every point, zero depths included;
  2. an unweighted least-squares line of ln k_T against 1/T over the four temperatures.
The paper's own per-temperature refits at 1650 and 1700 C (rates 2-8 times smaller than the baseline, Section 2.5) are
a cross-check, not fit data.  The 1627-1700 C depths are the independent points of check V6, so a model run with this
law can no longer be checked against them; it is a sensitivity case only.

Output: kd_arrhenius.json (E, k_ref at 2073.15 K, the per-temperature constants and the OxParams settings).
Run: python runs/final_review/kd_arrhenius.py
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DATA = ROOT / "literature" / "notes" / "oxidation_furnace_validation_points.csv"
R = 8.314462618
T_REF = 2073.15  # K, 1800 C, the temperature of the baseline fit


def load():
    rows = [r for r in csv.DictReader(l for l in open(DATA, encoding="utf-8") if not l.startswith("#"))]
    pts = {}
    for r in rows:
        if r["source"] == "Shugart2014" and r["quantity"].lower().startswith("sic-depleted"):
            pts.setdefault(float(r["T_C"]), []).append((float(r["t_min"]), float(r["value_um"])))
    return pts


def main():
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from uhtc import data as dt

    pts = load()
    k_base = dt.K_DEP * 1e12 * 60.0  # um^2/min
    per_T = {}
    for T_C in sorted(pts):
        t = np.array([p[0] for p in pts[T_C]])
        x = np.array([p[1] for p in pts[T_C]])
        per_T[T_C] = dict(n=len(t), n_zero=int(np.sum(x == 0)), k_um2_min=float(np.sum(x ** 2 * t) / np.sum(t * t)))
    TK = np.array([T + 273.15 for T in per_T])
    lk = np.log([v["k_um2_min"] for v in per_T.values()])
    A = np.vstack([np.ones_like(TK), 1.0 / TK]).T
    (c0, c1), *_ = np.linalg.lstsq(A, lk, rcond=None)
    E = -c1 * R
    k_ref = math.exp(c0 + c1 / T_REF)
    fit = lambda T_C: k_ref * math.exp(-E / R * (1.0 / (T_C + 273.15) - 1.0 / T_REF))
    resid = {T: math.log(v["k_um2_min"] / fit(T)) for T, v in per_T.items()}
    out = dict(
        E_kJ_mol=E / 1e3, k_ref_um2_min=k_ref, T_ref_K=T_REF, k_base_um2_min=k_base,
        kd_mult=k_ref / k_base, per_T=per_T, ln_residual=resid,
        fit_at={str(T): fit(T) for T in (1627, 1650, 1700, 1800, 1927, 2000, 2150)},
        ratio_base_to_fit={str(T): k_base / fit(T) for T in (1650, 1700)},
        cross_check="paper Section 2.5: per-temperature refits at 1650 and 1700 C give rates 2-8 times smaller than the "
                    "baseline (incubation fits 60 and 269 um^2/min, literature/notes/HighT_kinetics_alternatives.md row 6)",
        oxparams=dict(kd_E=E, kd_T_ref=T_REF, kd_mult=k_ref / k_base),
    )
    (HERE / "kd_arrhenius.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("E_kJ_mol", "k_ref_um2_min", "k_base_um2_min", "kd_mult", "per_T",
                                          "fit_at", "ratio_base_to_fit")}, indent=1))


if __name__ == "__main__":
    main()
