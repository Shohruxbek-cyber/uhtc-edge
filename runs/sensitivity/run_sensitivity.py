"""Phase 4, Q4: which inputs control the time to failure?

Base cases (chosen from the envelope maps so that the baseline time to failure is finite and well inside the 1800 s
cap; see base_cases.json).  For each base case:
  OAT  one-at-a-time low/high values of every input over its literature range (INPUTS.md); tornado of the change in
       time to failure (and of the failure mode);
  LHS  200-sample Latin hypercube over the continuous inputs (uniform; log-uniform for the multipliers on kinetics and
       mass transfer) with the catalysis bound drawn as a 50/50 binary; Spearman rank correlation and partial rank
       correlation (PRCC) of log10(time to failure) with each input; runs that survive the cap are censored at 1800 s.
Usage: python runs/sensitivity/run_sensitivity.py [oat] [lhs]
Output: runs/sensitivity/oat_<base>.jsonl/.csv, lhs_<base>.jsonl/.csv (+ lhs_design_<base>.csv).
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common  # noqa: E402

HERE = Path(__file__).resolve().parent

# name: (where, key, low, base, high, scale)   where in {case, ox}; scale 'lin' or 'log' (LHS sampling)
PARAMS = {
    "emissivity, glass": ("ox", "eps_glass", 0.70, 0.88, 0.95, "lin"),
    "emissivity, oxide (after glass loss)": ("ox", "eps_oxide", 0.50, 0.75, 0.85, "lin"),
    "substrate conductivity (x)": ("case", "k_mult", 0.70, 1.0, 1.05, "lin"),
    "heat capacity (x)": ("case", "cp_mult", 0.97, 1.0, 1.03, "lin"),
    "heating correlation (x)": ("case", "heat_mult", 0.85, 1.0, 1.15, "lin"),
    "oxide growth rate kp (x)": ("ox", "kp_mult", 1 / 3, 1.0, 3.0, "log"),
    "SiC depletion rate (x)": ("ox", "kd_mult", 0.5, 1.0, 2.0, "log"),
    "depletion onset T_dep (K)": ("ox", "T_dep", 1875.15, 1900.15, 1925.15, "lin"),
    "passive-to-active offset (decades)": ("ox", "pa_offset", 0.3, 0.65, 1.0, "lin"),
    "mass-transfer coefficient (x)": ("case", "mtf", 0.5, 1.0, 2.0, "log"),
    "glass conductivity (W/m/K)": ("ox", "k_glass", 0.7, 1.37, 2.8, "lin"),
    "ZrO2 scale conductivity (W/m/K)": ("ox", "k_ZrO2", 1.0, 1.5, 2.5, "lin"),
    "temperature limit (K)": ("ox", "T_limit", 2289.15, 2422.15, 2543.0, "lin"),  # upper: ZrB2-SiC eutectic (D-FAIL-5)
    "SiC fraction": ("case", "f_SiC", 0.1, 0.2, 0.3, "discrete"),
    "nose radius (mm)": ("case", "R_n", None, None, None, "radius"),
    "O2 pressure exponent n_p": ("ox", "n_p", 0.0, 0.0, 0.5, "lin"),
    "SiO2 from depleted layer lost": ("ox", "f_SiO_lost", 0.0, 0.0, 1.0, "lin"),
    "evaporation coefficient": ("ox", "alpha_evap", 0.1, 1.0, 1.0, "log"),
    "catalysis bound": ("case", "catalysis", "nc", "fc", "fc", "binary"),
}


def base_cases():
    return json.loads((HERE / "base_cases.json").read_text(encoding="utf-8"))


def make_case(base, overrides, cid):
    c = dict(mach=base["mach"], alt=base["alt"], R_n=base["R_n"], catalysis=base["catalysis"], f_SiC=0.2,
             stop="first", id=cid, ox={})
    for name, val in overrides.items():
        where, key = PARAMS[name][0], PARAMS[name][1]
        if where == "ox":
            c["ox"][key] = val
        else:
            c[key] = val
    return c


def oat_cases(bname, base):
    cs = [make_case(base, {}, f"{bname}_base")]
    for name, (where, key, lo, b, hi, sc) in PARAMS.items():
        if sc == "radius":
            lo, hi = base["R_n"] * 2 / 3, base["R_n"] * 1.5
        for tag, v in (("low", lo), ("high", hi)):
            if v is None or (sc != "radius" and v == b):
                continue
            cs.append(make_case(base, {name: v}, f"{bname}_{key}_{tag}"))
    return cs


def lhs_design(n, seed=20261102):
    rng = np.random.default_rng(seed)
    names = [k for k, p in PARAMS.items() if p[5] in ("lin", "log") and p[2] != p[4]]
    d = len(names)
    u = (rng.permuted(np.tile(np.arange(n), (d, 1)), axis=1).T + rng.random((n, d))) / n
    X = {}
    for j, name in enumerate(names):
        _, _, lo, b, hi, sc = PARAMS[name]
        X[name] = lo * (hi / lo) ** u[:, j] if sc == "log" else lo + (hi - lo) * u[:, j]
    X["catalysis bound"] = np.where(rng.random(n) < 0.5, "fc", "nc")
    return X


def lhs_cases(bname, base, n=200):
    X = lhs_design(n)
    with open(HERE / f"lhs_design_{bname}.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["case_id"] + list(X))
        for i in range(n):
            w.writerow([f"{bname}_lhs{i:03d}"] + [X[k][i] for k in X])
    return [make_case(base, {k: (float(X[k][i]) if k != "catalysis bound" else str(X[k][i])) for k in X},
                      f"{bname}_lhs{i:03d}") for i in range(n)]


if __name__ == "__main__":
    what = [a for a in sys.argv[1:] if a in ("oat", "lhs")] or ["oat", "lhs"]
    for bname, base in base_cases().items():
        if "oat" in what:
            common.run_all(oat_cases(bname, base), HERE / f"oat_{bname}.jsonl")
        if "lhs" in what:
            common.run_all(lhs_cases(bname, base), HERE / f"lhs_{bname}.jsonl")
