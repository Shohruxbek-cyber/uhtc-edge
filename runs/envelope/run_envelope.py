"""Phase 4, Q1-Q3: constant-condition exposures over the flight envelope.

Grid: Mach 5-12 (step 0.5) x altitude 20-40 km (step 1 km).  Sets:
  A  R_n = 1, 3, 10 mm; catalysis fully catalytic (fc) and non-catalytic (nc); 20 vol% SiC        (1890 cases)
  B  R_n = 3 mm; fc and nc; 10 and 30 vol% SiC                                                       (1260 cases)
  C  R_n = 1, 3, 10 mm; catalysis switched by the oxidation state (D-HEAT-4); 20 vol% SiC          ( 945 cases)
  D  R_n = 1, 3, 10 mm; fc; 20 vol% SiC; passive-to-active line of Balat-Pichelin et al. 2015 (D-OX-12) ( 945 cases)
Each case: 1D stagnation-line model, start at 300 K, exposure up to the cap of 1800 s (D-FAIL-2); the time to failure
is the first of the three criteria.  Sets A, B, D stop at the first failure (revision: with the permeation criterion a
thin, replenished glass stays protective in hot cases, so "run until all three criteria fired" made most hot runs
last the full 1800 s); set C continues after the first failure until all three times are known (or T_s > 3400 K),
because its purpose is what happens after the glass is lost.
Usage: python runs/envelope/run_envelope.py [A|B|C|D ...]    (default: A B C D).  Restartable.
Output: runs/envelope/set_<X>.jsonl and .csv, meta.json.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common  # noqa: E402

HERE = Path(__file__).resolve().parent
MACH = np.round(np.arange(5.0, 12.0 + 1e-9, 0.5), 2)
ALT = np.arange(20e3, 40e3 + 1, 1e3)


def cases(set_name):
    out = []
    if set_name == "A":
        combos = [(R, cat, 0.2) for R in (0.001, 0.003, 0.010) for cat in ("fc", "nc")]
    elif set_name == "B":
        combos = [(0.003, cat, f) for cat in ("fc", "nc") for f in (0.1, 0.3)]
    elif set_name in ("C", "D"):
        combos = [(R, "switch" if set_name == "C" else "fc", 0.2) for R in (0.001, 0.003, 0.010)]
    else:
        raise ValueError(set_name)
    for R, cat, f in combos:
        for M in MACH:
            for z in ALT:
                out.append(dict(id=f"{set_name}_R{R * 1e3:.0f}_{cat}_SiC{f:.2f}_M{M:.1f}_z{z / 1e3:.0f}",
                                mach=float(M), alt=float(z), R_n=R, catalysis=cat, f_SiC=f))
                if set_name == "D":
                    out[-1]["ox"] = dict(pa_law="balat2015")
                if set_name != "C":
                    out[-1]["stop"] = "first"
    return out


if __name__ == "__main__":
    sets = [a for a in sys.argv[1:] if a in ("A", "B", "C", "D")] or ["A", "B", "C", "D"]
    common.write_meta(HERE / "meta.json", grid_mach=MACH.tolist(), grid_alt_m=ALT.tolist(), sets=sets,
                      model="1D, depth 20 mm adiabatic back, start 300 K")
    for s in sets:
        common.run_all(cases(s), HERE / f"set_{s}.jsonl")
