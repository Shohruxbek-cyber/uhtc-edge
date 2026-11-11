"""Task 2.3: mass-transfer conductance with a Lewis-number correction, as one-at-a-time cases on the 1D base cases.

Eq. 6 uses g = heat-transfer conductance (Lewis number 1), while Eq. 2 uses Le = 1.4.  The correction factors tested are
Le^0.52 = 1.191 (the exponent of the Fay-Riddell bracket, Eq. 2) and Le^(2/3) = 1.251 (Chilton-Colburn).  The baseline is
not changed; these runs only place the correction inside the x0.5-x2 band of the OAT set.
Output: runs/sensitivity/oat_lewis.jsonl/.csv
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common  # noqa: E402

HERE = Path(__file__).resolve().parent
LE = 1.4
FACTORS = {"Le052": LE ** 0.52, "Le23": LE ** (2 / 3)}

if __name__ == "__main__":
    bases = json.loads((HERE / "base_cases.json").read_text(encoding="utf-8"))
    cases = [dict(mach=b["mach"], alt=b["alt"], R_n=b["R_n"], catalysis=b["catalysis"], f_SiC=0.2, stop="first",
                  id=f"{bn}_{tag}", mtf=v, ox={}) for bn, b in bases.items() for tag, v in FACTORS.items()]
    common.run_all(cases, HERE / "oat_lewis.jsonl")
