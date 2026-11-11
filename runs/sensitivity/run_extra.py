"""Phase 5 / revision additions to Q4: one-at-a-time cases that run_sensitivity.py does not include.
  k of the SiC-depleted layer 10 and 60 W/m/K (baseline 30);
  atomic-oxygen case of the passive-to-active line (offset 0.65 + 1.65 = 2.3 decades, INPUTS section 8);
  revision (Parts 2.3, 3.1, 3.3, 4):
  glass loss at a fixed threshold of 0.1, 0.5 and 1 um instead of criterion (a) (D-FAIL-4; the first version of this
  script set only d_fail, which criterion (a) ignores - those cases changed nothing and were replaced);
  O2 permeability of the glass low / high (D-FAIL-4); growth law above 2073 K: E = 250 / 354 kJ/mol and k_p x0.5 / x3
  (D-OX-11); incubation of SiC depletion (D-OX-10); passive-to-active line of Balat-Pichelin et al. 2015 (D-OX-12).
The Lewis number (1.4 -> 1.0) enters only the non-catalytic bound and moves it towards the catalytic one, so its effect
is bounded by the catalysis-bound row of the OAT set and is not run separately.
Case ids are <base>_<key>_<tag> with tag low/high (figure convention: low = blue, high = orange).
Base cases and run settings as in run_sensitivity.py (base_cases.json).  Output: runs/sensitivity/oat_extra.jsonl/.csv
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common  # noqa: E402

HERE = Path(__file__).resolve().parent
EXTRA = {"k_dep_low": {"k_dep": 10.0}, "k_dep_high": {"k_dep": 60.0}, "pa_offset_atomicO": {"pa_offset": 2.3},
         "dfix0.1_high": {"glass_criterion": "fixed", "d_fail": 0.1e-6},
         "dfix0.5_high": {"glass_criterion": "fixed", "d_fail": 0.5e-6},
         "dfix1.0_high": {"glass_criterion": "fixed", "d_fail": 1.0e-6},
         "perm_low": {"perm_case": "low"}, "perm_high": {"perm_case": "high"},
         "kpE_low": {"kp_E_hi": 250e3}, "kpE_high": {"kp_E_hi": 354e3},
         "kphi_low": {"kp_mult_hi": 0.5}, "kphi_high": {"kp_mult_hi": 3.0},
         "incub_high": {"incubation": True}, "paBalat_high": {"pa_law": "balat2015"}}

if __name__ == "__main__":
    bases = json.loads((HERE / "base_cases.json").read_text(encoding="utf-8"))
    cases = []
    for bn, b in bases.items():
        for tag, ox in EXTRA.items():
            cases.append(dict(mach=b["mach"], alt=b["alt"], R_n=b["R_n"], catalysis=b["catalysis"], f_SiC=0.2,
                              stop="first", id=f"{bn}_{tag}", ox=ox))
    common.run_all(cases, HERE / "oat_extra.jsonl")
