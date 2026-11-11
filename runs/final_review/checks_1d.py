"""Final review (8 Nov 2026): short 1D checks on the base cases B1 and B2 asked for by the reviewers.

Base cases as in runs/sensitivity (base_cases.json; f_SiC 0.2, stop at the first failure, 1800 s cap).
  hot_slope      B1 with E = 138 kJ/mol above 2073 K (the local slope of Shugart's constants between 1700 and 1800 C),
                 anchored at the fit value at 2073 K (R1-08)
  pore_silica    B1 with the pores of the ZrO2 skeleton filled with silica only (OxParams.pore_fill, R1-11)
  slab_50mm      B1 and B2 with a 50 mm instead of a 20 mm slab (R3-06)
  numerics       B1 at the V2 'finest' settings (dx0 0.5 um, growth 1.04, dT_max 2.5 K, oxide step x0.25) (R3-10)
  o2_budget      B1: O2 consumed by the advancing ZrO2 front against the O2 the outer glass passes by Henry's-law
                 permeation, Pi p_O2,w / delta_g, while the glass exists (R1-03)
Writes runs/final_review/checks_1d.json and checks_1d.md.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "runs"))
from uhtc import coupled as cp, materials as mt  # noqa: E402
from uhtc.oxidation import LayerOxidation, OxParams  # noqa: E402
from uhtc.thermo import MOLAR_MASS  # noqa: E402
import common  # noqa: E402

OUT = Path(__file__).resolve().parent
BASE = json.loads((ROOT / "runs" / "sensitivity" / "base_cases.json").read_text(encoding="utf-8"))


def case(b, **kw):
    c = dict(mach=BASE[b]["mach"], alt=BASE[b]["alt"], R_n=BASE[b]["R_n"], catalysis=BASE[b]["catalysis"], f_SiC=0.2,
             stop="first", id=b, ox={})
    c.update(kw)
    return c


def brief(r):
    return {k: r[k] for k in ("t_fail", "mode", "t_glass", "t_active", "t_limit", "T_s_max", "error") if k in r}


def numerics_b1():
    """B1 at production and at the V2 'finest' settings (same refinement as validation/V2_numerics/run_v2.py)."""
    out = {}
    b = BASE["B1"]
    for name, st in (("baseline", dict(dx0=2e-6, growth=1.08, dT_max=10.0, ox=0.05)),
                     ("finest", dict(dx0=0.5e-6, growth=1.04, dT_max=2.5, ox=0.0125))):
        ox = LayerOxidation(OxParams(f_SiC=0.2))
        orig = ox.rate_limit
        ox.rate_limit = lambda s, T, f, _o=orig, _k=st["ox"] / 0.05: _o(s, T, f) * _k
        state = common.stag_state(b["mach"], b["alt"])
        cond = cp.ConstantConditions(b["mach"], b["alt"], b["R_n"], state=state, catalysis=b["catalysis"])
        r = cp.run_1d(cond, mt.composite(0.2), ox, common.T_CAP, dx0=st["dx0"], growth=st["growth"],
                      dT_max=st["dT_max"], depth=0.02, stop_on=common.stop_first)
        tf, mode = cp.first_failure(r.events)
        out[name] = dict(t_fail=tf, mode=mode, T_s_max=float(np.max(r.T_s)), steps=len(r.t))
    out["rel_diff_t_fail"] = out["baseline"]["t_fail"] / out["finest"]["t_fail"] - 1.0
    return out


def o2_budget_b1():
    """Ratio of the O2 demand of the front to the permeation capacity of the outer glass along B1."""
    b = BASE["B1"]
    ox = LayerOxidation(OxParams(f_SiC=0.2))
    state = common.stag_state(b["mach"], b["alt"])
    cond = cp.ConstantConditions(b["mach"], b["alt"], b["R_n"], state=state, catalysis=b["catalysis"])
    r = cp.run_1d(cond, mt.composite(0.2), ox, common.T_CAP, depth=0.02, stop_on=common.stop_first, record_every=1.0)
    nu = (1 - 0.2) * ox.cZ * 2.5 + 0.2 * ox.cS * 1.5  # mol O2 per m^3 converted: ZrB2 + 5/2 O2, SiC + 3/2 O2
    rows = []
    for i in range(1, len(r.t)):
        s0, s1 = r.states[i - 1], r.states[i]
        dg = ox.glass_thickness(s1)
        if dg <= 0.0 or r.t[i] < 10.0:
            continue
        v = (s1["x_r"] - s0["x_r"]) / (r.t[i] - r.t[i - 1])
        flow = cond.at(r.t[i])
        perm = ox.o2_permeability(r.T_s[i], ox.x_B(s1)) * flow.p_O2 / dg
        j_gas = flow.g_mass * math.log1p(0.2315) / MOLAR_MASS["O2"]
        rows.append((float(r.t[i]), v * nu, perm, j_gas))
    a = np.array(rows)
    ratio = a[:, 1] / a[:, 2]
    tg = r.events["t_glass_loss"]
    return dict(t_glass=tg, t_range=[float(a[0, 0]), float(a[-1, 0])],
                demand_mol_m2s=[float(a[:, 1].min()), float(a[:, 1].max())],
                permeation_mol_m2s=[float(a[:, 2].min()), float(a[:, 2].max())],
                demand_over_permeation=[float(ratio.min()), float(np.median(ratio)), float(ratio.max())],
                gas_supply_over_demand=[float((a[:, 3] / a[:, 1]).min()), float((a[:, 3] / a[:, 1]).max())])


def main():
    res = {}
    res["B1_base"] = brief(common.run_case(case("B1")))
    res["B2_base"] = brief(common.run_case(case("B2")))
    res["hot_slope"] = brief(common.run_case(case("B1", ox=dict(kp_E_hi=138e3))))
    res["pore_silica"] = brief(common.run_case(case("B1", ox=dict(pore_fill="silica"))))
    res["slab_50mm_B1"] = brief(common.run_case(case("B1", depth=0.05)))
    res["slab_50mm_B2"] = brief(common.run_case(case("B2", depth=0.05)))
    res["numerics_B1"] = numerics_b1()
    res["o2_budget_B1"] = o2_budget_b1()
    (OUT / "checks_1d.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    lines = ["# Final review: 1D checks on B1 and B2 (generated by runs/final_review/checks_1d.py)", ""]
    for k, v in res.items():
        lines.append(f"- **{k}**: {json.dumps(v)}")
    (OUT / "checks_1d.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
