"""Revision Part 1.2-1.3: 2D production maps of the nose tip and the uncertainty band of the 30-minute boundary.

Geometry and model: axisymmetric sphere-cone nose tip, half-angle 5 deg, 30 mm flank, adiabatic back plane
(D-2D-2/4), baseline heating distribution BASE_DIST (D-2D-7), production mesh and 2 s step cap (D-NUM-3), 20 vol% SiC,
1800 s cap, stopped at the first failure.

Stages (python run_envelope_2d.py <stage>; every stage is restartable, results in runs/envelope_2d/<stage>.jsonl/.csv):
  grid     Mach 5-12 (step 0.5) x 20-40 km (step 1) x R_n 1/3/10 mm x fully/non-catalytic bound (1890 runs)
  refine   for every (R_n, bound, altitude) column: the midpoint (+0.25) between the highest surviving and the lowest
           failing grid Mach number
  band     fully catalytic bound, altitudes 20/25/30/35/40 km, three radii: the 30-minute boundary re-located for the
           low/high heating distribution, low/high stagnation heating (x0.85/1.15), low/high emissivities, and the
           two combined corners "worst" and "best" (adaptive search in steps of 0.5, then +0.25; D-2D-10)
  glass    as band, for the glass-loss criterion variants (fixed thresholds 0.1/0.5/1 um; D-FAIL-4)
  active   as band, with the passive-to-active line of Balat-Pichelin et al. 2015 for criterion (b) (D-OX-12)
  kinetics as band, for the growth-law extrapolation band, the incubation option and the permeability band
           (D-OX-10, D-OX-11, D-FAIL-4)
  flank    fully catalytic bound, five altitudes, R_n 1 and 3 mm: boundary for flank heating x0.5 and x2 and for the
           sharp-cone running-length correlation (D-2D-12)
  oxidant  as band, sensitivity case with an effective oxidant pressure that counts O atoms as O2 (D-OX-14)
  transport  as band, final review 8 Nov 2026 (R1-05, R1-06, R3-03, R4-02): the inputs of the glass balance and the scale
           resistance that the band left out: mass-transfer conductance x1.25 (Lewis-number factor Le^(2/3)) and x2,
           evaporation coefficient 0.1 and 0.005, and the ZrO2-scale conductivity at the ends of its range (1.0, 2.5)
  active also runs "pa_fe": Eq. 8 refitted with study fixed effects (common slope, R1-01)
  external as band, external review 9 Nov 2026: temperature-dependent depletion rate (A6, run R2; Arrhenius law of
           runs/final_review/kd_arrhenius.py, extrapolated above 2073 K, and held at its 2073 K value above it) and
           the converted-depth gate of criterion (a) at 0.5 and 2 um instead of 1 um (A9, run R5)
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common  # noqa: E402

HERE = Path(__file__).resolve().parent
BASE_DIST = "lees"  # sourced baseline (D-2D-7); "tauber" = draft distribution, kept as a comparison
MACHS = [5.0 + 0.5 * i for i in range(15)]
ALTS = [20e3 + 1e3 * i for i in range(21)]
RADII = [0.001, 0.003, 0.010]
BAND_ALTS = [20e3, 25e3, 30e3, 35e3, 40e3]
BAND_VARIANTS = {  # name: overrides of the baseline case (D-2D-10)
    "heat_lo": dict(heat_dist="lees_lo"), "heat_hi": dict(heat_dist="lees_hi"),
    "q0_lo": dict(heat_mult=0.85), "q0_hi": dict(heat_mult=1.15),
    "eps_lo": dict(ox=dict(eps_glass=0.70, eps_oxide=0.50)), "eps_hi": dict(ox=dict(eps_glass=0.95, eps_oxide=0.85)),
    "worst": dict(heat_dist="lees_hi", heat_mult=1.15, ox=dict(eps_glass=0.70, eps_oxide=0.50)),
    "best": dict(heat_dist="lees_lo", heat_mult=0.85, ox=dict(eps_glass=0.95, eps_oxide=0.85)),
}
ACTIVE_VARIANTS = {"pa_balat": dict(ox=dict(pa_law="balat2015")),  # D-OX-12
                   "pa_fe": dict(ox=dict(pa_law="jacobson_fe")),  # final review R1-01: study fixed-effects slope
                   # third pass (R4p3-03): the refitted line with the low end of the offset range of Table 1
                   "pa_fe_off0.3": dict(ox=dict(pa_law="jacobson_fe", pa_offset=0.3))}
TRANSPORT_VARIANTS = {  # final review: inputs of the glass balance and the scale resistance not in the band
    "mtf_x1.25": dict(mtf=1.25), "mtf_x2": dict(mtf=2.0),
    "alpha_0.1": dict(ox=dict(alpha_evap=0.1)), "alpha_0.005": dict(ox=dict(alpha_evap=0.005)),
    "kZrO2_1.0": dict(ox=dict(k_ZrO2=1.0)), "kZrO2_2.5": dict(ox=dict(k_ZrO2=2.5)),
}
KINETICS_VARIANTS = {  # D-OX-11 band of the extrapolated growth law and the incubation option D-OX-10
    "kp_E250": dict(ox=dict(kp_E_hi=250e3)), "kp_E354": dict(ox=dict(kp_E_hi=354e3)),
    "kp_x0.5": dict(ox=dict(kp_mult_hi=0.5)), "kp_x3": dict(ox=dict(kp_mult_hi=3.0)),
    "incubation": dict(ox=dict(incubation=True)), "perm_low": dict(ox=dict(perm_case="low")),
    "perm_high": dict(ox=dict(perm_case="high")),
    # final review, second pass (R4p2-06): the whole growth law, inside its fit range too, x1/3 and x3 (Table 1 range)
    "kp_all_x0.33": dict(ox=dict(kp_mult=1.0 / 3.0)), "kp_all_x3": dict(ox=dict(kp_mult=3.0)),
}
FLOW_VARIANTS = {  # version 3, Task 1: glass flow on (uhtc.glassflow, D-FLOW-3 to D-FLOW-7), baseline and sensitivity
    "flow": dict(glass_flow="on"),
    "flow_alpha_0.005": dict(glass_flow="on", ox=dict(alpha_evap=0.005)),  # compare transport_alpha_0.005 (flow off)
    "flow_F320": dict(glass_flow="on", flow_visc="F320"),
    "flow_tau_x0.5": dict(glass_flow="on", flow_tau_mult=0.5), "flow_tau_x2": dict(glass_flow="on", flow_tau_mult=2.0),
    "flow_no_pressure": dict(glass_flow="on", flow_pressure=False),
}
FLANKINT_VARIANTS = {"flank_integral": dict(heat_dist="integral")}  # version 3, Task 3 (D-HEAT-7; V12 item 3), all 15 columns
GEOMETRY_VARIANTS = {  # version 3, Task 4.1: boundary re-search for the geometry cases, 1 mm tip only, five band altitudes
    "theta_10": dict(theta_deg=10.0), "flank_15mm": dict(L_flank=0.015), "flank_60mm": dict(L_flank=0.06),
}
FLANK_VARIANTS = {"flank_x0.5": dict(flank_mult=0.5), "flank_x2": dict(flank_mult=2.0),  # D-2D-12
                  "sharp_cone": dict(heat_dist="tauber")}
OXIDANT_VARIANTS = {"o_as_o2": dict(oxidant="O_as_O2")}  # D-OX-14, sensitivity case only
GLASS_VARIANTS = {"dfail_0.1": dict(ox=dict(glass_criterion="fixed", d_fail=0.1e-6)),
                  "dfail_0.5": dict(ox=dict(glass_criterion="fixed", d_fail=0.5e-6)),
                  "dfail_1.0": dict(ox=dict(glass_criterion="fixed", d_fail=1.0e-6))}
_KD_FILE = Path(__file__).resolve().parents[1] / "final_review" / "kd_arrhenius.json"
_KD = json.loads(_KD_FILE.read_text(encoding="utf-8"))["oxparams"] if _KD_FILE.exists() else None
EXTERNAL_VARIANTS = {  # external review 9 Nov 2026 (runs R2 and R5)
    "gate_0.5": dict(ox=dict(x_ox_min=0.5e-6)), "gate_2": dict(ox=dict(x_ox_min=2e-6)),
}
if _KD is not None:
    EXTERNAL_VARIANTS = {"kd_arrh": dict(ox=dict(kd_E=_KD["kd_E"], kd_mult=_KD["kd_mult"])),
                         "kd_arrh_cap": dict(ox=dict(kd_E=_KD["kd_E"], kd_mult=_KD["kd_mult"], kd_T_max=_KD["kd_T_ref"])),
                         **EXTERNAL_VARIANTS}


def case(M, z, R, cat, tag="base", **over):
    c = dict(mach=M, alt=z, R_n=R, catalysis=cat, f_SiC=0.2, geometry="nose", heat_dist=BASE_DIST, stop="first",
             id=f"{tag}_R{R * 1e3:g}_{cat}_M{M:g}_z{z / 1e3:g}")
    ox = dict(over.pop("ox", {}))
    c.update(over)
    if ox:
        c["ox"] = ox
    return c


def survives(r):
    return r["mode"] == "none" and not r.get("error")


def columns(rows):
    """{(R, cat, z): (highest surviving Mach, lowest failing Mach)} from result rows."""
    out = {}
    for r in rows:
        k = (r["R_n"], r["catalysis"], r["alt"])
        lo, hi = out.get(k, (None, None))
        if survives(r):
            lo = r["mach"] if lo is None else max(lo, r["mach"])
        else:
            hi = r["mach"] if hi is None else min(hi, r["mach"])
        out[k] = (lo, hi)
    return out


def boundary_search(tag, over, start, out):
    """Adaptive 30-minute boundary search for a variant: starting from the baseline boundary in each column, step by
    0.5 Mach until survival and failure are bracketed, then one midpoint (resolution 0.25)."""
    start = {k: v for k, v in start.items()}
    for _ in range(12):  # rounds; each round runs the next Mach number of every unbracketed column in parallel
        done = common.load(out) if out.exists() else []
        cols = columns(done)
        todo = []
        for (R, cat, z), (lo0, hi0) in start.items():
            lo, hi = cols.get((R, cat, z), (None, None))
            if lo is not None and hi is not None and hi - lo <= 0.5:
                if hi - lo == 0.5:
                    todo.append(case(lo + 0.25, z, R, cat, tag, **json.loads(json.dumps(over))))
                continue
            if lo is None and hi is None:
                M = lo0 if lo0 is not None else 5.0
            elif lo is None:
                M = hi - 0.5
            else:
                M = lo + 0.5
            if 5.0 <= M <= 12.0:
                todo.append(case(M, z, R, cat, tag, **json.loads(json.dumps(over))))
        todo = [c for c in todo if c["id"] not in {r["case_id"] for r in done}]
        if not todo:
            return
        common.run_all(todo, out, runner=common.run_case_2d)


if __name__ == "__main__":
    stage = sys.argv[1]
    if stage == "grid":
        cases = [case(M, z, R, cat) for R in RADII for cat in ("fc", "nc") for z in ALTS for M in MACHS]
        common.run_all(cases, HERE / "grid.jsonl", runner=common.run_case_2d)
    elif stage == "refine":
        cols = columns(common.load(HERE / "grid.jsonl"))
        cases = [case(lo + 0.25, z, R, cat, "refine") for (R, cat, z), (lo, hi) in sorted(cols.items())
                 if lo is not None and hi is not None and hi - lo == 0.5]
        common.run_all(cases, HERE / "refine.jsonl", runner=common.run_case_2d)
    elif stage == "flowgrid":  # version 3, Task 1: fully catalytic maps of the three radii with glass flow on, per radius
        radii = [float(a) * 1e-3 for a in sys.argv[2:]] or RADII
        for R in radii:
            cases = [case(M, z, R, "fc", "flowgrid", glass_flow="on") for z in ALTS for M in MACHS]
            common.run_all(cases, HERE / f"flowgrid_R{R * 1e3:g}.jsonl", runner=common.run_case_2d)
            cols = columns(common.load(HERE / f"flowgrid_R{R * 1e3:g}.jsonl"))
            cases = [case(lo + 0.25, z, R, "fc", "flowrefine", glass_flow="on") for (R_, cat, z), (lo, hi) in sorted(cols.items())
                     if lo is not None and hi is not None and hi - lo == 0.5]
            common.run_all(cases, HERE / f"flowrefine_R{R * 1e3:g}.jsonl", runner=common.run_case_2d)
    elif stage == "lpgrid":  # journal revision, ASK 5: the full 2D grid (both bounds) with the local wall pressure at every
        # surface node (uhtc.coupled.run_2d, local_pressure=True), failure criteria at the stagnation node as in the baseline;
        # a check of all 126 columns, not a new baseline. Per radius, each followed by its midpoint refinement.
        radii = [float(a) * 1e-3 for a in sys.argv[2:]] or RADII
        for R in radii:
            cases = [case(M, z, R, cat, "lpgrid", local_pressure=True) for cat in ("fc", "nc") for z in ALTS for M in MACHS]
            common.run_all(cases, HERE / f"lpgrid_R{R * 1e3:g}.jsonl", runner=common.run_case_2d)
            cols = columns(common.load(HERE / f"lpgrid_R{R * 1e3:g}.jsonl"))
            cases = [case(lo + 0.25, z, R_, cat, "lprefine", local_pressure=True) for (R_, cat, z), (lo, hi) in sorted(cols.items())
                     if lo is not None and hi is not None and hi - lo == 0.5]
            common.run_all(cases, HERE / f"lprefine_R{R * 1e3:g}.jsonl", runner=common.run_case_2d)
    elif stage in ("band", "glass", "active", "kinetics", "flank", "oxidant", "transport", "external", "flow", "flankint", "geometry"):
        cols = columns(common.load(HERE / "grid.jsonl"))
        variants = dict(band=BAND_VARIANTS, glass=GLASS_VARIANTS, active=ACTIVE_VARIANTS, kinetics=KINETICS_VARIANTS,
                        flank=FLANK_VARIANTS, oxidant=OXIDANT_VARIANTS, transport=TRANSPORT_VARIANTS,
                        external=EXTERNAL_VARIANTS, flow=FLOW_VARIANTS, flankint=FLANKINT_VARIANTS,
                        geometry=GEOMETRY_VARIANTS)[stage]
        radii = RADII[:2] if stage == "flank" else RADII[:1] if stage == "geometry" else RADII
        for name, over in variants.items():
            start = {(R, "fc", z): cols[(R, "fc", z)] for R in radii for z in BAND_ALTS}
            boundary_search(name, over, start, HERE / f"{stage}_{name}.jsonl")
    else:
        raise SystemExit(__doc__)
