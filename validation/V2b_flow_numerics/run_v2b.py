"""V2b: mesh and step convergence of the 2D model with glass flow on (plan: validation/V2b_plan.md, committed before this run).

Usage: python validation/V2b_flow_numerics/run_v2b.py   -> v2b.jsonl/.csv, report.md
"""
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "runs"))
import common  # noqa: E402

FINE = dict(n_arc=60, n_flank=80, n_eta=45, eta1=0.0015)  # as runs/envelope_2d/check_mesh.py
MACHS = (8.5, 9.0, 9.5, 10.0)
SETTINGS = (("prod", common.MESH_2D, 2.0), ("fine", FINE, 0.5))


def main():
    cases = [dict(mach=M, alt=30e3, R_n=1e-3, catalysis="fc", f_SiC=0.2, geometry="nose", heat_dist="lees", stop="first",
                  t_cap=1800.0, glass_flow="on", mesh=mesh, dt_max=dtm, id=f"M{M:g}_{tag}") for M in MACHS for tag, mesh, dtm in SETTINGS]
    out = HERE / "v2b.jsonl"
    common.run_all(cases, out, runner=common.run_case_2d)
    rows = {r["case_id"]: r for r in common.load(out)}
    res, ok_all = [], True
    for M in MACHS:
        a, b = rows[f"M{M:g}_prod"], rows[f"M{M:g}_fine"]
        dT = abs(a["T_s_max"] / b["T_s_max"] - 1.0)
        same_mode = a["mode"] == b["mode"] and not a.get("error") and not b.get("error")
        fails = a["t_fail"] == a["t_fail"] and b["t_fail"] == b["t_fail"]
        dt_rel = abs(a["t_fail"] / b["t_fail"] - 1.0) if fails else math.nan
        ok = same_mode and dT <= 0.005 and (not fails or dt_rel <= 0.02)
        ok_all &= ok
        res.append(dict(mach=M, mode_prod=a["mode"], mode_fine=b["mode"], T_prod=a["T_s_max"], T_fine=b["T_s_max"], dT_rel=dT,
                        t_prod=a["t_fail"], t_fine=b["t_fail"], dt_rel=dt_rel, ok=ok,
                        audit=max(a.get("flow_audit_max", 0.0), b.get("flow_audit_max", 0.0))))
    (HERE / "v2b_results.json").write_text(json.dumps(dict(passed=ok_all, results=res), indent=1, default=float), encoding="utf-8")
    L = ["# V2b: 2D convergence with glass flow on", "", "Plan: `validation/V2b_plan.md` (committed before this run). R_n = 1 mm, 30 km, "
         "fully catalytic, Lees distribution, flow on (LA, shear x1, pressure on). Production mesh and 2 s cap against the refined "
         "mesh and 0.5 s cap. Criterion: peak temperature within 0.5 %, same first failure, failure time within 2 %.", "",
         "| Mach | Mode (prod / fine) | Peak T prod / fine (K) | Difference | Failure time prod / fine (s) | Difference | Film audit | Meets |",
         "|---|---|---|---|---|---|---|---|"]
    for r in res:
        tt = (f"{r['t_prod']:.1f} / {r['t_fine']:.1f}", f"{100 * r['dt_rel']:.2f} %") if r["dt_rel"] == r["dt_rel"] else ("none", "")
        L.append(f"| {r['mach']:g} | {r['mode_prod']} / {r['mode_fine']} | {r['T_prod']:.1f} / {r['T_fine']:.1f} | {100 * r['dT_rel']:.3f} % | "
                 f"{tt[0]} | {tt[1]} | {r['audit']:.1e} | {'yes' if r['ok'] else 'NO'} |")
    L += ["", f"**V2b {'PASSES' if ok_all else 'FAILS'}.**"]
    (HERE / "report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
