"""Convergence of the 2D production settings (revision, before the 2D maps): mesh (n_arc, n_flank, n_eta, first cell)
and time-step cap, at three conditions spanning the radii, run to 600 s with the production model.
Criterion (set before the run): peak stagnation temperature within 0.5 % and failure time within 2 % of the refined
settings.  Output: runs/envelope_2d/check_mesh.jsonl/.csv and a printed table.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common  # noqa: E402

HERE = Path(__file__).resolve().parent
FINE = dict(n_arc=60, n_flank=80, n_eta=45, eta1=0.0015)
CELLS = [(9.0, 30e3, 0.001), (9.5, 35e3, 0.010), (10.0, 40e3, 0.003)]
# failing cells (coarse 2D map, Phase 4): glass loss after 690 s and 522 s, temperature limit after 48 s
FAIL_CELLS = [(9.5, 30e3, 0.001), (11.0, 40e3, 0.003), (10.0, 30e3, 0.001)]

if __name__ == "__main__":
    dist = sys.argv[1] if len(sys.argv) > 1 else "tauber"
    failing = "failing" in sys.argv[2:]
    if failing:  # second check: failure times, production step cap 2 s and the candidate 5 s, 1800 s exposures
        CELLS = FAIL_CELLS
    cases = []
    for M, z, R in CELLS:
        variants = ((("prod", common.MESH_2D, 2.0), ("dt5", common.MESH_2D, 5.0), ("fine", FINE, 0.5)) if failing
                    else (("prod", common.MESH_2D, 2.0), ("fine", FINE, 0.5)))
        for tag, mesh, dtm in variants:
            cases.append(dict(mach=M, alt=z, R_n=R, catalysis="fc", f_SiC=0.2, t_cap=600.0, stop="all",
                              heat_dist=dist, mesh=mesh, dt_max=dtm, id=f"{dist}_M{M}_z{z / 1e3:.0f}_R{R * 1e3:.0f}_{tag}"))
            if failing:
                cases[-1].update(t_cap=1800.0, stop="first")
    out = HERE / (f"check_mesh_{dist}" + ("_failing" if failing else "") + ".jsonl")
    common.run_all(cases, out, runner=common.run_case_2d)
    rows = {r["case_id"]: r for r in common.load(out)}
    worst_T = worst_t = 0.0
    for M, z, R in CELLS:
        for tag in (("prod", "dt5") if failing else ("prod",)):
            a, b = rows[f"{dist}_M{M}_z{z / 1e3:.0f}_R{R * 1e3:.0f}_{tag}"], rows[f"{dist}_M{M}_z{z / 1e3:.0f}_R{R * 1e3:.0f}_fine"]
            dT = abs(a["T_s_max"] / b["T_s_max"] - 1)
            ta, tb = a["t_fail"], b["t_fail"]
            dt = abs(ta / tb - 1) if (ta == ta and tb == tb) else (0.0 if (ta != ta and tb != tb) else float("inf"))
            if tag == "prod":
                worst_T, worst_t = max(worst_T, dT), max(worst_t, dt)
            print(f"[{tag}] M{M} z{z / 1e3:.0f} R{R * 1e3:.0f}: T_max {a['T_s_max']:.1f} vs {b['T_s_max']:.1f} "
                  f"({100 * dT:.2f} %), t_fail {ta} vs {tb} ({100 * dt:.2f} %), modes {a['mode']}/{b['mode']}, "
                  f"cpu {a['cpu_s']}/{b['cpu_s']} s")
    print(f"worst: T {100 * worst_T:.2f} % (crit 0.5 %), t_fail {100 * worst_t:.2f} % (crit 2 %) ->",
          "PASS" if worst_T <= 0.005 and worst_t <= 0.02 else "FAIL")
