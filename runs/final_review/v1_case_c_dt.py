"""Final review (8 Nov 2026, R3-01): V1 case C (2D sphere and cylinder under a uniform flux) at a time step of
0.25 ms instead of 1 ms, so that the observed spatial order of the 2D solver is not masked by the time error.
Uses the functions of validation/V1_conduction/run_v1_analytic.py unchanged.  Writes v1_case_c_dt.json."""
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("v1", ROOT / "validation" / "V1_conduction" / "run_v1_analytic.py")
v1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v1)

T_END, F, A = 0.12, 2e6, 3e-3  # as in V1 case C (Fo = 0.22)
res = {}
for axi, exact, name in ((True, v1.exact_sphere_surface(T_END, F, A), "sphere"), (False, v1.exact_cyl_surface(T_END, F, A), "cylinder")):
    for dt in (1e-3, 2.5e-4):
        errs = []
        for n_arc, n_eta in ((20, 10), (40, 20), (80, 40)):
            Tm, spread, r, Tbar = v1.run_sphere2d(axi, n_arc, n_eta, dt, T_END, F=F, a=A)
            errs.append(abs(Tm - exact) / (exact - v1.T0))
        res[f"{name}_dt{dt:g}"] = dict(errors=errs, order=[float(x) for x in v1.order(errs)])
(Path(__file__).resolve().parent / "v1_case_c_dt.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
print(json.dumps(res, indent=1))
