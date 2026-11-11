"""One-command reproduction of the paper's results.

    python reproduce.py <target> [<target> ...]        (on systems with make: make <target>)
    python reproduce.py list                           lists the targets with their cost

Saved outputs are in the repository, so the cheap targets redraw and recheck everything from them. Run scripts skip the cases
already stored in their .jsonl file: to recompute a stage, move that file away first. Costs were measured on 4 workers of an
8-thread Windows laptop; see README.md.
"""
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PY = sys.executable
FIGS = ["fig1", "fig2", "fig3", "fig4", "fig5", "fig6", "fig7", "fig8", "figS1", "figS2", "figS3", "figS4"]
VAL = {"v1": "V1_conduction/run_v1_analytic.py", "v2": "V2_numerics/run_v2.py", "v3": "V3_radiative_equilibrium/run_v3.py",
       "v4": "V4_heating/run_v4.py", "v5": "V5_heating_distribution/run_v5.py", "v6": "V6_oxidation/run_v6.py",
       "v7": "V7_arcjet/run_v7.py", "v8": "V8_measured_flux/run_v8.py", "v9": "V9_plasmatron/run_v9.py"}
AUDIT = ["extract_numbers", "dump_trace", "link_numbers", "numbers_extra", "crosscheck", "v3_numbers", "build_numbers_csv",
         "numbers_extended"]
STAGES_2D = ["grid", "refine", "band", "glass", "active", "kinetics", "flank", "oxidant"]


def py(*args):
    return [PY, *args]


TARGETS = {  # name: (cost, list of commands or of other targets)
    "setup": ("under 1 min, needs internet", [py("literature/janaf/fetch_janaf.py")]),
    "test": ("under 1 min", [py("-m", "pytest", "-q", "tests")]),
    "table1": ("seconds", [py("runs/analysis/make_tables.py", "table1")]),
    "table2": ("seconds", [py("runs/analysis/make_tables.py", "table2")]),
    "table3": ("seconds", [py("runs/analysis/make_tables.py", "table3")]),
    "trajectory": ("seconds", [py("runs/trajectory/make_x43a.py")]),
    "maps1d": ("about 35 min", [py("runs/envelope/run_envelope.py", "A", "B", "C", "D")]),
    "histories": ("about 1 min", [py("runs/histories/run_histories.py"), py("runs/histories/run_jump.py"), py("runs/trajectory/run_x43a.py")]),
    "sensitivity": ("about 4 min", [py("runs/sensitivity/run_sensitivity.py"), py("runs/sensitivity/run_extra.py"), py("runs/sensitivity/run_lewis.py")]),
    "q5": ("about 20 min", [py("runs/q5_2d/run_q5.py"), py("runs/q5_2d/run_q5_sens.py")]),
    "maps2d": ("about 19 h", [py("runs/envelope_2d/run_envelope_2d.py", s) for s in STAGES_2D] + [py("runs/envelope_2d/run_front.py")]),
    "estimates": ("about 10 min", [py("audit/tools/glass_flow_estimate.py"), py("audit/tools/thermal_stress_screen.py"),
                                   py("audit/tools/edge_oxygen.py"), py("audit/tools/flank_ratio.py")]),
    "analysis": ("about 1 min", [py("runs/analysis/summarize.py"), py("runs/analysis/analyze_2d.py")]),
    "trace": ("under 1 min", [py("runs/analysis/trace_numbers.py")]),
    "numbers": ("under 1 min", [py(f"audit/tools/{t}.py") for t in AUDIT]),
    "sample": ("about 25 min", [py("runs/analysis/sample_rerun.py")]),
    "figures": ("about 1 min", ["@figures"]),
    "validation": ("about 13 min", list(VAL)),
    "vvu": ("about 2.5 h with 7 workers (post hoc validation uncertainty, validation/VVU; the committed run took 2 h 51 min, and 4 runs that do not terminate are stopped after 30 min)", [py("validation/VVU/run_vvu.py", "gate"), py("validation/VVU/run_vvu.py", "time"), py("validation/VVU/run_vvu.py", "run", "100"), py("validation/VVU/analyse_vvu.py"), py("validation/VVU/per_quantity.py")]),
    # version 3 (optional modules; the paper's maps are those with every module off, see "regression")
    "regression": ("about 13 min (79 saved cases, modules off, bit-exact)", [py("runs/analysis/regression_v3.py")]),
    "v3checks": ("about 30 min (V10, V2b, V9b, V11, V12)", [py("validation/V10_glassflow/run_v10.py"), py("validation/V2b_flow_numerics/run_v2b.py"),
                 py("validation/V9b_glassflow/beta_estimate.py"), py("validation/V9b_glassflow/run_v9b.py"), py("validation/V11_thermoelastic/run_v11.py"),
                 py("validation/V11_thermoelastic/run_v11_quadratic.py"), py("validation/V12_flank/run_v12.py")]),
    "v3runs": ("about 9 h (glass flow: band columns, regression with flow on, three full maps; flank and geometry columns; stress)",
               [py("runs/envelope_2d/run_envelope_2d.py", "flow"), py("runs/analysis/regression_v3.py", "--glass-flow", "on"),
                py("runs/envelope_2d/run_envelope_2d.py", "flowgrid", "1", "3", "10"), py("runs/envelope_2d/run_envelope_2d.py", "flankint"),
                py("runs/envelope_2d/run_envelope_2d.py", "geometry"), py("runs/thermoelastic/run_stress.py", "all"), py("runs/v3/frozen_shock.py"),
                py("runs/v3/b2o3_advection.py")]),
    "v3analysis": ("about 1 min", [py("runs/analysis/analyze_v3.py"), py("runs/analysis/analyze_flowgrid.py"), py("runs/analysis/flow_numbers.py"),
                   py("figures/make_figures_v3.py")]),
    "check": ("about 5 min", ["setup", "test", "figures", "table2", "table3", "trace", "numbers"]),
    "all": ("about 20 h", ["setup", "test", "trajectory", "maps1d", "histories", "sensitivity", "q5", "maps2d", "validation", "estimates",
                           "analysis", "figures", "trace", "numbers"]),
}
for k, v in VAL.items():
    TARGETS[k] = ("minutes", [py(f"validation/{v}")])
for f in FIGS:
    TARGETS[f] = ("seconds", [f"@{f}"])
# a copy of the repository without the manuscript (paper/) and its number audit has no "numbers" target
if not (ROOT / "audit" / "tools" / "build_numbers_csv.py").exists():
    del TARGETS["numbers"]
    for name in ("check", "all"):
        TARGETS[name] = (TARGETS[name][0], [s for s in TARGETS[name][1] if s != "numbers"])


def draw(names):
    subprocess.run(py("figures/make_figures.py", *names), cwd=ROOT, check=True)
    if not names:  # version-3 Supplement figures S5-S7
        subprocess.run(py("figures/make_figures_v3.py"), cwd=ROOT, check=True)
    if not (ROOT / "paper").is_dir():  # figures stay in figures/out/ when the manuscript is not in the repository
        return
    (ROOT / "paper" / "figures").mkdir(exist_ok=True)
    for p in [*(ROOT / "figures" / "out").glob("*.pdf"), *(ROOT / "figures" / "out").glob("*.png")]:
        shutil.copy2(p, ROOT / "paper" / "figures" / p.name)


def run(target):
    if target not in TARGETS:
        raise SystemExit(f"unknown target {target!r}; try: python reproduce.py list")
    t0 = time.time()
    for step in TARGETS[target][1]:
        if isinstance(step, str) and step.startswith("@"):
            draw([] if step == "@figures" else [step[1:]])
        elif isinstance(step, str):
            run(step)
        else:
            print("+", " ".join(Path(a).name if a == PY else a for a in step), flush=True)
            subprocess.run(step, cwd=ROOT, check=True, env={**os.environ, "PYTHONUTF8": "1"})
    print(f"[{target}: {time.time() - t0:.0f} s]", flush=True)


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] == "list":
        print(__doc__)
        for k, (cost, _) in TARGETS.items():
            print(f"  {k:12s} {cost}")
    else:
        for t in sys.argv[1:]:
            run(t)
