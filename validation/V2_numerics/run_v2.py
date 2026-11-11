"""V2: convergence of the full coupled 1D model (conduction + surface balance + oxide scale) on production
cases.  Refines the first-cell size and growth ratio, the surface-temperature step limit dT_max and the oxide step
limit (5 % -> 2.5 %), and the slab depth, and compares: surface temperature at 10, 60, 300 s, time to failure,
oxide thickness at the end.  Pass criterion: baseline settings within 0.5 % (temperatures, times) and 2 % (thickness)
of the finest settings.  Baseline (production) settings: dx0 = 2 um, growth 1.08, dT_max = 10 K, depth 20 mm.
"""
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from uhtc import coupled as cp, materials as mt  # noqa: E402
import uhtc.oxidation as oxm  # noqa: E402

CASES = [(8.5, 30e3, 0.003, "fc"), (10.0, 30e3, 0.003, "nc"), (9.0, 25e3, 0.001, "fc")]
SETTINGS = {
    "baseline": dict(dx0=2e-6, growth=1.08, dT_max=10.0, depth=0.02, ox=0.05),
    "fine mesh": dict(dx0=0.5e-6, growth=1.04, dT_max=10.0, depth=0.02, ox=0.05),
    "fine time": dict(dx0=2e-6, growth=1.08, dT_max=2.5, depth=0.02, ox=0.0125),
    "finest": dict(dx0=0.5e-6, growth=1.04, dT_max=2.5, depth=0.02, ox=0.0125),
    "deeper slab": dict(dx0=2e-6, growth=1.08, dT_max=10.0, depth=0.05, ox=0.05),
}


def run(case, st):
    M, z, R, cat = case
    ox = oxm.LayerOxidation(oxm.OxParams(f_SiC=0.2))
    orig = ox.rate_limit
    ox.rate_limit = lambda s, T, f, _o=orig, _k=st["ox"] / 0.05: _o(s, T, f) * _k
    cond = cp.ConstantConditions(M, z, R, catalysis=cat)
    t0 = time.time()
    r = cp.run_1d(cond, mt.composite(0.2), ox, 600.0, dx0=st["dx0"], growth=st["growth"], dT_max=st["dT_max"],
                  depth=st["depth"], record_every=0.0)
    tf, mode = cp.first_failure(r.events)
    L = ox.layers(r.states[-1])
    return dict(T10=np.interp(10, r.t, r.T_s), T60=np.interp(60, r.t, r.T_s), T300=np.interp(300, r.t, r.T_s),
                tf=tf, mode=mode, oxide=(L["ZrO2"] + L["glass"]) * 1e6, steps=len(r.t), cpu=time.time() - t0)


def main():
    lines = ["# V2: convergence of the coupled 1D model on production cases", "",
             "Settings and criterion: `run_v2.py` header. Each case runs 600 s (no stop at failure).", ""]
    allok = True
    for case in CASES:
        res = {k: run(case, v) for k, v in SETTINGS.items()}
        ref = res["finest"]
        lines += [f"## Mach {case[0]}, {case[1] / 1e3:.0f} km, R_n = {case[2] * 1e3:.0f} mm, {case[3]}", "",
                  "| Settings | T_s(10 s) | T_s(60 s) | T_s(300 s) | first failure (s, mode) | oxide at 600 s (µm) | steps | CPU (s) |",
                  "|---|---|---|---|---|---|---|---|"]
        for k, v in res.items():
            lines.append(f"| {k} | {v['T10']:.1f} | {v['T60']:.1f} | {v['T300']:.1f} | {v['tf']:.2f} ({v['mode']}) | "
                         f"{v['oxide']:.2f} | {v['steps']} | {v['cpu']:.1f} |")
        b = res["baseline"]
        dev = max(abs(b[k] / ref[k] - 1) for k in ("T10", "T60", "T300"))
        dtf = abs(b["tf"] / ref["tf"] - 1) if b["tf"] == b["tf"] and ref["tf"] == ref["tf"] else (0.0 if b["tf"] != b["tf"] and ref["tf"] != ref["tf"] else 1.0)
        dox = abs(b["oxide"] / ref["oxide"] - 1)
        ok = dev <= 0.005 and dtf <= 0.005 and dox <= 0.02
        allok &= ok
        lines += ["", f"Baseline vs finest: temperatures {dev * 100:.2f} %, failure time {dtf * 100:.2f} %, oxide {dox * 100:.2f} % "
                  f"→ **{'PASS' if ok else 'FAIL'}**.", ""]
        print(lines[-2], flush=True)
    lines.append(f"Overall V2: **{'PASS' if allok else 'FAIL'}**.")
    (HERE / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
