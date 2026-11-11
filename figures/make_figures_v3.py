"""Version-3 figures (Supplement), drawn from saved outputs with the style of make_figures.py (no model runs here).

Usage: python figures/make_figures_v3.py [figS5 figS6 figS7]     (default: all)
Output: figures/out/<name>.pdf and .png
  figS5_flank_models  three flank-heating models: ratios to the baseline at the map corners, and the measured cone of V5
                      (validation/V12_flank/v12_results.json)
  figS6_glass_flow    fully catalytic 2D boundary with and without glass flow, and the glass-loss time with flow relative to
                      without at the old lowest failing Mach (runs/analysis/flowgrid_boundaries.csv)
  figS7_stress        largest tensile stress within 2 R_n relative to the bend strength, cold start at the five conditions and
                      the X-43A path of the 1 mm tip (runs/thermoelastic/stress.jsonl)
"""
import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_figures as MF  # noqa: E402  (style, colours, save)
from matplotlib.ticker import NullFormatter  # noqa: E402

plt = MF.plt
ROOT = MF.ROOT
C, INK2 = MF.C, MF.INK2
RADII = (1.0, 3.0, 10.0)


def figS5_flank_models():
    d = json.loads((ROOT / "validation/V12_flank/v12_results.json").read_text(encoding="utf-8"))
    fig, ax = plt.subplots(1, 2, figsize=(MF.W2, 2.4))
    ov = [r for r in d["overlay"] if r["point"] != "tangency"]
    for j, (m, lab) in enumerate((("integral", "Integral"), ("sharpcone", "Sharp cone"))):
        for R, mk in ((1.0, "o"), (10.0, "s")):
            sel = [r for r in ov if r["R_mm"] == R]
            ax[0].scatter([r["mach"] + (0.12 if j else -0.12) for r in sel], [r[m] / r["baseline"] for r in sel], marker=mk, s=14,
                          facecolors="none" if R == 10.0 else C[j], edgecolors=C[j], lw=0.8,
                          label=f"{lab}, {R:g} mm")
    ax[0].axhspan(0.5, 2.0, color=MF.BAND, zorder=0, lw=0)
    ax[0].axhline(1.0, color=INK2, lw=0.8)
    ax[0].set_yscale("log")
    ax[0].set_ylim(0.28, 3.6)  # room for the legend below the data
    ax[0].set_yticks([0.5, 1, 2, 3])
    ax[0].set_yticklabels(["0.5", "1", "2", "3"])
    ax[0].yaxis.set_minor_formatter(NullFormatter())
    ax[0].set_xticks([6, 9, 12])
    ax[0].set_xlabel("Mach number (20, 30, 40 km)")
    ax[0].set_ylabel("Flank heating relative\nto the baseline")
    ax[0].legend(loc="lower right", fontsize=7, frameon=False, handletextpad=0.2, borderaxespad=0.2)
    an = sorted(d["anchor"], key=lambda a: a["s_over_R"])
    x = [a["s_over_R"] for a in an]
    ax[1].scatter(x, [a["measured"] for a in an], marker="D", s=12, color=MF.INK, label="Measured (V5)", zorder=3)
    for j, (m, lab, ls) in enumerate((("baseline", "Baseline", "-"), ("integral", "Integral method", "--"),
                                      ("sharpcone", "Sharp-cone correlation", ":"))):
        ax[1].scatter(x, [a[m] for a in an], marker="o", s=10, facecolors="none", edgecolors=[INK2, C[0], C[1]][j], lw=0.8, label=lab)
    ax[1].set_yscale("log")
    ax[1].set_yticks([0.05, 0.1, 0.2])
    ax[1].set_yticklabels(["0.05", "0.1", "0.2"])
    ax[1].yaxis.set_minor_formatter(NullFormatter())
    ax[1].set_xlabel("Surface distance $s/R_\\mathrm{n}$")
    ax[1].set_ylabel("$q/q_\\mathrm{s}$")
    ax[1].legend(loc="upper right", fontsize=7, frameon=False, handletextpad=0.2)
    for a, t in zip(ax, ("(a)", "(b)")):
        MF.grid(a)
        a.set_title(t, loc="left")
    fig.tight_layout()
    MF.save(fig, "figS5_flank_models")


def figS6_glass_flow():
    rows = list(csv.DictReader(open(ROOT / "runs/analysis/flowgrid_boundaries.csv", encoding="utf-8")))
    fig, ax = plt.subplots(1, 2, figsize=(MF.W2, 2.5))
    for i, R in enumerate(RADII):
        rr = sorted((r for r in rows if float(r["R_n_mm"]) == R), key=lambda r: float(r["alt_km"]))
        z = [float(r["alt_km"]) for r in rr]
        ax[0].step([float(r["M_s_base"]) for r in rr], z, where="mid", color=C[i], lw=1.4, label=f"{R:g} mm")
        ax[0].step([float(r["M_s_flow"]) for r in rr], z, where="mid", color=C[i], lw=1.0, ls="--")
        g = [(float(r["alt_km"]), float(r["t_fail_flow_at_old"]) / float(r["t_fail_base"])) for r in rr
             if r["t_fail_flow_at_old"] not in ("", "None") and r["mode_base"] == "glass" and r["mode_flow_at_old"] == "glass"]
        ax[1].plot([b for a, b in g], [a for a, b in g], marker="o", ms=3, color=C[i], lw=1.0, label=f"{R:g} mm")
    ax[0].set_xlabel("30-minute boundary (Mach)")
    ax[0].set_ylabel("Altitude (km)")
    from matplotlib.lines import Line2D
    h, lab = ax[0].get_legend_handles_labels()
    h += [Line2D([], [], color=INK2, lw=1.4), Line2D([], [], color=INK2, lw=1.0, ls="--")]
    lab += ["no flow", "flow"]
    ax[0].legend(h, lab, loc="lower right", fontsize=7, frameon=False, ncol=1, handlelength=1.8, borderaxespad=0.2)
    ax[1].axvline(1.0, color=INK2, lw=0.8)
    ax[1].set_xlabel("Glass-loss time with flow\nrelative to without")
    ax[1].set_xlim(0.6, 1.02)
    ax[1].legend(loc="lower left", fontsize=7, frameon=False)
    for a, t in zip(ax, ("(a)", "(b)")):
        MF.grid(a)
        a.set_title(t, loc="left")
        a.set_ylim(19.5, 40.5)
    fig.tight_layout()
    MF.save(fig, "figS6_glass_flow")


def figS7_stress():
    rows = [json.loads(l) for l in open(ROOT / "runs/thermoelastic/stress.jsonl", encoding="utf-8") if l.strip()]
    cases = [(8.5, 20.0, 1.0), (9.25, 30.0, 1.0), (10.25, 40.0, 1.0), (9.75, 30.0, 3.0), (10.25, 30.0, 10.0)]
    fig, axs = plt.subplots(1, 2, figsize=(MF.W2, 2.6), sharey=True)
    for i, (M, z, R) in enumerate(cases):
        rr = {r["t_s"]: r for r in rows if r["kind"] == "case" and (r["mach"], r["alt_km"], r["R_mm"]) == (M, z, R)}
        rr = [rr[t] for t in sorted(rr)]  # one row per snapshot time (the mesh-check rows repeat two of them)
        t = [r["t_s"] for r in rr]
        axs[0].plot(t, [r["free_TPSX"]["ratio_max_2R"] for r in rr], color=C[i % len(C)], lw=1.2, marker="o", ms=2.5,
                    label=f"Mach {M:g}, {z:g} km, {R:g} mm")
        axs[0].plot(t, [r["free_Neuman"]["ratio_max_2R"] for r in rr], color=C[i % len(C)], lw=0.8, ls="--")
        axs[1].plot(t, [r["restrained_TPSX"]["ratio_max_body"] for r in rr], color=C[i % len(C)], lw=1.2, marker="o", ms=2.5)
    xr = sorted((r for r in rows if r["kind"] == "x43a"), key=lambda r: r["t_s"])
    axs[0].plot([r["t_s"] for r in xr], [r["free_TPSX"]["ratio_max_2R"] for r in xr], color=MF.INK, lw=1.0, ls=":", marker="s", ms=2.5,
                label="X-43A path, 1 mm")
    axs[0].set_ylabel("Largest tensile stress\n/ bend strength")
    # legend below both panels, clear of every curve (journal revision; inside panel (a) it covered the 10 mm curve)
    h, lab = axs[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=7, frameon=False, handlelength=1.6,
               columnspacing=1.2)
    for a, t in zip(axs, ("(a) within $2R_\\mathrm{n}$, free tip", "(b) whole body, restrained tip")):
        a.set_xscale("log")
        a.set_xlabel("Time after the start (s)")
        a.set_ylim(bottom=0)
        MF.grid(a)
        a.set_title(t, loc="left")
    fig.tight_layout()
    MF.save(fig, "figS7_stress")


FIGS = dict(figS5=figS5_flank_models, figS6=figS6_glass_flow, figS7=figS7_stress)

if __name__ == "__main__":
    for k in sys.argv[1:] or FIGS:
        FIGS[k]()
    import shutil
    if (ROOT / "paper/figures").is_dir():  # the copies go to the manuscript only when it is in the repository
        for k in sys.argv[1:] or FIGS:
            name = FIGS[k].__name__
            for ext in ("pdf", "png"):
                shutil.copyfile(MF.OUT / f"{name}.{ext}", ROOT / "paper/figures" / f"{name}.{ext}")
