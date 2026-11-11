"""Draw every figure of the paper from saved outputs (no model runs here).

Usage: python figures/make_figures.py [fig1 fig2 ...]     (default: all)
Output: figures/out/<name>.pdf and .png (and copies used by paper/figures/).
Inputs: validation/*/ results, runs/*/ results (see each function).
Colours: categorical slots in fixed order (blue, orange, aqua, yellow; validated for colour-vision deficiency with
the dataviz validator - aqua and yellow need direct labels or line styles, which every figure provides);
sequential maps use one blue ramp, light = short time to failure, dark = long.
Conventions (STYLE_RULES.md, Section G), the same in every figure: fully catalytic solid, non-catalytic dashed; thick
lines the 2D nose-tip model, thin lines the 1D stagnation-line model; failure times are dotted vertical lines; censored
values (no failure within the 1800 s cap, or none up to Mach 12) are drawn as open triangles pointing the way the true
value lies. Figure text is at least 7 pt at the printed width (the \\includegraphics widths in paper/main.tex); only
the sub- and superscripts of symbols and units are smaller (mathtext scales them by 0.7).
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap, to_rgb  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures" / "out"
OUT.mkdir(parents=True, exist_ok=True)
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#6250d6", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
BAND = "#e4e2dd"  # shading of a tolerance or uncertainty band: visible on white in print
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
LS = ["-", "--", "-.", ":"]
W1, W2 = 3.4, 5.3  # drawn at the printed width (one column, or 5.3 in for the wide and the S figures): sizes below are printed sizes, none under 7 pt
TS_LABEL = "Stagnation surface\ntemperature $T_\\mathrm{s}$ (K)"  # the text's term for T_s, one label in every figure

plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "legend.fontsize": 7, "hatch.linewidth": 0.9,
                     "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.edgecolor": INK2, "axes.linewidth": 0.6,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.labelcolor": INK, "text.color": INK,
                     "lines.linewidth": 1.4, "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
                     "savefig.dpi": 300, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42})  # TrueType instead of Type 3 fonts in the PDFs (final review R8p2-18)


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", metadata={"CreationDate": None})  # no time stamp: reruns are byte-identical
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", dpi=200)
    plt.close(fig)
    print("wrote", name)


def grid(ax):
    ax.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)


def tint(col, f=0.55):
    """Lighter tint of a colour (mixed with white), for a second quantity drawn in the same colour."""
    return tuple(1 - f * (1 - v) for v in to_rgb(col))


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(l for l in f if not l.startswith("#")))


def fnum(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return math.nan
    return v


# ---------------------------------------------------------------------------------------------------- Fig. 1
def fig1_schematic():
    """Model schematic: sharp nose with heating distribution (left) and the layer stack with energy terms (right)."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(W2, 2.5), gridspec_kw=dict(width_ratios=[0.86, 1.14], wspace=0.04))
    # left: sphere-cone nose
    Rn, th = 1.0, np.radians(8)
    phi_t = np.pi / 2 - th
    phi = np.linspace(0, phi_t, 60)
    xa, ya = Rn * (1 - np.cos(phi)), Rn * np.sin(phi)
    L = 3.3
    xf = np.r_[xa[-1], xa[-1] + L * np.cos(th)]
    yf = np.r_[ya[-1], ya[-1] + L * np.sin(th)]
    for sgn in (1, -1):
        a.plot(xa, sgn * ya, color=INK, lw=1.2)
        a.plot(xf, sgn * yf, color=INK, lw=1.2)
    a.fill_between(np.r_[xa, xf[1:]], -np.r_[ya, yf[1:]], np.r_[ya, yf[1:]], color="#e9e8e4")
    a.plot([xf[1], xf[1]], [-yf[1], yf[1]], color=INK2, lw=0.8, ls="--")
    a.text(xf[1] + 0.12, 0, "adiabatic back plane", fontsize=7, va="center", ha="left", rotation=90, color=INK2)
    for p_ in np.linspace(0, np.radians(70), 6):
        x0, y0 = Rn * (1 - np.cos(p_)), Rn * np.sin(p_)
        ln = 0.9 * np.cos(p_) + 0.15
        a.annotate("", xy=(x0, y0), xytext=(x0 - ln * np.cos(p_), y0 + ln * np.sin(p_)),
                   arrowprops=dict(arrowstyle="->", color=C[1], lw=0.9))
    # the three parts of g (Section 2.2): cos to 70 deg, linear bridge to the tangency, Lees on the flank
    a.text(-1.5, 3.02, r"heating: $g=\cos\varphi$ to 70°,", color=C[1], fontsize=7)
    a.text(-1.5, 2.55, "linear bridge to the tangency,", fontsize=7, color=C[1])
    a.text(-1.5, 2.08, "Lees local similarity on the flank", fontsize=7, color=C[1])
    a.annotate("", xy=(0, 0), xytext=(-1.5, 0), arrowprops=dict(arrowstyle="->", color=INK2, lw=0.8))
    a.text(-1.5, -0.42, r"$M$", fontsize=7, color=INK2)
    a.plot([Rn, Rn], [0, 0], "o", ms=2, color=INK2)
    a.annotate(r"$R_\mathrm{n}$", xy=(Rn * (1 - np.cos(1.0)), Rn * np.sin(1.0)), xytext=(Rn + 0.1, 0.42),
               arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6), fontsize=7, color=INK2, va="center")
    a.annotate("1D stagnation-line model", xy=(0.6, 0.0), xytext=(-1.5, -2.2), fontsize=7, color=C[0],
               arrowprops=dict(arrowstyle="->", color=C[0], lw=0.8))
    a.plot([0, 2.0], [0, 0], color=C[0], lw=2.0)
    a.set_xlim(-1.7, 4.9)
    a.set_ylim(-2.5, 3.45)
    a.set_aspect("equal")
    a.axis("off")
    a.set_title("(a) stagnation region, heating", loc="left")
    # right: layer stack (no semicolons in the artwork, rule 4)
    layers = [("outer glass\n(B$_2$O$_3$–SiO$_2$)", "#cde2fb", 0.62, "evaporation:\nB$_2$O$_3$(g), SiO(g), SiO$_2$(g)"),
              ("ZrO$_2$ skeleton + glass", "#f3d2c4", 0.72, r"growth $\mathrm{d}(x_\mathrm{r}^2)/\mathrm{d}t=k_\mathrm{p}$," "\nO$_2$ inward"),
              ("SiC-depleted ZrB$_2$", "#d9d8d4", 0.62, r"$\mathrm{d}(\delta_\mathrm{d}^2)/\mathrm{d}t=k_\mathrm{d}$"),
              ("ZrB$_2$–SiC substrate\n(1D or 2D conduction)", "#b9c6d6", 0.95, "scale: thermal resistance\nin series, clocks (a)–(c)")]
    y = 3.05
    # gas state in words: no nested subscripts (p_O2 printed at 3.4 pt)
    b.text(0.03, 3.42, "gas: convective heating,\nO$_2$ pressure", fontsize=7, va="center")
    for txt, col, h, note in layers:
        b.add_patch(plt.Rectangle((0, y - h), 1.6, h, color=col, lw=0))
        b.text(0.04, y - h / 2, txt, fontsize=7, va="center")
        b.text(1.66, y - h / 2, note, fontsize=7, va="center", color=C[4] if "evap" in note else INK2)
        y -= h
    b.annotate("", xy=(1.86, 3.08), xytext=(1.86, 3.75), arrowprops=dict(arrowstyle="->", color=C[1]))
    b.text(1.9, 3.42, r"$q_\mathrm{conv}$", fontsize=8, color=C[1], va="center")
    b.annotate("", xy=(2.58, 3.75), xytext=(2.58, 3.08), arrowprops=dict(arrowstyle="->", color=C[2]))
    b.text(2.62, 3.42, r"$ϵ\sigma T_\mathrm{s}^4$", fontsize=8, color=C[2], va="center")  # lunate epsilon, as in Eq. (4)
    b.set_xlim(0, 3.3)
    b.set_ylim(0.1, 3.8)
    b.axis("off")
    b.set_title("(b) surface balance, layers (not to scale)", loc="left")
    save(fig, "fig1_model")


# ---------------------------------------------------------------------------------------------------- Fig. S1
def fig2_heating_checks():
    """Supplementary Figure S1, panels in the order of the checks: (a) V3, (b) V4, (c) V5."""
    v4 = read_csv(ROOT / "validation" / "V4_heating" / "v4_results.csv")
    v3 = read_csv(ROOT / "validation" / "V3_radiative_equilibrium" / "v3_results.csv")
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(5.1, 3.05))
    # (a) V3: velocity at 2149 C against Kolodziej. Open (published) symbols drawn larger and on top, so that the
    # filled model symbols inside them stay visible where the two agree.
    for src, lab, mk in (("Table1_FCHW", "tabulated, fully catalytic", "s"), ("Fig1_FCHW", "digitised, fully catalytic", "o"),
                         ("Fig1_NCHW", "digitised, non-catalytic", "^")):
        rows = [r for r in v3 if r["source"] == src]
        z = [fnum(r["alt_km"]) for r in rows]
        col = C[1] if "non" in lab else C[0]
        a.plot([fnum(r["V_model_SG_kms"]) for r in rows], z, mk, ms=2.6, color=col, zorder=3)
        a.plot([fnum(r["V_kolodziej_kms"]) for r in rows], z, mk, ms=6.2, mfc="none", mew=0.9, color=col, zorder=4, label=lab)
    a.set_xlabel("Velocity at 2422 K (km s$^{-1}$)")  # journal revision: kelvin
    a.set_ylabel("Altitude (km)")
    a.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(-0.3, -0.3), title="open: Kolodziej\nfilled: this model",
             title_fontsize=7, labelspacing=0.25, borderaxespad=0.1, alignment="left")
    grid(a)
    a.set_title("(a) V3", loc="left")
    # (b) V4: Sutton-Graves against Fay-Riddell (method names in plain text: no 5.7 pt subscripts)
    for i, Tw in enumerate(("300", "1111", "2000")):
        rows = [r for r in v4 if r["T_w"] == Tw]
        b.plot([fnum(r["H_s_MJkg"]) for r in rows], [fnum(r["SG_over_FR"]) for r in rows], "o", ms=4.5,
               color=C[i], mfc="none" if i == 1 else C[i], label=rf"$T_\mathrm{{s}}$ = {Tw} K")
    b.axhspan(0.85, 1.15, color=BAND, zorder=0)
    b.axhline(1.0, color=INK2, lw=0.7)
    b.set_xlabel("$H_\\mathrm{s}$ (MJ kg$^{-1}$)")
    b.set_ylabel("Heating ratio,\nSutton–Graves / Fay–Riddell")
    b.set_ylim(0.8, 1.2)
    b.text(6.6, 1.16, "±15 % band", fontsize=7, color=INK2, ha="right", va="bottom")
    b.legend(frameon=False, loc="lower left", labelspacing=0.2, borderaxespad=0.2)  # inside, in the empty lower left
    grid(b)
    b.set_title("(b) V4", loc="left")
    # (c) V5: distribution around a 15 deg sphere-cone (Cleary 1969, measured) vs the baseline and its band
    cur = read_csv(ROOT / "validation" / "V5_heating_distribution" / "v5_curve.csv")
    s = np.array([fnum(r["s_over_R"]) for r in cur])
    q = np.array([fnum(r["q_over_q0"]) for r in cur])
    beyond = s > np.radians(70.0)
    hb = c.fill_between(s[beyond], 0.66 * q[beyond], 1.5 * q[beyond], color=BAND, lw=0, zorder=0)
    hl, = c.plot(s, q, color=INK, lw=1.0)
    v5b = [r for r in read_csv(ROOT / "validation" / "V5_heating_distribution" / "v5_results.csv")
           if r["case"].startswith("Cleary")]
    hm = []
    for i, (re, mk) in enumerate((("0.4e6", "o"), ("1.2e6", "s"), ("1.8e6", "^"))):
        rr = [r for r in v5b if r["Re_per_ft"] == re]
        per_m = float(re) / 0.3048 / 1e6  # unit Reynolds number per foot (Cleary) given per metre
        h, = c.plot([fnum(r["s_over_R"]) for r in rr], [fnum(r["measured"]) for r in rr], mk, ms=3.6, color=C[i],
                    mfc="none" if i else C[i], lw=0, label=f"{per_m:.1f} million")
        hm.append(h)
    c.set_yscale("log")
    c.set_xlim(0, 17)
    c.set_ylim(1e-2, 1.3)
    c.set_xlabel(r"$s/R_\mathrm{n}$")
    c.set_ylabel("Heating distribution $g$")
    leg1 = c.legend([hl, hb], ["baseline", "band ×0.66–1.5"], frameon=False, fontsize=7, loc="lower left", labelspacing=0.2,
                    borderaxespad=0.2, handletextpad=0.4)
    c.add_artist(leg1)
    c.legend(handles=hm, frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(-0.3, -0.3),
             title="Cleary, measured, Reynolds\nnumber per metre:", title_fontsize=7, labelspacing=0.25, borderaxespad=0.1,
             handletextpad=0.3, alignment="left")
    grid(c)
    c.set_title("(c) V5, 15° cone", loc="left")
    fig.tight_layout(w_pad=1.2)
    save(fig, "figS1_heating_checks")


# ---------------------------------------------------------------------------------------------------- Fig. 2
def fig3_furnace():
    """Check V6: model against measured furnace thicknesses. Depleted-layer points have their own marker; the two
    depleted-layer measurements of zero thickness (not on log axes) are drawn at the left edge; glass points are
    reported but not counted in V6."""
    rows = read_csv(ROOT / "validation" / "V6_oxidation" / "v6_results.csv")
    fig, ax = plt.subplots(figsize=(W1, 2.45))
    is_dep = lambda r: r["quantity"].startswith("SiC-depleted")  # noqa: E731
    pos = lambda r: fnum(r["measured_um"]) > 0  # noqa: E731
    groups = [("calibration set (Shugart)", "s", C[0], "none", lambda r: r["set"] == "calibration"),
              ("independent: oxide, ZrO$_2$ layer, recession", "o", C[1], C[1], lambda r: r["set"] == "independent" and not is_dep(r)),
              ("independent: SiC-depleted layer", "v", C[1], C[1], lambda r: r["set"] == "independent" and is_dep(r)),
              ("glass thickness (excluded from V6)", "^", C[2], "none", lambda r: r["set"].startswith("glass"))]
    for lab, mk, col, mfc, sel_ in groups:
        sel = [r for r in rows if sel_(r) and pos(r)]
        ax.loglog([fnum(r["measured_um"]) for r in sel], [fnum(r["model_um"]) for r in sel], mk, ms=4.5 if mk != "v" else 5.0,
                  color=col, mfc=mfc, label=lab, lw=0)
    # incubation option (D-OX-10): only the points whose model value changes with the option are drawn
    inc = [r for r in rows if r["set"] == "independent" and pos(r)
           and abs(fnum(r["model_inc_um"]) - fnum(r["model_um"])) > 1e-6 and fnum(r["model_inc_um"]) > 0]
    ax.loglog([fnum(r["measured_um"]) for r in inc], [fnum(r["model_inc_um"]) for r in inc], "D", ms=3.2,
              color=C[3], mfc="none", lw=0, label="with incubation option (changed points)")
    # depleted layers measured as zero: no place on a log axis, drawn at the left edge with the model value
    zero = [r for r in rows if r["set"] == "independent" and is_dep(r) and not pos(r)]
    ax.plot([3.0] * len(zero), [fnum(r["model_um"]) for r in zero], "<", ms=5.0, color=C[1], mfc="white", mew=0.9, lw=0,
            clip_on=False, zorder=5, label="SiC-depleted, measured as zero (left edge)")
    x = np.array([1, 5000])
    ax.plot(x, x, color=INK, lw=0.8)
    ax.plot(x, 3 * x, color=INK2, lw=0.7, ls=":")
    ax.plot(x, x / 3, color=INK2, lw=0.7, ls=":")
    ax.text(4.5, 4.5 * 4.6, "×3", fontsize=7, color=INK2)  # above the dotted line, not on it (R8p2-04)
    ax.text(300, 300 / 4.3, "÷3", fontsize=7, color=INK2)
    ax.text(1050, 1500, "1:1", fontsize=7, color=INK, ha="right")
    ax.set_xlim(3, 2000)
    ax.set_ylim(3, 2000)
    ax.set_xlabel("Measured thickness (µm)")
    ax.set_ylabel("Model thickness (µm)")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.42, -0.2), ncol=1, handletextpad=0.3,
              labelspacing=0.2)  # below the axis label: clear of the data; one column keeps the column width
    grid(ax)
    save(fig, "fig2_furnace")


# ---------------------------------------------------------------------------------------------------- Fig. 3
def fig4_arcjet():
    res = json.load(open(ROOT / "validation" / "V7_arcjet" / "v7_results.json", encoding="utf-8"))
    meas = [r for r in read_csv(ROOT / "literature" / "digitized" / "Savino2010_TOAEJ_ArcJetUHTC_fig8.csv") if r["series"] == "hemisphere"]
    t = np.array([fnum(r["t_s"]) for r in meas])
    T = np.array([fnum(r["T_K"]) for r in meas])
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(W2, 2.5), gridspec_kw=dict(width_ratios=[1.15, 1]))
    sel = T > 1285
    ax.plot(t[sel], T[sel], ".", color=INK, ms=2.5, label="measured")
    # Journal revision (Phase 8): 1D and 2D differ by shade and marker as well as by line width, so that the panel reads in
    # greyscale print. 1D: light tint, thin line, open circles. 2D: full colour, thick line, no markers.
    show = [("1d-fc", C[1], "-", "1D, fully catalytic"), ("1d-nc", C[0], "--", "1D, non-catalytic"),
            ("2d-fc", C[1], "-", "2D, fully catalytic"), ("2d-nc", C[0], "--", "2D, non-catalytic")]
    for k, col, ls, lab in show:
        if k in res:
            if "2d" in k:
                ax.plot(res[k]["t"], res[k]["Ts"], color=col, ls=ls, lw=1.7, label=lab)
            else:
                tt = np.asarray(res[k]["t"])
                every = [int(i) for i in np.searchsorted(tt, np.arange(60, 560, 60))]
                ax.plot(tt, res[k]["Ts"], color=tint(col, 0.6), ls=ls, lw=0.9, marker="o", ms=3.0, mfc="white",
                        mec=tint(col, 0.6), mew=0.8, markevery=every, label=lab)
    ax.set_xlim(40, 560)
    ax.set_ylim(1200, 3000)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(TS_LABEL)
    ax.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=2)  # below the axis label: clear of the data
    grid(ax)
    ax.set_title("(a) V7: SPES arc jet", loc="left")
    # (b) V8: measured cold-wall flux; model minus measured surface temperature (1D, both catalysis bounds)
    rb = json.load(open(ROOT / "validation" / "V8_measured_flux" / "v8_results.json", encoding="utf-8"))
    cases = [(f"Ames, {h} MJ kg$^{{-1}}$", lambda c, h=h: rb[f"V8_1_H{h}_{c}"]["passive"]["T_max"], 1630.0)
             for h in (15, 30, 45)]
    cases += [("HIT, condition 1", lambda c: rb[f"V8_2_{c}"]["c1"]["T_end"], 1960.0),
              ("HIT, condition 2", lambda c: rb[f"V8_2_{c}"]["c2"]["T_end"], 1900.0)]
    bx.axvspan(-80, 80, color="#f0efec", zorder=0)
    bx.axvline(0, color=INK2, lw=0.7)
    for j, (cat, col, mk, lab) in enumerate((("fc", C[1], "o", "fully catalytic"), ("nc", C[0], "s", "non-catalytic"))):
        y = np.arange(len(cases)) + (0.12 if j else -0.12)
        bx.plot([f(cat) - 273.15 - Tm for _, f, Tm in cases], y, mk, ms=4.5, color=col, mfc=col if j == 0 else "none",
                lw=0, label=lab)
    bx.set_yticks(range(len(cases)))
    bx.set_yticklabels([c[0] for c in cases], fontsize=7)
    bx.invert_yaxis()
    bx.set_xlabel(r"Model − measured $T_\mathrm{s}$ (K)")
    bx.set_ylim(len(cases) + 0.5, -0.5)
    bx.text(90, len(cases) + 0.05, "±80 K band", fontsize=7, color=INK2, ha="left", va="center")  # empty bottom row
    bx.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.45, -0.3), ncol=2)  # below the axis label: clear of the band
    grid(bx)
    bx.set_title("(b) V8: measured flux", loc="left")
    fig.tight_layout(w_pad=1.0)
    save(fig, "fig3_arcjet")


# ---------------------------------------------------------------------------------------------------- data helpers
def load_jsonl(path):
    return [json.loads(l) for l in open(path, encoding="utf-8") if l.strip()]


def x43a_path():
    a = np.loadtxt(ROOT / "runs" / "trajectory" / "x43a_mach10.csv", delimiter=",")
    return a[:, 0], a[:, 1] / 1e3, a[:, 2]


T_LIMIT = 2422.15
BINS = [0, 10, 30, 100, 300, 1000, 1800, 1e9]
BIN_LABELS = ["<10 s", "10–30 s", "30–100 s", "100–300 s", "300–1000 s", "1000–1800 s", "≥1800 s"]
BIN_TICKS = [5, 20, 65, 200, 650, 1400, 5e8]  # mid-bin under the uniform colour-bar spacing (the last bin is [1800, 1e9])


# ---------------------------------------------------------------------------------------------------- Fig. 4
def fig5_temperature():
    """Q1: stagnation temperature histories - constant condition (Mach 9, 30 km) and the X-43A flight."""
    H = read_csv(ROOT / "runs" / "histories" / "histories.csv")
    E = read_csv(ROOT / "runs" / "histories" / "events.csv")
    fig = plt.figure(figsize=(W2, 2.55))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.9, 2.5], hspace=0.08, wspace=0.62)  # room for the (b) y-label
    a = fig.add_subplot(gs[:, 0])
    for i, R in enumerate(("1.0", "3.0", "10.0")):
        for cat, ls in (("fc", "-"), ("nc", "--")):
            rows = [r for r in H if r["mach"] == "9.0" and r["R_n_mm"] == R and r["catalysis"] == cat]
            a.plot([fnum(r["t"]) for r in rows], [fnum(r["T_s"]) for r in rows], color=C[i], ls=ls, lw=1.2,
                   label=rf"$R_\mathrm{{n}}$ = {float(R):g} mm" if cat == "fc" else None)
    a.axhline(T_LIMIT, color=INK2, lw=0.8, ls=":")
    a.text(2.3, T_LIMIT + 30, "limit (c)", ha="left", fontsize=7, color=INK2)  # short, so it ends before the curves cross the line
    Ttr = fnum([e for e in E if e["mach"] == "9.0"][0]["T_transition"])
    a.axhline(Ttr, color=C[4], lw=0.8, ls=":")
    a.text(2.3, Ttr + 45, rf"SiC active above {Ttr:.0f} K", ha="left", fontsize=7, color=C[4])  # clear of its own line
    a.set_xscale("log")
    a.set_xlim(2, 1800)
    a.set_ylim(300, 3100)
    a.set_xlabel("Time (s)")
    a.set_ylabel(TS_LABEL)
    a.legend(frameon=True, framealpha=1, edgecolor="none", loc="lower right", title="solid: fully catalytic\ndashed: non-catalytic",
             title_fontsize=7, labelspacing=0.2, borderaxespad=0.1)  # opaque: clear of the 10 mm curve (R8p2-04)
    grid(a)
    a.set_title("(a) Mach 9, 30 km, constant", loc="left")
    tt, zz, MM = x43a_path()
    b0 = fig.add_subplot(gs[0, 1])
    b0.plot(tt, MM, color=INK2, lw=1.0)
    b0.set_ylabel("Mach")
    b0.set_xlim(30, 640)
    b0.tick_params(labelbottom=False)
    grid(b0)
    b0r = b0.twinx()  # labelled tick axis only for altitude would be a dual axis: avoided, altitude given in text
    b0r.set_axis_off()
    b0.set_title("(b) X-43A Mach 10 flight path", loc="left")
    b1 = fig.add_subplot(gs[1, 1], sharex=b0)
    Q = [r for r in read_csv(ROOT / "runs" / "q5_2d" / "q5_histories.csv") if r["kind"] == "x43a" and r["catalysis"] == "fc"]
    for i, R in enumerate(("1.0", "3.0", "10.0")):
        for model, ls, lw in (("2d", "-", 1.7), ("1d", "-", 0.7)):
            rows = [r for r in Q if r["R_n_mm"] == R and r["model"] == model]
            b1.plot([fnum(r["t"]) + tt[0] for r in rows], [fnum(r["T_s"]) for r in rows], color=C[i], lw=lw, ls=ls,
                    label=rf"$R_\mathrm{{n}}$ = {float(R):g} mm" if model == "2d" else None)
    b1.axhline(T_LIMIT, color=INK2, lw=0.8, ls=":")
    b1.set_xlabel("Flight time (s)")
    b1.set_ylabel(TS_LABEL)
    b1.set_ylim(300, 3100)
    b1.text(630, 2330, "fully catalytic\nthick: 2D\nthin: 1D", ha="right", va="top", fontsize=7, color=INK2)
    grid(b1)
    save(fig, "fig4_temperature")


# ---------------------------------------------------------------------------------------------------- Fig. 5
def fig6_layers():
    """Q2: glass, ZrO2 and SiC-depleted thicknesses vs time for three illustrative cases (top), and the glass thickness
    against the critical thickness of criterion (a) and the fixed sensitivity thresholds (bottom; Part 2.4)."""
    H = read_csv(ROOT / "runs" / "histories" / "histories.csv")
    E = read_csv(ROOT / "runs" / "histories" / "events.csv")
    # each title names its catalysis bound: the first two cases are fully catalytic, the third non-catalytic
    cases = [("8.0", "1.0", "fc", "Mach 8, 30 km, 1 mm\nfully catalytic"), ("9.0", "3.0", "fc", "Mach 9, 30 km, 3 mm\nfully catalytic"),
             ("10.0", "10.0", "nc", "Mach 10, 35 km, 10 mm\nnon-catalytic")]
    fig, axs = plt.subplots(2, 3, figsize=(W2, 3.55), sharex=True, gridspec_kw=dict(height_ratios=[1, 0.9]))
    for c_i, (M, R, cat, title) in enumerate(cases):
        ax, bx = axs[0, c_i], axs[1, c_i]
        ls = "-" if cat == "fc" else "--"  # fully catalytic solid, non-catalytic dashed, as in every figure
        rows = [r for r in H if r["mach"] == M and r["R_n_mm"] == R and r["catalysis"] == cat]
        t = np.array([fnum(r["t"]) for r in rows])
        for key, col, lab in (("glass_um", C[0], "outer glass"), ("ZrO2_um", C[1], "ZrO$_2$ layer"), ("dep_um", C[2], "SiC-depleted layer")):
            ax.plot(t, [fnum(r[key]) for r in rows], color=col, ls=ls, lw=1.3, label=lab)
        ev = [e for e in E if e["mach"] == M and e["R_n_mm"] == R and e["catalysis"] == cat][0]
        marks = []
        for k, name, col in (("t_limit", "limit (c)", INK2), ("t_glass", "glass loss (a)", C[0]), ("t_active", "SiC active (b)", C[4])):
            v = fnum(ev[k])
            if v == v:
                for a_ in (ax, bx):
                    a_.axvline(v, color=col, lw=0.7, ls=":")
                marks.append((v, name, col))
        marks.sort()
        for n_, (v, name, col) in enumerate(marks):  # criterion names above the axes, staggered so that close times stay legible
            close = n_ > 0 and v - marks[n_ - 1][0] < 600.0
            ha = "left" if v < 500 else ("right" if v > 1300 else "center")
            ax.text(v, 1.03 + (0.12 if close else 0.0), name, transform=ax.get_xaxis_transform(), fontsize=7, color=col,
                    ha=ha, va="bottom")
        ax.set_title(title, loc="left", fontsize=7.5, pad=24)
        ax.set_yscale("symlog", linthresh=1.0)
        ax.set_ylim(0, 2000)
        grid(ax)
        # bottom: glass thickness vs critical thickness of criterion (a); after the glass is lost it has no thickness
        # to draw on a log axis, so the line ends there (no floor value)
        g = np.array([fnum(r["glass_um"]) for r in rows])
        g = np.where(g > 0, g, np.nan)
        dc = np.array([fnum(r.get("dcrit_nm")) for r in rows]) * 1e-3
        bx.plot(t, g, color=C[0], lw=1.3, ls=ls)
        bx.plot(t, dc, color=INK, lw=1.0, ls=ls)
        for th in (0.1, 0.5, 1.0):
            bx.axhline(th, color=INK2, lw=0.6, ls=":")
        bx.set_yscale("log")
        bx.set_ylim(1e-5, 3e3)
        bx.set_xlim(0, 1800)
        bx.set_xticks([0, 500, 1000, 1500])
        bx.set_xlabel("Time (s)")
        grid(bx)
    axs[0, 0].set_ylabel("Thickness (µm)")
    axs[1, 0].set_ylabel("Thickness (µm)")
    for a_ in list(axs[0, 1:]) + list(axs[1, 1:]):
        a_.tick_params(labelleft=False)
    fig.tight_layout(rect=(0, 0.13, 1, 1), h_pad=0.4, w_pad=0.8)
    h = [plt.Line2D([], [], color=C[0], lw=1.3), plt.Line2D([], [], color=C[1], lw=1.3), plt.Line2D([], [], color=C[2], lw=1.3),
         plt.Line2D([], [], color=INK, lw=1.0), plt.Line2D([], [], color=INK2, lw=0.6, ls=":")]
    fig.legend(h, ["outer glass", "ZrO$_2$ layer", "SiC-depleted layer", r"critical thickness $\delta_\mathrm{c}$", "fixed thresholds 0.1, 0.5, 1 µm"],
               frameon=False, loc="lower center", ncol=3, columnspacing=1.2, handlelength=1.8,
               title="solid: fully catalytic, dashed: non-catalytic", title_fontsize=7)
    save(fig, "fig5_layers")


# ---------------------------------------------------------------------------------------------------- Fig. 6 and S2
def envelope_grid(rows, R, cat, f=0.2):
    sel = [r for r in rows if abs(r["R_n"] - R) < 1e-9 and r["catalysis"] == cat and abs(r["f_SiC"] - f) < 1e-9]
    Ms = sorted({r["mach"] for r in sel})
    zs = sorted({r["alt"] for r in sel})
    tf = np.full((len(zs), len(Ms)), np.nan)
    mode = np.full((len(zs), len(Ms)), "", dtype=object)
    for r in sel:
        i, j = zs.index(r["alt"]), Ms.index(r["mach"])
        t = r["t_fail"]
        tf[i, j] = t if t == t else 1e9
        mode[i, j] = r["mode"]
    return np.array(Ms), np.array(zs) / 1e3, tf, mode


def draw_map(ax, Ms, zs, tf, mode, show_path=True, hatch=True, seq=None):
    """Binned map of the time to first failure. Cells where glass loss (a) comes first carry one diagonal stroke each
    (hatching), dark on light cells and white on dark cells, so that it stays visible at print size; SiC active
    oxidation first (b) would carry the opposite diagonal. hatch=False leaves the failure mode to its own figure
    (journal revision: Figure 6 shows the time, Figure 7 the mode)."""
    cmap = ListedColormap(seq or SEQ)
    norm = BoundaryNorm(BINS, cmap.N)
    dM, dz = Ms[1] - Ms[0], zs[1] - zs[0]
    xe = np.r_[Ms - dM / 2, Ms[-1] + dM / 2]
    ye = np.r_[zs - dz / 2, zs[-1] + dz / 2]
    pc = ax.pcolormesh(xe, ye, tf, cmap=cmap, norm=norm, edgecolors="face", linewidth=0.2)  # no white cell grid: it hid the hatching
    segs, cols = {"glass": [], "active": []}, {"glass": [], "active": []}
    for i in range(len(zs) if hatch else 0):
        for j in range(len(Ms)):
            m = mode[i, j]
            if m in segs:
                x0, x1 = xe[j] + 0.18 * dM, xe[j + 1] - 0.18 * dM
                y0, y1 = ye[i] + 0.15 * dz, ye[i + 1] - 0.15 * dz
                segs[m].append([(x0, y0), (x1, y1)] if m == "glass" else [(x0, y1), (x1, y0)])
                k = int(np.searchsorted(BINS, tf[i, j], side="right")) - 1
                cols[m].append(INK if k <= 3 else "white")
    for m in segs:
        if segs[m]:
            ax.add_collection(LineCollection(segs[m], colors=cols[m], linewidths=0.9, capstyle="round", zorder=2))
    if show_path:
        tt, za, MM = x43a_path()
        ax.plot(MM, za, color=C[1], lw=1.2)
    ax.set_xlim(xe[0], xe[-1])
    ax.set_ylim(ye[0], ye[-1])
    return pc


def mode_handles(rows):
    n_act = sum(1 for r in rows if r.get("mode") == "active")
    h = [Patch(facecolor=SEQ[4], edgecolor="white", lw=0, hatch="///", label="glass loss first (a)")]
    if n_act:  # never in the baseline maps: said in the caption instead of a legend entry for an absent pattern
        h.append(Patch(facecolor=SEQ[4], edgecolor="white", lw=0, hatch="\\\\\\", label="SiC active oxidation first (b)"))
    h.append(Patch(facecolor=SEQ[4], edgecolor="none", label="temperature limit first (c)"))
    return h


def map_figure(width):
    """2 x 3 maps with the colour bar in its own axes, clear of the panel titles (K8b, R8-02)."""
    fig, axs = plt.subplots(2, 3, figsize=(width, 3.45), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.1, right=0.835, top=0.93, bottom=0.13, wspace=0.13, hspace=0.26)
    cax = fig.add_axes([0.875, 0.30, 0.02, 0.52])
    return fig, axs, cax


def map_colorbar(fig, pc, cax):
    cb = fig.colorbar(pc, cax=cax, orientation="vertical", ticks=BIN_TICKS)
    cb.ax.set_yticklabels(BIN_LABELS, fontsize=7)
    cb.set_label("Time to first failure", fontsize=7)


def map_titles(ax, R, cat, r_i, c_i):
    ax.set_title(rf"{R * 1e3:g} mm, {'fully catalytic' if cat == 'fc' else 'non-catalytic'}", loc="left", fontsize=7.5)
    if r_i == 1:
        ax.set_xlabel("Mach number")
    if c_i == 0:
        ax.set_ylabel("Altitude (km)")


def figA1_envelope_1d():
    """Appendix: the 1D stagnation-line maps (conservative bound, Part 1.4)."""
    rows = load_jsonl(ROOT / "runs" / "envelope" / "set_A.jsonl")
    fig, axs, cax = map_figure(4.8)
    for r_i, cat in enumerate(("fc", "nc")):
        for c_i, R in enumerate((0.001, 0.003, 0.010)):
            ax = axs[r_i, c_i]
            pc = draw_map(ax, *envelope_grid(rows, R, cat))
            map_titles(ax, R, cat, r_i, c_i)
    map_colorbar(fig, pc, cax)
    handles = mode_handles(rows) + [plt.Line2D([], [], color=C[1], lw=1.2, label="X-43A path")]
    # below the axis labels, in its own row: the tight bounding box extends the page to hold it
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, fontsize=7, bbox_to_anchor=(0.47, 0.02), columnspacing=1.0,
               handlelength=1.5)
    save(fig, "figS2_envelope_1d")


def boundary_mach(cols, key):
    """Mach number at which the 30-minute boundary is drawn on a map: the edge of the last passing cell, midway between
    the highest passing and the lowest failing point (runs/analysis/analyze_2d.boundary); 4.75 / 12.25 when the column
    fails at Mach 5 / passes at 12 (outside the map). Tables 4 and 5 list the highest passing Mach number instead."""
    import analyze_2d as A2
    b = A2.boundary(cols.get(key, []))
    if b is None:
        return np.nan
    if b["M_s"] is None:
        return 4.75
    if b["M_f"] is None:
        return 12.25
    return 0.5 * (b["M_s"] + b["M_f"])


def step_line(ax, Mb, zs, **kw):
    x = np.repeat(Mb, 2)
    y = np.ravel([[z - 0.5, z + 0.5] for z in zs])
    ax.plot(x, y, path_effects=[pe.withStroke(linewidth=kw.pop("halo", 2.6), foreground=kw.pop("halo_fg", "white"))],
            **kw)


def band_bar(ax, lo, hi, z):
    """Uncertainty band of the 2D boundary at one altitude: black bar with a white rim and end caps, not clipped at
    20 and 40 km. An end beyond the map (a corner that passes at Mach 12 or fails at Mach 5) is open: a triangle
    pointing outwards instead of a cap, so the bar is not read as ending at the map edge."""
    halo = [pe.withStroke(linewidth=2.6, foreground="white")]
    kw = dict(color=INK, path_effects=halo, clip_on=False, zorder=5)
    ax.plot([lo, hi], [z, z], lw=1.1, solid_capstyle="butt", **kw)
    for x, open_, mk in ((lo, lo <= 4.75, "<"), (hi, hi >= 12.25, ">")):
        ax.plot([x], [z], marker=mk if open_ else "|", ms=3.6 if open_ else 4.8, mew=1.1, ls="none", **kw)


B2D = dict(color=INK, lw=2.0, halo=3.3)  # 2D boundary: thick black line, white rim
B1D = dict(color=INK, lw=0.8, halo=1.7)  # 1D boundary: thin black line, white rim (same colour: thick/thin is the 2D/1D code)


def fig7_envelope():
    """Headline (revision, Parts 1.2-1.4, 3.3): 2D nose-tip maps of the time to first failure (grid), with the 2D
    30-minute boundary (grid + refinement), the 1D boundary (conservative bound), the uncertainty band of the boundary
    (fully catalytic; worst and best corners, D-2D-10) and the line beyond which the surface exceeded 2073 K, the end
    of the high-temperature growth-law fit (extrapolated kinetics, D-OX-11)."""
    sys.path.insert(0, str(ROOT / "runs" / "analysis"))
    import analyze_2d as A2
    E2 = ROOT / "runs" / "envelope_2d"
    grid2 = [r for r in A2.load(E2 / "grid.jsonl") if not r.get("error")]
    ref = [r for r in A2.load(E2 / "refine.jsonl") if not r.get("error")]
    cols2 = A2.columns(grid2 + ref)
    cols1 = A2.columns([r for r in A2.load(ROOT / "runs" / "envelope" / "set_A.jsonl") if abs(r["f_SiC"] - 0.2) < 1e-9])
    corner = {n: A2.columns(A2.load(E2 / f"band_{n}.jsonl")) for n in ("worst", "best")}
    # highest surviving points whose reaction front (where k_p is evaluated) exceeded 2073 K (D-OX-13)
    front_hot = {}
    if (E2 / "front.jsonl").exists():
        for r in A2.load(E2 / "front.jsonl"):
            if r["case_id"].startswith("front_s") and not r.get("error") and r["T_f_max"] > A2.T_FIT_END:
                front_hot.setdefault((r["R_n"], r["catalysis"]), []).append((r["alt"] / 1e3, r["mach"]))
    # Journal revision (Phase 8): the former six-layer map is split. Figure 6 shows the time to first failure with the 2D and 1D
    # boundaries, the band of the 2D boundary and the X-43A path. Figure 7 (fig7_failure_mode) shows the failure mode and the
    # regions above 2073 K. Same saved outputs, no new runs.
    fig, axs, cax = map_figure(4.85)
    for r_i, cat in enumerate(("fc", "nc")):
        for c_i, R in enumerate((0.001, 0.003, 0.010)):
            ax = axs[r_i, c_i]
            Ms, zs, tf, mode = envelope_grid(grid2, R, cat)
            if len(Ms) < 2 or len(zs) < 2:  # map not run yet (preview while the 2D stages are running)
                continue
            pc = draw_map(ax, Ms, zs, tf, mode, hatch=False, seq=SEQ[:-1] + ["#f1f1f0"])  # passing cells neutral (referee)
            # both boundaries black with a white rim: thick = 2D, thin = 1D (dashes are reserved for the non-catalytic bound)
            step_line(ax, [boundary_mach(cols1, (R, cat, z * 1e3)) for z in zs], zs, **B1D)
            step_line(ax, [boundary_mach(cols2, (R, cat, z * 1e3)) for z in zs], zs, **B2D)
            if cat == "fc" and any(corner[n] for n in corner):
                for z in (20, 25, 30, 35, 40):
                    lo = boundary_mach(corner["worst"], (R, cat, z * 1e3))
                    hi = boundary_mach(corner["best"], (R, cat, z * 1e3))
                    if lo == lo and hi == hi:
                        band_bar(ax, lo, hi, z)
            map_titles(ax, R, cat, r_i, c_i)
    map_colorbar(fig, pc, cax)
    handles = [
        plt.Line2D([], [], color=INK, lw=B2D["lw"], label="30-minute boundary, 2D (thick)",
                   path_effects=[pe.withStroke(linewidth=B2D["halo"], foreground="white")]),
        plt.Line2D([], [], color=INK, lw=B1D["lw"], label="30-minute boundary, 1D (thin)",
                   path_effects=[pe.withStroke(linewidth=B1D["halo"], foreground="white")]),
        plt.Line2D([], [], color=INK, lw=1.1, marker="|", ms=5, mec=INK, label="heating and emissivity band",
                   path_effects=[pe.withStroke(linewidth=2.6, foreground="white")]),
        plt.Line2D([], [], color=C[1], lw=1.2, label="X-43A path")]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, fontsize=7, bbox_to_anchor=(0.47, 0.02), columnspacing=2.0,
               handlelength=1.6)
    save(fig, "fig6_envelope")
    fig7_failure_mode(grid2, cols2, front_hot)


MODE_COLOURS = {"pass": "#f1f1f0", "glass": SEQ[5], "limit": C[3], "active": C[5]}  # light, dark, mid: distinct in greyscale


def fig7_failure_mode(grid2, cols2, front_hot):
    """Journal revision, Figure 7: the first failure mode of the 2D maps (glass loss, temperature limit, no failure within
    30 minutes), with the 2D 30-minute boundary, the line beyond which the surface exceeds 2073 K (end of the growth-law
    fit) and the highest passing conditions whose reaction front exceeds 2073 K. From the saved grid and front runs."""
    sys.path.insert(0, str(ROOT / "runs" / "analysis"))
    import analyze_2d as A2
    fig, axs = plt.subplots(2, 3, figsize=(4.85, 3.2), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.93, bottom=0.14, wspace=0.13, hspace=0.26)
    keys = ["pass", "glass", "limit", "active"]
    cmap = ListedColormap([MODE_COLOURS[k] for k in keys])
    seen = set()
    for r_i, cat in enumerate(("fc", "nc")):
        for c_i, R in enumerate((0.001, 0.003, 0.010)):
            ax = axs[r_i, c_i]
            Ms, zs, tf, mode = envelope_grid(grid2, R, cat)
            if len(Ms) < 2 or len(zs) < 2:
                continue
            code = np.zeros(tf.shape)
            for i in range(len(zs)):
                for j in range(len(Ms)):
                    m = "pass" if tf[i, j] >= 1800 else {"glass": "glass", "active": "active"}.get(mode[i, j], "limit")
                    code[i, j] = keys.index(m)
                    seen.add(m)
            dM, dz = Ms[1] - Ms[0], zs[1] - zs[0]
            xe = np.r_[Ms - dM / 2, Ms[-1] + dM / 2]
            ye = np.r_[zs - dz / 2, zs[-1] + dz / 2]
            ax.pcolormesh(xe, ye, code, cmap=cmap, vmin=-0.5, vmax=len(keys) - 0.5, edgecolors="face", linewidth=0.2)
            Tm = np.full(tf.shape, np.nan)
            for r in grid2:
                if abs(r["R_n"] - R) < 1e-9 and r["catalysis"] == cat:
                    Tm[list(zs).index(r["alt"] / 1e3), list(Ms).index(r["mach"])] = r["T_s_max"]
            if np.nanmax(Tm) > A2.T_FIT_END > np.nanmin(Tm):
                cs = ax.contour(Ms, zs, Tm, levels=[A2.T_FIT_END], colors=[INK], linewidths=1.1, linestyles=":")
                for coll in getattr(cs, "collections", [cs]):
                    coll.set_path_effects([pe.withStroke(linewidth=2.4, foreground="white")])
            step_line(ax, [boundary_mach(cols2, (R, cat, z * 1e3)) for z in zs], zs, **B2D)
            for z, Mz in front_hot.get((R, cat), []):
                ax.plot([Mz], [z], marker="o", ms=3.4, mfc="white", mec=INK, mew=0.8, ls="none", zorder=6)
            ax.set_xlim(xe[0], xe[-1])
            ax.set_ylim(ye[0], ye[-1])
            map_titles(ax, R, cat, r_i, c_i)
    labels = {"pass": "no failure within 30 minutes", "glass": "glass loss first (a)", "limit": "temperature limit first (c)",
              "active": "SiC active oxidation first (b)"}
    handles = [Patch(facecolor=MODE_COLOURS[k], edgecolor=INK2, lw=0.4, label=labels[k]) for k in keys if k in seen]
    handles += [plt.Line2D([], [], color=INK, lw=B2D["lw"], label="30-minute boundary, 2D",
                           path_effects=[pe.withStroke(linewidth=B2D["halo"], foreground="white")]),
                plt.Line2D([], [], color=INK, lw=1.1, ls=":", label="surface above 2073 K to the right",
                           path_effects=[pe.withStroke(linewidth=2.4, foreground="white")])]
    if any(front_hot.values()):
        handles.append(plt.Line2D([], [], marker="o", ms=3.4, mfc="white", mec=INK, mew=0.8, ls="none",
                                  label="reaction front above 2073 K"))
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, fontsize=7, bbox_to_anchor=(0.54, 0.02), columnspacing=2.0,
               handlelength=1.6)
    save(fig, "fig7_failure_mode")


# ---------------------------------------------------------------------------------------------------- Fig. 7
CRITERION = ("dfix1.0", "dfix0.5", "dfix0.1", "paBalat")  # changes of the failure criterion, not input ranges (R4-23)


def fig8_sensitivity():
    """Q4: one-at-a-time effects for B1 and B2 as the ratio of the time to failure to that of the base case (log axis);
    the rank correlations are in Supplementary Table S5. The thirteen largest effects of each case are drawn: input
    ranges first, changes of the failure criterion below them in grey. Cases with no failure within the 1800 s cap are
    drawn at the cap with an open-ended triangle (the effect is at least that large)."""
    base = json.loads((ROOT / "runs" / "sensitivity" / "base_cases.json").read_text(encoding="utf-8"))
    names = list(base)
    fig, axs = plt.subplots(1, len(names), figsize=(W2, 3.0))
    sys.path.insert(0, str(ROOT / "runs" / "sensitivity"))
    import run_sensitivity as rs
    labels = {f"{p[1]}": n.replace(" (x)", "").replace(" (K)", "").replace(" (W/m/K)", "").replace(" (mm)", "").replace(" (decades)", "")
              for n, p in rs.PARAMS.items()}
    # the terms of Table 1 and the text (rule 26)
    labels.update({"k_dep": "conductivity, depleted layer", "dfix0.1": "fixed threshold 0.1 µm", "dfix0.5": "fixed threshold 0.5 µm",
                   "dfix1.0": "fixed threshold 1 µm", "perm": "O$_2$ permeability", "kpE": r"$k_\mathrm{p}$ above 2073 K, $E$",
                   "kphi": r"$k_\mathrm{p}$ above 2073 K, level", "incub": "incubation", "paBalat": "Balat-Pichelin line",
                   "kp_mult": r"ZrO$_2$-layer $k_\mathrm{p}$", "kd_mult": r"SiC depletion rate $k_\mathrm{d}$",
                   "heat_mult": "heating constant $K$", "mtf": "mass-transfer conductance", "eps_glass": "emissivity, glass",
                   "eps_oxide": "emissivity, oxide", "alpha_evap": "evaporation coefficient", "n_p": r"O$_2$-pressure exponent $n_\mathrm{p}$",
                   "f_SiO_lost": "SiO$_2$ of depleted layer lost", "T_limit": "temperature limit", "f_SiC": "SiC fraction",
                   "k_mult": "substrate conductivity", "R_n": "nose radius", "catalysis": "catalysis bound", "k_glass": "glass conductivity",
                   "k_ZrO2": "ZrO$_2$-layer conductivity", "cp_mult": "specific heat capacity", "T_dep": "depletion onset",
                   "pa_offset": "passive-to-active offset"})
    mode_name = {"glass": "glass loss", "temperature": "temperature limit", "active": "SiC active oxidation"}
    fx = ROOT / "runs" / "sensitivity" / "oat_extra.jsonl"
    extra = load_jsonl(fx) if fx.exists() else []
    CRIT_COL = "#9b9a95"
    for ax, bn in zip(axs, names):
        rows = load_jsonl(ROOT / "runs" / "sensitivity" / f"oat_{bn}.jsonl")
        rows += [r for r in extra if r["case_id"].startswith(f"{bn}_") and not r["case_id"].endswith("atomicO")]
        b = [r for r in rows if r["case_id"] == f"{bn}_base"][0]
        tb = b["t_fail"] if b["t_fail"] == b["t_fail"] else 1800.0
        eff, cens = {}, set()
        for r in rows:
            if r["case_id"] == f"{bn}_base":
                continue
            key, tag = r["case_id"][len(bn) + 1:].rsplit("_", 1)
            t = r["t_fail"] if r["t_fail"] == r["t_fail"] else 1800.0
            if r["t_fail"] != r["t_fail"]:
                cens.add((key, tag))
            eff.setdefault(key, {})[tag] = t / tb
        mag = lambda k: max(abs(math.log10(v)) for v in eff[k].values())  # noqa: E731
        order = sorted(eff, key=mag)[-13:]  # the thirteen largest effects, as before
        inputs = sorted((k for k in order if k not in CRITERION), key=mag, reverse=True)
        crit = sorted((k for k in order if k in CRITERION), key=mag, reverse=True)
        ypos, y = {}, 0.0
        for k in inputs:
            ypos[k] = y
            y += 1.0
        y_sep = y - 0.5 + 0.35
        if crit:
            y += 1.3  # gap that holds the block heading
        for k in crit:
            ypos[k] = y
            y += 1.0
        ratios = []
        for k, y0 in ypos.items():
            if k in CRITERION:
                bars = [("high", CRIT_COL, y0)]
            else:
                bars = [("low", C[0], y0 - 0.17), ("high", C[1], y0 + 0.17)]
            for tag, col, yb in bars:
                if tag not in eff[k]:
                    continue
                rr = eff[k][tag]
                ratios.append(rr)
                ax.barh(yb, rr - 1.0, left=1.0, height=0.32 if k not in CRITERION else 0.42, color=col, align="center")
                if (k, tag) in cens:  # no failure within the cap: the bar is a lower bound
                    ax.plot([rr], [yb], ">", ms=4.6, color=col, mec=INK, mew=0.5, clip_on=False, zorder=4)
        if crit:
            ax.axhline(y_sep, color=INK2, lw=0.5)
            ax.text(0.99, y_sep + 0.6, "changes of criterion", transform=ax.get_yaxis_transform(), ha="right", va="center",
                    fontsize=7, color=INK2, style="italic", zorder=6,
                    bbox=dict(facecolor="white", edgecolor="none", pad=0.6))  # not crossed by the ratio-1 line (R8p2-04)
        ax.set_xscale("log")
        lo, hi = min(ratios + [1.0]) / 1.3, max(ratios + [1.0]) * 1.3
        ax.set_xlim(lo, hi)
        cand = (0.03, 0.1, 0.3, 1, 3) if hi / lo > 20 else (0.25, 0.5, 1, 2)  # spaced so the labels never touch
        ax.xaxis.set_major_locator(FixedLocator([v for v in cand if lo <= v <= hi]))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.set_yticks(list(ypos.values()))
        ax.set_yticklabels([labels.get(k, k) for k in ypos], fontsize=7)
        ax.set_ylim(y - 0.4, -0.6)
        ax.axvline(1.0, color=INK, lw=0.7)
        ax.set_xlabel("Time to failure,\nratio to the base case")
        tshow = int(tb) if tb < 200 else int(tb // 10 * 10)  # rounded down as in the text: at least 200 s to tens
        ax.set_title(f"{bn}: {tshow} s ({mode_name.get(b['mode'], b['mode'])})", loc="right", fontsize=7.5)  # over the labels, inside the width
        grid(ax)
    fig.tight_layout(rect=(0, 0.1, 1, 1), w_pad=0.6)
    fig.legend([Patch(color=C[0]), Patch(color=C[1]), Patch(color=CRIT_COL),
                plt.Line2D([], [], marker=">", ms=4.6, color=INK2, mec=INK, mew=0.5, ls="none")],
               ["low value of the input", "high value of the input", "change of the failure criterion",
                "no failure within 1800 s: effect at least as shown"], frameon=False, loc="lower center", ncol=2, columnspacing=1.5)
    save(fig, "fig7_sensitivity")


# ---------------------------------------------------------------------------------------------------- Fig. 8
def fig9_q5():
    """Q5: lateral conduction. (a) stagnation temperature, 1D vs 2D, at Mach 9 / 30 km with the baseline and the sharp-cone
    flank distributions; (b) peak stagnation temperature 2D vs 1D for all constant-condition cases (600 s); (c) gain of
    the 30-minute boundary from 1D to 2D over altitude, from the full maps (runs/analysis/boundaries_2d.csv): the
    quantity of Tables 4 and 5, highest passing 2D minus highest passing 1D Mach number. Columns where the 2D map has
    no failure up to Mach 12 give only a lower bound, drawn as an open triangle pointing right."""
    R = read_csv(ROOT / "runs" / "q5_2d" / "q5_results.csv")
    H = read_csv(ROOT / "runs" / "q5_2d" / "q5_histories.csv")
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(W2, 2.5), gridspec_kw=dict(width_ratios=[1.1, 1.0, 0.95]))
    for i, Rn in enumerate(("1.0", "3.0")):
        for model, dist, lw, mk, lab in (("1d", "lees", 0.7, None, "1D"), ("2d", "lees", 1.7, None, "2D"),
                                         ("2d", "tauber", 1.0, "o", "2D, sharp-cone correlation")):
            rows = [r for r in H if r["kind"] == "const" and r["mach"] == "9.0" and r["R_n_mm"] == Rn
                    and r["theta_deg"] == "5.0" and r["model"] == model and r["heat_dist"] == dist]
            a.plot([fnum(r["t"]) for r in rows], [fnum(r["T_s"]) for r in rows], color=C[i], lw=lw, marker=mk, ms=2, markevery=0.1,
                   label=f"{float(Rn):g} mm, {lab}")
    a.axhline(T_LIMIT, color=INK2, lw=0.8, ls=":")
    a.text(1.1, T_LIMIT + 40, "limit (c)", fontsize=7, color=INK2, ha="left", va="bottom")
    a.set_xscale("log")
    a.set_xlim(1, 600)
    a.set_ylim(500, 2800)
    a.set_xlabel("Time (s)")
    a.set_ylabel(TS_LABEL)
    a.set_title("(a) Mach 9, 30 km", loc="left")
    grid(a)
    for dist, th, mk, lab in (("lees", "5.0", "o", "baseline, 5° cone"), ("tauber", "5.0", "s", "sharp-cone correlation, 5°"),
                              ("lees", "10.0", "^", "baseline, 10° cone")):
        for i, Rn in enumerate(("1.0", "3.0")):
            rows = [r for r in R if r["kind"] == "const" and r["heat_dist"] == dist and r["theta_deg"] == th
                    and r["R_n_mm"] == Rn]
            if rows:
                b.plot([fnum(r["1d_T_max"]) for r in rows], [fnum(r["2d_T_max"]) for r in rows], mk, color=C[i], ms=4.5,
                       mfc=C[i] if dist == "lees" and th == "5.0" else "none",
                       label=f"{lab}" if i == 0 else None)
    lim = [1400, 3000]
    b.plot(lim, lim, color=INK, lw=0.7)
    b.text(2560, 2700, "1:1", fontsize=7, color=INK, ha="right", va="bottom")
    b.set_xlim(lim)
    b.set_ylim(lim)
    b.set_xlabel(r"Peak $T_\mathrm{s}$, 1D (K)")
    b.set_ylabel(r"Peak $T_\mathrm{s}$, 2D (K)")
    b.set_title("(b) peaks over 600 s", loc="left")
    grid(b)
    B2 = read_csv(ROOT / "runs" / "analysis" / "boundaries_2d.csv")
    for i, Rn in enumerate(("1.0", "3.0", "10.0")):
        for cat, ls in (("fc", "-"), ("nc", "--")):
            rows = sorted((r for r in B2 if r["R_n_mm"] == Rn and r["catalysis"] == cat), key=lambda r: fnum(r["alt_km"]))
            z = np.array([fnum(r["alt_km"]) for r in rows])
            gain = np.array([fnum(r["M_s_2d"]) - fnum(r["M_s_1d"]) for r in rows])  # as in Tables 4 and 5
            open_ = np.array([fnum(r["M_f_2d"]) != fnum(r["M_f_2d"]) for r in rows])  # 2D passes at Mach 12: lower bound only
            # step profile: each altitude column is a vertical segment at its value (z +- 0.5 km), like the map boundaries
            c.plot(np.repeat(np.where(open_, np.nan, gain), 2), np.ravel([[zz - 0.5, zz + 0.5] for zz in z]), color=C[i], ls=ls, lw=1.2)
            if open_.any():
                c.plot(gain[open_], z[open_], ">", ms=4.0, mfc="white", mec=C[i], mew=0.9, ls="none", zorder=4, clip_on=False)
    c.axvline(0, color=INK, lw=0.7)
    c.set_xlabel("2D minus 1D\nboundary (Mach)")
    c.set_ylabel("Altitude (km)")
    c.set_ylim(19.5, 40.5)
    c.set_xlim(-0.25, 4.9)
    c.set_title("(c) 1D to 2D shift", loc="left")
    grid(c)
    fig.tight_layout(rect=(0, 0.2, 1, 1), w_pad=0.4)
    L = plt.Line2D
    h = [L([], [], color=C[0], lw=1.2), L([], [], color=C[1], lw=1.2), L([], [], color=C[2], lw=1.2), L([], [], color=INK2, lw=1.2),
         L([], [], color=INK2, lw=0.7), L([], [], color=INK2, lw=1.7), L([], [], color=INK2, lw=1.0, marker="o", ms=2),
         L([], [], color=INK2, lw=1.2, ls="--"),
         L([], [], color=INK2, marker="o", ls="none", ms=4), L([], [], color=INK2, marker="s", mfc="none", ls="none", ms=4),
         L([], [], color=INK2, marker="^", mfc="none", ls="none", ms=4), L([], [], color=INK2, marker=">", mfc="white", mew=0.9, ls="none", ms=4)]
    lab = [r"$R_\mathrm{n}$ = 1 mm", r"$R_\mathrm{n}$ = 3 mm", r"$R_\mathrm{n}$ = 10 mm, (c) only", "fully catalytic",
           "1D (thin)", "2D (thick)", "2D, sharp-cone", "non-catalytic, (c) only",
           "baseline, 5°", "sharp-cone, 5°", "baseline, 10°", "lower bound: 2D passes Mach 12"]
    fig.legend(h, lab, frameon=False, loc="lower center", ncol=3, columnspacing=1.4, handlelength=1.8, labelspacing=0.3)
    save(fig, "fig8_q5_2d")


# ---------------------------------------------------------------------------------------------------- Fig. S3
def fig10_variants():
    """(a) SiC fraction: survival boundary (t_fail >= 1800 s) and 60 s contour, R_n = 3 mm; (b) temperature jump after glass loss."""
    A = load_jsonl(ROOT / "runs" / "envelope" / "set_A.jsonl")
    B = load_jsonl(ROOT / "runs" / "envelope" / "set_B.jsonl")
    Cs = load_jsonl(ROOT / "runs" / "envelope" / "set_C.jsonl")
    fig, (a, b) = plt.subplots(1, 2, figsize=(5.2, 2.8))

    def contour(ax, rows, R, cat, f, col, ls, lab):
        Ms, zs, tf, _ = envelope_grid(rows, R, cat, f)
        lt = np.log10(np.minimum(tf, 1e4))
        # one line width for both (thick/thin means 2D/1D elsewhere): the 60 s contour is a light tint of the colour
        ax.contour(Ms, zs, lt, levels=[math.log10(60.0), math.log10(1799.0)], colors=[tint(col), col], linestyles=[ls, ls],
                   linewidths=[1.4, 1.4])
        ax.plot([], [], color=col, ls=ls, lw=1.4, label=lab)
    for f, col, ls in ((0.1, C[0], "-"), (0.2, C[1], "-"), (0.3, C[2], "-")):
        contour(a, A + B, 0.003, "fc", f, col, ls, f"{int(f * 100)} vol%")
    a.plot([], [], color=INK2, lw=1.4, label="dark: 30-minute boundary")
    a.plot([], [], color=tint(INK2), lw=1.4, label="light: 60 s contour")
    a.set_title(r"(a) SiC fraction, 3 mm, fully catalytic", loc="left", fontsize=7.5)
    a.set_xlabel("Mach number")
    a.set_xlim(7.5, 10.5)
    a.set_ylabel("Altitude (km)")
    a.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.45, -0.25), ncol=2, columnspacing=1.0, handlelength=1.5,
             labelspacing=0.25)  # below: clear of the contours; colour column, then shade column
    grid(a)
    J = read_csv(ROOT / "runs" / "histories" / "jump_histories.csv")
    b.remove()
    b = fig.add_subplot(1, 2, 2)
    for cat, col, ls, lab in (("fc", C[1], "-", "fully catalytic"), ("nc", C[0], "--", "non-catalytic"),
                              ("switch", C[3], "-.", "switch at glass loss")):
        rows = [r for r in J if r["mach"] == "10.5" and r["alt_km"] == "38.0" and r["catalysis"] == cat]
        b.plot([fnum(r["t"]) for r in rows], [fnum(r["T_s"]) for r in rows], color=col, ls=ls, lw=1.5 if cat == "switch" else 1.3, label=lab)
    sw = [r for r in J if r["mach"] == "10.5" and r["alt_km"] == "38.0" and r["catalysis"] == "switch"]
    tf_ = fnum(sw[0]["t_fail"])  # first failure of the switched wall
    tg = min(fnum(r["t"]) for r in sw if fnum(r["t"]) > 30 and fnum(r["glass_um"]) <= 1e-6)  # glass gone (2 s records)
    b.axvline(tg, color=INK2, lw=0.7, ls=":")  # failure times are dotted vertical lines in every figure
    if abs(tg - tf_) < 20:
        b.text(tg + 30, 2340, f"limit (c) and glass loss (a)\nafter about {tf_ // 10 * 10:.0f} s", fontsize=7, color=INK2,
               ha="left", va="top")
    else:
        b.text(tg - 15, 2700, "glass lost (a)", fontsize=7, color=INK2, ha="right")
    b.axhline(T_LIMIT, color=INK2, lw=0.8, ls=":")
    b.text(400, T_LIMIT + 10, "limit (c)", fontsize=7, color=INK2)  # between the curves, clear of the rising lines
    b.set_xlim(0, 1800)
    b.set_ylim(2000, 2750)
    b.set_xlabel("Time (s)")
    b.set_ylabel(TS_LABEL)
    b.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.45, -0.25), ncol=1, handlelength=2.2, labelspacing=0.25)  # below: clear of the glass-loss line
    b.set_title(r"(b) Mach 10.5, 38 km, 10 mm", loc="left", fontsize=7.5)
    grid(b)
    fig.tight_layout(w_pad=1.5)
    save(fig, "figS3_variants")


# ---------------------------------------------------------------------------------------------------- Fig. S4
VARIANT_LABELS = [  # (stage, name, label) in display order
    ("band", "worst", "worst corner"), ("band", "best", "best corner"),
    ("band", "q0_hi", "stagnation heating ×1.15"), ("band", "q0_lo", "stagnation heating ×0.85"),
    # the band of the heating beyond 70 deg (Section 4.3); "flank heating x0.5 / x2" are the different runs of Table 5
    ("band", "heat_hi", "heating beyond 70° ×1.5"), ("band", "heat_lo", "heating beyond 70° ×0.66"),
    ("band", "eps_lo", "emissivity 0.70 / 0.50"), ("band", "eps_hi", "emissivity 0.95 / 0.85"),
    ("glass", "dfail_1.0", "glass loss at fixed 1 µm"), ("glass", "dfail_0.5", "glass loss at fixed 0.5 µm"),
    ("glass", "dfail_0.1", "glass loss at fixed 0.1 µm"),
    # external review (stage "external"); the first element is the colour group of the row
    ("glass", "gate_0.5", "gate of criterion (a), 0.5 µm"), ("glass", "gate_2", "gate of criterion (a), 2 µm"),
    ("kinetics", "perm_low", "O$_2$ permeability low"), ("kinetics", "perm_high", "O$_2$ permeability high"),
    ("kinetics", "kp_x0.5", r"$k_\mathrm{p}$ above 2073 K ×0.5"), ("kinetics", "kp_x3", r"$k_\mathrm{p}$ above 2073 K ×3"),
    ("kinetics", "kp_E250", r"$k_\mathrm{p}$ above 2073 K, $E$ = 250 kJ mol$^{-1}$"),
    ("kinetics", "kp_E354", r"$k_\mathrm{p}$ above 2073 K, $E$ = 354 kJ mol$^{-1}$"),
    ("kinetics", "kp_all_x0.33", r"$k_\mathrm{p}$ ×1/3 at all temperatures"),
    ("kinetics", "kp_all_x3", r"$k_\mathrm{p}$ ×3 at all temperatures"),
    ("kinetics", "incubation", "incubation of SiC depletion"),
    ("kinetics", "kd_arrh", r"$k_\mathrm{d}$ Arrhenius, extrapolated"),
    ("kinetics", "kd_arrh_cap", r"$k_\mathrm{d}$ Arrhenius, held above 2073 K"),
    ("transport", "mtf_x1.25", "mass-transfer conductance ×1.25"), ("transport", "mtf_x2", "mass-transfer conductance ×2"),
    ("transport", "alpha_0.1", "evaporation coefficient 0.1"), ("transport", "alpha_0.005", "evaporation coefficient 0.005"),
    ("transport", "kZrO2_1.0", "ZrO$_2$-layer conductivity 1.0 W m$^{-1}$ K$^{-1}$"),
    ("transport", "kZrO2_2.5", "ZrO$_2$-layer conductivity 2.5 W m$^{-1}$ K$^{-1}$"),
    ("active", "pa_balat", "active line of Balat-Pichelin et al."),
    ("active", "pa_fe", "active line refitted with study intercepts"),
    ("active", "pa_fe_off0.3", "refitted line, offset 0.3 decade"),
]


def fig11_variants_2d():
    """Shift of the 2D 30-minute boundary (fully catalytic) for each variant, over the five band altitudes
    (runs/analysis/variants_2d.csv): bar = range over 20-40 km, dot = 30 km."""
    V = read_csv(ROOT / "runs" / "analysis" / "variants_2d.csv")
    have = [(st, n, lab) for st, n, lab in VARIANT_LABELS if any(r["variant"] == n for r in V)]
    fig, axs = plt.subplots(1, 3, figsize=(W2, 0.18 * len(have) + 0.9), sharey=True, sharex=True)  # fits one Supplement page
    grp_col = {"band": C[0], "glass": C[2], "kinetics": C[1], "transport": C[3], "active": C[4]}
    for ax, Rn in zip(axs, ("1.0", "3.0", "10.0")):
        for i, (st, n, lab) in enumerate(have):
            rows = [r for r in V if r["variant"] == n and r["R_n_mm"] == Rn]
            sh = [fnum(r["shift"]) for r in rows]
            sh = [x for x in sh if x == x]
            if not sh:
                continue
            ax.plot([min(sh), max(sh)], [i, i], color=grp_col[st], lw=2.4, solid_capstyle="butt", alpha=0.55)
            z30 = [fnum(r["shift"]) for r in rows if r["alt_km"] == "30.0"]
            if z30 and z30[0] == z30[0]:
                ax.plot(z30[0], i, "o", color=grp_col[st], ms=3.6)
            open_ = [r for r in rows if r["M_f"] in ("", "None")]
            if open_:
                ax.text(max(sh) + 0.25, i, "≥", fontsize=7, ha="left", va="center", color=INK2)  # clear of the bar end and dot
        ax.axvline(0, color=INK, lw=0.7)
        ax.set_title(rf"$R_\mathrm{{n}}$ = {float(Rn):g} mm", loc="left", fontsize=7.5)
        grid(ax)
    axs[0].set_yticks(range(len(have)))
    axs[0].set_yticklabels([lab for _, _, lab in have], fontsize=7)
    from matplotlib.ticker import MultipleLocator
    for ax in axs:  # readable shifts: labels every Mach number, minor ticks at the 0.25 grid step (R8p2-08)
        ax.xaxis.set_major_locator(MultipleLocator(1.0))
        ax.xaxis.set_minor_locator(MultipleLocator(0.25))
    axs[0].invert_yaxis()
    fig.supxlabel("Shift of the 30-minute boundary (Mach number)", fontsize=7.5, y=0.02)
    fig.tight_layout()
    save(fig, "figS4_variants_2d")


FIGS = {"fig1": fig1_schematic, "fig2": fig3_furnace, "fig3": fig4_arcjet, "fig4": fig5_temperature, "fig5": fig6_layers,
        "fig6": fig7_envelope, "fig7": fig8_sensitivity, "fig8": fig9_q5, "figS1": fig2_heating_checks, "figS2": figA1_envelope_1d,
        "figS3": fig10_variants, "figS4": fig11_variants_2d}

if __name__ == "__main__":
    names = sys.argv[1:] or list(FIGS)
    for n in names:
        FIGS[n]()
