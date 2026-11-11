"""Revision (Parts 1.2-1.4, 2.3, 3.3): 30-minute boundaries of the 2D nose-tip maps, their uncertainty band, the shifts
from the criterion / kinetics / transition-line variants, the 1D comparison, and the cells that rely on the extrapolated
high-temperature growth law.  Every number about the 2D maps quoted in the paper comes from the files written here.

Inputs : runs/envelope_2d/{grid,refine,band_*,glass_*,active_*,kinetics_*}.jsonl, runs/envelope/set_A.jsonl (1D)
Outputs: runs/analysis/boundaries_2d.csv   one row per (R_n, bound, altitude): 2D and 1D boundary, first mode beyond it
         runs/analysis/variants_2d.csv     one row per (variant, R_n, altitude), fully catalytic bound
         runs/analysis/envelope_2d.json    summary used by the figures and trace_numbers.py
         runs/analysis/envelope_2d.md      the same as readable tables

Boundary of a column: M_s = highest Mach number that survives the 1800 s cap below the lowest failing Mach number M_f
(all points are first-failure runs); resolution M_f - M_s (0.25 where refined, 0.5 on the plain grid).  '<5' if Mach 5
fails, '>12' if nothing up to Mach 12 fails.  A column is flagged non-monotone if a point above M_f survives.
Extrapolated kinetics (D-OX-11): a point "uses the extrapolated law" if its surface exceeded 2073 K (the end of the
high-temperature fit) at any time before it stopped; the surface temperature bounds the reaction-front temperature
from above, so the marking is conservative.
"""
import csv
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
E2 = ROOT / "runs" / "envelope_2d"
OUT = ROOT / "runs" / "analysis"
sys.path.insert(0, str(ROOT / "runs" / "envelope_2d"))
import run_envelope_2d as R2  # noqa: E402

T_FIT_END = 2073.15
CAP = 1800.0


def load(p):
    p = Path(p)
    if not p.exists():
        return []
    lines = p.read_text(encoding="utf-8").splitlines()
    out = []
    for i, l in enumerate(lines):
        if not l.strip():
            continue
        try:
            out.append(json.loads(l))
        except json.JSONDecodeError:
            if i == len(lines) - 1:  # last line still being written by a running stage
                continue
            raise
    return out


def ok(r):
    return not r.get("error")


def survives(r):
    return r["mode"] == "none"


def boundary(rows):
    """(M_s, M_f, mode at M_f, t_fail at M_f, nonmonotone, extrapolated at M_s, extrapolated at M_f)."""
    rows = sorted((r for r in rows if ok(r)), key=lambda r: r["mach"])
    if not rows:
        return None
    fails = [r for r in rows if not survives(r)]
    rf = fails[0] if fails else None
    M_f = rf["mach"] if rf else math.inf
    below = [r for r in rows if survives(r) and r["mach"] < M_f]
    rs = below[-1] if below else None
    nonmono = any(survives(r) and r["mach"] > M_f for r in rows)
    return dict(M_s=rs["mach"] if rs else None, M_f=None if rf is None else M_f,
                mode_f=rf["mode"] if rf else "", t_fail_f=rf["t_fail"] if rf else None,
                T_max_s=rs["T_s_max"] if rs else None, T_max_f=rf["T_s_max"] if rf else None,
                extrap_s=bool(rs and rs["T_s_max"] > T_FIT_END), extrap_f=bool(rf and rf["T_s_max"] > T_FIT_END),
                nonmono=nonmono, n=len(rows))


def fmt(b):
    if b is None:
        return "n/a"
    if b["M_s"] is None:
        return "<5"
    if b["M_f"] is None:
        return ">12"
    return f"{b['M_s']:g}"


def columns(rows):
    out = {}
    for r in rows:
        out.setdefault((r["R_n"], r["catalysis"], r["alt"]), []).append(r)
    return out


def main():
    grid = load(E2 / "grid.jsonl")
    ref = load(E2 / "refine.jsonl")
    base_cols = columns(grid + ref)
    one_d = columns([r for r in load(ROOT / "runs" / "envelope" / "set_A.jsonl") if abs(r["f_SiC"] - 0.2) < 1e-9])
    errors = [r["case_id"] for r in grid + ref if not ok(r)]
    rows_b, summ = [], dict(n_grid=len(grid), n_refine=len(ref), errors=errors, columns={})
    for key in sorted(base_cols):
        R, cat, z = key
        b2 = boundary(base_cols[key])
        b1 = boundary(one_d.get(key, []))
        g = [r for r in grid if (r["R_n"], r["catalysis"], r["alt"]) == key and ok(r)]
        n_ext = sum(1 for r in g if r["T_s_max"] > T_FIT_END)
        rows_b.append(dict(R_n_mm=R * 1e3, catalysis=cat, alt_km=z / 1e3, M_s_2d=b2["M_s"], M_f_2d=b2["M_f"],
                           mode_2d=b2["mode_f"], t_fail_2d=b2["t_fail_f"], T_max_s_2d=b2["T_max_s"],
                           extrap_s_2d=b2["extrap_s"], nonmono_2d=b2["nonmono"],
                           M_s_1d=b1["M_s"] if b1 else None, M_f_1d=b1["M_f"] if b1 else None,
                           mode_1d=b1["mode_f"] if b1 else "", n_grid_extrap=n_ext, n_grid=len(g)))
        summ["columns"][f"R{R * 1e3:g}_{cat}_z{z / 1e3:g}"] = dict(b2d=b2, b1d=b1, label_2d=fmt(b2),
                                                                 label_1d=fmt(b1))
    # modes of all failing grid points, per map
    maps = {}
    for r in grid:
        if not ok(r):
            continue
        k = f"R{r['R_n'] * 1e3:g}_{r['catalysis']}"
        m = maps.setdefault(k, dict(n=0, survive=0, glass=0, active=0, temperature=0, extrap=0))
        m["n"] += 1
        m["survive" if survives(r) else r["mode"]] += 1
        m["extrap"] += r["T_s_max"] > T_FIT_END
    summ["maps"] = maps
    # variants (fully catalytic, band altitudes)
    rows_v = []
    stages = dict(band=R2.BAND_VARIANTS, glass=R2.GLASS_VARIANTS, active=R2.ACTIVE_VARIANTS,
                  kinetics=R2.KINETICS_VARIANTS, flank=R2.FLANK_VARIANTS, oxidant=R2.OXIDANT_VARIANTS,
                  transport=R2.TRANSPORT_VARIANTS, external=R2.EXTERNAL_VARIANTS)
    summ["variants"] = {}
    for stage, variants in stages.items():
        for name in variants:
            vrows = load(E2 / f"{stage}_{name}.jsonl")
            if not vrows:
                continue
            vc = columns(vrows)
            for R in R2.RADII:
                for z in R2.BAND_ALTS:
                    key = (R, "fc", z)
                    # the variant search only runs near the boundary; points of the variant alone define it
                    if key not in vc:  # the flank stage has no 10 mm columns
                        continue
                    bv = boundary(vc.get(key, []))
                    b0 = boundary(base_cols.get(key, []))
                    shift = None
                    if bv and b0 and bv["M_s"] is not None and b0["M_s"] is not None:
                        shift = bv["M_s"] - b0["M_s"]
                    rows_v.append(dict(stage=stage, variant=name, R_n_mm=R * 1e3, alt_km=z / 1e3,
                                       M_s_base=b0["M_s"] if b0 else None, M_s=bv["M_s"] if bv else None,
                                       M_f=bv["M_f"] if bv else None, shift=shift,
                                       mode=bv["mode_f"] if bv else "", t_fail=bv["t_fail_f"] if bv else None,
                                       nonmono=bv["nonmono"] if bv else None, n=bv["n"] if bv else 0))
                    summ["variants"][f"{name}_R{R * 1e3:g}_z{z / 1e3:g}"] = rows_v[-1]
    # band = range spanned by the corners, and by all band variants (D-2D-10)
    band = {}
    for R in R2.RADII:
        for z in R2.BAND_ALTS:
            k = f"R{R * 1e3:g}_z{z / 1e3:g}"
            vals = {r["variant"]: r["M_s"] for r in rows_v
                    if r["stage"] == "band" and r["R_n_mm"] == R * 1e3 and r["alt_km"] == z / 1e3}
            base = summ["columns"].get(f"R{R * 1e3:g}_fc_z{z / 1e3:g}", {}).get("b2d", {}) or {}
            band[k] = dict(base=base.get("M_s"), worst=vals.get("worst"), best=vals.get("best"),
                           oat={v: vals.get(v) for v in ("heat_lo", "heat_hi", "q0_lo", "q0_hi", "eps_lo", "eps_hi")})
    summ["band"] = band
    for name, data in (("boundaries_2d.csv", rows_b), ("variants_2d.csv", rows_v)):
        if data:
            with open(OUT / name, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=list(data[0]))
                w.writeheader()
                w.writerows(data)
    (OUT / "envelope_2d.json").write_text(json.dumps(summ, indent=1, default=str), encoding="utf-8")
    # readable tables
    L = ["# 2D nose-tip maps: 30-minute boundaries (generated by runs/analysis/analyze_2d.py)", "",
         f"Grid points: {len(grid)}, refinement points: {len(ref)}, errors: {len(errors)}.", ""]
    for k, m in sorted(maps.items()):
        L.append(f"- {k}: {m['survive']}/{m['n']} survive 1800 s; first failures: glass {m['glass']}, active "
                 f"{m['active']}, temperature {m['temperature']}; points whose surface exceeded 2073 K: {m['extrap']}")
    L += ["", "## Highest surviving Mach number (2D | 1D) per altitude", ""]
    hdr = ["alt (km)"] + [f"R{R * 1e3:g} {c}" for R in R2.RADII for c in ("fc", "nc")]
    L += ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    for z in R2.ALTS:
        cells = [f"{z / 1e3:g}"]
        for R in R2.RADII:
            for c in ("fc", "nc"):
                s = summ["columns"].get(f"R{R * 1e3:g}_{c}_z{z / 1e3:g}")
                if not s:
                    cells.append("")
                    continue
                b2 = s["b2d"]
                tag = {"glass": "g", "active": "a", "temperature": "T", "": ""}[b2["mode_f"]] if b2 else ""
                ext = "*" if b2 and b2["extrap_s"] else ""
                nm = "!" if b2 and b2["nonmono"] else ""
                cells.append(f"{s['label_2d']}{tag}{ext}{nm} | {s['label_1d']}".replace(" | ", " / "))
        L.append("| " + " | ".join(cells) + " |")
    L += ["", "Suffix: first failure just beyond the boundary g = glass loss, a = SiC active, T = temperature limit; "
          "* = the surviving boundary point exceeded 2073 K (extrapolated growth law); ! = non-monotone column.", ""]
    L += ["## Band and variants (fully catalytic): highest surviving Mach number", ""]
    names = [n for st in stages.values() for n in st]
    L += ["| R (mm) | alt (km) | base | " + " | ".join(names) + " |", "|" + "---|" * (3 + len(names))]
    for R in R2.RADII:
        for z in R2.BAND_ALTS:
            b0 = summ["columns"].get(f"R{R * 1e3:g}_fc_z{z / 1e3:g}", {}).get("label_2d", "")
            cells = [f"{R * 1e3:g}", f"{z / 1e3:g}", b0]
            for n in names:
                v = summ["variants"].get(f"{n}_R{R * 1e3:g}_z{z / 1e3:g}")
                if v is None:
                    cells.append("")
                else:
                    cells.append("<5" if v["M_s"] is None else (">12" if v["M_f"] is None else f"{v['M_s']:g}")
                                 + {"glass": "g", "active": "a", "temperature": "T", "": ""}[v["mode"]])
            L.append("| " + " | ".join(cells) + " |")
    (OUT / "envelope_2d.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
