"""Tables of the paper from saved outputs.  Usage: python runs/analysis/make_tables.py table1|table2|table3

table1  the inputs varied in Q4 with their low, baseline and high values (the full table is in the paper and INPUTS.md)
table2  outcome of each validation check, read from its report in validation/
table3  gain of the 30-minute boundary from 1D to 2D for the flank-heating cases (Table 5 of the paper) -> runs/analysis/table3_gain.csv
"""
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "runs"))
sys.path.insert(0, str(ROOT / "runs" / "sensitivity"))

REPORTS = [("V1", "V1_conduction/report_part_a.md"), ("V2", "V2_numerics/report.md"), ("V3", "V3_radiative_equilibrium/report.md"),
           ("V4", "V4_heating/report.md"), ("V5", "V5_heating_distribution/report.md"), ("V6", "V6_oxidation/report.md"),
           ("V7", "V7_arcjet/report.md"), ("V8", "V8_measured_flux/report.md"), ("V9", "V9_plasmatron/report.md")]


def table1():
    import run_sensitivity as rs
    print(f"{'input':42s} {'low':>10s} {'base':>10s} {'high':>10s}")
    for name, (_, _, lo, b, hi, _) in rs.PARAMS.items():
        print(f"{name:42s} {lo!s:>10.10s} {b!s:>10.10s} {hi!s:>10.10s}")


def table2():
    for v, rel in REPORTS:
        text = (ROOT / "validation" / rel).read_text(encoding="utf-8")
        hits = [l.strip() for l in text.splitlines() if re.search(r"\*\*[^*]*(PASS|FAIL)[^*]*\*\*|Result:|Overall|inside the band", l)]
        print(f"{v}: " + (" / ".join(h[:120] for h in hits[-2:] if v == "V8" or h == hits[-1]) if hits else "see " + rel))


def table3():
    b = {(r["R_n_mm"], r["alt_km"]): r for r in csv.DictReader(open(ROOT / "runs/analysis/boundaries_2d.csv", encoding="utf-8")) if r["catalysis"] == "fc"}
    v = list(csv.DictReader(open(ROOT / "runs/analysis/variants_2d.csv", encoding="utf-8")))
    names = ("flank_x0.5", "flank_x2", "sharp_cone")
    rows = []
    for (R, z), r in sorted(b.items(), key=lambda kv: (float(kv[0][0]), float(kv[0][1]))):
        if R not in ("1.0", "3.0") or z not in ("20.0", "25.0", "30.0", "35.0", "40.0"):
            continue
        m1 = float(r["M_s_1d"])
        row = dict(R_n_mm=R, alt_km=z, baseline=float(r["M_s_2d"]) - m1)
        for n in names:
            row[n] = [float(x["M_s"]) - m1 for x in v if x["variant"] == n and x["R_n_mm"] == R and x["alt_km"] == z][0]
        rows.append(row)
    with open(ROOT / "runs/analysis/table3_gain.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r)


if __name__ == "__main__":
    {"table1": table1, "table2": table2, "table3": table3}[sys.argv[1]]()
