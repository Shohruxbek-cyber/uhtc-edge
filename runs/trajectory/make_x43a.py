"""Build the X-43A flight 3 (Mach 10) trajectory used in the study from the curves extracted from
Karlgaard et al. (2005), "Hyper-X Mach 10 Trajectory Reconstruction", AIAA 2005-5920 (NTRS 20050232845;
literature/pdf/Karlgaard2005_HyperX_Mach10_reconstruction.pdf):
  * boost: Fig. 9(a) altitude and 9(b) Mach number, solid line = best-estimated trajectory (BET), PDF page 12;
  * descent: Fig. 11(a) and 11(b), solid line (BET), PDF page 15.
The curves were extracted from the PDF vector paths (runs/trajectory/digitize_vec2.py; axis calibration from
the tick-label positions, residual <= 0.05 s and <= 0.001 axis units) -> literature/digitized/Karlgaard2005_*.csv.
DIGITIZED values.  Between the end of Fig. 9 (88.8 s) and the start of Fig. 11 (111.8 s) the vehicle separated
(88.35 s, Table 3: 109,653 ft, Mach 9.702) and flew its engine test near Mach 9.7; this 23 s gap is filled by
linear interpolation (ASSUMPTION, DECISIONS.md D-TRAJ-1).  Altitude: ft x 0.3048 (exact) -> m.

Output: runs/trajectory/x43a_mach10.csv  (t_s, altitude_m, mach) on a 0.5 s grid, Mach >= 3 only
(the correlations are for hypersonic flow).
"""
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "literature" / "digitized"


def load(name, scale=1.0):
    a = np.loadtxt(D / f"Karlgaard2005_{name}.csv", delimiter=",")
    t, v = a[:, 0], a[:, 1] * scale
    # bin to 0.25 s and average (the vector path is dense and single-valued)
    b = np.round(t * 4) / 4
    tb = np.unique(b)
    return tb, np.array([v[b == x].mean() for x in tb])


def main():
    ta1, a1 = load("boost_alt", 1e4 * 0.3048)
    tm1, m1 = load("boost_mach")
    ta2, a2 = load("desc_alt", 1e4 * 0.3048)
    tm2, m2 = load("desc_mach")
    # separation point from Table 3 (BET column)
    t_sep, alt_sep, m_sep = 88.35, 109653 * 0.3048, 9.702
    ta = np.r_[ta1[ta1 < t_sep], t_sep, ta2]
    al = np.r_[a1[ta1 < t_sep], alt_sep, a2]
    tm = np.r_[tm1[tm1 < t_sep], t_sep, tm2]
    ma = np.r_[m1[tm1 < t_sep], m_sep, m2]
    t = np.arange(np.ceil(max(ta[0], tm[0])), np.floor(min(ta[-1], tm[-1])) + 1e-9, 0.5)
    alt = np.interp(t, ta, al)
    mach = np.interp(t, tm, ma)
    keep = mach >= 3.0
    i0, i1 = np.argmax(keep), len(keep) - np.argmax(keep[::-1]) - 1
    out = ROOT / "runs" / "trajectory" / "x43a_mach10.csv"
    np.savetxt(out, np.c_[t[i0:i1 + 1], alt[i0:i1 + 1], mach[i0:i1 + 1]], delimiter=",", fmt="%.2f,%.1f,%.4f",
               header="t_s,altitude_m,mach  (X-43A flight 3; Karlgaard et al. 2005 Figs. 9 and 11, DIGITIZED; see make_x43a.py)")
    print(out, i1 - i0 + 1, "rows; t", t[i0], "-", t[i1], "s; max Mach", mach.max().round(3),
          "at", t[np.argmax(mach)], "s, alt", alt[np.argmax(mach)].round(0), "m")


if __name__ == "__main__":
    main()
