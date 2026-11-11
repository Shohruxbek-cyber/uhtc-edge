"""Literature data used by the model, with their sources (see INPUTS.md for the full provenance table).

Every array here is copied from the cited location; fits are computed from these arrays at import time so that
the numbers in the code can be traced back to the page they came from.
"""
from __future__ import annotations

import numpy as np

# ------------------------------------------------------------------------------------------------ properties
# NASA Ames TPSX database, material id 28 "Zirconium Diboride (ZrB2/SiC)", ZrB2 + 20 vol% SiC, ManLabs data
# (AFML-TR-68-190); literature/notes/TPSX2006_ZrB2-20SiC_ManLabs.md.  Thermal conductivity (property 3), +/-5 %
# to 2073 K, +/-8 % above.  At 1273.15 K the two methods (78.2 cut-bar, 82.0 laser flash) are averaged.
TPSX_K_T = np.array([373.15, 473.15, 573.15, 673.15, 773.15, 873.15, 973.15, 1073.15, 1173.15, 1273.15,
                     1373.15, 1473.15, 1573.15, 1673.15, 1773.15, 1873.15, 1973.15, 2073.15, 2173.15, 2273.15])
TPSX_K = np.array([98.7, 92.0, 88.3, 85.4, 84.5, 83.3, 81.6, 80.8, 79.5, 80.1,
                   79.5, 77.4, 75.7, 74.5, 73.2, 72.4, 71.5, 70.7, 69.9, 69.0])
# Density (property 1), same record: 5560 kg/m3 at 293 K, 5460 at 1273 K, 5310 at 2273 K (INDEX_properties 2.3)
TPSX_RHO_T = np.array([293.0, 1273.15, 2273.15])
TPSX_RHO = np.array([5560.0, 5460.0, 5310.0])

RHO_ZrB2 = 6090.0  # kg/m3, theoretical density, Shugart 2014 PhD thesis, PDF p. 16 ("theoretical density (6.09g/cm3)")
RHO_SiC = 3210.0  # kg/m3, TPSX K&K single crystal; PubChem lists 3.16-3.23 (literature/notes/PubChem2026_*.md)
RHO_SiO2_GLASS = 2200.0  # kg/m3, PubChem/HSDB (Merck Index), amorphous silica
RHO_B2O3_GLASS = 1800.0  # kg/m3, PubChem/HSDB (Merck Index), amorphous B2O3 ("1.8 g/cu cm (amorphous)")
ZrO2_EXPANSION = 0.169  # volume expansion ZrB2 -> ZrO2, Campbell et al. 2021 (AIAA 2021-3168) p. 7 ("16.9% for ZrO2")
K_SiO2_GLASS = 1.37  # W/m/K at 298 K, PubChem/HSDB citing Ullmann ("Thermal conductivity = 1.37 W/m-deg K at 298")
K_ZrO2_SCALE = 1.5  # W/m/K, porous zirconia, within 1.0-1.9 at 1316 C (Zhu & Miller 2002, NASA TM-2002-211481)

# Composition scaling of k (UNSURE, see INPUTS.md): 10 vol% from Van Laningham 2012 digitized ratios (10.7 vs 21.9
# vol% SiC: 67.6/58.7, 63.0/53.3, 56.9/47.5 -> 1.15-1.2); 30 vol%: 0.92, interpolated in the same Van Laningham data
# between 21.9 and 48.7 vol% SiC (Supplementary Section S2).  Zimmermann et al. 2008 (abstract: 56 W/m/K at 1675 K) against
# TPSX 20 % (74.5 W/m/K at 1673 K) would give 0.75; that ratio is not used (final review R1-20).
K_SCALE_SiC = {0.1: 1.17, 0.15: 1.085, 0.2: 1.0, 0.3: 0.92}  # 15 vol% (V7 material): linear between 10 and 20; 0.3: D-TH-4

# ------------------------------------------------------------------------------------------------ oxidation
# Shugart 2014 (UVA PhD), Table 9, PDF p. 75: ZrO2-layer parabolic constants (um^2/h) of ZrB2-30 vol% SiC.
# 1300-1500 C: box furnace, ambient (stagnant) air; 1650-1800 C: resistive heating in flowing O2 (900 sccm).
SHUGART_KP_LOW_T = np.array([1300.0, 1400.0, 1500.0]) + 273.15
SHUGART_KP_LOW = np.array([295.0, 280.0, 390.0])
SHUGART_KP_HIGH_T = np.array([1650.0, 1700.0, 1800.0]) + 273.15
SHUGART_KP_HIGH = np.array([8500.0, 28000.0, 42000.0])
# SiC-depletion onset 1627 +/- 25 C (Shugart 2014 pp. 79-87; Opila & Halbig 2001 depletion at 1627 C)
T_DEP_ONSET = 1627.0 + 273.15
# SiC-depleted layer at 1800 C (flowing O2): 43+/-10, 74+/-6, 104+/-8 um at 5, 15, 20 min (Shugart Fig. 44 caption,
# PDF p. 85); Shugart reports parabolic growth with no clear temperature dependence.
SHUGART_DEP_T_MIN = np.array([5.0, 15.0, 20.0])
SHUGART_DEP_UM = np.array([43.0, 74.0, 104.0])
# Opila & Halbig 2001 (Ceram. Eng. Sci. Proc. 22(3)), Table 1, PDF p. 6: measured recession kp'' (mm^2/min) of
# ZrB2-20 vol% SiC in stagnant air at 1327, 1627, 1927 C (used for the low-temperature slope and in V6).
OPILA_KPR_T = np.array([1327.0, 1627.0, 1927.0]) + 273.15
OPILA_KPR = np.array([7.9e-5, 1.3e-3, 2.2e-2])

# ------------------------------------------------------------------------------------------------ SiC transition
# Active-to-passive transition points of SiC in molecular O2 (45 points, 1473-1924 K) DIGITIZED from the
# literature compilation of Jacobson, Myers & Harder (2011 poster, NTRS); CSV in
# literature/notes/Jacobson2011_MoO8_poster_transition_points_DIGITIZED.csv (Rosner & Allendorf atomic-O data excluded).
# Passive-to-active lies 0.3-1 decade lower in p_O2 (Jacobson et al. 2012, Fig. 9, DIGITIZED).
PA_OFFSET_DECADES = 0.65  # midpoint of 0.3-1 (ASSUMPTION within the sourced range)


def _arrhenius(T, k):
    """Least-squares fit ln k = ln A - E/(R T); returns (A, E [J/mol])."""
    slope, intercept = np.polyfit(1.0 / np.asarray(T), np.log(np.asarray(k)), 1)
    return float(np.exp(intercept)), float(-slope * 8.314462618)


UM2H = 1e-12 / 3600.0  # um^2/h -> m^2/s
KP_LOW_A, KP_LOW_E = _arrhenius(SHUGART_KP_LOW_T, SHUGART_KP_LOW * UM2H)
KP_HIGH_A, KP_HIGH_E = _arrhenius(SHUGART_KP_HIGH_T, SHUGART_KP_HIGH * UM2H)
# below 1573 K the low-T law is continued with the activation energy of the Opila & Halbig recession constants
# between 1327 and 1627 C (the Shugart low-T data are nearly flat in T and must not be extrapolated downwards)
E_OPILA_LOW = float(-np.polyfit(1.0 / OPILA_KPR_T[:2], np.log(OPILA_KPR[:2]), 1)[0] * 8.314462618)
# depletion: x^2 = k_d t through the origin (no incubation: conservative)
_t = SHUGART_DEP_T_MIN * 60.0
_x2 = (SHUGART_DEP_UM * 1e-6) ** 2
K_DEP = float(np.sum(_x2 * _t) / np.sum(_t * _t))  # m^2/s


def _transition_fit():
    import csv
    from pathlib import Path
    f = Path(__file__).resolve().parents[2] / "literature" / "notes" / "Jacobson2011_MoO8_poster_transition_points_DIGITIZED.csv"
    rows = [r for r in csv.reader(l for l in open(f, encoding="utf-8") if not l.startswith("#"))][1:]
    T = np.array([float(r[2]) for r in rows])
    lp = np.array([float(r[4]) for r in rows])
    mol = np.array([r[0] != "RosnerAllendorf1970" for r in rows])
    b, a = np.polyfit(1.0 / T[mol], lp[mol], 1)
    rms = float(np.sqrt(np.mean((lp[mol] - (a + b / T[mol])) ** 2)))
    return float(a), float(b), rms, int(mol.sum())


TR_A, TR_B, TR_RMS, TR_N = _transition_fit()  # log10 p_O2,AP [bar] = TR_A + TR_B / T


def _transition_fit_fe():
    """Sensitivity case of the final review (R1-01): the same molecular-O2 points fitted with one common slope and one
    intercept per study (study fixed effects), so that offsets between studies do not flatten the slope.  The overall
    intercept makes the line pass through the centroid of all points, as the pooled fit does.  Returns (a, b, rms, n)."""
    import csv
    from pathlib import Path
    f = Path(__file__).resolve().parents[2] / "literature" / "notes" / "Jacobson2011_MoO8_poster_transition_points_DIGITIZED.csv"
    rows = [r for r in csv.reader(l for l in open(f, encoding="utf-8") if not l.startswith("#"))][1:]
    rows = [r for r in rows if r[0] != "RosnerAllendorf1970"]
    src = np.array([r[0] for r in rows])
    x = 1.0 / np.array([float(r[2]) for r in rows])
    y = np.array([float(r[4]) for r in rows])
    xd, yd = x.copy(), y.copy()
    for s in set(src):  # within-study deviations
        m = src == s
        xd[m] -= x[m].mean()
        yd[m] -= y[m].mean()
    b = float(np.sum(xd * yd) / np.sum(xd * xd))
    a = float(y.mean() - b * x.mean())
    res = yd - b * xd  # residuals about each study's own line
    return a, b, float(np.sqrt(np.mean(res ** 2))), len(rows)


TR_FE_A, TR_FE_B, TR_FE_RMS, TR_FE_N = _transition_fit_fe()

# ------------------------------------------------------------------ O2 permeability of the glass (criterion (a), D-FAIL-4)
# Molecular O2 permeability K (mol m^-1 s^-1 Pa^-1) as K = A exp(-Theta/T); flux through a film of thickness d is
# K p_O2 / d (Henry's law).  Sources and derivations: literature/notes/GlassOxygenTransport_silica.md,
# GlassOxygenTransport_borosilicate.md and GlassOxygenTransport_SiC-kp_protective-thickness.md; INPUTS.md section 8.
#   SiO2  base: implied by the CVD-SiC parabolic constant (measured 1473-1773 K on SiC);
#         low: Kajihara et al. 2004 D*S (eqs. 11-12, 1073-1373 K); high: Deal-Grove Si, dry O2 (Jacobson et al. 2005 Tab. II)
#   B2O3  base: O self-diffusion in the melt (Tokuda, via Schlichting 1982 Bild 3, 942-1010 K, E = 138 kJ/mol) times the
#         O2-equivalent density (mapping checked on silica at 750 C); high: Luthra's value as printed in Campbell et
#         al. 2021 eq. 7 (sign-typo reading, UNSURE); low = base.
# Mixed glass: linear in the B2O3 mole fraction (Campbell et al. 2021 eqs. 10-11; reproduces Schlichting Bild 4 within
# x1.5 at 750 C).  All values above ~1800 K are extrapolations.
O2_PERM = {
    "base": dict(SiO2=(2.15e-14, 14169.0), B2O3=(7.7e-7, 16600.0)),
    "low": dict(SiO2=(2.83e-16, 8471.0), B2O3=(7.7e-7, 16600.0)),
    "high": dict(SiO2=(6.50e-14, 14914.0), B2O3=(0.15 / 101325.0, 16000.0)),
}
Y_O2_AIR = 0.2095 * 0.0319988 / 0.0289644  # mass fraction of oxygen in air counted as O2 (D-OX-1)

# ------------------------------------------------------------------ incubation of SiC depletion (model option, D-OX-10)
# Shugart 2014 Fig. 40 (PDF p. 81) and text pp. 79-80: ZrB2-30SiC in flowing O2 (~1 atm).  Onset brackets: 1627 C
# <= 40 min, 1650 C 10-20 min, 1700 C 10-15 min, 1800 C <= 5 min (literature/digitized/Shugart2014_fig40_...csv,
# literature/notes/SiC_depletion_incubation.md).  The values below are the bracket midpoints (0-40 -> 20, 10-20 -> 15,
# 10-15 -> 12.5, 0-5 -> 2.5 min): our construction, not a published law; held constant beyond the end points.
INC_T = np.array([1627.0, 1650.0, 1700.0, 1800.0]) + 273.15
INC_MIN = np.array([20.0, 15.0, 12.5, 2.5])

# upper end of the high-temperature k_p fit (Shugart 1800 C); beyond it the law is extrapolated (D-OX-11)
KP_FIT_TMAX = 1800.0 + 273.15

# passive-to-active transition measured on SiC-based composites (IXV C/SiC, "HER") in air plasma and a solar furnace,
# p_O2 = A exp(-B/T) with p_O2 in Pa, fitted at 42-1050 Pa O2 (Balat-Pichelin et al. 2015, PDF pp. 4-5; HAL hal-02313822;
# literature/notes/CHECK_2_jump_and_active_oxidation_T_p.md).  Sensitivity case of criterion (b) only (D-OX-12).
PA_BALAT_A, PA_BALAT_B = 1.0e15, 55850.0
