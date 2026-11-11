"""U.S. Standard Atmosphere 1976, 0-86 km (temperature, pressure, density, speed of sound, viscosity).

Coded from the published document (NOAA-S/T 76-1562, NASA-TM-X-74335, NTRS 19770009539;
literature/pdf/USSA1976.pdf), not from a third-party package:
  * defining constants: Table 2 (doc p. 2) with the corrections of the errata bound in the same PDF
    (p. 242 of the PDF: R* = 8.31432e3 N m/(kmol K) as on doc p. 3; r0 = 6356.766 km; S = 110.4 K);
  * sea-level composition and molecular weights: Table 3 (doc p. 3); M0 = 28.9644 kg/kmol (doc p. 9);
  * reference levels and gradients of the molecular-scale temperature: Table 4 (doc p. 3).
Below 80 km (Table 8: M/M0 = 1 at 79 km) the molecular-scale temperature T_M equals the kinetic temperature T (M = M0); the study
uses 0-60 km only, so the small 80-86 km molecular-weight correction (Table 8) is not applied and
altitudes above 79 km are refused.
Verified against Table I (geometric altitude, metric units) in tests/test_atmosphere.py.
"""
from __future__ import annotations

import numpy as np

R_STAR = 8.31432e3  # N m / (kmol K)  (doc p. 3; Table 2 misprints the exponent, see errata)
G0 = 9.80665  # m/s^2 (Table 2)
R0_EARTH = 6356.766e3  # m (Table 2 corrected by the errata: 6.356766e3 km)
P0 = 1.013250e5  # Pa (Table 2)
T0 = 288.15  # K (Table 2)
GAMMA = 1.40  # (Table 2)
SUTH_BETA = 1.458e-6  # kg / (s m K^1/2) (Table 2)
SUTH_S = 110.4  # K (Table 2 corrected by the errata)

# Table 3: (molecular weight kg/kmol, fractional volume)
COMPOSITION = {
    "N2": (28.0134, 0.78084), "O2": (31.9988, 0.209476), "Ar": (39.948, 0.00934),
    "CO2": (44.00995, 0.000314), "Ne": (20.183, 0.00001818), "He": (4.0026, 0.00000524),
    "Kr": (83.80, 0.00000114), "Xe": (131.30, 0.000000087), "CH4": (16.04303, 0.000002),
    "H2": (2.01594, 0.0000005),
}
M0 = 28.9644  # kg/kmol, sea-level mean molecular weight (doc p. 9, eq. 21; constant below 80 km, Table 8)
M0_TABLE3 = sum(m * f for m, f in COMPOSITION.values()) / sum(f for _, f in COMPOSITION.values())  # check value
R_AIR = R_STAR / M0  # J/(kg K)

# Table 4: geopotential base heights (km') and molecular-scale temperature gradients (K/km')
H_B = np.array([0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.8520]) * 1e3  # m'
L_B = np.array([-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0]) * 1e-3  # K/m'


def _base_values():
    """Temperature and pressure at the base of each layer, integrated from sea level."""
    T_b, P_b = [T0], [P0]
    for i in range(len(L_B)):
        dH = H_B[i + 1] - H_B[i]
        T1 = T_b[-1] + L_B[i] * dH
        if L_B[i] == 0.0:
            P1 = P_b[-1] * np.exp(-G0 * M0 * dH / (R_STAR * T_b[-1]))
        else:
            P1 = P_b[-1] * (T_b[-1] / T1) ** (G0 * M0 / (R_STAR * L_B[i]))
        T_b.append(T1)
        P_b.append(P1)
    return np.array(T_b), np.array(P_b)


T_BASE, P_BASE = _base_values()


def geopotential(z):
    """Geopotential height H (m') from geometric altitude z (m) (doc eq. 18)."""
    z = np.asarray(z, dtype=float)
    return R0_EARTH * z / (R0_EARTH + z)


def state(z):
    """Return (T [K], p [Pa], rho [kg/m^3]) at geometric altitude z [m] (scalar or array)."""
    z = np.asarray(z, dtype=float)
    if np.any(z < -5e3) or np.any(z > 79e3):
        raise ValueError("atmosphere model used only for -5 km <= z <= 79 km")
    H = geopotential(z)
    b = np.clip(np.searchsorted(H_B, H, side="right") - 1, 0, len(L_B) - 1)
    Tb, Pb, Lb, Hb = T_BASE[b], P_BASE[b], L_B[b], H_B[b]
    T = Tb + Lb * (H - Hb)
    with np.errstate(divide="ignore", invalid="ignore"):
        p_grad = Pb * (Tb / T) ** (G0 * M0 / (R_STAR * np.where(Lb == 0, 1.0, Lb)))
    p_iso = Pb * np.exp(-G0 * M0 * (H - Hb) / (R_STAR * Tb))
    p = np.where(Lb == 0.0, p_iso, p_grad)
    rho = p * M0 / (R_STAR * T)
    if T.ndim == 0:
        return float(T), float(p), float(rho)
    return T, p, rho


def speed_of_sound(z):
    T, _, _ = state(z)
    return np.sqrt(GAMMA * R_AIR * np.asarray(T))


def viscosity(T):
    """Sutherland's law of the Standard (doc eq. 51), Pa s; used here also above 1000 K, where it is
    an extrapolation (stated in INPUTS.md)."""
    T = np.asarray(T, dtype=float)
    return SUTH_BETA * T ** 1.5 / (T + SUTH_S)
