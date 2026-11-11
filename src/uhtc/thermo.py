"""Thermochemistry from the NIST-JANAF tables (NIST SRD 13, public domain; literature/janaf/*.txt).

Each table gives, per temperature, Cp, S, -[G-H(Tr)]/T, H-H(Tr) (kJ/mol), delta-f H, delta-f G (kJ/mol)
and log Kf, with Tr = 298.15 K and a standard pressure of 1 bar.  Values are interpolated linearly in T
(tables at 100 K steps).  "Absolute" enthalpies use the formation basis: H(T) = delta-f H(298.15) +
[H - H(Tr)](T), so that elements in their reference state have H(298.15) = 0.

Uses: enthalpy of undissociated air (hot-wall correction), equilibrium air (edge state, dissociation
enthalpy h_D), vapor pressures over the glass (B2O3, SiO, SiO2) and reaction enthalpies.  Argon is
treated as an ideal monatomic gas (cp = 5/2 R), exact for an ideal gas, because no Ar table is used.
"""
from __future__ import annotations

from contextlib import contextmanager
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

R = 8.314462618  # J/(mol K), CODATA 2018
P_STD = 1.0e5  # Pa, JANAF standard state (1 bar)
JANAF_DIR = Path(__file__).resolve().parents[2] / "literature" / "janaf"

# standard atomic weights, abridged (CIAAW, https://www.ciaaw.org/abridged-atomic-weights.htm, read 2 Nov 2026)
ATOMIC_WEIGHT = {"N": 14.007, "O": 15.999, "Si": 28.085, "C": 12.011, "B": 10.81, "Zr": 91.222, "Ar": 39.95}
_FORMULA = {"N2": {"N": 2}, "O2": {"O": 2}, "N": {"N": 1}, "O": {"O": 1}, "NO": {"N": 1, "O": 1},
            "Ar": {"Ar": 1}, "B2O3": {"B": 2, "O": 3}, "SiO2": {"Si": 1, "O": 2}, "SiO": {"Si": 1, "O": 1},
            "SiC": {"Si": 1, "C": 1}, "CO": {"C": 1, "O": 1}, "CO2": {"C": 1, "O": 2},
            "ZrB2": {"Zr": 1, "B": 2}, "ZrO2": {"Zr": 1, "O": 2}}
MOLAR_MASS = {sp: 1e-3 * sum(n * ATOMIC_WEIGHT[el] for el, n in f.items()) for sp, f in _FORMULA.items()}  # kg/mol

FILES = {"N2": "N2_ref", "O2": "O2_ref", "N": "N_g", "O": "O_g", "NO": "NO_g", "CO": "CO_g",
         "CO2": "CO2_g", "SiO": "SiO_g", "SiO2(g)": "SiO2_g", "B2O3(g)": "B2O3_g",
         "B2O3(l)": "B2O3_l", "B2O3(cr,l)": "B2O3_cr-l", "SiO2(l)": "SiO2_l", "SiO2(cr,l)": "SiO2_cr-l",
         "SiC": "SiC_cr_beta", "ZrB2": "ZrB2_cr-l", "ZrO2": "ZrO2_cr-l"}


# Opt-in extension of the tables above their end (6000 K for the air species), for ground-test edge states only
# (PrescribedConditions): Cp is held at its last tabulated value, H - H(Tr) grows with that Cp and delta-f G is
# extended linearly with the slope of the last table interval.  Neutral species only (no ionisation).  The error of
# the extension is checked in tests/test_thermo.py by extending tables cut at 5000 K to 6000 K.  0 = strict.
_T_EXT_MAX = 0.0


@contextmanager
def extended_tables(T_max: float = 9000.0):
    """Allow table look-ups up to T_max (K) inside the with-block; outside it the tables stay strict."""
    global _T_EXT_MAX
    old = _T_EXT_MAX
    _T_EXT_MAX = float(T_max)
    try:
        yield
    finally:
        _T_EXT_MAX = old


class Table:
    def __init__(self, path: Path):
        rows = []
        with open(path, encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()
        self.title = lines[0]
        for line in lines[2:]:
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            try:
                vals = [float(x) for x in parts[:7]]
            except ValueError:  # TRANSITION rows or INFINITE entries
                continue
            rows.append(vals)
        a = np.array(rows)
        # transition rows appear twice at the same T (both phases); keep the higher-T phase value
        T = a[:, 0]
        keep = np.r_[T[1:] != T[:-1], True]
        a = a[keep]
        self.T, self.Cp, self.S = a[:, 0], a[:, 1], a[:, 2]
        self.dH = a[:, 4] * 1e3  # J/mol, H - H(298.15)
        self.dfH = a[:, 5] * 1e3
        self.dfG = a[:, 6] * 1e3
        self.dfH298 = float(np.interp(298.15, self.T, self.dfH))

    def _chk(self, T):
        T = np.asarray(T, float)
        if np.any(T < self.T[0]) or np.any(T > max(self.T[-1], _T_EXT_MAX)):
            raise ValueError(f"{self.title}: T outside {self.T[0]}-{self.T[-1]} K")
        return T

    def _above(self, T):
        """Excess of T over the table end (0 inside the table)."""
        return np.maximum(T - self.T[-1], 0.0)

    def H(self, T):
        """Absolute (formation-basis) molar enthalpy, J/mol."""
        T = self._chk(T)
        return self.dfH298 + np.interp(T, self.T, self.dH) + self.Cp[-1] * self._above(T)

    def G_f(self, T):
        """Gibbs energy of formation from the elements, J/mol."""
        T = self._chk(T)
        slope = (self.dfG[-1] - self.dfG[-2]) / (self.T[-1] - self.T[-2])
        return np.interp(T, self.T, self.dfG) + slope * self._above(T)

    def cp(self, T):
        return np.interp(self._chk(T), self.T, self.Cp)


@lru_cache(maxsize=None)
def table(species: str) -> Table:
    f = JANAF_DIR / f"{FILES[species]}.txt"
    if not f.exists():
        raise FileNotFoundError(f"{f} is missing: run 'python literature/janaf/fetch_janaf.py' first (NIST-JANAF tables are not redistributed)")
    return Table(f)


def H_species(sp, T):
    if sp == "Ar":
        return 2.5 * R * (np.asarray(T, float) - 298.15)
    return table(sp).H(T)


# ---------------------------------------------------------------------------------- air
# mole fractions of undissociated dry air: N2, O2, Ar of US 1976 Table 3, CO2 and traces lumped into Ar
X_AIR = {"N2": 0.78084, "O2": 0.209476, "Ar": 1.0 - 0.78084 - 0.209476}
M_AIR = sum(X_AIR[s] * MOLAR_MASS[s] for s in X_AIR)
N_OVER_O = 2 * X_AIR["N2"] / (2 * X_AIR["O2"])
AR_OVER_O = X_AIR["Ar"] / (2 * X_AIR["O2"])
SPECIES = ("N2", "O2", "NO", "N", "O", "Ar")


def h_air_frozen(T):
    """Specific enthalpy of undissociated air, J/kg (formation basis: 0 at 298.15 K)."""
    T = np.asarray(T, float)
    return sum(X_AIR[s] * H_species(s, T) for s in X_AIR) / M_AIR


def cp_air_frozen(T, dT=0.5):
    return (h_air_frozen(T + dT) - h_air_frozen(T - dT)) / (2 * dT)


def T_air_frozen(h):
    return brentq(lambda T: h_air_frozen(T) - h, 100.0, 5999.0, xtol=1e-6)


def _K(sp, T):
    """Equilibrium constant of formation of sp from its elements (bar units)."""
    return np.exp(-table(sp).G_f(T) / (R * T))


def equilibrium_air(T, p):
    """Mole fractions of N2, O2, NO, N, O, Ar for air at T (K) and p (Pa) in chemical equilibrium
    (neutral species only; ionisation is negligible below about 9000 K edge temperatures)."""
    P = p / P_STD
    KO, KN, KNO = _K("O", T), _K("N", T), _K("NO", T)

    def solve_given_a(a):
        # a = sqrt(p_O2) (bar^0.5); element ratio N/O fixes b = sqrt(p_N2)
        pO2, pO = a * a, KO * a
        nO = 2 * pO2 + pO

        def f(b):
            return (2 * b * b + KN * b + KNO * a * b) - N_OVER_O * (nO + KNO * a * b)
        bhi = 1.0
        while f(bhi) < 0:
            bhi *= 2
        b = brentq(f, 0.0, bhi, xtol=1e-16, rtol=1e-14)
        pN2, pN, pNO = b * b, KN * b, KNO * a * b
        pAr = AR_OVER_O * (2 * pO2 + pO + pNO)
        return np.array([pN2, pO2, pNO, pN, pO, pAr])

    def total(la):
        return np.log(solve_given_a(np.exp(la)).sum() / P)
    la = brentq(total, -60.0, 10.0, xtol=1e-13)
    pp = solve_given_a(np.exp(la))
    return dict(zip(SPECIES, pp / pp.sum()))


def mixture_props(x: dict, T):
    """(molar mass kg/mol, specific enthalpy J/kg, chemical enthalpy h_D J/kg) of a mixture."""
    M = sum(x[s] * MOLAR_MASS[s] for s in x)
    h = sum(x[s] * H_species(s, T) for s in x) / M
    hD = sum(x[s] * (table(s).dfH298 if s != "Ar" else 0.0) for s in x) / M
    return M, h, hD


def h_air_equilibrium(T, p):
    return mixture_props(equilibrium_air(T, p), T)[1]


def T_air_equilibrium(h, p, T_lo=150.0, T_hi=5999.0):
    """Temperature of equilibrium air with specific enthalpy h (J/kg) at pressure p (Pa)."""
    return brentq(lambda T: h_air_equilibrium(T, p) - h, T_lo, T_hi, xtol=1e-4)


# ---------------------------------------------------------------------------------- condensed / vapor
def dG_reaction(T, reactants: dict, products: dict):
    """Delta-r G (J/mol) from the delta-f G columns; keys are table names, values stoichiometry.
    Elements in their reference state (O2, N2) have delta-f G = 0 by construction."""
    T = np.asarray(T, float)
    g = sum(n * table(s).G_f(T) for s, n in products.items())
    g -= sum(n * table(s).G_f(T) for s, n in reactants.items())
    return g


def dH_reaction(T, reactants: dict, products: dict):
    T = np.asarray(T, float)
    h = sum(n * table(s).H(T) for s, n in products.items())
    h -= sum(n * table(s).H(T) for s, n in reactants.items())
    return h


def p_B2O3(T, activity=1.0, liquid="B2O3(l)"):
    """Equilibrium partial pressure (Pa) of B2O3(g) over B2O3 liquid of the given activity."""
    K = np.exp(-dG_reaction(T, {liquid: 1}, {"B2O3(g)": 1}) / (R * np.asarray(T, float)))
    return activity * K * P_STD


def p_SiO_over_SiO2(T, p_O2, activity=1.0, liquid="SiO2(l)"):
    """SiO2(l) = SiO(g) + 1/2 O2(g): equilibrium p_SiO (Pa) for a given oxygen partial pressure (Pa)."""
    K = np.exp(-dG_reaction(T, {liquid: 1}, {"SiO": 1}) / (R * np.asarray(T, float)))
    return activity * K / np.sqrt(np.maximum(p_O2, 1e-30) / P_STD) * P_STD


def p_SiO2_gas(T, activity=1.0, liquid="SiO2(l)"):
    K = np.exp(-dG_reaction(T, {liquid: 1}, {"SiO2(g)": 1}) / (R * np.asarray(T, float)))
    return activity * K * P_STD


def L_vap_B2O3(T):
    """Enthalpy of vaporisation of liquid B2O3, J/kg."""
    return dH_reaction(T, {"B2O3(l)": 1}, {"B2O3(g)": 1}) / MOLAR_MASS["B2O3"]
