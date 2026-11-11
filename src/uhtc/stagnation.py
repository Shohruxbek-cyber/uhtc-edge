"""Free stream -> stagnation-point state (total enthalpy, stagnation pressure, boundary-layer edge,
oxygen partial pressure at the wall).

* Free stream from the US 1976 atmosphere; V = M a_inf with a_inf = sqrt(gamma R T) (gamma = 1.4).
* Total enthalpy H_s = h(T_inf) + V^2/2 (energy conservation, exact for any gas), with h of undissociated
  air from JANAF (uhtc.thermo).
* Stagnation pressure, two options:
    'ideal'       Rayleigh pitot formula for a calorically perfect gas, gamma = 1.4 (Hamilton 1982 eq. 10
                  gives the ideal-gas shock density ratio used in the same way);
    'equilibrium' normal shock in equilibrium air (Rankine-Hugoniot with uhtc.thermo.equilibrium_air),
                  then p_s = p_inf + rho_inf V^2 (1 - eps/2), eps = rho_inf/rho_2 (Tauber 1989, p. 6).
  The production runs use 'equilibrium'; the ideal value is kept for the comparison in V4.
* Edge state: equilibrium air at (p_s, H_s): T_e, rho_e, mu_e (Sutherland, US 1976 eq. 51, extrapolated),
  dissociation (chemical) enthalpy h_D.
* Oxygen at the wall: p_O2,w = x_O2 p_s with x_O2 = 0.2095, i.e. all oxygen counted as O2
  (DECISIONS.md D-OX-1: atoms are assumed to recombine at the wall or to be equivalent to O2 for oxidation).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from . import atmosphere as atm
from . import thermo as th

GAMMA = 1.4
X_O2 = th.X_AIR["O2"]


@dataclass
class StagnationState:
    mach: float
    altitude: float
    V: float
    T_inf: float
    p_inf: float
    rho_inf: float
    H_s: float  # J/kg (formation basis)
    p_s: float  # Pa
    p_s_ideal: float  # Pa
    eps: float  # rho_inf / rho_2 across the normal shock
    T_e: float  # K, equilibrium edge temperature at (p_s, H_s)
    rho_e: float
    mu_e: float
    h_D: float  # J/kg, chemical (dissociation) enthalpy at the edge
    x_e: dict
    p_O2_w: float  # Pa


def pitot_ideal(M, p_inf, gamma=GAMMA):
    """Rayleigh pitot formula: stagnation pressure behind a normal shock, ideal gas."""
    a = ((gamma + 1) ** 2 * M ** 2 / (4 * gamma * M ** 2 - 2 * (gamma - 1))) ** (gamma / (gamma - 1))
    b = (1 - gamma + 2 * gamma * M ** 2) / (gamma + 1)
    return p_inf * a * b


def normal_shock_equilibrium(V, p1, rho1, h1):
    """Equilibrium-air normal shock. Returns (eps, p2, h2, T2)."""
    def resid(eps):
        p2 = p1 + rho1 * V ** 2 * (1 - eps)
        h2 = h1 + 0.5 * V ** 2 * (1 - eps ** 2)
        T2 = th.T_air_equilibrium(h2, p2)
        M2, _, _ = th.mixture_props(th.equilibrium_air(T2, p2), T2)
        rho2 = p2 * M2 / (th.R * T2)
        return rho1 / rho2 - eps
    eps = brentq(resid, 0.02, 0.9, xtol=1e-7)
    p2 = p1 + rho1 * V ** 2 * (1 - eps)
    h2 = h1 + 0.5 * V ** 2 * (1 - eps ** 2)
    return eps, p2, h2, th.T_air_equilibrium(h2, p2)


def stagnation_state(mach, altitude, pressure_model="equilibrium"):
    T, _, _ = atm.state(altitude)
    return stagnation_state_V(mach * np.sqrt(GAMMA * atm.R_AIR * T), altitude, pressure_model)


def stagnation_state_V(V, altitude, pressure_model="equilibrium"):
    """As stagnation_state, for a given flight velocity V (m/s)."""
    T, p, rho = atm.state(altitude)
    mach = V / np.sqrt(GAMMA * atm.R_AIR * T)
    h_inf = float(th.h_air_frozen(T))
    H_s = h_inf + 0.5 * V ** 2
    p_ideal = float(pitot_ideal(mach, p))
    eps, _, _, _ = normal_shock_equilibrium(V, p, rho, h_inf)
    p_eq = p + rho * V ** 2 * (1 - eps / 2)
    p_s = p_eq if pressure_model == "equilibrium" else p_ideal
    T_e = th.T_air_equilibrium(H_s, p_s)
    x = th.equilibrium_air(T_e, p_s)
    M_e, _, h_D = th.mixture_props(x, T_e)
    rho_e = p_s * M_e / (th.R * T_e)
    return StagnationState(mach, altitude, V, T, p, rho, H_s, p_s, p_ideal, eps, T_e, rho_e,
                           float(atm.viscosity(T_e)), h_D, x, X_O2 * p_s)
