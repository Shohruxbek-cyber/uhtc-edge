"""Evaporation from the glass surface: Hertz-Knudsen-Langmuir kinetic limit in series with
boundary-layer diffusion.

Kinetic (free-evaporation) limit, from the author's fusion-wall code
(package ife_wall, module materials, class Vapor, method hkl_flux):
    j_HK = alpha p_eq(T) sqrt(M / (2 pi R T))                      [kg/(m^2 s)]
Diffusion limit through the boundary layer (film model with the heat/mass-transfer analogy,
DECISIONS.md D-HEAT-3): with the mass-transfer conductance g = rho_e u_e C_M = K_SG sqrt(p_s/R_n)
(the Sutton-Graves heat-transfer coefficient), the vapour mass fraction at the wall
y_w = (p_eq/p_s)(M_v/M_mix) and zero at the edge,
    j_D = g ln(1 + B),   B = y_w / (1 - y_w)                        (blowing form; -> g y_w for small y_w)
Series combination (the two resistances act on the same driving partial pressure, exact for the
linearised forms):  j = j_HK j_D / (j_HK + j_D).
When p_eq approaches p_s the liquid boils; y_w is capped at 0.99 and the rate is then set by the energy
balance in reality - not reached in the cases reported (flagged by the caller).
"""
from __future__ import annotations

import numpy as np

R = 8.314462618  # J/(mol K)
M_AIR = 0.028964  # kg/mol (US 1976 M0 = 28.9644 kg/kmol)


def hkl_flux(p_eq, T, M, alpha=1.0):
    """Hertz-Knudsen-Langmuir free-evaporation flux, kg/(m^2 s) (ife_wall.materials.Vapor.hkl_flux)."""
    return alpha * p_eq * np.sqrt(M / (2.0 * np.pi * R * T))


def diffusion_flux(p_eq, p_s, M, g):
    """Boundary-layer diffusion limit, kg/(m^2 s); g = mass-transfer conductance kg/(m^2 s)."""
    x = np.clip(p_eq / p_s, 0.0, 0.99)
    y = x * M / (x * M + (1 - x) * M_AIR)
    y = np.minimum(y, 0.99)
    return g * np.log1p(y / (1.0 - y))


def evaporation_flux(p_eq, T, p_s, M, g, alpha=1.0):
    jk = hkl_flux(p_eq, T, M, alpha)
    jd = diffusion_flux(p_eq, p_s, M, g)
    with np.errstate(invalid="ignore", divide="ignore"):
        j = np.where(jk + jd > 0, jk * jd / (jk + jd), 0.0)
    return j
