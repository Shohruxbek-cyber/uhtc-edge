"""Stagnation-point convective heating and its distribution around a blunt nose.

Correlations (sources in INPUTS.md, section "Heating"):
  * Sutton & Graves (1971), NASA TR R-376, eq. (33) with the air coefficient of Table II (p. 39):
        q = K sqrt(p_s / R_n) (H_s - h_w),   K = 0.1113 kg s^-1 m^-3/2 atm^-1/2,  p_s in atm.
    Derived for an equilibrium boundary layer (wall at 300 and 1111 K, H_s 2.3-116.2 MJ/kg; abstract),
    i.e. effectively a fully catalytic wall.  The factor (H_s - h_w) is the hot-wall correction.
  * Fay & Riddell (1958) as reproduced by Hamilton (1982), NASA TM-84580, eq. (6):
        q = 0.76 Pr^-0.6 (rho_w mu_w)^0.1 (rho_e mu_e)^0.4 [1 + (Le^0.52 - 1) h_D/H_s] (H_s - h_w) sqrt(du_e/dx)
    with the Newtonian velocity gradient of Tauber (1989), NASA TP-2914, eq. (39a):
        du_e/dx = (1/R_n) sqrt(2 (p_s - p_inf) / rho_s).
  * Tauber (1989) eq. (40) (sphere, simplified) and eq. (41) (unswept cylinder):
        q_sph = 1.83e-4 sqrt(rho_inf/R_n) V^3 (1 - h_w/H_s),   q_cyl = 1.29e-4 sqrt(rho_inf/R_n) V^3 (1 - h_w/h_aw).
  * Distribution on a hemisphere or a cylinder normal to the stream: q/q_s = cos(theta), which Tauber
    (1989, p. 8, Figs. 5-6) shows to follow Lees' and Kemp et al.'s solutions up to about 70 deg.
  * Flank of a sharp cone / wedge (laminar): Tauber (1989) eq. (46) (cone) and eq. (48) (wedge, flat plate).

Catalysis bounds (DECISIONS.md D-HEAT-2): the fully catalytic value is the Sutton-Graves result with the wall
enthalpy of equilibrium air at (T_w, p_s).  The non-catalytic bound follows Kolodziej (1997, NASA TM-112204,
doc p. 3): "noncatalytic hot wall (NCHW) heating ... is obtained by setting Le = 0" in the Fay-Riddell bracket,
i.e. the bracket [1 + (Le^0.52 - 1) h_D/H_s] becomes (1 - h_D/H_s):
        q_nc = q_fc (1 - h_D/H_s) / [1 + (Le^0.52 - 1) h_D/H_s],
with the wall enthalpy of undissociated (frozen) air at T_w.
"""
from __future__ import annotations

import numpy as np

ATM = 101325.0  # Pa
K_SG_AIR = 0.1113  # kg / (s m^1.5 atm^0.5), Sutton & Graves 1971, Table II (air: 0.2320 O2 - 0.7680 N2)
C_TAUBER_SPHERE = 1.83e-4  # SI, Tauber 1989 eq. (40)
C_TAUBER_CYL = 1.29e-4  # SI, Tauber 1989 eq. (41) (sweep angle 0)
C_TAUBER_CONE = 4.03e-5  # SI, Tauber 1989 eq. (46)
C_TAUBER_WEDGE = 2.42e-5  # SI, Tauber 1989 eq. (48)
FR_COEF = 0.76  # Hamilton 1982 eq. (6)


def sutton_graves(p_s, R_n, H_s, h_w, K=K_SG_AIR):
    """Stagnation-point heat flux (W/m^2) on a sphere of nose radius R_n (m); p_s in Pa, enthalpies J/kg."""
    return K * np.sqrt(np.asarray(p_s) / ATM / R_n) * (np.asarray(H_s) - np.asarray(h_w))


def sutton_graves_cylinder(p_s, R_n, H_s, h_w, K=K_SG_AIR):
    """2D stagnation line of an unswept cylinder: the sphere value times C_cyl/C_sph of Tauber (1989)
    eqs. (40)-(41) (= 0.705, close to 1/sqrt(2) for the Newtonian velocity gradients of sphere and cylinder)."""
    return sutton_graves(p_s, R_n, H_s, h_w, K) * (C_TAUBER_CYL / C_TAUBER_SPHERE)


def newtonian_velocity_gradient(R_n, p_s, p_inf, rho_s):
    """du_e/dx at the stagnation point, 1/s (Tauber 1989 eq. 39a)."""
    return np.sqrt(2.0 * (p_s - p_inf) / rho_s) / R_n


def fay_riddell(rho_w, mu_w, rho_e, mu_e, H_s, h_w, h_D, due_dx, Pr=0.71, Le=1.4):
    """Fay-Riddell equilibrium stagnation-point heat flux, W/m^2 (Hamilton 1982 eq. 6). SI units."""
    bracket = 1.0 + (Le ** 0.52 - 1.0) * h_D / H_s
    return (FR_COEF * Pr ** -0.6 * (rho_w * mu_w) ** 0.1 * (rho_e * mu_e) ** 0.4 * bracket
            * (H_s - h_w) * np.sqrt(due_dx))


def noncatalytic_factor(H_s, h_D, Le=1.4):
    """Ratio q_noncatalytic / q_fully_catalytic (Fay-Riddell bracket with Le = 0 over the bracket with Le),
    Kolodziej 1997 p. 3; see module doc."""
    return (1.0 - h_D / H_s) / (1.0 + (Le ** 0.52 - 1.0) * h_D / H_s)


def tauber_sphere(rho_inf, V, R_n, h_w, H_s):
    return C_TAUBER_SPHERE * np.sqrt(rho_inf / R_n) * V ** 3 * (1.0 - h_w / H_s)


def tauber_cone(rho_inf, V, x, half_angle, h_w, h_aw):
    """Laminar heating on a sharp cone at distance x (m) from the tip (Tauber 1989 eq. 46)."""
    return (C_TAUBER_CONE * np.sqrt(rho_inf * np.cos(half_angle) / x) * V ** 3.2 * np.sin(half_angle)
            * (1.0 - h_w / h_aw))


def tauber_wedge(rho_inf, V, x, half_angle, h_w, h_aw):
    """Laminar heating on a wedge / inclined flat plate at distance x (m) (Tauber 1989 eq. 48)."""
    return (C_TAUBER_WEDGE * np.sqrt(rho_inf * np.cos(half_angle) / x) * V ** 3.2 * np.sin(half_angle)
            * (1.0 - h_w / h_aw))


def cosine_distribution(theta):
    """q/q_s on a hemisphere or a cylinder normal to the stream at body angle theta (rad) from the
    stagnation point (Tauber 1989 p. 8 and Figs. 5-6; valid to about 70 deg)."""
    return np.maximum(np.cos(theta), 0.0)
