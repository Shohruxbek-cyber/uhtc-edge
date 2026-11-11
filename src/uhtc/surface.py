"""Surface energy balance with a thin oxide scale in series (shared by the 1D and 2D solvers).

Energy balance at the outer surface (temperature T_s), per unit area:
    q_conv(T_s) - eps(T_s) sigma T_s^4 + q_chem(T_s) - (T_s - T_sub)/R = 0
    q_conv(T_s)                              hot-wall convective flux supplied by the caller (uhtc.coupled:
                                              Phi (H_s - h_w(T_s)) with the catalysis bound applied)
    q_chem                                   reaction heat minus evaporation enthalpy (from the oxidation
                                              model), with its derivative
    R = R_scale (+ half-cell of the 1D mesh) thermal resistance between the outer surface and the
                                              temperature node of the solid, K m^2/W.
The scale's heat capacity is neglected (quasi-steady scale: rho c delta << that of the heated substrate;
DECISIONS.md D-TH-3).  Returns the flux into the solid node and its derivative with respect to the node
temperature (implicit-function theorem), for the Newton iteration of the conduction solver.
"""
from __future__ import annotations

SIGMA = 5.670374419e-8  # W/(m^2 K^4), CODATA 2018


def _F(Ts, conv_fn, eps_fn, chem_fn):
    e, de = eps_fn(Ts)
    qc, dqc = conv_fn(Ts)
    f = qc - e * SIGMA * Ts ** 4
    df = dqc - 4 * e * SIGMA * Ts ** 3 - de * SIGMA * Ts ** 4
    if chem_fn is not None:
        c, dc = chem_fn(Ts)
        f += c
        df += dc
    return f, df


def balance(T_node, R, conv_fn, eps_fn, chem_fn=None, T_bounds=(150.0, 5900.0)):
    """Solve for T_s.  Returns (q_into_solid, dq/dT_node, T_s).

    conv_fn(T) -> (q_conv, dq_conv/dT); eps_fn(T) -> (eps, d eps/dT);
    chem_fn(T) -> (q_chem, dq_chem/dT) or None."""
    if R <= 0.0:
        f, df = _F(T_node, conv_fn, eps_fn, chem_fn)
        return f, df, T_node
    # G(Ts) = F(Ts) - (Ts - T_node)/R is strictly decreasing: safeguarded Newton on a bracket
    lo, hi = T_bounds
    Ts = T_node
    for _ in range(100):
        f, df = _F(Ts, conv_fn, eps_fn, chem_fn)
        g = f - (Ts - T_node) / R
        dg = df - 1.0 / R
        if g > 0:
            lo = max(lo, Ts)
        else:
            hi = min(hi, Ts)
        step = -g / dg
        Tn = Ts + step
        if not (lo < Tn < hi):
            Tn = 0.5 * (lo + hi)
        if abs(Tn - Ts) < 1e-11 * max(Ts, 1.0):  # tight: q = (Ts - T_node)/R with R ~ 1e-7 m2K/W
            Ts = Tn
            break
        Ts = Tn
    f, df = _F(Ts, conv_fn, eps_fn, chem_fn)
    q = (Ts - T_node) / R
    dTs = (1.0 / R) / (-df + 1.0 / R)
    return q, (dTs - 1.0) / R, Ts


def constant_eps(e):
    return lambda T: (e, 0.0)


def radiative_equilibrium(conv_fn, eps_fn, chem_fn=None):
    """Surface temperature with no conduction into the solid (R -> infinity limit)."""
    lo, hi = 150.0, 5900.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        f, _ = _F(mid, conv_fn, eps_fn, chem_fn)
        if f > 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-6:
            break
    return 0.5 * (lo + hi)
