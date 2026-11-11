"""Thin-film flow of the outer glass over the 2D surface (version 3, Task 1). Off unless run_2d gets a FlowParams.

Lubrication approximation on the wetted coordinate s of the surface nodes, flow from the stagnation point towards the flank:

    q(s) = tau_w(s) delta^2 / (2 mu(T_s)) + (-dp/ds) delta^3 / (3 mu(T_s))          (volume flux per unit width)
    d(delta)/dt + (1/r) d(r q)/ds = S - w_evap

Operator splitting inside one oxide step: LayerOxidation.advance applies the source S and the evaporation w_evap node by node
(one source term, not duplicated), then film_step moves glass between nodes over the same dt. The flux form on the lumped surface
weights of the conduction solver (Solver2D.w_surf, the integral of r ds per node, 2 pi dropped) conserves glass volume exactly;
glass leaving the last flank node is "shed" and tallied. First-order upwind in the flow direction; explicit sub-steps under a CFL
limit, or a linearised implicit upwind sweep when more than max_explicit sub-steps would be needed (positive and conservative
either way). The composition of each node is kept: a node keeps its x_B when it loses glass, and glass arriving on a bare node
takes the composition of the node it came from (D-FLOW-6).

Decisions (DECISIONS.md): D-FLOW-3 wall shear (Reynolds analogy, Pr^(2/3), band x0.5-2), D-FLOW-4 pressure gradient
(modified Newtonian), D-FLOW-5 viscosity at the local surface temperature (Mauro et al. 2022, LA and F320), D-FLOW-6 composition
not advected, D-FLOW-7 numerics.
"""
import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

R_GAS = 8.314462618
# Mauro et al. 2022, Table II (fibre-draw fits, 1973-2573 K): log10 mu [Pa s] = log10 mu_inf + E / (ln 10 R T)
VISC = {"LA": (587.4e3, -8.60), "F320": (495.8e3, -6.80)}
PR = 0.71  # Prandtl number of the Reynolds analogy (as the Supplement S9 estimate)
GAMMA = 1.4  # isentropic expansion from the stagnation point (as the S9 estimate and the Lees distribution)


@dataclass(frozen=True)
class FlowParams:
    law: str = "LA"  # viscosity law, "LA" (natural silica, baseline) or "F320" (fluorine-doped, low-viscosity case)
    tau_mult: float = 1.0  # band on the wall shear, x0.5-2 (D-FLOW-3)
    pressure: bool = True  # pressure-gradient term (sensitivity: off, D-FLOW-4)
    cfl: float = 0.9
    max_explicit: int = 10_000  # above this many explicit sub-steps per oxide step: linearised implicit sweep (D-FLOW-7)


def viscosity(T, law="LA"):
    """Silica-melt viscosity (Pa s) at T (K)."""
    E, l0 = VISC[law]
    return 10.0 ** (l0 + E / (math.log(10.0) * R_GAS * np.asarray(T, float)))


@dataclass
class Surface:
    """Surface nodes of the 2D solver, ordered from the stagnation point."""
    s: np.ndarray  # wetted coordinate of the nodes, m
    r: np.ndarray  # body radius of the nodes (1 for a planar body), m
    w: np.ndarray  # lumped surface weight of each node, integral of r ds (Solver2D.w_surf), m^2 per radian
    phi: np.ndarray  # body angle of the nodes, rad
    r_face: np.ndarray  # body radius at the n-1 internal faces and the outflow face at the end of the flank (n values)
    R_n: float
    axisymmetric: bool


def surface_from_solver(sol):
    sn = sol.nodes[sol.surf]
    axi = sol.g.axisymmetric
    r = sn[:, 1].copy() if axi else np.ones(len(sol.surf))
    r_face = np.r_[0.5 * (r[:-1] + r[1:]), r[-1]]
    return Surface(np.asarray(sol.s, float), r, np.asarray(sol.w_surf, float), np.asarray(sol.phi, float), r_face, sol.g.R_n, axi)


def pressure_ratio(phi, eps):
    """Modified-Newtonian wall pressure over stagnation pressure, as the Lees distribution and local_pressure use."""
    return (1.0 - eps) * np.cos(phi) ** 2 + eps


def drivers(surf: Surface, flow, gshape, rho_e, p_inf, params: FlowParams):
    """Wall shear tau (Pa) at the faces and G = -dp/ds (Pa/m) at the faces. Faces 0..n-2 lie between nodes i and i+1, face
    n-1 is the outflow face at the end of the flank (no pressure gradient there)."""
    eps = min(max(p_inf / flow.p_s, 0.0), 1.0)
    pr = pressure_ratio(surf.phi, eps)
    ue = np.sqrt(np.maximum(2.0 * GAMMA / (GAMMA - 1.0) * flow.p_s / rho_e * (1.0 - pr ** ((GAMMA - 1.0) / GAMMA)), 0.0))
    tau_node = params.tau_mult * flow.Phi * np.asarray(gshape, float) * PR ** (2.0 / 3.0) * ue  # Reynolds analogy
    tau = np.r_[0.5 * (tau_node[:-1] + tau_node[1:]), tau_node[-1]]
    if params.pressure:
        p = flow.p_s * pr
        G = np.r_[-np.diff(p) / np.diff(surf.s), 0.0]
    else:
        G = np.zeros(len(surf.s))
    return np.maximum(tau, 0.0), np.maximum(G, 0.0)


def film_step(delta, surf: Surface, tau, G, mu_node, dt, params: FlowParams):
    """Move the outer glass (volume per area delta, m, at the nodes) over dt. Returns (new delta, shed volume per radian over
    the step, number of explicit sub-steps or -1 for the implicit sweep). Upwind: the face flux uses the thickness and the
    viscosity of the node upstream of it."""
    d = np.array(delta, float)
    w, rf = surf.w, surf.r_face
    a = rf * tau / (2.0 * mu_node)  # face flux F = a d^2 + b d^3 of the upstream node
    b = rf * G / (3.0 * mu_node)
    shed, t, n_sub = 0.0, 0.0, 0
    while t < dt * (1.0 - 1e-12):
        c = 2.0 * a * d + 3.0 * b * d * d  # dF/dd
        with np.errstate(divide="ignore"):
            lim = np.where(c > 0, params.cfl * w / np.where(c > 0, c, 1.0), np.inf)
        h = min(dt - t, float(lim.min()))
        if n_sub >= params.max_explicit or (n_sub == 0 and dt / max(h, 1e-300) > params.max_explicit):
            return _implicit(np.array(delta, float), surf, a, b, dt)
        F = a * d * d + b * d ** 3
        div = F.copy()
        div[1:] -= F[:-1]
        d = d - h * div / w
        np.maximum(d, 0.0, out=d)  # rounding only: the CFL limit keeps every node non-negative
        shed += h * F[-1]
        t += h
        n_sub += 1
    return d, shed, n_sub


def _implicit(d0, surf: Surface, a, b, dt):
    """Linearised implicit upwind sweep: F_i = k_i d_i(new) with k_i = a_i d_i(old) + b_i d_i(old)^2. Forward sweep from the
    stagnation point; positive and conservative for any dt."""
    k = a * d0 + b * d0 * d0
    d = np.empty_like(d0)
    inflow = 0.0
    for i in range(len(d0)):
        d[i] = (surf.w[i] * d0[i] / dt + inflow) / (surf.w[i] / dt + k[i])
        inflow = k[i] * d[i]
    return d, dt * k[-1] * d[-1], -1


@lru_cache(maxsize=256)
def _stag(mach, altitude):
    from . import stagnation as st
    return st.stagnation_state(mach, altitude)


def edge_state(cond, flow):
    """(rho_e, p_inf) of the current flight state: from the condition's own stagnation state when it has one."""
    s = getattr(cond, "state", None)
    if s is None or abs(s.mach - flow.mach) > 1e-9 or abs(s.altitude - flow.altitude) > 1e-6:
        if flow.mach <= 0:
            raise ValueError("glass flow needs a flight state (mach > 0)")
        s = _stag(round(float(flow.mach), 6), round(float(flow.altitude), 3))
    return s.rho_e, s.p_inf


class FlowLedger:
    """Glass volume bookkeeping of the outer film over the whole surface (per radian for an axisymmetric body): cumulative
    outer-glass source, evaporated and shed volume; the audit compares them with the film left on the surface."""

    def __init__(self):
        self.source = 0.0
        self.evap = 0.0
        self.shed = 0.0
        self.n_explicit = 0
        self.n_implicit = 0
        self.max_sub = 0

    def audit(self, film):
        """Relative closure error of film + evaporated + shed against the cumulative source."""
        return abs(film + self.evap + self.shed - self.source) / max(self.source, 1e-300)
