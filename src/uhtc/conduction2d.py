"""Quasi-2D conduction in the nose region: linear triangles (P1 finite elements), axisymmetric
(sphere-cone nose tip) or planar (cylinder-wedge wing leading edge).

Geometry (half body, y >= 0): nose arc of radius R_n centred on (R_n, 0), from the stagnation point
(body angle phi = 0) to the tangency point phi_t = pi/2 - theta_c, then a straight flank at half-angle
theta_c over a wetted length L_f.  Every surface point is joined to the axis along its inward normal
(arc normals meet at the centre, flank normals meet the axis), which gives a structured (s, eta) grid
whose cells are split into triangles; the back boundary is the inward normal from the end of the flank
(a plane perpendicular to the flank).  Boundary conditions: surface = shared surface balance
(uhtc.surface, with an optional oxide-scale resistance per node); axis = symmetry; back plane =
adiabatic (DECISIONS.md D-2D-2).

Heating distribution (DECISIONS.md D-2D-1): q = Phi(s) (H_s - h_w(T)), Phi(s) = Phi_s g(s);
on the arc g = cos(phi) (Tauber 1989 p. 8, Figs. 5-6); on the flank g = q_flank(x = s)/q_s,cold with the
laminar sharp-cone (axisymmetric, Tauber eq. 46) or wedge (planar, Tauber eq. 48) correlation, running
length taken as the wetted distance s from the stagnation point.

Time integration: backward Euler on nodal enthalpy (lumped, r-weighted capacity), Newton iteration with
the conductivity lagged by one iterate (Picard on k); sparse direct solve.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import brentq

from . import heating as ht


@dataclass
class NoseGeometry:
    R_n: float  # m
    theta_c: float  # rad, cone / wedge half-angle
    L_flank: float  # m, wetted length of the flank
    axisymmetric: bool = True


def _eta_grid(n_eta, eta1):
    if eta1 * n_eta >= 1.0:
        return np.linspace(0.0, 1.0, n_eta + 1)
    r = brentq(lambda q: eta1 * (q ** n_eta - 1) / (q - 1) - 1.0, 1.0 + 1e-12, 3.0)
    eta = np.r_[0.0, np.cumsum(eta1 * r ** np.arange(n_eta))]
    eta[-1] = 1.0
    return eta


def _assemble(xs, ys, xe, ye, eta, scale):
    """Structured (surface index i, depth index j) grid between surface points and inward end points ->
    merged nodes, counter-clockwise triangles (degenerate ones dropped), surface node indices."""
    ns, ne = len(xs), len(eta)
    X = xs[:, None] + eta[None, :] * (xe - xs)[:, None]
    Y = ys[:, None] + eta[None, :] * (ye - ys)[:, None]
    pts = np.c_[X.ravel(), Y.ravel()]
    idx = np.arange(ns * ne).reshape(ns, ne)
    key = np.round(pts / (scale * 1e-9)).astype(np.int64)
    _, uniq, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inv = inv.ravel()
    nodes = pts[uniq]
    a, b = inv[idx[:-1, :-1]].ravel(), inv[idx[1:, :-1]].ravel()
    c, d = inv[idx[1:, 1:]].ravel(), inv[idx[:-1, 1:]].ravel()
    tri = np.r_[np.c_[a, b, c], np.c_[a, c, d]]
    P = nodes[tri]
    area2 = (P[:, 1, 0] - P[:, 0, 0]) * (P[:, 2, 1] - P[:, 0, 1]) - (P[:, 2, 0] - P[:, 0, 0]) * (P[:, 1, 1] - P[:, 0, 1])
    keep = np.abs(area2) > 1e-14 * scale ** 2
    tri, area2 = tri[keep], area2[keep]
    flip = area2 < 0
    tri[flip] = tri[flip][:, [0, 2, 1]]
    return nodes, tri, inv[idx[:, 0]]


def build_mesh(g: NoseGeometry, n_arc=60, n_flank=60, n_eta=40, eta1=0.004, flank_growth=1.04):
    """Return nodes (N,2) [x, y], triangles (M,3), surface node indices (ordered from the stagnation
    point), their wetted coordinate s and body angle phi, and the number of arc cells."""
    phi_t = np.pi / 2 - g.theta_c
    u = np.linspace(0, 1, n_arc + 1)
    phi = phi_t * (u ** 1.15)  # mild clustering towards the stagnation point
    s_arc = g.R_n * phi
    ds0 = s_arc[-1] - s_arc[-2]
    seg = np.array([ds0 * flank_growth ** i for i in range(n_flank)])
    seg = seg * g.L_flank / np.sum(seg)
    s_fl = np.r_[0.0, np.cumsum(seg)]
    xt, yt = g.R_n * (1 - np.cos(phi_t)), g.R_n * np.sin(phi_t)
    xs = np.r_[g.R_n * (1 - np.cos(phi)), xt + s_fl[1:] * np.cos(g.theta_c)]
    ys = np.r_[g.R_n * np.sin(phi), yt + s_fl[1:] * np.sin(g.theta_c)]
    s = np.r_[s_arc, s_arc[-1] + s_fl[1:]]
    phis = np.r_[phi, np.full(n_flank, phi_t)]
    # inward end points: centre for the arc, axis intersection along the normal for the flank
    xe = np.where(np.arange(len(xs)) <= n_arc, g.R_n, xs + ys * np.tan(g.theta_c))
    ye = np.zeros_like(xs)
    nodes, tri, surf = _assemble(xs, ys, xe, ye, _eta_grid(n_eta, eta1), g.R_n)
    return nodes, tri, surf, s, phis, n_arc


def build_sphere_mesh(R, n_arc=80, n_eta=40, eta1=0.01):
    """Full sphere (axisymmetric) or full disc (planar) of radius R: meridian phi in [0, pi], all normals
    to the centre.  Verification only (exact solutions exist for a uniform surface flux)."""
    phi = np.linspace(0.0, np.pi, n_arc + 1)
    xs, ys = R * (1 - np.cos(phi)), R * np.sin(phi)
    ys[-1] = 0.0
    nodes, tri, surf = _assemble(xs, ys, np.full_like(xs, R), np.zeros_like(xs), _eta_grid(n_eta, eta1), R)
    return nodes, tri, surf, R * phi, phi, n_arc


class Solver2D:
    def __init__(self, geom: NoseGeometry, material, T_init, mesh_kw=None, mesh=None):
        """mesh: optional prebuilt (nodes, tri, surf, s, phi, n_arc), e.g. build_sphere_mesh (verification)."""
        self.g = geom
        self.mat = material
        nodes, tri, surf, s, phi, n_arc = mesh if mesh is not None else build_mesh(geom, **(mesh_kw or {}))
        self.nodes, self.tri, self.surf, self.s, self.phi = nodes, tri, surf, s, phi
        self.on_flank = np.arange(len(s)) > n_arc
        n = len(nodes)
        P = nodes[tri]
        x, y = P[:, :, 0], P[:, :, 1]
        A = 0.5 * ((x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0]))
        # shape-function gradients
        b = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], 1) / (2 * A[:, None])
        c = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], 1) / (2 * A[:, None])
        if geom.axisymmetric:
            rc = y.mean(1)
            w_el = A * rc  # int r dA (2 pi dropped everywhere)
            lump = A[:, None] / 12 * (2 * y + y.sum(1, keepdims=True) - y)  # int N_i r dA
        else:
            w_el = A
            lump = np.repeat(A[:, None] / 3, 3, 1)
        self.A, self.w_el = A, w_el
        self.Ke_geo = (b[:, :, None] * b[:, None, :] + c[:, :, None] * c[:, None, :]) * w_el[:, None, None]
        self.vol = np.zeros(n)
        np.add.at(self.vol, tri, lump)
        rows = np.repeat(tri, 3, axis=1).ravel()
        cols = np.tile(tri, (1, 3)).ravel()
        self._rows, self._cols = rows, cols
        # surface edge weights (int N_i r ds lumped to nodes)
        sn = nodes[surf]
        L = np.hypot(np.diff(sn[:, 0]), np.diff(sn[:, 1]))
        if geom.axisymmetric:
            r = sn[:, 1]
            wa = L * (2 * r[:-1] + r[1:]) / 6
            wb = L * (r[:-1] + 2 * r[1:]) / 6
        else:
            wa = wb = L / 2
        self.w_surf = np.zeros(len(surf))
        self.w_surf[:-1] += wa
        self.w_surf[1:] += wb
        T0 = float(T_init)
        self.m = self.vol * float(material.rho(np.array(T0)))
        self.H = np.full(n, float(material.h_of_T(T0)))
        self.t = 0.0
        self.energy_in = 0.0
        self.E0 = float(np.sum(self.m * self.H))

    @property
    def T(self):
        return self.mat.T_of_h(self.H)

    def heating_shape(self, V, rho_inf, q_s_cold, dist="tauber", **kw):
        """g(s) = q(s)/q_s for the named heating distribution (HEAT_DISTRIBUTIONS)."""
        if dist not in HEAT_DISTRIBUTIONS:
            raise ValueError(f"unknown heating distribution {dist!r}; known: {sorted(HEAT_DISTRIBUTIONS)}")
        return HEAT_DISTRIBUTIONS[dist](self, V=V, rho_inf=rho_inf, q_s_cold=q_s_cold, **kw)

    def _shape_tauber(self, V, rho_inf, q_s_cold, **kw):
        """Draft baseline: cos(phi) on the arc; on the flank the cold-wall laminar cone (axisymmetric) or wedge
        (planar) correlation at running length s, divided by the cold-wall stagnation value q_s_cold
        computed with the same correlation family (Tauber 1989 eq. 40 sphere / eq. 41 cylinder)."""
        gs = ht.cosine_distribution(self.phi)
        if np.any(self.on_flank):
            fn = ht.tauber_cone if self.g.axisymmetric else ht.tauber_wedge
            q_fl = fn(rho_inf, V, self.s[self.on_flank], self.g.theta_c, 0.0, 1.0)  # cold wall
            gs = gs.copy()
            gs[self.on_flank] = q_fl / q_s_cold
        return gs

    def step(self, dt, surf_flux, max_iter=60, tol=1e-4):
        """surf_flux(T_surface_nodes) -> (q, dq/dT) arrays (flux into the solid, per unit area)."""
        H0 = self.H.copy()
        H = H0.copy()
        n = len(H)
        for it in range(max_iter):
            try:
                T = self.mat.T_of_h(H)
            except ValueError:
                return False, np.inf
            dTdH = self.mat.dTdh(H)
            kel = self.mat.k(T[self.tri].mean(1))
            Ke = self.Ke_geo * kel[:, None, None]
            K = sp.csr_matrix((Ke.ravel(), (self._rows, self._cols)), shape=(n, n))
            q, dq = surf_flux(T[self.surf])
            F = np.zeros(n)
            np.add.at(F, self.surf, q * self.w_surf)
            R = self.m * (H - H0) / dt + K @ T - F
            Jd = np.zeros(n)
            np.add.at(Jd, self.surf, -dq * self.w_surf)
            J = (K + sp.diags(Jd)) @ sp.diags(dTdH) + sp.diags(self.m / dt)
            dH = spla.spsolve(J.tocsc(), -R)
            dT = np.abs(dH * dTdH)
            if dT.max() > 200:
                dH *= 200 / dT.max()
            H = H + dH
            if np.abs(dH * dTdH).max() < tol and it > 0:
                break
        else:
            return False, np.inf
        Tn = self.mat.T_of_h(H)
        q, _ = surf_flux(Tn[self.surf])
        self.energy_in += float(np.sum(q * self.w_surf)) * dt
        dTmax = float(np.abs(Tn - self.mat.T_of_h(H0)).max())
        self.H = H
        self.t += dt
        return True, dTmax

    def energy_residual(self):
        d = float(np.sum(self.m * self.H)) - self.E0
        return (self.energy_in - d) / max(abs(self.energy_in), 1e-300)


# ------------------------------------------------------------------------------------ heating distributions (D-2D-7)
GAMMA_SL = 1.4  # isentropic exponent of the surface flow (perfect gas, as in the Crabtree et al. 1965 charts)
_DIG = Path(__file__).resolve().parents[2] / "literature" / "digitized"


def lees_distribution(g: NoseGeometry, s, mach, V, H_s, p_s, p_inf, n_fine=6000):
    """Laminar local-similarity heating distribution q(s)/q_s over the arc and the flank in one integral
    (Lees 1956 / Kemp, Rose & Detra 1959 as reproduced by Cleary 1969 NASA TN D-5450 eqs. 2, 4, 5 and DeJarnette &
    Tai 1969 NASA CR-101678 eqs. 10-11; Crabtree et al. 1965 ARC R&M 3637 eqs. 3.4-3.5 for the pressure-gradient factor
    a = 0.332 (1 + 0.1853 beta^1/2)), with modified-Newtonian pressure p/p_s = (1 - eps) cos^2(psi) + eps,
    eps = p_inf/p_s, psi = angle between the surface normal and the stream (phi on the arc, 90 deg - theta_c on the
    flank), and isentropic surface flow from the stagnation point (u_e/V = [(2 H_s/V^2)(1 - (p/p_s)^((g-1)/g))]^1/2):

        q/q_s = (p/p_s)(u_e/V)(r/R)^j F(beta) / ( [2 (1 + j)]^1/2 G0 I^1/2 F(beta_0) ),
        I = int_0^x (p/p_s)(u_e/V)(r/R)^2j dx,  beta = 2 (d(u_e/V)/dx) I / ((p/p_s)(u_e/V)^2 (r/R)^2j),
        G0 = [d(u_e/V)/dx]_0^1/2,  F = 1 + 0.1853 beta^1/2,  beta_0 = 1/2 (j = 1) or 1 (j = 0),  x = s/R,

    j = 1 axisymmetric (nose tip), 0 planar (wing edge).  The normalisation gives q/q_s -> 1 at the stagnation point."""
    j = 1 if g.axisymmetric else 0
    phi_t = np.pi / 2 - g.theta_c
    x_t = phi_t
    x_max = float(np.max(s)) / g.R_n * 1.0001 + 1e-9
    x = x_max * np.linspace(0.0, 1.0, n_fine) ** 2  # fine near the stagnation point (integrand ~ x^3 there)
    on_arc = x <= x_t
    psi = np.where(on_arc, x, phi_t)
    r = np.where(on_arc, np.sin(x), np.sin(phi_t) + (x - x_t) * np.sin(g.theta_c))
    eps = p_inf / p_s
    pr = (1.0 - eps) * np.cos(psi) ** 2 + eps
    e = (GAMMA_SL - 1.0) / GAMMA_SL
    k = 2.0 * H_s / V ** 2
    u = np.sqrt(np.maximum(k * (1.0 - pr ** e), 0.0))
    D0 = np.sqrt(k * e * (1.0 - eps))  # d(u_e/V)/dx at the stagnation point (small-x limit of the expression for u)
    G0 = np.sqrt(D0)
    w = pr * u * r ** (2 * j)
    I = np.r_[0.0, np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(x))]
    du = np.gradient(u, x)
    with np.errstate(divide="ignore", invalid="ignore"):
        q = pr * u * r ** j / (np.sqrt(2.0 * (1 + j)) * G0 * np.sqrt(I))
        beta = np.where(I > 0, 2.0 * du * I / (pr * u ** 2 * r ** (2 * j)), 1.0 / (1 + j))
    beta = np.clip(np.nan_to_num(beta, nan=0.0), 0.0, None)
    F = (1.0 + 0.1853 * np.sqrt(beta)) / (1.0 + 0.1853 * np.sqrt(1.0 / (1 + j)))
    q = q * F
    q[0] = 1.0
    return np.interp(np.asarray(s) / g.R_n, x, q)


def _kemp_shoulder(phi, axisymmetric):
    """Upper-bound shoulder values between 70 deg and the tangency: Kemp-Rose-Detra / Lees curves of Tauber 1989
    Fig. 5 (hemisphere-cylinder) or Fig. 6 (cylinder), DIGITIZED (literature/digitized/Tauber1989_fig5/6.csv),
    joined linearly to cos(70 deg) at 70 deg."""
    f = _DIG / ("Tauber1989_fig5.csv" if axisymmetric else "Tauber1989_fig6.csv")
    import csv
    rows = [r for r in csv.DictReader(l for l in open(f, encoding="utf-8") if not l.startswith("#"))]
    key = [k for k in rows[0] if k.startswith("q_over_qs")][0]
    a = np.r_[70.0, [float(r["body_angle_deg"]) for r in rows]]
    v = np.r_[np.cos(np.radians(70.0)), [float(r[key]) for r in rows]]
    return np.interp(np.degrees(phi), a, v)


PHI_COS = np.radians(70.0)  # cos(phi) is supported to about 70 deg (Tauber 1989 p. 8)


def baseline_distribution(g: NoseGeometry, s, phi, on_flank, mach, V, H_s, p_s, p_inf):
    """Baseline heating distribution (D-2D-7): cos(phi) on the arc up to 70 deg (Tauber 1989 p. 8, Figs. 5-6);
    from the tangency on, the Lees / Kemp-Rose-Detra local-similarity result (lees_distribution; reproduces the
    Crabtree et al. 1965 Fig. 12 hemisphere-cone charts to within 10 %); between 70 deg and the tangency, where the
    sources disagree (0.25-2 x cos), a linear bridge in phi between the two."""
    lees = lees_distribution(g, s, mach, V, H_s, p_s, p_inf)
    phi_t = np.pi / 2 - g.theta_c
    lees_t = float(lees_distribution(g, np.array([g.R_n * phi_t]), mach, V, H_s, p_s, p_inf)[0])
    w = np.clip((phi - PHI_COS) / (phi_t - PHI_COS), 0.0, 1.0)
    arc = np.where(phi <= PHI_COS, np.cos(phi), (1.0 - w) * np.cos(PHI_COS) + w * lees_t)
    return np.where(on_flank, lees, arc)


def _shape_lees(self, V, rho_inf, q_s_cold, flow=None, factor=None, **kw):
    from . import atmosphere as atm
    p_inf = float(atm.state(flow.altitude)[1])
    gs = baseline_distribution(self.g, self.s, self.phi, self.on_flank, flow.mach, V, flow.H_s, flow.p_s, p_inf)
    if factor is not None:  # band cases (D-2D-7): scale the distribution beyond 70 deg (arc and flank)
        lo, hi = factor
        beyond = self.phi > np.radians(70.0)
        if lo is not None:
            gs = np.where(beyond, lo * gs, gs)
        if hi is not None:
            kemp = np.where(self.on_flank, 0.0, _kemp_shoulder(self.phi, self.g.axisymmetric))
            gs = np.where(beyond, np.maximum(hi * gs, kemp), gs)
    return gs


def integral_distribution(g: NoseGeometry, s, mach, V, H_s, p_s, p_inf, n_fine=400, T_wall=300.0, perfect_gas=False):
    """Version 3, Task 3 (flank model "integral"; D-HEAT-7): laminar engineering integral method of Zoby et al. as printed by
    Riley 1992 (NASA TM-107630, Eqs. 3.3.1-3.3.2): q = 0.22 Re_theta^-1 (rho*/rho_e)(mu*/mu_e) rho_e u_e (h_aw - h_w) Pr_w^-0.6
    with theta = 0.664 (int rho* mu* u_e r^2j ds)^1/2 / (rho_e u_e r^j), which reduces to
        q = (0.22/0.664) rho* mu* u_e r^j (h_aw - h_w) Pr_w^-0.6 / (int_0^s rho* mu* u_e r^2j ds)^1/2,
    marched along s. Reference state (Eckert, as Crabtree et al. 1965 Eqs. 1.3-1.4): h* = h_e + 0.5 (h_w - h_e) + 0.22 (h_aw - h_e),
    h_aw = h_e + Pr^1/2 u_e^2/2, rho* and mu* of equilibrium air at (p, h*) (Sutherland viscosity). Edge state option (a) of the
    prompt: modified-Newtonian pressure and isentropic expansion from the stagnation point, as lees_distribution. Cold wall at
    T_wall; the stagnation limit (u_e = beta s, r = s^j) normalises q/q_s -> 1. perfect_gas (wind-tunnel anchor of V12 only):
    ideal air with h = cp T on an absolute basis (H_s = cp T_0) and the Sutherland viscosity, instead of equilibrium air."""
    from . import atmosphere as atm, thermo as th
    j = 1 if g.axisymmetric else 0
    phi_t = np.pi / 2 - g.theta_c
    x_max = float(np.max(s)) / g.R_n * 1.0001 + 1e-9
    x = x_max * np.linspace(0.0, 1.0, n_fine) ** 2
    on_arc = x <= phi_t
    psi = np.where(on_arc, x, phi_t)
    r = np.where(on_arc, np.sin(x), np.sin(phi_t) + (x - phi_t) * np.sin(g.theta_c))
    eps = p_inf / p_s
    pr = (1.0 - eps) * np.cos(psi) ** 2 + eps
    e = (GAMMA_SL - 1.0) / GAMMA_SL
    k = 2.0 * H_s / V ** 2
    u = V * np.sqrt(np.maximum(k * (1.0 - pr ** e), 0.0))  # u_e, m/s
    p = p_s * pr
    CP_PG, R_AIR = 1004.5, 287.05
    h_w = CP_PG * T_wall if perfect_gas else float(th.h_air_frozen(T_wall))
    Pr = 0.71
    h_e = H_s - 0.5 * u * u
    h_aw = h_e + np.sqrt(Pr) * 0.5 * u * u
    h_s = h_e + 0.5 * (h_w - h_e) + 0.22 * (h_aw - h_e)
    rm = np.empty_like(x)  # rho* mu*
    for i in range(len(x)):
        if perfect_gas:
            Ts = float(h_s[i]) / CP_PG
            rm[i] = float(p[i]) / (R_AIR * Ts) * float(atm.viscosity(Ts))
            continue
        Ts = th.T_air_equilibrium(float(h_s[i]), float(p[i]))
        M, _, _ = th.mixture_props(th.equilibrium_air(Ts, float(p[i])), Ts)
        rm[i] = float(p[i]) * M / (th.R * Ts) * float(atm.viscosity(Ts))
    xs = x * g.R_n
    w = rm * u * (r * g.R_n) ** (2 * j)
    I = np.r_[0.0, np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(xs))]
    beta = V * np.sqrt(k * e * (1.0 - eps)) / g.R_n  # du_e/ds at the stagnation point
    q0 = np.sqrt(2.0 * (1 + j) * rm[0] * beta) * (H_s - h_w)
    with np.errstate(divide="ignore", invalid="ignore"):
        q = rm * u * (r * g.R_n) ** j * (h_aw - h_w) / np.sqrt(I) / q0
    q[0] = 1.0
    return np.interp(np.asarray(s) / g.R_n, x, q)


def _shape_integral(self, V, rho_inf, q_s_cold, flow=None, **kw):
    """Flank model "integral": as the baseline (cos(phi) to 70 deg, linear bridge to the tangency value) with the integral method
    in place of the Lees local similarity beyond 70 deg and on the flank."""
    from . import atmosphere as atm
    p_inf = float(atm.state(flow.altitude)[1])
    phi_t = np.pi / 2 - self.g.theta_c
    qq = integral_distribution(self.g, np.r_[self.s, self.g.R_n * phi_t], flow.mach, V, flow.H_s, flow.p_s, p_inf)
    q, q_t = qq[:-1], float(qq[-1])
    w = np.clip((self.phi - PHI_COS) / (phi_t - PHI_COS), 0.0, 1.0)
    arc = np.where(self.phi <= PHI_COS, np.cos(self.phi), (1.0 - w) * np.cos(PHI_COS) + w * q_t)
    return np.where(self.on_flank, q, arc)


HEAT_DISTRIBUTIONS = {
    "integral": lambda sol, **kw: _shape_integral(sol, **kw),  # version 3, Task 3 (D-HEAT-7)
    "tauber": Solver2D._shape_tauber,  # draft baseline, kept as a comparison case
    "lees": lambda sol, **kw: _shape_lees(sol, **kw),
    "lees_lo": lambda sol, **kw: _shape_lees(sol, factor=(0.66, None), **kw),
    "lees_hi": lambda sol, **kw: _shape_lees(sol, factor=(None, 1.5), **kw),
}
