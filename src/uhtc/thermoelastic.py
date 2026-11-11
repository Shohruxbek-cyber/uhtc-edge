"""Axisymmetric, quasi-static, linear, uncoupled thermoelasticity on the triangle meshes of the 2D model (version 3, Task 2).

The temperature field of a 2D run drives the stress; nothing feeds back to the thermal solve. Linear (constant-strain)
axisymmetric triangles: strains (eps_r, eps_z, eps_theta, gamma_rz) with the hoop strain u_r/r at the element centroid,
isotropic D(E, nu), thermal strain alpha(T) (T - T_ref) in r, z and theta (stress-free at T_ref = 300 K, the cold start of
every run), E and alpha evaluated at the element temperature (mean of its three nodes). Integrals per unit radian (2 pi
dropped), as in uhtc.conduction2d. Mesh convention of uhtc.conduction2d: nodes[:, 0] = x (axial), nodes[:, 1] = y (radius).

Boundary conditions (D-MECH-1): u_r = 0 on the axis always; "free": u_z = 0 at one node only (removes the rigid axial
translation, does not change the stress), the free-standing tip; "restrained": u_z = 0 on the back plane, the mounted tip.
Verified by V11 (validation/V11_plan.md).
"""
from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


@dataclass
class Props:
    """Temperature-dependent properties: E(T) Pa and alpha(T) 1/K as callables (vectorised) or constants; nu constant."""
    E: object
    alpha: object
    nu: float

    def e(self, T):
        return self.E(T) if callable(self.E) else np.full_like(T, float(self.E))

    def a(self, T):
        return self.alpha(T) if callable(self.alpha) else np.full_like(T, float(self.alpha))


def _geometry(nodes, tri):
    r = nodes[tri][:, :, 1]
    z = nodes[tri][:, :, 0]
    M = np.stack([np.ones_like(r), r, z], axis=2)  # rows: nodes, columns: 1, r, z
    Minv = np.linalg.inv(M)
    dNdr, dNdz = Minv[:, 1, :], Minv[:, 2, :]
    area = 0.5 * np.abs(np.linalg.det(M))
    rc = r.mean(1)
    return dNdr, dNdz, area, rc


def _B(dNdr, dNdz, rc):
    m = len(rc)
    B = np.zeros((m, 4, 6))
    B[:, 0, 0::2] = dNdr
    B[:, 1, 1::2] = dNdz
    B[:, 2, 0::2] = (1.0 / 3.0) / rc[:, None]
    B[:, 3, 0::2] = dNdz
    B[:, 3, 1::2] = dNdr
    return B


def _D(E, nu):
    lam = E / ((1 + nu) * (1 - 2 * nu))
    D = np.zeros((len(E), 4, 4))
    for i in range(3):
        for j in range(3):
            D[:, i, j] = lam * (1 - nu if i == j else nu)
    D[:, 3, 3] = lam * (1 - 2 * nu) / 2
    return D


def solve(nodes, tri, T_nodes, props: Props, T_ref=300.0, bc="free", back_nodes=None, axis_tol=None, z_fix_node=None):
    """Displacements (N, 2) [u_x axial, u_y radial] and element stresses (M, 4) [sigma_r, sigma_z, sigma_theta, tau_rz], Pa,
    with the element centroids (M, 2) [x, y]."""
    nodes = np.asarray(nodes, float)
    tri = np.asarray(tri)
    n = len(nodes)
    dNdr, dNdz, area, rc = _geometry(nodes, tri)
    B = _B(dNdr, dNdz, rc)
    Te = np.asarray(T_nodes, float)[tri].mean(1)
    E, al = props.e(Te), props.a(Te)
    D = _D(E, props.nu)
    w = area * rc
    eps0 = (al * (Te - T_ref))[:, None] * np.array([1.0, 1.0, 1.0, 0.0])
    Ke = np.einsum("mki,mkl,mlj->mij", B, D, B) * w[:, None, None]
    fe = np.einsum("mki,mkl,ml->mi", B, D, eps0) * w[:, None]
    dof = np.empty((len(tri), 6), int)  # dof 2k = u_r of node k, 2k+1 = u_z
    dof[:, 0::2] = 2 * tri
    dof[:, 1::2] = 2 * tri + 1
    rows = np.repeat(dof, 6, axis=1).ravel()
    cols = np.tile(dof, (1, 6)).ravel()
    K = coo_matrix((Ke.ravel(), (rows, cols)), shape=(2 * n, 2 * n)).tocsr()
    f = np.zeros(2 * n)
    np.add.at(f, dof.ravel(), fe.ravel())
    tol = axis_tol if axis_tol is not None else 1e-9 * max(np.ptp(nodes[:, 1]), 1e-30)
    fixed = list(2 * np.where(np.abs(nodes[:, 1]) <= tol)[0])  # u_r = 0 on the axis
    if bc == "free":
        k = int(z_fix_node) if z_fix_node is not None else int(np.argmin(nodes[:, 0]))
        fixed.append(2 * k + 1)
    elif bc == "restrained":
        if back_nodes is None or not len(back_nodes):
            raise ValueError("restrained: back_nodes needed")
        fixed += list(2 * np.asarray(back_nodes) + 1)
    else:
        raise ValueError(bc)
    fixed = np.unique(fixed)
    free = np.setdiff1d(np.arange(2 * n), fixed)
    u = np.zeros(2 * n)
    u[free] = spsolve(K[free][:, free], f[free])
    ue = u[dof]
    strain = np.einsum("mij,mj->mi", B, ue)
    sig = np.einsum("mij,mj->mi", D, strain - eps0)
    U = np.c_[u[1::2], u[0::2]]  # [u_x (axial), u_y (radial)]
    cent = nodes[tri].mean(1)
    return U, sig, cent


# 6-point Gauss rule for triangles (degree 4): area coordinates and weights (sum 1)
_GA, _GB, _GW1 = 0.816847572980459, 0.091576213509771, 0.109951743655322
_GC, _GD, _GW2 = 0.108103018168070, 0.445948490915965, 0.223381589678011
GAUSS6 = np.array([[_GA, _GB, _GB], [_GB, _GA, _GB], [_GB, _GB, _GA], [_GC, _GD, _GD], [_GD, _GC, _GD], [_GD, _GD, _GC]])
GW6 = np.array([_GW1] * 3 + [_GW2] * 3)


def quadratic_mesh(nodes, tri):
    """Straight-sided 6-node triangles from 3-node ones: mid-side nodes added once per edge. Returns (nodes6, tri6) with
    tri6 columns (v0, v1, v2, m01, m12, m20)."""
    nodes = np.asarray(nodes, float)
    edges = {}
    extra = []
    tri6 = np.empty((len(tri), 6), int)
    tri6[:, :3] = tri
    for k, (a, b, c) in enumerate(np.asarray(tri)):
        for col, (i, j) in zip((3, 4, 5), ((a, b), (b, c), (c, a))):
            key = (min(i, j), max(i, j))
            if key not in edges:
                edges[key] = len(nodes) + len(extra)
                extra.append(0.5 * (nodes[i] + nodes[j]))
            tri6[k, col] = edges[key]
    return np.vstack([nodes, np.array(extra)]), tri6


def solve_quadratic(nodes, tri, T_nodes, props: Props, T_ref=300.0, bc="free", back_nodes=None, z_fix_node=None):
    """As solve, with 6-node straight-sided axisymmetric triangles (the rule of the V11 plan when the linear solve misses the
    5 % mesh check). Mid-side temperatures are interpolated linearly (the thermal solution is linear on the 3-node mesh).
    Returns displacements of the 6-node mesh, stresses (M*6, 4) at the Gauss points, their positions (M*6, 2) [x, y], and
    a dict with the volume (per radian) and the temperature of every Gauss point and the 6-node mesh."""
    nodes = np.asarray(nodes, float)
    tri = np.asarray(tri)
    T_nodes = np.asarray(T_nodes, float)
    n3 = len(nodes)
    nodes6, tri6 = quadratic_mesh(nodes, tri)
    T6 = np.r_[T_nodes, np.zeros(len(nodes6) - n3)]
    for col, (i, j) in zip((3, 4, 5), ((0, 1), (1, 2), (2, 0))):
        T6[tri6[:, col]] = 0.5 * (T_nodes[tri[:, i]] + T_nodes[tri[:, j]])
    n = len(nodes6)
    dNdr3, dNdz3, area, _ = _geometry(nodes, tri)  # dL_i/dr, dL_i/dz (rows: elements, cols: vertices)
    rv, zv = nodes[tri][:, :, 1], nodes[tri][:, :, 0]
    m = len(tri)
    K_rows, K_cols, K_vals = [], [], []
    f = np.zeros(2 * n)
    dof = np.empty((m, 12), int)
    dof[:, 0::2] = 2 * tri6
    dof[:, 1::2] = 2 * tri6 + 1
    gp_sig_parts, gp_info = [], []
    Bs, Ds, eps0s, wts = [], [], [], []
    for g, wgt in zip(GAUSS6, GW6):
        L1, L2, L3 = g
        N = np.array([L1 * (2 * L1 - 1), L2 * (2 * L2 - 1), L3 * (2 * L3 - 1), 4 * L1 * L2, 4 * L2 * L3, 4 * L3 * L1])
        dNdL = np.array([[4 * L1 - 1, 0, 0], [0, 4 * L2 - 1, 0], [0, 0, 4 * L3 - 1], [4 * L2, 4 * L1, 0], [0, 4 * L3, 4 * L2],
                         [4 * L3, 0, 4 * L1]])  # (6 shape, 3 area coords)
        dNdr = dNdr3 @ dNdL.T  # (m, 6)
        dNdz = dNdz3 @ dNdL.T
        r = rv @ g
        z = zv @ g
        Tg = T6[tri6] @ N
        B = np.zeros((m, 4, 12))
        B[:, 0, 0::2] = dNdr
        B[:, 1, 1::2] = dNdz
        B[:, 2, 0::2] = N[None, :] / r[:, None]
        B[:, 3, 0::2] = dNdz
        B[:, 3, 1::2] = dNdr
        E, al = props.e(Tg), props.a(Tg)
        D = _D(E, props.nu)
        w = wgt * area * r
        e0 = (al * (Tg - T_ref))[:, None] * np.array([1.0, 1.0, 1.0, 0.0])
        Ke = np.einsum("mki,mkl,mlj->mij", B, D, B) * w[:, None, None]
        fe = np.einsum("mki,mkl,ml->mi", B, D, e0) * w[:, None]
        K_rows.append(np.repeat(dof, 12, axis=1).ravel())
        K_cols.append(np.tile(dof, (1, 12)).ravel())
        K_vals.append(Ke.ravel())
        np.add.at(f, dof.ravel(), fe.ravel())
        Bs.append(B), Ds.append(D), eps0s.append(e0), wts.append(w)
        gp_info.append((np.c_[z, r], w, Tg))
    K = coo_matrix((np.concatenate(K_vals), (np.concatenate(K_rows), np.concatenate(K_cols))), shape=(2 * n, 2 * n)).tocsr()
    tol = 1e-9 * max(np.ptp(nodes6[:, 1]), 1e-30)
    fixed = list(2 * np.where(np.abs(nodes6[:, 1]) <= tol)[0])
    if bc == "free":
        k = int(z_fix_node) if z_fix_node is not None else int(np.argmin(nodes6[:, 0]))
        fixed.append(2 * k + 1)
    elif bc == "restrained":
        bn = set(int(b) for b in back_nodes)
        on_back = [tri6[e, c] for e in range(m) for c, (i, j) in zip((3, 4, 5), ((0, 1), (1, 2), (2, 0)))
                   if int(tri6[e, i]) in bn and int(tri6[e, j]) in bn]
        fixed += list(2 * np.asarray(sorted(bn | set(int(v) for v in on_back))) + 1)
    else:
        raise ValueError(bc)
    fixed = np.unique(fixed)
    free = np.setdiff1d(np.arange(2 * n), fixed)
    u = np.zeros(2 * n)
    u[free] = spsolve(K[free][:, free], f[free])
    ue = u[dof]
    sig = np.concatenate([np.einsum("mij,mj->mi", D, np.einsum("mij,mj->mi", B, ue) - e0) for B, D, e0 in zip(Bs, Ds, eps0s)])
    pos = np.concatenate([p for p, _, _ in gp_info])
    info = dict(vol=np.concatenate([w for _, w, _ in gp_info]), T=np.concatenate([t for _, _, t in gp_info]), nodes6=nodes6, tri6=tri6)
    U = np.c_[u[1::2], u[0::2]]
    return U, sig, pos, info


def principal_max(sig):
    """Largest principal stress of each element: the in-plane eigenvalues of (sigma_r, sigma_z, tau_rz) and the hoop stress."""
    sr, sz, st, trz = sig.T
    c, rad = 0.5 * (sr + sz), np.sqrt((0.5 * (sr - sz)) ** 2 + trz ** 2)
    return np.maximum(c + rad, st)


def back_plane_nodes(sol, rel_tol=1e-6):
    """Nodes of the back plane of a uhtc.conduction2d nose mesh: the line from the last surface node along the inward normal of
    the flank to the axis."""
    nodes = sol.nodes
    p0 = nodes[sol.surf[-1]]
    th = sol.g.theta_c
    pe = np.array([p0[0] + p0[1] * np.tan(th), 0.0])  # end point on the axis (build_mesh)
    d = pe - p0
    L = np.hypot(*d)
    rel = nodes - p0
    dist = np.abs(rel[:, 0] * d[1] - rel[:, 1] * d[0]) / L
    t = (rel @ d) / L ** 2
    return np.where((dist <= rel_tol * L) & (t >= -1e-9) & (t <= 1 + 1e-9))[0]
