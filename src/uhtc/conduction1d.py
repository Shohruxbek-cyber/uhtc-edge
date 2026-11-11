"""1D transient conduction along the stagnation line (enthalpy form, backward Euler, Newton).

Adapted from the author's fusion-wall solver (package ife_wall, module thermal)
(classes Layer, Mesh, EnergyLedger, ThermalSolver; functions _geometric_cells and step) and
src/ife_wall/materials.py (class Material: tabulated enthalpy h(T), T(h), dT/dh).  What was kept:
the Lagrangian mass coordinate (cell areal mass m_i, thickness m_i/rho(T_i)), the enthalpy unknown,
backward Euler with Newton iteration and a tridiagonal (banded) solve, the geometric mesh, the
harmonic-mean face conductances and the energy ledger checked against the change in stored enthalpy.
What changed: no volumetric sources, melting or evaporation of the solid (the ceramic does not melt
below the limits used here); the front boundary is a nonlinear surface balance q(T_0, t) supplied by
the caller with its derivative (used in the Newton Jacobian); geometry factors allow a planar slab
(default) or a spherical shell (cells at radius r, used only in tests).

The front-boundary callable receives the temperature of the first cell centre and returns
(q_into_solid [W/m^2], dq/dT_0 [W/m^2/K]).  The surface temperature and the thin oxide scale are
handled by the caller (uhtc.coupled), which folds the half-cell resistance dx_0/(2 k_0) into the
surface balance.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np

try:  # numba is optional; the pure-Python fallback gives identical results, only slower
    from numba import njit
except ImportError:  # pragma: no cover
    def njit(*a, **k):
        return (lambda f: f) if not a or not callable(a[0]) else a[0]

T_REF = 298.15  # K, enthalpy reference


@njit(cache=True)
def thomas(lower, diag, upper, rhs):
    """Tridiagonal solve (lower[i] multiplies x[i-1] in row i, upper[i] multiplies x[i+1])."""
    n = diag.shape[0]
    c = np.empty(n)
    d = np.empty(n)
    c[0] = upper[0] / diag[0]
    d[0] = rhs[0] / diag[0]
    for i in range(1, n):
        den = diag[i] - lower[i] * c[i - 1]
        c[i] = upper[i] / den if i < n - 1 else 0.0
        d[i] = (rhs[i] - lower[i] * d[i - 1]) / den
    x = np.empty(n)
    x[n - 1] = d[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = d[i] - c[i] * x[i + 1]
    return x


@dataclass
class Material:
    """Tabulated material (adapted from ife_wall.materials.Material, without melting)."""

    name: str
    rho: Callable  # kg/m^3 of T
    k: Callable  # W/(m K) of T
    cp: Callable  # J/(kg K) of T
    T_min: float = 200.0
    T_max: float = 4000.0
    dT_grid: float = 1.0
    sources: dict = field(default_factory=dict)

    def __post_init__(self):
        grid = np.arange(self.T_min, self.T_max + self.dT_grid, self.dT_grid)
        mid = 0.5 * (grid[1:] + grid[:-1])
        h = np.concatenate([[0.0], np.cumsum(self.cp(mid) * np.diff(grid))])
        h -= np.interp(T_REF, grid, h)
        if np.any(np.diff(h) <= 0):
            raise ValueError(f"{self.name}: enthalpy table not strictly increasing")
        self._T, self._h = grid, h
        self._dTdh = np.diff(grid) / np.diff(h)
        rho_tab, k_tab = np.asarray(self.rho(grid), float), np.asarray(self.k(grid), float)
        self.rho_fn, self.k_fn = self.rho, self.k
        self.rho = lambda T, _g=grid, _v=rho_tab: np.interp(T, _g, _v)
        self.k = lambda T, _g=grid, _v=k_tab: np.interp(T, _g, _v)

    def h_of_T(self, T):
        T = np.asarray(T, dtype=float)
        if np.any(T < self.T_min) or np.any(T > self.T_max):
            raise ValueError(f"{self.name}: T outside {self.T_min}-{self.T_max} K")
        return np.interp(T, self._T, self._h)

    def T_of_h(self, h):
        h = np.asarray(h, dtype=float)
        if np.any(h < self._h[0]) or np.any(h > self._h[-1]):
            raise ValueError(f"{self.name}: enthalpy outside table")
        return np.interp(h, self._h, self._T)

    def dTdh(self, h):
        i = np.clip(np.searchsorted(self._h, h, side="right") - 1, 0, len(self._dTdh) - 1)
        return self._dTdh[i]


def geometric_cells(L, dx0, r=1.1, dx_max=np.inf):
    """Cell widths filling length L with geometric growth (from ife_wall.thermal._geometric_cells)."""
    cells, x = [], 0.0
    while x < L - 1e-15 * L:
        d = min(dx0 * r ** len(cells), dx_max, L - x)
        cells.append(d)
        x += d
    if len(cells) > 1 and cells[-1] < 0.5 * cells[-2]:
        sliver = cells.pop()
        cells[-1] += sliver
    return np.array(cells)


@dataclass
class Ledger:
    front_in: float = 0.0  # J/m^2 (per unit front area) entering through the front boundary
    back_out: float = 0.0
    stored_initial: float = 0.0
    stored_final: float = 0.0

    def residual(self):
        d = self.stored_final - self.stored_initial
        res = self.front_in - self.back_out - d
        return res, res / max(abs(self.front_in), abs(d), 1e-300)


class Solver1D:
    """Planar (geometry='planar') or spherical-shell (geometry='sphere', outer radius R_out) conduction.

    Areas are normalised by the front area, so fluxes and the ledger are per unit front area."""

    def __init__(self, material: Material, thickness, T_init, dx0=2e-6, growth=1.08, dx_max=np.inf,
                 front: Optional[Callable[[float, float], tuple]] = None,
                 back: str = "adiabatic", h_back=0.0, T_back=300.0, geometry="planar"):
        self.mat = material
        dx = geometric_cells(thickness, dx0, growth, dx_max)
        self.geometry = geometry
        self.L = float(thickness)
        faces = np.r_[0.0, np.cumsum(dx)]  # depth of faces from the surface (initial temperature)
        if geometry == "planar":
            self.area_face = np.ones(len(faces))
            vol = dx
        elif geometry == "sphere":
            R = self.L
            r = R - faces
            self.area_face = (r / R) ** 2
            vol = (r[:-1] ** 3 - r[1:] ** 3) / (3 * R ** 2)
        else:
            raise ValueError(geometry)
        self.depth_c = 0.5 * (faces[:-1] + faces[1:])
        self.faces = faces
        T0 = float(T_init)
        self.vol = vol  # m (volume per unit front area) at T_init
        self.m = vol * float(material.rho(np.array(T0)))  # kg per unit front area (fixed)
        self.H = np.full(len(dx), float(material.h_of_T(T0)))
        self.front = front if front is not None else (lambda T0, t: (0.0, 0.0))
        self.back, self.h_back, self.T_back = back, h_back, T_back
        self.t = 0.0
        self.ledger = Ledger(stored_initial=self.stored())

    # --- helpers -------------------------------------------------------------------------
    def T_of(self, H):
        return self.mat.T_of_h(H)

    @property
    def T(self):
        return self.T_of(self.H)

    def stored(self):
        return float(np.sum(self.m * self.H))

    def _conductances(self, T):
        """Face conductances G_{i+1/2} (W/K per unit front area); geometry from the initial mesh
        (thermal expansion ~1 % neglected in the spacing; mass, hence heat capacity, is exact)."""
        k = self.mat.k(T)
        # half-cell resistances in series (harmonic combination), times the face area
        half = np.diff(self.faces) / 2
        r_half_lo = half[:-1] / k[:-1]
        r_half_hi = half[1:] / k[1:]
        G = self.area_face[1:-1] / (r_half_lo + r_half_hi)
        G_back = 0.0
        if self.back == "convective":
            G_back = self.area_face[-1] / (half[-1] / k[-1] + 1.0 / self.h_back)
        return G, G_back

    def half_cell_resistance(self):
        """dx_0 / (2 k_0), K m^2/W, for the caller's surface balance."""
        return float((self.faces[1] - self.faces[0]) / 2 / self.mat.k(self.T[:1])[0])

    # --- one step ---------------------------------------------------------------------------
    def step(self, dt, max_iter=50, tol_T=1e-6):
        """Backward-Euler step.  Returns (accepted, max |dT|).  Front flux and its derivative are
        re-evaluated at every Newton iterate (fully implicit surface balance)."""
        t1 = self.t + dt
        H0 = self.H.copy()
        H = H0.copy()
        n = len(H)
        polish = False
        for it in range(max_iter):
            try:
                T = self.T_of(H)
            except ValueError:
                return False, np.inf
            dTdH = self.mat.dTdh(H)
            G, G_back = self._conductances(T)
            q, dq = self.front(float(T[0]), t1)
            q_last, dq_last, T0_last = q, dq, float(T[0])
            flux = np.zeros(n)
            flux[:-1] -= G * (T[:-1] - T[1:])
            flux[1:] += G * (T[:-1] - T[1:])
            flux[-1] -= G_back * (T[-1] - self.T_back)
            flux[0] += q * self.area_face[0]
            R = self.m * (H - H0) / dt - flux
            diag = self.m / dt + np.r_[G, 0.0] * dTdH + np.r_[0.0, G] * dTdH
            diag[-1] += G_back * dTdH[-1]
            diag[0] -= dq * self.area_face[0] * dTdH[0]
            upper = np.r_[-G * dTdH[1:], 0.0]  # row i, column i+1
            lower = np.r_[0.0, -G * dTdH[:-1]]  # row i, column i-1
            dH = thomas(lower, diag, upper, -R)
            dT_est = np.abs(dH * dTdH)
            if dT_est.max() > 200.0:
                dH *= 200.0 / dT_est.max()
            H = H + dH
            if polish:
                break
            if np.abs(dH * dTdH).max() < tol_T and it > 0:
                polish = True  # one more Newton iteration: the stiff surface balance (dq/dT ~ 1/R_half)
                # turns a 1e-6 K update into a flux error of ~10 W/m^2; the extra quadratic step removes it
        else:
            return False, np.inf
        try:
            T_new = self.T_of(H)
        except ValueError:
            return False, np.inf
        T_old = self.T_of(H0)
        G, G_back = self._conductances(T_new)
        # the front flux exactly as the last (linear) Newton system used it, so the ledger matches the solve
        q = q_last + dq_last * (float(T_new[0]) - T0_last)
        self.ledger.front_in += q * self.area_face[0] * dt
        self.ledger.back_out += G_back * (T_new[-1] - self.T_back) * dt
        self.H = H
        self.t = t1
        return True, float(np.abs(T_new - T_old).max())

    def finish(self):
        self.ledger.stored_final = self.stored()
        return self.ledger
