"""Oxide-scale model for ZrB2-SiC at the stagnation point (engineering, regime-based).

Layers (from the surface): outer glass (B2O3 + SiO2) / ZrO2 skeleton with glass in its pores / SiC-depleted ZrB2 /
unreacted ZrB2-SiC (Eakins 2011 Fig. 1; Shugart 2014).  State per unit area:
    n_B, n_S   mol/m^2 of B2O3 and SiO2 in the outer glass
    x_r        m, depth of composite converted to ZrO2 (+ glass); the ZrO2 layer is taken to occupy this depth
    d_dep      m, thickness of the SiC-depleted layer below it
plus the event times of the three failure criteria.

Kinetics (DECISIONS.md D-OX-2...D-OX-8; data and fits in uhtc.data; values in INPUTS.md section 7), evaluated at
the reaction-front temperature T_f (the temperature of the first solid node below the scale):
  * ZrO2 layer:  d(x_r^2)/dt = k_p(T_f) * kp_mult * (p_O2 / p_ref)^n_p
      T_f <  T_dep (1900 K): Arrhenius fit to Shugart Table 9 (1573-1773 K, air), continued below 1573 K with the
                             activation energy of the Opila & Halbig recession constants (1600-1900 K);
      T_f >= T_dep:          Arrhenius fit to Shugart Table 9 (1923-2073 K, flowing O2), extrapolated above.
  * SiC-depleted layer, T_f >= T_dep:  d(d_dep^2)/dt = k_d (Shugart 1800 C data, no T dependence; no incubation).
  * Glass production: the B2O3 of the converted ZrB2 and the SiO2 of the converted SiC (and of the SiC removed from
    the depleted layer, assumed re-oxidised in the scale; fraction f_SiO_lost escapes as gas instead) minus the glass
    that fills the pores of the ZrO2 skeleton; the rest joins the outer layer.
  * Evaporation from the outer glass at the surface temperature T_s: B2O3(g) with activity x_B2O3, SiO(g) from
    SiO2(l) = SiO + 1/2 O2 at p_O2,w and SiO2(g), activity x_SiO2 (ideal solution), each by Hertz-Knudsen-Langmuir in
    series with boundary-layer diffusion (uhtc.evaporation; fast scalar copy below, checked in tests).
  * Surface heat terms: oxidation enthalpies (ZrB2 + 5/2 O2 -> ZrO2 + B2O3(l); SiC + 3/2 O2 -> SiO2(l) + CO) minus the
    vaporisation enthalpies, from JANAF.
Thermal resistance of the scale: d_g/k_glass + x_r/k_ZrO2 + d_dep/k_dep.  Emissivity: eps_glass while the outer glass
is at least d_eps thick, falling linearly to eps_oxide as it disappears.

Failure criteria (DECISIONS.md D-FAIL-1, D-FAIL-4):
  (a) no protective glass: more than x_ox_min (1 um) of composite has been converted but the outer glass is thinner
      than the critical thickness d_c - either lost after forming or evaporating as fast as it forms.  Baseline
      (glass_criterion = 'permeation'): the film stops limiting the oxygen supply when its permeation capacity
      K(T, x_B) p_O2,w / d exceeds the largest O2 flux the boundary layer can deliver, J_gas = g ln(1 + Y_O2)/M_O2, i.e.
      d_c = K(T_s, x_B) p_O2,w / J_gas, with K the molecular O2 permeability of the B2O3-SiO2 glass (data.O2_PERM,
      linear in the B2O3 mole fraction).  Draft rule, kept as a sensitivity case ('fixed'): d_c = d_fail (0.1 um);
  (b) SiC active oxidation: p_O2,w below the passive-to-active line at T_s:
        log10 p_PA[bar] = TR_A + TR_B/T_s - offset   (fit to the literature active-to-passive points, offset 0.65)
  (c) temperature limit: T_s >= T_limit.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from . import data as dt
from . import thermo as th

R = th.R
T_CHEM_MIN, T_CHEM_MAX = 300.0, 2990.0  # K, range of the chemistry tables (JANAF B2O3(l) ends at 3000 K)
M = th.MOLAR_MASS
M_AIR = 0.028964


@dataclass
class OxParams:
    f_SiC: float = 0.2
    # kinetics
    kp_mult: float = 1.0
    kp_E_hi: Optional[float] = None  # J/mol: activation energy above the fit range (D-OX-11 band; None = fit value)
    kp_mult_hi: float = 1.0  # multiplier on k_p above the fit range (D-OX-11 band)
    p_ref: float = 101325.0  # Pa: O2 pressure of the high-T calibration data (flowing O2, 1 atm)
    n_p: float = 0.0  # pressure exponent (D-OX-6; sensitivity 0.5)
    T_dep: float = dt.T_DEP_ONSET
    kd_mult: float = 1.0
    # activation energy of the depletion rate (J/mol; 0 = baseline, no temperature dependence).  Sensitivity case of the
    # external review (A6): k_d(T) = K_DEP kd_mult exp(-kd_E/R (1/T - 1/kd_T_ref)) (runs/final_review/kd_arrhenius.py)
    kd_E: float = 0.0
    kd_T_ref: float = 1800.0 + 273.15
    kd_T_max: Optional[float] = None  # K: above it the Arrhenius factor is held at its value there (no extrapolation)
    incubation: bool = False  # option: SiC depletion starts after an incubation time (D-OX-10, Shugart Fig. 40)
    f_SiO_lost: float = 0.0
    # pore glass: 'production' (baseline: the pores keep glass of the production composition) or 'silica' (sensitivity,
    # final review R1-11: the pores fill with silica only, so all B2O3 and less SiO2 joins the outer glass)
    pore_fill: str = "production"
    # evaporation
    alpha_evap: float = 1.0
    # thermal / optical
    k_glass: float = dt.K_SiO2_GLASS
    k_ZrO2: float = dt.K_ZrO2_SCALE
    k_dep: float = 30.0  # W/m/K, ASSUMPTION (porous ZrB2 left after SiC removal), range 10-60
    eps_glass: float = 0.88
    eps_oxide: float = 0.75
    d_eps: float = 1e-6  # m
    # SiC passive-to-active transition
    pa_offset: float = dt.PA_OFFSET_DECADES
    pa_law: str = "jacobson"  # transition line: 'jacobson' (baseline, Eq. pa) or 'balat2015' (sensitivity, D-OX-12)
    # failure
    glass_criterion: str = "permeation"  # (a): 'permeation' (D-FAIL-4) or 'fixed' (threshold d_fail; draft, D-FAIL-1)
    perm_case: str = "base"  # O2 permeability law of the glass: 'base' | 'low' | 'high' (data.O2_PERM)
    d_fail: float = 0.1e-6  # m, fixed threshold (glass_criterion = 'fixed')
    x_ox_min: float = 1e-6  # m, converted depth after which a missing glass counts as failure (a)
    T_limit: float = 2149.0 + 273.15  # K (Kolodziej 1997 reusable leading-edge limit; D-FAIL-1)
    x_r0: float = 0.05e-6  # m, initial converted depth (native oxide) to start the parabolic law
    sources: dict = field(default_factory=dict)


@dataclass
class _Tables:
    T: np.ndarray
    lnKB: np.ndarray  # B2O3(l) = B2O3(g), K in bar
    lnKSiO: np.ndarray  # SiO2(l) = SiO(g) + 1/2 O2, K in bar^1.5
    lnKS2: np.ndarray  # SiO2(l) = SiO2(g)
    dHB: np.ndarray  # J/mol
    dHSiO: np.ndarray
    dHS2: np.ndarray
    hZ: np.ndarray  # heat released per mol ZrB2 oxidised, J/mol
    hS: np.ndarray  # per mol SiC


_TAB = None
_NT = int(T_CHEM_MAX - T_CHEM_MIN) + 1


def tables():
    global _TAB
    if _TAB is None:
        T = np.arange(T_CHEM_MIN, T_CHEM_MAX + 0.5, 1.0)
        RT = R * T
        _TAB = _Tables(
            T=T,
            lnKB=-th.dG_reaction(T, {"B2O3(l)": 1}, {"B2O3(g)": 1}) / RT,
            lnKSiO=-th.dG_reaction(T, {"SiO2(l)": 1}, {"SiO": 1}) / RT,  # O2 has delta-f G = 0
            lnKS2=-th.dG_reaction(T, {"SiO2(l)": 1}, {"SiO2(g)": 1}) / RT,
            dHB=th.dH_reaction(T, {"B2O3(l)": 1}, {"B2O3(g)": 1}),
            dHSiO=th.dH_reaction(T, {"SiO2(l)": 1}, {"SiO": 1, "O2": 0.5}),
            dHS2=th.dH_reaction(T, {"SiO2(l)": 1}, {"SiO2(g)": 1}),
            hZ=-th.dH_reaction(T, {"ZrB2": 1, "O2": 2.5}, {"ZrO2": 1, "B2O3(l)": 1}),
            hS=-th.dH_reaction(T, {"SiC": 1, "O2": 1.5}, {"SiO2(l)": 1, "CO": 1}),
        )
        keys = ("lnKB", "lnKSiO", "lnKS2", "dHB", "dHSiO", "dHS2", "hZ", "hS")
        _TAB.arr = np.array([getattr(_TAB, k) for k in keys], dtype=np.float64)  # rows in this order, for _chem_pair
        for k in keys:
            setattr(_TAB, k, [float(v) for v in getattr(_TAB, k)])
    return _TAB


try:  # numba is optional (as in conduction1d); the fallback runs the same arithmetic in Python, only slower
    from numba import njit
except ImportError:  # pragma: no cover
    def njit(*a, **k):
        return a[0] if (len(a) == 1 and callable(a[0])) else (lambda f: f)


@njit(cache=True)
def _lk(arr, row, T, tmin, tmax, nt):
    """Same linear interpolation as _lookup, on row `row` of the table array."""
    x = min(max(T, tmin), tmax) - tmin
    i = int(x)
    if i >= nt - 1:
        i = nt - 2
    w = x - i
    return arr[row, i] * (1.0 - w) + arr[row, i + 1] * w


@njit(cache=True)
def _evs(p_eq, T, p_s, Mv, g, alpha, R_, M_air):
    """Same as evap_scalar."""
    if p_eq <= 0.0:
        return 0.0
    jk = alpha * p_eq * math.sqrt(Mv / (2.0 * math.pi * R_ * T))
    x = min(p_eq / p_s, 0.99)
    y = min(x * Mv / (x * Mv + (1 - x) * M_air), 0.99)
    jd = g * math.log1p(y / (1.0 - y))
    return jk * jd / (jk + jd)


@njit(cache=True)
def _chem_f(T, xB, glass, cZv, cSv, p_O2, p_s, g, alpha, arr, c):
    """q_ox - q_evap at surface temperature T: the arithmetic of LayerOxidation.chem/evap (kept there as the
    reference implementation, chem_reference) compiled for speed.  c = (tmin, tmax, nt, P_STD, R, M_B2O3, M_SiO,
    M_SiO2, M_air)."""
    tmin, tmax, nt, pst, R_, mB, mSiO, mS2, mair = c[0], c[1], int(c[2]), c[3], c[4], c[5], c[6], c[7], c[8]
    q_ox = cZv * _lk(arr, 6, T, tmin, tmax, nt) + cSv * _lk(arr, 7, T, tmin, tmax, nt)
    if not glass:
        return q_ox
    pB = xB * math.exp(_lk(arr, 0, T, tmin, tmax, nt)) * pst
    pSiO = (1 - xB) * math.exp(_lk(arr, 1, T, tmin, tmax, nt)) / math.sqrt(max(p_O2, 1e-30) / pst) * pst
    pS2 = (1 - xB) * math.exp(_lk(arr, 2, T, tmin, tmax, nt)) * pst
    Te = min(max(T, tmin), tmax)
    jB = _evs(pB, Te, p_s, mB, g, alpha, R_, mair)
    jSiO = _evs(pSiO, Te, p_s, mSiO, g, alpha, R_, mair)
    jS2 = _evs(pS2, Te, p_s, mS2, g, alpha, R_, mair)
    q = (jB / mB) * _lk(arr, 3, T, tmin, tmax, nt) + (jSiO / mSiO) * _lk(arr, 4, T, tmin, tmax, nt) \
        + (jS2 / mS2) * _lk(arr, 5, T, tmin, tmax, nt)
    return q_ox - q


_CHEM_CONSTS = np.array([T_CHEM_MIN, T_CHEM_MAX, float(_NT), th.P_STD, R, M["B2O3"], M["SiO"], M["SiO2"], M_AIR])


def _lookup(arr, T):
    """Linear interpolation on the 1 K chemistry grid (arr is a Python list for speed)."""
    x = min(max(T, T_CHEM_MIN), T_CHEM_MAX) - T_CHEM_MIN
    i = int(x)
    if i >= _NT - 1:
        i = _NT - 2
    w = x - i
    return arr[i] * (1.0 - w) + arr[i + 1] * w


def evap_scalar(p_eq, T, p_s, Mv, g, alpha):
    """Scalar copy of uhtc.evaporation.evaporation_flux (HKL in series with boundary-layer diffusion)."""
    if p_eq <= 0.0:
        return 0.0
    jk = alpha * p_eq * math.sqrt(Mv / (2.0 * math.pi * R * T))
    x = min(p_eq / p_s, 0.99)
    y = min(x * Mv / (x * Mv + (1 - x) * M_AIR), 0.99)
    jd = g * math.log1p(y / (1.0 - y))
    return jk * jd / (jk + jd)


class LayerOxidation:
    def __init__(self, p: OxParams = None):
        self.p = p or OxParams()
        p = self.p
        self.cZ = dt.RHO_ZrB2 / M["ZrB2"]  # mol/m3 of pure ZrB2
        self.cS = dt.RHO_SiC / M["SiC"]
        self.VB = M["B2O3"] / dt.RHO_B2O3_GLASS  # m3/mol
        self.VS = M["SiO2"] / dt.RHO_SiO2_GLASS
        # per m^3 of composite converted: ZrO2 volume (16.9 % expansion of the ZrB2 fraction) and glass volume
        self.v_ZrO2 = (1 - p.f_SiC) * (1 + dt.ZrO2_EXPANSION)
        self.v_glass = (1 - p.f_SiC) * self.cZ * self.VB + p.f_SiC * self.cS * self.VS  # front SiC always -> SiO2
        self.pore = max(1.0 - self.v_ZrO2, 0.0)  # glass held in the skeleton per m^3 converted
        self.frac_out = max(self.v_glass - self.pore, 0.0) / self.v_glass
        self.frac_out_B = self.frac_out_S = self.frac_out
        if p.pore_fill == "silica":
            v_SiO2 = p.f_SiC * self.cS * self.VS
            self.frac_out_B, self.frac_out_S = 1.0, max(v_SiO2 - self.pore, 0.0) / v_SiO2
        tables()

    # ---------------------------------------------------------------- state
    def init_state(self, T0):
        return dict(n_B=0.0, n_S=0.0, x_r=self.p.x_r0, d_dep=0.0, t=0.0, glass_seen=False,
                    t_glass_loss=math.nan, t_active=math.nan, t_limit=math.nan, Ts_max=T0, Tf_max=T0,
                    m_B_evap=0.0, m_S_evap=0.0)

    def glass_thickness(self, s):
        return s["n_B"] * self.VB + s["n_S"] * self.VS

    def x_B(self, s):
        n = s["n_B"] + s["n_S"]
        return s["n_B"] / n if n > 0 else 0.0

    def o2_permeability(self, T, xB):
        """Molecular O2 permeability of the B2O3-SiO2 glass, mol m^-1 s^-1 Pa^-1 (data.O2_PERM; linear in x_B)."""
        law = dt.O2_PERM[self.p.perm_case]
        kS = law["SiO2"][0] * math.exp(-law["SiO2"][1] / T)
        kB = law["B2O3"][0] * math.exp(-law["B2O3"][1] / T)
        return xB * kB + (1.0 - xB) * kS

    def delta_crit(self, s, T, flow):
        """Critical glass thickness of criterion (a) (m): below it the film no longer limits the O2 supply."""
        if self.p.glass_criterion == "fixed":
            return self.p.d_fail
        j_gas = flow.g_mass * math.log1p(dt.Y_O2_AIR) / M["O2"]  # mol/m^2/s, largest boundary-layer O2 supply
        n = s["n_B"] + s["n_S"]
        xB = s["n_B"] / n if n > 0 else s.get("_xB_last", 0.0)  # no glass left: composition of the last glass
        return self.o2_permeability(T, xB) * max(flow.p_O2, 0.0) / max(j_gas, 1e-300)

    def has_glass(self, s):
        """Glass present and protective by criterion (a) at the last advanced state (critical thickness lagged one
        step); the catalysis switch (D-HEAT-4) uses the same test as criterion (a)."""
        return not self.glass_lost(s)

    def glass_lost(self, s):
        dc = s.get("_dcrit", self.p.d_fail)
        return (self.glass_thickness(s) <= dc) and (s["x_r"] - self.p.x_r0 > self.p.x_ox_min)

    def resistance(self, s):
        p = self.p
        return self.glass_thickness(s) / p.k_glass + s["x_r"] / p.k_ZrO2 + s["d_dep"] / p.k_dep

    def emissivity(self, s):
        p = self.p
        w = min(self.glass_thickness(s) / p.d_eps, 1.0)
        e = p.eps_oxide + (p.eps_glass - p.eps_oxide) * w
        return lambda T: (e, 0.0)

    # ---------------------------------------------------------------- rates
    def kp(self, T, p_O2):
        p = self.p
        if T >= p.T_dep:
            k = dt.KP_HIGH_A * math.exp(-dt.KP_HIGH_E / (R * T))
            if T > dt.KP_FIT_TMAX and (p.kp_E_hi is not None or p.kp_mult_hi != 1.0):
                # extrapolation band above the fit range (D-OX-11): other activation energy anchored at the fit's value
                # at its upper end, and/or a multiplier
                if p.kp_E_hi is not None:
                    kT = dt.KP_HIGH_A * math.exp(-dt.KP_HIGH_E / (R * dt.KP_FIT_TMAX))
                    k = kT * math.exp(-p.kp_E_hi / R * (1.0 / T - 1.0 / dt.KP_FIT_TMAX))
                k *= p.kp_mult_hi
        else:
            T1 = dt.SHUGART_KP_LOW_T[0]
            if T >= T1:
                k = dt.KP_LOW_A * math.exp(-dt.KP_LOW_E / (R * T))
            else:
                k1 = dt.KP_LOW_A * math.exp(-dt.KP_LOW_E / (R * T1))
                k = k1 * math.exp(-dt.E_OPILA_LOW / R * (1.0 / T - 1.0 / T1))
        if p.n_p:
            k *= (max(p_O2, 1e-6) / p.p_ref) ** p.n_p
        return k * p.kp_mult

    def kd(self, T):
        if T < self.p.T_dep:
            return 0.0
        k = dt.K_DEP * self.p.kd_mult
        if self.p.kd_E:
            Ta = T if self.p.kd_T_max is None else min(T, self.p.kd_T_max)
            k *= math.exp(-self.p.kd_E / R * (1.0 / Ta - 1.0 / self.p.kd_T_ref))
        return k

    def log10_p_PA(self, T):
        """log10 of the passive-to-active transition O2 pressure (bar) at T."""
        if self.p.pa_law == "balat2015":  # p_O2 [Pa] = A exp(-B/T) (D-OX-12 sensitivity)
            return math.log10(dt.PA_BALAT_A) - dt.PA_BALAT_B / (T * math.log(10.0)) - 5.0
        if self.p.pa_law == "jacobson_fe":  # same points, study fixed-effects slope (final review R1-01, sensitivity)
            return dt.TR_FE_A + dt.TR_FE_B / T - self.p.pa_offset
        return dt.TR_A + dt.TR_B / T - self.p.pa_offset

    def T_transition(self, p_O2):
        """Surface temperature at which SiC turns active for an O2 pressure p_O2 (Pa)."""
        if self.p.pa_law == "balat2015":
            return dt.PA_BALAT_B / math.log(dt.PA_BALAT_A / p_O2)
        if self.p.pa_law == "jacobson_fe":
            return dt.TR_FE_B / (math.log10(p_O2 / 1e5) - dt.TR_FE_A + self.p.pa_offset)
        return dt.TR_B / (math.log10(p_O2 / 1e5) - dt.TR_A + self.p.pa_offset)

    def evap(self, s, T, flow):
        """(j_B2O3, j_SiO2-equivalent) kg/m^2/s and the vaporisation heat (W/m^2) at surface temperature T."""
        n = s["n_B"] + s["n_S"]
        if n <= 0.0:
            return 0.0, 0.0, 0.0
        xB = s["n_B"] / n
        t = tables()
        a = self.p.alpha_evap
        pB = xB * math.exp(_lookup(t.lnKB, T)) * th.P_STD
        pSiO = (1 - xB) * math.exp(_lookup(t.lnKSiO, T)) / math.sqrt(max(flow.p_O2, 1e-30) / th.P_STD) * th.P_STD
        pS2 = (1 - xB) * math.exp(_lookup(t.lnKS2, T)) * th.P_STD
        Te = min(max(T, T_CHEM_MIN), T_CHEM_MAX)
        jB = evap_scalar(pB, Te, flow.p_s, M["B2O3"], flow.g_mass, a)
        jSiO = evap_scalar(pSiO, Te, flow.p_s, M["SiO"], flow.g_mass, a)
        jS2 = evap_scalar(pS2, Te, flow.p_s, M["SiO2"], flow.g_mass, a)
        q = (jB / M["B2O3"]) * _lookup(t.dHB, T) + (jSiO / M["SiO"]) * _lookup(t.dHSiO, T) \
            + (jS2 / M["SiO2"]) * _lookup(t.dHS2, T)
        nS = jSiO / M["SiO"] + jS2 / M["SiO2"]
        return jB, nS * M["SiO2"], q

    def front_velocity(self, s, T_f, flow):
        return self.kp(T_f, flow.p_O2) / (2.0 * s["x_r"])

    def chem(self, s, flow):
        """Surface heat terms as a function of the surface temperature, (q, dq/dT); the oxidation heat uses the
        front velocity of the current state (lagged), the evaporation term is evaluated at T.  Compiled version of
        chem_reference (same arithmetic; tests/test_oxidation.py compares the two)."""
        t = tables()
        v = s.get("_v_last", 0.0)
        cZv = (1 - self.p.f_SiC) * self.cZ * v
        cSv = self.p.f_SiC * self.cS * v
        n = s["n_B"] + s["n_S"]
        xB = s["n_B"] / n if n > 0.0 else 0.0
        args = (xB, n > 0.0, cZv, cSv, float(flow.p_O2), float(flow.p_s), float(flow.g_mass),
                float(self.p.alpha_evap), t.arr, _CHEM_CONSTS)

        def fd(T, d=1.0):
            f0 = _chem_f(T, *args)
            return f0, (_chem_f(T + d, *args) - f0) / d
        return fd

    def chem_reference(self, s, flow):
        """Pure-Python reference of chem (the original implementation)."""
        t = tables()
        v = s.get("_v_last", 0.0)
        cZv = (1 - self.p.f_SiC) * self.cZ * v
        cSv = self.p.f_SiC * self.cS * v

        def f(T):
            q_ox = cZv * _lookup(t.hZ, T) + cSv * _lookup(t.hS, T)
            return q_ox - self.evap(s, T, flow)[2]

        def fd(T, d=1.0):
            # forward difference: the derivative only steers the Newton iterations of the surface balance and of the
            # conduction step (the converged solution does not depend on it)
            f0 = f(T)
            return f0, (f(T + d) - f0) / d
        return fd

    # ---------------------------------------------------------------- advance
    def advance(self, s, dt_, T_s, T_f, flow):
        s = dict(s)
        p = self.p
        k = self.kp(T_f, flow.p_O2)
        x_new = math.sqrt(s["x_r"] ** 2 + k * dt_)
        dx = x_new - s["x_r"]
        s["_v_last"] = dx / dt_
        s["x_r"] = x_new
        PB = dx * (1 - p.f_SiC) * self.cZ * self.frac_out_B  # mol/m^2 produced over the step
        PS = dx * p.f_SiC * self.cS * self.frac_out_S  # SiC converted at the front always forms SiO2 glass
        kd = self.kd(T_f)
        if kd > 0.0 and p.incubation and s.get("inc", 0.0) < 1.0:
            # incubation clock (additivity rule, our assumption): depletion starts when sum dt / t_inc(T_f) reaches 1
            t_inc = 60.0 * float(np.interp(T_f, dt.INC_T, dt.INC_MIN))
            s["inc"] = s.get("inc", 0.0) + dt_ / t_inc
            kd = 0.0 if s["inc"] < 1.0 else kd * min((s["inc"] - 1.0) * t_inc / dt_, 1.0)
        if kd > 0.0:
            d_new = math.sqrt(s["d_dep"] ** 2 + kd * dt_)
            PS += (d_new - s["d_dep"]) * p.f_SiC * self.cS * (1 - p.f_SiO_lost)
            s["d_dep"] = d_new
        # Evaporation, backward Euler.  The flux from the liquid film does not depend on its thickness, only on its
        # composition: j_i = a_i x_i M_i with a_i the molar flux at unit activity (linearised at the composition of the
        # current glass mixed with this step's production).  With A = n_B + P_B, B = n_S + P_S and c_i = dt a_i,
        #   n_B' = A N/(N + c_B),  n_S' = B N/(N + c_S),  N = n_B' + n_S' solving A/(N + c_B) + B/(N + c_S) = 1,
        # or N = 0 when no positive root exists (the film evaporates completely within the step).  Exact for the
        # linearised law and unconditionally stable (B2O3 leaves on a millisecond time scale).
        nB0, nS0 = s["n_B"], s["n_S"]
        g0 = nB0 + nS0
        A, B = nB0 + PB, nS0 + PS
        nB1 = nS1 = 0.0
        if A + B > 0:
            xm = A / (A + B)
            jB, jS, _ = self.evap(dict(s, n_B=A, n_S=B), T_s, flow)
            aB = jB / (M["B2O3"] * xm) if xm > 0 else 0.0
            aS = jS / (M["SiO2"] * (1 - xm)) if xm < 1 else 0.0
            cB, cS = dt_ * aB, dt_ * aS

            def h(N):
                return (A / (N + cB) if A > 0 else 0.0) + (B / (N + cS) if B > 0 else 0.0) - 1.0
            if (A > 0 and cB == 0.0) or (B > 0 and cS == 0.0) or h(0.0) > 0.0:
                lo, hi = 0.0, A + B  # h(lo) > 0 (or +inf), h(hi) <= 0
                N = hi
                for _ in range(200):  # safeguarded Newton on a monotone function
                    hv = h(N)
                    dh = -(A / (N + cB) ** 2 if A > 0 else 0.0) - (B / (N + cS) ** 2 if B > 0 else 0.0)
                    if hv > 0:
                        lo = N
                    else:
                        hi = N
                    Nn = N - hv / dh if dh != 0 else 0.5 * (lo + hi)
                    if not (lo < Nn < hi):
                        Nn = 0.5 * (lo + hi)
                    if abs(Nn - N) <= 1e-14 * max(A + B, 1e-300):
                        N = Nn
                        break
                    N = Nn
                nB1 = A * N / (N + cB) if A > 0 else 0.0
                nS1 = B * N / (N + cS) if B > 0 else 0.0
        s["m_B_evap"] += (nB0 + PB - nB1) * M["B2O3"]
        s["m_S_evap"] += (nS0 + PS - nS1) * M["SiO2"]
        s["n_B"], s["n_S"] = nB1, nS1
        s["_dn_dt"] = ((nB1 + nS1) - g0) / dt_
        t0 = s["t"]
        s["t"] += dt_
        # failure times are interpolated linearly inside the step from the start-of-step values (the step may be long
        # compared with the crossing); dg0 / Ts0 are the glass thickness and surface temperature at the step start
        dg0 = s.get("_dg_prev", 0.0)
        Ts0 = s.get("_Ts_prev", T_s)

        def cross(v0, v1, lev):
            return t0 + dt_ * min(max((lev - v0) / (v1 - v0), 0.0), 1.0) if v1 != v0 else s["t"]
        dg = self.glass_thickness(s)
        dc = self.delta_crit(s, T_s, flow)  # critical thickness of criterion (a) at the end of the step
        if s["n_B"] + s["n_S"] > 0:
            s["_xB_last"] = s["n_B"] / (s["n_B"] + s["n_S"])
        dc0 = s.get("_dcrit", dc)
        s["_dcrit"] = dc
        s["_dcrit_min"] = min(s.get("_dcrit_min", dc), dc)
        s["_dcrit_max"] = max(s.get("_dcrit_max", dc), dc)
        if dg > dc:
            s["glass_seen"] = True
        elif s["x_r"] - p.x_r0 > p.x_ox_min and math.isnan(s["t_glass_loss"]):
            # the surface has oxidised but carries no protective glass: lost after forming, or evaporated as fast as
            # formed; crossing of dg - dc through zero interpolated inside the step
            s["t_glass_loss"] = cross(dg0 - dc0, dg - dc, 0.0) if dg0 > dc0 else s["t"]
        lp = math.log10(max(flow.p_O2, 1e-30) / 1e5)
        if math.isnan(s["t_active"]) and lp < self.log10_p_PA(T_s):
            g0_, g1_ = lp - self.log10_p_PA(Ts0), lp - self.log10_p_PA(T_s)
            s["t_active"] = cross(g0_, g1_, 0.0) if g0_ >= 0 else s["t"]
        if math.isnan(s["t_limit"]) and T_s >= p.T_limit:
            s["t_limit"] = cross(Ts0, T_s, p.T_limit) if Ts0 < p.T_limit else s["t"]
        s["_dg_prev"], s["_Ts_prev"] = dg, T_s
        s["Ts_max"] = max(s["Ts_max"], T_s)
        s["Tf_max"] = max(s["Tf_max"], T_f)
        return s

    def rate_limit(self, s, T_s, flow):
        """Largest step keeping the relative change of the converted depth and of the glass inventory below ~5 %."""
        v = self.kp(min(s["Tf_max"] + 50.0, T_CHEM_MAX + 500.0), flow.p_O2) / (2.0 * s["x_r"])
        lim = 0.05 * s["x_r"] / max(v, 1e-30)
        n = s["n_B"] + s["n_S"]
        dndt = s.get("_dn_dt", 0.0)
        if n > 0 and dndt < 0 and self.glass_thickness(s) < 5e-6:
            # while a thin glass (< 5 um) is thinning, at most 10 % of it (but not less than 0.05 um) per step, so that
            # the time of glass loss is resolved; the evaporation update itself is unconditionally stable
            lim = min(lim, max(0.10 * n, 0.05e-6 / self.VS) / -dndt)
        return max(lim, 1e-4)

    def events(self, s):
        return {k: s[k] for k in ("t_glass_loss", "t_active", "t_limit")}

    def layers(self, s):
        """Thicknesses (m): outer glass, ZrO2 layer, depleted layer, total affected depth below the original surface,
        and the outward shift of the outer surface (glass + ZrO2 swelling)."""
        dg = self.glass_thickness(s)
        return dict(glass=dg, ZrO2=s["x_r"], depleted=s["d_dep"], consumed=s["x_r"] + s["d_dep"],
                    surface_shift=dg, x_B2O3=self.x_B(s))
