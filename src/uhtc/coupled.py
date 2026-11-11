"""Coupled time integration: stagnation-line conduction (1D) + surface balance + oxide-scale model.

Operator splitting with a short step: in each step the conduction problem is solved implicitly with the
surface balance (uhtc.surface) using the oxide state frozen at the start of the step (scale resistance,
emissivity, chemical heat terms); then the oxide state is advanced over the same step with the new surface
temperature.  The step is limited so that the surface temperature changes by at most dT_max and the scale
thicknesses by at most a fraction f_ox per step (DECISIONS.md D-NUM-1), which keeps the splitting error small
(checked by halving both limits in V2).

Flow conditions come from a `Conditions` object: constant (sweeps) or interpolated along a trajectory.
Wall enthalpy: equilibrium air for the fully catalytic bound, frozen air for the non-catalytic bound
(DECISIONS.md D-FLOW-2), tabulated on a (T, log p) grid once and interpolated.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Callable, Optional

import numpy as np

from . import heating as ht
from . import stagnation as st
from . import surface as sf
from . import thermo as th
from .conduction1d import Solver1D

CACHE = Path(__file__).resolve().parents[2] / "runs" / "_cache"


# ------------------------------------------------------------------------------------ wall enthalpy
@lru_cache(maxsize=1)
def _heq_table():
    """h_eq(T, p) of equilibrium air on T = 200..5000 K (10 K) x log10 p = 2..6.6 (0.1)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / "h_eq_air_table_v1.npz"
    T = np.arange(200.0, 5000.0 + 1e-9, 10.0)
    lp = np.arange(2.0, 6.6 + 1e-9, 0.1)
    if f.exists():
        d = np.load(f)
        if d["T"].shape == T.shape and d["lp"].shape == lp.shape:
            return d["T"], d["lp"], d["h"]
    h = np.empty((len(T), len(lp)))
    for j, l in enumerate(lp):
        for i, t in enumerate(T):
            h[i, j] = th.h_air_equilibrium(t, 10 ** l) if t > 1200 else th.h_air_frozen(t)
    np.savez(f, T=T, lp=lp, h=h)
    return T, lp, h


@lru_cache(maxsize=1)
def _hfr_table():
    T = np.arange(200.0, 5000.0 + 1e-9, 10.0)
    return T, th.h_air_frozen(T)


class WallEnthalpy:
    """h_w(T) and dh_w/dT at a fixed pressure for the chosen catalysis bound."""

    def __init__(self, p, catalytic=True):
        if catalytic:
            T, lp, h = _heq_table()
            l = np.clip(np.log10(p), lp[0], lp[-1])
            j = min(int((l - lp[0]) / 0.1), len(lp) - 2)
            w = (l - lp[j]) / (lp[j + 1] - lp[j])
            self.T, self.h = T, (1 - w) * h[:, j] + w * h[:, j + 1]
        else:
            self.T, self.h = _hfr_table()
        self.dh = np.gradient(self.h, self.T)
        # uniform 10 K grid: index arithmetic instead of np.interp per scalar call (same linear interpolation)
        self._T0, self._dT, self._n = float(self.T[0]), float(self.T[1] - self.T[0]), len(self.T)
        self._Tl, self._hl, self._dhl = [float(x) for x in self.T], [float(x) for x in self.h], [float(x) for x in self.dh]

    def __call__(self, T):
        Tl, hl, dhl = self._Tl, self._hl, self._dhl
        if T <= Tl[0]:
            return hl[0], dhl[0]
        if T >= Tl[-1]:
            return hl[-1], dhl[-1]
        j = min(int((T - self._T0) / self._dT), self._n - 2)
        if T < Tl[j]:  # guard against round-off in the index
            j -= 1
        x0, x1 = Tl[j], Tl[j + 1]
        return ((hl[j + 1] - hl[j]) / (x1 - x0) * (T - x0) + hl[j],
                (dhl[j + 1] - dhl[j]) / (x1 - x0) * (T - x0) + dhl[j])


def convective_flux(flow, catalytic, hw_eq, hw_fr):
    """Hot-wall convective flux q(T_w) and dq/dT_w for the chosen catalysis bound (D-HEAT-2).

    Fully catalytic: Phi (H_s - h_eq(T_w)).  Non-catalytic: Phi f_nc (H_s - h_frozen(T_w)) with the Fay-Riddell
    Le = 0 factor f_nc, capped at the fully catalytic value: the multiplicative cold-wall form would otherwise exceed
    the catalytic flux near convective blocking (wall temperature approaching the edge temperature), which is
    impossible (a catalytic wall receives at least the non-catalytic heat; DECISIONS D-HEAT-2)."""
    Phi, H = flow.Phi, flow.H_s
    if catalytic:
        def conv(T):
            h, dh = hw_eq(T)
            return Phi * (H - h), -Phi * dh
        return conv
    f = flow.f_nc

    def conv(T):
        h1, d1 = hw_eq(T)
        h2, d2 = hw_fr(T)
        q_fc, q_nc = Phi * (H - h1), Phi * f * (H - h2)
        if q_nc <= q_fc:
            return q_nc, -Phi * f * d2
        return q_fc, -Phi * d1
    return conv


# ------------------------------------------------------------------------------------ conditions
@dataclass
class FlowState:
    H_s: float
    p_s: float
    p_O2: float
    h_D: float
    Phi: float  # kg/(m^2 s): K_SG sqrt(p_s/R_n) x geometry factor (catalysis applied separately)
    f_nc: float  # non-catalytic factor (Fay-Riddell bracket with Le = 0, D-HEAT-2)
    g_mass: float  # mass-transfer conductance kg/(m^2 s) (D-HEAT-3)
    mach: float = 0.0
    altitude: float = 0.0


class Conditions:
    """Flow at the stagnation point as a function of time.

    catalysis: 'fc' fully catalytic bound, 'nc' non-catalytic bound, 'switch' non-catalytic while the outer glass
    exists and fully catalytic after it is lost (D-HEAT-4); also accepted: catalytic=True/False."""

    def __init__(self, R_n, catalytic=True, cylinder=False, mass_transfer_factor=1.0, catalysis=None,
                 heat_mult=1.0):
        self.catalysis = catalysis or ("fc" if catalytic else "nc")
        if self.catalysis not in ("fc", "nc", "switch"):
            raise ValueError(self.catalysis)
        self.catalytic = self.catalysis == "fc"
        self.R_n, self.cylinder = R_n, cylinder
        self.mtf = mass_transfer_factor
        self.heat_mult = heat_mult

    def _flow(self, s: st.StagnationState) -> FlowState:
        geo = (ht.C_TAUBER_CYL / ht.C_TAUBER_SPHERE) if self.cylinder else 1.0
        Phi = ht.K_SG_AIR * np.sqrt(s.p_s / ht.ATM / self.R_n) * geo
        f = float(ht.noncatalytic_factor(s.H_s, s.h_D))
        return FlowState(s.H_s, s.p_s, s.p_O2_w, s.h_D, Phi * self.heat_mult, f, Phi * self.mtf, s.mach, s.altitude)

    def at(self, t) -> FlowState:  # pragma: no cover - overridden
        raise NotImplementedError

    def next_change(self, t):
        return np.inf


class ConstantConditions(Conditions):
    def __init__(self, mach, altitude, R_n, state=None, **kw):
        super().__init__(R_n, **kw)
        self.state = state if state is not None else st.stagnation_state(mach, altitude)
        self.flow = self._flow(self.state)

    def at(self, t):
        return self.flow


class TrajectoryConditions(Conditions):
    """Stagnation conditions along a trajectory (t, altitude, Mach), evaluated on the trajectory's own time
    grid (decimated to dt_eval) and interpolated linearly in time."""

    def __init__(self, t, alt, mach, R_n, dt_eval=2.0, state_fn=None, **kw):
        super().__init__(R_n, **kw)
        te = np.arange(t[0], t[-1] + 1e-9, dt_eval)
        fn = state_fn or st.stagnation_state
        flows = [self._flow(fn(round(float(np.interp(x, t, mach)), 4), round(float(np.interp(x, t, alt)), 1)))
                 for x in te]
        self.te = te
        self.arr = {k: np.array([getattr(f, k) for f in flows]) for k in FlowState.__dataclass_fields__}

    def at(self, t):
        return FlowState(**{k: float(np.interp(t, self.te, v)) for k, v in self.arr.items()})


T_EDGE_MAX = 9000.0  # K, upper limit of the extended air tables for ground-test edge states (ionisation neglected)


class PrescribedConditions(Conditions):
    """Piecewise-constant stagnation conditions given directly as total enthalpy and stagnation pressure (ground
    tests: arc-jet steps).  steps = [(t_start, H_s [J/kg], p_s [Pa]), ...] or, for a measured cold-wall heat flux
    on the same geometry, [(t_start, H_s, p_s, q_cw [W/m^2]), ...]: then Phi = q_cw / (H_s - h(300 K)) replaces the
    Sutton-Graves value (V8).  x_O2 = O2 mole fraction of the test gas with all oxygen counted as O2 (D-OX-1).  The
    edge dissociation enthalpy h_D comes from equilibrium air at (p_s, H_s)."""

    def __init__(self, steps, R_n, x_O2=0.2095, **kw):
        super().__init__(R_n, **kw)
        self.t0 = np.array([s[0] for s in steps], float)
        self.flows = []
        self.edge_T = []  # equilibrium edge temperature of each step, K
        self.edge_held = []  # True where H lies above the extended tables and the edge state is held at T_EDGE_MAX
        for st in steps:
            _, H, p = st[:3]
            # arc-jet edge states can lie beyond the 6000 K end of the JANAF tables: the tables are extended to
            # T_EDGE_MAX (thermo.extended_tables; neutral species, cp held at 6000 K).  Above T_EDGE_MAX the state is
            # held at T_EDGE_MAX, where the neutral air is fully dissociated, so h_D is then close to its neutral limit
            # (ionisation not included).  h_D enters only the non-catalytic factor.  Before the final review the state
            # was clipped at 5998 K (V8, review/FINDINGS.csv).
            with th.extended_tables(T_EDGE_MAX):
                h_max = float(th.h_air_equilibrium(T_EDGE_MAX - 1.0, p))
                Te = th.T_air_equilibrium(min(H, h_max), p, T_hi=T_EDGE_MAX - 1.0)
                _, _, hD = th.mixture_props(th.equilibrium_air(Te, p), Te)
            self.edge_T.append(float(Te))
            self.edge_held.append(bool(H > h_max))
            geo = (ht.C_TAUBER_CYL / ht.C_TAUBER_SPHERE) if self.cylinder else 1.0
            Phi = ht.K_SG_AIR * np.sqrt(p / ht.ATM / self.R_n) * geo
            if len(st) > 3:  # measured cold-wall flux (fully catalytic calorimeter at ~300 K)
                Phi = st[3] / (H - float(th.h_air_frozen(300.0)))
            self.flows.append(FlowState(H, p, x_O2 * p, hD, Phi * self.heat_mult,
                                        float(ht.noncatalytic_factor(H, hD)), Phi * self.mtf))

    def at(self, t):
        i = int(np.clip(np.searchsorted(self.t0, t, side="right") - 1, 0, len(self.t0) - 1))
        return self.flows[i]

    def next_change(self, t):
        later = self.t0[self.t0 > t + 1e-9]
        return float(later[0]) if len(later) else np.inf


# ------------------------------------------------------------------------------------ oxide model interface
class NoOxide:
    """Inert surface with a fixed emissivity (used by tests and by the thermal-only comparison runs)."""

    def __init__(self, eps=0.8):
        self.eps = eps

    def has_glass(self, s):
        return True  # inert surface: never "loses its glass", so the catalysis switch keeps the non-catalytic bound

    def init_state(self, T0):
        return {}

    def resistance(self, s):
        return 0.0

    def emissivity(self, s):
        return lambda T: (self.eps, 0.0)

    def chem(self, s, flow):
        return None

    def advance(self, s, dt, T_s, T_sub, flow):
        return s

    def rate_limit(self, s, T_s, flow):
        return np.inf

    def events(self, s):
        return {}


# ------------------------------------------------------------------------------------ driver
@dataclass
class RunResult:
    t: np.ndarray
    T_s: np.ndarray
    T_sub: np.ndarray
    q_conv: np.ndarray
    states: list
    events: dict
    stopped: str
    energy_residual: float
    extra: dict = field(default_factory=dict)


def run_1d(cond: Conditions, material, oxide, t_end, T_init=300.0, depth=0.02, dx0=2e-6, growth=1.08,
           dT_max=10.0, f_ox=0.05, dt_max=2.0, dt0=1e-4, stop_on: Optional[Callable] = None,
           record_every=0.0, back="adiabatic"):
    """Integrate to t_end (or until stop_on(events, T_s) returns a reason string)."""
    s_ox = oxide.init_state(T_init)
    solver = Solver1D(material, depth, T_init, dx0=dx0, growth=growth, back=back)
    cache = {}

    def make_front(flow, s_ox):
        cat = cond.catalysis == "fc" or (cond.catalysis == "switch" and not oxide.has_glass(s_ox))

        def hw_tab(c):
            h = cache.get(("hw", flow.p_s, c))
            if h is None:
                h = WallEnthalpy(flow.p_s, c)
                if len(cache) > 8:
                    cache.clear()
                cache[("hw", flow.p_s, c)] = h
            return h
        conv = convective_flux(flow, cat, hw_tab(True), hw_tab(False))
        R_sc = oxide.resistance(s_ox)
        eps_fn = oxide.emissivity(s_ox)
        chem_fn = oxide.chem(s_ox, flow)
        info = {}

        def front(T0, t):
            R = R_sc + solver.half_cell_resistance()
            q, dq, Ts = sf.balance(T0, R, conv, eps_fn, chem_fn)
            info["Ts"] = Ts
            info["q_conv"] = conv(Ts)[0]
            return q, dq
        return front, info

    ts, Tss, Tsub, qc, states = [], [], [], [], []
    dt = dt0
    last_rec = -np.inf
    stopped = "t_end"
    Ts_prev = T_init
    while solver.t < t_end * (1 - 1e-12):
        flow = cond.at(solver.t)
        front, info = make_front(flow, s_ox)
        solver.front = front
        lim = min(dt_max, oxide.rate_limit(s_ox, Ts_prev, flow), cond.next_change(solver.t) - solver.t, t_end - solver.t)
        dt = min(dt, lim)
        H_save, t_save = solver.H.copy(), solver.t
        led_save = (solver.ledger.front_in, solver.ledger.back_out)
        ok, dT = solver.step(dt)
        if ok:
            Ts_new = info["Ts"]
            if abs(Ts_new - Ts_prev) > 2 * dT_max and dt > 1e-6:
                ok = False
        if not ok:
            solver.H, solver.t = H_save, t_save
            solver.ledger.front_in, solver.ledger.back_out = led_save
            dt *= 0.3
            if dt < 1e-9:
                raise RuntimeError(f"time step underflow at t = {solver.t}")
            continue
        T_sub = float(solver.T[0])
        s_ox = oxide.advance(s_ox, dt, Ts_new, T_sub, flow)
        dTs = abs(Ts_new - Ts_prev)
        Ts_prev = Ts_new
        reason = stop_on(oxide.events(s_ox), Ts_new, solver.t) if stop_on is not None else None
        if solver.t - last_rec >= record_every or solver.t >= t_end * (1 - 1e-12) or reason:
            ts.append(solver.t)
            Tss.append(Ts_new)
            Tsub.append(T_sub)
            qc.append(info["q_conv"])
            states.append(dict(s_ox))
            last_rec = solver.t
        if reason:
            stopped = reason
            break
        dt *= min(1.5, max(0.5, 0.8 * dT_max / max(dTs, 1e-9)))
    led = solver.finish()
    return RunResult(np.array(ts), np.array(Tss), np.array(Tsub), np.array(qc), states, oxide.events(s_ox),
                     stopped, led.residual()[1], {"solver": solver})


def first_failure(events):
    """(time, mode) of the earliest failure among the stored event times, or (nan, 'none')."""
    names = {"t_glass_loss": "glass", "t_active": "active", "t_limit": "temperature"}
    best = (np.nan, "none")
    for k, m in names.items():
        v = events.get(k, np.nan)
        if v == v and (best[0] != best[0] or v < best[0]):
            best = (float(v), m)
    return best


def stop_at_first_failure(events, T_s, t):
    tf, mode = first_failure(events)
    return mode if tf == tf else None


# ------------------------------------------------------------------------------------ 2D driver (Q5)
def run_2d(cond: Conditions, geom, material, oxide, t_end, T_init=300.0, mesh_kw=None, dT_max=10.0, dt_max=2.0,
           dt0=1e-4, stop_on: Optional[Callable] = None, record_every=0.0, flank_mult=1.0, heat_dist="tauber",
           local_pressure=False, glass_flow=None):
    """Quasi-2D nose region (uhtc.conduction2d) with the same surface model as run_1d, applied node by node.

    Local heating Phi(s) = Phi_s g(s) (cosine on the arc, flank correlation, D-2D-1); local mass-transfer
    conductance scaled the same way; the O2 pressure of the stagnation point is used on the whole surface
    (D-2D-3).  Each surface node has its own oxide state; the failure clocks are those of the stagnation node.
    flank_mult scales g(s) on the flank only (sensitivity case of D-2D-1); heat_dist names the heating
    distribution (uhtc.conduction2d.HEAT_DISTRIBUTIONS).

    local_pressure (check of the external review, item A4; off in every production run): the oxide model of each
    node sees the local wall pressure p_s (p/p_s) and O2 pressure p_O2 (p/p_s) instead of the stagnation values, with
    the modified-Newtonian ratio p/p_s = (1 - eps) cos^2(phi) + eps, eps = p_inf/p_s, that the Lees distribution
    uses (phi = body angle, 90 deg - theta_c on the flank).  The convective heating is unchanged.

    glass_flow (version 3, Task 1; None in every run of the paper): a uhtc.glassflow.FlowParams switches on the thin-film flow
    of the outer glass between surface nodes after each oxide step (uhtc.glassflow, D-FLOW-3 to D-FLOW-7)."""
    from .conduction2d import Solver2D

    sol = Solver2D(geom, material, T_init, mesh_kw=mesh_kw)
    ns = len(sol.surf)
    states = [oxide.init_state(T_init) for _ in range(ns)]
    gf = None
    if glass_flow is not None:
        from . import glassflow
        gf = dict(mod=glassflow, surf=glassflow.surface_from_solver(sol), ledger=glassflow.FlowLedger(), audit_max=0.0)
    flow0 = cond.at(0.0)
    shape_cache = {}
    pr_cache = {}

    def pratio_for(flow):
        """Local-to-stagnation wall-pressure ratio of every surface node (ones unless local_pressure)."""
        if not local_pressure or flow.mach <= 0:
            return np.ones(ns)
        key = (round(flow.mach, 3), round(flow.altitude, 0))
        if key not in pr_cache:
            from . import atmosphere as atm
            _, p_inf, _ = atm.state(flow.altitude)
            eps = min(float(p_inf) / flow.p_s, 1.0)
            pr_cache[key] = (1.0 - eps) * np.cos(sol.phi) ** 2 + eps
            if len(pr_cache) > 64:
                pr_cache.pop(next(iter(pr_cache)))
        return pr_cache[key]

    def gshape_for(flow):
        """Heating distribution for the current flight state (flank law depends weakly on V and rho_inf)."""
        key = (round(flow.mach, 3), round(flow.altitude, 0))
        if key not in shape_cache:
            if flow.mach > 0:
                from . import atmosphere as atm
                Tinf, _, rho_inf = atm.state(flow.altitude)
                V = flow.mach * np.sqrt(atm.GAMMA * atm.R_AIR * Tinf)
                q_s_cold = (ht.tauber_sphere(rho_inf, V, geom.R_n, 0.0, 1.0) if geom.axisymmetric
                            else ht.C_TAUBER_CYL * np.sqrt(rho_inf / geom.R_n) * V ** 3)
                g = sol.heating_shape(V, rho_inf, q_s_cold, dist=heat_dist, flow=flow).copy()
                g[sol.on_flank] *= flank_mult
                shape_cache[key] = g
            else:  # ground test given by (H_s, p_s) only: cosine on the arc, no flank heating law
                g = ht.cosine_distribution(sol.phi).copy()
                g[sol.on_flank] = 0.0
                shape_cache[key] = g
            if len(shape_cache) > 64:
                shape_cache.pop(next(iter(shape_cache)))
        return shape_cache[key]
    gshape = gshape_for(flow0)
    hw_tabs = {}

    def hw_for(p, cat):
        k = (round(p, 3), cat)
        if k not in hw_tabs:
            hw_tabs[k] = WallEnthalpy(p, cat)
        return hw_tabs[k]

    def surf_flux_factory(flow):
        conv_fns, eps_fns, chem_fns, Rs = [], [], [], []
        gshape = gshape_for(flow)
        pr = pratio_for(flow)
        for i in range(ns):
            s_i = states[i]
            cat = cond.catalysis == "fc" or (cond.catalysis == "switch" and not oxide.has_glass(s_i))
            fl_loc = FlowState(flow.H_s, flow.p_s, flow.p_O2, flow.h_D, flow.Phi * gshape[i], flow.f_nc,
                               flow.g_mass, flow.mach, flow.altitude)
            conv = convective_flux(fl_loc, cat, hw_for(flow.p_s, True), hw_for(flow.p_s, False))
            fl_i = FlowState(flow.H_s, flow.p_s * pr[i], flow.p_O2 * pr[i], flow.h_D, flow.Phi * gshape[i],
                             flow.f_nc, flow.g_mass * gshape[i], flow.mach, flow.altitude)
            conv_fns.append(conv)
            eps_fns.append(oxide.emissivity(s_i))
            chem_fns.append(oxide.chem(s_i, fl_i))
            Rs.append(oxide.resistance(s_i))
        out = {"Ts": np.zeros(ns)}

        def surf_flux(Tn):
            q = np.empty(ns)
            dq = np.empty(ns)
            for i in range(ns):
                q[i], dq[i], out["Ts"][i] = sf.balance(float(Tn[i]), Rs[i], conv_fns[i], eps_fns[i], chem_fns[i])
            return q, dq
        return surf_flux, out

    ts, Ts_stag, Tmax_body, recs = [], [], [], []
    dt, last_rec, stopped = dt0, -np.inf, "t_end"
    Ts_prev = np.full(ns, T_init)
    while sol.t < t_end * (1 - 1e-12):
        flow = cond.at(sol.t)
        surf_flux, out = surf_flux_factory(flow)
        lim = min(dt_max, oxide.rate_limit(states[0], Ts_prev[0], flow), cond.next_change(sol.t) - sol.t, t_end - sol.t)
        dt = min(dt, lim)
        H_save, t_save, E_save = sol.H.copy(), sol.t, sol.energy_in
        ok, _ = sol.step(dt, surf_flux)
        if ok:
            Ts_new = out["Ts"].copy()
            if np.abs(Ts_new - Ts_prev).max() > 2 * dT_max and dt > 1e-6:
                ok = False
        if not ok:
            sol.H, sol.t, sol.energy_in = H_save, t_save, E_save
            dt *= 0.3
            if dt < 1e-9:
                raise RuntimeError(f"2D time step underflow at t = {sol.t}")
            continue
        Tn = sol.T[sol.surf]
        gshape = gshape_for(flow)
        pr = pratio_for(flow)
        before = [dict(states[i]) for i in range(ns)] if gf is not None else None
        for i in range(ns):
            fl_i = FlowState(flow.H_s, flow.p_s * pr[i], flow.p_O2 * pr[i], flow.h_D, flow.Phi * gshape[i],
                             flow.f_nc, flow.g_mass * gshape[i], flow.mach, flow.altitude)
            states[i] = oxide.advance(states[i], dt, float(Ts_new[i]), float(Tn[i]), fl_i)
        if gf is not None:
            _glass_flow_step(gf, glass_flow, cond, flow, gshape, oxide, before, states, Ts_new, dt)
        dTs = np.abs(Ts_new - Ts_prev).max()
        Ts_prev = Ts_new
        reason = stop_on(oxide.events(states[0]), Ts_new[0], sol.t) if stop_on is not None else None
        if sol.t - last_rec >= record_every or sol.t >= t_end * (1 - 1e-12) or reason:
            if gf is not None:  # mass audit of the film at every output time
                film = float(np.sum(gf["surf"].w * np.array([oxide.glass_thickness(s_) for s_ in states])))
                gf["audit_max"] = max(gf["audit_max"], gf["ledger"].audit(film))
            ts.append(sol.t)
            Ts_stag.append(Ts_new[0])
            Tmax_body.append(Ts_new.max())
            recs.append(dict(states[0]))
            last_rec = sol.t
        if reason:
            stopped = reason
            break
        dt *= min(1.5, max(0.5, 0.8 * dT_max / max(dTs, 1e-9)))
    extra = {"solver": sol, "states": states, "gshape": gshape, "Ts_surface": Ts_prev}
    if gf is not None:
        extra["glass_flow"] = dict(ledger=gf["ledger"], audit_max=gf["audit_max"], surface=gf["surf"])
    return RunResult(np.array(ts), np.array(Ts_stag), np.array(Tmax_body), np.array([]), recs,
                     oxide.events(states[0]), stopped, sol.energy_residual(), extra)


def _glass_flow_step(gf, params, cond, flow, gshape, oxide, before, states, Ts, dt):
    """Film flow over one oxide step (uhtc.glassflow): ledger, flow between nodes, composition, criterion (a) after the flow."""
    from .oxidation import M as MM
    g = gf["mod"]
    surf, led = gf["surf"], gf["ledger"]
    n = len(states)
    v0 = np.array([oxide.glass_thickness(b) for b in before])
    v1 = np.array([oxide.glass_thickness(s_) for s_ in states])
    ev = np.array([(s_["m_B_evap"] - b["m_B_evap"]) / MM["B2O3"] * oxide.VB + (s_["m_S_evap"] - b["m_S_evap"]) / MM["SiO2"] * oxide.VS
                   for s_, b in zip(states, before)])
    led.source += float(np.sum(surf.w * (v1 - v0 + ev)))
    led.evap += float(np.sum(surf.w * ev))
    rho_e, p_inf = g.edge_state(cond, flow)
    tau, G = g.drivers(surf, flow, gshape, rho_e, p_inf, params)
    mu = g.viscosity(np.asarray(Ts, float), params.law)
    d, shed, nsub = g.film_step(v1, surf, tau, G, mu, dt, params)
    led.shed += shed
    if nsub < 0:
        led.n_implicit += 1
    else:
        led.n_explicit += 1
        led.max_sub = max(led.max_sub, nsub)
    p = oxide.p
    for i in range(n):
        s_ = dict(states[i])
        nb, nsv = s_["n_B"], s_["n_S"]
        if v1[i] > 0.0:
            f = d[i] / v1[i]
            s_["n_B"], s_["n_S"] = nb * f, nsv * f
        elif d[i] > 0.0:  # glass arriving on a bare node takes the composition of the node it came from (D-FLOW-6)
            up = states[i - 1] if i > 0 else s_
            nu = up["n_B"] + up["n_S"]
            x = up["n_B"] / nu if nu > 0 else up.get("_xB_last", 0.0)
            ntot = d[i] / (x * oxide.VB + (1.0 - x) * oxide.VS)
            s_["n_B"], s_["n_S"] = x * ntot, (1.0 - x) * ntot
        b = before[i]
        s_["_dn_dt"] = ((s_["n_B"] + s_["n_S"]) - (b["n_B"] + b["n_S"])) / dt
        dg = oxide.glass_thickness(s_)
        dc = s_["_dcrit"]
        if dg > dc:
            s_["glass_seen"] = True
        elif s_["x_r"] - p.x_r0 > p.x_ox_min and math.isnan(s_["t_glass_loss"]):
            # criterion (a) after the flow, interpolated inside the step as in LayerOxidation.advance
            dg0, dc0 = b.get("_dg_prev", 0.0), b.get("_dcrit", dc)
            t0 = s_["t"] - dt
            if dg0 > dc0 and (dg - dc) != (dg0 - dc0):
                s_["t_glass_loss"] = t0 + dt * min(max((0.0 - (dg0 - dc0)) / ((dg - dc) - (dg0 - dc0)), 0.0), 1.0)
            else:
                s_["t_glass_loss"] = s_["t"]
        s_["_dg_prev"] = dg
        states[i] = s_
