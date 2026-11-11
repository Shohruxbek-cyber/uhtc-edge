"""Coupled driver without oxidation: approach to radiative equilibrium and energy bookkeeping."""
import numpy as np
import pytest

from uhtc import coupled as cp, surface as sf
from uhtc.conduction1d import Material


def const_mat():
    f = lambda v: (lambda T: np.full_like(np.asarray(T, float), v))
    return Material("const", rho=f(5000.0), k=f(50.0), cp=f(600.0), T_min=100, T_max=6000)


def test_thin_slab_reaches_radiative_equilibrium():
    cond = cp.ConstantConditions(8.0, 30e3, 0.003, catalytic=True)
    flow = cond.at(0.0)
    hw = cp.WallEnthalpy(flow.p_s, True)
    conv = lambda T: (flow.Phi * (flow.H_s - hw(T)[0]), -flow.Phi * hw(T)[1])
    T_re = sf.radiative_equilibrium(conv, sf.constant_eps(0.8))
    r = cp.run_1d(cond, const_mat(), cp.NoOxide(0.8), t_end=400.0, depth=0.002, dx0=1e-5, dt_max=5.0)
    assert r.T_s[-1] == pytest.approx(T_re, abs=0.5)
    assert abs(r.energy_residual) < 1e-9
    assert np.all(np.diff(r.T_s) >= -1e-4)  # monotonic heating (round-off ~1e-6 K at equilibrium)


@pytest.mark.parametrize("M,z", [(8.0, 30e3), (10.0, 35e3), (6.0, 25e3)])
def test_noncatalytic_never_exceeds_catalytic(M, z):
    cond = cp.ConstantConditions(M, z, 0.003, catalysis="nc")
    fl = cond.at(0.0)
    eq, fr = cp.WallEnthalpy(fl.p_s, True), cp.WallEnthalpy(fl.p_s, False)
    q_nc = cp.convective_flux(fl, False, eq, fr)
    q_fc = cp.convective_flux(fl, True, eq, fr)
    for T in np.linspace(300.0, 0.98 * cond.state.T_e, 60):
        assert q_nc(T)[0] <= q_fc(T)[0] + 1e-9


def test_flank_mult_scales_flank_heating_only():
    from uhtc.conduction2d import NoseGeometry
    g = NoseGeometry(R_n=0.001, theta_c=np.radians(5.0), L_flank=0.01, axisymmetric=True)
    cond = cp.ConstantConditions(8.0, 30e3, 0.001, catalysis="fc")
    mesh = dict(n_arc=8, n_flank=6, n_eta=5, eta1=0.01)
    shapes = {}
    for fm in (1.0, 2.0):
        r = cp.run_2d(cond, g, const_mat(), cp.NoOxide(0.8), t_end=1e-3, mesh_kw=mesh, flank_mult=fm)
        shapes[fm], fl = r.extra["gshape"], r.extra["solver"].on_flank
    assert np.any(fl) and np.any(~fl)
    assert np.allclose(shapes[2.0][fl], 2.0 * shapes[1.0][fl])
    assert np.array_equal(shapes[2.0][~fl], shapes[1.0][~fl])


def test_fixed_glass_threshold_reproduces_draft_result():
    """glass_criterion='fixed', d_fail=0.1 um is the draft criterion (a): base case B1 (Mach 8, 35 km, 1 mm, fully
    catalytic) lost its glass after 710.862 s in the draft runs (runs/sensitivity/oat_B1.csv)."""
    from uhtc import materials as mt
    from uhtc.oxidation import LayerOxidation, OxParams
    ox = LayerOxidation(OxParams(f_SiC=0.2, glass_criterion="fixed", d_fail=0.1e-6))
    cond = cp.ConstantConditions(8.0, 35e3, 0.001, catalysis="fc")
    r = cp.run_1d(cond, mt.composite(0.2), ox, 1800.0, depth=0.02,
                  stop_on=lambda e, T, t: "first" if any(v == v for v in e.values()) else None)
    assert r.events["t_glass_loss"] == pytest.approx(710.862, abs=0.01)


def test_ground_test_edge_state_beyond_6000K():
    """Arc-jet steps above the 6000 K table end get their own edge state (not the clipped 5998 K one); above
    T_EDGE_MAX the state is held there and flagged."""
    from uhtc import thermo as th
    p = 0.042 * 101325.0
    cond = cp.PrescribedConditions([(0.0, 15e6, p), (10.0, 30e6, p), (20.0, 45e6, p)], R_n=0.01)
    T15, T30, T45 = cond.edge_T
    assert T15 < 6000.0 < T30 < T45 == pytest.approx(cp.T_EDGE_MAX - 1.0)
    assert cond.edge_held == [False, False, True]
    _, _, hD_clip = th.mixture_props(th.equilibrium_air(5998.0, p), 5998.0)
    assert cond.flows[1].h_D > 1.15 * hD_clip


def test_local_pressure_option_scales_oxide_pressures():
    """External review A4: with local_pressure=True each node's oxide model sees p/p_s = (1 - eps) cos^2(phi) + eps
    times the stagnation wall and O2 pressures (eps = p_inf/p_s); the default keeps the stagnation values everywhere."""
    from uhtc import atmosphere as atm, materials as mt
    from uhtc.conduction2d import NoseGeometry
    from uhtc.oxidation import LayerOxidation, OxParams
    g = NoseGeometry(R_n=0.001, theta_c=np.radians(5.0), L_flank=0.01, axisymmetric=True)
    cond = cp.ConstantConditions(8.0, 30e3, 0.001, catalysis="fc")
    fl0 = cond.at(0.0)
    eps = atm.state(30e3)[1] / fl0.p_s
    mesh = dict(n_arc=8, n_flank=6, n_eta=5, eta1=0.01)
    seen = {}
    for lp in (False, True):
        ox = LayerOxidation(OxParams(f_SiC=0.2))
        calls = []
        orig = ox.advance
        ox.advance = lambda s, dt_, Ts, Tf, fl, _o=orig, _c=calls: (_c.append((fl.p_s, fl.p_O2)), _o(s, dt_, Ts, Tf, fl))[1]
        r = cp.run_2d(cond, g, mt.composite(0.2), ox, t_end=1e-3, mesh_kw=mesh, local_pressure=lp)
        ns = len(r.extra["solver"].surf)
        seen[lp] = np.array(calls[-ns:])
        phi = r.extra["solver"].phi
    assert np.allclose(seen[False][:, 0], fl0.p_s) and np.allclose(seen[False][:, 1], fl0.p_O2)
    pr = (1.0 - eps) * np.cos(phi) ** 2 + eps
    assert pr[0] == pytest.approx(1.0)
    assert np.allclose(seen[True][:, 0], fl0.p_s * pr) and np.allclose(seen[True][:, 1], fl0.p_O2 * pr)
    assert seen[True][-1, 1] < 0.02 * fl0.p_O2  # flank: under 2 % of the stagnation value
