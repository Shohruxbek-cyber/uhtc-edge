"""Oxidation model and material properties."""
import math

import numpy as np
import pytest

from uhtc import data as dt, evaporation as ev, materials as mt, thermo as th
from uhtc.coupled import FlowState
from uhtc.oxidation import LayerOxidation, OxParams, evap_scalar


def flow(p_s=1e5, g=1.0):
    return FlowState(H_s=3e6, p_s=p_s, p_O2=0.2095 * p_s, h_D=0.0, Phi=1.0, f_nc=1.0, g_mass=g)


def test_material_tables():
    m = mt.composite(0.2)
    assert float(m.k(np.array(1673.15))) == pytest.approx(74.5, abs=1e-3)  # TPSX table value (1 K grid)
    assert float(m.k(np.array(3000.0))) == pytest.approx(69.0)  # held constant above 2273 K
    # JANAF rule of mixtures vs TPSX drop calorimetry 717 / 806 J/kg/K at 1273 / 2273 K (within ~3 %)
    assert float(m.cp(np.array(1273.15))) == pytest.approx(717.0, rel=0.04)
    assert float(m.cp(np.array(2273.15))) == pytest.approx(806.0, rel=0.04)
    assert mt.rho_rt(0.2) == pytest.approx(0.8 * 6090 + 0.2 * 3210)
    with pytest.raises(ValueError):
        mt.composite(0.25)


def test_kinetic_fits_reproduce_data():
    R = th.R
    kl = dt.KP_LOW_A * np.exp(-dt.KP_LOW_E / (R * dt.SHUGART_KP_LOW_T)) / dt.UM2H
    kh = dt.KP_HIGH_A * np.exp(-dt.KP_HIGH_E / (R * dt.SHUGART_KP_HIGH_T)) / dt.UM2H
    assert np.all(np.abs(np.log(kl / dt.SHUGART_KP_LOW)) < np.log(1.25))
    assert np.all(np.abs(np.log(kh / dt.SHUGART_KP_HIGH)) < np.log(1.6))  # data R^2 0.21-0.52
    assert dt.KP_HIGH_E / 1e3 == pytest.approx(326, abs=2)
    assert dt.TR_N == 45


def test_parabolic_growth_exact_at_constant_T():
    ox = LayerOxidation(OxParams(f_SiC=0.3))
    s = ox.init_state(300.0)
    f = flow()
    T = 1750.0
    for _ in range(100):
        s = ox.advance(s, 60.0, T, T, f)
    expect = math.sqrt(ox.p.x_r0 ** 2 + ox.kp(T, f.p_O2) * 6000.0)
    assert s["x_r"] == pytest.approx(expect, rel=1e-12)


def test_element_balance_B_and_Si():
    """Boron and silicon: converted + depleted = in glass + in pores + evaporated."""
    p = OxParams(f_SiC=0.2)
    ox = LayerOxidation(p)
    s = ox.init_state(300.0)
    f = flow(p_s=2e4, g=2.0)
    for T in np.r_[np.linspace(1500, 2400, 40), np.full(40, 2400.0)]:
        s = ox.advance(s, 2.0, T, T - 30.0, f)
    conv = s["x_r"] - p.x_r0
    B_made = conv * (1 - p.f_SiC) * ox.cZ  # mol B2O3 per m^2
    B_pores = B_made * (1 - ox.frac_out)
    B_out = s["n_B"] + s["m_B_evap"] / th.MOLAR_MASS["B2O3"]
    assert B_out + B_pores == pytest.approx(B_made, rel=1e-10)
    S_made = (conv * ox.frac_out + s["d_dep"] * (1 - p.f_SiO_lost)) * p.f_SiC * ox.cS
    S_out = s["n_S"] + s["m_S_evap"] / th.MOLAR_MASS["SiO2"]
    assert S_out == pytest.approx(S_made, rel=1e-10)
    assert s["m_B_evap"] > 0 and s["m_S_evap"] > 0


@pytest.mark.parametrize("p_eq", [1.0, 1e3, 3e4])
def test_scalar_evaporation_matches_module(p_eq):
    a = evap_scalar(p_eq, 2200.0, 5e4, 0.0696, 1.5, 1.0)
    b = float(ev.evaporation_flux(p_eq, 2200.0, 5e4, 0.0696, 1.5, 1.0))
    assert a == pytest.approx(b, rel=1e-12)


def test_transition_temperature_consistent():
    ox = LayerOxidation()
    for pO2 in (2e2, 2e3, 2e4):
        T = ox.T_transition(pO2)
        assert ox.log10_p_PA(T) == pytest.approx(math.log10(pO2 / 1e5), abs=1e-12)
    assert ox.T_transition(2e4) > ox.T_transition(2e3)  # higher O2 pressure keeps SiC passive to higher T


def test_failure_flags():
    ox = LayerOxidation(OxParams(T_limit=2000.0))
    s = ox.init_state(300.0)
    f = flow(p_s=5e3, g=3.0)
    for _ in range(200):
        s = ox.advance(s, 1.0, 2600.0, 2550.0, f)
    ev_ = ox.events(s)
    assert ev_["t_limit"] == pytest.approx(1.0)
    assert ev_["t_active"] == pytest.approx(1.0)  # p_O2 ~ 1 kPa: active at 2600 K
    assert not s["glass_seen"]  # evaporates as fast as it forms at 2600 K and 1 kPa O2
    assert not math.isnan(ev_["t_glass_loss"])  # ... which counts as criterion (a)


def test_f_SiO_lost_applies_to_depleted_layer_only():
    """f_SiO_lost (DECISIONS D-OX-7) removes only the SiC of the depleted layer from the glass budget."""
    import dataclasses
    f = flow(p_s=2e4, g=2.0)
    res = {}
    for lost in (0.0, 1.0):
        ox = LayerOxidation(OxParams(f_SiC=0.2, f_SiO_lost=lost))
        s = ox.init_state(300.0)
        for _ in range(50):
            s = ox.advance(s, 1.0, 1700.0, 1700.0, f)  # below T_dep: no depleted layer
        res[lost] = s["n_S"] + s["m_S_evap"] / th.MOLAR_MASS["SiO2"]
    assert res[0.0] == pytest.approx(res[1.0], rel=1e-12)


@pytest.mark.parametrize("nB,nS,v", [(0.0, 0.0, 0.0), (1e-3, 2e-3, 1e-7), (1e-6, 5e-2, 3e-8), (0.2, 1e-4, 0.0)])
def test_compiled_chem_matches_reference(nB, nS, v):
    ox = LayerOxidation(OxParams(f_SiC=0.2))
    s = ox.init_state(300.0)
    s.update(n_B=nB, n_S=nS, _v_last=v)
    for fl in (flow(2e4, 0.05), flow(2e5, 1.0)):
        fast, ref = ox.chem(s, fl), ox.chem_reference(s, fl)
        for T in (250.0, 900.0, 1650.5, 2321.7, 2989.9, 3200.0):
            a, b = fast(T), ref(T)
            assert a[0] == pytest.approx(b[0], rel=1e-12, abs=1e-9)
            assert a[1] == pytest.approx(b[1], rel=1e-9, abs=1e-9)


def test_wall_enthalpy_lookup_matches_interp():
    from uhtc.coupled import WallEnthalpy
    for cat in (True, False):
        w = WallEnthalpy(3e4, cat)
        for T in (100.0, 200.0, 305.0, 1234.5, 2000.0, 4999.0, 5000.0, 6000.0):
            h, dh = w(T)
            assert h == pytest.approx(float(np.interp(T, w.T, w.h)), rel=1e-13, abs=1e-6)
            assert dh == pytest.approx(float(np.interp(T, w.T, w.dh)), rel=1e-13, abs=1e-9)


def test_glass_criterion_permeation_and_fixed():
    fl = flow(2e4, 1.0)
    fixed = LayerOxidation(OxParams(f_SiC=0.2, glass_criterion="fixed", d_fail=0.3e-6))
    s = fixed.init_state(300.0)
    s.update(n_B=1e-3, n_S=1e-3)
    assert fixed.delta_crit(s, 2000.0, fl) == 0.3e-6
    ox = LayerOxidation(OxParams(f_SiC=0.2))
    silica, boria = dict(s, n_B=0.0, n_S=1e-3), dict(s, n_B=1e-3, n_S=0.0)
    # pure silica: the film limits the O2 supply down to sub-nanometre thickness; B2O3 is far more permeable
    assert ox.delta_crit(silica, 2300.0, fl) < 1e-9
    assert ox.delta_crit(boria, 1500.0, fl) > 1e3 * ox.delta_crit(silica, 1500.0, fl)
    # d_c = K p_O2 / J_gas: proportional to p_O2 / g, increasing with T
    fl2 = flow(4e4, 1.0)
    assert ox.delta_crit(silica, 2000.0, fl2) == pytest.approx(2.0 * ox.delta_crit(silica, 2000.0, fl), rel=1e-12)
    assert ox.delta_crit(silica, 2200.0, fl) > ox.delta_crit(silica, 2000.0, fl)
    lo, hi = (LayerOxidation(OxParams(f_SiC=0.2, perm_case=c)) for c in ("low", "high"))
    assert lo.delta_crit(silica, 2000.0, fl) < ox.delta_crit(silica, 2000.0, fl) < hi.delta_crit(silica, 2000.0, fl)


def test_kp_extrapolation_band_continuous_at_fit_end():
    base = LayerOxidation(OxParams(f_SiC=0.2))
    for kw in (dict(kp_E_hi=250e3), dict(kp_E_hi=354e3), dict(kp_mult_hi=0.5), dict(kp_mult_hi=3.0)):
        ox = LayerOxidation(OxParams(f_SiC=0.2, **kw))
        assert ox.kp(dt.KP_FIT_TMAX, 1e4) == pytest.approx(base.kp(dt.KP_FIT_TMAX, 1e4), rel=1e-12)  # unchanged up to the fit end
        assert ox.kp(2000.0, 1e4) == base.kp(2000.0, 1e4)
    lo, hi = LayerOxidation(OxParams(f_SiC=0.2, kp_E_hi=250e3)), LayerOxidation(OxParams(f_SiC=0.2, kp_E_hi=354e3))
    assert lo.kp(2400.0, 1e4) < base.kp(2400.0, 1e4) < hi.kp(2400.0, 1e4)


def test_pore_fill_silica_option():
    """Sensitivity option of the final review (R1-11): with the pores filled by silica only, all B2O3 and a smaller share
    of the SiO2 join the outer glass; the baseline keeps the production composition in the pores."""
    from uhtc.oxidation import LayerOxidation, OxParams
    base = LayerOxidation(OxParams(f_SiC=0.2))
    sil = LayerOxidation(OxParams(f_SiC=0.2, pore_fill="silica"))
    assert base.frac_out_B == base.frac_out_S == base.frac_out
    assert sil.frac_out_B == 1.0
    v_SiO2 = 0.2 * sil.cS * sil.VS
    assert sil.frac_out_S == pytest.approx((v_SiO2 - sil.pore) / v_SiO2)
    assert sil.frac_out_S < base.frac_out_S


def test_kd_arrhenius_option():
    """Sensitivity option of the external review (A6): k_d(T) = K_DEP kd_mult exp(-E/R (1/T - 1/T_ref)); kd_E = 0 is the
    baseline (no temperature dependence), and no depletion below T_dep in either case."""
    base = LayerOxidation(OxParams(f_SiC=0.2))
    arr = LayerOxidation(OxParams(f_SiC=0.2, kd_E=343e3, kd_mult=0.736))
    assert base.kd(2000.0) == base.kd(2300.0) == dt.K_DEP
    assert arr.kd(2073.15) == pytest.approx(0.736 * dt.K_DEP, rel=1e-12)
    assert arr.kd(1950.0) < arr.kd(2073.15) < arr.kd(2300.0)
    assert base.kd(dt.T_DEP_ONSET - 1.0) == arr.kd(dt.T_DEP_ONSET - 1.0) == 0.0
