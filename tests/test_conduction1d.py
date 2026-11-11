"""1D conduction solver against analytic solutions (constant properties)."""
import numpy as np
import pytest
from scipy.special import erf, erfc

from uhtc.conduction1d import Material, Solver1D

RHO, K, CP = 5000.0, 50.0, 600.0
ALPHA = K / (RHO * CP)


def const_mat():
    return Material("const", rho=lambda T: np.full_like(np.asarray(T, float), RHO),
                    k=lambda T: np.full_like(np.asarray(T, float), K),
                    cp=lambda T: np.full_like(np.asarray(T, float), CP), T_min=100, T_max=5000)


def run(s, t_end, dt0=1e-4, growth=1.15, dt_max=0.05):
    dt = dt0
    while s.t < t_end - 1e-12:
        dt = min(dt, t_end - s.t, dt_max)
        ok, _ = s.step(dt)
        assert ok
        dt *= growth
    return s


def test_constant_flux_semi_infinite():
    q, T0, t = 2.0e6, 300.0, 2.0
    s = Solver1D(const_mat(), thickness=0.05, T_init=T0, dx0=5e-6, growth=1.06,
                 front=lambda T, tt: (q, 0.0))
    run(s, t, dt0=1e-5, dt_max=0.01)
    x = s.depth_c
    exact = T0 + 2 * q / K * np.sqrt(ALPHA * t / np.pi) * np.exp(-x ** 2 / (4 * ALPHA * t)) \
        - q * x / K * erfc(x / (2 * np.sqrt(ALPHA * t)))
    err = np.abs(s.T - exact).max()
    rise = exact[0] - T0
    assert err / rise < 5e-3
    led = s.finish()
    assert abs(led.residual()[1]) < 1e-9


def test_step_temperature_semi_infinite():
    Ts, T0, t = 1800.0, 300.0, 1.0
    s = Solver1D(const_mat(), thickness=0.05, T_init=T0, dx0=5e-6, growth=1.06)
    Rh = s.half_cell_resistance()
    s.front = lambda T, tt: ((Ts - T) / Rh, -1.0 / Rh)
    run(s, t, dt0=1e-6, dt_max=0.01, growth=1.1)
    x = s.depth_c
    exact = Ts + (T0 - Ts) * erf(x / (2 * np.sqrt(ALPHA * t)))
    assert np.abs(s.T - exact).max() / (Ts - T0) < 5e-3


def test_adiabatic_slab_energy():
    # constant flux into a thin adiabatic slab: mean temperature rises as q t / (rho c L)
    q, L = 1e5, 0.002
    s = Solver1D(const_mat(), thickness=L, T_init=300.0, dx0=1e-5, growth=1.1,
                 front=lambda T, tt: (q, 0.0))
    run(s, 5.0)
    Tmean = np.sum(s.m * s.T) / np.sum(s.m)
    assert Tmean - 300.0 == pytest.approx(q * 5.0 / (RHO * CP * L), rel=1e-6)
