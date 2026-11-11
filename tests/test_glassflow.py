"""Unit tests of the glass-flow module (version 3, Task 1): conservation, positivity, integrator switch, drivers, viscosity."""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

from uhtc import glassflow as gf
from uhtc.coupled import FlowState

ROOT = Path(__file__).resolve().parents[1]


def _surface(n=60, axisymmetric=True):
    s = np.linspace(0.0, 5e-3, n)
    r = 1e-3 * np.sin(np.minimum(s / 1e-3, np.pi / 2)) + 1e-6 if axisymmetric else np.ones(n)
    ds = np.diff(s)
    w = np.zeros(n)
    w[:-1] += ds * (2 * r[:-1] + r[1:]) / 6
    w[1:] += ds * (r[:-1] + 2 * r[1:]) / 6
    return gf.Surface(s, r, w, np.minimum(s / 1e-3, np.pi / 2), np.r_[0.5 * (r[:-1] + r[1:]), r[-1]], 1e-3, axisymmetric)


@pytest.mark.parametrize("max_explicit", [10_000, 0])
def test_volume_conservation_and_positivity(max_explicit):
    rng = np.random.default_rng(1)
    surf = _surface()
    d0 = rng.uniform(0, 60e-6, len(surf.s))
    d0[5] = 0.0
    tau, G = rng.uniform(0, 3000, len(surf.s)), rng.uniform(0, 5e6, len(surf.s))
    mu = rng.uniform(1e2, 1e5, len(surf.s))
    d, shed, nsub = gf.film_step(d0, surf, tau, G, mu, 50.0, gf.FlowParams(max_explicit=max_explicit))
    assert (nsub < 0) == (max_explicit == 0)
    assert np.all(d >= 0.0)
    before, after = float(np.sum(surf.w * d0)), float(np.sum(surf.w * d)) + shed
    assert abs(after - before) <= 1e-12 * before
    assert shed >= 0.0


def test_no_glass_no_flow():
    surf = _surface()
    d, shed, _ = gf.film_step(np.zeros(len(surf.s)), surf, np.full(len(surf.s), 1e3), np.zeros(len(surf.s)),
                              np.full(len(surf.s), 1e3), 10.0, gf.FlowParams())
    assert not d.any() and shed == 0.0


def test_drivers_vanish_at_the_stagnation_point_and_point_outward():
    surf = _surface()
    flow = FlowState(10e6, 50e3, 10e3, 1e6, 2.0, 0.5, 2.0, 9.0, 30e3)
    tau_node_g = np.cos(surf.phi)
    tau, G = gf.drivers(surf, flow, tau_node_g, rho_e=0.02, p_inf=1.2e3, params=gf.FlowParams())
    assert np.all(tau >= 0) and np.all(G >= 0)
    assert tau[0] < tau[3]  # shear grows away from the stagnation point (u_e = 0 there)
    _, G_off = gf.drivers(surf, flow, tau_node_g, 0.02, 1.2e3, gf.FlowParams(pressure=False))
    assert not G_off.any()


def test_viscosity_matches_the_supplement_estimate():
    sys.path.insert(0, str(ROOT / "audit" / "tools"))
    sys.path.insert(0, str(ROOT / "runs"))
    import glass_flow_estimate as est
    for T in (2000.0, 2300.0):
        for law in ("LA", "F320"):
            assert gf.viscosity(T, law) == pytest.approx(est.mu(T, law), rel=1e-12)


def test_shear_fan_short():
    """Shortened V10 check 1 (shear-driven expansion fan) on a coarse grid: within 3 %."""
    n, L, mu, tau, d0 = 401, 10e-3, 1e4, 1000.0, 20e-6
    s = np.linspace(0, L, n)
    w = np.full(n, L / (n - 1))
    w[0] = w[-1] = w[1] / 2
    surf = gf.Surface(s, np.ones(n), w, np.zeros(n), np.ones(n), 1.0, False)
    i = n // 4
    t = 4 * mu * s[i] / (tau * d0)
    d, _, _ = gf.film_step(np.full(n, d0), surf, np.full(n, tau), np.zeros(n), np.full(n, mu), t, gf.FlowParams())
    assert d[i] == pytest.approx(mu * s[i] / (tau * t), rel=0.03)


def test_flow_off_is_the_default_of_run_2d():
    import inspect
    from uhtc import coupled
    assert inspect.signature(coupled.run_2d).parameters["glass_flow"].default is None
