"""2D P1 finite-element conduction: geometry weights, energy conservation, 1D short-time limit."""
import numpy as np
import pytest

from uhtc.conduction1d import Material
from uhtc.conduction2d import NoseGeometry, Solver2D

RHO, K, CP = 5000.0, 50.0, 600.0


def const_mat():
    return Material("const", rho=lambda T: np.full_like(np.asarray(T, float), RHO),
                    k=lambda T: np.full_like(np.asarray(T, float), K),
                    cp=lambda T: np.full_like(np.asarray(T, float), CP), T_min=100, T_max=5000)


def exact_sphere_cone(Rn, th, Lf):
    """Volume and wetted area of the half-body as swept (axisymmetric, without 2 pi) - exact."""
    phit = np.pi / 2 - th
    # arc area: int y ds = Rn^2 (1 - cos phit)
    A_arc = Rn ** 2 * (1 - np.cos(phit))
    yt = Rn * np.sin(phit)
    A_fl = Lf * (yt + 0.5 * Lf * np.sin(th))
    return A_arc + A_fl


@pytest.mark.parametrize("axi", [True, False])
def test_surface_weights_and_energy(axi):
    g = NoseGeometry(R_n=0.003, theta_c=np.radians(8), L_flank=0.02, axisymmetric=axi)
    s = Solver2D(g, const_mat(), 300.0, mesh_kw=dict(n_arc=30, n_flank=30, n_eta=25))
    if axi:
        assert s.w_surf.sum() == pytest.approx(exact_sphere_cone(0.003, np.radians(8), 0.02), rel=1e-3)
    else:
        L_arc = 0.003 * (np.pi / 2 - np.radians(8))
        assert s.w_surf.sum() == pytest.approx(L_arc + 0.02, rel=1e-4)  # chords vs arc
    q = 1e5
    flux = lambda T: (np.full(len(T), q), np.zeros(len(T)))
    for _ in range(20):
        ok, _ = s.step(0.05, flux)
        assert ok
    assert abs(s.energy_residual()) < 1e-8
    Tmean = np.sum(s.m * s.T) / np.sum(s.m)
    assert Tmean - 300 == pytest.approx(q * s.w_surf.sum() * 1.0 / (RHO * CP * s.vol.sum()), rel=1e-8)


def test_short_time_matches_1d():
    # uniform flux; at t where the penetration depth << R_n the stagnation node follows the 1D solution
    Rn = 0.01
    g = NoseGeometry(R_n=Rn, theta_c=np.radians(8), L_flank=0.02, axisymmetric=True)
    s = Solver2D(g, const_mat(), 300.0, mesh_kw=dict(n_arc=40, n_flank=20, n_eta=60, eta1=0.0005))
    q, t_end = 2e6, 0.02
    flux = lambda T: (np.full(len(T), q), np.zeros(len(T)))
    dt = 1e-6
    while s.t < t_end - 1e-15:
        dt = min(dt * 1.1, t_end - s.t)
        assert s.step(dt, flux)[0]
    exact = 300 + 2 * q * np.sqrt(t_end / (np.pi * K * RHO * CP))
    Ts = s.T[s.surf[0]]
    assert Ts == pytest.approx(exact, rel=0.02)


def test_lees_distribution_normalised_and_reproduces_crabtree_charts():
    """Local-similarity distribution: q/q_s -> 1 at the stagnation point, and on a 5-15 deg hemisphere-cone (perfect gas,
    Newtonian pressure) it reproduces the Crabtree et al. 1965 Fig. 12 charts (DIGITIZED) within 12 %."""
    import csv
    import math
    from pathlib import Path
    from uhtc.conduction2d import NoseGeometry, lees_distribution
    g = 1.4

    def pitot(M):
        return ((g + 1) ** 2 * M ** 2 / (4 * g * M ** 2 - 2 * (g - 1))) ** (g / (g - 1)) * (1 - g + 2 * g * M ** 2) / (g + 1)
    f = Path(__file__).resolve().parents[1] / "literature" / "digitized" / "Crabtree1965_fig12.csv"
    rows = [r for r in csv.DictReader(l for l in open(f, encoding="utf-8") if not l.startswith("#"))]
    worst = 0.0
    for th in (5, 10, 15):
        geo = NoseGeometry(R_n=1.0, theta_c=math.radians(th), L_flank=10.0, axisymmetric=True)
        for M in (5, 8, 10, 14):
            pts = [(float(r["s_over_R"]), float(r["q_over_qs"])) for r in rows if r["theta_c_deg"] == str(th)
                   and r["Mach"] == str(M)]
            V = float(M)
            q = lees_distribution(geo, np.array([p[0] for p in pts]), M, V, 0.5 * V ** 2 * (1 + 2 / ((g - 1) * M ** 2)),
                                  pitot(M), 1.0)
            worst = max(worst, max(abs(a / b[1] - 1) for a, b in zip(q, pts)))
            assert lees_distribution(geo, np.array([1e-4]), M, V, 0.5 * V ** 2, pitot(M), 1.0)[0] == pytest.approx(1, abs=2e-3)
    assert worst < 0.12
