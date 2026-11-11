"""Version 3, Task 2: fast unit tests of the thermoelastic module (the full verification is check V11)."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "validation" / "V11_thermoelastic"))
import run_v11 as v  # noqa: E402
from uhtc.thermoelastic import principal_max, solve, solve_quadratic  # noqa: E402


def _patch(solver):
    nodes, tri = v.grid(4, 3, 1e-3, 0.7e-3, jitter=0.3)
    dT = 800.0
    k = int(np.argmin(nodes[:, 0] + nodes[:, 1]))
    out = solver(nodes, tri, np.full(len(nodes), 300.0 + dT), v.P, bc="free", z_fix_node=k)
    return out[1], v.E * v.AL * dT


def test_patch_linear_free_expansion_is_stress_free():
    sig, scale = _patch(solve)
    assert np.abs(sig).max() < 1e-9 * scale


def test_patch_quadratic_free_expansion_is_stress_free():
    sig, scale = _patch(solve_quadratic)
    assert np.abs(sig).max() < 1e-9 * scale


def test_uniform_temperature_at_reference_gives_zero_stress():
    nodes, tri = v.grid(3, 3, 1e-3, 1e-3)
    _, sig, _, _ = solve_quadratic(nodes, tri, np.full(len(nodes), 300.0), v.P, bc="free", z_fix_node=0)
    assert np.abs(sig).max() < 1e-12 * v.E  # round-off only (a few micro-pascals against E of order 1e11 Pa)


def test_principal_max_of_a_uniaxial_state():
    sig = np.array([[2.0, -1.0, 0.5, 0.0]])  # one stress state in the module's component order
    s1 = principal_max(sig)
    assert s1.shape == (1,)
    assert s1[0] >= 2.0 - 1e-12
