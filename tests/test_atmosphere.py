"""US Standard Atmosphere 1976: check against Table I of the published document (geometric altitude,
metric units; literature/pdf/USSA1976.pdf, PDF pages 76, 78, 80 = doc pages 61, 63, 65)."""
import numpy as np
import pytest

from uhtc import atmosphere as atm

# (z m, T K, p Pa, rho kg/m^3) read from Table I (p in mb x 100)
TABLE_I = [
    (20000.0, 216.650, 5529.3, 8.8910e-2),
    (30000.0, 226.509, 1197.0, 1.8410e-2),
    (40000.0, 250.350, 287.14, 3.9957e-3),
]


@pytest.mark.parametrize("z,T,p,rho", TABLE_I)
def test_table_I(z, T, p, rho):
    Tm, pm, rhom = atm.state(z)
    assert Tm == pytest.approx(T, abs=0.002)
    assert pm == pytest.approx(p, rel=2e-4)  # table gives 5 significant figures
    assert rhom == pytest.approx(rho, rel=2e-4)


def test_sea_level_and_M0():
    T, p, rho = atm.state(0.0)
    assert T == pytest.approx(288.15)
    assert p == pytest.approx(101325.0)
    # doc p. 9: M0 = 28.9644 kg/kmol from Table 3 (our sum over Table 3 gives 28.96451)
    assert atm.M0_TABLE3 == pytest.approx(atm.M0, abs=2e-4)
    assert rho == pytest.approx(1.2250, rel=1e-4)  # Table II (geopotential), H = 0 row, PDF p. 67


def test_geopotential_row():
    # Table I: Z = 30000 m <-> H = 29859 m'
    assert atm.geopotential(30000.0) == pytest.approx(29859.0, abs=0.6)


def test_array_and_range():
    z = np.array([20e3, 30e3, 40e3])
    T, p, rho = atm.state(z)
    assert np.all(np.diff(p) < 0)
    with pytest.raises(ValueError):
        atm.state(90e3)
