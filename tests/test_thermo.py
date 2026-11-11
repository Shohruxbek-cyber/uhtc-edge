"""JANAF parsing and thermochemistry checks."""
import math
import re

import numpy as np
import pytest

from uhtc import thermo as th


@pytest.mark.parametrize("sp", sorted(th.FILES))
def test_dfG_consistent_with_logKf(sp):
    """Every parsed delta-f G must equal -R T ln10 log Kf (catches rows whose minus signs NIST dropped:
    literature/janaf/README.md lists 8 such rows; the tab-based parser skips them)."""
    t = th.table(sp)
    path = th.JANAF_DIR / f"{th.FILES[sp]}.txt"
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines()[2:]:
        f = re.split(r"\s+", line.strip())
        if len(f) == 8:
            try:
                rows.append([float(x) for x in f])
            except ValueError:
                pass
    checked = 0
    for T, *_, dfG, logKf in rows:
        if T < 200 or T not in t.T or logKf == 0.0:
            continue
        i = int(np.where(t.T == T)[0][0])
        expect = -th.R * T * math.log(10) * logKf
        assert t.dfG[i] == pytest.approx(expect, rel=2e-3, abs=50.0), (sp, T)
        checked += 1
    if sp not in ("N2", "O2"):
        assert checked > 5


def test_quirk_rows_are_skipped():
    t = th.table("B2O3(l)")
    assert 800.0 not in t.T and 700.0 in t.T and 900.0 in t.T


def test_air_properties():
    assert th.M_AIR == pytest.approx(0.028965, rel=2e-4)
    assert th.cp_air_frozen(300.0) == pytest.approx(1005.0, rel=0.01)  # dry air near room temperature
    x = th.equilibrium_air(300.0, 1e5)
    assert x["O"] < 1e-30 and x["O2"] == pytest.approx(0.2095, rel=1e-3)


def test_tables_strict_outside_extension():
    """Above the table end the look-ups raise unless they sit inside thermo.extended_tables (ground tests only)."""
    with pytest.raises(ValueError):
        th.table("N").H(6500.0)
    with th.extended_tables(9000.0):
        assert th.table("N").H(6000.0) == pytest.approx(th.table("N").H(6000.0 - 1e-9), rel=1e-12)
        assert th.table("N").H(6500.0) > th.table("N").H(6000.0)
    with pytest.raises(ValueError):
        th.table("N").H(6500.0)


def test_table_extension_backtest(monkeypatch):
    """Tables cut at 5000 K and extended to 6000 K reproduce equilibrium air at 6000 K: enthalpy and h_D within 1 %,
    the non-catalytic factor within 0.002 (the quantity the extension feeds; review/FINDINGS.csv, V8)."""
    import copy

    from uhtc import heating as ht
    full = {sp: th.table(sp) for sp in ("N2", "O2", "N", "O", "NO")}
    cuts = {}
    for sp, t in full.items():
        c, k = copy.copy(t), t.T <= 5000.0
        c.T, c.Cp, c.S, c.dH, c.dfH, c.dfG = t.T[k], t.Cp[k], t.S[k], t.dH[k], t.dfH[k], t.dfG[k]
        cuts[sp] = c
    for p in (4.26e3, 18e3):
        _, h0, hD0 = th.mixture_props(th.equilibrium_air(6000.0, p), 6000.0)
        monkeypatch.setattr(th, "table", lambda sp: cuts[sp])
        with th.extended_tables(9000.0):
            _, h1, hD1 = th.mixture_props(th.equilibrium_air(6000.0, p), 6000.0)
        monkeypatch.setattr(th, "table", lambda sp: full[sp])
        assert h1 == pytest.approx(h0, rel=0.01)
        assert hD1 == pytest.approx(hD0, rel=0.01)
        assert abs(ht.noncatalytic_factor(h1, hD1) - ht.noncatalytic_factor(h0, hD0)) < 0.002
