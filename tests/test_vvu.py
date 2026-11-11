"""Tests of the post hoc validation-uncertainty analysis (validation/VVU): the per-quantity statistics reproduce the text, and
the baseline sample of the runner reproduces the stored check results (the reproducibility gate of VVU_plan.md, Section 7)."""
import importlib.util
import json
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
VVU = ROOT / "validation" / "VVU"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, VVU / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_per_quantity_reproduces_v6_text():
    pq = _load("per_quantity")
    _, summary = pq.v6()
    assert summary["oxide_type"] == "15 of 18"
    assert summary["gm_oxide_type"] == 1.09
    assert summary["depleted"] == "9 of 19"


def test_gate_record_passed():
    gate = json.loads((VVU / "vvu_gate.json").read_text())
    assert gate["passed"] and not gate["problems"]
    assert gate["r7"] == pytest.approx(2.0, abs=0.01)


def test_gate_v9_baseline_reproduces_stored_values():
    """The fast part of the gate: the non-catalytic V9 baseline through the VVU wrapper equals the stored V9 results."""
    pytest.importorskip("uhtc")
    rv = _load("run_vvu")
    base = rv.designs("V9", 1)["base"][0]
    out = rv.run_v9_case("nc", {**base, "gamma_cu": None})
    stored = json.loads((ROOT / "validation/V9_plasmatron/v9_results.json").read_text())
    for s in ("1.8", "1.5", "1.6"):
        ref = stored[f"nc_{s}"]
        assert out[s]["T_low"] == pytest.approx(ref["T_low"], rel=1e-9)
        assert out[s]["T_ref"] == pytest.approx(ref["T_ref"], rel=1e-9)
        assert math.isnan(ref["t_glass"]) == math.isnan(out[s]["t_glass"])


def test_calorimeter_ratio_limits():
    """A fully catalytic calorimeter reads the fully catalytic flux; a weakly catalytic one reads less."""
    rv = _load("run_vvu")
    r1, f_nc, phi1 = rv.calorimeter_ratio(18.6e6, 1.0e4 + 75.0, 202e4, 1.0, 350.0)
    r01, _, phi01 = rv.calorimeter_ratio(18.6e6, 1.0e4 + 75.0, 202e4, 0.01, 350.0)
    assert 0.0 < f_nc < r01 < r1 <= 1.0
    assert phi1 > phi01


def test_s10_tables_match_results():
    """Every cell of the two tables of Supplementary Section S10 is rebuilt from vvu_results.json and per_quantity.json."""
    import subprocess
    import sys
    if not (VVU / "check_tables.py").exists() or not (ROOT / "paper" / "main.tex").exists():
        pytest.skip("the manuscript is not in this copy of the repository")
    r = subprocess.run([sys.executable, str(VVU / "check_tables.py")], capture_output=True, text=True, encoding="utf-8")
    assert r.returncode == 0, r.stdout


def test_censored_glass_loss_is_a_lower_bound():
    """A V9 glass-loss time kept to the end of the run is reported as E > ln(time left / D), never as a value."""
    res = json.loads((VVU / "vvu_results.json").read_text())["results"]
    r = next(x for x in res if x["quantity"] == "sample 1.5: glass loss after the switch" and x["bound"] == "fc")
    assert r["E"] != r["E"] and r["E_lower"] == pytest.approx(math.log(123.0 / 68.0))
    assert "cannot be applied" in r["label"]
