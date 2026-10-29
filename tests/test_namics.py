"""Optional cross-check against Namics.  Runs only if the environment variable
NAMICS points to a namics binary (github.com/leermakers/Namics)."""

import os
import pathlib
import sys

import pytest

NAMICS = os.environ.get("NAMICS")
pytestmark = pytest.mark.skipif(not NAMICS, reason="set NAMICS=/path/to/namics to run")

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "comparison" / "namics"))


@pytest.mark.parametrize("case", __import__("compare").all_cases(), ids=lambda c: c["tag"])
def test_agrees_with_namics(case):
    import compare

    r = compare.compare(NAMICS, **case)
    assert r["max_dphi"] < 1e-8
    assert r["theta_sf1d"] == pytest.approx(r["theta_namics"], rel=1e-8)
    assert r["F_sf1d"] == pytest.approx(r["F_namics"], rel=1e-7, abs=1e-9)
