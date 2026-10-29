"""Checks against exact limits and thermodynamic identities."""

import numpy as np
import pytest
from scipy.optimize import brentq

from sf1d import SF1D, critical_chis_transfer_matrix, gap_scan, parse_sequence

AB = dict(chi={"AS": 0.8, "BS": 0.4}, phib=1e-4)


def test_parse_sequence():
    assert "".join(parse_sequence("A2B3")) == "AABBB"
    with pytest.raises(ValueError):
        parse_sequence("S10")  # S is reserved for the solvent


def test_defaults_without_surface_parameters():
    m = SF1D("A10", phib=1e-3, M=20, walls=2).solve()  # no chis/chis2 given: inert walls
    assert np.all(m.chis == 0) and np.all(m.chis2 == 0)


def test_monomer_adsorption_is_langmuir():
    phib, chis = 1e-3, 1.5
    m = SF1D("A1", chis={"A": chis}, phib=phib, M=10).solve()
    exact = phib * np.exp(chis) / (1 - phib + phib * np.exp(chis))
    assert m.phi[0, 0] == pytest.approx(exact, rel=1e-9)


def test_monomer_interfacial_tension():
    phib, chis = 1e-3, 1.5
    m = SF1D("A1", chis={"A": chis}, phib=phib, M=10).solve()
    assert m.free_energy() == pytest.approx(-np.log(1 - phib + phib * np.exp(chis)), rel=1e-9)


def test_critical_adsorption_energy():
    assert critical_chis_transfer_matrix() == pytest.approx(np.log(6 / 5), abs=1e-4)


def test_zero_excess_point_approaches_critical_value():
    def zero_excess(N):
        exc = lambda c: SF1D(f"A{N}", chis={"A": c}, phib=1e-6, M=120).solve().analyse()[0]["theta_excess"]
        return brentq(exc, 0.0, 0.6, xtol=1e-5)

    c50, c200 = zero_excess(50), zero_excess(200)
    assert c50 < c200 < np.log(6 / 5)
    assert np.log(6 / 5) - c200 < 0.003


def test_bulk_recovery_and_incompressibility():
    m = SF1D("A200", chis={"A": 1.0}, chi={"AS": 0.5}, phib=1e-3, M=80).solve()
    assert m.phi[0, -1] == pytest.approx(1e-3, rel=1e-4)
    assert np.max(np.abs(m.phi.sum(axis=0) - 1)) < 1e-9


@pytest.mark.parametrize("walls", [1, 2])
def test_free_energy_derivative_wrt_surface_affinity(walls):
    """Stationarity: dOmega/dchis_A = -phi_A in the wall layer(s)."""
    om = lambda c: SF1D("A10B90", chis={"A": c}, M=40, walls=walls, **AB).solve()
    h = 1e-4
    fd = (om(2 + h).free_energy() - om(2 - h).free_energy()) / (2 * h)
    m = om(2.0)
    exact = -(m.phi[0, 0] + (m.phi[0, -1] if walls == 2 else 0.0))
    assert fd == pytest.approx(exact, rel=1e-6)


def test_free_energy_derivative_wrt_chi():
    """Stationarity: the total derivative dOmega/dchi_AS (re-solving) equals the
    explicit derivative at frozen fields (Hellmann-Feynman)."""
    om = lambda c: SF1D("A10B90", chis={"A": 2.0}, chi={"AS": c, "BS": 0.4}, phib=1e-4, M=40).solve()
    h = 1e-4
    total = (om(0.8 + h).free_energy() - om(0.8 - h).free_energy()) / (2 * h)

    m = om(0.8)
    iA, iS = m.idx["A"], m.idx["S"]

    def frozen(c):
        m.chi[iA, iS] = m.chi[iS, iA] = c
        return m.free_energy()

    explicit = (frozen(0.8 + h) - frozen(0.8 - h)) / (2 * h)
    assert total == pytest.approx(explicit, rel=1e-6)


def test_critical_adsorption_shifts_with_solvent_quality():
    """Scheutjens-Fleer convention: chis_c = ln(6/5) - lambda1*chi for long dilute chains."""
    exc = lambda c: SF1D("A200", chis={"A": c}, chi={"AS": 0.5}, phib=1e-6, M=120).solve().analyse()[0]["theta_excess"]
    assert brentq(exc, -0.2, 0.6, xtol=1e-5) == pytest.approx(np.log(6 / 5) - 0.5 / 6, abs=0.003)


def test_distant_walls_and_restricted_equilibrium():
    gamma = SF1D("A10B90", chis={"A": 2.0}, M=80, **AB).solve().free_energy()
    full = SF1D("A10B90", chis={"A": 2.0}, M=120, walls=2, **AB).solve()
    assert full.free_energy() == pytest.approx(2 * gamma, rel=1e-9)
    restricted = SF1D("A10B90", chis={"A": 2.0}, M=120, walls=2,
                      theta=full.analyse()[0]["theta_total"], **AB).solve()
    assert restricted.free_energy() == pytest.approx(full.free_energy(), rel=1e-9)
    assert np.max(np.abs(restricted.phi - full.phi)) < 1e-10


def test_long_chain_restricted_equilibrium_without_overflow():
    """Rescaled propagators: N=1000, restricted at the full-equilibrium amount."""
    kw = dict(chis={"A": 1.0}, chi={"AS": 0.5}, phib=1e-4, M=30, walls=2)
    with np.errstate(over="raise"):
        full = SF1D("A1000", **kw).solve()
        restricted = SF1D("A1000", theta=full.analyse()[0]["theta_total"], **kw).solve()
    assert restricted.free_energy() == pytest.approx(full.free_energy(), rel=1e-9)


def test_decomposition_and_bridges():
    m = SF1D("A10B90", chis={"A": 2.0}, M=12, walls=2, **AB).solve()
    s, p = m.analyse()
    assert np.all(p["loops"] > -1e-12) and np.all(p["tails"] > -1e-12)
    assert s["fraction_trains"] + s["fraction_loops"] + s["fraction_tails"] == pytest.approx(1.0)
    assert 0.0 < s["bridging_fraction_of_adsorbed_chains"] < 1.0


def test_homopolymer_gap_scan_matches_scheutjens_fleer_trends():
    """Full equilibrium: attraction.  Restricted: bridging minimum, then repulsion."""
    res = gap_scan("A200", {"AS": 0.5}, {"A": 1.0}, None, 1e-4, np.arange(3, 31), verbose=False)
    full, restr = res["full"], res["restricted"]
    assert np.all(full["dF"] <= 1e-8)
    assert restr["dF"].min() < 0 and restr["dF"][-1] > 0
