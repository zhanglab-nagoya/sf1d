# Changelog

## 0.2.0

**Breaking: χs convention.** `chis` is now the Silberberg exchange parameter of the
standard Scheutjens–Fleer model (wall slot empty in the mixing term), as in Namics and
the literature. In 0.1.0 the wall slot was filled with bulk, so whenever χ ≠ 0 the same
`chis` gave weaker adsorption. Athermal results are unchanged; to reproduce a 0.1.0
result, subtract λ₁(Σ_y χ_xy φ_y^b − Σ_y χ_Sy φ_y^b) from each `chis`. The critical
adsorption energy is now χs_c = ln(6/5) − χ/6.

- Cross-check against Namics (`comparison/namics`, optional `tests/test_namics.py`):
  agreement to ~1e-11 in profiles, amounts, tensions and ΔF(D), full and restricted.
- Rescaled propagators: no overflow of the chain partition function for long, strongly
  adsorbed chains in restricted equilibrium.
- The Newton solver now fails cleanly (instead of raising a LinAlgError) when the
  propagators overflow at a trial point, so gap scans fall back to a cold start.
- Fixed a crash when neither `chis` nor `chis2` is given.
- Examples regenerated under the new convention.

## 0.1.0

First release: one- and two-wall SF lattice SCF for linear copolymers, train/loop/tail
decomposition, bridging, grand potential, gap scans in full and restricted equilibrium.
