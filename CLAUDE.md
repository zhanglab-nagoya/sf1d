# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`sf1d` is a small, tested 1D Scheutjens–Fleer self-consistent field (SCF) solver for
polymer adsorption from solution onto one flat wall or into the slit between two
parallel walls. It's a research/physics tool: two source files implement the whole
model, dependencies are only NumPy and SciPy (matplotlib optional for plots).

## Commands

```bash
pip install -e ".[plot,test]"   # editable install; plot = matplotlib, test = pytest
pytest                          # run the full test suite (6 Namics cross-check tests skip without NAMICS set)
pytest tests/test_sf1d.py::test_decomposition_and_bridges   # run a single test
bash examples/run_all.sh        # regenerate everything in examples/output/ (~1 min)
NAMICS=/path/to/namics/binary pytest   # also runs tests/test_namics.py against a built Namics binary
```

There is no linter/formatter configured in this repo.

CLI entry point (`sf1d = "sf1d.cli:main"` in pyproject.toml):

```bash
sf1d --seq A10B90 --chis A=2.0 --chi AS=0.8 --chi BS=0.4 --phib 1e-4 --M 80 --out run --plot
sf1d --seq A10B90 --chis A=2.0 --chi AS=0.8 --chi BS=0.4 --phib 1e-4 --scan 4 60 --mode both --out gap --plot
```

## Architecture

Everything lives under `src/sf1d/`:

- **`core.py`** — the entire model: `SF1D` class (lattice setup, propagators,
  Newton solver, free energy, adsorbed-chain analysis) plus the two-wall
  `gap_scan` driver and a `critical_chis_transfer_matrix` reference calculation.
  This is the only file to touch for anything physics/numerics-related.
- **`cli.py`** — argument parsing and glue: builds an `SF1D` (or calls `gap_scan`),
  writes `*_summary.json` / `*_profiles.csv` / `*_scan_{full,restricted}.csv`,
  and the optional matplotlib plots. No physics here.
- **`__init__.py`** — public API re-exports (`SF1D`, `gap_scan`, `parse_sequence`,
  `critical_chis_transfer_matrix`).

### Model structure (read `docs/theory.md` before changing `core.py`)

- All species (polymer segment types + solvent `"S"`) live on a simple cubic
  lattice, one layer per `z = 1..M`, solved via Bragg–Williams mean field within
  each layer. One wall uses a mirror boundary at `z = M+1`; two walls put a
  second wall there with its own `chis2`.
- The unknowns are segment potentials `u_x(z)`; `SF1D.residual` combines the
  incompressibility constraint (as `ln(sum phi)`, chosen to keep Newton steps
  bounded) with the requirement that all species see the same Lagrange field.
  `SF1D.solve` is a damped Newton method with a forward-difference Jacobian
  (batched across perturbations) and backtracking line search — **not**
  Picard/Anderson mixing, which is unstable here because `phi` depends on `u`
  roughly as `exp(N*u)` for long chains.
- Chain statistics use forward/backward propagators (`propagators`) built per
  segment from free-segment weights `G = exp(-u)`; propagators are rescaled at
  every step with the scale factors carried as logarithms (`lf`/`lb`), so long,
  strongly-adsorbed chains can't overflow. `densities` turns these into volume
  fractions without ever forming the raw partition function. `analyse()` reruns
  propagators with wall layers masked out to separate free vs. adsorbed chains,
  then derives the train/loop/tail decomposition and (two walls) bridging
  chains by inclusion–exclusion.
- `chis` is the **Silberberg adsorption exchange parameter** (matches Namics
  and the literature): the mixing term's wall-neighbour slot is conceptually
  empty, not bulk-filled. Internally `SF1D.effective_chis` converts it to an
  equivalent bulk-filled-slot value before use — a deliberate implementation
  detail, not the physical definition. This changed in 0.2.0 (breaking); don't
  assume a `chis` number predating that is still correct when `chi != 0`.
- `gap_scan` compresses a slit from wide to narrow gaps, warm-starting each
  `SF1D` solve from the previous gap's potentials via `remap` (falls back to a
  cold start if Newton fails from the warm start). Both full equilibrium
  (chains exchange with the bulk) and restricted equilibrium (`theta` fixed at
  its full-equilibrium value at the widest gap) can be computed in one scan.

### Tests (`tests/test_sf1d.py`, `tests/test_namics.py`)

`test_sf1d.py` validates the physics itself, not just code paths: exact
Langmuir-type results for a monomer, the critical adsorption energy
`chis_c = ln(6/5) - chi/6`, stationarity of the free energy via finite
differences w.r.t. `chis` and `chi`, agreement between full and restricted
equilibrium at large separation, overflow-free long-chain restricted
equilibrium, and qualitative reproduction of Scheutjens & Fleer (1985)
gap-scan trends. When touching `core.py`, these identities are the right way
to check correctness — prefer extending them over ad hoc scripts.

`test_namics.py` cross-checks `sf1d` against [Namics](https://github.com/leermakers/Namics)
(the Wageningen group's SF code) on 6 cases; skipped unless `NAMICS` points
to a built Namics binary. See `comparison/namics/README.md` for how to build
Namics and what the comparison caught (it's what found the 0.1.0 `chis` bug).

## Docs

- `README.md` — overview, install, quick start, example gallery, scope/limitations.
- `CHANGELOG.md` — version history; read before assuming a parameter's meaning
  hasn't changed across versions (`chis`'s definition changed in 0.2.0).
- `docs/theory.md` — the physics: lattice conventions, segment potentials, chain
  statistics, the Newton solver, free energy, train/loop/tail decomposition, limitations.
- `docs/usage.md` — full CLI/Python API reference, output file formats, practical
  notes on choosing segment size and interpreting `chi`/`chis`.
- `comparison/namics/README.md` — Namics cross-check: parameter mapping, results, caveats.
