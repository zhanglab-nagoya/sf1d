# Plan: deriving `chi` and `chis` from MD

For PSMA-*stat*-PMMA and PSMA-*block*-PMMA (PSMA = `ST` = stearyl
methacrylate, the "buoy"; PMMA = `MMA` = methyl methacrylate, the "anchor")
adsorbing from hexadecane onto Fe2O3 (`substrate/fe2o3.data`). Flat geometry
throughout — curved/probe geometry is explicitly out of scope for this plan.

## Step 0 — segment size and chain mapping

Pick the lattice segment volume `a^3` = hexadecane molecular volume, from a
pure hexadecane NPT box (`a^3 = M/(rho * N_A)`). Work out how many lattice
segments one ST or MMA repeat unit maps to by matching molar volumes (ST's
C18 tail is much bigger than one hexadecane molecule, so likely >1 segment
per repeat unit). This ratio sets `N` in `sf1d`'s `--seq` for a given real
degree of polymerization.

## Step 1 — `chi_AS`, `chi_BS` (segment-solvent)

- **Primary, cheap:** Hildebrand solubility parameters from short pure-liquid
  NPT runs — hexadecane, and melts of short ST/MMA oligomer fragments (built
  from `monomers/`). `chi_xy = V_ref (delta_x - delta_y)^2 / RT` with
  `V_ref = a^3` from step 0.
- **Refinement:** alchemical transfer free energy (TI/BAR) of a repeat-unit
  fragment decoupled in hexadecane vs. in its own melt — catches specific
  ester-alkane effects the solubility-parameter route averages away.
- **Independent check:** single-chain R_g of the real sequence in explicit
  hexadecane, matched against a bulk-only (`chis=0`, large `M`) `sf1d` run
  tuned over `chi`.

## Step 2 — `chis_A`, `chis_B` (segment-Fe2O3, relative to solvent)

- **Primary:** dilute slab MD of ST-analog/MMA-analog fragments (and
  hexadecane separately) near the Fe2O3 surface. Read `chis_x` off the
  contact-vs-bulk density ratio, correcting for the solvent's own layering:

      chis_x = ln[phi_contact,x / phi_bulk,x] - ln[phi_contact,S / phi_bulk,S]

  This is exactly the dilute-limit relation `sf1d` itself assumes
  (`phi(wall)/phi_bulk ~= exp(chis)`, see
  `tests/test_sf1d.py::test_monomer_adsorption_is_langmuir`), and — since
  `sf1d` 0.2.0 uses the standard Silberberg convention — matches Namics and
  the literature directly, with no translation needed.
- **Cross-check:** umbrella-sampling PMF along the surface normal gives
  `chis_x = [dG_ads(x) - dG_ads(S)] / kT` directly.
- **Literature sanity check:** published SF-theory `chis` values for
  PMMA/alkyl-methacrylate-like segments on oxide surfaces are now directly
  comparable (the convention match is new as of sf1d 0.2.0) — use as a prior
  before trusting new numbers.
- **Comb/branched refinement (Namics, not sf1d):** `sf1d` only gives one
  `chis` per whole repeat unit, lumping ST's ester backbone with its long
  C18 tail. [Namics](https://github.com/leermakers/Namics) supports
  branched/comb architectures, so a ST/MMA comb model (shared ester-type
  backbone bead + a side-chain branch — none for MMA, a long alkyl branch
  for ST) lets `chis` be derived per *chemical moiety* instead: one short
  ester-analog PMF shared by both monomers, one separate alkane-tail PMF
  specific to ST. This is a genuinely more faithful model, not just a
  validation pass, and is the main reason to reach for Namics over `sf1d`
  in this plan. See `../comparison/namics/` for the parameter mapping
  between the two codes (now exact, since sf1d 0.2.0).
- **Whole-chain consistency check:** run the real copolymer adsorbing from
  dilute solution, compare adsorbed amount/train fraction against `sf1d`
  predictions using the fragment-derived `chis`; a mismatch points to
  chain-connectivity/cooperativity effects the isolated-monomer estimate
  misses.

## Step 3 — validate together

Feed derived `chi`/`chis` into `sf1d` for the actual PSMA-stat-PMMA (random
sequence — pass a raw array to `SF1D()`, not the block-string CLI syntax,
since `SF1D.__init__` accepts `np.asarray(seq)` directly) and
PSMA-block-PMMA sequences, compare predicted `phi(z)`/trains/loops/tails/
adsorbed amount against full-chain MD slab profiles. Bulk-`chi` issues show
up as good-solvent exponent mismatches (`docs/theory.md` limitations);
surface-`chis` issues show up as irreversible/kinetic trapping
(strongly-binding polar anchors not equilibrating).

## Division of labor

Do the fitting/consistency loop (steps 1 and 3, and the primary route in
step 2) in `sf1d` — its Python API is a direct in-process call, far faster
to iterate than Namics' file-based I/O. Reach for Namics specifically for
the comb/branched moiety-level `chis` in step 2, which `sf1d` structurally
can't do.

## Still open

- `forcefield/oplsaa.ff/` is needed (via `-DLOPLS`) to grompp
  `solvent/hexadecane.itp` and the PAMA monomer `.ff` files — this
  machine's stock GROMACS 2026.0 `oplsaa.ff` does not have that patch.
- Namics itself isn't built on this machine yet
  (`tests/test_namics.py` skips without `NAMICS=/path/to/namics/binary`).
