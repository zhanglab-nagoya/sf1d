# Hexadecane solvent

Extracted unchanged from `zhanglab-nagoya/2026-01-26-generate-aa-pama`,
`in/100_solvent/` (default branch `main`).

Coordinates (`hexadecane-raw.pdb` → `hexadecane.pdb`/`hexadecane.gro`, built
with `rdkit` via `create.py`) have been unchanged since creation.
`hexadecane_charged.itp` is the original LigParGen (self-contained
`opls_8XX` types) topology; it is superseded by `hexadecane.itp`, rebuilt
(2026-10-28) to use standard typed OPLS-AA atoms (`opls_135`/`136` for
CH3/CH2, `opls_140A`/`140B` for their hydrogens) per
[wesbarnett/lopls](https://github.com/wesbarnett/lopls), so bonded terms
resolve from `oplsaa.ff/ffbonded.itp` by atom class instead of listing
explicit constants. Both `.itp` files are kept for provenance.

**`hexadecane.itp` requires `define = -DLOPLS` in any `.mdp` that grompps
it** — without it `opls_140A`/`opls_140B` don't exist and `grompp` fails
fast (by design, not a bug). This means the GROMACS `oplsaa.ff` used must
have the `#ifdef LOPLS` blocks from the lopls repo patched in (nonbonded:
`opls_135`/`136`, `opls_140A`/`140B`, `142`, `144`; bonded: the 14
LOPLS-tagged dihedrals) — source repo stages this locally at
`in/000_polyply/oplsaa.ff` (gitignored there, not copied here). Check
whether the same patch is needed against this machine's GROMACS 2026.0
`oplsaa.ff` before running anything with `hexadecane.itp`.

Charges are exactly neutral by construction
(2×CH3×(−0.222) + 14×CH2×(−0.148) + 6×H(CH3)×0.074 + 28×H(CH2)×0.074 = 0),
so `renormalize_charges.py` (kept here, used for the PAMA monomers'
charge normalization, not for this file) is not needed for hexadecane.

Not copied from the source repo: the `.mdp`/`packmol`/`.top` template files
(`em.mdp`, `nvt.mdp`, `npt.mdp`, `packmol.inp.template`,
`solvent.top.template`) — those are this specific pipeline's box-building
steps rather than the hexadecane model itself; regenerate or re-fetch them from the source repo if actually running an
NPT/CED calculation to derive `chi`.
