# Patched oplsaa.ff (L-OPLS)

Copied from the local checkout at `../../../2026-10-28-generate-aa-pama/in/000_polyply/oplsaa.ff/`
(a clone of `zhanglab-nagoya/2026-01-26-generate-aa-pama`; that path is
gitignored *there*, so it couldn't be fetched via the GitHub API — this is
a local filesystem copy, not a GitHub extraction like the rest of
`forcefield/`).

Base: the stock `oplsaa.ff` shipped with this machine's GROMACS 2026.0
(`~/opt/gromacs-2026.0/share/gromacs/top/oplsaa.ff/`), confirmed identical
file-for-file except `ffbonded.itp` and `ffnonbonded.itp`, which carry the
`#ifdef LOPLS` blocks from [wesbarnett/lopls](https://github.com/wesbarnett/lopls):

- `ffnonbonded.itp`: `opls_135`/`opls_136` (CH3/CH2 carbon, LOPLS-specific
  partial charges), `opls_140A`/`opls_140B` (CH3-H/CH2-H, LOPLS-specific
  epsilon split), `opls_142`/`opls_144` (LOPLS alkene C/H) — all gated
  behind `#ifdef LOPLS`/`#ifndef LOPLS`, so default (non-LOPLS) behavior is
  unchanged when `-DLOPLS` isn't passed to `grompp`.
- `ffbonded.itp`: 14 dihedral entries gated the same way (alkane backbone
  `CT CT CT CT`, ester `C_2 OS CT CT`, alcohol/ether `CT CT CT O*`, etc.),
  each with an L-OPLS-refit alternative under `#ifdef LOPLS`.

Confirms what `forcefield/solvent/README.md` already states:
`forcefield/solvent/hexadecane.itp` and the PAMA monomer `.ff` files in
`forcefield/monomers/` need `-DLOPLS` at grompp time, and **this machine's
stock `oplsaa.ff` does not have the patch** — point `grompp`'s include path
at this directory (or patch the installed one the same way) rather than the
stock GROMACS force field.
