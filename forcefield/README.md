# MD reference data: PSMA / PMMA on Fe2O3 in hexadecane

See `chi_chis_from_md.md` for the plan this data feeds into: deriving
`sf1d`'s `chi`/`chis` parameters from MD.

Force field and structure data for the system behind `sf1d`'s parameters:
PSMA-*stat*-PMMA and PSMA-*block*-PMMA chains (PSMA = stearyl methacrylate
"buoy" = `ST`, PMMA = methyl methacrylate "anchor" = `MMA`) adsorbing from
hexadecane onto a Fe2O3 substrate. This is reference input for atomistic MD
used to *derive* `sf1d`'s `chi`/`chis` parameters; it is not consumed by
`sf1d` itself.

```
forcefield/
├── monomers/     — ST/MMA topology (GROMACS OPLS-AA) + 3D structure + chain-end caps
├── solvent/      — hexadecane topology, coordinates, build script
├── substrate/    — Fe2O3 slab (LAMMPS data file)
├── oplsaa.ff/    — patched stock GROMACS force field (adds L-OPLS #ifdef blocks)
├── custom-ffnonbonded.itp, PAMA_ST_MMA.oplsaa.LigParGen_links.ff, citations.bib
└── README.md (this file)
```

## Provenance

**Monomer topology + links**, from `jotelha/polyply_1.0`, branch
`2026-11-12-PAMA-parameters`, under `polyply/data/oplsaaLigParGen/`
(originally pulled from `2026-10-28-PAMA-parameters`; updated 2026-11-16 —
see "Update 2026-11-12/16" below):

| File here | Source |
|---|---|
| `monomers/PAMA_ST.oplsaa.LigParGen.ff` | `PAMA_ST.oplsaa.LigParGen.ff` |
| `monomers/PAMA_MMA.oplsaa.LigParGen.ff` | `PAMA_MMA.oplsaa.LigParGen.ff` |
| `monomers/CH3.ff`, `monomers/CH3n.ff` | unchanged since 2026-10-29 (generic PAMA chain-end methyl caps) |
| `PAMA_ST_MMA.oplsaa.LigParGen_links.ff` | filtered from `PAMA.oplsaa.LigParGen_links.ff`, unchanged since 2026-10-29 (last touched upstream 2026-03-16, before either PAMA-parameters branch existed) |
| `custom-ffnonbonded.itp` | from `../../etc/custom-ffnonbonded.itp` |
| `citations.bib` | filtered from `citations.bib` (dropped two entries for an unrelated battery-polymer monomer) |

**Monomer 3D structures**, from the same fork, `etc/PAMA-monomers/`:

| File here | Source |
|---|---|
| `monomers/st.pdb` | `etc/PAMA-monomers/st.pdb` (unchanged; LigParGen-generated geometry, job `ST_C8C92E`) |
| `monomers/mma.pdb` | `etc/PAMA-monomers/mma.pdb` (unchanged; job `MMA_7EA60E`) |

**Solvent and substrate**, from `zhanglab-nagoya/2026-01-26-generate-aa-pama`
(the actual indentation/friction-modifier MD pipeline this data feeds into):
see `solvent/README.md` and `substrate/README.md` for full provenance,
including the L-OPLS caveat on `hexadecane.itp` and the ~11 MB size of
`substrate/fe2o3.data`.

**Patched oplsaa.ff**, resolving that L-OPLS caveat: `oplsaa.ff/` is a copy
of the same source repo's locally-staged, gitignored
`in/000_polyply/oplsaa.ff/` (confirmed identical to this machine's stock
GROMACS 2026.0 `oplsaa.ff` except for the L-OPLS `#ifdef` patches in
`ffbonded.itp`/`ffnonbonded.itp`) — see `oplsaa.ff/README.md`. Needed
because this machine's installed `oplsaa.ff` does **not** have that patch.

The source repo covers six PAMA side-chain variants (`EH`, `LA`, `MMA`, `OC`,
`ST`, `TD`) plus PEO/PEPP/PMAP/PTMA/PVBP monomers for other projects. Its
`PAMA.oplsaa.LigParGen_links.ff` is the combinatorial backbone-bond/angle/
dihedral table across *all* of those residue names (1752 `[ link ]` blocks).
`PAMA_ST_MMA.oplsaa.LigParGen_links.ff` here keeps only the 288 blocks whose
`resname`s are a subset of `{ST, MMA, CH3, CH3n}` — i.e. exactly what's
needed to connect ST/MMA backbone segments and cap the chain ends, dropping
everything about the other four side chains. See the source repo's
`etc/generate_pama_links.py` if the link set ever needs regenerating (e.g.
to add a third monomer).

`custom-ffnonbonded.itp` is kept whole (254 lines, covers all six monomers)
rather than filtered — it's small and filtering by atom-type ID risked
silently dropping a type referenced only via an angle/dihedral line.

## Update 2026-11-12/16: backbone/LJ audit fixes, L-OPLS reverted

Pulled in from the source repo's `2026-11-12-PAMA-parameters` branch
(superseding `2026-10-28-PAMA-parameters`, checked into this repo on a
dated copy, `2026-11-16-sf1d`): `monomers/PAMA_ST.oplsaa.LigParGen.ff`,
`monomers/PAMA_MMA.oplsaa.LigParGen.ff`, `custom-ffnonbonded.itp`, and
`citations.bib` all updated; `CH3.ff`/`CH3n.ff`/the links file/the two
monomer `.pdb` geometries are unchanged. Per the source repo's
`etc/report.md` sections 6-10:

- **L-OPLS reverted** (section 7): the 2026-10-28 retyping of ST's long
  alkyl tail to L-OPLS-specific charges/LJ/dihedrals is reverted back to
  plain OPLS-AA, so the whole ST repeat unit is now consistently standard
  OPLS-AA with library charges throughout — no L-OPLS content remains in
  any monomer. The `citations.bib` L-OPLS entries are kept, but only as a
  record of the reverted attempt (see that file's own comment).
- **sp2→sp3 backbone residuals fixed** (sections 8-10): a systematic
  audit of every bond/angle/dihedral/LJ value touching the backbone
  carbons and ester region against stock OPLS-AA turned up several more
  terms LigParGen had parametrized against the free (sp2, vinyl) monomer
  rather than the polymerized (sp3) one — stale `CM`-flavored bond/angle
  values at the alpha carbon's two other substituents, two mistyped
  `CM-HC`-flavored Cbeta-H bonds, several internally-inconsistent
  per-molecule dihedral fits around the alpha carbon, a wrong generic-
  ether LJ on the ester ether oxygen (should be the ester-specific
  `opls_467`), and (caught separately) a bug in the L-OPLS-revert script
  that used the wrong stock dihedral curve for the raw `CT-CT-CT-CT` tail
  torsion. All now corrected and tagged in-line (`; sp3 backbone fix
  (was ...)`, `; corrected ... (was ...)`, `; OPLS-AA (reverted from
  L-OPLS...)`). Net charge confirmed exactly 0 for both ST and MMA after
  the update (re-checked directly in this repo, not just taken on
  faith).
- Everything else (atom numbering, residue/atom names, `[ pairs ]`,
  overall monomer structure) is unchanged — this is a parameter
  correction, not a re-derivation.

Since ST's alkyl tail no longer uses L-OPLS, the only remaining reason
this repo needs `oplsaa.ff/`'s L-OPLS `#ifdef` patch (and every `mdp/
*.mdp`'s `-DLOPLS` define) is hexadecane's own L-OPLS dependency
(`solvent/README.md`) — unrelated to this update, not touched here.

## What's in each monomer

- **ST** = stearyl methacrylate (long C18 alkyl ester side chain, standard
  OPLS-AA as of the 2026-11-12/16 update above; atom types
  `opls_1135`–`opls_11xx`). This is the oil-soluble "buoy" monomer.
- **MMA** = methyl methacrylate (short polar ester side chain; atom types
  `opls_1084`–`opls_1098`). This is the surface-affine "anchor" monomer.
- **CH3** / **CH3n** = generic methyl caps for the chain head and tail
  (standard OPLS-AA types `opls_135`/`opls_140`, already in stock GROMACS
  `oplsaa.ff` — not redefined here).

Each monomer `.ff` file is a GROMACS `[ moleculetype ]` fragment (atoms,
bonds, angles, dihedrals, pairs) with its `[ atomtypes ]` block already
stripped (per the source repo's convention — atom types live only in
`custom-ffnonbonded.itp` plus stock `oplsaa.ff`). Backbone sp3
hybridization fixes (methacrylate C=C → saturated backbone) were already
applied upstream; see that repo's `etc/report.md` for the full rationale
if force constants look surprising.

`monomers/st.pdb` and `monomers/mma.pdb` are standalone single-monomer 3D
structures (LigParGen output geometry, atom naming matches the `.ff` files)
— useful as a starting geometry for building oligomer melts for the
solubility-parameter route in the chi/chis-from-MD plan, or for visual
sanity-checking atom types against the topology.

## Using with polyply

Point polyply at this directory (or symlink it alongside a standard
`oplsaa.ff/`) and give it a sequence of `ST`/`MMA` residues framed by
`CH3n` on **both** ends (not `CH3`) to build a linear PSMA/PMMA chain —
statistical (random residue order) or block (contiguous runs), per the
project's compositions of interest. `custom-ffnonbonded.itp` must be on
the GROMACS include path alongside the stock `oplsaa.ff/ffnonbonded.itp`.

`CH3.ff`'s plain literature OPLS charges (C -0.3029, H 0.0694 ×3) sum to
-0.0947, not 0 — only `CH3n.ff`'s charges (same values, uniformly shifted
by +0.023675/atom) sum to exactly 0. `PAMA_ST_MMA.oplsaa.LigParGen_links.ff`
bonds `CH3`/`CH3n` identically to either backbone attachment carbon, so
they're structurally interchangeable, but capping with `CH3` on even one
end leaves the whole molecule charged (confirmed live 2026-11-12: every
oligomer/fragment this pipeline built was off by exactly -0.0947 e, and
every LAMMPS run warned "System is not charge neutral", until both ends
were switched to `CH3n`). `pipeline/rules/step0_step1_bulk.smk` and
`step2_dilute.smk` already build every `-seq` string as `CH3n:1 ... CH3n:1`
for this reason — `CH3.ff` is kept here only for provenance, never for
capping.

## Citations

See `citations.bib`: OPLS-AA / LigParGen (Jorgensen & Tirado-Rives 2005;
Dodda et al. 2017 ×2), and polyply itself (Grunewald et al. 2022). The
L-OPLS entries (Siu et al. 2012; Pluhackova et al. 2015) are kept only as
a record of the reverted 2026-10-28 L-OPLS retyping attempt (see "Update
2026-11-12/16" above) — no monomer here actually uses L-OPLS.
