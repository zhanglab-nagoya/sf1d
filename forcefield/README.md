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
├── monomers/     — ST/MMA topology (GROMACS OPLS-AA/L-OPLS) + 3D structure + chain-end caps
├── solvent/      — hexadecane topology, coordinates, build script
├── substrate/    — Fe2O3 slab (LAMMPS data file)
├── oplsaa.ff/    — patched stock GROMACS force field (adds L-OPLS #ifdef blocks)
├── custom-ffnonbonded.itp, PAMA_ST_MMA.oplsaa.LigParGen_links.ff, citations.bib
└── README.md (this file)
```

## Provenance

**Monomer topology + links**, from `jotelha/polyply_1.0`, branch
`2026-10-28-PAMA-parameters`, under `polyply/data/oplsaaLigParGen/`:

| File here | Source |
|---|---|
| `monomers/PAMA_ST.oplsaa.LigParGen.ff` | `PAMA_ST.oplsaa.LigParGen.ff` (unchanged) |
| `monomers/PAMA_MMA.oplsaa.LigParGen.ff` | `PAMA_MMA.oplsaa.LigParGen.ff` (unchanged) |
| `monomers/CH3.ff`, `monomers/CH3n.ff` | unchanged (generic PAMA chain-end methyl caps) |
| `PAMA_ST_MMA.oplsaa.LigParGen_links.ff` | filtered from `PAMA.oplsaa.LigParGen_links.ff` |
| `custom-ffnonbonded.itp` | from `../../etc/custom-ffnonbonded.itp` (unchanged) |
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

## What's in each monomer

- **ST** = stearyl methacrylate (long C18 alkyl ester side chain via L-OPLS;
  atom types `opls_1135`–`opls_11xx`). This is the oil-soluble "buoy" monomer.
- **MMA** = methyl methacrylate (short polar ester side chain; atom types
  `opls_1084`–`opls_1098`). This is the surface-affine "anchor" monomer.
- **CH3** / **CH3n** = generic methyl caps for the chain head and tail
  (standard OPLS-AA types `opls_135`/`opls_140`, already in stock GROMACS
  `oplsaa.ff` — not redefined here).

Each monomer `.ff` file is a GROMACS `[ moleculetype ]` fragment (atoms,
bonds, angles, dihedrals, pairs) with its `[ atomtypes ]` block already
stripped (per the source repo's convention — atom types live only in
`custom-ffnonbonded.itp` plus stock `oplsaa.ff`). Backbone sp3 hybridization
fixes (methacrylate C=C → saturated backbone) and L-OPLS alkyl retyping were
already applied upstream; see that repo's `etc/report.md` for the full
rationale if force constants look surprising.

`monomers/st.pdb` and `monomers/mma.pdb` are standalone single-monomer 3D
structures (LigParGen output geometry, atom naming matches the `.ff` files)
— useful as a starting geometry for building oligomer melts for the
solubility-parameter route in the chi/chis-from-MD plan, or for visual
sanity-checking atom types against the topology.

## Using with polyply

Point polyply at this directory (or symlink it alongside a standard
`oplsaa.ff/`) and give it a sequence of `ST`/`MMA` residues framed by `CH3`
(head) and `CH3n` (tail) to build a linear PSMA/PMMA chain — statistical
(random residue order) or block (contiguous runs), per the project's
compositions of interest. `custom-ffnonbonded.itp` must be on the GROMACS
include path alongside the stock `oplsaa.ff/ffnonbonded.itp`.

## Citations

See `citations.bib`: OPLS-AA / LigParGen (Jorgensen & Tirado-Rives 2005;
Dodda et al. 2017 ×2), L-OPLS for the alkyl side chains (Siu et al. 2012;
Pluhackova et al. 2015), and polyply itself (Grunewald et al. 2022).
