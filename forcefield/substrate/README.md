# Fe2O3 substrate

## Generated from scratch (current pipeline, as of 2026-11-05)

`pipeline/rules/step2_substrate.smk` builds the substrate the pipeline
actually uses (`results/substrate/substrate_minimized.data`, pointed to by
`config.yaml`'s `slab.substrate_data`) from the files in this directory,
instead of reading the pre-built `fe2o3.data` below. Same methodology as
`../2026-11-05-pama-in-nanogap`'s `299_substrate_creation` stage: cut
hematite (α-Fe2O3) from a bulk CIF, tile/orthogonalize to a target
footprint, assign INTERFACE force field charges + LJ pair coefficients,
and relax (CG minimization, lateral box/relax) under that force field —
see that repo's `CLAUDE.md` ("Substrate creation (stage 299) pitfalls") for
the two geometry bugs this guards against (a non-rectangular surface cell
from a naive hexagonal-to-orthorhombic cut, and a stoichiometry-breaking
z-trim) that don't announce themselves as crashes.

- `9015964.cif`: hematite (α-Fe2O3), from the Crystallography Open
  Database (COD), <https://www.crystallography.net/cod/9015964.html>.
  Original data from the American Mineralogist Crystal Structure Database.
  COD entries may be used within the scientific community with attribution
  to the original journal article (see the CIF header) — Finger, L. W. &
  Hazen, R. M.
- `write_coeffs.lmp.template` / `minimize.lmp.template`: Jinja2-rendered
  LAMMPS scripts (geometry/FF knobs come from `config.yaml`'s
  `slab.substrate_creation` block) that inject INTERFACE pair coefficients
  and then CG-minimize the tiled slab.

### Force field: INTERFACE (Kanhaiya et al. 2023)

> Kanhaiya, Krishan, et al. "Accurate force fields for atomistic
> simulations of oxides, hydroxides, and organic hybrid materials up to the
> micrometer scale." *Journal of Chemical Theory and Computation* 19.22
> (2023): 8293-8322.

Table 2 of the paper gives **four** parameter variants per compound
(LJ form x mixing rule); only the **OPLS, CVFF** row (12-6 LJ, **geometric**
mixing) is compatible with OPLS-AA, which is what the rest of this
pipeline's organic force field uses — the CHARMM/AMBER row's ε/r_min values
are tuned for arithmetic mixing instead and are not interchangeable. For
α-Fe2O3 (that row):

| Atom | charge (e) | LJ ε (kcal/mol) | LJ σ (Å) |
|------|-----------:|----------------:|---------:|
| Fe   | +1.74      | 0.3             | 1.853    |
| O    | −1.16      | 0.2             | 3.127    |

(σ = 2⁻¹ᐟ⁶·r_min from the paper's r_min values, Fe 2.08 Å / O 3.51 Å.)
`pair_style lj/cut/coul/long`, `pair_modify mix geometric`. As an
independent cross-check, the paper's Table 6 reports the α-Fe2O3 (0001)
surface energy under IFF-CVFF/OPLS as 1.89 J/m² against an experimental
best estimate of 1.9 ± 0.1 J/m² — and `(0,0,1)` is exactly the Miller
index this pipeline cuts (`slab.substrate_creation.surface_miller`).

### Sizing: deliberately much smaller than the nanogap repo's slab

`config.yaml`'s `target_nm: [5.0, 5.0, 3.0]` (vs. the nanogap repo's
production `[16.0, 15.5, 5.0]`, 126720 atoms) gives ~8400 atoms here —
plenty for a flat patch under a single dilute adsorbing fragment (step 2's
actual use case, not a confined single-polymer shear system with a huge
solvent column). Confirmed live 2026-11-05: 8400 atoms (3360 Fe, 5040 O,
exactly 2:3 stoichiometric, 0 undercoordinated atoms), CG-converged in 297
iterations (~18 s). Revisit if step2_chis.smk's real adsorption rules turn
out to need a bigger footprint than this.

---

## Removed: `fe2o3.data` (superseded, dropped from git 2026-11-05)

Was extracted unchanged from `zhanglab-nagoya/2026-01-26-generate-aa-pama`,
`in/300_substrate_preparation/` (default branch `main`): a LAMMPS data file
(`write_data`, LAMMPS 22Jul2025, `units real`, `atom_style full`,
`pair_style lj/cut/coul/long`), 126,720 atoms, 2 atom types (Fe mass
55.845, O mass 15.9994), box 161.5 × 157.4 × 69.8 Å. Created and minimized
by Ji Ce (zhanglab), 2026-05-12; per the source repo, minimization
inputs/outputs lived on an SMB share
(`zhanglab5.local/public/JiC/1-Research/2026-04-21-sub-minimize`), never in
that git repo or here — this opaque provenance is exactly why the pipeline
now builds its own substrate from scratch instead, and why the ~11 MB file
itself was removed (`git rm`) rather than kept now that nothing reads it
(`slab.substrate_data` points at the generated `substrate_minimized.data`
above). Pull it back from the source repo's path above if a future diff
against it is ever needed.

`sort.tcl` (VMD/TopoTools) is the companion atom-reordering script from the
source repo, kept alongside for provenance (not run here) — small enough
that there's no reason to drop it along with `fe2o3.data`.

Not copied: `fe2o3_probe.data` (`in/325_probe_preparation/`, 22,452-atom
Fe2O3 probe tip for indentation, same provenance/author).
