"""Step 2: chis_A, chis_B from dilute ST-analog/MMA-analog adsorption near
the Fe2O3 substrate (results/substrate/substrate_minimized.data, built
from scratch by step2_substrate.smk -- see forcefield/substrate/README.md).

STATUS: design only, not implemented. Unlike step0/step1 (plain GROMACS
liquid boxes), this needs a GROMACS-build-then-LAMMPS-convert-then-merge
pipeline, because the substrate is a LAMMPS data file
(pair_style lj/cut/coul/long) while the organic fragment + hexadecane side
is naturally built and parametrized in GROMACS (OPLS-AA/L-OPLS, PME).
zhanglab-nagoya/2026-01-26-generate-aa-pama already solves exactly this
conversion (in/050_conversion, in/310_polymer_preparation,
in/315_solvent_preparation, in/320_polymer_in_solvent_on_substrate) via a
GROMACS singlepoint run -> VMD/TopoTools atom-order canonicalization
(sort.tcl) -> LAMMPS write_data -> merge.tcl.template. Per the decision to
keep this pipeline standalone (not touching that repo), the same general
approach needs reimplementing here -- that merge/convert step is nontrivial
and was deliberately NOT fabricated blind; treat the stage names below as
the intended shape, to fill in and test incrementally, not as working rules.

Intended stages (one per species: ST-analog fragment, MMA-analog fragment,
and hexadecane itself as the reference for the correction term):

1. build_dilute_fragment_box(species)  - DONE, see rules/step2_dilute.smk:
   one fragment molecule + many hexadecane molecules, em/nvt/npt
   equilibrated in bulk first (reuse step0/step1's mdp/ and polyply
   machinery). Produces results/slab/{species}/npt.gro, equilibrated to
   ~748 kg/m^3 (confirmed live 2026-11-05 for ST and MMA), matching
   hexadecane's own reference density as expected.
2. convert_to_lammps(species)          - DONE, see rules/step2_convert.smk
   and scripts/gmx_to_lammps.py: unwraps PBC (gmx trjconv -pbc whole) and
   converts via InterMol's Python API (not its CLI -- needs the LOPLS
   #define, which the CLI has no way to pass) to a LAMMPS input+data
   file pair, both already pair_style lj/cut/coul/long like the
   substrate. Cross-checked against a GROMACS single-point energy
   (gmx mdrun -rerun): ~0.5% agreement for ST, residual expected from
   GROMACS PME vs. LAMMPS PPPM grid/cutoff differences, not a bug.
3. merge_with_substrate(species)       - DONE, see rules/step2_merge.smk
   and scripts/merge_slab.py: concatenates the converted organic system
   with results/substrate/substrate_minimized.data into one fully
   3D-periodic sandwich (liquid wraps around in z to meet the
   substrate's own bottom face again -- no vacuum-padded free surface,
   no kspace_modify slab correction needed), offsetting atom/molecule/
   type IDs and relying on geometric mixing (both force fields' own
   native rule) for every cross pair term. Smoke-tested for both
   species: the merged.data + merged_pair_coeffs.lmp pair loads and
   evaluates cleanly in LAMMPS (finite energies, no errors).
4. minimize / nvt / npt slab(species)  - LAMMPS equilibration near the
   surface.
5. production_slab(species)            - longer LAMMPS run, dump a
   z-resolved density profile for the fragment and for hexadecane.
6. compute_chis                        - python: read the two density
   profiles (fragment's and hexadecane's), take the contact-layer vs. bulk
   ratio for each, and combine per forcefield/chi_chis_from_md.md step 2:

       chis_x = ln[phi_contact,x/phi_bulk,x] - ln[phi_contact,S/phi_bulk,S]

Do not add placeholder rules with `touch()`/`false` here the way
step0_step1_bulk.smk does for vacuum_ref -- that pattern is fine for a
single well-scoped missing piece, but this whole stage is substantial new
work (its own box-building, conversion, and equilibration protocol) that
deserves a real design pass, informed by the working conversion scripts in
the generate-aa-pama repo, before any code is written here.
"""
