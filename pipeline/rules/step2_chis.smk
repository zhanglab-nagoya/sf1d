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
   species: the merged.data + merged_coeffs.lmp pair loads and
   evaluates cleanly in LAMMPS (finite energies, no errors).
4. minimize / nvt equilibrate_merged(species) - DONE, originally one
   continuous LAMMPS session (rules/step2_equilibrate.smk,
   lmp/equilibrate_merge.lmp.template) since `bond_style hybrid
   harmonic`'s coefficients did not survive a write_data/read_data or
   write_restart/read_restart round trip reliably between separate `lmp`
   invocations (confirmed live 2026-11-06). Rewritten 2026-11-08 into
   separate per-stage rules (rules/step2_equilibrate.smk,
   lmp/minimize.lmp.template + lmp/nvt.lmp.template + others) once
   that was fixed at the source instead -- see
   rules/step2_equilibrate.smk's module docstring. Verified for both
   species (original monolithic version): minimization converges cleanly
   on the energy-tolerance criterion in ~220 steps (no NaN/blow-up), and
   the NVT stage holds the mobile liquid at the target temperature
   (confirmed via a group-scoped `compute mobtemp mobile temp` -- the
   plain `temp` thermo column reads low because it normalizes by the
   frozen substrate's DOF too, even though those atoms carry zero kinetic
   energy; the thermostat itself was never affected, only that display
   column).
5. production_slab(species)            - IN PROGRESS (see STATUS below
   for what the first real run found and rules/step2_equilibrate.smk for
   the since-rewritten z-equilibration/fragment-placement design meant to
   fix it). Dumps a z-resolved number-density profile for the fragment group (molecule
   1:3) and the hexadecane group (molecule 4:max_org_molid) via
   `fix ave/chunk`, 1 Angstrom bins. config.yaml's slab.production_ns
   (10 ns as of 2026-11-06) was sized from a real flow2-ai throughput
   benchmark (profiles/flow2-ai/, job 102982.jobm: 6.8 ns/day/species on
   48 OpenMP threads) and submitted there directly (jobs 103064.jobm ST,
   103065.jobm MMA, 2026-11-06, ~35 h/species, ~530 points/species) --
   not yet complete. Syntax-verified locally first with a short run
   before committing to that cost. See config.yaml's slab.production_ns
   comment for why 10 ns is an EXPLORATORY length, not a guaranteed-
   converged one (single fragment molecule, might not revisit both
   contact and bulk if it adsorbs strongly -- check
   fragment_density.profile's later windows for continued drift, not
   just trust the run finished).
6. compute_chis                        - WRITTEN, see rule compute_chis /
   scripts/compute_chis.py below -- but it reveals a real sampling
   failure in stage 5's production data, not a usable chis_x yet (see
   STATUS below).

   Reads the two density profiles per species, classifies bins into a
   substrate-contact region (per face, validated against the solvent's
   own density there -- see compute_chis.py's module docstring for why a
   face can be invalid: a vacuum/free-surface edge, not a real substrate
   contact) and a bulk region, builds a running cumulative time-average
   over the retained windows (per-window ratios are unusable -- a single
   fragment molecule's instantaneous density is usually exactly zero in
   one of the two regions), and combines per
   forcefield/chi_chis_from_md.md step 2:

       chis_x = ln[phi_contact,x/phi_bulk,x] - ln[phi_contact,S/phi_bulk,S]

STATUS: run 2026-11-08 against both species' real production data
(103064.jobm ST, 103065.jobm MMA, both finished cleanly on flow2-ai, exit
status 0, 34.5h/36.4h real time, profiles pulled back locally since
flow2-ai has no conda env for this pipeline yet). The code runs and the
region classification is correct (independently confirmed by tracing the
fragment's own populated bins across all 199 windows for both species),
but the result is NOT a usable chis_x: for BOTH species, the single
fragment molecule was already at one fixed location (near the box's
periodic z-seam) by the start of production and never moved from it
across the full 10 ns -- see scripts/merge_slab.py's periodic-sandwich
design, and rules/step2_equilibrate.smk's module docstring for the
root cause (an under-equilibrated box z-length) and the since-rewritten
fix. For ST that location is on the INVALID face
(a real ~10 A vacuum gap formed there instead of the intended liquid-
substrate contact -- the fragment is stuck at a liquid/vacuum free
surface, not the Fe2O3 substrate), giving phi_contact=phi_bulk=0 (chis_A
= nan). For MMA that location IS a valid substrate contact, but the
fragment never once visited the bulk region either, giving phi_bulk=0
(chis_B = +inf). This is exactly the single-fragment convergence risk
flagged by config.yaml's slab.production_ns comment and
forcefield/chi_chis_from_md.md's umbrella-sampling-PMF fallback -- 10 ns
was not enough for either species' fragment to explore both states even
once. results/chis.yaml is NOT committed/trustworthy as of this run; see
the pipeline's own project memory for what to try next (longer run,
multiple independent replicas, starting the fragment in the bulk instead
of wherever the dilute-box building happened to leave it, or switching
straight to the PMF route) before spending more flow2-ai points on this.
"""

SLAB_RES2 = "results/slab"  # same dir as step2_equilibrate.smk's SLAB_RES


rule compute_chis:
    input:
        fragment_density=expand(f"{SLAB_RES2}/{{species}}/fragment_density.profile", species=["ST", "MMA"]),
        hexadecane_density=expand(f"{SLAB_RES2}/{{species}}/hexadecane_density.profile", species=["ST", "MMA"]),
        merged_data=expand(f"{SLAB_RES2}/{{species}}/merged.data", species=["ST", "MMA"]),
        segment_volume="results/segment_volume.yaml",
    output:
        chis="results/chis.yaml",
        diagnostics="results/chis_diagnostics.yaml",
    params:
        species=["ST", "MMA"],
        n_sub_types=len(config["slab"]["substrate_creation"]["masses"]),
        contact_width_a=config["slab"]["chis_contact_width_a"],
        bulk_margin_a=config["slab"]["chis_bulk_margin_a"],
        min_contact_fraction=config["slab"]["chis_min_contact_density_fraction"],
        discard_fraction=config["slab"]["chis_discard_fraction"],
        convergence_tol=config["slab"]["chis_convergence_tol"],
    conda:
        "../envs/analysis.yaml"
    localrule: True
    script:
        "../scripts/compute_chis.py"
