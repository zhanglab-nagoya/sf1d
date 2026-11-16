"""Step 2 chis stage 3: merge the equilibrated solvent-only box
(results/slab/solvent/converted.lmp, step2_convert.smk) onto the Fe2O3
substrate (results/substrate/substrate_minimized.data,
step2_substrate.smk) into one periodic LAMMPS data file.
SPECIES-INDEPENDENT since 2026-11-13's solvent/fragment split -- see
step2_dilute.smk's module docstring for the full root-cause history of why
the fragment no longer lives in this box at all (it used to be centered
here via scripts/center_fragment.py, replacing an even earlier post-
equilibration scripts/insert_fragment.py design; both are gone now --
equilibrate solvent+substrate first, insert the fragment only afterward
into the converged result, see rules/step2_equilibrate.smk).

See scripts/merge_slab.py's own docstring for merge_solvent_with_substrate's
two deliberate design choices: a fully 3D-periodic substrate/liquid
sandwich (liquid's top wraps around to the substrate's own bottom face)
instead of a vacuum-padded free surface -- avoiding LAMMPS's
kspace_modify slab correction and the large vacuum padding it would need,
consistent with the plain 3D PME/PPPM treatment used everywhere else in
this pipeline -- and relying on geometric mixing (already each force
field's own native combination rule) for every organic<->substrate cross
pair term, so no explicit cross pair_coeffs are ever written.

Smoke-tested live 2026-11-06 for ST (solvent+fragment still combined at
that point): the merged.data + merged_coeffs.lmp pair loads and evaluates
cleanly in LAMMPS (`read_data` + `run 0`, pair_style lj/cut/coul/long 11.0
11.0, kspace_style pppm 1.0e-4, pair_modify mix geometric). E_vdwl comes
out large and positive (~1e5 kcal/mol) -- this looks alarming but is NOT
an interface clash: the bare substrate alone (no liquid at all) already
has E_vdwl ~+104000 kcal/mol under this same force field (confirmed
live), since for this ionic force field the LJ channel mostly encodes
short-range Pauli repulsion between ions while the actual cohesion is
carried by the (hugely negative) Coulomb terms -- see config.yaml's
slab.merge_minimize comment for the actual (small) interface strain this
merge introduces, and step2_equilibrate.smk for how it gets resolved.

STATUS: re-run and validated for both species under the 2026-11-13
solvent/fragment split (same merge_slab.py logic, now on a solvent-only
box) as part of the 2026-11-14 end-to-end verification -- see
rules/step2_equilibrate.smk's module docstring for that verification's
results.
"""

SLAB_RES = "results/slab"  # shared with step2_dilute.smk/step2_convert.smk
SOLVENT_RES = f"{SLAB_RES}/solvent"


rule merge_solvent_with_substrate:
    input:
        organic_data=f"{SOLVENT_RES}/converted.lmp",
        organic_input=f"{SOLVENT_RES}/converted.input",
        substrate_data=config["slab"]["substrate_data"],
    output:
        data=f"{SOLVENT_RES}/merged.data",
        coeffs=f"{SOLVENT_RES}/merged_coeffs.lmp",
    params:
        gap=config["slab"]["merge_gap_angstrom"],
    localrule: True
    shell:
        "python scripts/merge_slab.py "
        "--organic-data {input.organic_data} --organic-input {input.organic_input} "
        "--substrate-data {input.substrate_data} --gap {params.gap} "
        "-o {output.data} --out-coeffs {output.coeffs}"
