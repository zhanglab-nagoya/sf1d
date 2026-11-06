"""Step 2 chis stage 3: merge each species' converted organic dilute box
(results/slab/{species}/converted.lmp, step2_convert.smk) onto the Fe2O3
substrate (results/substrate/substrate_minimized.data, step2_substrate.smk)
into one periodic LAMMPS data file.

See scripts/merge_slab.py's own docstring for the two deliberate design
choices: a fully 3D-periodic substrate/liquid sandwich (liquid's top
wraps around to the substrate's own bottom face) instead of a
vacuum-padded free surface -- avoiding LAMMPS's kspace_modify slab
correction and the large vacuum padding it would need, consistent with
the plain 3D PME/PPPM treatment used everywhere else in this pipeline --
and relying on geometric mixing (already each force field's own native
combination rule) for every organic<->substrate cross pair term, so no
explicit cross pair_coeffs are ever written.

Smoke-tested live 2026-11-06 for ST: the merged.data + merged_pair_coeffs.lmp
pair loads and evaluates cleanly in LAMMPS (`read_data` + `run 0`, pair_style
lj/cut/coul/long 11.0 11.0, kspace_style pppm 1.0e-4, pair_modify mix
geometric) -- E_vdwl is large and positive (~1e5 kcal/mol), expected for a
freshly-placed, unrelaxed interface and exactly what stage 4's energy
minimization (not yet written) needs to resolve, not a sign of a broken
merge.

STATUS: scaffolded and smoke-tested for ST; not yet run for MMA, not
wired into an actual Snakemake-level LAMMPS energy check here (that
belongs to stage 4's minimize rule).
"""

SLAB_RES = "results/slab"  # shared with step2_dilute.smk/step2_convert.smk


rule merge_with_substrate:
    input:
        organic_data=f"{SLAB_RES}/{{species}}/converted.lmp",
        organic_input=f"{SLAB_RES}/{{species}}/converted.input",
        substrate_data=config["slab"]["substrate_data"],
    output:
        data=f"{SLAB_RES}/{{species}}/merged.data",
        pair_coeffs=f"{SLAB_RES}/{{species}}/merged_pair_coeffs.lmp",
    params:
        gap=config["slab"]["merge_gap_angstrom"],
    wildcard_constraints:
        species="ST|MMA",
    localrule: True
    shell:
        "python scripts/merge_slab.py "
        "--organic-data {input.organic_data} --organic-input {input.organic_input} "
        "--substrate-data {input.substrate_data} --gap {params.gap} "
        "-o {output.data} --out-pair-coeffs {output.pair_coeffs}"
