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
geometric). E_vdwl comes out large and positive (~1e5 kcal/mol) -- this
looks alarming but is NOT an interface clash: the bare substrate alone
(no liquid at all) already has E_vdwl ~+104000 kcal/mol under this same
force field (confirmed live), since for this ionic force field the LJ
channel mostly encodes short-range Pauli repulsion between ions while the
actual cohesion is carried by the (hugely negative) Coulomb terms -- see
config.yaml's slab.merge_minimize comment for the actual (small) interface
strain this merge introduces, and step2_equilibrate.smk for how it gets
resolved.

STATUS: DONE. Smoke-tested for both ST and MMA, and now exercised
end-to-end by stage 4's minimize+NVT run (rules/step2_equilibrate.smk) for
both species.
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
