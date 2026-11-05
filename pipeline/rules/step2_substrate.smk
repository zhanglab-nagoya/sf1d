"""Step 2 substrate: build the Fe2O3 substrate from scratch, instead of
reading a pre-built, opaque-provenance .data file (forcefield/substrate/
fe2o3.data -- see forcefield/substrate/README.md: externally prepared,
minimization inputs/outputs never in any git history).

Same methodology as ../2026-11-05-pama-in-nanogap's 299_substrate_creation
stage: cut hematite (alpha-Fe2O3) from a bulk CIF (COD 9015964), tile to a
target footprint, orthogonalize, assign INTERFACE force field charges + LJ
pair coefficients (OPLS-AA-compatible combination rules), and relax (CG
minimization, lateral box/relax) under that force field. See config.yaml's
slab.substrate_creation block for the geometry/FF knobs (deliberately a
much smaller footprint/thickness than that repo's nanogap-shear-scale
slab -- this one only needs to support a dilute fragment, not a confined
single-polymer system) and forcefield/substrate/README.md for the full
force-field provenance and the geometry pitfalls (non-rectangular surface
cells, stoichiometric z-trimming) build_slab.py handles.

STATUS: runs end-to-end (substrate_minimized.data target), confirmed
2026-11-05. Result at this repo's (smaller) target_nm: 8400 atoms (3360
Fe, 5040 O, exactly 2:3 stoichiometric -- 0 undercoordinated atoms per
orthogonalize_slab's sanity check), footprint 50.3 x 52.3 A pre-relaxation
growing to 51.8 x 53.8 A (~2.9%) post-minimization, CG-converged on energy
tolerance in 297 iterations (~18 s wall time). The LAMMPS steps need
`source ~/opt/lammps-22Jul2025/env.sh` (or equivalent) before launching
snakemake, for LD_LIBRARY_PATH (liblammps.so.0), not just PATH.
"""

SUB_DIR = "results/substrate"
_SUB = config["slab"]["substrate_creation"]
_CIF = _SUB["cif_file"]


rule copy_cif:
    input:
        f"{FF_DIR}/substrate/{_CIF}",
    output:
        f"{SUB_DIR}/{_CIF}",
    localrule: True
    shell:
        "cp {input} {output}"


rule build_slab:
    """Cut, tile, and z-trim the bulk CIF. Pure geometry, no force field."""
    input:
        f"{SUB_DIR}/{_CIF}",
    output:
        f"{SUB_DIR}/slab_triclinic.xyz",
    params:
        miller=_SUB["surface_miller"],
        layers=_SUB["layers"],
        target_nm=_SUB["target_nm"],
        vacuum=_SUB["vacuum"],
    conda:
        "../envs/ase.yaml"
    localrule: True
    shell:
        "python scripts/build_slab.py {input} -o {output} "
        "--miller {params.miller} --layers {params.layers} "
        "--target-nm {params.target_nm} --vacuum {params.vacuum}"


rule orthogonalize_slab:
    input:
        f"{SUB_DIR}/slab_triclinic.xyz",
    output:
        f"{SUB_DIR}/slab_orthorhombic.xyz",
    conda:
        "../envs/ase.yaml"
    localrule: True
    shell:
        "python scripts/orthogonalize_slab.py {input} -o {output}"


rule write_lammps_charge_data:
    """Assign INTERFACE FF charges and write a LAMMPS data file (no pair
    coefficients yet -- see lmp_write_coeffs)."""
    input:
        f"{SUB_DIR}/slab_orthorhombic.xyz",
    output:
        f"{SUB_DIR}/slab_charged.data",
    params:
        charge_fe=_SUB["charges"]["Fe"],
        charge_o=_SUB["charges"]["O"],
    conda:
        "../envs/ase.yaml"
    localrule: True
    shell:
        "python scripts/write_lammps_charge_data.py {input} -o {output} "
        "--charge-fe {params.charge_fe} --charge-o {params.charge_o}"


rule render_write_coeffs_lmp:
    input:
        f"{FF_DIR}/substrate/write_coeffs.lmp.template",
    output:
        f"{SUB_DIR}/write_coeffs.lmp",
    params:
        cutoff=_SUB["cutoff"],
        kspace_accuracy=_SUB["kspace_accuracy"],
        mass_fe=_SUB["masses"]["Fe"],
        mass_o=_SUB["masses"]["O"],
        eps_fe=_SUB["pair_coeffs"]["Fe"]["epsilon"],
        sigma_fe=_SUB["pair_coeffs"]["Fe"]["sigma"],
        eps_o=_SUB["pair_coeffs"]["O"]["epsilon"],
        sigma_o=_SUB["pair_coeffs"]["O"]["sigma"],
    localrule: True
    script:
        "../scripts/render_template.py"


rule lmp_write_coeffs:
    """Inject masses + INTERFACE pair coefficients into the charge-only
    data file. LAMMPS needs LD_LIBRARY_PATH for liblammps.so.0 as well as
    PATH -- `source ~/opt/lammps-22Jul2025/env.sh` before running snakemake,
    same as GROMACS/packmol need PATH set (see pipeline/README.md)."""
    input:
        lmp=f"{SUB_DIR}/write_coeffs.lmp",
        data=f"{SUB_DIR}/slab_charged.data",
    output:
        f"{SUB_DIR}/slab_with_coeffs.data",
    shell:
        'cd "$(dirname {input.lmp})" && lmp -in "$(basename {input.lmp})"'


rule render_minimize_lmp:
    input:
        f"{FF_DIR}/substrate/minimize.lmp.template",
    output:
        f"{SUB_DIR}/minimize.lmp",
    params:
        cutoff=_SUB["cutoff"],
        kspace_accuracy=_SUB["kspace_accuracy"],
        etol=_SUB["minimize"]["etol"],
        ftol=_SUB["minimize"]["ftol"],
        maxiter=_SUB["minimize"]["maxiter"],
        maxeval=_SUB["minimize"]["maxeval"],
    localrule: True
    script:
        "../scripts/render_template.py"


rule lmp_minimize:
    """CG-minimize under the INTERFACE FF with the lateral box relaxed to
    zero lateral stress (z fixed) -- a 3D-periodic, bulk-like relaxation of
    the slab cell, not yet a 2D-periodic/slab-corrected treatment (that
    only matters once this substrate is placed in an actual confined
    adsorption geometry, downstream in step2_chis.smk)."""
    input:
        lmp=f"{SUB_DIR}/minimize.lmp",
        data=f"{SUB_DIR}/slab_with_coeffs.data",
    output:
        data=f"{SUB_DIR}/substrate_minimized.data",
        nocoeff=f"{SUB_DIR}/substrate_minimized_nocoeff.data",
        restart=f"{SUB_DIR}/substrate_minimized.restart",
    shell:
        'cd "$(dirname {input.lmp})" && lmp -in "$(basename {input.lmp})"'
