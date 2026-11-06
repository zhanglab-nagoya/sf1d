"""Step 2 chis stage 4: LAMMPS minimize + NVT equilibration of the merged
substrate/liquid slab (results/slab/{species}/merged.data,
step2_merge.smk), in one continuous LAMMPS session -- see
lmp/equilibrate_merge.lmp.template's module docstring for why this is
NOT split into separate minimize/NVT rules: `bond_style hybrid
harmonic`'s state didn't survive a write_data/read_data OR
write_restart/read_restart round trip reliably between two separate
`lmp` invocations (confirmed live 2026-11-06), and keeping both stages in
one session sidesteps that whole problem.

Both stages freeze the substrate (already relaxed on its own terms by
step2_substrate.smk and not re-relaxed or thermally perturbed here --
this pipeline treats it as a static wall, matching sf1d's own SCF
model).

STATUS: DONE. Run end-to-end for both ST and MMA 2026-11-06: minimization
converges cleanly on the energy-tolerance criterion in ~220 steps (no
NaN/blow-up), and NVT holds the mobile liquid at the target temperature
(see lmp/equilibrate_merge.lmp.template's comment on why the thermo
`temp` column itself reads low -- a display-only artifact of counting the
frozen substrate's DOF, not a thermostat problem).
"""

SLAB_RES = "results/slab"  # shared with step2_dilute.smk/convert/merge
N_SUB_TYPES = len(config["slab"]["substrate_creation"]["masses"])  # Fe, O


def n_org_types(merged_data_path):
    with open(merged_data_path) as fh:
        for line in fh:
            if "atom types" in line:
                return int(line.split()[0]) - N_SUB_TYPES
    raise ValueError(f"no 'atom types' line found in {merged_data_path}")


rule render_equilibrate_merge_lmp:
    input:
        "lmp/equilibrate_merge.lmp.template",
    output:
        f"{SLAB_RES}/{{species}}/equilibrate_merge.lmp",
    params:
        n_org_types=lambda wc, input: n_org_types(f"{SLAB_RES}/{wc.species}/merged.data"),
        cutoff=config["slab"]["merge_cutoff_angstrom"],
        kspace_accuracy=config["slab"]["merge_kspace_accuracy"],
        etol=config["slab"]["merge_minimize"]["etol"],
        ftol=config["slab"]["merge_minimize"]["ftol"],
        maxiter=config["slab"]["merge_minimize"]["maxiter"],
        maxeval=config["slab"]["merge_minimize"]["maxeval"],
        temperature_K=config["temperature_K"],
        seed=config["slab"]["merge_nvt_seed"],
        timestep_fs=config["slab"]["merge_nvt_timestep_fs"],
        nsteps=lambda wc: int(config["slab"]["merge_nvt_ps"] * 1000 / config["slab"]["merge_nvt_timestep_fs"]),
    wildcard_constraints:
        species="ST|MMA",
    localrule: True
    script:
        "../scripts/render_template.py"


rule lmp_equilibrate_merged:
    input:
        lmp=f"{SLAB_RES}/{{species}}/equilibrate_merge.lmp",
        data=f"{SLAB_RES}/{{species}}/merged.data",
        pair_coeffs=f"{SLAB_RES}/{{species}}/merged_pair_coeffs.lmp",
    output:
        minimized_data=f"{SLAB_RES}/{{species}}/minimized.data",
        minimized_restart=f"{SLAB_RES}/{{species}}/minimized.restart",
        nvt_data=f"{SLAB_RES}/{{species}}/nvt.data",
        nvt_restart=f"{SLAB_RES}/{{species}}/nvt.restart",
    shell:
        'cd "$(dirname {input.lmp})" && lmp -in "$(basename {input.lmp})"'
