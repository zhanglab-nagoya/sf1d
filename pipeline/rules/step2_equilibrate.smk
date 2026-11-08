"""Step 2 chis stages 4+5: LAMMPS minimize + z-direction density
equilibration + fragment placement + restrained re-equilibration + NVT +
production of the merged substrate/liquid slab
(results/slab/{species}/merged.data, step2_merge.smk) -- one Snakemake
rule per stage, each its own small `lmp/*.lmp.template` file sharing a
common preamble (`lmp/_common_setup.lmp.template`, `{% include %}`d by
every stage template).

Rewritten 2026-11-08 from one continuous ~250-line LAMMPS session into
this per-stage shape. The original reason for one session was that
`write_data`/`write_restart` were found (2026-11-06) to silently drop
`Bond Coeffs`/`Angle Coeffs`/`Dihedral Coeffs` sections for this system's
`bond_style hybrid harmonic` setup, losing those coefficients across a
`read_data`/`read_restart` round trip between separate `lmp` invocations.
That's now fixed at the source instead of worked around:
`scripts/merge_slab.py`'s `merged_coeffs.lmp` sidecar (already the home
of `pair_coeff` commands, since InterMol never embeds those in a data
file at all) now also carries explicit `bond_coeff`/`angle_coeff`/
`dihedral_coeff` commands, one per row in merged.data's own Bond/Angle/
Dihedral Coeffs sections (only 8/9/16 *distinct* (style, params) tuples
for ST, confirmed live -- scripts/dedup_lammps_bonded_types.py already
collapses what InterMol originally wrote as one type per bond/angle/
dihedral *instance* down to the true count). Every stage template
`include`s that sidecar right after its own `read_data`, so coefficients
are re-supplied explicitly and deterministically every time, independent
of whatever write_data did or didn't preserve -- the actual blocker, not
"keep everything in one session".

Every stage freezes the substrate (already relaxed on its own terms by
step2_substrate.smk and not re-relaxed or thermally perturbed here --
this pipeline treats it as a static wall, matching sf1d's own SCF
model).

The old `displace_atoms` + blind minimize approach for placing the
fragment (a rigid translation straight into whatever hexadecane already
occupied that spot, confirmed live 2026-11-08 to produce a real
~3e13 kcal/mol transient that the minimizer happened to resolve) is
replaced by scripts/insert_fragment.py: it moves the fragment the same
way, but explicitly finds and removes the hexadecane molecules it would
otherwise overlap first, so the one remaining minimize
(`lmp_insert_minimize`) only has to relax a much smaller, localized
disturbance.

Stage 5's production run is long (slab.production_ns, 2026-11-06 default
10 ns/species) and genuinely needs real parallelism -- confirmed live via
a throughput benchmark on flow2-ai (profiles/flow2-ai/, job 102982.jobm):
48 OpenMP threads (1 NUMA, t1-share) gave 17.4x over this workstation's
serial baseline (0.111 s/step -> 0.00636 s/step), 6.8 ns/day/species. Every
LAMMPS rule below takes the same `omp_threads`/`-sf omp -pk omp N` pattern
(not just production), defaulting to 1 locally and meant to be driven up
per-rule via the `ompthreads` resource on a cluster profile.

STATUS: the z-direction equilibration + fragment placement + restrained
re-equilibration design (previously stages 4b-4d of one session) was
smoke-tested end-to-end 2026-11-08 on ST's real merged.data (tiny step
counts) with the OLD monolithic template and `fix press/berendsen` (not
`fix npt`/Nose-Hoover, which never converges on this system regardless of
target -- confirmed live on two independent attempts, see
lmp/npt_z.lmp.template's comment): zero LAMMPS errors, temperature stable
~295-305 K through every stage, fragment ends up sitting in the bulk
region it was placed in, not the periodic seam. This per-stage rule split
and scripts/insert_fragment.py are new as of this rewrite and still need
their own smoke test before trusting them at the real production_ns or on
flow2-ai.
"""

SLAB_RES = "results/slab"  # shared with step2_dilute.smk/convert/merge
N_SUB_TYPES = len(config["slab"]["substrate_creation"]["masses"])  # Fe, O
LMP_COMMON_TEMPLATE = "lmp/_common_setup.lmp.template"


def n_org_types(merged_data_path):
    with open(merged_data_path) as fh:
        for line in fh:
            if "atom types" in line:
                return int(line.split()[0]) - N_SUB_TYPES
    raise ValueError(f"no 'atom types' line found in {merged_data_path}")


def max_org_molid(merged_data_path):
    """merge_slab.py never renumbers organic molecule-IDs and gives the
    substrate the single id (max organic molid + 1) -- so the largest
    molid anywhere in the file, minus 1, is exactly the organic side's
    own max, with no need to separately locate the substrate's atoms."""
    max_molid = 0
    in_atoms = False
    with open(merged_data_path) as fh:
        for line in fh:
            stripped = line.strip()
            if stripped == "Atoms" or stripped.startswith("Atoms "):
                in_atoms = True
                continue
            if in_atoms:
                fields = stripped.split()
                if not fields:
                    continue
                if not fields[0].isdigit():
                    break
                max_molid = max(max_molid, int(fields[1]))
    if max_molid == 0:
        raise ValueError(f"no Atoms section found in {merged_data_path}")
    return max_molid - 1


def merged_data_path(species):
    return f"{SLAB_RES}/{species}/merged.data"


# --- stage 4a: minimize --------------------------------------------------

rule render_slab_minimize_lmp:
    input:
        "lmp/minimize.lmp.template",
        LMP_COMMON_TEMPLATE,
    output:
        f"{SLAB_RES}/{{species}}/minimize.lmp",
    params:
        n_org_types=lambda wc, input: n_org_types(merged_data_path(wc.species)),
        max_org_molid=lambda wc, input: max_org_molid(merged_data_path(wc.species)),
        cutoff=config["slab"]["merge_cutoff_angstrom"],
        kspace_accuracy=config["slab"]["merge_kspace_accuracy"],
        timestep_fs=config["slab"]["merge_nvt_timestep_fs"],
        temperature_K=config["temperature_K"],
        read_data="merged.data",
        write_prefix="minimized",
        velocity_create=True,
        seed=config["slab"]["merge_nvt_seed"],
        etol=config["slab"]["merge_minimize"]["etol"],
        ftol=config["slab"]["merge_minimize"]["ftol"],
        maxiter=config["slab"]["merge_minimize"]["maxiter"],
        maxeval=config["slab"]["merge_minimize"]["maxeval"],
    wildcard_constraints:
        species="ST|MMA",
    localrule: True
    script:
        "../scripts/render_template.py"


rule lmp_slab_minimize:
    input:
        lmp=f"{SLAB_RES}/{{species}}/minimize.lmp",
        data=f"{SLAB_RES}/{{species}}/merged.data",
        coeffs=f"{SLAB_RES}/{{species}}/merged_coeffs.lmp",
    output:
        data=f"{SLAB_RES}/{{species}}/minimized.data",
        restart=f"{SLAB_RES}/{{species}}/minimized.restart",
    params:
        omp_threads=lambda wc, resources: resources.get("ompthreads", 1),
    shell:
        'cd "$(dirname {input.lmp})" && '
        'LMP_BIN="$(command -v lmp_mpi || command -v lmp)" && '
        '"$LMP_BIN" -sf omp -pk omp {params.omp_threads} -in "$(basename {input.lmp})"'


# --- stage 4b: z-direction density equilibration -------------------------

rule render_npt_z_lmp:
    input:
        "lmp/npt_z.lmp.template",
        LMP_COMMON_TEMPLATE,
    output:
        f"{SLAB_RES}/{{species}}/npt_z.lmp",
    params:
        n_org_types=lambda wc, input: n_org_types(merged_data_path(wc.species)),
        max_org_molid=lambda wc, input: max_org_molid(merged_data_path(wc.species)),
        cutoff=config["slab"]["merge_cutoff_angstrom"],
        kspace_accuracy=config["slab"]["merge_kspace_accuracy"],
        timestep_fs=config["slab"]["merge_nvt_timestep_fs"],
        temperature_K=config["temperature_K"],
        read_data="minimized.data",
        pressure_atm=config["slab"]["npt_pressure_atm"],
        pressure_damp_fs=config["slab"]["npt_pressure_damp_fs"],
        modulus_atm=config["slab"]["npt_modulus_atm"],
        npt_z_nsteps=lambda wc: int(config["slab"]["npt_z_ps"] * 1000 / config["slab"]["merge_nvt_timestep_fs"]),
    wildcard_constraints:
        species="ST|MMA",
    localrule: True
    script:
        "../scripts/render_template.py"


rule lmp_npt_z:
    input:
        lmp=f"{SLAB_RES}/{{species}}/npt_z.lmp",
        data=f"{SLAB_RES}/{{species}}/minimized.data",
        coeffs=f"{SLAB_RES}/{{species}}/merged_coeffs.lmp",
    output:
        data=f"{SLAB_RES}/{{species}}/npt_equilibrated.data",
        restart=f"{SLAB_RES}/{{species}}/npt_equilibrated.restart",
    params:
        omp_threads=lambda wc, resources: resources.get("ompthreads", 1),
    shell:
        'cd "$(dirname {input.lmp})" && '
        'LMP_BIN="$(command -v lmp_mpi || command -v lmp)" && '
        '"$LMP_BIN" -sf omp -pk omp {params.omp_threads} -in "$(basename {input.lmp})"'


# --- stage 4c: move the fragment to the bulk-arc midpoint, drop whatever
# hexadecane it would otherwise overlap (Python, not LAMMPS -- see
# scripts/insert_fragment.py) ---------------------------------------------

rule insert_fragment:
    input:
        data=f"{SLAB_RES}/{{species}}/npt_equilibrated.data",
        script="scripts/insert_fragment.py",
    output:
        data=f"{SLAB_RES}/{{species}}/inserted.data",
    params:
        n_sub_types=N_SUB_TYPES,
        clash_angstrom=config["slab"]["fragment_insert_clash_angstrom"],
    wildcard_constraints:
        species="ST|MMA",
    conda:
        "../envs/analysis.yaml"
    localrule: True
    shell:
        "python {input.script} {input.data} -o {output.data} "
        "--n-sub-types {params.n_sub_types} --clash-angstrom {params.clash_angstrom}"


# --- stage 4c-ii: relax whatever local clash remains after insertion -----

rule render_insert_minimize_lmp:
    input:
        "lmp/minimize.lmp.template",
        LMP_COMMON_TEMPLATE,
    output:
        f"{SLAB_RES}/{{species}}/insert_minimize.lmp",
    params:
        n_org_types=lambda wc, input: n_org_types(merged_data_path(wc.species)),
        max_org_molid=lambda wc, input: max_org_molid(merged_data_path(wc.species)),
        cutoff=config["slab"]["merge_cutoff_angstrom"],
        kspace_accuracy=config["slab"]["merge_kspace_accuracy"],
        timestep_fs=config["slab"]["merge_nvt_timestep_fs"],
        temperature_K=config["temperature_K"],
        read_data="inserted.data",
        write_prefix="insert_minimized",
        velocity_create=False,
        seed=config["slab"]["merge_nvt_seed"],
        etol=config["slab"]["merge_minimize"]["etol"],
        ftol=config["slab"]["merge_minimize"]["ftol"],
        maxiter=config["slab"]["merge_minimize"]["maxiter"],
        maxeval=config["slab"]["merge_minimize"]["maxeval"],
    wildcard_constraints:
        species="ST|MMA",
    localrule: True
    script:
        "../scripts/render_template.py"


rule lmp_insert_minimize:
    input:
        lmp=f"{SLAB_RES}/{{species}}/insert_minimize.lmp",
        data=f"{SLAB_RES}/{{species}}/inserted.data",
        coeffs=f"{SLAB_RES}/{{species}}/merged_coeffs.lmp",
    output:
        data=f"{SLAB_RES}/{{species}}/insert_minimized.data",
        restart=f"{SLAB_RES}/{{species}}/insert_minimized.restart",
    params:
        omp_threads=lambda wc, resources: resources.get("ompthreads", 1),
    shell:
        'cd "$(dirname {input.lmp})" && '
        'LMP_BIN="$(command -v lmp_mpi || command -v lmp)" && '
        '"$LMP_BIN" -sf omp -pk omp {params.omp_threads} -in "$(basename {input.lmp})"'


# --- stage 4d: short restrained re-equilibration --------------------------

rule render_restrained_npt_lmp:
    input:
        "lmp/restrained_npt.lmp.template",
        LMP_COMMON_TEMPLATE,
    output:
        f"{SLAB_RES}/{{species}}/restrained_npt.lmp",
    params:
        n_org_types=lambda wc, input: n_org_types(merged_data_path(wc.species)),
        max_org_molid=lambda wc, input: max_org_molid(merged_data_path(wc.species)),
        cutoff=config["slab"]["merge_cutoff_angstrom"],
        kspace_accuracy=config["slab"]["merge_kspace_accuracy"],
        timestep_fs=config["slab"]["merge_nvt_timestep_fs"],
        temperature_K=config["temperature_K"],
        read_data="insert_minimized.data",
        pressure_atm=config["slab"]["npt_pressure_atm"],
        pressure_damp_fs=config["slab"]["npt_pressure_damp_fs"],
        modulus_atm=config["slab"]["npt_modulus_atm"],
        restraint_k=config["slab"]["fragment_restraint_k_kcal_per_mol_A2"],
        restrained_npt_nsteps=lambda wc: int(
            config["slab"]["restrained_npt_ps"] * 1000 / config["slab"]["merge_nvt_timestep_fs"]
        ),
    wildcard_constraints:
        species="ST|MMA",
    localrule: True
    script:
        "../scripts/render_template.py"


rule lmp_restrained_npt:
    input:
        lmp=f"{SLAB_RES}/{{species}}/restrained_npt.lmp",
        data=f"{SLAB_RES}/{{species}}/insert_minimized.data",
        coeffs=f"{SLAB_RES}/{{species}}/merged_coeffs.lmp",
    output:
        data=f"{SLAB_RES}/{{species}}/restrained.data",
        restart=f"{SLAB_RES}/{{species}}/restrained.restart",
    params:
        omp_threads=lambda wc, resources: resources.get("ompthreads", 1),
    shell:
        'cd "$(dirname {input.lmp})" && '
        'LMP_BIN="$(command -v lmp_mpi || command -v lmp)" && '
        '"$LMP_BIN" -sf omp -pk omp {params.omp_threads} -in "$(basename {input.lmp})"'


# --- stage 4e: final NVT, restraint released ------------------------------

rule render_final_nvt_lmp:
    input:
        "lmp/final_nvt.lmp.template",
        LMP_COMMON_TEMPLATE,
    output:
        f"{SLAB_RES}/{{species}}/final_nvt.lmp",
    params:
        n_org_types=lambda wc, input: n_org_types(merged_data_path(wc.species)),
        max_org_molid=lambda wc, input: max_org_molid(merged_data_path(wc.species)),
        cutoff=config["slab"]["merge_cutoff_angstrom"],
        kspace_accuracy=config["slab"]["merge_kspace_accuracy"],
        timestep_fs=config["slab"]["merge_nvt_timestep_fs"],
        temperature_K=config["temperature_K"],
        read_data="restrained.data",
        nsteps=lambda wc: int(config["slab"]["merge_nvt_ps"] * 1000 / config["slab"]["merge_nvt_timestep_fs"]),
    wildcard_constraints:
        species="ST|MMA",
    localrule: True
    script:
        "../scripts/render_template.py"


rule lmp_final_nvt:
    input:
        lmp=f"{SLAB_RES}/{{species}}/final_nvt.lmp",
        data=f"{SLAB_RES}/{{species}}/restrained.data",
        coeffs=f"{SLAB_RES}/{{species}}/merged_coeffs.lmp",
    output:
        data=f"{SLAB_RES}/{{species}}/nvt.data",
        restart=f"{SLAB_RES}/{{species}}/nvt.restart",
    params:
        omp_threads=lambda wc, resources: resources.get("ompthreads", 1),
    shell:
        'cd "$(dirname {input.lmp})" && '
        'LMP_BIN="$(command -v lmp_mpi || command -v lmp)" && '
        '"$LMP_BIN" -sf omp -pk omp {params.omp_threads} -in "$(basename {input.lmp})"'


# --- stage 5: production, with z-resolved density profiles ---------------

rule render_production_lmp:
    input:
        "lmp/production.lmp.template",
        LMP_COMMON_TEMPLATE,
    output:
        f"{SLAB_RES}/{{species}}/production.lmp",
    params:
        n_org_types=lambda wc, input: n_org_types(merged_data_path(wc.species)),
        max_org_molid=lambda wc, input: max_org_molid(merged_data_path(wc.species)),
        cutoff=config["slab"]["merge_cutoff_angstrom"],
        kspace_accuracy=config["slab"]["merge_kspace_accuracy"],
        timestep_fs=config["slab"]["merge_nvt_timestep_fs"],
        temperature_K=config["temperature_K"],
        read_data="nvt.data",
        production_nsteps=lambda wc: int(
            config["slab"]["production_ns"] * 1.0e6 / config["slab"]["merge_nvt_timestep_fs"]
        ),
        bin_width_angstrom=config["slab"]["density_bin_angstrom"],
        chunk_nevery=config["slab"]["density_chunk_nevery"],
        chunk_nrepeat=config["slab"]["density_chunk_nrepeat"],
        chunk_nfreq=config["slab"]["density_chunk_nfreq"],
    wildcard_constraints:
        species="ST|MMA",
    localrule: True
    script:
        "../scripts/render_template.py"


rule lmp_production:
    input:
        lmp=f"{SLAB_RES}/{{species}}/production.lmp",
        data=f"{SLAB_RES}/{{species}}/nvt.data",
        coeffs=f"{SLAB_RES}/{{species}}/merged_coeffs.lmp",
    output:
        data=f"{SLAB_RES}/{{species}}/production.data",
        restart=f"{SLAB_RES}/{{species}}/production.restart",
        fragment_density=f"{SLAB_RES}/{{species}}/fragment_density.profile",
        hexadecane_density=f"{SLAB_RES}/{{species}}/hexadecane_density.profile",
    params:
        omp_threads=lambda wc, resources: resources.get("ompthreads", 1),
    shell:
        'cd "$(dirname {input.lmp})" && '
        'LMP_BIN="$(command -v lmp_mpi || command -v lmp)" && '
        '"$LMP_BIN" -sf omp -pk omp {params.omp_threads} -in "$(basename {input.lmp})"'
