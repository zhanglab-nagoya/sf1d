"""Step 2 chis stages 4+5: LAMMPS minimize + NVT equilibration + production
of the merged substrate/liquid slab (results/slab/{species}/merged.data,
step2_merge.smk), in one continuous LAMMPS session -- see
lmp/equilibrate_merge.lmp.template's module docstring for why this is
NOT split into separate minimize/NVT/production rules: `bond_style hybrid
harmonic`'s state didn't survive a write_data/read_data OR
write_restart/read_restart round trip reliably between two separate
`lmp` invocations (confirmed live 2026-11-06), and keeping every stage in
one session sidesteps that whole problem.

Every stage freezes the substrate (already relaxed on its own terms by
step2_substrate.smk and not re-relaxed or thermally perturbed here --
this pipeline treats it as a static wall, matching sf1d's own SCF
model).

Stage 5's production run is long (slab.production_ns, 2026-11-06 default
10 ns/species) and genuinely needs real parallelism -- confirmed live via
a throughput benchmark on flow2-ai (profiles/flow2-ai/, job 102982.jobm):
48 OpenMP threads (1 NUMA, t1-share) gave 17.4x over this workstation's
serial baseline (0.111 s/step -> 0.00636 s/step), 6.8 ns/day/species. The
`omp_threads` param below defaults to 1 locally (matching the originally
verified local invocation) and is meant to be driven up via the
`ompthreads` resource on a cluster profile; `-sf omp -pk omp N` applies
the OPENMP package's threaded pair/kspace/bond/angle/dihedral styles
(`Last active /omp style is ...` in the log confirms it actually switched
styles, not just set an unused env var).

STATUS: stage 4 DONE, run end-to-end for both ST and MMA 2026-11-06:
minimization converges cleanly on the energy-tolerance criterion in ~220
steps (no NaN/blow-up), and NVT holds the mobile liquid at the target
temperature (see lmp/equilibrate_merge.lmp.template's comment on why the
thermo `temp` column itself reads low -- a display-only artifact of
counting the frozen substrate's DOF, not a thermostat problem). Stage 5
is scaffolded (density-profile chunking added to the same template) but
not yet run for the real production_ns -- that run is expensive enough
(~35 h/species on flow2-ai) that it should go through the cluster
profile, not this workstation.
"""

SLAB_RES = "results/slab"  # shared with step2_dilute.smk/convert/merge
N_SUB_TYPES = len(config["slab"]["substrate_creation"]["masses"])  # Fe, O


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


rule render_equilibrate_merge_lmp:
    input:
        "lmp/equilibrate_merge.lmp.template",
    output:
        f"{SLAB_RES}/{{species}}/equilibrate_merge.lmp",
    params:
        n_org_types=lambda wc, input: n_org_types(f"{SLAB_RES}/{wc.species}/merged.data"),
        max_org_molid=lambda wc, input: max_org_molid(f"{SLAB_RES}/{wc.species}/merged.data"),
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
        production_data=f"{SLAB_RES}/{{species}}/production.data",
        production_restart=f"{SLAB_RES}/{{species}}/production.restart",
        fragment_density=f"{SLAB_RES}/{{species}}/fragment_density.profile",
        hexadecane_density=f"{SLAB_RES}/{{species}}/hexadecane_density.profile",
    params:
        omp_threads=lambda wc, resources: resources.get("ompthreads", 1),
    shell:
        'cd "$(dirname {input.lmp})" && '
        'LMP_BIN="$(command -v lmp_mpi || command -v lmp)" && '
        '"$LMP_BIN" -sf omp -pk omp {params.omp_threads} -in "$(basename {input.lmp})"'
