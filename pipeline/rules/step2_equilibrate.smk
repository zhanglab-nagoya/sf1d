"""Step 2 chis stages 4+5: LAMMPS minimize + z-NPT equilibration +
fragment placement + restrained re-equilibration + NVT + production of
the merged substrate/liquid slab (results/slab/{species}/merged.data,
step2_merge.smk), in one continuous LAMMPS session -- see
lmp/equilibrate_merge.lmp.template's module docstring for why this is
NOT split into separate rules: `bond_style hybrid harmonic`'s state
didn't survive a write_data/read_data OR write_restart/read_restart
round trip reliably between two separate `lmp` invocations (confirmed
live 2026-11-06), and keeping every stage in one session sidesteps that
whole problem -- including the fragment's repositioning, done with
LAMMPS's own xcm()/displace_atoms rather than a second Python pass.

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

STATUS: rewritten 2026-11-08 after the first real 10 ns/species
production run (on flow2-ai) revealed a real ~10 A vacuum gap for ST and
an asymmetric density profile for MMA -- traced to the old fixed-volume
NVT stage never actually being at this merged system's true equilibrium
density (merge_slab.py sized the box's z-length from a single raw,
unequilibrated snapshot), and to the single fragment molecule -- started
wherever step2_dilute.smk's GROMACS packing happened to leave it, not a
deliberate placement -- getting swept to the box's periodic seam by that
same under-equilibrated transient within the first 20 ps of the old NVT.
Stages 4b-4d below (z-only barostat equilibration, deliberate fragment
placement at the equilibrated bulk-arc midpoint, short restrained
re-equilibration) are new. The barostat went through two real failures
before landing on `fix press/berendsen` (see
lmp/equilibrate_merge.lmp.template's module docstring) -- `fix npt`
(Nose-Hoover) never converges on this system regardless of target,
confirmed live 2026-11-08 on two independent attempts. Smoke-tested
locally end-to-end 2026-11-08 on ST's real merged.data (tiny step counts
-- 2 ps NPT-z, 1 ps restrained, 1 ps final NVT, 1 ps production): zero
LAMMPS errors, temperature stable ~295-305 K through every stage, and
the fragment -- deliberately placed at the equilibrated bulk-arc
midpoint (z~68.75) -- was still sitting in that same bulk region
(z~65-83, still mobile, not stuck) by the end, nowhere near the
periodic seam. Not yet run at the real production_ns or on flow2-ai;
that's the next step, pending the user's go-ahead given the point
budget.
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


def substrate_z_range(merged_data_path):
    """(sub_zmin, sub_zmax): the substrate's own static z-extent, read
    directly from merged.data's Atoms section. Unaffected by any later
    dynamics -- the substrate is always frozen -- so this is passed into
    the template as a fixed param; stage 4c computes the equilibrated
    liquid arc's midpoint from these plus the running session's OWN
    zhi/zlo (which can differ from merged.data's after stage 4b's
    barostat), not from a second Python pass over a checkpoint file."""
    n_org = n_org_types(merged_data_path)
    zs = []
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
                if int(fields[2]) > n_org:
                    zs.append(float(fields[6]))
    if not zs:
        raise ValueError(f"no substrate atoms (type > {n_org}) found in {merged_data_path}")
    return min(zs), max(zs)


rule render_equilibrate_merge_lmp:
    input:
        "lmp/equilibrate_merge.lmp.template",
    output:
        f"{SLAB_RES}/{{species}}/equilibrate_merge.lmp",
    params:
        n_org_types=lambda wc, input: n_org_types(f"{SLAB_RES}/{wc.species}/merged.data"),
        max_org_molid=lambda wc, input: max_org_molid(f"{SLAB_RES}/{wc.species}/merged.data"),
        sub_zmin=lambda wc, input: substrate_z_range(f"{SLAB_RES}/{wc.species}/merged.data")[0],
        sub_zmax=lambda wc, input: substrate_z_range(f"{SLAB_RES}/{wc.species}/merged.data")[1],
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
        pressure_atm=config["slab"]["npt_pressure_atm"],
        pressure_damp_fs=config["slab"]["npt_pressure_damp_fs"],
        modulus_atm=config["slab"]["npt_modulus_atm"],
        npt_z_nsteps=lambda wc: int(config["slab"]["npt_z_ps"] * 1000 / config["slab"]["merge_nvt_timestep_fs"]),
        restrained_npt_nsteps=lambda wc: int(
            config["slab"]["restrained_npt_ps"] * 1000 / config["slab"]["merge_nvt_timestep_fs"]
        ),
        restraint_k=config["slab"]["fragment_restraint_k_kcal_per_mol_A2"],
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
        npt_equilibrated_data=f"{SLAB_RES}/{{species}}/npt_equilibrated.data",
        npt_equilibrated_restart=f"{SLAB_RES}/{{species}}/npt_equilibrated.restart",
        restrained_data=f"{SLAB_RES}/{{species}}/restrained.data",
        restrained_restart=f"{SLAB_RES}/{{species}}/restrained.restart",
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
