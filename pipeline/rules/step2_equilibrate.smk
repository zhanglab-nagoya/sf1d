"""Step 2 chis stages 4+5: LAMMPS minimize + z-direction density
equilibration of the solvent+substrate system, fragment insertion, short
restrained relaxation, and production -- rewritten 2026-11-13 to decouple
the fragment from the solvent/substrate equilibration entirely (see the
"Decouple fragment insertion from solvent/substrate equilibration" plan
for the full root-cause history this STATUS section condenses).

One Snakemake rule per stage, each its own small `lmp/*.lmp.template` file
sharing a common preamble (`lmp/_common_setup.lmp.template`, `{% include %}`d
by every stage template).

Two halves, split the same way step2_dilute.smk/step2_convert.smk/
step2_merge.smk already are:

1. SPECIES-INDEPENDENT (results/slab/solvent/, no `{species}` wildcard):
   lmp_slab_minimize -> lmp_nvt -> close_liquid_gap -> lmp_gap_minimize ->
   lmp_settle, run on the solvent+substrate system alone (no fragment
   present anywhere in this half). `_common_setup.lmp.template`'s
   fragment/hexadecane molid-range groups are never defined here --
   "mobile" IS all of hexadecane. Built and run exactly ONCE, reused by
   both ST and MMA.

   This used to carry a single restrained fragment through the whole
   chain (`spring/self`, released only at production). Confirmed live
   2026-11-12/13 across three independently-tuned barostat variants
   (isotropic-average vs. Pzz-only pressure formula; naive vs.
   mass-conservation-based initial box sizing) that they ALL converged to
   the same wrong ~271 A/~1500 atm plateau regardless of starting
   conditions -- the real mechanism: `change_box hexadecane z ... remap`
   rescales hexadecane proportionally from a fixed zlo every cycle, but
   the restrained fragment never moves, so hexadecane downstream of it
   gets pulled progressively further away over a long (100 ps/400-cycle)
   run, opening a real, permanent, growing cavity right at the fragment's
   location (confirmed live: density crashed from ~150-270 to 0 within
   ~10 A of the fragment's own edge in the final npt_equilibrated.data).
   Because `hexadecane_volume` is read from that group's own `zmax-zmin`
   span, the empty cavity got folded into the barostat's own pressure
   signal, letting it settle on a self-consistent but wrong fixed point.
   Removing the fragment from this stage entirely removes the anchor the
   mismatch depended on.

   `lmp_npt_z`/the z-direction barostat (both the original whole-box-
   pressure version and a later hand-rolled, liquid-only-pressure
   Berendsen rescale) was removed entirely 2026-11-13, replaced by
   `lmp_settle` (plain, box-fixed NVT, no pressure feedback at all -- see
   config.yaml's slab.settle_ps comment for the full history). This
   wasn't a tuning fix, it was dropping the barostat: direct empirical
   testing established the hand-rolled version's convergence rate is
   governed by the liquid's own real compressibility (an effective ~120
   ps relaxation time) regardless of clamp/modulus tuning, so closing the
   ~17,000-20,000 atm gap left by close_liquid_gap.py's rigid-body
   translation would cost roughly 1 ns, not the ~100 ps budgeted --  and
   researching how ../2026-07-28-fe2o3-cg solves the identical
   frozen-substrate-plus-solvent problem confirmed this is a structural
   dead end there too, not a tuning gap specific to our own
   implementation. Its fix -- a deterministic box resize to a reference-
   density target followed by a barostat-free settle, no pressure
   feedback anywhere -- is what `close_liquid_gap.py` (resize) +
   `lmp_gap_minimize` (local-overlap relax) + `lmp_settle` (heating
   dissipation) already amounts to once the barostat step between the
   last two is simply removed.

2. SPECIES-SPECIFIC (results/slab/{species}/): `insert_fragment`
   (scripts/insert_fragment.py, recreated 2026-11-13 in a new design --
   see its own docstring) combines the converged
   results/slab/solvent/npt_equilibrated.data with the species' own
   standalone results/slab/{species}/fragment_converted.lmp
   (step2_convert.smk's convert_fragment_to_lammps), finding the bulk
   liquid arc's midpoint, translating the fragment there, and deleting
   whichever hexadecane molecules it clashes with (whole molecules only,
   never partial -- would leave a dangling bond). `lmp_insert_minimize`
   (lmp/minimize.lmp.template, `velocity_create_group="fragment"` so
   hexadecane's own already-equilibrated velocities are untouched)
   relaxes the resulting local clash. `lmp_insert_relax`
   (lmp/nvt.lmp.template, `has_fragment=True`, box fixed, short --
   slab.insert_relax_ps) lets the local disturbance thermalize with the
   fragment restrained (`spring/self`, same mechanism as before),
   released only at the start of `lmp_production`.

   Deliberately NO barostat stage after insertion: the box is already at
   its correct equilibrium size from stage 1 above, and a barostat here
   would reopen the same structural problem removing `lmp_npt_z` from
   stage 1 exists to avoid (see this file's own `settle_ps` comment and
   config.yaml's slab.settle_ps comment for the full history -- it isn't
   specific to the restrained-fragment case, a barostat on this frozen-
   substrate system is unreliable regardless of what else is present).
   If the short NVT alone doesn't fully relax the local clash-deletion
   disturbance (watch the hexadecane density profile right around the
   fragment's own position in relaxed.data/production.data -- the same
   check that originally surfaced the old cavity bug), the fallback is a
   longer `insert_relax_ps`, not a barostat stage of any kind.

Every stage freezes the substrate (already relaxed on its own terms by
step2_substrate.smk and not re-relaxed or thermally perturbed here --
this pipeline treats it as a static wall, matching sf1d's own SCF model).

Every `merged_coeffs.lmp` (solvent/, and per-species once insert_fragment
has run) carries explicit `pair_coeff`/`bond_coeff`/`angle_coeff`/
`dihedral_coeff` commands (see
scripts/merge_slab.py's own docstring for why: `write_data` silently drops
hybrid-style Bond/Angle/Dihedral Coeffs sections across a
write_data/read_data round trip, confirmed live 2026-11-06). Every stage
template `include`s the relevant sidecar right after its own `read_data`.

NVT (thermalize at fixed volume) runs BEFORE close_liquid_gap's
deterministic resize and lmp_settle's post-resize settle in the
solvent-only chain -- thermalizing before resizing keeps the gap-closing
math working from a properly-equilibrated (not freshly-minimized, still
cold) starting configuration.

Production is long (slab.production_ns, default 10 ns/species) and
genuinely needs real parallelism -- confirmed live via a throughput
benchmark on flow2-ai (profiles/flow2-ai/, job 102982.jobm): 48 OpenMP
threads (1 NUMA, t1-share) gave 17.4x over this workstation's serial
baseline, 6.8 ns/day/species. Every LAMMPS rule below takes the same
`omp_threads`/`-sf omp -pk omp N` pattern (not just production), defaulting
to 1 locally and meant to be driven up per-rule via the `ompthreads`
resource on a cluster profile.

STATUS: the solvent/fragment split above has not yet been run end-to-end
(2026-11-13). Prior to this split, the per-stage rule chain (fragment
restrained throughout, as described in point 1 above) ran successfully
for real (production_ns=10, both species) on flow2-ai, completing
2026-11-11 (ST job 108984.jobm, MMA job 108985.jobm) -- but that run and
subsequent local verification surfaced the fragment-restraint cavity bug
this split fixes, plus earlier, now-fixed bugs along the way: atom types
were never deduplicated (scripts/dedup_lammps_types.py now also collapses
atom types, not just bond/angle/dihedral types); the z-barostat was
originally driven by LAMMPS's whole-box `compute pressure` (hard-requires
group "all", so it always included the frozen substrate's own enormous
internal virial); and -- discovered verifying the solvent/fragment split
itself, 2026-11-13 -- even a liquid-only-pressure hand-rolled barostat
replacement for that was itself a structural dead end on this
frozen-substrate system, not a tuning problem (see this file's own
`lmp_settle` comment above for the full history). `lmp_npt_z` and
`lmp/_liquid_barostat.lmp.template`/`lmp/npt_z.lmp.template` are removed
entirely; `lmp_settle` (plain NVT, no barostat) replaces that stage.
Local verification of the solvent-only chain through `lmp_settle` is
in progress as of this writing. Both species' full production needs
rerunning once this split is verified end-to-end (see the plan's own
Verification section for the checklist).
"""

SLAB_RES = "results/slab"
SOLVENT_RES = f"{SLAB_RES}/solvent"  # species-independent, shared by ST/MMA
N_SUB_TYPES = len(config["slab"]["substrate_creation"]["masses"])  # Fe, O
LMP_COMMON_TEMPLATE = "lmp/_common_setup.lmp.template"
SOLVENT_MERGED_DATA = f"{SOLVENT_RES}/merged.data"


def n_org_types(data_path):
    with open(data_path) as fh:
        for line in fh:
            if "atom types" in line:
                return int(line.split()[0]) - N_SUB_TYPES
    raise ValueError(f"no 'atom types' line found in {data_path}")


def max_org_molid(data_path):
    """merge_slab.py never renumbers organic molecule-IDs and gives the
    substrate the single id (max organic molid + 1) -- so the largest
    molid anywhere in a solvent-only file, minus 1, is exactly the
    hexadecane side's own max, with no need to separately locate the
    substrate's atoms. NOT valid on inserted.data (which has the
    fragment's own molids beyond substrate's too) -- use
    insert_molid_info() for that instead."""
    max_molid = 0
    in_atoms = False
    with open(data_path) as fh:
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
        raise ValueError(f"no Atoms section found in {data_path}")
    return max_molid - 1


def insert_molid_info(inserted_data_path):
    """inserted.data's molid layout (scripts/insert_fragment.py): hexadecane
    keeps its own [1..H] unchanged, substrate keeps its existing [H+1], the
    fragment's 3 residues (cap/core/cap, see build_fragment_itp) are
    appended last at [H+2, H+3, H+4]. Returns (hex_max_molid, frag_lo,
    frag_hi) derived from the file's own overall max molid M = H+4."""
    overall_max = 0
    in_atoms = False
    with open(inserted_data_path) as fh:
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
                overall_max = max(overall_max, int(fields[1]))
    if overall_max == 0:
        raise ValueError(f"no Atoms section found in {inserted_data_path}")
    frag_hi = overall_max
    frag_lo = overall_max - 2
    hex_max_molid = overall_max - 4
    return hex_max_molid, frag_lo, frag_hi


def inserted_data_path(species):
    return f"{SLAB_RES}/{species}/inserted.data"


# ===========================================================================
# Part 1: species-independent solvent+substrate equilibration
# ===========================================================================

# --- stage 4a: minimize --------------------------------------------------

rule render_slab_minimize_lmp:
    input:
        "lmp/minimize.lmp.template",
        LMP_COMMON_TEMPLATE,
        # merged.data is only read ad hoc by the params lambdas below
        # (n_org_types) -- but it MUST be declared here too, not just as a
        # param side-effect, or Snakemake's DAG doesn't know this rule
        # depends on it and can schedule it before merge_solvent_with_
        # substrate has actually written the file it's about to read.
        SOLVENT_MERGED_DATA,
    output:
        f"{SOLVENT_RES}/minimize.lmp",
    params:
        n_org_types=lambda wc, input: n_org_types(SOLVENT_MERGED_DATA),
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
        soft_cutoff=config["slab"]["soft_premin"]["soft_cutoff"],
        soft_a=config["slab"]["soft_premin"]["soft_a"],
        soft_etol=config["slab"]["soft_premin"]["soft_etol"],
        soft_ftol=config["slab"]["soft_premin"]["soft_ftol"],
        soft_maxiter=config["slab"]["soft_premin"]["soft_maxiter"],
        soft_maxeval=config["slab"]["soft_premin"]["soft_maxeval"],
    localrule: True
    script:
        "../scripts/render_template.py"


rule lmp_slab_minimize:
    input:
        lmp=f"{SOLVENT_RES}/minimize.lmp",
        data=SOLVENT_MERGED_DATA,
        coeffs=f"{SOLVENT_RES}/merged_coeffs.lmp",
    output:
        data=f"{SOLVENT_RES}/minimized.data",
        restart=f"{SOLVENT_RES}/minimized.restart",
    params:
        omp_threads=lambda wc, resources: resources.get("ompthreads", 1),
    shell:
        'cd "$(dirname {input.lmp})" && '
        'LMP_BIN="$(command -v lmp_mpi || command -v lmp)" && '
        '"$LMP_BIN" -sf omp -pk omp {params.omp_threads} -in "$(basename {input.lmp})"'


# --- stage 4b: thermalize (NVT), box fixed, no fragment -----------------

rule render_nvt_lmp:
    input:
        "lmp/nvt.lmp.template",
        LMP_COMMON_TEMPLATE,
        SOLVENT_MERGED_DATA,
    output:
        f"{SOLVENT_RES}/nvt.lmp",
    params:
        n_org_types=lambda wc, input: n_org_types(SOLVENT_MERGED_DATA),
        cutoff=config["slab"]["merge_cutoff_angstrom"],
        kspace_accuracy=config["slab"]["merge_kspace_accuracy"],
        timestep_fs=config["slab"]["merge_nvt_timestep_fs"],
        temperature_K=config["temperature_K"],
        read_data="minimized.data",
        write_prefix="nvt",
        has_fragment=False,
        nsteps=lambda wc: int(config["slab"]["merge_nvt_ps"] * 1000 / config["slab"]["merge_nvt_timestep_fs"]),
    localrule: True
    script:
        "../scripts/render_template.py"


rule lmp_nvt:
    input:
        lmp=f"{SOLVENT_RES}/nvt.lmp",
        data=f"{SOLVENT_RES}/minimized.data",
        coeffs=f"{SOLVENT_RES}/merged_coeffs.lmp",
    output:
        data=f"{SOLVENT_RES}/nvt.data",
        restart=f"{SOLVENT_RES}/nvt.restart",
    params:
        omp_threads=lambda wc, resources: resources.get("ompthreads", 1),
    shell:
        'cd "$(dirname {input.lmp})" && '
        'LMP_BIN="$(command -v lmp_mpi || command -v lmp)" && '
        '"$LMP_BIN" -sf omp -pk omp {params.omp_threads} -in "$(basename {input.lmp})"'


# --- stage 4b-ii: deterministically close the vacuum-like gap NVT leaves
# in the liquid, before NPT's barostat ever sees it (scripts/close_liquid_
# gap.py) -- same pipeline position as ../2026-10-28-generate-aa-pama's own
# shrink_box_zhi, between its NVT and NPT stages ------------------------

rule close_liquid_gap:
    input:
        data=f"{SOLVENT_RES}/nvt.data",
        segment_volume="results/segment_volume.yaml",
        script="scripts/close_liquid_gap.py",
    output:
        data=f"{SOLVENT_RES}/nvt_shrunk.data",
    params:
        n_sub_types=N_SUB_TYPES,
    localrule: True
    shell:
        "python {input.script} {input.data} -o {output.data} "
        "--n-sub-types {params.n_sub_types} --segment-volume {input.segment_volume}"


# --- stage 4b-iii: relax whatever local clash the deterministic gap
# closure leaves behind -- it translates each molecule as a rigid body
# (zero internal distortion), but two independently-shifted, long/floppy
# hexadecane chains can still land with overlapping conformations close
# to each other ------------------------------------------------------------

rule render_gap_minimize_lmp:
    input:
        "lmp/minimize.lmp.template",
        LMP_COMMON_TEMPLATE,
        SOLVENT_MERGED_DATA,
    output:
        f"{SOLVENT_RES}/gap_minimize.lmp",
    params:
        n_org_types=lambda wc, input: n_org_types(SOLVENT_MERGED_DATA),
        cutoff=config["slab"]["merge_cutoff_angstrom"],
        kspace_accuracy=config["slab"]["merge_kspace_accuracy"],
        timestep_fs=config["slab"]["merge_nvt_timestep_fs"],
        temperature_K=config["temperature_K"],
        read_data="nvt_shrunk.data",
        write_prefix="nvt_shrunk_minimized",
        velocity_create=False,
        seed=config["slab"]["merge_nvt_seed"],
        etol=config["slab"]["merge_minimize"]["etol"],
        ftol=config["slab"]["merge_minimize"]["ftol"],
        maxiter=config["slab"]["merge_minimize"]["maxiter"],
        maxeval=config["slab"]["merge_minimize"]["maxeval"],
        soft_cutoff=config["slab"]["soft_premin"]["soft_cutoff"],
        soft_a=config["slab"]["soft_premin"]["soft_a"],
        soft_etol=config["slab"]["soft_premin"]["soft_etol"],
        soft_ftol=config["slab"]["soft_premin"]["soft_ftol"],
        soft_maxiter=config["slab"]["soft_premin"]["soft_maxiter"],
        soft_maxeval=config["slab"]["soft_premin"]["soft_maxeval"],
    localrule: True
    script:
        "../scripts/render_template.py"


rule lmp_gap_minimize:
    input:
        lmp=f"{SOLVENT_RES}/gap_minimize.lmp",
        data=f"{SOLVENT_RES}/nvt_shrunk.data",
        coeffs=f"{SOLVENT_RES}/merged_coeffs.lmp",
    output:
        data=f"{SOLVENT_RES}/nvt_shrunk_minimized.data",
        restart=f"{SOLVENT_RES}/nvt_shrunk_minimized.restart",
    params:
        omp_threads=lambda wc, resources: resources.get("ompthreads", 1),
    shell:
        'cd "$(dirname {input.lmp})" && '
        'LMP_BIN="$(command -v lmp_mpi || command -v lmp)" && '
        '"$LMP_BIN" -sf omp -pk omp {params.omp_threads} -in "$(basename {input.lmp})"'


# --- stage 4c: barostat-free settle at the already-correct density ------
# (see config.yaml's slab.settle_ps comment for the full root-cause
# history of why this replaces a z-direction barostat here: a barostat
# acting on a system with a large frozen, immobile substrate component is
# fundamentally unreliable, not just hard to tune -- confirmed by
# ../2026-07-28-fe2o3-cg's own identical finding on the identical
# substrate+solvent pairing. scripts/close_liquid_gap.py already resized
# the box deterministically to the solvent's own reference bulk density
# (mass-conservation based), and lmp_gap_minimize's soft-core pre-relax
# already resolved the local overlap that rigid-body translation leaves
# -- this stage is just a short, fixed-volume NVT to dissipate whatever
# compression heating remains, reusing lmp/nvt.lmp.template's
# has_fragment=False branch (box fixed, no deform, no barostat) rather
# than a dedicated template.

rule render_settle_lmp:
    input:
        "lmp/nvt.lmp.template",
        LMP_COMMON_TEMPLATE,
        SOLVENT_MERGED_DATA,
    output:
        f"{SOLVENT_RES}/settle.lmp",
    params:
        n_org_types=lambda wc, input: n_org_types(SOLVENT_MERGED_DATA),
        cutoff=config["slab"]["merge_cutoff_angstrom"],
        kspace_accuracy=config["slab"]["merge_kspace_accuracy"],
        timestep_fs=config["slab"]["merge_nvt_timestep_fs"],
        temperature_K=config["temperature_K"],
        read_data="nvt_shrunk_minimized.data",
        write_prefix="npt_equilibrated",
        has_fragment=False,
        nsteps=lambda wc: int(config["slab"]["settle_ps"] * 1000 / config["slab"]["merge_nvt_timestep_fs"]),
    localrule: True
    script:
        "../scripts/render_template.py"


rule lmp_settle:
    input:
        lmp=f"{SOLVENT_RES}/settle.lmp",
        data=f"{SOLVENT_RES}/nvt_shrunk_minimized.data",
        coeffs=f"{SOLVENT_RES}/merged_coeffs.lmp",
    output:
        data=f"{SOLVENT_RES}/npt_equilibrated.data",
        restart=f"{SOLVENT_RES}/npt_equilibrated.restart",
    params:
        omp_threads=lambda wc, resources: resources.get("ompthreads", 1),
    shell:
        'cd "$(dirname {input.lmp})" && '
        'LMP_BIN="$(command -v lmp_mpi || command -v lmp)" && '
        '"$LMP_BIN" -sf omp -pk omp {params.omp_threads} -in "$(basename {input.lmp})"'


# ===========================================================================
# Part 2: species-specific fragment insertion + short relax + production
# ===========================================================================

rule insert_fragment:
    input:
        solvent_data=f"{SOLVENT_RES}/npt_equilibrated.data",
        solvent_coeffs=f"{SOLVENT_RES}/merged_coeffs.lmp",
        fragment_data=f"{SLAB_RES}/{{species}}/fragment_converted.lmp",
        fragment_input=f"{SLAB_RES}/{{species}}/fragment_converted.input",
        script="scripts/insert_fragment.py",
    output:
        data=f"{SLAB_RES}/{{species}}/inserted.data",
        # Literally named merged_coeffs.lmp, not inserted_coeffs.lmp --
        # lmp/_common_setup.lmp.template hardcodes `include
        # merged_coeffs.lmp` (not parameterized), and every downstream
        # LAMMPS run below `cd`s into this species' own directory, so
        # reusing that exact filename here (now meaning "hexadecane +
        # fragment + substrate coeffs", not "hexadecane + substrate") is
        # simpler than parameterizing the include.
        coeffs=f"{SLAB_RES}/{{species}}/merged_coeffs.lmp",
    params:
        n_sub_types=N_SUB_TYPES,
        clash_angstrom=config["slab"]["insert_clash_angstrom"],
    wildcard_constraints:
        species="ST|MMA",
    localrule: True
    shell:
        "python {input.script} {input.solvent_data} "
        "--solvent-coeffs {input.solvent_coeffs} "
        "--fragment-data {input.fragment_data} --fragment-input {input.fragment_input} "
        "--n-sub-types {params.n_sub_types} --clash-angstrom {params.clash_angstrom} "
        "-o {output.data} --out-coeffs {output.coeffs}"


rule render_insert_minimize_lmp:
    input:
        "lmp/minimize.lmp.template",
        LMP_COMMON_TEMPLATE,
        inserted_data_path("{species}"),
    output:
        f"{SLAB_RES}/{{species}}/insert_minimize.lmp",
    params:
        n_org_types=lambda wc, input: n_org_types(inserted_data_path(wc.species)),
        max_org_molid=lambda wc, input: insert_molid_info(inserted_data_path(wc.species))[0],
        frag_molid_lo=lambda wc, input: insert_molid_info(inserted_data_path(wc.species))[1],
        frag_molid_hi=lambda wc, input: insert_molid_info(inserted_data_path(wc.species))[2],
        cutoff=config["slab"]["merge_cutoff_angstrom"],
        kspace_accuracy=config["slab"]["merge_kspace_accuracy"],
        timestep_fs=config["slab"]["merge_nvt_timestep_fs"],
        temperature_K=config["temperature_K"],
        read_data="inserted.data",
        write_prefix="inserted_minimized",
        velocity_create=True,
        velocity_create_group="fragment",
        seed=config["slab"]["merge_nvt_seed"],
        etol=config["slab"]["merge_minimize"]["etol"],
        ftol=config["slab"]["merge_minimize"]["ftol"],
        maxiter=config["slab"]["merge_minimize"]["maxiter"],
        maxeval=config["slab"]["merge_minimize"]["maxeval"],
        soft_cutoff=config["slab"]["soft_premin"]["soft_cutoff"],
        soft_a=config["slab"]["soft_premin"]["soft_a"],
        soft_etol=config["slab"]["soft_premin"]["soft_etol"],
        soft_ftol=config["slab"]["soft_premin"]["soft_ftol"],
        soft_maxiter=config["slab"]["soft_premin"]["soft_maxiter"],
        soft_maxeval=config["slab"]["soft_premin"]["soft_maxeval"],
    wildcard_constraints:
        species="ST|MMA",
    localrule: True
    script:
        "../scripts/render_template.py"


rule lmp_insert_minimize:
    input:
        lmp=f"{SLAB_RES}/{{species}}/insert_minimize.lmp",
        data=inserted_data_path("{species}"),
        coeffs=f"{SLAB_RES}/{{species}}/merged_coeffs.lmp",
    output:
        data=f"{SLAB_RES}/{{species}}/inserted_minimized.data",
        restart=f"{SLAB_RES}/{{species}}/inserted_minimized.restart",
    params:
        omp_threads=lambda wc, resources: resources.get("ompthreads", 1),
    wildcard_constraints:
        species="ST|MMA",
    shell:
        'cd "$(dirname {input.lmp})" && '
        'LMP_BIN="$(command -v lmp_mpi || command -v lmp)" && '
        '"$LMP_BIN" -sf omp -pk omp {params.omp_threads} -in "$(basename {input.lmp})"'


# --- short, box-fixed, restrained relaxation of the local insertion
# disturbance -- see this file's module docstring for why there is
# deliberately NO barostat stage here --------------------------------

rule render_insert_relax_lmp:
    input:
        "lmp/nvt.lmp.template",
        LMP_COMMON_TEMPLATE,
        inserted_data_path("{species}"),
    output:
        f"{SLAB_RES}/{{species}}/insert_relax.lmp",
    params:
        n_org_types=lambda wc, input: n_org_types(inserted_data_path(wc.species)),
        max_org_molid=lambda wc, input: insert_molid_info(inserted_data_path(wc.species))[0],
        frag_molid_lo=lambda wc, input: insert_molid_info(inserted_data_path(wc.species))[1],
        frag_molid_hi=lambda wc, input: insert_molid_info(inserted_data_path(wc.species))[2],
        cutoff=config["slab"]["merge_cutoff_angstrom"],
        kspace_accuracy=config["slab"]["merge_kspace_accuracy"],
        timestep_fs=config["slab"]["merge_nvt_timestep_fs"],
        temperature_K=config["temperature_K"],
        read_data="inserted_minimized.data",
        write_prefix="relaxed",
        has_fragment=True,
        restraint_k=config["slab"]["fragment_restraint_k_kcal_per_mol_A2"],
        nsteps=lambda wc: int(config["slab"]["insert_relax_ps"] * 1000 / config["slab"]["merge_nvt_timestep_fs"]),
    wildcard_constraints:
        species="ST|MMA",
    localrule: True
    script:
        "../scripts/render_template.py"


rule lmp_insert_relax:
    input:
        lmp=f"{SLAB_RES}/{{species}}/insert_relax.lmp",
        data=f"{SLAB_RES}/{{species}}/inserted_minimized.data",
        coeffs=f"{SLAB_RES}/{{species}}/merged_coeffs.lmp",
    output:
        data=f"{SLAB_RES}/{{species}}/relaxed.data",
        restart=f"{SLAB_RES}/{{species}}/relaxed.restart",
    params:
        omp_threads=lambda wc, resources: resources.get("ompthreads", 1),
    wildcard_constraints:
        species="ST|MMA",
    shell:
        'cd "$(dirname {input.lmp})" && '
        'LMP_BIN="$(command -v lmp_mpi || command -v lmp)" && '
        '"$LMP_BIN" -sf omp -pk omp {params.omp_threads} -in "$(basename {input.lmp})"'


# --- stage 5: production, with z-resolved density profiles ---------------

rule render_production_lmp:
    input:
        "lmp/production.lmp.template",
        LMP_COMMON_TEMPLATE,
        inserted_data_path("{species}"),
    output:
        f"{SLAB_RES}/{{species}}/production.lmp",
    params:
        n_org_types=lambda wc, input: n_org_types(inserted_data_path(wc.species)),
        max_org_molid=lambda wc, input: insert_molid_info(inserted_data_path(wc.species))[0],
        frag_molid_lo=lambda wc, input: insert_molid_info(inserted_data_path(wc.species))[1],
        frag_molid_hi=lambda wc, input: insert_molid_info(inserted_data_path(wc.species))[2],
        cutoff=config["slab"]["merge_cutoff_angstrom"],
        kspace_accuracy=config["slab"]["merge_kspace_accuracy"],
        timestep_fs=config["slab"]["merge_nvt_timestep_fs"],
        temperature_K=config["temperature_K"],
        read_data="relaxed.data",
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
        data=f"{SLAB_RES}/{{species}}/relaxed.data",
        coeffs=f"{SLAB_RES}/{{species}}/merged_coeffs.lmp",
    output:
        data=f"{SLAB_RES}/{{species}}/production.data",
        restart=f"{SLAB_RES}/{{species}}/production.restart",
        fragment_density=f"{SLAB_RES}/{{species}}/fragment_density.profile",
        hexadecane_density=f"{SLAB_RES}/{{species}}/hexadecane_density.profile",
    params:
        omp_threads=lambda wc, resources: resources.get("ompthreads", 1),
    wildcard_constraints:
        species="ST|MMA",
    shell:
        'cd "$(dirname {input.lmp})" && '
        'LMP_BIN="$(command -v lmp_mpi || command -v lmp)" && '
        '"$LMP_BIN" -sf omp -pk omp {params.omp_threads} -in "$(basename {input.lmp})"'
