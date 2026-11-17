"""Step 0 (segment size) + Step 1 (chi_AS, chi_BS): bulk hexadecane and
capped ST/MMA oligomer "melt" boxes, GROMACS em -> nvt -> squeeze -> npt,
then read off density (step0) and cohesive energy density / solubility
parameter (step1).

STATUS: runs end-to-end through both segment_volume (step0) and chi.yaml
(step1) targets. The nvt->squeeze->npt chain and run lengths were
validated via standalone gmx_mpi benchmarking (not yet re-confirmed
through this exact rule chain) on 2026-11-05: HEX and ST10 both
reach/overshoot their real liquid density within the squeeze phase and
relax to a stable plateau within ~30 ps of switching back to ref_p=1 bar
(HEX: ~748 kg/m^3 vs. a ~770 kg/m^3 literature reference -- a normal
few-percent OPLS/L-OPLS deviation, not a red flag).

solubility_chi's cohesive-energy-density numbers, current as of the
2026-11-16 fix below: delta_HEX=15.10, delta_ST=11.22, delta_MMA=13.04
MPa^0.5, chi_AS=1.83, chi_BS=0.51.

**2026-11-05 run was methodologically broken, not just noisy:**
delta_HEX=14.06, delta_ST=11.44, delta_MMA=14.12, chi_AS=0.83, chi_BS~0.
chi_BS coming out ~0 (MMA's melt looking exactly as cohesive as plain
hexadecane) turned out to be a real artifact, confirmed live 2026-11-16
via energy decomposition (gmx energy LJ-14/Coulomb-14/LJ-(SR)/
Coulomb-(SR)/Coul.-recip.): an UNRESTRAINED floppy 10-mer collapses onto
itself over vacuum_nvt's 100 ps, self-solvating its own polar groups
intramolecularly almost as effectively as the real liquid phase does
intermolecularly -- MMA10's vacuum Coulomb(SR) (-1304 kJ/mol) nearly
equaled the liquid's per-molecule Coulomb(SR) (-1319 kJ/mol), collapsing
the vacuum/liquid energy gap the whole CED method depends on. ST10 showed
the same mechanism, just diluted by its large nonpolar C18 side chain
(LJ(SR) only ~50% cancelled, vs. MMA's near-total Coulomb(SR)
cancellation). Fixed by adding a one-sided pull restraint (see
vacuum_nvt's docstring and mdp/vacuum_nvt.mdp) that stops the chain
folding below roughly an extended conformation without rigidly freezing
it. Two bugs surfaced getting this restraint to actually work, both
confirmed live: grompp's `-n` suppresses its own auto-generated "System"
index group (had to be written explicitly into the custom .ndx), and
GROMACS's `flat-bottom` pull type restrains ABOVE the target distance,
not below it -- `flat-bottom-high` is the one-sided-below type actually
wanted here.

chi_AS still coming out noticeably larger than chi_BS (the polar-ester
MMA binding hexadecane better than the alkyl-rich ST does) remains the
opposite of naive "like dissolves like" intuition, and the fix above
didn't change that ordering (if anything, sharpened it). Not yet
resolved -- directionally plausible (ST10's comb architecture, a C18 side
chain hanging off every backbone unit, packs less efficiently than linear
hexadecane than naive "alkyl similarity" reasoning would suggest) but
could also still be under-sampling (bulk_npt's 200 ps production window
may be too short for ST10's floppier chains to reach a conformationally
representative ensemble even though density itself plateaus much
faster). This is exactly the known limitation forcefield/chi_chis_from_md.md's
step 1 already calls out ("primary, cheap" solubility-parameter route,
with alchemical TI/BAR and single-chain R_g matching as listed
refinements) -- not a reason to block on this before moving on, but a
reason to sanity-check these specific numbers again at step 3's
full-chain MD comparison before trusting them in isolation.
"""

RES = "results/bulk"


rule build_itp_oligomer:
    """Capped ST/MMA oligomer topology via polyply gen_params."""
    input:
        ff=[
            f"{FF_DIR}/monomers/PAMA_{{species}}.oplsaa.LigParGen.ff",
            f"{FF_DIR}/monomers/CH3.ff",
            f"{FF_DIR}/monomers/CH3n.ff",
            f"{FF_DIR}/PAMA_ST_MMA.oplsaa.LigParGen_links.ff",
        ],
    output:
        itp=f"{RES}/{{species}}/{{species}}{N_OLIGO}.itp",
    wildcard_constraints:
        species="ST|MMA",
    conda:
        "../envs/polyply.yaml"
    localrule: True
    params:
        # polyply's -seq takes "resname:n_blocks" tokens (current PyPI
        # release, 1.8.0) -- not a flat space-repeated list of names.
        # CH3n on BOTH ends, not CH3 -- confirmed live 2026-11-12: CH3.ff's
        # plain literature OPLS charges (C -0.3029, H 0.0694 x3) sum to
        # -0.0947, while CH3n.ff's are the same values uniformly shifted by
        # +0.023675/atom to sum to exactly 0 (PAMA_ST_MMA.oplsaa.LigParGen_
        # links.ff bonds CH3 and CH3n identically to either backbone
        # attachment carbon, so they're structurally interchangeable -- CH3
        # itself is never used again, kept only as the other half of that
        # now-unused correction pair). This is why every oligomer/fragment
        # built here carried a net charge of exactly -0.0947 e and every
        # LAMMPS run warned "System is not charge neutral" regardless of
        # species.
        seq=lambda wc: f"CH3n:1 {wc.species}:{N_OLIGO} CH3n:1",
        name=lambda wc: f"{wc.species}{N_OLIGO}",
    shell:
        "polyply gen_params -name {params.name} -f {input.ff} "
        "-o {output.itp} -seq {params.seq}"


rule bulk_itp_hex:
    """Hexadecane is already a complete itp (forcefield/solvent/hexadecane.itp)."""
    input:
        f"{FF_DIR}/solvent/hexadecane.itp",
    output:
        f"{RES}/HEX/HEX.itp",
    localrule: True
    shell:
        "cp {input} {output}"


rule bulk_top:
    """Written AFTER bulk_pack (not before -- see that rule) so the
    [ molecules ] count always matches what the packer actually placed, not
    what config.yaml asked for. gmx insert-molecules' success count was
    non-deterministic run to run (confirmed live 2026-11-04: the exact same
    -nmol/-box/-try settings placed all 40 ST10 copies once and only 39 the
    next time) -- deriving n_chains from box0.gro's real atom count instead
    of trusting the request avoided a grompp "number of coordinates ...
    does not match topology" failure every time the random placement came
    up even one molecule short. bulk_pack now uses packmol instead, which is
    far more reliable but still not guaranteed deterministic (see its
    docstring), so this safety net stays."""
    input:
        itp=lambda wc: f"{RES}/{wc.species}/{molname(wc.species)}.itp",
        ff_itp=f"{FF_DIR}/oplsaa.ff/forcefield.itp",
        nb=f"{FF_DIR}/custom-ffnonbonded.itp",
        box=f"{RES}/{{species}}/box0.gro",
        single=f"{RES}/{{species}}/single_em.gro",
    output:
        top=f"{RES}/{{species}}/topol.top",
    localrule: True
    params:
        n_chains_requested=lambda wc: config["bulk"]["n_chains"][wc.species],
        mol=lambda wc: molname(wc.species),
    run:
        import os
        import sys

        def count_atoms(gro_path):
            with open(gro_path) as fh:
                fh.readline()
                return int(fh.readline())

        n_atoms_per_mol = count_atoms(input.single)
        n_chains = count_atoms(input.box) // n_atoms_per_mol
        if n_chains != params.n_chains_requested:
            print(
                f"bulk_top: gmx insert-molecules placed {n_chains} of "
                f"{params.n_chains_requested} requested {params.mol} "
                "molecules -- using the actual count.",
                file=sys.stderr,
            )

        ff_itp_rel = os.path.relpath(input.ff_itp, os.path.dirname(output.top))
        nb_rel = os.path.relpath(input.nb, os.path.dirname(output.top))
        itp_rel = os.path.basename(input.itp)
        with open(output.top, "w") as fh:
            fh.write(f'#include "{ff_itp_rel}"\n')
            fh.write(f'#include "{nb_rel}"\n')
            fh.write(f'#include "{itp_rel}"\n\n')
            fh.write("[ system ]\n")
            fh.write(f"{params.mol} bulk box\n\n")
            fh.write("[ molecules ]\n")
            fh.write(f"{params.mol}    {n_chains}\n")


rule single_top:
    """n=1 topology, for resolving one molecule's real geometry before
    replicating it (see single_mol_em / bulk_pack below)."""
    input:
        itp=lambda wc: f"{RES}/{wc.species}/{molname(wc.species)}.itp",
        ff_itp=f"{FF_DIR}/oplsaa.ff/forcefield.itp",
        nb=f"{FF_DIR}/custom-ffnonbonded.itp",
    output:
        top=f"{RES}/{{species}}/single.top",
    localrule: True
    params:
        mol=lambda wc: molname(wc.species),
    run:
        import os

        ff_itp_rel = os.path.relpath(input.ff_itp, os.path.dirname(output.top))
        nb_rel = os.path.relpath(input.nb, os.path.dirname(output.top))
        itp_rel = os.path.basename(input.itp)
        with open(output.top, "w") as fh:
            fh.write(f'#include "{ff_itp_rel}"\n')
            fh.write(f'#include "{nb_rel}"\n')
            fh.write(f'#include "{itp_rel}"\n\n')
            fh.write("[ system ]\n")
            fh.write(f"single {params.mol}\n\n")
            fh.write("[ molecules ]\n")
            fh.write(f"{params.mol}    1\n")


rule single_mol_raw:
    """One isolated ST10/MMA10 molecule via polyply gen_coords, in a big
    empty box. See single_mol_em below for why this needs fixing up before
    it's usable."""
    input:
        top=f"{RES}/{{species}}/single.top",
        itp=lambda wc: f"{RES}/{wc.species}/{molname(wc.species)}.itp",
    output:
        gro=f"{RES}/{{species}}/single0.gro",
    wildcard_constraints:
        species="ST|MMA",
    conda:
        "../envs/polyply.yaml"
    localrule: True
    shell:
        "polyply gen_coords -p {input.top} -o {output.gro} -box 5.0 5.0 5.0"


rule single_mol_raw_hex:
    """HEX already has a correctly-built single-molecule structure
    (forcefield/solvent/hexadecane.itp's companion .gro, from LigParGen) --
    center it in a box instead of asking polyply to invent one."""
    input:
        gro=f"{FF_DIR}/solvent/hexadecane.gro",
    output:
        gro=f"{RES}/HEX/single0.gro",
    shell:
        "gmx_mpi editconf -f {input.gro} -o {output.gro} -c -box 5.0 5.0 5.0"


rule single_mol_em:
    """Resolve the single molecule's real geometry before replicating it.
    confirmed live 2026-11-04: polyply's gen_coords (both the PyPI release
    and the jotelha/polyply_1.0 "2026-10-28-PAMA-parameters" fork) builds de
    novo 3D structure for these all-atom OPLS oligomers/hexadecane with
    severely wrong local geometry -- bond-length-scale distances off by
    2-3x, hundreds of sub-0.05 nm atom clashes per molecule -- regardless of
    whether it's one isolated molecule or many packed together. Plain
    unconstrained steepest descent with a small step (mdp/em.mdp) resolves
    a single molecule cleanly (verified: ST10 -> 4.4e3 kJ/mol, MMA10 ->
    3.2e3 kJ/mol, both physically sane). This EM'd single molecule, not the
    raw single0.gro, is what bulk_pack replicates via gmx insert-molecules."""
    input:
        gro=f"{RES}/{{species}}/single0.gro",
        top=f"{RES}/{{species}}/single.top",
        mdp="mdp/em.mdp",
    output:
        gro=f"{RES}/{{species}}/single_em.gro",
        tpr=f"{RES}/{{species}}/single_em.tpr",
    params:
        prefix=f"{RES}/{{species}}/single_em",
        mdp_rendered=f"{RES}/{{species}}/single_em.rendered.mdp",
    threads: 8
    run:
        render_mdp(input.mdp, params.mdp_rendered, em_nsteps=config["bulk"]["em_nsteps"])
        shell(
            "gmx_mpi grompp -f {params.mdp_rendered} -c {input.gro} -p {input.top} "
            "-o {output.tpr} -maxwarn 2"
        )
        shell("gmx_mpi mdrun -deffnm {params.prefix} -ntomp {threads}")


rule bulk_pack:
    """Replicate the EM-resolved single molecule n_chains times via packmol
    (optimization-based packing) -- NOT polyply gen_coords on the full
    system directly (confirmed live 2026-11-04: compounds its own
    per-molecule distortion across hundreds of copies into catastrophic
    overlaps, potential energies up to ~1e19 kJ/mol), and NOT gmx
    insert-molecules either (confirmed live 2026-11-04: its random-reject
    placement is non-deterministic -- the exact same -nmol/-box/-try
    settings placed all 40 ST10 copies in one run and only 39 in the next,
    silently writing a box0.gro with fewer atoms than topol.top's molecule
    count declares, which grompp then rejects). packmol's optimization-based
    approach reliably hits the exact requested count every time, verified
    directly for both the easy case (HEX, <6 s) and the hard case (ST10's
    long, floppy C18 tail, ~2.3 min at this box size with the 2026-10-29-era
    monomer parameters).

    Needs a PDB round-trip: packmol doesn't read/write .gro, so the
    EM-resolved single molecule is converted to PDB first, and packmol's
    PDB output is converted back to .gro with the real box vectors
    reattached (a bare PDB has no periodic box info).

    Tolerates packmol's own STOP 173 ("ended without perfect packing"):
    confirmed live 2026-11-16, with the updated (2026-11-12-PAMA-parameters)
    ST monomer, GENCAN failed to reach a zero-overlap solution twice in a
    row at nloop=1000 (~90 min each, residual function value ~215 across
    26720 atoms) despite the single molecule's own geometry (extent/Rg)
    being unchanged in any obvious way from the pre-update case that always
    packed cleanly -- GENCAN's starting point is random and no seed is
    fixed, so this reads as this specific packing problem being
    close-to-but-not-quite within nloop=1000's margin, not a geometry bug.
    packmol's own STOP-173 message explicitly recommends using the result
    anyway ("still very likely that the current point fits your needs...
    minimize the energy of the solution found"), and bulk_em right after
    this rule is already designed (see mdp/em.mdp's own docstring) to
    resolve far worse residual overlaps than this (gmx insert-molecules'
    clashes, bond-length-scale distances 3x too short) with room to spare
    -- so retrying packmol for another ~90 min bought nothing a free EM
    step doesn't already give for free. On STOP 173, use packmol's own
    box-constraint-forced output (`packed.pdb_FORCED`) instead of the
    nominal `packed.pdb` (which may not even respect the box's periodic
    boundary in this case) and let bulk_em do the rest. Any other packmol
    exit code still fails the rule -- this carve-out is specific to the
    documented 173 case, not a blanket "ignore packmol errors"."""
    input:
        gro=f"{RES}/{{species}}/single_em.gro",
    output:
        gro=f"{RES}/{{species}}/box0.gro",
    params:
        n_chains=lambda wc: config["bulk"]["n_chains"][wc.species],
        box_nm=config["bulk"]["box_nm"],
        pdb_in=f"{RES}/{{species}}/single_em.pdb",
        packed_pdb=f"{RES}/{{species}}/packed.pdb",
        inp_file=f"{RES}/{{species}}/packmol.inp",
    run:
        import subprocess
        import sys

        box_ang = [b * 10 for b in params.box_nm]  # nm -> Angstrom (packmol/PDB units)

        shell("gmx_mpi editconf -f {input.gro} -o {params.pdb_in}")

        with open(params.inp_file, "w") as fh:
            fh.write("tolerance 2.0\n")
            fh.write("filetype pdb\n")
            # Default nloop (200 for a single structure type) isn't always
            # enough -- confirmed live 2026-11-04: the exact same ST10
            # inputs/box converged in one run and hit "Maximum number of
            # GENCAN loops achieved" (exit 173) in the next, since packmol's
            # GENCAN optimizer starts from a random configuration each time
            # (no seed fixed below). 1000 gives real margin without being
            # the reason to add complexity like a fixed seed + retry loop.
            fh.write("nloop 1000\n")
            fh.write(f"output {params.packed_pdb}\n\n")
            fh.write(f"structure {params.pdb_in}\n")
            fh.write(f"  number {params.n_chains}\n")
            fh.write(
                f"  inside box 0. 0. 0. {box_ang[0]} {box_ang[1]} {box_ang[2]}\n"
            )
            fh.write("end structure\n")

        with open(params.inp_file) as fh:
            result = subprocess.run(["packmol"], stdin=fh)

        pdb_for_editconf = params.packed_pdb
        if result.returncode == 173:
            # packmol names its box-constraint-forced fallback by appending
            # "_FORCED" to the whole output filename (packed.pdb_FORCED, not
            # packed_FORCED.pdb) -- gmx editconf identifies format purely by
            # the trailing extension, so that name is rejected outright
            # ("does not have a recognizable extension"); give it a
            # same-content, properly-suffixed copy instead of renaming the
            # original (keep the raw packmol output around for inspection).
            import shutil

            forced_pdb_raw = params.packed_pdb + "_FORCED"
            forced_pdb = str(params.packed_pdb).replace(".pdb", "_forced.pdb")
            shutil.copyfile(forced_pdb_raw, forced_pdb)
            print(
                f"bulk_pack: packmol ended without perfect packing (STOP 173) "
                f"for {wildcards.species} -- using its box-constraint-forced "
                f"best-effort solution {forced_pdb_raw} (copied to {forced_pdb} "
                f"for gmx) instead of {params.packed_pdb}; bulk_em resolves the "
                "residual overlaps from here.",
                file=sys.stderr,
            )
            pdb_for_editconf = forced_pdb
        elif result.returncode != 0:
            raise RuntimeError(
                f"packmol failed with exit code {result.returncode} "
                f"(not the tolerated STOP-173 'imperfect packing' case) -- see "
                f"{params.inp_file} and the console output above."
            )

        shell(
            "gmx_mpi editconf -f {pdb_for_editconf} -o {output.gro} "
            "-box {params.box_nm[0]} {params.box_nm[1]} {params.box_nm[2]}"
        )


rule bulk_em:
    input:
        gro=f"{RES}/{{species}}/box0.gro",
        top=f"{RES}/{{species}}/topol.top",
        mdp="mdp/em.mdp",
    output:
        gro=f"{RES}/{{species}}/em.gro",
        tpr=f"{RES}/{{species}}/em.tpr",
    params:
        prefix=f"{RES}/{{species}}/em",
        mdp_rendered=f"{RES}/{{species}}/em.rendered.mdp",
    threads: 8
    run:
        render_mdp(input.mdp, params.mdp_rendered, em_nsteps=config["bulk"]["em_nsteps"])
        shell(
            "gmx_mpi grompp -f {params.mdp_rendered} -c {input.gro} -p {input.top} "
            "-o {output.tpr} -maxwarn 2"
        )
        shell("gmx_mpi mdrun -deffnm {params.prefix} -ntomp {threads}")


rule bulk_nvt:
    input:
        gro=f"{RES}/{{species}}/em.gro",
        top=f"{RES}/{{species}}/topol.top",
        mdp="mdp/nvt.mdp",
    output:
        gro=f"{RES}/{{species}}/nvt.gro",
        tpr=f"{RES}/{{species}}/nvt.tpr",
    params:
        prefix=f"{RES}/{{species}}/nvt",
        mdp_rendered=f"{RES}/{{species}}/nvt.rendered.mdp",
        nsteps=lambda wc: int(config["bulk"]["nvt_ps"] / 0.002),
    threads: 8
    run:
        render_mdp(
            input.mdp,
            params.mdp_rendered,
            nvt_nsteps=params.nsteps,
            temperature_K=config["temperature_K"],
        )
        shell(
            "gmx_mpi grompp -f {params.mdp_rendered} -c {input.gro} -r {input.gro} "
            "-p {input.top} -o {output.tpr} -maxwarn 2"
        )
        shell("gmx_mpi mdrun -deffnm {params.prefix} -ntomp {threads}")


rule bulk_squeeze:
    """High-pressure pre-compression, between NVT and the real NPT -- see
    config.yaml's bulk.squeeze_pressure_bar/squeeze_ps comment for why this
    step exists (going straight from the loosely-packed box0.gro to a
    normal ref_p=1 bar NPT doesn't work: confirmed live 2026-11-05, density
    only crept up linearly and nowhere near the real target over 500 ps)."""
    input:
        gro=f"{RES}/{{species}}/nvt.gro",
        top=f"{RES}/{{species}}/topol.top",
        mdp="mdp/squeeze.mdp",
    output:
        gro=f"{RES}/{{species}}/squeeze.gro",
        tpr=f"{RES}/{{species}}/squeeze.tpr",
    params:
        prefix=f"{RES}/{{species}}/squeeze",
        mdp_rendered=f"{RES}/{{species}}/squeeze.rendered.mdp",
        nsteps=lambda wc: int(config["bulk"]["squeeze_ps"] / 0.002),
    threads: 8
    run:
        render_mdp(
            input.mdp,
            params.mdp_rendered,
            squeeze_nsteps=params.nsteps,
            temperature_K=config["temperature_K"],
            squeeze_pressure_bar=config["bulk"]["squeeze_pressure_bar"],
        )
        shell(
            "gmx_mpi grompp -f {params.mdp_rendered} -c {input.gro} -r {input.gro} "
            "-p {input.top} -o {output.tpr} -maxwarn 2"
        )
        shell("gmx_mpi mdrun -deffnm {params.prefix} -ntomp {threads}")


rule bulk_npt:
    input:
        gro=f"{RES}/{{species}}/squeeze.gro",
        top=f"{RES}/{{species}}/topol.top",
        mdp="mdp/npt.mdp",
    output:
        gro=f"{RES}/{{species}}/npt.gro",
        edr=f"{RES}/{{species}}/npt.edr",
        tpr=f"{RES}/{{species}}/npt.tpr",
    params:
        prefix=f"{RES}/{{species}}/npt",
        mdp_rendered=f"{RES}/{{species}}/npt.rendered.mdp",
        nsteps=lambda wc: int(config["bulk"]["npt_ps"] / 0.002),
    threads: 8
    run:
        render_mdp(
            input.mdp,
            params.mdp_rendered,
            npt_nsteps=params.nsteps,
            temperature_K=config["temperature_K"],
            pressure_bar=config["pressure_bar"],
        )
        shell(
            "gmx_mpi grompp -f {params.mdp_rendered} -c {input.gro} -r {input.gro} "
            "-p {input.top} -o {output.tpr} -maxwarn 2"
        )
        shell("gmx_mpi mdrun -deffnm {params.prefix} -ntomp {threads}")


rule bulk_density:
    """Step 0 input: average mass density -> a^3 (segment volume mapping).
    -b skips bulk_npt's own equilibration window (relaxing the squeeze
    phase's overshoot back down to the real ref_p) -- see config.yaml's
    bulk.npt_discard_ps comment."""
    input:
        edr=f"{RES}/{{species}}/npt.edr",
    output:
        xvg=f"{RES}/{{species}}/density.xvg",
    params:
        discard_ps=config["bulk"]["npt_discard_ps"],
    shell:
        'echo "Density" | gmx_mpi energy -f {input.edr} -b {params.discard_ps} -o {output.xvg}'


rule bulk_potential:
    """Step 1 input: average liquid-phase potential energy (same production
    window as bulk_density) -- U_total/n_chains feeds the cohesive-energy-
    density split in scripts/solubility_chi.py."""
    input:
        edr=f"{RES}/{{species}}/npt.edr",
    output:
        xvg=f"{RES}/{{species}}/liquid_potential.xvg",
    params:
        discard_ps=config["bulk"]["npt_discard_ps"],
    shell:
        'echo "Potential" | gmx_mpi energy -f {input.edr} -b {params.discard_ps} -o {output.xvg}'


rule segment_volume:
    """Step 0: a^3 from HEX density, and ST/MMA-per-lattice-segment ratio
    from the oligomer melt densities + known molar masses."""
    input:
        density=expand(f"{RES}/{{species}}/density.xvg", species=BULK_SPECIES),
    output:
        "results/segment_volume.yaml",
    params:
        species=BULK_SPECIES,
        n_oligo=N_OLIGO,
    conda:
        "../envs/analysis.yaml"
    localrule: True
    script:
        "../scripts/segment_volume.py"


rule vacuum_nvt:
    """Thermal average of the single molecule's intramolecular-only
    potential energy at temperature_K -- the vacuum-reference half of the
    cohesive-energy-density split (bulk_potential is the liquid half).
    Reuses mdp/nvt.mdp's own nonbonded settings VERBATIM (same PME, same
    1.1 nm cutoffs as the liquid-phase runs, via mdp/vacuum_nvt.mdp) on
    single_em.gro (already cleanly resolved -- see mdp/em.mdp), not a
    separate "plain cutoff" vacuum treatment -- confirmed live 2026-11-05
    that mismatch matters: an earlier attempt gave ST10/MMA10 *lower*
    solubility parameters than plain hexadecane, which is physically
    backwards for ester-containing oligomers, because the different
    nonbonded method biased the subtraction. Single-molecule periodic
    images at this 5 nm box's separation are negligible, so
    same-method-as-the-liquid is what actually matters, not strict
    non-periodicity.

    Confirmed live 2026-11-16: an UNRESTRAINED single oligomer in vacuum
    folds onto itself over this 100 ps run, nearly reproducing the liquid
    phase's own intermolecular cohesive energy *intramolecularly*
    (self-solvating its ester groups) -- MMA10's vacuum Coulomb(SR)
    (-1304 kJ/mol) came out almost equal to the liquid's per-molecule
    Coulomb(SR) (-1319 kJ/mol), collapsing the vacuum/liquid gap this whole
    method depends on and driving chi_BS to ~0. ST10 showed the same
    mechanism (less severely, since its C18 side tail dilutes the effect).
    mdp/vacuum_nvt.mdp adds a one-sided (flat-bottom) pull restraint
    between the oligomer's two terminal cap atoms, sized from its own
    backbone bond-path length (backbone_bond_path_length in common.smk) so
    it can't collapse below roughly an extended conformation, but is free
    to fluctuate above that -- not a rigid freeze, just removing the
    folding pathway that was corrupting the vacuum reference state.

    Confirmed live 2026-11-16: ST10's single_em.gro is ALREADY close to
    fully extended (2.78 nm end-to-end, vs. a ~2.75 nm backbone-bond-path
    estimate) -- grompp refuses to even set up the pull code in the
    shared 5 nm single-molecule box (single_mol_raw's box_nm, also used by
    bulk_pack) because 2.78 nm already exceeds its own minimum-image
    safety limit (0.49 * box size = 2.45 nm), regardless of the restraint
    target. Re-boxes into a bigger, vacuum_nvt-only box (not touching
    single_em.gro or anything downstream of it) rather than growing the
    shared box, which would needlessly cascade into bulk_pack/bulk_top for
    all three species."""
    input:
        gro=f"{RES}/{{species}}/single_em.gro",
        top=f"{RES}/{{species}}/single.top",
        itp=lambda wc: f"{RES}/{wc.species}/{molname(wc.species)}.itp",
        mdp="mdp/vacuum_nvt.mdp",
    output:
        box_gro=f"{RES}/{{species}}/vacuum_box.gro",
        gro=f"{RES}/{{species}}/vacuum_nvt.gro",
        edr=f"{RES}/{{species}}/vacuum_nvt.edr",
        tpr=f"{RES}/{{species}}/vacuum_nvt.tpr",
        ndx=f"{RES}/{{species}}/vacuum_pull.ndx",
    params:
        prefix=f"{RES}/{{species}}/vacuum_nvt",
        mdp_rendered=f"{RES}/{{species}}/vacuum_nvt.rendered.mdp",
        nsteps=lambda wc: int(config["bulk"]["vacuum_nvt_ps"] / 0.002),
        extend_fraction=config["bulk"]["vacuum_extend_fraction"],
        pull_k=config["bulk"]["vacuum_pull_k_kj_per_mol_nm2"],
        projected_bond_length_nm=0.125,  # tetrahedral sp3 zigzag projection
    threads: 4
    run:
        def count_atoms(gro_path):
            with open(gro_path) as fh:
                fh.readline()
                return int(fh.readline())

        n_atoms = count_atoms(input.gro)
        n_bonds = backbone_bond_path_length(input.itp, n_atoms)
        full_extension_nm = n_bonds * params.projected_bond_length_nm
        pull_target_nm = params.extend_fraction * full_extension_nm

        # Box must satisfy grompp's own minimum-image safety check
        # (pull-group distance < 0.49 * box size) for BOTH the restraint
        # target and whatever the starting structure's own end-to-end
        # distance already is -- 1.2x margin on top of the tighter of
        # the two 0.49-scaled requirements, floor at the original 5 nm.
        box_edge_nm = max(5.0, full_extension_nm / 0.49 * 1.2)
        shell(
            "gmx_mpi editconf -f {input.gro} -o {output.box_gro} "
            "-box {box_edge_nm} {box_edge_nm} {box_edge_nm} -c"
        )

        with open(output.ndx, "w") as fh:
            # grompp -n suppresses its own auto-generated default groups
            # (confirmed live 2026-11-16: "System" came back missing for
            # tc-grps once a custom index file was supplied) -- so System
            # must be written here explicitly, not just the two pull groups.
            system_atoms = " ".join(str(i) for i in range(1, n_atoms + 1))
            fh.write(f"[ System ]\n{system_atoms}\n\n")
            fh.write(f"[ ChainEnd1 ]\n1\n\n[ ChainEnd2 ]\n{n_atoms}\n")

        render_mdp(
            input.mdp,
            params.mdp_rendered,
            nvt_nsteps=params.nsteps,
            temperature_K=config["temperature_K"],
            pull_target_nm=f"{pull_target_nm:.4f}",
            pull_k=params.pull_k,
        )
        shell(
            "gmx_mpi grompp -f {params.mdp_rendered} -c {output.box_gro} -r {output.box_gro} "
            "-p {input.top} -n {output.ndx} -o {output.tpr} -maxwarn 2"
        )
        shell("gmx_mpi mdrun -deffnm {params.prefix} -ntomp {threads}")


rule vacuum_potential:
    """The vacuum-side half of the cohesive-energy-density split -- see
    bulk_potential for the liquid-side half."""
    input:
        edr=f"{RES}/{{species}}/vacuum_nvt.edr",
    output:
        xvg=f"{RES}/{{species}}/vacuum_potential.xvg",
    params:
        discard_ps=config["bulk"]["vacuum_discard_ps"],
    shell:
        'echo "Potential" | gmx_mpi energy -f {input.edr} -b {params.discard_ps} -o {output.xvg}'


rule solubility_chi:
    """Step 1: cohesive energy density -> Hildebrand solubility parameters
    -> chi_AS, chi_BS."""
    input:
        segment_volume="results/segment_volume.yaml",
        liquid_potential=expand(
            f"{RES}/{{species}}/liquid_potential.xvg", species=BULK_SPECIES
        ),
        vacuum_potential=expand(
            f"{RES}/{{species}}/vacuum_potential.xvg", species=BULK_SPECIES
        ),
        topol=expand(f"{RES}/{{species}}/topol.top", species=BULK_SPECIES),
        itp=expand(
            f"{RES}/{{species}}/{{molname}}.itp",
            zip,
            species=BULK_SPECIES,
            molname=[molname(sp) for sp in BULK_SPECIES],
        ),
    output:
        "results/chi.yaml",
    params:
        species=BULK_SPECIES,
    conda:
        "../envs/analysis.yaml"
    localrule: True
    script:
        "../scripts/solubility_chi.py"
