"""Step 2 dilute box: one ST-analog/MMA-analog fragment diluted in a normal
periodic hexadecane box, built and equilibrated entirely in GROMACS -- no
substrate here at all. The resulting snapshot is converted to LAMMPS and
merged onto results/substrate/substrate_minimized.data afterward (see
step2_merge.smk), not built next to the substrate directly: GROMACS has no
Fe2O3/INTERFACE force field, so there is nothing for it to equilibrate
against anyway. This mirrors step2_chis.smk's docstring's own stage 1/2
split (build_dilute_fragment_box, then convert_to_lammps).

The fragment is deliberately NOT the ST10/MMA10 10-mer oligomer step0/1
uses for bulk CED -- see config.yaml's slab.fragment_length comment for
why a single capped repeat unit (fragment_length=1) is the right analog
for a per-segment chis, not the whole oligomer. Reuses mdp/em.mdp,
mdp/nvt.mdp, mdp/npt.mdp verbatim (same force field, same settings) and
results/bulk/HEX/single_em.gro directly (hexadecane's EM'd single-molecule
structure is shared, identical regardless of which stage needs it -- no
reason to rebuild it here).

Packs close to hexadecane's own measured equilibrium density (n_hex in
config.yaml, sized from results/segment_volume.yaml's measured HEX
density) with only a modest z-margin, not bulk_pack's full
oversize-then-squeeze -- that whole problem was about packmol reliably
placing hundreds of ST10's long, floppy C18-tailed copies. But "modest
margin" turned out not to be zero margin: confirmed live 2026-11-05,
packing dead-on at the real density was too tight for packmol's optimizer
even for one small fragment (ST barely converged at loop 839/1000 after
~87 min; MMA exhausted all 1000 loops without converging) -- the same
"packing near the true density is the hard case" lesson from bulk_pack's
own box_nm comment, just at a smaller scale. slab.volume_margin
(config.yaml) now carries a ~12% margin instead; dilute_npt's z-only
relaxation easily absorbs that small a correction.

The box's lateral (x,y) dimensions are NOT sized independently at all --
dilute_pack reads them straight from results/substrate/
substrate_minimized.data (step2_substrate.smk), so this organic box ends
up commensurate with the Fe2O3 substrate's own relaxed footprint *before*
any MD here runs, rather than merging two differently-sized periodic
boxes later (step2_convert.smk/step2_chis.smk) and having to stretch one
to fit -- which would strain either the solvent's density or the
substrate's own crystal lattice for no good reason. Only z is free:
dilute_npt (mdp/npt_semiiso.mdp, semiisotropic pressure coupling with
compressibility_xy=0) relaxes z only, keeping x,y frozen at whatever
dilute_pack set them to. mdp/em.mdp and mdp/nvt.mdp are still reused
verbatim from step0/1 (no coupling-type concern at constant volume).

STATUS: runs end-to-end for both species with the substrate-commensurate
box above -- confirmed live 2026-11-06: lateral box frozen exactly at the
substrate's own relaxed footprint (5.1811 x 5.3844 nm, both species, as
intended), z relaxes to ~5.42-5.44 nm, density ~743-745 kg/m^3 (ST 743.1,
MMA 745.1), close to hexadecane's own ~750 kg/m^3 reference as expected.
"""

SLAB_RES = "results/slab"
FRAG_LEN = config["slab"]["fragment_length"]
FRAG_SPECIES = ["ST", "MMA"]


def frag_molname(species):
    return f"{species}{FRAG_LEN}"


rule build_fragment_itp:
    """Single capped repeat unit, CH3-{species}-CH3n -- see this file's
    docstring for why fragment_length=1, not oligomer_length=10."""
    input:
        ff=[
            f"{FF_DIR}/monomers/PAMA_{{species}}.oplsaa.LigParGen.ff",
            f"{FF_DIR}/monomers/CH3.ff",
            f"{FF_DIR}/monomers/CH3n.ff",
            f"{FF_DIR}/PAMA_ST_MMA.oplsaa.LigParGen_links.ff",
        ],
    output:
        itp=f"{SLAB_RES}/{{species}}/{{species}}{FRAG_LEN}.itp",
    wildcard_constraints:
        species="ST|MMA",
    conda:
        "../envs/polyply.yaml"
    localrule: True
    params:
        seq=lambda wc: f"CH3:1 {wc.species}:{FRAG_LEN} CH3n:1",
        name=lambda wc: f"{wc.species}{FRAG_LEN}",
    shell:
        "polyply gen_params -name {params.name} -f {input.ff} "
        "-o {output.itp} -seq {params.seq}"


rule fragment_single_top:
    """n=1 topology, for resolving the fragment's own geometry before
    mixing it into the hexadecane box (same reasoning as step0/1's
    single_top -- polyply's gen_coords needs a clean single-molecule EM
    pass regardless of what it's eventually mixed with)."""
    input:
        itp=lambda wc: f"{SLAB_RES}/{wc.species}/{frag_molname(wc.species)}.itp",
        ff_itp=f"{FF_DIR}/oplsaa.ff/forcefield.itp",
        nb=f"{FF_DIR}/custom-ffnonbonded.itp",
    output:
        top=f"{SLAB_RES}/{{species}}/single.top",
    localrule: True
    params:
        mol=lambda wc: frag_molname(wc.species),
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


rule fragment_single_mol_raw:
    input:
        top=f"{SLAB_RES}/{{species}}/single.top",
        itp=lambda wc: f"{SLAB_RES}/{wc.species}/{frag_molname(wc.species)}.itp",
    output:
        gro=f"{SLAB_RES}/{{species}}/single0.gro",
    wildcard_constraints:
        species="ST|MMA",
    conda:
        "../envs/polyply.yaml"
    localrule: True
    shell:
        "polyply gen_coords -p {input.top} -o {output.gro} -box 5.0 5.0 5.0"


rule fragment_single_mol_em:
    """Resolve the fragment's real geometry before mixing -- polyply's
    gen_coords gives wrong local geometry regardless of molecule size (see
    step0_step1_bulk.smk's single_mol_em docstring); plain unconstrained
    steepest descent with mdp/em.mdp's small step fixes it the same way."""
    input:
        gro=f"{SLAB_RES}/{{species}}/single0.gro",
        top=f"{SLAB_RES}/{{species}}/single.top",
        mdp="mdp/em.mdp",
    output:
        gro=f"{SLAB_RES}/{{species}}/single_em.gro",
        tpr=f"{SLAB_RES}/{{species}}/single_em.tpr",
    params:
        prefix=f"{SLAB_RES}/{{species}}/single_em",
        mdp_rendered=f"{SLAB_RES}/{{species}}/single_em.rendered.mdp",
    threads: 4
    run:
        render_mdp(input.mdp, params.mdp_rendered, em_nsteps=config["slab"]["em_nsteps"])
        shell(
            "gmx_mpi grompp -f {params.mdp_rendered} -c {input.gro} -p {input.top} "
            "-o {output.tpr} -maxwarn 2"
        )
        shell("gmx_mpi mdrun -deffnm {params.prefix} -ntomp {threads}")


rule dilute_pack:
    """Pack 1 EM'd fragment + n_hex copies of hexadecane's own EM'd single
    molecule (results/bulk/HEX/single_em.gro -- shared with step0/1, not
    rebuilt here) directly near the real liquid density via packmol. See
    bulk_pack's docstring for the PDB round-trip packmol needs (.gro <->
    .pdb, box vectors reattached after).

    Lateral (x,y) box dims come straight from the substrate's own relaxed
    footprint (see this file's module docstring for why); only z is
    sized here, from n_hex/HEX's measured molar mass+density
    (results/segment_volume.yaml) plus slab.volume_margin headroom for
    packmol -- same role bulk_pack's box_nm oversize played, just
    confined to one dimension instead of three since x,y are fixed."""
    input:
        frag_gro=f"{SLAB_RES}/{{species}}/single_em.gro",
        hex_gro="results/bulk/HEX/single_em.gro",
        substrate_data=config["slab"]["substrate_data"],
        segment_volume="results/segment_volume.yaml",
    output:
        gro=f"{SLAB_RES}/{{species}}/box0.gro",
    params:
        n_hex=config["slab"]["n_hex"],
        volume_margin=config["slab"]["volume_margin"],
        frag_pdb=f"{SLAB_RES}/{{species}}/single_em.pdb",
        hex_pdb=f"{SLAB_RES}/{{species}}/hex_single_em.pdb",
        packed_pdb=f"{SLAB_RES}/{{species}}/packed.pdb",
        inp_file=f"{SLAB_RES}/{{species}}/packmol.inp",
    run:
        import yaml

        N_A = 6.02214076e23  # 1/mol

        def substrate_xy_nm(data_path):
            xlo = xhi = ylo = yhi = None
            with open(data_path) as fh:
                for line in fh:
                    if "xlo xhi" in line:
                        xlo, xhi = map(float, line.split()[:2])
                    elif "ylo yhi" in line:
                        ylo, yhi = map(float, line.split()[:2])
            return (xhi - xlo) / 10.0, (yhi - ylo) / 10.0  # Angstrom -> nm

        x_nm, y_nm = substrate_xy_nm(input.substrate_data)

        with open(input.segment_volume) as fh:
            seg = yaml.safe_load(fh)
        m_hex = seg["hexadecane_molar_mass_g_per_mol"]
        rho_hex = seg["hexadecane_density_kg_per_m3"]

        # n_hex HEX molecules' own volume at their real measured density
        # (the 1 fragment's volume is a negligible correction here --
        # this is only sizing the packing box, not the production result).
        volume_nm3 = (params.n_hex * m_hex / (N_A * 1000.0) / rho_hex) * 1e27
        z_nm = volume_nm3 / (x_nm * y_nm) * params.volume_margin

        shell("gmx_mpi editconf -f {input.frag_gro} -o {params.frag_pdb}")
        shell("gmx_mpi editconf -f {input.hex_gro} -o {params.hex_pdb}")

        box_ang = [x_nm * 10, y_nm * 10, z_nm * 10]
        with open(params.inp_file, "w") as fh:
            fh.write("tolerance 2.0\n")
            fh.write("filetype pdb\n")
            fh.write("nloop 1000\n")
            fh.write(f"output {params.packed_pdb}\n\n")
            fh.write(f"structure {params.frag_pdb}\n")
            fh.write("  number 1\n")
            fh.write(
                f"  inside box 0. 0. 0. {box_ang[0]} {box_ang[1]} {box_ang[2]}\n"
            )
            fh.write("end structure\n\n")
            fh.write(f"structure {params.hex_pdb}\n")
            fh.write(f"  number {params.n_hex}\n")
            fh.write(
                f"  inside box 0. 0. 0. {box_ang[0]} {box_ang[1]} {box_ang[2]}\n"
            )
            fh.write("end structure\n")

        shell("packmol < {params.inp_file}")
        shell(
            f"gmx_mpi editconf -f {params.packed_pdb} -o {output.gro} "
            f"-box {x_nm:.4f} {y_nm:.4f} {z_nm:.4f}"
        )


rule dilute_top:
    """Written after dilute_pack, deriving both molecule counts from the
    packed box's real atom count -- same non-determinism safety net as
    bulk_top (see its docstring), applied to both species present here."""
    input:
        frag_itp=lambda wc: f"{SLAB_RES}/{wc.species}/{frag_molname(wc.species)}.itp",
        hex_itp="results/bulk/HEX/HEX.itp",
        ff_itp=f"{FF_DIR}/oplsaa.ff/forcefield.itp",
        nb=f"{FF_DIR}/custom-ffnonbonded.itp",
        box=f"{SLAB_RES}/{{species}}/box0.gro",
        frag_single=f"{SLAB_RES}/{{species}}/single_em.gro",
        hex_single="results/bulk/HEX/single_em.gro",
    output:
        top=f"{SLAB_RES}/{{species}}/topol.top",
    localrule: True
    params:
        frag_mol=lambda wc: frag_molname(wc.species),
        n_hex_requested=config["slab"]["n_hex"],
    run:
        import os
        import sys

        def count_atoms(gro_path):
            with open(gro_path) as fh:
                fh.readline()
                return int(fh.readline())

        n_frag_atoms = count_atoms(input.frag_single)
        n_hex_atoms = count_atoms(input.hex_single)
        n_box_atoms = count_atoms(input.box)
        n_hex = (n_box_atoms - n_frag_atoms) // n_hex_atoms
        if n_hex != params.n_hex_requested:
            print(
                f"dilute_top: packmol placed {n_hex} of "
                f"{params.n_hex_requested} requested HEX molecules -- "
                "using the actual count.",
                file=sys.stderr,
            )

        ff_itp_rel = os.path.relpath(input.ff_itp, os.path.dirname(output.top))
        nb_rel = os.path.relpath(input.nb, os.path.dirname(output.top))
        frag_itp_rel = os.path.relpath(input.frag_itp, os.path.dirname(output.top))
        hex_itp_rel = os.path.relpath(input.hex_itp, os.path.dirname(output.top))
        with open(output.top, "w") as fh:
            fh.write(f'#include "{ff_itp_rel}"\n')
            fh.write(f'#include "{nb_rel}"\n')
            fh.write(f'#include "{frag_itp_rel}"\n')
            fh.write(f'#include "{hex_itp_rel}"\n\n')
            fh.write("[ system ]\n")
            fh.write(f"{params.frag_mol} dilute in HEX\n\n")
            fh.write("[ molecules ]\n")
            fh.write(f"{params.frag_mol}    1\n")
            fh.write(f"HEX    {n_hex}\n")


rule dilute_em:
    input:
        gro=f"{SLAB_RES}/{{species}}/box0.gro",
        top=f"{SLAB_RES}/{{species}}/topol.top",
        mdp="mdp/em.mdp",
    output:
        gro=f"{SLAB_RES}/{{species}}/em.gro",
        tpr=f"{SLAB_RES}/{{species}}/em.tpr",
    params:
        prefix=f"{SLAB_RES}/{{species}}/em",
        mdp_rendered=f"{SLAB_RES}/{{species}}/em.rendered.mdp",
    threads: 8
    run:
        render_mdp(input.mdp, params.mdp_rendered, em_nsteps=config["slab"]["em_nsteps"])
        shell(
            "gmx_mpi grompp -f {params.mdp_rendered} -c {input.gro} -p {input.top} "
            "-o {output.tpr} -maxwarn 2"
        )
        shell("gmx_mpi mdrun -deffnm {params.prefix} -ntomp {threads}")


rule dilute_nvt:
    input:
        gro=f"{SLAB_RES}/{{species}}/em.gro",
        top=f"{SLAB_RES}/{{species}}/topol.top",
        mdp="mdp/nvt.mdp",
    output:
        gro=f"{SLAB_RES}/{{species}}/nvt.gro",
        tpr=f"{SLAB_RES}/{{species}}/nvt.tpr",
    params:
        prefix=f"{SLAB_RES}/{{species}}/nvt",
        mdp_rendered=f"{SLAB_RES}/{{species}}/nvt.rendered.mdp",
        nsteps=lambda wc: int(config["slab"]["nvt_ps"] / 0.002),
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


rule dilute_npt:
    """Packed directly near the real density in z (unlike bulk_pack), so
    no squeeze phase needed -- this is the "small perturbation" regime a
    normal ref_p=1 bar barostat handles fine (see config.yaml's bulk.
    squeeze_pressure_bar comment for why that was NOT true for bulk_pack's
    deliberately oversized boxes). Uses mdp/npt_semiiso.mdp, not the
    shared isotropic mdp/npt.mdp step0/1 use: x,y must stay exactly at
    whatever dilute_pack set them to (the substrate's own footprint, see
    dilute_pack's docstring), so only z is allowed to relax."""
    input:
        gro=f"{SLAB_RES}/{{species}}/nvt.gro",
        top=f"{SLAB_RES}/{{species}}/topol.top",
        mdp="mdp/npt_semiiso.mdp",
    output:
        gro=f"{SLAB_RES}/{{species}}/npt.gro",
        edr=f"{SLAB_RES}/{{species}}/npt.edr",
        tpr=f"{SLAB_RES}/{{species}}/npt.tpr",
    params:
        prefix=f"{SLAB_RES}/{{species}}/npt",
        mdp_rendered=f"{SLAB_RES}/{{species}}/npt.rendered.mdp",
        nsteps=lambda wc: int(config["slab"]["npt_ps"] / 0.002),
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
