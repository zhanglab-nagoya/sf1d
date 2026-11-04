"""Step 0 (segment size) + Step 1 (chi_AS, chi_BS): bulk hexadecane and
capped ST/MMA oligomer "melt" boxes, GROMACS em -> nvt -> npt, then read off
density (step0) and cohesive energy density / solubility parameter (step1).

STATUS: scaffolded, not yet run end-to-end on this machine. The GROMACS
settings (mdp/*.mdp) and box sizes (config.yaml: bulk.*) are reasonable
defaults, not tuned — review run lengths before trusting any density/CED
number this produces. The vacuum single-molecule reference needed for a
correct cohesive-energy-density split (see solubility_chi rule below) is the
one piece that's a real TODO, not just "needs longer runs": the shell
command there must be reviewed, not just executed, before use.
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
        seq=lambda wc: " ".join(["CH3"] + [wc.species] * N_OLIGO + ["CH3n"]),
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
    input:
        itp=lambda wc: f"{RES}/{wc.species}/{molname(wc.species)}.itp",
        ff_itp=f"{FF_DIR}/oplsaa.ff/forcefield.itp",
        nb=f"{FF_DIR}/custom-ffnonbonded.itp",
    output:
        top=f"{RES}/{{species}}/topol.top",
    localrule: True
    params:
        n_chains=lambda wc: config["bulk"]["n_chains"][wc.species],
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
            fh.write(f"{params.mol} bulk box\n\n")
            fh.write("[ molecules ]\n")
            fh.write(f"{params.mol}    {params.n_chains}\n")


rule bulk_pack:
    """Build the box from scratch with polyply (no seed coordinates needed)."""
    input:
        top=f"{RES}/{{species}}/topol.top",
        itp=lambda wc: f"{RES}/{wc.species}/{molname(wc.species)}.itp",
    output:
        gro=f"{RES}/{{species}}/box0.gro",
    conda:
        "../envs/polyply.yaml"
    localrule: True
    params:
        box=config["bulk"]["box_nm"],
    shell:
        "polyply gen_coords -p {input.top} -o {output.gro} "
        "-box {params.box[0]} {params.box[1]} {params.box[2]}"


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
    shell:
        """
        python -c "import sys; sys.path.insert(0,'rules'); from common import render_mdp; \
render_mdp('{input.mdp}', '{params.mdp_rendered}', em_nsteps={config[bulk][em_nsteps]})"
        gmx_mpi grompp -f {params.mdp_rendered} -c {input.gro} -p {input.top} \
            -o {output.tpr} -maxwarn 2
        gmx_mpi mdrun -deffnm {params.prefix} -ntmpi 1
        """


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
    shell:
        """
        python -c "import sys; sys.path.insert(0,'rules'); from common import render_mdp; \
render_mdp('{input.mdp}', '{params.mdp_rendered}', nvt_nsteps={params.nsteps}, \
temperature_K={config[temperature_K]})"
        gmx_mpi grompp -f {params.mdp_rendered} -c {input.gro} -r {input.gro} \
            -p {input.top} -o {output.tpr} -maxwarn 2
        gmx_mpi mdrun -deffnm {params.prefix} -ntmpi 1
        """


rule bulk_npt:
    input:
        gro=f"{RES}/{{species}}/nvt.gro",
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
    shell:
        """
        python -c "import sys; sys.path.insert(0,'rules'); from common import render_mdp; \
render_mdp('{input.mdp}', '{params.mdp_rendered}', npt_nsteps={params.nsteps}, \
temperature_K={config[temperature_K]}, pressure_bar={config[pressure_bar]})"
        gmx_mpi grompp -f {params.mdp_rendered} -c {input.gro} -r {input.gro} \
            -p {input.top} -o {output.tpr} -maxwarn 2
        gmx_mpi mdrun -deffnm {params.prefix} -ntmpi 1
        """


rule bulk_density:
    """Step 0 input: average mass density -> a^3 (segment volume mapping)."""
    input:
        edr=f"{RES}/{{species}}/npt.edr",
    output:
        xvg=f"{RES}/{{species}}/density.xvg",
    shell:
        'echo "Density" | gmx_mpi energy -f {input.edr} -o {output.xvg}'


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


rule vacuum_ref:
    """TODO (review before trusting): single isolated molecule, same force
    field, large box, no PBC images close enough to interact — gives the
    intramolecular reference energy needed to split liquid U_total into
    U_intra + U_inter for the cohesive-energy-density (step 1) calculation.
    Not yet wired up to bulk_em/nvt/npt; needs its own minimal mdp (vacuum
    electrostatics, e.g. plain cutoff, not PME) and its own small rule chain
    mirroring bulk_* above but with n_chains=1 and a big (~10 nm) box."""
    output:
        touch("results/bulk/{species}/vacuum_ref.TODO"),
    localrule: True
    shell:
        "echo 'vacuum reference not implemented — see rule docstring' && false"


rule solubility_chi:
    """Step 1: cohesive energy density -> Hildebrand solubility parameters
    -> chi_AS, chi_BS. Depends on vacuum_ref (not implemented yet)."""
    input:
        segment_volume="results/segment_volume.yaml",
        vacuum_ref=expand(
            "results/bulk/{species}/vacuum_ref.TODO", species=BULK_SPECIES
        ),
    output:
        "results/chi.yaml",
    conda:
        "../envs/analysis.yaml"
    localrule: True
    script:
        "../scripts/solubility_chi.py"
