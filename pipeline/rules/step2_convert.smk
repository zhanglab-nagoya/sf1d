"""Step 2 chis stage 2: convert each species' equilibrated dilute-box
GROMACS snapshot (results/slab/{species}/npt.gro, see step2_dilute.smk)
into a LAMMPS input+data file, ready to merge onto the Fe2O3 substrate
(step2_merge.smk, not yet written).

Goes through InterMol's Python API directly (scripts/gmx_to_lammps.py),
not the `intermol-convert` CLI: the CLI has no way to pass GROMACS
preprocessor #defines, and this force field's L-OPLS atom types
(opls_140A/opls_140B, see forcefield/oplsaa.ff/ffnonbonded.itp) only exist
behind `#ifdef LOPLS` -- the same define every mdp/*.mdp file sets for
grompp. Without it, InterMol's own GROMACS parser (which does its own
#ifdef/#else resolution, independent of grompp) silently drops that whole
block and crashes with `KeyError: 0` on atom.epsilon[0] for every
hydrogen using those two types -- confirmed live 2026-11-05 (see
scripts/gmx_to_lammps.py's own docstring for the full diagnosis).

Cross-checked against GROMACS: a single-point potential energy
comparison (gmx mdrun -rerun on the unwrapped structure vs. a LAMMPS
`run 0` on the converted files) agreed to ~0.5% for ST (33296.8 kJ/mol
GROMACS vs. 7997.62 kcal/mol = 33474.0 kJ/mol LAMMPS) -- the residual is
expected from the different long-range electrostatics grids/cutoffs
(GROMACS PME vs. LAMMPS PPPM, not tuned to match exactly here) rather
than a conversion bug.

convert_to_lammps runs `module load intermol` itself, scoped to its own
shell subprocess via `bash -lc` -- confirmed live 2026-11-05 this must NOT
be done globally before invoking snakemake (unlike gmx_mpi/packmol/lmp's
plain-PATH convention): `module load intermol` points PYTHONPATH at a
spack-built numpy for Python 3.14, and if that leaks into the whole
snakemake process (python 3.13, its own conda env), Snakemake's own
numpy/pandas imports break with an ABI mismatch
("No module named 'numpy._core._multiarray_umath'") when it tries to
serialize job parameters for its own run metadata -- nothing to do with
this rule's actual work.

STATUS: scaffolded and validated (energy cross-check above) for ST; not
yet run for MMA.
"""

SLAB_RES = "results/slab"  # shared with step2_dilute.smk


rule gmx_trjconv_whole:
    """Unwrap the equilibrated NPT snapshot across PBC -- InterMol's own
    GROMACS parser does not rebuild whole molecules from a wrapped .gro
    the way GROMACS itself does internally, so a wrapped structure would
    hand LAMMPS molecules split across box images."""
    input:
        gro=f"{SLAB_RES}/{{species}}/npt.gro",
        tpr=f"{SLAB_RES}/{{species}}/npt.tpr",
    output:
        gro=f"{SLAB_RES}/{{species}}/unwrapped.gro",
    wildcard_constraints:
        species="ST|MMA",
    localrule: True
    shell:
        "echo 0 | gmx_mpi trjconv -f {input.gro} -s {input.tpr} -pbc whole -o {output.gro}"


rule convert_to_lammps:
    """See scripts/gmx_to_lammps.py's docstring for why this goes through
    InterMol's Python API rather than its CLI (the --define LOPLS below
    is the whole reason).

    Also deduplicates bond/angle/dihedral types right after conversion
    (scripts/dedup_lammps_bonded_types.py) -- InterMol gives every
    individual bond/angle/dihedral its own unique type (e.g. 14773 "bond
    types" for 14773 bonds, confirmed live 2026-11-06 to collapse to just
    8 actually-distinct parameter sets), which breaks LAMMPS's own
    restart-file serialization for `bond_style hybrid` once the type
    count gets this large (needed downstream by
    step2_equilibrate.smk's minimize->NVT handoff) -- see that script's
    docstring."""
    input:
        top=f"{SLAB_RES}/{{species}}/topol.top",
        gro=f"{SLAB_RES}/{{species}}/unwrapped.gro",
    output:
        lmp_input=f"{SLAB_RES}/{{species}}/converted.input",
        lmp_data=f"{SLAB_RES}/{{species}}/converted.lmp",
    params:
        oname=f"{SLAB_RES}/{{species}}/converted",
        raw_lmp=f"{SLAB_RES}/{{species}}/converted_raw.lmp",
    localrule: True
    shell:
        "bash -lc 'module load intermol >/dev/null 2>&1 && "
        "python scripts/gmx_to_lammps.py {input.top} {input.gro} "
        "-o {params.oname} --define LOPLS' && "
        "mv {output.lmp_data} {params.raw_lmp} && "
        "python scripts/dedup_lammps_bonded_types.py {params.raw_lmp} -o {output.lmp_data}"
