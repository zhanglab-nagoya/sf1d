"""Convert a GROMACS .gro/.top into a LAMMPS input+data file pair via
InterMol's Python API directly, instead of the `intermol-convert` CLI --
the CLI has no way to pass GROMACS preprocessor #defines, and this force
field's L-OPLS atom types (opls_140A/opls_140B, see
forcefield/oplsaa.ff/ffnonbonded.itp) only exist behind `#ifdef LOPLS`,
the same define every mdp/*.mdp file sets for grompp. Without it,
InterMol's own GROMACS parser (which does its own #ifdef/#else resolution,
independent of grompp) silently drops the whole #ifdef LOPLS block and
later crashes with `KeyError: 0` on `atom.epsilon[0]` for every hydrogen
using those two types -- confirmed live 2026-11-05 debugging
rules/step2_dilute.smk's npt.gro conversion.
"""

import argparse

import intermol.gromacs as gmx
import intermol.lammps as lmp


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("top")
    ap.add_argument("gro")
    ap.add_argument(
        "-o", "--oname", required=True,
        help="output prefix (writes <oname>.input, referencing <oname>.input's own data section)",
    )
    ap.add_argument(
        "--define", action="append", default=[],
        help="GROMACS preprocessor define, e.g. LOPLS (repeatable)",
    )
    args = ap.parse_args()

    defines = {name: "" for name in args.define}
    system = gmx.load(top_file=args.top, gro_file=args.gro, defines=defines)
    lmp.save(f"{args.oname}.input", system)


if __name__ == "__main__":
    main()
