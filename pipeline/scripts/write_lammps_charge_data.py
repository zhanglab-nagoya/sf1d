#!/usr/bin/env python
"""Assign per-element point charges and write a LAMMPS (atom_style full) data file.

Charges (and, downstream, the LJ pair coefficients applied by LAMMPS in a
later rule) come from the INTERFACE force field for Fe2O3: Kanhaiya,
Krishan, et al. "Accurate force fields for atomistic simulations of oxides,
hydroxides, and organic hybrid materials up to the micrometer scale."
Journal of Chemical Theory and Computation 19.22 (2023): 8293-8322.
"""
import argparse

from ase.io import read, write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Orthorhombic slab (xyz)")
    parser.add_argument("-o", "--output", required=True, help="LAMMPS data file (atom_style full)")
    parser.add_argument("--charge-fe", type=float, required=True, help="Fe point charge (e)")
    parser.add_argument("--charge-o", type=float, required=True, help="O point charge (e)")
    args = parser.parse_args()

    atoms = read(args.input)

    charge_by_symbol = {"Fe": args.charge_fe, "O": args.charge_o}
    charges = []
    for atom in atoms:
        if atom.symbol not in charge_by_symbol:
            raise ValueError(f"Unexpected element {atom.symbol!r} in {args.input}")
        charges.append(charge_by_symbol[atom.symbol])
    atoms.set_initial_charges(charges)

    write(args.output, atoms, format="lammps-data", atom_style="full")

    n_fe = sum(1 for a in atoms if a.symbol == "Fe")
    n_o = sum(1 for a in atoms if a.symbol == "O")
    print(f"Wrote {args.output}: {len(atoms)} atoms ({n_fe} Fe, {n_o} O)")


if __name__ == "__main__":
    main()
