#!/usr/bin/env python
"""Snap a slab cell to exactly orthorhombic (floating-point cleanup) and wrap atoms.

build_slab.py already constructs a crystallographically exact orthogonal
cell (not an approximation), but repeated centering/tiling can leave tiny
(~1e-10 A) off-diagonal residues; this drops them (keeping only the
diagonal lattice lengths) and wraps all atom positions back into the
resulting box. Also reports any undercoordinated atoms (a quick sanity
check on the cut/wrap, not a correctness guarantee) using a simple
distance-cutoff neighbor count.
"""
import argparse
from collections import Counter

from ase.io import read, write
from ase.neighborlist import neighbor_list


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Tiled/trimmed slab (xyz)")
    parser.add_argument("-o", "--output", required=True, help="Orthorhombic slab (xyz)")
    parser.add_argument(
        "--coord-cutoff", type=float, default=3.0,
        help="Cutoff (Angstrom) for the coordination sanity check (Fe-O nearest-neighbour scale)",
    )
    parser.add_argument(
        "--min-coord", type=int, default=2,
        help="Flag atoms with fewer neighbors than this within --coord-cutoff",
    )
    args = parser.parse_args()

    slab = read(args.input)

    Lx, Ly, Lz = slab.cell[0][0], slab.cell[1][1], slab.cell[2][2]
    slab.set_cell([[Lx, 0, 0], [0, Ly, 0], [0, 0, Lz]])

    positions = slab.get_positions()
    positions[:, 0] %= Lx
    positions[:, 1] %= Ly
    positions[:, 2] %= Lz
    slab.set_positions(positions)

    write(args.output, slab)
    print(f"Wrote {args.output}")
    print(f"Cell (A):\n{slab.get_cell()}")

    i, _, _ = neighbor_list("ijd", slab, args.coord_cutoff)
    cnt = Counter(i)
    under = [atom.index for atom in slab if cnt.get(atom.index, 0) < args.min_coord]
    print(
        f"Possibly undercoordinated atoms (< {args.min_coord} neighbors within "
        f"{args.coord_cutoff} A): {len(under)} (first 20: {under[:20]})"
    )


if __name__ == "__main__":
    main()
