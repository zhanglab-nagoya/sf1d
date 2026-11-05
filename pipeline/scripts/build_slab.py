#!/usr/bin/env python
"""Build a tiled, z-trimmed slab from a bulk crystal structure (e.g. a CIF).

Cuts the bulk crystal along the given Miller plane, orthogonalizes the
resulting surface cell EXACTLY (see _orthogonalize_cell below -- needed
because many low-index cuts, including hematite's (0001), have a hexagonal
in-plane cell with a=b and a 120 degree angle between them, not a
rectangular one), tiles the orthogonal cell in-plane to reach the target
footprint, trims a centered slice in z to the target thickness (avoiding
cutting through an atomic plane), and adds z vacuum for 2D periodicity. No
force field information (charges, LJ, bonds) is assigned here -- see
write_lammps_charge_data.py for that.
"""
import argparse

import numpy as np
from ase.build import cut, surface
from ase.io import read, write


def _orthogonalize_cell(slab_cell, tol=1e-3, max_k=12):
    """Return (orthogonal_cell, k) with an exactly rectangular a/b in-plane cell.

    `ase.build.surface` returns the PRIMITIVE surface cell, which for many
    low-index cuts (e.g. a hexagonal/rhombohedral (0001) surface, as here)
    has a_vec along x but b_vec skewed by some angle other than 90 degrees.
    Dropping that skew and wrapping atom positions into a rectangular box
    (the naive fix) only reproduces the true periodic structure if the
    accumulated x-skew over the chosen number of b-repeats happens to land
    on an exact multiple of a_vec's length -- otherwise it silently
    introduces a stacking fault at the wrap seam.

    Instead, find the smallest integer k such that k*b_vec's x-component is
    an integer multiple (c) of a_vec's length, so that
    `new_b_vec = k*b_vec + (-c)*a_vec` is exactly perpendicular to a_vec (y
    component only). `ase.build.cut` then constructs the resulting
    `k`-times-larger cell (spanned by a_vec and new_b_vec) from the original
    atoms -- a real orthogonal repeating unit, not an approximation.
    """
    a_vec = slab_cell.get_cell()[0]
    b_vec = slab_cell.get_cell()[1]
    if abs(a_vec[1]) > tol or abs(a_vec[2]) > tol:
        raise RuntimeError(f"Expected a_vec along x; got {a_vec}")
    a_len = a_vec[0]

    skew_ratio = b_vec[0] / a_len
    for k in range(1, max_k + 1):
        c = k * skew_ratio
        if abs(c - round(c)) < tol:
            c = int(round(c))
            break
    else:
        raise RuntimeError(
            f"Could not find an integer combination (up to k={max_k}) that "
            f"orthogonalizes this surface cell (b_x/a_len = {skew_ratio:.6f}); "
            "this cut may not reduce to a rectangular supercell this way."
        )

    ortho_cell = cut(slab_cell, a=(1, 0, 0), b=(-c, k, 0), c=(0, 0, 1))
    new_b_vec = ortho_cell.get_cell()[1]
    if abs(new_b_vec[0]) > tol * a_len:
        raise RuntimeError(f"Orthogonalization failed: new b_vec not along y: {new_b_vec}")
    return ortho_cell, k


def _stoichiometric_z_trim(slab_tiled, target_thickness, plane_tol=0.3):
    """Trim to a z-window of whole atomic planes with exact Fe:O = 2:3.

    Groups atoms into z-planes (atoms within `plane_tol` Angstrom of a
    plane's starting z are the same plane -- true atomic planes are flat to
    well under this for a rigid, untilted cut), then searches all contiguous
    plane windows for the one closest to `target_thickness` whose cumulative
    Fe and O counts satisfy Fe*3 == O*2 exactly (a stoichiometric,
    charge-neutral Fe2O3 slice). Only Fe/O elements are expected.
    """
    positions = slab_tiled.positions
    symbols = slab_tiled.get_chemical_symbols()
    order = np.argsort(positions[:, 2])

    planes = []  # (z_start, {symbol: count})
    plane_start = None
    counts = {}
    for idx in order:
        zi = positions[idx, 2]
        sym = symbols[idx]
        if plane_start is None or zi - plane_start > plane_tol:
            if plane_start is not None:
                planes.append((plane_start, counts))
            plane_start = zi
            counts = {}
        counts[sym] = counts.get(sym, 0) + 1
    if plane_start is not None:
        planes.append((plane_start, counts))

    n_planes = len(planes)
    best = None  # (|span - target|, i, j)
    for i in range(n_planes):
        fe = o = 0
        for j in range(i, n_planes):
            fe += planes[j][1].get("Fe", 0)
            o += planes[j][1].get("O", 0)
            if fe == 0 and o == 0:
                continue
            if fe * 3 == o * 2:
                span = planes[j][0] - planes[i][0]
                score = abs(span - target_thickness)
                if best is None or score < best[0]:
                    best = (score, i, j)

    if best is None:
        raise RuntimeError(
            "No contiguous window of atomic planes is exactly Fe:O = 2:3 "
            "stoichiometric; cannot build a charge-neutral slab this way."
        )
    _, i, j = best
    zlo = planes[i][0] - plane_tol / 2
    zhi = planes[j][0] + plane_tol / 2
    mask = (positions[:, 2] >= zlo) & (positions[:, 2] <= zhi)
    slab = slab_tiled[mask]

    fe_n = sum(1 for s in slab.get_chemical_symbols() if s == "Fe")
    o_n = sum(1 for s in slab.get_chemical_symbols() if s == "O")
    actual_thickness = planes[j][0] - planes[i][0]
    print(
        f"Stoichiometric z-trim: planes [{i}:{j}] of {n_planes}, "
        f"thickness {actual_thickness:.3f} A (target {target_thickness:.3f} A), "
        f"Fe={fe_n} O={o_n} (ratio {fe_n / o_n:.6f}, expect 0.666667)"
    )
    return slab


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("structure", help="Bulk crystal structure (e.g. CIF)")
    parser.add_argument("-o", "--output", required=True, help="Output tiled/trimmed slab (xyz)")
    parser.add_argument("--miller", type=int, nargs=3, default=[0, 0, 1], help="Surface Miller indices")
    parser.add_argument("--layers", type=int, default=12, help="Initial number of atomic layers to cut")
    parser.add_argument(
        "--target-nm", type=float, nargs=3, required=True,
        help="Target footprint and thickness (x y z, nm) before any later relaxation",
    )
    parser.add_argument("--vacuum", type=float, default=10.0, help="Vacuum padding added in z (Angstrom)")
    args = parser.parse_args()

    target = np.array(args.target_nm) * 10.0  # nm -> Angstrom

    bulk = read(args.structure)

    # Build a primitive slab cell (no large xy tiling yet) with enough layers
    # that the local stacking is present.
    slab_cell = surface(bulk, tuple(args.miller), layers=args.layers, vacuum=0.0)
    slab_cell.center(axis=2)

    ortho_cell, k = _orthogonalize_cell(slab_cell)

    a_vec = ortho_cell.get_cell()[0]
    b_vec = ortho_cell.get_cell()[1]
    a_len = a_vec[0]
    b_len = b_vec[1]

    nx = int(np.ceil(target[0] / a_len))
    ny = int(np.ceil(target[1] / b_len))
    print(
        f"Orthogonalized primitive surface cell: {a_len:.3f} x {b_len:.3f} A "
        f"({k}x the raw (0,0,1) primitive cell's area)"
    )
    print(f"Repeating orthogonal cell: nx={nx}, ny={ny}")

    slab_tiled = ortho_cell.repeat((nx, ny, 1))

    # z-trim to (close to) the target thickness, picking a window of whole
    # atomic planes that is EXACTLY Fe:O = 2:3 stoichiometric -- not just one
    # that avoids slicing an atom in half. Corundum's (0001) stacking
    # alternates Fe- and O-only planes, so a plane-aligned cut can still
    # leave a non-stoichiometric (charged) slab if it starts/ends at the
    # wrong point in that alternation; searching for an exactly-neutral
    # window avoids having to hardcode the stacking period.
    slab = _stoichiometric_z_trim(slab_tiled, target[2])

    slab.center(vacuum=args.vacuum, axis=2)
    slab.set_pbc((True, True, False))

    write(args.output, slab)
    print(f"Wrote {args.output}")
    print(f"Cell (A):\n{slab.get_cell()}")
    print(f"Atom count: {len(slab)}")


if __name__ == "__main__":
    main()
