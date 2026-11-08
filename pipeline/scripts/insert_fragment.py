"""Move the single adsorbing fragment (molecule 1-3, see
lmp/_common_setup.lmp.template's comment) to the middle of the
z-equilibrated liquid arc, removing whatever hexadecane molecules it
would otherwise overlap.

Replaces the earlier approach (LAMMPS's own `displace_atoms` followed by
a blind minimize) -- confirmed live 2026-11-08 that a rigid translation
straight into occupied solvent produces a real ~3e13 kcal/mol transient
(a near-complete atomic overlap) that the minimizer happens to resolve,
but that's luck, not design. This script instead explicitly finds every
hexadecane atom within a clash distance of the (already-translated)
fragment and drops those molecules whole (never a partial molecule --
that would leave a dangling bond) before any LAMMPS run touches the
result, so the one remaining minimize (lmp_insert_minimize) only has to
relax a much smaller, localized disturbance.

All fragment-placement geometry lives here now, not split with a LAMMPS
`variable` block the way an earlier draft had it -- see
liquid_arc_target()'s own docstring for the gap-based method this uses
to find the bulk-liquid arc's midpoint, robust to the z-only barostat
(lmp/npt_z.lmp.template) shifting the box's zlo into what used to be
the substrate's own range.

Run as a plain script (not a snakemake `script:` block), matching
merge_slab.py's own style, so it's independently testable from the
shell:

    python scripts/insert_fragment.py npt_equilibrated.data -o inserted.data \\
        --n-sub-types 2 --clash-angstrom 4.0
"""

import argparse
import re

import numpy as np
from scipy.spatial import cKDTree

SECTION_NAMES = [
    "Masses", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs",
    "Dihedral Coeffs", "Improper Coeffs", "Atoms", "Velocities",
    "Bonds", "Angles", "Dihedrals", "Impropers",
]
COUNT_RE = re.compile(r"^\s*(\d+)\s+(atoms|bonds|angles|dihedrals|impropers|"
                       r"atom types|bond types|angle types|dihedral types|"
                       r"improper types)\s*$")
BOX_RE = re.compile(r"^\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+(xlo xhi|ylo yhi|zlo zhi)\s*$")


def read_data_sections(path):
    """Same small parser merge_slab.py already uses: (header counts dict,
    box (xlo,xhi,ylo,yhi,zlo,zhi) dict, {section name: [raw lines]})."""
    with open(path) as fh:
        lines = fh.readlines()

    header = {}
    box = {}
    i, n = 0, len(lines)
    while i < n:
        line = lines[i].rstrip("\n")
        stripped = line.split("#")[0].strip()
        if stripped in SECTION_NAMES:
            break
        m = COUNT_RE.match(line)
        if m:
            header[m.group(2)] = int(m.group(1))
        m = BOX_RE.match(line)
        if m:
            key = m.group(3).split()[0]  # "xlo", "ylo", "zlo"
            box[key] = (float(m.group(1)), float(m.group(2)))
        i += 1

    sections = {}
    current = None
    for line in lines[i:]:
        stripped = line.rstrip("\n").split("#")[0].strip()
        if stripped in SECTION_NAMES:
            current = stripped
            sections[current] = []
            continue
        if current is None or line.strip() == "":
            continue
        sections[current].append(line.rstrip("\n"))

    return header, box, sections


def liquid_arc_target(sections, n_org_types, zlo, zhi):
    """(arc_length, target_z): the bulk-liquid arc's length and its
    absolute-z midpoint. NOT a naive max(z)-min(z) on the substrate's own
    atoms -- confirmed live 2026-11-08 that after the z-only barostat
    (lmp/npt_z.lmp.template) shrinks the box, the new zlo can land
    *inside* where the substrate used to sit (the substrate itself is
    never dilated, but the box boundary moving around it is enough to
    push some of its atoms outside [zlo, zhi), where LAMMPS's own
    periodic wrap folds them back in near the opposite edge) -- a plain
    min/max over those atoms' z then spans almost the entire box, not
    the substrate's real ~30 A thickness.

    Instead: fold every substrate atom's z into the box-relative
    coordinate s = (z - zlo) mod L and sort. The substrate is one dense,
    lattice-packed block (confirmed live 2026-11-08, internal gaps
    ~0.84 A for both species), so by far the largest gap between
    consecutive sorted s values -- treated circularly, including the
    wrap from the last value back to the first -- is the liquid region
    itself, found this way regardless of whether the substrate happens
    to straddle the s=0/L seam. That gap's own size is exactly the
    liquid arc's length, and its midpoint converts straight back to an
    absolute target z."""
    L = zhi - zlo
    s = sorted((float(line.split()[6]) - zlo) % L for line in sections["Atoms"]
               if int(line.split()[2]) > n_org_types)
    n = len(s)
    gaps = [(s[(i + 1) % n] - s[i]) % L for i in range(n)]
    k = max(range(n), key=lambda i: gaps[i])
    arc_length = gaps[k]
    target_s = (s[k] + arc_length / 2.0) % L
    target_z = zlo + target_s
    return arc_length, target_z


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("data", help="z-equilibrated LAMMPS data file (e.g. npt_equilibrated.data)")
    ap.add_argument("-o", "--out-data", required=True)
    ap.add_argument("--n-sub-types", type=int, default=2, help="substrate atom types (Fe, O)")
    ap.add_argument("--clash-angstrom", type=float, default=4.0,
                     help="hexadecane atoms within this distance of a (translated) fragment atom "
                          "get their whole molecule removed")
    args = ap.parse_args()

    header, box, sections = read_data_sections(args.data)
    n_atom_types = header["atom types"]
    n_org_types = n_atom_types - args.n_sub_types

    atoms = [line.split() for line in sections["Atoms"]]
    # atom_style full: id mol type q x y z
    atom_by_id = {int(f[0]): f for f in atoms}
    masses = {int(line.split()[0]): float(line.split()[1]) for line in sections["Masses"]}

    lx = box["xlo"][1] - box["xlo"][0]
    ly = box["ylo"][1] - box["ylo"][0]
    lz = box["zlo"][1] - box["zlo"][0]
    # zlo need not be 0 here (unlike merge_slab.py's own output): the
    # upstream z-only barostat (lmp/npt_z.lmp.template) can shift zlo as
    # it rescales the box -- confirmed live 2026-11-08. Everything below
    # is already expressed in absolute box coordinates (the cKDTree
    # origin, target_z), so this needs no special-casing.
    arc_length, target_z = liquid_arc_target(sections, n_org_types, box["zlo"][0], box["zlo"][1])

    frag_ids = [int(f[0]) for f in atoms if int(f[1]) <= 3]
    frag_mass = np.array([masses[int(atom_by_id[i][2])] for i in frag_ids])
    frag_z = np.array([float(atom_by_id[i][6]) for i in frag_ids])
    frag_cz = float(np.sum(frag_mass * frag_z) / np.sum(frag_mass))
    dz = target_z - frag_cz

    for i in frag_ids:
        f = atom_by_id[i]
        f[6] = f"{float(f[6]) + dz:.7f}"

    # --- find hexadecane molecules clashing with the (now-moved) fragment ---
    # molid > 3 alone is NOT enough to mean "hexadecane" -- the substrate
    # also has molid > 3 (a single id, max_org_molid + 1, see
    # merge_slab.py) and must be excluded by atom type too, or its whole
    # 8400 atoms get flagged as "clashing" and deleted (confirmed live
    # 2026-11-08: a first pass without the type check reported removing
    # "3 hexadecane molecules (8500 atoms)" -- impossible for a
    # 50-atom-per-molecule species, and 8400 + 2*50 = 8500 exactly).
    hex_atoms = [f for f in atoms if int(f[1]) > 3 and int(f[2]) <= n_org_types]
    hex_ids = np.array([int(f[0]) for f in hex_atoms])
    hex_molids = np.array([int(f[1]) for f in hex_atoms])
    hex_xyz = np.array([[float(f[4]), float(f[5]), float(f[6])] for f in hex_atoms])
    # cKDTree's periodic boxsize needs every coordinate folded into [0, L);
    # merge_slab.py's box has a small negative xlo/ylo, so shift into that frame.
    origin = np.array([box["xlo"][0], box["ylo"][0], box["zlo"][0]])
    boxsize = np.array([lx, ly, lz])
    hex_frac = np.mod(hex_xyz - origin, boxsize)

    frag_xyz = np.array([[float(atom_by_id[i][4]), float(atom_by_id[i][5]), float(atom_by_id[i][6])] for i in frag_ids])
    frag_frac = np.mod(frag_xyz - origin, boxsize)

    tree = cKDTree(hex_frac, boxsize=boxsize)
    clashing = tree.query_ball_point(frag_frac, r=args.clash_angstrom)
    clashing_hex_idx = sorted({idx for group in clashing for idx in group})
    removed_molids = set(hex_molids[clashing_hex_idx])
    removed_atom_ids = {int(hex_ids[i]) for i in range(len(hex_ids)) if hex_molids[i] in removed_molids}

    kept_atoms = [f for f in atoms if int(f[0]) not in removed_atom_ids]
    kept_vel = [
        line for line in sections.get("Velocities", [])
        if int(line.split()[0]) not in removed_atom_ids
    ]

    def filter_topology(section_name, n_atom_cols):
        kept = []
        for line in sections.get(section_name, []):
            fields = line.split()
            if any(int(fields[2 + k]) in removed_atom_ids for k in range(n_atom_cols)):
                continue
            kept.append(line)
        return kept

    kept_bonds = filter_topology("Bonds", 2)
    kept_angles = filter_topology("Angles", 3)
    kept_dihedrals = filter_topology("Dihedrals", 4)

    with open(args.out_data, "w") as fh:
        fh.write("Fragment moved to bulk-arc midpoint, overlapping hexadecane removed (insert_fragment.py)\n\n")
        fh.write(f"{len(kept_atoms)} atoms\n")
        fh.write(f"{len(kept_bonds)} bonds\n")
        fh.write(f"{len(kept_angles)} angles\n")
        fh.write(f"{len(kept_dihedrals)} dihedrals\n")
        fh.write(f"{header.get('impropers', 0)} impropers\n\n")
        fh.write(f"{n_atom_types} atom types\n")
        for key in ("bond types", "angle types", "dihedral types"):
            if key in header:
                fh.write(f"{header[key]} {key}\n")
        fh.write("\n")
        fh.write(f"{box['xlo'][0]:.7f} {box['xlo'][1]:.7f} xlo xhi\n")
        fh.write(f"{box['ylo'][0]:.7f} {box['ylo'][1]:.7f} ylo yhi\n")
        fh.write(f"{box['zlo'][0]:.7f} {box['zlo'][1]:.7f} zlo zhi\n\n")

        fh.write("Masses\n\n")
        fh.write("\n".join(sections["Masses"]) + "\n\n")

        for coeff_section in ("Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs"):
            if coeff_section in sections:
                fh.write(f"{coeff_section}\n\n")
                fh.write("\n".join(sections[coeff_section]) + "\n\n")

        fh.write("Atoms # full\n\n")
        fh.write("\n".join(" ".join(f) for f in kept_atoms) + "\n\n")

        if kept_vel:
            fh.write("Velocities\n\n")
            fh.write("\n".join(kept_vel) + "\n\n")

        for name, kept in (("Bonds", kept_bonds), ("Angles", kept_angles), ("Dihedrals", kept_dihedrals)):
            if kept:
                fh.write(f"{name}\n\n")
                fh.write("\n".join(kept) + "\n\n")

    print(
        f"fragment centroid z: {frag_cz:.3f} -> target {target_z:.3f} (dz={dz:.3f}); "
        f"removed {len(removed_molids)} hexadecane molecules "
        f"({len(removed_atom_ids)} atoms) within {args.clash_angstrom} A"
    )


if __name__ == "__main__":
    main()
