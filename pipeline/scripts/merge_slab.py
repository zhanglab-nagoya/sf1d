"""Merge a dilute organic (fragment+hexadecane) LAMMPS box (converted via
InterMol, see gmx_to_lammps.py) onto the Fe2O3 substrate
(results/substrate/substrate_minimized.data) into one periodic LAMMPS
data file, for step 2's slab-adsorption production.

Design choices, both deliberate:

1. Fully 3D-periodic sandwich, not a vacuum-padded free liquid surface.
   The organic box is already laterally commensurate with the substrate's
   own relaxed footprint (step2_dilute.smk's dilute_pack reads it
   straight from substrate_minimized.data), so x,y need no further work
   here. In z, rather than adding a large vacuum region above the liquid
   (which LAMMPS's kspace_modify slab correction would need, typically
   >=2-3x the filled region's own extent -- hundreds of Angstrom here,
   wasteful for this deliberately small system), this stacks substrate
   and liquid so the box wraps around periodically: liquid's top face
   meets the substrate's own bottom face again via the z periodic
   boundary. This is the same "periodic wall" mental model sf1d's own
   SCF model already uses (one wall with a mirror boundary, or a slit
   between two walls) and the same plain 3D PME/PPPM treatment GROMACS
   has used everywhere else in this pipeline -- no slab correction, no
   extra vacuum bookkeeping.

2. Pair coefficients are NOT written into the merged data file. Neither
   the organic side (converted.lmp has no "Pair Coeffs" section --
   InterMol writes `pair_coeff` commands straight into converted.input
   instead) nor this merge embeds them; both the organic diagonal
   pair_coeff values (parsed from converted.input) and the substrate's
   (parsed from substrate_minimized.data's own "Pair Coeffs" section,
   offset by the organic atom-type count) are written to a small sidecar
   file instead (merged_pair_coeffs.lmp), for the next stage's run
   scripts to `include` alongside a single, consistent
   `pair_style`/`kspace_style`/`pair_modify mix geometric` -- geometric
   mixing is already the combination rule BOTH force fields use (OPLS-AA
   for the organic side, confirmed by converted.input's own
   `pair_modify mix geometric` line; INTERFACE/OPLS-CVFF for the
   substrate, see forcefield/substrate/README.md), so no explicit
   organic<->substrate cross terms are needed at all.

Atom IDs, types, and molecule IDs are offset/renumbered so the two sides
never collide: organic keeps its own IDs/types unchanged (1..n_org_atoms,
1..n_org_types), substrate atoms are renumbered contiguously starting at
n_org_atoms+1 and its 2 types offset to n_org_types+1/+2; substrate's own
molecule-id (0 for every atom) is replaced by a single new id
(max organic molid + 1) so it doesn't collide with a real organic
molecule id. Bonds/Angles/Dihedrals (organic-only -- the substrate is a
purely nonbonded ionic lattice, no bonds/angles/dihedrals sections at
all) are copied verbatim, since organic atom IDs/types are unchanged.
"""

import argparse
import re


def read_data_sections(path):
    """Parse a LAMMPS data file into (header, box, sections) where header
    is a dict of the count lines (e.g. header["atoms"] = 15074), box is
    (xlo, xhi, ylo, yhi, zlo, zhi), and sections maps each section name
    (e.g. "Atoms", "Masses") to its list of raw data lines (no header/
    blank lines)."""
    with open(path) as fh:
        lines = fh.readlines()

    header = {}
    box = {}
    count_re = re.compile(r"^\s*(\d+)\s+(atoms|bonds|angles|dihedrals|impropers|"
                           r"atom types|bond types|angle types|dihedral types|"
                           r"improper types)\s*$")
    box_re = re.compile(r"^\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+(xlo xhi|ylo yhi|zlo zhi)\s*$")

    i = 0
    n = len(lines)
    # Header block: everything before the first section keyword.
    section_names = [
        "Masses", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs",
        "Dihedral Coeffs", "Improper Coeffs", "Atoms", "Velocities",
        "Bonds", "Angles", "Dihedrals", "Impropers",
    ]
    while i < n:
        line = lines[i].rstrip("\n")
        stripped = line.split("#")[0].strip()
        if stripped in section_names:
            break
        m = count_re.match(line)
        if m:
            header[m.group(2)] = int(m.group(1))
        m = box_re.match(line)
        if m:
            key = m.group(3).split()[0]  # "xlo", "ylo", "zlo"
            box[key] = (float(m.group(1)), float(m.group(2)))
        i += 1

    sections = {}
    current = None
    for line in lines[i:]:
        stripped = line.rstrip("\n").split("#")[0].strip()
        if stripped in section_names:
            current = stripped
            sections[current] = []
            continue
        if current is None:
            continue
        if line.strip() == "":
            continue
        sections[current].append(line.rstrip("\n"))

    return header, box, sections


def read_pair_coeffs_from_input(path):
    """Parse `pair_coeff N N eps sigma` lines from an InterMol-written
    LAMMPS input script (organic side only has diagonal entries)."""
    coeffs = {}
    pat = re.compile(r"^\s*pair_coeff\s+(\d+)\s+(\1)\s+([-\d.eE+]+)\s+([-\d.eE+]+)")
    with open(path) as fh:
        for line in fh:
            m = pat.match(line)
            if m:
                coeffs[int(m.group(1))] = (float(m.group(3)), float(m.group(4)))
    return coeffs


def read_pair_coeffs_from_data(sections):
    coeffs = {}
    for line in sections.get("Pair Coeffs", []):
        fields = line.split()
        coeffs[int(fields[0])] = (float(fields[1]), float(fields[2]))
    return coeffs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organic-data", required=True)
    ap.add_argument("--organic-input", required=True, help="converted.input, for its pair_coeff lines")
    ap.add_argument("--substrate-data", required=True)
    ap.add_argument("--gap", type=float, default=3.0, help="buffer gap in Angstrom at each substrate/liquid/periodic-wrap interface")
    ap.add_argument("-o", "--out-data", required=True)
    ap.add_argument("--out-pair-coeffs", required=True)
    args = ap.parse_args()

    org_header, org_box, org_sections = read_data_sections(args.organic_data)
    sub_header, sub_box, sub_sections = read_data_sections(args.substrate_data)

    n_org_atoms = org_header["atoms"]
    n_org_types = org_header["atom types"]
    n_sub_types = sub_header["atom types"]

    org_pair_coeffs = read_pair_coeffs_from_input(args.organic_input)
    sub_pair_coeffs = read_pair_coeffs_from_data(sub_sections)
    assert len(org_pair_coeffs) == n_org_types
    assert len(sub_pair_coeffs) == n_sub_types

    # --- geometry -----------------------------------------------------
    org_zlo = org_box["zlo"][0]
    org_lz = org_box["zlo"][1] - org_box["zlo"][0]
    org_xlo = org_box["xlo"][0]
    org_ylo = org_box["ylo"][0]
    sub_xlo = sub_box["xlo"][0]
    sub_ylo = sub_box["ylo"][0]
    lx = org_box["xlo"][1] - org_box["xlo"][0]
    ly = org_box["ylo"][1] - org_box["ylo"][0]

    sub_atom_lines = sub_sections["Atoms"]
    sub_zs = [float(l.split()[6]) for l in sub_atom_lines]
    sub_zmin, sub_zmax = min(sub_zs), max(sub_zs)
    sub_thickness = sub_zmax - sub_zmin

    g = args.gap
    sub_shift_z = g - sub_zmin

    org_shift_x = sub_xlo - org_xlo
    org_shift_y = sub_ylo - org_ylo

    # --- atoms ----------------------------------------------------------
    org_molids = [int(l.split()[1]) for l in org_sections["Atoms"]]
    max_org_molid = max(org_molids)
    sub_molid = max_org_molid + 1

    # Substrate atoms keep their own x,y (the reference frame for the
    # merged box's lateral origin) and shift only in z.
    new_sub_atoms = []
    for new_id, line in enumerate(sub_atom_lines, start=n_org_atoms + 1):
        fields = line.split()
        old_type = int(fields[2])
        q = float(fields[3])
        x, y, z = float(fields[4]), float(fields[5]), float(fields[6])
        new_type = old_type + n_org_types
        new_sub_atoms.append(
            f"{new_id} {sub_molid} {new_type} {q:.6f} "
            f"{x:.7f} {y:.7f} {(z + sub_shift_z):.7f}"
        )

    # gmx_trjconv_whole (-pbc whole, see step2_convert.smk) deliberately
    # writes some molecules with coordinates outside the organic box's
    # own [zlo, zhi) so no molecule is split across that box's periodic
    # boundary. The merged box has a DIFFERENT z length, so first undo
    # that per-molecule (never per-atom, which would re-split bonds):
    # shift each whole molecule by an integer multiple of the ORIGINAL
    # org_lz so its centroid lands back in [org_zlo, org_zlo+org_lz).
    from collections import defaultdict

    molecule_atoms = defaultdict(list)
    for line in org_sections["Atoms"]:
        fields = line.split()
        molecule_atoms[fields[1]].append(fields)

    for molid, atoms in molecule_atoms.items():
        centroid_z = sum(float(f[6]) for f in atoms) / len(atoms)
        k = (centroid_z - org_zlo) // org_lz
        if k != 0:
            for f in atoms:
                f[6] = f"{float(f[6]) - k * org_lz:.7f}"

    all_org_atoms = [f for atoms in molecule_atoms.values() for f in atoms]

    # Individual atoms can still dangle beyond their own molecule's
    # wrapped centroid (a floppy hexadecane chain's far end, say) -- use
    # the ACTUAL occupied z-range after wrapping, not the nominal box
    # length, so no atom can end up placed below the substrate's top
    # surface (which would be a real steric clash against a rigid
    # crystal, unlike dipping past a periodic liquid-liquid boundary,
    # which is harmless).
    org_zs_wrapped = [float(f[6]) for f in all_org_atoms]
    org_occupied_min = min(org_zs_wrapped)
    org_occupied_thickness = max(org_zs_wrapped) - org_occupied_min

    liquid_start = g + sub_thickness + g
    org_shift_z = liquid_start - org_occupied_min
    box_zhi = liquid_start + org_occupied_thickness + g

    # Shift organic atoms (x,y to align with the substrate's own lateral
    # origin, z to sit above the substrate -- see module docstring).
    shifted_org_atoms = []
    for fields in all_org_atoms:
        atom_id, molid, atype, q = fields[0], fields[1], fields[2], fields[3]
        x, y, z = float(fields[4]), float(fields[5]), float(fields[6])
        shifted_org_atoms.append(
            f"{atom_id} {molid} {atype} {q} "
            f"{(x + org_shift_x):.7f} {(y + org_shift_y):.7f} {(z + org_shift_z):.7f}"
        )
    out_atoms = shifted_org_atoms + new_sub_atoms

    # --- velocities -------------------------------------------------
    out_vel = list(org_sections.get("Velocities", []))
    for new_id, line in enumerate(sub_sections.get("Velocities", []), start=n_org_atoms + 1):
        fields = line.split()
        out_vel.append(f"{new_id} {fields[1]} {fields[2]} {fields[3]}")

    # --- masses -------------------------------------------------------
    out_masses = list(org_sections["Masses"])
    for line in sub_sections["Masses"]:
        fields = line.split()
        new_type = int(fields[0]) + n_org_types
        out_masses.append(f"{new_type} {fields[1]}")

    # --- combined pair coeffs sidecar --------------------------------
    combined_pair_coeffs = dict(org_pair_coeffs)
    for old_type, (eps, sigma) in sub_pair_coeffs.items():
        combined_pair_coeffs[old_type + n_org_types] = (eps, sigma)

    with open(args.out_pair_coeffs, "w") as fh:
        for t in sorted(combined_pair_coeffs):
            eps, sigma = combined_pair_coeffs[t]
            fh.write(f"pair_coeff {t} {t} {eps:.7f} {sigma:.7f}\n")

    # --- write merged data file --------------------------------------
    n_atoms = n_org_atoms + sub_header["atoms"]
    n_atom_types = n_org_types + n_sub_types

    with open(args.out_data, "w") as fh:
        fh.write("Merged organic dilute slab + Fe2O3 substrate (merge_slab.py)\n\n")
        fh.write(f"{n_atoms} atoms\n")
        fh.write(f"{org_header.get('bonds', 0)} bonds\n")
        fh.write(f"{org_header.get('angles', 0)} angles\n")
        fh.write(f"{org_header.get('dihedrals', 0)} dihedrals\n")
        fh.write(f"{org_header.get('impropers', 0)} impropers\n\n")
        fh.write(f"{n_atom_types} atom types\n")
        if "bond types" in org_header:
            fh.write(f"{org_header['bond types']} bond types\n")
        if "angle types" in org_header:
            fh.write(f"{org_header['angle types']} angle types\n")
        if "dihedral types" in org_header:
            fh.write(f"{org_header['dihedral types']} dihedral types\n")
        fh.write("\n")
        fh.write(f"{sub_xlo:.7f} {sub_xlo + lx:.7f} xlo xhi\n")
        fh.write(f"{sub_ylo:.7f} {sub_ylo + ly:.7f} ylo yhi\n")
        fh.write(f"0.0 {box_zhi:.7f} zlo zhi\n\n")

        fh.write("Masses\n\n")
        fh.write("\n".join(out_masses) + "\n\n")

        for coeff_section in ("Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs"):
            if coeff_section in org_sections:
                fh.write(f"{coeff_section}\n\n")
                fh.write("\n".join(org_sections[coeff_section]) + "\n\n")

        fh.write("Atoms # full\n\n")
        fh.write("\n".join(out_atoms) + "\n\n")

        if out_vel:
            fh.write("Velocities\n\n")
            fh.write("\n".join(out_vel) + "\n\n")

        for topo_section in ("Bonds", "Angles", "Dihedrals", "Impropers"):
            if topo_section in org_sections:
                fh.write(f"{topo_section}\n\n")
                fh.write("\n".join(org_sections[topo_section]) + "\n\n")


if __name__ == "__main__":
    main()
