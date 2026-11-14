"""Insert a single fragment molecule into an already-equilibrated
solvent+substrate LAMMPS system, combining two INDEPENDENTLY converted
LAMMPS files rather than relocating an already-present molecule within one
file (replacing both the earlier `insert_fragment.py`, which rigid-
translated a fragment already present since dilute_pack through
already-equilibrated hexadecane and left a real, never-filled vacancy at
its origin -- confirmed live 2026-11-12 -- and the `center_fragment.py`
design that replaced it, which turned out to hide a different problem:
restraining the fragment through the whole z-barostat run made it a fixed
anchor the barostat's own proportional rescale stretched hexadecane away
from every cycle, opening a permanent cavity over a long run -- see the
"Decouple fragment insertion from solvent/substrate equilibration" plan
for the full root-cause history).

Reuses two established patterns from this repo rather than inventing new
ones:

- Type/atom/molid offsetting and a coeffs-sidecar, same idea as
  scripts/merge_slab.py's organic+substrate combination (`read_data_sections`,
  `read_pair_coeffs_from_input`, offset-and-append, no fragment<->hexadecane
  cross `pair_coeff` terms needed -- `pair_modify mix geometric` already
  handles that, same justification merge_slab.py's own docstring gives for
  organic<->substrate).
- Clash detection and whole-molecule deletion (never partial -- that would
  leave a dangling bond), same cKDTree + periodic-boxsize approach the
  deleted `insert_fragment.py` used, including its `liquid_arc_target()`
  gap-detection method to find the bulk liquid arc's midpoint as the
  insertion point (robust to substrate atoms wrapping across the box
  boundary).

Ordering constraint this script alone has to satisfy: lmp/_common_setup.
lmp.template HARDCODES `group substrate type {{ n_org_types+1 }}:
{{ n_org_types+2 }}` -- substrate must always be the LAST two atom types.
The solvent file already has substrate's types right after hexadecane's
own (merge_slab.py's own layout); inserting the fragment's atom types
means slotting them BETWEEN hexadecane's and substrate's, not after:
hexadecane keeps its own atom ids/types/molids completely unchanged,
fragment's own (already deduped) atom types get fresh ids
[n_hex_types+1 .. n_hex_types+n_frag_types], and substrate's atom TYPE
column (only -- not its atom ids or molid) shifts by +n_frag_types.
Fragment's own bond/angle/dihedral types have no such constraint
(substrate has none at all) -- plain offset-and-append after hexadecane's
own max bond/angle/dihedral type count.

Atom/bond/angle/dihedral IDs are never renumbered to stay contiguous after
deletion (confirmed from this script's own predecessor: LAMMPS's
`read_data` only needs unique IDs matching the declared count, not a
contiguous 1..N range) -- surviving solvent IDs are left exactly as they
were; only fragment's own new atoms/bonds/angles/dihedrals get fresh IDs
starting after the solvent file's own original maximum (computed before
any deletion, so it is a safe upper bound regardless of what got removed).

Run as a plain script (not a snakemake `script:` block), matching
merge_slab.py's/the deleted insert_fragment.py's own style:

    python scripts/insert_fragment.py results/slab/solvent/npt_equilibrated.data \\
        --solvent-coeffs results/slab/solvent/merged_coeffs.lmp \\
        --fragment-data results/slab/ST/fragment_converted.lmp \\
        --fragment-input results/slab/ST/fragment_converted.input \\
        --n-sub-types 2 --clash-angstrom 4.0 \\
        -o results/slab/ST/inserted.data --out-coeffs results/slab/ST/merged_coeffs.lmp
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
    """Same small parser merge_slab.py/the deleted insert_fragment.py
    already use: (header counts dict, box (xlo,xhi,ylo,yhi,zlo,zhi) dict,
    {section name: [raw lines]})."""
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


def read_pair_coeffs_from_input(path):
    """Parse `pair_coeff N N eps sigma` lines from an InterMol-written
    LAMMPS input script (organic side only has diagonal entries) --
    identical to merge_slab.py's own helper."""
    coeffs = {}
    pat = re.compile(r"^\s*pair_coeff\s+(\d+)\s+(\1)\s+([-\d.eE+]+)\s+([-\d.eE+]+)")
    with open(path) as fh:
        for line in fh:
            m = pat.match(line)
            if m:
                coeffs[int(m.group(1))] = (float(m.group(3)), float(m.group(4)))
    return coeffs


def read_pair_coeffs_from_commands(path):
    """Same as read_pair_coeffs_from_input, but for a coeffs sidecar
    file (e.g. solvent/merged_coeffs.lmp) that also has bond_coeff/
    angle_coeff/dihedral_coeff lines mixed in -- the pair_coeff regex
    alone already ignores those."""
    return read_pair_coeffs_from_input(path)


def read_bonded_coeffs_from_commands(path, command):
    """{type: 'style params...'} for bond_coeff/angle_coeff/dihedral_coeff
    lines in a coeffs sidecar file."""
    coeffs = {}
    pat = re.compile(rf"^\s*{command}\s+(\d+)\s+(.*)$")
    with open(path) as fh:
        for line in fh:
            m = pat.match(line)
            if m:
                coeffs[int(m.group(1))] = m.group(2).rstrip()
    return coeffs


def liquid_arc_target(sections, n_hex_types, zlo, zhi):
    """(target_x, target_y, target_z): the lateral box centre and the
    bulk-liquid arc's z-midpoint. Folds every substrate atom's z into the
    box-relative coordinate s = (z - zlo) mod L and finds the largest gap
    between consecutive sorted values (circularly) -- the substrate is one
    dense, lattice-packed block, so that gap is the liquid region itself,
    found this way regardless of whether the substrate happens to straddle
    the s=0/L seam. Same method the deleted insert_fragment.py used."""
    L = zhi - zlo
    s = sorted((float(line.split()[6]) - zlo) % L for line in sections["Atoms"]
               if int(line.split()[2]) > n_hex_types)
    n = len(s)
    gaps = [(s[(i + 1) % n] - s[i]) % L for i in range(n)]
    k = max(range(n), key=lambda i: gaps[i])
    arc_length = gaps[k]
    target_s = (s[k] + arc_length / 2.0) % L
    target_z = zlo + target_s
    return target_z


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("solvent_data", help="results/slab/solvent/npt_equilibrated.data")
    ap.add_argument("--solvent-coeffs", required=True, help="results/slab/solvent/merged_coeffs.lmp")
    ap.add_argument("--fragment-data", required=True, help="results/slab/{species}/fragment_converted.lmp")
    ap.add_argument("--fragment-input", required=True, help="results/slab/{species}/fragment_converted.input")
    ap.add_argument("--n-sub-types", type=int, default=2, help="substrate atom types (Fe, O)")
    ap.add_argument("--clash-angstrom", type=float, default=4.0,
                     help="hexadecane atoms within this distance of a (translated) fragment atom "
                          "get their whole molecule removed")
    ap.add_argument("-o", "--out-data", required=True)
    ap.add_argument("--out-coeffs", required=True)
    args = ap.parse_args()

    sv_header, sv_box, sv_sections = read_data_sections(args.solvent_data)
    fr_header, fr_box, fr_sections = read_data_sections(args.fragment_data)

    n_sub_types = args.n_sub_types
    n_hex_types = sv_header["atom types"] - n_sub_types
    n_frag_types = fr_header["atom types"]
    n_hex_bond_types = sv_header.get("bond types", 0)
    n_hex_angle_types = sv_header.get("angle types", 0)
    n_hex_dihedral_types = sv_header.get("dihedral types", 0)

    # --- geometry: bulk liquid arc midpoint + lateral box centre --------
    zlo, zhi = sv_box["zlo"]
    target_z = liquid_arc_target(sv_sections, n_hex_types, zlo, zhi)
    target_x = (sv_box["xlo"][0] + sv_box["xlo"][1]) / 2.0
    target_y = (sv_box["ylo"][0] + sv_box["ylo"][1]) / 2.0

    # --- translate fragment (mass-weighted centroid -> target) ----------
    fr_atoms = [line.split() for line in fr_sections["Atoms"]]
    fr_masses = {int(line.split()[0]): float(line.split()[1]) for line in fr_sections["Masses"]}
    fr_mass = np.array([fr_masses[int(f[2])] for f in fr_atoms])
    fr_xyz = np.array([[float(f[4]), float(f[5]), float(f[6])] for f in fr_atoms])
    frag_centroid = (fr_mass[:, None] * fr_xyz).sum(axis=0) / fr_mass.sum()
    shift = np.array([target_x, target_y, target_z]) - frag_centroid
    fr_xyz_new = fr_xyz + shift
    for f, xyz in zip(fr_atoms, fr_xyz_new):
        f[4], f[5], f[6] = (f"{v:.7f}" for v in xyz)

    # --- find hexadecane molecules clashing with the (now-translated)
    # fragment, drop them whole (never partial -- would leave a dangling
    # bond) --------------------------------------------------------------
    sv_atoms = [line.split() for line in sv_sections["Atoms"]]
    lx = sv_box["xlo"][1] - sv_box["xlo"][0]
    ly = sv_box["ylo"][1] - sv_box["ylo"][0]
    lz = zhi - zlo
    origin = np.array([sv_box["xlo"][0], sv_box["ylo"][0], zlo])
    boxsize = np.array([lx, ly, lz])

    hex_atoms = [f for f in sv_atoms if int(f[2]) <= n_hex_types]
    hex_ids = np.array([int(f[0]) for f in hex_atoms])
    hex_molids = np.array([int(f[1]) for f in hex_atoms])
    hex_xyz = np.array([[float(f[4]), float(f[5]), float(f[6])] for f in hex_atoms])
    hex_frac = np.mod(hex_xyz - origin, boxsize)
    frag_frac = np.mod(fr_xyz_new - origin, boxsize)

    tree = cKDTree(hex_frac, boxsize=boxsize)
    clashing = tree.query_ball_point(frag_frac, r=args.clash_angstrom)
    clashing_hex_idx = sorted({idx for group in clashing for idx in group})
    removed_molids = set(hex_molids[clashing_hex_idx])
    removed_atom_ids = {int(hex_ids[i]) for i in range(len(hex_ids)) if hex_molids[i] in removed_molids}

    kept_sv_atoms = [f for f in sv_atoms if int(f[0]) not in removed_atom_ids]
    kept_sv_vel = [
        line for line in sv_sections.get("Velocities", [])
        if int(line.split()[0]) not in removed_atom_ids
    ]

    def filter_topology(sections, section_name, n_atom_cols):
        kept = []
        for line in sections.get(section_name, []):
            fields = line.split()
            if any(int(fields[2 + k]) in removed_atom_ids for k in range(n_atom_cols)):
                continue
            kept.append(fields)
        return kept

    kept_sv_bonds = filter_topology(sv_sections, "Bonds", 2)
    kept_sv_angles = filter_topology(sv_sections, "Angles", 3)
    kept_sv_dihedrals = filter_topology(sv_sections, "Dihedrals", 4)

    # --- substrate's atom TYPE column shifts by +n_frag_types (ids/molid
    # untouched) -- hexadecane's own atoms are left completely unchanged --
    for f in kept_sv_atoms:
        t = int(f[2])
        if t > n_hex_types:
            f[2] = str(t + n_frag_types)

    # --- fragment: fresh atom/bond/angle/dihedral IDs after the solvent
    # file's own ORIGINAL maximum (computed before deletion -- a safe
    # upper bound regardless of what got removed), fresh atom types
    # [n_hex_types+1 .. n_hex_types+n_frag_types], fresh molids placed
    # after solvent's own max molid (substrate's, the largest thing in
    # the file) ------------------------------------------------------------
    max_atom_id = max(int(f[0]) for f in sv_atoms)
    max_molid = max(int(f[1]) for f in sv_atoms)
    max_bond_id = max((int(line.split()[0]) for line in sv_sections.get("Bonds", [])), default=0)
    max_angle_id = max((int(line.split()[0]) for line in sv_sections.get("Angles", [])), default=0)
    max_dihedral_id = max((int(line.split()[0]) for line in sv_sections.get("Dihedrals", [])), default=0)

    frag_molids_present = sorted({int(f[1]) for f in fr_atoms})
    molid_map = {old: max_molid + i + 1 for i, old in enumerate(frag_molids_present)}

    id_map = {}
    new_fr_atoms = []
    for i, f in enumerate(fr_atoms):
        old_id = int(f[0])
        new_id = max_atom_id + i + 1
        id_map[old_id] = new_id
        new_type = int(f[2]) + n_hex_types
        new_molid = molid_map[int(f[1])]
        # Image flags (0 0 0, fragment is fresh/never-wrapped) appended to
        # match npt_equilibrated.data's own Atoms column count -- LAMMPS
        # requires a uniform column count across the whole Atoms section,
        # and solvent/substrate atoms (wrapped by real MD) all carry them.
        new_fr_atoms.append(
            [str(new_id), str(new_molid), str(new_type), f[3], f[4], f[5], f[6], "0", "0", "0"]
        )

    def remap_fragment_topology(section_name, n_atom_cols, type_offset, id_start):
        out = []
        for i, line in enumerate(fr_sections.get(section_name, [])):
            fields = line.split()
            new_id = id_start + i + 1
            new_type = int(fields[1]) + type_offset
            atoms = [str(id_map[int(fields[2 + k])]) for k in range(n_atom_cols)]
            out.append([str(new_id), str(new_type)] + atoms)
        return out

    new_fr_bonds = remap_fragment_topology("Bonds", 2, n_hex_bond_types, max_bond_id)
    new_fr_angles = remap_fragment_topology("Angles", 3, n_hex_angle_types, max_angle_id)
    new_fr_dihedrals = remap_fragment_topology("Dihedrals", 4, n_hex_dihedral_types, max_dihedral_id)

    out_atoms = [" ".join(f) for f in kept_sv_atoms] + [" ".join(f) for f in new_fr_atoms]
    out_bonds = [" ".join(f) for f in kept_sv_bonds] + [" ".join(f) for f in new_fr_bonds]
    out_angles = [" ".join(f) for f in kept_sv_angles] + [" ".join(f) for f in new_fr_angles]
    out_dihedrals = [" ".join(f) for f in kept_sv_dihedrals] + [" ".join(f) for f in new_fr_dihedrals]

    # --- masses: order within the section doesn't matter to LAMMPS (it's
    # parsed as {type: mass} regardless of line order) -- just emit every
    # type [1 .. n_total] exactly once -------------------------------------
    out_masses = []
    for line in sv_sections["Masses"]:
        fields = line.split()
        t = int(fields[0])
        if t <= n_hex_types:
            out_masses.append(line)
        else:
            out_masses.append(f"{t + n_frag_types} {fields[1]}")
    for line in fr_sections["Masses"]:
        fields = line.split()
        out_masses.append(f"{int(fields[0]) + n_hex_types} {fields[1]}")

    # --- coeffs sidecar: hexadecane's own (unchanged) + substrate's
    # (type-shifted) pair_coeff, hexadecane's own (unchanged) bond/angle/
    # dihedral_coeff, plus fragment's own new ones (offset, appended) ------
    sv_pair_coeffs = read_pair_coeffs_from_commands(args.solvent_coeffs)
    sv_bond_coeffs = read_bonded_coeffs_from_commands(args.solvent_coeffs, "bond_coeff")
    sv_angle_coeffs = read_bonded_coeffs_from_commands(args.solvent_coeffs, "angle_coeff")
    sv_dihedral_coeffs = read_bonded_coeffs_from_commands(args.solvent_coeffs, "dihedral_coeff")

    combined_pair_coeffs = {}
    for t, (eps, sigma) in sv_pair_coeffs.items():
        new_t = t if t <= n_hex_types else t + n_frag_types
        combined_pair_coeffs[new_t] = (eps, sigma)
    fr_pair_coeffs = read_pair_coeffs_from_input(args.fragment_input)
    for t, (eps, sigma) in fr_pair_coeffs.items():
        combined_pair_coeffs[t + n_hex_types] = (eps, sigma)

    fr_bond_coeffs = {
        int(line.split()[0]) + n_hex_bond_types: " ".join(line.split()[1:])
        for line in fr_sections.get("Bond Coeffs", [])
    }
    fr_angle_coeffs = {
        int(line.split()[0]) + n_hex_angle_types: " ".join(line.split()[1:])
        for line in fr_sections.get("Angle Coeffs", [])
    }
    fr_dihedral_coeffs = {
        int(line.split()[0]) + n_hex_dihedral_types: " ".join(line.split()[1:])
        for line in fr_sections.get("Dihedral Coeffs", [])
    }

    with open(args.out_coeffs, "w") as fh:
        for t in sorted(combined_pair_coeffs):
            eps, sigma = combined_pair_coeffs[t]
            fh.write(f"pair_coeff {t} {t} {eps:.7f} {sigma:.7f}\n")
        for command, sv_coeffs, fr_coeffs in (
            ("bond_coeff", sv_bond_coeffs, fr_bond_coeffs),
            ("angle_coeff", sv_angle_coeffs, fr_angle_coeffs),
            ("dihedral_coeff", sv_dihedral_coeffs, fr_dihedral_coeffs),
        ):
            for t in sorted(sv_coeffs):
                fh.write(f"{command} {t} {sv_coeffs[t]}\n")
            for t in sorted(fr_coeffs):
                fh.write(f"{command} {t} {fr_coeffs[t]}\n")

    # --- write combined data file -----------------------------------------
    n_atom_types = n_hex_types + n_frag_types + n_sub_types
    n_bond_types = n_hex_bond_types + fr_header.get("bond types", 0)
    n_angle_types = n_hex_angle_types + fr_header.get("angle types", 0)
    n_dihedral_types = n_hex_dihedral_types + fr_header.get("dihedral types", 0)

    with open(args.out_data, "w") as fh:
        fh.write("Fragment inserted into equilibrated solvent+substrate (insert_fragment.py)\n\n")
        fh.write(f"{len(out_atoms)} atoms\n")
        fh.write(f"{len(out_bonds)} bonds\n")
        fh.write(f"{len(out_angles)} angles\n")
        fh.write(f"{len(out_dihedrals)} dihedrals\n")
        fh.write(f"{sv_header.get('impropers', 0)} impropers\n\n")
        fh.write(f"{n_atom_types} atom types\n")
        if n_bond_types:
            fh.write(f"{n_bond_types} bond types\n")
        if n_angle_types:
            fh.write(f"{n_angle_types} angle types\n")
        if n_dihedral_types:
            fh.write(f"{n_dihedral_types} dihedral types\n")
        fh.write("\n")
        fh.write(f"{sv_box['xlo'][0]:.7f} {sv_box['xlo'][1]:.7f} xlo xhi\n")
        fh.write(f"{sv_box['ylo'][0]:.7f} {sv_box['ylo'][1]:.7f} ylo yhi\n")
        fh.write(f"{sv_box['zlo'][0]:.7f} {sv_box['zlo'][1]:.7f} zlo zhi\n\n")

        fh.write("Masses\n\n")
        fh.write("\n".join(out_masses) + "\n\n")

        fh.write("Atoms # full\n\n")
        fh.write("\n".join(out_atoms) + "\n\n")

        if kept_sv_vel:
            # LAMMPS' Velocities section, if present, must give exactly one
            # line per atom in the file -- a partial section (hexadecane +
            # substrate only) makes it misread the next section header as a
            # velocity line ("Incorrect format in Velocities section...:
            # Bonds", confirmed live 2026-11-14). Fragment atoms get plain
            # zero placeholders: lmp_insert_minimize's `minimize` never
            # reads velocities, and its trailing `velocity fragment create`
            # overwrites these before any dynamics (lmp/minimize.lmp.template)
            # ever uses them -- so this placeholder is never actually acted on.
            frag_vel = [f"{f[0]} 0 0 0" for f in new_fr_atoms]
            fh.write("Velocities\n\n")
            fh.write("\n".join(kept_sv_vel + frag_vel) + "\n\n")

        for name, kept in (("Bonds", out_bonds), ("Angles", out_angles), ("Dihedrals", out_dihedrals)):
            if kept:
                fh.write(f"{name}\n\n")
                fh.write("\n".join(kept) + "\n\n")

    print(
        f"fragment centroid -> target ({target_x:.3f}, {target_y:.3f}, {target_z:.3f}); "
        f"removed {len(removed_molids)} hexadecane molecules ({len(removed_atom_ids)} atoms) "
        f"within {args.clash_angstrom} A; wrote {len(out_atoms)} atoms total "
        f"({len(fr_atoms)} fragment, {n_frag_types} new atom types)"
    )


if __name__ == "__main__":
    main()
