"""Deduplicate bond/angle/dihedral/atom "types" in an InterMol-written
LAMMPS data file (+ companion .input script for atom types).

InterMol's LAMMPS writer gives every individual bond/angle/dihedral its
own unique type (one row per *instance*, not per distinct parameter
set) -- confirmed live 2026-11-06: converted.lmp's 14773 "bond types"
collapse to just 8 actually-distinct (style, params) tuples, and its
28942 "angle types" to 9. This is wasteful on its own, but it also
breaks something concrete downstream: LAMMPS's restart-file
serialization for `bond_style hybrid` truncates the last several dozen
bond types' coefficients ("is not set" on read_restart) once the type
count gets this large -- a LAMMPS-side limitation, not something to work
around file-format-by-file-format. Collapsing to the true, small number
of distinct types sidesteps it entirely and is also just the normal,
intended way to use `bond_style hybrid`.

This small type count is also what makes it cheap to re-supply
bond/angle/dihedral coefficients explicitly as `bond_coeff`/
`angle_coeff`/`dihedral_coeff` commands in every one of
rules/step2_equilibrate.smk's now-separate per-stage LAMMPS rules
(scripts/merge_slab.py's `merged_coeffs.lmp` sidecar, added 2026-11-08) --
without this dedup there would be thousands of coefficient lines to carry
across every `read_data`/`include`, not 33.

InterMol does the exact same one-type-per-instance thing for ATOM types,
just split across two files: `Masses` lives in the data file, but the
organic side's `pair_coeff` (epsilon/sigma) lines live in the companion
`.input` script instead (InterMol never writes a "Pair Coeffs" section at
all -- see merge_slab.py's own docstring point 2). Confirmed live
2026-11-11: ST's converted.input has 71 raw atom types, with many
identical (mass, epsilon, sigma) triples (e.g. types 1/4/7/8/14/17/20 are
all the same (0.066, 3.5) carbon). `--pair-coeffs-input`/
`--pair-coeffs-output` turn this pass on, rewriting both files together.
"""

import argparse
import re

SECTION_RE = re.compile(
    r"^(Masses|Pair Coeffs|Bond Coeffs|Angle Coeffs|Dihedral Coeffs|"
    r"Improper Coeffs|Atoms|Velocities|Bonds|Angles|Dihedrals|Impropers)"
)
COUNT_RE = re.compile(
    r"^(\d+)\s+(atoms|bonds|angles|dihedrals|impropers|atom types|"
    r"bond types|angle types|dihedral types|improper types)\s*$"
)


def parse_sections(lines):
    header_end = 0
    for i, line in enumerate(lines):
        if SECTION_RE.match(line.split("#")[0].strip()):
            header_end = i
            break
    header = lines[:header_end]

    sections = {}
    current = None
    for line in lines[header_end:]:
        stripped = line.rstrip("\n").split("#")[0].strip()
        m = SECTION_RE.match(stripped)
        if m:
            current = line.rstrip("\n")
            sections[current] = []
            continue
        if current is None or line.strip() == "":
            continue
        sections[current].append(line.rstrip("\n"))
    return header, sections


def dedup(coeffs_key, topo_key, sections, topo_type_col=1):
    """Collapse <coeffs_key>'s rows into the distinct (style, params)
    tuples they actually represent, and remap <topo_key>'s type column
    (index topo_type_col) accordingly. Returns the new type count."""
    if coeffs_key not in sections:
        return None

    old_to_params = {}
    for line in sections[coeffs_key]:
        fields = line.split()
        old_to_params[fields[0]] = tuple(fields[1:])

    new_id_of_params = {}
    old_to_new = {}
    for old_id, params in old_to_params.items():
        if params not in new_id_of_params:
            new_id_of_params[params] = str(len(new_id_of_params) + 1)
        old_to_new[old_id] = new_id_of_params[params]

    sections[coeffs_key] = [
        f"{new_id} {' '.join(params)}"
        for params, new_id in sorted(new_id_of_params.items(), key=lambda kv: int(kv[1]))
    ]

    if topo_key in sections:
        new_topo = []
        for line in sections[topo_key]:
            fields = line.split()
            fields[topo_type_col] = old_to_new[fields[topo_type_col]]
            new_topo.append(" ".join(fields))
        sections[topo_key] = new_topo

    return len(new_id_of_params)


def read_pair_coeffs_from_input(path):
    """Parse `pair_coeff N N eps sigma` lines from an InterMol-written
    LAMMPS input script (organic side only has diagonal entries) --
    same small regex scripts/merge_slab.py's own
    read_pair_coeffs_from_input uses, duplicated here rather than
    imported (this repo's convention of small self-contained per-script
    helpers)."""
    coeffs = {}
    pat = re.compile(r"^\s*pair_coeff\s+(\d+)\s+(\1)\s+([-\d.eE+]+)\s+([-\d.eE+]+)")
    with open(path) as fh:
        for line in fh:
            m = pat.match(line)
            if m:
                coeffs[m.group(1)] = (m.group(3), m.group(4))
    return coeffs


def dedup_atom_types(sections, pair_coeffs, precision=10):
    """Collapse atom types sharing an identical (mass, epsilon, sigma)
    triple -- the same InterMol one-type-per-instance pattern as bonded
    types, just split across the data file's "Masses" section and the
    companion .input script's diagonal pair_coeff lines. Remaps
    "Atoms"'s type column (index 2: `id mol type q x y z nx ny nz`).
    Returns (new_atom_type_count, old_to_new) or None if there's no
    "Masses" section."""
    if "Masses" not in sections:
        return None

    old_to_mass = {}
    for line in sections["Masses"]:
        fields = line.split()
        old_to_mass[fields[0]] = fields[1]

    def key_of(old_id):
        eps, sigma = pair_coeffs[old_id]
        return (
            round(float(old_to_mass[old_id]), precision),
            round(float(eps), precision),
            round(float(sigma), precision),
        )

    new_id_of_key = {}
    old_to_new = {}
    for old_id in old_to_mass:
        key = key_of(old_id)
        if key not in new_id_of_key:
            new_id_of_key[key] = str(len(new_id_of_key) + 1)
        old_to_new[old_id] = new_id_of_key[key]

    mass_of_new_id = {new_id: old_to_mass[old_id] for old_id, new_id in old_to_new.items()}
    sections["Masses"] = [
        f"{new_id} {mass_of_new_id[new_id]}"
        for new_id in sorted(new_id_of_key.values(), key=int)
    ]

    if "Atoms" in sections:
        new_atoms = []
        for line in sections["Atoms"]:
            fields = line.split()
            fields[2] = old_to_new[fields[2]]
            new_atoms.append(" ".join(fields))
        sections["Atoms"] = new_atoms

    return len(new_id_of_key), old_to_new


def rewrite_pair_coeffs_input(input_path, output_path, old_to_new):
    """Rewrite a converted.input-style script's diagonal `pair_coeff N N
    eps sigma` lines to the deduplicated/remapped ids from
    dedup_atom_types, dropping now-redundant duplicates. Any other line
    (including a non-diagonal pair_coeff, not expected to occur -- see
    this script's module docstring) is copied through unchanged."""
    pat = re.compile(r"^(\s*pair_coeff\s+)(\d+)\s+(\d+)(\s+.*)$")
    seen_new_ids = set()
    out_lines = []
    with open(input_path) as fh:
        for line in fh:
            m = pat.match(line)
            if m and m.group(2) == m.group(3):
                new_id = old_to_new[m.group(2)]
                if new_id in seen_new_ids:
                    continue
                seen_new_ids.add(new_id)
                out_lines.append(f"{m.group(1)}{new_id} {new_id}{m.group(4)}\n")
            else:
                out_lines.append(line)
    with open(output_path, "w") as fh:
        fh.writelines(out_lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument(
        "--pair-coeffs-input",
        help="converted.input, for its diagonal pair_coeff lines -- "
        "turns on atom-type dedup when given",
    )
    ap.add_argument(
        "--pair-coeffs-output",
        help="path to write the deduplicated/remapped companion .input file "
        "(required with --pair-coeffs-input)",
    )
    args = ap.parse_args()
    if args.pair_coeffs_input and not args.pair_coeffs_output:
        ap.error("--pair-coeffs-output is required with --pair-coeffs-input")

    with open(args.input) as fh:
        lines = fh.readlines()

    header, sections = parse_sections(lines)

    new_counts = {}
    for coeffs_key, topo_key, count_name in (
        ("Bond Coeffs", "Bonds", "bond types"),
        ("Angle Coeffs", "Angles", "angle types"),
        ("Dihedral Coeffs", "Dihedrals", "dihedral types"),
        ("Improper Coeffs", "Impropers", "improper types"),
    ):
        n = dedup(coeffs_key, topo_key, sections)
        if n is not None:
            new_counts[count_name] = n

    if args.pair_coeffs_input:
        pair_coeffs = read_pair_coeffs_from_input(args.pair_coeffs_input)
        result = dedup_atom_types(sections, pair_coeffs)
        if result is not None:
            n_atom_types, old_to_new_atoms = result
            new_counts["atom types"] = n_atom_types
            rewrite_pair_coeffs_input(
                args.pair_coeffs_input, args.pair_coeffs_output, old_to_new_atoms
            )

    new_header = []
    for line in header:
        m = COUNT_RE.match(line.rstrip("\n"))
        if m and m.group(2) in new_counts:
            new_header.append(f"{new_counts[m.group(2)]} {m.group(2)}\n")
        else:
            new_header.append(line)

    with open(args.output, "w") as fh:
        fh.writelines(new_header)
        for section_header, body_lines in sections.items():
            fh.write(f"{section_header}\n\n")
            fh.write("\n".join(body_lines))
            fh.write("\n\n")


if __name__ == "__main__":
    main()
