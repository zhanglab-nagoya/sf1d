"""Deduplicate bond/angle/dihedral "types" in an InterMol-written LAMMPS
data file.

InterMol's LAMMPS writer gives every individual bond/angle/dihedral its
own unique type (one row per *instance*, not per distinct parameter
set) -- confirmed live 2026-11-06: converted.lmp's 14773 "bond types"
collapse to just 8 actually-distinct (style, params) tuples, and its
28942 "angle types" to 9. This is wasteful on its own, but it also
breaks something concrete downstream: LAMMPS's restart-file
serialization for `bond_style hybrid` (needed to carry per-type coeffs
between step2_equilibrate.smk's minimize and NVT stages, since
`write_data` silently drops Bond/Angle/Dihedral Coeffs sections for this
many-type hybrid setup) truncates the last several dozen bond types'
coefficients ("is not set" on read_restart) once the type count gets
this large -- a LAMMPS-side limitation, not something to work around
file-format-by-file-format. Collapsing to the true, small number of
distinct types sidesteps it entirely and is also just the normal,
intended way to use `bond_style hybrid`.
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


def dedup(coeffs_key, topo_key, sections):
    """Collapse <coeffs_key>'s rows into the distinct (style, params)
    tuples they actually represent, and remap <topo_key>'s type column
    (index 1) accordingly. Returns the new type count."""
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
            fields[1] = old_to_new[fields[1]]
            new_topo.append(" ".join(fields))
        sections[topo_key] = new_topo

    return len(new_id_of_params)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

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
