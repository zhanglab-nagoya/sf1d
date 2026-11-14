"""Deterministically close the vacuum-like gap that forms in the mobile
(hexadecane) liquid during NVT at fixed volume -- the same failure mode
../2026-10-28-generate-aa-pama's own `shrink_substrate_box.py` fixes
for its free-surface geometry ("the substrate.data box typically includes
vacuum padding above the slab... this makes the barostat place the polymer
too high"), just adapted to our fully 3D-periodic sandwich (no free surface,
no literal padding -- confirmed live 2026-11-12 on ST's own nvt.data: the
liquid's own thermal tail has thinned out and partly diffused clean through
the periodic boundary, landing stray atoms just below the substrate's own
z_min, 33 A away from where the bulk liquid actually sits).

Run between lmp_nvt and lmp_gap_minimize/lmp_settle
(rules/step2_equilibrate.smk) -- same pipeline position as the sibling
repo's own shrink_box_zhi, between its NVT (333) and NPT (335) stages.
This script is this pipeline's analog of ../2026-07-28-fe2o3-cg's own
`fix deform` step (in/530_.../compress.lmp.template): both size the box
deterministically from the solvent's own known reference bulk density,
no pressure feedback of any kind. A barostat was tried here instead
first (both an isotropic-average and a Pzz-only, liquid-only-pressure
hand-rolled Berendsen rescale) and confirmed live 2026-11-12/13 to never
reliably converge on this frozen-substrate system -- the box just kept
inflating without the pressure signal ever settling near target, and
even the corrected liquid-only version's convergence rate turned out to
be set by the liquid's own real compressibility (an effective ~120 ps
relaxation time) regardless of tuning, matching fe2o3-cg's own
independent conclusion that a barostat is a structural dead end here,
not a tuning problem. There is no barostat stage left in this pipeline
at all -- lmp_settle (plain, fixed-volume NVT) is the only step after
this one and lmp_gap_minimize.

Unlike shrink_substrate_box.py's literal "zhi = actual atom z_max" (valid
only because nothing at all exists past the real material in a free-surface
geometry, a hard vacuum edge), there is no sharp edge here to detect at
all -- confirmed live 2026-11-12 that the mobile population's own density
profile merely TAPERS continuously (a liquid film's natural leading-edge
falloff, not a step), so no density-threshold/profile-shape heuristic can
robustly tell "real tail" from "actual vacuum": tightening the threshold
made the residual gap WORSE, not better, and iterating the whole detection
on its own output never converges -- it just keeps finding "more gap" and
crushing the whole system, cycle after cycle.

So this sizes the box from MASS CONSERVATION instead of from the density
profile's shape at all: the mobile population's own total mass (every
mobile atom's own mass, directly from the Masses section -- not a
molar-mass-times-molecule-count estimate) divided by hexadecane's own
measured reference density (results/segment_volume.yaml's
hexadecane_density_kg_per_m3, the same value step2_dilute.smk's solvent_pack
already sizes the ORIGINAL packing box from) gives the volume this
material should occupy; dividing by the box's own (untouched) lateral
footprint gives the target arc length directly, with no profile-shape
ambiguity at all. Runs on the solvent-only system since 2026-11-13's
solvent/fragment split (rules/step2_dilute.smk's module docstring) -- no
fragment is present at this stage at all, so there is no correction term
to worry about here anymore.

Every mobile molecule is then translated as a RIGID BODY (never an
individual-atom rescale, which would mechanically compress whatever
fraction of a bond happens to point along z -- confirmed live this
actually matters: an earlier draft that rescaled atoms directly measurably
heated the system within the first 500 steps, then lost a bond once the
barostat's own first change_box compounded it) by an amount proportional
to its own distance from the substrate's own top edge -- the same affine
operation `_liquid_barostat.lmp.template`'s own `change_box ... remap`
already applies every cycle, just computed once from geometry instead of
many times from a slow/non-converging pressure signal. The substrate is
never touched.

Two independently-shifted, long/floppy hexadecane chains can still land
with overlapping conformations near each other even with zero intramolecular
distortion (confirmed live, down to 0.23 A between different molecules in
one real case) -- that's expected and is why this is always followed by
lmp_gap_minimize (rules/step2_equilibrate.smk), the same established
"move, then minimize whatever local clash remains" pattern the old
scripts/insert_fragment.py (now retired) + lmp_insert_minimize used.

Run as a plain script (matching merge_slab.py/center_fragment.py's own
style):

    python scripts/close_liquid_gap.py nvt.data -o nvt_shrunk.data \\
        --n-sub-types 2 --segment-volume results/segment_volume.yaml
"""

import argparse
import re

import yaml

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
    """Same small parser merge_slab.py/insert_fragment.py/center_fragment.py
    already use: (header counts dict, box dict, {section name: [raw lines]})."""
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
            key = m.group(3).split()[0]
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


def unwrap_z(z, zlo, lz, sub_zmax):
    """Fold a mobile atom's z into a single continuous coordinate starting
    at the substrate's own top (sub_zmax) -- an atom that has diffused
    across the periodic boundary and landed just below the substrate
    (z < sub_zmax, confirmed live on ST's own nvt.data: a handful of
    hexadecane atoms at z<3 while the substrate itself starts at z=3) is
    really just the far tail of the SAME liquid arc, so it gets +lz to
    place it logically right after the rest of that arc instead of before
    the substrate."""
    return z if z >= sub_zmax else z + lz


def arc_position(z, zlo, lz, sub_zmax):
    return unwrap_z(z, zlo, lz, sub_zmax) - sub_zmax


N_A = 6.02214076e23  # 1/mol -- same constant step2_dilute.smk's solvent_pack already uses


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("data", help="post-NVT LAMMPS data file (nvt.data)")
    ap.add_argument("-o", "--out-data", required=True)
    ap.add_argument("--n-sub-types", type=int, default=2, help="substrate atom types (Fe, O)")
    ap.add_argument("--segment-volume", required=True,
                     help="results/segment_volume.yaml, for hexadecane_density_kg_per_m3 -- "
                          "the same reference density step2_dilute.smk's solvent_pack already "
                          "sizes the original packing box from")
    args = ap.parse_args()

    header, box, sections = read_data_sections(args.data)
    n_atom_types = header["atom types"]
    n_org_types = n_atom_types - args.n_sub_types
    masses = {int(line.split()[0]): float(line.split()[1]) for line in sections["Masses"]}

    zlo, zhi = box["zlo"]
    lz = zhi - zlo
    lx = box["xlo"][1] - box["xlo"][0]
    ly = box["ylo"][1] - box["ylo"][0]

    atoms = [line.split() for line in sections["Atoms"]]
    mobile = [f for f in atoms if int(f[2]) <= n_org_types]
    substrate = [f for f in atoms if int(f[2]) > n_org_types]

    sub_z = [float(f[6]) for f in substrate]
    sub_zmin, sub_zmax = min(sub_z), max(sub_z)
    # The liquid's own occupied arc: box length minus the substrate's own
    # (never-moving) footprint -- NOT the same as lz, so a rescale factor
    # must be computed from THIS length, not lz, or it under/over-corrects
    # by however much of the box the substrate itself occupies (~30 A out
    # of ~106 A here -- confirmed live this was a real bug in an earlier
    # draft: rescaling by new_lz/lz instead of new_arc/arc dragged the
    # dense bulk right next to the substrate backward INTO it).
    arc_length = lz - (sub_zmax - sub_zmin)

    with open(args.segment_volume) as fh:
        seg = yaml.safe_load(fh)
    target_density_kg_per_m3 = seg["hexadecane_density_kg_per_m3"]

    # Sum each mobile atom's own mass (g/mol) -- not molar-mass times a
    # molecule count. Simpler than a per-molecule-type count, and exact
    # either way since `mobile` is pure hexadecane here (no fragment
    # present at this stage at all, see this file's module docstring).
    total_mass_g_per_mol = sum(masses[int(f[2])] for f in mobile)
    total_mass_kg = (total_mass_g_per_mol / N_A) / 1000.0

    # Angstrom^3 -> m^3 is 1e-30; volume_m3 = mass_kg / density_kg_per_m3.
    target_volume_angstrom3 = (total_mass_kg / target_density_kg_per_m3) / 1e-30
    new_arc_length = target_volume_angstrom3 / (lx * ly)

    # new_zhi = sub_zmax + new_arc_length, NOT zlo + (lz - (arc_length -
    # new_arc_length)) -- confirmed live 2026-11-12 that second form is
    # wrong whenever the kept arc (as it always has so far) sits entirely
    # within the un-wrapped part of the arc (before the wrap to zlo):
    # subtracting "arc_length - new_arc_length" from zhi silently credits
    # back the small, UNRELATED, always-present wrapped stub near zlo
    # (sub_zmin - zlo, part of arc_length's own total but nothing to do
    # with the far-end gap being closed here) as if it were available
    # headroom, under-shrinking zhi by exactly that amount.
    unwrapped_extent = zhi - sub_zmax
    if new_arc_length > unwrapped_extent:
        raise ValueError(
            f"target arc length {new_arc_length:.3f} extends past the un-wrapped "
            f"portion ({unwrapped_extent:.3f}) -- the gap-closing math above only "
            "handles a target entirely before the wrap to zlo, which has held for "
            "every real case seen so far; this run needs the formula generalized."
        )
    new_zhi = sub_zmax + new_arc_length
    new_lz = new_zhi - zlo
    # mu scales the LIQUID's own arc length only (anchored at sub_zmax, the
    # substrate's own top -- so an atom sitting right there, zero distance
    # into the arc, never moves) -- not the whole box length, which would
    # incorrectly also drag the dense bulk right next to the substrate
    # backward by however much the (untouched) substrate itself occupies.
    mu = new_arc_length / arc_length if arc_length > 0 else 1.0

    molecule_atoms = {}
    for f in mobile:
        molecule_atoms.setdefault(f[1], []).append(f)

    for mol_atoms in molecule_atoms.values():
        # Translate each molecule as a RIGID BODY, by its own centroid's
        # required shift -- never rescale individual atoms, which would
        # mechanically compress whatever fraction of each bond happens to
        # point along z. Confirmed live 2026-11-12 this matters: an
        # earlier draft scaled every atom's own distance into the arc by
        # mu directly; that shortens any bond with a z-component by the
        # same factor, injecting real strain -- confirmed by the elevated
        # c_mobtemp (642 K vs. 300 K target) within the first 500 steps,
        # then a lost bond once the barostat's own first change_box
        # compounded it. A pure per-molecule translation, by contrast,
        # preserves every molecule's internal geometry exactly (zero
        # distortion) -- the same safe pattern merge_slab.py/
        # center_fragment.py already use for their own whole-molecule
        # shifts, just computed from each molecule's own centroid here
        # instead of a single global offset.
        atom_unwrapped = [unwrap_z(float(f[6]), zlo, lz, sub_zmax) for f in mol_atoms]
        old_centroid_unwrapped = sum(atom_unwrapped) / len(atom_unwrapped)
        new_centroid_unwrapped = sub_zmax + (old_centroid_unwrapped - sub_zmax) * mu
        delta = new_centroid_unwrapped - old_centroid_unwrapped
        # Shift every atom's own UNWRAPPED position by the same delta
        # (preserving the molecule's internal geometry exactly), then
        # convert back to a valid wrapped coordinate in the NEW, shrunk
        # box -- using the OLD lz to unwrap above and the NEW lz to
        # re-wrap here is required, not a shortcut-able no-op: a molecule
        # that started past sub_zmin (wrapped under the OLD box) or ends
        # up past the new zhi (wrapped under the NEW box) needs a
        # DIFFERENT box length on each side of that conversion. The
        # wrap decision is made ONCE per molecule, from its own new
        # centroid, and applied identically to every atom in it -- same
        # periodic-image convention merge_slab.py/center_fragment.py
        # already use elsewhere (never per-atom, which would split bonds).
        wrap = new_centroid_unwrapped >= new_zhi
        for f, z_unwrapped in zip(mol_atoms, atom_unwrapped):
            new_z_unwrapped = z_unwrapped + delta
            final_z = new_z_unwrapped - new_lz if wrap else new_z_unwrapped
            f[6] = f"{final_z:.7f}"

    sections["Atoms"] = [" ".join(f) for f in (mobile + substrate)]
    with open(args.out_data, "w") as fh:
        fh.write("Liquid vacuum gap closed (scripts/close_liquid_gap.py)\n\n")
        fh.write(f"{header['atoms']} atoms\n")
        fh.write(f"{header['bonds']} bonds\n")
        fh.write(f"{header['angles']} angles\n")
        fh.write(f"{header['dihedrals']} dihedrals\n")
        fh.write(f"{header.get('impropers', 0)} impropers\n\n")
        fh.write(f"{n_atom_types} atom types\n")
        for key in ("bond types", "angle types", "dihedral types"):
            if key in header:
                fh.write(f"{header[key]} {key}\n")
        fh.write("\n")
        fh.write(f"{box['xlo'][0]:.7f} {box['xlo'][1]:.7f} xlo xhi\n")
        fh.write(f"{box['ylo'][0]:.7f} {box['ylo'][1]:.7f} ylo yhi\n")
        fh.write(f"{zlo:.7f} {new_zhi:.7f} zlo zhi\n\n")

        for section_name in ("Masses", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs"):
            if section_name in sections:
                fh.write(f"{section_name}\n\n")
                fh.write("\n".join(sections[section_name]) + "\n\n")

        fh.write("Atoms # full\n\n")
        fh.write("\n".join(sections["Atoms"]) + "\n\n")

        for section_name in ("Velocities", "Bonds", "Angles", "Dihedrals", "Impropers"):
            if section_name in sections:
                fh.write(f"{section_name}\n\n")
                fh.write("\n".join(sections[section_name]) + "\n\n")

    current_density = total_mass_kg / ((lx * ly * arc_length) * 1e-30)
    print(
        f"mobile population currently at {current_density:.1f} kg/m^3 over arc "
        f"length {arc_length:.3f}; target {target_density_kg_per_m3:.1f} kg/m^3 "
        f"-> target arc length {new_arc_length:.3f}; "
        f"shrinking Lz {lz:.3f} -> {new_lz:.3f} (mu={mu:.6f}); zhi {zhi:.3f} -> {new_zhi:.3f}"
    )


if __name__ == "__main__":
    main()
