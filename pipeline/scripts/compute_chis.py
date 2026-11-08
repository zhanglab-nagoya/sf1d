"""Step 2 chis stage 6: chis_A (ST), chis_B (MMA) from stage 5's z-resolved
density profiles (results/slab/{species}/{fragment,hexadecane}_density.profile),
per forcefield/chi_chis_from_md.md step 2:

    chis_x = ln[phi_contact,x / phi_bulk,x] - ln[phi_contact,S / phi_bulk,S]

x = the adsorbing fragment (ST-analog -> A, MMA-analog -> B), S = hexadecane
(the solvent). Because this is a ratio of the SAME species' density between
two regions, the volume-fraction-vs-number-density distinction cancels (both
regions share the same 1 Angstrom bin width and the same species' per-atom
volume), so the raw `density/number` column from `fix ave/chunk` can be used
directly as phi -- no mass/volume conversion needed.

Geometry: merged.data places the substrate as a single contiguous block
[sub_zmin, sub_zmax] (confirmed live 2026-11-08 by directly reading its
Atoms section, for both species -- NOT two separate periodic-wrapped
chunks as the box's zero/nonzero hexadecane-density regions alone might
suggest). Everything outside that range is a single contiguous "liquid
arc" running from sub_zmax up through the box's periodic wrap back down to
sub_zmin, parameterized by arc-length coordinate s = (z - sub_zmax) % Lz,
s in (0, arc_length). s=0 is the substrate's "upper" face; s=arc_length is
its "lower" face (approached by the liquid via the periodic wrap). This
generalizes cleanly over an asymmetry confirmed live 2026-11-08: ST's 10 ns
production run left a genuine ~10 A vacuum gap at the wrap-around (lower)
face -- the liquid never fully closed that side -- while MMA's did close
it and sees contact on both faces. Each face's validity (real liquid
contact vs. a vacuum/free-surface edge) is checked independently per
species against config.yaml's slab.chis_min_contact_density_fraction,
using hexadecane's own mean density there (hundreds of molecules, well
sampled, unlike the single fragment).

Because the fragment group is a single molecule, its per-window contact/
bulk density is very often exactly zero (simply absent from that region
in a given 50 ps window) -- ln(0) is undefined, so chis_x is NOT computed
per window. Instead this builds a cumulative (running) time-average over
the retained, time-ordered windows and evaluates chis_x on that running
average at every point -- the final value (full-average) is the reported
chis_x, and comparing the running estimate at the halfway point vs. the
end is the only convergence check available from a single run (see
config.yaml's slab.chis_convergence_tol and production_ns comments on why
this is a heuristic, not a guarantee, for a single-fragment system).

results/chis.yaml is written as a bare {"A": ..., "B": ...} dict (nothing
else at the top level) because scripts/sf1d_predict.py loads it and passes
it straight through as SF1D(..., chis=chis, ...), which iterates
chis.items() expecting exactly that shape (src/sf1d/core.py). All the
diagnostics (contact/bulk densities, valid faces, convergence check, the
running series) go to the separate results/chis_diagnostics.yaml instead.

Run as a snakemake `script:` -- `snakemake` is injected (input/output/params).
"""

import numpy as np
import yaml

SUBSTRATE_LABEL = {"ST": "A", "MMA": "B"}


def read_data_header_and_substrate_z(data_path, n_sub_types):
    """box (xlo,xhi,ylo,yhi,zlo,zhi) and the substrate's own [zmin,zmax],
    reading merged.data directly (atom_style full: id mol type q x y z).
    The substrate always has the LAST n_sub_types atom types (merge_slab.py
    offsets the organic side's types to 1..n_org_types and appends the
    substrate's at n_org_types+1..n_org_types+n_sub_types)."""
    box = {}
    n_atom_types = None
    atoms = []
    in_atoms = False
    with open(data_path) as fh:
        for line in fh:
            stripped = line.split("#")[0].strip()
            if not in_atoms:
                if stripped.endswith("atom types"):
                    n_atom_types = int(stripped.split()[0])
                for key in ("xlo xhi", "ylo yhi", "zlo zhi"):
                    if stripped.endswith(key):
                        lo, hi = stripped.split()[:2]
                        box[key.split()[0]] = (float(lo), float(hi))
                if stripped == "Atoms" or stripped.startswith("Atoms "):
                    in_atoms = True
                continue
            fields = stripped.split()
            if not fields:
                continue
            if not fields[0].isdigit():
                break
            atoms.append(fields)
    if n_atom_types is None or len(box) != 3:
        raise ValueError(f"missing header fields in {data_path}")

    n_org_types = n_atom_types - n_sub_types
    sub_zs = [float(f[6]) for f in atoms if int(f[2]) > n_org_types]
    if not sub_zs:
        raise ValueError(f"no substrate atoms (type > {n_org_types}) found in {data_path}")

    lz = box["zlo"][1] - box["zlo"][0]
    if abs(box["zlo"][0]) > 1e-6:
        raise ValueError(f"{data_path}: expected zlo==0, got {box['zlo'][0]} -- the arc-length math below assumes it")
    return lz, min(sub_zs), max(sub_zs)


def parse_density_profile(path):
    """fix ave/chunk density/number output -> (timesteps, coords, dens),
    dens shape (n_windows, n_bins), coords shape (n_bins,) -- bin centers,
    identical across windows (one shared `compute chunk/atom` for both
    species' fixes, lmp/production.lmp.template)."""
    with open(path) as fh:
        lines = fh.readlines()
    timesteps, coords, dens = [], None, []
    i, n = 0, len(lines)
    while i < n:
        line = lines[i].strip()
        if not line or line.startswith("#"):
            i += 1
            continue
        ts, nchunks, _total = line.split()
        nchunks = int(nchunks)
        block = lines[i + 1 : i + 1 + nchunks]
        c = np.array([float(l.split()[1]) for l in block])
        d = np.array([float(l.split()[3]) for l in block])
        if coords is None:
            coords = c
        timesteps.append(int(ts))
        dens.append(d)
        i += 1 + nchunks
    return np.array(timesteps), coords, np.array(dens)


def classify_bins(coords, lz, sub_zmin, sub_zmax, contact_width, bulk_margin):
    """Returns (upper_contact, lower_contact, bulk) boolean masks over
    coords, per the module docstring's arc-length parameterization."""
    in_substrate = (coords >= sub_zmin) & (coords <= sub_zmax)
    s = np.mod(coords - sub_zmax, lz)  # 0 at the upper face, arc_length at the lower
    arc_length = lz - (sub_zmax - sub_zmin)

    upper_contact = (~in_substrate) & (s <= contact_width)
    lower_contact = (~in_substrate) & (s >= arc_length - contact_width)
    bulk = (~in_substrate) & (s > bulk_margin) & (s < arc_length - bulk_margin)
    return upper_contact, lower_contact, bulk


def running_cumulative_mean(per_window_region_mean):
    """Cumulative mean over the time axis (axis 0), one value per window
    -- the last entry is the full-average estimate."""
    cs = np.cumsum(per_window_region_mean)
    counts = np.arange(1, len(per_window_region_mean) + 1)
    return cs / counts


def compute_species_chis(species, fragment_path, hexadecane_path, merged_data_path,
                          n_sub_types, a_angstrom, contact_width_a, bulk_margin_a,
                          min_contact_fraction, discard_fraction, convergence_tol):
    lz, sub_zmin, sub_zmax = read_data_header_and_substrate_z(merged_data_path, n_sub_types)

    frag_ts, coords, frag_dens = parse_density_profile(fragment_path)
    hex_ts, coords_hex, hex_dens = parse_density_profile(hexadecane_path)
    if not np.array_equal(coords, coords_hex) or not np.array_equal(frag_ts, hex_ts):
        raise ValueError(f"{species}: fragment/hexadecane profiles don't share the same bins/windows")

    n_windows = len(frag_ts)
    n_discard = int(round(discard_fraction * n_windows))
    frag_dens = frag_dens[n_discard:]
    hex_dens = hex_dens[n_discard:]
    n_retained = len(frag_dens)
    if n_retained < 2:
        raise ValueError(f"{species}: only {n_retained} windows left after discarding {n_discard}/{n_windows}")

    contact_width = contact_width_a * a_angstrom
    bulk_margin = bulk_margin_a * a_angstrom
    upper_contact, lower_contact, bulk = classify_bins(
        coords, lz, sub_zmin, sub_zmax, contact_width, bulk_margin
    )
    if not bulk.any():
        raise ValueError(f"{species}: no bins classified as bulk -- widen the liquid region or shrink chis_bulk_margin_a")

    hex_mean_profile = hex_dens.mean(axis=0)
    bulk_ref = float(np.median(hex_mean_profile[hex_mean_profile > 1e-4]))

    faces = {"upper": upper_contact, "lower": lower_contact}
    valid_face_masks = []
    face_diagnostics = {}
    for name, mask in faces.items():
        if not mask.any():
            face_diagnostics[name] = {"n_bins": 0, "valid": False}
            continue
        hex_face_density = float(hex_mean_profile[mask].mean())
        valid = hex_face_density >= min_contact_fraction * bulk_ref
        face_diagnostics[name] = {
            "n_bins": int(mask.sum()),
            "hexadecane_density": hex_face_density,
            "valid": bool(valid),
        }
        if valid:
            valid_face_masks.append(mask)

    if not valid_face_masks:
        raise ValueError(
            f"{species}: neither substrate face has a real liquid contact "
            f"(bulk_ref={bulk_ref:.4g}, see face diagnostics) -- the production "
            "run's geometry doesn't support a contact-layer estimate at all"
        )
    contact_mask = np.logical_or.reduce(valid_face_masks)

    frag_contact = frag_dens[:, contact_mask].mean(axis=1)
    frag_bulk = frag_dens[:, bulk].mean(axis=1)
    hex_contact = hex_dens[:, contact_mask].mean(axis=1)
    hex_bulk = hex_dens[:, bulk].mean(axis=1)

    frag_contact_cum = running_cumulative_mean(frag_contact)
    frag_bulk_cum = running_cumulative_mean(frag_bulk)
    hex_contact_cum = running_cumulative_mean(hex_contact)
    hex_bulk_cum = running_cumulative_mean(hex_bulk)

    with np.errstate(divide="ignore", invalid="ignore"):
        chis_running = (
            np.log(frag_contact_cum / frag_bulk_cum)
            - np.log(hex_contact_cum / hex_bulk_cum)
        )

    chis_final = float(chis_running[-1])
    mid = n_retained // 2
    chis_mid = float(chis_running[mid]) if np.isfinite(chis_running[mid]) else None
    converged = (
        chis_mid is not None
        and np.isfinite(chis_final)
        and abs(chis_final - chis_mid) <= convergence_tol
    )

    return chis_final, {
        "n_windows_total": int(n_windows),
        "n_windows_discarded": int(n_discard),
        "n_windows_retained": int(n_retained),
        "substrate_zmin_A": sub_zmin,
        "substrate_zmax_A": sub_zmax,
        "box_lz_A": lz,
        "bulk_reference_density": bulk_ref,
        "faces": face_diagnostics,
        "chis_x_final": chis_final if np.isfinite(chis_final) else None,
        "chis_x_at_midpoint": chis_mid,
        "converged_heuristic": converged,
        "phi_contact_fragment": float(frag_contact_cum[-1]),
        "phi_bulk_fragment": float(frag_bulk_cum[-1]),
        "phi_contact_hexadecane": float(hex_contact_cum[-1]),
        "phi_bulk_hexadecane": float(hex_bulk_cum[-1]),
    }


def main():
    params = snakemake.params
    with open(snakemake.input.segment_volume) as fh:
        seg_vol = yaml.safe_load(fh)
    a_angstrom = seg_vol["a3_nm3"] ** (1.0 / 3.0) * 10.0

    n_sub_types = params.n_sub_types

    chis = {}
    diagnostics = {}
    for species, fragment_path, hexadecane_path, merged_data_path in zip(
        params.species,
        snakemake.input.fragment_density,
        snakemake.input.hexadecane_density,
        snakemake.input.merged_data,
    ):
        chis_x, diag = compute_species_chis(
            species, fragment_path, hexadecane_path, merged_data_path,
            n_sub_types, a_angstrom,
            params.contact_width_a, params.bulk_margin_a,
            params.min_contact_fraction, params.discard_fraction,
            params.convergence_tol,
        )
        label = SUBSTRATE_LABEL[species]
        chis[label] = chis_x
        diagnostics[species] = diag

    with open(snakemake.output.chis, "w") as fh:
        yaml.safe_dump(chis, fh, sort_keys=False)
    with open(snakemake.output.diagnostics, "w") as fh:
        yaml.safe_dump(diagnostics, fh, sort_keys=False)


if __name__ == "__main__":
    main()
