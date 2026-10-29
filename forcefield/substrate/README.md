# Fe2O3 substrate

Extracted unchanged from `zhanglab-nagoya/2026-01-26-generate-aa-pama`,
`in/300_substrate_preparation/` (default branch `main`).

`fe2o3.data` is a LAMMPS data file (`write_data`, LAMMPS 22Jul2025,
`units real`, `atom_style full`, `pair_style lj/cut/coul/long`):
126,720 atoms, 2 atom types (Fe mass 55.845, O mass 15.9994), box
161.5 × 157.4 × 69.8 Å. Created and minimized by Ji Ce (zhanglab), 2026-05-12;
per the source repo, minimization inputs/outputs live on an SMB share
(`zhanglab5.local/public/JiC/1-Research/2026-04-21-sub-minimize`), not in
that git repo or here. Treat this as an externally-prepared, already-
equilibrated slab — not something to regenerate from the `.data` file alone.

`sort.tcl` (VMD/TopoTools) is the companion atom-reordering script from the
source repo, kept alongside for provenance; not run here.

Not copied: `fe2o3_probe.data` (`in/325_probe_preparation/`, 22,452-atom
Fe2O3 probe tip for indentation, same provenance/author) — only the flat
substrate was requested. Pull it in the same way
(`in/325_probe_preparation/fe2o3_probe.data`) if the indentation/probe
geometry is needed later.

## Note on repo size

This single file is ~11 MB, making `forcefield/` the large majority of this
repository's size (the rest of `sf1d` is well under 1 MB). It hasn't been
committed yet — flag if you'd rather keep large MD inputs like this out of
git history (e.g. reference by path to the source repo instead, or a
`.gitattributes`/LFS setup) before running `git add`.
