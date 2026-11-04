# chi/chis-from-MD pipeline

Snakemake implementation of the plan in `../forcefield/chi_chis_from_md.md`.
Standalone and self-contained inside `sf1d/` — it does not touch or extend
`zhanglab-nagoya/2026-01-26-generate-aa-pama` or the other related repos on
this machine (`2026-07-28-fe2o3-cg`, `2026-06-17-cg-parameter-optimization`),
even though those do overlapping AA-reference-characterization work for the
same physical system; that was a deliberate choice, not an oversight — see
`chi-chis-from-md-plan` in this project's Claude memory if that choice ever
needs revisiting.

## Status: what's real vs. design-only

| Stage | File | Status |
|---|---|---|
| Step 0 — segment size | `rules/step0_step1_bulk.smk` (`bulk_*`, `segment_volume`) | **Runnable.** Dry-run verified (`snakemake -n`), not yet executed. |
| Step 1 — `chi_AS`/`chi_BS`, cohesive-energy-density route | `rules/step0_step1_bulk.smk` (`vacuum_ref`, `solubility_chi`) | **Not implemented.** `vacuum_ref` needs a real vacuum single-molecule mdp/box; `solubility_chi` raises `NotImplementedError` on purpose until then. |
| Step 2 — `chis_A`/`chis_B`, dilute slab adsorption | `rules/step2_chis.smk` | **Design only.** Needs a GROMACS-build → LAMMPS-convert → merge-with-`fe2o3.data` pipeline; see the module docstring for the intended stages and why this wasn't fabricated blind. |
| Step 3 — validate against MD | `rules/step3_validate.smk` | `sf1d_predict` is real (pure Python, no MD) but blocked on step 2's `results/chis.yaml`; the MD side is design-only, same status as step 2. |

Nothing in this pipeline has been executed yet — only dry-run (`-n`) checked.
Running it for real means committing real wall-clock time on this (single,
shared) workstation; review `config.yaml`'s run lengths and box sizes first,
they're placeholders sized for "does this even work", not for production
statistics.

## Running it

```bash
cd pipeline
conda run -n snakemake snakemake --use-conda -c 4 segment_volume   # the one fully-wired target
```

`--use-conda` creates isolated environments from `envs/*.yaml` under
`.snakemake/conda/` (same pattern as the user's other dated project repos —
this never touches the persistent named conda envs like `rdkit` or
`psma-bi`). GROMACS itself is *not* conda-managed: rules call `gmx_mpi`
directly, assuming `~/opt/gromacs-2026.0/bin` is on `PATH` (it has no plain
`gmx`, only `gmx_mpi` — single-rank invocation works fine for everything
here) and using `../forcefield/oplsaa.ff` (the L-OPLS-patched copy) rather
than the stock installed force field, which lacks that patch.

## Design notes

- **Species for step0/step1:** `HEX` (hexadecane, already a complete itp —
  no polyply needed), `ST`, `MMA` (capped oligomers,
  `CH3-(monomer)_n-CH3n`, built with `polyply gen_params`/`gen_coords` from
  `forcefield/monomers/` + `forcefield/PAMA_ST_MMA.oplsaa.LigParGen_links.ff`).
  `oligomer_length` (config.yaml, default 10) is *not* the real chain
  length — it's just long enough that the end caps are a small fraction of
  the melt's mass, for a reasonable per-repeat-unit molar volume estimate.
- **Molar masses are read from the force field files, not hardcoded** —
  `scripts/segment_volume.py` sums the `[ atoms ]` mass column straight out
  of `forcefield/monomers/*.ff` and `forcefield/solvent/hexadecane.itp`.
- **polyply** isn't in any of this machine's named conda envs (checked);
  `envs/polyply.yaml` installs it from PyPI into its own isolated env rather
  than depending on `~/git/polyply_1.0/venv` (which does have a working
  install — useful for quickly checking CLI syntax by hand, e.g.
  `polyply gen_params --help`, but not a dependency this pipeline should
  have on an unrelated sibling repo's venv).
- **Why GROMACS for step0/1 but LAMMPS for step2:** the Fe2O3 substrate
  (`forcefield/substrate/fe2o3.data`) only exists as a LAMMPS data file
  (`pair_style lj/cut/coul/long`), prepared externally (Ji Ce, zhanglab) —
  not something to rebuild from scratch in GROMACS. Step 2 has to bridge
  the two, which is exactly what `rules/step2_chis.smk`'s docstring lays
  out as still needing real design work.
