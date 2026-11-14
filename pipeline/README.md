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
| Step 0 — segment size | `rules/step0_step1_bulk.smk` (`bulk_*`, `segment_volume`) | **Runs end-to-end.** `a3_nm3=0.5015 nm^3`, HEX density 749.8 kg/m^3 vs. a ~770 kg/m^3 literature reference (normal few-percent OPLS/L-OPLS deviation). |
| Step 1 — `chi_AS`/`chi_BS`, cohesive-energy-density route | `rules/step0_step1_bulk.smk` (`vacuum_nvt`, `solubility_chi`) | **Runs end-to-end**, but see that file's STATUS note on sampling-length caveats before trusting the numbers quantitatively (`chi_AS=0.83`, `chi_BS~0` as of 2026-11-05 — the chi_AS > chi_BS direction is counterintuitive and likely under-sampled, not a bug). |
| Step 2 substrate — Fe2O3 from scratch | `rules/step2_substrate.smk` | **Runs end-to-end.** Builds `results/substrate/substrate_minimized.data` (8400 atoms, same CIF->ASE->INTERFACE-FF->LAMMPS-CG-minimize method as `../2026-11-05-pama-in-nanogap`'s `299_substrate_creation`) instead of reading the legacy, opaque-provenance `forcefield/substrate/fe2o3.data`. |
| Step 2 chis — dilute slab adsorption | `rules/step2_dilute.smk` (solvent GROMACS build + fragment GROMACS build), `rules/step2_convert.smk` (GROMACS->LAMMPS), `rules/step2_merge.smk` (solvent+substrate merge), `rules/step2_equilibrate.smk` (equilibration + fragment insertion + production), `rules/step2_chis.smk` (contact-vs-bulk extraction) | **Partly real, redesigned 2026-11-13, not yet re-verified end-to-end.** Three real bugs were found and fixed across 2026-11-08/11/12 in the original design (fragment packed into the solvent box from the start, restrained through the whole barostat run): atom types never deduplicated; the z-barostat driven by LAMMPS's whole-box `compute pressure` (hard-requires group "all", so it always included the frozen substrate's own enormous virial) instead of a group-scoped one; `insert_fragment.py`'s rigid translation leaving a real vacancy at the fragment's origin (fixed by centering the fragment pre-merge instead, `center_fragment.py`). That last fix then surfaced a FOURTH, more fundamental bug: restraining the fragment through the whole 100 ps/400-cycle `lmp_npt_z` barostat run made it a fixed anchor the barostat's own `change_box ... remap` (anchored at a fixed `zlo`) stretched hexadecane away from every cycle, opening a real, permanent cavity at the fragment's location — confirmed across three independently-tuned barostat variants (isotropic vs. Pzz-only pressure formula, naive vs. mass-conservation box sizing) that all converged to the same wrong ~271 Å/~1500 atm plateau regardless of starting conditions. **2026-11-13 redesign**: decouple the fragment from the solvent/substrate equilibration entirely. `rules/step2_dilute.smk`'s `solvent_pack`/`solvent_top`/`solvent_em`/`solvent_nvt`/`solvent_npt` build a hexadecane-ONLY box (no fragment, no `{species}` wildcard — built once, shared by ST and MMA), merged onto the substrate (`merge_solvent_with_substrate`) and run through the full `lmp_slab_minimize -> lmp_nvt -> close_liquid_gap -> lmp_gap_minimize -> lmp_settle` chain with no fragment present anywhere, so there is no restrained anchor for a barostat to stretch away from — moot anyway since `lmp_npt_z`'s barostat (both the whole-box and later liquid-only-pressure hand-rolled versions) was removed entirely 2026-11-13 after confirming live it's a structural dead end on this frozen-substrate system regardless of tuning (same conclusion `../2026-07-28-fe2o3-cg` reached independently on the identical substrate+solvent pairing); `lmp_settle` (plain, box-fixed NVT, no pressure feedback) replaces it, since `close_liquid_gap`+`lmp_gap_minimize` already leave the box at its correct, deterministically-sized density. Only once that converges does each species' own fragment (built standalone via `convert_fragment_to_lammps`) get combined with the equilibrated solvent by a recreated `scripts/insert_fragment.py` (new design: combines two independently-converted LAMMPS files, clash-deletes whichever hexadecane molecules the translated fragment overlaps, reuses `merge_slab.py`'s type-offsetting pattern) — `lmp_insert_minimize` relaxes the local clash, a short restrained `lmp_insert_relax` (box fixed, no barostat — deliberately, to avoid reopening the same cavity mechanism) lets it thermalize, then `lmp_production` releases the restraint. `center_fragment.py`/the old `insert_fragment.py` design are deleted. See `rules/step2_equilibrate.smk`'s module docstring and the "Decouple fragment insertion from solvent/substrate equilibration" plan for the full history and the verification checklist this redesign still needs to clear (same checklist also doubles as a direct test of the cavity diagnosis above: if the fragment-free solvent system converges cleanly to ~1 atm, that confirms it). Stage 6 (`compute_chis`, contact-vs-bulk density ratio) is written and runs correctly, but no real production dataset has been clean enough to trust its output yet. |
| Step 3 — validate against MD | `rules/step3_validate.smk` | `sf1d_predict` is real (pure Python, no MD) but blocked on step 2 chis' `results/chis.yaml`; the MD side is design-only, same status as step 2 chis. |

## Running it

```bash
cd pipeline
conda run -n snakemake snakemake --use-conda -c 8 segment_volume      # step 0
conda run -n snakemake snakemake --use-conda -c 8 results/chi.yaml    # step 0 + step 1
source ~/opt/lammps-22Jul2025/env.sh  # LD_LIBRARY_PATH for liblammps.so.0
conda run -n snakemake snakemake --use-conda -c 4 results/substrate/substrate_minimized.data  # step 2 substrate
conda run -n snakemake snakemake --use-conda -c 8 results/slab/ST/npt.edr results/slab/MMA/npt.edr        # step 2 chis stage 1
conda run -n snakemake snakemake --use-conda -c 4 results/slab/ST/converted.lmp results/slab/MMA/converted.lmp  # step 2 chis stage 2
```

Stage 2 (`convert_to_lammps`) needs InterMol, available via `module load
intermol` -- but *not* loaded globally before invoking snakemake like the
other module-based tools: the rule loads it itself in its own shell
subprocess (`bash -lc 'module load intermol && ...'`), because
`module load intermol` points `PYTHONPATH` at a spack-built numpy for
Python 3.14, and leaking that into the whole snakemake process (its own
conda env, Python 3.13) breaks snakemake's own numpy/pandas imports --
confirmed live 2026-11-05 (see `rules/step2_convert.smk`'s docstring).

The bulk GROMACS rules need `threads: 8` worth of free cores (set via `-c`)
to get real throughput — without it Snakemake pins every `mdrun` to 1
thread regardless of the machine's core count (bitten this exact pipeline
once already; see `rules/step0_step1_bulk.smk`'s `bulk_squeeze`/`bulk_npt`
history).

`--use-conda` creates isolated environments from `envs/*.yaml` under
`.snakemake/conda/` (same pattern as the user's other dated project repos —
this never touches the persistent named conda envs like `rdkit` or
`psma-bi`). GROMACS itself is *not* conda-managed: rules call `gmx_mpi`
directly, assuming `~/opt/gromacs-2026.0/bin` is on `PATH` (it has no plain
`gmx`, only `gmx_mpi` — single-rank invocation works fine for everything
here) and using `../forcefield/oplsaa.ff` (the L-OPLS-patched copy) rather
than the stock installed force field, which lacks that patch.

### On flow2-ai (HPC)

`profiles/flow2-ai/` runs the same rules on Nagoya University's Flow 2.0
supercomputer (PBS Professional) instead of this workstation — see
`profiles/flow2-ai/README.md` for setup and scope. Only step 0/1's bulk
GROMACS rules go through PBS there today; step 2 has no rules yet.

```bash
snakemake --workflow-profile profiles/flow2-ai --software-deployment-method conda segment_volume
```

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
  (`results/substrate/substrate_minimized.data`, built by
  `rules/step2_substrate.smk`) is a LAMMPS data file (`pair_style
  lj/cut/coul/long`, INTERFACE force field) — GROMACS has no equivalent FF
  for this material. Step 2 chis has to bridge the two (GROMACS-built
  organic fragment + LAMMPS-built substrate), which is exactly what
  `rules/step2_chis.smk`'s docstring lays out as still needing real design
  work. `forcefield/substrate/fe2o3.data` (prepared externally, Ji Ce,
  zhanglab) is now legacy/unused — see `forcefield/substrate/README.md`.
