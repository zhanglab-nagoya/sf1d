# flow2-ai profile

Runs `pipeline/` on **Flow 2.0** (不老・弐), Nagoya University's PBS
Professional supercomputer, via its generative-AI login nodes
(`flow2-ai.cc.nagoya-u.ac.jp`, host `~/.ssh/config` entry `flow2-ai`).

## What's covered

Only step 0/1's bulk GROMACS rules (`bulk_em`, `bulk_nvt`, `bulk_npt`,
`bulk_density`) go through PBS here — the same scope `rule all` targets
today (see `../../README.md`). Step 2 (`rules/step2_chis.smk`, LAMMPS +
Fe2O3 slab) has no rules yet, so there's nothing to schedule for it; once
it's written, size its queue/walltime from a real run, not a guess (see the
note at the bottom of `config.yaml`).

## One-time setup on flow2-ai

Checked live on flow2-ai (2026-11-04): `gromacs/2026.3` and `lammps/29Aug2024`
are site modules (both need `oneapi impi` loaded first — done unconditionally
in `pbs-jobscript.sh`), but **`snakemake`, `polyply`, `intermol` and `packmol`
are not** — there's no module for any of them, and no conda env exists yet
under `~/.conda/envs`. Build one with the site's `miniforge3` module:

```bash
ssh flow2-ai
module load miniforge3
source $MINIFORGE3_DIR/etc/profile.d/conda.sh
conda create -n sf1d-pipeline -c conda-forge python=3.11 snakemake-minimal pandas numpy scipy pyyaml packmol
conda activate sf1d-pipeline
pip install polyply   # same as pipeline/envs/polyply.yaml, from PyPI
```

`intermol` isn't on PyPI in a form that just works (the local workstation
uses a patched jotelha fork built via spack — see
`/home/jotelha/CLAUDE.md`'s custom package overlay). It's only needed once
step 2 exists; don't install it speculatively now.

Then get the repo onto flow2-ai (it's a normal git clone — `forcefield/`,
including the ~11 MB `fe2o3.data`, is already committed):

```bash
git clone git@github.com:zhanglab-nagoya/sf1d.git ~/sf1d   # or scp/rsync if no outbound git+ssh from the login node
```

## Running

```bash
ssh flow2-ai
module load miniforge3 && source $MINIFORGE3_DIR/etc/profile.d/conda.sh && conda activate sf1d-pipeline
cd ~/sf1d/pipeline
snakemake -n --workflow-profile profiles/flow2-ai --software-deployment-method conda segment_volume   # dry-run first, always
snakemake --workflow-profile profiles/flow2-ai --software-deployment-method conda -c 2 segment_volume
qstat -u $USER            # watch jobs
tail -f logs/bulk_npt/*.out
```

## Design notes

- **PBS, not RCCS/Slurm.** flow2-ai has a plain, unwrapped PBS Professional
  (`qsub`/`qstat`/`qdel`) — unlike the sibling generate-aa-pama repo's `rccs`
  profile, which wraps RCCS's own `jsub`/`jdel`/`jobinfo`/`joblog`. So
  `submit.py` and `status.sh` here talk to `qsub`/`qstat -f` directly.
- **`localrule: True`** is set (in `../../rules/*.smk`, not here) on every
  rule that's cheap enough to run on the login node directly: topology
  generation (`build_itp_oligomer`, needs the `polyply` conda env, not
  GROMACS), file copies/writes (`bulk_itp_hex`, `bulk_top`), box packing
  (`bulk_pack`), and the pure-Python analysis rules (`segment_volume`,
  `vacuum_ref`, `solubility_chi`, `sf1d_predict`). Only the actual MD
  dynamics (`bulk_em/nvt/npt`) and the `gmx energy` extraction
  (`bulk_density`, needs the `gromacs` module) go through PBS.
- **Queue: `t1-share`** (CPU, NUMA-granularity — 1 of 8 NUMA domains per
  node, 48 cores each). The step0/1 bulk boxes are tiny pure-liquid/melt
  systems (`config.yaml`'s `bulk.box_nm`, generous 6×6×6 nm — not the full
  Fe2O3 substrate), so a shared, cheap queue is the right tier; see
  flow2-ai's own `~/CLAUDE.md` §7 for the points-accounting rates that make
  exclusive-node queues much more expensive for jobs this small.
- **Walltimes are unbenchmarked placeholders** (same caveat as
  `../../README.md` for the box sizes/run lengths themselves) — pad them if
  a real run gets killed; `status.sh`/`qstat -f <id>` reports `Exit_status`
  either way.
- **GROMACS→LAMMPS conversion for step 2**, when it's written: the sibling
  `2026-10-28-generate-aa-pama` repo's `workflow/050_conversion.smk` shows
  the working pattern (`gmx trjconv -pbc whole` → `intermol-convert
  --gro_in ... --lammps --lmppath $(which lmp) --gropath $(dirname $(which
  gmx))`, after sourcing `GMXRC`) — useful as a reference for what step 2
  needs to reproduce here, not something this pipeline calls directly (see
  `../../rules/step2_chis.smk`'s own docstring for why that stage is
  deliberately not fabricated blind yet).
