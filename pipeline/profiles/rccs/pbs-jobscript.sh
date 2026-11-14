#!/bin/bash
# properties = {properties}
#
# PBS Pro jobscript for the sf1d chi/chis-from-MD pipeline on RCCS
# (ccfep.ims.ac.jp, Institute for Molecular Science). Submitted via gxy's
# own jsub wrapper (profiles/rccs/submit.py), which injects the #PBS header
# itself -- this script carries no #PBS lines of its own, matching the
# pattern already proven working for 2026-02-02-pc-on-sio2's own RCCS
# profile (workflow/profiles/rccs/pbs-jobscript.sh).
#
# Modules are loaded unconditionally here (not per-rule via Snakemake's
# envmodules: directive, which this pipeline's rules don't use anywhere --
# same reasoning as profiles/flow2-ai/pbs-jobscript.sh). lammps/gromacs are
# RCCS's own system modules (no custom spack build needed -- this pipeline
# only uses standard pair_style lj/cut/coul/long + pppm + bond/angle/
# dihedral hybrid harmonic, unlike pc-on-sio2's `fix controller` need).
# packmol/intermol are gxy's private modules (RCCS has no system module for
# either -- see root CLAUDE.md's RCCS deployment section); included here for
# completeness even though today's PBS-submitted rules (the lmp_* dynamics
# stages) don't need them -- only this pipeline's localrule:True build/
# convert scripts do, and those run on the login node, not inside this job.
#
# lammps/2025-Jul22-u3 (RCCS system module) auto-loads intelmpi/2021.16,
# which conflicts with gromacs/2025.4's own openmpi/4.1.8/gcc13 dependency
# -- confirmed live 2026-11-14, `module load gromacs/2025.4` after lammps
# errors out ("Load of requirement openmpi/4.1.8/gcc13 failed") and gmx/
# gmx_mpi never end up on PATH. This is harmless for lmp_* rules (lammps
# itself loads and runs fine regardless) but means bulk_*/gmx_* rules
# cannot run under this jobscript as written. pc-on-sio2 avoids this by
# using its own private spack-built lammps/20250722, whose MPI stack is
# resolved to openmpi instead of intelmpi -- same fix here if a GROMACS
# rule ever needs to run on RCCS for real, swapping the system lammps
# module for an equivalent private build.

echo " "
echo "### PBS job environment:"
echo " "
echo "Job ID:            $PBS_JOBID"
echo "Job name:          $PBS_JOBNAME"
echo "Submit dir:        $PBS_O_WORKDIR"
echo "Nodes file:        $PBS_NODEFILE"
echo "Nodes:             $(sort -u "$PBS_NODEFILE" | wc -l)"
echo "Tasks (PBS_NP):    $PBS_NP"
echo " "

cd "$PBS_O_WORKDIR"
pwd

module use --append "$HOME/privatemodules/linux-rocky8-zen3"
module use --append "$HOME/privatemodules/linux-rocky8-zen"

module load lammps/2025-Jul22-u3
module load gromacs/2025.4
module load packmol/20.0.0
module load intermol/2026-01-28-use-gmx-mpi

echo " "
echo "### Loaded modules:"
echo " "
module list

echo " "
echo "### Run rule:"
echo " "

{exec_job}
