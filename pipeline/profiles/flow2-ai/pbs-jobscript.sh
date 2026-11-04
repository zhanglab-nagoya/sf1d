#!/bin/bash -l
# properties = {properties}
#
# PBS Professional jobscript for the sf1d chi/chis-from-MD pipeline on
# flow2-ai (Flow 2.0, Nagoya University ITC). "#!/bin/bash -l" is required:
# without a login shell, `module` and the compute-node HTTP(S) proxy are not
# initialized (flow2-ai's own ~/CLAUDE.md, section 6.4).
#
# Modules are loaded unconditionally here (not per-rule via Snakemake's
# envmodules: directive, which this pipeline's rules don't use) -- cheap,
# and keeps every rule's shell block free of site-specific module names.
# gromacs and lammps both need oneapi+impi loaded first (`show_module -k
# gromacs|lammps` on flow2-ai).

echo " "
echo "### PBS job environment:"
echo " "
echo "Job ID:            $PBS_JOBID"
echo "Job name:          $PBS_JOBNAME"
echo "Submit dir:        $PBS_O_WORKDIR"
echo "Nodes file:        $PBS_NODEFILE"
echo "OMP_NUM_THREADS:   $OMP_NUM_THREADS"
echo " "

cd "$PBS_O_WORKDIR"
pwd

echo " "
echo "### Load modules:"
echo " "
module purge
module load oneapi impi
module load gromacs lammps
module list

echo " "
echo "### Run rule:"
echo " "

{exec_job}
