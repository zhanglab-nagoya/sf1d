#!/usr/bin/env python3
"""Submit cmd for Snakemake's cluster-generic executor, targeting RCCS
(ccfep.ims.ac.jp) via its own `jsub` wrapper -- not plain PBS `qsub` (unlike
../flow2-ai/submit.py): jsub does not accept `-l` resource arguments on its
CLI, so the #PBS header has to be injected into the jobscript text itself
before handing it to jsub. Modeled directly on
2026-02-02-pc-on-sio2/workflow/profiles/rccs/submit.py (confirmed working
there), simplified for this pipeline's own resource convention: every
lmp_*/bulk_* rule here runs ONE process (`lmp -sf omp -pk omp {ompthreads}`
or plain `gmx`/`gmx_mpi` without MPI decomposition -- see
rules/step2_equilibrate.smk and rules/step0_step1_bulk.smk's shell blocks),
so `mpiprocs` is always 1 and `ompthreads` is the only thing that varies
per rule, not pc-on-sio2's tasks_per_node/mpi_tasks split.
"""

import sys
import subprocess
from pathlib import Path

from snakemake.utils import read_job_properties


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: submit.py <jobscript>")

    jobscript = Path(sys.argv[1])
    props = read_job_properties(str(jobscript))
    resources = props.get("resources", {})
    rule = props.get("rule", "unknown")
    jobid = props.get("jobid", "unknown")
    wildcards = props.get("wildcards", {})
    wc_str = ".".join(f"{k}={v}" for k, v in wildcards.items()) if wildcards else str(jobid)

    def res(name, default=None):
        return resources.get(name, default)

    nodes = res("nodes", 1)
    ncpus = res("ncpus", 1)
    ompthreads = res("ompthreads", ncpus)
    walltime = res("time", "01:00:00")

    select_line = f"#PBS -l select={nodes}:ncpus={ncpus}:mpiprocs=1:ompthreads={ompthreads}"

    logdir = Path(f"logs/{rule}")
    logdir.mkdir(parents=True, exist_ok=True)

    pbs_header_lines = [
        select_line,
        f"#PBS -l walltime={walltime}",
        f"#PBS -o logs/{rule}/{rule}-{wc_str}.out",
        "#PBS -j oe",
        "",
    ]

    lines = jobscript.read_text().splitlines()
    if not lines:
        raise ValueError("Jobscript is empty")
    if not lines[0].startswith("#!"):
        lines.insert(0, "#!/bin/bash")

    if any(line.startswith("#PBS") for line in lines):
        # Re-submission of an already-materialized jobscript (Snakemake
        # reuses the same temp file on retry) -- don't double-inject.
        new_lines = lines
    else:
        new_lines = [lines[0]] + pbs_header_lines + lines[1:]

    jobscript.write_text("\n".join(new_lines) + "\n")

    cmd = ["jsub", "-N", f"{rule}-{jobid}", str(jobscript)]

    print("Submitting:", " ".join(cmd), file=sys.stderr)

    # jsub prints the PBS job ID to stdout -- do not capture; Snakemake reads it.
    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
