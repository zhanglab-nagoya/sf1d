#!/usr/bin/env python3
"""Submit cmd for Snakemake's cluster-generic executor, targeting flow2-ai's
PBS Professional scheduler directly via its own `qsub` (no site wrapper --
unlike ../rccs/submit.py, which wraps RCCS's jsub; flow2-ai has plain PBS).

Builds the `-l select=...` resource string from whichever resource keys the
rule declared, following flow2-ai's own ~/CLAUDE.md section 6.2:
  - nnumas       -> shared NUMA-granularity queues (t1-share, t1-lm-small):
                    `select=1:nnumas=<N>:ompthreads=<M>`
  - ngpus        -> shared GPU-granularity queues (t2*-share/-mig, t2-dev):
                    `select=1:ngpus=<N>:ompthreads=<M>`
  - otherwise    -> exclusive whole-node queues (t1-single/small/..., t2-*):
                    `select=<nodes>[:mpiprocs=<M>][:ompthreads=<M>]` --
                    deliberately NOT adding ncpus/mem: PBS expands a bare
                    `select=<nodes>` to the full node automatically, and
                    specifying it yourself can be rejected as out of range.
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

    queue = res("queue", "t1-share")
    walltime = res("time", "01:00:00")

    if "nnumas" in resources:
        select = f"select=1:nnumas={res('nnumas', 1)}:ompthreads={res('ompthreads', 48)}"
    elif "ngpus" in resources:
        select = f"select=1:ngpus={res('ngpus', 1)}:ompthreads={res('ompthreads', 36)}"
    else:
        select = f"select={res('nodes', 1)}"
        mpiprocs = res("mpiprocs")
        ompthreads = res("ompthreads")
        if mpiprocs:
            select += f":mpiprocs={mpiprocs}"
        if ompthreads:
            select += f":ompthreads={ompthreads}"

    logdir = Path(f"logs/{rule}")
    logdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "qsub",
        "-q", queue,
        "-l", select,
        "-l", f"walltime={walltime}",
        "-N", f"smk-{rule}",
        "-o", f"logs/{rule}/{rule}-{wc_str}.out",
        "-j", "oe",
        str(jobscript),
    ]

    print("Submitting:", " ".join(cmd), file=sys.stderr)

    # qsub prints the new job ID ("12345.jobm") to stdout on success --
    # Snakemake reads that directly, so stdout must not be captured/altered.
    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
