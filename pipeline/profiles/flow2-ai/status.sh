#!/usr/bin/env bash
# Status-check cmd for Snakemake's cluster-generic executor on flow2-ai
# (plain PBS Professional -- no RCCS-style jsub/jobinfo/joblog wrappers
# here, so this is a direct `qstat -f` parse, simpler than
# ../rccs/status-jobinfo.sh).
#
# `-x` (job history) is required, not optional `qstat -f`: confirmed live
# 2026-11-09 that a short job (lmp_slab_minimize, ~5-10s real time) can
# already have dropped out of the *live* queue by the time Snakemake's
# first status poll fires, while `qstat -x -f` still returns its full
# record (job_state = F, Exit_status = 0) via PBS's retained job history.
# Without `-x`, that race made Snakemake treat a successfully-finished job
# as failed and delete its (correct) output.
#
# Usage: status.sh <jobid>   (jobid as printed by qsub, e.g. "12345.jobm")
# Prints exactly one of: running | success | failed
set -eu

jobid="$1"

output=$(qstat -x -f "$jobid" 2>/dev/null) || true

if [[ -z "$output" ]]; then
    # Even with -x, a job could in principle age out of history or the
    # jobid could be bogus -- keep this fallback, now just a safety net
    # rather than the common case it was before -x.
    echo failed
    exit 0
fi

job_state=$(echo "$output" | awk -F'= ' '/[[:space:]]job_state =/ {print $2}' | tr -d ' \n\r')

case "$job_state" in
    F)
        exit_status=$(echo "$output" | awk -F'= ' '/[[:space:]]Exit_status =/ {print $2}' | tr -d ' \n\r')
        if [[ "$exit_status" == "0" ]]; then
            echo success
        else
            echo failed
        fi
        ;;
    R|Q|H|E|B|W|S|"")
        # Running/queued/held/exiting/begun/waiting/suspended, or a state we
        # didn't parse -- never guess "failed" here, let Snakemake keep polling.
        echo running
        ;;
    *)
        echo running
        ;;
esac
