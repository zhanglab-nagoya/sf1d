#!/usr/bin/env bash
# Status-check cmd for Snakemake's cluster-generic executor on flow2-ai
# (plain PBS Professional -- no RCCS-style jsub/jobinfo/joblog wrappers
# here, so this is a direct `qstat -f` parse, simpler than
# ../rccs/status-jobinfo.sh).
#
# Usage: status.sh <jobid>   (jobid as printed by qsub, e.g. "12345.jobm")
# Prints exactly one of: running | success | failed
set -eu

jobid="$1"

output=$(qstat -f "$jobid" 2>/dev/null) || true

if [[ -z "$output" ]]; then
    # PBS normally keeps a finished job visible (job_state = F) for a while
    # before it drops out of qstat entirely. A job that has already vanished
    # by the time we first check it is treated as failed -- Snakemake has no
    # other way to recover its exit code at that point.
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
