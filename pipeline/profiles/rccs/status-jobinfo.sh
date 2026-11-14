#!/usr/bin/env bash
# Check status of a PBS Pro job for Snakemake cluster-generic executor.
# This system uses RCCS wrappers — there is no qstat command.
#
# Usage: status-jobinfo.sh <jobid>
# Prints one of: running | success | failed
#
# jobid may be the full PBS job ID (e.g. 96293.ccpbs1) or just the numeric
# part.  jobinfo and joblog both use the numeric part only.
#
# Two-source jobinfo strategy
# ────────────────────────────
# "jobinfo -c" is real-time: jobs appear and disappear promptly, and can
#   show "Exit" state briefly just after a job finishes.
# "jobinfo" (plain) is a stale cache (~1–2 min lag): newly submitted jobs
#   may not yet appear in -c but are visible here once the cache picked them
#   up; finished jobs still show "Run" until the cache refreshes.
#
# The script therefore checks -c first (authoritative), then falls back to
# plain jobinfo (catches the -c registration-lag window for new jobs, and
# also returns "running" for recently finished jobs before the cache clears).
#
# joblog lag
# ──────────
# joblog has a significant lag (potentially tens of minutes) after a job
# finishes before its record appears.  To avoid false "failed" reports during
# that window the script returns "running" whenever the job is absent from
# joblog, so Snakemake keeps retrying until the record appears.
#
# State machine
# ─────────────
#  jobinfo -c → Run/Queue/Hold                   → running
#  jobinfo -c → Exit, joblog → 0                 → success
#  jobinfo -c → Exit, joblog → non-zero          → failed
#  jobinfo -c → Exit, joblog → absent            → running  (joblog lag; retry)
#  jobinfo -c → absent, plain jobinfo → any Run  → running  (registration lag
#                                                             OR stale cache for
#                                                             recently finished job)
#  both absent, joblog → 0                        → success
#  both absent, joblog → non-zero                 → failed
#  both absent, joblog → absent                   → running  (conservative; retry)
set -eu

jobid="$1"

# Strip hostname suffix: "96293.ccpbs1" → "96293"
numeric_id="${jobid%%.*}"

# -----------------------------------------------------------------------
# Step 1a — real-time status from jobinfo -c
# -----------------------------------------------------------------------
status=$(jobinfo -c 2>/dev/null | awk -v jid="$numeric_id" '$2 == jid { print $4 }')

case "$status" in
    Run|Queue|Hold)
        echo running
        exit 0
        ;;
    Exit)
        # Job exited — fall through to joblog to read the exit code.
        ;;
    "")
        # Not found in real-time view.
        # ---------------------------------------------------------------
        # Step 1b — stale-cache fallback: catches newly submitted jobs
        # that haven't yet appeared in -c, and recently finished jobs
        # whose cache entry still shows "Run".
        # Any hit here means "running" — we cannot trust the stale status
        # for final-state detection.
        # ---------------------------------------------------------------
        stale=$(jobinfo 2>/dev/null | awk -v jid="$numeric_id" '$2 == jid { print $4 }')
        if [[ -n "$stale" ]]; then
            echo running
            exit 0
        fi
        # Not found in either view — fall through to joblog.
        ;;
esac

# -----------------------------------------------------------------------
# Step 2 — look up exit status in joblog
#
# joblog -d 8 covers 8 days (sufficient for the 168 h max walltime + buffer).
#
# joblog output format (|-delimited table):
#   +--------+------+
#   | Job ID | exit |
#   +--------+------+
#   |  96293 |    0 |
# -----------------------------------------------------------------------
exit_code=$(joblog -d 8 -o jobid,exit_status 2>/dev/null \
    | awk -v jid="$numeric_id" '
        /^\|/ {
            split($0, f, "|")
            id = f[2]; gsub(/[[:space:]]/, "", id)
            es = f[3]; gsub(/[[:space:]]/, "", es)
            if (id ~ /^[0-9]+$/ && id == jid) { print es; found = 1 }
        }
        END { if (!found) print "NOT_FOUND" }
    ')

case "$exit_code" in
    0)
        echo success
        ;;
    NOT_FOUND)
        # Job not in joblog yet (lag) — return "running" so Snakemake retries.
        echo running
        ;;
    *)
        # Explicit non-zero exit code recorded in joblog
        echo failed
        ;;
esac
