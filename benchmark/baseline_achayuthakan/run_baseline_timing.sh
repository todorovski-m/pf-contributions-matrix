#!/bin/zsh
# Runtime of the Achayuthakan et al. (2010) baseline against the port (brief T14h, step 4), under caffeinate -i, into
# benchmark/baseline_achayuthakan/results/<timestamp>/ (ignored by .gitignore). Run on a quiet machine (mains power, Spotlight and
# sync services paused), in the same session as the timing series. Extra options go to baseline_timing.py run
# (for example --cases case9,case118 or --budget 10).
# Usage (from anywhere): zsh benchmark/baseline_achayuthakan/run_baseline_timing.sh [options]
setopt ERR_EXIT NO_UNSET PIPE_FAIL
ROOT=${0:A:h:h:h}
cd $ROOT
PY=${PYTHON:-/Users/mirko/.pyenv/shims/python}
OUT=benchmark/baseline_achayuthakan/results/$(date -u +%Y%m%dT%H%M%SZ)
$PY -B benchmark/baseline_achayuthakan/baseline_timing.py estimate
print "results: $OUT"
print "git: $(git rev-parse HEAD)$([[ -n "$(git status --porcelain)" ]] && print ' (dirty)')"
export PFC_CAFFEINATE=1
caffeinate -i $PY -B benchmark/baseline_achayuthakan/baseline_timing.py run --out $OUT "$@"
print "done: $OUT/summary.md"
