#!/bin/zsh
# Run the timing series of the matrix method (Python port) against the old iterative method (Python) under
# caffeinate -i, and write benchmark/timing/results/<timestamp>/ (raw samples, environment record, memory, iteration
# counts, guard results, summary). See benchmark/timing/README.md.
#
# Usage (from the root of the matrix repository):
#     zsh benchmark/timing/run_timing.sh [--smoke] [any option of timing.py run, for example --cases case9,case118]
# Defaults: the 14 small cases and every large case found in benchmark/inputs_large/ (made by export_large.sh);
# budget 30 s per cell, 5 to 200 repetitions, 2 warm-up runs, time-out 600 s per single run of the old method.
# The synthetic cases are exported first with tests/export_synthetic.sh into a scratch folder (Octave is needed for
# that export only; it is outside every timed region). $PYTHON selects the interpreter (default: python).

setopt ERR_EXIT NO_UNSET PIPE_FAIL
ROOT=${0:A:h:h:h}
cd $ROOT
[[ -f apfc.m && -d benchmark/timing ]] || { print -u2 "run_timing.sh: cannot find the matrix repository"; exit 2 }
PY=${PYTHON:-python}
STAMP=$(date -u +%Y%m%dT%H%M%SZ)
OUT=benchmark/timing/results/$STAMP
mkdir -p $OUT
SYN=$(mktemp -d /tmp/run_timing_syn.XXXXXX)
zsh tests/export_synthetic.sh $SYN > $OUT/export_synthetic.log 2>&1 || { print -u2 "synthetic export failed, see $OUT/export_synthetic.log"; exit 1 }

CASES=(case1 case9 case14 case24_ieee_rts case30 case33 case118 tap_loop syn9_neg_pd syn9_neg_pg syn9_neg_qd syn9_neg_qg syn9_case3_tap syn9_carrier_tap)
for c in case300 case1354pegase case2869pegase case9241pegase case_ACTIVSg2000 case_ACTIVSg10k; do
  [[ -f benchmark/inputs_large/$c.mat ]] && CASES+=($c)
done
LIST=${(j:,:)CASES}
HAS_CASES=0
for a in "$@"; do [[ $a == --cases ]] && HAS_CASES=1; done
EXTRA=()
(( HAS_CASES )) || EXTRA=(--cases $LIST)

print "results: $OUT"
print "git: $(git rev-parse HEAD)$([[ -n "$(git status --porcelain)" ]] && print ' (dirty)')"
export PFC_CAFFEINATE=1
caffeinate -i $PY -B benchmark/timing/timing.py run --out $OUT --syn-inputs $SYN $EXTRA "$@"
$PY -B benchmark/timing/summarize.py $OUT
print "done: $OUT/summary.md"
