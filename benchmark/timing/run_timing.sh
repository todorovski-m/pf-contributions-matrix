#!/bin/zsh
# Run the timing series of the matrix method (Python port) against the old iterative method (Python) under
# caffeinate -i, and write benchmark/timing/results/<timestamp>/ (raw samples, environment record, memory, iteration
# counts, guard results, summary). See benchmark/timing/README.md.
#
# Usage (from the root of the matrix repository):
#     zsh benchmark/timing/run_timing.sh [--smoke] [any option of timing.py run, for example --cases case9,case118]
# Defaults: the (case, kind) combinations of cases.py (COMBOS: the 14 small cases with both kinds, and the large
# combinations that are usable, if benchmark/inputs_large/<case>.mat exists; made by export_large.sh); budget 30 s
# per cell, 5 to 200 repetitions, 2 warm-up runs; the slow-cell rule (a single run above 30 s: 3 repetitions, flagged
# slow); time-out 1800 s per single run (the old method of the large cases); the stall probe for predicted stalls
# (itermax 300, time-out 120 s for the guard run; --no-stall-probe restores the full run).
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

COMBOS=$($PY -B -c "import os, sys; sys.path.insert(0, 'benchmark/timing'); import cases as C; print(C.combos_text([(c, k) for c, k in C.COMBOS if c in C.SMALL or os.path.exists(os.path.join('benchmark', 'inputs_large', c + '.mat'))]))")
HAS_CASES=0; HAS_TIMEOUT=0; HAS_PROBE=0
for a in "$@"; do [[ $a == --cases || $a == --combos ]] && HAS_CASES=1; [[ $a == --timeout ]] && HAS_TIMEOUT=1; [[ $a == --stall-probe || $a == --no-stall-probe ]] && HAS_PROBE=1; done
EXTRA=()
(( HAS_CASES )) || EXTRA+=(--combos $COMBOS)
(( HAS_TIMEOUT )) || EXTRA+=(--timeout 1800)
(( HAS_PROBE )) || EXTRA+=(--stall-probe)          # T14f: default on; --no-stall-probe restores the full guard run

print "results: $OUT"
print "git: $(git rev-parse HEAD)$([[ -n "$(git status --porcelain)" ]] && print ' (dirty)')"
export PFC_CAFFEINATE=1
caffeinate -i $PY -B benchmark/timing/timing.py run --out $OUT --syn-inputs $SYN $EXTRA "$@"
$PY -B benchmark/timing/summarize.py $OUT
$PY -B benchmark/timing/plot_timing.py $OUT
print "done: $OUT/summary.md"
