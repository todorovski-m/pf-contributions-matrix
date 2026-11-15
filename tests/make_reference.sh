#!/bin/zsh
# Write reference_results/: the outputs of the current apfc.m and rpfc.m, run under Octave (with
# tests/octave_shims) on the six shared cases and tap_loop, with 15 significant digits as Matlab's
# writetable writes them (tests/format15.py).
#
# Usage (from the root of the matrix repository): zsh tests/make_reference.sh [output dir]
# Default output dir: reference_results. The run itself happens in a scratch folder under /tmp.
setopt ERR_EXIT NO_UNSET PIPE_FAIL
[[ -f apfc.m && -f rpfc.m ]] || { print -u2 "run this script from the root of the matrix repository"; exit 2 }
OUT=${1:-reference_results}
CASES=(case1 case9 case14 case30 case33 case118 tap_loop)
RAW=$(mktemp -d /tmp/make_reference.XXXXXX)
zsh tests/run_matrix.sh $RAW $CASES
python3 tests/format15.py $RAW $OUT
print "reference results written to $OUT (raw %.17g outputs in $RAW)"
