#!/bin/zsh
# Run an Octave test script from tests/unit/ in a fresh scratch copy under /tmp (with the stand-ins
# of tests/octave_shims), so that nothing is written inside the repository.
# Usage (from the root of the matrix repository): zsh tests/octave_in_scratch.sh <script name without .m>
# The exit status is that of Octave; the scripts print one "PASS <name>" or "FAIL <name>: ..." line per test.
setopt ERR_EXIT NO_UNSET PIPE_FAIL
fail() { print -u2 "octave_in_scratch.sh: ERROR: $*"; exit 2 }
[[ -f apfc.m && -d tests/unit ]] || fail "run this script from the root of the matrix repository"
(( $# == 1 )) || fail "usage: zsh tests/octave_in_scratch.sh <script>"
SCRATCH=$(mktemp -d /tmp/octave_in_scratch.XXXXXX)
cp apfc.m rpfc.m pfc_checks.m $SCRATCH/
mkdir -p $SCRATCH/cases $SCRATCH/results $SCRATCH/unit
cp cases/*.m $SCRATCH/cases/
for f in tests/cases_synthetic/*.m(N); do cp $f $SCRATCH/cases/; done
cp tests/unit/*.m $SCRATCH/unit/
cp tests/data/*.csv $SCRATCH/unit/
SHIMS=$PWD/tests/octave_shims
cd $SCRATCH && octave --no-gui --eval "addpath('$SHIMS'); addpath('$SCRATCH/unit'); $1"
