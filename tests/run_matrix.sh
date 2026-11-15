#!/bin/zsh
# Run apfc.m and rpfc.m of this working tree under Octave (with tests/octave_shims) in a
# fresh scratch copy under /tmp, and collect the CSV outputs.
#
# Usage (from the root of the matrix repository):
#     zsh tests/run_matrix.sh <output dir> [case ...]
# Cases are looked up in cases/ and tests/cases_synthetic/. Default case list: the six
# shared cases, case24_ieee_rts, tap_loop and all synthetic cases.

setopt ERR_EXIT NO_UNSET PIPE_FAIL
fail() { print -u2 "run_matrix.sh: ERROR: $*"; exit 1 }
[[ -f apfc.m && -f rpfc.m && -d cases ]] || fail "run this script from the root of the matrix repository"
(( $# >= 1 )) || fail "usage: zsh tests/run_matrix.sh <output dir> [case ...]"
OUT=$1; shift
CASES=("$@")
if (( ${#CASES} == 0 )); then
  CASES=(case1 case9 case14 case24_ieee_rts case30 case33 case118 tap_loop)
  for f in tests/cases_synthetic/*.m(N); do CASES+=(${f:t:r}); done
fi
mkdir -p $OUT
SCRATCH=$(mktemp -d /tmp/run_matrix.XXXXXX)
cp apfc.m rpfc.m pfc_checks.m $SCRATCH/
mkdir -p $SCRATCH/cases $SCRATCH/results
cp cases/*.m $SCRATCH/cases/
for f in tests/cases_synthetic/*.m(N); do cp $f $SCRATCH/cases/; done
SHIMS=$PWD/tests/octave_shims
LIST="{$(printf "'%s'," $CASES | sed 's/,$//')}"
print "scratch: $SCRATCH"
set +e
( cd $SCRATCH && octave --no-gui --eval "addpath('$SHIMS'); c = $LIST; for i = 1:numel(c), apfc(['cases/' c{i}]); rpfc(['cases/' c{i}]); end" ) > $SCRATCH/octave.log 2>&1
RC=$?
set -e
(( RC == 0 )) || { tail -20 $SCRATCH/octave.log >&2; fail "octave exited with code $RC (log: $SCRATCH/octave.log)" }
for c in $CASES; do
  for k in apfc rpfc; do
    [[ -f $SCRATCH/results/${c}_$k.csv ]] || fail "missing output $SCRATCH/results/${c}_$k.csv"
    cp $SCRATCH/results/${c}_$k.csv $OUT/
  done
done
print "wrote ${#CASES} cases (apfc and rpfc) to $OUT"
