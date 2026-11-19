#!/bin/zsh
# Run apfc.m and rpfc.m of this working tree under Octave (with tests/octave_shims) on the cases of
# benchmark/reconstruction/cases/, in a fresh scratch copy under /tmp, as tests/run_matrix.sh does for the cases of
# cases/ and tests/cases_synthetic/. Usage (from the matrix root): zsh benchmark/reconstruction/run_octave_cases.sh <out dir> [case ...]
setopt ERR_EXIT NO_UNSET PIPE_FAIL
fail() { print -u2 "run_octave_cases.sh: ERROR: $*"; exit 1 }
[[ -f apfc.m && -f rpfc.m && -d benchmark/reconstruction/cases ]] || fail "run this script from the root of the matrix repository"
(( $# >= 1 )) || fail "usage: zsh benchmark/reconstruction/run_octave_cases.sh <out dir> [case ...]"
OUT=$1; shift
CASES=("$@")
(( ${#CASES} )) || for f in benchmark/reconstruction/cases/*.m; do CASES+=(${f:t:r}); done
mkdir -p $OUT
SCRATCH=$(mktemp -d /tmp/run_octave_cases.XXXXXX)
cp apfc.m rpfc.m pfc_checks.m $SCRATCH/
mkdir -p $SCRATCH/cases $SCRATCH/results
cp benchmark/reconstruction/cases/*.m $SCRATCH/cases/
SHIMS=$PWD/tests/octave_shims
LIST="{$(printf "'%s'," $CASES | sed 's/,$//')}"
set +e
( cd $SCRATCH && octave --no-gui --eval "addpath('$SHIMS'); c = $LIST; for i = 1:numel(c), apfc(['cases/' c{i}]); rpfc(['cases/' c{i}]); end" ) > $SCRATCH/octave.log 2>&1
RC=$?
set -e
(( RC == 0 )) || { tail -20 $SCRATCH/octave.log >&2; fail "octave exited with code $RC (log: $SCRATCH/octave.log)" }
for c in $CASES; do for k in apfc rpfc; do
  [[ -f $SCRATCH/results/${c}_$k.csv ]] || fail "missing output $SCRATCH/results/${c}_$k.csv"
  cp $SCRATCH/results/${c}_$k.csv $OUT/
done; done
print "wrote ${#CASES} cases (apfc and rpfc) to $OUT; warnings in $SCRATCH/octave.log: $(grep -c -i warning $SCRATCH/octave.log)"
