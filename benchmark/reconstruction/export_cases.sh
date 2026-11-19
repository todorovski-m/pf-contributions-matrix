#!/bin/zsh
# Export the cases of benchmark/reconstruction/cases/*.m as inputs for the port and the old method, with the procedure
# of tests/export_synthetic.sh: for each case, a fresh scratch copy of the pinned submodule under /tmp,
# octave --eval "my_runpf('<case>')", and checks of the exit code, of a newly written results/results.mat and of
# success == 1. These cases are kept out of tests/cases_synthetic/, so that the default case lists of the timing
# series and of the test suites do not change.
# Usage (from the root of the matrix repository): zsh benchmark/reconstruction/export_cases.sh <output dir>
setopt ERR_EXIT NO_UNSET PIPE_FAIL
fail() { print -u2 "export_cases.sh: ERROR: $*"; exit 1 }
[[ -f apfc.m && -d benchmark/reconstruction/cases ]] || fail "run this script from the root of the matrix repository"
(( $# == 1 )) || fail "usage: zsh benchmark/reconstruction/export_cases.sh <output dir>"
OUT=$1
SUB=external/iterative-method
[[ -f $SUB/my_runpf.m ]] || fail "submodule $SUB is missing; run: git submodule update --init"
[[ -z "$(git -C $SUB status --porcelain)" ]] || fail "submodule $SUB has local modifications"
mkdir -p $OUT
for f in benchmark/reconstruction/cases/*.m; do
  c=${f:t:r}
  SCRATCH=$(mktemp -d /tmp/export_reconstruction.$c.XXXXXX)
  cp -R $SUB/ $SCRATCH/
  cp $f $SCRATCH/$c.m
  rm -f $SCRATCH/results/results.mat
  touch $SCRATCH/.start_marker
  set +e
  ( cd $SCRATCH && octave --no-gui --eval "my_runpf('$c')" ) > $SCRATCH/octave.log 2>&1
  RC=$?
  set -e
  (( RC == 0 )) || { tail -20 $SCRATCH/octave.log >&2; fail "octave exited with code $RC for $c (log: $SCRATCH/octave.log)" }
  MAT=$SCRATCH/results/results.mat
  [[ -f $MAT && $MAT -nt $SCRATCH/.start_marker ]] || fail "$c: $MAT was not newly written"
  python3 -B benchmark/manifest.py check $MAT || fail "$c: success check failed"
  cp $MAT $OUT/$c.mat
  print "$c: $(grep -h 'converged in\|Converged' $SCRATCH/octave.log | head -1) -> $OUT/$c.mat"
done
