#!/bin/zsh
# Export the synthetic test cases (tests/cases_synthetic/*.m) as inputs for the old method,
# with the same procedure as benchmark/export_inputs.sh: for each case, a fresh scratch copy
# of the pinned submodule under /tmp, octave --eval "my_runpf('<case>')", and checks of the
# exit code, of a newly written results/results.mat and of success == 1.
#
# Usage (from the root of the matrix repository):
#     zsh tests/export_synthetic.sh <output dir>
# Writes <output dir>/<case>.mat. Nothing is written inside the repository.

setopt ERR_EXIT NO_UNSET PIPE_FAIL
fail() { print -u2 "export_synthetic.sh: ERROR: $*"; exit 1 }
[[ -f apfc.m && -d tests/cases_synthetic ]] || fail "run this script from the root of the matrix repository"
(( $# == 1 )) || fail "usage: zsh tests/export_synthetic.sh <output dir>"
OUT=$1
SUB=external/iterative-method
[[ -f $SUB/my_runpf.m ]] || fail "submodule $SUB is missing; run: git submodule update --init"
[[ -z "$(git -C $SUB status --porcelain)" ]] || fail "submodule $SUB has local modifications"
mkdir -p $OUT
for f in tests/cases_synthetic/*.m; do
  c=${f:t:r}
  SCRATCH=$(mktemp -d /tmp/export_synthetic.$c.XXXXXX)
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
