#!/bin/zsh
# Export the large cases (MATPOWER 7.1 case files on the Octave path) as inputs for both methods, with the pipeline of
# benchmark/export_inputs.sh: for each case, a fresh scratch copy of the pinned submodule under /tmp, then
#   1. prepare_case.m (in the scratch copy): the buses are renumbered 1..nb in table order and the matrices are
#      stripped to the standard MATPOWER columns (bus 13, gen 21, branch 13); the mapping old -> new bus number is
#      saved as benchmark/inputs_large/<case>_busmap.csv. No MATPOWER file, my_runpf.m or submodule file is changed;
#   2. neutrality_check.m: the power flow of the prepared case, mapped back, against that of the original case
#      (benchmark/inputs_large/neutrality.csv);
#   3. octave --eval "my_runpf('prepared_<case>.mat')" (my_runpf.m unchanged), with checks of the exit code, of a
#      newly written results/results.mat and of success == 1.
# A case that does not converge is dropped and recorded; it is not tuned.
#
# Usage (from the root of the matrix repository):
#     zsh benchmark/timing/export_large.sh [case ...]
# Default cases (the ladder, then the stretch cases): case300 case_ACTIVSg500 case1354pegase case_ACTIVSg2000
# case2383wp case2869pegase case3120sp, case_ACTIVSg10k case9241pegase.
# Writes benchmark/inputs_large/<case>.mat (git-ignored; regenerate with this script) and
# benchmark/inputs_large/export_log.txt (one line per case: exported, or dropped with the reason).

setopt ERR_EXIT NO_UNSET PIPE_FAIL
fail() { print -u2 "export_large.sh: ERROR: $*"; exit 1 }
[[ -f apfc.m && -d benchmark ]] || fail "run this script from the root of the matrix repository"
SUB=external/iterative-method
[[ -f $SUB/my_runpf.m ]] || fail "submodule $SUB is missing; run: git submodule update --init"
[[ -z "$(git -C $SUB status --porcelain)" ]] || fail "submodule $SUB has local modifications"
command -v octave > /dev/null || fail "octave not found"
OUT=benchmark/inputs_large
mkdir -p $OUT
CASES=("$@")
(( ${#CASES} )) || CASES=(case300 case_ACTIVSg500 case1354pegase case_ACTIVSg2000 case2383wp case2869pegase case3120sp case_ACTIVSg10k case9241pegase)
TIMING=$PWD/benchmark/timing
NEU=$PWD/$OUT/neutrality.csv          # absolute: used inside the scratch folder
print "case,max_abs_vm,max_abs_va_deg,max_abs_pf,max_abs_qf,max_abs_pt,max_abs_qt,iterations_original,iterations_prepared,success_original,success_prepared,busmap_bijection,maps_back" > $NEU
LOG=$OUT/export_log.txt
print "# export_large.sh $(date -u +%Y-%m-%dT%H:%M:%SZ); submodule $(git -C $SUB rev-parse HEAD); $(octave --version | head -1)" > $LOG
for c in $CASES; do
  SCRATCH=$(mktemp -d /tmp/export_large.$c.XXXXXX)
  cp -R $SUB/ $SCRATCH/
  rm -f $SCRATCH/results/results.mat
  touch $SCRATCH/.start_marker
  set +e
  ( cd $SCRATCH && octave --no-gui --eval "addpath('$TIMING'); printf('case file: %s\n', which('$c')); prepare_case('$c', 'prepared_$c.mat', 'busmap.csv'); neutrality_check('$c', 'prepared_$c.mat', 'busmap.csv', '$NEU'); my_runpf('prepared_$c.mat')" ) > $SCRATCH/octave.log 2>&1
  RC=$?
  set -e
  MAT=$SCRATCH/results/results.mat
  if (( RC != 0 )) || [[ ! -f $MAT || ! $MAT -nt $SCRATCH/.start_marker ]]; then
    print "$c: dropped (octave exit $RC or no new results.mat; log $SCRATCH/octave.log)" | tee -a $LOG; continue
  fi
  if ! python3 -B benchmark/manifest.py check $MAT > /dev/null 2>&1; then
    print "$c: dropped (power flow did not converge: success != 1; log $SCRATCH/octave.log)" | tee -a $LOG; continue
  fi
  cp $MAT $OUT/$c.mat
  cp $SCRATCH/busmap.csv $OUT/${c}_busmap.csv
  print "$c: exported ($(grep -h 'case file:' $SCRATCH/octave.log); $(grep -h 'converged in\|Converged' $SCRATCH/octave.log | head -1))" | tee -a $LOG
done
