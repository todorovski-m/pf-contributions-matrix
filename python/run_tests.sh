#!/bin/zsh
# Run the tests of the Python port (matrix/python/tests) with pytest.
#
# Usage (from anywhere; the submodule external/iterative-method must be initialised):
#     zsh python/run_tests.sh            (from the root of the matrix repository)
# The Python interpreter is $PYTHON if set, else "python" (for this project: /Users/mirko/.pyenv/shims/python).
#
# If Octave is available, the synthetic exports (tests/export_synthetic.sh) and the Octave outputs of apfc.m and
# rpfc.m (tests/run_matrix.sh) are created in a scratch folder under /tmp and passed to the tests through
# PFC_SYN_INPUTS and PFC_OCTAVE_OUT; then every test runs. Without Octave, the tests that need these files are
# skipped, and their number is printed. Nothing is written inside the repository.
# Exit status: that of pytest (0 if every test passed).

setopt NO_UNSET PIPE_FAIL
ROOT=${0:A:h:h}
cd $ROOT || exit 2
[[ -f apfc.m && -d python/tests ]] || { print -u2 "run_tests.sh: cannot find the matrix repository at $ROOT"; exit 2 }
PY=${PYTHON:-python}
WORK=$(mktemp -d /tmp/python_run_tests.XXXXXX)
print "work folder: $WORK"
print "python: $(command -v $PY) ($($PY -c 'import sys; print(sys.version.split()[0])'))"

if command -v octave > /dev/null; then
  print "octave: $(octave --version | head -1)"
  zsh tests/export_synthetic.sh $WORK/syn_inputs > $WORK/export_synthetic.log 2>&1 \
    || { print -u2 "export of the synthetic inputs failed, see $WORK/export_synthetic.log"; exit 1 }
  zsh tests/run_matrix.sh $WORK/octave > $WORK/run_matrix.log 2>&1 \
    || { print -u2 "apfc.m and rpfc.m under Octave failed, see $WORK/run_matrix.log"; exit 1 }
  export PFC_SYN_INPUTS=$WORK/syn_inputs PFC_OCTAVE_OUT=$WORK/octave
  print "synthetic exports in $PFC_SYN_INPUTS, Octave outputs in $PFC_OCTAVE_OUT"
else
  print "octave: not available; the tests that need the synthetic exports or the Octave outputs are skipped"
  unset PFC_SYN_INPUTS PFC_OCTAVE_OUT
fi

( cd python && $PY -B -m pytest -p no:cacheprovider -rs tests ) 2>&1 | tee $WORK/pytest.log
RC=${pipestatus[1]}
NSKIP=$(grep -oE '[0-9]+ skipped' $WORK/pytest.log | tail -1 | grep -oE '[0-9]+'); NSKIP=${NSKIP:-0}
print "\nskipped tests: $NSKIP$([[ -z ${PFC_OCTAVE_OUT:-} ]] && print ' (Octave not available)')"
print "pytest exit status: $RC (log: $WORK/pytest.log)"
exit $RC
