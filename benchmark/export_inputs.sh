#!/bin/zsh
# Export the canonical benchmark inputs.
#
# For each case, run my_runpf.m of the pinned iterative method (external/iterative-method)
# in a fresh scratch copy under /tmp and copy the resulting results/results.mat to
# benchmark/inputs/<case>.mat. Also write benchmark/inputs/MANIFEST.json and
# benchmark/environment.txt.
#
# Usage (from the root of the matrix repository):
#     zsh benchmark/export_inputs.sh
# Optional: EXPORT_DIR=/some/dir zsh benchmark/export_inputs.sh writes the .mat files and
# MANIFEST.json there instead (used by the reproducibility check); environment.txt is then
# written into EXPORT_DIR as well.

setopt ERR_EXIT NO_UNSET PIPE_FAIL

fail() { print -u2 "export_inputs.sh: ERROR: $*"; exit 1 }

[[ -f apfc.m && -f rpfc.m && -d cases ]] || fail "run this script from the root of the matrix repository"
SUB=external/iterative-method
[[ -f $SUB/my_runpf.m && -f $SUB/power_flow.py ]] || fail "submodule $SUB is missing; run: git submodule update --init"
[[ -z "$(git -C $SUB status --porcelain)" ]] || fail "submodule $SUB has local modifications"
command -v octave >/dev/null || fail "octave not found"
command -v python3 >/dev/null || fail "python3 not found"

OUT=${EXPORT_DIR:-benchmark/inputs}
ENVFILE=${EXPORT_DIR:+$EXPORT_DIR/environment.txt}
ENVFILE=${ENVFILE:-benchmark/environment.txt}
mkdir -p $OUT

CASES=(case1 case9 case14 case24_ieee_rts case30 case33 case118 tap_loop)
SHARED=(case1 case9 case14 case30 case33 case118)

SUB_COMMIT=$(git -C $SUB rev-parse HEAD)
MATRIX_COMMIT=$(git rev-parse HEAD)
MATRIX_DIRTY=0; [[ -z "$(git status --porcelain)" ]] || MATRIX_DIRTY=1

# The six shared case files must be byte-identical in both repositories.
for c in $SHARED; do
  cmp -s $SUB/$c.m cases/$c.m || fail "$SUB/$c.m and cases/$c.m differ"
done

PAIRS=()
for c in $CASES; do
  SCRATCH=$(mktemp -d /tmp/export_inputs.$c.XXXXXX)
  cp -R $SUB/ $SCRATCH/
  if [[ $c == tap_loop ]]; then
    cp cases/tap_loop.m $SCRATCH/tap_loop.m
    CASE_FILE=cases/tap_loop.m
  else
    [[ -f $SUB/$c.m ]] || fail "case file $SUB/$c.m not found"
    CASE_FILE=$SUB/$c.m
  fi
  # remove the results.mat copied from the submodule, so that a new file must be written
  rm -f $SCRATCH/results/results.mat
  touch $SCRATCH/.start_marker
  print "== $c (scratch $SCRATCH)"
  set +e
  ( cd $SCRATCH && octave --no-gui --eval "my_runpf('$c')" ) > $SCRATCH/octave.log 2>&1
  RC=$?
  set -e
  (( RC == 0 )) || { tail -20 $SCRATCH/octave.log >&2; fail "octave exited with code $RC for $c (log: $SCRATCH/octave.log)" }
  MAT=$SCRATCH/results/results.mat
  [[ -f $MAT ]] || fail "$c: $MAT was not written"
  [[ $MAT -nt $SCRATCH/.start_marker ]] || fail "$c: $MAT is not newer than the start of the run"
  python3 benchmark/manifest.py check $MAT || fail "$c: success check failed"
  cp $MAT $OUT/$c.mat
  PAIRS+=("$c=$CASE_FILE")
done

python3 benchmark/manifest.py manifest $OUT $SUB_COMMIT $MATRIX_COMMIT $MATRIX_DIRTY $PAIRS

# environment record
{
  print "# Environment of the export (written by benchmark/export_inputs.sh)"
  print "date_utc: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  print "octave: $(octave --version 2>&1 | head -1)"
  print "matpower: $(octave --no-gui --eval 'mpver' 2>&1 | grep '^MATPOWER ' | tr -s ' ')"
  print "python3: $(command -v python3) $(python3 --version 2>&1)"
  python3 -c "
import importlib
for p in ['numpy', 'pandas', 'scipy', 'networkx', 'matplotlib']:
    print(f'{p}: {importlib.import_module(p).__version__}')"
  print "macos: $(sw_vers -productVersion) (build $(sw_vers -buildVersion))"
  print "cpu: $(sysctl -n machdep.cpu.brand_string)"
  print "cores (hw.ncpu / physical / logical): $(sysctl -n hw.ncpu) / $(sysctl -n hw.physicalcpu) / $(sysctl -n hw.logicalcpu)"
  print "iterative-method submodule commit: $SUB_COMMIT"
  print "matrix commit: $MATRIX_COMMIT (working tree dirty: $MATRIX_DIRTY)"
} > $ENVFILE

cat > $OUT/README.md <<'EOF'
# Canonical benchmark inputs

This folder holds one input file per test case, shared by the matrix method (this
repository) and the iterative method (the submodule `external/iterative-method`, pinned
at tag `benchmark-v1`).

## What the files are

- `<case>.mat`: the file `results/results.mat` written by `my_runpf.m` of the iterative
  method for that case. It contains the solved MATPOWER case (`bus`, `gen`, `branch`, ...)
  after the preprocessing of `my_runpf.m` (sign of `PT` and `QT` reversed so that they are
  the flows delivered at the to-end, branches with all four flows below 1e-4 removed,
  single flows below 1e-4 set to zero), and the derived quantities `PS`, `QS`, `QSF` and
  `QST`. The function `load_results` in `external/iterative-method/power_flow.py` turns a
  file into the DataFrames `bus`, `gen` and `branch` without calling Octave.
- `MANIFEST.json`: for each case, the SHA-256 of the case file used and of the exported
  `.mat` file, a content hash that ignores `et`, the number of buses, branches (after the
  deletion in `my_runpf.m`) and generators, the commits of the submodule and of this
  repository, and the date.
- `../environment.txt`: Octave, MATPOWER, Python and library versions, macOS version,
  CPU model and number of cores of the machine that produced the files.

## How they are produced

From the root of this repository:

    git submodule update --init
    zsh benchmark/export_inputs.sh

For each case the script copies the submodule into a fresh scratch directory under
`/tmp`, runs `octave --eval "my_runpf('<case>')"` there, stops with an error if Octave
fails, if `results/results.mat` is not newly written or if `success` is not 1, and copies
the file here. The case files of the six shared cases (case1, case9, case14, case30,
case33, case118) are taken from the submodule after checking that they are
byte-identical to `cases/`; `case24_ieee_rts` exists only in the submodule; `tap_loop` is
taken from `cases/tap_loop.m`.

## Variation between runs

The variable `et` (MATPOWER's elapsed time) varies between runs and must be ignored in
any comparison. The file hash `mat_sha256` also changes on every export, because Octave
writes the time of writing into the MAT-file header. Compare files with
`content_sha256_excluding_et` or with `benchmark/check_inputs.py`, which compares every
array except `et`.
EOF

print "Exported ${#CASES} cases to $OUT; environment written to $ENVFILE"
