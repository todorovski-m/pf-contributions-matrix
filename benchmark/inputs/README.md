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
