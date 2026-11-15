# Test data

`T0_D_bad_1.csv`: the cells of the reactive results in which the v1.0 code and the old iterative
method differ, from the T0 analysis (`reports/T0_scripts/out/D_bad_1.csv` of the revision work,
produced by `reports/T0_scripts/d_compare.py` with the Case-3 labels corrected). The rows with
case `case14` are the cells that depend on the negative reactive demand at bus 4.
`tests/compare_legacy.py` uses them as the set of expected differences against `paper_results/`.
