# Scaling with the number of sources

`sources_scaling.py` measures how the time of the matrix method (Python port) changes with the number k of sources,
and what one assembly and one factorization with k right-hand sides save against k separate single-source computations
(brief T14g, part A; reviewer comments R1.4 and R3). It uses the machinery of `benchmark/timing/` (measurement with
warm-up and repetitions, median and interquartile range, the quiet wait, the slow-cell rule, the environment record),
imported and not changed.

## Definitions

- **Source:** a column of the right-hand side B, which is also a column of the output table. APFC: the generators with
  positive output and the buses with negative demand. RPFC: the generators, the injecting bus shunts, the buses with
  negative demand and the terminal injections of the Case-3 branches. The self-term columns `E_i_j_s` of RPFC are not
  columns of B (they come from a separate solve) and are not counted.
- **Subsets:** the first k entries of one fixed permutation of the source columns (numpy `default_rng(20261119)`), so
  they are deterministic and nested; k = 1, 2, 4, 8, ... below the number of sources ns, and ns.
- **"one run":** one complete computation for k sources: input (`Case.from_frames`, `read_branch_params`), prepare,
  assemble, the solve of the k selected columns of B (the port's `solve_active` / `solve_reactive` on the system with B
  restricted to those columns; in RPFC without the self-term solve), and the table block of the k columns.
- **"repeated":** k such computations with one column each. Measured up to k = 16; above, it is reported as k times the
  measured single-column time and marked `extrapolated` (in `sources.csv` and with a cross in the figures).
- **Table block:** the public builders `build_active_table` and `build_reactive_table` always label every source, so the
  block of the k columns is formed with their formula (demand rows, from-end rows X, to-end rows a·X, and in RPFC the
  Case-3 rows of the selected Case-3 columns). The method code is not changed.

## Correctness before timing

```sh
python -B benchmark/sources_scaling/sources_scaling.py check --cases case9,case118 --k 1,2,4
```

For each case, kind and k, the k-column block must equal the stack of the k single-column blocks and the same columns
of the public table (`apfc`, `rpfc`), within 1e-12.

## The run (on a quiet machine, mains power, Spotlight and sync services paused)

From the root of the matrix repository:

```sh
python -B benchmark/sources_scaling/sources_scaling.py estimate
caffeinate -i /Users/mirko/.pyenv/shims/python -B benchmark/sources_scaling/sources_scaling.py run --out benchmark/sources_scaling/results/$(date -u +%Y%m%dT%H%M%SZ)
```

Defaults: the networks case118, case_ACTIVSg500 and case_ACTIVSg2000 (the large ones from `benchmark/inputs_large/`,
made by `benchmark/timing/export_large.sh`), both kinds, a budget of 20 s per cell with 5 to 200 repetitions and 2
warm-up runs, the slow-cell rule of the timing harness, and the quiet wait (at most 300 s per cell). `estimate` gives the
number of cells and the lower bound of the duration from these settings (80 cells, at least 1600 s); the parts that
depend on run times are named, not estimated.

**Outputs** (in the folder given by `--out`; `results/` is ignored by git): `sources.csv` (network, kind, k, mode, n,
median, min and IQR in nanoseconds, quiet, one-minute load, extrapolated), `summary.md` (with the ratio repeated / one
run), `environment.json`, and the figures `sources_time.pdf` (time against k, both modes, log-log) and
`sources_ratio.pdf` (repeated / one run against k). `--smoke` labels a test run "smoke, not for publication".

## Tests

```sh
python -B -m pytest -p no:cacheprovider benchmark/sources_scaling/tests
```

They check the subset selection, the k values, the correctness of the k-column computation on case9, the extrapolation
labels and the CSV schema (with a tiny smoke run in a temporary folder); no time is checked.
