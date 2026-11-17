# Timing harness: matrix method (Python port) and iterative method (Python)

This folder measures the run time of the matrix method (the Python port `python/pfc_matrix`) and of the earlier
iterative method (`external/iterative-method`, Python), both in Python and from the same input. It contains no
results. Octave and Matlab timings are not part of this study.

## How to run the series

On a quiet machine (see below), from the root of the matrix repository, with the submodule initialised:

```sh
zsh benchmark/timing/export_large.sh          # once: inputs of the large cases into benchmark/inputs_large/
PYTHON=/Users/mirko/.pyenv/shims/python zsh benchmark/timing/run_timing.sh
```

`run_timing.sh`:
- exports the synthetic cases into a scratch folder (`tests/export_synthetic.sh`; Octave is used only for that export);
- runs `timing.py run` under `caffeinate -i` on the 14 small cases and on every large case found in `benchmark/inputs_large/`;
- writes `benchmark/timing/results/<timestamp>/` (git-ignored) and summarises it with `summarize.py`.

Every option of `timing.py run` can be added, for example `--cases case9,case118`, `--budget 30`, `--timeout 600` or `--no-memory`; `--smoke` labels a test run as "smoke, not for publication".

`python -B -m pytest -p no:cacheprovider benchmark/timing/tests` tests the machinery.

## Files

| file | content |
|---|---|
| `timing.py` | measurement, statistics, time-out, guard, memory measurement, environment record, the series (`run`) |
| `cases.py` | the case lists, the loading through `load_results`, the import of the iterative method, the QSF correction |
| `summarize.py` | per-cell statistics, phase fractions, the ratio of the medians, log-log slopes against `n_e` |
| `run_timing.sh` | the full series under `caffeinate -i` |
| `export_large.sh`, `raw_counts.m`, `scan.py` | export and scan of the large cases |

A results folder contains:
- `samples.csv`: every raw sample in long format (case, method, kind, checks, repetition, nanoseconds, phase);
- `environment.json`: the environment record, with the per-case load averages, quiet flags and QSF corrections;
- `cells.csv`: the status of each cell;
- `guard.csv`: the agreement of the two methods before timing;
- `iterations.csv`: the iteration counts of the iterative method;
- `memory.csv`: the peak resident set size;
- `summary.md`, `summary.csv`, `summary_slopes.csv`: the outputs of `summarize.py`.

## Timed regions

Both methods start from the in-memory DataFrames `bus`, `gen` and `branch` of `load_results` (a solved power flow, exported by `my_runpf.m`). The power flow, `load_results`, Octave and any CSV writing are outside every timed region.

| method | inside the timed region | outside |
|---|---|---|
| port, APFC or RPFC | `read_branch_params` (reads BR_R, BR_X, BR_B and baseMVA from the exported file) and `Case.from_frames`, both included to be conservative against the port; then the public phase functions: (a) `prepare_*`, (b) `assemble_*` (with the checks), (c) `solve_*`, (d) `build_*_table`; with `checks` on, also `condest_A` | loading, writing the table |
| iterative method | `makeAPFC` or `makeRPFC(bus, gen, branch, None, itermax, False)`: their own preprocessing (`makeActiveFlows`, `makeReactiveFlows`, which copy the input), the iteration, and the result table | loading; the QSF correction; the iteration count (from an untimed run) |

- **Series:** the port is measured twice, with `checks=False` (the headline series, without `condest_A`) and with `checks=True`.
- **Phases:** the times of the port's phases are taken inside the same runs, with `time.perf_counter_ns` between the phase calls. Their sum is the run time except for the few calls between them.
- **Logging:** the iterative method's logging calls create no records (root logger at level WARNING, no handler), so it prints nothing. The text of its debug messages is still formatted inside its code, as written.
- **Iteration count:** it is read from an untimed run with an INFO collector. The iteration cap is `--itermax` (default 10000, as in `tests/compare_old.py`). A run that reaches the cap is "did not finish (iteration cap)"; for example the reactive part of `tap_loop`, which then has only the port timed.

## Statistics

- **Warm-up:** at least 2 warm-up runs per cell, not recorded.
- **Repetitions:** as many as fit in a budget of 30 s per cell, with at least 5 and at most 200.
- **Time-out:** a single run longer than 600 s ends the cell as "did not finish (time-out)".
- **Clock and garbage collector:** `time.perf_counter_ns`; `gc.collect()` before each cell, and the garbage collector disabled during a single run.
- **Reported per cell:** n, minimum, median, mean, standard deviation and interquartile range; every raw sample is kept.
- **Memory:** the peak resident set size of one run of each method (`resource.getrusage`, `ru_maxrss`), in a fresh subprocess per method and case, before and after the run. The unit is in the environment record (bytes on macOS, kibibytes on Linux).
- **Threads:** `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS` and `NUMEXPR_NUM_THREADS` are set to 1 before numpy is imported; the BLAS configuration of numpy is recorded.
- **Summary:** `summarize.py` gives the phase fractions of the port (median of a phase divided by the median of the run), the ratio old/port of the medians (port with `checks=False`), and, for each method and kind, the slope of a least-squares fit of log10(median) against log10(`n_e`), with the cases used.

## Factorisation facts

- The port factorises A once per method and run: one `scipy.sparse.linalg.splu` (SuperLU, COLAMD ordering), used for all right-hand sides. For RPFC, that one factorisation serves both the main solve and the self-term solve.
- `rpfc.m` under Matlab or Octave uses `\` twice, so it factorises twice.
- `condest_A` (series with `checks` on) makes its own factorisation.
- The iterative method has no matrix; its iteration counts are recorded per case.

## Guard before timing

Before a case is timed, both methods are run once in memory and compared, with the helper functions of `tests/compare_old.py` (unchanged). The case is timed only if:
- all rows map, and every cell agrees within 1e-10 when the iterative method finishes;
- the row sums equal the flows within 1e-10, for both methods;
- the port's phase composition gives the same table as `apfc`/`rpfc`.

Otherwise the case is flagged and not timed.

**Known issue of the iterative pipeline:** `my_runpf.m` computes QSF with B/2 instead of B/(2·tap²). Before the guard and the timing, QSF is therefore replaced in the branches with tap ≠ 1 and BR_B ≠ 0 by the Qsf that `rpfc.m` computes (`branch_shunt_flows` of the port). This happens outside the timed region, and `environment.json` lists the cases and rows.

**Cells not timed:** a cell whose dense result tables are estimated above half of the physical memory is skipped (`--mem-fraction`). So is a case that `load_results` cannot read, or whose input the port rejects. Each is recorded with its reason in `cells.csv`.

## Environment record and a quiet machine

`environment.json` records:
- the date and time, CPU model, cores, memory and OS version;
- the Python, numpy, scipy, pandas and networkx versions, the thread settings and the BLAS configuration;
- the git commits of `matrix` and of the submodule (and whether `matrix` has local changes);
- the power source and battery state (`pmset -g batt`), whether `caffeinate` is in use;
- the 1-, 5- and 15-minute load averages before and after each case.

A case is `quiet` if the 1-minute load average before it is below 0.25 times the number of cores. For a series meant for publication, use a machine on mains power with no other work running, and check that every case is quiet.

## Large cases

`export_large.sh` exports case300, case1354pegase, case2869pegase, case9241pegase, case_ACTIVSg2000 and case_ACTIVSg10k (MATPOWER 7.1) with the pipeline of `benchmark/export_inputs.sh`. A case that does not converge is dropped, not tuned. `scan.py` tabulates their properties. Not every exported case can be run by both methods:
- the port, like `apfc.m`, needs buses numbered 1..nb;
- `load_results` reads only the 13 standard bus columns;
- the methods' guards and the memory of the dense tables also apply.

`cells.csv` gives the reason for every case that is not timed.
