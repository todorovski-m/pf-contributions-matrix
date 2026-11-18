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
| `summarize.py` | per-cell statistics, phase fractions, the ratio of the medians, table sizes, local and large-end log-log slopes against `n_e` |
| `plot_timing.py` | figures of a results folder: time against `n_e` (log-log) and the ratio old/port against `n_e`, PDF and PNG |
| `old_stall_predictor.py` | static prediction of where the old method stalls (no old-method run) |
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
- **Summary:** `summarize.py` gives per cell the phase fractions of the port (median of a phase divided by the median of the run), the ratio old/port of the medians (port with `checks=False`; for every cell where both finished, cells flagged slow included), the dense table size in bytes and the number of table columns. For each method, checks setting and kind it gives the slopes of log10(median) against log10(`n_e`):
  - **local slopes** between consecutive cases with `n_e` >= 500, each with its two cases (`summary_local_slopes.csv`);
  - **one least-squares slope over the cases with `n_e` >= 500**, with the number of points;
  - **the earlier global slope over all cases**, kept in its own column labelled "all cases, dominated by fixed overhead": the small cases sit at a fixed overhead of a few milliseconds and do not describe the large end.

  A median whose cell is not quiet is marked. Local slopes between cases of nearly equal `n_e` reflect differences in structure as much as size and can even be negative; they are not growth rates.
- **Figures:** `plot_timing.py` draws time against `n_e` (log-log, checks off, both methods and kinds) and the ratio old/port against `n_e`. A cell where the old method did not finish is drawn at the run's time-out with one marker per failure kind (time-out, iteration cap, error); that position is a limit, not a measured time. `run_timing.sh` calls it after `summarize.py`.

## Factorisation facts

- The port factorises A once per method and run: one `scipy.sparse.linalg.splu` (SuperLU, COLAMD ordering), used for all right-hand sides. For RPFC, that one factorisation serves both the main solve and the self-term solve.
- `rpfc.m` under Matlab or Octave uses `\` twice, so it factorises twice.
- `condest_A` (series with `checks` on) makes its own factorisation.
- The iterative method has no matrix; its iteration counts are recorded per case.

## Guard before timing

Before a case is timed, both methods are run once in memory, with the row mapping of `correctness.py` (helper functions of `tests/compare_old.py`, unchanged).
- **Always (the port's own conditions):** the port's phase composition gives the same table as `apfc`/`rpfc`, and its row sums equal the flows within 1e-10.
- **When the old method finishes, also:** every port row maps, every mapped cell agrees within 1e-10, and the port row sums equal the old flows within 1e-10.
- **When the old method does not finish** (iteration cap, time-out or error; T14e): only the port's own conditions are checked, and the old status is recorded. An unfinished old table is not a reference, so its cells, its row errors and the port rows that cannot be mapped to it are not compared. `guard.csv` records which checks were applied (`guard_checks`).

A case that fails the guard is flagged and not timed. Before T14e, `case1354pegase` RPFC failed it on the unfinished old table: three port rows had no counterpart, and inductive branch-end shunt rows `Tf_<bus>` of buses with several such rows were paired by order of occurrence, which differs between the methods.

**Known issue of the iterative pipeline:** `my_runpf.m` computes QSF with B/2 instead of B/(2·tap²). Before the guard and the timing, QSF is therefore replaced in the branches with tap ≠ 1 and BR_B ≠ 0 by the Qsf that `rpfc.m` computes (`branch_shunt_flows` of the port). This happens outside the timed region, and `environment.json` lists the cases and rows.

**Cells not timed:** a cell whose dense result tables are estimated above half of the physical memory is skipped (`--mem-fraction`). So is a case that `load_results` cannot read, or whose input the port rejects. Each is recorded with its reason in `cells.csv`.

## Environment record and a quiet machine

`environment.json` records:
- the date and time, CPU model, cores, memory and OS version;
- the Python, numpy, scipy, pandas and networkx versions, the thread settings and the BLAS configuration;
- the git commits of `matrix` and of the submodule (and whether `matrix` has local changes);
- the power source and battery state (`pmset -g batt`), whether `caffeinate` is in use;
- the 1-, 5- and 15-minute load averages before and after each case.

**Quiet wait (T14e).** Before each cell, outside every timed region, the harness reads the one-minute load average every 10 s (`--quiet-poll`) and waits until it is below 0.25 times the number of logical cores, for at most 300 s (`--quiet-wait`).
- **Recorded per cell** (`cells.csv`): the wait time, the load at the start of the wait and the load when the cell starts.
- **When the wait expires:** the cell is still timed but flagged `quiet` False, and `environment.json` records the three processes with the highest CPU share at that moment (`ps`).
- **Per case:** a case is `quiet` unless the wait of one of its cells expired.
- **Why the threshold leaves room for the benchmark itself:** the benchmark runs in one thread, which by itself adds about 1 to the one-minute load average while it runs. On the 8 logical cores of the reference machine the threshold is 2.0, so the benchmark's own thread and the decay of its previous cell do not keep the machine from counting as quiet, while other work using more than about one core does.

For a series meant for publication, use a machine on mains power with no other work running, and check that every cell is quiet.

## Large cases

### Export: renumbering and stripping

`export_large.sh` exports case300, case_ACTIVSg500, case1354pegase, case_ACTIVSg2000, case2383wp, case2869pegase,
case3120sp, case_ACTIVSg10k and case9241pegase (MATPOWER 7.1) into `benchmark/inputs_large/<case>.mat` (git-ignored),
with the pipeline of `benchmark/export_inputs.sh`. A case that does not converge is dropped, not tuned; `export_log.txt`
records each case. Before the power flow, `prepare_case.m` changes two things in a scratch copy of the case:
- **Renumbering:** the buses are renumbered 1..nb in the order of the bus table, and the bus numbers of the generators,
  the branches and the area reference buses are mapped in the same way. The port, like `apfc.m`, needs buses numbered
  1..nb; several large cases have other numbers. The mapping (old bus, new bus) is saved as `<case>_busmap.csv`, so
  results can be mapped back.
- **Stripping:** the bus, generator and branch matrices are cut to the standard MATPOWER columns (13, 21 and 13).
  `load_results` reads only the standard columns, and the extra columns of some cases (for example market data) are
  not used by the power flow.

No MATPOWER file, no file of the submodule and `my_runpf.m` are changed.

**Neutrality check.** `neutrality_check.m` runs the power flow of the original case and of the prepared case with the
same options, maps the prepared one back, and writes to `benchmark/inputs_large/neutrality.csv` the largest absolute
differences of VM, VA, PF, QF, PT and QT, the iteration counts, the success flags, and whether the mapping is a
bijection onto 1..nb that maps every bus, generator and branch back. A prepared case is used only if all differences
are zero, both power flows converge in the same number of iterations, and both flags are 1. The test
`test_renumbering_of_the_export` (in `tests/test_timing.py`) checks the same on case9 with every bus number increased
by 100.

### Case and kind list

`cases.py` (`LARGE_COMBOS`) lists the (case, kind) combinations that the scan (`scan.py`) marks usable, in order of
size. A combination is left out if either method's input guard rejects it or its dense result tables do not fit:

| case | APFC | RPFC | reason for an exclusion |
|---|---|---|---|
| case300 | yes | yes | |
| case_ACTIVSg500 | yes | yes | |
| case1354pegase | yes | yes | |
| case_ACTIVSg2000 | yes | yes | |
| case2383wp | yes | yes | |
| case2869pegase | yes | yes | |
| case3120sp | no | yes | APFC: `pfc:apfc:negativeLoss` (branches that deliver active power at both ends; the data have branches with negative series resistance) |
| case_ACTIVSg10k | yes | no | RPFC: the estimated dense reactive tables are above the memory limit |
| case9241pegase | no | no | APFC: `pfc:apfc:spectralRadius` (a circulating 2-cycle of gain 1 between two buses); RPFC: the estimated dense reactive tables are above the memory limit |

The scan itself, with the counts per case, is produced by `scan.py`; `cells.csv` of a series gives the reason for
every cell that is not timed.

### Slow cells

A single run (warm-up or recorded) longer than 30 s (`--slow`) makes the cell **slow**: no further warm-up run is
made, exactly 3 repetitions are recorded (fewer only on a time-out), and the budget and the minimum number of
repetitions do not apply. Slow cells are flagged in `cells.csv` and marked "(slow)" in `summary.md`. In
`run_timing.sh` the time-out of a single run is 1800 s, so that the old method of the large cases can finish.

### The series for publication

On a quiet machine on mains power, from the root of the matrix repository:

```sh
zsh benchmark/timing/export_large.sh
PYTHON=/Users/mirko/.pyenv/shims/python zsh benchmark/timing/run_timing.sh
```

Without `--smoke`, `--combos` or `--cases`, this runs every combination of `cases.py` (the 14 small cases with both
kinds and every usable large combination whose export exists) with the defaults above. A run with `--smoke` is a test
of the machinery and is labelled "smoke, not for publication".

- **Order of the two commands:** run `export_large.sh` well before the series, not right before it. The export runs Octave and leaves a high load average, which the first series (T14e) started under. The series now waits up to 300 s per cell for a quiet machine (`--quiet-wait`), but a cell whose wait expires is still timed and flagged `quiet` False.
- **After the series:** check that every cell in `cells.csv` has `quiet` True. For a cell flagged False, `environment.json` names the three busiest processes at that moment.
- **Expected stalls:** the old-method runs that stall (`old_stall_predictor.py`) run until the iteration cap or the time-out of 1800 s in the guard, and are recorded as "did not finish".

### Where the old method stalls: a static predictor

The old method is a sweep over the buses, not a numerical iteration (T14d): it processes a bus once, when every branch that flows into it is known, and never revisits it. `old_stall_predictor.py` builds the old method's initial inflow and outflow sets with its own preprocessing (`makeActiveFlows`, `makeReactiveFlows`) and the sign rules of its set construction. It then finds the buses that can never be processed: those in cyclic strongly connected components of the dependency digraph, those with an orphan inflow (a branch that flows into a bus but out of none, for example a reactive branch with zero exported from-end flow), and every bus downstream of these.
- **Can say:** whether the old method can finish ("can finish") or stalls ("predicted stall"), and which buses stay unprocessed. It was validated against the observed runs of T14d (equal sets) and every combination of the first series, for both kinds.
- **Cannot say:** whether the old method raises an exception inside a sweep (the division by zero of four large combinations), how many sweeps it needs, or how long a sweep takes. A time-out of a combination predicted "can finish" is a matter of size.
- **Usage:** `python -B benchmark/timing/old_stall_predictor.py --combos case2869pegase:rpfc,case_ACTIVSg10k:apfc [--observed <results>/iterations.csv]`. It does not run the old method.

### Correctness before timing

`correctness.py` compares, outside the timing, the port with the old method in memory on the large combinations
(iteration cap 10000, a time-out per run and a total budget), and reports the row sums and the determinism of the
port's tables. `octave_large.sh` (with `octave_exported.m`) runs `apfc.m` or `rpfc.m` under Octave on an exported
solved power flow, for the comparison of the port with Octave. The old method is not a reference for the rows of
RPFC dead-end branches and of APFC branches fed from both ends (pattern (+,−)); `correctness.py` excludes those rows
and counts them.
