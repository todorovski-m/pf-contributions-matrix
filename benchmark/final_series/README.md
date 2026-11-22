# Selected result files of the final benchmark series

Copies of selected result files of the three final runs behind the benchmark results of the paper. Each run wrote its
files into a `results/` folder, which `.gitignore` excludes. The files here are copies, made on 2026-11-22 with
`cp -p`, byte-identical to the originals with one exception: `guard.csv`, `memory.csv` and `summary.csv` of the timing
series were written with CRLF line ends, and git stores them with LF line ends (`core.autocrlf=input`); their content
is otherwise identical. The commit, the date and the labels below are those recorded by each run in its
`environment.json`. All three runs are series runs (label `series`, `smoke` false) of a clean working tree.

Not copied: the raw samples, the per-run logs and the other intermediate files (for the timing series: `samples.csv`,
`cells.csv`, `iterations.csv`, `summary_*slopes*.csv`, `export_synthetic.log` and the PNG figures). Not used: the
superseded runs `timing-20261117T203418Z`, `timing-20261119T090534Z` and the sources run `20261120T081314Z`.

| folder | run | produced by | matrix commit | date (local) | files (all copies) |
|---|---|---|---|---|---|
| `timing-20261119T203432Z/` | timing series `20261119T203432Z` (`benchmark/timing/results/20261119T203432Z`) | `benchmark/timing/run_timing.sh` (`timing.py`, `summarize.py`, `plot_timing.py`) | `9ae0f0b` | 2026-11-19 21:35 +01:00 | `summary.csv`, `summary.md`, `environment.json`, `memory.csv`, `guard.csv`, `timing_cells.pdf`, `timing_loglog.pdf`, `timing_ratio.pdf` |
| `baseline-20261120T095859Z/` | baseline run `benchmark/baseline_achayuthakan/results/20261120T095859Z` | `benchmark/baseline_achayuthakan/run_baseline_timing.sh` | `5f12bbf` | 2026-11-20 10:59 +01:00 | `baseline_timing.csv`, `guard.csv`, `summary.md`, `environment.json` |
| `sources-20261120T082459Z/` | sources run `benchmark/sources_scaling/results/20261120T082459Z` | `benchmark/sources_scaling/sources_scaling.py` (subcommands `run` and `plot`) | `8189cd5` | 2026-11-20 09:25 +01:00 | `sources.csv`, `summary.md`, `environment.json`, `sources_ratio.pdf`, `sources_ratio.png`, `sources_time.pdf`, `sources_time.png` |

The meaning of the columns of each file is documented in the README of the folder that produced it
(`benchmark/timing/README.md`, `benchmark/baseline_achayuthakan/README.md`, `benchmark/sources_scaling/README.md`).
The times were measured on one machine (recorded in `environment.json`: CPU, cores, thread settings and library
versions).
