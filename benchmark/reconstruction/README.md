# Reconstruction errors, properties and phase-shifter checks

## Reconstruction errors and properties (`reconstruction.py`)

For case1, case9, case30, case118 and tap_loop, APFC and RPFC, with the code of this repository (rule Q3 version 2.0):
- **reconstruction errors:** the row sums of the contribution table against the flows of the rows, separately for the
  from-end rows (P^f, Q^f), the to-end rows (P^t, Q^t) and the demand rows; the largest absolute error (MW for APFC,
  Mvar for RPFC) and the same divided by the largest absolute flow of those rows; for the port and for Octave;
- **agreement** of the port with Octave (labels and the largest absolute cell difference);
- **residual** ||A X - B||_inf of the port's solve (and of the RPFC self-term solve);
- **properties:** n_e (order of A), table rows and columns, nnz(N), density nnz(N)/n_e², the condition estimate of A
  (the port's `condest_A`: a 1-norm estimate with `scipy.sparse.linalg.onenormest`, as Matlab's `condest`), and the
  peak resident memory of one port run with checks off.

**Input files used:**
- `benchmark/inputs/<case>.mat` (the exported power flows) for the port;
- the Octave outputs `<case>_apfc.csv`, `<case>_rpfc.csv` of `tests/run_matrix.sh` (apfc.m and rpfc.m of this repository), in a
  scratch folder;
- `memory.csv` of the timing series `benchmark/timing/results/20261119T090534Z` (git-ignored) for the peak memory, not rerun.
  That series ran the code of commit d0701fa, before rule Q3 version 2.0; the memory of these small cases does not depend
  on the shunt rule in any visible way, but it was measured with the earlier code.

**Outputs** (`tables/`; a folder named `results/` would be ignored by `.gitignore`): `reconstruction.csv`, `properties.csv`, and `reconstruction_paper_cases.tex` (a booktabs
`tabular` with the four absolute reconstruction errors of the port for case1, case9, case30 and case118; it needs
`booktabs` and `amsmath`).

```sh
zsh tests/run_matrix.sh /tmp/recon_octave case1 case9 case30 case118 tap_loop
python3 -B benchmark/reconstruction/reconstruction.py /tmp/recon_octave benchmark/timing/results/20261119T090534Z benchmark/reconstruction/tables
```

The script exits with status 1 if any reconstruction error is above 1e-10.

## Phase-shifter cases and checks

- `cases/`: `syn9_pst`, `syn9_pst_tap` (the shifted branch is a reactive Case-3 branch), `syn9_pst_carrier`,
  `syn9_pst_carrier_tap` (the shifted branch is a reactive carrier), made by `make_pst_cases.m` (headers by `add_header.py`;
  the carrier variants were found by `search_pst_carrier.m`). They are kept out of `tests/cases_synthetic/`, so that the
  default case lists of the test suites and of the timing series do not change.
- `export_cases.sh` exports them (procedure of `tests/export_synthetic.sh`); `run_octave_cases.sh` runs apfc.m and rpfc.m on them.
- `pst_check.py`: the port against Octave, the row sums, the flows of the phase-shifting branch from the 2×2 stamp and its
  shunt and series flows against the branch model, and the old method fed the line-charging shunts.
- `screening_rule_q3_v2.md`: the screening of 15 further MATPOWER cases under rule Q3 version 2.0 (brief T16).
