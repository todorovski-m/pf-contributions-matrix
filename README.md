# Active Power Flow Components (APFC) and Reactive Power Flow Components (RPFC)

MATLAB implementation of matrix-based power flow contribution methods (version 2.0).

## Overview

This package provides matrix-based algorithms for decomposing active and reactive power flows in power networks into contributions from individual generators and reactive sources. The methods use a proportional-sharing principle and handle meshed networks, transformers with tap ratios, shunt elements, and reactive circulating flows.

**Key Features:**

- **APFC (Active Power Flow Components)**: Allocates active power flows to generators
- **RPFC (Reactive Power Flow Components)**: Allocates reactive power flows to generators and reactive sources
- **Matrix Formulation**: All source contributions follow from one linear system, solved with a single factorization for all right-hand sides
- **Branch-Level Resolution**: Separates from-end and to-end contributions for every branch
- **Shunt Handling**: Attributes line-charging capacitances, transformer tap-model shunts and bus-mounted capacitor banks
- **Circulating Flow Support**: Reactive circulating flows lead to a cyclic flow digraph; the system is solved as long as the spectral radius of the propagation matrix is below 1, which is checked
- **Validity Checks and Diagnostics**: Acyclicity, spectral radius, attenuation range and inflow guards, with a diagnostics output (`diag`)

## Requirements

- **MATLAB** (R2018b or later)
- **MATPOWER** (version 7.x)

**Octave.** The code also runs in GNU Octave, which lacks the Matlab table functions `array2table`, `addvars` and `writetable`. The folder `tests/octave_shims/` contains minimal stand-ins for them. They are **for validation only**: they keep the data and the labels, and write CSV files with 17 significant digits.

**What was run.** Version 2.0 was validated under Octave 11.1.0 with MATPOWER 7.1 and the stand-ins. **Nothing was run in Matlab.** `tests/matlab_check.m` is provided for a check in Matlab (see "Tests").

### Installation

1. Install MATPOWER:

   ```matlab
   % Add MATPOWER to MATLAB path
   addpath('/path/to/matpower7.x');
   ```

2. Ensure this code folder is in your MATLAB path:

   ```matlab
   addpath('/path/to/code');
   ```

3. Under Octave only, also add the stand-ins:

   ```matlab
   addpath('/path/to/code/tests/octave_shims');
   ```

4. The tests compare with the earlier iterative method, which is included as a git submodule. Clone with `git clone --recurse-submodules`, or run `git submodule update --init` in an existing clone.

## File Structure

### Core Functions

- **`apfc.m`**: APFC computation
  - Input: a case name, a MATPOWER case struct, or a solved case struct (see "Options")
  - Output: APFC contribution table, the preprocessed case struct, and the diagnostics struct
  - Constructs the demand-as-branches model with a virtual sink bus
  - Builds the system matrix and solves for all source contributions

- **`rpfc.m`**: RPFC computation
  - Input and output as for `apfc.m`
  - Classifies branches into three reactive cases based on the signs of the terminal flows
  - Computes series and shunt flow components and the shunt self-contributions
  - Handles reactive circulating flows; the spectral radius of the propagation matrix is checked

- **`pfc_checks.m`**: options and validity checks used by both functions (see "Validity Checks")

- **`batch_solve.m`**: runs APFC and RPFC on the test cases and writes the CSV files to `results/`

### Test Cases (`cases/` directory)

Standard MATPOWER test networks:

- **`case1.m`**: Three-bus system (verification example)
- **`case9.m`**: IEEE 9-bus system
- **`case14.m`**: IEEE 14-bus system
- **`case24_ieee_rts.m`**: IEEE RTS 24-bus system (several generators at one bus)
- **`case30.m`**: IEEE 30-bus system (transformers and capacitor banks)
- **`case33.m`**: 33-bus distribution network
- **`case118.m`**: IEEE 118-bus system
- **`tap_loop.m`**: Parallel transformers with a reactive circulating flow

Synthetic variants of case9 are in `tests/cases_synthetic/`: negative active and reactive demand and generation, and a Case-3 branch with an inductive shunt (`syn9_case3_tap`).

### Results

- **`paper_results/`**: the outputs of the archived version 1.0 (without two stray files, see `RELEASE_NOTES_v2.0.md`), including its known labelling defect (the Case-3 source columns of `case9_rpfc.csv` and `case118_rpfc.csv` are named in a different order from their data). `reference_results/` supersedes it; see `RELEASE_NOTES_v2.0.md` for all differences.
- **`reference_results/`**: the outputs of version 2.0 for the six shared cases (case1, case9, case14, case30, case33, case118) and `tap_loop`. They were written by `tests/make_reference.sh` under Octave, with 15 significant digits in the format of Matlab's `writetable` (C format `%.15g`, negative zero as `-0`). Matlab writes the same format. Its values agree with the Octave values to about 1e-11 (for the version 1.0 code, the largest difference between the Matlab files in `paper_results/` and the Octave results was below 1e-11), so individual last digits, and values near zero (rounding noise below 1e-13), can differ.
- **`results/`**: the output folder of `apfc`, `rpfc` and `batch_solve` (not tracked).

## Usage

### Basic Usage

#### Compute APFC for a single case

```matlab
[APFC_tbl, mpc, diag] = apfc('cases/case1');
```

#### Compute RPFC for a single case

```matlab
[RPFC_tbl, mpc, diag] = rpfc('cases/case1');
```

#### Reuse one power flow for both functions

```matlab
r = runpf(loadcase('cases/case9'), mpoption('verbose', 0, 'out.all', 0));
opts = struct('solved', true, 'name', 'cases/case9');
[APFC_tbl, ~, dA] = apfc(r, opts);
[RPFC_tbl, ~, dR] = rpfc(r, opts);
```

#### Run batch computation

```matlab
batch_solve;
```

### Options

`apfc(casedata, opts)` and `rpfc(casedata, opts)` accept an optional struct `opts`. Unknown fields are an error.

| Field | Default | Meaning |
|---|---|---|
| `tol` | 1e-6 | Threshold for deletion and classification decisions, in the units of the case (MW, Mvar): removal of zero-flow branches and zero-output generators, and the sign tests of demands, generators, bus shunts and inductive branch shunts. |
| `tol_div` | 1e-9 | Used for three purposes: (1) the guards against division by a near-zero value (the series flow at the from-end before the self-contribution of a branch is computed, and the total inflow of a bus with an outgoing branch, see "Validity Checks"); (2) the threshold below which a branch shunt self-contribution is neglected; (3) the slack of the upper bound of the active attenuation, 1 + `tol_div`. |
| `write` | true | Write `results/<case>_apfc.csv` or `results/<case>_rpfc.csv`. A struct input writes a file only if `name` is given. |
| `solved` | false | `casedata` is a struct that already holds a solved power flow (for example the output of `runpf`); `runpf` is not called. |
| `name` | `''` | Case name for the output file when `casedata` is a struct. |

`casedata` is a case name, an unsolved MATPOWER case struct (it is solved with `runpf`), or, with `solved = true`, a solved struct. A power flow that does not converge is an error.

### Diagnostics (`diag`)

| Field | Meaning |
|---|---|
| `case` | case name (empty for a struct without `name`) |
| `n_bus`, `n_branch`, `n_demand`, `n_source`, `n_e` | numbers of buses, physical branches (after preprocessing), demand branches, source columns, and branches of the expanded network |
| `nnz_N`, `nnz_A` | non-zeros of the propagation matrix N and the system matrix A = I + N |
| `issparse` | struct with `N`, `A`, `B`: whether the matrices are sparse |
| `condest_A` | estimate of the 1-norm condition number of A (`condest`) |
| `attenuation_range` | smallest and largest attenuation factor of the physical branches |
| `removed` | numbers of removed branches and generators |
| `flipped_branches` | number of branches whose orientation was reversed |
| `converted` | numbers of converted elements (negative demand to source, negative generation to demand; in `rpfc` also inductive bus shunts, Case-2 and Case-3 branches with both terminals non-zero, inductive transformer shunts), and `zero_terminal_demands`, `zero_terminal_sources`: the branches with exactly one zero terminal that became one terminal demand or one terminal source (see "Terminal Classification"; `zero_terminal_sources` is always 0 in `apfc`) |
| `acyclic`, `spectral_radius` | acyclicity of the flow digraph (Kahn's topological sort) and the spectral radius of N (exactly 0 for an acyclic digraph) |
| `warnings` | identifiers of the warnings issued |
| `negative_series_flow` | `rpfc` only: `count` and `branches` (index in the returned `mpc.branch` and bus pair) of the branches with a negative series flow |
| `opts` | the options used |

### Output Format

Both functions return tables with:

- **Rows**: demands first (in APFC: `D_b`, `G_b`, then the terminal demands `C2_b` of branches with one zero terminal), then the from-ends and the to-ends of the branches; RPFC appends the rows of the Case-3 branches (their own terminal injections) at the end
- **Columns**: source contributions
- **Values**: Power contributions (MW for APFC, Mvar for RPFC)
- **Row sums**: Match solved power flows (verification check)

## Preprocessing Rules

| Rule | APFC | RPFC |
|---|---|---|
| Tolerance | `opts.tol` (default 1e-6) | `opts.tol` (default 1e-6) |
| Deletion (per quantity) | a branch is removed if both active terminal flows are below `tol` in magnitude; a generator if its active output is | a branch is removed if both reactive terminal flows are below `tol` in magnitude; a generator if its reactive output is |
| Orientation | a branch with a negative from-end flow is reversed (sign test at exactly 0), except a branch with one zero terminal | a branch with both terminal flows negative is reversed (sign tests at exactly 0), except a branch with one zero terminal |
| Demand | a bus with PD > `tol` is a demand (row `D_b`) | a bus with QD > `tol` is a demand (row `D_b`) |
| Negative demand | a bus with PD < -`tol` is a source (column `D_b`) | a bus with QD < -`tol` is a source (column `D_b`) |
| Negative generation | a generator with PG < -`tol` is a demand (row `G_b`) | a generator with QG < -`tol` is a demand (row `G_b`) |
| Reactive branch cases (after orientation) | - | Case 1, both terminal flows positive: a carrier. Case 2, from-end positive and to-end negative (consumes at both ends): replaced by two demands (rows `C2_b` at both buses). Case 3, from-end negative and to-end positive (injects at both ends): removed from the carriers, its two terminal injections are sources (columns `E_i_j_f`, `E_i_j_t`). The tests use the signs at exactly 0. |
| Bus shunts | - | Q = VM² · BS: a source (column `S_b`) if above `tol`, a demand (row `Bs_b`) if below -`tol` |
| Branch shunts | - | The shunt flows at both ends follow from the series current. An inductive shunt (below -`tol`) is a demand (row `Tf_b` or `Tt_b`), except on a Case-2 branch, whose two terminal demands already contain its shunt power. Capacitive branch shunts give the self-contribution columns `E_i_j_s`. |
| Attenuation | a = Pt / Pf | a = Qe / Qs, from the series flows (shunts removed) |

### Terminal Classification

The pattern of the signs of the two terminal flows decides how a branch enters the method. A terminal flow x has the sign + if x > `tol`, - if x < -`tol`, and 0 otherwise. The terminal flows are PF and the delivered PT (APFC), and QF and the delivered QT (RPFC). The pattern is taken before the orientation flip, on the branches that are not deleted:

| Pattern (from, to) | APFC | RPFC |
|---|---|---|
| (+,+) | carrier | carrier (Case 1) |
| (-,-) | reversed, then carrier | reversed, then carrier |
| (+,-) | error `pfc:apfc:nonPositiveAttenuation` | Case 2: two terminal demands (rows `C2_b`) |
| (-,+) | reversed, then error `pfc:apfc:nonPositiveAttenuation` | Case 3: two terminal sources (columns `E_i_j_f`, `E_i_j_t`) |
| (+,0) | one terminal demand PF at the from bus (row `C2_b`) | Case 2 with the zero terminal omitted: one demand QF at the from bus (row `C2_b`) |
| (0,-) | one terminal demand -PT at the to bus (row `C2_b`) | Case 2 with the zero terminal omitted: one demand -QT at the to bus (row `C2_b`) |
| (-,0) | error `pfc:apfc:negativeLoss` | Case 3 with the zero terminal omitted: one source -QF at the from bus (column `E_i_j_f`, row `L_{f.i-j}`) |
| (0,+) | error `pfc:apfc:negativeLoss` | Case 3 with the zero terminal omitted: one source QT at the to bus (column `E_i_j_t`, row `L_{t.i-j}`) |

- **Dead-end branches.** A branch with exactly one zero terminal typically connects a dead-end bus (no load, no generator, no other branch) and carries only its own loss (APFC) or its own charging or shunt power (RPFC). It is not a carrier: its non-zero terminal is a demand or a source at the other bus. The dead-end bus therefore has no outgoing branch, and the zero-inflow guard does not fire there.
- **Shunt rules.** A one-zero branch in RPFC follows the shunt rules of its case: a Case-2 member gets no transformer-shunt rows (`Tf_b`, `Tt_b`); a Case-3 member keeps them.
- **Robustness.** Such a branch is never reversed and its result does not depend on the sign of a residue at its zero end.
- **Row and column order:**
  - APFC: the terminal demands follow the loads and the generators with negative output; first the from-terminals, then the to-terminals.
  - RPFC: the Case-2 demands are first the from-terminals of all members, then their to-terminals. A Case-3 member contributes one column in its place, and its single row is appended with the other Case-3 rows.
- **Other branches:** every branch with both terminals non-zero keeps the exact-sign tests above.

### Difference from the earlier iterative method

The earlier iterative method (submodule `external/iterative-method`) prepares its input with `my_runpf.m`:
- It removes a branch only when all four terminal flows (active and reactive) are below 1e-4.
- It then sets every single flow below 1e-4 to zero.

The matrix code differs in two ways:
- It removes a branch per quantity, when both terminal flows of that quantity are below `opts.tol`.
- It does not set single flows to zero; it decides the signs of the flows at exactly 0.

In the six shared test cases (case1, case9, case14, case30, case33, case118), no flow has a magnitude between 1e-6 and 1e-4, so the two thresholds classify every flow in the same way. The only difference is the branch 7-8 of case14, which carries no active power: the matrix code removes it from the active network, while the iterative method keeps it with zero active flows. In all test cases, the contributions agree with the iterative method to 1e-10 (see "Tests and Validation").

## Validity Checks

The checks are in `pfc_checks.m` and use identifiers of the form `pfc:...`. They fall into three groups:

| Group | Condition | Identifier |
|---|---|---|
| Error (results not valid) | cyclic active digraph with spectral radius of N at or above 1 | `pfc:apfc:spectralRadius` |
| | cyclic reactive digraph with spectral radius of N at or above 1 (the sufficient condition for the Neumann series is not met) | `pfc:rpfc:spectralRadius` |
| | active attenuation at or below 0 (a defensive check; reachable only by a terminal flow equal to the tolerance, so the message asks to check the input) | `pfc:apfc:nonPositiveAttenuation` |
| | active power delivered at one end of a branch whose other end carries zero flow (negative loss; patterns (0,+) and (-,0), see "Terminal Classification") | `pfc:apfc:negativeLoss` |
| | total inflow below `tol_div` at a bus with an outgoing branch | `pfc:apfc:zeroInflow`, `pfc:rpfc:zeroInflow` |
| | to-end shunt self-term of a reactive carrier below `-tol_div` (it would not be propagated; cannot occur for a valid carrier, see `rpfc.m`) | `pfc:rpfc:negativeSelfTerm` |
| | invalid input, unknown option, unsolved or non-converged power flow | `pfc:options:...`, `pfc:input:...`, `pfc:runpf:notConverged` |
| Warning (can invalidate the results) | cyclic active digraph with spectral radius below 1 | `pfc:apfc:cyclic` |
| | active attenuation above 1 + `tol_div` | `pfc:apfc:attenuationAboveOne` |
| Information (in `diag` only, not printed) | cyclic reactive digraph with spectral radius below 1 (a reactive circulating flow) | `diag.acyclic`, `diag.spectral_radius` |
| | branches with a negative series flow (negative reactive attenuation) | `diag.negative_series_flow` |
| | condition number estimate of A | `diag.condest_A` |

What each check does:
- **Acyclicity.** Kahn's topological sort on the pattern of the upstream adjacency matrix decides whether the flow digraph is acyclic. Acyclic means N is nilpotent, so the spectral radius is exactly 0.
- **Spectral radius.** Only for a cyclic digraph: the largest absolute eigenvalue of N on its cyclic core, which is what remains after repeatedly removing branches without upstream or without downstream branches. `eig` is used up to 2000 core branches and `eigs` above.
- **Attenuation range (active).** Every a = Pt / Pf must lie in (0, 1 + `tol_div`].
- **Inflow guard.** The shares of the outgoing branches of a bus divide by its total inflow.
- **Negative series flow (reactive).** Information only, see below.

Not checked: a spectral radius close to (but below) 1, and a near-singular A. `diag.spectral_radius` and `diag.condest_A` report the values.

**Limitations.**
- **Negative loss (APFC).** A branch with one zero terminal that delivers active power at its other end, patterns (0,+) and (-,0), is an error (`pfc:apfc:negativeLoss`), not a source. Such a branch generates active power: in a MATPOWER case this means a negative series resistance in the data. The message names the branches.
- **Fed from both ends (APFC).** A branch whose two active terminal flows are both non-zero and draw power from both buses, pattern (+,-), remains the error `pfc:apfc:nonPositiveAttenuation`. It occurs where a branch carries almost no active power and its loss is larger than the transfer.

## Negative Reactive Attenuation in Standard Cases

Some branches carry positive reactive flow at both terminals, but the series element absorbs reactive power from both ends, so the series flow at the to-end is negative and the reactive attenuation Qe / Qs is negative. This occurs mainly on transformers with an off-nominal tap. Negative attenuation does not invalidate the results; the contributions agree with the earlier iterative method. But some contributions become negative. Generated by `python3 -B tests/negative_attenuation_table.py <rpfc outputs> benchmark/inputs`:

| case | carrier branches with negative attenuation | of which transformers | most negative attenuation | negative contribution cells | smallest negative contribution (Mvar) |
|---|---|---|---|---|---|
| case1 | 0 | 0 | - | 0 | - |
| case9 | 0 | 0 | - | 0 | - |
| case14 | 2 | 1 | -2.012 | 249 | -5.792 |
| case24_ieee_rts | 4 | 2 | -0.918 | 77 | -4.466 |
| case30 | 0 | 0 | - | 0 | - |
| case33 | 0 | 0 | - | 0 | - |
| case118 | 12 | 8 | -4.708 | 1084 | -130.175 |
| tap_loop | 1 | 1 | -0.327 | 4 | -1.658 |

"Negative contribution cells" excludes the own terminal injection of a Case-3 branch in its own row, which is negative by design. The branches are listed in `diag.negative_series_flow`.

## Row and Column Labels

| Label | Kind | Meaning |
|---|---|---|
| `G_b` | column | generator at bus b with positive output |
| `S_b` | column (RPFC) | bus shunt at bus b that injects reactive power |
| `D_b` | column | bus b with negative demand, used as a source |
| `E_i_j_f`, `E_i_j_t` | column (RPFC) | from-end and to-end injection of the Case-3 branch i-j |
| `E_i_j_s` | column (RPFC) | shunt self-contribution of branch i-j |
| `D_b` | row | demand at bus b |
| `G_b` | row | generator at bus b with negative output (a demand) |
| `Bs_b` | row (RPFC) | inductive bus shunt at bus b |
| `C2_b` | row (RPFC) | terminal demand at bus b of a Case-2 branch |
| `Tf_b`, `Tt_b` | row (RPFC) | inductive branch shunt at the from-end or to-end, at bus b |
| `L_{f.i-j}`, `L_{t.i-j}` | row | from-end and to-end of branch i-j (after orientation) |

`D_b` and `G_b` appear as rows and as columns. A row is a demand; a column is a source.

**Suffixes.** A label that would repeat gets the suffix `_2`, `_3`, ... on its second and later occurrences. Examples: parallel branches (`L_{f.1-2}`, `L_{f.1-2}_2`, with the matching to-end rows); several generators at one bus (`G_22`, `G_22_2`, ...); two Case-2 branch ends at one bus (`C2_49`, `C2_49_2`). A suffix that would collide with an existing label is skipped.

## Key Results from Test Cases

### Case 1 (Three-Bus System)

This verification case is the three-bus system used in the paper's worked
example. It demonstrates the demand-as-branches construction, the virtual sink
bus used by APFC, and the reactive preprocessing used by RPFC.

![Three-bus system](figures/simple_system.png)

Active power flows, attenuation factors, and inflow weights:

| ℓ   | from |  to | P<sup>f</sup><sub>ℓ</sub> | P<sup>t</sup><sub>ℓ</sub> | a<sup>P</sup><sub>ℓ</sub> |
| --- | ---: | --: | ------------------------: | ------------------------: | ------------------------: |
| 1   |    1 |   2 |                    38.608 |                    36.545 |                    0.9466 |
| 2   |    1 |   3 |                    46.119 |                    43.506 |                    0.9434 |
| 3   |    2 |   3 |                     6.545 |                     6.494 |                    0.9921 |
| 4   |    2 |   - |                        50 |                        50 |                         1 |
| 5   |    3 |   - |                        50 |                        50 |                         1 |

Reactive power flows and branch classification:

| ℓ   | from |  to | Q<sup>f</sup><sub>ℓ</sub> | Q<sup>t</sup><sub>ℓ</sub> | Q<sup>sf</sup><sub>ℓ</sub> | Q<sup>st</sup><sub>ℓ</sub> | Q<sup>s</sup><sub>ℓ</sub> | Q<sup>e</sup><sub>ℓ</sub> | Case |
| --- | ---: | --: | ------------------------: | ------------------------: | -------------------------: | -------------------------: | ------------------------: | ------------------------: | ---: |
| 1   |    1 |   2 |                    26.148 |                    22.378 |                      1.832 |                      1.430 |                    27.980 |                    20.949 |    1 |
| 2   |    1 |   3 |                    25.595 |                    19.945 |                      1.832 |                      1.422 |                    27.427 |                    18.523 |    1 |
| 3   |    2 |   3 |                    -2.622 |                     0.055 |                      1.430 |                      1.422 |                    -1.192 |                    -1.367 |    3 |

APFC result (MW):

| Row                         | G<sub>1</sub> | G<sub>2</sub> |  Total |
| --------------------------- | ------------: | ------------: | -----: |
| D<sub>2</sub>               |        32.315 |        17.685 |     50 |
| D<sub>3</sub>               |        47.703 |         2.297 |     50 |
| ℓ<sub>1-2</sub><sup>f</sup> |        38.608 |             0 | 38.608 |
| ℓ<sub>1-3</sub><sup>f</sup> |        46.119 |             0 | 46.119 |
| ℓ<sub>2-3</sub><sup>f</sup> |         4.230 |         2.315 |  6.545 |
| ℓ<sub>1-2</sub><sup>t</sup> |        36.545 |             0 | 36.545 |
| ℓ<sub>1-3</sub><sup>t</sup> |        43.506 |             0 | 43.506 |
| ℓ<sub>2-3</sub><sup>t</sup> |         4.197 |         2.297 |  6.494 |

RPFC result (Mvar):

| Row                         | G<sub>1</sub> | ℓ<sub>2-3</sub><sup>f</sup> | ℓ<sub>2-3</sub><sup>t</sup> | ℓ<sub>1-2</sub><sup>s</sup> | ℓ<sub>1-3</sub><sup>s</sup> |  Total |
| --------------------------- | ------------: | --------------------------: | --------------------------: | --------------------------: | --------------------------: | -----: |
| D<sub>2</sub>               |        15.662 |                       2.097 |                           0 |                       2.241 |                           0 |     20 |
| D<sub>3</sub>               |        17.286 |                           0 |                       0.055 |                           0 |                       2.659 |     20 |
| G<sub>2</sub>               |         3.915 |                       0.524 |                           0 |                       0.560 |                           0 |      5 |
| ℓ<sub>1-2</sub><sup>f</sup> |        26.148 |                           0 |                           0 |                           0 |                           0 | 26.148 |
| ℓ<sub>1-3</sub><sup>f</sup> |        25.595 |                           0 |                           0 |                           0 |                           0 | 25.595 |
| ℓ<sub>1-2</sub><sup>t</sup> |        19.577 |                           0 |                           0 |                       2.801 |                           0 | 22.378 |
| ℓ<sub>1-3</sub><sup>t</sup> |        17.286 |                           0 |                           0 |                           0 |                       2.659 | 19.945 |
| ℓ<sub>2-3</sub><sup>f</sup> |             0 |                      -2.622 |                           0 |                           0 |                           0 | -2.622 |
| ℓ<sub>2-3</sub><sup>t</sup> |             0 |                           0 |                       0.055 |                           0 |                           0 |  0.055 |

### Case 9 (IEEE 9-Bus)

- Illustrates meshed topology with multi-source loads
- Demonstrates line-charging dominance in reactive supply
- Reference: Generators supply residual reactive balance

### Case 30 (IEEE 30-Bus)

- Tests transformers with tap-ratio attribution
- Validates capacitor bank contributions
- Verifies from-end to-end reactive losses

### Case 118 (IEEE 118-Bus)

- Large network with 54 generators and 186 branches
- One factorization of the system matrix serves all source columns
- Contains carrier branches with negative reactive attenuation and negative contributions (see the table above)

### Case 24 (IEEE RTS 24-Bus)

- Several generators at one bus (unique column names `G_b`, `G_b_2`, ...)
- Case-2 transformers with an inductive shunt (no `Tf`/`Tt` demand for them, see "Preprocessing Rules")

### Tap Loop (Circulating Flows)

- Parallel transformers with a tap mismatch (taps 1.05 and 0.95)
- The reactive flow digraph is cyclic: `diag.acyclic` is false and `diag.spectral_radius` is 0.17937561102916488 (below 1, so A is invertible by the Neumann series)
- The earlier iterative method cannot resolve this case: each of the two buses waits for a contribution that only the other bus can provide, so no bus is ever processed
- Transformer A has a negative series flow at its to-end (`diag.negative_series_flow`)

## Tests and Validation

The row sums of the contribution matrices match the solved power flows, at every from-end and to-end, in every test case. The tests check this and more. They run under Octave with MATPOWER, in scratch folders under `/tmp`. They need zsh, the submodule `external/iterative-method`, and Python 3 with numpy, pandas, scipy, networkx and matplotlib (for the earlier iterative method).

```sh
zsh tests/run_tests.sh            # all tests, one PASS or FAIL line each
zsh tests/make_reference.sh       # regenerate reference_results/
```

`tests/run_tests.sh` runs:
- the unit tests in `tests/unit/` (options and input forms, the checks with crafted inputs, the diagnostics of all test cases, and the terminal classification: the nine sign patterns, the dead-end cases at the from side and at the to side, invariance to a perturbation of the zero terminal, and the equivalence of a dead-end branch with the terminal injection it represents);
- the comparison of the active dead ends with the earlier iterative method (`tests/dead_end_equivalence.py`);
- the comparison with the earlier iterative method on in-memory results (cells and row sums within 1e-10);
- the comparison with `paper_results/`, where only the documented differences of version 2.0 are allowed;
- the regeneration of `reference_results/`, which must match the committed files byte for byte.

The synthetic case `tests/cases_synthetic/syn9_carrier_tap.m` (case9 with tap 1.05 on line 8-9, which stays a reactive carrier with an off-nominal tap and line charging) is run by `apfc.m` and `rpfc.m` and by the Python port, but it is not compared with the earlier iterative method: that method's preprocessing (`my_runpf.m`) computes the from-end shunt flow of such a line with B/2 instead of B/(2·tap²), so it disagrees there (see "Python Port"). The case list of that comparison is given explicitly in `tests/run_tests.sh`.

The synthetic cases `tests/cases_synthetic/syn9_dead_*.m` add to case9 a dead-end bus attached to bus 5 by a branch with line charging, with an inductive susceptance, or with a bus shunt at the dead end, at the from side and at the to side. They are not in that case list either: the earlier iterative method is not a reference for reactive dead ends (it does not finish, or its row sums differ from the flows). Their active part is compared with it by `tests/dead_end_equivalence.py`.

**Matlab check.** In Matlab, from the root of the repository and without `tests/octave_shims` on the path, run `matlab_check('/some/output/folder')`. It runs both functions on case9, case118 and case24_ieee_rts, writes the CSV files with Matlab's `writetable`, and prints the label agreement and the maximum differences against the Octave results in `tests/octave_results/`. It has not been run, because Matlab was not available.

## Python Port

`python/pfc_matrix` is a Python port of `apfc.m`, `rpfc.m` and `pfc_checks.m`. It reads the same exported inputs as the earlier iterative method: the `results.mat` files written by `my_runpf.m` (`benchmark/inputs/`), loaded with `load_results` of `external/iterative-method/power_flow.py`. It does not run a power flow, and the package does not import the iterative method.

**Installation.** Python 3.10 with numpy, scipy and pandas; networkx, matplotlib and pytest for the tests. The versions used are in `python/requirements.txt`. No virtual environment is needed; with an existing Python, `python -m pip install -r python/requirements.txt` installs them.

**API.**

```python
import sys; sys.path.insert(0, "python")
from pfc_matrix import Case, read_branch_params, apfc, rpfc, write_csv, Options
bus, gen, branch, _ = load_results("benchmark/inputs/case9.mat")   # from external/iterative-method/power_flow.py
case = Case.from_frames(bus, gen, branch, "case9", params=read_branch_params("benchmark/inputs/case9.mat"))
table, info = apfc(case)                        # pandas DataFrame (index "Line") and dict
table_q, info_q = rpfc(case, Options(tol=1e-6, tol_div=1e-9))
write_csv(table_q, "case9_rpfc.csv")            # the CSV layout of apfc.m/rpfc.m, %.17g
```

- `Case.from_frames(bus, gen, branch, name, params=None)` copies the DataFrames and checks the columns and the bus numbering (1..nb). `params` holds BR_R, BR_X, BR_B and baseMVA, which `load_results` drops. `read_branch_params(path)` reads them from the same `.mat` file, and `Case.from_frames` checks that its rows align with the branch table. `apfc` does not need them; `rpfc` does, because it computes the branch shunt flows as `rpfc.m` does.
- `apfc(case, options=None, tol=None, tol_div=None)` and `rpfc(...)` return `(table, info)`. The table has the row labels, column labels, order and values of the CSV files of `apfc.m` and `rpfc.m`. `info` has every field of `diag` under the same name, plus `solver` and, for `rpfc`, `self_toend`.
- `info["removed"]["branches"]` counts the branches removed from the exported input. It can be smaller than in `diag` by the branches that `my_runpf.m` deleted before the export (all four flows below 1e-4).
- `Options(tol=1e-6, tol_div=1e-9, write=True, solved=False, name='', checks=True)` has the names and defaults of `opts`. `write` and `solved` have no effect: the port never writes files (use `write_csv`), and its input is always solved. The keyword arguments `tol` and `tol_div` win over the object. Unknown names and invalid values raise `PfcError` (`pfc:options:unknownField`, `pfc:options:value`).

**Phases.** Each method runs four functions in order, none of which prints, times or writes a file:
- `prepare_active` or `prepare_reactive`: preprocessing and the index and flow arrays;
- `assemble_active` or `assemble_reactive`: the sparse matrices N, A and B (and the self-term right-hand side), with the validity checks;
- `solve_active` or `solve_reactive`: one SuperLU factorisation of A, used for all right-hand sides;
- `build_active_table` or `build_reactive_table`: the output table and its labels.

**Checks and identifiers.** `pfc_matrix.checks` ports `pfc_checks.m`, with the same identifiers and messages. Errors raise `PfcError` (a `ValueError`; the identifier is in `.identifier`). Warnings are collected in `info["warnings"]` and issued with `warnings.warn` as `PfcWarning`. The identifiers:
- errors: `pfc:apfc:spectralRadius`, `pfc:rpfc:spectralRadius`, `pfc:apfc:nonPositiveAttenuation`, `pfc:apfc:negativeLoss`, `pfc:apfc:zeroInflow`, `pfc:rpfc:zeroInflow`, `pfc:rpfc:negativeSelfTerm`, `pfc:options:...`, `pfc:input:casedata` and `pfc:input:notSolved`;
- errors that exist only in Python: `pfc:input:columns`, `pfc:input:busNumbering` and `pfc:input:branchParams`;
- warnings: `pfc:apfc:cyclic` and `pfc:apfc:attenuationAboveOne`.

`condest_A` is estimated with `scipy.sparse.linalg.onenormest` (seeded, so that two calls give the same value), and may differ from Matlab's `condest`. `Options(checks=False)` switches off the attenuation check, the warning `pfc:apfc:cyclic` and `condest_A` (which is then `None`). The zero-inflow, spectral-radius, negative-self-term and negative-loss guards are always on.

**Terminal classification.** Phase (a) classifies the branches with one zero terminal as `apfc.m` and `rpfc.m` do (see "Terminal Classification"), and `info["converted"]` has `zero_terminal_demands` and `zero_terminal_sources`. The negative-loss error is raised at the start of phase (b), before the other checks.

**Tests.** `zsh python/run_tests.sh` runs the pytest tests (`python -B -m pytest -p no:cacheprovider`). If Octave is available, it first creates the synthetic exports and the Octave outputs in a scratch folder (`tests/export_synthetic.sh`, `tests/run_matrix.sh`), so that the port is compared with `apfc.m` and `rpfc.m` on every case. Without Octave, those tests are skipped and the runner prints their number. Set `PYTHON` to choose the interpreter.

**Known issue of the earlier iterative pipeline.** `my_runpf.m` computes the from-end shunt flow QSF of a branch with B/2 instead of B/(2·tap²), which matters only for a branch with both an off-nominal tap and line charging. The port therefore does not use the exported QSF and QST, but computes the shunt flows as `rpfc.m` does.

## Related Work

- [Todorovski, M. and Rajičić, D., "Decomposing total power flows into contributions originating from power sources in loads, generators, and network elements," _Electric Power Systems Research_, 2024.](https://doi.org/10.1016/j.epsr.2024.110433)
- [Enshaee, A. and Enshaee, P., "New reactive power flow tracing and loss allocation algorithms for power grids using matrix calculation," _International Journal of Electrical Power & Energy Systems_, 2017.](https://doi.org/10.1016/j.ijepes.2016.11.010)
- [Gutierrez, A., et al., "Power Flow Tracing in Power Systems with Active and Reactive Circulating Power Flows," _Electric Power Components and Systems_, 2023.](https://doi.org/10.1080/15325008.2023.2175928)
- [Achayuthakan, C., Dent, C.J., Bialek, J.W., and Ongsakul, W., "Electricity Tracing in Systems With and Without Circulating Flows: Physical Insights and Mathematical Proofs," _IEEE Transactions on Power Systems_, vol. 25, no. 2, pp. 1078–1087, May 2010.](https://doi.org/10.1109/TPWRS.2009.2037506)

## License

This project is licensed under the MIT License. See the LICENSE file for details.
