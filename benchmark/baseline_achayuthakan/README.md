# Independent baseline: Achayuthakan et al. (2010), active power

A bus-level implementation of the matrix tracing method of C. Achayuthakan, C. J. Dent, J. W. Bialek and W. Ongsakul,
"Electricity tracing in systems with and without circulating flows: physical insights and mathematical proofs",
IEEE Trans. Power Syst., vol. 25, no. 2, pp. 1078-1086, 2010, compared with APFC at demand level (brief T14h; reviewer
comment R1.4). The method code of the matrix method is not used or changed; `achayuthakan.py` imports nothing of it.

## Files

| file | content |
|---|---|
| `SPEC.md` | the published method in our notation, with page and equation numbers, and its ambiguities |
| `achayuthakan.py` | the baseline: network, distribution matrices, contributions; explicit inverse and sparse LU |
| `compare.py` | the comparison with APFC at demand level, matrix properties and peak memory (no timing) |
| `tables/` | `comparison.csv`, `comparison_paper_cases.tex` (booktabs `tabular`; needs `amsmath`), `summary.md` |
| `baseline_timing.py`, `run_baseline_timing.sh` | the runtime measurement (run on a quiet machine) |
| `tests/` | tests of the baseline and of the timing machinery |

## Decisions (stop condition 1 of T14h)

The paper derives the method for a lossless network and names two loss variants without choosing one (p. 1080). Both
are implemented, with the actual generations and demands:
- **net-downstream:** A_d from the receiving-end flows, nodal flow P = inflow (generation plus received flows), and
  eq. (14): the contribution of the generator at bus i to the demand at bus k is (P_G,i / P_i) [A_d^-1]_ik P_D,k.
  The contributions to a demand sum to the demand;
- **gross-upstream:** A_u from the sending-end flows, P = outflow (demand plus sent flows), and the upstream
  counterpart (P_D,k / P_k) [A_u^-1]_ki P_G,i. The contributions of a generator sum to its generation, so the demands
  carry the losses.

Branches fed from both ends (+,-) and branches with one zero end (+,0), (0,-) are not lines of the digraph but demands
at their own ends, as the Case-2 rows of `apfc.m`; the paper does not cover them (they occur in case_ACTIVSg2000 only).
Parallel lines are added. Several generators at a bus are one nodal generation (the APFC columns are summed per bus
for the comparison).

## Comparison with APFC (no timing)

```sh
python -B benchmark/baseline_achayuthakan/compare.py benchmark/baseline_achayuthakan/tables
```

Seven networks (case1, case9, case14, case30, case118, case_ACTIVSg500, case_ACTIVSg2000), the exported power flows of
the final timing series. The results are in `tables/`.
- **net-downstream equals APFC at demand level** within 4.3e-14 MW in every network: it uses the same loss treatment
  (the shares of the nodal inflow, with the received to-end flows).
- **gross-upstream differs** by up to 5.24 MW (9.3 % of a demand): it allocates the network losses to the demands. The
  allocated total equals the losses of the lines of its digraph up to the power-flow mismatch of the export, and
  differs from sum(PF - PT) in case_ACTIVSg2000 by the 0.626 MW of the four branches treated as terminal demands.
- Matrix orders: n_b (baseline, buses with a nonzero nodal flow) against n_e (APFC, branches and demands).

## Runtime (on a quiet machine)

```sh
cd /Users/mirko/Mirko/research/components_2026/matrix
python -B benchmark/baseline_achayuthakan/baseline_timing.py estimate
zsh benchmark/baseline_achayuthakan/run_baseline_timing.sh
```

For each network: the port (checks off) and the four baseline methods (net and gross, explicit inverse and sparse LU),
each as the complete computation from the power-flow quantities to the demand-level contribution table, with the
machinery of `benchmark/timing/` (imported, not changed): a guard (net against APFC within 1e-9 MW, inverse against
sparse within 1e-10), the quiet wait before each cell, warm-up and repetitions within a budget of 20 s per cell, the
median and the slow-cell rule, and the environment record. Output: `results/<timestamp>/` (ignored by git) with
`baseline_timing.csv`, `guard.csv`, `environment.json` and `summary.md` (with the ratio of each baseline method to the
port). Estimate from the settings: 35 cells, at least 700 s; plus up to 300 s of quiet wait per cell on a busy machine.

## Tests

```sh
python -B -m pytest -p no:cacheprovider benchmark/baseline_achayuthakan/tests
```

They reproduce T_d of the paper's Fig. 1 and Fig. 4 examples (eqs. (21), (26)) and A_d of eq. (11), check that net and
gross agree on those lossless examples, and check conservation, nonnegativity and the agreement of the two solution
methods on case9 and case118; the timing machinery is tested with a stub instead of measurement (no time is measured).
