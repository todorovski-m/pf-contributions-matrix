"""Table of negative reactive attenuation and negative contributions per case (for README.md).

Usage (from the root of the matrix repository):
    python3 -B tests/negative_attenuation_table.py <dir with <case>_rpfc.csv> <inputs dir> [case ...]

For each case:
  - carrier branches with a negative series reactive flow (negative attenuation Qe/Qs): computed from the
    series flows QF1 = QF + QSF and QT1 = QT - QST of the old iterative method (makeReactiveFlows of the
    pinned submodule; QSF and QST of the tapped or phase-shifting branches replaced by the line charging of
    rule Q3 version 2.0, compare_old.formula_shunts), branches that inject at both ends (Case 3: QF < 0, QT > 0) excluded, since they
    are sources, not carriers; the same set is reported by rpfc.m in diag.negative_series_flow;
  - the most negative attenuation QT1/QF1 among them;
  - the negative cells (below -1e-9) of the rpfc.m output, excluding the own terminal injection of a
    Case-3 branch in its own row (negative by design), and the smallest such value.
Adapted from the analysis scripts of the revision work (negative_series_flows.py, negative_cells.py).
"""
import sys
sys.dont_write_bytecode = True
import os, re, logging
import pandas as pd, scipy.io as sio
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_old as co

out, inputs = sys.argv[1], sys.argv[2]
cases = sys.argv[3:] or ["case1", "case9", "case14", "case24_ieee_rts", "case30", "case33", "case118", "tap_loop"]
pf, pc = co.import_old()
logging.getLogger().addHandler(logging.NullHandler())
print("| case | carrier branches with negative attenuation | of which transformers | most negative attenuation | "
      "negative contribution cells | smallest negative contribution (Mvar) |")
print("|---|---|---|---|---|---|")
for c in cases:
    bus, gen, branch, xy = pf.load_results(os.path.join(inputs, f"{c}.mat"))
    branch, _ = co.formula_shunts(branch, bus, os.path.join(inputs, f"{c}.mat"))     # rule Q3 version 2.0 (T16)
    raw = sio.loadmat(os.path.join(inputs, f"{c}.mat"))["branch"]
    fq = pc.makeReactiveFlows(bus, gen, branch)
    E = fq[fq.TYPE == "E"]
    neg = E[((E.QF1 < -1e-6) | (E.QT1 < -1e-6)) & ~((E.QF < 0) & (E.QT > 0))]
    ntr = sum(1 for i in neg.index if raw[int(i[1:]) - 1, 8] != 0)
    att = (neg.QT1 / neg.QF1).min() if len(neg) else None
    M = pd.read_csv(os.path.join(out, f"{c}_rpfc.csv"), float_precision="round_trip").set_index("Line")
    v = M.to_numpy(float); cells = []
    for i, r in enumerate(M.index):
        for j, col in enumerate(M.columns):
            if v[i, j] < -1e-9:
                mr = re.fullmatch(r"L_\{[ft]\.(\d+)-(\d+)\}(?:_\d+)?", r); mc = re.fullmatch(r"E_(\d+)_(\d+)_[ft](?:_\d+)?", col)
                if not (mr and mc and mr.groups() == mc.groups()):
                    cells.append(v[i, j])
    print(f"| {c} | {len(neg)} | {ntr} | {f'{att:.3f}' if att is not None else '-'} | {len(cells)} | "
          f"{f'{min(cells):.3f}' if cells else '-'} |")
