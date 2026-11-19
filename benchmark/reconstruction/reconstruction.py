"""Reconstruction errors and properties of the matrix method (brief T14g, part B; T17 item 3), with the code of this
working tree, for case1, case9, case30, case118 and tap_loop, APFC and RPFC.

Reconstruction error: the row sums of the contribution table against the flows of the rows, reported separately for
  - the from-end rows of the carriers (P^f in APFC, Q^f in RPFC; in RPFC with the from-end rows of the Case-3 branches),
  - the to-end rows (P^t, Q^t; in RPFC with the to-end rows of the Case-3 branches),
  - the demand rows (D_b, G_b, Bs_b, C2_b, Tf_b, Tt_b),
  as the largest absolute error (MW for APFC, Mvar for RPFC) and relative to the largest absolute flow of the same rows.
  For the port and for Octave (apfc.m, rpfc.m of this working tree, outputs of tests/run_matrix.sh), with the flows of
  the port's preprocessing; the Octave tables must have the port's labels and order.
Residual: ||A X - B||_inf of the port's solve (and ||A X_self - B_self||_inf for the RPFC self-terms).
Properties: n_e (order of A), the number of table columns, nnz(N), density nnz(N)/n_e^2, the condition estimate of A
  (the port's condest_A: scipy.sparse.linalg.onenormest of A and of its inverse, a 1-norm estimate, as Matlab's
  condest), and the peak resident memory of one port run with checks off, read from memory.csv of a timing series
  (not rerun).
Writes reconstruction.csv, properties.csv and reconstruction_paper_cases.tex (booktabs tabular for case1, case9, case30
and case118) into the output folder. Exits with status 1 if a reconstruction error is above 1e-10.

Usage (from the matrix root):
    python -B benchmark/reconstruction/reconstruction.py <octave dir> <timing results dir> <output dir>
Inputs: benchmark/inputs/<case>.mat (exported flows), <octave dir>/<case>_{apfc,rpfc}.csv, <timing results dir>/memory.csv.
"""
import csv
import os
import sys
import warnings

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
MATRIX = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(MATRIX, "benchmark", "timing"))
sys.path.insert(0, os.path.join(MATRIX, "python"))
import cases as C   # noqa: E402
from pfc_matrix import (Case, Options, read_branch_params, apfc, rpfc, prepare_active, assemble_active,  # noqa: E402
                        solve_active, prepare_reactive, assemble_reactive, solve_reactive)

CASES = ["case1", "case9", "case30", "case118", "tap_loop"]
PAPER = ["case1", "case9", "case30", "case118"]
TOL = 1e-10
g = lambda x: "" if x is None else f"{x:.17g}"


def groups(kind, prep):
    """Row index groups of the table: demand rows, from-end rows, to-end rows, and the flows of all rows."""
    nd, nl = prep.nd, prep.nl
    if kind == "apfc":
        flows = np.r_[prep.Pf[nl:], prep.Pf[:nl], prep.Pt[:nl]]
        dem, fr, to = np.arange(nd), np.arange(nd, nd + nl), np.arange(nd + nl, nd + 2 * nl)
    else:
        hf, ht = np.asarray(prep.case3_hasf, bool), np.asarray(prep.case3_hast, bool)
        n3f, n3t = int(hf.sum()), int(ht.sum())
        flows = np.r_[prep.Qf[nl:], prep.Qf[:nl], prep.Qt[:nl], -prep.case3_QF[hf], prep.case3_QT[ht]]
        base = nd + 2 * nl
        dem = np.arange(nd)
        fr = np.r_[np.arange(nd, nd + nl), np.arange(base, base + n3f)]
        to = np.r_[np.arange(nd + nl, nd + 2 * nl), np.arange(base + n3f, base + n3f + n3t)]
    return flows, {"demand": dem, "from": fr, "to": to}


def errors(rowsums, flows, grp):
    out = {}
    for name, idx in grp.items():
        if len(idx) == 0:
            out[name] = (None, None)
            continue
        e = float(np.abs(rowsums[idx] - flows[idx]).max())
        scale = float(np.abs(flows[idx]).max())
        out[name] = (e, e / scale if scale > 0 else None)
    return out


def main(octd, timing, outd):
    os.makedirs(outd, exist_ok=True)
    mem = {}
    with open(os.path.join(timing, "memory.csv")) as fh:
        for r in csv.DictReader(fh):
            if r["method"] == "port" and r["checks"] == "off" and r["status"] == "ok":
                mem[(r["case"], r["kind"])] = (int(r["ru_maxrss_after"]), r["unit"])
    rec, prop = [], []
    worst = 0.0
    for c in CASES:
        data = C.load(c, C.input_dirs())
        case = Case.from_frames(data["bus"], data["gen"], data["branch"], c, params=read_branch_params(data["path"]))
        for kind in ("apfc", "rpfc"):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                t, info = (apfc if kind == "apfc" else rpfc)(case, Options(checks=True))
                if kind == "apfc":
                    prep = prepare_active(case, 1e-6); sysm = assemble_active(prep); X = solve_active(sysm)
                    res = float(np.abs(sysm.A @ X - sysm.B.toarray()).max()) if X.size else 0.0
                    res_self = None
                else:
                    prep = prepare_reactive(case, 1e-6); sysm = assemble_reactive(prep); X, Xs = solve_reactive(sysm)
                    res = float(np.abs(sysm.A @ X - sysm.B.toarray()).max()) if X.size else 0.0
                    res_self = float(np.abs(sysm.A @ Xs - sysm.B_self.toarray()).max()) if Xs.size else 0.0
            flows, grp = groups(kind, prep)
            assert len(flows) == len(t)
            ep = errors(t.to_numpy().sum(axis=1), flows, grp)
            o = pd.read_csv(os.path.join(octd, f"{c}_{kind}.csv"), float_precision="round_trip", dtype={"Line": str}).set_index("Line")
            lab = list(o.index) == list(t.index) and list(o.columns) == list(t.columns)
            agree = float(np.abs(o.to_numpy() - t.to_numpy()).max()) if lab else None
            eo = errors(o.to_numpy().sum(axis=1), flows, grp) if lab else {k: (None, None) for k in grp}
            for side in (ep, eo):
                for e, _ in side.values():
                    if e is not None:
                        worst = max(worst, e)
            q = "P" if kind == "apfc" else "Q"
            unit = "MW" if kind == "apfc" else "Mvar"
            row = {"case": c, "kind": kind, "unit": unit, "quantity_from": f"{q}^f", "quantity_to": f"{q}^t",
                   "port_from_abs": ep["from"][0], "port_from_rel": ep["from"][1], "port_to_abs": ep["to"][0], "port_to_rel": ep["to"][1],
                   "port_demand_abs": ep["demand"][0], "port_demand_rel": ep["demand"][1],
                   "octave_from_abs": eo["from"][0], "octave_to_abs": eo["to"][0], "octave_demand_abs": eo["demand"][0],
                   "octave_labels_equal": lab, "port_octave_max_abs_diff": agree,
                   "residual_AX_minus_B_inf": res, "residual_self_inf": res_self,
                   "largest_abs_flow": float(np.abs(flows).max())}
            rec.append(row)
            m = mem.get((c, kind))
            prop.append({"case": c, "kind": kind, "n_e": info["n_e"], "table_rows": len(t), "table_columns": t.shape[1],
                         "nnz_N": info["nnz_N"], "density_nnz_N_over_n_e_squared": info["nnz_N"] / info["n_e"] ** 2,
                         "condest_A_1norm": info["condest_A"], "acyclic": info["acyclic"], "spectral_radius": info["spectral_radius"],
                         "peak_rss_port_checks_off": m[0] if m else None, "peak_rss_unit": m[1] if m else None})
    for name, rows in (("reconstruction.csv", rec), ("properties.csv", prop)):
        with open(os.path.join(outd, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            for r in rows:
                w.writerow({k: (g(v) if isinstance(v, float) else v) for k, v in r.items()})
    # booktabs table of the four paper cases: the four reconstruction errors (port), absolute
    tex = ["% Generated by benchmark/reconstruction/reconstruction.py (do not edit). Largest absolute error of the row sums",
           "% against the flows: P^f, P^t in MW (APFC), Q^f, Q^t in Mvar (RPFC); port, rule Q3 version 2.0.",
           "\\begin{tabular}{@{}lrrrr@{}}", "\\toprule",
           "Case & $P^f$ (MW) & $P^t$ (MW) & $Q^f$ (Mvar) & $Q^t$ (Mvar) \\\\", "\\midrule"]
    by = {(r["case"], r["kind"]): r for r in rec}
    fmt = lambda x: "0" if x == 0 else f"{x:.1e}".replace("e-", "\\text{e-}")
    for c in PAPER:
        a, r = by[(c, "apfc")], by[(c, "rpfc")]
        tex.append(f"{c} & ${fmt(a['port_from_abs'])}$ & ${fmt(a['port_to_abs'])}$ & ${fmt(r['port_from_abs'])}$ & ${fmt(r['port_to_abs'])}$ \\\\")
    tex += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(outd, "reconstruction_paper_cases.tex"), "w").write("\n".join(tex) + "\n")
    print(f"largest reconstruction error (port and Octave, all row groups): {worst:.17g}")
    return 0 if worst <= TOL else 1


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:4]))
