"""Static predictor of the old iterative method's stall (T14e). Nothing is run through makeAPFC or makeRPFC.

The old method (external/iterative-method/power_components.py) is a sweep over the buses, not a numerical iteration:
it processes a bus once, when every branch in its inflow set is known, and never revisits it (T14d). So it can
finish only if every bus can be processed. From the initial sets this predictor builds the old method's dependency
digraph and finds the buses that can never be processed:
  - edge u -> k for every branch that leaves bus u (in the outflow set of u) and enters bus k (in the inflow set of
    k) and is not known initially;
  - orphan inflow: a branch in the inflow set of some bus, in no outflow set, and not known initially (it can never
    become known);
  - never processed: the buses of cyclic strongly connected components, the buses with an orphan inflow, and every
    bus reachable from them.
Label "predicted stall" if that set is not empty, else "can finish".

The sets are built exactly as the old code builds them, from the old method's own preprocessing of the DataFrames
(makeActiveFlows / makeReactiveFlows, which only rearrange the input; the corrected QSF of cases.correct_qsf does not
change the signs used here):
  RPFC (makeRPFC, lines 640-657): at bus k, a branch entering k with QT > 0 is an inflow (D[k]), otherwise an outflow
       (C[k]); a branch leaving k with QF > 0 is an outflow, otherwise an inflow. Known initially: the generator
       branches (TYPE G, line 697) and the branches with TYPE E, QF < 0 and QT > 0 (lines 702-705).
  APFC (makeAPFC, lines 435-452): the same with PT, PF and the sets B (inflow) and A (outflow). Known initially: the
       generator branches only (line 492).
  Buses: the sorted bus numbers of F_BUS and T_BUS without the first, 0 (lines 410-413 and 613-616).

What it cannot say: whether the old method raises an exception inside a sweep (the division by zero at lines 557 and
807 of four large combinations), how many sweeps it needs, or how long a sweep takes; a time-out of a combination
predicted "can finish" is a matter of size, not of a stall.

Usage (from the matrix repository root):
    python -B benchmark/timing/old_stall_predictor.py --combos case300:rpfc,case_ACTIVSg10k:apfc [--observed iterations.csv]
"""
import argparse
import csv
import os
import sys

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cases as C   # noqa: E402


def initial_sets(data, kind, branch=None):
    """(buses, out_sets, in_sets, known) as the old code builds them before its loop."""
    _, pc, _ = C.import_old()
    br = data["branch"] if branch is None else branch
    flows = (pc.makeReactiveFlows if kind == "rpfc" else pc.makeActiveFlows)(data["bus"], data["gen"], br)
    xf, xt = ("QF", "QT") if kind == "rpfc" else ("PF", "PT")
    F = flows["F_BUS"].to_numpy(); T = flows["T_BUS"].to_numpy()
    XF = flows[xf].to_numpy(float); XT = flows[xt].to_numpy(float)
    ids = list(flows.index)
    buses = sorted(set(F.tolist()) | set(T.tolist()))[1:]
    out_sets = {k: set() for k in buses}; in_sets = {k: set() for k in buses}
    for i, b in enumerate(ids):
        t, f = T[i], F[i]
        if t in in_sets:                                  # the branch enters bus t (an in-edge of t)
            (in_sets if XT[i] > 0 else out_sets)[t].add(b)
        if f in out_sets:                                 # the branch leaves bus f (an out-edge of f)
            (out_sets if XF[i] > 0 else in_sets)[f].add(b)
    typ = flows["TYPE"].to_numpy()
    known = {b for b, ty in zip(ids, typ) if ty == "G"}
    if kind == "rpfc":
        known |= {b for b, ty, a, c in zip(ids, typ, XF, XT) if ty == "E" and a < 0 and c > 0}
    return buses, out_sets, in_sets, known


def blocking(buses, out_sets, in_sets, known):
    pos = {b: i for i, b in enumerate(buses)}
    up = {}
    for u in buses:
        for b in out_sets[u]:
            up.setdefault(b, []).append(u)
    rows, cols, orphan = [], [], {}
    for k in buses:
        for b in in_sets[k]:
            if b in known:
                continue
            if b not in up:
                orphan.setdefault(k, []).append(b)
            for u in up.get(b, []):
                rows.append(pos[u]); cols.append(pos[k])
    n = len(buses)
    G = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    ncomp, lab = connected_components(G, directed=True, connection="strong")
    size = np.bincount(lab, minlength=ncomp)
    selfloop = np.zeros(n, bool)
    for r, c in zip(rows, cols):
        if r == c:
            selfloop[r] = True
    cyc = (size[lab] > 1) | selfloop
    seeds = set(np.where(cyc)[0]) | {pos[k] for k in orphan}
    blocked, stack = set(seeds), list(seeds)
    while stack:
        u = stack.pop()
        for v in G.indices[G.indptr[u]:G.indptr[u + 1]]:
            if v not in blocked:
                blocked.add(v); stack.append(v)
    return {"buses": n, "cyclic_buses": int(cyc.sum()), "orphan_buses": len(orphan),
            "orphan_branches": sorted({b for v in orphan.values() for b in v}),
            "never_processed": sorted(buses[i] for i in blocked)}


def predict(data, kind, branch=None):
    r = blocking(*initial_sets(data, kind, branch))
    r["predicted"] = len(r["never_processed"])
    r["label"] = "predicted stall" if r["predicted"] else "can finish"
    return r


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--combos", required=True, help="case:kind list")
    ap.add_argument("--observed", help="iterations.csv of a series: the observed old status per case and kind")
    a = ap.parse_args(argv)
    obs = {}
    if a.observed:
        for r in csv.DictReader(open(a.observed)):
            obs[(r["case"], r["kind"])] = r["old_status"]
    dirs = C.input_dirs()
    print("| case | kind | buses | buses in cyclic components | buses with an orphan inflow | predicted never processed | prediction | observed old status |")
    print("|---|---|---|---|---|---|---|---|")
    for c, kinds in C.parse_combos(a.combos):
        data = C.load(c, dirs)
        br = C.correct_qsf(data)[0]
        for k in kinds:
            r = predict(data, k, br)
            print(f"| {c} | {k} | {r['buses']} | {r['cyclic_buses']} | {r['orphan_buses']} | {r['predicted']} | {r['label']} | {obs.get((c, k), '')} |")


if __name__ == "__main__":
    main()
