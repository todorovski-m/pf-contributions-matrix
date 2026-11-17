"""Dead-end branches (terminal classification, brief T15b): comparison of apfc.m with the old iterative method.

The old method is a reference for active dead ends where it finishes (report T15a): it treats the dead-end branch
as a branch whose zero end carries nothing, and the loss row of that branch is shared like a demand at the other
bus. apfc.m represents the same branch by one terminal demand (row C2_<b>). This script runs the old method's
makeAPFC (pinned submodule, in memory, itermax 10000) on the exported input of each case and compares its table
with the apfc.m output, using the unchanged helper functions of tests/compare_old.py:
    every row and column except those of the dead-end branch is mapped as in compare_old.py;
    the old rows E<k>_PF and E<k>_PT of the dead-end branch (old ID E<k>): the row at the zero end must be zero,
    and the row at the other end, as a demand (E<k>_PF if the dead end is at the to side; minus E<k>_PT if it is at
    the from side, because the old to-end row holds the delivered flow), is compared with the row C2_<b>.
The old method is not a reference for reactive dead ends (it does not finish, or its row sums differ from the
flows; report T15a), so RPFC is not compared here; the reactive results are checked by the equivalence tests of
tests/unit/test_terminal.m.

Usage (from the root of the matrix repository):
    python3 -B tests/dead_end_equivalence.py <dir with <case>_apfc.csv> <inputs dir> [--cases c1,c2,...]
Default cases: syn9_dead_charge_from, syn9_dead_charge_to, syn9_dead_shunt_to. Exit status 0 if every case passes.
"""
import sys
sys.dont_write_bytecode = True
import os, re, argparse, logging, importlib.util
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TOL = 1e-10
ITERMAX = 10000
DEFAULT_CASES = ["syn9_dead_charge_from", "syn9_dead_charge_to", "syn9_dead_shunt_to"]
DEAD_BUS = 10
g = lambda x: f"{x:.17g}"


def load_compare_old():
    spec = importlib.util.spec_from_file_location("compare_old", os.path.join(HERE, "compare_old.py"))
    co = importlib.util.module_from_spec(spec); spec.loader.exec_module(co)
    return co


def run_old_apfc(co, pf, pc, mat):
    root = logging.getLogger(); root.setLevel(logging.NOTSET)
    for h in list(root.handlers): root.removeHandler(h)
    root.addHandler(logging.NullHandler())
    h = co.Collect(); root.addHandler(h)
    bus, gen, branch, _ = pf.load_results(mat)
    _, A = pc.makeAPFC(bus, gen, branch, None, ITERMAX, False)
    root.removeHandler(h)
    m = [x for x in (re.search(r"APFC: nr = \d+, iters = (\d+), DONE = (\w+)", l) for l in h.lines) if x][-1]
    return A, branch, m.group(2) == "True", int(m.group(1))


def compare(co, A, branch, M):
    """Returns a dict with the numbers of one case."""
    bf = [int(x) for x in branch.F_BUS]
    k = [i for i in range(len(branch)) if DEAD_BUS in (int(branch.F_BUS.iloc[i]), int(branch.T_BUS.iloc[i]))]
    assert len(k) == 1, "one dead-end branch expected"
    k = k[0]; eid = f"E{k + 1}"
    dead_from = int(branch.F_BUS.iloc[k]) == DEAD_BUS
    other = int(branch.T_BUS.iloc[k]) if dead_from else int(branch.F_BUS.iloc[k])
    own_row, zero_row, sign = (f"{eid}_PT", f"{eid}_PF", -1.0) if dead_from else (f"{eid}_PF", f"{eid}_PT", 1.0)
    target = f"C2_{other}"
    mk = co.matrix_keys(M.Line, "apfc"); mpos = {kk: i for i, kk in enumerate(mk)}
    mcols = list(M.columns[1:])
    ocols = list(A.columns[:list(A.columns).index("F_BUS")])
    cm, unm_o, unm_m = co.colmap(ocols, mcols, "apfc", None)
    okeys = [None if (idx.startswith(eid + "_")) else co.old_row_key(idx, r, "apfc", bf) for idx, r in A.iterrows()]
    matched, mx, rsmax, unm_old = set(), 0.0, 0.0, []
    for j, ((kk, n), idx) in enumerate(zip(co.occ(okeys), A.index)):
        if kk in (None, "skip"):
            continue
        if (kk, n) not in mpos:
            unm_old.append(idx); continue
        i = mpos[(kk, n)]; matched.add(i)
        orow, mrow = A.iloc[j], M.iloc[i]
        rsmax = max(rsmax, abs(float(mrow.iloc[1:].astype(float).sum()) - (float(orow["Total"]) - float(orow["Error"]))))
        for oc, mc in cm.items():
            mx = max(mx, abs(float(orow[oc]) - float(mrow[mc].astype(float).sum())))
    # the dead-end branch
    it = [i for i, l in enumerate(M.Line) if l == target]
    own = float("inf"); zero = float("inf")
    if len(it) == 1:
        i = it[0]; matched.add(i)
        orow, mrow = A.loc[own_row], M.iloc[i]
        own = max(abs(sign * float(orow[oc]) - float(mrow[mc].astype(float).sum())) for oc, mc in cm.items())
        zero = max(abs(float(A.loc[zero_row, oc])) for oc in ocols)
    unm_mat = [M.Line[i] for i in range(len(M)) if i not in matched]
    ok = mx <= TOL and rsmax <= TOL and own <= TOL and zero <= TOL and not unm_mat and not unm_old and not unm_o and not unm_m
    return {"old_id": eid, "dead_end_side": "from" if dead_from else "to", "row": f"{own_row} ({'minus ' if sign < 0 else ''}old) -> {target}",
            "rows_mapped": len(matched), "rows": len(M), "max_other": mx, "old_rowsum": rsmax, "max_own": own,
            "zero_row": zero, "unmapped": unm_old + unm_mat + unm_o + unm_m, "ok": ok}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir"); ap.add_argument("inputs", nargs="+"); ap.add_argument("--cases")
    a = ap.parse_args()
    cases = a.cases.split(",") if a.cases else DEFAULT_CASES
    co = load_compare_old()
    pf, pc = co.import_old()
    print("| case | old DONE | old iterations | old ID of the dead-end branch | dead end at | dead-end row mapping | rows mapped (matrix) "
          "| max abs diff, other rows | max abs diff, dead-end row against C2 row | max abs value of the old row at the zero end "
          "| max abs old row sum - flow | unmapped | pass |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    fails = 0
    for c in cases:
        mat = next((os.path.join(d, f"{c}.mat") for d in a.inputs if os.path.exists(os.path.join(d, f"{c}.mat"))), None)
        A, branch, done, iters = run_old_apfc(co, pf, pc, mat)
        M = pd.read_csv(os.path.join(a.outdir, f"{c}_apfc.csv"), float_precision="round_trip")
        r = compare(co, A, branch, M)
        ok = r["ok"] and done
        fails += not ok
        print(f"| {c} | {done} | {iters} | {r['old_id']} | {r['dead_end_side']} side | {r['row']} | {r['rows_mapped']} / {r['rows']} | "
              f"{g(r['max_other'])} | {g(r['max_own'])} | {g(r['zero_row'])} | {g(r['old_rowsum'])} | {', '.join(r['unmapped']) or 'none'} | {ok} |")
    print(f"\nCases failing: {fails}")
    return 0 if fails == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
