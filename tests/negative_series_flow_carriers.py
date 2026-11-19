"""Write tests/data/negative_series_flow_carriers.csv: per case, the carrier branches with a negative series reactive
flow as the old iterative method sees them (makeReactiveFlows of the pinned submodule): QF1 = QF + QSF < -1e-6 or
QT1 = QT - QST < -1e-6, branches that inject at both ends (Case 3: QF < 0, QT > 0) excluded. QSF and QST of the tapped or
phase-shifting branches are the line charging of rule Q3 version 2.0 (compare_old.formula_shunts), as for the
comparison with the old method. The bus pair is the old method's (after its flip of branches with both flows
negative). Every case has one row without a branch (the coverage marker). rpfc.m and the port report the same set in
diag.negative_series_flow (tests/unit/test_diag.m, python/tests/test_checks.py).

Usage (from the root of the matrix repository):
    python3 -B tests/negative_series_flow_carriers.py <inputs dir> [<inputs dir> ...] > tests/data/negative_series_flow_carriers.csv
"""
import sys
sys.dont_write_bytecode = True
import os, logging
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_old as co

CASES = ["case1", "case118", "case14", "case24_ieee_rts", "case30", "case33", "case9", "syn9_neg_pd", "syn9_neg_pg",
         "syn9_neg_qd", "syn9_neg_qg", "tap_loop", "syn9_case3_tap"]
dirs = sys.argv[1:]
pf, pc = co.import_old()
logging.getLogger().addHandler(logging.NullHandler())
print("case,from_bus,to_bus")
for c in CASES:
    mat = next(os.path.join(d, f"{c}.mat") for d in dirs if os.path.exists(os.path.join(d, f"{c}.mat")))
    bus, gen, branch, xy = pf.load_results(mat)
    branch, _ = co.formula_shunts(branch, bus, mat)
    fq = pc.makeReactiveFlows(bus, gen, branch)
    E = fq[fq.TYPE == "E"]
    neg = E[((E.QF1 < -1e-6) | (E.QT1 < -1e-6)) & ~((E.QF < 0) & (E.QT > 0))]
    for _, r in neg.iterrows():
        print(f"{c},{int(r.F_BUS)},{int(r.T_BUS)}")
    print(f"{c},,")
