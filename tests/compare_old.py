"""CSV values are read exactly (pandas float_precision="round_trip"; the default parser can be off by
about 1e-13 relative for 17-digit values).

Compare the output of apfc.m/rpfc.m with the old iterative method (the reference for the corrected behaviour).

Usage (from the root of the matrix repository):
    python3 -B tests/compare_old.py <dir with <case>_apfc.csv, <case>_rpfc.csv> <inputs dir> [<inputs dir> ...]
        [--cases c1,c2,...] [--diff-csv FILE] [--rowsum-csv FILE]

For each case, the input <case>.mat is taken from the first inputs dir that contains it
(benchmark/inputs for the shared cases and tap_loop; the export of tests/export_synthetic.sh for
the synthetic cases). The old method's makeAPFC and makeRPFC are imported by path from the pinned
submodule external/iterative-method and run on load_results(<case>.mat), with QSF and QST of the tapped or
phase-shifting branches replaced by the line charging of rule Q3 version 2.0 (formula_shunts), with itermax = 10000 and
logging at INFO without a file handler. Their in-memory tables are compared, not the saved CSV files
(cleanRows() in the old code drops the rows of identical parallel branches).

Row mapping (old in-memory row -> matrix row):
    E<k>_PF / E<k>_PT, E<k>_QF / E<k>_QT (type E) -> L_{f.F-T} / L_{t.F-T}, F and T from the old row;
        parallel branches by order of appearance (matrix suffix _2, _3, ...)
    L<b>_PF, L<b>_QF (type L)            -> D_<b>     (demand)
    G<b>_PF, G<b>_QF (type L)            -> G_<b>     (generator with negative output, as a demand)
    S<b>_QF (type L)                     -> Bs_<b>    (inductive bus shunt)
    E<k>_QF_QF / E<k>_QT_QF (type L)     -> C2_<bus>  (branch consuming at both ends; from-end rows
                                            form the first half of the C2 block, to-end rows the second)
    E<k>.S<b>_QF (type L)                -> Tf_<b> if b is the from-bus of branch k in the input file,
                                            else Tt_<b> (inductive transformer shunt)
    rows of type G (sources) have no matrix row.
Column mapping: G<b> -> G_<b>; S<b> -> S_<b>; L<b> (negative demand as a source) -> D_<b>;
    E<k> (rpfc) -> sum of E_F_T_s[_n], E_F_T_f and E_F_T_t of the same branch.
Several generators at one bus: the old code names them G<b>.1, G<b>.2, ... (replaceDuplicates), the
matrix code G_<b>, G_<b>_2, ... (unique_labels). Both are reduced to the base name and paired by order
of appearance, for columns and for rows. Matrix demand-row labels with a uniqueness suffix
(C2_<b>_2, Tf_<b>_2, ...) are matched the same way.
Row sums: the matrix row sum is compared with the flow of the mapped old row (Total - Error).

The mapping follows reports/T0_scripts/d_compare.py of the T0 analysis, adapted to the in-memory tables.
"""
import sys
sys.dont_write_bytecode = True
import os, re, argparse, collections, logging, importlib.util
import numpy as np, pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
SUB = os.path.join(ROOT, "external", "iterative-method")
TOL = 1e-10
ITERMAX = 10000
DEFAULT_CASES = ["case1", "case9", "case14", "case24_ieee_rts", "case30", "case33", "case118", "tap_loop",
                 "syn9_neg_pd", "syn9_neg_pg", "syn9_neg_qd", "syn9_neg_qg", "syn9_case3_tap"]
g = lambda x: f"{x:.17g}"


def import_old():
    mods = {}
    for name in ["power_flow", "power_components"]:
        spec = importlib.util.spec_from_file_location(name, os.path.join(SUB, f"{name}.py"))
        mod = importlib.util.module_from_spec(spec); sys.modules[name] = mod; spec.loader.exec_module(mod)
        mods[name] = mod
    def _raise(*a, **k): raise RuntimeError("os.system called")
    os.system = _raise; mods["power_flow"].os.system = _raise
    return mods["power_flow"], mods["power_components"]


class Collect(logging.Handler):
    def __init__(self): super().__init__(logging.INFO); self.lines = []
    def emit(self, r): self.lines.append(r.getMessage())


def formula_shunts(branch, bus, mat):
    """A copy of the branch DataFrame of load_results in which QSF and QST are the line charging of the branch model
    of rule Q3 version 2.0 (ideal transformer with ratio tau = TAP e^(j SHIFT), then the pi circuit):
    QSF = (B/2) VM_f^2 / |TAP|^2 baseMVA and QST = (B/2) VM_t^2 baseMVA, on every branch with TAP != 0 and != 1 or
    SHIFT != 0. my_runpf.m exports the shunt flows of the real-ratio pi-equivalent of the tap instead (with B/2 at the
    from-end), which equal the line charging only on branches without a tap; the old method is run unchanged on the
    corrected input, so that both methods split the transformers alike. Returns (branch, rows replaced)."""
    import numpy as np
    import scipy.io as sio
    m = sio.loadmat(mat)
    bm = np.asarray(m["branch"], float); base = float(np.asarray(m["baseMVA"]).ravel()[0])
    assert len(bm) == len(branch)
    tap = np.abs(bm[:, 8]); tap[tap == 0] = 1.0
    rows = np.where((tap != 1) | (bm[:, 9] != 0))[0]
    out = branch.copy()
    if len(rows):
        vm = dict(zip(bus["BUS_I"].astype(int), bus["VM"].astype(float)))
        vf = np.array([vm[int(b)] for b in branch["F_BUS"]]); vt = np.array([vm[int(b)] for b in branch["T_BUS"]])
        half_b = bm[:, 4] / 2
        out.loc[out.index[rows], "QSF"] = (vf ** 2 * half_b / tap ** 2 * base)[rows]
        out.loc[out.index[rows], "QST"] = (vt ** 2 * half_b * base)[rows]
    return out, rows


def run_old(pf, pc, mat):
    root = logging.getLogger(); root.setLevel(logging.NOTSET)
    for h in list(root.handlers): root.removeHandler(h)
    root.addHandler(logging.NullHandler())
    h = Collect(); root.addHandler(h)
    bus, gen, branch, xy = pf.load_results(mat)
    branch, _ = formula_shunts(branch, bus, mat)       # rule Q3 version 2.0 (T16)
    fp, A = pc.makeAPFC(bus, gen, branch, None, ITERMAX, False)
    fq, R = pc.makeRPFC(bus, gen, branch, None, ITERMAX, False)
    root.removeHandler(h)
    st = {}
    for name in ["APFC", "RPFC"]:
        m = [re.search(rf"{name}: nr = \d+, iters = (\d+), DONE = (\w+)", l) for l in h.lines]
        m = [x for x in m if x][-1]
        st[name] = (m.group(2) == "True", int(m.group(1)))
    return A, R, fq, branch, st


def occ(keys):
    c = collections.Counter(); out = []
    for k in keys:
        c[k] += 1; out.append((k, c[k]))
    return out


def matrix_keys(labels, kind):
    """Keys (base label, occurrence) for matrix rows; C2 rows are split into from-end and to-end halves."""
    lab = list(labels)
    c2 = [i for i, l in enumerate(lab) if l.startswith("C2_")]
    for n, i in enumerate(c2):
        lab[i] = ("C2f_" if n < len(c2) // 2 else "C2t_") + lab[i][3:]
    return occ([base_label(l) for l in lab])


def base_label(l):
    """Remove a uniqueness suffix _n from a matrix row label or column name."""
    m = re.fullmatch(r"(L_\{[ft]\.\d+-\d+\}|(?:G|D|S|Bs|Tf|Tt|C2|C2f|C2t)_\d+)_\d+", l)
    return m.group(1) if m else l


def old_row_key(idx, r, kind, branch_from):
    q = "P" if kind == "apfc" else "Q"
    m = re.fullmatch(rf"E(\d+)_{q}([FT])", idx)
    if r.TYPE == "E" and m:
        return f"L_{{{m.group(2).lower()}.{int(r.F_BUS)}-{int(r.T_BUS)}}}"
    if r.TYPE != "L":
        return None
    if not idx.endswith(f"_{q}F"):
        return "skip"                       # the _PT/_QT copy of a demand row
    s = idx[:-3]
    if re.fullmatch(r"L\d+", s): return f"D_{int(r.F_BUS)}"
    if re.fullmatch(r"G\d+(\.\d+)?", s): return f"G_{int(r.F_BUS)}"
    if re.fullmatch(r"S\d+", s): return f"Bs_{int(r.F_BUS)}"
    m = re.fullmatch(r"E\d+_Q([FT])", s)
    if m: return ("C2f_" if m.group(1) == "F" else "C2t_") + str(int(r.F_BUS))
    m = re.fullmatch(r"E(\d+)\.S(\d+)", s)
    if m:
        k, b = int(m.group(1)), int(m.group(2))
        return ("Tf_" if b == branch_from[k - 1] else "Tt_") + str(b)
    return f"UNMAPPED:{idx}"


def colmap(old_cols, mcols, kind, fq):
    cm, unm = {}, []
    # source columns of generators, shunts and negative demands: pair by base name and order of appearance
    mocc = occ([base_label(x) for x in mcols]); mpos = {k: mcols[i] for i, k in enumerate(mocc)}
    obase = []
    for oc in old_cols:
        m = re.fullmatch(r"([GSL])(\d+)(?:\.\d+)?", oc)
        obase.append({"G": "G_", "S": "S_", "L": "D_"}[m.group(1)] + m.group(2) if m else None)
    oocc = occ(obase)
    if kind == "rpfc":
        E = fq[fq.TYPE == "E"]
        keys = {}
        cnt = collections.Counter()
        for i, r in E.iterrows():
            ft = (int(r.F_BUS), int(r.T_BUS)); case3 = r.QF < 0 and r.QT > 0
            if not case3: cnt[ft] += 1
            keys[i] = (ft, cnt[ft] if not case3 else 1, case3)
    for oc, key in zip(old_cols, oocc):
        if key[0] is not None:
            cand = [mpos[key]] if key in mpos else []
        elif kind == "rpfc" and re.fullmatch(r"E\d+", oc) and oc in keys:
            (f, t), n, _ = keys[oc]
            cand = [f"E_{f}_{t}_s" + ("" if n == 1 else f"_{n}"), f"E_{f}_{t}_f", f"E_{f}_{t}_t"]
        else:
            cand = []
        cand = [x for x in cand if x in mcols]
        if cand: cm[oc] = cand
        else: unm.append(oc)
    used = {x for v in cm.values() for x in v}
    return cm, unm, [x for x in mcols if x not in used]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir"); ap.add_argument("inputs", nargs="+")
    ap.add_argument("--cases"); ap.add_argument("--diff-csv"); ap.add_argument("--rowsum-csv")
    a = ap.parse_args()
    cases = a.cases.split(",") if a.cases else DEFAULT_CASES
    pf, pc = import_old()
    hdr = ("| case | file | old DONE | old iterations | rows mapped (old / matrix) | cells compared | max abs diff | cells above 1e-10 "
           "| max abs row sum - flow | unmapped old rows | unmapped matrix rows | old columns without counterpart | matrix columns without counterpart |")
    print(hdr); print("|" + "---|" * (hdr.count("|") - 1))
    diffs, rows_out, fails = [], [], 0
    for c in cases:
        mat = next((os.path.join(d, f"{c}.mat") for d in a.inputs if os.path.exists(os.path.join(d, f"{c}.mat"))), None)
        if mat is None:
            print(f"| {c} | | input not found | | | | | | | | | | |"); fails += 1; continue
        A, R, fq, branch, st = run_old(pf, pc, mat)
        branch_from = [int(x) for x in branch.F_BUS]
        for kind, T in [("apfc", A), ("rpfc", R)]:
            done, iters = st[kind.upper()]
            M = pd.read_csv(os.path.join(a.outdir, f"{c}_{kind}.csv"), float_precision="round_trip")
            mk = matrix_keys(M.Line, kind); mpos = {k: i for i, k in enumerate(mk)}
            mcols = list(M.columns[1:])
            src_end = list(T.columns).index("F_BUS")
            ocols = list(T.columns[:src_end])
            cm, unm_o, unm_m = colmap(ocols, mcols, kind, fq)
            okeys = [old_row_key(i, r, kind, branch_from) for i, r in T.iterrows()]
            ok_occ = occ(okeys)
            matched, nmapped, ncell, mx, nbad, rsmax = set(), 0, 0, 0.0, 0, 0.0
            unm_old = []
            for j, ((k, n), idx) in enumerate(zip(ok_occ, T.index)):
                if k in (None, "skip"): continue
                if (k, n) not in mpos:
                    unm_old.append(idx); continue
                i = mpos[(k, n)]; matched.add(i); nmapped += 1
                orow, mrow = T.iloc[j], M.iloc[i]
                flow = float(orow["Total"]) - float(orow["Error"])
                rs = float(mrow.iloc[1:].astype(float).sum()); rsd = abs(rs - flow); rsmax = max(rsmax, rsd)
                rows_out.append((c, kind, M.Line[i], idx, rs, flow, rsd))
                if kind == "rpfc" and not done:
                    continue                  # no reference for the cells
                for oc, mc in cm.items():
                    x = float(orow[oc]); y = float(mrow[mc].astype(float).sum()); d = abs(x - y)
                    ncell += 1; mx = max(mx, d)
                    if d > TOL:
                        nbad += 1; diffs.append((c, kind, M.Line[i], idx, oc, "+".join(mc), x, y, d))
            unm_mat = [M.Line[i] for i in range(len(M)) if i not in matched]
            ref = "" if (done or kind == "apfc") else " (no reference: old method did not finish)"
            fails += (nbad > 0) + (rsmax > TOL) + (len(unm_mat) > 0)
            cells = f"{ncell}" if ncell or done else "0" + ref
            print(f"| {c} | {kind} | {done} | {iters} | {nmapped} / {len(M)} | {cells} | {g(mx)} | {nbad} | {g(rsmax)} | "
                  f"{', '.join(unm_old) or 'none'} | {', '.join(unm_mat) or 'none'} | {', '.join(unm_o) or 'none'} | {', '.join(unm_m) or 'none'} |")
    print(f"\nChecks failed (cells above 1e-10, row sums above 1e-10, or unmapped matrix rows): {fails}")
    if a.diff_csv:
        pd.DataFrame(diffs, columns=["case", "file", "matrix row", "old row", "old col", "matrix cols", "old", "matrix", "absdiff"]).to_csv(a.diff_csv, index=False, float_format="%.17g")
    if a.rowsum_csv:
        pd.DataFrame(rows_out, columns=["case", "file", "matrix row", "old row", "matrix row sum", "flow", "absdiff"]).to_csv(a.rowsum_csv, index=False, float_format="%.17g")
    return 0 if fails == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
