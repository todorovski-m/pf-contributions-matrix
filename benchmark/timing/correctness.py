"""Correctness of the port on the exported (large) cases (brief T14c, B2), per (kind, case) combination:
    the port's row sums against the flows of the rows, and determinism (two runs, SHA-256 of the CSV of write_csv);
    the old method in memory (makeAPFC/makeRPFC of the pinned submodule, iteration cap 10000, corrected QSF where
    tap != 1 and BR_B != 0, as in the T14a harness), with a time-out per run and a total budget for all old runs,
    and the cell comparison with the port.
No timing is measured. Run from the root of the matrix repository:
    python -B benchmark/timing/correctness.py --runs apfc:case300,rpfc:case300,... [--timeout 1800] [--budget 10800]
        [--out-md FILE] [--csv-dir DIR]

Row mapping (old in-memory row -> port row). The unchanged helpers of tests/compare_old.py are used (old_row_key,
base_label, occ, colmap), except for the rows that the terminal classification creates:
    APFC: a branch with one zero terminal (patterns (+,0), (0,-); T15b) is one terminal demand C2_<b>; the old
        method keeps it as a branch, and its non-zero row (E<k>_PF for (+,0); minus E<k>_PT for (0,-), because the
        old to-end row holds the delivered flow) is compared with the C2 row, and its zero row must be zero.
        A branch fed from both ends (pattern (+,-); T14c) is two terminal demands; the old method is not a reference
        for its two rows (each is twice the flow, decision on the T14c report), so they are excluded and counted.
        The C2 rows of the port are identified by position (the order rule: from-terminals of the (+,0) and (+,-)
        branches, then to-terminals of the (0,-) and (+,-) branches, each in branch order).
    RPFC: the old method is not a reference for branches with one zero reactive terminal (T15a, T15b); their rows
        and their source columns are excluded and counted. The other Case-2 terminal demands are identified by
        position (from-terminals of all Case-2 members, then their to-terminals), the Case-3 rows by position at the
        end of the table.
A combination passes if the port's row sums equal the flows within 1e-10 (or the scale-appropriate value reported)
and, where the old method finishes, every compared cell is within 1e-10 and the old row sums equal the flows.
"""
import os, sys, time, hashlib, argparse, warnings, collections
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "..", "python"))
import numpy as np
import timing as T
import cases as C

TOL, CELL_TOL = 1e-6, 1e-10
g = lambda x: "" if x is None else (f"{x:.17g}" if isinstance(x, float) else str(x))


def sgn(x):
    return (x > TOL).astype(int) - (x < -TOL).astype(int)


def port_table(data, kind):
    from pfc_matrix import Case, read_branch_params, apfc, rpfc, prepare_active, prepare_reactive
    case = Case.from_frames(data["bus"], data["gen"], data["branch"], data["case"], params=read_branch_params(data["path"]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t, info = (apfc if kind == "apfc" else rpfc)(case)
    if kind == "apfc":
        p = prepare_active(case, TOL); nl = p.nl
        flows = np.r_[p.Pf[nl:], p.Pf[:nl], p.Pt[:nl]]
    else:
        p = prepare_reactive(case, TOL); nl = p.nl
        hf, ht = np.asarray(p.case3_hasf, bool), np.asarray(p.case3_hast, bool)
        flows = np.r_[p.Qf[nl:], p.Qf[:nl], p.Qt[:nl], -p.case3_QF[hf], p.case3_QT[ht]]
    return t, info, flows


def csv_hash(table, path):
    from pfc_matrix import write_csv
    write_csv(table, path)
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def classify(data, kind):
    """Per exported branch row k: kept after the deletion, the sign pattern by tol, and the bus numbers."""
    br = data["branch"]
    x = br["PF" if kind == "apfc" else "QF"].to_numpy(float); y = br["PT" if kind == "apfc" else "QT"].to_numpy(float)
    keep = ~((np.abs(x) < TOL) & (np.abs(y) < TOL))
    return keep, sgn(x), sgn(y), x, y, br["F_BUS"].to_numpy(int), br["T_BUS"].to_numpy(int)


def port_keys(data, kind, t):
    """Key of every port row (None: excluded) and the number of excluded rows, by the rules of the docstring."""
    co = C.import_old()[2]
    keep, sf, st, x, y, F, Tb = classify(data, kind)
    rows = list(t.index); keys = [None] * len(rows); excluded = 0
    special = set()                                  # positions handled here
    tail, c3order = [], []
    one = keep & (((sf != 0) & (st == 0)) | ((sf == 0) & (st != 0)))
    if kind == "apfc":
        c2 = keep & (sf > 0) & (st < 0)
        frm = [k for k in range(len(F)) if keep[k] and ((sf[k] > 0 and st[k] == 0) or c2[k])]
        to = [k for k in range(len(F)) if keep[k] and ((sf[k] == 0 and st[k] < 0) or c2[k])]
        pos = [i for i, r in enumerate(rows) if r.startswith("C2_")]
        assert len(pos) == len(frm) + len(to), "C2 rows do not match the terminal classification"
        for i, (k, side) in zip(pos, [(k, "f") for k in frm] + [(k, "t") for k in to]):
            special.add(i)
            if c2[k]:
                excluded += 1
            else:
                keys[i] = ("ZT", k)
    else:
        # Case-2 members (never flipped): (+,-) with both terminals non-zero, and the one-zero (+,0), (0,-)
        c2m = [k for k in range(len(F)) if keep[k] and ((sf[k] > 0 and st[k] < 0) or (sf[k] > 0 and st[k] == 0) or (sf[k] == 0 and st[k] < 0))]
        frm = [(k, F[k]) for k in c2m if sf[k] != 0]
        to = [(k, Tb[k]) for k in c2m if st[k] != 0]
        pos = [i for i, r in enumerate(rows) if r.startswith("C2_")]
        assert len(pos) == len(frm) + len(to), "C2 rows do not match the Case-2 members"
        cnt = collections.Counter()
        for i, (k, b), side in zip(pos, frm + to, ["f"] * len(frm) + ["t"] * len(to)):
            special.add(i)
            if one[k]:
                excluded += 1
            else:
                key = ("C2f_" if side == "f" else "C2t_") + str(b); cnt[key] += 1; keys[i] = (key, cnt[key])
        # Case-3 members, appended at the end: from-rows of the members with a from-end source, then the to-rows
        c3m = [k for k in range(len(F)) if keep[k] and ((sf[k] < 0 and st[k] > 0) or (sf[k] < 0 and st[k] == 0) or (sf[k] == 0 and st[k] > 0))]
        n3 = sum(sf[k] != 0 for k in c3m) + sum(st[k] != 0 for k in c3m)
        tail = list(range(len(rows) - n3, len(rows)))
        c3order = [k for k in c3m if sf[k] != 0] + [k for k in c3m if st[k] != 0]
        excluded += sum(bool(one[k]) for k in c3order)
    # every other row: the base label and its occurrence (as compare_old.matrix_keys, without its C2 halves)
    occ = collections.Counter()
    tpos = {i: n for n, i in enumerate(tail)}
    for i, r in enumerate(rows):
        if i in special:                              # C2 rows: keyed above (or excluded)
            continue
        if i in tpos and one[c3order[tpos[i]]]:       # Case-3 row of a branch with one zero terminal: excluded
            continue
        b = co.base_label(r); occ[b] += 1; keys[i] = (b, occ[b])
    return keys, excluded


def old_keys(data, kind, A):
    """Key of every old row, the sign to apply, and the rows that must be zero (APFC zero ends)."""
    co = C.import_old()[2]
    keep, sf, st, x, y, F, Tb = classify(data, kind)
    bf = [int(v) for v in data["branch"].F_BUS]
    q = "P" if kind == "apfc" else "Q"
    one = keep & (((sf != 0) & (st == 0)) | ((sf == 0) & (st != 0)))
    keys, signs, zero = [], [], []
    occ = collections.Counter(); cnt = collections.Counter()
    for idx, r in A.iterrows():
        k = None
        if idx.startswith("E"):
            num = idx[1:].split("_")[0].split(".")[0]
            k = int(num) - 1 if num.isdigit() else None
        key, s = None, 1.0
        if kind == "apfc" and k is not None and r.TYPE == "E" and keep[k] and (one[k] or (sf[k] > 0 and st[k] < 0)):
            if one[k]:
                own = f"E{k + 1}_PF" if sf[k] > 0 else f"E{k + 1}_PT"
                if idx == own:
                    key, s = ("ZT", k), (1.0 if sf[k] > 0 else -1.0)
                else:
                    zero.append(idx)
        elif kind == "rpfc" and k is not None and keep[k] and one[k]:
            key = None                                    # not a reference (one zero reactive terminal)
        else:
            base = co.old_row_key(idx, r, kind, bf)
            if base not in (None, "skip"):
                if base.startswith(("C2f_", "C2t_")):
                    cnt[base] += 1; key = (base, cnt[base])
                else:
                    occ[base] += 1; key = (base, occ[base])
        keys.append(key); signs.append(s)
    return keys, signs, zero


def compare(data, kind, A, fq, t):
    co = C.import_old()[2]
    pk, n_excl = port_keys(data, kind, t)
    ok_, sg_, zero = old_keys(data, kind, A)
    ppos = {k: i for i, k in enumerate(pk) if k is not None}
    mcols = list(t.columns); ocols = list(A.columns[:list(A.columns).index("F_BUS")])
    cm, unm_o, unm_m = co.colmap(ocols, mcols, kind, fq)
    excl_cols = []
    if kind == "rpfc":
        keep, sf, st, x, y, F, Tb = classify(data, kind)
        one = keep & (((sf != 0) & (st == 0)) | ((sf == 0) & (st != 0)))
        for k in np.where(one)[0]:
            oc = f"E{k + 1}"
            if oc in cm:
                excl_cols += cm.pop(oc)
            elif oc in unm_o:
                unm_o.remove(oc)
    unm_m = [c for c in unm_m if c not in excl_cols]
    M = t.to_numpy(); colpos = {c: j for j, c in enumerate(mcols)}
    matched, mx, rsmax, worst, rs_flow = set(), 0.0, 0.0, None, 0.0
    psum = M.sum(axis=1)
    for j, (key, s) in enumerate(zip(ok_, sg_)):
        if key is None or key not in ppos:
            continue
        i = ppos[key]; matched.add(i)
        orow = A.iloc[j]
        rsmax = max(rsmax, abs(float(orow["Error"])))       # the old row error: Total - flow
        rs_flow = max(rs_flow, abs(float(psum[i]) - s * (float(orow["Total"]) - float(orow["Error"]))))   # port row sum - old flow
        for oc, mc in cm.items():
            d = abs(s * float(orow[oc]) - float(sum(M[i, colpos[c]] for c in mc)))
            if d > mx:
                mx, worst = d, (A.index[j], t.index[i], oc)
    zmax = max([abs(float(A.loc[z, oc])) for z in zero for oc in ocols] or [0.0])
    unmapped = [t.index[i] for i, k in enumerate(pk) if k is not None and i not in matched]
    return {"max_cell": mx, "worst": worst, "old_row_error": rsmax, "port_rowsum_vs_old_flow": rs_flow, "zero_rows": zmax, "excluded_rows": n_excl,
            "excluded_columns": len(excl_cols), "unmapped_port_rows": unmapped, "unmapped_old_columns": unm_o,
            "unmapped_port_columns": unm_m}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True); ap.add_argument("--timeout", type=float, default=1800.0)
    ap.add_argument("--budget", type=float, default=10800.0); ap.add_argument("--out-md")
    ap.add_argument("--csv-dir", default="/tmp/correctness_csv"); ap.add_argument("--no-old", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.csv_dir, exist_ok=True)
    dirs = C.input_dirs(None, None)
    head = ("| kind | case | port rows x columns | port max abs (row sum - flow) | largest abs flow | two runs byte-identical | "
            "old method | old iterations | old elapsed s | max abs cell difference | at (old row, port row, old column) | "
            "max abs old row error (Total - flow) | max abs old value at an APFC zero end | excluded rows | excluded columns | "
            "unmapped port rows | pass |")
    lines = [head, "|" + "---|" * (head.count("|") - 1)]
    print(lines[0]); print(lines[1])
    spent = 0.0
    for kc in a.runs.split(","):
        kind, name = kc.split(":")
        data = C.load(name, dirs)
        t, info, flows = port_table(data, kind)
        rs = float(np.abs(t.to_numpy().sum(axis=1) - flows).max())
        path = os.path.join(a.csv_dir, f"{name}_{kind}.csv")
        h1 = csv_hash(t, path)
        t2, _, _ = port_table(data, kind)
        h2 = csv_hash(t2, path); del t2
        row = {"rows": f"{t.shape[0]} x {t.shape[1]}", "rs": rs, "fmax": float(np.abs(flows).max()), "det": h1 == h2}
        status, iters, el, cmpd = "not run (option --no-old)", None, None, None
        if not a.no_old:
            if spent >= a.budget:
                status = "not run (budget)"
            else:
                branch_old, _, _, _ = C.correct_qsf(data)
                t0 = time.monotonic()
                try:
                    fq, A, done, iters = T.old_with_iterations(data, kind, 10000, branch_old, min(a.timeout, a.budget - spent))
                    el = time.monotonic() - t0
                    status = "finished" if done else "did not finish (iteration cap)"
                    cmpd = compare(data, kind, A, fq, t)
                except T.RunTimeout:
                    el = time.monotonic() - t0
                    status = "did not finish (time-out)" if el >= a.timeout - 1 else "did not finish (budget)"
                spent += el
        ok = rs <= CELL_TOL and row["det"]
        if status == "finished":
            ok = ok and cmpd["max_cell"] <= CELL_TOL and cmpd["old_row_error"] <= CELL_TOL and cmpd["zero_rows"] <= CELL_TOL \
                and not cmpd["unmapped_port_rows"]
        c = cmpd or {}
        line = (f"| {kind} | {name} | {row['rows']} | {g(rs)} | {g(row['fmax'])} | {row['det']} | {status} | {g(iters)} | "
                f"{'' if el is None else f'{el:.0f}'} | {g(c.get('max_cell')) if status == 'finished' else ''} | "
                f"{c.get('worst') if status == 'finished' else ''} | {g(c.get('old_row_error')) if status == 'finished' else ''} | {g(c.get('zero_rows')) if status == 'finished' else ''} | "
                f"{c.get('excluded_rows', '')} | {c.get('excluded_columns', '')} | "
                f"{', '.join(c.get('unmapped_port_rows', [])[:5]) + (' ...' if len(c.get('unmapped_port_rows', [])) > 5 else '') if c else ''} | {ok} |")
        print(line, flush=True); lines.append(line)
        del t
    if a.out_md:
        open(a.out_md, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
