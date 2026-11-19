"""CSV values are read exactly (pandas float_precision="round_trip"; the default parser can be off by
about 1e-13 relative for 17-digit values).

Compare the output of the current apfc.m/rpfc.m with paper_results/ (the frozen reference of the submitted paper).

Usage (from the root of the matrix repository):
    python3 tests/compare_legacy.py <dir with new <case>_apfc.csv and <case>_rpfc.csv> [--diff-csv FILE] [--t0-bad FILE]

Rows are matched by label and order of appearance (parallel-branch suffixes _2, _3, ... are removed
before matching, because paper_results/ has no suffixes in apfc files). Columns are matched by name.

In paper_results/, the reactive source columns of branches that inject at both ends (E_i_j_f, E_i_j_t)
are named in interleaved order while the data are stored as all from-ends followed by all to-ends.
Before comparing, the committed files are relabelled to the true identity of each column
(function relabel_case3, adapted from reports/T0_scripts/d_compare.py of the T0 analysis).
The appended rows of those branches (last 2*n3 rows) were written with the interleaved names, so
they keep their names. Uniqueness suffixes of the current code are removed before matching (base()).

Expected differences (acceptance 2 of brief T12a): parallel-branch labels in apfc; the column names
and order in case9 and case118 rpfc (removed by the relabelling); in case14 rpfc, the cells listed in
the T0 file D_bad_1.csv (they depend on the negative demand at bus 4) and the new source column D_4.
"""
import sys, os, re, argparse, collections
import numpy as np, pandas as pd

TOL = 1e-10
CASES = ["case1", "case9", "case14", "case30", "case33", "case118", "tap_loop"]


def tapped(case):
    """True if the exported input has a branch with TAP != 0 and != 1 or SHIFT != 0. The RPFC of such a case changed
    with rule Q3 version 2.0 (T16: line charging only, instead of the real-ratio pi-equivalent of the tap of version 1.0,
    which paper_results/ holds), so it is not compared with paper_results/."""
    import scipy.io as sio
    bm = np.asarray(sio.loadmat(os.path.join("benchmark", "inputs", f"{case}.mat"))["branch"], float)
    tap = np.abs(bm[:, 8]); tap[tap == 0] = 1
    return bool(((tap != 1) | (bm[:, 9] != 0)).any())
g = lambda x: f"{x:.17g}"


def relabel_case3(M):
    cols = list(M.columns[1:])
    blk = [i for i, x in enumerate(cols) if re.fullmatch(r"E_\d+_\d+_[ft]", x)]
    n3 = len(blk) // 2
    if n3 <= 1:
        return M
    pairs = [re.match(r"E_(\d+_\d+)_f", cols[blk[2 * i]]).group(1) for i in range(n3)]
    true = [f"E_{p}_f" for p in pairs] + [f"E_{p}_t" for p in pairs]
    ren = {cols[blk[j]]: true[j] for j in range(2 * n3)}
    body = M.iloc[:-2 * n3].rename(columns=ren)
    tail = M.iloc[-2 * n3:]
    return pd.concat([body, tail.reindex(columns=body.columns)], ignore_index=True)


def base(label):
    """Remove the uniqueness suffixes _2, _3, ... of the current code (parallel branches, generators at one
    bus, demand rows C2_, Tf_, Tt_, G_, Bs_ at one bus); paper_results/ has none."""
    m = re.fullmatch(r"(L_\{[ft]\.\d+-\d+\}|(?:G|D|S|Bs|Tf|Tt|C2)_\d+)_\d+", label)
    return m.group(1) if m else label


def keyed(labels):
    c = collections.Counter(); out = []
    for l in labels:
        b = base(l); c[b] += 1; out.append((b, c[b]))
    return out


def compare(new, ref):
    """Return (rows only in new, rows only in ref, cols only in new, cols only in ref, list of cell diffs, max abs diff, n cells)."""
    kn, kr = keyed(list(new.Line)), keyed(list(ref.Line))
    pn = {k: i for i, k in enumerate(kn)}; pr = {k: i for i, k in enumerate(kr)}
    # columns: by base name and order of appearance (generator columns G_b, G_b_2, ... of the current code)
    cn = [f"{b}#{k}" for b, k in keyed(list(new.columns[1:]))]
    cr = [f"{b}#{k}" for b, k in keyed(list(ref.columns[1:]))]
    new = new.set_axis(["Line"] + cn, axis=1); ref = ref.set_axis(["Line"] + cr, axis=1)
    common_c = [c for c in cn if c in cr]
    diffs, mx, n = [], 0.0, 0
    for k in kn:
        if k not in pr: continue
        a, b = new.iloc[pn[k]], ref.iloc[pr[k]]
        for c in common_c:
            d = abs(float(a[c]) - float(b[c])); n += 1; mx = max(mx, d)
            if d > TOL: diffs.append((k, c.split("#")[0], float(a[c]), float(b[c]), d))
    return ([k for k in kn if k not in pr], [k for k in kr if k not in pn],
            [c.split("#")[0] for c in cn if c not in cr], [c.split("#")[0] for c in cr if c not in cn], diffs, mx, n,
            [(k, new.iloc[pn[k]].Line, ref.iloc[pr[k]].Line) for k in kn if k in pr and new.iloc[pn[k]].Line != ref.iloc[pr[k]].Line])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("newdir"); ap.add_argument("--diff-csv"); ap.add_argument("--t0-bad")
    a = ap.parse_args()
    expected14 = set()
    if a.t0_bad and os.path.exists(a.t0_bad):
        b = pd.read_csv(a.t0_bad)
        b = b[(b.case == "case14") & (b.file == "rpfc")]
        for r, c in zip(b["matrix row"], b["matrix cols"]):
            expected14.add((re.sub(r"^C2[ft]_", "C2_", r), c))
    print("| case | file | rows compared | cells compared | max abs diff | cells above 1e-10 | of which expected | rows only in new | rows only in reference | columns only in new | columns only in reference | labels changed |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    allrows = []; unexpected_total = 0
    not_compared = []
    for c in CASES:
        for k in ["apfc", "rpfc"]:
            if k == "rpfc" and tapped(c):
                not_compared.append(c)
                print(f"| {c} | {k} | not compared: rule Q3 version 2.0 (T16) changes the RPFC of cases with tapped branches | | | | | | | | | |")
                continue
            new = pd.read_csv(os.path.join(a.newdir, f"{c}_{k}.csv"), float_precision="round_trip")
            ref = pd.read_csv(os.path.join("paper_results", f"{c}_{k}.csv"), float_precision="round_trip")
            if k == "rpfc": ref = relabel_case3(ref)
            on, orf, cn, cr, diffs, mx, n, relab = compare(new, ref)
            exp = 0
            for (key, col, x, y, d) in diffs:
                is_exp = c == "case14" and k == "rpfc" and (key[0], col) in expected14
                exp += is_exp
                allrows.append((c, k, key[0], key[1], col, x, y, d, is_exp))
            # a new column is expected only for case14 rpfc (D_4)
            col_unexp = [x for x in cn if not (c == "case14" and k == "rpfc" and x == "D_4")]
            unexpected_total += (len(diffs) - exp) + len(on) + len(orf) + len(col_unexp) + len(cr)
            rel = len(relab)
            print(f"| {c} | {k} | {len(new) - len(on)} | {n} | {g(mx)} | {len(diffs)} | {exp} | {', '.join(x[0] for x in on) or 'none'} | {', '.join(x[0] for x in orf) or 'none'} | {', '.join(cn) or 'none'} | {', '.join(cr) or 'none'} | {rel} |")
    if expected14:
        found = {(r[2], r[4]) for r in allrows if r[0] == "case14" and r[1] == "rpfc"}
        print(f"\nT0 cells of case14 rpfc expected to differ: {len(expected14)}; found among the differences: {len(expected14 & found)}; "
              f"differences outside that set: {len(found - expected14)}")
    if not_compared:
        print(f"\nRPFC not compared with paper_results/ (cases with tapped branches, rule Q3 version 2.0, T16): {', '.join(not_compared)}")
    print(f"\nUnexpected differences (cells outside the expected set, unmatched rows, unexpected columns): {unexpected_total}")
    if a.diff_csv:
        pd.DataFrame(allrows, columns=["case", "file", "row", "occurrence", "column", "new", "paper_results", "absdiff", "expected"]).to_csv(a.diff_csv, index=False, float_format="%.17g")
    return 0 if unexpected_total == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
