"""CSV values are read exactly (pandas float_precision="round_trip"; the default parser can be off by
about 1e-13 relative for 17-digit values).

Compare two output folders of apfc.m/rpfc.m (for example the current code against the T12a end state).

Usage (from the root of the matrix repository):
    python3 tests/compare_runs.py <dir A (new)> <dir B (reference)> [--tol 1e-12]

Every <case>_apfc.csv and <case>_rpfc.csv present in both folders is compared by position: same
number of rows and columns, values within the tolerance. Labels (row labels and column names) are
compared as strings; a label that differs only by a uniqueness suffix (_2, _3, ...) added to the
reference label is reported as "suffix only", any other label change as "other".
"""
import sys, os, re, argparse
import numpy as np, pandas as pd

g = lambda x: f"{x:.17g}"


def suffix_only(new, ref):
    return new != ref and re.fullmatch(re.escape(ref) + r"_\d+", new) is not None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("new"); ap.add_argument("ref"); ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--skip-rpfc", default="", help="comma-separated cases whose rpfc file is not compared (for example "
                    "the cases with tapped or phase-shifting branches, whose RPFC rule Q3 version 2.0 changed, T16)")
    a = ap.parse_args()
    skip = {f"{c}_rpfc.csv" for c in a.skip_rpfc.split(",") if c}
    files = sorted(f for f in os.listdir(a.ref) if f.endswith(("_apfc.csv", "_rpfc.csv")) and os.path.exists(os.path.join(a.new, f))
                   and f not in skip)
    only_ref = sorted(f for f in os.listdir(a.ref) if f.endswith(".csv") and not os.path.exists(os.path.join(a.new, f)))
    print(f"| file | same shape | max abs diff | values within {a.tol:g} | row labels changed (suffix only / other) | column names changed (suffix only / other) |")
    print("|---|---|---|---|---|---|")
    fails, changes = 0, []
    for f in files:
        N = pd.read_csv(os.path.join(a.new, f), float_precision="round_trip"); R = pd.read_csv(os.path.join(a.ref, f), float_precision="round_trip")
        if N.shape != R.shape:
            print(f"| {f} | False ({N.shape} vs {R.shape}) | | False | | |"); fails += 1; continue
        d = np.abs(N.iloc[:, 1:].to_numpy(float) - R.iloc[:, 1:].to_numpy(float))
        mx = float(d.max()) if d.size else 0.0
        rl = [(x, y) for x, y in zip(N.Line, R.Line) if x != y]
        cl = [(x, y) for x, y in zip(N.columns[1:], R.columns[1:]) if x != y]
        rs, ro = sum(suffix_only(x, y) for x, y in rl), sum(not suffix_only(x, y) for x, y in rl)
        cs, co = sum(suffix_only(x, y) for x, y in cl), sum(not suffix_only(x, y) for x, y in cl)
        ok = mx <= a.tol and ro == 0 and co == 0
        fails += not ok
        changes += [(f, "row", y, x) for x, y in rl] + [(f, "column", y, x) for x, y in cl]
        print(f"| {f} | True | {g(mx)} | {mx <= a.tol} | {rs} / {ro} | {cs} / {co} |")
    if only_ref:
        print(f"\nFiles only in the reference folder: {only_ref}")
    if changes:
        print("\nLabel changes (file, kind, reference label -> new label):\n")
        for c in changes:
            print(f"- {c[0]}, {c[1]}: `{c[2]}` -> `{c[3]}`")
    if skip:
        print(f"\nNot compared (rpfc of the cases with tapped or phase-shifting branches; rule Q3 version 2.0, T16): {sorted(skip)}")
    print(f"\nFiles compared: {len(files)}; failed (values above tolerance, shape change, or non-suffix label change): {fails}")
    return 0 if fails == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
