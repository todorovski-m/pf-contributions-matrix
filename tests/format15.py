"""Rewrite the numbers of CSV files written by the Octave stand-in writetable (%.17g) with 15
significant digits, the precision of Matlab's writetable for double values (C format %.15g, for
example 8.84367205134871e-05; negative zero is written as -0). The header and the row labels are
copied unchanged. Checked against paper_results/ (written by Matlab from the v1.0 code): after this
conversion, the Octave outputs of the v1.0 code differ from them only in the last digit or in
values near zero (below 1e-13), never in the format.

Usage: python3 tests/format15.py <input dir> <output dir> [file ...]   (default: every *.csv)"""
import sys, os, glob


def fmt(tok):
    try:
        v = float(tok)
    except ValueError:
        return tok
    return "%.15g" % v      # also "-0" for negative zero, as Matlab writes it


def convert(src, dst):
    with open(src) as f:
        lines = f.read().splitlines()
    out = [lines[0]] + [",".join([l.split(",", 1)[0]] + [fmt(t) for t in l.split(",")[1:]]) for l in lines[1:]]
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    src_dir, dst_dir = sys.argv[1], sys.argv[2]
    files = sys.argv[3:] or sorted(os.path.basename(p) for p in glob.glob(os.path.join(src_dir, "*.csv")))
    os.makedirs(dst_dir, exist_ok=True)
    for f in files:
        convert(os.path.join(src_dir, f), os.path.join(dst_dir, f))
    print(f"converted {len(files)} files to {dst_dir}")
