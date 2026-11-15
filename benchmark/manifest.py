"""Helpers for benchmark/export_inputs.sh: check one exported file, and write MANIFEST.json.

Usage:
    python3 benchmark/manifest.py check <results.mat>
        Exit with status 1 unless the file loads and its variable `success` equals 1.
    python3 benchmark/manifest.py manifest <inputs dir> <submodule commit> <matrix commit> <matrix dirty> <case>=<case file> ...
        Write <inputs dir>/MANIFEST.json.

The content hash (`content_sha256_excluding_et`) is computed over every variable of the
.mat file except `et` (MATPOWER's elapsed time), in sorted order of the variable names,
from the name, dtype, shape and bytes of each array. It is reproducible between runs.
The file hash (`mat_sha256`) is not: the MAT-file header written by Octave contains the
time of writing, and `et` differs between runs.
"""
import sys, os, json, hashlib, datetime
import numpy as np
import scipy.io as sio

IGNORED = {"et"}


def load(path):
    m = sio.loadmat(path)
    return {k: v for k, v in m.items() if not k.startswith("__")}


def content_hash(path):
    h = hashlib.sha256()
    m = load(path)
    for k in sorted(m):
        if k in IGNORED:
            continue
        v = m[k]
        h.update(k.encode())
        h.update(str(v.dtype).encode())
        h.update(str(v.shape).encode())
        if v.dtype.kind in "biufc":
            h.update(np.ascontiguousarray(v).tobytes())
        else:
            h.update(repr(v.tolist()).encode())
    return h.hexdigest()


def sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def check(path):
    try:
        m = load(path)
    except Exception as e:
        print(f"ERROR: cannot load {path}: {e}", file=sys.stderr)
        return 1
    if "success" not in m:
        print(f"ERROR: {path} has no variable 'success'", file=sys.stderr)
        return 1
    s = float(np.asarray(m["success"]).ravel()[0])
    if s != 1:
        print(f"ERROR: {path}: success = {s!r}, expected 1 (Newton iteration did not converge)", file=sys.stderr)
        return 1
    return 0


def manifest(outdir, sub_commit, matrix_commit, matrix_dirty, pairs):
    cases = {}
    for p in pairs:
        case, case_file = p.split("=", 1)
        mat = os.path.join(outdir, f"{case}.mat")
        m = load(mat)
        cases[case] = {
            "case_file": os.path.relpath(case_file),
            "case_file_sha256": sha256(case_file),
            "mat_file": f"{case}.mat",
            "mat_sha256": sha256(mat),
            "content_sha256_excluding_et": content_hash(mat),
            "n_bus": int(m["bus"].shape[0]),
            "n_branch": int(m["branch"].shape[0]),
            "n_gen": int(m["gen"].shape[0]),
        }
    doc = {
        "description": "Canonical benchmark inputs: results/results.mat written by my_runpf.m of the iterative method, one file per case.",
        "iterative_method_submodule_commit": sub_commit,
        "matrix_commit": matrix_commit,
        "matrix_working_tree_dirty": matrix_dirty == "1",
        "date_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "ignored_variables": sorted(IGNORED),
        "note": "n_branch counts the branches in the exported file, after the deletion done by my_runpf.m. "
                "mat_sha256 changes on every export (Octave writes the time into the MAT-file header, and 'et' is an elapsed time); "
                "content_sha256_excluding_et is reproducible.",
        "cases": cases,
    }
    with open(os.path.join(outdir, "MANIFEST.json"), "w") as f:
        json.dump(doc, f, indent=2)
        f.write("\n")
    return 0


if __name__ == "__main__":
    if sys.argv[1] == "check":
        sys.exit(check(sys.argv[2]))
    elif sys.argv[1] == "manifest":
        sys.exit(manifest(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6:]))
    sys.exit(f"unknown command {sys.argv[1]}")
