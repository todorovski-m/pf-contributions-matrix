"""Checks of the canonical benchmark inputs (benchmark/inputs/*.mat).

Run from the root of the matrix repository. Each command prints Markdown tables.

    python3 benchmark/check_inputs.py repro <dir A> <dir B>
        Every array of every <case>.mat must be identical in A and B, except `et`.
    python3 benchmark/check_inputs.py compat [inputs dir]
        Load each file with load_results from the pinned submodule, with os.system replaced
        by a function that raises (so no Octave call can occur). Report shapes and warnings.
    python3 benchmark/check_inputs.py consistency [inputs dir]
        For each of the six shared cases, run runpf of the submodule (this calls Octave) in a
        fresh scratch copy under /tmp, in a separate process, and compare its DataFrames
        exactly with load_results(<inputs>/<case>.mat).
    python3 benchmark/check_inputs.py oldmethod [inputs dir]
        Run makeAPFC and makeRPFC of the submodule on each exported input, without Octave,
        with itermax=200, and report DONE and the number of iterations.

The submodule code is imported by path from external/iterative-method; nothing is copied.
"""
import sys
sys.dont_write_bytecode = True  # never write __pycache__ into the submodule
import os, re, pickle, shutil, subprocess, tempfile, warnings, logging, importlib.util
import numpy as np
import pandas as pd
import scipy.io as sio

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
SUB = os.path.join(ROOT, "external", "iterative-method")
INPUTS = os.path.join(ROOT, "benchmark", "inputs")
CASES = ["case1", "case9", "case14", "case24_ieee_rts", "case30", "case33", "case118", "tap_loop"]
SHARED = ["case1", "case9", "case14", "case30", "case33", "case118"]
IGNORED = {"et"}
g = lambda x: f"{x:.17g}"


def _raise(*a, **k):
    raise RuntimeError("os.system was called (an Octave call was attempted)")


def import_submodule(block_os_system=True):
    """Import power_flow and power_components from the submodule by path."""
    mods = {}
    for name in ["power_flow", "power_components"]:
        if name in sys.modules:
            mods[name] = sys.modules[name]
            continue
        spec = importlib.util.spec_from_file_location(name, os.path.join(SUB, f"{name}.py"))
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod  # power_components does "from power_flow import runpf"
        spec.loader.exec_module(mod)
        mods[name] = mod
    if block_os_system:
        os.system = _raise
        mods["power_flow"].os.system = _raise
    return mods["power_flow"], mods["power_components"]


def load_mat(path):
    return {k: v for k, v in sio.loadmat(path).items() if not k.startswith("__")}


def cmd_repro(a, b):
    print("| case | variables compared (excluding et) | all identical | differing variables | et in A | et in B |\n|---|---|---|---|---|---|")
    ok = True
    for c in CASES:
        x, y = load_mat(os.path.join(a, f"{c}.mat")), load_mat(os.path.join(b, f"{c}.mat"))
        keys = sorted((set(x) | set(y)) - IGNORED)
        bad = []
        for k in keys:
            if k not in x or k not in y:
                bad.append(f"{k} (missing)"); continue
            u, v = x[k], y[k]
            same = u.shape == v.shape and u.dtype == v.dtype and (
                np.array_equal(u, v, equal_nan=True) if u.dtype.kind in "biufc" else repr(u.tolist()) == repr(v.tolist()))
            if not same: bad.append(k)
        ok &= not bad
        print(f"| {c} | {len(keys)} | {not bad} | {', '.join(bad) or 'none'} | {g(float(x['et'].ravel()[0]))} | {g(float(y['et'].ravel()[0]))} |")
    print(f"\nAll cases identical except et: {ok}")
    return 0 if ok else 1


def cmd_compat(inputs):
    pf, _ = import_submodule()
    print("| case | bus shape | gen shape | branch shape | xy | warnings |\n|---|---|---|---|---|---|")
    for c in CASES:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            bus, gen, branch, xy = pf.load_results(os.path.join(inputs, f"{c}.mat"))
        msgs = sorted({str(m.message).splitlines()[0] for m in w})
        print(f"| {c} | {bus.shape} | {gen.shape} | {branch.shape} | {'None' if xy is None else xy.shape} | {msgs or 'none'} |")
    print(f"\nos.system in power_flow is the raising stand-in: {pf.os.system is _raise}")
    return 0


RUNPF_CHILD = r"""
import sys, pickle, importlib.util, os
spec = importlib.util.spec_from_file_location("power_flow", sys.argv[1])
pf = importlib.util.module_from_spec(spec); spec.loader.exec_module(pf)
out = pf.runpf(sys.argv[2])
pickle.dump(out, open(sys.argv[3], "wb"))
"""


def cmd_consistency(inputs):
    pf, _ = import_submodule()
    scratch = tempfile.mkdtemp(prefix="check_inputs_consistency.", dir="/tmp")
    shutil.copytree(SUB, scratch, dirs_exist_ok=True)
    print(f"Scratch copy of the submodule: {scratch}\n")
    print("| case | runpf exit status | bus equal | gen equal | branch equal | xy equal |\n|---|---|---|---|---|---|")
    ok = True
    for c in SHARED:
        pkl = os.path.join(scratch, f"{c}.pkl")
        r = subprocess.run([sys.executable, "-B", "-c", RUNPF_CHILD, os.path.join(SUB, "power_flow.py"), c, pkl],
                           cwd=scratch, capture_output=True, text=True)
        if r.returncode != 0:
            ok = False; print(f"| {c} | {r.returncode} | | | | |"); continue
        ref = pickle.load(open(pkl, "rb"))
        got = pf.load_results(os.path.join(inputs, f"{c}.mat"))
        eq = []
        for a, b in zip(got[:3], ref[:3]):
            try:
                pd.testing.assert_frame_equal(a, b, check_exact=True); eq.append(True)
            except AssertionError:
                eq.append(False)
        xy_eq = (got[3] is None and ref[3] is None) or (got[3] is not None and ref[3] is not None and np.array_equal(got[3], ref[3]))
        ok &= all(eq) and xy_eq
        print(f"| {c} | 0 | {eq[0]} | {eq[1]} | {eq[2]} | {xy_eq} |")
    t1 = "/tmp/iter_after/results/results.mat"
    if os.path.exists(t1):
        x, y = load_mat(t1), load_mat(os.path.join(inputs, "case118.mat"))
        keys = sorted((set(x) | set(y)) - IGNORED)
        same = all(k in x and k in y and x[k].shape == y[k].shape and
                   (np.array_equal(x[k], y[k]) if x[k].dtype.kind in "biufc" else repr(x[k].tolist()) == repr(y[k].tolist()))
                   for k in keys)
        print(f"\nAdditional: {t1} (last case of the T1 run, case118) against the exported case118.mat, "
              f"all {len(keys)} variables except et identical: {same}")
    print(f"\nAll shared cases equal: {ok}")
    return 0 if ok else 1


class Collect(logging.Handler):
    def __init__(self):
        super().__init__(level=logging.INFO); self.lines = []
    def emit(self, record):
        self.lines.append(record.getMessage())


def cmd_oldmethod(inputs):
    pf, pc = import_submodule()
    # same levels as in T0: root logger NOTSET, one handler at INFO; no file handler
    root = logging.getLogger(); root.setLevel(logging.NOTSET)
    for h in list(root.handlers): root.removeHandler(h)
    print("| case | APFC DONE | APFC iterations | APFC rows not done | RPFC DONE | RPFC iterations | RPFC rows not done | exception |\n|---|---|---|---|---|---|---|---|")
    notdone = []
    for c in CASES:
        h = Collect(); root.addHandler(h)
        res = {"APFC": ("", "", ""), "RPFC": ("", "", "")}; exc = "none"
        try:
            bus, gen, branch, xy = pf.load_results(os.path.join(inputs, f"{c}.mat"))
            bus_pos = None if xy is None else {int(i): [x, y] for i, x, y in xy}
            _, A = pc.makeAPFC(bus, gen, branch, bus_pos, 200, False)
            _, R = pc.makeRPFC(bus, gen, branch, bus_pos, 200, False)
            for name, T in [("APFC", A), ("RPFC", R)]:
                m = [re.search(rf"{name}: nr = \d+, iters = (\d+), DONE = (\w+)", l) for l in h.lines]
                m = [x for x in m if x][-1]
                nd = ~T["DONE"].astype(bool)
                res[name] = (m.group(2), m.group(1), str(int(nd.sum())))
                if nd.any():
                    notdone.append((c, name, list(T.index[nd])))
        except Exception as e:
            exc = f"{type(e).__name__}: {e}"
        root.removeHandler(h)
        print(f"| {c} | {res['APFC'][0]} | {res['APFC'][1]} | {res['APFC'][2]} | {res['RPFC'][0]} | {res['RPFC'][1]} | {res['RPFC'][2]} | {exc} |")
    for c, name, rows in notdone:
        print(f"\nRows not done in {name} of {c}: {', '.join(rows)}")
    return 0


if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "repro":
        sys.exit(cmd_repro(*args))
    inputs = args[0] if args else INPUTS
    sys.exit({"compat": cmd_compat, "consistency": cmd_consistency, "oldmethod": cmd_oldmethod}[cmd](inputs))
