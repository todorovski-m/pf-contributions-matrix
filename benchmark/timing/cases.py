"""Case list and loading for the timing harness.

Both methods take the same in-memory DataFrames of load_results (external/iterative-method/power_flow.py) of an
exported results.mat (my_runpf.m of the pinned submodule). Loading, the power flow and Octave are outside every
timed region. The old method is imported by path; os.system is replaced by a raising function, so an attempt to call
Octave fails instead of running. Before timing, the old method's input is aligned with rule Q3 version 2.0 (T16): in the
branches with a tap or a phase shift, the exported QSF and QST (the real-ratio pi-equivalent of my_runpf.m) are replaced
by the line charging that rpfc.m and the port use (correct_qsf).
"""
import importlib.util
import logging
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MATRIX = os.path.abspath(os.path.join(HERE, "..", ".."))
SUB = os.path.join(MATRIX, "external", "iterative-method")
sys.path.insert(0, os.path.join(MATRIX, "python"))

SMALL = ["case1", "case9", "case14", "case24_ieee_rts", "case30", "case33", "case118", "tap_loop",
         "syn9_neg_pd", "syn9_neg_pg", "syn9_neg_qd", "syn9_neg_qg", "syn9_case3_tap", "syn9_carrier_tap"]
# The large (case, kinds) combinations, in order of size. The nine cases of T14c (B1) are in by usability. Left out:
#   case3120sp apfc (pfc:apfc:negativeLoss: six branches with negative series resistance in the data);
#   case_ACTIVSg10k rpfc and case9241pegase rpfc (estimated dense reactive tables above the memory limit);
#   case9241pegase apfc (pfc:apfc:spectralRadius: a circulating 2-cycle of gain 1 between buses 880 and 6671).
# Of the 15 cases screened in T14f, only combinations that are usable, predicted "can finish" by old_stall_predictor.py,
# without a division-by-zero risk and within the memory limit are added: case_ACTIVSg200 (both kinds). The reasons for
# every other screened combination are in the README (table of the screened cases).
LARGE_COMBOS = [("case_ACTIVSg200", ("apfc", "rpfc")), ("case300", ("apfc", "rpfc")), ("case_ACTIVSg500", ("apfc", "rpfc")), ("case1354pegase", ("apfc", "rpfc")),
                ("case_ACTIVSg2000", ("apfc", "rpfc")), ("case2383wp", ("apfc", "rpfc")), ("case2869pegase", ("apfc", "rpfc")),
                ("case3120sp", ("rpfc",)), ("case_ACTIVSg10k", ("apfc",))]
LARGE = [c for c, _ in LARGE_COMBOS]
COMBOS = [(c, ("apfc", "rpfc")) for c in SMALL] + LARGE_COMBOS


def parse_combos(text):
    """'case9:apfc,case9:rpfc,case3120sp:rpfc' -> [(case, kinds)] in the order of first appearance."""
    out = {}
    for item in text.split(","):
        c, k = item.split(":")
        assert k in ("apfc", "rpfc"), item
        out.setdefault(c, []).append(k)
    return [(c, tuple(k)) for c, k in out.items()]


def combos_text(combos):
    return ",".join(f"{c}:{k}" for c, ks in combos for k in ks)


def input_dirs(syn_inputs=None, large_inputs=None):
    """Folders searched for <case>.mat, in this order: benchmark/inputs, the synthetic exports (argument or
    PFC_SYN_INPUTS), python/tests/cases (syn9_carrier_tap), the large inputs (argument or benchmark/inputs_large)."""
    d = [os.path.join(MATRIX, "benchmark", "inputs")]
    syn = syn_inputs or os.environ.get("PFC_SYN_INPUTS")
    if syn:
        d.append(syn)
    d.append(os.path.join(MATRIX, "python", "tests", "cases"))
    d.append(large_inputs or os.path.join(MATRIX, "benchmark", "inputs_large"))
    return d


def find_input(case, dirs):
    for d in dirs:
        p = os.path.join(d, f"{case}.mat")
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"{case}.mat not found in {dirs}")


_OLD = {}


def _raise_system(*a, **k):
    raise RuntimeError("os.system was called (an Octave call was attempted)")


def import_old():
    """(power_flow, power_components, compare_old) of the pinned submodule and tests/compare_old.py, imported by path
    without writing bytecode; os.system raises. The root logger gets a NullHandler and level WARNING, so that the
    logging calls of the old method create no records and print nothing."""
    if _OLD:
        return _OLD["pf"], _OLD["pc"], _OLD["co"]
    sys.dont_write_bytecode = True
    os.system = _raise_system
    root = logging.getLogger()
    for h in list(root.handlers):
        root.removeHandler(h)
    root.addHandler(logging.NullHandler())
    root.setLevel(logging.WARNING)
    mods = {}
    for name, path in (("power_flow", os.path.join(SUB, "power_flow.py")),
                       ("power_components", os.path.join(SUB, "power_components.py")),
                       ("compare_old", os.path.join(MATRIX, "tests", "compare_old.py"))):
        spec = importlib.util.spec_from_file_location(name, path)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
        mods[name] = mod
    mods["power_flow"].os.system = _raise_system
    for h in list(root.handlers):          # compare_old or the old code may have configured logging
        if not isinstance(h, logging.NullHandler):
            root.removeHandler(h)
    root.setLevel(logging.WARNING)
    _OLD.update(pf=mods["power_flow"], pc=mods["power_components"], co=mods["compare_old"])
    return _OLD["pf"], _OLD["pc"], _OLD["co"]


def load(case, dirs):
    """The DataFrames of load_results and the path of the exported file."""
    pf, _, _ = import_old()
    path = find_input(case, dirs)
    bus, gen, branch, _ = pf.load_results(path)
    return {"case": case, "path": path, "bus": bus, "gen": gen, "branch": branch}


def qsf_rows(data):
    """Rows of the branches with |TAP| != 1 (TAP = 0 means 1) or SHIFT != 0, read from the branch matrix of the exported
    .mat file: the branches whose exported QSF and QST (my_runpf.m: the real-ratio pi-equivalent of the tap, with B/2 at
    the from-end) differ from the line charging of rule Q3 version 2.0."""
    import numpy as np
    import scipy.io as sio
    bm = np.asarray(sio.loadmat(data["path"])["branch"], float)
    tap = np.abs(bm[:, 8]); tap[tap == 0] = 1.0
    return np.where((tap != 1) | (bm[:, 9] != 0))[0]


def correct_qsf(data):
    """A copy of the branch DataFrame in which QSF and QST of the tapped or phase-shifting branches are the line charging
    of rule Q3 version 2.0 (T16), QSF = (B/2) VM_f^2/|TAP|^2 baseMVA and QST = (B/2) VM_t^2 baseMVA, which rpfc.m and the
    port use (tests/compare_old.formula_shunts). Every other branch keeps its exported values (they equal the line
    charging there). The old method is run unchanged on this input. Returns (branch, rows, old, new) with old and new
    the [QSF, QST] pairs of the replaced rows."""
    co = import_old()[2]
    branch, rows = co.formula_shunts(data["branch"], data["bus"], data["path"])
    old = data["branch"][["QSF", "QST"]].to_numpy(float)[rows].tolist()
    new = branch[["QSF", "QST"]].to_numpy(float)[rows].tolist()
    return branch, rows, old, new
