"""Case list and loading for the timing harness.

Both methods take the same in-memory DataFrames of load_results (external/iterative-method/power_flow.py) of an
exported results.mat (my_runpf.m of the pinned submodule). Loading, the power flow and Octave are outside every
timed region. The old method is imported by path; os.system is replaced by a raising function, so an attempt to call
Octave fails instead of running. Before timing, the known issue of the old pipeline is corrected: in the branches with
tap != 1 and BR_B != 0, the exported QSF is replaced by the Qsf that rpfc.m computes (branch_shunt_flows of the port).
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
LARGE = ["case300", "case1354pegase", "case2869pegase", "case9241pegase", "case_ACTIVSg2000", "case_ACTIVSg10k"]


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


def qsf_rows(branch, params):
    """Rows of the branches with tap != 1 and BR_B != 0 (load_results has already set TAP = 0 to 1)."""
    import numpy as np
    tap = branch["TAP"].to_numpy(float)
    return np.where((tap != 1) & (params.frame["BR_B"].to_numpy(float) != 0))[0]


def correct_qsf(data):
    """A copy of the branch DataFrame in which QSF is replaced, in the branches with tap != 1 and BR_B != 0 only, by
    the Qsf that rpfc.m computes (pfc_matrix.reactive.branch_shunt_flows). Returns (branch, rows, old, new)."""
    from pfc_matrix import Case, read_branch_params
    from pfc_matrix.reactive import branch_shunt_flows
    params = read_branch_params(data["path"])
    rows = qsf_rows(data["branch"], params)
    branch = data["branch"].copy()
    if len(rows) == 0:
        return branch, rows, [], []
    case = Case.from_frames(data["bus"], data["gen"], data["branch"], data["case"], params=params)
    qsf, _ = branch_shunt_flows(case)
    old = branch["QSF"].to_numpy(float)[rows].tolist()
    branch.loc[branch.index[rows], "QSF"] = qsf[rows]
    return branch, rows, old, qsf[rows].tolist()
