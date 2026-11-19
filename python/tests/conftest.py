import os
import sys

sys.dont_write_bytecode = True
PYTHON_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
MATRIX_DIR = os.path.abspath(os.path.join(PYTHON_DIR, ".."))
if PYTHON_DIR not in sys.path:
    sys.path.insert(0, PYTHON_DIR)


import pytest  # noqa: E402


@pytest.fixture(autouse=True)
def crafted_shunt_flows(monkeypatch):
    """For crafted inputs only (branch table with the columns QSF_TEST and QST_TEST, made by test_rpfc.make_case):
    rule Q3 returns the crafted shunt flows, so that those tests keep checking the rules after Q3. Every other input
    uses the real rule Q3 (pfc_matrix.reactive.branch_shunt_flows)."""
    import pfc_matrix.reactive as reactive
    real = reactive.branch_shunt_flows

    def shunt_flows(case):
        br = case.branch
        if case.params is not None and "QSF_TEST" in br.columns and "QST_TEST" in br.columns:
            return br["QSF_TEST"].to_numpy(float).copy(), br["QST_TEST"].to_numpy(float).copy()
        return real(case)
    monkeypatch.setattr(reactive, "branch_shunt_flows", shunt_flows)
