"""Input object of the Python port: a solved case, built from the DataFrames that ``load_results`` of the
iterative method (``external/iterative-method/power_flow.py``) returns for an exported ``results.mat``.

The exported inputs have already been preprocessed by ``my_runpf.m`` of the iterative method:
- the to-end flows PT and QT are the flows delivered at the to-end (MATPOWER sign reversed);
- branches whose four flows are all below 1e-4 in magnitude are deleted;
- single flows below 1e-4 in magnitude are set to zero;
- PS, QS (bus shunts) and QSF, QST (branch shunt flows) are computed.
The package itself does not import the iterative method.

The reactive port (RPFC) computes the branch shunt flows as rpfc.m does, from BR_R, BR_X, TAP, VM, VA and
baseMVA. ``load_results`` drops BR_R, BR_X, BR_B and baseMVA; ``read_branch_params`` reads them from the same
exported ``.mat`` file, and ``Case.from_frames(..., params=...)`` receives them. The exported QSF and QST are not
used: QSF of ``my_runpf.m`` has the from-end line charging as B/2 instead of B/(2 tap^2).
"""
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
import scipy.io as sio

BUS_COLUMNS = ("BUS_I", "PD", "QD", "VM", "VA")
GEN_COLUMNS = ("GEN_BUS", "PG", "QG")
BRANCH_COLUMNS = ("F_BUS", "T_BUS", "PF", "QF", "PT", "QT")
PARAM_COLUMNS = ("F_BUS", "T_BUS", "BR_R", "BR_X", "BR_B")
# columns of the MATPOWER branch matrix (0-based)
_F_BUS, _T_BUS, _BR_R, _BR_X, _BR_B = 0, 1, 2, 3, 4


@dataclass(frozen=True)
class BranchParams:
    """Branch impedance data and system base, read from the exported .mat file (rows as in the branch table)."""
    frame: pd.DataFrame      # columns F_BUS, T_BUS, BR_R, BR_X, BR_B
    base_mva: float


def read_branch_params(path):
    """Read BR_R, BR_X, BR_B (with F_BUS and T_BUS, for the alignment check) and baseMVA from an exported
    results.mat file (written by my_runpf.m). The rows are those of the saved branch matrix, which are also the
    rows of the branch DataFrame of load_results."""
    m = sio.loadmat(path)
    if "branch" not in m or "baseMVA" not in m:
        raise ValueError(f"{path}: no branch matrix or no baseMVA")
    br = np.asarray(m["branch"], dtype=float)
    frame = pd.DataFrame({"F_BUS": br[:, _F_BUS].astype(int), "T_BUS": br[:, _T_BUS].astype(int),
                          "BR_R": br[:, _BR_R], "BR_X": br[:, _BR_X], "BR_B": br[:, _BR_B]})
    return BranchParams(frame=frame, base_mva=float(np.asarray(m["baseMVA"]).ravel()[0]))


@dataclass(frozen=True)
class Case:
    """A solved power flow case: bus, generator and branch tables (copies of the input DataFrames), and,
    optionally, the branch parameters needed by the reactive port."""
    name: str
    bus: pd.DataFrame
    gen: pd.DataFrame
    branch: pd.DataFrame
    params: Optional[BranchParams] = None

    @classmethod
    def from_frames(cls, bus, gen, branch, name="", params=None):
        """Build a Case from the DataFrames ``bus``, ``gen`` and ``branch`` of ``load_results``.

        The frames are copied. ``branch.PT`` and ``branch.QT`` must be the flows delivered at the to-end
        (as exported by ``my_runpf.m``). The buses must be numbered 1..nb in the order of the bus table:
        ``apfc.m`` uses the bus numbers directly as matrix indices, so other numberings are not supported.

        ``params`` (optional, a ``BranchParams``, for example from ``read_branch_params``) is needed only by
        ``rpfc``. Its rows must align with ``branch`` (same number of rows, same F_BUS and T_BUS); otherwise a
        ValueError is raised.
        """
        for frame, cols, what in ((bus, BUS_COLUMNS, "bus"), (gen, GEN_COLUMNS, "gen"), (branch, BRANCH_COLUMNS, "branch")):
            missing = [c for c in cols if c not in frame.columns]
            if missing:
                raise ValueError(f"{what} table lacks the columns {missing}")
        nb = len(bus)
        if not np.array_equal(bus["BUS_I"].to_numpy(int), np.arange(1, nb + 1)):
            raise ValueError("buses must be numbered 1..nb in table order (as assumed by apfc.m)")
        if params is not None:
            check_alignment(branch, params)
            params = BranchParams(frame=params.frame.copy(), base_mva=float(params.base_mva))
        return cls(name=name, bus=bus.copy(), gen=gen.copy(), branch=branch.copy(), params=params)


def check_alignment(branch, params):
    """Raise ValueError unless the rows of ``params.frame`` align with ``branch`` (row count, F_BUS, T_BUS)."""
    missing = [c for c in PARAM_COLUMNS if c not in params.frame.columns]
    if missing:
        raise ValueError(f"branch parameters lack the columns {missing}")
    if len(params.frame) != len(branch):
        raise ValueError(f"branch parameters have {len(params.frame)} rows, the branch table {len(branch)}")
    for c in ("F_BUS", "T_BUS"):
        if not np.array_equal(params.frame[c].to_numpy(int), branch[c].to_numpy(int)):
            raise ValueError(f"branch parameters and branch table differ in {c}")
