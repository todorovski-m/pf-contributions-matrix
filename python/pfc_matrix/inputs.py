"""Input object of the Python port: a solved case, built from the DataFrames that ``load_results`` of the
iterative method (``external/iterative-method/power_flow.py``) returns for an exported ``results.mat``.

The exported inputs have already been preprocessed by ``my_runpf.m`` of the iterative method:
- the to-end flows PT and QT are the flows delivered at the to-end (MATPOWER sign reversed);
- branches whose four flows are all below 1e-4 in magnitude are deleted;
- single flows below 1e-4 in magnitude are set to zero;
- PS, QS (bus shunts) and QSF, QST (branch shunt flows) are computed.
The package itself does not import the iterative method.
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd

BUS_COLUMNS = ("BUS_I", "PD", "QD", "VM", "VA")
GEN_COLUMNS = ("GEN_BUS", "PG", "QG")
BRANCH_COLUMNS = ("F_BUS", "T_BUS", "PF", "QF", "PT", "QT")


@dataclass(frozen=True)
class Case:
    """A solved power flow case: bus, generator and branch tables (copies of the input DataFrames)."""
    name: str
    bus: pd.DataFrame
    gen: pd.DataFrame
    branch: pd.DataFrame

    @classmethod
    def from_frames(cls, bus, gen, branch, name=""):
        """Build a Case from the DataFrames ``bus``, ``gen`` and ``branch`` of ``load_results``.

        The frames are copied. ``branch.PT`` and ``branch.QT`` must be the flows delivered at the to-end
        (as exported by ``my_runpf.m``). The buses must be numbered 1..nb in the order of the bus table:
        ``apfc.m`` uses the bus numbers directly as matrix indices, so other numberings are not supported.
        """
        for frame, cols, what in ((bus, BUS_COLUMNS, "bus"), (gen, GEN_COLUMNS, "gen"), (branch, BRANCH_COLUMNS, "branch")):
            missing = [c for c in cols if c not in frame.columns]
            if missing:
                raise ValueError(f"{what} table lacks the columns {missing}")
        nb = len(bus)
        if not np.array_equal(bus["BUS_I"].to_numpy(int), np.arange(1, nb + 1)):
            raise ValueError("buses must be numbered 1..nb in table order (as assumed by apfc.m)")
        return cls(name=name, bus=bus.copy(), gen=gen.copy(), branch=branch.copy())
