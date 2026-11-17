"""Reactive power flow contributions (RPFC): Python port of rpfc.m (branch revision of the matrix repository).

The computation is split into four phases, which can be called (and later timed) on their own:
    prepare_reactive(case, tol)           (a) shunt flows, conversions, deletion, orientation, Cases 1 to 3,
                                              demands, sources; index and flow arrays of the expanded network
    assemble_reactive(prep, tol, tol_div) (b) incidence matrices, attenuation, sparse M, N, A, B, the self-terms
                                              and the self-term right-hand side
    solve_reactive(system)                (c) one sparse LU factorisation of A, used for both solves
    build_reactive_table(prep, system, X, X_self)   (d) output table (pandas DataFrame) and labels
``rpfc(case, tol, tol_div)`` calls them in this order. No phase prints, times or writes files. The validity checks
of rpfc.m (module checks: inflow, series flow, cycle, and the guard pfc:rpfc:negativeSelfTerm) run inside phase (b);
the self-terms are also returned in ``info["self_toend"]``.

The rules are those of rpfc.m; the numbered rules Q1 to Q26 of the T13b report are referenced as [Qn].
"""
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .labels import disambiguate_parallel, unique_labels
from .options import resolve_options, PfcError
from .inputs import Case
from . import checks

PERMC_SPEC = "COLAMD"   # column ordering of scipy.sparse.linalg.splu (SuperLU)
REACTIVE_BUS_COLUMNS = ("QS",)
REACTIVE_BRANCH_COLUMNS = ("TAP",)


@dataclass
class ReactivePrep:
    """Result of phase (a): the expanded network (carriers followed by demand branches) and the sources."""
    nb: int
    nl: int                 # number of carriers (Case 1 after orientation)
    nd: int                 # number of demand branches
    ns: int                 # number of source columns (without the self-term columns)
    n_pre: int              # number of source columns before the Case-3 columns
    f: np.ndarray           # 0-based from-bus index of every branch of the expanded network (length ne)
    t: np.ndarray           # 0-based to-bus index (the sink for demand branches)
    Qf: np.ndarray
    Qt: np.ndarray
    Qsf: np.ndarray
    Qst: np.ndarray
    src: np.ndarray         # 0-based bus index of every source column
    Qsrc: np.ndarray
    src_names: list         # source column names before unique_labels
    d_labels: list          # demand row labels before unique_labels
    case3_fbus: np.ndarray
    case3_tbus: np.ndarray
    case3_QF: np.ndarray    # injection into the from-bus (positive)
    case3_QT: np.ndarray    # injection into the to-bus (positive)
    case3_hasf: np.ndarray  # the member has a from-end source (False for a (0,+) member, T15b)
    case3_hast: np.ndarray  # the member has a to-end source (False for a (-,0) member, T15b)
    removed: dict
    flipped_branches: int
    converted: dict

    @property
    def ne(self):
        return self.nl + self.nd

    @property
    def n_case3(self):
        return len(self.case3_fbus)


@dataclass
class ReactiveSystem:
    """Result of phase (b)."""
    M: sp.csc_matrix
    N: sp.csc_matrix
    A: sp.csc_matrix
    B: sp.csc_matrix
    B_self: sp.csc_matrix   # right-hand side of the self-term solve (columns: carriers with self_toend > tol_div)
    i_self: np.ndarray      # 0-based carrier indices of the columns of B_self
    aQ: np.ndarray
    Qs: np.ndarray
    Qe: np.ndarray
    self_toend: np.ndarray  # self-term of every carrier [Q18]
    negative_series_flow: list
    acyclic: bool = True    # [Q15] checks
    spectral_radius: float = 0.0
    warnings: list = field(default_factory=list)


def _require(frame, cols, what):
    missing = [c for c in cols if c not in frame.columns]
    if missing:
        raise PfcError("pfc:input:columns", f"{what} table lacks the columns {missing} (needed by rpfc)")


def branch_shunt_flows(case):
    """[Q3] Shunt flows at the from-end and the to-end of every branch of the input, computed as in rpfc.m from
    BR_R, BR_X, TAP, VM, VA and baseMVA: I = Yb/tap (Vf - Vt), Qs = Im(Vf conj(I)) baseMVA,
    Qe = Im(Vt conj(I)) baseMVA, Qsf = Qs - QF, Qst = QT - Qe (QT delivered)."""
    if case.params is None:
        raise PfcError("pfc:input:branchParams", "rpfc needs the branch parameters: build the Case with params=read_branch_params(path)")
    br, bus, p = case.branch, case.bus, case.params
    tap = br["TAP"].to_numpy(float).copy()
    tap[tap == 0] = 1.0
    V = bus["VM"].to_numpy(float) * np.exp(1j * bus["VA"].to_numpy(float) / 180 * np.pi)
    fb = br["F_BUS"].to_numpy(int) - 1
    tb = br["T_BUS"].to_numpy(int) - 1
    Yb = 1.0 / (p.frame["BR_R"].to_numpy(float) + 1j * p.frame["BR_X"].to_numpy(float))
    I = Yb / tap * (V[fb] - V[tb])
    Qs = np.imag(V[fb] * np.conj(I)) * p.base_mva
    Qe = np.imag(V[tb] * np.conj(I)) * p.base_mva
    return Qs - br["QF"].to_numpy(float), br["QT"].to_numpy(float) - Qe


def prepare_reactive(case, tol=1e-6):
    """Phase (a): rules [Q3] to [Q12] on the exported input."""
    _require(case.bus, REACTIVE_BUS_COLUMNS, "bus")
    _require(case.branch, REACTIVE_BRANCH_COLUMNS, "branch")
    bus, gen, br = case.bus, case.gen, case.branch
    nb = len(bus)
    F = br["F_BUS"].to_numpy(int); T = br["T_BUS"].to_numpy(int)
    QF = br["QF"].to_numpy(float); QT = br["QT"].to_numpy(float)   # [Q2] QT is the delivered flow
    Qsf, Qst = branch_shunt_flows(case)                               # [Q3]
    # Terminal classification (T15b), before the flip, sign by tol: zt = 1 (+,0), 2 (0,-) are Case-2 members and
    # zt = 3 (-,0), 4 (0,+) Case-3 members with the zero terminal omitted; zt = 0 keeps the exact-sign tests below
    s_f = (QF > tol).astype(int) - (QF < -tol).astype(int)
    s_t = (QT > tol).astype(int) - (QT < -tol).astype(int)
    zt = np.zeros(len(QF), dtype=int)
    zt[(s_f > 0) & (s_t == 0)] = 1
    zt[(s_f == 0) & (s_t < 0)] = 2
    zt[(s_f < 0) & (s_t == 0)] = 3
    zt[(s_f == 0) & (s_t > 0)] = 4
    # [Q4] Case-2 mask on all branches; [Q5] inductive branch-end shunts become demands (before the deletion)
    is_case2 = ((QF > 0) & (QT < 0) & (zt == 0)) | (zt == 1) | (zt == 2)
    i_neg_qsf = (Qsf < -tol) & ~is_case2
    extra_bus_qsf, extra_Q_qsf = F[i_neg_qsf], -Qsf[i_neg_qsf]
    Qsf = np.where(i_neg_qsf, 0.0, Qsf)
    i_neg_qst = (Qst < -tol) & ~is_case2
    extra_bus_qst, extra_Q_qst = T[i_neg_qst], -Qst[i_neg_qst]
    Qst = np.where(i_neg_qst, 0.0, Qst)
    # [Q6] deletion of branches and generators
    i_zero = (np.abs(QF) < tol) & (np.abs(QT) < tol)
    keep = ~i_zero
    F, T, QF, QT, Qsf, Qst, zt = F[keep], T[keep], QF[keep], QT[keep], Qsf[keep], Qst[keep], zt[keep]
    QG_all = gen["QG"].to_numpy(float)
    i_zero_gen = np.abs(QG_all) < tol
    gbus, QG = gen["GEN_BUS"].to_numpy(int)[~i_zero_gen], QG_all[~i_zero_gen]
    # [Q7] flip of the branches with both flows negative
    flip = (QF < 0) & (QT < 0) & (zt == 0)
    F, T = np.where(flip, T, F), np.where(flip, F, T)
    QF, QT = np.where(flip, -QT, QF), np.where(flip, -QF, QT)
    Qsf, Qst = np.where(flip, Qst, Qsf), np.where(flip, Qsf, Qst)
    # [Q8] Case 3: injects at both ends -> two sources (one for a member with a zero terminal)
    c3 = ((QF < 0) & (QT > 0) & (zt == 0)) | (zt == 3) | (zt == 4)
    case3_fbus, case3_tbus, case3_QF, case3_QT = F[c3], T[c3], -QF[c3], QT[c3]
    case3_hasf, case3_hast = zt[c3] != 4, zt[c3] != 3
    n_case3_zt = int((zt[c3] != 0).sum())
    F, T, QF, QT, Qsf, Qst, zt = F[~c3], T[~c3], QF[~c3], QT[~c3], Qsf[~c3], Qst[~c3], zt[~c3]
    # [Q9] Case 2: consumes at both ends -> two demands (all from-ends, then all to-ends; one for a member with a
    # zero terminal)
    c2 = ((QF > 0) & (QT < 0) & (zt == 0)) | (zt == 1) | (zt == 2)
    c2f, c2t = c2 & (zt != 2), c2 & (zt != 1)
    c2_bus = np.r_[F[c2f], T[c2t]]
    c2_Q = np.r_[QF[c2f], -QT[c2t]]
    n_case2_zt, n_case2_both = int((c2 & (zt != 0)).sum()), int((c2 & (zt == 0)).sum())
    F, T, QF, QT, Qsf, Qst = F[~c2], T[~c2], QF[~c2], QT[~c2], Qsf[~c2], Qst[~c2]
    nl = len(F)
    # [Q10] demands
    bus_i = bus["BUS_I"].to_numpy(int); QD = bus["QD"].to_numpy(float)
    Qsh = bus["QS"].to_numpy(float)      # VM^2 * BS (my_runpf.m computes BS .* VM.^2, the same product)
    i_qd = QD > tol
    i_qg_neg = QG < -tol
    i_qs_neg = Qsh < -tol
    d_bus = np.r_[bus_i[i_qd], gbus[i_qg_neg], bus_i[i_qs_neg], c2_bus, extra_bus_qsf, extra_bus_qst].astype(int)
    Qd = np.r_[QD[i_qd], -QG[i_qg_neg], -Qsh[i_qs_neg], c2_Q, extra_Q_qsf, extra_Q_qst]
    d_labels = ([f"D_{b}" for b in bus_i[i_qd]] + [f"G_{b}" for b in gbus[i_qg_neg]] +
                [f"Bs_{b}" for b in bus_i[i_qs_neg]] + [f"C2_{b}" for b in c2_bus] +
                [f"Tf_{b}" for b in extra_bus_qsf] + [f"Tt_{b}" for b in extra_bus_qst])
    nd = len(Qd)
    # [Q11] sources: G_b, S_b, D_b, then the Case-3 injections interleaved per branch (f_1, t_1, f_2, t_2, ...)
    i_qg_pos = QG > tol
    i_qs_pos = Qsh > tol
    i_qd_neg = QD < -tol
    case3_has = np.column_stack([case3_hasf, case3_hast]).ravel()
    case3_bus = np.column_stack([case3_fbus, case3_tbus]).ravel()[case3_has]
    case3_Q = np.column_stack([case3_QF, case3_QT]).ravel()[case3_has]
    src_bus = np.r_[gbus[i_qg_pos], bus_i[i_qs_pos], bus_i[i_qd_neg], case3_bus].astype(int)
    Qsrc = np.r_[QG[i_qg_pos], Qsh[i_qs_pos], -QD[i_qd_neg], case3_Q]
    n_pre = int(i_qg_pos.sum() + i_qs_pos.sum() + i_qd_neg.sum())
    src_names = ([f"G_{b}" for b in gbus[i_qg_pos]] + [f"S_{b}" for b in bus_i[i_qs_pos]] +
                 [f"D_{b}" for b in bus_i[i_qd_neg]])
    for a, b, hf, ht in zip(case3_fbus, case3_tbus, case3_hasf, case3_hast):
        src_names += ([f"E_{a}_{b}_f"] if hf else []) + ([f"E_{a}_{b}_t"] if ht else [])
    # [Q12] expanded network
    sink = nb + 1
    return ReactivePrep(
        nb=nb, nl=nl, nd=nd, ns=len(Qsrc), n_pre=n_pre,
        f=(np.r_[F, d_bus] - 1).astype(np.int64), t=(np.r_[T, np.full(nd, sink)] - 1).astype(np.int64),
        Qf=np.r_[QF, Qd], Qt=np.r_[QT, Qd], Qsf=np.r_[Qsf, np.zeros(nd)], Qst=np.r_[Qst, np.zeros(nd)],
        src=(src_bus - 1).astype(np.int64), Qsrc=Qsrc, src_names=src_names, d_labels=d_labels,
        case3_fbus=case3_fbus, case3_tbus=case3_tbus, case3_QF=case3_QF, case3_QT=case3_QT,
        case3_hasf=case3_hasf, case3_hast=case3_hast,
        removed={"branches": int(i_zero.sum()), "generators": int(i_zero_gen.sum())},
        flipped_branches=int(flip.sum()),
        converted={"negative_generation_to_demand": int(i_qg_neg.sum()),
                   "negative_demand_to_source": int(i_qd_neg.sum()),
                   "inductive_bus_shunt_to_demand": int(i_qs_neg.sum()),
                   "case2_branches_to_demands": n_case2_both,
                   "case3_branches_to_sources": int(c3.sum()) - n_case3_zt,
                   "inductive_from_end_shunt_to_demand": int(len(extra_bus_qsf)),
                   "inductive_to_end_shunt_to_demand": int(len(extra_bus_qst)),
                   "zero_terminal_demands": n_case2_zt,
                   "zero_terminal_sources": n_case3_zt},
    )


def assemble_reactive(prep, tol=1e-6, tol_div=1e-9):
    """Phase (b): rules [Q13], [Q14], the informational series-flow list of [Q15], [Q16], [Q18] and the
    right-hand side of [Q19]."""
    ne, nb, ns, nl = prep.ne, prep.nb, prep.ns, prep.nl
    rows = np.arange(ne)
    # [Q13] series flows and attenuation
    with np.errstate(divide="ignore", invalid="ignore"):
        Qs = prep.Qf + prep.Qsf
        Qe = prep.Qt - prep.Qst
        aQ = Qe / Qs
    # [Q14] incidence matrices, total inflow
    Cf = sp.csc_matrix((np.ones(ne), (rows, prep.f)), shape=(ne, nb + 1))
    Ct = sp.csc_matrix((np.ones(ne), (rows, prep.t)), shape=(ne, nb + 1))
    Cs = sp.csc_matrix((np.ones(ns), (np.arange(ns), prep.src)), shape=(ns, nb + 1))
    Qtot = Ct.T @ prep.Qt + Cs.T @ prep.Qsrc
    with np.errstate(divide="ignore"):
        h = 1.0 / Qtot
    hf = h[prep.f]
    # [Q15] checks in the order of rpfc.m: inflow guard (always); carriers with a negative series flow
    # (informational, labelled 'k (f-t)'); cycle and spectral radius after A; negative self-term after [Q18]
    labels = [f"{k + 1} ({prep.f[k] + 1}-{prep.t[k] + 1})" for k in range(nl)]
    checks.check_inflow(Qtot, prep.f, tol_div, "rpfc")
    _, i_ser = checks.check_series_flow(Qs[:nl], Qe[:nl], labels, tol)
    neg = [labels[k] for k in i_ser]
    # [Q16] system
    D = sp.diags(prep.Qf * hf, 0, shape=(ne, ne), format="csc")
    M = (Cf @ Ct.T).tocsc()
    N = (-(D @ (M @ sp.diags(aQ, 0, shape=(ne, ne), format="csc")))).tocsc()
    A = (sp.identity(ne, format="csc") + N).tocsc()
    acyclic, rho, w_cyc = checks.check_cycle(N, M, "rpfc")
    B = (D @ (Cf @ Cs.T) @ sp.diags(prep.Qsrc, 0, shape=(ns, ns), format="csc")).tocsc()
    # [Q18] self-terms of the carriers
    self_toend = np.zeros(nl)
    big = np.abs(Qs[:nl]) > tol_div
    self_toend[big] = (Qe[:nl][big] / Qs[:nl][big]) * prep.Qsf[:nl][big] + prep.Qst[:nl][big]
    checks.check_selfterm(self_toend, labels, tol_div)     # guard pfc:rpfc:negativeSelfTerm (always on)
    # [Q19] right-hand side of the self-term solve
    i_self = np.where(self_toend > tol_div)[0]
    n_self = len(i_self)
    Css = sp.csc_matrix((np.ones(n_self), (np.arange(n_self), prep.t[i_self])), shape=(n_self, nb + 1))
    B_self = (D @ (Cf @ Css.T) @ sp.diags(self_toend[i_self], 0, shape=(n_self, n_self), format="csc")).tocsc()
    return ReactiveSystem(M=M, N=N, A=A, B=B, B_self=B_self, i_self=i_self, aQ=aQ, Qs=Qs, Qe=Qe,
                          self_toend=self_toend, negative_series_flow=neg, acyclic=acyclic, spectral_radius=rho,
                          warnings=w_cyc)


def solve_reactive(system):
    """Phase (c): one sparse LU factorisation of A (SuperLU, COLAMD ordering), used for both solves:
    [Q17] X = A \\ B and [Q19] X_self = A \\ B_self, each for all columns at once."""
    ne = system.A.shape[0]
    if system.B.shape[1] == 0 and system.B_self.shape[1] == 0:
        return np.zeros((ne, 0)), np.zeros((ne, 0))
    lu = spla.splu(system.A, permc_spec=PERMC_SPEC)
    X = lu.solve(system.B.toarray()) if system.B.shape[1] else np.zeros((ne, 0))
    X_self = lu.solve(system.B_self.toarray()) if system.B_self.shape[1] else np.zeros((ne, 0))
    return X, X_self


def build_reactive_table(prep, system, X, X_self):
    """Phase (d): [Q17] to-end contributions, [Q19] and [Q20] the self-term block, [Q21] to [Q25] labels,
    row order and the appended Case-3 rows. Returns a DataFrame with index 'Line'."""
    nl, ne, ns = prep.nl, prep.ne, prep.ns
    aQ = system.aQ
    Xt = aQ[:, None] * X
    Xf_self = np.zeros((ne, nl)); Xt_self = np.zeros((ne, nl))
    for k, ii in enumerate(system.i_self):
        Xf_self[:, ii] = X_self[:, k]
        Xt_self[:, ii] = aQ * X_self[:, k]
    for i in range(nl):                                  # [Q20] direct self-term at the own to-end
        Xt_self[i, i] += system.self_toend[i]
    full_f = np.hstack([X, Xf_self]); full_t = np.hstack([Xt, Xt_self])
    # [Q21] column names: sources, then the self-term columns with the parallel scheme
    fb, tb = prep.f[:nl] + 1, prep.t[:nl] + 1
    self_names = [f"E_{a}_{b}_s" for a, b in zip(fb, tb)]
    self_names, _ = disambiguate_parallel(self_names, self_names)
    names = list(prep.src_names) + self_names
    # [Q22] carrier rows with the parallel scheme; data rows: demands, from-ends, to-ends
    lf = [f"L_{{f.{a}-{b}}}" for a, b in zip(fb, tb)]
    lt = [f"L_{{t.{a}-{b}}}" for a, b in zip(fb, tb)]
    lf, lt = disambiguate_parallel(lf, lt)
    data = np.vstack([full_f[nl:], full_f[:nl], full_t[:nl]])
    rows = unique_labels(prep.d_labels) + lf + lt          # [Q23]
    # [Q24] Case-3 rows: all from-end rows, then all to-end rows (a member with a zero terminal has only the row
    # of its non-zero terminal); each in the member's own column
    n3 = prep.n_case3
    if n3:
        hf, ht = np.asarray(prep.case3_hasf, bool), np.asarray(prep.case3_hast, bool)
        col_end = prep.n_pre + np.cumsum(hf.astype(int) + ht.astype(int))   # one past the member's last column
        colf, colt = col_end - 1 - ht.astype(int), col_end - 1
        i_f, i_t = np.where(hf)[0], np.where(ht)[0]
        c3f = np.zeros((len(i_f), ns + nl)); c3t = np.zeros((len(i_t), ns + nl))
        for r, i in enumerate(i_f):
            c3f[r, colf[i]] = -prep.case3_QF[i]
        for r, i in enumerate(i_t):
            c3t[r, colt[i]] = prep.case3_QT[i]
        data = np.vstack([data, c3f, c3t])
        rows += [f"L_{{f.{prep.case3_fbus[i]}-{prep.case3_tbus[i]}}}" for i in i_f]
        rows += [f"L_{{t.{prep.case3_fbus[i]}-{prep.case3_tbus[i]}}}" for i in i_t]
    # [Q25] final unique-label pass
    return pd.DataFrame(data, index=pd.Index(unique_labels(rows), name="Line"), columns=unique_labels(names))


def _nnz(m):
    """Number of non-zero values (as Matlab's nnz), not the number of stored entries."""
    return int(np.count_nonzero(m.data))


def rpfc(case, options=None, tol=None, tol_div=None):
    """RPFC of a solved case (built with branch parameters). ``options`` is an Options (or a dict, or None for
    the defaults); the keyword arguments tol and tol_div win over it. Returns (table, info): the table has the row labels,
    column labels, order and values of the CSV of rpfc.m; info holds the sizes and counts of diag of rpfc.m [Q26],
    the solver, and the self-terms (for the guard of T13c).
    info["removed"]["branches"] counts the branches removed from the exported input (rule R2 / Q6). It can be smaller
    than diag.removed.branches of the Matlab code by the branches that my_runpf.m of the iterative method had already
    deleted before the export (all four flows below 1e-4), for example branch 9-11 of case30; the exported file does
    not record them.
    """
    opts = resolve_options(options, tol, tol_div)
    if not isinstance(case, Case):
        raise PfcError("pfc:input:casedata", "case must be a pfc_matrix.Case (see Case.from_frames).")
    tol, tol_div = opts.tol, opts.tol_div
    prep = prepare_reactive(case, tol)
    system = assemble_reactive(prep, tol, tol_div)
    X, X_self = solve_reactive(system)
    table = build_reactive_table(prep, system, X, X_self)
    info = {
        "case": case.name, "n_bus": prep.nb, "n_branch": prep.nl, "n_demand": prep.nd, "n_source": prep.ns,
        "n_e": prep.ne, "nnz_N": _nnz(system.N), "nnz_A": _nnz(system.A),
        "issparse": {"N": sp.issparse(system.N), "A": sp.issparse(system.A), "B": sp.issparse(system.B)},
        "removed": prep.removed, "flipped_branches": prep.flipped_branches, "converted": prep.converted,
        "condest_A": checks.condest(system.A) if opts.checks else None,
        "attenuation_range": ([float(system.aQ[:prep.nl].min()), float(system.aQ[:prep.nl].max())] if prep.nl else []),
        "acyclic": bool(system.acyclic), "spectral_radius": float(system.spectral_radius),
        "warnings": list(system.warnings),
        "solver": f"scipy.sparse.linalg.splu (SuperLU), permc_spec={PERMC_SPEC}; one factorisation of A "
                  f"used for both solves (main and self-term)",
        "negative_series_flow": {"count": len(system.negative_series_flow),
                                 "branches": list(system.negative_series_flow)},
        "self_toend": [float(v) for v in system.self_toend],
        "opts": opts.as_dict(),
    }
    return table, info
