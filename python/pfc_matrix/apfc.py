"""Active power flow contributions (APFC): Python port of apfc.m (branch revision of the matrix repository).

The computation is split into four phases, which can be called (and later timed) on their own:
    prepare_active(case, tol)      (a) preprocessing; index and flow arrays
    assemble_active(prep)          (b) assembly of M, N, A and B (scipy.sparse)
    solve_active(system)           (c) sparse LU factorisation of A and solve for all columns of B
    build_active_table(prep, X)    (d) output table (pandas DataFrame) and labels
``apfc(case, options)`` calls them in this order. No phase prints, times or writes files; the checks of apfc.m
(module checks) run inside phase (b), raise PfcError or issue PfcWarning.

The rules are those of apfc.m; the numbered rules of the T13a report are referenced as [Rn].
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


@dataclass
class ActivePrep:
    """Result of phase (a): the expanded network (physical branches followed by demand branches)."""
    nb: int                 # number of buses; the virtual sink has the 0-based index nb
    nl: int                 # number of physical branches after deletion
    nd: int                 # number of demand branches
    ng: int                 # number of source columns
    f: np.ndarray           # 0-based from-bus index of every branch of the expanded network (length ne)
    t: np.ndarray           # 0-based to-bus index (the sink for demand branches)
    f_bus: np.ndarray       # bus numbers of the physical branches after orientation (for labels)
    t_bus: np.ndarray
    Pf: np.ndarray          # from-end flows (MW)
    Pt: np.ndarray          # to-end flows (MW)
    g: np.ndarray           # 0-based bus index of every source column
    Pg: np.ndarray          # source injections (MW)
    src_bus: np.ndarray     # bus numbers of the sources (for labels)
    n_gpos: int             # number of generator sources (the rest are negative demands)
    d_bus: np.ndarray       # bus numbers of the demands (for labels)
    nd_pd: int              # number of demands from PD (then nd_pg generators with negative output, then C2_b)
    removed: dict
    flipped_branches: int
    converted: dict
    nd_pg: int = 0          # number of demands from generators with negative output (the rest are terminal demands C2_b)
    negative_loss_to: list = field(default_factory=list)    # labels of the (0,+) branches (error in phase (b))
    negative_loss_from: list = field(default_factory=list)  # labels of the (-,0) branches (error in phase (b))

    @property
    def ne(self):
        return self.nl + self.nd


def terminal_sign(x, tol):
    """Sign of terminal flows by tol: +1 if x > tol, -1 if x < -tol, 0 otherwise."""
    return (x > tol).astype(int) - (x < -tol).astype(int)


@dataclass
class ActiveSystem:
    """Result of phase (b)."""
    M: sp.csc_matrix
    N: sp.csc_matrix
    A: sp.csc_matrix
    B: sp.csc_matrix
    a: np.ndarray           # attenuation factors Pt / Pf of the expanded network
    acyclic: bool = True    # [R20] checks: acyclicity of the flow digraph (Kahn on M)
    spectral_radius: float = 0.0
    warnings: list = field(default_factory=list)


def prepare_active(case, tol=1e-6):
    """Phase (a): preprocessing of apfc.m on the exported input, and the index and flow arrays."""
    bus, gen, br = case.bus, case.gen, case.branch
    nb = len(bus)
    # [R2] delete branches with both active flows below tol; [R3] delete generators with |PG| below tol
    pf, pt = br["PF"].to_numpy(float), br["PT"].to_numpy(float)
    del_br = (np.abs(pf) < tol) & (np.abs(pt) < tol)
    fb, tb = br["F_BUS"].to_numpy(int)[~del_br], br["T_BUS"].to_numpy(int)[~del_br]
    pf, pt = pf[~del_br], pt[~del_br]
    pg_all = gen["PG"].to_numpy(float)
    del_gen = np.abs(pg_all) < tol
    gbus, pg = gen["GEN_BUS"].to_numpy(int)[~del_gen], pg_all[~del_gen]
    # [R4] the to-end flows of the input are already the delivered flows (sign reversed by my_runpf.m)
    # Terminal classification (T15b), before the flip, sign by tol: (+,0) one terminal demand PF at the from bus,
    # (0,-) one terminal demand -PT at the to bus (rows C2_b); both leave the branch set. (0,+) and (-,0) are the
    # error pfc:apfc:negativeLoss, raised at the start of phase (b); labels 'k (f-t)' with k the 1-based index after
    # the deletion and f-t as in the input. Every other branch keeps the path below.
    s_f, s_t = terminal_sign(pf, tol), terminal_sign(pt, tol)
    lab_in = [f"{k + 1} ({a}-{b})" for k, (a, b) in enumerate(zip(fb, tb))]
    neg_to = [lab_in[k] for k in np.where((s_f == 0) & (s_t > 0))[0]]
    neg_from = [lab_in[k] for k in np.where((s_f < 0) & (s_t == 0))[0]]
    zt_from = (s_f > 0) & (s_t == 0)
    zt_to = (s_f == 0) & (s_t < 0)
    zt_bus = np.r_[fb[zt_from], tb[zt_to]]
    zt_P = np.r_[pf[zt_from], -pt[zt_to]]
    keep = ~(zt_from | zt_to)
    fb, tb, pf, pt = fb[keep], tb[keep], pf[keep], pt[keep]
    # [R5] reverse the branches with PF < 0: swap the ends, PF <- -PT, PT <- -PF
    flip = pf < 0
    fb2, tb2 = np.where(flip, tb, fb), np.where(flip, fb, tb)
    pf2, pt2 = np.where(flip, -pt, pf), np.where(flip, -pf, pt)
    # [R6] demands: buses with PD > tol, then generators with PG < -tol (as demands)
    bus_i, pd_ = bus["BUS_I"].to_numpy(int), bus["PD"].to_numpy(float)
    i_pd = pd_ > tol
    i_pgneg = pg < -tol
    d_bus = np.r_[bus_i[i_pd], gbus[i_pgneg], zt_bus].astype(int)     # then the terminal demands (C2_b)
    Pd = np.r_[pd_[i_pd], -pg[i_pgneg], zt_P]
    # [R7] sources: generators with PG > tol, then buses with PD < -tol
    i_gpos = pg > tol
    i_pdneg = pd_ < -tol
    src_bus = np.r_[gbus[i_gpos], bus_i[i_pdneg]]
    Pg = np.r_[pg[i_gpos], -pd_[i_pdneg]]
    nl, nd = len(pf2), len(Pd)
    # [R8] expanded network: physical branches, then demand branches from the demand bus to the sink
    f = np.r_[fb2, d_bus] - 1
    t = np.r_[tb2, np.full(nd, nb + 1)] - 1
    return ActivePrep(
        nb=nb, nl=nl, nd=nd, ng=len(Pg), f=f.astype(np.int64), t=t.astype(np.int64),
        f_bus=fb2, t_bus=tb2, Pf=np.r_[pf2, Pd], Pt=np.r_[pt2, Pd],
        g=(src_bus - 1).astype(np.int64), Pg=Pg, src_bus=src_bus, n_gpos=int(i_gpos.sum()),
        d_bus=d_bus, nd_pd=int(i_pd.sum()),
        removed={"branches": int(del_br.sum()), "generators": int(del_gen.sum())},
        flipped_branches=int(flip.sum()),
        converted={"negative_demand_to_source": int(i_pdneg.sum()),
                   "negative_generation_to_demand": int(i_pgneg.sum()),
                   "zero_terminal_demands": int(len(zt_bus)),
                   "zero_terminal_sources": 0},
        nd_pg=int(i_pgneg.sum()), negative_loss_to=neg_to, negative_loss_from=neg_from,
    )


def assemble_active(prep, tol_div=1e-9, run_checks=True):
    """Phase (b): incidence matrices, total inflows, attenuation, and the sparse matrices M, N, A, B, with the
    checks of apfc.m in its order [R20]: negative loss (always; raised by apfc.m at the classification, before all
    other checks), inflow guard (always), attenuation range (if run_checks), cycle and spectral radius (the error
    always, the warning pfc:apfc:cyclic if run_checks)."""
    checks.check_negative_loss(prep.negative_loss_to, prep.negative_loss_from)
    ne, nb, ng = prep.ne, prep.nb, prep.ng
    rows = np.arange(ne)
    # [R9] Cf, Ct: branch -> bus incidence (ne x nb+1, the last column is the sink); Cg: source -> bus
    Cf = sp.csc_matrix((np.ones(ne), (rows, prep.f)), shape=(ne, nb + 1))
    Ct = sp.csc_matrix((np.ones(ne), (rows, prep.t)), shape=(ne, nb + 1))
    Cg = sp.csc_matrix((np.ones(ng), (np.arange(ng), prep.g)), shape=(ng, nb + 1))
    # [R10] total inflow per bus, reciprocal h, attenuation a = Pt / Pf
    Ptot = Ct.T @ prep.Pt + Cg.T @ prep.Pg
    with np.errstate(divide="ignore", invalid="ignore"):
        h = 1.0 / Ptot
        a = prep.Pt / prep.Pf
    hf = h[prep.f]
    # [R20] checks: inflow guard, attenuation range of the physical branches (labels 'k (f-t)', 1-based)
    nl = prep.nl
    labels = [f"{k + 1} ({prep.f[k] + 1}-{prep.t[k] + 1})" for k in range(nl)]
    checks.check_inflow(Ptot, prep.f, tol_div, "apfc")
    w_att = checks.check_attenuation(a[:nl], labels, tol_div) if run_checks else []
    # [R11] M = Cf Ct' (M(i,j) = 1 if branch j ends at the from-bus of branch i); N = -diag(Pf.hf) M diag(a)
    M = (Cf @ Ct.T).tocsc()
    D = sp.diags(prep.Pf * hf, 0, shape=(ne, ne), format="csc")
    N = (-(D @ M @ sp.diags(a, 0, shape=(ne, ne), format="csc"))).tocsc()
    # [R12] A = I + N; [R13] B = diag(Pf.hf) S diag(Pg), S = Cf Cg'
    A = (sp.identity(ne, format="csc") + N).tocsc()
    acyclic, rho, w_cyc = checks.check_cycle(N, M, "apfc", warn_cyclic=run_checks)
    S = Cf @ Cg.T
    B = (D @ S @ sp.diags(prep.Pg, 0, shape=(ng, ng), format="csc")).tocsc()
    return ActiveSystem(M=M, N=N, A=A, B=B, a=a, acyclic=acyclic, spectral_radius=rho, warnings=w_att + w_cyc)


def solve_active(system):
    """Phase (c): sparse LU factorisation of A (SuperLU, COLAMD ordering) and solve for all columns of B.
    [R14] X = A \\ B (from-end contributions of every branch of the expanded network)."""
    if system.B.shape[1] == 0:
        return np.zeros(system.B.shape)
    lu = spla.splu(system.A, permc_spec=PERMC_SPEC)
    return lu.solve(system.B.toarray())


def build_active_table(prep, X, a):
    """Phase (d): [R15] to-end contributions X.a; [R16] rows: demands, from-ends, to-ends;
    [R17]-[R19] labels and unique-label rule. Returns a DataFrame with index 'Line'."""
    nl, nd = prep.nl, prep.nd
    Xt = X * a[:, None]
    data = np.vstack([X[nl:], X[:nl], Xt[:nl]])
    cols = [f"G_{b}" for b in prep.src_bus[:prep.n_gpos]] + [f"D_{b}" for b in prep.src_bus[prep.n_gpos:]]
    k1, k2 = prep.nd_pd, prep.nd_pd + prep.nd_pg
    rows = ([f"D_{b}" for b in prep.d_bus[:k1]] + [f"G_{b}" for b in prep.d_bus[k1:k2]] +
            [f"C2_{b}" for b in prep.d_bus[k2:]])
    lf = [f"L_{{f.{i}-{j}}}" for i, j in zip(prep.f_bus, prep.t_bus)]
    lt = [f"L_{{t.{i}-{j}}}" for i, j in zip(prep.f_bus, prep.t_bus)]
    lf, lt = disambiguate_parallel(lf, lt)
    rows = unique_labels(rows + lf + lt)
    cols = unique_labels(cols)
    return pd.DataFrame(data, index=pd.Index(rows, name="Line"), columns=cols)


def _nnz(m):
    """Number of non-zero values (as Matlab's nnz), not the number of stored entries."""
    return int(np.count_nonzero(m.data))


def apfc(case, options=None, tol=None, tol_div=None):
    """APFC of a solved case. ``options`` is an Options (or a dict, or None for the defaults); the keyword
    arguments tol and tol_div win over it. Returns (table, info): the table has the row labels, column labels,
    order and values of the CSV of apfc.m; info holds the fields of diag of apfc.m.
    info["removed"]["branches"] counts the branches removed from the exported input (rule R2 / Q6). It can be smaller
    than diag.removed.branches of the Matlab code by the branches that my_runpf.m of the iterative method had already
    deleted before the export (all four flows below 1e-4), for example branch 9-11 of case30; the exported file does
    not record them.
    """
    opts = resolve_options(options, tol, tol_div)
    if not isinstance(case, Case):
        raise PfcError("pfc:input:casedata", "case must be a pfc_matrix.Case (see Case.from_frames).")
    prep = prepare_active(case, opts.tol)
    system = assemble_active(prep, opts.tol_div, opts.checks)
    X = solve_active(system)
    table = build_active_table(prep, X, system.a)
    info = {
        "case": case.name, "n_bus": prep.nb, "n_branch": prep.nl, "n_demand": prep.nd, "n_source": prep.ng,
        "n_e": prep.ne, "nnz_N": _nnz(system.N), "nnz_A": _nnz(system.A),
        "issparse": {"N": sp.issparse(system.N), "A": sp.issparse(system.A), "B": sp.issparse(system.B)},
        "removed": prep.removed, "flipped_branches": prep.flipped_branches, "converted": prep.converted,
        "condest_A": checks.condest(system.A) if opts.checks else None,
        "attenuation_range": ([float(system.a[:prep.nl].min()), float(system.a[:prep.nl].max())] if prep.nl else []),
        "acyclic": bool(system.acyclic), "spectral_radius": float(system.spectral_radius),
        "warnings": list(system.warnings),
        "solver": f"scipy.sparse.linalg.splu (SuperLU), permc_spec={PERMC_SPEC}",
        "opts": opts.as_dict(),
    }
    return table, info
