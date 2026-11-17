"""Validity checks of the Python port: a port of pfc_checks.m (head of branch revision). The identifiers and the
message texts are those of the Matlab code. Errors raise PfcError; warnings are returned as identifiers (the caller
collects them in info["warnings"]) and are issued with warnings.warn(PfcWarning). Nothing is printed.

    kahn(M)                              acyclicity (Kahn's topological sort on the pattern of M)
    cyclic_core(P), spectral_radius(N, acyclic)
    check_cycle(N, M, kind)              pfc:apfc:cyclic (warning), pfc:apfc:spectralRadius, pfc:rpfc:spectralRadius
    check_attenuation(a, labels, slack)  pfc:apfc:nonPositiveAttenuation, pfc:apfc:attenuationAboveOne (warning)
    check_inflow(total, f, tol_div, kind)  pfc:apfc:zeroInflow, pfc:rpfc:zeroInflow
    check_series_flow(qs, qe, labels, tol, warn=False)  informational (pfc:rpfc:negativeSeriesFlow only with warn)
    check_selfterm(s, labels, tol_div)   pfc:rpfc:negativeSelfTerm
    condest(A)                           1-norm condition number estimate (scipy onenormest; Matlab: condest)
Indices are 0-based in Python; labels and messages use the 1-based numbers of the Matlab code.
"""
import warnings

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .options import PfcError, PfcWarning

EIG_DENSE_MAX = 2000      # eig on the dense cyclic core up to this size, eigs above (as pfc_checks.m)
CONDEST_SEED = 0          # seed of the random columns of onenormest (the global numpy state is restored)


def _warn(identifier, message):
    warnings.warn(PfcWarning(identifier, message), stacklevel=3)
    return identifier


def _pattern(M):
    P = sp.csc_matrix(M, copy=True)
    P.data = (P.data != 0).astype(float)
    P.eliminate_zeros()
    P.sort_indices()
    return P


def kahn(M):
    """[acyclic, order] = kahn(M) of pfc_checks.m: Kahn's algorithm on the pattern of M (edge j -> i for
    M(i,j) != 0), with the same stack discipline (last in, first out). ``order`` is 0-based."""
    P = _pattern(M)
    n = P.shape[0]
    indeg = np.asarray(P.sum(axis=1)).ravel().astype(int)
    stack = list(np.where(indeg == 0)[0])
    order = []
    while stack:
        j = stack.pop()
        order.append(int(j))
        for i in P.indices[P.indptr[j]:P.indptr[j + 1]]:      # rows of column j, ascending (as find)
            indeg[i] -= 1
            if indeg[i] == 0:
                stack.append(i)
    return len(order) == n, np.array(order, dtype=int)


def cyclic_core(P):
    """Nodes left after repeatedly removing nodes without incoming edges (in the remaining graph) and nodes
    without outgoing edges (0-based indices)."""
    P = _pattern(P)
    keep = np.ones(P.shape[0], dtype=bool)
    changed = True
    while changed:
        idx = np.where(keep)[0]
        Q = P[idx][:, idx]
        no_in = np.asarray(Q.sum(axis=1)).ravel() == 0
        no_out = np.asarray(Q.sum(axis=0)).ravel() == 0
        drop = no_in | no_out
        changed = bool(drop.any())
        keep[idx[drop]] = False
    return np.where(keep)[0]


def spectral_radius(N, acyclic=None):
    """0 exactly if acyclic; otherwise the largest absolute eigenvalue of the cyclic core of N: eig on the dense
    core up to EIG_DENSE_MAX nodes, eigs (largest magnitude) above."""
    if acyclic is None:
        acyclic, _ = kahn(N)
    if acyclic:
        return 0.0
    core = cyclic_core(N)
    if core.size == 0:
        return 0.0       # not reachable when acyclic is computed from the pattern of N (the core is then non-empty)
    Nc = sp.csc_matrix(N)[core][:, core]
    if core.size <= EIG_DENSE_MAX:
        return float(np.max(np.abs(np.linalg.eigvals(Nc.toarray()))))
    return float(np.abs(spla.eigs(Nc, k=1, which="LM", return_eigenvectors=False)[0]))


def check_cycle(N, M, kind, warn_cyclic=True):
    """[acyclic, rho, w] = check_cycle(N, M, kind) of pfc_checks.m. With warn_cyclic=False (Options.checks=False)
    the warning pfc:apfc:cyclic is not issued; the errors are always raised."""
    w = []
    acyclic, _ = kahn(M)
    rho = spectral_radius(N, acyclic)
    if acyclic:
        return acyclic, rho, w
    if kind == "apfc":
        if rho < 1:
            if warn_cyclic:
                w.append(_warn("pfc:apfc:cyclic",
                               "The active power flow digraph is cyclic (spectral radius of N = %.17g < 1); "
                               "A = I + N is invertible by the Neumann series, but the acyclicity theorem does not apply." % rho))
        else:
            raise PfcError("pfc:apfc:spectralRadius",
                           "The active power flow digraph is cyclic and the spectral radius of N is "
                           "%.17g >= 1: the sufficient condition rho(N) < 1 for the invertibility of A = I + N is not met." % rho)
    elif kind == "rpfc":
        if rho >= 1:
            raise PfcError("pfc:rpfc:spectralRadius",
                           "The reactive power flow digraph is cyclic and the spectral radius of N^Q is "
                           "%.17g >= 1: the sufficient condition rho(N^Q) < 1 for the invertibility of A^Q = I + N^Q "
                           "(Neumann series) is not met." % rho)
    else:
        raise PfcError("pfc:checks:kind", f"Unknown kind '{kind}'.")
    return acyclic, rho, w


def check_attenuation(a, labels, slack):
    """W = check_attenuation(a, labels, slack) of pfc_checks.m: error for a <= 0, warning for a > 1 + slack."""
    a = np.asarray(a, dtype=float)
    w = []
    bad = np.where(a <= 0)[0]
    if bad.size:
        raise PfcError("pfc:apfc:nonPositiveAttenuation",
                       "Active attenuation Pt/Pf <= 0 (terminal flow at the tolerance; check the input) for branch(es): %s."
                       % ", ".join(labels[k] for k in bad))
    hi = np.where(a > 1 + slack)[0]
    if hi.size:
        w.append(_warn("pfc:apfc:attenuationAboveOne",
                       "Active attenuation Pt/Pf exceeds 1 by more than %.3g for branch(es): %s."
                       % (slack, ", ".join(labels[k] for k in hi))))
    return w


def check_inflow(total, f, tol_div, kind):
    """check_inflow(total, f, tol_div, kind) of pfc_checks.m: error if a bus with an outgoing branch (0-based
    from-bus indices f) has a total inflow below tol_div."""
    total = np.asarray(total, dtype=float).ravel()
    has_out = np.zeros(total.size, dtype=bool)
    has_out[np.unique(np.asarray(f, dtype=int))] = True
    bad = np.where(has_out & (total < tol_div))[0]
    if bad.size:
        raise PfcError(f"pfc:{kind}:zeroInflow",
                       "Total inflow below %.3g at bus(es) %s, which have an outgoing branch: the shares cannot be computed."
                       % (tol_div, ", ".join(str(k + 1) for k in bad)))


def check_negative_loss(labels_to, labels_from):
    """check_negative_loss(labels_to, labels_from) of pfc_checks.m: error pfc:apfc:negativeLoss for active power
    delivered at one end of a branch whose other end carries zero flow: labels_to names the (0,+) branches,
    labels_from the (-,0) branches; when both occur, the two sentences are joined by one space."""
    msg = []
    if labels_to:
        msg.append("Active power is delivered at the to-end of a branch with zero from-end flow (negative loss) "
                   "for branch(es): %s." % ", ".join(labels_to))
    if labels_from:
        msg.append("Active power is delivered at the from-end of a branch with zero to-end flow (negative loss) "
                   "for branch(es): %s." % ", ".join(labels_from))
    if msg:
        raise PfcError("pfc:apfc:negativeLoss", " ".join(msg))


def check_series_flow(qs, qe, labels, tol, warn=False):
    """[w, idx] = check_series_flow(qs, qe, labels, tol, warn) of pfc_checks.m (idx 0-based)."""
    idx = np.where((np.asarray(qs) < -tol) | (np.asarray(qe) < -tol))[0]
    w = []
    if warn and idx.size:
        w.append(_warn("pfc:rpfc:negativeSeriesFlow",
                       "Negative series reactive flow (negative attenuation) in branch(es): %s."
                       % ", ".join(labels[k] for k in idx)))
    return w, idx


def check_selfterm(s, labels, tol_div):
    """check_selfterm(s, labels, tol_div) of pfc_checks.m: error if a to-end shunt self-term is below -tol_div."""
    bad = np.where(np.asarray(s, dtype=float).ravel() < -tol_div)[0]
    if bad.size:
        raise PfcError("pfc:rpfc:negativeSelfTerm",
                       "To-end shunt self-term below -%.3g (it would not be propagated, and the downstream row sums "
                       "would not equal the flows) for branch(es): %s." % (tol_div, ", ".join(labels[k] for k in bad)))


def condest(A):
    """Estimate of the 1-norm condition number of A: onenormest(A) * onenormest(A^-1), with A^-1 applied through a
    sparse LU factorisation of A (its own factorisation). onenormest draws random columns from the global numpy
    generator; it is seeded with CONDEST_SEED and the previous state is restored, so the value is reproducible.
    Matlab's condest is a different estimator; the two values may differ (compared within a factor of 10)."""
    A = sp.csc_matrix(A)
    n = A.shape[0]
    if n == 0:
        return 0.0
    lu = spla.splu(A, permc_spec="COLAMD")
    inv = spla.LinearOperator((n, n), dtype=float,
                              matvec=lambda x: lu.solve(np.asarray(x, dtype=float).ravel()),
                              rmatvec=lambda x: lu.solve(np.asarray(x, dtype=float).ravel(), trans="T"))
    state = np.random.get_state()
    try:
        np.random.seed(CONDEST_SEED)
        return float(spla.onenormest(A) * spla.onenormest(inv))
    finally:
        np.random.set_state(state)
