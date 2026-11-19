"""Independent baseline (brief T14h): the matrix tracing method of Achayuthakan, Dent, Bialek and Ongsakul,
IEEE Trans. Power Syst. 25(2), 2010, for active power at bus level. Specification: SPEC.md in this folder.

Only numpy, scipy and pandas are used; nothing of the matrix method (pfc_matrix, apfc.m) is imported, so the comparison
is independent. The input is the solved power flow as the DataFrames bus (BUS_I, PD), gen (GEN_BUS, PG) and branch
(F_BUS, T_BUS, PF, PT with PT the flow delivered at the to-end), as load_results of the iterative method returns them.

Network (decided for the lossy case, SPEC.md section 6 and the decision after stop condition 1 of T14h):
- signs by tol (default 1e-6) of PF and the delivered PT; (0,0) branches are left out;
- (+,+): a line from F to T, sending-end flow PF, receiving-end flow PT; (-,-): a line from T to F, sending-end
  flow -PT, receiving-end flow -PF;
- (+,-) (fed from both ends), (+,0) and (0,-) (one zero end): not lines of the digraph, but demands at their own ends
  (PF at F, -PT at T), as the Case-2 rows of apfc.m. The paper does not cover them;
- (-,+), (0,+) and (-,0) (power delivered at an end of a branch without a source): an error (not in the paper, and
  apfc.m rejects them too);
- parallel lines between the same ordered pair of buses are added (eq. (4) of the paper adds them);
- nodal generation: the sum of PG > tol of the generators at a bus; nodal demands: PD > tol of the bus and the
  terminal demands above; negative generation or demand is an error (not covered by the paper; not in the networks of T14h).

Variants:
- "net-downstream": A_d from the receiving-end flows, P = nodal inflow (P_G + received flows), eq. (14):
  c_ik = (P_G,i / P_i) [A_d^-1]_ik P_D,k (generator bus i to demand k at bus k). The contributions to a demand sum
  exactly to the demand (net: the losses stay with the generators).
- "gross-upstream": A_u from the sending-end flows, P = nodal outflow (P_D + sent flows), the upstream counterpart:
  c_ik = (P_D,k / P_k) [A_u^-1]_ki P_G,i. The contributions of a generator sum exactly to its generation, and the
  contributions to a demand sum to the demand plus its allocated losses (gross).
Each with two solution methods: "inverse" (dense A, explicit inverse, as published) and "sparse" (one sparse LU of A,
solved for the generator buses as several right-hand sides).
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.sparse.linalg as spla

TOL = 1e-6


@dataclass
class Network:
    buses: np.ndarray            # bus numbers of the nodes kept (nodal flow > 0)
    pg: np.ndarray               # nodal generation per node
    demands: list                # (label, node index, value)
    send: dict                   # (i, j) node indices -> sending-end flow of the lines from i to j (summed)
    recv: dict                   # (i, j) -> receiving-end flow
    branch_demands: int          # number of terminal demands from branches
    lines: int                   # number of branches used as lines


def build_network(bus, gen, branch, tol=TOL):
    s = lambda x: np.where(x > tol, 1, np.where(x < -tol, -1, 0))
    PF, PT = branch["PF"].to_numpy(float), branch["PT"].to_numpy(float)
    F, T = branch["F_BUS"].to_numpy(int), branch["T_BUS"].to_numpy(int)
    sf, st = s(PF), s(PT)
    bad = ((sf < 0) & (st > 0)) | ((sf == 0) & (st > 0)) | ((sf < 0) & (st == 0))
    if bad.any():
        raise ValueError(f"{int(bad.sum())} branch(es) deliver power at an end without a source (patterns (-,+), (0,+), (-,0))")
    if (bus["PD"].to_numpy(float) < -tol).any() or (gen["PG"].to_numpy(float) < -tol).any():
        raise ValueError("negative demand or generation is not covered by the method")
    all_buses = bus["BUS_I"].to_numpy(int)
    pos = {b: k for k, b in enumerate(all_buses)}
    n = len(all_buses)
    pg = np.zeros(n)
    for b, p in zip(gen["GEN_BUS"].to_numpy(int), gen["PG"].to_numpy(float)):
        if p > tol:
            pg[pos[b]] += p
    demands = [(f"D_{b}", pos[b], float(p)) for b, p in zip(all_buses, bus["PD"].to_numpy(float)) if p > tol]
    send, recv = {}, {}
    nbd = nl = 0
    for k in range(len(PF)):
        f, t = pos[F[k]], pos[T[k]]
        if sf[k] > 0 and st[k] > 0:
            key, se, re_ = (f, t), PF[k], PT[k]
        elif sf[k] < 0 and st[k] < 0:
            key, se, re_ = (t, f), -PT[k], -PF[k]
        else:
            if sf[k] > 0:                       # (+,-) and (+,0): demand at the from-end
                demands.append((f"C2_{F[k]}", f, float(PF[k]))); nbd += 1
            if st[k] < 0:                       # (+,-) and (0,-): demand at the to-end
                demands.append((f"C2_{T[k]}", t, float(-PT[k]))); nbd += 1
            continue
        send[key] = send.get(key, 0.0) + se
        recv[key] = recv.get(key, 0.0) + re_
        nl += 1
    return Network(all_buses, pg, demands, send, recv, nbd, nl), n


def _matrix(net, n, kind):
    """The distribution matrix (sparse, order n) and the nodal flow P.
    kind 'net': A_d with receiving-end flows, P = inflow; [A_d]_ij = -recv(i, j) / P_j for a line i -> j.
    kind 'gross': A_u with sending-end flows, P = outflow; [A_u]_ij = -send(j, i) / P_j for a line j -> i."""
    if kind == "net":
        P = net.pg.copy()
        for (i, j), v in net.recv.items():
            P[j] += v
        rows, cols, vals = zip(*[(i, j, -v / P[j]) for (i, j), v in net.recv.items()]) if net.recv else ((), (), ())
    else:
        P = np.zeros(n)
        for _, k, v in net.demands:
            P[k] += v
        for (i, j), v in net.send.items():
            P[i] += v
        rows, cols, vals = zip(*[(j, i, -v / P[i]) for (i, j), v in net.send.items()]) if net.send else ((), (), ())
    keep = P > TOL
    if not keep.all():
        # nodes with zero nodal flow are removed (p. 1080); they must not carry any line
        dead = set(np.where(~keep)[0])
        if any(i in dead or j in dead for i, j in net.recv):
            raise ValueError("a node with zero nodal flow carries a line")
    A = sp.csc_matrix((np.r_[np.ones(n), np.asarray(vals, float)],
                       (np.r_[np.arange(n), np.asarray(rows, int)], np.r_[np.arange(n), np.asarray(cols, int)])), shape=(n, n))
    idx = np.where(keep)[0]
    return A[idx][:, idx].tocsc(), P, idx


def contributions(bus, gen, branch, kind="net", method="sparse", tol=TOL, return_info=False):
    """Contribution of every generator bus (columns G_<bus>) to every demand (rows D_<bus> for the bus demands and
    C2_<bus> for the terminal demands of branches, in the order of the bus table, then of the branch table)."""
    net, n = build_network(bus, gen, branch, tol)
    A, P, idx = _matrix(net, n, kind)
    loc = {k: m for m, k in enumerate(idx)}
    gbus = [k for k in idx if net.pg[k] > tol]
    gcol = [loc[k] for k in gbus]
    m = len(idx)
    if method == "inverse":
        Ainv = np.linalg.inv(A.toarray())
        # net: rows of A_d^-1 at the generator buses; gross: columns of A_u^-1 at the generator buses
        G = Ainv[gcol, :].T if kind == "net" else Ainv[:, gcol]
    else:
        E = np.zeros((m, len(gcol))); E[gcol, np.arange(len(gcol))] = 1.0
        if kind == "net":
            G = spla.splu(A.T.tocsc()).solve(E)      # columns: rows of A_d^-1 at the generator buses
        else:
            G = spla.splu(A).solve(E)                # columns of A_u^-1 at the generator buses
    # G[m, g]: net -> [A_d^-1]_{gen g, node m}; gross -> [A_u^-1]_{node m, gen g}
    data = np.zeros((len(net.demands), len(gbus)))
    for r, (_, k, v) in enumerate(net.demands):
        mk = loc[k]
        if kind == "net":
            data[r] = (net.pg[gbus] / P[gbus]) * G[mk] * v
        else:
            data[r] = (v / P[k]) * G[mk] * net.pg[gbus]
    table = pd.DataFrame(data, index=pd.Index([d[0] for d in net.demands], name="Line"),
                         columns=[f"G_{net.buses[k]}" for k in gbus])
    if not return_info:
        return table
    nnz = int(A.nnz)
    info = {"kind": kind, "method": method, "n_b": m, "nnz_A": nnz, "density": nnz / m ** 2, "lines": net.lines,
            "branch_demands": net.branch_demands, "generator_buses": len(gbus), "demands": len(net.demands),
            "P": P[idx], "A": A, "pg": net.pg[idx], "demand_values": np.array([d[2] for d in net.demands]),
            "demand_nodes": np.array([loc[d[1]] for d in net.demands])}
    return table, info


def condest_1norm(A):
    """1-norm condition estimate of a sparse matrix: ||A||_1 times the onenormest of A^-1 (one sparse LU), as the
    matrix method's condest_A and Matlab's condest."""
    lu = spla.splu(sp.csc_matrix(A))
    n = A.shape[0]
    Ainv = spla.LinearOperator((n, n), matvec=lu.solve, rmatvec=lambda x: lu.solve(x, trans="T"), dtype=float)
    return float(spla.norm(sp.csc_matrix(A), 1) * spla.onenormest(Ainv))
