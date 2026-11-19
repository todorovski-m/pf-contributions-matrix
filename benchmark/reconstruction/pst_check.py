"""Phase-shifter check (brief T14g, part C) on the cases of benchmark/reconstruction/cases/ (syn9_pst, syn9_pst_tap, and the T16 cases syn9_pst_carrier, syn9_pst_carrier_tap).

For each case:
 (a) the port (python/pfc_matrix, on the exported flows) against Octave (apfc.m and rpfc.m of this working tree on the
     case file, run_octave_cases.sh): row and column labels equal, every cell within 1e-10;
 (b) the row sums of the port's tables against the flows, separately for the demand rows, the from-end rows (P^f or
     Q^f), the to-end rows (P^t or Q^t) and, for RPFC, the Case-3 rows; flows in MW or Mvar, as in the export;
 (c) an independent recomputation of the terminal flows of the phase-shifting branch from the 2x2 stamp and the
     exported bus voltages: tau = TAP e^(j SHIFT) (TAP 0 read as 1), ys = 1/(R + jX),
     Ytt = ys + jB/2, Yff = Ytt/|tau|^2, Yft = -ys/conj(tau), Ytf = -ys/tau, I_f = Yff V_f + Yft V_t,
     I_t = Ytf V_f + Ytt V_t, S_f = V_f conj(I_f) baseMVA, S_t = V_t conj(I_t) baseMVA. Compared with the exported
     PF, QF and the delivered PT, QT (= -P_t, -Q_t of MATPOWER). The complex power across the ideal transformer
     (V_f conj(I_f) against V_i conj(I_i), V_i = V_f/tau, I_i = conj(tau) I_f) is conserved. The port's from-end
     shunt flow (branch_shunt_flows, generation positive) equals (B/2)|V_f|^2/|tau|^2 baseMVA, that is
     -Q^sf in the brief's sign (consumption positive), and its to-end shunt flow (B/2)|V_t|^2 baseMVA; the series flows
     Q^s = Q^f + Qsf(port) = Q^f - Q^sf(brief) and Q^e = Q^t - Qst(port) equal Im(V_i conj(I_s)) and
     Im(V_t conj(I_s)) baseMVA of the series current I_s = ys (V_i - V_t);
 (d) the old method (pinned submodule, in memory, corrected QSF of the timing harness, iteration cap 10000) against
     the port where it finishes, with the row mapping of benchmark/timing/correctness.py.
Prints Markdown. Usage (from the matrix root):
    python -B benchmark/reconstruction/pst_check.py <export dir> <octave dir>
"""
import os
import sys
import warnings

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
MATRIX = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(MATRIX, "benchmark", "timing"))
sys.path.insert(0, os.path.join(MATRIX, "python"))
import cases as C            # noqa: E402
import timing as T           # noqa: E402
import correctness as K      # noqa: E402
from pfc_matrix import Case, read_branch_params, apfc, rpfc, Options, prepare_active, prepare_reactive  # noqa: E402
from pfc_matrix.reactive import branch_shunt_flows  # noqa: E402

CASES = ["syn9_pst", "syn9_pst_tap", "syn9_pst_carrier", "syn9_pst_carrier_tap"]
TOL = 1e-10

g = lambda x: f"{x:.17g}"


def load_case_file(name):
    """baseMVA and (from, to, R, X, B, TAP, SHIFT) of the phase-shifting branch (the one branch with SHIFT != 0), read
    from the case file through Octave-free parsing of the MATPOWER text (the branch matrix rows are plain numbers)."""
    txt = open(os.path.join(HERE, "cases", f"{name}.m")).read()
    base = float(txt.split("mpc.baseMVA = ")[1].split(";")[0])
    rows = txt.split("mpc.branch = [")[1].split("];")[0].strip().splitlines()
    for r in rows:
        v = [float(x) for x in r.replace(";", " ").split()]
        if v[9] != 0:
            return base, dict(F=int(v[0]), T=int(v[1]), R=v[2], X=v[3], B=v[4], TAP=v[8], SHIFT=v[9])
    raise ValueError("no branch with SHIFT != 0")


def main(exp, octd):
    out = []
    for name in CASES:
        data = C.load(name, [exp])
        case = Case.from_frames(data["bus"], data["gen"], data["branch"], data["case"], params=read_branch_params(data["path"]))
        out += [f"## {name}", ""]
        base, b = load_case_file(name)
        F_, T_ = b["F"], b["T"]
        out += [f"Branch {F_}-{T_}: R {b['R']}, X {b['X']}, B {b['B']} pu, TAP {b['TAP']}, SHIFT {b['SHIFT']} degrees; baseMVA {base}.", ""]
        # (a) and (b)
        out += ["| kind | (a) labels equal | (a) max abs diff port - Octave | (a) cells above 1e-10 | "
                "(b) demand rows | (b) from-end rows (P^f or Q^f) | (b) to-end rows (P^t or Q^t) | (b) Case-3 rows | rows |",
                "|---|---|---|---|---|---|---|---|---|"]
        fails = 0
        for kind in ("apfc", "rpfc"):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                t, _ = (apfc if kind == "apfc" else rpfc)(case, Options(checks=True))
            o = pd.read_csv(os.path.join(octd, f"{name}_{kind}.csv"), float_precision="round_trip", dtype={"Line": str}).set_index("Line")
            lab = list(t.index) == list(o.index) and list(t.columns) == list(o.columns)
            d = np.abs(t.to_numpy() - o.to_numpy()) if lab else np.array([np.inf])
            if kind == "apfc":
                p = prepare_active(case, 1e-6); nd, nl = p.nd, p.nl
                flows = np.r_[p.Pf[nl:], p.Pf[:nl], p.Pt[:nl]]; n3 = 0
            else:
                p = prepare_reactive(case, 1e-6); nd, nl = p.nd, p.nl
                hf, ht = np.asarray(p.case3_hasf, bool), np.asarray(p.case3_hast, bool)
                flows = np.r_[p.Qf[nl:], p.Qf[:nl], p.Qt[:nl], -p.case3_QF[hf], p.case3_QT[ht]]; n3 = int(hf.sum() + ht.sum())
            dev = np.abs(t.to_numpy().sum(axis=1) - flows)
            parts = [dev[:nd], dev[nd:nd + nl], dev[nd + nl:nd + 2 * nl], dev[nd + 2 * nl:]]
            mx = [float(x.max()) if x.size else None for x in parts]
            fails += (not lab) + int((d > TOL).sum()) + sum(1 for m in mx if m is not None and m > TOL)
            out.append(f"| {kind} | {lab} | {g(float(d.max()))} | {int((d > TOL).sum())} | "
                       + " | ".join("" if m is None else g(m) for m in mx) + f" | {len(t)} (Case-3 rows {n3}) |")
        # (c) the 2x2 stamp
        bus = data["bus"].set_index("BUS_I"); br = data["branch"]
        row = br[(br.F_BUS == F_) & (br.T_BUS == T_)].iloc[0]
        k = int(np.where(((br.F_BUS == F_) & (br.T_BUS == T_)).to_numpy())[0][0])
        Vf = bus.loc[F_, "VM"] * np.exp(1j * np.deg2rad(bus.loc[F_, "VA"]))
        Vt = bus.loc[T_, "VM"] * np.exp(1j * np.deg2rad(bus.loc[T_, "VA"]))
        tap = b["TAP"] if b["TAP"] != 0 else 1.0
        tau = tap * np.exp(1j * np.deg2rad(b["SHIFT"]))
        ys = 1 / (b["R"] + 1j * b["X"])
        Ytt = ys + 1j * b["B"] / 2; Yff = Ytt / abs(tau) ** 2; Yft = -ys / np.conj(tau); Ytf = -ys / tau
        If = Yff * Vf + Yft * Vt; It = Ytf * Vf + Ytt * Vt
        Sf = Vf * np.conj(If) * base; St = Vt * np.conj(It) * base
        Vi = Vf / tau; Ii = np.conj(tau) * If
        Is = ys * (Vi - Vt)
        qsf, qst = branch_shunt_flows(case)
        Qs_port = float(row.QF) + qsf[k]; Qe_port = float(row.QT) - qst[k]
        chk = [
            ("P^f: exported PF against Re S_f", float(row.PF), Sf.real),
            ("Q^f: exported QF against Im S_f", float(row.QF), Sf.imag),
            ("P^t (delivered): exported PT against -Re S_t", float(row.PT), -St.real),
            ("Q^t (delivered): exported QT against -Im S_t", float(row.QT), -St.imag),
            ("ideal transformer: V_f conj(I_f) against V_i conj(I_i), real part", (Vf * np.conj(If)).real * base, (Vi * np.conj(Ii)).real * base),
            ("ideal transformer: V_f conj(I_f) against V_i conj(I_i), imaginary part", (Vf * np.conj(If)).imag * base, (Vi * np.conj(Ii)).imag * base),
            ("from-end shunt: port Qsf against (B/2)|V_f|^2/|tau|^2 baseMVA (= -Q^sf of the brief)", qsf[k], b["B"] / 2 * abs(Vf) ** 2 / abs(tau) ** 2 * base),
            ("to-end shunt: port Qst against (B/2)|V_t|^2 baseMVA", qst[k], b["B"] / 2 * abs(Vt) ** 2 * base),
            ("Q^s = Q^f + Qsf(port) = Q^f - Q^sf(brief) against Im(V_i conj(I_s)) baseMVA", Qs_port, (Vi * np.conj(Is)).imag * base),
            ("Q^e = Q^t - Qst(port) against Im(V_t conj(I_s)) baseMVA", Qe_port, (Vt * np.conj(Is)).imag * base),
        ]
        out += ["", f"(c) Branch {F_}-{T_} from the 2x2 stamp and the exported voltages (V_f = {g(abs(Vf))} at {g(float(bus.loc[F_, 'VA']))} degrees, "
                f"V_t = {g(abs(Vt))} at {g(float(bus.loc[T_, 'VA']))} degrees; |tau| = {g(abs(tau))}):", "",
                "| quantity | value used by the method | recomputed | abs difference |", "|---|---|---|---|"]
        for q, a_, b_ in chk:
            dd = abs(a_ - b_); fails += dd > TOL
            out.append(f"| {q} | {g(a_)} | {g(b_)} | {g(dd)} |")
        # (d) the old method
        out += ["", "(d) The old method (in memory, corrected QSF, cap 10000) against the port:", "",
                "| kind | old method | iterations | max abs cell difference | unmapped port rows | port row sum - old flow | old row error |",
                "|---|---|---|---|---|---|---|"]
        branch_old = C.correct_qsf(data)[0]
        for kind in ("apfc", "rpfc"):
            tt, _, _ = K.port_table(data, kind)
            try:
                fq, A, done, it = T.old_with_iterations(data, kind, 10000, branch_old, 600)
                if done:
                    c = K.compare(data, kind, A, fq, tt)
                    fails += (c["max_cell"] > TOL) + bool(c["unmapped_port_rows"]) + (c["port_rowsum_vs_old_flow"] > TOL)
                    out.append(f"| {kind} | finished | {it} | {g(c['max_cell'])} | {c['unmapped_port_rows']} | {g(c['port_rowsum_vs_old_flow'])} | {g(c['old_row_error'])} |")
                else:
                    out.append(f"| {kind} | did not finish (iteration cap) | {it} | | | | |")
            except Exception as e:
                out.append(f"| {kind} | did not finish (error: {T.error_text(e)}) | | | | | |")
        out += ["", f"Checks above 1e-10 or failing in {name}: {fails}", ""]
    print("\n".join(out))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
