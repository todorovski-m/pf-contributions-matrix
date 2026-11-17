function [RPFC_tbl, mpc, info] = rpfc(casedata, opts)
%RPFC  Reactive Power Flow Components (RPFC) for a MATPOWER case.
%
%   [RPFC_TBL, MPC, DIAG] = RPFC(CASEDATA, OPTS) runs a power flow for
%   CASEDATA and computes the reactive power flow contributions from each
%   reactive source to every branch end and demand in the network.
%
%   CASEDATA is a case name (for example 'cases/case9'), an unsolved MATPOWER
%   case struct, or, with OPTS.solved = true, a struct that already holds a
%   solved power flow (for example the output of RUNPF); then RUNPF is not
%   called. OPTS is optional; see PFC_CHECKS('options') for the fields
%   (tol, tol_div, write, solved, name). The defaults reproduce the earlier
%   behaviour: RPFC(CASEDATA) runs RUNPF and writes results/<case>_rpfc.csv.
%
%   DIAG is a struct with sizes, the number of non-zeros of N and A, whether
%   N, A and B are sparse, an estimate of the 1-norm condition number of A
%   (CONDEST), the range of the attenuation factors, and the numbers of
%   removed, flipped and converted elements.
%
%   Checks (see PFC_CHECKS) fall into three groups:
%     errors, for conditions under which the results are not valid: the
%       reactive flow digraph is cyclic and the spectral radius of N^Q is not
%       below 1 (pfc:rpfc:spectralRadius); a bus with an outgoing branch has
%       (near) zero total inflow (pfc:rpfc:zeroInflow); a to-end shunt
%       self-term is below -tol_div (pfc:rpfc:negativeSelfTerm; cannot occur
%       for a valid carrier, see the comment at the self-terms);
%     warnings, for conditions that can invalidate the results; they are
%       listed in DIAG.warnings (rpfc currently issues none; APFC warns for a
%       cyclic active digraph and for an active attenuation above 1; a
%       spectral radius close to 1 and a near-singular A^Q are not yet
%       checked, see DIAG.spectral_radius and DIAG.condest_A);
%     information only, recorded in DIAG and never printed: a cyclic
%       reactive digraph with spectral radius below 1 (DIAG.acyclic,
%       DIAG.spectral_radius; circulating flows are allowed), and branches
%       with a negative series flow, that is a negative reactive attenuation
%       (DIAG.negative_series_flow.count and .branches). Negative attenuation
%       occurs in standard cases (case14, case118) and leads to negative
%       contributions in some cells; it does not invalidate the results.
%
%   The method uses a matrix-based formulation based on the proportionality
%   principle for reactive flow-component allocation.
%
%   Three branch cases for reactive power:
%     Case 1 (normal carrier): positive flow at both ends.
%     Case 2 (consumes from both ends): positive at from-end and negative
%       at to-end, handled as two loads.
%     Case 3 (injects at both ends): negative at from-end and positive at
%       to-end, treated as pure reactive sources.
%   Branches with exactly one zero terminal (sign s(x) = + if x > tol, - if
%   x < -tol, 0 otherwise, on QF and the delivered QT, before the flip; a
%   dead-end branch that carries only its own charging or shunt power) are
%   members of Case 2 or Case 3 with the zero terminal omitted:
%     (+,0) one demand QF at the from bus, (0,-) one demand -QT at the to bus
%       (rows C2_b, with the Case-2 shunt rules);
%     (-,0) one source -QF at the from bus (column E_i_j_f, row L_{f.i-j}),
%       (0,+) one source QT at the to bus (column E_i_j_t, row L_{t.i-j}),
%       with the Case-3 shunt rules.
%   Such a branch is never flipped and is not a carrier. DIAG.converted
%   counts them in zero_terminal_demands and zero_terminal_sources; the
%   counts case2_branches_to_demands and case3_branches_to_sources are those
%   of the branches with both terminals non-zero.
%
%   Attenuation uses branch series flows after removing local shunt effects
%   from the terminal measurements.
%
%   Shunt self-contributions at the to-end of each branch are computed in
%   post-processing and propagated downstream by the same proportionality.
%
%   Conversions of negative injections (as in the iterative method):
%     a generator with negative reactive output (QG < -tol) is a reactive
%       demand (row G_b);
%     a bus with negative reactive demand (QD < -tol) is a reactive source
%       (column D_b);
%     a bus shunt is a source if it injects (column S_b) and a demand if it
%       absorbs (row Bs_b);
%     an inductive transformer shunt (negative shunt flow at the from-end or
%       the to-end) is a demand (row Tf_b or Tt_b), except for a Case-2
%       branch: its two terminal demands (rows C2_b) already contain its
%       shunt power, so it gets no Tf_b or Tt_b row (as in the iterative
%       method, which converts Case-2 branches before it converts the
%       transformer shunts). Case-3 branches keep their Tf_b and Tt_b rows,
%       as in the iterative method.
%
%   Source columns, in this order: G_b (generators with QG > tol), S_b
%   (injecting bus shunts), D_b (buses with negative reactive demand),
%   E_i_j_f and E_i_j_t (the two end injections of each Case-3 branch,
%   interleaved per branch), then E_i_j_s (branch shunt self-contributions).
%   Labels are unique: a repeated label gets the suffix _2, _3, ... (for
%   example several generators at bus 22: G_22, G_22_2; two Case-2 branch
%   ends at bus 49: C2_49, C2_49_2).
%
%   Output RPFC_TBL is a table with contribution values and branch labels.
%   MPC is the modified MATPOWER case after power flow.
%

%% Options, load case and run MATPOWER power flow
if nargin < 1 || isempty(casedata)
    casedata = 'cases/case1';
end
if nargin < 2
    opts = struct();
end
opts = pfc_checks('options', opts);
tol = opts.tol;          % deletion and classification threshold
tol_div = opts.tol_div;  % division guard and negligible self-term threshold

define_constants;
[mpc, name] = load_solved(casedata, opts);

%% Normalize branch orientations
%  MATPOWER stores PT/QT as power from the to-bus into the branch; sign-correct to power leaving the branch into the to-bus.
mpc.branch(:, [PT QT]) = -mpc.branch(:, [PT QT]);

%% Calculate reactive power flows in branches' shunt elements
tap = mpc.branch(:, TAP);  % tap = 0 in MATPOWER means transmission line
tap(tap == 0) = 1; % set tap = 1 for transmission lines
V = mpc.bus(:, VM).*exp(1j*mpc.bus(:,VA)/180*pi); % complex bus voltages
fb = mpc.branch(:, F_BUS); % vector of "from" buses
tb = mpc.branch(:, T_BUS); % vector of "to" buses
Yb = 1./(mpc.branch(:, BR_R) + 1j * mpc.branch(:, BR_X)); % branch series admittance

I   = Yb./tap .* (V(fb) - V(tb));
Qs  = imag(V(fb).*conj(I)) * mpc.baseMVA;
Qe  = imag(V(tb).*conj(I)) * mpc.baseMVA;
Qsf = Qs - mpc.branch(:, QF);
Qst = mpc.branch(:, QT) - Qe;

% Branches that consume reactive power at both ends (Case 2, QF > 0 and
% delivered QT < 0) are later replaced by their two terminal demands (C2_b),
% which already contain their shunt power. As in the iterative method, they
% get no transformer-shunt demand (Tf_b, Tt_b). The test on the delivered
% terminal flows here gives the same Case-2 set as the classification below:
% the orientation flip never applies to a Case-2 branch.
% Branches with exactly one zero terminal (sign by tol on QF and the
% delivered QT, before the flip): zt = 1 (+,0) and zt = 2 (0,-) are Case-2
% members with the zero terminal omitted (one terminal demand); zt = 3 (-,0)
% and zt = 4 (0,+) are Case-3 members with the zero terminal omitted (one
% terminal source). zt = 0 for every other branch, which keeps the exact-sign
% tests below. Such a branch is never flipped.
sgn = @(x) (x > tol) - (x < -tol);
s_f = sgn(mpc.branch(:, QF));
s_t = sgn(mpc.branch(:, QT));
zt = zeros(size(mpc.branch, 1), 1);
zt(s_f > 0 & s_t == 0) = 1;
zt(s_f == 0 & s_t < 0) = 2;
zt(s_f < 0 & s_t == 0) = 3;
zt(s_f == 0 & s_t > 0) = 4;
is_case2_branch = ((mpc.branch(:, QF) > 0) & (mpc.branch(:, QT) < 0) & zt == 0) | zt == 1 | zt == 2;

% Handle negative Qsf (inductive shunt at from-end): convert to extra demand
i_neg_qsf       = Qsf < -tol & ~is_case2_branch;
extra_d_bus_qsf = mpc.branch(i_neg_qsf, F_BUS);
extra_Qd_qsf    = -Qsf(i_neg_qsf);
Qsf(i_neg_qsf) = 0;

% Handle negative Qst (inductive shunt at to-end): convert to extra demand
i_neg_qst       = Qst < -tol & ~is_case2_branch;
extra_d_bus_qst = mpc.branch(i_neg_qst, T_BUS);
extra_Qd_qst    = -Qst(i_neg_qst);
Qst(i_neg_qst) = 0;

%% Remove inactive branches and zero-power generators
i_zero = abs(mpc.branch(:, QF)) < tol & abs(mpc.branch(:, QT)) < tol;
mpc.branch = mpc.branch(~i_zero, :);
Qsf = Qsf(~i_zero);
Qst = Qst(~i_zero);
zt  = zt(~i_zero);

i_zero_gen = abs(mpc.gen(:, QG)) < tol;
mpc.gen = mpc.gen(~i_zero_gen, :);

%% Flip branches where both QF < 0 and QT < 0
i_flip = (mpc.branch(:, QF) < 0) & (mpc.branch(:, QT) < 0) & zt == 0;
if any(i_flip)
    mpc.branch(i_flip, [F_BUS T_BUS]) = mpc.branch(i_flip, [T_BUS F_BUS]);
    tmp = mpc.branch(i_flip, QF);
    mpc.branch(i_flip, QF) = -mpc.branch(i_flip, QT);
    mpc.branch(i_flip, QT) = -tmp;
    tmp = Qsf(i_flip);
    Qsf(i_flip) = Qst(i_flip);
    Qst(i_flip) = tmp;
end

%% Case 3: injects at both ends (pure source)
%  Remove from branch set; register as source columns. Members with one zero
%  terminal (zt = 3, 4) have only the source of their non-zero terminal: no
%  column and no appended row for the zero terminal.
i_case3    = ((mpc.branch(:, QF) < 0) & (mpc.branch(:, QT) > 0) & zt == 0) | zt == 3 | zt == 4;
case3_fbus = mpc.branch(i_case3, F_BUS);
case3_tbus = mpc.branch(i_case3, T_BUS);
case3_QF   = -mpc.branch(i_case3, QF);   % injection into f-bus (positive)
case3_QT   =  mpc.branch(i_case3, QT);   % injection into t-bus (positive)
case3_hasf = zt(i_case3) ~= 4;           % from-end source present
case3_hast = zt(i_case3) ~= 3;           % to-end source present
case3_flab = arrayfun(@(a,b) sprintf('L_{f.%d-%d}',a,b), ...
    case3_fbus(case3_hasf), case3_tbus(case3_hasf), 'UniformOutput', false);
case3_tlab = arrayfun(@(a,b) sprintf('L_{t.%d-%d}',a,b), ...
    case3_fbus(case3_hast), case3_tbus(case3_hast), 'UniformOutput', false);
n_case3    = sum(i_case3);
n_case3_zt = nnz(zt(i_case3) ~= 0);

mpc.branch = mpc.branch(~i_case3, :);
Qsf = Qsf(~i_case3);
Qst = Qst(~i_case3);
zt  = zt(~i_case3);

%% Case 2: consumes from both ends
%  Split into two load branches at f-bus and t-bus: first the from-terminals
%  of all members, then their to-terminals. Members with one zero terminal
%  (zt = 1, 2) have only the demand of their non-zero terminal.
i_case2 = ((mpc.branch(:, QF) > 0) & (mpc.branch(:, QT) < 0) & zt == 0) | zt == 1 | zt == 2;
c2_hasf = i_case2 & zt ~= 2;
c2_hast = i_case2 & zt ~= 1;
c2_extra_bus = [mpc.branch(c2_hasf, F_BUS); mpc.branch(c2_hast, T_BUS)];
c2_extra_Q   = [mpc.branch(c2_hasf, QF);   -mpc.branch(c2_hast, QT)   ];
n_case2_zt   = nnz(i_case2 & zt ~= 0);
n_case2_both = nnz(i_case2 & zt == 0);
mpc.branch = mpc.branch(~i_case2, :);
Qsf = Qsf(~i_case2);
Qst = Qst(~i_case2);

nl   = size(mpc.branch, 1);
nb   = size(mpc.bus, 1);
sink = nb + 1;
%% Reactive demand branches (all to sink)
% QD loads
i_qd     = mpc.bus(:, QD) > tol;
d_bus_qd = mpc.bus(i_qd, BUS_I);
Qd_load  = mpc.bus(i_qd, QD);

% Negative-QG generators -> reactive loads
i_qg_neg    = mpc.gen(:, QG) < -tol;
d_bus_qgneg = mpc.gen(i_qg_neg, GEN_BUS);
Qd_qgneg    = -mpc.gen(i_qg_neg, QG);

% Negative bus shunts -> reactive loads
Qsh = mpc.bus(:, VM).^2 .* mpc.bus(:, BS);
i_qs_neg    = Qsh < -tol;
d_bus_qsneg = mpc.bus(i_qs_neg, BUS_I);
Qd_qsneg    = -Qsh(i_qs_neg);

% Combine all demands (including extra demands from negative transformer shunts)
d_bus = [d_bus_qd; d_bus_qgneg; d_bus_qsneg; c2_extra_bus; extra_d_bus_qsf; extra_d_bus_qst];
Qd    = [Qd_load;  Qd_qgneg;    Qd_qsneg;    c2_extra_Q;   extra_Qd_qsf;    extra_Qd_qst   ];
nd    = numel(Qd);

%% Reactive sources
% Positive-QG generators
i_qg_pos  = mpc.gen(:, QG) > tol;
src_bus_g = mpc.gen(i_qg_pos, GEN_BUS);
Qg_src    = mpc.gen(i_qg_pos, QG);
ng        = numel(Qg_src);

% Positive bus shunts
i_qs_pos   = Qsh > tol;
src_bus_qs = mpc.bus(i_qs_pos, BUS_I);
Qs_src     = Qsh(i_qs_pos);
n_qs       = numel(Qs_src);

% Buses with negative reactive demand -> reactive sources (column D_b)
i_qd_neg   = mpc.bus(:, QD) < -tol;
src_bus_qd = mpc.bus(i_qd_neg, BUS_I);
Qd_src     = -mpc.bus(i_qd_neg, QD);
n_qd       = numel(Qd_src);

% Case 3 end-injections, interleaved per branch: [f_1; t_1; f_2; t_2; ...],
% the same order as the column names (E_i_j_f, E_i_j_t) and the appended rows
% (the terminal of a member with one zero terminal is omitted)
case3_bus = reshape([case3_fbus, case3_tbus].', [], 1);
case3_Q   = reshape([case3_QF,   case3_QT  ].', [], 1);
case3_has = reshape([case3_hasf, case3_hast].', [], 1);
case3_bus = case3_bus(case3_has);
case3_Q   = case3_Q(case3_has);
src_bus = [src_bus_g;       src_bus_qs;  src_bus_qd;  case3_bus];
Qsrc    = [Qg_src;          Qs_src;      Qd_src;      case3_Q  ];
ns      = numel(Qsrc);
n_pre   = ng + n_qs + n_qd;   % number of source columns before the Case-3 columns

%% Unified branch set
ne = nl + nd;
f  = [mpc.branch(:, F_BUS); d_bus            ];
t  = [mpc.branch(:, T_BUS); repmat(sink,nd,1)];
Qf  = [mpc.branch(:, QF);  Qd          ];
Qt  = [mpc.branch(:, QT);  Qd          ];
Qsf = [Qsf;                zeros(nd,1) ];
Qst = [Qst;                zeros(nd,1) ];

%% Series flows and attenuation factors
Qs  = Qf + Qsf;   % series entry flow (from-end shunt removed)
Qe  = Qt - Qst;   % series exit flow  (to-end shunt removed)
aQ  = Qe ./ Qs;   % attenuation factor

%% Incidence matrices and total inflow per bus
Cf   = sparse(1:ne, f,       1, ne, nb+1);
Ct   = sparse(1:ne, t,       1, ne, nb+1);
Cs   = sparse(1:ns, src_bus, 1, ns, nb+1);

Qtot = Ct' * Qt + Cs' * Qsrc;
h    = 1 ./ Qtot;
hf   = h(f);

%% Validity checks (see pfc_checks.m)
%  branch labels: index in the returned mpc.branch and the oriented bus pair
br_labels = arrayfun(@(k) sprintf('%d (%d-%d)', k, f(k), t(k)), (1:nl)', 'UniformOutput', false);
pfc_checks('inflow', Qtot, f, tol_div, 'rpfc');
[~, i_ser] = pfc_checks('series_flow', Qs(1:nl), Qe(1:nl), br_labels, tol);   % informational, no warning

%% Build the system matrix and right-hand side, then solve.
M   = Cf * Ct';
N   = -spdiags(Qf .* hf, 0, ne, ne) * (M * spdiags(aQ, 0, ne, ne));
A   = speye(ne) + N;

% acyclicity of the reactive flow digraph and spectral radius of N^Q
[acyclic, rho, w_cyc] = pfc_checks('cycle', N, M, 'rpfc');
S   = Cf * Cs';
B   = spdiags(Qf .* hf, 0, ne, ne) * S * spdiags(Qsrc, 0, ns, ns);
RPFCf = A \ B;
RPFCt = spdiags(aQ, 0, ne, ne) * RPFCf;

%% Post-processing: shunt self-contributions.
%  Compute branch to-end self terms and propagate them with the same
%  system matrix used above.
self_toend = zeros(nl, 1);
for i = 1:nl
    if abs(Qs(i)) > tol_div
        self_toend(i) = (Qe(i) / Qs(i)) * Qsf(i) + Qst(i);
    end
end
% Guard: a negative self-term would not be propagated (only terms above
% tol_div are), so the downstream row sums would not equal the flows.
% With Qs = Qf + Qsf and Qe = Qt - Qst the self-term equals
%     e^st = (Qe/Qs)*Qsf + Qst = Qt - aQ*Qf.
% A valid carrier has Qf > 0 and Qt > 0 (Case 2 and Case 3 removed, both-
% negative branches flipped) and Qsf, Qst >= 0 after the conversion of
% inductive shunts to demands. Then e^st > Qt > 0 if aQ < 0, and
% e^st = aQ*Qsf + Qst >= 0 if aQ >= 0, so the error cannot occur; it can only
% be reached through shunt flows in [-tol, 0), which are not converted.
pfc_checks('selfterm', self_toend, br_labels, tol_div);

i_self    = find(self_toend > tol_div);
n_self    = numel(i_self);
RPFCf_self  = zeros(ne, nl);
RPFCt_self = zeros(ne, nl);

if n_self > 0
    Css      = sparse(1:n_self, t(i_self), 1, n_self, nb+1);
    Ss       = Cf * Css';
    B_self = spdiags(Qf .* hf, 0, ne, ne) * Ss * spdiags(self_toend(i_self), 0, n_self, n_self);
    sol      = A \ B_self;
    for k = 1:n_self
        ii = i_self(k);
        RPFCf_self(:, ii)  = sol(:, k);
        RPFCt_self(:, ii) = aQ .* sol(:, k);
    end
end
% Add direct shunt self-contribution at each branch's own to-end
for i = 1:nl
    RPFCt_self(i,i) = RPFCt_self(i,i) + self_toend(i);
end

%% Assemble output table
% Column names
src_names = cell(1, ns);
for i = 1:ng
    src_names{i} = sprintf('G_%d', src_bus_g(i));
end
for i = 1:n_qs
    src_names{ng+i} = sprintf('S_%d', src_bus_qs(i));
end
for i = 1:n_qd
    src_names{ng+n_qs+i} = sprintf('D_%d', src_bus_qd(i));
end
% Case-3 columns: per member E_i_j_f (if present), then E_i_j_t (if present);
% case3_colf(i), case3_colt(i) are the column indices (0 if absent)
case3_colf = zeros(n_case3, 1);
case3_colt = zeros(n_case3, 1);
c = n_pre;
for i = 1:n_case3
    if case3_hasf(i)
        c = c + 1; case3_colf(i) = c;
        src_names{c} = sprintf('E_%d_%d_f', case3_fbus(i), case3_tbus(i));
    end
    if case3_hast(i)
        c = c + 1; case3_colt(i) = c;
        src_names{c} = sprintf('E_%d_%d_t', case3_fbus(i), case3_tbus(i));
    end
end
self_names = arrayfun(@(i) sprintf('E_%d_%d_s', f(i), t(i)), ...
    (1:nl)', 'UniformOutput', false);

% Disambiguate duplicate self-contribution names (parallel branches)
seen = containers.Map('KeyType','char','ValueType','int32');
for i = 1:nl
    nm = self_names{i};
    if isKey(seen, nm)
        seen(nm) = seen(nm) + 1;
        self_names{i} = sprintf('%s_%d', nm, seen(nm));
    else
        seen(nm) = 1;
    end
end

all_names  = [src_names, self_names'];

% Full matrices
RPFCf_full = [RPFCf, RPFCf_self];
RPFCt_full = [RPFCt, RPFCt_self];

% Row labels — each demand group gets a distinct prefix:
%   D_b  : reactive load (QD > 0) at bus b
%   G_b  : absorbing generator (QG < 0) at bus b
%   Bs_b : inductive bus shunt (BS < 0) at bus b
%   C2_b : Case 2 branch terminal at bus b
%   Tf_b : inductive transformer shunt, from-end at bus b
%   Tt_b : inductive transformer shunt, to-end at bus b
d_labels = [
    arrayfun(@(b) sprintf('D_%d',  b), d_bus_qd,        'UniformOutput', false);
    arrayfun(@(b) sprintf('G_%d',  b), d_bus_qgneg,     'UniformOutput', false);
    arrayfun(@(b) sprintf('Bs_%d', b), d_bus_qsneg,     'UniformOutput', false);
    arrayfun(@(b) sprintf('C2_%d', b), c2_extra_bus,    'UniformOutput', false);
    arrayfun(@(b) sprintf('Tf_%d', b), extra_d_bus_qsf, 'UniformOutput', false);
    arrayfun(@(b) sprintf('Tt_%d', b), extra_d_bus_qst, 'UniformOutput', false);
    ];
lf_labels = arrayfun(@(i) sprintf('L_{f.%d-%d}',f(i),t(i)), (1:nl)', 'UniformOutput', false);
lt_labels = arrayfun(@(i) sprintf('L_{t.%d-%d}',f(i),t(i)), (1:nl)', 'UniformOutput', false);

% Disambiguate duplicate row labels for parallel branches
seen_rows = containers.Map('KeyType','char','ValueType','int32');
for i = 1:nl
    nm = lf_labels{i};
    if isKey(seen_rows, nm)
        seen_rows(nm) = seen_rows(nm) + 1;
        lf_labels{i} = sprintf('%s_%d', nm, seen_rows(nm));
        lt_labels{i} = sprintf('%s_%d', lt_labels{i}, seen_rows(nm));
    else
        seen_rows(nm) = 1;
    end
end

% Assemble rows: demands first, then branch from-ends and to-ends
RPFC_data    = [RPFCf_full(nl+1:end,:); RPFCf_full(1:nl,:); RPFCt_full(1:nl,:)];
% Demand rows: two conversions at the same bus (C2_b, Tf_b, Tt_b, G_b, ...)
% get the suffix _2, _3, ... (the scheme of the parallel-branch rows)
d_labels = unique_labels(d_labels);
row_labels = [d_labels; lf_labels; lt_labels];

% Append Case 3 branch rows (each end sourced entirely by that branch's injection)
n_tot = ns + nl;
if n_case3 > 0
    % from-end rows of the members with a from-end source, then the to-end rows
    i_f = find(case3_hasf);
    i_t = find(case3_hast);
    c3f_data = zeros(numel(i_f), n_tot);
    c3t_data = zeros(numel(i_t), n_tot);
    for k = 1:numel(i_f)
        c3f_data(k, case3_colf(i_f(k))) = -case3_QF(i_f(k));
    end
    for k = 1:numel(i_t)
        c3t_data(k, case3_colt(i_t(k))) =  case3_QT(i_t(k));
    end
    RPFC_data    = [RPFC_data;    c3f_data; c3t_data   ];
    row_labels = [row_labels; case3_flab; case3_tlab];
end

%% Diagnostics
info = struct();
info.case = name;
info.n_bus = nb;
info.n_branch = nl;
info.n_demand = nd;
info.n_source = ns;
info.n_e = ne;
info.nnz_N = nnz(N);
info.nnz_A = nnz(A);
info.issparse = struct('N', issparse(N), 'A', issparse(A), 'B', issparse(B));
info.condest_A = condest(A);
if nl > 0
    info.attenuation_range = [min(aQ(1:nl)), max(aQ(1:nl))];
else
    info.attenuation_range = zeros(1, 0);
end
info.removed = struct('branches', nnz(i_zero), 'generators', nnz(i_zero_gen));
info.flipped_branches = nnz(i_flip);
info.converted = struct( ...
    'negative_generation_to_demand', numel(d_bus_qgneg), ...
    'negative_demand_to_source', n_qd, ...
    'inductive_bus_shunt_to_demand', numel(d_bus_qsneg), ...
    'case2_branches_to_demands', n_case2_both, ...
    'case3_branches_to_sources', n_case3 - n_case3_zt, ...
    'inductive_from_end_shunt_to_demand', numel(extra_d_bus_qsf), ...
    'inductive_to_end_shunt_to_demand', numel(extra_d_bus_qst), ...
    'zero_terminal_demands', n_case2_zt, ...
    'zero_terminal_sources', n_case3_zt);
info.acyclic = acyclic;
info.spectral_radius = rho;
info.negative_series_flow = struct('count', numel(i_ser), 'branches', {br_labels(i_ser).'});
info.warnings = w_cyc;
info.opts = opts;

%% Output table
% Unique names: several generators at one bus give G_b, G_b_2, ...; the final
% pass on all names and all rows leaves already unique labels unchanged.
all_names = unique_labels(all_names);
row_labels = unique_labels(row_labels);
RPFC_tbl = array2table(full(RPFC_data), 'VariableNames', all_names);   % full only for the table
RPFC_tbl = addvars(RPFC_tbl, row_labels, 'Before', 1, 'NewVariableNames', 'Line');

if opts.write && ~isempty(name)
    result = strrep(name, 'cases/', 'results/');
    writetable(RPFC_tbl, [result '_rpfc.csv']);
end
end

function [mpc, name] = load_solved(casedata, opts)
%LOAD_SOLVED  Load CASEDATA and run the power flow unless OPTS.solved.
define_constants;
if ischar(casedata)
    name = casedata;
    mpc = loadcase(casedata);
elseif isstruct(casedata)
    name = opts.name;
    if opts.solved
        mpc = casedata;          % keep all fields of the solved struct
    else
        mpc = loadcase(casedata);
    end
else
    error('pfc:input:casedata', 'CASEDATA must be a case name or a MATPOWER case struct.');
end
if opts.solved
    if size(mpc.branch, 2) < QT
        error('pfc:input:notSolved', ...
            'opts.solved is true, but the branch matrix has no power flow columns (PF, QF, PT, QT).');
    end
    if isfield(mpc, 'success') && ~mpc.success
        error('pfc:input:notConverged', 'opts.solved is true, but the power flow did not converge (success = 0).');
    end
else
    mpc = runpf(mpc, mpoption('verbose', 0, 'out.all', 0));
    if ~mpc.success
        error('pfc:runpf:notConverged', 'The MATPOWER power flow did not converge.');
    end
end
end

function lab = unique_labels(lab)
%UNIQUE_LABELS  Make labels unique: the second and later occurrences of a
%   label get the suffix _2, _3, ... (the scheme of the parallel-branch rows).
%   A suffix is skipped if the resulting label is already in use.
seen = containers.Map('KeyType', 'char', 'ValueType', 'int32');
used = containers.Map('KeyType', 'char', 'ValueType', 'logical');
for i = 1:numel(lab)
    used(lab{i}) = true;
end
for i = 1:numel(lab)
    nm = lab{i};
    if isKey(seen, nm)
        k = seen(nm) + 1;
        while isKey(used, sprintf('%s_%d', nm, k))
            k = k + 1;
        end
        seen(nm) = k;
        lab{i} = sprintf('%s_%d', nm, k);
        used(lab{i}) = true;
    else
        seen(nm) = 1;
    end
end
end
