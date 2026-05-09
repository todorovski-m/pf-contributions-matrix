function [RPFC_tbl, mpc] = rpfc(casedata)
%RPFC  Reactive Power Flow Components (RPFC) for a MATPOWER case.
%
%   [RPFC_TBL, MPC] = RPFC(CASEDATA) runs a power flow for CASEDATA and computes
%   the reactive power flow contributions from each reactive source to every
%   branch end and demand in the network.
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
%
%   Attenuation uses branch series flows after removing local shunt effects
%   from the terminal measurements.
%
%   Shunt self-contributions at the to-end of each branch are computed in
%   post-processing and propagated downstream by the same proportionality.
%
%   Output RPFC_TBL is a table with contribution values and branch labels.
%   MPC is the modified MATPOWER case after power flow.
%

%% Load case and run MATPOWER power flow
if nargin < 1
    casedata = 'cases/case1';
end

define_constants;
opt = mpoption('verbose', 0, 'out.all', 0);
mpc = loadcase(casedata);
mpc = runpf(mpc, opt);

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

% Handle negative Qsf (inductive shunt at from-end): convert to extra demand
i_neg_qsf       = Qsf < -1e-6;
extra_d_bus_qsf = mpc.branch(i_neg_qsf, F_BUS);
extra_Qd_qsf    = -Qsf(i_neg_qsf);
Qsf(i_neg_qsf) = 0;

% Handle negative Qst (inductive shunt at to-end): convert to extra demand
i_neg_qst       = Qst < -1e-6;
extra_d_bus_qst = mpc.branch(i_neg_qst, T_BUS);
extra_Qd_qst    = -Qst(i_neg_qst);
Qst(i_neg_qst) = 0;

%% Remove inactive branches and zero-power generators
i_zero = abs(mpc.branch(:, QF)) < 1e-6 & abs(mpc.branch(:, QT)) < 1e-6;
mpc.branch = mpc.branch(~i_zero, :);
Qsf = Qsf(~i_zero);
Qst = Qst(~i_zero);

i_zero_gen = abs(mpc.gen(:, QG)) < 1e-6;
mpc.gen = mpc.gen(~i_zero_gen, :);

%% Flip branches where both QF < 0 and QT < 0
i_flip = (mpc.branch(:, QF) < 0) & (mpc.branch(:, QT) < 0);
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
%  Remove from branch set; register as source columns.
i_case3    = (mpc.branch(:, QF) < 0) & (mpc.branch(:, QT) > 0);
case3_fbus = mpc.branch(i_case3, F_BUS);
case3_tbus = mpc.branch(i_case3, T_BUS);
case3_QF   = -mpc.branch(i_case3, QF);   % injection into f-bus (positive)
case3_QT   =  mpc.branch(i_case3, QT);   % injection into t-bus (positive)
case3_flab = arrayfun(@(a,b) sprintf('L_{f.%d-%d}',a,b), ...
    case3_fbus, case3_tbus, 'UniformOutput', false);
case3_tlab = arrayfun(@(a,b) sprintf('L_{t.%d-%d}',a,b), ...
    case3_fbus, case3_tbus, 'UniformOutput', false);
n_case3    = sum(i_case3);

mpc.branch = mpc.branch(~i_case3, :);
Qsf = Qsf(~i_case3);
Qst = Qst(~i_case3);

%% Case 2: consumes from both ends
%  Split into two load branches at f-bus and t-bus.
i_case2 = (mpc.branch(:, QF) > 0) & (mpc.branch(:, QT) < 0);
if any(i_case2)
    c2_extra_bus = [mpc.branch(i_case2, F_BUS); mpc.branch(i_case2, T_BUS)];
    c2_extra_Q   = [mpc.branch(i_case2, QF);   -mpc.branch(i_case2, QT)   ];
    mpc.branch = mpc.branch(~i_case2, :);
    Qsf = Qsf(~i_case2);
    Qst = Qst(~i_case2);
else
    c2_extra_bus = [];
    c2_extra_Q   = [];
end

nl   = size(mpc.branch, 1);
nb   = size(mpc.bus, 1);
sink = nb + 1;
%% Reactive demand branches (all to sink)
% QD loads
i_qd     = mpc.bus(:, QD) > 1e-6;
d_bus_qd = mpc.bus(i_qd, BUS_I);
Qd_load  = mpc.bus(i_qd, QD);

% Negative-QG generators -> reactive loads
i_qg_neg    = mpc.gen(:, QG) < -1e-6;
d_bus_qgneg = mpc.gen(i_qg_neg, GEN_BUS);
Qd_qgneg    = -mpc.gen(i_qg_neg, QG);

% Negative bus shunts -> reactive loads
Qsh = mpc.bus(:, VM).^2 .* mpc.bus(:, BS);
i_qs_neg    = Qsh < -1e-6;
d_bus_qsneg = mpc.bus(i_qs_neg, BUS_I);
Qd_qsneg    = -Qsh(i_qs_neg);

% Combine all demands (including extra demands from negative transformer shunts)
d_bus = [d_bus_qd; d_bus_qgneg; d_bus_qsneg; c2_extra_bus; extra_d_bus_qsf; extra_d_bus_qst];
Qd    = [Qd_load;  Qd_qgneg;    Qd_qsneg;    c2_extra_Q;   extra_Qd_qsf;    extra_Qd_qst   ];
nd    = numel(Qd);

%% Reactive sources
% Positive-QG generators
i_qg_pos  = mpc.gen(:, QG) > 1e-6;
src_bus_g = mpc.gen(i_qg_pos, GEN_BUS);
Qg_src    = mpc.gen(i_qg_pos, QG);
ng        = numel(Qg_src);

% Positive bus shunts
i_qs_pos   = Qsh > 1e-6;
src_bus_qs = mpc.bus(i_qs_pos, BUS_I);
Qs_src     = Qsh(i_qs_pos);
n_qs       = numel(Qs_src);

% Case 3 end-injections
src_bus = [src_bus_g;       src_bus_qs;  case3_fbus; case3_tbus];
Qsrc    = [Qg_src;          Qs_src;      case3_QF;   case3_QT  ];
ns      = numel(Qsrc);

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

%% Build the system matrix and right-hand side, then solve.
M   = Cf * Ct';
N   = -diag(Qf .* hf) * (M * diag(aQ));
A   = speye(ne) + N;
S   = Cf * Cs';
B   = diag(Qf .* hf) * S * diag(Qsrc);
RPFCf = A \ B;
RPFCt = diag(aQ) * RPFCf;

%% Post-processing: shunt self-contributions.
%  Compute branch to-end self terms and propagate them with the same
%  system matrix used above.
self_toend = zeros(nl, 1);
for i = 1:nl
    if abs(Qs(i)) > 1e-9
        self_toend(i) = (Qe(i) / Qs(i)) * Qsf(i) + Qst(i);
    end
end

i_self    = find(self_toend > 1e-9);
n_self    = numel(i_self);
RPFCf_self  = zeros(ne, nl);
RPFCt_self = zeros(ne, nl);

if n_self > 0
    Css      = sparse(1:n_self, t(i_self), 1, n_self, nb+1);
    Ss       = Cf * Css';
    B_self = diag(Qf .* hf) * Ss * diag(self_toend(i_self));
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
for i = 1:n_case3
    src_names{ng+n_qs+2*i-1} = sprintf('E_%d_%d_f', case3_fbus(i), case3_tbus(i));
    src_names{ng+n_qs+2*i  } = sprintf('E_%d_%d_t', case3_fbus(i), case3_tbus(i));
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
row_labels = [d_labels; lf_labels; lt_labels];

% Append Case 3 branch rows (each end sourced entirely by that branch's injection)
n_tot = ns + nl;
if n_case3 > 0
    c3f_data = zeros(n_case3, n_tot);
    c3t_data = zeros(n_case3, n_tot);
    for i = 1:n_case3
        c3f_data(i, ng+n_qs+2*i-1) = -case3_QF(i);
        c3t_data(i, ng+n_qs+2*i  ) =  case3_QT(i);
    end
    RPFC_data    = [RPFC_data;    c3f_data; c3t_data   ];
    row_labels = [row_labels; case3_flab; case3_tlab];
end

RPFC_tbl = array2table(RPFC_data, 'VariableNames', all_names);
RPFC_tbl = addvars(RPFC_tbl, row_labels, 'Before', 1, 'NewVariableNames', 'Line');

result = strrep(casedata, 'cases/', 'results/');
writetable(RPFC_tbl, [result '_rpfc.csv']);
