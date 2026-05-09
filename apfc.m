function [APFC_tbl, mpc] = apfc(casedata)
%APFC  Active Power Flow Components (APFC) for a MATPOWER case.
%
%   [APFC_TBL, MPC] = APFC(CASEDATA) runs a power flow for CASEDATA, constructs
%   an expanded network that includes demand-as-branches and a virtual sink
%   bus, builds the APFC allocation matrix, and computes generator
%   contributions through each branch (including demand "branches").
%
%   Output APFC_TBL is a table with contribution values and line labels.
%   MPC is the modified MATPOWER case after preprocessing.
%

%% Load case and run MATPOWER power flow
if nargin < 1
    casedata = 'cases/case1';
end

define_constants;
opt = mpoption('verbose', 0, 'out.all', 0);

mpc = loadcase(casedata);
mpc = runpf(mpc, opt);

%% Remove inactive network elements

% remove branches carrying zero power
i = abs(mpc.branch(:,PF)) < 1e-6 & abs(mpc.branch(:,PT)) < 1e-6;
mpc.branch = mpc.branch(~i,:);

% remove generators producing zero power
i = abs(mpc.gen(:,PG)) < 1e-6;
mpc.gen = mpc.gen(~i,:);

%% Normalize branch orientations
%  MATPOWER stores PT/QT as power from the to-bus into the branch; sign-correct to power leaving the branch into the to-bus.
mpc.branch(:, [PT QT]) = -mpc.branch(:, [PT QT]);

% flip branches that still have PF < 0
i = mpc.branch(:, PF) < 0;
if any(i)
    % swap endpoints
    mpc.branch(i, [F_BUS T_BUS]) = mpc.branch(i, [T_BUS F_BUS]);
    % flip flows
    mpc.branch(i, [PF QF PT QT]) = -mpc.branch(i, [PT QT PF QF]);
end

nb = size(mpc.bus, 1);  % number of buses
sink = nb + 1;          % index of virtual sink bus

%% Demand-as-branches construction
%  Each demand is converted to a branch whose "to" bus = sink
i = mpc.bus(:, PD) > 0;
d_bus = mpc.bus(i, BUS_I);   % buses with demand
Pd = mpc.bus(i, PD);         % demand values
nd = numel(Pd);              % number of demands

%% Branch data (including appended demand branches)
nl = size(mpc.branch, 1);   % number of physical lines
ne = nl + nd;               % total number of branches (lines + demands)

% assemble unified from/to bus vectors
f  = [mpc.branch(:, F_BUS); d_bus];
t  = [mpc.branch(:, T_BUS); repmat(sink, nd, 1)];

% assemble unified flows
Pf = [mpc.branch(:, PF); Pd];
Pt = [mpc.branch(:, PT); Pd];

%% Generator data
g_bus = mpc.gen(:, GEN_BUS);
Pg    = mpc.gen(:, PG);
ng    = numel(Pg);

%% Build incidence matrices
%  Cf, Ct: branch → bus incidence
%  Cg: generator → bus incidence
Cf = sparse(1:ne, f, 1, ne, nb+1);
Ct = sparse(1:ne, t, 1, ne, nb+1);
Cg = sparse(1:ng, g_bus, 1, ng, nb+1);

%% Compute total inflow per bus and, weight h and attenuation a
Ptot = Ct' * Pt + Cg' * Pg;
h = 1 ./ Ptot;
a = Pt ./ Pf;

%% Build APFC system matrix A
% Construct the parent-child adjacency matrix (ne x ne):
% M(i,j) = 1 if branch j ends at the "from" bus of branch i, 0 otherwise
M = Cf * Ct';

% Map reciprocal total inflows to the "from" bus of each branch
hf = h(f);

% Compute the contribution matrix N:
% each element N(i,j) = - Pf(i) * h(from_bus_i) * a(j) * M(i,j)
N = -diag(Pf .* hf) * M * diag(a);

% APFC system matrix: A = I + N
A = eye(ne) + N;

%% Build RHS matrix B for each generator
%  Each generator contributes only to branches that originate at its bus.
%  B(ell, k) = Pf(ell) * h(f(ell)) * Pg(k) if f(ell) = g_bus(k), else 0
S = Cf * Cg'; % maps generators to outgoing branches (ne x ng)
B = diag(Pf .* hf) * S * diag(Pg); % (ne x ng)

%% Solve contribution system
%  APFC_tbl = contributions at "from" side of each branch
%  APFCt    = contributions at "to" side (after losses)
% A = I + N, A^{-1} = \sum_{k=0}^{\infty}(-N)^k = I - N + N^2 - N^3 ...
APFC_tbl = A \ B;
APFCt = APFC_tbl .* a;

% reorder rows: demands, lines-from, lines-to
APFC_tbl = [APFC_tbl(nl+1:end,:); APFC_tbl(1:nl,:); APFCt(1:nl,:)];

%% Build output table with readable labels
gen_names = cell(1,ng);
for i = 1:ng
    gen_names{i} = sprintf('G_%i', g_bus(i));
end

line_names = cell(ne,1);
for i = 1:nd
    line_names{i} = sprintf('D_%d', d_bus(i));
end
for i = 1:nl
    line_names{nd + i}       = sprintf('L_{f.%d-%d}', f(i), t(i));
    line_names{nd + nl + i}  = sprintf('L_{t.%d-%d}', f(i), t(i));
end

APFC_tbl = array2table(APFC_tbl, 'VariableNames', gen_names);
APFC_tbl = addvars(APFC_tbl, line_names, 'Before', 1, 'NewVariableNames', 'Line');

result = strrep(casedata, 'cases/', 'results/');
writetable(APFC_tbl, [result '_apfc.csv']);
