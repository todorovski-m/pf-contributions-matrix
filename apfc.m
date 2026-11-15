function [APFC_tbl, mpc, info] = apfc(casedata, opts)
%APFC  Active Power Flow Components (APFC) for a MATPOWER case.
%
%   [APFC_TBL, MPC, DIAG] = APFC(CASEDATA, OPTS) runs a power flow for
%   CASEDATA, constructs an expanded network that includes demand-as-branches
%   and a virtual sink bus, builds the APFC allocation matrix, and computes
%   generator contributions through each branch (including demand "branches").
%
%   CASEDATA is a case name (for example 'cases/case9'), an unsolved MATPOWER
%   case struct, or, with OPTS.solved = true, a struct that already holds a
%   solved power flow (for example the output of RUNPF); then RUNPF is not
%   called. OPTS is optional; see PFC_CHECKS('options') for the fields
%   (tol, tol_div, write, solved, name). The defaults reproduce the earlier
%   behaviour: APFC(CASEDATA) runs RUNPF and writes results/<case>_apfc.csv.
%
%   DIAG is a struct with sizes, the number of non-zeros of N and A, whether
%   N, A and B are sparse, an estimate of the 1-norm condition number of A
%   (CONDEST), the range of the attenuation factors, and the numbers of
%   removed, flipped and converted elements.
%
%   Conversions of negative injections (as in the iterative method):
%     a bus with negative active demand (PD < -tol) is an active source
%       (column D_b);
%     a generator with negative active output (PG < -tol) is an active
%       demand (row G_b).
%
%   Source columns: G_b (generators with PG > 0), then D_b (buses with
%   negative active demand). Rows: D_b (demands), G_b (generators with
%   negative output), then the from-ends and the to-ends of the branches.
%
%   Output APFC_TBL is a table with contribution values and line labels.
%   MPC is the modified MATPOWER case after preprocessing.
%

%% Options, load case and run MATPOWER power flow
if nargin < 1 || isempty(casedata)
    casedata = 'cases/case1';
end
if nargin < 2
    opts = struct();
end
opts = pfc_checks('options', opts);
tol = opts.tol;

define_constants;
[mpc, name] = load_solved(casedata, opts);

%% Remove inactive network elements

% remove branches carrying zero power
i_del_br = abs(mpc.branch(:,PF)) < tol & abs(mpc.branch(:,PT)) < tol;
mpc.branch = mpc.branch(~i_del_br,:);

% remove generators producing zero power
i_del_gen = abs(mpc.gen(:,PG)) < tol;
mpc.gen = mpc.gen(~i_del_gen,:);

%% Normalize branch orientations
%  MATPOWER stores PT/QT as power from the to-bus into the branch; sign-correct to power leaving the branch into the to-bus.
mpc.branch(:, [PT QT]) = -mpc.branch(:, [PT QT]);

% flip branches that still have PF < 0
i_flip = mpc.branch(:, PF) < 0;
if any(i_flip)
    % swap endpoints
    mpc.branch(i_flip, [F_BUS T_BUS]) = mpc.branch(i_flip, [T_BUS F_BUS]);
    % flip flows
    mpc.branch(i_flip, [PF QF PT QT]) = -mpc.branch(i_flip, [PT QT PF QF]);
end

nb = size(mpc.bus, 1);  % number of buses
sink = nb + 1;          % index of virtual sink bus

%% Demand-as-branches construction
%  Each demand is converted to a branch whose "to" bus = sink
i = mpc.bus(:, PD) > tol;
d_bus_pd = mpc.bus(i, BUS_I);   % buses with demand
Pd_load  = mpc.bus(i, PD);      % demand values
% generators with negative output are demands (row G_b)
i = mpc.gen(:, PG) < -tol;
d_bus_pg = mpc.gen(i, GEN_BUS);
Pd_pgneg = -mpc.gen(i, PG);
d_bus = [d_bus_pd; d_bus_pg];
Pd    = [Pd_load;  Pd_pgneg];
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

%% Source data: generators with positive output, then buses with negative demand
i = mpc.gen(:, PG) > tol;
src_bus_g = mpc.gen(i, GEN_BUS);
Pg_src    = mpc.gen(i, PG);
i = mpc.bus(:, PD) < -tol;
src_bus_pd = mpc.bus(i, BUS_I);
Pd_src     = -mpc.bus(i, PD);
g_bus = [src_bus_g; src_bus_pd];
Pg    = [Pg_src;    Pd_src];
ng    = numel(Pg);
n_gpos = numel(Pg_src);

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

%% Validity checks (see pfc_checks.m)
%  branch labels: index in the returned mpc.branch and the oriented bus pair
br_labels = arrayfun(@(k) sprintf('%d (%d-%d)', k, f(k), t(k)), (1:nl)', 'UniformOutput', false);
pfc_checks('inflow', Ptot, f, opts.tol_div, 'apfc');
w_att = pfc_checks('attenuation', a(1:nl), br_labels, opts.tol_div);

%% Build APFC system matrix A
% Construct the parent-child adjacency matrix (ne x ne):
% M(i,j) = 1 if branch j ends at the "from" bus of branch i, 0 otherwise
M = Cf * Ct';

% Map reciprocal total inflows to the "from" bus of each branch
hf = h(f);

% Compute the contribution matrix N:
% each element N(i,j) = - Pf(i) * h(from_bus_i) * a(j) * M(i,j)
N = -spdiags(Pf .* hf, 0, ne, ne) * M * spdiags(a, 0, ne, ne);

% APFC system matrix: A = I + N
A = speye(ne) + N;

% acyclicity of the flow digraph and spectral radius of N
[acyclic, rho, w_cyc] = pfc_checks('cycle', N, M, 'apfc');

%% Build RHS matrix B for each generator
%  Each generator contributes only to branches that originate at its bus.
%  B(ell, k) = Pf(ell) * h(f(ell)) * Pg(k) if f(ell) = g_bus(k), else 0
S = Cf * Cg'; % maps generators to outgoing branches (ne x ng)
B = spdiags(Pf .* hf, 0, ne, ne) * S * spdiags(Pg, 0, ng, ng); % (ne x ng), sparse

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
for i = 1:n_gpos
    gen_names{i} = sprintf('G_%i', g_bus(i));
end
for i = n_gpos+1:ng
    gen_names{i} = sprintf('D_%i', g_bus(i));
end

line_names = cell(ne,1);
nd_pd = numel(d_bus_pd);
for i = 1:nd_pd
    line_names{i} = sprintf('D_%d', d_bus(i));
end
for i = nd_pd+1:nd
    line_names{i} = sprintf('G_%d', d_bus(i));
end
for i = 1:nl
    line_names{nd + i}       = sprintf('L_{f.%d-%d}', f(i), t(i));
    line_names{nd + nl + i}  = sprintf('L_{t.%d-%d}', f(i), t(i));
end

% Disambiguate duplicate row labels for parallel branches (same scheme as rpfc.m):
% the second and later branches with the same orientation get the suffix _2, _3, ...
seen_rows = containers.Map('KeyType','char','ValueType','int32');
for i = 1:nl
    nm = line_names{nd + i};
    if isKey(seen_rows, nm)
        seen_rows(nm) = seen_rows(nm) + 1;
        line_names{nd + i}      = sprintf('%s_%d', nm, seen_rows(nm));
        line_names{nd + nl + i} = sprintf('%s_%d', line_names{nd + nl + i}, seen_rows(nm));
    else
        seen_rows(nm) = 1;
    end
end

%% Diagnostics
info = struct();
info.case = name;
info.n_bus = nb;
info.n_branch = nl;
info.n_demand = nd;
info.n_source = ng;
info.n_e = ne;
info.nnz_N = nnz(N);
info.nnz_A = nnz(A);
info.issparse = struct('N', issparse(N), 'A', issparse(A), 'B', issparse(B));
info.condest_A = condest(A);
if nl > 0
    info.attenuation_range = [min(a(1:nl)), max(a(1:nl))];
else
    info.attenuation_range = zeros(1, 0);
end
info.removed = struct('branches', nnz(i_del_br), 'generators', nnz(i_del_gen));
info.flipped_branches = nnz(i_flip);
info.converted = struct('negative_demand_to_source', numel(src_bus_pd), ...
                        'negative_generation_to_demand', numel(d_bus_pg));
info.acyclic = acyclic;
info.spectral_radius = rho;
info.warnings = [w_att, w_cyc];
info.opts = opts;

%% Output table
APFC_tbl = array2table(full(APFC_tbl), 'VariableNames', gen_names);   % full only for the table
APFC_tbl = addvars(APFC_tbl, line_names, 'Before', 1, 'NewVariableNames', 'Line');

if opts.write && ~isempty(name)
    result = strrep(name, 'cases/', 'results/');
    writetable(APFC_tbl, [result '_apfc.csv']);
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
