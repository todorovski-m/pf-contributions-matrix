function varargout = pfc_checks(action, varargin)
%PFC_CHECKS  Options and validity checks shared by APFC and RPFC.
%
%   OPTS = PFC_CHECKS('options', OPTS) fills in the default options and
%   validates them. Fields (all optional):
%     tol      threshold for deletion and classification decisions
%              (default 1e-6, in the units of the MATPOWER case: MW, Mvar)
%     tol_div  guard against division by (near) zero and threshold for
%              negligible self-contribution terms (default 1e-9)
%     write    write results/<case>_apfc.csv or _rpfc.csv (default true;
%              skipped for struct input without a name)
%     solved   CASEDATA already holds a solved power flow, so RUNPF is not
%              called (default false)
%     name     case name used for the output file when CASEDATA is a struct
%              (default '': no file is written)
%
%   Checks (each can be unit tested with crafted inputs, see tests/unit/):
%
%   [ACYCLIC, ORDER] = PFC_CHECKS('acyclic', M)
%     Kahn's topological sort of the digraph with an edge j -> i for every
%     non-zero M(i,j) (branch j feeds the from-bus of branch i).
%   RHO = PFC_CHECKS('spectral_radius', N, ACYCLIC)
%     0 exactly if ACYCLIC; otherwise the largest absolute eigenvalue of the
%     cyclic core of N (the nodes left after repeatedly removing nodes
%     without incoming or without outgoing edges; the removed parts are
%     triangular blocks with zero eigenvalues). EIG is used for a core of at
%     most 2000 nodes, EIGS otherwise.
%   [ACYCLIC, RHO, W] = PFC_CHECKS('cycle', N, M, KIND)
%     KIND 'apfc': if cyclic, warning pfc:apfc:cyclic when RHO < 1 and error
%     pfc:apfc:spectralRadius otherwise. KIND 'rpfc': if cyclic, error
%     pfc:rpfc:spectralRadius when RHO >= 1 (no warning when RHO < 1,
%     reactive circulating flows are allowed).
%   W = PFC_CHECKS('attenuation', A, LABELS, SLACK)
%     Active attenuation A = Pt./Pf must lie in (0, 1 + SLACK]: error
%     pfc:apfc:nonPositiveAttenuation for A <= 0 (a branch fed from both
%     ends), warning pfc:apfc:attenuationAboveOne for A > 1 + SLACK.
%   PFC_CHECKS('inflow', TOTAL, F, TOL_DIV, KIND)
%     Error pfc:<KIND>:zeroInflow, naming the buses, if the total inflow
%     TOTAL(k) of a bus k with an outgoing branch (F == k) is below TOL_DIV.
%   [W, IDX] = PFC_CHECKS('series_flow', QS, QE, LABELS, TOL)
%     Warning pfc:rpfc:negativeSeriesFlow listing the branches whose series
%     flow at the from-end QS or at the to-end QE is below -TOL (negative
%     reactive attenuation).
%
%   W is a cell array of the warning identifiers that were issued.

switch action
    case 'options'
        varargout{1} = check_options(varargin{:});
    case 'acyclic'
        [varargout{1}, varargout{2}] = kahn(varargin{:});
    case 'spectral_radius'
        varargout{1} = spectral_radius(varargin{:});
    case 'cycle'
        [varargout{1}, varargout{2}, varargout{3}] = check_cycle(varargin{:});
    case 'attenuation'
        varargout{1} = check_attenuation(varargin{:});
    case 'inflow'
        check_inflow(varargin{:});
    case 'series_flow'
        [varargout{1}, varargout{2}] = check_series_flow(varargin{:});
    otherwise
        error('pfc:checks:unknownAction', 'pfc_checks: unknown action ''%s''.', action);
end
end

function opts = check_options(opts)
if nargin < 1 || isempty(opts)
    opts = struct();
end
if ~isstruct(opts)
    error('pfc:options:type', 'opts must be a struct.');
end
defaults = struct('tol', 1e-6, 'tol_div', 1e-9, 'write', true, 'solved', false, 'name', '');
known = fieldnames(defaults);
given = fieldnames(opts);
for k = 1:numel(given)
    if ~any(strcmp(given{k}, known))
        error('pfc:options:unknownField', 'Unknown option ''%s''. Known options: %s.', ...
            given{k}, strjoin(known', ', '));
    end
end
for k = 1:numel(known)
    if ~isfield(opts, known{k})
        opts.(known{k}) = defaults.(known{k});
    end
end
for k = {'tol', 'tol_div'}
    v = opts.(k{1});
    if ~(isnumeric(v) && isscalar(v) && isreal(v) && v > 0 && isfinite(v))
        error('pfc:options:value', 'opts.%s must be a positive finite scalar.', k{1});
    end
end
if ~ischar(opts.name)
    error('pfc:options:value', 'opts.name must be a character vector.');
end
opts.write = logical(opts.write);
opts.solved = logical(opts.solved);
end

function [acyclic, order] = kahn(M)
% Kahn's algorithm on the pattern of M (edge j -> i for M(i,j) ~= 0).
n = size(M, 1);
P = spones(sparse(M));
indeg = full(sum(P, 2));
stack = find(indeg == 0);
order = zeros(n, 1);
cnt = 0;
while ~isempty(stack)
    j = stack(end);
    stack(end) = [];
    cnt = cnt + 1;
    order(cnt) = j;
    nbr = find(P(:, j));
    for k = 1:numel(nbr)
        i = nbr(k);
        indeg(i) = indeg(i) - 1;
        if indeg(i) == 0
            stack(end + 1) = i; %#ok<AGROW>
        end
    end
end
acyclic = (cnt == n);
order = order(1:cnt);
end

function core = cyclic_core(P)
% Nodes left after repeatedly removing nodes without incoming edges (in the
% remaining graph) and nodes without outgoing edges.
keep = true(size(P, 1), 1);
changed = true;
while changed
    Q = P(keep, keep);
    idx = find(keep);
    no_in = full(sum(Q, 2)) == 0;
    no_out = full(sum(Q, 1)).' == 0;
    drop = no_in | no_out;
    changed = any(drop);
    keep(idx(drop)) = false;
end
core = find(keep);
end

function rho = spectral_radius(N, acyclic)
if nargin < 2
    acyclic = kahn(N);
end
if acyclic
    rho = 0;
    return
end
core = cyclic_core(spones(sparse(N)));
Nc = N(core, core);
if numel(core) <= 2000
    rho = max(abs(eig(full(Nc))));
else
    rho = abs(eigs(sparse(Nc), 1, 'lm'));
end
end

function [acyclic, rho, w] = check_cycle(N, M, kind)
w = {};
acyclic = kahn(M);
rho = spectral_radius(N, acyclic);
if acyclic
    return
end
switch kind
    case 'apfc'
        if rho < 1
            id = 'pfc:apfc:cyclic';
            warning(id, ['The active power flow digraph is cyclic (spectral radius of N = %.17g < 1); ' ...
                'A = I + N is invertible by the Neumann series, but the acyclicity theorem does not apply.'], rho);
            w{end + 1} = id;
        else
            error('pfc:apfc:spectralRadius', ['The active power flow digraph is cyclic and the spectral radius of N is ' ...
                '%.17g >= 1: the sufficient condition rho(N) < 1 for the invertibility of A = I + N is not met.'], rho);
        end
    case 'rpfc'
        if rho >= 1
            error('pfc:rpfc:spectralRadius', ['The reactive power flow digraph is cyclic and the spectral radius of N^Q is ' ...
                '%.17g >= 1: the sufficient condition rho(N^Q) < 1 for the invertibility of A^Q = I + N^Q ' ...
                '(Neumann series) is not met.'], rho);
        end
    otherwise
        error('pfc:checks:kind', 'Unknown kind ''%s''.', kind);
end
end

function w = check_attenuation(a, labels, slack)
w = {};
bad = find(a <= 0);
if ~isempty(bad)
    error('pfc:apfc:nonPositiveAttenuation', ['Active attenuation Pt/Pf <= 0 (branch fed from both ends, ' ...
        'not handled) for branch(es): %s.'], strjoin(labels(bad).', ', '));
end
hi = find(a > 1 + slack);
if ~isempty(hi)
    id = 'pfc:apfc:attenuationAboveOne';
    warning(id, 'Active attenuation Pt/Pf exceeds 1 by more than %.3g for branch(es): %s.', ...
        slack, strjoin(labels(hi).', ', '));
    w{end + 1} = id;
end
end

function check_inflow(total, f, tol_div, kind)
has_out = false(numel(total), 1);
has_out(unique(f)) = true;
bad = find(has_out & total(:) < tol_div);
if ~isempty(bad)
    error(['pfc:' kind ':zeroInflow'], ['Total inflow below %.3g at bus(es) %s, which have an outgoing ' ...
        'branch: the shares cannot be computed.'], tol_div, strjoin(arrayfun(@(k) sprintf('%d', k), bad.', ...
        'UniformOutput', false), ', '));
end
end

function [w, idx] = check_series_flow(qs, qe, labels, tol)
w = {};
idx = find(qs < -tol | qe < -tol);
if ~isempty(idx)
    id = 'pfc:rpfc:negativeSeriesFlow';
    warning(id, 'Negative series reactive flow (negative attenuation) in branch(es): %s.', ...
        strjoin(labels(idx).', ', '));
    w{end + 1} = id;
end
end
