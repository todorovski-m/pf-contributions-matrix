% Unit tests for pfc_checks.m with crafted inputs (brief T12b, Step 2, acceptance 3).
% Run with: zsh tests/octave_in_scratch.sh test_checks
warning('off', 'all');   % the tests read the identifiers with lastwarn; printing is not needed

% acyclic chain 1 -> 2 -> 3: M(i,j) ~= 0 means branch j feeds the from-bus of branch i
M = sparse([2 3], [1 2], 1, 3, 3);
[ac, ord] = pfc_checks('acyclic', M);
pfc_report('acyclic: chain is acyclic with a valid topological order', ac && isequal(ord(:)', [1 2 3]));
N = -0.7 * M;
pfc_report('spectral radius: exactly 0 for an acyclic digraph', pfc_checks('spectral_radius', N, ac) == 0);
lastwarn('', '');
[ac2, rho2, w2] = pfc_checks('cycle', N, M, 'apfc');
[~, id] = lastwarn;
pfc_report('cycle: acyclic active digraph gives no warning', ac2 && rho2 == 0 && isempty(w2) && isempty(id));

% cyclic active digraph (2-cycle) with spectral radius 0.5 < 1: warning
M2 = sparse([1 2], [2 1], 1, 2, 2);
lastwarn('', '');
[ac3, rho3, w3] = pfc_checks('cycle', -0.5 * M2, M2, 'apfc');
[~, id] = lastwarn;
pfc_report('cycle: cyclic active digraph with rho < 1 gives warning pfc:apfc:cyclic', ...
    ~ac3 && abs(rho3 - 0.5) < 1e-15 && strcmp(id, 'pfc:apfc:cyclic') && isequal(w3, {'pfc:apfc:cyclic'}), ...
    sprintf('acyclic=%d rho=%.17g id=%s', ac3, rho3, id));

% cyclic active digraph with spectral radius 1.2 >= 1: error
[ok, msg] = pfc_expect_error(@() pfc_checks('cycle', -1.2 * M2, M2, 'apfc'), 'pfc:apfc:spectralRadius');
pfc_report('cycle: cyclic active digraph with rho >= 1 is an error', ok, msg);
[ok, msg] = pfc_expect_error(@() pfc_checks('cycle', -1.0 * M2, M2, 'apfc'), 'pfc:apfc:spectralRadius');
pfc_report('cycle: rho = 1 exactly is an error', ok, msg);

% reactive: rho < 1 allowed without warning, rho >= 1 is an error naming the sufficient condition
lastwarn('', '');
[ac4, rho4, w4] = pfc_checks('cycle', -0.3 * M2, M2, 'rpfc');
[~, id] = lastwarn;
pfc_report('cycle: cyclic reactive digraph with rho < 1 gives no warning and no error', ~ac4 && isempty(w4) && isempty(id));
[ok, msg] = pfc_expect_error(@() pfc_checks('cycle', -1.5 * M2, M2, 'rpfc'), 'pfc:rpfc:spectralRadius');
pfc_report('cycle: cyclic reactive digraph with rho >= 1 is an error', ok && ~isempty(strfind(msg, 'rho(N^Q) < 1')), msg);

% cyclic core: a cycle with a long acyclic tail (Jordan chain) keeps the exact radius
n = 300;
Mt = sparse([2:n, 1], [1:n-1, 2], 1, n, n);   % chain 1 -> 2 -> ... -> n, plus 2 -> 1 (cycle between 1 and 2)
Nt = -0.4 * Mt;
rho_t = pfc_checks('spectral_radius', Nt, pfc_checks('acyclic', Mt));
pfc_report('spectral radius: cycle with a 298-node tail gives the radius of the cycle (0.4)', abs(rho_t - 0.4) < 1e-12, ...
    sprintf('rho = %.17g', rho_t));
% eigs path: 2500 independent 2-cycles of weight 0.3 (core larger than 2000 nodes)
k = 2500; i1 = 1:2:2*k; i2 = 2:2:2*k;
Me = sparse([i1 i2], [i2 i1], 1, 2*k, 2*k);
rho_e = pfc_checks('spectral_radius', -0.3 * Me, false);
pfc_report('spectral radius: eigs path for a 5000-node cyclic core', abs(rho_e - 0.3) < 1e-10, sprintf('rho = %.17g', rho_e));

% attenuation: a <= 0 is an error naming the branch, a > 1 + slack is a warning
lab = {'1 (1-2)'; '2 (2-3)'; '3 (3-4)'};
[ok, msg] = pfc_expect_error(@() pfc_checks('attenuation', [0.9; -0.1; 0.95], lab, 1e-9), 'pfc:apfc:nonPositiveAttenuation');
pfc_report('attenuation: a <= 0 is an error naming the branch', ok && ~isempty(strfind(msg, '2 (2-3)')), msg);
[ok, msg] = pfc_expect_error(@() pfc_checks('attenuation', [0.9; 0; 0.95], lab, 1e-9), 'pfc:apfc:nonPositiveAttenuation');
pfc_report('attenuation: a = 0 is an error', ok, msg);
lastwarn('', '');
w = pfc_checks('attenuation', [0.9; 1 + 1e-6; 1 + 1e-12], lab, 1e-9);
[m, id] = lastwarn;
pfc_report('attenuation: a > 1 + slack is a warning naming only that branch', strcmp(id, 'pfc:apfc:attenuationAboveOne') ...
    && ~isempty(strfind(m, '2 (2-3)')) && isempty(strfind(m, '3 (3-4)')), m);
lastwarn('', '');
w = pfc_checks('attenuation', [0.9; 1; 1 + 1e-12], lab, 1e-9);
[~, id] = lastwarn;
pfc_report('attenuation: values in (0, 1 + slack] give no warning', isempty(w) && isempty(id));

% zero inflow at a bus with an outgoing branch is an error naming the bus
[ok, msg] = pfc_expect_error(@() pfc_checks('inflow', [0; 5; 0], [1; 2], 1e-9, 'apfc'), 'pfc:apfc:zeroInflow');
pfc_report('inflow: zero inflow at a bus with an outgoing branch is an error naming the bus', ok && ~isempty(strfind(msg, 'bus(es) 1,')), msg);
[ok, msg] = pfc_expect_error(@() pfc_checks('inflow', [5; 0; 0], [1; 2], 1e-9, 'rpfc'), 'pfc:rpfc:zeroInflow');
pfc_report('inflow: reactive identifier pfc:rpfc:zeroInflow', ok, msg);
try
    pfc_checks('inflow', [5; 4; 0], [1; 2], 1e-9, 'apfc'); ok = true;
catch
    ok = false;
end
pfc_report('inflow: zero inflow at a bus without an outgoing branch (the sink) is not an error', ok);

% negative series flow: warning listing the branch
lastwarn('', '');
[w, idx] = pfc_checks('series_flow', [6.4; 12.4], [-2.1; 2.2], {'1 (2-1)'; '2 (1-2)'}, 1e-6);
[m, id] = lastwarn;
pfc_report('series flow: negative series flow gives warning pfc:rpfc:negativeSeriesFlow naming the branch', ...
    strcmp(id, 'pfc:rpfc:negativeSeriesFlow') && isequal(idx, 1) && ~isempty(strfind(m, '1 (2-1)')) && isempty(strfind(m, '2 (1-2)')), m);
lastwarn('', '');
[w, idx] = pfc_checks('series_flow', [6.4; 12.4], [2.1; 2.2], {'1 (2-1)'; '2 (1-2)'}, 1e-6);
pfc_report('series flow: positive series flows give no warning', isempty(w) && isempty(idx));
