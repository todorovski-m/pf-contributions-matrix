% Unit tests of the APFC branches whose two terminals both draw active power (brief T14c): the pattern (+,-) as two
% terminal demands on crafted solved structs and on the synthetic case syn9_both_ends, the error
% pfc:apfc:negativeLoss for (-,+) and its message, the counter converted.case2_branches_to_demands, the order of the
% terminal-demand rows, invariance to a perturbation of the terminal flows by +/-1e-12, the attenuation check at
% the (0,0) boundary, and the equivalence with two loads.
% Run with: zsh tests/octave_in_scratch.sh test_both_ends
warning('off', 'all');
define_constants;
TOL = 1e-10;
sopt = @(n) struct('solved', true, 'write', false, 'name', n);
has = @(t, lab) any(strcmp(pfc_tbl_rows(t), lab));
rowsum = @(t, lab) sum(pfc_tbl_row(t, lab));
q = [1 0.9];
base = runpf(loadcase('cases/case9'), mpoption('verbose', 0, 'out.all', 0));
[~, ~, d0] = apfc(base, sopt('base'));

% (+,-): two terminal demands, not a carrier
[t, ~, d] = apfc(pfc_crafted_branch(1, -0.1, q(1), q(2)), sopt('pm'));
pfc_report('(+,-): rows C2_4 = PF and C2_5 = -PT, not a carrier', has(t, 'C2_4') && has(t, 'C2_5') ...
    && abs(rowsum(t, 'C2_4') - 1) < TOL && abs(rowsum(t, 'C2_5') - 0.1) < TOL && ~has(t, 'L_{f.4-5}_2') ...
    && d.n_branch == d0.n_branch);
pfc_report('(+,-): counter case2_branches_to_demands = 1, zero_terminal_demands = 0', ...
    d.converted.case2_branches_to_demands == 1 && d.converted.zero_terminal_demands == 0 && d.converted.zero_terminal_sources == 0);
pfc_report('counter case2_branches_to_demands = 0 in case9', d0.converted.case2_branches_to_demands == 0);

% (-,+): error with the message of the brief, labels not reversed
[ok, msg] = pfc_expect_error(@() apfc(pfc_crafted_branch(-0.1, 1, q(1), q(2)), sopt('mp')), 'pfc:apfc:negativeLoss');
pfc_report('(-,+): error pfc:apfc:negativeLoss with the message of the brief', ok && ~isempty(strfind(msg, ...
    'Active power is delivered at both ends of a branch (negative loss) for branch(es): 10 (4-5).')), msg);
[ok, msg] = pfc_expect_error(@() pfc_checks('negative_loss', {'1 (1-2)'}, {'2 (3-4)'}, {'3 (5-6)'}), 'pfc:apfc:negativeLoss');
pfc_report('negative_loss check: three sentences in the order to-end, from-end, both ends', ok && ...
    ~isempty(strfind(msg, ['Active power is delivered at the to-end of a branch with zero from-end flow (negative loss) for branch(es): 1 (1-2). ' ...
    'Active power is delivered at the from-end of a branch with zero to-end flow (negative loss) for branch(es): 2 (3-4). ' ...
    'Active power is delivered at both ends of a branch (negative loss) for branch(es): 3 (5-6).'])), msg);
[ok, msg] = pfc_expect_error(@() pfc_checks('negative_loss', {}, {'2 (3-4)'}), 'pfc:apfc:negativeLoss');
pfc_report('negative_loss check: the two-argument form as before', ok, msg);

% order of the terminal demands: from-terminals first, then to-terminals, each in branch order
m = pfc_crafted_branch(0, -0.2, q(1), q(2));                       % branch 10: (0,-), demand at bus 5
r = m.branch(10, :); r([PF QF PT QT]) = [0.3 1 0.4 -0.9]; m.branch(11, :) = r;   % branch 11: (+,-) 4-5
[t, ~, d] = apfc(m, sopt('order'));
rl = pfc_tbl_rows(t); c2 = rl(strncmp(rl, 'C2_', 3))';
pfc_report('order of the terminal demands: C2_4 (from-terminal of 11), C2_5 (10), C2_5_2 (to-terminal of 11)', ...
    isequal(c2, {'C2_4', 'C2_5', 'C2_5_2'}) && abs(rowsum(t, 'C2_4') - 0.3) < TOL && abs(rowsum(t, 'C2_5') - 0.2) < TOL ...
    && abs(rowsum(t, 'C2_5_2') - 0.4) < TOL && d.converted.case2_branches_to_demands == 1 && d.converted.zero_terminal_demands == 1, ...
    strjoin(c2, ' '));

% the attenuation check stays reachable at the (0,0) boundary (terminals of magnitude exactly tol are not deleted)
[ok, msg] = pfc_expect_error(@() apfc(pfc_crafted_branch(1e-6, -1e-6, q(1), q(2)), sopt('b')), 'pfc:apfc:nonPositiveAttenuation');
pfc_report('attenuation check: reached by a (0,0) boundary branch with PF = tol, PT = -tol', ok, msg);

% the synthetic case
s0 = runpf(loadcase('cases/syn9_both_ends'), mpoption('verbose', 0, 'out.all', 0));
k = size(s0.branch, 1);
[t0, ~, d] = apfc(s0, sopt('syn'));
pfc_report('syn9_both_ends: branch 2-8 gives C2_2 and C2_8, counter 1', has(t0, 'C2_2') && has(t0, 'C2_8') ...
    && d.converted.case2_branches_to_demands == 1 && abs(rowsum(t0, 'C2_2') - s0.branch(k, PF)) < TOL ...
    && abs(rowsum(t0, 'C2_8') - s0.branch(k, PT)) < TOL);
[~, ~, d] = apfc('cases/syn9_both_ends', struct('write', false));
pfc_report('syn9_both_ends from the case name (runpf inside apfc)', d.converted.case2_branches_to_demands == 1);
[tq, ~, ~] = rpfc('cases/syn9_both_ends', struct('write', false));
cq = pfc_tbl_names(tq);
pfc_report('syn9_both_ends rpfc: the branch 2-8 is a reactive Case-3 branch (columns E_2_8_f and E_2_8_t)', ...
    any(strcmp(cq, 'E_2_8_f')) && any(strcmp(cq, 'E_2_8_t')));

% invariance to a perturbation of the terminal flows by +/-1e-12
mx = 0; msg = '';
for e = [1e-12, -1e-12]
    s1 = s0; s1.branch(k, [PF PT]) = s1.branch(k, [PF PT]) + e;
    [t1, ~] = apfc(s1, sopt('p'));
    [m1, g1] = pfc_map_diff(t0, t1); mx = max(mx, m1); msg = [msg g1]; %#ok<AGROW>
end
pfc_report(sprintf('invariance to +/-1e-12 at the terminal flows of 2-8 (largest difference %.3g)', mx), mx <= TOL, msg);

% equivalence with two loads (generators with negative output, rows G_2_2 and G_8)
rep = pfc_replace_both_ends(s0, k);
[tb, ~] = apfc(rep, sopt('rep'));
rb = pfc_tbl_rows(tb); g2 = rb(strncmp(rb, 'G_2', 3)); g8 = rb(strncmp(rb, 'G_8', 3));
[m_other, g1] = pfc_map_diff(t0, tb, {}, {}, {'C2_2', 'C2_8'}, [g2(:); g8(:)]);
[m_own, g2m] = pfc_map_diff(t0, tb, {'C2_2', g2{1}, 1; 'C2_8', g8{1}, 1}, {}, setdiff(pfc_tbl_rows(t0), {'C2_2', 'C2_8'}), ...
    setdiff(rb, [g2(:); g8(:)]));
fprintf('| syn9_both_ends | apfc | %.17g | %.17g |\n', m_other, m_own);
pfc_report('equivalence: syn9_both_ends apfc (other rows; the C2 rows against the load rows)', max(m_other, m_own) <= TOL, [g1 g2m]);
