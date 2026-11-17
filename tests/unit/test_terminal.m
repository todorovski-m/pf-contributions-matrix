% Unit tests of the terminal classification of branches with exactly one zero terminal (brief T15b):
% the nine sign patterns on crafted solved structs, the error pfc:apfc:negativeLoss, the counters, equal tables
% for the dead end at the from side and at the to side, invariance to a perturbation of the zero terminal by
% +/-1e-12, and the equivalence with the terminal injection that a dead-end branch represents.
% Run with: zsh tests/octave_in_scratch.sh test_terminal
warning('off', 'all');
define_constants;
TOL = 1e-10;
sopt = @(n) struct('solved', true, 'write', false, 'name', n);
has = @(t, lab) any(strcmp(pfc_tbl_rows(t), lab));
hascol = @(t, lab) any(strcmp(pfc_tbl_names(t), lab));

base = runpf(loadcase('cases/case9'), mpoption('verbose', 0, 'out.all', 0));
[~, ~, dpb] = apfc(base, sopt('base'));
[~, ~, dqb] = rpfc(base, sopt('base'));

%% the nine patterns, APFC (extra branch 4-5, index 10; reactive flows of a carrier)
q = [1 0.9];
[t, ~, d] = apfc(pfc_crafted_branch(1, 0.9, q(1), q(2)), sopt('pp'));
pfc_report('apfc pattern (+,+): carrier, no terminal demand', has(t, 'L_{f.4-5}_2') && d.n_branch == dpb.n_branch + 1 ...
    && d.converted.zero_terminal_demands == 0 && d.converted.zero_terminal_sources == 0);
[t, ~, d] = apfc(pfc_crafted_branch(-0.9, -1, q(1), q(2)), sopt('mm'));
pfc_report('apfc pattern (-,-): flipped, then carrier', has(t, 'L_{f.5-4}') && d.flipped_branches == dpb.flipped_branches + 1 ...
    && d.converted.zero_terminal_demands == 0);
[ok, msg] = pfc_expect_error(@() apfc(pfc_crafted_branch(1, -0.1, q(1), q(2)), sopt('pm')), 'pfc:apfc:nonPositiveAttenuation');
pfc_report('apfc pattern (+,-): error pfc:apfc:nonPositiveAttenuation as before', ok, msg);
[ok, msg] = pfc_expect_error(@() apfc(pfc_crafted_branch(-0.1, 1, q(1), q(2)), sopt('mp')), 'pfc:apfc:nonPositiveAttenuation');
pfc_report('apfc pattern (-,+): error pfc:apfc:nonPositiveAttenuation as before', ok, msg);
[t, ~, d] = apfc(pfc_crafted_branch(1, 0, q(1), q(2)), sopt('p0'));
[r, ~, x] = pfc_tbl(t);
pfc_report('apfc pattern (+,0): one terminal demand PF at the from bus (row C2_4), not a carrier', ...
    has(t, 'C2_4') && ~has(t, 'L_{f.4-5}_2') && abs(sum(x(strcmp(r, 'C2_4'), :)) - 1) < TOL ...
    && d.n_branch == dpb.n_branch && d.converted.zero_terminal_demands == 1 && d.converted.zero_terminal_sources == 0);
[t, ~, d] = apfc(pfc_crafted_branch(0, -1, q(1), q(2)), sopt('0m'));
[r, ~, x] = pfc_tbl(t);
pfc_report('apfc pattern (0,-): one terminal demand -PT at the to bus (row C2_5)', ...
    has(t, 'C2_5') && abs(sum(x(strcmp(r, 'C2_5'), :)) - 1) < TOL && d.converted.zero_terminal_demands == 1);
[ok, msg] = pfc_expect_error(@() apfc(pfc_crafted_branch(-1, 0, q(1), q(2)), sopt('m0')), 'pfc:apfc:negativeLoss');
pfc_report('apfc pattern (-,0): error pfc:apfc:negativeLoss naming the branch (from-end)', ...
    ok && ~isempty(strfind(msg, 'delivered at the from-end')) && ~isempty(strfind(msg, '10 (4-5)')), msg);
[ok, msg] = pfc_expect_error(@() apfc(pfc_crafted_branch(0, 1, q(1), q(2)), sopt('0p')), 'pfc:apfc:negativeLoss');
pfc_report('apfc pattern (0,+): error pfc:apfc:negativeLoss with the message of the brief', ok && ~isempty(strfind(msg, ...
    'Active power is delivered at the to-end of a branch with zero from-end flow (negative loss) for branch(es): 10 (4-5).')), msg);
[t, ~, d] = apfc(pfc_crafted_branch(0, 0, q(1), q(2)), sopt('00'));
pfc_report('apfc pattern (0,0): deleted as before', d.removed.branches == dpb.removed.branches + 1 && ~has(t, 'L_{f.4-5}_2'));
[ok, msg] = pfc_expect_error(@() pfc_checks('negative_loss', {'1 (1-2)'}, {'2 (3-4)'}), 'pfc:apfc:negativeLoss');
pfc_report('negative_loss check: both sub-patterns in one message', ok && ~isempty(strfind(msg, 'to-end of a branch')) ...
    && ~isempty(strfind(msg, 'from-end of a branch')) && ~isempty(strfind(msg, '1 (1-2)')) && ~isempty(strfind(msg, '2 (3-4)')), msg);
pfc_checks('negative_loss', {}, {});
pfc_report('negative_loss check: no error without such branches', true);

%% the nine patterns, RPFC (extra branch 4-5, index 10; active flows of a carrier)
p = [1 0.9];
[t, ~, d] = rpfc(pfc_crafted_branch(p(1), p(2), 1, 0.9), sopt('pp'));
pfc_report('rpfc pattern (+,+): carrier', has(t, 'L_{f.4-5}_2') && d.n_branch == dqb.n_branch + 1 ...
    && d.converted.zero_terminal_demands == 0 && d.converted.zero_terminal_sources == 0);
[t, ~, d] = rpfc(pfc_crafted_branch(p(1), p(2), -0.9, -1), sopt('mm'));
pfc_report('rpfc pattern (-,-): flipped, then carrier', has(t, 'L_{f.5-4}') && d.flipped_branches == dqb.flipped_branches + 1);
[t, ~, d] = rpfc(pfc_crafted_branch(p(1), p(2), 1, -0.5), sopt('pm'));
pfc_report('rpfc pattern (+,-): Case 2, two terminal demands as before', has(t, 'C2_4') && has(t, 'C2_5') ...
    && d.converted.case2_branches_to_demands == dqb.converted.case2_branches_to_demands + 1 && d.converted.zero_terminal_demands == 0);
[t, ~, d] = rpfc(pfc_crafted_branch(p(1), p(2), -1, 0.5), sopt('mp'));
pfc_report('rpfc pattern (-,+): Case 3, two terminal sources as before', hascol(t, 'E_4_5_f') && hascol(t, 'E_4_5_t') ...
    && d.converted.case3_branches_to_sources == dqb.converted.case3_branches_to_sources + 1 && d.converted.zero_terminal_sources == 0);
[t, ~, d] = rpfc(pfc_crafted_branch(p(1), p(2), 1, 0), sopt('p0'));
[r, ~, x] = pfc_tbl(t);
pfc_report('rpfc pattern (+,0): one demand QF at the from bus (C2_4 only), no transformer-shunt row (Case-2 rule)', ...
    has(t, 'C2_4') && ~has(t, 'C2_5') && abs(sum(x(strcmp(r, 'C2_4'), :)) - 1) < TOL && ~any(strncmp(r, 'Tf_', 3)) ...
    && d.converted.zero_terminal_demands == 1 && d.converted.case2_branches_to_demands == dqb.converted.case2_branches_to_demands);
[t, ~, d] = rpfc(pfc_crafted_branch(p(1), p(2), 0, -1), sopt('0m'));
[r, ~, x] = pfc_tbl(t);
pfc_report('rpfc pattern (0,-): one demand -QT at the to bus (C2_5 only)', has(t, 'C2_5') && ~has(t, 'C2_4') ...
    && abs(sum(x(strcmp(r, 'C2_5'), :)) - 1) < TOL && d.converted.zero_terminal_demands == 1);
[t, ~, d] = rpfc(pfc_crafted_branch(p(1), p(2), -1, 0), sopt('m0'));
[r, c, x] = pfc_tbl(t);
pfc_report('rpfc pattern (-,0): one source -QF at the from bus (column E_4_5_f only, its row holds QF)', ...
    hascol(t, 'E_4_5_f') && ~hascol(t, 'E_4_5_t') && abs(x(strcmp(r, 'L_{f.4-5}_2'), strcmp(c, 'E_4_5_f')) + 1) < TOL ...
    && d.converted.zero_terminal_sources == 1 && d.converted.case3_branches_to_sources == dqb.converted.case3_branches_to_sources);
[t, ~, d] = rpfc(pfc_crafted_branch(p(1), p(2), 0, 1), sopt('0p'));
[r, c, x] = pfc_tbl(t);
pfc_report('rpfc pattern (0,+): one source QT at the to bus (column E_4_5_t only, its row holds QT)', ...
    hascol(t, 'E_4_5_t') && ~hascol(t, 'E_4_5_f') && abs(x(strcmp(r, 'L_{t.4-5}_2'), strcmp(c, 'E_4_5_t')) - 1) < TOL ...
    && d.converted.zero_terminal_sources == 1);
[t, ~, d] = rpfc(pfc_crafted_branch(p(1), p(2), 0, 0), sopt('00'));
pfc_report('rpfc pattern (0,0): deleted as before', d.removed.branches == dqb.removed.branches + 1);

%% counters on the dead-end cases
cnt = {'syn9_dead_charge_from', 1, 0, 0, 1; 'syn9_dead_charge_to', 1, 0, 0, 1; 'syn9_dead_reactor_from', 1, 0, 1, 0; ...
       'syn9_dead_reactor_to', 1, 0, 1, 0; 'syn9_dead_shunt_to', 1, 0, 0, 0};
for k = 1:size(cnt, 1)
    [~, ~, d1] = apfc(['cases/' cnt{k, 1}], struct('write', false));
    [~, ~, d2] = rpfc(['cases/' cnt{k, 1}], struct('write', false));
    got = [d1.converted.zero_terminal_demands, d1.converted.zero_terminal_sources, ...
           d2.converted.zero_terminal_demands, d2.converted.zero_terminal_sources];
    pfc_report(sprintf('counters: %s (apfc %d, %d; rpfc %d, %d)', cnt{k, :}), isequal(got, [cnt{k, 2:5}]), mat2str(got));
end
c0 = {'case1', 'case9', 'case14', 'case24_ieee_rts', 'case30', 'case33', 'case118', 'tap_loop'};
z = 0;
for k = 1:numel(c0)
    [~, ~, d1] = apfc(['cases/' c0{k}], struct('write', false));
    [~, ~, d2] = rpfc(['cases/' c0{k}], struct('write', false));
    z = z + d1.converted.zero_terminal_demands + d1.converted.zero_terminal_sources + ...
        d2.converted.zero_terminal_demands + d2.converted.zero_terminal_sources;
end
pfc_report('counters: 0 in the shared cases, case24_ieee_rts and tap_loop', z == 0, sprintf('sum %d', z));

%% equal tables for the dead end at the from side and at the to side
pf = @(c) runpf(loadcase(['cases/' c]), mpoption('verbose', 0, 'out.all', 0));
for pair = {'charge', 'reactor'}
    sf = pf(['syn9_dead_' pair{1} '_from']); st = pf(['syn9_dead_' pair{1} '_to']);
    [ta, ~] = apfc(sf, sopt('f')); [tb, ~] = apfc(st, sopt('t'));
    [mx, msg] = pfc_map_diff(ta, tb);
    pfc_report(sprintf('from and to side equal: %s apfc (largest difference %.3g)', pair{1}, mx), mx <= TOL, msg);
    [ta, ~] = rpfc(sf, sopt('f')); [tb, ~] = rpfc(st, sopt('t'));
    if strcmp(pair{1}, 'charge')
        rp = {'L_{t.10-5}', 'L_{f.5-10}', -1}; cg = {'E_5_10_f', {'E_10_5_t'}};
    else
        rp = {}; cg = {};
    end
    [mx, msg] = pfc_map_diff(ta, tb, rp, cg);
    pfc_report(sprintf('from and to side equal: %s rpfc (largest difference %.3g)', pair{1}, mx), mx <= TOL, msg);
end

%% invariance to a perturbation of the zero terminal by +/-1e-12
inv = {'syn9_dead_charge_from', [PF QF]; 'syn9_dead_charge_to', [PT QT]; 'syn9_dead_reactor_from', [PF QF]; ...
       'syn9_dead_reactor_to', [PT QT]; 'syn9_dead_shunt_to', PT};
for k = 1:size(inv, 1)
    s0 = pf(inv{k, 1});
    [a0, ~] = apfc(s0, sopt('0')); [q0, ~] = rpfc(s0, sopt('0'));
    mx = 0; msg = '';
    for e = [1e-12, -1e-12]
        s1 = s0; s1.branch(end, inv{k, 2}) = e;
        [a1, ~] = apfc(s1, sopt('1')); [q1, ~] = rpfc(s1, sopt('1'));
        [m1, g1] = pfc_map_diff(a0, a1); [m2, g2] = pfc_map_diff(q0, q1);
        mx = max([mx, m1, m2]); msg = [msg, g1, g2]; %#ok<AGROW>
    end
    pfc_report(sprintf('invariance to +/-1e-12 at the zero terminal: %s (largest difference %.3g)', inv{k, 1}, mx), mx <= TOL, msg);
end

%% equivalence with the terminal injection (replacement by a generator at bus 5, or a bus shunt for the shunt case)
fprintf('| case | kind | largest difference of the other rows and columns | largest difference of the dead-end branch and its replacement |\n|---|---|---|---|\n');
eqv = {'syn9_dead_charge_from', 'L_{t.10-5}', 'E_10_5_t'; 'syn9_dead_charge_to', 'L_{f.5-10}', 'E_5_10_f'; ...
       'syn9_dead_reactor_from', '', ''; 'syn9_dead_reactor_to', '', ''; 'syn9_dead_shunt_to', '', ''};
for k = 1:size(eqv, 1)
    c = eqv{k, 1}; s0 = pf(c); shunt = strcmp(c, 'syn9_dead_shunt_to');
    rep = pfc_replace_dead_end(s0, 10, shunt);
    % APFC: row C2_5 of the original and row G_5 (the generator with negative output) of the replacement
    [ta, ~] = apfc(s0, sopt('o')); [tb, ~] = apfc(rep, sopt('r'));
    [m_other, g1] = pfc_map_diff(ta, tb, {}, {}, {'C2_5'}, {'G_5'});
    [m_own, g2] = pfc_map_diff(ta, tb, {'C2_5', 'G_5', 1}, {}, setdiff(pfc_tbl_rows(ta), {'C2_5'}), setdiff(pfc_tbl_rows(tb), {'G_5'}));
    fprintf('| %s | apfc | %.17g | %.17g |\n', c, m_other, m_own);
    pfc_report(sprintf('equivalence: %s apfc', c), max(m_other, m_own) <= TOL, [g1 g2]);
    % RPFC
    [ta, ~] = rpfc(s0, sopt('o')); [tb, ~] = rpfc(rep, sopt('r'));
    [rA, cA, xA] = pfc_tbl(ta);
    if shunt
        % the bus shunt S_10 and the self-term of the carrier 10-5 deliver at bus 5 what the bus shunt S_5 injects
        [m_other, g1] = pfc_map_diff(ta, tb, {}, {'S_5', {'S_10', 'E_10_5_s'}}, {'L_{f.10-5}', 'L_{t.10-5}'}, {});
        m_own = abs(sum(xA(strcmp(rA, 'L_{t.10-5}'), :)) - (rep.bus(5, BS) - s0.bus(5, BS)) * rep.bus(5, VM)^2); g2 = '';
    elseif isempty(eqv{k, 2})
        [m_other, g1] = pfc_map_diff(ta, tb, {}, {}, {'C2_5'}, {'G_5'});
        [m_own, g2] = pfc_map_diff(ta, tb, {'C2_5', 'G_5', 1}, {}, setdiff(rA, {'C2_5'}), setdiff(pfc_tbl_rows(tb), {'G_5'}));
    else
        % source: column of the dead-end branch against the column G_5 of the replacement generator; the dead-end
        % row holds the terminal value, which equals the generator output
        [m_other, g1] = pfc_map_diff(ta, tb, {}, {'G_5', {eqv{k, 3}}}, {eqv{k, 2}}, {});
        v = xA(strcmp(rA, eqv{k, 2}), strcmp(cA, eqv{k, 3}));
        m_own = abs(abs(v) - rep.gen(end, QG)); g2 = '';
    end
    fprintf('| %s | rpfc | %.17g | %.17g |\n', c, m_other, m_own);
    pfc_report(sprintf('equivalence: %s rpfc', c), max(m_other, m_own) <= TOL, [g1 g2]);
end
