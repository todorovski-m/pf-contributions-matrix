function [mx, msg] = pfc_map_diff(tA, tB, rowpairs, colgroups, exclA, exclB)
%PFC_MAP_DIFF  Largest absolute difference between two output tables after mapping rows and columns by label.
%   ROWPAIRS: cell array {labelA, labelB, sign; ...} (row labelB of B compared with sign * row labelA of A);
%   every other row is compared with the row of the same label. COLGROUPS: {labelB, {labelsA}; ...} (column labelB
%   of B compared with the sum of the columns labelsA of A); every other column with the column of the same label.
%   EXCLA, EXCLB: row labels left out of the comparison. Every row and column that is not left out must be mapped;
%   MX is Inf and MSG names the first unmapped label otherwise.
[rA, cA, dA] = pfc_tbl(tA);
[rB, cB, dB] = pfc_tbl(tB);
if nargin < 3 || isempty(rowpairs), rowpairs = cell(0, 3); end
if nargin < 4 || isempty(colgroups), colgroups = cell(0, 2); end
if nargin < 5, exclA = {}; end
if nargin < 6, exclB = {}; end
mx = Inf; msg = '';
% columns: for every column of B the columns of A
colsA = cell(1, numel(cB)); usedA = false(1, numel(cA));
for j = 1:numel(cB)
    g = find(strcmp(colgroups(:, 1), cB{j}));
    if isempty(g), lab = cB(j); else, lab = colgroups{g, 2}; end
    [ok, idx] = ismember(lab, cA);
    if ~all(ok), msg = sprintf('column %s of B has no counterpart in A', cB{j}); return; end
    colsA{j} = idx; usedA(idx) = true;
end
if ~all(usedA), msg = sprintf('column %s of A is not mapped', cA{find(~usedA, 1)}); return; end
% rows
mx = 0; usedB = ismember(rB, exclB);
for i = 1:numel(rA)
    if any(strcmp(exclA, rA{i})), continue; end
    p = find(strcmp(rowpairs(:, 1), rA{i}));
    if isempty(p), lb = rA{i}; s = 1; else, lb = rowpairs{p, 2}; s = rowpairs{p, 3}; end
    ib = find(strcmp(rB, lb));
    if numel(ib) ~= 1, mx = Inf; msg = sprintf('row %s of A has no counterpart %s in B', rA{i}, lb); return; end
    usedB(ib) = true;
    for j = 1:numel(cB)
        mx = max(mx, abs(s * sum(dA(i, colsA{j})) - dB(ib, j)));
    end
end
if ~all(usedB), mx = Inf; msg = sprintf('row %s of B is not mapped', rB{find(~usedB, 1)}); end
end
