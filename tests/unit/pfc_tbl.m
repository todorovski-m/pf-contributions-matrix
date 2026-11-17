function [rows, names, data] = pfc_tbl(t)
%PFC_TBL  Row labels, column names and values of an output table of apfc or rpfc, for a Matlab table and for
%   the stand-in of tests/octave_shims (a struct with fields rows, names, data).
if isstruct(t)
    rows = t.rows(:); names = t.names(:)'; data = full(t.data);
else
    rows = t.Line; names = t.Properties.VariableNames(2:end); data = t{:, 2:end};
end
end
