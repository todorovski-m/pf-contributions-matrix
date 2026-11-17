function v = pfc_tbl_row(t, lab)
%PFC_TBL_ROW  Values of the row with label LAB of an output table (see PFC_TBL).
[r, ~, x] = pfc_tbl(t);
v = x(strcmp(r, lab), :);
end
