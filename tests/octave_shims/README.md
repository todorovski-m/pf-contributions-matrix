# Octave stand-ins for validation only

`array2table.m`, `addvars.m` and `writetable.m` are minimal stand-ins for the Matlab
functions of the same names, which are not available in Octave without the `datatypes`
package. They exist only so that `apfc.m` and `rpfc.m` can be validated under Octave:

- `array2table` keeps the data matrix and the `VariableNames`;
- `addvars` stores the row-label column;
- `writetable` writes a CSV file with a header `Line,<names>` and every value as `%.17g`.

They are not a general replacement for the Matlab table type. In Matlab, do not put this
folder on the path; the built-in functions are used, and `writetable` writes at most 15
significant digits.

Usage under Octave, from the repository root:

    octave --eval "addpath('tests/octave_shims'); apfc('cases/case9'); rpfc('cases/case9')"
