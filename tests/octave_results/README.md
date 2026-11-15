# Octave results for the Matlab check

The outputs of `apfc.m` and `rpfc.m` for `case9`, `case118` and `case24_ieee_rts`, computed under
Octave 11.1.0 with MATPOWER 7.1 and the stand-ins of `tests/octave_shims/` (which write every value
with `%.17g`). They were produced by

    zsh tests/run_matrix.sh tests/octave_results case9 case118 case24_ieee_rts

from the code of the commit that added this folder. `tests/matlab_check.m` compares the Matlab
outputs against these files.
