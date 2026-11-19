# Octave results for the Matlab check

The outputs of `apfc.m` and `rpfc.m` for `case9`, `case118` and `case24_ieee_rts`, computed under
Octave 11.1.0 with MATPOWER 7.1 and the stand-ins of `tests/octave_shims/` (which write every value
with `%.17g`). They were produced by

    zsh tests/run_matrix.sh tests/octave_results case9 case118 case24_ieee_rts

from the code of the commit that added this folder, and again for rule Q3 version 2.0 (brief T16; the rpfc files of
the tapped cases case118 and case24_ieee_rts change, case9_rpfc only at the level of rounding).
`tests/matlab_check.m` compares the Matlab outputs against these files.
