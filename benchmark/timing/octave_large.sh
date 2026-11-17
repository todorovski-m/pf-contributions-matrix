#!/bin/zsh
# Run apfc.m or rpfc.m (this working tree) under Octave on the exported flows of large cases
# (benchmark/inputs_large/<case>.mat), through benchmark/timing/octave_exported.m: the solved struct is built from the
# exported data with PT and QT restored to the MATPOWER sign convention and run with opts.solved = true. One Octave
# process per case and kind, with a limit of LIMIT seconds (default 1200; perl alarm, since macOS has no timeout
# command). Writes <out dir>/<case>_<kind>.csv (writetable of tests/octave_shims, %.17g) and prints one status line
# per run. Compare with the port by reports/T13b_scripts/accept_octave.py or any CSV comparison.
# Usage (from the root of the matrix repository): zsh benchmark/timing/octave_large.sh <out dir> <kind>:<case> ...
setopt NO_UNSET
OUT=$1; shift
LIMIT=${LIMIT:-1200}
R=$PWD; T=$R/benchmark/timing
SCR=$(mktemp -d /tmp/octave_large.XXXXXX); mkdir -p $SCR/results $OUT
command cp -f apfc.m rpfc.m pfc_checks.m $SCR/
for kc in "$@"; do
  k=${kc%%:*}; c=${kc#*:}
  ( cd $SCR && perl -e "alarm $LIMIT; exec @ARGV" octave --no-gui --eval "addpath('$R/tests/octave_shims'); addpath('$T'); CNAME='$c'; KIND='$k'; INPUT='$R/benchmark/inputs_large/$c.mat'; octave_exported" ) > $SCR/$c.$k.log 2>&1
  rc=$?
  if (( rc == 0 )) && [[ -f $SCR/results/${c}_$k.csv ]]; then command cp -f $SCR/results/${c}_$k.csv $OUT/; print "$c $k: finished"
  elif (( rc == 142 )); then print "$c $k: not run to the end (time limit of $LIMIT s)"
  else print "$c $k: failed (exit $rc; $(grep -h 'error:' $SCR/$c.$k.log | head -1))"; fi
done
