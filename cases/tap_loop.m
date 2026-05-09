function mpc = tap_loop
% TAP_LOOP  Two-bus test case with transformer tap mismatch and reactive
%           circulating flow.
%
%   Synthetic network designed to produce a reactive circulating power
%   flow. Two parallel transformers connect bus 1 (110 kV slack) to bus 2
%   (10 kV) with intentionally mismatched tap ratios (1.05 and 0.95
%   pu). The mismatch drives around 6 Mvar around the loop 1-2-1.
%   The generator at bus 1 supplies 5 Mvar to the load and 1.264 Mvar reactive losses.
%
%   All branch shunts are zero (no line charging).
%   Active load is 10 MW at bus 2; reactive load is 5 Mvar.

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
    1	3	0	0	0	0	1	1	0	110	1	1.1	0.9;
    2	1	10	5	0	0	1	0.986683801	-1.54965234	10	1	1.1	0.9;
    ];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
    1	10.0689473	6.26403427	99	-99	1	100	1	99	0	0	0	0	0	0	0	0	0	0	0	0;
    ];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax	Pf	Qf	Pt	Qt
mpc.branch = [
    1	2	0.03	0.55	0	0	0	0	1.05	0	1	-360	360	4.2871	-6.1112	-4.2687	6.4492;
    1	2	0.03	0.55	0	0	0	0	0.95	0	1	-360	360	5.7818	12.3753	-5.7313	-11.4492;
    ];
