# Specification of the baseline: Achayuthakan, Dent, Bialek, Ongsakul (2010), active power

Source: C. Achayuthakan, C. J. Dent, J. W. Bialek and W. Ongsakul, "Electricity tracing in systems with and without
circulating flows: physical insights and mathematical proofs", IEEE Trans. Power Syst., vol. 25, no. 2, pp. 1078-1086,
May 2010 (the "paper" below; page numbers are the journal's). This file states the published method in our own words
and notation, for the independent baseline of brief T14h. It is restricted to active power at bus level.

## 1. Network data and flows

- **Nodes and branches** (p. 1079, Sec. II-A): a network of n_bus nodes and m branches; per node the nodal flow
  P_i, the nodal generation P_G,i and the nodal demand P_D,i; per branch the flow F_l. Every branch is oriented in the
  direction of its flow (incidence matrix B, eq. (1): +1 at the start, -1 at the end of each line; B = B_u + B_d,
  eq. (2)).
- **The paper's derivation is for a lossless network** (p. 1079, Sec. II, first paragraph: "the lossless downstream
  version of the methodology will be considered"). There one flow per line exists, and the nodal flow is equally the
  sum of the outflows or of the inflows (eq. (6)):
  P = P_D + F_d·1 (outflows) or P = P_G + F_dᵀ·1 (inflows),
  with F_d = -B_dᵀ diag(F) B_u (eq. (4)), the n_bus × n_bus matrix whose (i, j) element is the flow in line i-j
  towards node j.
- **Losses** (p. 1080, last paragraph of Sec. II-A): "The influence of transmission losses can be accounted for by
  considering gross and net network flows [2]. To obtain the net flows, the vector F used to calculate A_d in (9) must
  be replaced by the vector of flows at the receiving end of each line. To obtain the gross flows, the flow vector F
  used to calculate A_u in (9) must be replaced by the vector of flows at the sending end of each line." In words: the
  net variant uses receiving-end flows, the gross variant sending-end flows. Nothing else about losses is specified in
  the paper; it refers to [2] (Bialek 1996) for the gross and net versions.
- **Nodes with zero nodal flow** are removed from the digraph, so that diag(P)⁻¹ exists (p. 1080, after eq. (8)).

## 2. The distribution matrices and the linear system

- **Downstream distribution matrix** (eqs. (7), (9), (10), p. 1080): A_d = I + B_dᵀ diag(F) B_u (diag P)⁻¹, of order
  n_bus, with elements [A_d]_ii = 1 and [A_d]_ij = -|P_ij| / P_j for j in α_i^d (the nodes supplied directly from
  node i), 0 otherwise. Here P_ij is the flow in the line linking nodes i and j (eq. (10)). The system is
  P_D = A_d P (eq. (7)), so P = A_d⁻¹ P_D (eq. (13)).
- **Upstream distribution matrix** (eqs. (8), (9), (10)): A_u = I + B_uᵀ diag(F) B_d (diag P)⁻¹, with
  [A_u]_ij = -|P_ji| / P_j for j in α_i^u (the nodes supplying node i directly); P_G = A_u P, P = A_u⁻¹ P_G (eq. (13)).
- **D_d = I - A_d** (eq. (16)): [D_d]_ij is the share of P_j supplied directly from node i through line i-j
  (p. 1081). Its column sums are at most 1, and strictly below 1 in a column of a node with generation (p. 1081,
  p. 1084).
- **T_d = (A_d)⁻¹ = I + Σ_{n≥1} D_dⁿ** (eqs. (20), (24), p. 1083): [T_d]_ij is the share of the nodal flow at node j
  supplied from node i through all paths. The paper computes T_d by inverting A_d (the matrix power series is for
  insight only, p. 1083).
- **Dimension:** n_bus × n_bus (the network's buses; buses with zero nodal flow removed). The number of nonzeros of
  D_d equals the number of directed node pairs connected by a line, with parallel lines merged into one pair.

## 3. Generator contributions to loads (the output)

- **Eq. (14), p. 1080** (downstream version): the proportion of the nodal flow P_i coming from the local generator is
  P_G,i / P_i, and
  P_G,i = (P_G,i / P_i) Σ_k [A_d⁻¹]_ik P_D,k.
  The contribution of the generator at node i to the demand at node k is therefore
  c_ik = (P_G,i / P_i) [A_d⁻¹]_ik P_D,k,
  where [A_d⁻¹]_ik is the share of the nodal flow at node k supplied from node i (T_d), so [A_d⁻¹]_ik P_D,k is the
  part of the demand at k that comes from node i, and P_G,i / P_i the part of that which is the local generator.
- **Eq. (15)**: the same for the branch flows, P_ji = (P_ji / P_i) Σ_k [A_d⁻¹]_ik P_D,k for j in α_i^u (not needed
  here; the output of the baseline is at bus level: generators (per bus) to demands (per bus)).
- **Upstream version** (stated as symmetric, p. 1080 and p. 1081): with A_u, P = A_u⁻¹ P_G, and the share of the
  demand at node k in its nodal flow, P_D,k / P_k, gives the contribution of each generator to each demand. The paper
  does not write this formula out.
- **Output:** a matrix of contributions, generator buses × demand buses (bus level, not branch level).

## 4. Assumptions and invertibility

- **Assumptions** (p. 1084, end of Sec. III; p. 1085): (a) all column sums of D_d are at most 1; (b) for at least one
  column K the column sum is strictly below 1 (a node with generation); (c) all elements of D_d are nonnegative.
  These follow from the definition of D_d and the presence of at least one source.
- **Theorem (Sec. IV-A, eqs. (28)-(32), p. 1085):** under (a)-(c), every eigenvalue of D_d has magnitude below 1.
- **Corollaries (Sec. IV-B, IV-C):** the matrix power series (20) converges, and A_d = I - D_d is invertible (eqs.
  (33), (34)). This holds with and without circulating flows (cycles in the digraph of flows).
- **Non-invertible example (Sec. IV-D, eq. (35), Fig. 8(c)):** a pure circulating flow with no generation and no
  demand (all column sums equal to 1) gives a singular A_d.
- **Conditioning (Sec. IV-E):** A_d becomes ill-conditioned when the generations are small compared with the nodal
  flows (an eigenvalue of A_d close to 0).
- The paper does not treat negative generation or negative demand, nor a line whose two ends both take power from the
  network or both inject (a line "fed from both ends"); its flow digraph assumes one direction per line.

## 5. The two solution variants of the baseline (brief T14h, step 2)

- (a) **As published:** form A_d (dense, n_bus × n_bus), invert it explicitly, apply eq. (14).
- (b) **Sparse:** factorize the sparse A_d once (LU) and solve for the right-hand sides needed by eq. (14) (the columns
  of T_d at the demand buses, or the rows at the generator buses through the transpose), without forming the inverse.

## 6. Ambiguities

Each item says whether it changes the numbers on the seven comparison networks (case1, case9, case14, case30, case118,
case_ACTIVSg500, case_ACTIVSg2000). The network facts are in `reports/T14h_scripts/out/step1_facts.md`: every network
has losses (1.3 % to 4.9 % of the generation); case_ACTIVSg2000 has 2 branches fed from both ends (pattern (+,-)),
2 branches with one zero end, 539 parallel branch pairs and 9 buses with several generators; case118 has 7 and
case_ACTIVSg500 12 parallel pairs.

1. **Loss variant: net (receiving-end flows) or gross (sending-end flows).** *Changes the numbers in every network.*
   The paper gives both and does not say which one the baseline is. It also ties them to different matrices: net to
   A_d, gross to A_u (p. 1080: "net flows ... used to calculate A_d"; "gross flows ... used to calculate A_u").
2. **Downstream (A_d, eq. (14)) or upstream (A_u) version.** *Changes the numbers in lossy networks.* The paper derives
   eq. (14) for the downstream version and only states that the upstream version is symmetric. With losses, the two
   give different contributions; with item 1 they pair naturally as net-downstream and gross-upstream, but the paper
   does not state this pairing as the method.
3. **Nodal flow P under losses: inflows or outflows (eq. (6)).** *Changes the numbers in lossy networks.* In a lossy
   network the two sums of eq. (6) differ when F is one value per line. The natural choice for net flows is
   P = P_G + F_dᵀ·1 with receiving-end flows (the inflow of the node), and for gross flows P = P_D + F_d·1 with
   sending-end flows (the outflow); the paper does not say.
4. **Generation and demand in the lossy versions.** *Changes the numbers.* In Bialek's net and gross versions [2], the
   generations are reduced by allocated losses (net) or the demands increased (gross). The paper does not say whether
   eq. (14) uses the actual P_G, P_D or the netted or grossed ones; with actual values, the contributions to a demand
   do not sum exactly to the demand in one of the variants.
5. **Lines fed from both ends (pattern (+,-) of the active flows) and lines with one zero end.** *Changes the numbers
   in case_ACTIVSg2000 only.* The paper's digraph has one direction per line and no such lines.
6. **Parallel lines.** *Does not change the numbers* if their flows are added into one pair (i, j) (eq. (4) adds them:
   F_d = -B_dᵀ diag(F) B_u sums the lines between the same ordered pair). Opposite flows on parallel lines would give
   both (i, j) and (j, i); none occur in these networks except as part of item 5.
7. **Several generators at one bus.** *Does not change the bus-level numbers.* The paper works with nodal generation
   P_G,i; the APFC table has one column per generator, so the comparison must sum APFC's generator columns per bus.
8. **Negative generation or demand.** *Does not occur* in the seven networks (step-1 facts), so it changes nothing here.

**Decision needed (stop condition 1 of brief T14h):** items 1 to 4 are number-changing in every network and are not
determined by the paper. A proposal for the decision: implement the paper's two loss variants as stated on p. 1080,
both with eq. (14) and its upstream counterpart, and report both:
- **net-downstream:** A_d from receiving-end flows, P = nodal inflow (P_G + F_dᵀ·1), actual P_G and P_D, eq. (14);
- **gross-upstream:** A_u from sending-end flows, P = nodal outflow (P_D + F_d·1), actual P_G and P_D, the upstream
  counterpart of eq. (14) (contribution of generator i to demand k = (P_D,k / P_k) [A_u⁻¹]_ki P_G,i);
and treat a line fed from both ends (item 5) as two flows out of the line, one into each end (as APFC's terminal
demands), stating that the paper does not cover it.
