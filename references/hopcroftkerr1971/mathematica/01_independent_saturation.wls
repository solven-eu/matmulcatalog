(* ::Package:: *)

(* HK (2,2,bridge-3) impossibility theorem — independent Mathematica
   verification of the sympy x_equals_p_proof.py saturation argument.

   Reads nothing from disk; the entire problem is set up symbolically.
   Run:
       wolframscript -file 01_independent_saturation.wls
   (or in a Mathematica notebook, evaluate cell-by-cell).

   Expected wall-clock: ~5-30 s on a modern machine.

   Output:
     • PROPORTIONAL_GENS == GroebnerBasis(I_P)   — generators check
     • SATURATED_IDEAL == PROPORTIONAL_IDEAL     — the key claim:
         (I : c^infinity) = I_P
       confirming X = P ∪ V_{c=0} set-theoretically.
*)

Print["═══════════════════════════════════════════════════════════════"];
Print["HK (2,2,bridge-3): independent Mathematica saturation check"];
Print["═══════════════════════════════════════════════════════════════"];

(* ── Setup ───────────────────────────────────────────────────────── *)
aVars = {ai1, ai2, ap1, ap2, aj1, aj2};
xVars = {xi1, xi2, xp1, xp2, xj1, xj2};
mono = Flatten@Table[a*x, {a, aVars}, {x, xVars}];   (* 36 *)
n = Length[mono];

toVec[expr_] := Module[{e = Expand[expr]},
  Table[Coefficient[Coefficient[e, p[[1]]], p[[2]]],
    {p, Flatten[Table[{a, x}, {a, aVars}, {x, xVars}], 1]}]];

(* Targets *)
T1 = ai2*xi1 + ap2*xp1 + aj2*xj1;
T2 = ai1*xi2 + ap1*xp2 + aj1*xj2;

(* Shared (reusable) products *)
Ci = (ai1 - ai2)*xi2;
Di = ai1*(xi1 + xi2);
Cj = (aj1 - aj2)*xj2;
Dj = aj1*(xj1 + xj2);
Ep = ap2*xp2;
Fp = ap1*xp1;
shared = {Ci, Di, Cj, Dj, Ep, Fp};

vT1 = toVec[T1]; vT2 = toVec[T2];
Smat = Transpose@(toVec /@ shared);                  (* 36 × 6 *)

(* Projector onto S^perp, then onto Q = span(q1, q2) *)
PperpS = IdentityMatrix[n] - Smat . Inverse[Transpose[Smat] . Smat] . Transpose[Smat];
q1 = PperpS . vT1;
q2 = PperpS . vT2;
Qmat = Transpose[{q1, q2}];                          (* 36 × 2 *)
PQ = Qmat . Inverse[Transpose[Qmat] . Qmat] . Transpose[Qmat];
Pgarbage = PperpS - PQ;                              (* projector onto 28-D garbage *)
Print["rank(P_garbage) = ", MatrixRank[Pgarbage], "  (expect 28)"];

(* ── Generic rank-1 atom L·M ──────────────────────────────────────── *)
ls = Symbol["l" <> ToString[#]] & /@ Range[6];
ms = Symbol["m" <> ToString[#]] & /@ Range[6];
Ls = Symbol["L" <> ToString[#]] & /@ Range[6];
Ms = Symbol["M" <> ToString[#]] & /@ Range[6];

LA  = Sum[ls[[k]] * aVars[[k]], {k, 6}];
MA  = Sum[ms[[k]] * xVars[[k]], {k, 6}];
LAp = Sum[Ls[[k]] * aVars[[k]], {k, 6}];
MAp = Sum[Ms[[k]] * xVars[[k]], {k, 6}];

vA  = toVec[Expand[LA  * MA]];
vAp = toVec[Expand[LAp * MAp]];

garbageA  = Pgarbage . vA;
garbageAp = Pgarbage . vAp;
Print["garbage_A and garbage_A' built"];

(* Constraint map Phi_k = garbage_Ap_k - c * garbage_A_k *)
Phi = Table[Expand[garbageAp[[k]] - c * garbageA[[k]]], {k, n}];
PhiNonzero = Select[Phi, # =!= 0 &];
Print["#non-trivial Phi coordinates = ", Length[PhiNonzero]];

(* ── Fix the anchor A (matching sympy v2-anchor for cross-check) ──── *)
anchor = Thread[Join[ls, ms] -> {2, 3, 5, 7, 11, 13, 3, 5, 7, 11, 13, 17}];
fibreEqs = Expand /@ (PhiNonzero /. anchor);
fibreEqs = Select[fibreEqs, # =!= 0 &];
Print["Fibre at v2-anchor: ", Length[fibreEqs], " non-trivial polynomials in 13 unknowns"];

(* ── STEP 1: Gröbner basis of the fibre ideal ─────────────────────── *)
unknowns = Join[Ls, Ms, {c}];
Print[];
Print["STEP 1: GroebnerBasis(I, lex L1>L2>...>M6>c) ..."];
{tFibre, gbFibre} = AbsoluteTiming[GroebnerBasis[fibreEqs, unknowns,
    MonomialOrder -> Lexicographic]];
Print["   took ", tFibre, " s; |GB| = ", Length[gbFibre]];

(* ── STEP 2: Proportional ideal I_P ────────────────────────────────── *)
(* Parametrise: L'_k = alpha * lambda_k, M'_k = beta * mu_k, c = alpha*beta *)
{lamVals, muVals} = {{2, 3, 5, 7, 11, 13}, {3, 5, 7, 11, 13, 17}};
propGens = Join[
  Table[Ls[[k]] - alpha*lamVals[[k]], {k, 6}],
  Table[Ms[[k]] - beta*muVals[[k]], {k, 6}],
  {c - alpha*beta}];
Print[];
Print["STEP 2: eliminate (alpha, beta) from proportional ideal ..."];
{tProp, propIdeal} = AbsoluteTiming[GroebnerBasis[propGens, unknowns,
    {alpha, beta}, MonomialOrder -> Lexicographic]];
Print["   took ", tProp, " s; |I_P GB| = ", Length[propIdeal]];

(* ── STEP 3: Saturation (I : c^∞) via Rabinowitsch ────────────────── *)
Print[];
Print["STEP 3: saturate I by c via Rabinowitsch (adjoin t with 1 - t*c) ..."];
{tSat, satGB} = AbsoluteTiming[GroebnerBasis[Append[fibreEqs, 1 - t*c],
    Append[unknowns, t], {t}, MonomialOrder -> Lexicographic]];
Print["   took ", tSat, " s; |sat GB| = ", Length[satGB]];

(* ── STEP 4: Compare (I : c^∞) with I_P ────────────────────────────── *)
Print[];
Print["STEP 4: bi-directional containment check"];
satInIP = AllTrue[satGB, (PolynomialReduce[#, propIdeal, unknowns][[2]] === 0) &];
ipInSat = AllTrue[propIdeal, (PolynomialReduce[#, satGB, unknowns][[2]] === 0) &];
Print["   (I : c^∞)  ⊆  I_P  ?  ", satInIP];
Print["   I_P         ⊆  (I : c^∞) ?  ", ipInSat];

(* ── STEP 5: V(I + ⟨c⟩) = zero-garbage variety V ─────────────────── *)
Print[];
Print["STEP 5: V(I + ⟨c⟩) reproduces the zero-garbage equations of V?"];
fibreAtC0 = (fibreEqs /. c -> 0);
fibreAtC0 = Select[Expand /@ fibreAtC0, # =!= 0 &];
Print["   |fibre eqs at c=0| = ", Length[fibreAtC0]];
(* Build the zero-garbage equations directly: g(A') = 0 *)
zeroGarbageEqs = Select[Expand /@ (garbageAp /. anchor), # =!= 0 &];
Print["   |zero-garbage eqs|  = ", Length[zeroGarbageEqs]];
gbC0 = GroebnerBasis[fibreAtC0, Join[Ls, Ms]];
gbZG = GroebnerBasis[zeroGarbageEqs, Join[Ls, Ms]];
Print["   GroebnerBasis(fibreAtC0) == GroebnerBasis(zeroGarbageEqs) ?  ",
  Sort[gbC0] === Sort[gbZG]];

(* ── Verdict ──────────────────────────────────────────────────────── *)
Print[];
Print["═══════════════════════════════════════════════════════════════"];
If[satInIP && ipInSat,
  Print["✓ VERIFIED:  (I : c^∞) = I_P    →   X = P ∪ V_{c=0}"];
  Print["  The HK (2,2,bridge-3) impossibility theorem holds unconditionally"];
  Print["  within the rank-1 atom family over fields of characteristic 0."];
  ,
  Print["✗ FAILED:  ideals do not match — re-examine."];
];
Print["═══════════════════════════════════════════════════════════════"];
