(* ::Package:: *)

(* HK (2,2,bridge-3) impossibility — primary-decomposition cross-check.

   Sympy can't really do primary decomposition; Mathematica's Solve +
   GroebnerBasis machinery (plus the PolynomialIdeals` package, if
   available) is more capable. We re-derive the fibre ideal I at the
   v2 anchor and ask: into how many irreducible components does V(I)
   actually split?

   Expected (per sympy x_equals_p_proof.py):
     V(I)  =  V(I_P)  ∪  V(I + ⟨c⟩)  =  P-fibre  ∪  (V × {0})
   where V is the 8-component zero-garbage variety from
   zero_garbage_variety.py. So total irreducible components = 1 + 8 = 9
   (or 2 if we treat V as one reducible component).

   Run:
       wolframscript -file 03_primary_decomposition.wls

   Expected wall-clock: ~30 s – a few minutes. Solve over the lex GB
   of a 13-variable, ~30-equation system can be slow.
*)

Print["═══════════════════════════════════════════════════════════════"];
Print["HK (2,2,bridge-3): primary decomposition of V(I)"];
Print["═══════════════════════════════════════════════════════════════"];

(* ── Setup (duplicated from 01) ──────────────────────────────────── *)
aVars = {ai1, ai2, ap1, ap2, aj1, aj2};
xVars = {xi1, xi2, xp1, xp2, xj1, xj2};
mono = Flatten@Table[a*x, {a, aVars}, {x, xVars}];
n = Length[mono];

toVec[expr_] := Module[{e = Expand[expr]},
  Table[Coefficient[Coefficient[e, p[[1]]], p[[2]]],
    {p, Flatten[Table[{a, x}, {a, aVars}, {x, xVars}], 1]}]];

T1 = ai2*xi1 + ap2*xp1 + aj2*xj1;
T2 = ai1*xi2 + ap1*xp2 + aj1*xj2;

Ci = (ai1 - ai2)*xi2;  Di = ai1*(xi1 + xi2);
Cj = (aj1 - aj2)*xj2;  Dj = aj1*(xj1 + xj2);
Ep = ap2*xp2;          Fp = ap1*xp1;
shared = {Ci, Di, Cj, Dj, Ep, Fp};

vT1 = toVec[T1]; vT2 = toVec[T2];
Smat = Transpose@(toVec /@ shared);
PperpS = IdentityMatrix[n] - Smat . Inverse[Transpose[Smat] . Smat] . Transpose[Smat];
q1 = PperpS . vT1; q2 = PperpS . vT2;
Qmat = Transpose[{q1, q2}];
PQ = Qmat . Inverse[Transpose[Qmat] . Qmat] . Transpose[Qmat];
Pgarbage = PperpS - PQ;

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

Phi = Table[Expand[garbageAp[[k]] - c * garbageA[[k]]], {k, n}];
PhiNonzero = Select[Phi, # =!= 0 &];

anchor = Thread[Join[ls, ms] -> {2, 3, 5, 7, 11, 13, 3, 5, 7, 11, 13, 17}];
fibreEqs = Expand /@ (PhiNonzero /. anchor);
fibreEqs = Select[fibreEqs, # =!= 0 &];
unknowns = Join[Ls, Ms, {c}];
Print["Fibre ideal I has ", Length[fibreEqs], " generators in ",
  Length[unknowns], " unknowns"];

(* ── STEP 1: Try the PolynomialIdeals` package if available ───────── *)
Print[];
Print["STEP 1: Attempting PolynomialIdeals`PrimaryDecomposition ..."];
hasPD = Quiet[Check[Needs["PolynomialIdeals`"]; True, False]];
If[hasPD && ValueQ[PolynomialIdeals`PrimaryDecomposition],
  Print["   PolynomialIdeals` available — running PrimaryDecomposition"];
  {tPD, pd} = AbsoluteTiming[
    PolynomialIdeals`PrimaryDecomposition[
      PolynomialIdeals`PolynomialIdeal[fibreEqs, unknowns]]];
  Print["   took ", tPD, " s; |components| = ", Length[pd]];
  Do[Print["   component ", k, ": ",
      PolynomialIdeals`Dimension[pd[[k]]], "-dim, ",
      Length[PolynomialIdeals`GeneratingSet[pd[[k]]]], " gens"],
    {k, Length[pd]}];
  ,
  Print["   PolynomialIdeals`PrimaryDecomposition not available — ",
    "skipping; falling back to Solve."];
];

(* ── STEP 2: Solve the fibre system explicitly ───────────────────── *)
Print[];
Print["STEP 2: Solve[fibreEqs == 0, unknowns] (this may take a while) ..."];
solveEqs = Thread[fibreEqs == 0];
{tSolve, sols} = AbsoluteTiming[
  TimeConstrained[Solve[solveEqs, unknowns], 600,
    {"TIMEOUT"}]];
If[sols === {"TIMEOUT"},
  Print["   Solve timed out at 600 s — falling back to Reduce."];
  {tRed, sols} = AbsoluteTiming[
    TimeConstrained[Reduce[solveEqs, unknowns], 600,
      {"TIMEOUT"}]];
  Print["   Reduce took ", tRed, " s"];
  ,
  Print["   Solve took ", tSolve, " s; ", Length[sols], " branches"];
];

(* ── STEP 3: Classify branches ───────────────────────────────────── *)
Print[];
Print["STEP 3: Classify branches as proportional / zero-garbage / other"];

(* Proportional locus P: L' = alpha * lambda, M' = beta * mu,
   c = alpha * beta. With anchor lambda = (2,3,5,7,11,13),
   mu = (3,5,7,11,13,17), this means Ls and Ms are each a scalar
   multiple of the anchor vector. *)
{lamVals, muVals} = {{2, 3, 5, 7, 11, 13}, {3, 5, 7, 11, 13, 17}};
isProportional[sol_] := Module[{Lvals, Mvals, cVal, aFit, bFit},
  Lvals = Ls /. sol;
  Mvals = Ms /. sol;
  cVal = c /. sol;
  (* Check: does there exist alpha such that Lvals == alpha*lamVals,
     and beta such that Mvals == beta*muVals, and cVal == alpha*beta? *)
  Quiet[
    aFit = Solve[Lvals == alpha*lamVals, alpha];
    bFit = Solve[Mvals == beta*muVals, beta];
    Length[aFit] >= 1 && Length[bFit] >= 1 &&
      Simplify[(cVal /. sol) - (alpha /. aFit[[1]])*(beta /. bFit[[1]])] === 0
  ]];

isZeroGarbage[sol_] := Module[{},
  Quiet[Simplify[(c /. sol)] === 0]];

If[ListQ[sols] && Length[sols] > 0 && Head[sols[[1]]] === List,
  classified = Table[
    {k, isProportional[sols[[k]]], isZeroGarbage[sols[[k]]]},
    {k, Length[sols]}];
  Print["   Total branches:        ", Length[sols]];
  Print["   On proportional locus: ",
    Count[classified, {_, True, _}]];
  Print["   On c=0 (zero-garbage): ",
    Count[classified, {_, _, True}]];
  Print["   Other (unexpected):    ",
    Count[classified, {_, False, False}]];
  others = Select[classified, #[[2]] === False && #[[3]] === False &];
  If[Length[others] > 0,
    Print["   ⚠  unexpected branches at indices: ",
      others[[All, 1]]];
    Do[Print["     branch ", o[[1]], ": ", sols[[o[[1]]]]],
      {o, Take[others, UpTo[3]]}];
  ];
  ,
  Print["   sols is not a list-of-rules — Reduce output structure:"];
  Print["   ", Short[sols, 5]];
];

(* ── STEP 4: Cross-check expected count ──────────────────────────── *)
Print[];
Print["═══════════════════════════════════════════════════════════════"];
Print["Expected (from sympy):"];
Print["  V(I) = P̄ ∪ X_{c=0}"];
Print["  P̄  = proportional locus (1 irreducible component, dim 14)"];
Print["  X_{c=0} = R^12 × V × {0}, where V is the zero-garbage variety"];
Print["          with 8 irreducible components"];
Print["  ⇒  total irreducibles = 1 (P̄) + 8 (V slices) = 9"];
Print["     OR 2 if V is treated as one reducible component."];
Print["═══════════════════════════════════════════════════════════════"];
If[hasPD && ValueQ[PolynomialIdeals`PrimaryDecomposition] &&
    MatchQ[pd, _List] && (Length[pd] === 9 || Length[pd] === 2),
  Print["✓ Mathematica primary decomposition matches sympy prediction"];
  ,
  Print["⊙ Inspect the branch classification above; if all branches are"];
  Print["  proportional or zero-garbage, sympy's decomposition holds."];
];
Print["═══════════════════════════════════════════════════════════════"];
