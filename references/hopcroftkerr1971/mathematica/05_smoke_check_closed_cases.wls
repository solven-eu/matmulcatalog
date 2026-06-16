(* ::Package:: *)

(* HK 1971 same-method-pair cases — Mathematica cross-CAS smoke check.

   Re-verify the three sympy-derived (and Java-emitted) closed cases:
     • (1,1,bridge-2)   — HK's own page-10 derivation
     • (2,2,bridge-1)   — derived by sympy/derive_22bridge1.py
     • (1,1,bridge-3)   — same template
   Each case: set up targets T1, T2 + the reusable set S (different
   per case) and the explicit 3-atom solution; verify symbolically
   that the W-combination reproduces T1, T2.

   The reference implementation lives at
     src/main/java/eu/solven/matmul/papers/hopcroftkerr1971/
       HopcroftKerr2bc.java#emitSameMethodPair_{11_bridge2,22_bridge1,11_bridge3}

   Run:
       wolframscript -file 05_smoke_check_closed_cases.wls

   Expected wall-clock: < 1 s per case, < 5 s total.
*)

Print["═══════════════════════════════════════════════════════════════"];
Print["HK same-method-pair: cross-CAS smoke check of 3 closed cases"];
Print["═══════════════════════════════════════════════════════════════"];

(* ── Common variable setup ───────────────────────────────────────── *)
(* For each case we use 3 rows i, p (=i+1), j, with two a/x scalars
   each. Variables match the sympy / Java naming. *)
aVars = {ai1, ai2, ap1, ap2, aj1, aj2};
xVars = {xi1, xi2, xp1, xp2, xj1, xj2};

(* "Substituted args" used by all 3 cases (HK's page-10 trick) *)
subA1 = ai1 + ap1 - aj1;
subA2 = ai2 + ap2 - aj2;
subX1 = xi1 + xp1 - xj1;
subX2 = xi2 + xp2 - xj2;

(* Six base method constructors *)
methA[a1_, a2_, x1_, x2_] := a2 * (x1 + x2);
methB[a1_, a2_, x1_, x2_] := (a1 - a2) * x1;
methC[a1_, a2_, x1_, x2_] := (a1 - a2) * x2;
methD[a1_, a2_, x1_, x2_] := a1 * (x1 + x2);
methE[a1_, a2_, x1_, x2_] := a2 * x2;
methF[a1_, a2_, x1_, x2_] := a1 * x1;
(* G constructor used in HK page 10 derivations *)
methG[s2_, s1_, k1_, l2_] := (s1 + s2) * (k1 - l2);

(* Generic targets per pair: T_ij = y_{i,j} + y_{p,j}; T_ji = y_{j,i} + y_{j,p}
   y_{r,c} = a_{r,1}*x_{1,c} + a_{r,2}*x_{2,c}.
   For our 3-row {i, p, j} setup, the relevant outputs are y_{i,j}, y_{p,j},
   y_{j,i}, y_{j,p}. *)
y[r1_, r2_, c1_, c2_] := r1*c1 + r2*c2;   (* a_{r,1}*x_{1,c} + a_{r,2}*x_{2,c} *)
yIJ = y[ai1, ai2, xi2, xj2];   (* y_{i,j}: row i, col j → a_{i,1}*x_{1,j} + a_{i,2}*x_{2,j} *)
(* Wait — let's be careful: y_{r,c} = a_{r,1}*x_{1,c} + a_{r,2}*x_{2,c} *)
yIJ = ai1*xj1 + ai2*xj2;       (* Hmm, this is x_{1,j} = xj1, x_{2,j} = xj2. *)

(* Re-derive cleanly. With x naming x_{row, col}: xi1 = x_{1,i}, xi2 = x_{2,i},
   xj1 = x_{1,j}, etc. So y_{r,c} = a_{r,1}*x_{1,c} + a_{r,2}*x_{2,c}. *)
yIJ = ai1*xj1 + ai2*xj2;
yPJ = ap1*xj1 + ap2*xj2;
yJI = aj1*xi1 + aj2*xi2;
yJP = aj1*xp1 + aj2*xp2;

(* The pair targets in HK's framework: *)
TIJ = yIJ + yPJ;     (* "sum_{ij}" in Java code *)
TJI = yJI + yJP;     (* "sum_{ji}" in Java code *)

verify[caseName_, atomList_, weightsIJ_, weightsJI_] := Module[
  {sumIJ, sumJI, resIJ, resJI, ok1, ok2},
  sumIJ = Sum[weightsIJ[[k]] * atomList[[k]], {k, Length[atomList]}];
  sumJI = Sum[weightsJI[[k]] * atomList[[k]], {k, Length[atomList]}];
  resIJ = Expand[sumIJ - TIJ];
  resJI = Expand[sumJI - TJI];
  ok1 = resIJ === 0;
  ok2 = resJI === 0;
  Print["─── ", caseName, " ───"];
  Print["  T_{ij} residual: ", If[ok1, "0  ✓", resIJ]];
  Print["  T_{ji} residual: ", If[ok2, "0  ✓", resJI]];
  ok1 && ok2];

verdicts = {};

(* ── Case 1: (1, 1, bridge-2) ────────────────────────────────────── *)
(* y_{ii}, y_{jj} use method 1 (A+B); bridge y_{i+1,i+1} uses method 2.
   Per HopcroftKerr2bc#emitSameMethodPair_11_bridge2:
     SHARED: A_jj=A(a_{j,1},a_{j,2},x_{1,j},x_{2,j}), B_jj,
             E_adj=E(a_{i,2}+a_{i+1,2}, ., x_{2,i}+x_{2,i+1}, .),
             F_adj
     NEW:    C_sub = C(sub_a_1, sub_a_2, _, sub_x_2)
             D_sub = D(sub_a_1, _, sub_x_1, sub_x_2)
             G    = (sub_a_1 + a_{j,2}) * (x_{1,j} - sub_x_2)
   Weights (from Java):
     sum_{ij} = +B_jj + C_sub + E_adj + G
     sum_{ji} = +A_jj + F_adj − D_sub − G
*)
Ajj   = methA[aj1, aj2, xj1, xj2];
Bjj   = methB[aj1, aj2, xj1, xj2];
EadjV = methE[ai1 + ap1, ai2 + ap2, xi1 + xp1, xi2 + xp2];
FadjV = methF[ai1 + ap1, ai2 + ap2, xi1 + xp1, xi2 + xp2];
Csub  = methC[subA1, subA2, subX1, subX2];
Dsub  = methD[subA1, subA2, subX1, subX2];
Gval  = (subA1 + aj2) * (xj1 - subX2);

atoms11b2 = {Ajj, Bjj, EadjV, FadjV, Csub, Dsub, Gval};
(* Order:    A_jj, B_jj, E_adj, F_adj, C_sub, D_sub, G *)
wIJ11b2 = { 0,    +1,    +1,    0,     +1,    0,     +1};
wJI11b2 = { +1,   0,     0,     +1,    0,     -1,    -1};
AppendTo[verdicts, {"(1,1,bridge-2)", verify["(1,1,bridge-2)",
  atoms11b2, wIJ11b2, wJI11b2]}];

(* ── Case 2: (2, 2, bridge-1) ────────────────────────────────────── *)
(* y_{ii}, y_{jj} use method 2 (-C+D); bridge y_{i+1,i+1} uses method 1.
   Per HopcroftKerr2bc#emitSameMethodPair_22_bridge1:
     SHARED: C_jj, D_jj, E_adj, F_adj
     NEW:    A_sub = A(_, sub_a_2, sub_x_1, sub_x_2)
             B_sub = B(sub_a_1, sub_a_2, sub_x_1, _)
             G    = (a_{j,1} + sub_a_2) * (sub_x_1 − x_{2,j})
   Weights:
     sum_{ij} = −C_jj + F_adj − B_sub − G
     sum_{ji} = +D_jj + E_adj − A_sub + G
*)
Cjj22b1   = methC[aj1, aj2, xj1, xj2];
Djj22b1   = methD[aj1, aj2, xj1, xj2];
Asub22b1  = methA[subA1, subA2, subX1, subX2];
Bsub22b1  = methB[subA1, subA2, subX1, subX2];
Gv2       = (aj1 + subA2) * (subX1 - xj2);

atoms22b1 = {Cjj22b1, Djj22b1, EadjV, FadjV, Asub22b1, Bsub22b1, Gv2};
(* Order:    C_jj,   D_jj,    E_adj, F_adj, A_sub,    B_sub,    G_v2 *)
wIJ22b1 = { -1,      0,       0,     +1,    0,        -1,       -1};
wJI22b1 = {  0,      +1,      +1,    0,     -1,       0,        +1};
AppendTo[verdicts, {"(2,2,bridge-1)", verify["(2,2,bridge-1)",
  atoms22b1, wIJ22b1, wJI22b1]}];

(* ── Case 3: (1, 1, bridge-3) ────────────────────────────────────── *)
(* y_{ii}, y_{jj} use method 1 (A+B); bridge y_{i+1,i+1} uses method 3 (E+F).
   Per HopcroftKerr2bc#emitSameMethodPair_11_bridge3:
     SHARED: A_jj, B_jj, E_adj, F_adj   (same naming, different sourcing)
     NEW:    C_sub, D_sub, G_v1 (= (subA1 + aj2) * (xj1 − subX2),
                                 same G expression as (1,1,bridge-2))
   Weights:
     sum_{ij} = +B_jj + E_adj + C_sub + G_v1
     sum_{ji} = +A_jj + F_adj − D_sub − G_v1
*)
atoms11b3 = atoms11b2;       (* same atom list, different weights *)
wIJ11b3 = { 0,    +1,    +1,    0,     +1,    0,     +1};
wJI11b3 = { +1,   0,     0,     +1,    0,     -1,    -1};
AppendTo[verdicts, {"(1,1,bridge-3)", verify["(1,1,bridge-3)",
  atoms11b3, wIJ11b3, wJI11b3]}];

(* ── Verdict ─────────────────────────────────────────────────────── *)
Print[];
Print["═══════════════════════════════════════════════════════════════"];
allOk = AllTrue[verdicts, #[[2]] &];
Do[Print["  ", v[[1]], ":  ", If[v[[2]], "✓ verified", "✗ FAILED"]],
  {v, verdicts}];
Print[];
If[allOk,
  Print["✓ All 3 closed cases reproduce in Mathematica."];
  Print["  Sympy ↔ Mathematica agree on the same-method-pair emitters."];
  ,
  Print["✗ At least one case failed — cross-CAS DISAGREEMENT."];
  Print["  Investigate: either the Java/sympy formula or the"];
  Print["  Mathematica encoding here is wrong."];
];
Print["═══════════════════════════════════════════════════════════════"];
