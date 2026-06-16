(* ::Package:: *)

(* HK (2,2,bridge-3) impossibility theorem — does it extend to F₂?

   Re-runs the saturation argument of 01_independent_saturation.wls
   over the field of characteristic 2 (Modulus -> 2). The sympy proof
   works over characteristic 0; AlphaTensor 22-style F₂-escapes (e.g.
   ⟨4,4,4⟩ = 47/F₂ vs 48/C) make it worth checking whether F₂ admits
   a solution the characteristic-0 analysis missed.

   Run:
       wolframscript -file 02_characteristic_2_check.wls

   Expected wall-clock: ~10-60 s (GroebnerBasis with Modulus -> 2
   is sometimes faster, sometimes slower than over Q).

   Verdict:
     ✓  (I : c^∞) = I_P  over F₂      ⇒  impossibility extends to F₂
     ✗  saturated ideal strictly larger ⇒  candidate F₂ solutions exist
*)

Print["═══════════════════════════════════════════════════════════════"];
Print["HK (2,2,bridge-3): characteristic-2 saturation check"];
Print["═══════════════════════════════════════════════════════════════"];

(* ── Setup (duplicated from 01) ──────────────────────────────────── *)
aVars = {ai1, ai2, ap1, ap2, aj1, aj2};
xVars = {xi1, xi2, xp1, xp2, xj1, xj2};
mono = Flatten@Table[a*x, {a, aVars}, {x, xVars}];   (* 36 *)
n = Length[mono];

toVec[expr_] := Module[{e = Expand[expr]},
  Table[Coefficient[Coefficient[e, p[[1]]], p[[2]]],
    {p, Flatten[Table[{a, x}, {a, aVars}, {x, xVars}], 1]}]];

T1 = ai2*xi1 + ap2*xp1 + aj2*xj1;
T2 = ai1*xi2 + ap1*xp2 + aj1*xj2;

Ci = (ai1 - ai2)*xi2;
Di = ai1*(xi1 + xi2);
Cj = (aj1 - aj2)*xj2;
Dj = aj1*(xj1 + xj2);
Ep = ap2*xp2;
Fp = ap1*xp1;
shared = {Ci, Di, Cj, Dj, Ep, Fp};

vT1 = toVec[T1]; vT2 = toVec[T2];
Smat = Transpose@(toVec /@ shared);                  (* 36 × 6 *)

(* NOTE on projectors: over F₂ the rational-arithmetic projector
   (which inverts S^T·S in Q) is not well-defined. We build the
   projector over Q (rationals), use it to derive the SYMBOLIC garbage
   constraints, then reduce the resulting polynomials mod 2. The
   projector itself is a tool to extract the garbage-vanishing ideal;
   the polynomials are then handled over F₂. *)
PperpS = IdentityMatrix[n] - Smat . Inverse[Transpose[Smat] . Smat] . Transpose[Smat];
q1 = PperpS . vT1;
q2 = PperpS . vT2;
Qmat = Transpose[{q1, q2}];                          (* 36 × 2 *)
PQ = Qmat . Inverse[Transpose[Qmat] . Qmat] . Transpose[Qmat];
Pgarbage = PperpS - PQ;
Print["rank(P_garbage) over Q = ", MatrixRank[Pgarbage], "  (expect 28)"];

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

Phi = Table[Expand[garbageAp[[k]] - c * garbageA[[k]]], {k, n}];
PhiNonzero = Select[Phi, # =!= 0 &];
Print["#non-trivial Phi coordinates = ", Length[PhiNonzero]];

(* ── Fix the v2-anchor ────────────────────────────────────────────── *)
anchor = Thread[Join[ls, ms] -> {2, 3, 5, 7, 11, 13, 3, 5, 7, 11, 13, 17}];
fibreEqs = Expand /@ (PhiNonzero /. anchor);
fibreEqs = Select[fibreEqs, # =!= 0 &];

(* Clear rational denominators introduced by the projector so that
   Modulus -> 2 sees well-defined integer coefficients. *)
clearDen[p_] := Module[{lcm = LCM @@ Denominator /@ Cases[p, _Rational, {0, Infinity}]},
  If[lcm === Infinity || lcm === 0, p, Expand[lcm*p]]];
fibreEqsZ = clearDen /@ fibreEqs;
Print["Fibre at v2-anchor (Z-cleared): ", Length[fibreEqsZ], " polynomials"];

(* ── STEP 1: Gröbner basis over F₂ ────────────────────────────────── *)
unknowns = Join[Ls, Ms, {c}];
Print[];
Print["STEP 1: GroebnerBasis(I, Modulus -> 2) ..."];
{tFibre, gbFibre} = AbsoluteTiming[GroebnerBasis[fibreEqsZ, unknowns,
    Modulus -> 2, MonomialOrder -> Lexicographic]];
Print["   took ", tFibre, " s; |GB| = ", Length[gbFibre]];

(* ── STEP 2: Proportional ideal I_P over F₂ ──────────────────────── *)
{lamVals, muVals} = {{2, 3, 5, 7, 11, 13}, {3, 5, 7, 11, 13, 17}};
(* mod 2: lamVals -> {0,1,1,1,1,1}, muVals -> {1,1,1,1,1,1}.
   Use the SAME parametrisation so the comparison is meaningful. *)
propGens = Join[
  Table[Ls[[k]] - alpha*lamVals[[k]], {k, 6}],
  Table[Ms[[k]] - beta*muVals[[k]], {k, 6}],
  {c - alpha*beta}];
Print[];
Print["STEP 2: eliminate (alpha, beta) from I_P (Modulus -> 2) ..."];
{tProp, propIdeal} = AbsoluteTiming[GroebnerBasis[propGens, unknowns,
    {alpha, beta}, Modulus -> 2, MonomialOrder -> Lexicographic]];
Print["   took ", tProp, " s; |I_P GB| = ", Length[propIdeal]];

(* ── STEP 3: Saturation (I : c^∞) via Rabinowitsch, mod 2 ────────── *)
Print[];
Print["STEP 3: saturate I by c via Rabinowitsch (Modulus -> 2) ..."];
{tSat, satGB} = AbsoluteTiming[GroebnerBasis[Append[fibreEqsZ, 1 - t*c],
    Append[unknowns, t], {t}, Modulus -> 2,
    MonomialOrder -> Lexicographic]];
Print["   took ", tSat, " s; |sat GB| = ", Length[satGB]];

(* ── STEP 4: Compare (I : c^∞) with I_P over F₂ ──────────────────── *)
Print[];
Print["STEP 4: bi-directional containment check (mod 2)"];
reduceMod2[p_, basis_, vs_] := PolynomialReduce[p, basis, vs, Modulus -> 2][[2]];
satInIP = AllTrue[satGB, (reduceMod2[#, propIdeal, unknowns] === 0) &];
ipInSat = AllTrue[propIdeal, (reduceMod2[#, satGB, unknowns] === 0) &];
Print["   (I : c^∞)  ⊆  I_P  ?  ", satInIP];
Print["   I_P         ⊆  (I : c^∞) ?  ", ipInSat];

(* ── STEP 5: If they differ, surface the witness generators ──────── *)
If[!(satInIP && ipInSat),
  Print[];
  Print["DIVERGENCE — listing witness generators:"];
  extraInSat = Select[satGB,
    reduceMod2[#, propIdeal, unknowns] =!= 0 &];
  extraInIP = Select[propIdeal,
    reduceMod2[#, satGB, unknowns] =!= 0 &];
  Print["   sat-only generators (", Length[extraInSat], "):"];
  Do[Print["     ", g], {g, Take[extraInSat, UpTo[5]]}];
  Print["   I_P-only generators (", Length[extraInIP], "):"];
  Do[Print["     ", g], {g, Take[extraInIP, UpTo[5]]}];
];

(* ── Verdict ──────────────────────────────────────────────────────── *)
Print[];
Print["═══════════════════════════════════════════════════════════════"];
If[satInIP && ipInSat,
  Print["✓ VERIFIED over F₂:  (I : c^∞) = I_P"];
  Print["  HK (2,2,bridge-3) impossibility extends to characteristic 2."];
  Print["  No AlphaTensor-style F₂ escape exists for this fibre."];
  ,
  Print["✗ DIVERGENCE over F₂:  (I : c^∞) ≠ I_P"];
  Print["  The F₂ saturation strictly contains extra polynomials —"];
  Print["  candidate F₂-specific solutions may exist; investigate the"];
  Print["  witness generators listed above and search for explicit"];
  Print["  rank-1 atoms over F₂ that satisfy them."];
];
Print["═══════════════════════════════════════════════════════════════"];
