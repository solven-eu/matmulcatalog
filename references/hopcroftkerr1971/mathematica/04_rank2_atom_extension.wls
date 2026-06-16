(* ::Package:: *)

(* HK (2,2,bridge-3) — exploratory rank-2 atom extension.

   The impossibility theorem covers rank-1 atoms L·M. Here we set up
   the system allowing each "atom" to be a rank-2 sum
       L₁·M₁ + L₂·M₂
   (still 1 "product slot" in HK's framework if interpreted loosely,
   though arguably 2 mults). 3 atoms × 2 (L+M) × 2 (rank) × 6 coeffs
   = 72 free parameters.

   Question: with 72 parameters per the 3-atom team, does the system
   admit a solution to (W·Q-components = I₂, W·garbage = 0)?

   We do NOT attempt a full symbolic decision (the system is too
   large). We do:
     • set up Q (2-D) and the 28-D garbage space, as in 01;
     • parametrise 3 rank-2 atoms with 72 unknowns;
     • compute the symbolic (α, β) projection and the garbage of each;
     • at a few random (Z-valued) anchor points for the 72 atom
       parameters, compute the 3×2 Q-projection matrix and the 3×28
       garbage matrix, and report:
         - rank of the garbage matrix
         - dimension of its left null space
         - whether the Q-projection on that null space spans 2-D Q
     • report the dimensions and any feasibility witness.

   Run:
       wolframscript -file 04_rank2_atom_extension.wls

   Expected wall-clock: 30 s – a few minutes.
*)

Print["═══════════════════════════════════════════════════════════════"];
Print["HK (2,2,bridge-3): rank-2 atom extension (exploratory)"];
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
Print["dim Q = 2, dim garbage = ", MatrixRank[Pgarbage], " (expect 28)"];

(* ── Rank-2 atom: build one atom from 24 parameters ───────────────── *)
(* Atom_k = L_{k,1} * M_{k,1} + L_{k,2} * M_{k,2}                       *)
(* Each L_{k,r} has 6 coeffs in a-vars; each M_{k,r} has 6 in x-vars.   *)
(* Total per atom: 6+6+6+6 = 24. Three atoms: 72.                       *)

(* Symbolic parameters: l[k,r,i], m[k,r,i] for k∈{1,2,3}, r∈{1,2}, i∈1..6 *)
makeAtomVec[k_] := Module[{L1, M1, L2, M2},
  L1 = Sum[l[k, 1, i] * aVars[[i]], {i, 6}];
  M1 = Sum[m[k, 1, i] * xVars[[i]], {i, 6}];
  L2 = Sum[l[k, 2, i] * aVars[[i]], {i, 6}];
  M2 = Sum[m[k, 2, i] * xVars[[i]], {i, 6}];
  toVec[Expand[L1*M1 + L2*M2]]];

Print[];
Print["Building 3 symbolic rank-2 atom vectors (each 36-D) ..."];
{tBuild, atomVecs} = AbsoluteTiming[Table[makeAtomVec[k], {k, 3}]];
Print["   took ", tBuild, " s"];

(* Symbolic Q-projection and garbage of each atom *)
Print[];
Print["Computing symbolic Q-projection α, β per atom ..."];
{tProj, projData} = AbsoluteTiming[Table[
  {Qmat\[Transpose] . (PQ . atomVecs[[k]]),    (* unused (PQ . v ∈ image, projecting back gives same) *)
   Inverse[Transpose[Qmat] . Qmat] . (Transpose[Qmat] . atomVecs[[k]]),  (* (α, β) coords *)
   Pgarbage . atomVecs[[k]]                    (* garbage, 36-D *)
   }, {k, 3}]];
Print["   took ", tProj, " s"];

(* ── Numeric Jacobian / rank check at random parameter points ────── *)
Print[];
Print["═══════════════════════════════════════════════════════════════"];
Print["NUMERIC TANGENT-RANK CHECK at random rank-2 atom anchors"];
Print["═══════════════════════════════════════════════════════════════"];

allParams = Flatten@Table[
  {l[k, r, i], m[k, r, i]}, {k, 3}, {r, 2}, {i, 6}];
Print["Total parameters in 3 rank-2 atoms: ", Length[allParams], "  (expect 72)"];

SeedRandom[20260530];
nTrials = 5;
Do[
  randVals = Thread[allParams -> RandomInteger[{-5, 5}, Length[allParams]]];

  (* Q-projection matrix: 3 atoms × 2 (α, β) → 3×2 *)
  Qproj = Table[projData[[k, 2]] /. randVals, {k, 3}];

  (* Garbage matrix: 3 atoms × 36 → 3 × 36 (only 28 effective dims) *)
  Gmat = Table[Pgarbage . atomVecs[[k]] /. randVals, {k, 3}];

  rGarbage = MatrixRank[Gmat];
  nullDim = 3 - rGarbage;       (* left-null-space dim *)
  rQ = MatrixRank[Qproj];

  Print[];
  Print["Trial ", trial, ":"];
  Print["   rank(garbage matrix, 3×36) = ", rGarbage,
    "   ⇒ left-null-space dim = ", nullDim];
  Print["   rank(Q-projection,  3×2)  = ", rQ, "   (need 2 to span Q)"];

  (* Feasibility: do the W ∈ left-null-space rows project onto a 2-D
     subspace of Q? *)
  If[nullDim >= 1,
    Wbasis = NullSpace[Transpose[Gmat]];   (* W rows in R^3 *)
    Wproj = Wbasis . Qproj;                (* nullDim × 2 *)
    rWProj = MatrixRank[Wproj];
    Print["   image of null space in Q: dim = ", rWProj];
    If[rWProj === 2,
      Print["   ⚠  FEASIBLE at this anchor — span both q1 and q2!"];
      ,
      Print["   ✓ infeasible at this anchor — null-space image is ",
        rWProj, "-D ⊂ Q, cannot span 2-D Q"];
    ];
    ,
    Print["   ✓ infeasible — left null space is trivial"];
  ];
  , {trial, nTrials}];

(* ── Verdict ─────────────────────────────────────────────────────── *)
Print[];
Print["═══════════════════════════════════════════════════════════════"];
Print["NOTES:"];
Print["  • This is exploratory; a single FEASIBLE anchor does NOT"];
Print["    prove existence (the system has free residuals)."];
Print["  • If feasibility holds at multiple random anchors, follow up"];
Print["    by solving the linear system W·Qproj = I₂ ∧ W·garbage = 0"];
Print["    symbolically and reading off explicit L, M coefficients."];
Print["  • If infeasibility holds at every random anchor, the rank-2"];
Print["    extension likely also fails — but a clean proof would need"];
Print["    a Gröbner-style analysis of the augmented variety."];
Print["═══════════════════════════════════════════════════════════════"];
