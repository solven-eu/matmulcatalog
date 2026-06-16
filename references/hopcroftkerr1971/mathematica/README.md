# Mathematica cross-CAS verification scaffold for HK (2,2,bridge-3)

This directory contains Wolfram Language scripts that re-derive and
extend the sympy results in
[`../sympy/`](../sympy/) and the consolidated writeup in
[`../README.md`](../README.md).

The intent is **cross-CAS sanity**: if Mathematica disagrees with
sympy on any of the steps below, that's a major red flag — the
impossibility theorem rests on Gröbner-basis computations, and an
independent CAS reproducing the result is the strongest
non-paper-and-pencil check we have.

## The scripts

| Script | Purpose |
|--------|---------|
| [`01_independent_saturation.wls`](01_independent_saturation.wls) | Re-prove `(I : c^∞) = I_P` over ℚ (the core saturation step from `sympy/x_equals_p_proof.py`). Establishes `X = P̄ ∪ V_{c=0}` set-theoretically. |
| [`02_characteristic_2_check.wls`](02_characteristic_2_check.wls) | Re-run the saturation argument with `Modulus -> 2`. Does the impossibility extend to F₂, or does an AlphaTensor-style F₂ escape exist? |
| [`03_primary_decomposition.wls`](03_primary_decomposition.wls) | Primary-decomposition cross-check. Uses `PolynomialIdeals`PrimaryDecomposition` (if available) or falls back to `Solve`/`Reduce`. Expected: 1 (P̄) + 8 (zero-garbage `V` slices) = 9 irreducible components. |
| [`04_rank2_atom_extension.wls`](04_rank2_atom_extension.wls) | Exploratory: extend each atom from rank 1 to rank 2 (sum of two rank-1 tensors, 24 params per atom × 3 atoms = 72 total). Numeric tangent-rank check at several random anchors; report whether the rank-2 system admits feasible weights. |
| [`05_smoke_check_closed_cases.wls`](05_smoke_check_closed_cases.wls) | Cross-CAS smoke check of the three already-closed cases `(1,1,bridge-2)`, `(2,2,bridge-1)`, `(1,1,bridge-3)` against their Java emitters (`HopcroftKerr2bc.emitSameMethodPair_*`). Symbolic verification that the 3-atom solution + reusables reproduce `T1, T2`. |

## How to run

```bash
wolframscript -file 01_independent_saturation.wls
wolframscript -file 02_characteristic_2_check.wls
wolframscript -file 03_primary_decomposition.wls
wolframscript -file 04_rank2_atom_extension.wls
wolframscript -file 05_smoke_check_closed_cases.wls
```

Each script is **self-contained**: no `Get[]` of any other script,
no project imports. The setup block (12 variables, 36 monomials,
projectors onto S⊥ and Q, garbage projector) is duplicated across
scripts to keep them independent.

## Expected wall-clock

| Script | Typical runtime |
|--------|----------------:|
| 01 | 5–30 s |
| 02 | 10–60 s (`GroebnerBasis[…, Modulus -> 2]` can be erratic) |
| 03 | 30 s – a few minutes (`Solve`/`Reduce` over 13-var lex GB) |
| 04 | 30 s – a few minutes (5 random anchors, symbolic projection per atom) |
| 05 | < 5 s total (all 3 cases) |

## Dependencies

- **Mathematica 12+** (or any Wolfram Engine version supporting
  `GroebnerBasis[…, Modulus -> 2]`, `MonomialOrder -> Lexicographic`,
  and `PolynomialReduce[…, Modulus -> 2]`). Tested intent: 13.x.
- Script 03 will try to load the (built-in) `PolynomialIdeals`
  package; if it is unavailable in your distribution, the script
  prints a notice and falls back to `Solve`/`Reduce`.
- No external Mathematica packages required beyond the standard
  distribution.

## Cross-CAS verification policy

These scripts are deliberately written **without consulting the
sympy outputs** beyond the variable-naming conventions and the
v2-anchor `(2, 3, 5, 7, 11, 13, 3, 5, 7, 11, 13, 17)`. The Gröbner
bases and saturations are computed afresh, and the verdicts are
echoed at the end of each script with `✓` / `✗`.

**If Mathematica disagrees with sympy on any of:**

- `(I : c^∞) = I_P` (script 01) — the central saturation argument;
- the same equality holding over F₂ (script 02) — the
  characteristic-2 extension;
- the count of irreducible components of `V(I)` (script 03);
- the rank-1 closed-case formulas (script 05);

then the sympy result needs to be revisited. The current state
(2026-05) is that sympy passes all of these; a re-run in
Mathematica is meant to either reinforce confidence or surface a
genuine bug.

Script 04 is **exploratory** and not a cross-check: it asks a new
question (does the rank-2 atom family admit a closure?) that sympy
has not been asked. Its output is informational only — a clean
proof in either direction would require follow-up work.
