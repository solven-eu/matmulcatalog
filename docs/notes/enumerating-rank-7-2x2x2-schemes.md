# Enumerating rank-7 ⟨2,2,2⟩ schemes — three layers of coverage

A discussion-note (2026-06-01) on what "enumerate all rank-7 ⟨2,2,2⟩
schemes" means in our cost-minimisation context, and what tractable
layers of coverage exist.

## The framing question

For cost minimisation at a target like ⟨17,17,17⟩, we want the
*cheapest* rank-7 ⟨2,2,2⟩ outer scheme (at some allocation). The
question we're really asking: **have we exhausted the discrete orbit
of rank-7 ⟨2,2,2⟩ schemes that produce qualitatively different shape
multisets at unbalanced allocations?**

This note records three increasingly-thorough answers.

## What "equivalent" means in our cost perspective

Two rank-7 ⟨2,2,2⟩ schemes are **cost-equivalent for our purposes** if
they produce the same per-product shape multiset at every allocation.
This is *strictly stronger* than algebraic equivalence: de Groote
1978 says all rank-7 ⟨2,2,2⟩ schemes are GL-equivalent (single
9-dim orbit under GL₂(K)³), but most of that orbit collapses to a
small finite number of discrete shape-multiset equivalence classes.

**Strassen and its 8 axis-flip variants are NOT cost-equivalent at
unbalanced allocations**. At ⟨17,17,17⟩ with allocation (9,8)³:
- Strassen mask=0 (canonical): cost 2940
- Strassen mask=1..7: cost 2944 with 3 distinct multisets

The reason: axis-flip permutes which block each product touches. At
unequal block sizes, the `min(U_view, W_view)` reduction in
`Recombination.processAdditions` produces different effective shapes
per product depending on the relabeling. Strassen's full discrete
orbit at (9,8)³ is `{2940, 2944, 2944, 2944}` — four distinct
multisets, not one.

Same for Winograd: mask=0 gives 2940, mask=1 gives 2930 (the win),
other masks give other values. The 8 axis-flip variants of Winograd
span yet another discrete orbit, *disjoint* from Strassen's.

So **each canonical scheme's axis-flip orbit is a potentially distinct
discrete equivalence class**. Coverage of cost-relevant variants
requires enumerating the orbits of multiple canonical schemes.

## Layer 1 — discrete orbit of catalog schemes (~free)

**Approach**: take every rank-7 ⟨2,2,2⟩ scheme already in the catalog
(Strassen, Winograd 1971, FMM-Lille, AlphaTensor-Z, Probert-Fischer,
Pan 7-product variant if any, etc.). Apply the full cheap orbit
(`SymmetryTransforms.fullCheapOrbit` = S₃ × axis-flip, ≤ 48 variants
each). Score every (scheme, mask, allocation) at unbalanced cubic
targets via `AnalyticalMaskSearch`.

**Expected discoveries**:
- Confirm that Winograd's 2930 win at ⟨17,17,17⟩ is reproducible
  via Layer-1 search alone (we know it is — the cousin hunt found it).
- Surface analogous wins at ⟨19,19,19⟩, ⟨21,21,21⟩, etc. — if any
  catalog scheme's orbit produces a multiset that beats the current
  search prediction at these sizes.
- Identify schemes whose orbit DOES intersect with another's
  (cost-redundant) vs schemes whose orbit is unique.

**Cost**: implementation is trivial (sweep over `AnalyticalMaskSearch`
results). Running is sub-second per target shape.

**Limit**: Layer 1 covers what's in our catalog. If a new discrete
normal form exists that we never imported, Layer 1 misses it.

## Layer 2 — discrete normal-form classification

**Sources** (revised 2026-06-01 after literature audit):
- **de Groote 1978**: structural theorem; all rank-7 are GL-equivalent
- **Burichenko 2014** (arXiv:1408.6273): Strassen's discrete symmetry
  group has order 36 (72 extended). Useful for sizing the discrete
  orbit of a specific scheme, but not a normal-form enumeration.
- **Ikenmeyer-Moosbauer 2025** (arXiv:2503.05467, "Strassen's algorithm
  via orbit flip graphs"): uses an order-6 group action (S_3) on flip
  graphs. Most recent and most likely to give a constructive procedure.
  → Dig in via task #109.
- **Schwartz-Zwecher 2025 §2.20** (kin-row reduction): another angle
  on normal-form representation

(Earlier versions of this note cited "Heun 1994" and "Smirnov 2017"
as Layer-2 sources. Both turned out to be unverifiable —
see `heun-1994.md` for the negative finding. The actual literature
on rank-7 ⟨2,2,2⟩ classification is dominated by de Groote's
GL-uniqueness theorem; "discrete normal forms" only arise under
specific smaller groups and have no canonical published enumeration
beyond Burichenko's computation of Strassen's stabilizer.)

**Approach**: review these papers, extract the normal-form
representatives as Java/JSON schemes, register each in the catalog,
then run Layer 1 over the augmented catalog.

**Expected payoff**: if any of Heun/Smirnov's normal forms is NOT
GL-orbit-equivalent (at the discrete level) to our existing schemes,
it would unlock new cost minima at some target shapes. Whether this
yields new wins is empirical — most published catalogs likely
re-derive each other's schemes, but the existence of multiple
publication threads suggests at least some don't overlap.

**Cost**: 1-2 days. Paper-reading time dominates. Each normal form is
a small JSON to author.

**Limit**: Layer 2 covers what's published. Anything not in the
literature is missed.

## Layer 3 — full discrete SAT enumeration

**Approach**: encode the rank-7 ⟨2,2,2⟩ bilinear identity as a SAT/SMT
problem.

**Variables**:
- `u[i][j][k] ∈ {0,1}` for each block-position (i,j) ∈ {0,1}² × each
  product k ∈ {0..6}: "does product k's U-factor touch block (i,j)?"
- Similarly `v[j][l][k]`, `w[i][l][k]`
- (Optional) Coefficient sign `±1` per nonzero support

**Constraints**:
- Bilinear identity holds: for each output position (i, l, j),
  `Σ_k (u[i'][j'][k] · v[j'][l'][k] · w[i][l][k] · sign_factors) = δ_...`
- Symmetry-breaking lemmata to avoid counting orbit-equivalent
  solutions multiple times

**Output**: complete set of discrete support patterns that admit a
valid rank-7 ⟨2,2,2⟩ algorithm.

**Cost**: weeks. Encoding correctness, symmetry breaking, and verifying
solutions are all hard. Modern SAT solvers can likely handle the
search space at rank 7 / dim 2×2×2, but tuning is non-trivial.

**Limit**: even Layer 3 may not surface new shape-multiset
equivalence classes beyond what Layer 2 provides — most "novel"
solutions might be GL-discrete-orbit-equivalent to known ones. The
marginal value over Layer 2 is uncertain.

## What the literature says directly (informal recall)

- de Groote 1978: 9-dimensional GL orbit, single
  algebraic-equivalence class
- Heun 1994: finitely many discrete normal forms (the exact count
  needs paper verification)
- Smirnov 2017 / Smirnov 2013: published catalogs of explicit ⟨2,2,2⟩
  schemes — some are new normal forms, some re-derive prior work
- Sedoglavic 2017 (HAL): discusses the discrete vs continuous split
  in the orbit

## Recommended order of operations

1. **Layer 1 immediately** — sweep all catalog rank-7 ⟨2,2,2⟩ schemes
   × axis-flip orbit × allocation grid at n ∈ {17, 19, 21, 23, 25,
   27}. Identify new wins. ~1 hour.
2. **Layer 2 next** — review Heun 1994 + Smirnov 2017 specifically
   for ⟨2,2,2⟩=7 normal forms. Import any missing. Re-run Layer 1.
   ~2 days.
3. **Layer 3 only if Layer 2 leaves unexplained gaps** — encoding
   complexity makes this a research project, not routine work.

## Cross-refs

- `AnalyticalMaskSearch` (task #106) — the fast scoring primitive
  used by all three layers
- `TestAnalyticalMaskSearch` — cross-check vs brute-force
- `feedback_prefer_winograd_over_strassen.md` — discovered that the
  discrete orbit of Strassen does NOT contain the cheap support
  pattern that Winograd's discrete orbit contains
- task #105 — greedy shape-multiset search, the integrated version
  of Layer 1 in `BlockSplitSearch`
- task #102 — DisjointSumSearch (τ-theorem family) addresses the
  *cross-cardinality* gap (5-disjoint+TA pair vs 7-product), which
  is orthogonal to the rank-7 ⟨2,2,2⟩ classification question
