# Peeling zeros: input-side vs output-side reduction

Two distinct sources of "free" multiplication savings when we pad a
matmul target and run a fast scheme on the padded version.

## The setup

Suppose we want to multiply two `n × n` matrices, but the fast scheme
we'd like to use (Strassen ⟨2,2,2⟩=7, Laderman ⟨3,3,3⟩=23, …) has a
block structure that doesn't divide `n`. The standard trick is to
**pad** the inputs to a size that the block structure does divide
(e.g. `17 → 18 = 2·9` for Strassen), running the fast scheme on the
padded matrices, and discarding the padding from the output.

After padding:

- Some **input** rows / columns are entirely zero (the padding).
- Correspondingly, some **output** rows / columns are entirely zero
  (they receive only contributions from the zero input padding).
- We throw away the zero output entries at the end to recover the
  original-size result.

A "sub-product" in this context is one of the rank-1 atoms of the fast
scheme — for Strassen, one of M₁..M₇ — applied to some block layout of
the padded inputs.

## Two distinct optimizations

### Input-side reduction (already implemented)

> **"This sub-product's inputs are zero in some rows/columns, so we
> can use a smaller-shape scheme for it."**

If sub-product `M_k = (linear combo of A blocks) · (linear combo of B
blocks)` reads A only in its first `n_A` rows (because all other rows
are padding zeros), and B only in its first `n_B` columns for similar
reasons, we can replace the nominal ⟨n, n, n⟩ scheme with a smaller
⟨n_A, n, n_B⟩ scheme. The output is then padded back with the
appropriate zeros before being added into C.

Implemented in [`Recombination.processAdditions`](../src/main/java/eu/solven/matmul/catalog/Recombination.java).
Look for `subA = min(uu[0], ww[0])` — `uu[0]` is the input-A extent,
`ww[0]` is the input-A-side extent of the sub-product.

### Output-side reduction (IMPLEMENTED — see "What's implemented" below)

> **Status (2026-06)**: this section originally described the output-side
> case as missing. It is now implemented — the peel-aware
> `Recombination` overload shrinks each sub-product by the elementwise-min
> of input- and output-side effective extents, and `BlockSplitSearch`
> enumerates the `(alloc, peel)` patterns. The canonical γ5 case
> (⟨3,3,3⟩ via padded ⟨4,4,4⟩) is closed and tested. The genuine
> *leftover* is the ⟨17,17,17⟩-class odd-cubic gap, which is **not** a
> peel — see the end of this file. The explanation below is kept for the
> intuition.

> **"This sub-product's outputs go to positions we'll throw away
> anyway, so we can use a smaller-shape scheme for it."**

The dual case, equally real, not currently captured.

Islam (2009 MSc, Ch. 4) spells out the canonical example. When
Strassen computes a ⟨4,4,4⟩ recursion to land on a ⟨3,3,3⟩ target via
padding:

> *"γ5 = (Ã₂,₁ − Ã₁,₁)(B̃₁,₁ + B̃₁,₂). There is no zero row or column
> in any of the terms in this product. However, γ5 is used only to
> compute C̃₂,₂ = γ2 − γ3 + γ5 + γ7, and we know that C̃₂,₂ has only
> one non-zero term, so we only need one term in γ5: we reduce the
> cost from 7 to 2 multiplications."*

The input side of γ5 has no zeros (no help). The output side does:
C̃₂,₂ is the bottom-right block of the result, and after peeling the
padding back to ⟨3,3,3⟩, only one entry of C̃₂,₂ survives. So most
of γ5 is wasted work — Strassen's ⟨2,2,2⟩=7 collapses to a single
multiplication for the surviving entry.

Same idea generalises: for each sub-product, walk its W matrix to
find which output entries survive the post-peel, then shrink the
sub-product's nominal shape to match.

## Why it matters numerically

The motivating case at the time of writing is ⟨17,17,17⟩ (`R(⟨17⟩)`
over Q):

- **Our pipeline (input-side only)**: pad to ⟨18,18,18⟩, Strassen
  ⟨2,2,2⟩=7 on `[9,8]³` allocation, sum of 7 sub-products of mixed
  ⟨9,9,9⟩, ⟨8,8,9⟩-permutations and ⟨8,9,9⟩-permutations →
  **rank 2940**.
- **FMM-Lille's reported best**: same outer Strassen, but with
  Islam-style output reduction shrinking 2 of the ⟨9,9,9⟩
  sub-products to ⟨8,8,8⟩ (each saves 150 mults via the output-mask
  trick) and TA-pair-fusing the remaining diagonal ⟨9,9,9⟩ pair →
  **rank 2934**.

Six multiplications, all from one missing trick. Same pattern likely
explains the small +Δ at every odd cubic in our FMM gap report
(⟨21⟩, ⟨23⟩, ⟨25⟩, ⟨27⟩, ⟨29⟩, ⟨31⟩) — they all rely on padding to
the next-even and would all benefit from the output-side reduction.

## What's implemented (task #86 — peel scope COMPLETE; residual re-scoped below)

`Recombination.recombineWithAllocation(base, sota, allocA, allocB, allocC, peelA, peelB, peelC)`
is the peel-aware overload. It computes effective per-block sizes
`effA[i] = allocA[i] − peelA[i]` (and similar for B, C) and passes those
to `processAdditions`. Each sub-product's resulting shape is the
elementwise-min of input-side and output-side effective extents — i.e.
both directions of the Islam reduction in one pass.

**What this closes**: the canonical Islam γ5 case for ⟨3,3,3⟩ via
padded ⟨4,4,4⟩ Strassen: the sub-product whose W column lives entirely
in the peeled corner block C̃₂,₂ collapses from ⟨2,2,2⟩=7 to ⟨1,2,1⟩=2
multiplications. Test: `TestRecombinationOutputPeel.canonical_islam_gamma5_case_3x3x3_via_padded_strassen`.

**What this DOES NOT close**: the ⟨17,17,17⟩=2934 case. The
`[9,9]³+peel=[0,1]³` invocation is mathematically equivalent to direct
`[9,8]³` because none of the 7 standard Strassen products has its W
entirely inside the peeled (1,1)-block. FMM-Lille's 2934 recipe uses
a *different* Strassen product↔block mapping (with 2 ⟨9,9,9⟩'s on the
diagonal pair-fused + 2 ⟨8,8,8⟩'s in the corner). That's an outer-
template substitution, not a peel — see follow-up task on
non-permutation invertible S_{X,Y,Z} transforms.

## Search-layer enumeration (IMPLEMENTED)

The caller no longer has to pass peel arrays by hand: the
`maxPadding` variant of `BlockSplitSearch.findBestMultiBaseSplit`
enumerates over-allocations summing to `target+1 .. target+maxPadding`
and peels the excess off the last block via the peel-aware
`Recombination` overload (see the `maxPadding` Javadoc in
`BlockSplitSearch`). `maxPadding = base.n − 1` captures the dominant
pattern; the search-space cost is `≤ (maxPadding+1)^3` per base.

## Implementation sketch (task #86)

For each rank-`r` sub-product emitted by `Recombination.recombineWithAllocation`:

1. Look at the W column for this product (the `n·p`-shaped vector of
   per-output-entry coefficients).
2. Combine with the peel mask `peelN[i] ∈ {0, 1}` and `peelP[j] ∈ {0, 1}`
   (which output entries survive the unpad).
3. The effective output extent of this sub-product is the bounding
   box of `{(i, j) : W[i, j] ≠ 0 AND peelN[i] = 1 AND peelP[j] = 1}`.
4. If that box is strictly smaller than the nominal sub-shape,
   recurse with the smaller shape via `sota.getRank(...)`.

Cross-references: see
[`docs/diagnostics/17x17x17_pair_fuse.md`](diagnostics/17x17x17_pair_fuse.md)
for the per-product cost accounting that shows the savings land at
exactly 2934, and
[`src/test/java/eu/solven/matmul/search/TestPairFusingDiagonal17.java`](../src/test/java/eu/solven/matmul/search/TestPairFusingDiagonal17.java)
test (8) for the rejected lead (DIS09 permutation-only S_{X,Y,Z} orbit
— tested, doesn't close the gap on its own).

## Residual (re-scoped out of task #86)

**Task #86's peel scope is complete** (both input- and output-side
zero-peel, search-enumerated via `BlockSplitSearch`'s `maxPadding`,
canonical γ5 closed and tested). What remains is **not a peel** and is
re-filed as a separate "pair-fusion + template substitution" item:

The odd-cubic gaps (`⟨17,17,17⟩=2934`, and the analogous +Δ at ⟨21⟩,
⟨23⟩, ⟨25⟩, ⟨27⟩, ⟨29⟩, ⟨31⟩) need two things peel cannot provide:

1. **Pair-fusion on *non-uniform* sub-shapes.** The Pan-TA primitive
   exists (`PairFusedRecombination`) and is wired into `BlockSplitSearch`
   (`ofPairFused` / `predictBalancedCubic`), but only for the **uniform
   cubic** case `⟨k,k,k⟩`. The ⟨17³⟩ optimum sits at the non-uniform
   `(9,8)³` allocation (mixed ⟨9,9,9⟩/⟨8,9,9⟩/⟨8,8,9⟩ sub-shapes), which
   that MVP explicitly does not cover. Generalising pair-fusion to
   non-uniform diagonals is step one.
2. **A non-permutation invertible `S_{X,Y,Z}` template substitution** to
   place ⟨8,8,8⟩ in the corner and expose the diagonal ⟨9,9,9⟩ pair for
   fusion. The *permutation-only* `S_{X,Y,Z}` orbit is already searched
   (and rejected — see `TestPairFusingDiagonal17` test 8); the missing
   piece is the non-permutation family.

Caveat: for ⟨17³⟩ *specifically*, the gap is already beaten by a
**different** mechanism — the Winograd-cousin axis-flip lands at
**2930 < 2934** (`solven_winograd_cousin_axflip1`, §axis-flip orbit).
So this residual matters more for the *general* odd-cubic family than
for ⟨17³⟩ alone, and it is a `KnownTauIdentities`/`DisjointSumSearch`
generalisation, not a peeling task.
