# Q⟨17,17,17⟩ search: 2940 vs FMM-Lille 2934 — diagnosis

**Status:** root cause identified; **not** a bug in `PairedSubProducts.applyPairing`. No fix landed.

**Repro:** `mvn -ntp test -Dtest=TestPairFusingDiagonal17`
(`src/test/java/eu/solven/matmul/search/TestPairFusingDiagonal17.java`)

## What FMM-Lille publishes

Per https://fmm.univ-lille.fr/17x17x17.html, the best known rank for
`Q⟨17,17,17⟩` is **2934**, with recipe:

```
⟨17,17,17⟩ = ⟨9,9,8⟩ + ⟨8,8,8⟩ + ⟨8,8,8⟩ + ⟨8,9,9⟩ + ⟨9,8,9⟩
             + TA(⟨9,9,9⟩, ⟨9,9,9⟩)
           = 430 + 336 + 336 + 430 + 430 + 972
           = 2934
```

The `TA(⟨9,9,9⟩,⟨9,9,9⟩) = 9³ + 3·9² = 972` term is Pan's pair-product
(simultaneous compute of two cyclically-equivalent sub-products).

## What our search returns

`BlockSplitSearch.findBestMultiBaseSplit(17,17,17, pool, sotaQ, balancedOnly)`
returns **rank 2940** via Strassen `⟨2,2,2⟩` recombination on `[9,8]³`.

Diagnostic trace (sub-product shapes + rank lookup from the catalog,
Q-field, via `FormulaAwareSota`):

| base mult | sub-shape  | rank |
|-----------|------------|------|
| M1        | ⟨9,9,9⟩   | 486  |
| M2        | ⟨8,9,9⟩   | 430  |
| M3        | ⟨9,9,8⟩   | 430  |
| M4        | ⟨8,8,9⟩   | 388  |
| M5        | ⟨9,8,8⟩   | 388  |
| M6        | ⟨8,9,8⟩   | 388  |
| M7        | ⟨9,8,9⟩   | 430  |
| **total** |            | **2940** |

`applyPairing` on these shapes returns **2940** (zero savings) because the
only profitable cyclic pairings would be between cyclically-equivalent
shapes, and:

```
pairCost(8,9,9) = 8·9·9 + 8·9 + 9·9 + 9·8 = 873; 2·R(⟨8,9,9⟩) = 860 → loss 13
pairCost(8,8,9) = 8·8·9 + 8·8 + 8·9 + 9·8 = 784; 2·R(⟨8,8,9⟩) = 776 → loss 8
```

So pair-fusing the (5,5) cyclic candidates would make things *worse*. And
there's only ONE ⟨9,9,9⟩ sub-product (M1) — no pair to fuse.

## Why our recipe differs from FMM's

The two recipes disagree on the shape multiset of the 7 Strassen sub-products:

- **Ours (`[9,8]³` Strassen):** one ⟨9,9,9⟩ + four perms-of-⟨8,9,9⟩ (rank 430)
  + three perms-of-⟨8,8,9⟩ (rank 388). Multiset:
  `{⟨9,9,9⟩·1, perm-⟨8,9,9⟩·4, perm-⟨8,8,9⟩·3}` →
  `486 + 4·430 + 3·388 = 2940`.

- **FMM's (per recipe):** two ⟨9,9,9⟩ + three perms-of-⟨8,9,9⟩ + two ⟨8,8,8⟩.
  Multiset: `{⟨9,9,9⟩·2, perm-⟨8,9,9⟩·3, ⟨8,8,8⟩·2}` →
  `2·486 + 3·430 + 2·336 = 2934`.

Net trade: FMM swaps **3·⟨8,8,9⟩-perms (3·388 = 1164)** for
**1·⟨9,9,9⟩ + 2·⟨8,8,8⟩ (486 + 2·336 = 1158)** — saving 6 mults.

## Root cause (per task taxonomy): mostly (b), partially (a)

**Primary cause — (b):** Our Strassen recombination doesn't produce the
sub-shape distribution FMM uses. We exercised every ⟨2,2,2⟩=7 scheme in
the catalog (Strassen, Winograd, fmm-lille variant) under `[9,8]³`,
`[8,9]³`, and mixed `[9,8]/[8,9]/[9,8]` allocations:

| scheme + allocation                 | shapes (sorted multiset)                            | total |
|-------------------------------------|-----------------------------------------------------|-------|
| Strassen [9,8]³                     | ⟨9,9,9⟩·1, perm-⟨8,9,9⟩·3+1=4, perm-⟨8,8,9⟩·3       | 2940  |
| Strassen [8,9]³                     | same as above (mirrored)                            | 2940  |
| Strassen mixed [9,8]/[8,9]/[9,8]    | ⟨9,9,9⟩·2, perm-⟨8,9,9⟩·3, ⟨8,8,9⟩·1, ⟨9,8,8⟩·1, ⟨8,8,8⟩·1 | 2944  |
| Winograd [9,8]³                     | ⟨9,9,9⟩·2, perm-⟨8,9,9⟩·1, perm-⟨8,8,9⟩·4           | 2954  |
| fmm-lille (Strassen) [9,8]³         | same as Strassen                                    | 2940  |
| **FMM-Lille reported recipe**       | **⟨9,9,9⟩·2, perm-⟨8,9,9⟩·3, ⟨8,8,8⟩·2**            | **2934** |

None of these vanilla rank-7 schemes / allocations produces FMM's
`{⟨9,9,9⟩·2, perm-⟨8,9,9⟩·3, ⟨8,8,8⟩·2}` distribution. The closest is
the mixed allocation (2 ⟨9,9,9⟩, 1 ⟨8,8,8⟩, total 2944) — still 10
above FMM.

The recipe is mathematically consistent (5 solos + 1 paired = 6
"operations" producing 7 sub-product instances; total mults
`5·R + pairCost(9,9,9)`), but it is **not the output of any vanilla
⟨2,2,2⟩=7 recombination on a uniform 2-axis split**. It must be
either:

1. A Pan ⟨2,2,2;7⟩-symmetric / "Lambda-algorithm" outer scheme
   (a 14-mult joint scheme for two ⟨2,2,2⟩'s that, when applied to a
   single ⟨2,2,2⟩ via aggregation/disaggregation, yields a different
   sub-product distribution than vanilla Strassen);
2. A constructive recipe with the 7 base products written down by hand,
   exploiting structure not captured by `processAdditions`'s
   max-over-nonzeros heuristic (which is intentionally conservative —
   it ignores the *algebraic* cancellation possible between sub-products
   that share rows/columns).

**Secondary cause — none in (a):** `PairedSubProducts.cyclicallyEquivalent`
*does* recognise all 9 pairs of cyclic-equivalents on our [9,8]³ sub-shape
list (3 pairs in {⟨8,9,9⟩,⟨9,9,8⟩,⟨9,8,9⟩} and 3 pairs in
{⟨8,8,9⟩,⟨9,8,8⟩,⟨8,9,8⟩}). It correctly refuses to apply Pan-pairing to
them because `pairCost > 2·R` for those small thick shapes. Predicate is
correct.

**Not (c):** `pairCost(9,9,9) = 972` is correct.

**Not (d):** No allocation is being mis-scored — every allocation
score matches what `Recombination.recombineWithAllocation` computes.

## What would close the 6-mult gap

To recover the 2934 bound from our framework, we would need one of:

1. **A specific Strassen-equivalent ⟨2,2,2⟩=7 variant** whose
   `processAdditions` on `[9,8]³` actually produces the
   `{⟨9,9,9⟩·2, perm-⟨8,9,9⟩·3, ⟨8,8,8⟩·2}` distribution. None of
   the variants currently on disk does this; the symmetry group of
   Strassen ⟨2,2,2⟩=7 modulo permutations of products is small but it's
   plausible such a variant exists. **Add a candidate variant + re-test.**

2. **A non-uniform recombination layer** that maps the 7 base products
   directly to the 7 declared FMM sub-shapes, bypassing
   `processAdditions`. This would amount to inlining the FMM recipe as
   a Java constant for ⟨17,17,17⟩ (no search), citing FMM-Lille — the
   "register the bound as cited but not yet derived" path from
   `CLAUDE.md`. **Add to `docs/cited-bounds.json` per the policy.**

3. **A Lambda-algorithm / Pan-symmetric ⟨2,2,2;7⟩ outer scheme** that
   computes two ⟨2,2,2⟩'s in 14-7-7=14 (rather than 2·7=14) base
   multiplications and re-distributes the work to fewer mults via the
   trilinear-aggregation correction. The recipe `TA(⟨9,9,9⟩,⟨9,9,9⟩)`
   strongly suggests this is the underlying construction. **Wire
   `PanPairProduct.PairScheme` into `BlockSplitSearch` as an outer base
   alongside Strassen.**

Of these (2) is the lowest-risk way to bring the catalog into line with
FMM-Lille while (3) is the right long-term fix. Both are out of scope
for the present diagnostic — the test file documents the current
behaviour for the next-step author.

## Follow-up leads (raised in chat)

Two ideas surfaced that bear directly on the gap and are worth tracking
separately because they apply far beyond the 17³ case.

### Lead 1 — Specialized "with-zeros" schemes

Some published schemes are **conditionally** rank-r when a specific
input entry (or block) is zero. Example: a `⟨3,3,3⟩` scheme that
assumes `A₁,₃ = 0` may achieve rank 21 instead of 23 over a NC ring.
When our recombination ends up with a sub-product whose input block
has known zeros (after `processAdditions`-min reduction or after a
padding/peeling decision), we currently look up a **dense** rank for
that shape — we miss the cheaper "with-zeros" specialised scheme.

Where this could bite the 17³ search: padding the outer split to
`[9,9]³` then peeling produces sub-products with one zero block-row
or block-column. A specialised rank-with-zeros variant could trim
each of those.

**Action:** sweep the catalog (FMM-Lille, Perminov ZT, Smirnov 2013)
for schemes explicitly tagged as "with zero entries" / "conditional"
/ "specialized" and import them under a new
`src/main/resources/schemes/specialized/` tree, with metadata
`"zero_blocks": [[axis, idx], …]` so the resolver can match.

### Lead 2 — Islam-MS Chapter 4 / DIS09 "S_{X,Y,Z} transform" gap

**Islam (2009 MSc) Ch. 4** describes an alternative to peel/pad that
exploits **known zeros in the OUTPUT** of a sub-product. The example
(verbatim, with Strassen padded to compute ⟨3,3,3⟩ via ⟨4,4,4⟩):

> γ5 = (Ã₂,₁ − Ã₁,₁)(B̃₁,₁ + B̃₁,₂) — no zero rows / columns in either
> factor. *However*, γ5 is used only to compute C̃₂,₂ = γ2 − γ3 + γ5
> + γ7, and we know that C̃₂,₂ has only one non-zero term (because of
> the padding), so we only need ONE entry of γ5: we reduce its cost
> from 7 to 2 multiplications.

The trick: even if `processAdditions`-on-U and -V both say the
sub-product needs full 2×2 inputs, the **W matrix tells us which
output entries are downstream-used**. When the consumer C-block is
mostly zero (due to padding), only that one nonzero C-entry needs to
be computed — so the sub-product's output dimension shrinks
correspondingly, which in turn shrinks the sub-product's mult count.

Our current `processAdditions` looks at W's nonzero pattern (returns
max alloc over nonzero positions), but it does *not* propagate
"the consumer C-block is structurally smaller than the alloc says".
For padded targets specifically (e.g. ⟨17,17,17⟩ via ⟨18,18,18⟩
Strassen-recursive with row/col 18 = padding), this leaves savings
on the table.

**DIS09 §3 (eq. 12)** introduces the same idea via S_{X,Y,Z}(U,V,W)
transforms with X, Y, Z invertible (their search uses permutation
matrices only, but the framework allows arbitrary invertibles —
Probert-Fischer-style basis-changes). Our codebase has
`SymmetryTransforms.s3Orbit` which only does the **S₃ slot-permutation**
orbit (6 transforms) — NOT the per-axis-permutation orbit that DIS09
uses. For Strassen ⟨2,2,2⟩ with a 2×2 row-permutation on the A-axis
the resulting U is row-swapped, which under `[9,8]³` allocation
genuinely shifts which sub-products land at which shape. **We do not
currently enumerate these "internal" permutation variants** in
`BlockSplitSearch` — only the S₃ slot orbit.

**Action items** (separate task — out of scope here):

a. Extend `SymmetryTransforms` with `internalPermOrbit(alg)` that
   enumerates `S_{X,Y,Z}(U,V,W)` for `(X,Y,Z) ∈ Sₐ × S_b × S_c`,
   deduplicated by signature. For Strassen this yields up to
   `2! × 2! × 2! = 8` variants; for Laderman up to `6³ = 216`.
   Feed these into `BlockSplitSearch.findBestMultiBaseSplit` as
   extra pool entries.

   **Already-verified result (see test (8) in `TestPairFusingDiagonal17`):**
   the 8 internal-permutation variants of Strassen on `[9,8]³` split
   into 2 dedup classes by multiset:

   | swap (X,Y,Z) | sub-shape multiset [⟨888⟩,⟨889⟩-c,⟨899⟩-c,⟨999⟩] | rank |
   |---|---|---|
   | (I,I,I) and (J,J,J) | `[0, 3, 3, 1]` | 2940 |
   | the other 6        | `[1, 2, 2, 2]` | 2944 |
   | FMM target          | `[2, 0, 3, 2]` | **2934** |

   **No permutation-only Strassen variant** reproduces FMM's
   distribution. So the gap is NOT closed by extending the symmetry
   orbit to S_{X,Y,Z}-perm — it requires either a *non-permutation*
   invertible basis change (Probert-Fischer 1980) or the Islam-MS
   output-reduction below.

b. Implement the **Islam-MS Ch. 4 output-driven reduction** in
   `Recombination.recombineWithAllocation`: after computing
   `(uu, vv, ww)`, additionally inspect which entries of `base.W[*,r]`
   are nonzero AND map to a non-degenerate target C-block — if the
   union of those C-block positions does NOT fill the row/col extent
   captured by `ww`, shrink `subA` / `subC` accordingly. This is a
   strictly finer bound than `processAdditions`-max. The test
   `TestPairFusingDiagonal17` will be the canonical regression
   target — under the Islam reduction the [9,8]³ Strassen
   recombination should drop from 2940 toward (but possibly not all
   the way to) 2934.

c. **Probert-Fischer 1980** (in `references/papers/`) appears to be
   the original source for the basis-change strategy DIS09 §3
   formalises. Reading their decomposition theorems for square
   matmul might reveal the *non-permutation* basis changes that
   reproduce FMM's `{⟨9,9,9⟩·2, ⟨8,8,8⟩·2, ⟨8,9,9⟩-perm·3}`
   distribution exactly, which the permutation-only orbit above
   doesn't reach.

### Why Lead 2 likely closes the 17³ gap

Quick arithmetic with the Islam-style output-reduction lens:

- FMM's claimed 2·⟨8,8,8⟩ products correspond to a Strassen sub-product
  whose 9×9 output block is structurally restricted to its 8×8 sub-block
  (because of padding to ⟨18,18,18⟩ → peel to ⟨17,17,17⟩). The
  *input* maxes might both say 9×9, but the *output* mask says only
  the top-left 8×8 entries are needed → sub-shape collapses from
  ⟨9,9,9⟩=486 to ⟨8,8,8⟩=336, saving 150 mults *per* such product.
  Two of these gives 300 mults saved.
- Meanwhile, the FMM recipe still has 2·⟨9,9,9⟩ → 972 (vs our 1·486
  + 3·388 = 1650 for the 4 thick-cyclic sub-products).
- Net: −300 (output-reduction savings) +972 (added ⟨9,9,9⟩) −486
  (removed our ⟨9,9,9⟩) −3·388 (removed our 3 ⟨8,8,9⟩-cyclic)
  +3·430 (added 3 ⟨8,9,9⟩-cyclic) = exactly the trade that lands at
  2934.

So both leads are real, both are out of scope for this diagnostic, and
both should land before the catalog claims an independent
`Q⟨17,17,17⟩ ≤ 2934`. In the meantime, register the FMM bound in
`docs/cited-bounds.json` per `CLAUDE.md` policy.

## Test artefacts

- `src/test/java/eu/solven/matmul/search/TestPairFusingDiagonal17.java`
  — prints the full trace shown above; passes (no assertion on rank;
  just structural assertions). Re-run after any change to
  `PairedSubProducts`, `Recombination.processAdditions`, or the
  ⟨2,2,2⟩=7 scheme catalog to see the impact on the 17³ gap.
