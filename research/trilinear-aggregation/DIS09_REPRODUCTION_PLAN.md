# DIS09 reproduction plan

**Goal**: rebuild DIS09's per-format rank table using our modern catalog
of small-matrix algorithms (AlphaTensor, AlphaEvolve, etc.), to see how
much DIS09's results can be improved purely by swapping in better base
patterns — before any new composition techniques.

**Status: planning only — no code yet.** This doc captures what DIS09
does, what we already have, and what's missing.

---

## What DIS09 actually computed

DIS09 ([10]) tabulates the best rank they could find for each cubic
`⟨n,n,n⟩` with `n ∈ [2, 30]`, using a **fixed library of base
patterns** combined via **recursive composition**. They report both
non-commutative (Table 3) and commutative (Table 4) results.

The "Algorithm" column of each table names the WINNING base pattern
for that `n`. The patterns DIS09 uses:

### Cubic base algorithms used as recursion bases

| pattern | base format | base rank | notes |
|---|---|---|---|
| Strassen | `⟨2,2,2⟩` | 7 | Strassen 1969 |
| Winograd | `⟨2,2,2⟩` | 7 | Winograd 1971 variant (different sparseness from Strassen) |
| Winograd2 | `⟨2,2,2⟩` | 7 | DIS09's symmetry transform of Winograd |
| Laderman | `⟨3,3,3⟩` | 23 | Laderman 1976 |
| Makarov | `⟨5,5,5⟩` | 100 | Makarov 1970 (non-commutative) |
| Makarov333 | `⟨3,3,3⟩` | 22 | Makarov 1986 (commutative) |
| Waksman (cmt) | `⟨n,n,n⟩` family | `n(n²+2n-1)/2` for odd n | Waksman 1970 closed-form |

### Non-cubic "thin" base patterns

| pattern | base format | semantics |
|---|---|---|
| mul121 | `⟨1,2,1⟩` | 1×2 by 2×1 = 1×1 — actually `1·2·1 = 2` mults trivially. The "mul121" name refers to a recursion pattern where DIS09 uses this thin shape as a composition rule (not as a rank-saving base). |
| Hopcroft332 | `⟨3,3,2⟩` | Hopcroft-Kerr-style ⟨3,3,2⟩=15 used as a base shape for compositions |

### TA (trilinear aggregation)

Pan's trilinear aggregation framework, instantiated by DIS09 with
their own parameter choices. Approximates `(n³ + lower-order terms)`
for various n. **This is the family we already mostly cover** via
our `BlockSplitSearch` + `Compose.blockSplitCubic`.

---

## What we have implemented

✅ **Strassen base recursion** — `Compose.kroneckerGeneral` +
`Compositions.strassen3()` etc. Verified up to `⟨32,32,32⟩`.

✅ **Block-split via Strassen outer + mixed-shape inner** —
`Compose.blockSplitCubic` + `Recombination.constructWithAllocation`.
This generalizes Sedoglavic's `⟨7,7,7⟩=250` recipe and recovers it
exactly. Materialized 23 verified schemes for non-cubic gaps.

✅ **Per-field rank lookup + recursive propagation** —
`BlockSplitSearch` walks the catalog, applies Sedoglavic's closed-form
identity, propagates derived bounds back into the lookup so larger n
can benefit from smaller-n improvements.

✅ **Rosowski commutative formulas** — `RosowskiBound` covers Thm 2/3
(bilinear commutative) and Thm 4/5 (non-bilinear, non-commutative
via transpose involution). Used to populate derived bounds.

## What we don't have

❌ **Laderman base recursion** (`Recombination.constructWithAllocation`
with Laderman ⟨3,3,3⟩=23 as outer base instead of Strassen ⟨2,2,2⟩=7).
Should be a few lines — just pass `Laderman23.get()` instead of the
Strassen scheme. But the result quality depends on Laderman's specific
sparsity vs Strassen's, and on the per-mult sub-shape distribution
under min-reduction.

❌ **Winograd / Winograd2 base recursion**. Same as Laderman — we have
no `Winograd7.java`; we'd need to encode it. Different addition
profile from Strassen but same rank.

❌ **Makarov base recursion** for `⟨5n, 5n, 5n⟩`. Requires Makarov's
`⟨5,5,5⟩=100` algorithm explicitly encoded.

❌ **mul121 pattern as a composition rule**. DIS09 uses this for
n ∈ {9, 12, 15, 19, 22, 25, 27, 29} (commutative) and {9, 12, 27, 28}
(non-commutative). The pattern reduces some sub-products via thin
matrix-vector intermediate steps. Needs study of DIS09 Section 3.

❌ **Hopcroft332 pattern** — DIS09 uses for `n ∈ {17, 18, 23, 24}`
(commutative). Recursion via a ⟨3,3,2⟩-shaped sub-problem.

❌ **Multi-base search** — for each target n, try ALL bases and pick
the min. Our current `BlockSplitSearch` only uses Strassen as outer.
The optimal base for medium n is often Laderman or Winograd.

---

## Reproduction strategy

### Phase A: encode missing base algorithms

1. Hardcode `Winograd7.java` (variant 1) and `Winograd7Sym.java`
   (DIS09's Winograd2 = symmetric transform). Export as scheme JSON.
2. Hardcode `Makarov5.java` for `⟨5,5,5⟩=100`. We don't have the
   explicit factor matrices on disk; would need to transcribe from
   Makarov 1970 (Russian) or from an English source.
3. Verify all three with `Verifier.isExactNonCubic`.

### Phase B: multi-base block-split search

Extend `BlockSplitSearch.findBestSplit` to accept a list of outer
bases (currently fixed to Strassen). For each cubic target:
- Try each base
- For each base, try each allocation
- Pick the min rank
- Record which (base, allocation) won

This reproduces DIS09's logic for the cubic-base patterns
(Strassen, Winograd, Laderman, Makarov).

### Phase C: mul121 / Hopcroft332 patterns

Read DIS09 Section 3 carefully to understand:
- How mul121 reduces sub-product cost
- How Hopcroft332 enables thin-shape composition

Both are non-cubic outer bases (⟨1,2,1⟩, ⟨3,3,2⟩) applied with
specific allocation patterns. May require extending
`recombineWithAllocation` to handle non-cubic bases (currently
Strassen is the only base tested in practice).

### Phase D: post-DIS09 base swaps

Replace DIS09's bases with our modern catalog entries where they're
better:
- ⟨4,4,4⟩=47 (AlphaTensor F₂) → use as base for `⟨4n,4n,4n⟩` over F₂
- ⟨4,4,4⟩=48 (AlphaEvolve C) → use as base for `⟨4n,4n,4n⟩` over C
- ⟨5,5,5⟩=93 (AlphaEvolve R) → use as base for `⟨5n,5n,5n⟩` over R
- ⟨3,3,3⟩=23 (Laderman) — already in our catalog

For each cubic target n, the best base + allocation likely changes
when the catalog updates.

### Phase E: compare to DIS09 table

Emit a side-by-side comparison: DIS09 vs us, per (format, field).
Highlight where modern catalog improves DIS09's results. This becomes
a new "DIS09 vs SOTA" Markdown table or Pages section.

---

## Open questions — RESOLVED (2026-05-28)

1. **Makarov 1970's `⟨5,5,5⟩=100`** — **SKIPPED**. The catalog already
   has much better `⟨5,5,5⟩` results (Perminov / AlphaEvolve /
   Moosbauer at r=93, AT-Z at r=98). For modern DIS09 reproduction
   we'd use r=93 as the `⟨5,5,5⟩` base, making the historical
   Makarov 100 a non-target. Russian original PDF is still
   inaccessible but no longer load-bearing.

2. **What is `mul121`?** — It is literally `⟨1,2,1⟩=2` (the trivial
   1×2 × 2×1 dot product) used as a Kronecker base for
   axis-splitting compositions. Per DIS09 §3 the algorithm is:
   ```
   U = [[1,0],[0,1]]    V = [[1,0],[0,1]]    W = [[1,1]]
   ```
   Used to express `⟨n, m, p⟩ = ⟨n, m₁, p⟩ + ⟨n, m₂, p⟩` (middle
   axis split). Same idea for `mul211 = ⟨2,1,1⟩=2` (first-axis
   split) and `mul112 = ⟨1,1,2⟩=2` (third-axis split).

3. **What is `Hopcroft332`?** — It is just `⟨3,3,2⟩=15`
   (Hopcroft-Kerr 1971), used as a Kronecker base. Same for
   `Hopcroft233 = ⟨2,3,3⟩=15` and `Hopcroft323 = ⟨3,2,3⟩=15` (all
   the same algorithm under axis permutation; the three names refer
   to which axis carries the "thin" 2). **We already have this in
   the catalog.**

4. **Non-cubic outer bases in our framework** — Should work:
   `Recombination.recombineWithAllocation` and the helper
   `constructFromResult` take any `NonCubicBilinearAlgorithm` as base.
   The per-axis math (`processAdditions`, `embedFactor`, cumulative
   block offsets) is shape-generic. **Untested for non-cubic bases**
   but no design barrier identified.

## Revised implementation strategy

Most of DIS09's "missing pieces" are already in place. Concrete work:

### Phase A (small): trivial axis-split bases

Encode `Mul121.java`, `Mul211.java`, `Mul112.java` (10 lines each) +
the corresponding scheme JSON files at
`section2/mul121_1x2x1_r2.json` etc. These join the catalog as
candidate outer bases for the multi-base search.

### Phase B: multi-base search

Generalise `BlockSplitSearch.findBestSplit` to iterate over a
**configurable pool of outer bases** (today: Strassen only). Default
pool: { Strassen `⟨2,2,2⟩=7`, Laderman `⟨3,3,3⟩=23`, Hopcroft-Kerr
`⟨2,3,3⟩=15` and axis-perms, mul121/211/112, plus modern AlphaTensor /
AlphaEvolve cubic schemes }. Per target, pick the (base, allocation)
giving min total rank.

### Phase C: symmetry transforms (Winograd2-style)

DIS09 §3 generates per-base "symmetry transforms" via
`S_{X,Y,Z}(U,V,W)` where X, Y, Z are axis-permutation matrices. This
produces equivalent algorithms with different support patterns (and
hence different sub-shape distributions under min-reduction). Cheap
to add — multiply factor matrices by permutation matrices.

### Phase D: bake into derived-bounds

After A–C, regenerate `docs/derived-bounds.json` with the enriched
search. Many cubic and non-cubic rows should improve. Expected: lots
of `DIS09 (Strassen base)` rows in cited-bounds.json get superseded by
`derived: BlockSplitSearch (best base = AlphaEvolve, allocation=...)`
rows.

### Phase E: explicit DIS09 vs SOTA table

Generate a side-by-side comparison page in Pages (or as a section
header in `docs/index.html`): for each cubic `n ∈ [2, 30]`, show
DIS09's value (from cited-bounds.json), our best derived value, and
the improvement.

Phase A unlocks B; B alone delivers most of the value. C is polish.
D + E are presentation.
