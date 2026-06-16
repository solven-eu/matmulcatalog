# Strassen vs Winograd at ⟨17,17,17⟩ — why fewer additions = fewer multiplications

A focused note on a structural observation that emerged from the
2026-06-01 cousin-hunt result: **at unbalanced cubic sizes (like n=17),
choosing the rank-7 outer scheme with fewer additions is not just an
addition-count optimisation — it's a multiplication-count optimisation**.

## The observation

At ⟨17,17,17⟩ via (9,8)³ block-split recombination:

| Scheme | Multiplications | Distribution of sub-shapes |
|---|---|---|
| Strassen 1969 (18 additions internally) | 2940 | 1·⟨9,9,9⟩ + 3·⟨8,9,9⟩-cyc + 3·⟨8,8,9⟩-cyc |
| Winograd-Strassen 1971 (15 additions, axis-flip mask=1) | **2930** | 1·⟨9,9,9⟩ + 4·⟨8,9,9⟩-cyc + 1·⟨8,8,9⟩ + 1·⟨8,8,8⟩ |

**Strassen and Winograd are GL-equivalent at ⟨2,2,2⟩=7** (de Groote 1978) — they both compute 2×2 matmul in 7 multiplications. But their **internal product structure differs**, and at unbalanced allocations like (9,8)³ this structural difference translates into a 10-multiplication gap on ⟨17,17,17⟩.

The mechanism connecting "fewer additions" to "fewer multiplications" is the topic of this note.

## Padding at unbalanced sizes

When we use Strassen at (9,8)³, block A₁₁ is 9×9 and block A₂₂ is 8×8. To form `A₁₁ + A₂₂` (Strassen's M₁ left operand), we have to either:
- **Pad A₂₂ up to 9×9 with zero rows/columns** (and similarly pad B₂₂)
- **Restrict to a common 8×8 sub-block of both** (peel A₁₁ down to its top-left 8×8)

Most recombination frameworks (ours included) implicitly pad. The padding adds zero rows/columns to the smaller operand, then performs the matmul at the larger size. The result is a 9×9 product where many of the entries reduce to zero due to the input padding — but the multiplication ITSELF was performed at the larger 9×9 shape, costing R(⟨9,9,9⟩) = 486 multiplications instead of the smaller R(⟨8,8,8⟩) = 336.

This is where **wasted padding** becomes a structural cost.

## Strassen has NO single-block products

Look at Strassen's 7 products explicitly:
```
M₁ = (A₁₁ + A₂₂) · (B₁₁ + B₂₂)
M₂ = (A₂₁ + A₂₂) · B₁₁
M₃ = A₁₁ · (B₁₂ − B₂₂)
M₄ = A₂₂ · (B₂₁ − B₁₁)
M₅ = (A₁₁ + A₁₂) · B₂₂
M₆ = (A₂₁ − A₁₁) · (B₁₁ + B₁₂)
M₇ = (A₁₂ − A₂₂) · (B₂₁ + B₂₂)
```

**Every Strassen product is at least one of "sum on A side" or "sum on B side"**. M₂ uses only B₁₁ on right but sums A₂₁+A₂₂ on left. M₃ uses only A₁₁ on left but sums B₁₂−B₂₂ on right. There is no Strassen product `A_ij · B_kl` involving single blocks only.

Consequence at unbalanced (9,8)³:
- Every Strassen product has at least one operand built from multiple blocks of different sizes
- Every Strassen product therefore involves padding on at least one side
- The padding inflates the effective sub-shape to the max of the combined block sizes
- The minimum sub-shape Strassen can produce at (9,8)³ has at least one dimension = 9 (the max of {8,9})
- **Strassen at (9,8)³ cannot produce a `⟨8,8,8⟩` sub-product** — no combination of its operations isolates the corner 8×8 blocks

## Winograd HAS single-block products

Compare to Winograd's 7 products (using Strassen-Winograd 1971's form):
```
M₁' = (A₂₁ + A₂₂)              · B₁₂
M₂' = A₂₂                       · (B₂₁ + B₂₂ − B₁₁ − B₁₂)
M₃' = (A₂₁ + A₂₂ − A₁₁ − A₁₂)   · B₁₂
M₄' = (A₁₂ − A₂₁ − A₂₂ + A₁₁)   · (B₂₁ + B₂₂ − B₁₂)
M₅' = A₁₁                       · B₁₁           ← single-block × single-block
M₆' = A₁₂                       · B₂₁           ← single-block × single-block
M₇' = (A₂₁ − A₁₁ + A₁₂ − A₂₂)   · B₂₂
```

**M₅' and M₆' are single-block × single-block products** — no sums, no padding. Whatever blocks they reference go through with their natural sizes intact.

When we apply axis-flip orbit mask=1 (swap A's row labelling, so what was "block 1" becomes "block 2" on the A-row axis), one of these single-block products now references the corner 8×8 block on both inputs. This lands a sub-product at exactly `⟨8,8,8⟩` shape → R = 336.

The result: Winograd-axflip-mask-1 has shape distribution `1·⟨9,9,9⟩ + 4·⟨8,9,9⟩-cyc + 1·⟨8,8,9⟩ + 1·⟨8,8,8⟩` = 2930.

## Why "fewer additions" → "fewer multiplications"

Winograd's structural identity: **15 additions vs Strassen's 18**. The traditional reading is "fewer additions saves a tiny linear-algebra cost". But there's a deeper reading.

Each addition operation in the 7-product identity is a "sum two A blocks" or "sum two B blocks". At unbalanced allocations, each such addition INFLATES the effective sub-shape (since the sum requires padding). Conversely, each single-block product can stay at its native size without padding.

Roughly: **a 7-product identity with fewer additions has more single-block products. Single-block products avoid padding-induced shape inflation. At unbalanced allocations like (9,8)³, this directly reduces the total multiplication count.**

The correlation isn't accidental — it's structural:
- Strassen 1969: 18 additions, 0 single-block products → all products padded → 2940 at ⟨17,17,17⟩
- Winograd 1971: 15 additions, 2 single-block products → 2 corners avoid padding → 2930 at ⟨17,17,17⟩ (after axis-flip)

## Practical recommendation

**When picking a single rank-7 outer scheme for `⟨2,2,2⟩=7` recombination at any non-power-of-2 cubic target, prefer Winograd over Strassen.** Winograd:
- Has at most 15 additions internally
- Has at least 2 single-block products
- Produces tighter shape distributions at unbalanced allocations
- Is GL-equivalent to Strassen for balanced/cubic cases (no loss when it doesn't help)

The traditional "Strassen as default" comes from balanced-power-of-2 contexts where the additions matter as constant-time overhead. **In our research-catalog context (mixed shapes, unbalanced allocations), the addition count is a structural hint** that tells us which scheme has more "structural slack" to land sub-products at smaller shapes.

## Axis-flip orbits as structural slack discoverers

The axis-flip orbit (8 mask variants per scheme) doesn't change the algebraic identity — it just relabels which block is "1" vs "2" on each axis. At balanced allocations, all 8 masks give the same shape distribution.

At unbalanced (9,8)³, the masks differ:
- For Strassen, all 8 masks land at the same 1+3+3 distribution → no help (Strassen has no single-block products to "route")
- For Winograd, 18 of the 48 cheap-orbit variants (some masks × some S₃ relabelings) land at the 1+1+4+1 distribution → axis-flip routes its single-block products toward the corner shapes

The orbit-hunt is therefore a **structural-slack discovery mechanism**: it finds which mask variant best aligns the scheme's single-block products with the target's unbalanced allocation.

## Generalisation to bigger cubic shapes

If this pattern holds:
- For ⟨n,n,n⟩ at any odd n with (⌈n/2⌉, ⌊n/2⌋) split, Winograd should beat Strassen by routing single-block products to the smaller cubic shape `⟨⌊n/2⌋, ⌊n/2⌋, ⌊n/2⌋⟩`
- The savings should scale roughly with `R(⟨⌊n/2⌋⟩) − R(⟨other shapes⟩)` — bigger relative gap → bigger Winograd advantage

This is testable: run the cousin-hunt at n ∈ {19, 21, 23, 25, 27} and see if Winograd-axflip-mask-1 consistently beats Strassen by routing single-block products to cubic ⟨⌊n/2⌋⟩-class.

## Connection to greedy ω-tracking

A greedy scheme-selection algorithm would:

1. For target ⟨n,m,p⟩ with allocation (n₀,n₁)×(m₀,m₁)×(p₀,p₁), enumerate the 8-mask × pool of rank-7 ⟨2,2,2⟩ schemes
2. For each (scheme, mask), compute the shape multiset its 7 products produce at this allocation
3. Score each multiset by `Σ R(shape)` using the current catalog ranks
4. Pick the lowest-scoring (scheme, mask) pair

The structural hint: when the pool is too large to fully enumerate, **prioritise schemes with low addition counts** — they're more likely to have single-block products that route well at unbalanced allocations.

This is the greedy-ω framing: at each decomposition level, route products toward sub-shapes with low ω_eff, prioritising the scheme that has the most "structural slack" to land its products in low-ω territory.

## TL;DR

| Question | Answer |
|---|---|
| Why does Winograd beat Strassen at ⟨17,17,17⟩? | Winograd has 2 single-block products (M₅', M₆'); Strassen has 0. At unbalanced (9,8)³, axis-flip routes Winograd's single-block products to the corner 8×8 blocks, producing a `⟨8,8,8⟩=336` sub-product (cheaper than any Strassen alternative) |
| Why "fewer additions" = "fewer multiplications"? | Fewer additions in the rank-7 identity = more single-block products = less padding-induced shape inflation at unbalanced allocations |
| Should we default to Winograd over Strassen? | Yes, when targeting non-balanced cubic sizes. Strassen offers no advantage and a known disadvantage in shape routing. |
| Can axis-flip orbits help generally? | Yes, for any scheme with single-block products — the orbit chooses which corner the single-block product lands in. Adds structural slack at unbalanced allocations. |

## References & cross-refs

- `TestStrassenVsWinogradAt17` — side-by-side test producing the comparison table
- `TestStrassenCousinHunt` — cousin-orbit enumeration that found the 1+1+4+1 distribution
- Task #105 — pre-index schemes by shape multiset at canonical allocations (the systematic version of this)
- de Groote 1978 — proves all rank-7 ⟨2,2,2⟩ algorithms are GL-equivalent. This is the abstract algebraic-equivalence; structural differences (single-block product count, addition count) emerge at unbalanced allocations.
- DIS09 §6 ⟨3,3,3⟩=25 — the γ5 zero-peel mechanism partially recovers from padding overhead. Complementary to the "avoid padding via single-block products" approach here.
