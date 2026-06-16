# Bud structure — capacity, value, and the flip-graph lifecycle

*(2026-06-13; companion to `SERENDIPITY_RANK_TRADEOFF.md`,
`SERENDIPITOUS_PARTIAL_PRODUCT.md`, `PROJECTION_MARGIN_TRADEOFF.md`.
Grounds the in-Java flip engine — `eu.solven.matmul.search.flip`
(`FlipScheme` / `FlipObjectives` / `FlipGraphWalk` / `MetaFlipWalk`) — and the
σ-pricing in `catalog.SerendipitousBudProduct`. Status: the capacity bound and
the lifecycle claims are theorems / direct facts about the move set; the
"where useful bases live" remarks are heuristic.)*

## What a bud is — and what it is NOT

A bilinear scheme of rank `r` for `⟨n,m,p⟩` is a set of `r` rank-one triples
`(u_l, v_l, w_l)`, with the three vectors living in the flattened sizes of the
matrices they read:

| axis | reads | dimension `d` |
| --- | --- | --- |
| U | A (n×m) | `d_U = n·m` |
| V | B (m×p) | `d_V = m·p` |
| W | C (n×p) | `d_W = n·p` |

A **bud** is a class of ≥2 products whose vectors on **one** axis are equal up
to scaling (`FlipScheme.signClasses` for the ternary/sign case;
`SerendipitousBudProduct.independentClassSizes` for the general
direction-partition). `budScore` (`FlipObjectives.budScore`) sums all class
sizes ≥2 over the three axes *independently* — one product may be counted on
all three.

The load-bearing distinction, easy to get wrong:

- A **bud** = a pair sharing **one** slot. It contributes to `budScore` and is
  the pivot a *flip* needs. **It does nothing to rank.**
- A **double-bud** = a pair sharing **two** slots (up to sign). This is the
  *only* thing `FlipScheme.mergeOnce` can collapse → **rank −1**.

Buds are cheap and abundant; double-buds are the rare coincidence rank descent
feeds on. This is why **a rank-optimal scheme can be bud-rich**: a Perminov-ZT
`⟨2,9,10⟩` seed carries `budScore = 167` while being rank-tight, because none
of those 167 single-slot shares happens to have a *second* aligned slot. There
is no contradiction between "optimal rank" and "many buds" — only between
"optimal rank" and "many *double*-buds".

## The tight axis and the capacity bound

The W-flattening of the matmul tensor has rank `n·p`, so `span{w_l} ≥ n·p`,
so the products must use at least `d_W = n·p` **distinct** W-directions
(parallel vectors add no dimension). Same per axis. The **tight axis** is the
one with the largest `d` — the axis most forced to spread products across
distinct directions, hence with the least room to pile them into shared
classes.

Counting the W-axis: with `r` products needing ≥ `d_W` distinct directions, at
most `r − d_W` directions can hold a *second* product, and each pairs up 2
products, so at most `2·(r − d_W)` products sit in classes of size ≥2:

```
budScore_axis  ≤  min( r,  2·(r − d_axis) )
```

In words, on the V-axis of `⟨2,9,10⟩` (`d_V = 90`, `r = 140`): you must use
≥90 distinct V-directions, so at most `140 − 90 = 50` of them can double up,
pairing at most `100` products → `budScore_V ≤ 100`.

Consequences:

- **As `r → d_axis`** (an optimal scheme on its tight axis) the capacity → 0:
  a rank-minimal scheme has almost no room for buds on its tight axis, and
  ample room on its slack axes.
- Capacity is a *ceiling on the count*, not a prediction. `⟨2,9,10⟩` has
  capacity ≈ U:140 / V:100 / W:140 (~380) yet HK71 realizes `budScore = 0` —
  rank-minimal Euclidean-recursive bases land at bud-sterile vertices despite
  enormous room. Capacity ≠ realizability.

## Count vs value — the same buds, two scores

`budScore` counts; **σ prices**. A size-`k` class on axis `ax`, fused against a
partner inner `I = ⟨n₂,m₂,p₂⟩`, realizes as one enlarged block instead of `k`
inner copies, saving

```
σ_ax(k) = k·R(I) − R(k-enlarged along ax)
```

(U enlarges `p₂`, V enlarges `n₂`, W enlarges `m₂`;
`SerendipitousBudProduct.costOf` / `FlipObjectives.serendipitySavingByAxis`).
The catch: `σ` is often **0**. For inner `⟨2,2,2⟩`, a size-2 bud gives
`σ(2) = 2·7 − R⟨2,2,4⟩ = 14 − 14 = 0` — a real bud that saves no rank.

So "many buds" and "useful buds" are **the same buds evaluated against a
partner**, not two kinds of bud. We collect *all* buds (`findBuds` over all 6
orderings); `σ` is the price tag read off afterward.

### σ is relational, not an intrinsic property of the scheme

This is load-bearing and easy to misread. The **bud structure** — how many
classes, what sizes, on which axes — *is* intrinsic to a scheme (it is just
counting shared directions). The **value** `σ` puts on that structure is not;
it depends on two things outside the scheme:

1. **A chosen partner inner `I = ⟨n₂,m₂,p₂⟩`.** There is no single σ — there is
   a σ *per partner you intend to Kronecker with*. `FlipObjectives.serendipitous`
   takes the inner as a parameter; `selfSerendipitous` merely fixes the
   convention `I` = the scheme's own shape (what the `self_serendipity_savings`
   JSON stamp uses) — a *choice*, not an intrinsic figure.
2. **The current catalog ranks `R(I)` and `R(k-enlarged)`.** Both are
   `FieldAwareLookup.findRank` lookups — today's SOTA for *other* shapes — so
   **σ moves when the catalog improves anywhere.** A bud worth 0 today can
   become worth >0 if a better rank lands for the enlarged shape (or worse if
   `R(I)` drops).

So a bud's worth is a *relation* (scheme × partner × current catalog), never a
label on the scheme. This is why we say "useful *against this inner*", and why
`self_serendipity_savings` must be read as "value against itself, at today's
catalog", not as a fixed property.

Canonical evidence, the two `bud-bases/section10/2x9x10-*` artifacts:

| file | seed | regime | budScore | σ-savings [U,V,W] |
| --- | --- | --- | --- | --- |
| `…-r142-…` | HK71 r140 (budScore 0) | `weighted w-rank=10 w-bud=1 split-prob=0.01 --no-full-reduce` | 0→**28** | **[0,0,2]** |
| `…-r144-…` | Perminov-ZT r144 (budScore 167) | rank-first lexicographic | 167→**170** | **[0,0,0]** |

170 buds worth **0**; 28 buds worth **2**. `budScore` is a misleading proxy for
the real currency — optimize `FlipObjectives.serendipitous(inner)` directly
when chasing a win (the `ProbeFlipBudHarvest` lesson). The pricing tradeoff —
when one rank of slack out-buys its `R(I)` cost — is formalized separately in
`SERENDIPITY_RANK_TRADEOFF.md` (`δ_max = G / (R(I) − s)`).

## The flip-graph lifecycle: where surviving buds come from

Three moves, and only these, act on a fixed-format scheme:

| move | rank | bud effect |
| --- | --- | --- |
| **flip** (`FlipScheme.flip`) | 0 | rewrites gauge on a pair that already shares one slot; can create/destroy *other* single-slot shares as a side effect |
| **merge** (`mergeOnce`, inside `reduce()`) | −1 | **consumes** a double-bud |
| **split** (`split`, plus-transition) | +1 | **creates** a double-bud |

Reading off the dynamics:

1. **A zero-class vertex is frozen.** With `budScore = 0` there is no pivot, so
   `randomFlip` returns false (the HK71 rigidity). The only exits are *split*
   (rank +1) or *extend* (`MetaFlipWalk` / `FlipScheme.extendAxis`, a new
   format whose naive slice is bud-rich by construction). So **split is the
   bootstrap for a sterile vertex only** — not the general bud source.

2. **Once any bud exists, flips churn buds for free.** A flip rewrites two
   vectors and can make a third product newly parallel (size-2 → size-3) or
   break an existing share. Most bud structure on a bud-rich seed (Perminov)
   comes with the seed and is then reshuffled by flips — *no splits involved*.
   "10 buds" does **not** mean "10 splits".

3. **Split injects +4 budScore, not +1.** Splitting product `l` on its U-vector
   makes the new product a clone of `l` on **both** V and W → a size-2 V-class
   *and* a size-2 W-class at once. But that pair shares two slots — it is a
   *double-bud*, and `reduce()` would `mergeOnce` it straight back to the
   original product (split∘merge = identity). This is why `--no-full-reduce` is
   load-bearing.

4. **The surviving bud comes from the degrade step.** The stable macro is:

   ```
   split          rank +1   create a double-bud (shares V and W)   ← mergeable, would undo
   flip one slot  rank  0   break the V-share, keep the W-share    ← now a single-bud
                                                                     no longer mergeable → SURVIVES
   ```

   You never run the merge — you *suppress* it (`fullReduce=false`) and break
   the double-bud down to a harmless single-bud before `reduce()` can cash it
   back. The +1 rank stays paid **because** you refused the merge. That is the
   whole source of a bud-base's extra rank.

**The real antagonism**, stated correctly: a *double-bud* is a fork — cash it
for rank (merge) **or** keep it as structure (degrade to a single-bud) — never
both. Ordinary single-slot buds do not threaten rank at all, which is exactly
why bud-rich optimal schemes (Perminov) exist.

## Dense → optimal, and the 8→7 lesson

From the naive scheme (`r = n·m·p`, all unit vectors, maximally bud-rich) rank
descent is the cycle *flip · flip · … until a pair shares a second slot → merge
(rank −1)*. Each rank drop is one double-bud being cashed; the descent passes
through bud-rich intermediates but strips them on the way down, so the
endpoint is bud-poor.

`⟨2,2,2⟩` 8→7 is **one merge**. The ~12k ternary flips before it are a *blind*
hunt for the single rare gauge in which a pair shares two slots, over the most
restrictive move set (sign-only, `coefCap=1`). The lesson generalizes:

```
rank-descent cost ≈ (sparsity of merge-enabling configs) × (blindness of the search)
```

Random walks pay the full product; a **directed** generator that scores each
candidate flip by how close it brings *some* pair to a two-slot share
(merge-distance gradient) cuts the second factor. The 1-step-lookahead greedy
already built does this and is enough for easy cases but stalls at sharp optima
(`⟨6,8,9⟩=296`, the V-rigid `⟨3,2,3⟩`) — there the next levers are 2-step/beam
lookahead or a math-first existence oracle (`StructuredWAls`), not a better
random move-picker.

## Directed bud generation — the principled selectors

Instead of `randomFlip`'s uniform (class, pair, variant), the theory above
gives three targeting criteria:

1. **σ-axis priority.** For most shapes only one axis pays
   (`⟨6,8,9⟩`: σ_W only; reversed `⟨3,2,3⟩⊗⟨2,4,3⟩`: σ_V only). Direct *splits*
   onto the axes that feed the paying bud-axis (split U/V to grow W-buds, never
   split W). `serendipitySavingByAxis` already says which axis pays.
2. **Merge-distance gradient** for rank descent (above).
3. **Goal-directed via ALS.** `StructuredWAls` answers "does a rank-`r` scheme
   with *this* class profile exist?" constructively — the existence oracle a
   directed walk can aim at instead of wandering.

## Factorability gates serendipity (the prime effect)

`σ` itself imposes no divisibility requirement: a size-`k` bud enlarges one axis
*by* `k`, producing `⟨n₂,m₂,k·p₂⟩` — you create a multiple, you don't need a
pre-existing divisor, and `k` is whatever the base's structure happens to be.

The factorability constraint lives one level up, in the Kronecker decomposition
serendipity refines. The serendipitous product has shape `⟨n₁n₂, m₁m₂, p₁p₂⟩`;
to *address a target* `⟨N,M,P⟩` this way each dimension must **factor** as
`N = n₁·n₂`, etc. A **prime** dimension factors only as `1·N`, so it admits no
non-trivial base⊗inner split on that axis — serendipity (and ordinary Kronecker
recursion) has nothing to work with, and the shape must fall back to a direct
construction or be padded up to a composite (wasting the slack). Smooth
composite dimensions (powers of 2/3) have many factorizations → many
recursion/serendipity routes → cheap recursive upper bounds (`⟨2,2,2⟩=7`
cascades to `⟨4,4,4⟩, ⟨8,8,8⟩, …` for free; prime cubes `⟨5,5,5⟩`, `⟨7,7,7⟩`
need bespoke schemes).

This is a property of the **construction toolbox**, not a number-theoretic law
about rank: primes don't get the recursive discount, so our *achievable upper
bounds* for prime dimensions tend to sit closer to direct/padded, while the
true lower bounds favour neither. Slogan: recursion loves composites, and
serendipity inherits the taste. A gap-map / targeted-search program should
therefore weight candidate targets by factorability, not by raw σ-gap alone.

## Empirical gap map over the Q catalog (≤32) and what it does NOT mean

`docs.explore.ScanSerendipityGapMap` prices `σ_ax(k) = k·R(I) − R(k-enlarged)`
for every known shape from `findRank` alone (no scheme generation). Over Q,
shapes ≤32, k≤4: 5,456 known shapes, 16,368 priceable σ-candidates, 32,736
dropped as unpriceable. Findings (2026-06-13):

- **Absolute σ grows with shape size, as expected.** Largest single-bud saving:
  inner ⟨8,32,32⟩ with a size-4 V-bud → ⟨32,32,32⟩, σ = `4·4608 − 15096 = 3336`.
- **Every top row is a V-bud — a catalog-boundary artifact, not a result.** A
  V-bud enlarges the *smallest* axis (`n → k·n`), keeping the enlarged shape
  in-catalog; a U-bud enlarges the *largest* axis and leaves the ≤32 catalog
  (⟨32,32,128⟩ doesn't exist) → unpriceable → dropped. The scan is structurally
  blind to large-axis buds on big shapes; this is the bulk of the 2:1
  unpriceable:priced ratio, NOT evidence that U/W buds are worthless.
- **The fractional saving `σ/(k·R(I))` caps at ~25–28% and is scale-invariant.**
  The fractional ranking is a flat wall of `⟨3,m,p⟩ → ⟨12,m,p⟩` at ~25%
  (`R⟨12,m,p⟩ ≈ 3·R⟨3,m,p⟩`: ×4 in a dimension costs only ×3 in rank). Nothing
  in the catalog beats ~28% sub-additivity on any single enlargement. So the
  *rate* of structural opportunity does not grow with scale — only the raw
  point-count does. That ceiling is a usable bound on serendipity's leverage
  (feeds `δ_max`, `SERENDIPITY_RANK_TRADEOFF.md`).

**Load-bearing caveat — σ here measures sub-additivity ALREADY BOOKED in the
catalog.** `findRank` returns the catalog's *best* rank, which for ⟨12,…⟩,
⟨32,…⟩ is itself typically a recursive/Kronecker construction. So a big σ
usually means "the catalog already exploits this fusion", and rebuilding it as
base⊗inner would merely *re-derive* the known rank, not beat it. **A high σ
flags a bud-favorable inner; it is necessary, not sufficient, for a win.** The
actual win test is downstream: does (low-rank base with the bud) ⊗ (inner)
produce a *product shape* whose rank beats the catalog's current entry — which
for an in-catalog improvement also requires the product shape to be ≤32 (the
high-σ rows have large inners → products >32 → coverage extension, not
improvement). That verified test is `docs.SerendipitousSweep` (builds +
`Verifier.isExactNonCubic` + compares to `trueSota = min(findRank, plain-Kron)`).

> The sweep's SOTA oracle was a silent casualty of the 2026-06 filename rename:
> it parsed rank from a `_m\d+` filename token that the `-r{rank}-` convention
> matched 0/14,556 files → empty bar → no win could ever fire. Fixed to read
> SOTA from `findRank` (content-driven, stub-aware); guard
> `TestSerendipitousSweep.sota_oracle_is_content_driven`.

**Empirical result of the repaired sweep** (2026-06-13, Q, baseCap=12,
secondCap=8, targetCap=16 — the range where build+`isExactNonCubic` is
tractable): **1265 bud-rich bases, 0 wins below SOTA.** The "already-banked"
caveat is confirmed in the feasible range: every serendipitous product our
existing bases can form and verify is already matched by the catalog. Note the
verify cost `(nm)(mp)(np)·r` makes the large-target regime (≤32) intractable to
*confirm* — a candidate `pred < sota` at a near-⟨…,32⟩ product hangs the sweep —
so the genuinely-unexplored estate is large products, reachable only by
prediction (the gap map) + a cheaper verification than full Brent-equation
expansion.

## Base-agnostic win-potential screen — and why serendipity self-improvement is closed (Q≤32)

`docs.explore.ScanSerendipityWinPotential` is the operational layer: for every
target `T` and factorization `T = base ⊗ inner`, it asks *what bud profile a base
of shape ⟨n₁,m₁,p₁⟩ would need to beat SOTA(T)* — pre-certifying the hit so a
directed meta-flip search has a known goal. A rank-`r₁` base with `b` size-`k`
buds gives `r_s = r₁·R(inner) − b·σ(k)`, so the win needs
`b* = ⌊(r₁·R(inner) − SOTA)/σ(k)⌋ + 1` buds, feasible iff `b*·k ≤ r₁` and within
the capacity bound.

**Naive run (2026-06-13, Q≤32):** 21,136 "goals" over 1,434 targets — and
**every one assumes a bud at the base's OPTIMAL rank.** They are dominated by
"base ⟨2,2,2⟩ r₁=7 with a size-4 V-bud", which is *impossible*: ⟨2,2,2⟩=7 is
unique up to symmetry (de Groote), bud count is a symmetry invariant, and
Strassen has zero V-buds → no rank-7 ⟨2,2,2⟩ has one. The capacity bound
(`4 ≤ 2(7−4)=6`) does NOT catch this — rigidity is stronger than capacity.

**Charging the bud's rank cost kills all of them.** At optimal rank schemes are
bud-sterile, so a bud must be *created*, costing rank. Charge even the optimistic
1 rank/bud: `r_s_real = plainKron + b·(R(inner) − σ)`. By the **28%
sub-additivity ceiling**, `σ(k) ≤ ~0.28·k·R(inner) < R(inner)` (it only reaches
`≈R(inner)` at the k=4 ceiling), so `R(inner) − σ ≥ 0` and one extra rank pushes
`r_s_real ≥ plainKron ≥ SOTA`. **Result: 0 of 21,136 goals survive.**

**Verdict (provable form of the empirical 0):** serendipity *self*-improvement
is closed over Q≤32. The σ-density `s = σ/(k−1)` never exceeds the rank
exchange-rate `R(inner)` (the ceiling guarantees it), so by the `δ_max`
analysis (`SERENDIPITY_RANK_TRADEOFF.md`) no rank relaxation ever pays — a win
requires a bud *at optimal rank*, and the only test of that (`SerendipitousSweep`
over existing bases) finds 0 (rigidity at optimum). The naive screen's hundreds
of "wins" are all the δ=0 illusion on bud-sterile-at-optimum bases. The one
residual crack: a base that is genuinely bud-rich *at* its optimal rank, beyond
`SerendipitousSweep`'s caps — a targeted higher-cap sweep would close it.

## What this means for catalog targets

A near-optimal base is the *worst* place to start a bud search: low rank ⇒
tight-axis capacity near 0 ⇒ sterile (HK71). Useful serendipity bases live at
`r = R_opt + δ`, paid for by splits and kept by refusing the merge — but only
where `σ > 0` against the intended partner, which for border-ish shapes
(`⟨2,m,p⟩`) is often almost nowhere (the `⟨2,9,10⟩` σ ≈ 0 finding). The shapes
worth the rank are those where serendipity *is* the SOTA mechanism
(`⟨6,8,9⟩` / `⟨8,9,9⟩` families). See `SerendipityCeiling` for the certified
per-axis brackets and `SERENDIPITY_RANK_TRADEOFF.md` for the relaxation window.
