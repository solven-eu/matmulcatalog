# Serendipity vs rank — the relaxation tradeoff, formalized

*(2026-06-12; companion to `SERENDIPITOUS_PARTIAL_PRODUCT.md`,
`BUD_STRUCTURE_THEORY.md` (bud capacity / count-vs-σ / flip lifecycle), and
`SerendipityCeiling`. Status: the inequalities below are theorems within
the stated scope; the leverage hypothesis at the end is a conjecture.)*

## Setup

Fix an outer base shape `⟨n,m,p⟩`, an inner shape `I = ⟨n₂,m₂,p₂⟩`, and a
field (everything here is per-field; catalog prices are the Q-chain unless
stated). The serendipitous product (Smith 2002 eq. 69) prices a base scheme
`B` of rank `r` with bud structure as

```
cost(B) = r·R(I) − Σ_axes Σ_classes σ_ax(k_c)
```

where `σ_ax(k) = k·R(I) − R(k-enlarged along ax)` is the **fusion saving** of
a size-`k` bud class on axis `ax` (U-buds enlarge `p₂`, V-buds `n₂`, W-buds
`m₂`). All σ are catalog-priced, so every quantity below moves when the
catalog improves.

Define the **structure value function**

```
S*(r) = max { Σ σ(k_c) : a rank-r base scheme realizes the class profile }
```

— the best priced bud structure *realizable* at rank `r` (per axis or summed;
the per-axis version is what `SerendipityCeiling` brackets). `S*` is unknown
in general; we hold two handles on it:

- **Lower handle**: any concrete scheme (catalog: `S*(20) ≥ 4` for `⟨2,4,3⟩`
  against inner `⟨3,2,3⟩` — the four W-pairs of the record base).
- **Upper handle**: the certified ceiling `S_ceil(r)` (knapsack over the
  σ-table under the spanning bound `#classes ≥ flat` and the contracted-axis
  divisibility bump — `SerendipityCeiling.maxSavings`).

The product cost at relaxed rank `r₀+δ` is `C(δ) = (r₀+δ)·R(I) − S*(r₀+δ)`,
so **relaxation by δ pays iff**

```
S*(r₀+δ) − S*(r₀) > δ·R(I)            (RELAX)
```

— the marginal structure must out-buy rank at the exchange rate `R(I)` per
unit. This is the user's intuition made exact: a very tight-rank base may
also be structure-tight (small `S*(r₀)`), and one rank of slack may unlock a
disproportionate structure jump.

## The ceiling slope bounds the window

`S_ceil` grows at most linearly: converting singletons into fused classes,
each extra product spent on a size-`k` class yields at most

```
s := max_k σ(k)/(k−1)        (the σ-DENSITY of the inner)
```

savings per product. Hence `S_ceil(r₀+δ) ≤ S_ceil(r₀) + δ·s`, and combining
with (RELAX) and `S*(r₀+δ) ≤ S_ceil(r₀+δ)`:

```
δ · (R(I) − s)  <  S_ceil(r₀) − S*(r₀)  =: G(r₀)   (the FRUSTRATION gap)
```

**Theorem (relaxation window).** If `s < R(I)` (always true when fusion
saves less than a full inner block per extra product — every case priced so
far), rank relaxation can only pay within

```
δ < δ_max = G(r₀) / (R(I) − s).
```

Beyond `δ_max`, even a base that *hits the ceiling* loses to the incumbent.
The window is wide exactly when (a) the incumbent is far below its own
ceiling (large frustration `G`) and (b) the inner has high σ-density `s`.

**Worked instance — `⟨2,4,3⟩` ⊗ `⟨3,2,3⟩` (W axis, the only paying axis):**
σ_W = {2:1, 3:5, 4:5} → `s = 5/2`; `R(I) = 15`; `S_ceil(20) = 31`,
`S*(20) ≥ 4` → `G ≤ 27` → `δ_max < 27/12.5 ≈ 2.16`. So **only rank 21–22
could ever pay**, no matter how clever the construction. Rank 21 needs
`σ ≥ 20` (≈ 4 W-triples); the tied-W ALS found none (weak evidence — see
power caveat below) and rank 22 would need `σ ≥ 35` = the exact ceiling
(7 triples at the minimum class count) — extremal, and triples already look
obstructed at 21. The honest status: rank 20 richer-than-4-pairs is
*undecided-leaning-no*, the relaxation window is *provably* ≤ 2.

## Diminishing returns and the maximum effect

The user's "bud-structure has a maximum effect" is the statement
`S*(r) ≤ S_ceil(r)` plus the knapsack saturating: once every class is at the
best density size `k*`, extra rank adds at most `s` per unit while costing
`R(I)` — so `C(δ)` is eventually strictly increasing in δ with slope
`≥ R(I) − s > 0`. The cost curve is (weakly) U-shaped in δ; the minimum sits
within the window above.

## Leverage hypothesis (conjecture, testable)

Define the **serendipity leverage** of an inner shape as `λ(I) = s/R(I)`
(dimensionless, in `[0,1)`). The window `δ_max = G/(R(I)(1−λ))` widens as
`λ → 1`. Empirically λ grows with inner size — fusion digs deeper into
sub-additivity as shapes leave the `⟨2,m,p⟩` border (e.g. λ = 0 for U/V axes
of `⟨3,2,3⟩` since R⟨3,2,3k⟩ = 15k is exactly additive, vs λ = 1/6 on its W
axis; bigger inners should do better — *this is the user's "more and more
true when shape size increases"*). The conjecture: for inner families of
growing size, `λ(I)` is non-decreasing and bounded away from 0 on the
contracted axis, so relaxation windows widen with scale and rank-tight bases
are increasingly NOT serendipity-optimal. Test: tabulate λ per axis for all
catalog inners (cheap — pure `findRank` arithmetic) and regress against
inner volume.

## Decision power caveat (load-bearing)

Cold-start tied-slot ALS has **no power at tight rank** (the known-feasible
4-pair `⟨2,4,3⟩` control solves 0/300 cold; unconstrained rank-20 solves
0/100; over-rank 21/23 solve only ~2/100). Therefore: a cold ALS "no" at
tight rank is UNINFORMATIVE, not evidence. `ScanSerendipityFrontier` runs an
unconstrained control per (shape, rank) and labels verdicts accordingly
(`EVIDENCE` vs `UNINFORMATIVE`). Positive results (`SOLVED`) are always
valid (constructive, spot-checked; rationalize before catalog entry).
Certified NOs need fixed-structure SAT/Gröbner. Warm local repairs (keep
U/V + known classes, tie one more pair) retain power near the incumbent —
the 66-pairing sweep around the `⟨2,4,3⟩` record is genuine local evidence.
