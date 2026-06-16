# Multiset frontier: a multiset-first reframing of recombination search

Design note (2026-06-04). Consolidates the recombination/allocation/orbit/TA
algorithms under one objective, after a discussion clarifying what is actually
being minimised. Companion to `CLOSURE_OPTIMALITY_AND_PERF.md`.

## The core object is the *multiset*, not the allocation

A plain recombination of an outer base `⟨bn,bm,bp⟩` at a target `⟨n,m,p⟩`
produces a **multiset** of sub-product shapes, and its rank is exactly
```
cost(M) = Σ_{shape ∈ M} R*(shape)        (R* = best/closure rank of the sub-shape)
```
Nothing else about the base survives into the rank. So **the multiset is a
complete rank invariant**: any two routes (different base, ordering, allocation)
that land on the same multiset have identical rank. The multiset is the decision
variable; everything else is a *generator* of multisets.

Caveats (state them honestly — CLAUDE.md optimality discipline):
- The "complete invariant" claim holds for **plain additive recombination
  only**. With cross-product sharing (TA/disjoint-sum/cancellation) the
  *structure* matters, not just the bag — see the TA-aware cost below, which is
  the controlled exception.
- Multiset equivalence is for **rank**, not **additions** (different bases reach
  the same multiset with different addition structure). ⇒ two-pass (below).

## Three levels (this removes the "vaguely similar algorithms" confusion)

| level | object | operation |
| --- | --- | --- |
| **generators** | (base support pattern) × (ordered allocation) | *produce* a multiset |
| **cost** | a multiset `M` | `cost(M)` = TA-aware sum (below) — a closure lookup, not a min |
| **objective** | the target | `min over { distinct reachable M }` |

The allocation/base/ordering are **not** an inner min of a multiset's cost —
they enumerate *which multisets are reachable*; the cost of each is a (TA-aware)
sum; the answer is the min over the reachable set.

### Axis-flip ≡ allocation ordering (do NOT double-count)

For a 2-part base (`⟨2,2,2⟩`), **axis-flip on an axis = swapping the 2-part
allocation order on that axis** — provably the same multiset (both-block
products are flip-invariant ⟨max,max,max⟩; single-block products have their size
`a↔b` swapped identically by either operation). So:
```
multiset = f( support pattern,  per-axis part SIZES,  per-axis part ORDERING )
                                 └─ "allocation" ─┘    └─ "axis-flip" ─┘
```
Sizes × ordering are orthogonal coordinates, but **ordering *is* the flip** — a
fully *ordered* allocation already contains it. Enumerate **ordered allocations**
and flips come for free; or **sorted allocations × flips** — equivalent. Doing
both is double-counting. (The current closure does sorted-alloc + canonical base
= neither ordering ⇒ it gets ⟨17³⟩=2940; the cousin hunt adds flips ⇒ 2930.)
For bases with ≥3 parts/axis, "flip" (Z₂ reverse) is only a subgroup of the full
part-permutation `S_k` — use **AxisPermute**, not just AxisFlip.

### Cross-base dedup (per-target, not per-base)

Today each base runs its own allocation search, blind to the others. The
frontier is **per-target**: take the **union** of reachable multisets across all
bases (each contributing its ordered-allocation multisets), **dedup by canonical
multiset** (a bag reachable from two bases is costed once), then min. The base
stops being the search axis and becomes one generator among many.

## TA-aware cost (makes "seek TA-friendly multisets" automatic)

Trilinear aggregation (Pan) fuses a shape with a cyclic **rotation** of itself:
one `⟨a,b,c⟩` **+** one `⟨b,c,a⟩` computed jointly at
```
pairCost(a,b,c) = abc + ab + bc + ca        (rectangular; cubic = n³+3n²)
```
beating `R*(⟨a,b,c⟩) + R*(⟨b,c,a⟩)` **only when** `R* > pairCost/2`, i.e. when
the single-shape scheme is weak (near-naive). So the cost is:
```
cost(M) = Σ over shape-groups of min( plain sum,  best cyclic-pair fusion )
```
implemented by `PairedSubProducts.applyPairing(multiset, R*)` (greedy
max-savings matching). With this cost, the frontier **automatically**:
- discovers FMM-style "concentrate repeats on a weak cubic + TA" (e.g. FMM's
  `⟨17³⟩=2934` uses `TA(⟨9,9,9⟩,⟨9,9,9⟩)=972` because its ⟨9,9,9⟩ was weak), and
- prefers the cousin "spread onto strong shapes, no TA" when cheaper (our
  `⟨17³⟩=2930`, single ⟨9,9,9⟩=486 — TA there only ties),

and picks the min — no manual hunt.

**Pairing correctness (fixed 2026-06-04, `panPairable`).** A Pan pair computes
one `⟨a,b,c⟩` **+** one `⟨b,c,a⟩` (a shape + a non-trivial cyclic rotation). The
key subtlety: **transpose is a free isotropy** (`⟨a,b,c⟩` ↦ `⟨c,b,a⟩` via
`Cᵀ=BᵀAᵀ`, zero arithmetic), so each product may be transposed before matching.
Consequence for two *identical* copies:
- cubic `⟨n,n,n⟩` → pairable;
- exactly two equal dims `⟨8,9,9⟩` → **pairable**, because `⟨8,9,9⟩ᵀ = ⟨9,9,8⟩`
  *is* the rotation;
- all distinct `⟨8,9,10⟩` → **not** pairable (`⟨10,9,8⟩` is still not a rotation;
  and cyclic rotation is *not* a free isotropy of a bilinear product — it would
  move the output into an input slot).

The earlier `cyclicallyEquivalent` (pure rotation, includes identity) over-paired
*all* identical shapes — wrong for the all-distinct case. `applyPairing` now uses
`panPairable` (tries both transposes × the two non-trivial rotations);
`cyclicallyEquivalent` is retained only as the pure-rotation utility.
**Still missing:** the 3-way cyclic *triple* aggregation
(`⟨a,b,c⟩+⟨b,c,a⟩+⟨c,a,b⟩` jointly) — `applyPairing` does 2-way only; adding it
needs the verified rectangular 3-way rank formula (Pan 1982 / DIS09 /
`PanTrilinearAggregation`), tracked as a follow-up rather than guessed.

## Two-pass: rank then additions

- **Phase 1** — minimise **rank** over distinct multisets (above). Tier:
  *optimal-within-scope* (rule set + atom set + 2-way TA), and only a *bound*
  where the allocation/orbit enumeration is itself bounded.
- **Phase 2** — among rank-optimal **routes** (a min-rank multiset may be reached
  by several bases/orderings, each with different addition structure), minimise
  **additions** via the SLP minimiser (`SchemeAdditiveComplexity`, a *bound*).
  Equal multiset ⇒ equal rank but not equal additions, so additions can only be
  decided after fixing the route.

## What exists vs what's new

- `AllocationOptimizer` / `AssignmentOptimizer` — "min over allocations for ONE
  base" (allocation-centric; inner cost is the multiset sum). Reusable as a
  generator.
- Analytical mask search (#106) — multiset per axis-flip orbit member of a 2×2×2
  base. The orbit-coverage generator.
- Multiset signature cache (#124/#125), distribution-indexed search (#105) — the
  dedup-by-multiset machinery.
- `PairedSubProducts.applyPairing` — the TA-aware cost (2-way pairs).
- **New:** a `MultisetFrontier(target)` that unifies these — generators (base
  support × ordered allocation, across all bases) → canonical-multiset dedup →
  TA-aware cost → min; records the cheapest-additions route for materialisation.

## The ⟨2,2,2⟩ base/support enumeration (the open input)

Step-1 generators need the set of distinct rank-7 `⟨2,2,2⟩` **support patterns**
(only the support — not the coefficients — affects the multiset). By de Groote
1978 all rank-7 decompositions form one continuous `GL₂³` orbit; the relevant
*discrete* set is the distinct sparse normal forms and their `Z₂³×S₃` orbits.
Grounded in **Burichenko 2014** (arXiv:1408.6273; isotropy of Strassen ⟨2,2,2⟩=7
has order 36 — tells us how many distinct orbit members exist and why some
collapse to the same multiset) and **Heun 1994** (rank-7 normal forms). This is
also task #168.
