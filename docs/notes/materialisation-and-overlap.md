# Materialisation vs evaluation — and the non-overlap property

A methodology note. Two intertwined points, both load-bearing for
reading our sweep outputs honestly:

1. **What materialisation actually adds, beyond writing JSON.**
2. **Why our generic recombination cannot produce overlapping
   multiplications** — and why that matters for how we should be
   compared to hand-crafted schemes.

## Materialisation is *almost* purely a serialisation step

Our search produces two outputs per shape:

- **Evaluate mode** — a *predicted* rank `r̂ = Σ_k sota.getRank(sub-shape_k)`,
  with a strategy descriptor (outer base, allocation, peel) and the
  list of sub-shapes pulled from the catalog. No (U, V, W).
- **Materialise mode** — concrete factor matrices instantiating the
  strategy via `Recombination.constructWithAllocation`, written to disk
  as a `…_r{rank}_a{additions}.json`.

What does materialise mode *uniquely* contribute? Three things, none
of them are structural to the closure loop:

| Need | Does materialise solve it? | Could it be done another way? |
| --- | --- | --- |
| Random-matmul spot-check (`Verifier.passesRandomMatmulSpotCheck`) | yes — needs concrete (U, V, W) | yes — verify on a later pass over new schemes only, or full-catalog sanity sweeps. Not needed at write-time. |
| Lineage record | yes — stamped into the JSON | yes — lineage is recursive over sub-lineages; can be carried as a tree of `(strategy, alloc, sub-lineages)` without materialised matrices. |
| Addition count for the filename | yes — counted on the materialised (U, V, W) | yes (in principle) — additions are a recursive formula over sub-additions plus the strategy's outer additions (which are computable from the strategy template alone). |
| Reusability as a sub-product in future closure rounds | yes — `lookup.find` returns `Optional<NonCubicBilinearAlgorithm>` and `constructWithAllocation` needs the inner (U, V, W) to embed | **only if** the outer round itself materialises. If both rounds stay in evaluate mode and propagate `(rank, additions, lineage)` tuples, no materialisation is needed at the boundary. |

So the closure loop as a *pure search* — find the best rank chain for
each shape — does not need materialisation. The current code couples
them because we *also* want a registered catalog scheme at the end.
That coupling is a deliberate UX choice, not an algorithmic
requirement.

### Practical consequence

A future `--mode=closure-evaluate-only` is feasible and would let us
iterate much faster: round N propagates `(rank, additions, lineage)`
tuples; round N+1 reads them; only the final accepted frontier gets
materialised + verified + written, in one batch pass at the end. We
don't have this mode today — task to file if we ever care about
sweep throughput more than per-round disk artefacts.

## Multiplications never overlap in generic recombination

`Recombination.constructWithAllocation(outer, lookup, sota, allocA,
allocB, allocC)` does **block-substitution**. For each of the
`r_outer` products of the outer scheme, it instantiates the sub-scheme
matching that product's sub-shape and wires its `r_sub_k` columns into
distinct block-positions of (A, B, C). The resulting (U, V, W) has
exactly `Σ_k r_sub_k` columns. **No two columns compute the same
scalar bilinear form.**

This is structural, not coincidental:

- Outer products are indexed by `k ∈ [0, r_outer)`, and product `k1`'s
  columns are wired to a different `(A-block, B-block, C-block)`
  signature than product `k2`'s — by the bilinear form of the outer
  scheme.
- Even when two outer products use the same sub-shape **and** the same
  sub-scheme (common with symmetric allocations), they remain
  algebraically distinct because they read different input blocks and
  write to different output blocks.

So:

> **The predicted rank `r̂ = Σ_k r_sub_k` is the exact rank of the
> materialised scheme.** No hidden discovery happens during
> materialise. If `evaluate` says 2930 and `materialise` says 2930,
> that's arithmetic, not luck.

(Materialise *can* report a different rank only if there's a bug —
e.g. an overzealous peel that violated the recombination identity, or
a stale catalog lookup. SchemeSweep guards this with a
`predicted != actual` check at the call site.)

### Sharing exists, but in *separate constructors*

Several catalog-aware constructors deliberately produce schemes with
fewer than `Σ_k r_sub_k` columns by fusing across outer products:

- `PairFusedRecombination` (cyclic pair-fusion, task #99) — re-uses a
  sub-product across two outer products with matching structural
  positions.
- Pan TA (`PanTrilinearAggregation`) — aggregation identities sum
  several outer products through a shared bilinear core.
- Kin-row reductions (FMM-Lille ⟨17,17,17⟩=2934 → 2931) — collapse
  multiple disjoint-sum components into a shared bilinear form.
- SZ 2025 disjoint-sum recipes — structural pre-conditions that
  enable cross-component sharing.

These are *separate code paths* invoked when the strategy is tagged
for them. Generic recombination never falls into them silently.

### Warning to future-us: this property is current-state, not eternal

The non-overlap claim depends on `constructWithAllocation` doing
*pure* block-substitution. If we ever introduce something that
"salts" the materialisation — e.g. random sign-changes per
sub-block to enrich the symmetry orbit, a post-pass that fuses
duplicate columns after construction, or any coefficient-mixing
between outer products that go through the same sub-shape — the
"predicted == actual" invariant no longer holds and overlapping
multiplications can appear.

If such salting is ever added, two things change at once:

1. The `predicted != actual` guard at the SchemeSweep call site stops
   being a bug-catcher and starts firing legitimately (actual could
   be *lower* than predicted, which is a feature).
2. The "evaluate is exact" claim downstream of search-only modes
   becomes a lower bound, not a tight rank.

Search-only pipelines must then either (a) keep using the unsalted
`evaluate` path for rank propagation and apply salt only at the
final materialise stage, or (b) re-derive predictions to account for
the salting. Pick one and document it; don't let the search loop
believe its own salted prediction.

### Why this matters methodologically

Most hand-crafted small-matmul schemes in the literature are *built*
around discovering shared multiplications. Strassen's `⟨2,2,2⟩=7`
beats the naïve 8 because two of his 7 products cleverly do double
duty across multiple output cells. Laderman's `⟨3,3,3⟩=23` is the same
story at larger scale. Smirnov, Pan, AlphaTensor — all of them are, in
one way or another, in the business of finding non-obvious shared
bilinear cores.

**Our generic search is *not* in that business.** It searches over:

- which **outer base** to use (Strassen, Winograd, AlphaTensor 2×2×2,
  Smirnov 3×3, …),
- which **allocation** (`a + a' = n`, etc.) to split along,
- which **peel** to apply (output-side zero suppression),
- and **axis-flip orbits** of each base.

…and assembles the result by *embedding sub-products from the
catalog*. The cleverness comes from picking the right composition,
not from inventing a new bilinear identity. The non-overlap property
makes this an honest description of what the algorithm does.

This has two consequences worth keeping in mind:

1. **Our search can never beat a hand-crafted scheme by re-discovering
   its sharing pattern accidentally** — there's no mechanism for it.
   The wins come from composing better cataloged building blocks.
2. **When we close a gap that a hand-crafted recipe also closes** —
   e.g. ⟨17,17,17⟩=2930 via Strassen-recombination at (9,8)³ vs
   FMM-Lille's disjoint-sum recipe — those are *different* schemes
   landing at the same rank, not the same scheme found by a different
   route. The catalog should keep both, with distinct lineages.

The closure search and the hand-crafted-recipe space are two
complementary directions on the same SOTA frontier — neither subsumes
the other.

## Cross-refs

- `Recombination.constructWithAllocation` — the block-substitution
  step; non-overlap is structural to it.
- `PairFusedRecombination`, `PanTrilinearAggregation`,
  `KinRowReduction` — the sharing-aware constructors (separate paths).
- `SchemeSweep` — couples evaluate + materialise per round; this note
  argues the coupling is a UX choice, not an algorithmic requirement.
- `docs/notes/why-axis-flip-helps.md` — companion note on why
  axis-flip orbits matter despite the asymmetry from mixed-size sums.
