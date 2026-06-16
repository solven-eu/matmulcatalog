# Closure: global optimality & performance

Design notes (2026-06-04) for making the catalog closure (a) provably optimal
within its rule set and (b) fast, and for understanding what "global optimal"
can and cannot mean here.

## 1. The closure is a monotone DP over the shape-size DAG

For a fixed **atom set** (imported/base schemes) and a fixed **rule set**
(Kronecker product, single-axis concat on n/m/p, recombination-with-allocation),
the best achievable rank is a dynamic program:

```
r*(⟨n,m,p⟩) = min(
    atom_rank(⟨n,m,p⟩)                                 if an atom exists,
    min over Kron factorisations  r*(A)·r*(B),
    min over axis splits          r*(part1) + r*(part2),     # +p / +n / +m
    min over (base, allocation)   Σ r*(sub-block)            # recombination
)
```

Every decomposition consults **strictly smaller** sub-shapes:

- Kronecker: factors are strictly smaller in every axis.
- Concat: the split axis shrinks; the other two are unchanged → strictly
  smaller in one axis (so smaller `n·m·p`, and `max-dim` non-increasing).
- Recombination: each block is strictly smaller in every split axis.

So the dependency graph is a **DAG**, ordered by `(max-dim, then n·m·p)`. A single
bottom-up pass in that order computes the exact least fixpoint — i.e. the
**global optimum over the rule set**. No multi-round iteration is needed *if*
the order is respected and each shape's inner search (the recombination
allocation B&B) is run to optimality.

**Why the current full-range run still needs care:** it processes bottom-up but
(a) the recombination B&B is bounded (100k-node anytime, may miss the best
allocation), and (b) a few decompositions can consult a same-`max-dim` sibling
(e.g. concat on a non-max axis), so within a `max-dim` band a small local
fixpoint pass is still safer than a single sweep.

### ⇒ Process one dimension at a time, ascending

`for d = 2,3,4,…,32: close all shapes with max-dim == d` (a fresh JVM each),
reading dims `< d` from disk. This is the clean realisation of the DP order:

- **Correctness**: dims `< d` are final before `d` starts, so every `d`-shape
  sees optimal sub-ranks → no need to re-iterate lower dims.
- **Memory**: each JVM only ever builds matrices up to dim `d`; nothing
  accumulates across dimensions (the OOM in the single full run came from
  retaining built schemes across the whole range).
- **Parallelism**: within a band, all `d`-shapes are independent (they depend
  only on `< d`), so the band parallelises trivially and safely.
- **Incrementality**: when a new atom lands at dim `k`, only bands `≥ k` need
  re-closing (dirty-set propagation), not the whole catalog.

## 2. Performance levers

1. **Rank-DP before materialise.** The fixpoint is over *numbers*; it is cheap.
   Build factor matrices only for the **catalog-elected minimum** per shape, not
   for every candidate explored nor every intermediate win across rounds. (The
   first closure run's cost and OOM were Phase-2 building/retaining far more than
   the final electees.)
2. **Verify once, then trust the lineage.** Above `MATERIALISE_MAX_DIM` we store
   lineage-only stubs; rebuilding+spot-checking a ⟨30,32,32⟩ on every run is the
   dominant cost. Verify a stub when first written; on later runs, if the lineage
   and atom hashes are unchanged, skip the rebuild.
3. **Parallelise within a dimension band** (see §1) — near-linear speedup, no
   cross-shape contention because lower dims are immutable during the band.
4. **Incremental dirty-set closure** — recompute only shapes whose dependencies
   changed since the last run, keyed by (atom-set hash per sub-shape).
5. **Memoise sub-shape strategies** by `(shape, base, alloc-pattern, peel-pattern)`
   (partly done, #124/#125) so repeated recombination probes are O(1).
6. **Bound smarter, not harder.** Use the Kron/concat result as the B&B upper
   bound *before* the recombination allocation search (already wired, #145) so
   the B&B prunes aggressively and the 100k cap rarely binds.

## 3. What "global optimal" can and cannot mean

- **Optimal over the rule set + atom set** — achievable, as above. This is the
  honest target for the closure.
- **Optimal over all bilinear algorithms** — NOT achievable: that is the matmul
  tensor-rank problem (open/►hard). We can only ever be optimal relative to the
  atoms we hold and the compositions we know.

### The real gap to FMM is a *missing rule*, not a weak search

The shapes where FMM beats us are FMM applying a **composition technique our
rule set does not fully replicate**.

**What that technique is — be honest about the uncertainty.** We have NOT
inspected FMM's construction trace for the gap shapes, so we should not assert a
specific mechanism. An earlier draft of this note guessed "τ-theorem /
disjoint-sum (Schönhage)"; that was inference from the large-n recursive flavour,
not evidence. Per a conversation with **Perminov** (relayed by the user,
2026-06-04), the relevant **"serendipitous product" is a recursive *composition*
à la DIS09 (Drevet–Islam–Schost), NOT Schönhage's τ-theorem.** Treat that as the
working hypothesis; τ-theorem is at most a separate, unconfirmed possibility.
(Optimality/honesty discipline, CLAUDE.md: label hypotheses as hypotheses.)

So: the gap is a richer **serendipitous / DIS09-style composition** than our
current Kron + single-axis-concat + recombination. Adding it as a first-class
decomposition in `findBestStrategy` is what would close those gaps. Until then
our DP is globally optimal over a *strictly weaker* rule set than FMM's — which
is exactly why a bottom-up closure converges yet still trails FMM at dim ≥17
(post-closure 2026-06-04: every shape ≤16 is matched; the 2292 remaining gaps all
start at dim 17, only 6 of them cubic, ⟨17,17,17⟩ being gap 1).

**Priority ladder toward closing the gap:**
1. Per-dimension ascending closure (this note) — extracts the full optimum of the
   *current* rules, cheaply and safely.
2. Implement the **serendipitous / DIS09-style composition** rule (verify the
   exact mechanism first — read DIS09 / ask Perminov) and add it to
   `findBestStrategy` — the single highest-leverage addition for matching FMM.
   (τ-theorem #159/#160 stays a separate, lower-confidence track.)
3. Peel-position coverage beyond tail-only (#144) and multi-axis simultaneous
   splits — smaller incremental gains.

## 4. Additive complexity is a separate, parallel optimum

Rank (multiplications) and additions are independent objectives. The additive
optimum is its own NP-hard LSP-minimisation (`LinearCircuitMinimizer`, #190),
run *after* the rank-optimal scheme is chosen. Keep the two pipelines separate;
don't let addition counts influence rank election (and vice-versa) unless we
explicitly want a Pareto front.
