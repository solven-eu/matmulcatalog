# τ-theorem disjoint-sum identity search — scoping

**Status:** scoping (2026-06-06). Motivated by the FMM-vs-us comparison: of the
70 maxDim≤17 formats where FMM beats us, the bulk are `⟨a,b,17⟩` built by a
**τ-theorem disjoint-sum (N-term) identity** we cannot currently synthesise.

## Canonical worked example — ⟨6,17,17⟩

FMM (https://fmm.univ-lille.fr/6x17x17.html) computes **⟨6,17,17⟩ = 1106** as a
single trilinear-trace identity over **15** ⟨3,·,·⟩ blocks:

```
8·⟨3,6,6⟩=80  +  3·⟨3,6,5⟩=68  +  3·⟨3,5,6⟩=68  +  1·⟨3,5,5⟩=58  =  1106
```

- Axis partition: `n=6 → 2·3`, `m=17 → 6+6+5`, `p=17 → 6+6+5`.
- The **full block grid** is `2·3·3 = 18` blocks (sum = 1300). The identity uses
  **15** — the trilinear-trace ("disjoint sum") shares computation so 3 blocks
  drop. That 18→15 saving is the entire win.

**We already hold every block at FMM's exact rank** (⟨3,5,5⟩=58, ⟨3,5,6⟩=68,
⟨3,6,6⟩=80). So this is **not a block-coverage gap** — it is purely a
**combination-identity** gap. Our search reaches only **1154** (+48) via an
additive stitch `perminov⟨5,6,8⟩ +ₘ ⟨6,8,12⟩ +ₚ naive⟨6,1,9⟩ +ₘ ⟨6,9,16⟩`.

**Honesty caveat** (per `feedback-dont-overclaim-cousin-wins`): we have *not*
proven that no Strassen-cousin / recombination reaches 1106. The claim is only
that FMM's *published* recipe is a τ-identity our engine cannot generate.
⟨6,17,17⟩ is the **cheapest acceptance test** for any τ-identity search (all
blocks ⟨3,5..6,5..6⟩, target maxDim 17 → fully materialisable + verifiable).

## What exists today

- **`KnownTauIdentities`** — hand-extracted published identities (Schönhage τ,
  Pan borrow-and-correct), wrapped as candidates in `BlockSplitSearch`
  (#160). Static; does **not** contain the ⟨6,17,17⟩ identity, and hand-curation
  doesn't scale to the 70-format gap.
- **`DisjointSumSearch`** — beam search by ω_eff over shape-multisets covering
  the target. **Crucial limitation (its own javadoc):** it enforces only the
  **area-cover** necessary conditions (`Σnᵢmᵢ ≥ nm`, `Σmᵢpᵢ ≥ mp`,
  `Σnᵢpᵢ ≥ np`), which are **far from sufficient** → over-optimistic, NOT
  constructive (returns an unrealisable ⟨17³⟩≈2395, ⟨4,4,4⟩=28). So it can
  *propose* the multiset but cannot *certify* it.
- **`Lineage.DisjointSum(children, taLegs)`** — the lineage node already exists
  (renderers/parse/replay wired), so a certified identity has a home to record.

**The missing capability = constructive certification**: given a candidate
shape-multiset (+ axis partition), prove there is an actual trilinear embedding
realising the target at the summed rank — and emit the explicit factor matrices.

## Proposed search — the "grid-minus-k" family first

This sub-family covers ⟨6,17,17⟩ and most of the 70 `⟨a,b,17⟩` gaps, and is far
more tractable than the general τ-theorem:

1. **Enumerate axis partitions** of `(n,m,p)` into part sizes that are catalog
   atoms (e.g. `17 = 6+6+5`, `17 = 5+6+6`, …; `6 = 3+3` or single recursion
   factor 3). Bound the number of parts per axis (≤4) to keep it finite.
2. **Build the full block grid** = Cartesian product of the three partitions →
   the naive additive decomposition `Σ ⟨nᵢ,mⱼ,pₖ⟩` (rank = Σ block ranks; this
   is the honest UB baseline, e.g. 1300 for ⟨6,17,17⟩).
3. **Search for a trilinear-trace identity that drops blocks.** The target
   trilinear form is `T = Σ_{i,k} (Σ_j A_{ij}B_{jk}) C_{ki}` (trace form). For a
   grid, `T = Σ_blocks T_block`. Look for a subset S of blocks + correction
   terms such that `Σ_{S} T_block ≡ T` as trilinear forms (Strassen's disjoint-
   sum / Schönhage τ "the joint computation of the kept blocks embeds the
   target"). This is an **exact linear-algebra feasibility** problem over the
   trilinear coefficient space — the "TRUE constraint is algebraic" that
   `DisjointSumSearch` currently skips.
4. **Construct + verify.** When a feasible S is found, assemble the explicit
   U/V/W from the kept blocks' factor matrices + the identity's combination
   coefficients, and run `Verifier.isExactNonCubic`. Only verified results are
   reported (constructive UB); a passing-area-cover-but-uncertified multiset is
   a **heuristic, never a rank** (optimality discipline).
5. **Lineage:** emit `Lineage.DisjointSum(children = kept blocks, taLegs = the
   shared/dropped legs)` so the win is replayable like any other scheme.

## Acceptance / milestones

- **M0 (test fixture):** reproduce FMM's ⟨6,17,17⟩=1106 from our held
  ⟨3,·,·⟩ blocks — verified exact. This is the go/no-go for the certifier.
- **M1:** sweep the 70 maxDim≤17 FMM-better formats; report how many the
  grid-minus-k certifier closes (vs the area-cover heuristic's optimistic count,
  which must be reported separately and labelled non-constructive).
- **M2:** general τ-theorem (arbitrary sub-tensor embeddings, Pan TA pairs) —
  larger; defer until M0/M1 prove the certifier.

## Cost & honesty

- The certifier is the expensive part (trilinear feasibility per candidate
  subset); cap with a node budget and **log what was dropped** (no silent
  truncation). Report `optimality: bound` for certified results;
  area-cover-only predictions are **not** emitted as catalog ranks.
- Ties into ROADMAP **#196 (Unified MultisetFrontier, TA-aware)** and **#102
  (DisjointSumSearch feasibility)** — this doc is the concrete plan for the
  *certification* half they both presuppose.

## Cross-refs

- Memory `project-17x17x17-search-gap` (first concrete τ observation).
- `references/CLOSURE_OPTIMALITY_AND_PERF.md` (optimality tiers).
- `KnownTauIdentities` javadoc (the two identity families: plain disjoint sum
  vs Pan borrow-and-correct).
