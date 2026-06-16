# Future work

Open research/engineering items deliberately deferred. Each entry: what
it is, why it's open, what closing it would change, and what would be
needed to close it.

---

## Engineering: catalog integrity

### Broken-lineage stubs (filename shape ≠ replayed shape)

**What.** Many `derived_recursive-*` stubs replay (via `LineageReplayer`)
to a DIFFERENT shape than their filename declares — e.g.
`derived_recursive-24x27x27_m9000_a940584.json` → ⟨18,27,36⟩,
`27x27x32_m12330` → ⟨18,27,48⟩, `21x21x30_m7040` → ⟨15,21,42⟩.
`ProjectFmmGaps`'s scatter logged **18+ distinct** such stubs in the
first minutes of one pass. They're skipped safely (resolveParentUncached
falls through to the next candidate, or returns null), so they don't
corrupt results — but each is (a) a catalog-integrity defect and (b) a
*lost projection parent*: a correctly-stored ⟨24,27,27⟩=9000 could be a
strong projection source we're currently throwing away.

**Root cause (pinned 2026-06-08).** Pre-2026-06-04 `materialise`/concat
path recorded concat operands as the **sorted** shape ref
(`canon()+"-direct"`) with **no `OrientAs`** wrapper — e.g.
`ConcatCols(Atom("9x10x12-direct"), Atom("9x10x12-direct"))` for a
declared ⟨10,12,18⟩, which replays to ⟨9,10,24⟩ (operands ⟨9,10,12⟩,
concat p: 12+12). The forward bug was **fixed on 2026-06-04** (`1c1e31b8f`):
`materialise` now wraps a reoriented operand in `Lineage.OrientAs(leaf,
n,m,p)`, so new concat stubs reconstruct their declared shape. The ~150
phantom stubs are all pre-fix (first written 2026-05-31, `72b811b08`).

**Mitigations shipped (2026-06-08).**
1. `EnrichSchemeMetrics --revalidate-stubs` re-replays every stub
   (bypassing the projection_margin fast path) and stamps
   `corrupted:true` on shape-mismatch → manifest `corrupted_files` →
   the lookup skips them at index time. Run it + regen the manifest to
   clear the existing rot.
2. **Write-boundary guard** in `RecursiveMaterialiser.writeToDisk`:
   refuses to persist a stub whose lineage doesn't `replaysConsistently`
   to its declared shape+rank — so any future emitter regression is
   caught at write time, not discovered weeks later. (The projection
   paths already discarded un-replayable results; this extends it to
   the compose path.)

**Still open.** (a) A non-flaky regression test for the guard / the
OrientAs round-trip (deferred — needs a clean catalog, i.e. run after
the revalidate flags the broken stubs). (b) Confirm `LineageTrackingLookup`'s
sorted-`canon()+"-direct"` leaves are orientation-safe in their only use
(RecombinationN leaves); fix to native shape if any recombination stub
turns up mis-shaped.

---

## Engineering: orientation-canonical content hash

**What.** `SchemeIO.contentHash(alg)` is orientation-SENSITIVE: a scheme and
its axis-permutations have different matrices → different hashes. So a content
hash identifies *a scheme in a specific orientation*, not *a scheme up to
orientation*. Hash-refs (`{n}x{m}x{p}@{hash}`) therefore depend on which
orientation was hashed at emit time vs at resolve time. The two disagreed: the
registration/`findByHash` convention hashed the STORED form, but the
`SerendipitousProduct` emitter hashed the base AFTER orienting it to the
product frame (e.g. ⟨6,5,5⟩ while the file is stored ⟨5,5,6⟩). The ref
`6x5x5@4eb0fa` never matched a stored-form scan → fell back to the rank-best
bud-FREE sibling → `productViaBuds` threw `NoSuchElement` → the write-guard
discarded the result → ⟨18,20,30⟩/⟨24,24,27⟩ were unfillable.

**Mitigation shipped (2026-06-09).** `findByHash` now matches EITHER the stored
or the (n,m,p)-oriented content hash (backward-compatible; stored still wins).
This unblocked the last two MISSING shapes (now filled at the better
serendipitous rank, one beating FMM).

**Principled fix (open).** Hash the scheme in its **sorted-shape canonical
orientation**, so all six orientations share ONE identity and hash-refs are
orientation-proof by construction — removing the check-both workaround and any
future emit/resolve orientation mismatch. Cost: re-stamp every existing `hash`
field + hash-ref (a migration); the `4`-char short-hash filenames/refs would
also change. Defer until a hash-ref scheme actually needs the stronger
guarantee — the check-both fix covers correctness today.

---

## Engineering: pin the recombination assignment (task-#91 residual)

**What.** `RecombinationN(base, allocA, allocB, allocC, leaves)` records the base
(now hash-pinned, 2026-06-09) and the allocation *sizes*, but NOT the
**assignment** — which allocation block maps to which base axis-index. The
allocation optimiser picks the rank-minimising assignment (e.g. for ⟨3,21,21⟩ via
base ⟨2,3,3⟩ with `alloc=[2,1]`, it puts the size-1 block on the A-index that makes
more products collapse to the cheap naive ⟨1,7,7⟩=49 instead of ⟨2,7,7⟩=76,
reaching r=978). Replay re-derives the assignment positionally → maps every
product to ⟨2,7,7⟩ → a valid but higher-rank scheme → the write-guard discards it.
Observed rate ≈ 0.4 % of recombination wins (2/458 in T2); safe (no phantom), but
a missed win.

**Why it's the same family as task #91.** Both are "the lineage under-specifies the
recombination". The base-orientation half is fixed (hash-pinned base + precise
`AxisPermute`); this is the *assignment* half — also affected by the base's internal
axis-orientation the optimiser may have used.

**Fix.** Have `recombineWithAllocation` return the chosen per-product assignment
(the allocation→base-index permutation), add it to the `RecombinationN` node, and
have `applyRecomb` reconstruct with that fixed assignment instead of re-optimising.
The `[replay-diag]` WARN already surfaces any build-vs-replay rank mismatch.

---

## Search: cascade ordering, and the (still-open) allocation symmetry

**Cascade is the real gate on the ⟨24,*⟩ tail (measured 2026-06-09).** ⟨24,26,29⟩
via base ⟨2,3,3⟩ has its OPTIMUM at the *balanced* split (`[12,12]|[9,9,8]|[10,10,9]`
= 9576 = FMM), which the B&B incumbent seed finds at **node 0** (verified with
`--stagnation=1`). So the budget / lower-bound were never the blocker for it. T2
missed it only because it's the TOP gap → processed FIRST, before its sub-shapes
(⟨12,9,9⟩, ⟨12,8,10⟩, …) were improved by later passes. **Fix = a second / cascade
pass** (or dependency-ordered processing): re-running now lands 9576. (Correction:
an earlier note here blamed an under-pruned B&B; the 2 M nodes that run consumed were
spent *proving* the seed optimal, not *finding* it.)

**Allocation B&B perf (done, partial).** `AllocationOptimizer` now precomputes each
axis's per-product effective-dim vector and dedups identical vectors
(`uniqueAxis`), and early-exits the inner cost sum once it passes the incumbent —
an exact speedup (no per-node `axisDims`/`shapesAt` recompute). BUT measured: the
dim-vector dedup yields **no space reduction** in practice (per-axis vectors are
distinct even for the symmetric ⟨2,2,2⟩ base).

**Open: real allocation symmetry.** The genuine symmetry (the user's ⟨2,3,3⟩
`3!`-orderings + cross-axis allocB↔allocC for axis-symmetric bases) needs the
**base's index-stabiliser** to canonicalise allocations modulo the group — not the
per-axis dim-vector dedup. Compute the stabiliser per base, enumerate allocations up
to it. Would shrink the heavy-base spaces ~6–36× when it applies.

---

## Search: composite reduction sweep + serendipitous concat (cross-seam merging)

Two related mechanisms surfaced in discussion (2026-06-12), parked for later
investigation. Both are *reduction* passes over composed schemes — the
machinery overlaps; the seeds differ.

### A. Composite reduction sweep (the FMM "kin-row" mechanism)

**What.** After every Kron/concat/cascade materialisation, run a pass that
hunts proportional or linearly-dependent product pairs (`uᵢ⊗vᵢ ∥ uⱼ⊗vⱼ` up to
scaling, possibly after a GL change of basis) and merges them, followed by DCE.

**Why it matters.** Our own gap analysis concluded FMM's headline numbers like
⟨27,28,28⟩=10413 (vs our 10442 raw composition) come from a *global* row
reduction over a big composed scheme, with all local mechanisms falsified
(see the FMM-projection-gap notes). Composed schemes are exactly where such
redundancies accumulate. This is the one mechanism in the 17–32 WORSE band
that our closure verifiably lacks and that demonstrably produces wins.

### B. Serendipitous concat — cross-seam u-collision merging

**The key observation.** In a concat along `p`, the two sides share exactly
one axis — the A-operand. So the only "currency" for sparing products across
the seam is **u-collisions**: a side-1 product and a side-2 product with
proportional u-vectors. That is a bud (Perminov Def 2.9) *straddling the
seam*. `v`/`w` cannot collide across the seam (disjoint p-axis supports).

**Payoff 1 — cross-seam buds for downstream Kron (bookkeeping only).** Run
the bud recognizer on the *concatenated* scheme: cross-seam u-pairs appear as
`⟨1,1,2⟩` blocks and improve the serendipitous cost `r_s` when the concat
result is later tensored. Relative u-alignment between sides is invariant
under global sandwiching but changes when an *inequivalent variant* is
swapped in for one side → "which variant of S₂ to concat with S₁" is a
search axis; natural metaflip objective `cross_seam_bud_score` (walk side 2
against a frozen side 1). Hooks into #159/#196 bud machinery.

**Payoff 2 — lowering the concat's own rank.** Merge a u-collision pair:
`(u·A)(v₁·B) + (u·A)(v₂·B) → (u·A)((v₁+v₂)·B)` — saves 1 product but injects
a structured garbage term (`P₂` leaking into block-1 cells) that other
products must cancel. A merge is net-positive iff cancellation reuses
existing products. **HK71 is exactly this mechanism, solved to closure for
m=2** (the `max(p,n)/2` saving = the cells that pair across would-be seams;
the bridge identities are the cancellation web). For m ≥ 3, do NOT try to
generalise the HK combinatorics — treat it as **drop-and-repair on a concat
seed**: merge k collision pairs, compute the residual tensor exactly (it
lives in `u ⊗ (cross-block v) ⊗ w`, small and structured), and ask a small
SAT/linear solve whether surviving products + freed budget cover it. Sits
inside the Phase-2 SAT envelope, unlike from-scratch search.

**Where savings can exist.** Near-square concats only — the wide regime
telescopes losslessly (HK boundary-term arithmetic: `max(p,n)/2` is paid
once globally, per-piece under concat only when a piece is narrower than p).

**Benchmarks (known answers).**
- Smoke test: ⟨3,3,3⟩=23 vs concat ⟨3,3,2⟩+⟨3,3,1⟩=15+9=24 — one merge.
- Canonical: Smirnov ⟨3,3,6⟩=40 vs self-concat 2·⟨3,3,3⟩=46 — six merges;
  self-concat is the best case (every u-form collides with its twin).
- Inventory: the WORSE-band list (shapes where a published rank beats our
  concat prediction).

**Cheapest first step (do this before any solver).** A `docs.explore` scan
driver (e.g. `ScanConcatSeamCollisions`): for each WORSE-band shape, take our
best concat decomposition, bucket products by canonical u-direction across
the seam, report the disjoint collision count. That count is a **hard
ceiling** on the merge saving (each merge saves exactly 1) — it ranks targets
for free, and near-zero counts falsify the cheap version early (savings would
then require non-parallel, HK-bridge-style 3-way identities instead).

**Relation A↔B.** Same reduction machinery; concat seeds are the
highest-yield composites because operand-sharing *guarantees* collisions,
while Kron composites share multiplicatively and rarely leave parallel pairs
(their savings route through buds — #159). Adjacent but distinct from the
τ-identity certifier (#198), which targets disjoint-sum sharing, not seam
merging.

---

## Engineering: closure scheduling

### Frontier-multipass projection closure (shape-dependency graph)

**What.** Today's projection sweep (`ProjectFmmGaps` / `RecursiveMaterialiser.projectScatter`)
runs **flat multipass**: every pass re-scans the whole target list. The
intra-run redundancy is real — within a pass parents are static, so a
shape already at its projection-closure can only move when one of its
**parents** improves (a *cascade*). The fix is to treat projection as a
**dependency graph over shapes**: a node `⟨n,m,p⟩` depends on its parent
shapes `⟨n+δ,m+δ,p+δ⟩` (the projection edges). After pass *N* produces
wins at some shapes, only the **frontier** — the projection-children
*downstream* of those wins — needs re-evaluation in pass *N+1*, instead
of rescanning all 2061 shapes.

This is exactly the **impact-graph** concept already in the repo
(`evaluate-up, materialize-down`): a win materialised at a node
invalidates the cached results of its downward neighbours. Generalised,
the same graph carries the **composition strategies** too — Kron / concat
/ recombination / serendipitous edges, not just projection — so a single
impact-graph could schedule *all* closure operators by dirty-propagation
rather than blind full passes.

**Why it's open / deferred (user 2026-06-08).** Tricky to code correctly:
the edge set is per-strategy (projection edges are `±δ` per axis;
Kron/concat edges are factor-pair products), cascades can re-enter, and
the dirty-set must survive across passes without masking a genuine retry
when **code** changes (a still-behind-FMM shape must stay retryable when a
new construction lands — see the "don't foreclose future wins" rule). The
safe scope is **intra-run only**: across runs we deliberately retry the
full WORSE list with the latest code.

**What closing it would change.** Pass 2+ would touch only the cascade
frontier (handful of shapes) instead of all 2061 — a multipass closure
would cost ≈ pass-1 + ε rather than k × pass-1. No effect on *what* is
found, only the cost of finding cascades.

**What's needed.** A per-shape "projection-closed @ rank R against
parent-state H" marker (H = hash of the parent shapes' current ranks),
plus a dirty-propagation scheduler keyed on the projection (and later
composition) edge set. Margin-prune + parent-cache (shipped 2026-06-08)
already make each *necessary* attempt cheap; this only removes the
*unnecessary* re-attempts within a run.

---

## Mathematical derivation gaps

### Hopcroft-Kerr 1971: missing same-method-pair derivations

**What.** HK page 10 derives a same-method-pair formula
`(y_{ii}=method-1, y_{jj}=method-1, bridge y_{i+1,i+1}=method-2)`
— shorthand `(1,1,bridge-2)`. It claims "the other cases follow by
symmetry" but never derives them. There are 6 such cases:

- `(1,1,bridge-2)` — paper-explicit, implemented in `emitSameMethodPair_11_bridge2`
- `(2,2,bridge-1)` — implemented in `emitSameMethodPair_22_bridge1` (sympy-derived)
- `(1,1,bridge-3)` — implemented in `emitSameMethodPair_11_bridge3` (sympy-derived)
- `(2,2,bridge-3)` — **OPEN**. Current code looks up products (`E_adj`, `F_adj`) that the underlying Lemma-2 Case (2,3) bridge pair doesn't actually emit.
- `(3,3,bridge-1)`, `(3,3,bridge-2)` — also undefined; don't arise in Case-1 odd-n alternating coloring so paper non-issue, but Case 2 (even p) Step 3 boundary handling reaches them.

**Why it's open.** Sympy enumeration over the natural atom catalog
(shift-form products `c_i·a_i + c_p·a_{i+1} + c_j·a_j` with ±1 coeffs,
combined with method-style products A/B/C/D/E/F) returned no
3-new-product solution for `(2,2,bridge-3)`. Whether one exists in a
wider catalog (rational/larger coefficients, alternate shift bases) is
genuinely unknown — paper is silent.

**Closing would unlock.**

- HK Case 1 (odd p) at HK-formula rank `(3pn+n)/2` for all `p ≥ 7` (currently capped at `p ∈ {3, 5}`).
- HK Case 2 (even p) Steps 1-3, which targets the 3 FMM-gap shapes
  `⟨2,10,15⟩`, `⟨2,10,16⟩`, `⟨2,12,16⟩`.
- Best-case rank improvement vs current catalog: **1-4 mults per shape** on those three. FMM-Lille / Perminov already publish 234/249/298 (1-4 above HK formula), so the marginal gain is tight.

**What's needed to close.**

- Either: a careful manual analysis of HK's intended "symmetry" claim — find an explicit isomorphism mapping the `(1,1,bridge-2)` derivation onto `(2,2,bridge-3)`. The methods 1↔2 and 1↔3 swaps under `x_1↔x_2` and `a_1↔a_2` don't cleanly combine, but a layered substitution might. 4-8 hours of careful algebra + sympy verification; ~15% expected success.
- Or: extended brute-force enumeration with coefficients ∈ {-2,-1,0,1,2} and asymmetric shifts on `sa_*` vs `sx_*`. Hours of compute; ~10% expected success.
- Or: prove impossibility — show no 3-product solution exists in any reasonable atom space, accept +1 mult per bridge-3 pair.

**Workaround in current code.** `HopcroftKerr2bcAsymmetric.build(p, n)`
dispatches to `buildOdd(p, n)` (HK-optimal) for `p ∈ {3, 5}`, else to
`buildNaiveDCE(p, n)` (sub-optimal but verified).

---

## Catalog/engineering gaps

### Recursive materialiser for derived bounds (task #32)

`search.findBestStrategy` returns a *rank prediction* for shapes like
`⟨24,24,24⟩=7200` without writing a scheme to disk. FMM-Lille writes
them. Closing this means making the strategy chooser callable as a
materialiser — recursively dispatch each picked sub-strategy through
`Compose.kroneckerGeneral` / `Compose.concatRight` /
`Recombination.constructWithAllocation`, write the result via
`SchemeIO.write` (cache auto-syncs). Expect ~30+ new schemes in the
(24, 32) range.

### Lineage back-fill (task #44)

`MaterializeViaPanPair`, `MaterializeViaPanPair14`,
`MaterializeSolvenStrassen{New,21,777}`, `MaterializeRosowskiAlgorithm1`,
`MaterializeClosureLoop` — none emit a `Lineage.Node` yet. Pattern
established by `MaterializeTriple18` is mechanical; each materialiser
needs to construct the right tree. Once #32 lands, the recursive
materialiser can emit lineage automatically.

### Lineage on direct catalog hits is shallow

When `RecursiveMaterialiser.materialise(n, m, p)` finds the shape in
the catalog directly (no recursion needed), it returns
`new Lineage.Leaf(canon-direct)` — a flat string referencing the
canonical key. If the catalog leaf has its OWN `lineage` field
(because it was itself composed earlier), that deeper history is
lost. As a result, deep composition trees collapse to shallow ones:
the new ⟨18,18,18⟩=3200 we materialised is
`ConcatRight(18x18x9-direct, 18x18x9-direct)` rather than
`ConcatRight(KronProduct(fmm-lille_3x3x6_r40, …), …)`.

To fix: change the direct-hit branch in `materialise()` to
(a) look up the actual filename via `FieldAwareLookup`,
(b) parse the file's `lineage` field if present, and
(c) return that as the `Result`'s lineage instead of a flat Leaf.
Should also extract the actual filename rather than the canonical
sorted key — `2x9x18-direct` lineage entries are misleading when the
catalog file is at a different orientation.

Affects all schemes the recursive sweep writes; the more iterations
the sweep runs, the more lineage history is silently lost.

### Commutative asymptotic exponent ω_c

**Conceptual gap I had wrong earlier.** Rosowski 2019's central
contribution is **not** a commutative→NC bridge. It's that
commutative schemes recurse over commutative rings, so a commutative
cubic `⟨n,n,n⟩` at rank `r` gives `ω_c ≤ log_n(r)` — a *separate*
asymptotic exponent from NC `ω`. Example: Rosowski Thm 2/3
`⟨3,3,3⟩_c = 21` → `ω_c ≤ log_3(21) ≈ 2.771`, *better* than
Strassen's NC `log_2(7) ≈ 2.807`.

Islam 2009 (MSc thesis) discusses this asymptotic-from-commutative
framework — exact treatment TBD (digging through
`references/papers/islam_2009_msc_optim_matmul.pdf` deferred).

**Implications for the catalog/SPA:**
1. Add a commutative `ω_c`-history axis to the ⟨∞,∞,∞⟩ section —
   tracked in ROADMAP as "Commutative ω-history axis".
2. The `commutative: true` flag is enough for asymptotic
   aggregation — no conversion primitive needed.
3. Commutative-only schemes contribute to `ω_c` progress even
   though they don't lift to NC matmul over NC rings.

### Commutative → NC conversion primitives (low priority)

Separate question: turning a commutative `⟨a,b,c⟩` scheme into a NC
`⟨a,b,c⟩` (or `⟨a,2b,c⟩`) scheme. Known patterns:

1. **Trivial 2× lift**: `R_nc(⟨a,b,c⟩) ≤ 2 · R_c(⟨a,b,c⟩)`
   (compute `AB` and `BA` independently). Weak; rarely tight.
2. **Symmetrization** (`⟨a,b,c⟩_c → ⟨a,2b,c⟩_nc`): doubles the b-axis;
   no one-line formula I've seen.
3. **Rosowski Thm 4 transpose involution**: `R_nc(⟨n,n,n⟩) ≤
   n(n²+3n+1)/2`. Standalone NC construction (not a conversion);
   `papers.rosowski2019.RosowskiAlgorithm1` materialises it for
   `⟨n,3,3⟩` (31 files on disk).

Demoted from "needed" to "nice-to-have" — the commutative SOTA
contributes via `ω_c` directly, without needing to be projected
into NC.

### Commutative recursive materialiser

`RecursiveMaterialiser` currently uses a non-commutative leaf pool
(Strassen + Laderman). A parallel pass should run with a
**commutative** leaf pool: `{Waksman 1970, Rosowski Thm 2/3,
Makarov 1986, Islam 2009}`. Compositions (Kron, concat,
constructWithAllocation) preserve commutativity, so this pass writes
commutative-only schemes covering the (a, b, c) grid with their
own `c_F` improvements. Schemes are tagged `"commutative": true` and
stay separate from the NC catalog (commutative algorithms don't lift
to recursive matmul over NC rings, per RANK_KNOWLEDGE.md §1.2bis).

**Missing pieces**:
- `Waksman 1970` generator + materialiser (`MaterializeWaksman1970`)
  exist; 31 files on disk under `waksman-1970_{n}x{n}x{n}_r{r}_a{adds}_commutative.json`.
- `Rosowski Algorithm 1` materialiser exists
  (`papers.rosowski2019.RosowskiAlgorithm1`); 31 files on disk.
- **Need**: Islam 2009 materialiser (currently only the bound formula
  is registered); Rosowski Thm 2/3 materialiser (above);
  Makarov 1986 materialiser (above).
- **Need**: A `FieldAwareLookup` configured to filter by
  `"commutative": true` (so the commutative sweep can find commutative
  leaves).
- **Need**: A commutative-aware `BlockSplitSearch.defaultPool()`
  variant returning the commutative leaf set.

The expected output is a parallel set of `composed-recursive-cmt_*`
scheme files documenting the commutative SOTA across shapes — useful
for context (showing how much commutativity buys vs. NC), even though
they don't feed back into the NC recursion.

### Materialise Rosowski Thm 2/3 ⟨3,3,3⟩ = 21 (commutative)

`docs/derived-bounds.json` already lists `⟨3,3,3⟩ = 21` attributed to
Rosowski 2019 — but only as a rank claim. The explicit factor matrices
have not been emitted to `src/main/resources/schemes/section3/`.

The Rosowski 2019 paper (arXiv:1904.07683) gives Thm 2 + Thm 3 as
commutative bilinear constructions for `⟨n,n,n⟩` beating Waksman's
generic formula for small `n` (Thm 3 specifically targets `⟨3,3,3⟩ =
21`). Closing this means writing a `papers.rosowski2019.RosowskiThm2`
or `…Thm3` constructor + `MaterializeRosowskiThm2`, mirroring the
pattern in `MaterializeWaksman1970`.

Distinct from `RosowskiAlgorithm1` (already materialised, 31 files)
which targets `⟨n,3,3⟩` via the transpose involution and is
non-bilinear.

### Materialise Makarov 1987 ⟨5,5,5⟩ = 100 (non-commutative)

Separate from the Makarov ⟨3,3,3⟩ above (different paper).
Constructive recipe (Islam 2009 §3.2.5):

| Sub-product | Rank | Source | Status |
|---|---:|---|---|
| ⟨3,3,3⟩ | 23 | Laderman 1976 | ✓ on disk |
| ⟨2,3,3⟩, ⟨3,3,2⟩, ⟨3,2,3⟩ | 3·15 = 45 | Hopcroft-Musinski 1973 | ✓ via orientAs from ⟨2,3,3⟩=15 |
| ⟨2,2,2⟩ | 7 | Strassen 1969 | ✓ |
| ⟨2,2,3⟩ | 11 | Strassen recursion | ✓ |
| ⟨3,3,2⟩ *with one zero entry* | 14 | modified HM | **blocker** |
| **Total** | **100** | | |

6 of 7 sub-products are already in our catalog. The blocker is the
"⟨3,3,2⟩ with one zero matrix entry → 14 mults" specialisation — a
non-generic ⟨3,3,2⟩ scheme that exploits the known zero. Closing:

1. Locate the explicit construction in Hopcroft-Musinski 1973 (SIAM J.
   Comput. 2(3)), Makarov 1987, or Islam 2009 §3.2.5; OR derive it
   via sympy by symbolic elimination of one zero from HM ⟨3,3,2⟩=15.
2. Implement `MaterializeMakarov555` that glues the 7 sub-products at
   the right partitions of ⟨5,5,5⟩.

`ω` implication: log_5(100) = 2.862 — historical baseline, not SOTA
(AlphaEvolve ⟨5,5,5⟩=93 over Z gives 2.817; AlphaTensor F₂=96 gives
2.835). Materialising preserves chronology per CLAUDE.md.

### HK n > 2p decomposition (task #49)

Recursive split of `⟨p, 2, n⟩` when `n > 2p` into chunks. Doable but
the catalog's `ConcatSplitSearch` already finds these splits and they
don't beat known catalog ranks, so the engineering payoff is small.

---

## Conventions

If an item makes it to "closed" status, move it out of this file. If
new gaps emerge during work, append them here rather than letting them
evaporate.
