# Rectangular trilinear aggregation (rectangular Pan TA) — scoping (b)

*(2026-06-13 scoping → 2026-06-14 IMPLEMENTED. Status: DONE — the construction is
derived, exact-verified, and wired into decomposition + materialisation. See
"INTEGRATION — DONE" below. Do NOT ship a TA construction unless it is
exact-verified and replays bit-exactly — that discipline is now enforced at the
write boundary.)*

## Why

FMM closes prime-dimension shapes by **peeling** (a prime axis can't be
Kronecker-tiled) into a `⟨1,2,2⟩` block grid and **fusing the two off-diagonal
cross-blocks**, which are cyclic rotations, via Pan trilinear aggregation (TA).
Canonical case (`analyse-fmm-html-not-tau`, this session):

```
⟨26,29,29⟩ = 11693 = ⟨26,26,26⟩ 8658 + ⟨26,3,3⟩ 175 + TA(⟨26,3,26⟩, ⟨26,26,3⟩) 2860
```

The two cross-blocks `⟨26,3,26⟩` and `⟨26,26,3⟩` are cyclic rotations of
`{26,3,26}`; naive `1504+1504 = 3008`, TA-fused `2860` (−148). CORRECTION (2026-06-14): our `⟨26,3,3⟩=159` is the **commutative** Rosowski
scheme — INVALID as an NC leaf. The non-commutative corner is **175** (exactly
FMM's). So TA gives `8658 + 2860 + 175 = 11693` = **FMM exactly (a TIE, not a
beat)** — but still improves OUR catalog (prior 11808 → 11693). Lesson: compare
NC-to-NC; a commutative rank is not a valid leaf (cf. the ⟨17,17,17⟩=2868
commutative-leaf incident). The `(a)` map shows ~136 TA-class holdouts.

## What we have vs need

- **Have (cubic only):** `PanTrilinearAggregation.build(n)` →
  `PanTrilinearAggregationBuilder` — a verified port of Islam's Magma square TA,
  parameterized by ONE side `nn` with `(n+1)`-padding. It computes the cyclic
  triple of `⟨n,n,n⟩` jointly. `cubicBound(n)` is its rank.
- **Need:** the **rectangular** generalization — fuse a cyclic-rotation PAIR
  `⟨a,b,c⟩ + ⟨c,a,b⟩` (or the relevant pair of the triple) for unequal `a,b,c`,
  with (1) an exact **rank** and (2) an exact **factor-matrix construction**.
  The square builder does not generalize mechanically — the padding and the
  aggregation forms assume `a=b=c`.

## Open derivation (do this math-first, then verify)

1. Write Pan's TA bilinear identity for a disjoint pair of cyclic-rotation
   products with unequal dims; confirm the cross-term cancellation that yields
   the saving, and read off the **rank** as a function of `(a,b,c)`.
2. Confirm against the known instance: the pair `⟨26,3,26⟩,⟨26,26,3⟩` must come
   out at FMM's **2860** (or better).
3. Only THEN build factor matrices; **exact-verify** (`isExactNonCubic` /
   spot-check) at small dims before trusting large ones.

Sources to consult (ask the user for PDFs per the PDF policy): Pan 1978 TA,
Hopcroft–Kerr 1971 (pairwise rotation aggregation), Islam 2009 (our cubic port).

## Integration design (once the construction exists)

Per the user (2026-06-13), TA belongs **inside the recombination evaluator** as
a cost-gated leaf-fusion, NOT a one-off:

- Carriers: add `⟨1,2,2⟩` (+ other multi-axis naive splits) to the pool — a
  symmetric peel emits a disjoint cyclic-rotation cross-pair among its leaves.
- In `Recombination.recombineWithAllocation` (the `Σ sota.getRank(shapes[r])`
  leaf-sum), detect a **disjoint cyclic-rotation leaf-pair** and price it at the
  **fused TA rank** instead of the sum — **only when it beats the sum**
  (cost-gated), and **exact-verify** the fused block before keeping it.
- This is finite Pan TA, NOT the asymptotic τ-theorem
  ([[feedback_we_do_not_leverage_tau_theorem]] — finite TA is a real,
  implementable mechanism; the asymptotic τ is the pointless one).

## Construction progress (2026-06-13)

Verified building blocks in `RectangularTrilinearAggregation` (+
`TestRectangularTrilinearAggregation`, 5/5):

1. **Identity (3) atom** — `aggregatePair(a,b,c, a',b',c')` returns the 4-term
   `(a+a')(b+b')(c+c') − a'(b+b')c − a·b'(c+c') − (a+a')·b·c'`; bit-exact verified.
2. **Disjoint assembly** — `aggregateDisjoint(p1,…,p2,…)` embeds the two
   cross-blocks in a combined space, pairs term-by-term, aggregates each via (3);
   verified to compute `p1 ⊕ p2` exactly. Currently `4·r` terms (a LOSS).

**Key insight — uniting ≡ a two-factor merge.** Pan's uniting identity (4) sums
correction terms that share TWO factor-directions and differ in the third:
`u⊗v⊗w₁ + u⊗v⊗w₂ = u⊗v⊗(w₁+w₂)` — exactly the doubly-proportional merge we
already have (`FlipScheme.reduce` / the bud machinery). So the construction is
**aggregateDisjoint → merge**, and the merge is exact/safe.

**KIN-PAIRING — DERIVED & VERIFIED (2026-06-14).** For the DISJOINT cross-pair
the kin-pairing is **just the rotation** `T1(i,j,k) ↔ T2(k,i,j)` — **no shift**
(Pan's cubic shift was only to get the even/odd perfect matching within ONE LA;
the disjoint sum pairs P1→P2 bijectively for free). The three correction
families then unite naturally because each fixes two indices:
- `a'=A2(k,i)`, `c=C1(i,k)` fix `(i,k)` → corr1 unites over `j` → `n·p` terms
- `a=A1(i,j)`, `b'=B2(i,j)` fix `(i,j)` → corr2 unites over `k` → `n·r` terms
- `b=B1(j,k)`, `c'=C2(k,j)` fix `(j,k)` → corr3 unites over `i` → `r·p` terms

giving the rank formula

```
r_TA(n,r,p) = nrp + np + nr + rp        (fuse ⟨n,r,p⟩ ⊕ ⟨p,n,r⟩)
```

`⟨26,3,26⟩`: `2028 + 676 + 78 + 78 = 2860` = FMM's TA term, exactly. Implemented
in `RectangularTrilinearAggregation.{fusedRank,build}`; `build` is DETERMINISTIC
in `(n,r,p)` and verified **bit-exact** against the disjoint-sum tensor
(`TestRectangularTrilinearAggregation`, 6/6). Coefficients are ±1 (integer).

**Replay (user requirement):** `build` is a pure function of `(n,r,p)`, so the
lineage node need only record `(n,r,p)` + the two cross-block leaf refs (their
placement comes from the parent recombination's allocations) — bit-exact by
re-running `build`.

**INTEGRATION — DONE (2026-06-14).** TA is wired into decomposition (not a
top-level strategy), driven via `SchemeSweep --base=1x2x2`:

- **Scoring** (`BlockSplitSearch.findBestStrategy`): the default unbalanced /
  no-peel recombination path used the pairing-BLIND allocation optimizer, so the
  Pan saving was invisible (it even missed the `[N,s]` peel allocation, returning a
  balanced split). Fix: alongside the optimizer, run the pairing-aware mask sweep
  (`findBestMultiBaseSplit`, which prices cyclic pairs via
  `PairedSubProducts.applyPairing` = `fusedRank`) **over the width-1 CARRIER bases
  only** — zero cost when the pool has none (the common Strassen/Laderman case;
  cubic same-shape pairs stay on the balanced `PairFusedRecombination` path),
  bounded by `PAIRING_SWEEP_COMBO_BUDGET`. Keep whichever recombination is cheaper.
- **Materialisation** (`RecursiveMaterialiser.tryBuildTaPeel`): a recombination
  scored with pairing was being BUILT by gluing the cross-pair independently (→
  un-fused rank, rejected by the improve gate). The symmetric `⟨1,2,2⟩` peel now
  routes through `RectangularTrilinearAggregation.buildPeeledViaTa`, emitting a
  `Lineage.PeeledViaTa(N,s,cube,corner)`. The leaves are pinned to the
  NON-COMMUTATIVE best (`FieldAwareLookup.findFilesNonCommutative` — the corner
  must be ⟨26,3,3⟩=175, NOT the commutative 159 that the rank-best resolveShape
  would pick); the alg is built BY replaying the recorded tree, so build ≡ replay
  and the write-boundary `replaysConsistently` check cannot reject it on a
  leaf-resolution mismatch.

Verified end-to-end: `⟨26,29,29⟩ = 8658 + 2860 + 175 = 11693` (FMM tie; prior
11808) and `⟨28,31,31⟩ = 14043`. Sweep over all `⟨n,m,m⟩` (n<m≤32): only these two
beat both the incumbent AND concat — TA is a narrow win, not a broad one. Guard:
`TestTaPeelDecomposition`.

## GENERALISED to ANY naïve-grid base (2026-06-14)

Per the user, TA must be a **generic** mechanism over any base, not the ⟨1,2,2⟩
special case — FMM closes most ⟨n,m,m⟩ gaps with **larger block grids and MULTIPLE
TA fusions**, e.g. (from `fmm.univ-lille.fr/22x28x28.html`):

```
⟨22,28,28⟩ : grid 22→[12,10], 28→[10,9,9]², a ⟨2,3,3⟩ naïve grid
           = ⟨12,10,9⟩+⟨12,10,10⟩+⟨10,10,10⟩+… + 2× TA(⟨10,9,10⟩,⟨10,10,9⟩)
```

Implemented as `Recombination.constructWithTaFusion(base, subResolver, sota,
allocA,allocB,allocC)` → `TaFusedConstruction(alg, fusedPairs)`. For a **naïve-grid**
base (`isNaiveGrid`: `r==n·m·p` and every product is a single block A(i,j)·B(j,l)→C(i,l),
coeff +1), it greedily matches **disjoint cyclic-rotation single-block product pairs**
(distinct A/B/C blocks) and embeds the fused TA block `build(n,r,p)` onto the two
products' global block positions (`embedTaPair` — generalising `buildPeeledViaTa`'s
`taU/taV/taW`); unpaired products embed normally. The fused leaves need NO scheme (TA
replaces them); only UNPAIRED leaves are resolved (lazily, via a caller `SubResolver`
so STUBS replay).

- **Lineage:** `Lineage.RecombinationTaN(base, allocA, allocB, allocC, leaves)` — base
  is the naïve grid (`Atom "naive-NxMxP"`), `leaves` are the **NC-pinned** unpaired
  schemes (one per distinct sorted shape; the fused part is integer ±1, field-neutral).
  Replay (`LineageReplayer.replayRecombinationTa`) re-runs `constructWithTaFusion`,
  re-deriving pairs deterministically and resolving leaves by sorted-shape key. Build
  is done BY REPLAY so build ≡ replay (write-boundary `replaysConsistently` validated).
- **Scoring:** the pairing-aware mask sweep (`BlockSplitSearch.findBestStrategy`) now
  runs over the **naïve-grid** bases in the pool (was: width-1 carriers), bounded by
  `PAIRING_SWEEP_COMBO_BUDGET` — zero cost on the Strassen/Laderman default pool.
- **Driving:** `SchemeSweep --base=2x3x3` (etc.) — `userBasePool` synthesises the naïve
  grid and its S₃ orbit (the ⟨2,1,2⟩/⟨2,2,1⟩ isotropies are kept; the search picks the
  orientation that yields a beneficial cyclic pair).

Guards: `TestConstructWithTaFusion` (bit-exact ⟨2,3,3⟩-via-⟨1,2,2⟩, 2-fusion ⟨22,28,28⟩
over a ⟨2,3,3⟩ grid, `isNaiveGrid`). `PeeledViaTa` remains (older stubs replay via it)
but the materialiser now emits the generic `RecombinationTaN`.

**OPEN (search tuning, not mechanism):** finding FMM's *specific* unbalanced grid
allocations (e.g. `[12,10]/[10,9,9]²`) needs a higher combo budget or FMM-page-guided
allocations — the default budget degrades a ⟨2,3,3⟩-at-⟨22,28,28⟩ enumeration (≈2.6M
combos) to balanced, missing the FMM split. And our catalog's sub-schemes (⟨12,10,9⟩
etc.) must reach FMM's leaf ranks for the total to win.

## Honesty / non-goals

A TA win is only written after the fused block is built and **exact-verified**
(`passesRandomMatmulSpotCheck`) AND replays bit-exactly (`replaysConsistently` at
the write boundary). A `⟨1,2,2⟩` carrier alone (no fusion) just reproduces the
naive 4-block sum that concat+suminner already give — no saving is assumed without
a verified fused block. The general rectangular cross-pair (a non-`⟨1,2,2⟩` base
producing a cyclic pair embedded among many leaves) is still NOT materialised —
only the symmetric `⟨N,N+s,N+s⟩` peel that `buildPeeledViaTa` covers.
