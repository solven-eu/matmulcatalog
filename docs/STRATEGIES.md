# Strategies for finding matrix-multiplication schemes

Reading guide: for any FMM-better gap row, this file tells you which
strategy *would* close it, which strategies are *automatically* tried
during search, and which require manual invocation.

The codebase has accumulated a sizeable set of construction
strategies, lookup augmentations, and search drivers. Without an
index they're easy to confuse. This file is the index.

---

## 1. Construction primitives

These are direct algorithms — given inputs, produce a
`NonCubicBilinearAlgorithm` (or `PairScheme` / `NonBilinearAlgorithm`
variant).

| File | Strategy | Shapes | Field | Status |
|---|---|---|---|---|
| `Strassen7` (legacy) | Strassen ⟨2,2,2⟩=7 | exact | all | ✓ |
| `Compositions.atStrassen` | Strassen × AT-F2 ⟨4,4,4⟩=47 | composed | F2 | ✓ |
| `Compositions.strassenSquared` | Strassen² ⟨4,4,4⟩=49 | exact | all | ✓ |
| `HopcroftKerr2bc.buildSquareOdd(n)` | HK ⟨n,2,n⟩ Case 1 | odd n ∈ [3, 15+] | all | ✓ |
| `HopcroftKerr2bc.buildSquareEven(n)` | HK ⟨n,2,n⟩ Case 2 | even n ∈ [2, 16+] | all | ✓ |
| `HopcroftKerr2bc.buildSquare(n)` | dispatcher | any n | all | ✓ |
| (planned) HK asymmetric `⟨p,2,n⟩` | requires `LemmaOneAugmentation` + band-restricted square + back-sub | p ≤ n ≤ 2p | Q (denominators) | bedrock done (Lemma 1 + Vandermonde); rest pending (#48) |
| `Waksman1970.build(n)` | commutative ⟨n,n,n⟩ via (3n²+n)/2 formula | any cubic | commutative-only | ✓ |
| `Rosowski21.build()` | Rosowski Corollary 1 `⟨3,3,3⟩=21` commutative | exact | commutative | ✓ |
| `RosowskiAlgorithm1.build(n)` | Rosowski Algorithm 1 `⟨n,3,3⟩=6n+3` commutative | exact | commutative | ✓ |
| `PanPairProduct.build(a, b, c)` | Pan trilinear-aggregation pair: joint `⟨a,b,c⟩ + ⟨b,c,a⟩` at rank `abc+ab+bc+ca` | any cyclic-pair | all | ✓ |
| `PanTrilinearAggregation` | Pan TA closed-form rank for cubic ⟨n,n,n⟩ NC | any cubic | all | ✓ (formula only; constructor pending #44 lineage) |
| `HopcroftKerrBound.forShape(a,b,c)` | HK formula `(3ac+max(a,c))/2` for `⟨a,2,c⟩` | any with one axis = 2 | NC | ✓ (formula only) |
| `WaksmanBound.forShape(a,b,c)` | Waksman closed-form for commutative `⟨a,b,c⟩` | any (b axis structured) | commutative | ✓ (formula only) |
| `LemmaOneAugmentation.build(p, n)` | n × p integer matrix; first p rows = identity, every cyclic p-window non-singular | any p ≤ n | independent of scheme field | ✓ (Vandermonde) |
| `Strassen-Winograd ⟨2,2,2⟩=7, a=15` | Winograd 1971 variant of Strassen at 15 adds | exact | all | pending #52 |
| `IslamWaksmanBound` | Islam 2009 Prop 3 formula | rectangular | commutative | pending #51 |

## 2. Composition operators

Build a bigger scheme from smaller ones. These are deterministic given
inputs — no search. Either auto-invoked by a search driver or called
explicitly by a `Materialize*` script.

| File | Operator | Resulting shape | Rank | Notes |
|---|---|---|---|---|
| `Compose.kroneckerGeneral(outer, inner)` | Kronecker product (tensor product of the bilinear forms; Sedoglavic Lemma 13) | ⟨n₁n₂, m₁m₂, p₁p₂⟩ | `r₁·r₂` | Rank-preserving but never strictly improving over factor ranks |
| `Compose.kronecker(outer, inner)` | cubic Kronecker convenience wrapper | cubic | `r₁·r₂` | |
| `Compose.chain` / `chainGeneral` | repeated Kronecker | iterated cubic | `∏ r_k` | Strassen^k etc. |
| `Compose.concatRight(left, right)` | block-column concat on p axis | ⟨n, m, p₁+p₂⟩ | `r₁+r₂` | Cost-additive; for "narrow" shapes with one axis = 2, often within 1-2 mults of HK formula |
| `Compose.concatBelow(top, bottom)` | block-row concat on n axis | ⟨n₁+n₂, m, p⟩ | `r₁+r₂` | Symmetric of concatRight |
| `Recombination.recombineWithAllocation(base, sota, alloc)` | DIS09-style allocation-aware composition; addition-processing shrinks sub-products | varies (depends on alloc) | `Σ R(sub-shape)` over the allocation | The "Strassen[a,b]³ → 7·R(⟨a,a,a⟩)" pattern lives here |
| `Recombination.constructWithAllocation` | constructive variant of the above (returns the actual scheme) | as above | as above | What `MaterializeSolvenStrassenNew` calls |
| `RecombinationWithPair.constructWithPairing(base, lookup, alloc, pairing)` | as above, but fuses N cyclic-equivalent sub-product *pairs* through `PanPairProduct` | as above | `Σ R(unpaired) + Σ PanPair(paired)` | Closes `⟨14,14,14⟩=1719` etc. |
| `ConcatSplitSearch.materialise(split, lookup)` | materialise a concat decision | as concatRight/concatBelow | r₁+r₂ | |

## 3. SOTA resolvers (lookup augmentation)

A `SotaResolver` returns "best known rank for shape ⟨a,b,c⟩". Multiple
implementations exist; choosing the right one for a given search is
load-bearing.

| File | Resolver | Returns |
|---|---|---|
| `FieldAwareLookup` | catalog-only — walks `src/main/resources/schemes/` and returns the smallest rank for the requested field per the inclusion chain Z ⊂ Q ⊂ R ⊂ C; F₂ separate | `Optional<NonCubicBilinearAlgorithm>` |
| anonymous in `Recombination.constructWithAllocation` | trivial wrapper: catalog hit or `Integer.MAX_VALUE/100` | int |
| `FormulaAwareSota(lookup, includeCommutative)` | catalog + Pan TA + HK ⟨·,2,·⟩ formula + (optional) Waksman | int (best across catalog + formulas) |
| `RecursiveClosureSota(lookup, pool, balancedOnly, useFormulaBounds)` | memoised recursive: each query consults catalog, formulas, AND `findBestStrategy` (which itself recurses) — DIS09's T[a,b,c] fill | int (closure-minimum) |

## 4. Search drivers (what gets tried automatically)

| File / method | Strategies enumerated | Sota source | Notes |
|---|---|---|---|
| `BlockSplitSearch.findBestSplit(n)` | u+v splits of cubic ⟨n,n,n⟩ via Strassen | caller-provided | sums 7 sub-products at 3 cubic + 3 mixed shapes |
| `BlockSplitSearch.findBestSplitNonCubic(n,m,p, strassen, sota)` | Strassen ⟨2,2,2⟩ allocations on each axis | caller | non-cubic targets |
| `BlockSplitSearch.findBestMultiBaseSplit(...)` | every (base, allocation) pair in caller's `NamedBase` pool. For each (base, alloc) tuple also tries `PairedSubProducts.applyPairing` on the resulting sub-shape list, takes the min vs unpaired. | caller | `defaultPool() == rootPool()` returns 8 entries (see below) |
| `BlockSplitSearch.findBestStrategy(n,m,p, pool, sota, balancedOnly)` | min of `findBestMultiBaseSplit` AND `ConcatSplitSearch.findBest` | caller | "the headline picker" |
| `ConcatSplitSearch.findBest(n,m,p, sota)` | every 1-axis split on n and on p | caller | |
| **(NOT implemented)** systematic Kronecker enumeration | every `Kronecker(catalog_i, catalog_j)` matching shape | n/a | Task #55 — would let recursive closure auto-discover composite ranks |
| **(NOT in search)** PanPairProduct fusion | only tried inside `RecombinationWithPair` when caller specifies an explicit pairing | n/a | Task #50 — wire into `findBestStrategy` pairwise on cyclic sub-shapes |

### 4.1 Two populations: outer templates vs sub-shape leaves

A recombination search has **two** distinct uses of catalog schemes:

1. **Outer-template bases** — the "block layout" we recurse on
   (Strassen ⟨2,2,2⟩'s 7-product schedule, Laderman ⟨3,3,3⟩'s 23, etc.).
   These are the candidates in `BlockSplitSearch.defaultPool()` /
   `rootPool()`. Adding a base here means we'll try block-decomposing
   ANY target through it.

2. **Sub-shape SOTA leaves** — the cost of each sub-product is read via
   `FormulaAwareSota.getRank(sub_n, sub_m, sub_p)`, which consults the
   FULL catalog (every JSON in `src/main/resources/schemes/section*/`)
   plus Pan TA + HK formulas. So every leaf cost leverages everything.

The Kronecker path implicitly bridges the two: any `(a,b,c) ⊗ (d,e,f)`
factorization is enumerated and both halves go through the SOTA. But
the "outer template" pool is small by design (enumeration grows
combinatorially with base size).

`BlockSplitSearch.rootPool()` currently lists **8 entries**, generated
by S₃-orbit expansion (one representative per distinct shape) of these
historical-root NC Z-arithmetic schemes:

- Strassen ⟨2,2,2⟩=7 (Strassen 1969) — 1 entry, shape-symmetric.
- AT ⟨2,2,3⟩=11 — 3 entries: ⟨2,2,3⟩, ⟨2,3,2⟩, ⟨3,2,2⟩.
- AT ⟨2,3,3⟩=15 — 3 entries: ⟨2,3,3⟩, ⟨3,2,3⟩, ⟨3,3,2⟩.
- Laderman ⟨3,3,3⟩=23 (Laderman 1976) — 1 entry, shape-symmetric.

An `extendedPool()` covering every Leaf-tagged catalog scheme is
in flight (task #79); it's an opt-in slow pass for thorough audits.

## 5. Materialization drivers

`Materialize*` scripts under `src/test/java/eu/solven/matmul/research/`
walk specific lists of target shapes, invoke a construction+composition
chain, verify with `Verifier.passesRandomMatmulSpotCheck`, and write
JSON files to disk. They are explicitly invoked via
`mvn exec:java -Dexec.mainClass=...` — they do *not* run automatically
during build.

| Script | Target shapes | Strategy used |
|---|---|---|
| `MaterializeSolvenStrassenNew` | cubic n ∈ {14, 17, 19, 21, 23} | `Recombination.constructWithAllocation` with hardcoded allocs |
| `MaterializeSolvenStrassen21` | `⟨21,21,21⟩` only | as above, refined alloc |
| `MaterializeSolvenStrassen777` | cubic n ∈ {7, 8, 9, 10, 11, 12} | as above for the leaves used by larger composes |
| `MaterializeRosowskiAlgorithm1` | `⟨n,3,3⟩` for various n | direct `RosowskiAlgorithm1.build(n)` |
| `MaterializeViaPanPair14` | `⟨14,14,14⟩` only | `RecombinationWithPair.constructWithPairing` with 3-pair / 1-solo decision |
| `MaterializeViaPanPair` | cubic n ∈ {14, 20, 22} where Pan pair is profitable | as above, auto-skips break-even cases |
| `SchemeSweep --mode=closure` | cubic n ∈ [4, 32] | iterated `findBestStrategy` with dependency-ordered shape sweep |
| `MaterialiseGaps` | non-cubic shapes ≤ MAX_DIM | hardcoded Strassen alloc per shape; mostly historical |

## 6. Where each FMM-better gap sits

As of latest audit (cubic + rectangular, up to MAX_DIM=32):

| Shape | ours | Sedoglavic | Strategy that closes it |
|---|--:|--:|---|
| ⟨18,18,18⟩ | 3402 | 3200 | Unknown rank-128 factor + two rank-5 factors (see type-polynomial in BORDER_RANKS.md analysis); not reachable by Strassen[9,9]³ (break-even). Task #55 (Kronecker enum) might find it. |
| ⟨23,23,23⟩ | 6738 | 6678 | DIS09-style multi-base; might be reproducible if we add Pan TA constructive recipe (task #44). |
| ⟨21,21,21⟩ | 5258 | 5240 | as above |
| ⟨19,19,19⟩ | 4044 | 4030 | as above |
| ⟨17,17,17⟩ | 2940 | 2930 | as above |
| ⟨2,12,16⟩ | 298 | 296 | HK asymmetric ⟨12,2,16⟩ (task #48) |
| ⟨2,10,16⟩ | 249 | 248 | HK asymmetric ⟨10,2,16⟩ (task #48) |
| ⟨2,10,15⟩ | 234 | 233 | HK asymmetric ⟨10,2,15⟩ (task #48) |
| ⟨9,11,15⟩ | 958 | 956 | unclear; possibly multi-base recombination with non-default pool |

## 7. What's NOT yet a strategy

Cataloguing the *absences* helps decide what to build next.

- **Pan pair-product fused inside the recombination search** (task #50): currently the search uses `findBestStrategy`, which picks the best of multi-base + concat. Pan pair isn't in the candidate list; `MaterializeViaPanPair` is a manual invocation.
- **Output-side zero peeling (Islam Ch. 4 / γ5 reduction)** (task #86): `Recombination.processAdditions` looks at INPUT-side zeros (sub-product whose A reads only the first `n_A` rows shrinks to ⟨n_A, n, n⟩) but doesn't propagate OUTPUT-side zeros from padding/peel. Explained in [PEELING_ZEROS.md](PEELING_ZEROS.md); closes ⟨17,17,17⟩=2934 and probably the whole odd-cubic gap family.
- **Systematic Kronecker product enumeration** (task #55): `Compose.kroneckerGeneral` exists, but no driver iterates `(scheme_i, scheme_j)` pairs against shape factorisations.
- **Pan trilinear-aggregation constructor** (task #44 cousin): we have `PanTrilinearAggregation.cubicBound` as a closed-form rank, but no constructive recipe. DIS09's compositions use Pan TA at Layer 3 and the published 250 / 143 / etc. catalog rows are reachable only via this path.
- **Asymmetric HK** (task #48): Lemma 1 done; band-restricted square + back-sub remain.
- **Lemma 3 sequences for n > 2k+2** (task #47 partial): handled implicitly via the asymmetric augmentation; the Lemma 3 boundary case never hits for pure square n=2k+2.
- **AlphaTensor-style RL search** (intentionally not pursued): outside the scope of "small-matrix SOTA tracking"; we ingest AT results without re-discovering them.

## 8. Recommended next step for closing each gap

1. **3 narrow ⟨2,a,b⟩ gaps**: finish task #48 (asymmetric HK band-restricted). Bedrock done.
2. **5 cubic large-n gaps**: task #55 (Kronecker enumeration) most likely to find improvements; then task #44 (Pan TA constructor) for the remaining.
3. **⟨18,18,18⟩=3200** specifically: type-polynomial analysis suggests a 3-factor Kronecker product; enumeration should find it if the factor schemes exist in our catalog.

---

If a strategy isn't in this index, it isn't tried. Keep this file
up-to-date as you add new ones.
