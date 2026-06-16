# Search strategy — finding new matmul rank upper bounds

How this repo discovered the new bounds tracked in
[NEW_BOUNDS.md](../../research/NEW_BOUNDS.md). Written so future runs (with a
richer catalog, more bases, or a tighter pipeline) can reproduce the
methodology and add new findings without re-deriving the framework.

---

## Pipeline overview

For each target `⟨n,m,p⟩` we want to find a (rank, scheme) pair where
the rank improves on the literature baseline. The pipeline is:

```
target ⟨n,m,p⟩
       │
       ▼
┌──────────────────────────────┐
│ Layer 1: pick OUTER BASE     │  Strassen ⟨2,2,2⟩=7,
│                              │  Laderman ⟨3,3,3⟩=23,
│                              │  Hopcroft-Kerr ⟨2,3,3⟩=15,
│                              │  mul211/mul121/mul112 (axis-splits)
└──────────────┬───────────────┘
               ▼
┌──────────────────────────────┐
│ Layer 2: min-reduction       │  Per-rank-1 term, compute the
│   (DIS09 §2)                 │  smallest sub-shape that the
│                              │  outer base's column support needs
└──────────────┬───────────────┘
               ▼
┌──────────────────────────────┐
│ Layer 3: S₃ symmetry         │  Try all 6 axis orientations of each
│                              │  base; different orientations give
│                              │  different sub-shape multisets
└──────────────┬───────────────┘
               ▼
┌──────────────────────────────┐
│ Layer 4: paired sub-products │  Pan's `(abc + ab + bc + ca)` for
│   (DIS09 §3, Pan 1980)       │  cyclic-pair shapes — replaces
│                              │  `2·R(⟨a,b,c⟩)` when smaller
└──────────────┬───────────────┘
               ▼
┌──────────────────────────────┐
│ Layer 5: Pan TA closed form  │  Direct bound for cubic n:
│   (DIS09 §3 appendix)        │   - n even: (n³+12n²+11n)/3
│                              │   - n odd : (n³+15n²+14n−6)/3
└──────────────┬───────────────┘
               ▼
┌──────────────────────────────┐
│ Layer 6: enumeration         │  - balanced (DIS09 default)
│                              │  - non-balanced (full search)
│                              │  Cap per-base combo count → auto-
│                              │  fallback to balanced for Laderman
│                              │  on large n (avoid combinatorial
│                              │  blow-up).
└──────────────┬───────────────┘
               ▼
        rank candidate
```

The SOTA resolver (`FieldAwareLookup` for explicit schemes; formulas
for Pan TA / Waksman / Rosowski) provides inner sub-product ranks to
the recombiner. The output is a **derived bound** — an upper bound on
`R(⟨n,m,p⟩)` plus the recipe to construct the algorithm.

---

## Where each layer lives in the code

| Layer | Implementation |
|---|---|
| 1 — outer-base pool | `BlockSplitSearch.NamedBase` + the buildPool() helpers in test code |
| 2 — min-reduction | `Recombination.processAdditions` (returns sub-shape from base column support × axis allocation) |
| 3 — S₃ orbit | `SymmetryTransforms.s3Orbit(base)` returns 1–6 oriented bases |
| 4 — paired sub-products | `PairedSubProducts.applyPairing(shapes, sota)` — greedy max-savings matching |
| 5 — Pan TA closed form | `PanTrilinearAggregation.cubicBound(n)` — injected into `BlockSplitSearch.loadCatalogBestRanksForField` so recursive lookups see it |
| 6 — enumeration | `BlockSplitSearch.findBestMultiBaseSplit(n, m, p, pool, sota, balancedOnly, maxCombosPerBase)` |

Driver: `TestDIS09FullScan` (NC, R/C/F2 fields), `TestDIS09Table4Commutative`.

---

## SOTA resolver — what feeds the recombiner's lookups

The resolver returns a rank for any `⟨a,b,c⟩` sub-shape. It's the min over:

1. **Catalog scheme files** — `src/main/resources/schemes/section*/`
   filtered per field-class with the chain `R ⊆ R/Q/Z`, `C ⊇ R`,
   `F₂` isolated. Implementation:
   `BlockSplitSearch.loadCatalogBestRanksForField(field)`.
2. **Pan TA formula** for cubic `⟨n,n,n⟩` over any non-commutative
   ring (per DIS09 appendix).
3. **Waksman formula** for commutative `⟨a,b,c⟩` (only in the
   commutative scan, via `WaksmanBound`).
4. **Rosowski formulas** (`RosowskiBound.bestCommutativeBound`,
   `RosowskiBound.nonBilinearRankBound`).

**Field discipline is load-bearing** — every catalog file is
classified by `FieldAwareLookup.classifyFilenameField` (F₂ /
C / Z-class) and only files in the requested field's fallback chain
are used. Skipping this leads to false wins (e.g. claiming
`⟨16,16,16⟩=2209` over R when it's actually F₂-only).

---

## Materialisation — from derived bound to verified scheme JSON

A derived bound is just a number plus a recipe. To produce an actual
verified scheme:

```java
NonCubicBilinearAlgorithm scheme = Recombination.constructWithAllocation(
        outerBase, FieldAwareLookup, allocA, allocB, allocC);
boolean ok = Verifier.passesRandomMatmulSpotCheck(scheme);
SchemeIO.write(scheme, outFile, /* z2 */ false);
```

Verification is via the **randomised matmul spot-check** —
`O(samples · r · (nm + mp + np))` instead of the tensor-residual
verifier's `O(n²m²p²·r)`. Runs 5 random `(A, B)` pairs through the
algorithm and compares to naive `A·B`. ~1000× faster, equally
definitive.

For symbolic debugging (when a scheme fails the spot-check): use
`Verifier.symbolicDiscrepancies(alg, maxReport, tolerance)` which
reports specific bilinear-term coefficient errors — pinpoints which
inner sub-scheme or outer placement is wrong.

---

## Known limitations / wins on the table

### Layer 4 isn't constructive yet
`PairedSubProducts` computes a savings-aware rank ESTIMATE but doesn't
build the explicit `(α+β)·(γ+δ) − (α−β)·(γ−δ)` factor matrices.
Materialisation thus uses the unpaired construction, which is
slightly worse than the paired bound. Fix: write the constructive
version of `applyPairing` that emits factor matrices for paired
products.

### Catalog gaps
Several inner sub-shapes are looked up via formulas only — there's no
explicit scheme on disk. When constructing larger composites, the
construction can't proceed past a shape with no concrete scheme.
Mitigation: recursive materialiser (ROADMAP item) that builds
intermediate schemes on demand by recursing into smaller sub-shapes.

### Broken-import denylist
27 AT-Z composed-improvement schemes (dims 9–12) fail standalone
spot-check — our Python import has a convention bug. Listed in
`FieldAwareLookup.KNOWN_BROKEN_FILES` and skipped by lookup +
`GenerateCatalogManifest`. **Workaround**: fall back to verified
alternatives. **Fix**: re-fetch from FMM or Perminov status.json
(ROADMAP item — both catalogs re-host the same AT results with proper
attribution).

---

## How to re-run the search

```bash
# (1) Refresh catalog from on-disk schemes (regenerates docs/catalog.json)
mvn -q -o exec:java -Dexec.classpathScope=test \
    -Dexec.mainClass=io.cormoran.strassen.v3.catalog.GenerateCatalogManifest

# (2) Regenerate derived bounds (multi-base + symmetry + Pan TA + paired)
mvn -q -o exec:java -Dexec.classpathScope=test \
    -Dexec.mainClass=io.cormoran.strassen.v3.catalog.GenerateDerivedBounds

# (3) Compare against DIS09 — see wins/ties/losses per field
mvn -q -o test -Dtest=TestDIS09FullScan        # NC, R/C/F2 fields
mvn -q -o test -Dtest=TestDIS09Table4Commutative  # commutative

# (4) Materialise a specific new bound (n=19, 21, 23 as concrete examples)
mvn -q -o exec:java -Dexec.classpathScope=test \
    -Dexec.mainClass=io.cormoran.strassen.v3.catalog.MaterializeSolvenStrassenNew

# (5) Verify every published scheme (catches data-quality regressions)
mvn -q -o test -Dtest=TestPublishedCatalogVerifies
```

---

## How to extend the search

Most impactful improvements, ranked by expected payoff:

1. **Bigger base pool** — add ⟨3,4,3⟩=29, ⟨4,4,3⟩=38, AT-Z direct
   discoveries (⟨3,4,5⟩=47, ⟨4,4,5⟩=63, ⟨4,5,5⟩=76, ⟨5,5,5⟩=96/F₂,
   ⟨4,4,4⟩=47/F₂) as outer bases. Each enriches the per-shape coverage.
2. **Fix the AT-Z import** — re-fetch from FMM/Perminov, drop the
   denylist. Restores ~27 sub-shape options that today fall back to
   worse alternatives.
3. **Implement Layer 4 constructively** — `PairedSubProducts` currently
   only computes a rank estimate. The constructive version emits the
   explicit paired factor matrices, unlocking the 1720 / 4030 / 5258
   / 6731 ranks (vs the current materialised 1750 / 4053 / 5276 / 6738).
4. **Recursive materialiser** — top-down: when an inner sub-shape
   lacks an explicit scheme, recursively construct it. Unblocks
   materialisation for any derived bound.
5. **Larger n via Pan TA recursion** — Pan TA at large n is a closed
   form; combine with Strassen/Laderman recursion for `n > 30`.
6. **Look at non-cubic targets** — current scan focuses on cubic
   `⟨n,n,n⟩` (DIS09's table format). Many non-cubic targets are also
   tracked by the catalog and could see similar wins from the same
   pipeline.
7. **Commutative-specific bases** — Waksman family at small n, plus
   Rosowski Algorithm 1 / Corollary 1 (now in catalog as non-bilinear
   schemes) as outer bases for the commutative scan.

Each item is a focused PR. Most don't require new theory — just
plumbing more inputs through the existing layers.

---

## Methodology for adding a new finding

When you see a new bound from a paper / catalog / a fresh run:

1. **Reproduce** the rank claim — re-run the relevant test
   (`TestDIS09FullScan`, `TestDIS09Table4Commutative`, or a custom
   probe). Confirm the SOTA resolver returns the rank.
2. **Materialise** — construct the explicit scheme via
   `Recombination.constructWithAllocation` and verify it.
3. **Add to catalog** — write the scheme JSON to
   `src/main/resources/schemes/section{maxDim}/`, regenerate the
   manifest, ensure `TestPublishedCatalogVerifies` still passes.
4. **Cite** — add an entry to `NEW_BOUNDS.md` (verified column),
   update `DISCOVERIES_PENDING_ANALYSIS.md` if there's anything still
   provisional, and link to the originating commit/PR.
5. **Cross-check field** — make sure the rank actually applies in the
   claimed field (the `R` filter has been a source of subtle bugs;
   see ROADMAP "Re-fetch broken AT-Z composed schemes").

If steps 1 + 2 disagree (derived bound says X, materialised gives Y),
the gap is information — usually an inner sub-scheme that the SOTA
resolver knew about but isn't physically on disk. Treat as a
materialiser improvement (item 4 above).
