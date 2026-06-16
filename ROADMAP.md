# Roadmap

Open work items, organized by area. **Status legend**: 🟢 ready / quick-win
| 🟡 planned, scope clear | 🔴 research-grade, multi-day | ⚪ exploratory.

When picking up an item, move its status to `🚧 in-progress @<initials>` and
strike it through when done.

---

## Catalog & schemes

### On-demand scheme fetcher (FMM + Perminov) 🟡
**Status:** identified 2026-05-29. A CLI / Java tool that takes a
target `⟨n,m,p⟩` and returns the best (rank, scheme) by:

1. Querying [Perminov status.json](https://github.com/dronperminov/FastMatrixMultiplication/blob/master/schemes/status.json)
   for `ranks` at the sorted format. Pick the lowest rank across the
   Q/Z/ZT fields valid for the requested field-chain (see "Practical
   fallback rule" in RANK_KNOWLEDGE.md §1.2bis).
2. Querying [FMM Université de Lille](https://fmm.univ-lille.fr/{n}x{m}x{p}.html)
   for the page-listed rank — scraping or by precomputed index.
3. Pick the lower of the two; download the actual factor matrices
   (Perminov: fetch the `path` JSON from GitHub; FMM: download
   `{n}x{m}x{p}_raw.mpl.bz2`, bunzip, parse Maple). Save to our
   `src/main/resources/schemes/section{maxDim}/` with proper
   `attribution_for_rank`.

Today the "recursive materialiser" (next item) fails as soon as it
hits an inner sub-shape not in our local catalog — this fetcher
fills those gaps on demand.

### GitHub Action — periodic Perminov status.json sync 🟡
**Status:** identified 2026-05-28. The Perminov GitHub repo now
publishes
[`schemes/status.json`](https://github.com/dronperminov/FastMatrixMultiplication/blob/master/schemes/status.json)
— a synthetic index listing every scheme with `(rank, complexity)`,
the `source` sub-folder (which encodes TRUE provenance:
`classic/`, `alpha_tensor/`, `alpha_evolve/`, `a_60_addition/`,
`jakobmoosbauer_symmetric_flips/`, `meta_flip_graph/`, `tensor/`,
`fmm_add_reduction/`), and the JSON path in their repo.

Plan: a GH Action that runs weekly to (a) fetch status.json, (b)
diff against our `src/main/resources/schemes/`, (c) import any new
schemes preserving the `source` attribution (NOT tagging everything as
"Perminov" as the current importer does), (d) open a PR with the
delta + a summary. Closes the gap that drives the
`TestCatalogYearMetadata` worklist (~2,200 Perminov entries
currently lack proper attribution).

### Re-fetch broken AT-Z composed schemes from FMM / Perminov 🟡
**Status:** identified 2026-05-29. The AT-Z composed-improvement
schemes for dimensions ≥9 (27 files, all flagged in
`FieldAwareLookup.KNOWN_BROKEN_FILES`) don't spot-check after import.
Root cause is the Python import script (`tools/import_alphatensor.py`
or similar) — AT 2022's published factor matrices are correct, so our
JSON encoding has a convention bug specific to those higher-dimension
files.

Simpler path than debugging the Python: re-fetch the same schemes via
the on-demand fetcher (FMM `/{n}x{m}x{p}.html` or Perminov
`status.json`). Both catalogs re-host AT's results with correct
attribution. Once re-imported, drop the denylist entries.

Affected files: 27 listed in `FieldAwareLookup.KNOWN_BROKEN_FILES` —
mostly AT-Z files with dimensions 9, 10, 11, 12 in their format.

### AT bulk-import provenance audit ✅
Done 2026-05-28. 113 `alphatensor-*.json` files tagged with
`discovery: true|false` + `attribution_for_rank` reflecting the true
originator (per AT Extended Data Table 1). 77 confirmed AT
discoveries; 36 rediscoveries attributed to Hopcroft-Kerr 1971,
Smirnov 2013, Sedoglavic-Smirnov 2021, Strassen 1969 or Laderman 1976.

### Materialise DIS09-equality entries as explicit schemes 🟡
DIS09 published bounds in Tables 3 & 4 but **no factor matrices, and
the code generator described in §4 of the paper appears lost** — the
URL `csd.uwo.ca/~mislam63/` referenced in the paper is dead. For each
cubic ⟨n,n,n⟩ where this repo's pipeline matches or beats DIS09 (n=4..30),
materialise the explicit scheme via
`Recombination.constructWithAllocation` and save as a JSON file, so the
catalog has a verified scheme rather than just a bound claim. Affects
all 27 NC entries (15 wins + 12 ties) and 25 commutative entries.

### Synthesise Waksman 1970 commutative schemes 🟡
Waksman's algorithm is a closed-form family (rank `b(ac+a+c-1)/2` for
b even) but published as a construction, not as factor matrices.
[fmm-lille](https://fmm.univ-lille.fr/) hosts only non-commutative
schemes. Mezzarobba 2007 tabulated ranks but didn't publish matrices.
To make Waksman ⟨n,n,n⟩ a first-class catalog entry, write a
`WaksmanGenerator.java` that emits the explicit (quadratic-algorithm)
U/V/W matrices from the formula. This unlocks recursive use of Waksman
as an outer base in commutative compositions — currently we only have
the BOUND, not the scheme. Requires extending SchemeIO for quadratic
algorithms (U_l and V_l may each contain both A and B entries).

### Stabilizer annotations 🟡
**Status:** detector implemented (`StabilizerDetector.detectZ3`),
verified on Strassen (Z/3) and Laderman (trivial). NOT yet:
- run across the whole catalog (2,353 schemes) to populate the
  `"stabilizer"` field in scheme JSON
- surfaced in `COVERAGE.md` notes column
- surfaced in `docs/index.html` browser as a filter + column
- extended beyond Z/3: Z/2 (transpose-like), S_3 (full slot symmetric
  group), non-cubic-specific stabilizers (e.g. swap A ↔ C^T when n = p)
- proper sign-flip handling (current detector uses a canonical sign
  flip; some algorithms need per-column independent signs)

Plan: add a one-shot `ComputeStabilizers` main that loads each cubic
scheme, runs the detector, and edits the JSON in place to add
`"stabilizer": "Z/3"` / `"trivial"`. Then extend
`GenerateCatalogManifest` to include the field and `catalog.js` to
filter on it. See `SYMMETRIES.md` for the underlying group structure.

### Lower bounds catalog 🟢→🟡
**Status:** seeded `docs/lower-bounds.json` with 11 well-cited entries
covering ⟨2,2,2⟩ through ⟨5,5,5⟩. Pages view (`docs/catalog.js`)
cross-references the LB and shows a "gap" column per scheme row
(rank − LB; `0 (tight)` when an algorithm matches the LB).

NOT yet:
- broader coverage: ⟨2,n,m⟩ for `n,m ≤ 8` (Hopcroft-Kerr formulas often
  apply), ⟨3,3,n⟩, ⟨3,n,n⟩, ⟨4,4,n⟩
- URL for Wang 2026 (placeholder note; need confirmed arXiv link)
- per-field LBs vs cross-field (e.g. distinguish `R_R`, `R_C`, `R_{F₂}`,
  border-rank `R̲`)
- direct integration with `COVERAGE.md` and `BlockSplitSearch` so that
  a Phase-1 formula prediction is automatically tagged as "still above
  LB" vs "would match LB"

### fmm-lille biblio integration 🟢→🟡
**Status:** scraped + cross-referenced.
- `eu.solven.matmul.docs.verify.SyncReferenceCatalogs --fmm` →
  `references/fmm-lille-biblio.json` and
  `references/fmm-lille-catalog.json` (5,426 rows).
- `tools/compare_fmm_lille.py` → `references/fmm-lille-discrepancies.md`:
  **171 formats where our schemes are sub-optimal vs fmm-lille**, 463
  matched, 6 where we beat fmm-lille (F₂-specific from AlphaTensor), and
  **4,786 formats fmm-lille knows but we don't have a scheme for**.
- `docs/catalog.json` now embeds an `fmm_lille` cross-ref per scheme
  (best_rank + details_url + references).

NOT yet:
- per-format detail-page scrape (`--details` flag, ~45 min) to extract
  the full algorithm-variant list per format → richer reference linking
- per-entry annotation "what does this paper contribute" for the ~80–100
  novel-algorithm papers
- update `REFERENCES.md` with curated 30–40 highest-impact entries from
  the new biblio data
- catalog UI: sub-bibliography view (filter algorithms → list only their
  referenced papers)
- triage the 171 sub-optimal cases: are they fixable by importing a
  specific fmm-lille scheme, by composition, or by SAT search?
- pull stabilizer column from fmm-lille's catalog presentation page if
  we don't compute it ourselves (see "Stabilizer annotations" above)

### Basis-change support in SchemeIO 🟡
**Status:** none — current `SchemeIO` only handles bilinear schemes
without preprocessing/postprocessing transforms.

**Motivation:** algorithms like Karstadt-Schwartz 2017 reduce
⟨2,2,2⟩=7 to **12 additions** (vs Strassen's 15-18) by applying
invertible linear transforms `A' = P·A·Q⁻¹`, `B' = Q·B·R⁻¹` before
the algorithm and `C = P⁻¹·C'·R` after. The transforms amortize for
free over recursive matmul on large inputs but the inner addition
count is the headline claim.

**Required changes:**
- Extend scheme JSON: add optional `basis_change: {P, Q, R}` block.
- `SchemeIO.read/write` round-trips P, Q, R.
- `Verifier.isExactNonCubic` applies the transforms when present.
- Import Karstadt-Schwartz and similar basis-changed schemes.

Without this, claims like "12-add Strassen" can't be represented; the
underlying algorithm would appear with ~20+ adds because the
basis-change work bakes into the per-product linear combinations.

### Reproducible split: `known/` schemes vs `ours/` schemes 🟡
**Status:** partial — sources are distinguishable by filename prefix
(`alphatensor_*`, `perminov_*`, `fmm-lille_*`, `kauers_2026-*`, … =
imported; `derived_recursive-*`, `derived_kron_concat-*`, `solven-*` =
ours) but they live intermixed under `src/main/resources/schemes/sectionN/`.

**Goal:** mirror the same directory structure under two roots so the
provenance is structural, not just a filename convention:
- `schemes/known/sectionN/…` — schemes imported from the literature /
  community catalogs (the *inputs* we don't derive). Treated as
  read-only ground truth.
- `schemes/ours/sectionN/…` — everything our search/closure produces
  (recombination, Kronecker, serendipitous, projection wins, stubs).

**Why:** **droppable + rebuildable.** We should be able to delete the
entire `ours/` tree and regenerate it deterministically from `known/`
alone, by re-running the closure/sweep — and get back the same catalog
(content hashes make this checkable). This is the reproducibility
guarantee: the historical record (`known/`) is the seed; our
contribution (`ours/`) is a pure function of it plus the search rules.
It also makes "what did *we* actually find vs re-import" answerable by
`ls`, and keeps the [[feedback_discovery_is_multidimensional]] /
discovery-vs-rediscovery accounting honest.

**Required changes:**
- `FieldAwareLookup` / `SchemeSweep` / manifest generator read from BOTH
  roots (known ∪ ours), with `ours/` as the write target (staging
  promotes into `ours/`, never `known/`).
- A `--rebuild-ours` driver: wipe `ours/`, run the full closure from
  `known/`, diff the regenerated tree against the committed `ours/` by
  content hash to certify reproducibility.
- Migration: classify every current `sectionN/*.json` by source into
  `known/` vs `ours/` and move it (filename prefix + lineage atom-ness
  is the signal; `atom:true` ⇒ known, composed ⇒ ours).
- Keep the manifest (`docs/catalog.json`) merging both, tagging each row
  with its root so the SPA can filter "ours only".

Related: this is the structural form of the
[[project_research_grade_catalog]] direction and the #199 3-catalog
comparison; the content-hash work already gives the rebuild-equality
check its primitive.

### Catalog expansion 🟡
Schemes we know exist in the literature but haven't imported:
- **Smirnov's `{-1,0,+1}` extended catalog** — beyond what Perminov
  ships (covers many small non-cubic at lower ranks)
- **Kauers & Wood 2025** — recent improvements per fmm-lille (`⟨2,5,7⟩`,
  `⟨2,5,8⟩`, `⟨2,6,7⟩`, `⟨2,6,8⟩`, `⟨2,7,7⟩` …)
- **Holtz et al. 2025** — "Alternative bases for matrix-multiplication algorithms"
- **Perminov 2026** — newer entries (`⟨2,4,11⟩`, `⟨3,5,9⟩`, `⟨3,5,10⟩`)
  in fmm-lille's "most recent" pane
- 4,786 "missing" formats per the discrepancy report — many will be
  filled in once we import the relevant Smirnov / Heun / Sedoglavic /
  Kauers-Wood schemes.

### Schwartz-Zwecher 2025 — n=34..50 import deferred 🟡
The supplemental zip
(`https://www.cs.huji.ac.il/~odedsc/papers/trilinear_aggregation_algorithms_decomposed-2025-07-29.zip`)
ships explicit factor matrices for `Q⟨n,n,n⟩:r` for n ∈
{20, 22, …, 50}. We've imported n=20..32 (within the current
`CatalogLimits.MAX_DIM=32` envelope, file sizes 6–35 MB each).
**Imports for n=34..50 are deferred** because:
- our derivation procedure (`PanTrilinearAggregation` / kin-row
  unification) does not yet construct these schemes from first
  principles, so importing the explicit matrices would commit us
  to carrying very large (~100s of MB at n=44, the headline
  result with rank 36110) opaque blobs in the repo;
- our random-spot-check verifier handles them, but the rest of
  the catalog tooling (sweep drivers, lineage cache, JSON
  formatter) hasn't been benchmarked at these sizes.

When we extend our own constructor to support n ≥ 34 (Section §3
of arXiv:2508.01748: kin-row unification + Pan's aggregation
tables), the n=34..50 SZ2025 schemes can either be derived
internally or imported as anchors.

We may **make an exception for n=44** (the asymptotic-optimal
base case, `ω₀ ≈ 2.773203`) once we have a clear way to verify
it at acceptable cost and the catalog tooling tolerates a
single very-large scheme. Until then, the .npz file lives under
`references/schwartz-zwecher-2025/trilinear_aggregation_algorithms_decomposed/`
and the bound is registered without the scheme via
`docs/cited-bounds.json`.

### Schwartz-Zwecher 2025 — kin-row unification constructor 🟡
Derivation follow-up: implement the construction described in
arXiv:2508.01748 §3 (TA-New25 family) and §4 (TA-New25b) so we
can SYNTHESISE these schemes from Pan's aggregation tables plus
a chosen `⟨2,2,2;7⟩`-algorithm (cf. Theorem 2.22). Inputs:
- Pan's aggregation tables (`PanTrilinearAggregation` already
  emits the bound; needs an emitter for the explicit `U`/`V`/`W`);
- Hadas-Schwartz 2023 implicit-canceling transformation;
- Theorem 2.22's family of `⟨2,2,2;7⟩` algorithms parameterised
  by invertible `K_U, K_V ∈ F²ˣ²`.

Once that lands, the imported SZ2025 schemes can be flagged as
re-derivable, and we can drop them in favour of computed ones —
per the project's "prefer derivable over imported" rule. Tracker:
the imported schemes' JSONs carry
`"derivation_task": "TBD-SZ2025-kin-row-constructor"`.

### HK 2bc constructive closure ✅ (2026-06-11/12 — one residue open)
**Status:** SOLVED. The constructive procedure exists and shipped:
**465 integer (ℤ) schemes** over `3 ≤ p ≤ 32, p ≤ n ≤ 32` in
`schemes/constructed/`, **459 at the exact HK formula**
`⌈(3np + max(n,p))/2⌉`, quadruple-gate + independently certified
(465/465 exact). The old gap table closed with strict catalog wins:
`⟨2,10,15⟩=233` (was 234 everywhere), `⟨2,10,16⟩=248`, `⟨2,12,16⟩=296`
(was 298), plus 197+22 further strict improvements. Key unlocks:
arc-interior bridge selection; `(2,2,b3)`/`(1,1,b3)` solved over the
TRUE 12-product reusable set (the published impossibility's S had 6);
repaired Case-2 Step 3 (circulant matching + Z-pairs); chained
augmentation for `n > 2p−1` (DP over achieved segment ranks); and
all-unimodular Lemma-1 matrices by Euclidean recursion → integer
back-substitution → ℤ schemes. Recipe + proofs:
[`research/hopcroft-kerr-2np/CONSTRUCTIVE_METHOD.md`](research/hopcroft-kerr-2np/CONSTRUCTIVE_METHOD.md);
paper section `paper/sections/hk71.tex`; historical impossibility doc
(status-bannered): [`references/hopcroftkerr1971/README.md`](references/hopcroftkerr1971/README.md).

**Open residue (task #9):** the six `g ≥ 6` circulant shapes sit at
+1..+3 over the formula (still below every catalog) — provably
formula-impossible in-framework: arc-sum placement + the
`(3,3,bridge-1/2)` impossibility THEOREM over arbitrary local rank-1
atoms (9-product reusable set, `sympy/derive_33bridge_general.py`).
Resume route: the 12-product reusable set (finite quadratic system
after structural substitution — Gröbner), non-local atoms, or
4-product trades. Lower bounds: parked by explicit decision.

### Deferred `_reduced` schemes 🟡
**Status:** 167 `_reduced` schemes now verify end-to-end after the
nested-`_fresh` fix. NOT yet:
- ~50 deeper-nested schemes that still hit `ArrayIndexOutOfBounds` —
  the W-axis fresh expansion handles one level of nesting but some files
  have multi-level CSE chains we don't yet expand.
- `_reduced` writers — we only READ this format; SchemeIO.write always
  emits dense or our new sparse, never `_reduced`.

---

## Search & discovery

### Projection operator (project a larger scheme down + DCE) 🟢 — IMPLEMENTED
**Status:** implemented 2026-06-05 (#159). `Compose.project` + `ProjectionSearch`
(cheap survivor-count, build/verify only the winner) + a `Lineage.Project`
node (rendered/parsed/replayed/field-inferred) + a `RecursiveMaterialiser.tryProjection`
hook + a **downward projection closure** (`SchemeSweep` Phase 3,
`--projection-only` for a standalone run). Parents resolve via the disk lookup
(≤16) or stub-replay (>16); every emitted `Project` node is replay-checked
before persist. Tests: `TestProjection` (operator, search, lineage round-trip) —
all green.

**Operator is done & fast on the ≤16 band.** A 2026-06-05 probe at ⟨25³⟩←⟨26³⟩
established that the **17–32 headline win is gated by two PRE-EXISTING issues, not
by the operator:**
1. **Best-rank parents are non-replayable stubs.** The best ⟨26³⟩ we hold
   (`dis09_Q-26x26x26_m8658`) has a lineage referencing the **`DIS09Lemma4(n=26)`**
   parametric constructor, which `LineageReplayer` can't reconstruct (it only knows
   the small named bases). So `findWithSource`/replay can't surface it → projection
   never sees the good 8658 cube. This is **task #68** (stub-replay verification).
   *(Already fixed one such gap: `naive-NxMxP` trivial-axis leaf refs now replay.)*
2. **Per-target replay cost.** Resolving a parent re-replays a large stub to dense
   `double[][]` (a single ⟨26³⟩ replay is multi-minute), and the closure does this
   per target × up to 7 δ≤1 parents — far too slow. *(The combo enumeration itself
   is now cheap: `ProjectionSearch` precomputes each product's index-support once
   per parent, so survivor-counting is O(rank·tiny) not O(rank·keepArea·combos);
   `resolveParent` now logs replay timing so a run never stalls silently.)*

**Next step — parent-centric amortized closure (the practical 17–32 unlock):**
invert the pass — for each *good cube/near-cube parent*, replay it **once** and
project it **down to all δ≤1 children**, instead of re-replaying all parents per
target. Combined with making `DIS09Lemma4`-style stubs replayable (#68), this is
what actually lands FMM's projection wins on 17–32.

**Original identification:** 2026-06-05 from a US-vs-FMM diff. FMM beats us on **2530
shapes, entirely in the maxDim 17–32 band**, and reading its per-format pages
shows the construction is **overwhelmingly `projection`** (Perminov draft
Def 2.8 / meta-flip-graph `Project`), almost always projecting the **next
near-cube** down:
- `⟨21,22,22⟩=5476 = projection [[1,0],[0]] of ⟨22,22,22⟩=5566`
- `⟨27,28,28⟩=10413 = projection of ⟨28,28,28⟩=10556`
- `⟨29,30,30⟩=12588 = projection of ⟨30,30,30⟩=12710`
- `⟨25,25,25⟩=8359 = projection [[1,15],[15]] of ⟨26,26,26⟩=8658`
- `⟨9,17,17⟩=1600 = projection of ⟨9,17,18⟩` (rank unchanged)

**Mechanism:** to get `⟨n,m,p⟩` from `⟨n+1,m,p⟩` (etc.), zero one input slice
and don't compute the corresponding output slice; **DCE** the products that
fed *only* the dropped index. Rank drops by however many products were
localized there. The bracket notation `[[i,j],[k]]` selects *which* index to
drop per axis — and FMM picks the index that maximises the DCE, so the operator
must try **each** index, not just the last.

**Why it's the unlock:** we're at parity with FMM on the **even cubes**
(⟨26³⟩, ⟨28³⟩, ⟨30³⟩, …) — they're absent from the FMM-better list — so FMM's
edge is purely that it *projects* those good cubes onto every nearby non-cube
and odd cube, which we never do. Wiring this closes most of the 17–32 gap.

**Plan:**
- `Compose.project(scheme, dropA[], dropB[], dropC[])` → drop the chosen index
  rows/cols from U/V/W, DCE products feeding no kept output, verify.
- Search: for target `⟨n,m,p⟩`, try projecting catalog schemes of
  `⟨n+a,m+b,p+c⟩` (small a,b,c ≥ 0, prioritise the near-cube `⟨M,M,M⟩`,
  M=max), choosing the best drop-index per axis; keep min verified < SOTA.
- Lineage: a `Project(child, dropA, dropB, dropC)` node (replayable).
- Honesty: a projected scheme is exact; rank is a **bound** (best found),
  proven only if the parent + drop choice are jointly optimal (they're not).

**It is NOT a new primitive — it is our peel/DCE generalised.** Projecting
`⟨n+1,m,p⟩→⟨n,m,p⟩` = zero one input slice + skip the matching output + DCE the
products localised to it = padding+DCE run in reverse. We already have the DCE
machinery: output-zero masks (#86, "Islam γ5") + recombination tail-peel. Two
generalisations turn it into FMM's projection:
1. **source = any larger *catalog* scheme** (e.g. `⟨26³⟩=8658`), not just a
   recursion-padded small base pattern;
2. **drop the *best* index per axis**, not just the tail — exactly the pending
   **#144** (peel positions beyond tail-only); FMM's `[[1,15],[15]]` = "drop
   index 15," chosen to maximise DCE.

**Concrete unlock now:** we already hold every even cube at *exactly* FMM's rank
(⟨24³⟩=7000, ⟨26³⟩=8658, ⟨28³⟩=10550, ⟨30³⟩=12688, ⟨32³⟩=15096), and FMM builds
the odd cubes + near-cube non-cubes by projecting these down. So the sources are
in hand; projecting them is the bulk of closing the 17–32 gap.

**Multipass implication:** projection propagates cube quality *downward*
(better `⟨26³⟩` ⇒ better `⟨25³⟩`, `⟨25,26,26⟩`, …). Our closure only composes
*upward* today; add a downward projection pass and iterate to fixpoint — improve
cubes → project to neighbours → re-evaluate → repeat.

Complements serendipitous (multiplicative bud cases) and disjoint-sum/τ
(additive). Together these are FMM's construction taxonomy: author-year leaves
+ `⊗` + `+` + **projection** (= generalised peel/DCE).

### τ-theorem disjoint-sum identity synthesis (the additive FMM gap) 🟡 — RE-SCOPED
**Status:** scoped 2026-06-06, then **largely reframed the same day** (see below).
`references/TAU_DISJOINT_SUM_SEARCH.md`. The FMM-vs-us comparison (Q, NC) shows
**70 formats at maxDim≤17 where FMM wins**, almost all `⟨a,b,17⟩`, presented by FMM
as a **"Trace·Mul" N-term disjoint sum**. **Canonical case — ⟨6,17,17⟩=1106:**
`8·⟨3,6,6⟩=80 + 3·⟨3,6,5⟩=68 + 3·⟨3,5,6⟩=68 + 1·⟨3,5,5⟩=58` (15 blocks).

**⚠ REFRAMING (2026-06-06):** ⟨6,17,17⟩=1106 turned out to be **plain recombination
of Hopcroft–Kerr ⟨2,3,3⟩=15** at allocation n=(3,3), m=(6,6,5), p=(6,6,5) — the
15 blocks ARE HK's 15 products (full grid is 2·3·3=18; HK uses 15). Confirmed:
`SchemeSweep --shape=6x17x17 --base=2x3x3` → **1106** (beats our prior 1154).
So FMM's "disjoint sum / trace identity" is just **how they write a recombination
from a small outer base** — NOT a mechanism we lack. Our default closure missed it
only because ⟨2,3,3⟩=15 wasn't tried as an OUTER base at the (6,6,5) allocation.
Same lesson as ⟨17³⟩=2930: **the gap is pool-completeness, not τ-synthesis.**

**New primary plan (cheap): enrich the outer recombination pool** with the small
efficient bases (⟨2,3,3⟩=15, ⟨2,2,3⟩=11, ⟨3,3,3⟩=23, Smirnov bases, …) in all
axis orientations + keep unbalanced allocations, then re-sweep the 70 formats — most
should fall out as recombinations. The bespoke τ-certifier (below) is demoted to
the residue of genuinely-non-recombination cases (not yet observed). Use
`--base=<NxMxP>` to test a candidate base on any shape first.

**(Demoted) constructive τ-certifier**, for any true residue:
**We hold every block at FMM's exact rank** — so this was never block coverage.

What exists: `KnownTauIdentities` (hand-extracted, static — lacks this one) and
`DisjointSumSearch` (enumerates multisets but enforces only the **loose area-cover**
bounds → over-optimistic, non-constructive). **Missing = constructive
certification:** given a candidate multiset + axis partition, prove a trilinear
embedding realises the target at the summed rank, and emit verified U/V/W
(`Lineage.DisjointSum`). Start with the **grid-minus-k** family (covers
⟨6,17,17⟩ + most `⟨a,b,17⟩`); M0 acceptance test = reproduce ⟨6,17,17⟩=1106
verified. Honesty: certified results are `bound`; area-cover-only predictions are
heuristics, **never** emitted as ranks. Ties #196, #102.

### Incremental sweep on a single new scheme 🟡
**Status:** requested 2026-06-05. Once the project stabilises we won't want to
re-run full closures all the time. When **one new scheme `S = ⟨n,m,p⟩=r`** lands
(an import, an ALS/SAT find, a manual construction), only a bounded neighbourhood
can change — so sweep just that neighbourhood instead of the whole catalog:

- **Upward (S as a building block):** re-sweep the band of shapes **≥ `⟨n,m,p⟩`**
  (the shapes that could use `S` as a Kron/concat/recombination/serendipitous
  factor), *excluding projection* — projection never consumes a smaller scheme.
- **Downward (S as a projection parent):** re-sweep the band of shapes
  **≤ `⟨n,m,p⟩`** reachable by projecting `S` down (within `PROJECTION_MAX_DELTA`),
  considering **only `S` itself** as the projection parent — no need to re-project
  the rest of the catalog, which hasn't changed.

So the work is: `{shapes that can be built from S, upward, no projection}` ∪
`{shapes S can project to, downward, S-only}`. Everything else is provably
unaffected by adding `S`. Wire as a `SchemeSweep --incremental=NxMxP` (or
`--since-file=...`) mode that derives these two bands from the new scheme's shape
and runs the existing search + projection phases restricted to them. Big
wall-clock win vs a full sweep; keeps the catalog fresh on every single addition.

### Multi-scheme store per format — retain bud-diverse variants 🔴
**Status:** identified 2026-06-05 (serendipitous-product work, #159). Today
the catalog/closure keeps **one best-rank scheme per format**, which throws
away the bud-rich alternatives the **serendipitous product** needs. A scheme's
*bud structure* (terms sharing a `u`/`v`/`w` vector → elementary `⟨1,1,k⟩` /
`⟨k,1,1⟩` / `⟨1,k,1⟩` blocks; see `references/SERENDIPITOUS_PARTIAL_PRODUCT.md`
and `paper/sections/strategies.tex` §serendipitous) determines how much a
serendipitous product `T₁⊗⟨n₂,m₂,p₂⟩` saves over the naive Kronecker.

The catch: the **rank-optimal** base often has *fewer* buds than a slightly
higher-rank one, so a sub-optimal scheme can yield a *better* serendipitous
product. Clean condition: a base at rank `r₁+1` with one extra `U`-bud beats
the rank-`r₁` base iff the per-bud saving `s > r₂/2` (where `r₂ = R(⟨n₂,m₂,p₂⟩)`).
Perminov confirms (draft p.9): `⟨3,4,13⟩` has optimal rank 117, yet a **rank-123**
scheme was used because its bud structure gives a better product.

**Proposal:** keep *multiple* schemes per format — at minimum the min-rank one
plus bud-rich higher-rank variants, tagged by their (`U`,`V`,`W`)-bud
multiset — and let `SerendipitousSearch` try all of them, keeping the best
*final* product (not the best base). This changes the closure's core invariant
(per-shape state becomes a *scheme set*, not a single best rank). Phasing:
- **Phase 1** (done/quick): serendipitous as an upgraded-Kronecker candidate on
  the *current best* base per shape — catches divisible+bud-rich-optimal cases
  (e.g. AT `⟨2,3,4⟩=20 → ⟨4,8,12⟩=272`). `SerendipitousBudProduct` /
  `SerendipitousSearch` exist; needs a `Lineage.SerendipitousProduct` node +
  `findBestStrategy` hook.
- **Phase 2** (this item): the multi-scheme store → captures the `101-3 > 100-1`
  wins.
- **Phase 3** (research): a *bud-rich scheme generator* (flip-graph-style, à la
  Perminov) to actively produce bud-diverse schemes rather than only exploit
  those already present. Without it our reach is bounded by catalog/composed
  schemes' existing buds.

Related: the cheap, intrinsic `buds` per-scheme field (bud multiset + `has_buds`
flag) proposed alongside #159 would make this store's tagging free.

**Bud-alignment under isotropy (search enhancement, user 2026-06-06).** The
Kronecker product *multiplies* bud proportionalities: a U-bud of size `a` in
T₁ aligned with a U-bud of size `b` in T₂ (same axis) gives a U-bud of size
`a·b`; crossed (U vs V) it gives only `max(a,b)`. The **only** realignment
lever is the discrete **S₃ axis permutation** — a `GL` basis change preserves
proportionality, so it leaves the bud partition invariant. Two actionable
consequences for `SerendipitousSearch`:
1. **Asymmetric product** `T₁⊗⟨n₂,m₂,p₂⟩` (only T₁'s buds used): **try all S₃
   orientations of the bud-provider** — a U-bud enlarges to `⟨n₂,m₂,k·p₂⟩`, a
   V-bud to `⟨k·n₂,m₂,p₂⟩`, whose `R` differ, so orientation changes `r_s`
   directly. (We generate the orbit already; just exploit it here.)
2. **Recursive/full `T₁⊗T₂`**: pick orientations that **align both factors'
   bud axes** to maximise the `a·b` reinforcement → a bud-richer base for the
   next multiplication.
Storage is unaffected (one orbit representative suffices); this is purely a
search-time orientation choice.

**Derived vs non-derived axis (user 2026-06-06).** Independently of buds, keep at
least **two representatives per format: one *non-derived*** (Perminov meta-flip /
ALS / imported — `atom:true`) **and one *derived*** (our recombination/composition
— `atom:false`), because "what each methodology can *reach*" is itself the
research signal (derived vs non-derived reachability), even when one ties or
beats the other on rank.
- **Manifest half — DONE (2026-06-06):** `shaveByBestAdditions` now splits its
  shave-group by `atom`, so a derived scheme that rediscovers an imported rank is
  no longer shaved as a worse-adds dup — one best-of-each-class survives per
  (format, field, rank). The `atom` flag (catalog.json) lets the SPA/comparator
  segment them.
- **Materialisation half — PENDING:** on disk we only *write* a derived scheme
  when it strictly beats the existing entry (`improveExisting`), so shapes where
  an import is already optimal have **no derived representative at all**. To make
  the comparison exist everywhere, add a "materialise the best derived scheme even
  on a tie" mode (closure/materialise) so each format carries a derived
  representative alongside its imported one. This is the on-disk side of "register
  2 schemes per shape".

### Phase 1.6 — `⟨2,3,3⟩ Z/2 r=14` SAT 🟡
Encoder ready (non-cubic refactor landed; `TestSatNonCubic#denseZ2_223_r11_isSat`
validates). Need to actually run `Phase16Runner 14` to either find a
new `R_{F₂}(⟨2,3,3⟩) ≤ 14` algorithm or certify UNSAT.
Estimated: minutes-to-hours with kissat.

### Phase 2 — `⟨3,3,3⟩ Z/2 r=22` SAT 🔴
Blocked on encoder scaling (kissat handles n=3 over-rank but Laderman
r=23 needs BreakID + embedding cubes). See
[`SOLVING_STRATEGIES.md`](paper/theory/SOLVING_STRATEGIES.md) §10.3.

### Constructive recombination — asymmetric allocations 🟡
**Status:** `Recombination.construct` lands the symmetric case (all base
multiplications share a sub-format). NOT yet:
- the asymmetric case (different sub-formats per base mult, exploiting
  AlphaTensor's `min(u1, w1)` reductions)
- catalog-driven sub-algorithm lookup that picks the lowest-rank
  available scheme per sub-format

### Padded Kronecker / trilinear aggregation 🟡
**Status:** rank-counting via `Recombination.recombine` works
(Sedoglavic 2017's `min` reduction). The brute-force `FindBestPaddedCompositions`
turned out to be the wrong abstraction; killed and superseded by the
targeted recipe approach below.

**Phase 1 (search):** `BlockSplitSearch.findBestSplit` (Java) evaluates
the closed-form identity over all splits via catalog rank lookups —
cheap, formula-only. `main` prints the per-target table.
**Phase 2 (construct):** `Compose.blockSplitCubic` builds the actual
factor matrices for a chosen split; validated on uniform splits
(`⟨4,4,4⟩=56`, `⟨6,6,6⟩=184`), both verify under
`Verifier.isExactNonCubic`.

NOT yet:
- **`NonCubicBilinearAlgorithm.permuteAxes(σ)`** — prerequisite for
  mixed splits. Given a scheme for `⟨n,m,p⟩` and `σ ∈ S_3`, returns the
  scheme for `σ(n,m,p)` via mechanical factor-matrix reshuffling.
  Without this, `blockSplitCubic` for `u ≠ v` fails because sub-product
  shapes like `⟨4,4,3⟩` and `⟨3,4,4⟩` (same canonical form, different
  ordering) can't both be served from a single catalog entry.
  ~50 lines, deterministic.
- **`section1/` and thin trivial schemes**: `⟨1,m,p⟩` / `⟨n,1,p⟩` /
  `⟨n,m,1⟩` need explicit entries (rank = `n·m`) so the sub-problem
  lookup is grounded for splits like `n = (n-1) + 1`.
- ~~**Constructive Sedoglavic-exact**~~ ✅ delivered via
  `Recombination.constructWithAllocation(Strassen, lookup, [u, v]³)`.
  Sedoglavic's saving = "Strassen-outer × min-shape inner" — not a
  separate algebraic identity. `⟨7,7,7⟩=250` verified end-to-end in
  `TestStrassenRecombination777`. Required a one-line bounds-check
  fix in `embedFactor` to skip padded-zero positions in non-uniform
  allocations.
- **Save the constructive 250-mult `⟨7,7,7⟩` to disk** — currently
  computed on demand; should land as
  `section7/composed-sedoglavic-strassen-4-3_7x7x7_r250.json` so it
  shows up in the catalog browser and Phase-1 lookup.
- **Generalise**: same recipe applies to any `⟨n,n,n⟩` with `n = u+v`
  by forcing `[u, v]³` allocation. A `BlockSplitSearch`-style enumerator
  that tries every `(u, v)` split, picks the min rank, and constructs
  the winning algorithm gives a full Phase-1+Phase-2 pipeline for
  cubic targets.

### Z/3-equivariant ALS — `r=22` exploration 🔴
ALS is stable (Z3Als fixed-triple fix landed). Next: actually push ALS
at `r=22` for `⟨3,3,3⟩` with SAT warm-starts (the A2 strategy in
[`SOLVING_STRATEGIES.md`](paper/theory/SOLVING_STRATEGIES.md) §10).

---

## Verification & infrastructure

### Big-format full verification 🟡
**Status:** sampled-verified `⟨16,16,16⟩` and `⟨32,32,32⟩` compositions
(50k / 10k random positions, all match). NOT yet:
- a once-off full-residual run on `⟨16,16,16⟩` (≈ minutes-to-hours;
  would give us a single conclusive 0-residual for the biggest
  reasonable target)

### Composed schemes for the larger sections 🟢
Currently the disk catalog stops at `section16/` and the only `⟨8,8,8⟩`
entries are composed. We could persist a curated set of large composed
schemes (e.g. `Strassen⁵`, `AlphaTensor² ⊗ Strassen`) to disk if their
file sizes don't blow up (sparse format helps a lot).

### ω-history intermediate entries 🟢
[`OMEGA_HISTORY.md`](paper/theory/OMEGA_HISTORY.md) §1.1 lists six refinements
between Coppersmith–Winograd 1990 and Williams 2024 (Stothers 2010,
Vassilevska Williams 2012, Le Gall 2014, Alman–Vassilevska Williams
2020, Duan–Wu–Zhou 2022). Add proper `[N]` entries in
[`REFERENCES.md`](REFERENCES.md).

### Feasible ω track 🟢
Schwartz & Zwecher 2025 is at `[19]` in references and §1.3 of
[`OMEGA_HISTORY.md`](paper/theory/OMEGA_HISTORY.md). Round out with Pan 1982
(`O(n^{2.773372})`), Strassen recursion at `O(n^{2.807})`, the
crossover analysis at `n_0 = 28`.

### Spotless ↔ JSON formatter alignment 🟢
Spotless handles Java + Markdown. JSON goes through
`ReformatSchemes` (custom matrix-friendly formatter). Consider wiring
`ReformatSchemes` into a Maven goal so it runs alongside
`spotless:apply` automatically (currently manual).

---

## Documentation

### Per-format detail pages 🟡
[`COVERAGE.md`](generated/COVERAGE.md) shows best-per-field per format. The
catalog browser ([`docs/index.html`](docs/index.html)) supports
filtering. Some users may want a per-format Markdown deep-dive showing
**all** entries for a given format chronologically — partly delivered
by the "history" column added today; could expand into a per-format
section with `<details>` collapsible blocks.

### Cross-link audit 🟢
Many MDs still reference each other via plain paths
(`SOLVING_STRATEGIES.md` § N) but don't always link explicitly.
Pass through and add explicit `[label](file.md)` everywhere.

### Specialised-input schemes (e.g. "with one zero entry") 🟡

Some classical results — Makarov 1970 ⟨5,5,5⟩=100 explicitly relies on
"⟨3,3,2⟩ Hopcroft-Musinski *with one matrix entry known to be zero* →
14 mults" — a non-generic scheme that exploits a structural zero in
one input. Our current catalog only stores **generic** ⟨n,m,p⟩
schemes (no precondition on inputs).

Consider a separate catalog category for `specialised-*` schemes:
- New scheme metadata: `"input_constraints": "{a31=0}"` (or similar).
- New filename convention: `*_3x3x2_r14_a*_zero-a31.json`.
- `Verifier` must respect the constraint (set the relevant entry to
  zero before the spot check).
- `FieldAwareLookup.find(n,m,p)` should NOT return these by default
  (they're not drop-in replacements for the generic ⟨n,m,p⟩).
- Add to a new `specialised-bounds.json` data file.

Practical motivation: closes the Makarov ⟨5,5,5⟩ materialisation
blocker (see FUTURE_WORK), plus enables future specialised-input
results (e.g. structured / sparse / symmetric input families).

### Chronology of small-matmul enumeration: Probert-Fischer → Smith → Islam → DIS09 🟡

The lineage of "enumerate upper bounds for small matmul by recursive
decomposition rules" runs through four key works that we should
explicitly chronologically anchor in REFERENCES.md and any future
intro material:

1. **Probert & Fischer 1976** — first systematic table of upper bounds
   for ⟨m,n,p⟩ via Permutation / Additive / Multiplicative / Zero-Padding
   rules. Manual search, up to ⟨40,40,40⟩. Their decomposition example
   for ⟨12,12,12⟩=1125 (M(⟨3,2,3⟩=15, A(⟨4,4,4⟩=49, ⟨4,2,4⟩=26))) is the
   archetype for what `BlockSplitSearch`/`ConcatSplitSearch` do today.
2. **Smith 2002** — computer search, expanded scope, full rectangular
   coverage to ⟨11,11,11⟩, squares to ⟨28,28,28⟩.
3. **Islam 2009** — MSc thesis synthesizing prior techniques (Strassen,
   Pan, Laderman, Makarov, Hopcroft-Kerr, Hopcroft-Musinski, Waksman,
   Probert-Fischer, Smith). Rich groundwork; chapter 3 contains
   explicit algorithm constructions (e.g. Makarov ⟨3,3,3⟩=22) and the
   decomposition framework used by DIS09.
4. **DIS09** (DeGroote-Ibarra-Sahni 2009 republication) — the table
   we benchmark against today.

Our `BlockSplitSearch`/`ConcatSplitSearch`/`PairFusedRecombination`
are direct descendants of this lineage but use better resolvers
(formula-aware SOTA) and find improvements unreachable by hand
search. ROADMAP item: add a "Historical chronology" section to
`REFERENCES.md` or a new `HISTORY.md` that traces this lineage so
new contributors understand the genealogy.

Islam 2009 is rich enough that a separate `references/islam2009/`
directory (paralleling `references/hopcroftkerr1971/`) with chapter
notes + transcribed algorithms would be worth standing up. The
Makarov ⟨3,3,3⟩=22 transcription attempt at this writing failed
verification (residual 4.9 against the formula given in §3.3.1 of
the thesis), suggesting either a typo in the thesis or transcription
error; cross-checking against a second source (Heun 1994, or the
original Russian Makarov paper) would resolve.

### Per-paper notes 🔴
For each reference in [`REFERENCES.md`](REFERENCES.md), expand the
annotation to summarize what we'd cite the paper FOR (algorithm
contribution, lower bound, methodology). Some entries are already
rich (e.g. AlphaEvolve); others are stubs.

### Commutative ω-history axis 🔴
We currently render an NC `ω`-history tab (`docs/omega-history.json`)
in the SPA — one ⟨∞,∞,∞⟩ row per published improvement. The
**commutative** counterpart `ω_c` is missing.

Commutative schemes lift to commutative-matmul recursion (entries
must commute, e.g., scalars over a commutative ring), so a
commutative cubic `⟨n,n,n⟩` at rank `r` gives `ω_c ≤ log_n(r)`.
Examples:
- Rosowski Thm 2/3 `⟨3,3,3⟩_c = 21` → `ω_c ≤ log_3(21) ≈ 2.771`
- Waksman 1970 `⟨n,n,n⟩_c = n(n²+n+ceil)/2` → asymptotic stays at 3
- Islam 2009 MSc thesis discusses the commutative-asymptotic
  framework explicitly (see `references/papers/islam_2009_msc_optim_matmul.pdf`)

Build:
1. `docs/omega-history-commutative.json` — analogous schema with
   commutative-only entries (Waksman 1970, Rosowski 2019,
   Makarov 1986, anything else from the literature).
2. Update the SPA so the ⟨∞,∞,∞⟩ section has two parallel
   sub-tabs (NC `ω` vs commutative `ω_c`), or one chart with
   two series and a legend.
3. Backfill: write a small generator that walks the catalog
   for commutative cubic schemes (`"commutative": true`), computes
   `log_n(rank)`, and emits an entry per (n, rank) tuple.

### Normalize REFERENCES.md / bibliography 🔴
Two-section structure:
1. **First section** — every reference on a single line in a
   standard bibliographic format (author, year, title, venue, DOI).
   Skimmable; the index.
2. **Second section** (keep the existing richer content) — full
   annotations, BibTeX, attribution notes.

For each entry, normalize its access tagging:
- `📄 Local PDF` — archived under `references/papers/`
- `🔗 External PDF` — open-access link (arXiv, institutional repo)
- `🔒 External PayWall` — link to publisher's gated page; we cite
  but cannot redistribute. Springer/Elsevier journals fall here.
This makes it obvious which sources a contributor can dig into
themselves vs. which require external access.

### Tutorial 🔴
A "getting started" tutorial that walks a new contributor through:
1. Loading a scheme.
2. Verifying it.
3. Composing it with another scheme.
4. Submitting a new scheme (verify, save, regenerate manifest).

### SymPy/marimo notebook for independent scheme verification + lineage replay 🔴
Build an executable notebook (sympy-based) that:
- Loads any scheme JSON in the catalog.
- Verifies it independently of the Java `Verifier` — symbolic expansion
  in sympy, exact rational arithmetic, cross-check against the
  matmul-tensor target. Catches Java-side bugs (compose arithmetic,
  rounding, etc.) by providing a totally separate stack.
- For schemes with a `lineage` field, replays the construction step
  by step: load each leaf, apply the operation (Kron / concat / recomb /
  pair-fused), show the intermediate matrix at each step, verify that
  the final result matches what's stored on disk.
- Renders the lineage tree visually (one cell per node) with
  collapsible sub-trees so a reader can drill from "this 32×32×32
  scheme of rank 14238" down to "Strassen × Sedoglavic-⟨11,11,11⟩".

**Recommended platform**: marimo (newer-than-jupyter notebook with
reactive cells, reproducible across kernels, exports to a single
self-contained HTML — perfect for catalog-browser linking). Plain
jupyter works too but is less embeddable.

**Why it's high-value**:
- Independent cross-check beyond the Java suite — particularly
  important now that the recursive materialiser is writing thousands
  of composed schemes.
- Pedagogical: a researcher reading the catalog can see exactly how
  ⟨14,14,14⟩=1719 is built without diving into Java source.
- Catalog-browser integration: each scheme row gets a "↗ verify in
  notebook" link that opens its rendered notebook.

**Effort**: ~1-2 days for a polished notebook with examples for the
main composition strategies; ongoing if we want notebook cells
auto-generated from lineage JSON.

---

## Publication

### Publishable article in LaTeX 🟡
**Status:** opened 2026-06-02. The catalog + frontier-closure search
have reached a state where a coherent write-up is possible. The
target shape is a short research paper (≈ 20–30 pages) covering:

- **Why a fresh catalog now.** The 2022–2026 small-format burst
  (AlphaTensor, AlphaEvolve, FMM-Lille, Perminov, Schwartz-Zwecher
  2025) has scattered results across notations, fields, and
  attribution conventions. A unified machine-checkable catalog is
  missing.
- **Field discipline.** Every rank claim names its field. The
  ⟨4,4,4⟩ landscape (47/F₂, 48/C, 49/R) is the canonical why.
  Commutative results form a separate axis (Waksman, Rosowski,
  Makarov 1986).
- **Composition strategies that the catalog uses.** Each gets a
  dedicated subsection: Strassen vs Winograd (structural
  difference, axis-flip relevance), axis-flip orbit and its
  asymmetry under mixed-size sums, block-additive (axis-concat)
  decomposition, multiplicative (Kronecker) composition,
  recombination = padded outer + sub-products, output peeling
  (DIS09 γ5), pair fusion (Pan TA, cyclic), disjoint-sum / τ-theorem
  (SZ 2025).
- **The non-overlap property** (see `docs/notes/materialisation-and-overlap.md`)
  — our generic recombination cannot rediscover sharing
  accidentally, which is the methodological boundary against
  hand-crafted schemes.
- **Frontier-closure search.** Two-phase algorithm (search + lazy
  materialise), in-memory overlay propagation between rounds,
  verification as a separate pass.
- **Comparison tables.** Analogous to DIS09 Tables 3 and 4 but
  refreshed and per-field, with separate Commutative columns. Each
  row carries an `attribution_for_rank` to avoid Cluster-A
  misattribution to bulk importers.

**Initial drop:** `paper/` skeleton with abstract, intro, strategy
section, tables stubs. The skeleton can be edited live; LaTeX
compiles from `paper/article.tex`.

**Tooling follow-up (separate task):** scripts to regenerate the
comparison tables from `docs/catalog.json` + `docs/cited-bounds.json`
+ `docs/derived-bounds.json` so the numbers in the paper track the
catalog automatically. Needs the SOTA definition clarified first —
"best known rank per (n, m, p, field, commutativity)" is the working
definition, but corner cases (border rank, mixed-field, partially
verified) need decisions written down.

**Why "publishable" and not just "documentation":** the catalog is
genuinely new wrt scope (covers F₂ alongside Q/Z/R/C, separates
commutative axis, includes Schwartz-Zwecher 2025), and the
frontier-closure-via-overlay algorithm has not been written up
elsewhere as far as the literature shows. The non-overlap-property
observation is also worth a section in its own right — it draws a
useful line between this catalog and the hand-crafted-discovery
literature it does NOT pretend to replace.

## How to contribute to this roadmap

Open a PR that:
- Either marks an item as 🚧 + adds your initials, OR
- Strikes it through with a link to the merged PR that completed it, OR
- Adds a new entry under the most relevant heading with a status emoji.
