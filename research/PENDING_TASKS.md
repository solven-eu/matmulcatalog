# Pending tasks

Snapshot of open work (pending + in-progress), grouped by theme. IDs match the
internal task tracker. Newly-surfaced items from recent sessions (no tracker ID
yet) are tagged **(new)**. Keep this in sync when items close.

## Search & discovery (rank/addition improvements)

- **#159** — Constructive serendipitous / Schönhage borrow-and-correct engine.
  *Largely advanced*: `--strategies=serendipitous`, `SerendipitousBudProduct`,
  `LineageVerifier`, 396 wins promoted (bands 16–32). Remaining: fold serendipity
  into the closure loop; tighten `serendipitousCost` to buildable ranks.
- **#196** — Unified MultisetFrontier search (multiset-first, cross-base, TA-aware):
  one cost+argmin engine shared by recombination **and** serendipity.
- **#198** — τ-identity certifier (grid-minus-k disjoint-sum), M0 = ⟨6,17,17⟩=1106.
  The structural mechanism FMM uses for the big near-cube gaps we don't close.
- **#144** — Enumerate peel positions beyond tail-only (γ5 coverage gap).
- **#169** — Proper Z₃ search (ALS-based) — Perminov's Z₃ gap is a research opening.
- **#197** — Border rank as a SEPARATE capability (not needed for serendipitous).
- **Composite reduction sweep + serendipitous concat (new)** — two reduction
  passes over composed schemes: (A) global proportional-product merge + DCE after
  Kron/concat/cascade (the FMM "kin-row" mechanism behind ⟨27,28,28⟩=10413);
  (B) cross-seam u-collision merging on concats (HK71 generalised via
  drop-and-repair; benchmark Smirnov ⟨3,3,6⟩=40 vs self-concat 46). First step:
  `ScanConcatSeamCollisions` collision-count scan (free ceiling on merge savings).
  Full write-up: `research/FUTURE_WORK.md` § "composite reduction sweep +
  serendipitous concat".
- **Projection scope for non-cubic targets (new)** — the projection closure default
  is cubic-only, so non-cubic shapes like ⟨29,30,30⟩ are never projected from a
  better parent (we hold ⟨30,30,30⟩=12688 < FMM's 12710 but carry ⟨29,30,30⟩=13902
  vs FMM's projected 12588). Widen projection to non-cubic shapes / run a
  projection pass over the non-cubic band.
- **BudEffectTable (new)** — savings-based bud scoring: tabulate
  `saving(type,k,inner)=k·R(inner)−R(enlarge)` from the catalog; replace the
  structural `budScore` with realizable serendipity value; gate base generation.
- **Incremental dirty-tracking for closure re-sweeps (new)** — replay only the
  shapes whose dependencies improved, instead of re-sweeping everything each round.

## Catalog curation, provenance & attribution

- **#121** *(in progress)* — Audit FMM-Lille-attributed schemes for upstream
  Smirnov / Sedoglavic / etc. attribution.
- **#199** — 3-catalog comparison (ours vs FMM vs Perminov) + `solven_discovery`
  flag. *(us-vs-FMM exists in `GenerateFmmGapReport`: 1888 FMM-better, 3218 tie-or-beat.)*
- **Prefer-derived-at-tie policy (new)** — when a derived (lineage-bearing) scheme
  *ties or beats* a raw lineage-less import, elect the derived one so attribution
  flows through the lineage (e.g. ⟨3,6,6⟩=80 = 2×⟨3,3,6⟩, not an FMM "discovery").
  Must NOT blanket-discard FMM (1888 formats still strictly better).
- **known/ vs ours/ droppable-rebuildable split (new — in ROADMAP)** — mirror the
  directory structure so `ours/` can be wiped and regenerated from `known/`
  deterministically; `--rebuild-ours` driver + content-hash equality check.
- **Linz ⟨3,3,3⟩ catalog evaluation (new — REFERENCES [85])** — mine the JKU-Linz
  matrix-multiplication page for bud-rich ⟨3,3,3⟩ bases / orbit diversity;
  incremental value over our existing ⟨3,3,3⟩ entries is undetermined.
- **#177** — Lineage regeneration: Kron-detection for "composed-X-Y" entries.
- **#200** — Stamp `atom:true|false` into individual scheme JSONs (from lineage).
- **#91** — Record S_{X,Y,Z} transform for axis-flip / permutation orbit-derived
  materialised schemes.

## Specific shapes / constructions

- **#126** — Deep-dive ⟨17,17,17⟩: our LRP import 2931 vs FMM published 2934.
- **#127** *(in progress)* — Materialise ⟨17,17,17⟩=2930 via Strassen-recombination at (9,8)³.
- **#112** — HK ⟨2,2,bridge-3⟩ constructive closure (extend impossibility proof to
  rank-2+ atoms, or find a positive construction).
- **#113** — HK2bc: n > 2p decomposition.
- **#168** *(largely done)* — ⟨2,2,2⟩=7 isotropy collapse + wording note on
  shape-minimisation goal. Resolved by exact GL₂(ℚ)³-orbit enumeration of
  recombination multisets: **40** distinct canonical multisets over ℚ (all
  rank-7 ⟨2,2,2⟩, seed-independent — an independent check of de Groote's
  single-orbit theorem), **6** ternary-realisable. Code: `RecombinationMultisetOrbit`
  (base-agnostic, e.g. ⟨2,3,3⟩), `Ternary2x2x2Orbit` (ternary + exact symbolic
  via per-axis decoupling); paper §"Recombination multisets of a base"
  (`paper/sections/multisets.tex`). Registered the 2 uncovered ternary bases in
  `curated/section2/` (`solven_orbit_*`), incl. a novel `3·⟨n₁,n₁,n₁⟩`
  triple-TA base. Remaining: certify the ternary count (6) beyond
  seed/alphabet-stability; the blind brute force is infeasible (39k⁷).

## Additions / SLP

- **#190** *(in progress)* — Additive-complexity evaluator: (near-)minimal additions.
- **#189** — Express additions as a straight-line-program `schedule[]` array (CSE);
  derive `scheduled_additions` from it.

## Verification & infrastructure

- **#68** *(in progress → effectively done)* — Stub-replay verification: superseded
  by the compositional `LineageVerifier` (verify leaves + trust operators). Close
  once wired everywhere it's needed.
- **Legacy stub-lineage frame corruption (new)** — `EnrichSchemeMetrics` (which
  replays every stub to stamp `verified`/`additions`/`buds`/`projection_margin`)
  exposed that ~3300 legacy `derived_recursive-*` stubs encode their concat
  operands in a *permuted axis frame*: an OLDER `RecursiveMaterialiser` named
  direct leaves from the **sorted** shape (`canon()`), so e.g.
  `ConcatCols(2x12x15, …)` for a ⟨12,15,17⟩ result really means a ⟨12,15,2⟩ left
  operand written rotated. The current emitter is already correct
  (`naive-%dx%dx%d` / `atomFromFilename`), so this is **stale on-disk data**, not
  a live bug. `LineageReplayer` was made robust to it (try exact-frame first, then
  orient operands to satisfy the concat precondition; `orientAs`/concat are
  correctness-preserving; final result oriented to the file's declared shape +
  verified). Result: replay errors **3297 → 621**, μ coverage **8056 → 10728 /
  11429 (93.9 %)**, and **0** recovered files stamp `verified:false` (all 8084
  recovered `derived_recursive` verify true). **Residual 621** are deeply-nested
  legacy trees where per-node orientation can't *globally* reconcile the frame
  (`robustConcat: no compatible orientation`, 538) + `construct: missing
  sub-algorithm` (59) + safety-net rejects (22) + 1 stray
  `ConcurrentModificationException` (a remaining race in a shared static cache,
  down from 4 after sharing one replayer in enrich). To close fully: either
  re-materialise these shapes with the current (correct) emitter, or add a global
  frame-inference solver to the replayer. Not blocking — the catalog carries μ for
  94 % of schemes and the gap is logged honestly.
  *RESOLVED (OrientAs):* the root cause was the materialiser recording sub-shape
  refs by `canon()` (sorted dims) — orientation-lossy — so even fresh lineage
  wasn't standalone-replayable for deep concat trees. Fixed by a new
  `Lineage.OrientAs(child,n,m,p)` node that records the exact reorientation by
  target dims (wired through replayer/verifier/parser/field-infer/bud-infer);
  `RecursiveMaterialiser` emits it at the disk-hit + memoized-reuse sites, so
  `replay(lineage).shape == alg.shape` holds by construction. Re-sweep + promote +
  enrich then stamped 315 (vs 32) with 0 new corrupted (vs 273); pruning the
  now-dominated stubs took **corrupted 328 → 17** and μ coverage to **99.1%**. The
  17 residual are the heaviest shapes (⟨28,30,32⟩ etc.) that OOM at -Xmx12g during
  the spot-check expansion — recoverable with more heap or the sparse backing.
  *Historical mitigation (still in place):* corrupted stubs are stamped `"corrupted": true`
  (+ `corrupted_reason`) by `EnrichSchemeMetrics` (auto-detected on replay
  failure, **self-healing** — cleared when a future replay succeeds). The flag is
  forwarded to `docs/catalog.json` and **`FieldAwareLookup` treats corrupted
  schemes as ABSENT for search gating** (`manifestCorrupted()`), so a corrupted
  stub's unverifiable rank claim can no longer shadow a fresh, verifiable
  discovery (the dual of the phantom-win bug). They remain listed in the
  catalog/SPA, flagged "rank claimed, not currently reproducible".
- **#128** *(in progress)* — Split `SchemeSweep` closure into search / materialise /
  verify stages.
- **Sparsity in `NonCubicBilinearAlgorithm` (new)** — large schemes are ~98% zeros
  ({−1,0,1}); a sparse backing (or sparse view) cuts memory + speeds sparse-exact
  verify/expand. Not on any critical path now that verify is compositional.

## Docs, references, paper

- **#129** *(in progress)* — Publishable article LaTeX (sections + tables).
- **#130** — Paper tables: tooling to regenerate from catalog + bounds JSON.
- **#194** — MermaidJS lineage graph from lineage JSON.
- **#195** — Tag each REFERENCES.md entry with its result/technique class.
- **#77** *(in progress)* — SPA ω-history rendering polish.

## Catalog closure guard (new)

- **One-level closure guard** *(test shipped)* — `TestCatalogOneLevelClosure`
  asserts no catalog scheme (char-0, non-commutative) is strictly beaten by a
  single concat (direct-sum) or Kronecker composition of catalogued sub-shapes,
  up to `GUARD_MAX_DIM=17` (the dense + first-stub band, currently fully closed).
  Encodes the "derive what we can; an import is representative only when our own
  derivation isn't equal-or-better" rule (the recurring ⟨7,7,9⟩-style report).
- **18–32 one-level-closure frontier (new)** — ~700 char-0 NC shapes (all
  maxDim ≥ 18) are strictly beaten by a one-level derivation the *bounded*
  recursive search doesn't reach (it doesn't compose arbitrary catalog-best
  halves). Needs a **direct one-level-composition closer** (apply the known
  concat/Kron of catalog bests, verify, write) rather than the search. Raising
  `GUARD_MAX_DIM` tracks progress. (Note: many earlier "2371 dominated" were
  illusory — they mixed commutative/F₂ sub-ranks that don't lift to NC.)
