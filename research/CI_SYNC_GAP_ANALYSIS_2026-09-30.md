# Why the 2026-09 upstream refresh (KGP LITA) was not picked up by CI — root-cause analysis

**Symptom (2026-09-30).** A manual resync showed our catalog *worse* than
FMM-Lille on **529** shapes (114 815 multiplications in total; ⟨31,31,31⟩
14519 vs 13433, ⟨32,32,32⟩ 15079 vs 14197, …) and worse than
FMM ∪ Perminov on 1139, while the committed `references/fmm-cross-check.md`
still said "WORSE: 6" (its 2026-07-16 state) and the weekly upstream scan
had reported "Strict gaps: 0" on 2026-09-28. Every regular CI job was green.

## What actually happened upstream

Khoruzhii, Serafin, Gelß & Pokutta kept improving LITA in
[khoruzhii/lita](https://github.com/khoruzhii/lita): commits "Add LITA3"
(2026-08-01) … "Improve LITA7 odd" (2026-09-24), publishing explicit rational
`⟨N,N,N⟩` schemes for **13 ≤ N ≤ 32** (13³ = 1379 … 32³ = 14197, ω ≈ 2.7587).
FMM-Lille's index adopted them cube by cube — the committed FMM digest shows
⟨31,31,31⟩ at 14768 (2026-06-28) → 14519 (09-16) → 14216 (09-24) → **13433**
(09-29) — and gives every non-square shape of a family the cube's rank by
`Proj` (509 of the 529 "WORSE" rows are such inherited members; 20 are the
cube heads). Our `LitaTaConstruction` ports their *June* generator (N ≥ 19),
so every cube ≥ 13 fell behind and the whole 17–32 band followed.

## Layer by layer: what CI did and did not do

| Job (schedule) | Designed to | Did it? | Why not |
|---|---|---|---|
| `sync-reference-catalogs` (every 3 h) | refresh `references/catalogs/{fmm-lille,perminov}-catalog.json`; import **Perminov** scheme files | **Yes** — each FMM index move landed within hours; Perminov import ran ("0 imported, 2597 skipped-existing") | The only auto-importer is Perminov's, and it lists his own `schemes/results/*` only. He did mirror the KGP cubes 13³–16³ (`schemes/known/lita/`, 2026-09-18 and 09-26) — a `known/<sub>` folder no job lists — and his `status.json` stops at 16, so 17³–32³ are not there at all. FMM hosts Maple artifacts we have no CI importer for, and FMM's *citation* of a third-party repo is followed by nobody. |
| `regenerate-catalog` (daily + on push) | regenerate `docs/catalog.json` etc. | **Yes** — the manifest carried `external_best_rank` = FMM's 13433 every day, so the SPA *displayed* the gap | It regenerates only the three `docs/*.json`; the human-facing gap reports (`references/fmm-cross-check.md`, `references/reference-comparison.md`, `docs/comparison/fmm-gap-report.md`) were **not run by any workflow**, so the committed ones rotted at their July state. Nobody looks at 9 588 manifest rows for `external_best_rank < rank`. |
| `scan-upstream-sources` (weekly) | the **alert**: list ranks that strictly beat ours; "opens a tracking issue if any new gap appears" | **No** — "Strict gaps: 0" every week, including 2026-09-28 when the digest already had 14216 vs our 14519 | (a) `loadLocalRanks()` read "ours" from **filenames** with `[-_]NxMxP_[rm]R` — the legacy `{source}-{shape}_m{rank}` pattern. The 2026-06 catalog rename to `{n}x{m}x{p}-r{rank}-{note}-{hash7}` left it matching **7 of ~11 000** files, so every upstream shape became "MISSING" (a category nobody acts on), never a "GAP". No crash, a clean green run, an empty report — the canonical silent regression, and a direct violation of CLAUDE.md's "read metadata from content, never the filename". (b) The promised issue-opening step **did not exist**; the report only went to a run artifact and the job summary. |
| `verify-catalog` (daily) | verify every scheme on disk | Yes | Verification is about correctness, not competitiveness — by design. |
| (none) | *close* gaps (projection / sweeps) | — | No scheduled closure job exists; closing needs compute and judgment (the `/fmm-gap` playbook). Even a working scan could only have alerted. |

So the chain had **three independent breaks**: the alert was blind
(regex), the alert would have gone nowhere (no issue step), and the
committed reports that a human would read were never refreshed. The data
(digests, manifest) was correct throughout — the failure was entirely in
*surfacing*.

## The Perminov side — a second, separate failure (621 shapes)

Against Perminov's serendipitous 17–32 band (`perminov-serendipitous-catalog.json`,
971 formats, each naming the exact base file of its `s1 ⊗ˢ s2` recipe) we were
worse on **621** formats (11 465 multiplications) — again with a green sync job
("0 imported, 2597 skipped-existing"). Three causes, all silent:

1. **The importer's idempotence key was `(shape, rank)`.** Correct for rank
   results; wrong for `schemes/results/serendipitous_base/`, where several
   *content-distinct* schemes share a `(shape, rank)` (six `2x3x11_m55_*`) and
   differ in **bud structure**. At most one variant was imported — not
   necessarily the bud-rich one a recipe needs. Of 281 upstream bases we held 37.
   (Same lesson as the ⟨2,4,4⟩=26 representatives in issue #8: schemes at one
   `(shape, rank)` are not interchangeable.)
2. **The engine could not build the recipes anyway.** Most of that band is a
   *degenerate* product `base ⊗ˢ ⟨1,b,c⟩` — the second factor has a unit axis, the
   base's buds on that axis fuse into `⟨k,b,c⟩` blocks. Prediction priced the
   unit-axis inner (naive fallback) and the replayer resolved it, but the
   materialiser's build-time resolver returned "unavailable" for `⟨1,b,c⟩`
   (no catalog file: it *is* the naive scheme) and dropped every such candidate
   as unbuildable. `⟨13,20,21⟩ = ⟨13,4,3⟩:123 ⊗ˢ ⟨1,5,7⟩ = 3165` was unreachable
   (catalog: 3291) with the base on disk.

3. **The engine fused single-type buds only.** With every base on disk and the
   unit-axis inner buildable, 45 formats still sat above Perminov's rank. 41 of
   those recipes carry a *combined* bud in the published structure column
   (`references/perminov-serendipitous-17-32.json`): a `⟨1,2,2⟩` / `⟨2,1,2⟩` /
   `⟨2,2,1⟩` grid of four terms that fuses into ONE doubly-enlarged inner —
   `⟨9,9,18⟩ = 15·⟨6,3,3⟩ + 2·⟨12,3,3⟩ + ⟨6,6,6⟩:153 = 913`, where single-type
   buds give 920. `SerendipitousBudProduct`'s own class doc had it as "a
   follow-up". Two more needed a better *choice* of single-type buds than any of
   the six type orderings (a U-triple must give up a term to that term's only
   W-partner), and the last two were the recipes of two rational bases our exact
   verifier wrongly rejected (fix 14).

## The ≤ 16 band — a third-party dataset no importer listed (Witteveen 2026)

Merlijn S. Witteveen published `MerlijnW70/fmm-schemes` release 1.0 on
**2026-09-24**: ternary-integer schemes including a new best over every ring at
⟨11,13,15⟩ = 1364 (was 1371) and ⟨11,14,14⟩ = 1373 (was 1376); release 1.2
(09-30) added ⟨7,11,15⟩ = 772 (was 777). FMM-Lille cites the dataset for the
first two; Perminov mirrored release 1.1 on 09-29 08:40Z. Our 09:07Z sync ran
green, imported **two** of the 25 mirrored files, and credited both to Perminov.
Three more silent causes, all in `ImportPerminovSchemes`:

1. **Scope.** It lists `schemes/results/*` only — by design, since
   `schemes/known/<sub>/` is other people's work. But nothing else lists a *new*
   `known/<sub>`: `MerlijnW70_fmm_schemes/` (20 files, the two records among
   them) and `lita/` (the KGP cubes 13³–16³ as JSON) were mapped nowhere and
   imported by no job. The rank digest did move (`status.json` said 1364) — into
   reports no workflow regenerated (see above).
2. **The `(shape, rank)` key was blind to the coefficient class.** Perminov filed
   the five ⟨2,p,n⟩ ternary schemes under his own `results/ZT/`. Three were at a
   rank we already held as a *non-ternary* integer scheme — ⟨2,13,15⟩ = 304,
   ⟨2,13,16⟩ = 324, ⟨2,15,16⟩ = 374 — so "already have (shape, rank)" skipped
   them without downloading. Four of Perminov's own ZT records had been skipped
   the same way (⟨2,11,13⟩ = 221, ⟨2,11,14⟩ = 238, ⟨7,12,16⟩ = 878,
   ⟨9,12,13⟩ = 878). The key was also parsed from `-perminov_` *filenames*.
3. **Attribution by directory.** The two that did get in (⟨2,12,15⟩ = 280,
   ⟨2,14,16⟩ = 348) were stamped `source: "Perminov 2023"` because they sat
   under `results/` — five days after their author had published them.

## A third hider: the manual artifact audit never expired

`references/fmm-artifact-audit.json` (2026-07-09) marks FMM index ranks that
their published artifact did not back, and `FmmCrossCheck` moves those rows out
of WORSE ("upstream-unverified", deliberately excluded from every gap-closing
target list). The entries carried no record of *which* index value was audited,
so when FMM's index later moved — twelve of them to the KGP cube ranks, which
we hold and verified — the rows stayed hidden: ⟨27,28,28⟩ sat at 10442 vs 9847,
never targeted, and blocked the projection cascade to ⟨26,28,28⟩, ⟨25,28,28⟩, ….

## Fixes shipped with this analysis (branch `resync-perminov-fmm`)

1. **Scan made content-driven** — `ScanUpstreamSources.loadLocalRanks(Path)`
   reads `n`/`m`/`fields[]`/`commutative` from the JSON, excludes
   commutative-only and non-Q schemes (so an F₂ or Waksman rank can never
   "close" an NC char-0 gap). Regression guard
   `TestScanUpstreamSources` (new-convention names, no-token names,
   commutative/F₂ exclusion, and "the real catalog yields > 4000 shapes" — a
   regex-blind loader yields ~7).
2. **Scan alerts** — `scan-upstream-sources.yml` now sums the per-source
   "Strict gaps" counts and opens (or comments on) a single `upstream-gap`
   labelled issue with the report when the sum is > 0.
3. **Reports regenerated daily** — `regenerate-catalog.yml` also runs
   `FmmCrossCheck`, `CompareReferenceCatalogs`, `GenerateFmmGapReport`,
   commits the three reports, and prints the WORSE / BETTER headline in the
   job summary.
4. **KGP cubes imported** — `ImportKhoruzhiiLita` + the Python-free
   `io.NpzSchemeReader` (20 cubes, size-aware exact/spot-check gate; see
   REFERENCES.md [82] update). Re-runnable by hand when the repo moves.
5. **Projection closure repaired** — the scatter used an ambiguously-oriented
   parent ref (`NxMxP@hash` with two equal axes) and aborted the whole
   529-target campaign on the first predict/build divergence; it now pins the
   native `shape@hash` + exact-perm `Transpose`, skips a divergent candidate,
   and re-throws the collected divergences at the end.
6. **Perminov importer keyed by path for serendipitous bases** —
   `ImportPerminovSchemes` imports every content-distinct
   `serendipitous_base/` file into `bud-bases/` (hash-stamped reaction bases);
   the CI sync job therefore picks up new bases from now on. First run: 242
   imported (2 upstream files rejected by the exact verifier).
7. **Unit-axis inner in the serendipitous build** — one fallback
   (`⟨1,b,c⟩` → naive) in `RecursiveMaterialiser.trySerendipitous`; guard
   `TestSweepSpotsSota.compute_pipeline_reaches_13x20x21_3165_via_unit_axis_inner`.
8. **Audit entries expire** — each carries the `index_rank` it was audited
   against; `FmmCrossCheck` honours it only while the digest still shows that
   rank and lists the expired ones in the report.
9. **Witteveen's dataset imported, and pulled on a schedule** — the 26 schemes of
   release 1.2 are in `known/` (REFERENCES.md [90]); `ImportWitteveenSchemes`
   lists the *origin* repo and imports what the catalog lacks, as a step of
   `sync-reference-catalogs.yml` (every 3 h). Automated imports go through the
   contributed-scheme gate and are stamped `discovery: "TBD"`. It warns when
   the upstream listing matches no scheme file at all (a moved layout must not
   read as "nothing new").
10. **Class-aware, content-driven idempotence key** — `KnownSchemeKeys`:
    `(sorted shape, rank, ZT ⊂ Z ⊂ Q)` read from each file's `n` / `m` /
    `fields[]` / `zt` / coefficients, any source, any orientation. A ternary
    upstream scheme is no longer "already held" because an integer one is. First
    run: the four skipped Perminov ZT records imported. Guard `TestKnownSchemeKeys`.
11. **Unmapped `known/<sub>` folders are reported** — the importer names every
    sub-folder `PerminovKnownAttribution` does not map (log WARN, a GitHub
    `::warning::` annotation, the job summary). `MerlijnW70_fmm_schemes` and
    `lita` are now mapped (credited to their authors, pulled from their origin
    repos). Guard `TestImportPerminovSchemes`.
12. **Attribution of third-party files under `results/`** —
    `PerminovKnownAttribution.THIRD_PARTY_IN_RESULTS` (the five Witteveen paths;
    four confirmed by content hash against the origin files). ⟨2,14,16⟩ = 348 is
    re-attributed; the duplicate ⟨2,12,15⟩ = 280 file is removed in favour of the
    origin import. This one stays manual by nature: `status.json` names the path,
    not the author.

13. **Combined buds in the serendipitous product** — `GridBud` (`⟨1,b,c⟩`,
    `⟨a,1,c⟩`, `⟨a,b,1⟩`) detected, priced and built
    (`SerendipitousBudProduct.buildGridBlock`); the decomposition is now the
    cheapest of a structural family (grid strategies × type orderings × a
    fewest-options-first greedy), shared by search and replay. Math, table of
    mechanisms and the honesty tier (bound) in
    `references/SERENDIPITOUS_PARTIAL_PRODUCT.md` §6; guards in
    `TestSerendipitousGridBud` (three grid types against Strassen, scale
    absorption, ⟨9,9,18⟩ ≤ 913 built and exactly verified, four priced recipes).
14. **The exact verifier reads large denominators** — two published rational
    serendipitous bases (`2x4x6_m39_…_Q`, `2x5x7_m56_…_Q`) were rejected on every
    sync run, a `[FAIL]` line in a green job: `SymbolicVerifier` recovered
    denominators by scanning `d ≤ 1024`, so `3866/3705` read as irrational. It
    now falls back to continued fractions (den ≤ 10⁶, relative 1e-12) and scales
    numerators in `BigInteger` once the common denominator outgrows a double
    (⟨2,5,7⟩:56's is ~2·10¹⁷). Sound by construction: the integer identity is
    checked on the recovered fractions. Both bases are imported; the importer
    reports 0 failed. Guard `TestExactVerifierLargeDenominators`.

## Follow-ups (not in this PR)

- **CI-able KGP importer**: give `ImportKhoruzhiiLita` a `--download` mode
  (list `repos/khoruzhii/lita/contents/schemes`, import new/better
  `(shape, rank)` like `ImportWitteveenSchemes`) and run it in the sync job.
  With Perminov and Witteveen that would make three origin repos on a schedule;
  a fourth should become a small registry rather than a fourth class.
- **FMM Maple artifact importer in Java** (`MapleSchemeParser` exists for the
  17³ case; `tools/import_fmm_maple.py` needs Python) so a WORSE row backed by
  an FMM artifact can be pulled as a reaction base automatically.
- **Scheduled closure**: a weekly `ProjectFmmGaps --passes=2` over the WORSE
  list (projection is the cheap operator; it closed the inherited 509 here)
  with the wins committed like the Perminov import.
- **Guard the guard**: the scan's output should fail the job (or at least
  warn loudly) when the local-rank map is implausibly small — the June
  regression would have been caught the first Monday.
- **A cost-aware bud-partition search.** The serendipitous decomposition is the
  cheapest of a structural family, not an optimum over bud partitions (a weighted
  set packing). An optimiser that reads the ranks would need the chosen
  decomposition recorded in the `SerendipitousProduct` lineage node to stay
  replayable. Three-axis buds (`⟨a,b,c⟩`, all ≥ 2) are not searched either.
