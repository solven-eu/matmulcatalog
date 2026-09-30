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
| `sync-reference-catalogs` (every 3 h) | refresh `references/catalogs/{fmm-lille,perminov}-catalog.json`; import **Perminov** scheme files | **Yes** — each FMM index move landed within hours; Perminov import ran ("0 imported, 2597 skipped-existing") | The only auto-importer is Perminov's. KGP's cubes are not in Perminov's repo (his `status.json` is ≤ 16; the 17–32 band is a separate cited-bound digest), FMM hosts Maple artifacts we have no CI importer for, and FMM's *citation* of a third-party repo is followed by nobody. |
| `regenerate-catalog` (daily + on push) | regenerate `docs/catalog.json` etc. | **Yes** — the manifest carried `external_best_rank` = FMM's 13433 every day, so the SPA *displayed* the gap | It regenerates only the three `docs/*.json`; the human-facing gap reports (`references/fmm-cross-check.md`, `references/reference-comparison.md`, `docs/comparison/fmm-gap-report.md`) were **not run by any workflow**, so the committed ones rotted at their July state. Nobody looks at 9 588 manifest rows for `external_best_rank < rank`. |
| `scan-upstream-sources` (weekly) | the **alert**: list ranks that strictly beat ours; "opens a tracking issue if any new gap appears" | **No** — "Strict gaps: 0" every week, including 2026-09-28 when the digest already had 14216 vs our 14519 | (a) `loadLocalRanks()` read "ours" from **filenames** with `[-_]NxMxP_[rm]R` — the legacy `{source}-{shape}_m{rank}` pattern. The 2026-06 catalog rename to `{n}x{m}x{p}-r{rank}-{note}-{hash7}` left it matching **7 of ~11 000** files, so every upstream shape became "MISSING" (a category nobody acts on), never a "GAP". No crash, a clean green run, an empty report — the canonical silent regression, and a direct violation of CLAUDE.md's "read metadata from content, never the filename". (b) The promised issue-opening step **did not exist**; the report only went to a run artifact and the job summary. |
| `verify-catalog` (daily) | verify every scheme on disk | Yes | Verification is about correctness, not competitiveness — by design. |
| (none) | *close* gaps (projection / sweeps) | — | No scheduled closure job exists; closing needs compute and judgment (the `/fmm-gap` playbook). Even a working scan could only have alerted. |

So the chain had **three independent breaks**: the alert was blind
(regex), the alert would have gone nowhere (no issue step), and the
committed reports that a human would read were never refreshed. The data
(digests, manifest) was correct throughout — the failure was entirely in
*surfacing*.

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

## Follow-ups (not in this PR)

- **CI-able KGP importer**: give `ImportKhoruzhiiLita` a `--download` mode
  (list `repos/khoruzhii/lita/contents/schemes`, import new/better
  `(shape, rank)` like `ImportPerminovSchemes`) and run it in the sync job.
  Generalise: a small registry of "third-party scheme repos FMM cites"
  (KGP, `MerlijnW70/fmm-schemes`, …).
- **FMM Maple artifact importer in Java** (`MapleSchemeParser` exists for the
  17³ case; `tools/import_fmm_maple.py` needs Python) so a WORSE row backed by
  an FMM artifact can be pulled as a reaction base automatically.
- **Scheduled closure**: a weekly `ProjectFmmGaps --passes=2` over the WORSE
  list (projection is the cheap operator; it closed the inherited 509 here)
  with the wins committed like the Perminov import.
- **Guard the guard**: the scan's output should fail the job (or at least
  warn loudly) when the local-rank map is implausibly small — the June
  regression would have been caught the first Monday.
