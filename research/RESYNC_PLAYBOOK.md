# Resync playbook — every upstream source, one pass

The catalog tracks six upstream sources. "Resync" means: refresh what they
publish, import every scheme we lack through the exact gates, close what is still
above them with our own operators, regenerate the committed artifacts, and read the
comparison reports until nothing upstream is below us (or the residue is audited).
The 3-hourly CI job (`sync-reference-catalogs.yml`) runs steps 1–5 on its own; the
weekly scan (`scan-upstream-sources.yml`) runs step 8 and opens an `upstream-gap`
issue when a strict gap exists. A manual resync (PR #11 of 2026-09-30/10-12, this
one of 2026-10-13) is the same sequence plus the closure step.

## The sources

| source | what it publishes | channel | attribution |
|---|---|---|---|
| **FMM-Lille** (`fmm.univ-lille.fr`) | the index of best ranks per format (≤ 32), per-format pages with recipes, tensor artifacts `{shape}_tensor.mpl.bz2` | digest: `SyncReferenceCatalogs --fmm`; artifacts by hand: `ImportFmmTensorArtifacts` (parser `io.MapleTensorParser`); unverified index rows: `references/fmm-artifact-audit.json` | aggregator — `GenerateCatalogManifest` resolves the per-format citation |
| **Sedoglavic `fmm_digest`** | the same index as a JSON digest (`references/sedoglavic-fmm-sota.json`) | `ScanUpstreamSources --refresh` (both layouts; a digest that reads as 0 entries warns, never "0 gaps") | — |
| **Perminov** (`dronperminov/FastMatrixMultiplication`) | `status.json` (≤ 16, per Q/Z/ZT), `schemes/results/*` (his own explicit schemes incl. the serendipitous reaction bases), `schemes/known/<sub>/` (mirrors of others), the 17–32 serendipitous band (`references/perminov-serendipitous-17-32.json`, static) | `SyncReferenceCatalogs --perminov`; `ImportPerminovSchemes` (results/ only; keyed `(shape, rank, ZT⊂Z⊂Q)` from content; warns on an unmapped `known/<sub>`) | `PerminovKnownAttribution` (sub-folder → author; third-party files under `results/` listed explicitly) |
| **Witteveen 2026** (`MerlijnW70/fmm-schemes`, CC BY 4.0) | explicit ZT schemes ≤ 16 | `UpstreamSchemeRepos.WITTEVEEN` → `ImportUpstreamSchemeRepos` | `Witteveen 2026`, REFERENCES [90] |
| **Cloudroom Math 2026** (`matyasstoch/cloudroom-math`, Apache-2.0) | explicit schemes, mostly 17–32 | `UpstreamSchemeRepos.CLOUDROOM` → `ImportUpstreamSchemeRepos` | `Cloudroom Math 2026`, REFERENCES [91], issue #12 |
| **KGP LITA** (`khoruzhii/lita`, no licence file) | the cubes 13³–32³ as `.npz` archives | `ImportKhoruzhiiLita --download --execute` (archives vendored under `external/khoruzhii-lita/`, stubs on `TA_lita_npz(n=N,r=R)`) | `Khoruzhii 2026`, REFERENCES [82] |

Contributions that arrive as a GitHub issue (issues #7, #8, #12) go through
`ImportContributedSchemes --meta=references/contributions/<x>.meta.json` and get a
record under `references/contributions/`.

## The pass

```bash
# 0. environment (this repo's JDK/Maven setup); every long step: background + heap dump flags
M="mvn -q -ntp exec:java -Dexec.mainClass="

# 1. digests
$M eu.solven.matmul.docs.verify.SyncReferenceCatalogs -Dexec.args="--perminov"
$M eu.solven.matmul.docs.verify.SyncReferenceCatalogs -Dexec.args="--fmm"

# 2. explicit schemes we lack (each importer: content-keyed idempotence, exact gate, dry-run available)
$M eu.solven.matmul.docs.migrate.ImportPerminovSchemes            # results/*; warns on new known/<sub>
$M eu.solven.matmul.docs.migrate.ImportUpstreamSchemeRepos        # every registry repo
$M eu.solven.matmul.docs.migrate.ImportKhoruzhiiLita -Dexec.args="--download --execute"

# 3. where are we still above an upstream rank?
$M eu.solven.matmul.docs.verify.FmmCrossCheck                     # references/fmm-cross-check.md: WORSE / UNVERIFIED / BETTER
$M eu.solven.matmul.docs.verify.CompareReferenceCatalogs          # references/reference-comparison.md: vs min(FMM, Perminov)
#    (Perminov's 17–32 band: ProbeFindRank over the 971 formats of perminov-serendipitous-17-32.json)

# 4. close the WORSE rows with our operators (loop until no progress)
$M eu.solven.matmul.docs.explore.ProjectFmmGaps -Dexec.args="--passes=2"      # projection scatter
$M eu.solven.matmul.docs.SchemeSweep -Dexec.args="--mode=materialize --field=Q --config=simple \
     --strategies=concat,kron,recomb --shape-file=target/fmm-worse.txt --threads=4"
$M eu.solven.matmul.docs.migrate.MaterialiseSerendipitousWins -Dexec.args="--shapes=..."   # serendipitous products
#    an FMM "serendipitous" recipe we cannot price: the base REPRESENTATIVE is the variable —
#    walk it (persists a reaction base when it builds below the catalog), then materialise:
$M eu.solven.matmul.docs.explore.ProbeSerendipWalkGen -Dexec.args="16 20 28  4 5 7  4 4 4"
#    a row that survives: FMM artifact audit (download the artifact, read its rank → fmm-artifact-audit.json;
#    read the page recipe too — the index may be the recipe re-priced with today's ingredients)

# 5. regenerate the committed artifacts (CI order)
$M eu.solven.matmul.docs.generate.GenerateCatalogManifest
$M eu.solven.matmul.docs.generate.GenerateCitedBounds
$M eu.solven.matmul.docs.generate.GenerateDerivedBounds
$M eu.solven.matmul.docs.verify.FmmCrossCheck
$M eu.solven.matmul.docs.verify.CompareReferenceCatalogs
$M eu.solven.matmul.docs.generate.GenerateFmmGapReport

# 6. audit automated imports (discovery TBD → verdict) once the source's own table has been read
$M eu.solven.matmul.docs.migrate.RestampContributedProvenance -Dexec.args="--meta=references/contributions/<x>.meta.json --execute"

# 7. guards, then commit everything the steps wrote (schemes, external/, digests, docs/*.json, reports)
mvn -q test -Dtest='TestSweepSpotsSota,TestFmmArtifactAudit,TestScanUpstreamSources'

# 8. the alert view (what the weekly job sees)
$M eu.solven.matmul.docs.explore.ScanUpstreamSources -Dexec.args="--refresh --out generated/UPSTREAM_SCAN.md"
```

"Done" = `FmmCrossCheck` WORSE 0, every UNVERIFIED row audited at its current index
value, `CompareReferenceCatalogs` worse rows ⊆ the audited ones, the Perminov band
worse 0, the scan's strict gaps 0 — and a `research/DISCOVERIES_PENDING_ANALYSIS.md`
entry with the numbers.

## Adding a source

- A GitHub repo of Perminov-style JSON files: one `UpstreamSchemeRepos` entry (path
  pattern, label, `source`, citation), a REFERENCES entry, a `docs/catalog.js`
  source → reference mapping, and — if Perminov mirrors it — a
  `PerminovKnownAttribution` case. The importer's unmapped-subtree warning is how
  such a mirror announces itself.
- Any other format: a parser in `io/` and an importer in `docs.migrate` that ends in
  `ImportContributedSchemes.gate` / `write` (the `.npz` and Maple-tensor importers are
  the templates).
- Rank claims without matrices: `GenerateCitedBounds` (`docs/cited-bounds.json`).

## What silently broke before, and the guard that now exists

| failure | guard |
|---|---|
| scan read "ours" from filenames (7 of 11k matched) → "Strict gaps: 0" for 3 months | `TestScanUpstreamSources` (content-driven, > 4000 shapes) |
| scan parsed Sedoglavic's old map layout only → 0 entries, "0 gaps" (2026-10-09 layout change) | both layouts + 0-entries warning; `sedoglavic_digest_is_read_in_both_layouts` |
| importer keyed by `(shape, rank)` from filenames → ZT records and content-distinct bases skipped | `KnownSchemeKeys` + `TestKnownSchemeKeys` |
| new `known/<sub>` mirrors invisible (Witteveen, LITA, Cloudroom) | unmapped-subtree warning + `TestImportPerminovSchemes` |
| audit entries never expired → rows hidden after the index moved | `index_rank` on every entry, `FmmArtifactAudit` + `TestFmmArtifactAudit` |
| exact verifier rejected denominators > 1024 → two bases "FAIL" every run | `TestExactVerifierLargeDenominators` |
| re-stamp "succeeded" without writing (`updateFields` never clobbers) | `TestRestampContributedProvenance` |
| a stub's recorded rank drifted below what its lineage replays (legacy/bare refs follow the catalog-best, serendipitous nodes re-derive, a pinned stub pins a lineage not matrices) — two false claims found 2026-10-14, nothing re-read `m` after write | `VerifySchemes --check=rank` (daily schedule), `TestVerifySchemes.rank_check_flags_a_stub_recorded_below_its_replay`; repair with `RestampStubRanks` |
