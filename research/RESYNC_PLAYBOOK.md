# Resync playbook — every upstream source, one pass

The catalog tracks six upstream sources. "Resync" means: refresh what they
publish, import every scheme we lack through the exact gates, close what is still
above them with our own operators, regenerate the committed artifacts, and read the
comparison reports until nothing upstream is below us (or the residue is audited).
The 3-hourly CI job (`sync-reference-catalogs.yml`) runs steps 1, 2 and the manifest
half of 5 on its own; the daily `regenerate-catalog.yml` runs 3 and 5; the weekly scan
(`scan-upstream-sources.yml`) runs step 8, opens an `upstream-gap` issue when a strict
gap exists and uploads the **closure input** (`upstream-gaps.json` + `upstream-gaps.txt`,
see step 8) with the run's artifact. **Step 4 — the closure — runs in no CI job**: it is
CPU-heavy, writes scheme files and needs a human read of what it produced. A gap the scan
detects therefore stays open until someone runs step 4 from the gap list (issue #23,
2026-11-09: 44 FMM rows re-priced upstream by `ProjOfAlongMode`, detected on time,
closed by nothing). A manual resync (PR #11 of 2026-09-30/10-12, this one of 2026-10-13)
is the same sequence plus the closure step.

## The sources

| source | what it publishes | channel | attribution |
|---|---|---|---|
| **FMM-Lille** (`fmm.univ-lille.fr`) | the index of best ranks per format (≤ 32), per-format pages with recipes, tensor artifacts `{shape}_tensor.mpl.bz2` | digest: `SyncReferenceCatalogs --fmm`; artifacts by hand: `ImportFmmTensorArtifacts` (parser `io.MapleTensorParser`); unverified index rows: `references/fmm-artifact-audit.json` | aggregator — `GenerateCatalogManifest` resolves the per-format citation |
| **Sedoglavic `fmm_digest`** | the same index as a JSON digest (`references/sedoglavic-fmm-sota.json`) | `ScanUpstreamSources --refresh` (both layouts; a digest that reads as 0 entries warns, never "0 gaps") | — |
| **Perminov** (`dronperminov/FastMatrixMultiplication`) | `status.json` (≤ 16, per Q/Z/ZT), `schemes/results/*` (his own explicit schemes incl. the serendipitous reaction bases), `schemes/known/<sub>/` (mirrors of others), the 17–32 serendipitous band (`references/perminov-serendipitous-17-32.json`, static) | `SyncReferenceCatalogs --perminov`; `ImportPerminovSchemes` (results/ only; keyed `(shape, rank, ZT⊂Z⊂Q)` from content; warns on an unmapped `known/<sub>`) | `PerminovKnownAttribution` (sub-folder → author; third-party files under `results/` listed explicitly) |
| **Witteveen 2026** (`MerlijnW70/fmm-schemes`, CC BY 4.0) | explicit ZT schemes ≤ 16 | `UpstreamSchemeRepos.WITTEVEEN` → `ImportUpstreamSchemeRepos` | `Witteveen 2026`, REFERENCES [90] |
| **Cloudroom Math 2026** (`matyasstoch/cloudroom-math`, Apache-2.0) | explicit schemes, mostly 17–32 | `UpstreamSchemeRepos.CLOUDROOM` → `ImportUpstreamSchemeRepos` | `Cloudroom Math 2026`, REFERENCES [91], issue #12 |
| **mizopl / we*star 2026** (`mizopl/we-star`, MIT) | explicit ZT/Z rank schemes ≤ 16 (`zt_rank_schemes/`) and 48 addition-reduced SLPs in Perminov's `_reduced` layout (`slp_addition_reduced/`); a Perminov-format `status.json` | `UpstreamSchemeRepos.MIZOPL` → `ImportUpstreamSchemeRepos` (a program is replayed with `SlpJson` and attached to the scheme we hold with the same matrices as `min_additions` + `slp`, or written as a new scheme); `ScanUpstreamSources` reads the `status.json` on the Q / Z / ZT axes | `mizopl 2026`, REFERENCES [92], issue #20; Perminov mirrors four under `known/mizopl/` (`PerminovKnownAttribution` → SKIP_FRESH_IMPORT) |
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
#    Input: the shape list the weekly scan uploaded (download it, or write it locally with step 8's
#    --shapes-out). Leading token NxMxP, then ours / upstream / gap / kind / sources — the kind
#    (from FMM's recipe cell) says what the row needs: closure-candidate:{projection,kron,block-sum,
#    serendipitous} → the operators below; base-candidate → import a reaction base (fmm-gap §3),
#    sweeping cannot reach it; upstream-unverified rows are not in the file.
gh run download <run_id> -n upstream-scan-<run_id> -D target/upstream-scan
cp target/upstream-scan/upstream-gaps.txt target/fmm-worse.txt
$M eu.solven.matmul.docs.explore.ProjectFmmGaps -Dexec.args="--shapes-file=target/fmm-worse.txt --passes=2"   # coordinate projection scatter
#    an FMM row that reads `ProjOfAlongMode(⟨parent⟩, mode)` is NOT a coordinate projection: it is Lemma 1 of
#    fmm.univ-lille.fr/projalongmode.html — drop one index of p / m / n by SUBSTITUTING the slice's variables,
#    which kills m / n / p products from ANY representative (the coordinate scatter above saw 0 on every axis
#    of the Kronecker stubs behind issue #23: ⟨9,18,18⟩:1600 → ⟨9,17,18⟩:1591 → ⟨9,17,17⟩:1574 upstream).
#    `catalog.SubstitutionProjection` + the `SubstProject` lineage node implement it; this driver replays it over
#    every one-larger parent on disk, in every axis order (p-then-m gives ⟨9,17,17⟩=1570, below FMM's 1574):
$M eu.solven.matmul.docs.migrate.MaterialiseSubstitutionProjections -Dexec.args="--shapes-file=target/fmm-worse.txt --execute"
#    then the concat/kron tier below, so block-sum rows (⟨9,17,k⟩ = ⟨9,17,18⟩ + ⟨9,17,k−18⟩, ⟨17,17,26⟩ = ⟨9,17,17⟩ +
#    ⟨17,17,17⟩) follow from the new parents; iterate the pair once more. Coordinate margins still add on top
#    where a factor REPRESENTATIVE carries a common rank-1 direction — rebased factors as margin bases:
$M eu.solven.matmul.docs.migrate.MaterialiseVectorProjections -Dexec.args="--shapes-file=target/fmm-worse.txt --execute"
$M eu.solven.matmul.docs.SchemeSweep -Dexec.args="--mode=materialize --field=Q --config=simple \
     --strategies=concat,kron,recomb --shape-file=target/fmm-worse.txt --threads=4"
#    `simple` is the cubic root pool (six cubic bases, no ⟨2,3,3⟩ / ⟨2,4,4⟩): its fixpoint says nothing
#    about rectangular bases. Run the rectangular tier, then the default pool (no --config =
#    defaultLarge: rectangular ≤ 5, cubic ≤ 8, PLUS every file the catalog's lineages already stand on,
#    whatever its dimension — pool = cap + proven-useful list, 2026-10-21, `TestProvenOuterBases`):
$M eu.solven.matmul.docs.SchemeSweep -Dexec.args="--mode=materialize --field=Q --config=rectangular \
     --baseShape=rectangular --strategies=recomb --shape-file=target/fmm-worse.txt --threads=4"
$M eu.solven.matmul.docs.SchemeSweep -Dexec.args="--mode=materialize --field=Q \
     --strategies=recomb --shape-file=target/fmm-worse.txt --threads=4"
#    a base newly registered in the root pool re-prices every stub standing on that shape (issue #8 put
#    two ⟨2,4,4⟩ reps in; the 360 ⟨2,4,4⟩-outer stubs were only ever optimised against one — re-sweeping
#    them gave 129 records in 20 s, 2026-10-21): list the shapes, then
#    grep -rlE '"lineage_compact": "R\[(2x4x4|4x2x4|4x4x2);' src/main/resources/schemes | xargs grep -h '"n": \[' \
#      | sed -E 's/[^0-9,]//g' | sort -u | awk -F, '{print $1"x"$2"x"$3}' > target/rerun-shapes.txt
$M eu.solven.matmul.docs.SchemeSweep -Dexec.args="--mode=materialize --field=Q --config=rectangular \
     --baseFilter=2,4,4 --strategies=recomb --maxNodes=60000000 --shape-file=target/rerun-shapes.txt --threads=12"
$M eu.solven.matmul.docs.migrate.MaterialiseSerendipitousWins -Dexec.args="--shapes=..."   # serendipitous products
#    an FMM "serendipitous" recipe we cannot price: the base REPRESENTATIVE is the variable —
#    walk it (persists a reaction base when it builds below the catalog), then materialise:
$M eu.solven.matmul.docs.explore.ProbeSerendipWalkGen -Dexec.args="16 20 28  4 5 7  4 4 4"
#    a row that survives: FMM artifact audit (download the artifact, read its rank → fmm-artifact-audit.json;
#    read the page recipe too — the index may be the recipe re-priced with today's ingredients)
#    stubs whose leaves are bare / legacy-stem refs replay LOWER once ingredients improved — restamp them
#    (2026-10-18: 510 loose stubs, 114 of them below the catalog-best their shape reported):
$M eu.solven.matmul.docs.migrate.RestampStubRanks -Dexec.args="--direction=all --execute"

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

# 8. the alert view (what the weekly job sees) + the closure input it hands over (step 4's shape file)
$M eu.solven.matmul.docs.explore.ScanUpstreamSources -Dexec.args="--refresh --out generated/UPSTREAM_SCAN.md \
     --gaps-out generated/upstream-gaps.json --shapes-out generated/upstream-gaps.txt"
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
| the lineage-ref audit kept its own grammar: 17 859 "UNRESOLVABLE" refs that all replayed (canonical stems, `TA_lita(n=N)`, `… :: CANONICAL` labels) and 4 "dangling" pins that were F₂ files under a C lookup | `AuditLineageRefs` delegates to `LineageReplayer.classifyLeaf` under the file's `fields[]`; `TestAuditLineageRefs` |
| the materialiser emitted canonical filename stems as UNPINNED leaves (419 distinct), so stub ranks moved with every import — 510 loose stubs on 2026-10-18, 114 hiding a record | `Lineage.atomFromFilename` pins `shape@hash7` (`TestLineageAtomFromFilename`); `RestampStubRanks --direction=all --execute` after every resync |
| the lineage JSON parser read a Recombination's `leaves` before its `base` while the writer renders base first — a shared-base back-reference became a literal `Atom("@ref?:L0")` on every read → write round trip (39 files after the 2026-10-19 restamp, 1 older) | `SchemeIO.readLineage` resolves back-references on demand; `TestLineageJsonSharedNodes` |
| the weekly scan's `fetch` wrote the download over the existing file without truncating (`BodyHandlers.ofFile(dst)` = CREATE+WRITE): a longer destination kept its old tail and the digest stopped parsing — would have crashed the job the day upstream shrinks the file | `TRUNCATE_EXISTING`; `TestScanUpstreamSources.fetch_truncates_a_longer_existing_file` |
| the replayer resolved a bare shape ref with `findFile`, which ranks commutative-only entries too (Rosowski ⟨2,2,n⟩ = 3n+1 below the NC best; unreplayable, so a warning-and-fallback hid it) | `FieldAwareLookup.findNcFile` in `resolveShape`, `classifyLeaf`, the cycle checks; `TestReplayLeafFallback.bare_shape_ref_resolves_to_the_non_commutative_best` |
| the default recombination pool capped non-cubic bases at dim 5 while committed stubs stood on dim-6/7 outer bases that only one-off `--base=` runs had offered (⟨3,4,6⟩, ⟨2,5,6⟩, ⟨3,4,7⟩, ⟨2,6,6⟩, …): the closure could not re-derive what the catalog held, so its "no improvement" fixpoint over those shapes was a statement about the pool; and the playbook's closure step ran `--config=simple`, six cubic bases | pool = cap + proven-useful list: every file the lineages name as a recombination base is admitted whatever its dimension (`ProvenOuterBases`, read from the catalog, not curated); `TestProvenOuterBases` (every outer base of every committed lineage is offered by the default pool); step 4 runs the rectangular tier and the default pool |
| the weekly scan crashed (digest download kept an old tail) while its `if: always()` steps posted the stale committed `generated/UPSTREAM_SCAN.md` from 2026-09-30 — issue #16 received "4151 strict gaps" twice (2026-10-26, 2026-11-02) for a 1-gap reality | the report file is deleted before the run and untracked; the summary and issue steps run only on `steps.scan.outcome == success` |
| the scan compared the char-0 (Q) rank only, so a ZT or Z record we lacked at a shape whose Q rank we matched was invisible — the ten mizopl ZT records of issue #20 (and Perminov's own Z ranks below ours) never alerted | Z and ZT axes per Perminov-layout source, from `fields[]` and `SchemeIO.ternaryOf` (the catalog's one ZT definition); `TestScanUpstreamSources.local_z_and_zt_axes_follow_fields_and_coefficients` |
| `ImportContributedSchemes.write` rewrote an identical existing file from its matrices ("re-stamping metadata only") — every stamped key (`min_additions`, `slp`, provenance) wiped on a re-import, the normal path of an addition-reduced program over held matrices | the existing file is kept and only missing keys are added; `TestImportContributedSchemes.write_over_an_identical_existing_file_keeps_its_addition_metadata` |
| a stored `slp` was write-only: nothing could read a program back, so a wrong schedule (or count) would have sat in the catalog unchallenged | `SlpJson` reads both layouts and `SchemeSlp.reconstructs` replays them; the importer refuses a program that does not reconstruct its matrices; `TestSlpJson` |
