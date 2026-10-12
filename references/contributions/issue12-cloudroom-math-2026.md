# Issue #12 — Cloudroom Math Campaign (2026): 27 explicit schemes below the catalog — all loaded, `discovery: true`

**Source.** GitHub issue
[solven-eu/matmulcatalog#12](https://github.com/solven-eu/matmulcatalog/issues/12)
(opened 2026-10-09 by [@matyasstoch](https://github.com/matyasstoch)), pointing at
the dataset [matyasstoch/cloudroom-math](https://github.com/matyasstoch/cloudroom-math)
(release 2026-10-08, **Apache-2.0**, `CITATION.cff` author "Cloudroom Math
Campaign"): 27 explicit schemes in Perminov's JSON layout, with per-result
manuscripts, a verifier (`verification/check_all.sh`) and a baselines table
(`verification/baselines.md`) naming the previous best of *every* source per
format. The search, the code and the write-ups were done by AI agents (Claude
Opus 5.5, orchestrated in Cloudroom); every scheme was verified exactly before
release. REFERENCES.md [91]. Processed 2026-10-11/12 (this note is the record).

**What was contributed.** Three result families:

- `⟨4,6,15⟩ = 261` (was 263): a serendipitous product of Perminov's `⟨2,3,15⟩:75`
  base with `⟨2,2,1⟩` where an exact change of basis *inside* two bud classes
  exposes two **hidden 2×3 grids**, so two blocks use `R(⟨2,4,3⟩) = 20` instead
  of 24. This is `SerendipitousBudProduct.GridBud` (combined buds, 2026-09-30)
  one step further: our grids must be visible in the given basis; theirs are
  found after an exact re-basing of the free factors — a well-defined engine
  follow-up ("hidden grids").
- `⟨9,11,11⟩ = 710` (was 711): a mod-3 flip-graph walk on Perminov's rational 711
  (denominators 4 and 8 reduce mod 3) and a core-restricted exact Hensel lift of
  the 48 changed terms.
- 25 formats in the 17–32 band: serendipitous products of existing schemes with
  best known rational block ranks, six with a hidden grid. No new small scheme.

## Outcome in one table

All 27 pass the contributed-scheme gate (`ImportContributedSchemes`: exact
BigInteger symbolic proof, random spot-check, F₂/F₃ reductions from the
coefficients, canonical-disk round-trip) and are loaded as
`schemes/known/section{maxDim}/…-cloudroom_2026-*.json` (PR #11, 2026-10-11).
"Catalog" columns are the issue's own (this catalog at its 2026-10-05 data
commit); the last column is the gate's verdict against the catalog **at import
time** (2026-10-11), after the combined-bud serendipitous passes of 2026-10-01/02.

| shape | catalog `rank` (10-05) | catalog `external_best` (10-05) | FMM page (10-08) | **new** | coefficients (issue) | fields stamped here | at import (10-11) |
|---|---:|---:|---:|---:|---|---|---|
| ⟨4,6,15⟩ | 263 | 263 | 263 | **261** | Z[1/2] | F3 Q R C | **below** (263) |
| ⟨4,6,27⟩ | 470 | 469 | 470 | **468** | Z | F2 F3 Z Q R C | tie (468) |
| ⟨4,9,22⟩ | 550 | 549 | 550 | **548** | Z | F2 F3 Z Q R C | **below** (549) |
| ⟨4,9,30⟩ | 748 | 746 | 748 | **742** | Z[1/2] | F3 Q R C | **below** (744) |
| ⟨4,12,30⟩ | 959 | 957 | 960 | **954** | Z[1/2] | F3 Q R C | **below** (957) |
| ⟨4,20,26⟩ | 1331 | 1328 | 1349 | **1326** | Z[1/2] | F3 Q R C | tie (1326) |
| ⟨6,9,22⟩ | 745 | 741 | 756 | **737** | Z[1/2] | F3 Q R C | **below** (741) |
| ⟨6,10,24⟩ | 912 | 910 | 912 | **900** | Z[1/2] | F3 Q R C | tie (900) |
| ⟨6,10,32⟩ | 1234 | 1234 | 1234 | **1231** | Z[1/2] | F3 Q R C | tie (1231) |
| ⟨6,14,20⟩ | 1079 | 1079 | 1079 | **1073** | Z[1/2] | F3 Q R C | tie (1073) |
| ⟨6,14,24⟩ | 1240 | 1240 | 1264 | **1232** | Z[1/2] | F3 Q R C | tie (1232) |
| ⟨6,16,26⟩ | 1548 | 1528 | 1554 | **1527** | Z | F2 F3 Z Q R C | tie (1527) |
| ⟨7,15,20⟩ | 1357 | 1347 | 1357 | **1346** | Z | F2 F3 Z Q R C | tie (1346) |
| ⟨8,8,24⟩ | 985 | 982 | 985 | **981** | Q | Q R C | tie (981) |
| ⟨8,9,22⟩ | 1000 | 997 | 1005 | **994** | Z | F2 F3 Z Q R C | **below** (997) |
| ⟨8,10,26⟩ | 1328 | 1328 | 1333 | **1326** | Z[1/2] | F3 Q R C | tie (1326) |
| ⟨8,10,28⟩ | 1418 | 1413 | 1429 | **1410** | Z[1/2] | F3 Q R C | tie (1410) |
| ⟨8,11,32⟩ | 1769 | 1754 | 1771 | **1752** | Q | Q R C | tie (1752) |
| ⟨8,13,24⟩ | 1556 | 1544 | 1568 | **1541** | Z | F2 F3 Z Q R C | tie (1541) |
| ⟨8,13,28⟩ | 1800 | 1800 | 1800 | **1796** | Z[1/2] | F3 Q R C | tie (1796) |
| ⟨8,13,32⟩ | 2065 | 2064 | 2065 | **2056** | Q | Q R C | tie (2056) |
| ⟨8,14,24⟩ | 1667 | 1666 | 1682 | **1655** | Q | Q R C | tie (1655) |
| ⟨8,14,28⟩ | 1946 | 1940 | 1946 | **1928** | Z[1/2] | F3 Q R C | tie (1928) |
| ⟨9,11,11⟩ | 711 | 711 | 715 | **710** | Z[1/2] | F3 Q R C | **below** (711) |
| ⟨9,14,20⟩ | 1568 | 1564 | 1575 | **1563** | Z | F2 F3 Z Q R C | tie (1563) |
| ⟨9,14,28⟩ | 2147 | 2133 | 2153 | **2121** | Z[1/2] | F3 Q R C | tie (2121) |
| ⟨12,14,21⟩ | 2133 | 2126 | 2133 | **2121** | Z[1/2] | F3 Q R C | tie (2121) |

**Discovery status: `true` for all 27.** On the release date each rank was below
every published rank — FMM-Lille's live page, this catalog's `rank` and
`external_best_rank` (the latter already folds in Perminov's ≤ 16 tables and his
2026 serendipitous bounds) — per the issue's table and the dataset's
`verification/baselines.md`, re-checked here. Stamped by
`RestampContributedProvenance --meta=references/contributions/cloudroom-math-2026.meta.json`
(2026-10-12); the per-file `discovery_note` carries the claim and the audit.

**The 20 ties are two independent searches landing on the same rank.** Between
2026-09-30 and 10-02 the combined-bud serendipitous passes (PR #11, branch
`resync-perminov-fmm`, merged 2026-10-12) reached exactly these 20 ranks here —
e.g. ⟨6,14,20⟩=1073, ⟨6,14,24⟩=1232, ⟨8,13,28⟩=1796, ⟨9,14,28⟩=2121 — from the
same mechanism (grids of bud classes against best block ranks), while the
dataset's 2026-10-05 comparison baseline was `master`, which did not have them.
The release of 2026-10-08 is the earlier *publication*; the catalog's own stubs
stay as independent derivations with their lineage. Where both exist, the
explicit import and the stub are different schemes at the same rank.

**Ranks are ℚ-field claims** (the `fields` column is what the gate proved: the
seven integer files reduce to F₂ and F₃; the sixteen `Z[1/2]` files to F₃ only — a
denominator 2 is not invertible mod 2; the four `Q` files (denominators divisible
by 3 as well) to neither). Optimality tier of every rank: **bound**.

## Reproduce

```bash
# gate (dry-run), then load — the verdicts above are the gate's log lines
git clone https://github.com/matyasstoch/cloudroom-math /tmp/cloudroom
mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.ImportContributedSchemes \
    -Dexec.args="/tmp/cloudroom/schemes --meta=references/contributions/cloudroom-math-2026.meta.json --field=Q"
# audit verdict (discovery TBD → true), content-selected by source
mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.RestampContributedProvenance \
    -Dexec.args="--meta=references/contributions/cloudroom-math-2026.meta.json --execute"
```

Guard: `TestSweepSpotsSota.retains_issue12_cloudroom_records` (≤ on five of the
ranks). Perminov mirrors two of the files under `schemes/known/cloudroom-math/`
(mapped in `PerminovKnownAttribution`, `SKIP_FRESH_IMPORT`: we pull from the
origin). A scheduled pull from the origin repo, like `ImportWitteveenSchemes`,
is a follow-up; so is the hidden-grid re-basing in the engine.
