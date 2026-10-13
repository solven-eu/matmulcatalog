# FMM artifact audit — index rank vs published tensor artifact (2026-07-08)

**Method:** for every row of the cross-check `WORSE` table (44 rows,
digest of 2026-07-08), download FMM-Lille's published tensor artifact
(`https://fmm.univ-lille.fr/{shape}_tensor.mpl.bz2`) and count its actual
rank (triads, via `tools/import_fmm_maple.py`'s parser — analysis only,
nothing lands in the catalog trees). Compare: FMM index rank (what the
cross-check charges us against) vs our catalog rank vs the artifact rank.

**Why this matters:** the cross-check treats FMM's *index* as ground truth.
This audit shows the index is NOT always backed by a published scheme —
14 of 44 rows are **phantom gaps** (their artifact is no better than our
catalog, sometimes far worse), and ~5 more rest on recipe-only placeholder
files with no tensor data at all. Re-check whenever the digest refreshes.

## Classes

### PHANTOM — index unattained; artifact ≥ ours (14 rows, NOT real gaps)

| shape | index | ours | FMM artifact | note |
| --- | ---: | ---: | ---: | --- |
| ⟨2,12,18⟩ | 333 | 334 | 334 | index = HK71 formula; our task-#9 theorems say locally impossible |
| ⟨2,14,21⟩ | 452 | 453 | 453 | ditto |
| ⟨2,16,24⟩ | 588 | 590 | 592 | artifact WORSE than ours |
| ⟨2,18,27⟩ | 743 | 745 | 747 | artifact worse than ours |
| ⟨2,20,30⟩ | 915 | 918 | 920 | artifact worse than ours |
| ⟨2,24,30⟩ | 1095 | 1096 | 1104 | artifact worse than ours |
| ⟨21,28,30⟩ | 9473 | 9477 | 9782 | artifact 305 above index |
| ⟨21,28,32⟩ | 10134 | 10143 | 10477 | |
| ⟨22,23,23⟩ | 6435 | 6476 | 6501 | |
| ⟨22,28,28⟩ | 9383 | 9411 | 9552 | |
| ⟨22,30,30⟩ | 10534 | 10555 | 10705 | |
| ⟨27,28,28⟩ | 10413 | 10442 | 10442 | artifact == ours exactly |
| ⟨27,28,29⟩ | 11169 | 11198 | 11198 | artifact == ours exactly |
| ⟨27,28,30⟩ | 11561 | 11590 | 11590 | artifact == ours exactly |

### RECIPE-ONLY — placeholder artifact, no tensor data (≥3 rows)

`⟨16,20,28⟩`, `⟨16,20,29⟩`, `⟨20,24,25⟩`: the `.mpl` contains only symbolic
matrix declarations ending `Tensor:=Tensor:` (~10 KB). The index rank rests
on the page's construction recipe (serendipitous formulas — arithmetic
verifies against FMM's own ingredient ranks, so the ranks are plausibly
valid, but no explicit scheme is published; our engine can't execute the
span-compressed fusion, see the 2026-07-08 span-bud log entry).

### PARAMETRIC artifacts — resolved 2026-07-09, all PHANTOM (5 more rows)

These artifacts carry SYMBOLIC coefficients (`-a61/a41`,
`(a11*a32-a12*a31)/(a21*a32-a22*a31)` …) — parametric scheme families the
numeric importer rejects. Rank counted structurally (`Matrix(` count − 3)/3,
calibrated exact on ⟨2,12,18⟩=334:

| shape | index | ours | FMM artifact | note |
| --- | ---: | ---: | ---: | --- |
| ⟨27,28,31⟩ | 12072 | 12101 | 12117 | artifact worse than ours |
| ⟨27,28,32⟩ | 12339 | 12368 | 12394 | artifact worse than ours |
| ⟨27,29,31⟩ | 12890 | 12919 | 12942 | artifact worse than ours |
| ⟨27,29,32⟩ | 13187 | 13201 | 13234 | artifact worse than ours |
| ⟨4,18,30⟩ | 1394 | 1395 | 1397 | artifact worse than ours |

**Updated tally: 19 PHANTOM, 3 recipe-only placeholders, 22 attains-index.**

### CORRECTION (2026-07-09, user catch): the two "base-level mismatches" were OUR importer's mirror-preference artifacts, NOT FMM errors

`tools/import_fmm_maple.py` prefers the Perminov GitHub mirror and falls back
to FMM only on 404. The claimed base-level mismatches ⟨4,9,10⟩ "250 vs 255"
and ⟨4,7,8⟩ "161 vs 164" were counts of PERMINOV'S mirror copies; FMM's own
artifacts, re-downloaded directly, are **250 and 161 — exactly their index
ranks** (⟨4,9,10⟩ page correctly credits Perminov Feb 2026). Retract both
from the report. Provenance re-check of all 14 non-parametric PHANTOM rows:
every one 404s on the mirror, so their audit counts came from FMM's own
files — the phantom classifications STAND (the 5 parametric rows were
curl'd from FMM directly). Residual note: Perminov's own
`known/tensor/4x9x10_tensor.mpl` (255) lags the 250 FMM credits him with —
an upstream-of-upstream curiosity, not an FMM issue.

### REVISION (2026-07-09): "phantom" ≠ "index wrong" — two are now REPRODUCED

The INDEX-page `<td content="description">` recipes (not audited above — the
detail pages were) back several "phantom" rows. Two are now reproduced by our
own engine after the ambiguous-orientation fix (⟨7,7,5⟩-style dims-repeat
masking in budBasesAt): **⟨21,28,30⟩ = 9473** (exact tie; their artifact 9782
is just stale) and **⟨4,18,30⟩ = 1394** (tie; beats their artifact 1397). So
the correct report language for FMM is "artifact stale / missing", not "index
over-claim", wherever an index recipe verifies. The ⟨2,·,·⟩ HK-formula six
and ⟨16,20,28⟩ remain UNVERIFIED index claims: for ⟨16,20,28⟩ the recipe
`(⟨4,5,7⟩−17)⊗⟨4,4,4⟩+8·⟨4,4,8⟩` is term-count inconsistent (17 terms cannot
fit 8 size-2 blocks; slice-rank analysis of their published ⟨4,5,7⟩:104 —
7×U2+1×U3, every class slice-rank-tight — proves no partition prices 4944;
best expressible is 4989).

### ATTAINS-INDEX — real gaps, explicit schemes exist (≈25 rows)

Artifact rank ≤ index (sometimes better, e.g. ⟨8,27,30⟩ artifact 3736 vs
index 3744). These are genuine: FMM publishes explicit factor matrices at
ranks our engines currently cannot derive (absorbing-pad / overlap
constructions, per the 2026-07-08 gap-tree census). Full rows in
`results2.txt` of the audit run; key examples:
⟨3,29,29⟩ 1840 (ours 1843), ⟨17,17,19⟩ 3266 (3267), ⟨19,19,22⟩ 4536 (4538),
⟨12,20,23⟩ 3184 (3224), ⟨20,23,23⟩ 5906 (5945), ⟨20,20,27⟩ 6006 (6022).

## Action items

1. Cross-check honesty: annotate/split PHANTOM rows in `FmmCrossCheck`
   output (index-only over-claims are not actionable gaps and overstate
   our deficit). Consider a periodic artifact audit driver in
   `docs.verify` replacing the ad-hoc shell run.
2. Report the 14 phantom rows upstream (Sedoglavic) — stale artifacts or
   optimistic index entries.
3. The 25 real rows are the true frontier: dominated by absorbing-pad /
   overlap constructions (analyze the now-downloaded artifacts to extract
   the device — they are explicit!).

## Expiry rule (2026-09-30)

An audit entry is a statement about **one index value** ("FMM's index said R and
their artifact did not back it"). Each entry in `fmm-artifact-audit.json` now
carries `index_rank` — the FMM index rank it was audited against — and
`FmmCrossCheck` honours the entry **only while the digest still shows that rank**.
When the index moves, the row returns to the ordinary WORSE list (and is named in
the report's "expired artifact-audit entries" line) until someone re-audits it.

Why: twelve of the July entries (⟨27,28,x⟩, ⟨27,29,x⟩, ⟨22,23,23⟩, ⟨22,28,28⟩,
⟨22,30,30⟩, ⟨16,20,28/29⟩) kept hiding rows whose FMM index had since dropped — most
of them to the Khoruzhii–Serafin–Gelß–Pokutta LITA cube ranks, which we hold and
verified. ⟨27,28,28⟩ sat at 10442 vs 9847, was never a closure target, and blocked
the projection cascade to ⟨26,28,28⟩ / ⟨25,28,28⟩. As of 2026-09-30 only
⟨20,24,25⟩ (placeholder), ⟨21,28,32⟩ and ⟨22,30,30⟩ (index_only) are still valid at
their audited index value. When adding an entry, always record `index_rank`.

## Re-audit 2026-09-30 — the 13 rows left after the resync

After the KGP / Witteveen imports and the closure campaigns the cross-check had 13
WORSE rows. Each FMM artifact (`{shape}_tensor.mpl.bz2`) was downloaded and its rank
read twice — the `i=1..R` of the file's own verification line and the `Triad(` count
(they agree on every file).

| shape | FMM index | ours | FMM artifact | verdict |
| --- | ---: | ---: | ---: | --- |
| ⟨6,14,25⟩ | 1320 | 1322 | 1346 | index only — artifact 24 above ours |
| ⟨7,7,17⟩ | 572 | 573 | 573 | index only — artifact == ours (a parametric family, free symbol `l`) |
| ⟨7,8,17⟩ | 635 | 636 | 637 | index only — artifact above ours |
| ⟨7,9,19⟩ | 779 | 780 | 780 | index only — artifact == ours |
| ⟨7,9,25⟩ | 1020 | 1022 | 1022 | index only — artifact == ours |
| ⟨9,19,25⟩ | 2541 | 2542 | 2542 | index only — artifact == ours |
| ⟨11,11,17⟩ | 1313 | 1314 | 1314 | index only — artifact == ours |
| ⟨11,13,13⟩ | 1204 | 1205 | 1205 | index only — artifact == ours |
| ⟨7,11,30⟩ | 1493 | 1495 | **1493** | **backed** — imported (`ImportFmmTensorArtifacts`) |
| ⟨13,19,29⟩ | 4248 | 4249 | **4248** | **backed** — imported |
| ⟨8,27,30⟩ | 3744 | 3750 | 3736 (**invalid**) | parses to 3736 triads that do **not** compute ⟨8,27,30⟩ — 95 of 240 output cells wrong on random input, exact term (0,9,18) absent |
| ⟨16,20,28⟩ | 4936 | 4974 | placeholder | recipe-only; the page's recipe gives 4944 |
| ⟨16,20,29⟩ | 5256 | 5272 | placeholder | recipe-only; the page's recipe gives 5264 |

So of the 13: **8 are not gaps** (FMM's index is one or two below what FMM itself
publishes; the per-shape page still states the artifact's rank), **2 were real and
are closed by importing the artifact**, **1 rests on an artifact that does not
verify** (`artifact_invalid`, a new audit class — the July audit decoded a valid
3744-triad artifact for this format; the current 3736-triad file is broken), and
**2 are recipe-only**. The eight, the one and the two are recorded in `fmm-artifact-audit.json` with the index value audited, so they
return to WORSE the day the index moves again.

The two recipe-only rows are a real open gap against the *recipe*, not against the
index: FMM's page builds ⟨16,20,28⟩ as `(⟨4,5,7⟩:104 − 17) ⊗ ⟨4,4,4⟩:48 + 8·⟨4,4,8⟩:96 = 4944`,
and our engine reaches 4974. Seventeen base terms are absorbed into eight doubled
blocks — more than any partition into class buds allows (see the 2026-07-08
span-compressed-bud entry in `research/DISCOVERIES_PENDING_ANALYSIS.md`); the
combined buds added on 2026-09-30 do not reach it either.

The pages and the digest index disagree on eight rows, which is itself worth
reporting upstream: the digest's `rank` for those formats is below the rank in the
page title and in the artifact.

## Re-audit 2026-10-13 — the 13 upstream-unverified rows, read from the pages

The 2026-10-13 resync (PR #14) left exactly the 13 audited rows above FMM's index.
This pass asked, for each one, what FMM's per-format page *constructs today* (recipe
text and, where it exists, the tensor artifact), and whether our operators can reach
it. Pages were re-fetched on 2026-10-13; the 2026-09-30 artifact readings are
unchanged.

| shape | FMM index | ours (before → after) | what the page constructs | reading |
| --- | ---: | ---: | --- | --- |
| ⟨16,20,28⟩ | 4936 | 4967 → **4960** | serendipitous `(⟨4,5,7⟩:104 − 17) ⊗ ⟨4,4,4⟩:48 + 8·⟨4,4,8⟩:96 = 4944`; artifact placeholder | the index is the recipe re-priced with today's ⟨4,4,8⟩ = 95 (87·48 + 8·95 = 4936) — no hidden construction; the recipe itself does not balance under bud semantics (see below) |
| ⟨22,30,30⟩ | 10534 | 10555 | artifact 10705 = 3·⟨10,15,15⟩:1435 + 4·⟨12,15,15⟩:1600 (7-block Strassen recombination) | our 10555 is the *same* recombination with ⟨10,15,15⟩ = 1385 (FMM's own index value); the index sits 21 below anything the page or the artifact describes |
| ⟨16,20,29⟩ | 5256 | 5272 | `⟨16,20,1⟩:320 + ⟨16,20,28⟩:4944 = 5264` | follows ⟨16,20,28⟩ (index = 320 + 4936); our 5272 is `R[4x4x4; 4,4,4,4 \| 5,5,5,5 \| 8,7,7,7]`, and 320 + 4960 = 5280 does not beat it |
| ⟨21,28,32⟩ | 10134 | 10143 | artifact 10477 (15-block recombination of ⟨9–12, 9–10, 9–12⟩ leaves) | we beat the published construction by 334; the index is 343 below the page — unexplained |
| ⟨20,24,25⟩ | 6466 | 6474 | serendipitous `(⟨5,12,5⟩:204 − 42) ⊗ ⟨4,2,5⟩:32 + ⟨16,2,5⟩:126 + 2·⟨4,6,5⟩:90 + 16·⟨4,4,5⟩:61`; artifact placeholder | balanced bud arithmetic (σ = 2 + 2·6 + 16·3 = 62) — needs a ⟨5,12,5⟩:204 representative with one V-quad, two W-triples and sixteen W-pairs; every published representative (FMM's own artifact, Perminov's two) has eighteen pairs and nothing larger → 6474 |
| ⟨8,27,30⟩ | 3744 | 3750 | serendipitous `(⟨4,9,10⟩:250 − 12) ⊗ ⟨2,3,3⟩:15 + 6·⟨4,3,3⟩:29`; artifact 3736, **invalid** | balanced (six V-pairs, σ = 1 each) — needs a ⟨4,9,10⟩:250 representative with six B-side pairs; FMM's own ⟨4,9,10⟩ artifact and Perminov's file carry six *U*-pairs (⟨2,3,6⟩ = 30, σ = 0) → 3750 |
| ⟨7,9,25⟩ | 1020 | 1022 | recombination summing to 1022 | page = ours |
| ⟨7,7,17⟩ | 572 | 573 | `2·⟨3,4,8⟩:73 + ⟨3,4,9⟩:83 + ⟨3,3,9⟩:63 + ⟨3,4,8⟩:73 + 2·⟨4,4,9⟩:104 = 573` | page = ours; every ingredient is at FMM's index value, so the 572 is not ingredient drift |
| ⟨7,8,17⟩ | 635 | 636 | recombination summing to 637 | we beat the page by 1 |
| ⟨7,9,19⟩ | 779 | 780 | recombination summing to 780 | page = ours |
| ⟨9,19,25⟩ | 2541 | 2542 | `⟨9,19,7⟩:780 + ⟨9,19,18⟩:1762 = 2542` | page = ours |
| ⟨11,11,17⟩ | 1313 | 1314 | recombination summing to 1314 | page = ours |
| ⟨11,13,13⟩ | 1204 | 1205 | recombination summing to 1205 | page = ours |

**Two kinds of index/page disagreement.** For ⟨16,20,28⟩ and ⟨16,20,29⟩ the digest
is simply the page's recipe re-priced with today's ingredient ranks (the page text is
frozen at generation time; `⟨4,4,8⟩` went 96 → 95). Every other disagreement is
genuine: the digest rank is below anything the page constructs and below the
artifact, with no ingredient drift to explain it (⟨7,7,17⟩'s ingredients are all at
index value). Those remain upstream claims.

**What closed, and how.** The serendipitous recipes hinge on the *representative*
of the base, not on its rank. A metaflip walk over the base at fixed rank
(`ProbeSerendipWalkGen`: objective = the structural family's product cost against
the recipe's inner, tie-broken by a per-axis σ-potential) found a ⟨4,5,7⟩:104
representative with eight U-pairs, two U-triples, four V-pairs, six W-pairs and two
⟨2,1,2⟩ grids — `bud-bases/section7/4x5x7-r104-flipwalk_serendip-4c05684.json`,
attributed as a reaction base only. Against ⟨4,4,4⟩ the family prices it at
`80·48 + 6·95 + 2·95 + 2·180 = 4960` (built, exact), and the same base gives four
further records: ⟨8,15,21⟩ = 1548 (FMM 1558), ⟨8,25,28⟩ = 3302 (3319),
⟨12,21,30⟩ = 4146 (4160), ⟨16,21,30⟩ = 5572 (5592). The same walk on ⟨5,5,12⟩:204
(eight seeds, 200k steps each) never left the eighteen-pair plateau, and on
⟨4,9,10⟩:250 never produced a B-side pair — those two rows stay open as
representative-search problems, not as engine gaps.

**⟨16,20,28⟩'s recipe does not balance.** `(104 − 17) ⊗ 48 + 8·⟨4,4,8⟩` realises
seventeen base products with sixteen product slots (8 doubled blocks). Under bud
semantics a class of k same-U products is a ⟨1,1,k⟩ sub-tensor of rank k, so eight
⟨4,4,8⟩ blocks cover sixteen products; a seventeenth can only be absorbed if some
class's V- or W-forms are linearly dependent — which would rewrite the base at rank
103, contradicting ⟨4,5,7⟩ = 104 being FMM's own best. Either the construction uses
a device beyond buds (a product-level fusion of the kind decoded in the ⟨8,27,30⟩
artifact on 2026-07-09) or the printed accounting is the generator's, not a built
scheme's — the artifact is a placeholder, so this cannot be settled from public data.
With seven pairs and one triple priced honestly the recipe would give
`87·48 + 7·95 + ⟨4,4,12⟩:141 = 4982`, above our 4960.
