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
