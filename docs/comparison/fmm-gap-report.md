# Catalog vs FMM-Lille digest — rank gaps

Pure data diff of `references/catalogs/fmm-lille-catalog.json` (FMM published NC ranks) vs our on-disk best (`FieldAwareLookup("R").findRank`). Non-trivial formats only.

- formats compared (we have a scheme): 5426
- **FMM strictly better than us: 13** (import/derive targets)
- we tie or beat FMM: 5413
- FMM has it, we have NOTHING: 0
- SKIPPED synthesized/unverified FMM bounds (empty page / HK-formula): 0

### split of the 13 FMM-better gaps by FMM provenance
- with a published REFERENCE (real import/cite targets): **0**
- NO reference = FMM-derived composition (re-derivable by us): **13**

## (A) FMM-better WITH a reference — published results to import (sorted by gap)

| shape | FMM | ours | gap | reference |
| --- | ---: | ---: | ---: | --- |

## (B) FMM-better with NO reference — FMM-derived, re-derivable by us (sorted by gap)

| shape | FMM | ours | gap |
| --- | ---: | ---: | ---: |
| ⟨16,20,28⟩ | 4936 | 4974 | 38 |
| ⟨22,30,30⟩ | 10534 | 10555 | 21 |
| ⟨16,20,29⟩ | 5256 | 5272 | 16 |
| ⟨21,28,32⟩ | 10134 | 10143 | 9 |
| ⟨20,24,25⟩ | 6466 | 6474 | 8 |
| ⟨8,27,30⟩ | 3744 | 3750 | 6 |
| ⟨7,9,25⟩ | 1020 | 1022 | 2 |
| ⟨7,7,17⟩ | 572 | 573 | 1 |
| ⟨7,8,17⟩ | 635 | 636 | 1 |
| ⟨7,9,19⟩ | 779 | 780 | 1 |
| ⟨9,19,25⟩ | 2541 | 2542 | 1 |
| ⟨11,11,17⟩ | 1313 | 1314 | 1 |
| ⟨11,13,13⟩ | 1204 | 1205 | 1 |

## FMM lists a shape we have no scheme for (first 200 by size)

| shape | FMM |
| --- | ---: |
