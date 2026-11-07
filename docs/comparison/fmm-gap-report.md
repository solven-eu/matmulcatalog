# Catalog vs FMM-Lille digest — rank gaps

Pure data diff of `references/catalogs/fmm-lille-catalog.json` (FMM published NC ranks) vs our on-disk best (`FieldAwareLookup("R").findRank`). Non-trivial formats only.

- formats compared (we have a scheme): 5426
- **FMM strictly better than us: 57** (import/derive targets)
- we tie or beat FMM: 5369
- FMM has it, we have NOTHING: 0
- SKIPPED synthesized/unverified FMM bounds (empty page / HK-formula): 0

### split of the 57 FMM-better gaps by FMM provenance
- with a published REFERENCE (real import/cite targets): **0**
- NO reference = FMM-derived composition (re-derivable by us): **57**

## (A) FMM-better WITH a reference — published results to import (sorted by gap)

| shape | FMM | ours | gap | reference |
| --- | ---: | ---: | ---: | --- |

## (B) FMM-better with NO reference — FMM-derived, re-derivable by us (sorted by gap)

| shape | FMM | ours | gap |
| --- | ---: | ---: | ---: |
| ⟨9,17,17⟩ | 1574 | 1600 | 26 |
| ⟨17,17,26⟩ | 4297 | 4323 | 26 |
| ⟨16,20,28⟩ | 4936 | 4960 | 24 |
| ⟨6,21,23⟩ | 1779 | 1800 | 21 |
| ⟨9,21,23⟩ | 2499 | 2520 | 21 |
| ⟨11,18,21⟩ | 2499 | 2520 | 21 |
| ⟨15,21,23⟩ | 4139 | 4160 | 21 |
| ⟨21,23,28⟩ | 7382 | 7403 | 21 |
| ⟨21,23,30⟩ | 7894 | 7915 | 21 |
| ⟨21,23,31⟩ | 8102 | 8123 | 21 |
| ⟨22,30,30⟩ | 10534 | 10555 | 21 |
| ⟨16,20,29⟩ | 5256 | 5272 | 16 |
| ⟨15,17,21⟩ | 3145 | 3160 | 15 |
| ⟨12,17,18⟩ | 2148 | 2160 | 12 |
| ⟨12,17,21⟩ | 2508 | 2520 | 12 |
| ⟨21,23,27⟩ | 7174 | 7184 | 10 |
| ⟨9,17,18⟩ | 1591 | 1600 | 9 |
| ⟨9,17,19⟩ | 1744 | 1753 | 9 |
| ⟨9,17,20⟩ | 1829 | 1838 | 9 |
| ⟨9,17,21⟩ | 1933 | 1942 | 9 |
| ⟨9,17,22⟩ | 2020 | 2029 | 9 |
| ⟨9,17,23⟩ | 2112 | 2121 | 9 |
| ⟨9,17,24⟩ | 2151 | 2160 | 9 |
| ⟨9,17,25⟩ | 2292 | 2301 | 9 |
| ⟨9,17,26⟩ | 2373 | 2382 | 9 |
| ⟨9,17,27⟩ | 2477 | 2486 | 9 |
| ⟨9,17,32⟩ | 2933 | 2942 | 9 |
| ⟨9,22,23⟩ | 2706 | 2715 | 9 |
| ⟨9,23,23⟩ | 2821 | 2830 | 9 |
| ⟨11,17,18⟩ | 2059 | 2068 | 9 |
| ⟨13,17,18⟩ | 2395 | 2404 | 9 |
| ⟨13,21,23⟩ | 3755 | 3764 | 9 |
| ⟨17,18,27⟩ | 4622 | 4631 | 9 |
| ⟨21,28,32⟩ | 10134 | 10143 | 9 |
| ⟨22,23,31⟩ | 8331 | 8340 | 9 |
| ⟨23,23,32⟩ | 8830 | 8839 | 9 |
| ⟨20,24,25⟩ | 6466 | 6474 | 8 |
| ⟨6,15,23⟩ | 1265 | 1272 | 7 |
| ⟨15,17,18⟩ | 2705 | 2712 | 7 |
| ⟨8,27,30⟩ | 3744 | 3750 | 6 |
| ⟨9,23,29⟩ | 3551 | 3557 | 6 |
| ⟨13,14,14⟩ | 1579 | 1582 | 3 |
| ⟨14,14,27⟩ | 3172 | 3175 | 3 |
| ⟨15,16,16⟩ | 2220 | 2223 | 3 |
| ⟨16,16,31⟩ | 4456 | 4459 | 3 |
| ⟨7,9,25⟩ | 1020 | 1022 | 2 |
| ⟨12,23,24⟩ | 3788 | 3790 | 2 |
| ⟨14,14,25⟩ | 2960 | 2962 | 2 |
| ⟨7,7,17⟩ | 572 | 573 | 1 |
| ⟨7,8,17⟩ | 635 | 636 | 1 |
| ⟨7,9,19⟩ | 779 | 780 | 1 |
| ⟨8,15,23⟩ | 1713 | 1714 | 1 |
| ⟨9,17,28⟩ | 2579 | 2580 | 1 |
| ⟨9,19,25⟩ | 2541 | 2542 | 1 |
| ⟨11,11,17⟩ | 1313 | 1314 | 1 |
| ⟨11,13,13⟩ | 1204 | 1205 | 1 |
| ⟨21,22,31⟩ | 7619 | 7620 | 1 |

## FMM lists a shape we have no scheme for (first 200 by size)

| shape | FMM |
| --- | ---: |
