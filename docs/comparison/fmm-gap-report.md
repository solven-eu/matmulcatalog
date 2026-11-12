# Catalog vs FMM-Lille digest — rank gaps

Pure data diff of `references/catalogs/fmm-lille-catalog.json` (FMM published NC ranks) vs our on-disk best (`FieldAwareLookup("R").findRank`). Non-trivial formats only.

- formats compared (we have a scheme): 5426
- **FMM strictly better than us: 53** (import/derive targets)
- we tie or beat FMM: 5373
- FMM has it, we have NOTHING: 0
- SKIPPED synthesized/unverified FMM bounds (empty page / HK-formula): 0

### split of the 53 FMM-better gaps by FMM provenance
- with a published REFERENCE (real import/cite targets): **7**
- NO reference = FMM-derived composition (re-derivable by us): **46**

## (A) FMM-better WITH a reference — published results to import (sorted by gap)

| shape | FMM | ours | gap | reference |
| --- | ---: | ---: | ---: | --- |
| ⟨14,14,15⟩ | 1715 | 1789 | 74 | Sidebottom:2026aa |
| ⟨14,18,19⟩ | 2829 | 2890 | 61 | Shukla:2026aa |
| ⟨10,19,24⟩ | 2723 | 2754 | 31 | Shukla:2026aa |
| ⟨15,16,16⟩ | 2202 | 2219 | 17 | Sidebottom:2026aa |
| ⟨13,14,14⟩ | 1565 | 1578 | 13 | Sidebottom:2026aa |
| ⟨13,15,15⟩ | 1763 | 1776 | 13 | Sidebottom:2026aa |
| ⟨14,15,15⟩ | 1842 | 1848 | 6 | Sidebottom:2026aa |

## (B) FMM-better with NO reference — FMM-derived, re-derivable by us (sorted by gap)

| shape | FMM | ours | gap |
| --- | ---: | ---: | ---: |
| ⟨14,14,29⟩ | 3308 | 3382 | 74 |
| ⟨14,14,30⟩ | 3430 | 3487 | 57 |
| ⟨14,15,29⟩ | 3557 | 3592 | 35 |
| ⟨9,27,29⟩ | 4053 | 4080 | 27 |
| ⟨9,24,29⟩ | 3576 | 3600 | 24 |
| ⟨14,15,28⟩ | 3430 | 3454 | 24 |
| ⟨16,20,28⟩ | 4936 | 4960 | 24 |
| ⟨6,21,29⟩ | 2179 | 2200 | 21 |
| ⟨9,21,29⟩ | 3139 | 3160 | 21 |
| ⟨22,30,30⟩ | 10534 | 10555 | 21 |
| ⟨23,29,30⟩ | 10721 | 10742 | 21 |
| ⟨14,14,19⟩ | 2247 | 2265 | 18 |
| ⟨16,16,31⟩ | 4438 | 4455 | 17 |
| ⟨16,20,29⟩ | 5256 | 5272 | 16 |
| ⟨14,14,17⟩ | 2016 | 2031 | 15 |
| ⟨14,14,31⟩ | 3609 | 3624 | 15 |
| ⟨14,14,25⟩ | 2947 | 2960 | 13 |
| ⟨14,14,27⟩ | 3158 | 3171 | 13 |
| ⟨9,26,29⟩ | 3982 | 3993 | 11 |
| ⟨6,23,29⟩ | 2450 | 2460 | 10 |
| ⟨6,24,29⟩ | 2496 | 2506 | 10 |
| ⟨24,29,30⟩ | 10824 | 10834 | 10 |
| ⟨6,26,29⟩ | 2767 | 2776 | 9 |
| ⟨9,17,30⟩ | 2711 | 2720 | 9 |
| ⟨21,28,32⟩ | 10134 | 10143 | 9 |
| ⟨9,22,29⟩ | 3400 | 3408 | 8 |
| ⟨20,24,25⟩ | 6466 | 6474 | 8 |
| ⟨22,29,31⟩ | 10701 | 10709 | 8 |
| ⟨26,29,32⟩ | 12318 | 12326 | 8 |
| ⟨8,18,29⟩ | 2517 | 2523 | 6 |
| ⟨8,27,30⟩ | 3744 | 3750 | 6 |
| ⟨9,25,29⟩ | 3837 | 3843 | 6 |
| ⟨14,15,19⟩ | 2399 | 2405 | 6 |
| ⟨6,22,29⟩ | 2353 | 2358 | 5 |
| ⟨22,28,29⟩ | 9654 | 9659 | 5 |
| ⟨15,16,17⟩ | 2442 | 2446 | 4 |
| ⟨7,9,25⟩ | 1020 | 1022 | 2 |
| ⟨14,15,26⟩ | 3253 | 3255 | 2 |
| ⟨7,7,17⟩ | 572 | 573 | 1 |
| ⟨7,8,17⟩ | 635 | 636 | 1 |
| ⟨7,9,19⟩ | 779 | 780 | 1 |
| ⟨9,19,25⟩ | 2541 | 2542 | 1 |
| ⟨9,23,28⟩ | 3393 | 3394 | 1 |
| ⟨11,11,17⟩ | 1313 | 1314 | 1 |
| ⟨11,13,13⟩ | 1204 | 1205 | 1 |
| ⟨14,14,21⟩ | 2491 | 2492 | 1 |

## FMM lists a shape we have no scheme for (first 200 by size)

| shape | FMM |
| --- | ---: |
