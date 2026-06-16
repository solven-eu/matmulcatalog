# Attribution audit — Smirnov 2017 "Several Bilinear Algorithms" ⟨3, P, Q⟩

Reference: [REFERENCES.md §[70]](../REFERENCES.md#70-smirnov2017).

This file lists, for each rank claim of Smirnov 2017, the on-disk
scheme JSONs in `src/main/resources/schemes/` that hit that rank, and
checks whether their `attribution_for_rank` field correctly points
back to **Smirnov 2017** (or an even earlier source). It does **not**
modify any scheme JSON — it is purely a tracking document.

Convention used in the table below:
- `OK` — the JSON already credits Smirnov 2017 (or an explicitly
  earlier source that *also* established the bound).
- `MISSING` — the JSON has no `attribution_for_rank` field at all.
- `WRONG` — the JSON credits a *later* re-importer (AlphaTensor 2022,
  Perminov, AlphaEvolve 2025, meta-flip-graph, Moosbauer flips, …)
  for a rank Smirnov 2017 already published in January 2017.
- `SUPERSEDED-OK` — the JSON's rank is **lower** than Smirnov 2017's
  claim (a later result genuinely improved the bound), and the
  attribution to the later source is therefore *correct* — Smirnov
  2017 should NOT be credited for the strictly-better bound.

Per project policy ("Don't bypass safety in git operations" / "Don't
be selective about improvements only" / "Distinguish discoveries from
re-discoveries"), the `WRONG`/`MISSING` rows are the ones that should
be re-tagged in a follow-up PR — **this PR does NOT modify any scheme
JSON**.

## Shape `R⟨3,4,6⟩` — Smirnov 2017 claims rank 56

| file | rank | current `attribution_for_rank` | verdict |
| --- | --: | --- | --- |
| `section6/alphaevolve_3x4x6_r54_a538_0.5xZ.json` | 54 | (AlphaEvolve `source`) | SUPERSEDED-OK (54 < 56) |
| `section6/perminov-ZT_3x4x6_r54_a700.json` | 54 | — (assumed AlphaEvolve / Perminov reduction) | SUPERSEDED-OK |
| `section6/perminov-ZT_3x4x6_r56_a359.json` | 56 | `"Moosbauer (symmetric flips)"` | **WRONG** — should credit Smirnov 2017 |
| `section6/perminov-cr209_fv56_cn359_ZT_reduced_3x4x6_r56_a359.json` | 56 | `"Moosbauer (symmetric flips)"` | **WRONG** — should credit Smirnov 2017 |
| `section6/perminov-ZT_3x4x6_r57_a405.json` | 57 | — | SUPERSEDED-OK (Smirnov gives 56, this row is rank-57 from elsewhere; out of scope) |

## Shape `R⟨3,5,5⟩` — Smirnov 2017 claims rank 58

| file | rank | current `attribution_for_rank` | verdict |
| --- | --: | --- | --- |
| `section5/alphatensor-Z_3x5x5_r58_a369.json` | 58 (over ℤ) | `"Sedoglavic-Smirnov 2021"` | **WRONG** — Smirnov 2017 is 4 years earlier than 2021. Should credit Smirnov 2017. |
| `section5/alphatensor-F2_3x5x5_r58_a413.json` | 58 (over F₂) | `"Sedoglavic-Smirnov 2021"` | **WRONG** — Smirnov 2017 is 4 years earlier. Should credit Smirnov 2017 (the ℤ scheme reduces mod 2). |
| `section5/perminov-c351_ZT_3x5x5_r58_a351.json` | 58 | `"AlphaTensor 2022"` | **WRONG** — Smirnov 2017 predates AlphaTensor 2022 by 5 years. Should credit Smirnov 2017. |
| `section5/perminov-cr221_cn357_ZT_reduced_3x5x5_r58_a357.json` | 58 | `"AlphaTensor 2022"` | **WRONG** — same as row above. |

Note: it is *possible* that the actual ⟨3,5,5⟩=58 scheme as encoded by
Sedoglavic–Smirnov 2021 / AlphaTensor 2022 is structurally distinct
from the one in Smirnov 2017 (different factor matrices realising the
same rank). The attribution policy still credits the **earliest known
publication of the bound**, which is Smirnov 2017. A separate scheme-
equivalence check (S₃-orbit / flip-graph distance) could disambiguate
whether the on-disk tensor IS the Smirnov 2017 one, but that does not
change the attribution for the *rank*.

## Shape `R⟨3,4,7⟩` — Smirnov 2017 claims rank 66

| file | rank | current `attribution_for_rank` | verdict |
| --- | --: | --- | --- |
| `section7/alphaevolve_3x4x7_r63_a588_0.5xC.json` | 63 (over ℂ) | (AlphaEvolve `source`) | SUPERSEDED-OK (63 < 66) |
| `section7/perminov-ZT_3x4x7_r64_a454.json` | 64 | `"meta-flip-graph search"` | SUPERSEDED-OK (64 < 66) |
| `section7/perminov-cr249_cn446_ZT_reduced_3x4x7_r64_a446.json` | 64 | `"meta-flip-graph search"` | SUPERSEDED-OK |

**Observation**: no on-disk scheme is at rank 66 — every encoded
⟨3,4,7⟩ scheme is strictly better. So Smirnov 2017's ⟨3,4,7⟩=66
contribution is **not currently represented in the JSON catalog** (it
has been superseded outright in our schemes). Cited-bound only —
candidate for `docs/cited-bounds.json`.

## Shape `R⟨3,4,8⟩` — Smirnov 2017 claims rank 75

| file | rank | current `attribution_for_rank` | verdict |
| --- | --: | --- | --- |
| `section8/perminov-ZT_3x4x8_r73_a976.json` | 73 | `"meta-flip-graph search"` | SUPERSEDED-OK (73 < 75) |
| `section8/perminov-Z_3x4x8_r73_a1902.json` | 73 | `"meta-flip-graph search"` | SUPERSEDED-OK |
| `section8/alphaevolve_3x4x8_r74_a461_Z.json` | 74 | (AlphaEvolve `source`) | SUPERSEDED-OK (74 < 75) |
| `section8/perminov-ZT_3x4x8_r74_a461.json` | 74 | `"AlphaEvolve 2025"` | SUPERSEDED-OK |
| `section8/perminov-cr267_cn461_ZT_reduced_3x4x8_r74_a461.json` | 74 | `"AlphaEvolve 2025"` | SUPERSEDED-OK |

**Observation**: no on-disk scheme is at rank 75. Smirnov 2017's
⟨3,4,8⟩=75 contribution is fully superseded in the JSON catalog.
Cited-bound only.

## Shape `R⟨3,5,7⟩` — Smirnov 2017 claims rank 82

| file | rank | current `attribution_for_rank` | verdict |
| --- | --: | --- | --- |
| `section7/perminov-ZT_3x5x7_r79_a520.json` | 79 | `"meta-flip-graph search"` | SUPERSEDED-OK (79 < 82) |
| `section7/perminov-Z_3x5x7_r79_a1927.json` | 79 | `"meta-flip-graph search"` | SUPERSEDED-OK |
| `section7/alphaevolve_3x5x7_r80_a652_Z.json` | 80 | (AlphaEvolve `source`) | SUPERSEDED-OK |
| `section7/perminov-ZT_3x5x7_r81_a712.json` | 81 | `"Perminov (FastMatrixMultiplication)"` | SUPERSEDED-OK |
| `section7/perminov-ZT_3x5x7_r82_a811.json` | 82 | `"Perminov (FastMatrixMultiplication)"` | **WRONG** — should credit Smirnov 2017 |
| `section7/perminov-ZT_3x5x7_r83_a499.json` | 83 | `"Perminov (FastMatrixMultiplication)"` | out of scope (rank 83 > 82) |
| `section7/perminov-cr303_cn494_ZT_reduced_3x5x7_r83_a494.json` | 83 | `"Perminov (FastMatrixMultiplication)"` | out of scope |

## Summary

| shape | rank | on-disk schemes at the Smirnov-2017 rank | WRONG | MISSING | SUPERSEDED-OK | OK |
| --- | --: | --: | --: | --: | --: | --: |
| ⟨3,4,6⟩ | 56 | 2 | 2 | 0 | 3 | 0 |
| ⟨3,5,5⟩ | 58 | 4 | 4 | 0 | 0 | 0 |
| ⟨3,4,7⟩ | 66 | 0 | 0 | 0 | 3 | 0 |
| ⟨3,4,8⟩ | 75 | 0 | 0 | 0 | 5 | 0 |
| ⟨3,5,7⟩ | 82 | 1 | 1 | 0 | 4 | 0 |
| **total** | | **7** | **7** | **0** | **15** | **0** |

**7 of 7** scheme files at a Smirnov-2017 rank currently credit a
later source — they should be re-tagged in a follow-up PR. None of
the files are missing an `attribution_for_rank` field outright; they
all have one, just pointing to the wrong (later) source.

## Action items (for a follow-up PR, NOT this one)

1. Re-tag the 7 `WRONG` rows above to `"attribution_for_rank":
   "Smirnov 2017"` (referencing [REFERENCES.md §[70]](../REFERENCES.md#70-smirnov2017)).
2. Add ⟨3,4,7⟩=66, ⟨3,4,8⟩=75 (and ⟨3,5,7⟩=82 if no on-disk scheme
   captures it cleanly) to `docs/cited-bounds.json` via
   `GenerateCitedBounds.java` so the historical Smirnov-2017 bound is
   recorded even where modern schemes have superseded it.
3. Once a local PDF is archived to
   `references/papers/smirnov_2017_3pq_several_bilinear.pdf`, cross-
   check the actual factor matrices to confirm the on-disk Z-tensors
   for ⟨3,5,5⟩=58 and ⟨3,5,7⟩=82 ARE (or ARE NOT) the Smirnov-2017
   tensors. If they match, also add a `tensor_identity_check_2017`
   metadata tag.

## Caveat

This audit relies on the user-supplied abstract of Smirnov 2017
(ranks 56 / 58 / 66 / 75 / 82 for the five shapes). The local PDF
could not be auto-downloaded (ResearchGate and Academia.edu both
return HTTP 403 to unauthenticated curl/WebFetch). Action item #3
above is the safety net: once the PDF is in-tree, the factor matrices
should be re-checked. Until then, the rank attributions in this audit
are based on the published claim, not on a side-by-side scheme match.
