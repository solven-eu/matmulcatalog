# Curated improvement examples (verified)

Source for the paper's *illustrative* examples (Section "Comparison: structure
and conventions", `ssec:no-static-ranks`) — a small, hand-vetted set of
**non-commutative**, **integer (ℤ)** cases where this work's own derivation beats
the external state of the art, **one per construction mechanism**. Not
exhaustive: the live, complete comparison is the SPA / `docs/comparison/*.json`.

## Selection + verification rules
- **Non-commutative only**, **field ℤ** (so it also reduces over F₂/F₃ and lifts
  to Q/R/C; no rational coefficients).
- Strictly beats the external best over R = `min(FMM-Lille, Perminov)`, where
  **Perminov now includes his serendipitous 17–32 catalog** (registered as 971
  cited-bounds in `docs/cited-bounds.json`, sourced from
  `references/perminov-serendipitous-17-32.json`).
- **Each re-verified now** by replaying its lineage stub and running
  `Verifier.isExactNonCubic` (`docs.verify.VerifyOneScheme`). Times are real.

> **Why this is the right discipline (a live example):** before registering
> Perminov's 17–32 catalog, ⟨14,20,22⟩=3580 read as a clean serendipitous *win*
> (FMM lists 3682). Registering Perminov revealed he **independently reaches the
> same 3580** — so it is a *tie*, not a beat. Simultaneously ⟨6,22,27⟩=2209
> turned from "beats FMM only" into "beats FMM **and** Perminov". Static rank
> tables would have mis-stated both. This is why the paper points to the live
> catalog instead of freezing a table.

## Verified — beats the current external best (one per mechanism)

| Mechanism | Shape | Our rank | FMM | Perminov | beats | verified |
| --- | --- | ---: | ---: | ---: | --- | --- |
| Coordinate projection | ⟨6,22,27⟩ | **2209** | 2236 | 2211 | **FMM + Perminov** | exact ✓ (29 s) |
| Recombination | ⟨9,10,17⟩ | **980** | 988 | — | FMM (Perminov n/a <18) | exact ✓ (2.3 s) |
| Concatenation | ⟨6,14,22⟩ | **1176** | 1185 | — | FMM (Perminov n/a here) | exact ✓ (4.0 s) |
| Serendipitous product | ⟨14,20,22⟩ | 3580 | 3682 | **3580** | FMM; **ties** Perminov | exact ✓ (132 s) |

- **Projection ⟨6,22,27⟩=2209 is the headline**: an integer scheme beating *both*
  FMM-Lille (2236) and Perminov's serendipitous result (2211). Lineage:
  `⟨6,22,28⟩ ↓ (peel one p-column)`.
- **Serendipitous ⟨14,20,22⟩=3580** is a *convergent* result: our `2x4x11 ⊗ˢ ⟨7,5,2⟩`
  and Perminov independently reach 3580, both beating FMM's 3682 — an honest
  "our method matches the SOTA-setter and beats FMM" example rather than a
  unique win.
- Lineages: recombination `R[alphatensor_Z ⟨2,3,3⟩=15 ^ABC→BAB; 3,3,3 | 5,5 | 6,5,6]`;
  concat `6x2x22 +m 6x12x22`.

Files (lineage stubs):
`derived/section27/6x22x27-r2209-derived-b966997.json`,
`derived/section17/9x10x17-r980-derived-3ec9766.json`,
`derived/section22/6x14x22-r1176-derived-76fe5e7.json`,
`derived/section22/14x20x22-r3580-derived-2760225.json`.

## Further Perminov improvements (candidates, verification pending)

Beyond ⟨6,22,27⟩, our work beats Perminov's serendipitous 17–32 ranks on more
shapes (per `references/perminov-serendipitous-17-32.json`), e.g.:

| Shape | Our rank | Perminov | FMM | gap vs Perminov | mechanism |
| --- | ---: | ---: | ---: | ---: | --- |
| ⟨16,26,28⟩ | 6488 | 6545 | 6552 | 57 | projection |
| ⟨10,22,25⟩ | 3316 | 3332 | 3343 | 16 | concat |
| ⟨20,20,26⟩ | 5793 | 5796 | 5838 | 3 | serendipitous |

(These are larger shapes — confirm with `VerifyOneScheme` before quoting.)

## Reaching SOTA by our own construction (a tie, not a beat)

| Mechanism | Shape | Our rank | FMM | Perminov | note |
| --- | --- | ---: | ---: | ---: | --- |
| HK71 constructive | ⟨2,10,15⟩ | 233 | 233 | 233 | exact ✓ (0.06 s) |

`hk71-constructive(10,15)` matches the SOTA bound of 233 rather than beating it,
but is the explicit verified scheme we hold (the Perminov entry we carry is 234).
File: `constructed/section15/2x10x15-r233-hk71-2a4c9ea.json`.

## Honest caveats (load-bearing for the paper)
- **Snapshots, not permanent.** As the ⟨14,20,22⟩ flip shows, FMM/Perminov move;
  re-sync (`SyncReferenceCatalogs`, register any new Perminov file) and re-verify
  before quoting any number.
- **No Kronecker example:** the only Kron "win" (⟨16,16,16⟩=2209 via 4×4×4=47 ⊗
  4×4×4=47) is **F₂-only** (47 is AlphaTensor's F₂ rank; over R the best 4×4×4 is
  48 → 2304), so it is excluded as not R/ℤ-valid.
- **Comparison digest gap (open):** `references/catalogs/perminov-catalog.json` (which the
  us-vs-catalogs comparison reads) still stops at maxdim 16. The serendipitous
  17–32 data is registered as cited-bounds and is in the catalog/SPA, but
  `GenerateSourceComparison` will not show Perminov for 17–32 until
  `SyncReferenceCatalogs` is extended to merge the serendipitous file into the
  digest.
