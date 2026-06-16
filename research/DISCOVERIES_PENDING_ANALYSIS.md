# Discoveries pending analysis

Scratch file for experimental results that should NOT be lost if a
process is killed mid-investigation. Once a finding here is digested
into the catalog (cited-bounds.json / derived-bounds.json / a
materialised scheme JSON), move it from here to its permanent home.

---

## 2026-06-12 — Kaporin 2024 C⟨4,4,4⟩=48 predates AlphaEvolve 2025 (attribution audit)

Imported Kaporin's explicit complex `(4,4,4;48)` scheme from *Semi-analytical
solution of Brent equations*, Doklady Mathematics 518(1):29–34 (**2024**)
— REFERENCES.md [87], scheme
`schemes/known/section4/4x4x4-r48-kaporin_2024-8a17320.json`. Extracted from the
author's `test444r48.for` and **re-verified exactly** (`Verifier.isExactComplex`,
residual ~4.7e-15). Coefficients are complex floats (numerical, not rationals);
the scheme **differs from AlphaEvolve's C=48 original**.

**Pending analysis (discovery="TBD"):**
1. **Chronology / attribution.** This is an explicit C⟨4,4,4⟩=48 published a year
   *before* AlphaEvolve 2025 — so AlphaEvolve is NOT the originator of 48 over C.
   The canonical field-discipline example ("48/C (AlphaEvolve 2025)") in
   `CLAUDE.md` and the SPA legend should be revisited to credit Kaporin 2024 (and
   the Li–Zhang–Ke 2023 existence conjecture, arXiv:2310.11686). **Not yet edited
   CLAUDE.md** — load-bearing doc, surface to user first.
2. **Rationalisation.** Does Kaporin's numerical 48 rationalise to Q/R like the
   Dumas–Pernet–Sedoglavic 2025 rationalisation of AlphaEvolve's 48? Open.
3. **C⟨2,4,5⟩=32** (also in the paper, below-HK) is **not yet loaded** — the .for
   only carries the (4,4,4;48). Would need the (2,4,5;32) coefficients from the
   author.

---

## 2026-06-08 — FMM ⟨27,28,28⟩=10413 is NOT a stronger projection (kin-row gap is global)

Investigated why the FMM cross-check lists `⟨27,28,28⟩` as WORSE
(ours 11718 vs FMM 10413). Findings (all from
`ProbeStructuredProjection.java`, `AnalyzeFmmKinRows.java`, and the
downloaded `references/fmm-lille/27x28x28/27x28x28_raw.mpl`):

1. **Our coordinate projection already equals FMM's published scheme.**
   FMM's `27x28x28_raw.mpl` has **10442 products** — exactly our
   exhaustive coordinate-drop projection of the held ⟨28,28,28⟩=10556
   PanTA cube (μ=114 → 10442, verified). FMM's `[[1,0],[0]]` projection
   is **not** stronger than ours.
2. **The headline 10413 is `10442 − 29` "kin-row unification"**, which is
   NOT realised in FMM's own published raw scheme. We attempted to
   reproduce it and **falsified every local mechanism**:
   - stronger/linear projection (quotient an aggregation direction
     `e_a±e_d`): kills 0 extra products;
   - pairwise rank-1 merge (two products sharing 2 of 3 factor
     directions u/v/w): **0** on every axis-pair;
   - per-left-factor slice re-decomposition: the 27 size-3 shared-left
     groups (27 = the reduced dim) all have **full-rank** slices → **0**
     savings.
   ⟹ the 29-product gap is a **global tensor-rank reduction** (or the
   headline comes from a different construction the `_raw` export does
   not realise). Not closeable by any cheap local op.

**Actionable takeaway (NOT a missing projection capability):** most of
the 2083 "WORSE" FMM gaps are simply that the **downward projection
closure has not been run** over the 17–32 band — our operator reaches
FMM's *real* published ranks (e.g. ⟨27,28,28⟩ 11718 → 10442) once the
closure runs. Tracked by the `ProjectFmmGaps` projection sweep. The last
~0.3% (kin-row unification) is parked as open research; do not treat it
as a projection deficit.

Repro: `mvn -q exec:java -Dexec.mainClass=eu.solven.matmul.docs.explore.AnalyzeFmmKinRows`
and `…ProbeStructuredProjection -Dexec.args="28"`.

---

## 2026-06-04 — SchemeSweep band 11-16 (evaluate, stagnation cap) — 5 composition-beats-catalog candidates

First completed sweep of the band `11 ≤ maxDim ≤ 16` (515 shapes, R,
config=simple) after wiring balance-first ordering + a stagnation cap
into `AllocationOptimizer` / `BlockSplitSearch` (so large-base targets
terminate instead of timing out). The composition-vs-catalog section
reports **better: 5 / tie: 300 / worse: 210**. The 5 "better" rows —
single-base, non-recursive recombination strictly below the on-disk
catalog SOTA (NC, R):

| shape | composition | catalog | Δ |
| --- | ---: | ---: | ---: |
| ⟨2,10,15⟩ | 233 | 234 | −1 |
| ⟨2,10,16⟩ | 248 | 249 | −1 |
| ⟨2,12,16⟩ | 296 | 298 | −2 |
| ⟨3,3,14⟩ | 95 | 98 | −3 |
| ⟨3,3,15⟩ | 103 | 105 | −2 |

Repro:
```
mvn -q -o exec:java -Dexec.mainClass=eu.solven.matmul.docs.SchemeSweep \
  -Dexec.args="--mode=evaluate --field=R --band=11-16 --config=simple \
  --stagnation=1000000 --out=target/sweep-band11-16.md"
```
(see the "## Composition vs catalog" section of the output md.)

Status: **[3 false positives + 2 REAL WINS] (verified 2026-06-04).**

UPDATE 2026-06-04 (after explicit build+Verifier, `VerifyConcat3314`):
the 5 split into two groups, NOT all false:

- **⟨2,10,15⟩ / ⟨2,10,16⟩ / ⟨2,12,16⟩ — FALSE POSITIVES.** Came from the
  HK1971 formula PREDICTOR (unverified). Fixed: `ConstructiveMethod.Prediction`
  now carries a `verified` flag; `HopcroftKerr1971Method` sets it false;
  `findBestStrategy` no longer elects unverified predictions
  (`TestUnverifiedBoundsIgnored`). After the fix these shapes elect a real
  `concat-p` at 234/249/298 = catalog. No win.
- **⟨3,3,14⟩=95 and ⟨3,3,15⟩=103 — REAL & VERIFIED, but RE-DISCOVERIES of
  FMM-Lille (NOT novel).** FMM-Lille catalog AND Sedoglavic FMM-SOTA both
  already list ⟨3,3,14⟩=95 and ⟨3,3,15⟩=103 (see `references/fmm-lille-catalog.json`
  entries + `references/sedoglavic-fmm-sota.json`). Our on-disk best was the
  worse Perminov import (98 / 105) — a **catalog-completeness gap**, not a
  discovery. attribution_for_rank = FMM-Lille; our construction is the obvious
  column-block concat (below). ⚠ The FMM cross-check (task #180) should have
  flagged these FMM-better ranks and did not — investigate (subset of shapes? /
  ran before FMM listed them?).
  Plain column-block concat of on-disk NC atoms:
  - ⟨3,3,14⟩ = ⟨3,3,2⟩(=⟨2,3,3⟩=15, alphatensor_Z) ⊞ ⟨3,3,12⟩(=80, fmm_lille)
    = **95** < catalog 98 (Perminov). `isExactNonCubic`=true.
  - ⟨3,3,15⟩ = ⟨3,3,3⟩(=23, Laderman) ⊞ ⟨3,3,12⟩(=80) = **103** < catalog 105.
  The recursive materialiser MISSED these because it could not orient
  ⟨2,3,3⟩→⟨3,3,2⟩ and fell back to naive ⟨3,3,2⟩=18 (18+80=98). This implies
  a SYSTEMIC materialiser orientation gap — likely more such concat wins exist
  catalog-wide.
  TODO: (1) check FMM/Perminov for ⟨3,3,14⟩/⟨3,3,15⟩ < 95/103 to set
  attribution (discovery vs re-discovery); (2) materialise + register the two
  concat schemes (lineage = concat ⊞ of the named atoms); (3) fix the
  materialiser/atom-lookup orientation so concat sub-shapes orient through the
  6-fold matmul symmetry.

--- original (pre-build) note below, kept for history ---

Verified via `VerifyBand11to16Candidates` (materialise each with the
recursive materialiser, `balancedOnly=false`, no disk write, then
`Verifier.isExactNonCubic`):

| shape | evalComp | best strategy | catalog | materialised+verified | fromDisk |
| --- | ---: | --- | ---: | ---: | :--: |
| ⟨2,10,15⟩ | 233 | `HK1971` | 234 | **234 ✓** | yes |
| ⟨2,10,16⟩ | 248 | `HK1971` | 249 | **249 ✓** | yes |
| ⟨2,12,16⟩ | 296 | `HK1971` | 298 | **298 ✓** | yes |
| ⟨3,3,14⟩  | 95  | `concat-p[2,12]` | 98 | **98 ✓** | yes |
| ⟨3,3,15⟩  | 103 | `concat-p[3,12]` | 105 | **105 ✓** | yes |

Root cause — **predictor-tier ranks that are not buildable**:
- The ⟨2,·,·⟩ trio: best strategy is the **HK1971 formula PREDICTOR**
  (upper-rank tier), which claims 233/248/296 but yields no factor
  matrices; the materialiser can't realise it, so the verified scheme is
  the **Perminov-imported** ⟨2,b,c⟩ at 234/249/298 (task #111).
- ⟨3,3,14⟩/⟨3,3,15⟩: `concat-p` splits off ⟨3,3,12⟩ (real,
  `fmm_lille-3x3x12_m80`, buildable) **+ ⟨3,3,2⟩ (= ⟨2,3,3⟩), which has
  NO on-disk scheme** — its resolver rank is again an HK1971 prediction.
  So the concat prediction 95/103 isn't materialisable → falls back to
  the real ⟨3,3,14⟩=98 / ⟨3,3,15⟩=105.

Consequence (real bug, follow-up): **`SchemeSweep.catalogComparison`'s
"better" count is unsound when the SOTA resolver contains predictor-tier
ranks** — it compares a *predicted* composition cost against a *buildable*
catalog scheme. Fix: count "better" only when the composition is
buildable (sub-shape ranks come from constructors / on-disk schemes, or
the result passes `Verifier`). Until then, treat evaluate "better" rows
as candidates requiring materialisation, exactly as done here.

Open sub-question worth a look: HK1971 *predicts* ⟨2,10,15⟩=233 — one
below the Perminov import 234. Hopcroft–Kerr 1971 is constructive, so
either the predictor formula is optimistic for this format, or there is a
genuine HK construction at 233 we have not materialised. Settle by making
the HK constructor emit factor matrices for ⟨2,10,15⟩ and `Verifier`-ing
the result; only then is it a real improvement.

---

## 2026-06-01 — Schwartz-Zwecher 2025 import

Imported 7 explicit cubic schemes from the supplemental .npz of
arXiv:2508.01748 ("Towards Faster Feasible Matrix Multiplication by
Trilinear Aggregation", Schwartz & Zwecher, Hebrew University):

- `Q⟨20,20,20⟩:r=4378` — `section20/schwartz-zwecher-2025_20x20x20_r4378_a521506_Q.json` (6.1 MB)
- `Q⟨22,22,22⟩:r=5596` — `section22/...` (8.6 MB)
- `Q⟨24,24,24⟩:r=7020` — `section24/...` (12 MB)
- `Q⟨26,26,26⟩:r=8666` — `section26/...` (16 MB)
- `Q⟨28,28,28⟩:r=10550` — `section28/...` (21 MB) **discovery**
- `Q⟨30,30,30⟩:r=12688` — `section30/...` (27 MB) **discovery**
- `Q⟨32,32,32⟩:r=15096` — `section32/...` (35 MB) **discovery**

Status of pending analysis:

- **n=20..26 are NOT new discoveries.** At those n, DIS09's Pan TA
  bound is slightly better (n=20: 4340 vs SZ 4378; etc.). The .npz
  files for those n are provided in the supplemental zip for
  completeness — the TA-New25 construction produces the same scheme
  shape for every even n≠16, but only beats Pan starting at n=28.
  Their `attribution_for_rank` is set to "Pan 1982" and
  `discovery: false`.
- **n=28, 30, 32 ARE new discoveries** matching the paper's Table 1.
- **n=34..50 deferred.** See ROADMAP entry "Schwartz-Zwecher 2025 —
  n=34..50 import deferred" — sizes grow rapidly (n=44 alone would be
  hundreds of MB), and per project policy we prefer to derive these
  ourselves via the kin-row unification constructor (TBD).
- **Derivation follow-up tracker.** Every imported JSON carries
  `"derivation_task": "TBD-SZ2025-kin-row-constructor"`. Implementing
  the constructor (paper §3, Theorem 2.22 + 3.4) is filed in ROADMAP
  under "Schwartz-Zwecher 2025 — kin-row unification constructor".

Verification: all 7 schemes pass `Verifier.passesRandomMatmulSpotCheck`
(5 random A,B sample residual < 1e-9). `SymbolicVerifier.verifyBilinear`
is intractable at this scale (`O(n⁶·r)` BigInteger multiplies), so we
rely on the fast random-input verifier. The float→Q rationalisation
in the import script
(`tools/sz2025-import/import_schwartz_zwecher_2025.py`) cross-checks
against the algorithm by re-running the spot check on the snapped
matrices — all 7 schemes pass after snapping.

Supplemental data: `references/schwartz-zwecher-2025/trilinear_aggregation_algorithms_decomposed/`
(retained on-disk; .npz files plus README and sample.py).

Cited-bounds: 13 SZ entries added to `docs/cited-bounds.json` via
`GenerateCitedBounds.java` (covers n=28..50 + n=60).

---

## 2026-05-28 — DIS09 ⟨n,n,n⟩ comparison — PER-FIELD CORRECTED

**Important field-discipline note**: DIS09 Table 3 is the
**non-commutative** upper bound, applicable to any non-commutative
ring (equivalently: it must hold over a generic ring, so it's the
strict upper bound for R/Q/Z). DIS09 Table 4 is the **commutative**
version — different bounds, NOT compared in this report.

For each target field we apply the same SOTA resolver pipeline
(catalog ∪ Pan TA) but filter the catalog to the schemes valid for
that field:

- **R**: integer / rational / real schemes (Z, Q, R, ZT-tagged).
  Excludes F₂-only (AlphaTensor, ATf2-composed), complex-only
  (AlphaEvolve 0.5xC, AE2 composed), and commutative-only
  (Waksman, Rosowski, Makarov86).
- **C**: R-friendly ∪ complex-only (AlphaEvolve `_0.5xC`).
- **F2**: F₂-native (AlphaTensor `_F2_`, ATf2-composed).
  Note: integer schemes also reduce to F₂ mod 2, but this is NOT yet
  included in the F₂ filter — see "Known gaps" below.

**Repro commands**:

```bash
mvn -q -o test -Dtest=TestDIS09FullScan#dis09_full_comparison_R_4_to_30_unbalanced
mvn -q -o test -Dtest=TestDIS09FullScan#dis09_full_comparison_C_4_to_30_unbalanced
mvn -q -o test -Dtest=TestDIS09FullScan#dis09_full_comparison_F2_4_to_30_unbalanced
```

**Pool** (S₃-orbit-expanded, 21 entries):
Strassen `⟨2,2,2⟩=7` ×6, Laderman `⟨3,3,3⟩=23` ×6, mul211/mul121/mul112 ×3 each.

### Field R (R/Q/Z, non-commutative) vs DIS09 Table 3 — 15 wins, 12 ties, 0 losses

| n  | DIS09 | ours | source                            |
|----|-------|------|-----------------------------------|
| 4  | 49    | 49   | = Strassen[2,2]³                  |
| 5  | 100   | 93   | AlphaEvolve Z (Perminov)          |
| 6  | 161   | 153  | Perminov ZT direct                |
| 7  | 258   | 250  | Sedoglavic Strassen[3,4]³         |
| 8  | 343   | 343  | = Strassen[4,4]³ (TIE)            |
| 9  | 522   | 486  | Perminov ZT direct                |
| 10 | 700   | 651  | Strassen[5,5]³                    |
| 11 | 923   | 873  | Strassen[5,6]³                    |
| 12 | 1125  | 1068 | Perminov ZT direct                |
| 13 | 1450  | 1421 | Perminov Q direct                 |
| 14 | 1728  | **1720** | Strassen[7,7]³ + paired ⭐    |
| 15 | 2108  | 2058 | Perminov ZT direct                |
| 16 | 2401  | 2304 | fmm-lille Q direct                |
| 17 | 2972  | 2955 | Strassen[8,9]³ + paired ⭐        |
| 18 | 3306  | 3306 | = Pan TA bound (TIE)              |
| 19 | 4073  | **4030** | Strassen[9,10]³ + paired ⭐    |
| 20 | 4340  | 4340 | = Pan TA bound (TIE)              |
| 21 | 5365  | **5258** | Strassen[9,12]/[9,12]/[10,11] non-balanced ⭐⭐⭐ |
| 22 | 5566  | 5566 | = Pan TA bound                    |
| 23 | 6806  | **6731** | Strassen[11,12]³ + paired ⭐  |
| 24 | 7000  | 7000 | = Pan TA bound                    |
| 25 | 8448  | 8448 | = Pan TA bound                    |
| 26 | 8658  | 8658 | = Pan TA bound                    |
| 27 | 10330 | 10330| = Pan TA bound                    |
| 28 | 10556 | 10556| = Pan TA bound                    |
| 29 | 12468 | 12468| = Pan TA bound                    |
| 30 | 12710 | 12710| = Pan TA bound                    |

Aggregate: 113,825 vs DIS09 114,466 (−0.56%). All ⭐ rows are
**genuine new bounds** attributable to this repo (n=14, 19, 21, 23).

### Field C (complex) — 17 wins, 10 ties, 0 losses

| n | DIS09 | ours | source                          |
|---|-------|------|---------------------------------|
| 4 | 49    | **48**  | AlphaEvolve 0.5×C            |
| 5 | 100   | 93   | (same as R)                     |
| 7 | 258   | **249** | Strassen[3,4]³ + paired      |
| 8 | 343   | **336** | Strassen[4,4]³ + paired      |
| (others same as R)               |
Aggregate: 113,816 (−0.57% vs DIS09).

**Note**: Most wins over R also apply over C since R ⊂ C. The
additional C-specific wins (n=4: 48, n=7: 249, n=8: 336) come from
AlphaEvolve's `⟨4,4,4⟩=48` lifting into the recursion.

### Field F2 — 10 wins, 12 ties, 5 losses (current state)

| n  | DIS09 | ours | source                            |
|----|-------|------|-----------------------------------|
| 4  | 49    | **47**  | AlphaTensor F₂ direct          |
| 5  | 100   | 96   | AlphaTensor F₂                    |
| 7  | 258   | 248  | Strassen[3,4]³ + paired           |
| 8  | 343   | 329  | Strassen[4,4]³                    |
| 9  | 522   | 513  | Strassen[4,5]³                    |
| 10 | 700   | 672  | Strassen[5,5]³                    |
| 11 | 923   | 895  | Laderman[4,3,4]/[3,4,4]/[4,4,3]   |
| 12 | 1125  | 1081 | Laderman[4,4,4]³                  |
| 13 | 1450  | **1487** (+37) | LOSS — F₂ catalog gap   |
| 14 | 1728  | **1750** (+22) | LOSS — F₂ catalog gap   |
| 15 | 2108  | **2208** (+100) | LOSS — F₂ catalog gap  |
| 16 | 2401  | 2209 | AlphaTensor F₂² (47²)             |
| 19 | 4073  | **4178** (+105) | LOSS — F₂ catalog gap  |
| 21 | 5365  | **5388** (+23) | LOSS — F₂ catalog gap   |
| (Pan TA matches DIS09 from n=18 onward except above)            |

Aggregate: 114,408 (−0.05% vs DIS09 — essentially break-even).

**The F₂ losses are a filter issue**, not a real algorithmic gap.
Integer-coefficient schemes (Z, ZT) reduce to F₂ mod 2 and would
re-introduce the wins seen in R. The current F2 filter only matches
F₂-tagged files. Fix is a one-line filter change — see "Known gaps".

### Known gaps / things to fix

1. **F2 filter missing integer schemes**. `loadCatalogBestRanksForField("F2")`
   currently includes only files tagged `f2`/`z2`. Any Z/ZT/Q-integer
   scheme reduces mod 2 and should also count. Concrete fix:
   ```java
   case "F2" -> loadCatalogBestRanksFiltered(name ->
       isF2(name) || isRfriendly(name));
   ```
   But: Q schemes with denominators (½, ⅓) DON'T reduce. Need to read
   JSON `.field` metadata and check whether coefficients are integral.
2. **JSON metadata not consulted**. Filter relies only on filename
   conventions. 2,318 of 2,381 schemes have no `field` key in JSON.
   Would need a metadata audit pass to canonicalise.
3. **Commutative comparison absent**. DIS09 Table 4 (commutative
   ⟨n,n,n⟩) was not reproduced. Catalog has Waksman, Rosowski 2019,
   and Makarov 1986 commutative schemes — would need a separate scan
   pipeline.

---

## 2026-05-28 — DIS09 Table 4 (COMMUTATIVE) comparison

**13 wins, 9 ties, 5 losses** vs DIS09 Table 4. Total ours = 109,872
vs DIS09 = 110,213 (**−0.31%**).

**Commutative SOTA**: min of
- NC catalog (R-filtered, any NC scheme is also a commutative upper bound)
- Waksman 1970 family (`WaksmanBound.forShape` over all 3 axis orientations)
- Rosowski 2019 commutative bilinear (`RosowskiBound.bestCommutativeBound`,
  all 6 axis perms; includes the `⟨3,3,3⟩=21` formula — non-bilinear
  scheme is in `references/rosowski-algorithms.md`)
- Pan TA closed form (NC ≥ commutative as upper bound)

**Repro**:
```bash
mvn -q -o test -Dtest=TestDIS09Table4Commutative
```

### Highlights

| n  | DIS09 Cmt | ours | source                                        |
|----|-----------|------|-----------------------------------------------|
| 5  | 93        | **85**  | mul211 + [2,3]/[5]/[5] — Rosowski ⟨5,3,5⟩=51 + Waksman ⟨5,2,5⟩=34 ⭐ |
| 7  | 235       | **217** | mul211 + [2,5]/[7]/[7]                       |
| 9  | 472       | **441** | mul211 + [2,7]/[9]/[9]                       |
| 11 | 825       | **781** | mul211 + [2,9]/[11]/[11]                     |
| 13 | 1318      | **1261**| Strassen[6,7]³                               |
| 14 | 1525      | **1519**| Strassen[7,7]³                               |
| 15 | 1941      | **1891**| Strassen[7,8]³                               |
| 17 | 2762      | **2673**| Strassen[8,9]³                               |
| 19 | 3757      | **3673**| Strassen[9,10]³                              |
| 21 | 4938      | **4861**| Strassen[10,11]³                             |
| 23 | 6382      | **6315**| Strassen[11,12]³                             |
| 25 | 8083      | **7993**| Strassen[12,13]³                             |
| 27 | 9994      | **9985**| Strassen[13,14]³                             |

### Losses (where DIS09 used Hopcroft332 or mul121 + commutative tricks we don't replicate)

| n  | DIS09 | ours | gap |
|----|-------|------|-----|
| 18 | 3060 (Hopcroft332) | 3087 | +27  |
| 20 | 4158 (Strassen)    | 4165 | +7   |
| 22 | 5440 (mul121)      | 5467 | +27  |
| 24 | 6900 (Hopcroft332) | 7000 | +100 |
| 29 | 12109 (mul121)     | 12237 | +128 |

These losses come from DIS09 using **Hopcroft-Kerr ⟨3,3,2⟩=15** as a
recursive base (which we don't have in the multi-base pool) and from
more aggressive mul121 recursion. The Hopcroft332 base would close
n=18 and n=24.

### Wins note

The dominant pattern is **mul211 + axis-split + Rosowski formula** —
DIS09 didn't have Rosowski 2019 (paper appeared 10 years after DIS09),
so any sub-product shape where Rosowski's `n(lm+l+m−1)/2` beats
Waksman gives us a free win. Specifically the `(b odd)` Rosowski
branch is tighter than Waksman's `(b−1)/2 + ac` for many shapes.

---

### Genuinely new bounds vs DIS09 Table 3 (all R/Q/Z)

These come from this repo on top of DIS09's framework:

- **n=14 R(⟨14,14,14⟩) ≤ 1720** via Strassen[7,7]³ + Pan paired sub-products
- **n=19 R(⟨19,19,19⟩) ≤ 4030** via Strassen[9,10]³ + paired
- **n=21 R(⟨21,21,21⟩) ≤ 5258** via non-balanced Strassen[9,12]/[9,12]/[10,11]
- **n=23 R(⟨23,23,23⟩) ≤ 6731** via Strassen[11,12]³ + paired

Emitted to `docs/derived-bounds.json` tagged `solven-strassen [30]`.

---

## 2026-05-28 — AT bulk-import provenance audit (TBD)

**Problem**: every `alphatensor-F2_*.json` and `alphatensor-Z_*.json`
scheme was imported with `source: "AlphaTensor"`, but the Fawzi 2022
Nature paper only highlights a subset as actual discoveries. The
rest match known bounds (Smirnov 2013, Strassen-recursion derivatives,
Pan, Hopcroft-Kerr) that AlphaTensor reproduced but did not discover.

**Impact**: catalog rows tagged "AlphaTensor 2022" for non-discovery
ranks misattribute the historical record. Concrete example: AT-Z
`⟨2,3,3⟩=15` is genuinely Hopcroft-Kerr 1971 (proven tight that year);
AT found the same rank computationally — not an improvement.

**Fix required**:

1. Read AlphaTensor Supplementary Table 2 — colour codes mark
   improvements vs matches.
2. For each `alphatensor-*.json`:
   - If improvement → add `"discovery": true`.
   - Else → add `"discovery": false, "attribution_for_rank": "..."`
     pointing to the actual originator (Smirnov 2013, Strassen 1969 +
     Sedoglavic 2017 recursion, Hopcroft-Kerr 1971, Pan, etc.).
   - If unsure → `"discovery": "TBD"` with a `"notes"` field
     describing what to investigate.
3. Update `GenerateCitedBounds.attributionLabel` and Pages display
   logic to prefer `attribution_for_rank` when `discovery: false`.
4. Apply the same audit to other bulk imports (`perminov-*` from
   FastMatrixMultiplication may have similar issues — that catalog
   aggregates schemes from many original authors).

**Out of scope for now** — flagged in CLAUDE.md as a project
convention. Track as a follow-up PR.

---

## How to use this file

- **Add to top** of the relevant section as new findings appear.
- Each entry should be self-contained: claim, repro command,
  breakdown, status.
- When a finding is digested into the permanent catalog
  (`docs/cited-bounds.json` / `docs/derived-bounds.json` / a
  materialised scheme JSON), DELETE the entry from here with a brief
  note in the commit message.
- If a finding becomes superseded by a later result, mark it
  `[SUPERSEDED by …]` rather than deleting — keeps history.

## Serendipitous sweep (2026-06-05, #159 Phase 1)

**Run:** `SerendipitousSweep --baseCap=8 --secondCap=6 --targetCap=18 --field=Q`
(205 bud-rich bases of 2439 catalog files).

**Result: no genuinely new wins in range.** The one flagged hit —
`⟨18,18,18⟩=3200 via fmm_lille-3x6x6_m80 ⊗ ⟨6,3,3⟩` (verified) — is a
**re-derivation** of a rank we already hold (`composed_kron_concat-18x18x18_m3200`,
`composed_recursive-18x18x18_m3200`). It was flagged only because
`FieldAwareLookup.findWithSource(18,18,18)` returns the dense `dis09_Q` **3306**
instead of the on-disk **3200** — the 3200 files are *stubs* (maxDim 18 > 16) and
the lookup skips them.

**Takeaways:**
- The serendipitous engine correctly re-derives a known catalog best (good validation).
- **Bug to fix:** the SOTA comparator under-reports for stub shapes → spurious "wins".
  Future sweeps should compare against the stub rank (parse `_m{R}` from stub
  filenames or expand), not just `findWithSource`.
- Real new wins likely need (a) wider caps and (b) the Phase-2 multi-scheme store
  (bud-rich suboptimal bases), since the rank-optimal catalog base often has fewer buds.

### Correction (2026-06-05): the ⟨18,18,18⟩=3200 flag was a FALSE POSITIVE

3200 is a **plain Kronecker** product `⟨3,3,6⟩=40 ⊗ ⟨6,6,3⟩=80 = 3200`
(FMM `18x18x18.html` states exactly this; both bases on disk). Serendipitous
adds **nothing** for it. The sweep flagged it because its keep-criterion is too
weak: it compares `r_s` only against (a) the per-base naive Kronecker — whose
`R(⟨6,3,3⟩)` came back > 40 from a lookup orientation-miss — and (b) the
stub-skipping lookup SOTA 3306. Neither is the true target best.

**Net: zero genuine serendipitous wins in the sweep.** Required fix before any
sweep result is trustworthy: compare `r_s` against the *true* target SOTA,
computed consistently — including plain Kronecker over **all** factorisations
and orientations, plus stub ranks (parse `_m{R}`). Only then does a reported
"win" mean `r_s` < everything else we can already build.

---

## 2026-06-06 — Bud-base factory (serendipitous product) — 6 verified wins over catalog

First run of `BudBaseFactory` (serendipitous product from bud-rich
bases, cost via `SerendipitousBudProduct.serendipitousCost`, then
build+`Verifier.passesRandomMatmulSpotCheck`). It predicted **29**
serendipitous results below the current catalog; the **8 smallest were
built and spot-checked**, of which **6 CONFIRMED** (valid scheme,
strictly below catalog) and 2 unconfirmed (see caveat).

| shape | serendipitous | catalog (was) | Δ | base ⊗ inner |
| --- | ---: | ---: | ---: | --- |
| ⟨6,9,21⟩ | **709** | 715 | −6 | ⟨2,3,7⟩ ⊗ ⟨3,3,3⟩ |
| ⟨6,12,24⟩ | **1029** | 1040 | −11 | ⟨2,4,8⟩ ⊗ ⟨3,3,3⟩ |
| ⟨6,16,24⟩ | **1383** | 1404 | −21 | ⟨2,4,8⟩ ⊗ ⟨3,4,3⟩ |
| ⟨6,20,15⟩ | **1136** | 1141 | −5 | ⟨2,5,5⟩ ⊗ ⟨3,4,3⟩ |
| ⟨4,18,24⟩ | **1119** | 1120 | −1 | ⟨2,6,6⟩ ⊗ ⟨2,3,4⟩ |
| ⟨4,21,28⟩ | **1518** | 1520 | −2 | ⟨2,7,7⟩ ⊗ ⟨2,3,4⟩ |

All wins are **rectangular, flat-axis bud-bases** (a "2" axis), as the
serendipitous theory predicts — the inner enlargement is sub-additive
there, unlike the cube. These are the first ranks the bud-base lever
found that the recombination/Kronecker pool missed.

Repro:
```
mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.research.BudBaseFactory
```

CAVEAT (optimality discipline — predictions are UPPER BOUNDS, verify
each): `serendipitousCost` sums `findRank` of the enlarged inners
(which may be a stub rank), while the actual build uses `findWithSource`
(materialised) — so a predicted cost can be optimistic vs the buildable
rank. A bisection (`BudConstructionDebug`) showed the bud-block
construction itself is sound for U/V/W types (all verify), so the 2
unconfirmed cases are cost-optimism / orientation edges, NOT a uniform
construction bug. The build+verify gate is therefore mandatory; the 6
above passed it.

Status: **[6 REAL WINS, verified 2026-06-06] — not yet materialised to
disk JSON.** Next: materialise + register the 6 as scheme files; verify
the remaining 21 predictions; debug the cost-vs-build gap.

---

## Projection gap: non-cubic neighbours of strong cubes (2026-06-07)

**Claim.** FMM beats us at ⟨29,30,30⟩ (FMM 12588, ours 13902). FMM's 12588
is a **projection [[1,0],[0]] of a ⟨30,30,30⟩=12710** scheme (per their
page). We hold the *same* best cube FMM does — `Schwartz:2025aa` =
`schwartz_zwecher_2025-30x30x30_m12688` (the digest lists FMM's own best 30³
as **12688**; the 12710 they cite is just the projection-parent for this
format, not their best cube). So we do NOT have a better cube — we have the
**same** one — but we never project it (or a 12710-class parent) down to the
non-cubic neighbour, so we carry a recombination 13902 instead of ≤~12588.

**Correction (2026-06-07):** an earlier draft of this note wrongly said "we
hold a *better* cube (12688 < 12710)". The FMM *digest*
(`references/fmm-lille-catalog.json`) reflects ⟨29,30,30⟩=12588 and
⟨30,30,30⟩=12688 correctly — nothing stale there; the gap is ours to close
by projecting, not a digest error.

**Why we don't find it — two root causes:**
1. **Scope.** The closure's default scope is *cubic-only*
   (`SchemeSweep` CLOSURE: `non-cubic off`). So non-cubic targets like
   ⟨29,30,30⟩ are never handed to `projectInto`, even though
   `PROJECTION_MAX_DELTA=1` and a delta-1 parent (⟨30,30,30⟩) exists.
2. **Heavy dense projection.** The SZ cubes store their matrices in an
   external `source_data_file` and densify to ~270 MB of `double[][]`
   (rank≈12.7k, dim 30). A single `--shape=29x30x30 --strategies=projection`
   probe ran > 6 min without completing (load + densify + project).
   Projection over the large-dim band is impractical with the current
   dense operator → ties to the **sparsity-in-NonCubicBilinearAlgorithm**
   ROADMAP item (a sparse/lazy projection would make this cheap).

**Repro:** `--mode=materialize --field=Q --shape=29x30x30
--strategies=projection --improve=true` (heavy; needs sparse projection
to be practical).

**Status: open gap, representative of many FMM-better non-cubic shapes
that are projections of cubes we already match/beat.** Fix path: (a) widen
projection to non-cubic shapes; (b) sparse/lazy projection so large-dim
parents are affordable. NOT a blanket "discard FMM" — 1888 formats remain
genuinely FMM-better (see `target/fmm-gap-report.md`).

## VerifyAllSchemes memory: throttle heavy dense schemes (2026-06-07)

**Symptom.** A full `VerifyAllSchemes` over the grown catalog (11328 files)
GC-thrashed near a 5 GB cap (5421/11328 in 38 min, then ~4 h projected) and
an earlier `-Xmx10g` run left a 10 GB OOM dump.

**Cause.** Large *dense* schemes (SZ n≥24 / FMM cubes) densify to hundreds
of MB of `double[][]`; all 15 worker threads loading one at once blows the
heap. (The earlier 10 GB OOM was pre-stub-fix — it densified the big STUBS
too; the compositional `LineageVerifier` stub path fixed that half.)

**Fix (applied).** Size-gated concurrency: schemes with
`weight = r·(nm+mp+np) > 10M` are throttled to `nThreads/4` concurrent
(`Semaphore`), small schemes (the majority) run unthrottled. Peak heap is
now bounded (~3 heavy × ~270 MB). CI `verify-catalog` heap bumped to 6 GB,
timeout to 120 min, with heap-dump-on-OOM + artifact upload.

**Status: fix committed to working tree; full local re-verify NOT re-run
(70+ min) — CI will exercise it. If still slow, next lever is descending-
weight scheduling or skipping unchanged raw imports (change-tracking).**
