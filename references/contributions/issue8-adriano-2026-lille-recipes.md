# Issue #8 — Marcos Adriano (2026): 22 FMM-Lille "recipe" ranks as explicit schemes — re-derived here, 13 of them beaten

**Source.** GitHub issue
[solven-eu/matmulcatalog#8](https://github.com/solven-eu/matmulcatalog/issues/8)
(opened 2026-09-24 by [@marcosadriano27](https://github.com/marcosadriano27),
"22 explicit schemes at FMM-Lille table ranks (the table's own recipes,
executed)"), attachment `lille-recipes-22-schemes.zip`
(<https://github.com/user-attachments/files/32619676/lille-recipes-22-schemes.zip>,
5.2 MB zipped / **321 MB** unzipped: 9 sparse integer files, 13 dense files with
quoted rationals, 1.3–49 MB each). The contributor states the work was done
*with the help of AI tools*. Processed 2026-09-30 (this note is the record).

**What was contributed.** For 22 shapes where FMM-Lille lists a rank but no
public explicit scheme reaches it, the table's own `construction` recipe
(a recombination over unequal blocks — outer ⟨2,4,4⟩:26 Hopcroft–Kerr for 14
shapes, Strassen ⟨2,2,2⟩:7 for 8 — with every outer product computed at the
smallest size its two sides need) executed exactly as stated, as explicit
factor matrices. Not new ranks (FMM lists them); explicit witnesses where none
was public. Pieces from FMM-Lille files and, for ⟨7,12,16⟩:876 / ⟨7,15,16⟩:1081,
from Perminov's repo at `b28490c` (MIT).

## Outcome in one table

Our engine **re-derives all 22** as lineage stubs (no explicit matrices
committed — see "Why not import"), and lands **strictly below the recipe on 13
of the 14 ⟨2,4,4⟩-outer shapes** (below FMM's own index too), ties the
remaining 9. Ranks are ℚ-field claims (fields per stub below); optimality tier:
**bound** (recombination B&B under a node cap — `[node cap]` in the probe log
means the optimizer stopped at its budget, so *lower may exist*).

| shape | FMM recipe (issue) | FMM index (digest 2026-09-29) | ours before | **ours now** | Δ vs recipe | route (our stub's lineage) | stub fields |
|---|---:|---:|---:|---:|---:|---|---|
| ⟨6,14,25⟩ | 1322 | 1320 | 1324 | **1322** | 0 | R[⟨2,4,4⟩ AT-Z rep; 3,3 \| 4,4,3,3 \| 7,6,6,6] | F3,Q,R,C |
| ⟨10,19,31⟩ | 3532 | 3532 | 3535 | **3492** | **−40** | R[⟨2,4,4⟩ Perminov-ZT rep; 5,5 \| 5,5,5,4 \| 8,7,8,8] | F2,F3,Z,Q,R,C |
| ⟨10,22,25⟩ | 3308 | 3308 | 3316 | **3288** | **−20** | R[⟨2,4,4⟩ P-ZT; 5,5 \| 6,6,6,4 \| 7,6,6,6] | F2,F3,Z,Q,R,C |
| ⟨10,22,29⟩ | 3795 | 3795 | 3800 | **3772** | **−23** | R[⟨2,4,4⟩ P-ZT; 5,5 \| 6,6,6,4 \| 8,7,7,7] | F2,F3,Z,Q,R,C |
| ⟨10,23,23⟩ | 3183 | 3183 | 3186 | **3146** | **−37** | R[⟨2,4,4⟩ P-ZT; 5,5 \| 6,6,6,5 \| 6,5,6,6] | F2,F3,Z,Q,R,C |
| ⟨10,23,26⟩ | 3571 | 3571 | 3588 | **3534** | **−37** | R[⟨2,4,4⟩ P-ZT; 5,5 \| 6,6,6,5 \| 7,5,7,7] | F2,F3,Z,Q,R,C |
| ⟨10,23,27⟩ | 3688 | 3688 | 3691 | **3648** | **−40** | R[⟨2,4,4⟩ P-ZT; 5,5 \| 6,6,6,5 \| 7,6,7,7] | F2,F3,Z,Q,R,C |
| ⟨10,23,30⟩ | 4076 | 4076 | 4087 | **4036** | **−40** | R[⟨2,4,4⟩ P-ZT; 5,5 \| 6,6,6,5 \| 8,6,8,8] | F2,F3,Z,Q,R,C |
| ⟨10,23,31⟩ | 4193 | 4193 | 4196 | **4150** | **−43** | R[⟨2,4,4⟩ P-ZT; 5,5 \| 6,6,6,5 \| 8,7,8,8] | F2,F3,Z,Q,R,C |
| ⟨10,25,26⟩ | 3864 | 3864 | 3884 | **3850** | **−14** | R[⟨2,4,4⟩ P-ZT; 5,5 \| 6,6,7,6 \| 7,5,7,7] | F2,F3,Z,Q,R,C |
| ⟨10,26,26⟩ | 4028 | 4028 | 4060 | **4014** | **−14** | R[⟨2,4,4⟩ P-ZT; 5,5 \| 7,7,7,5 \| 7,7,6,6] | F2,F3,Z,Q,R,C |
| ⟨10,26,27⟩ | 4172 | 4172 | 4194 | **4144** | **−28** | R[⟨2,4,4⟩ P-ZT; 5,5 \| 7,7,7,5 \| 7,6,7,7] | F2,F3,Z,Q,R,C |
| ⟨10,26,29⟩ | 4482 | 4482 | 4487 | **4482** | 0 | R[⟨2,4,4⟩ AT-Z rep; 5,5 \| 6,6,7,7 \| 8,7,7,7] | F2,F3,Z,Q,R,C |
| ⟨10,27,27⟩ | 4322 | 4322 | 4326 | **4276** | **−46** | R[⟨2,4,4⟩ P-ZT; 5,5 \| 7,7,7,6 \| 7,6,7,7] | F2,F3,Z,Q,R,C |
| ⟨13,23,31⟩ | 5396 | 5396 | 5398 | **5396** | 0 | R[Strassen; 7,6 \| 12,11 \| 16,15] over ⟨7,12,16⟩:876 | Q,R,C |
| ⟨13,23,32⟩ | 5552 | 5552 | 5556 | **5552** | 0 | R[Strassen; 7,6 \| 12,11 \| 16,16] | Q,R,C |
| ⟨13,24,31⟩ | 5540 | 5540 | 5544 | **5540** | 0 | R[Strassen; 7,6 \| 12,12 \| 16,15] | Q,R,C |
| ⟨13,25,32⟩ | 6008 | 6008 | 6012 | **6008** | 0 | R[Strassen; 7,6 \| 13,12 \| 16,16] | Q,R,C |
| ⟨13,29,32⟩ | 6910 | 6910 | 6914 | **6910** | 0 | R[Strassen; 7,6 \| 15,14 \| 16,16] over ⟨7,15,16⟩:1081 | Q,R,C |
| ⟨13,30,31⟩ | 6913 | 6913 | 6919 | **6913** | 0 | R[Strassen; 6,7 \| 15,15 \| 16,15] | Q,R,C |
| ⟨13,31,32⟩ | 7322 | 7322 | 7326 | **7322** | 0 | R[Strassen; 7,6 \| 16,15 \| 16,16] | Q,R,C |
| ⟨15,31,32⟩ | 8185 | 8185 | 8189 | **8185** | 0 | R[Strassen; 7,8 \| 16,15 \| 16,16] | Q,R,C |

Totals: the contribution saved 170 multiplications over our catalog; the
re-derivation saves **170 + 382 = 552**. ⟨6,14,25⟩ stays +2 over FMM's *index*
(1320, which moved below the recipe's 1322 after the issue was written; that
1320 is index-only — no artifact). The Strassen-outer stubs are `[Q,R,C]`
because a leaf carries `1/3` (no F₃) and another `1/2` (no F₂) — the
contributor's files for those say `F3,Q,R,C`; ours are stamped from the actual
leaves and are the conservative intersection.

## What we verified

1. **The contributed explicit files** — every one through the
   `ImportContributedSchemes` gate (dry-run, nothing written): exact BigInteger
   symbolic proof of the bilinear identity, random spot-check, fields from the
   coefficients; 30 s–6 min per file, the dense rank-8185 ⟨15,31,32⟩ being the
   slowest. Result: **22/22 exact** (run 2026-09-30; the per-file verdict lines
   are reproduced by re-running the gate on the attachment — dry-run, no
   `--execute`). Coefficient check: the 13 "rational" files
   are `±1/2` (⟨13,·,·⟩, ⟨15,31,32⟩, ⟨6,14,25⟩ → F₃ yes / F₂ no) or `±1/3, ±2/3`
   (⟨10,19,31⟩, ⟨10,22,29⟩, ⟨10,23,30⟩, ⟨10,23,31⟩ → F₂ yes / F₃ no) — the
   issue text says "1/2" for all 13, the files' own `fields` headers agree with
   *our* coefficient-driven result, not with the issue text.
2. **Re-derivation probe** — `docs.explore.ProbeIssue8Recipes` (read-only):
   per target, our leaf cost for every recipe piece (all 22 recipes re-sum
   exactly with our leaves — no piece gap), the exact allocation optimizer
   over every ⟨2,4,4⟩=26 / ⟨2,2,2⟩=7 scheme we hold in every orientation
   (12 / 48 variants), and the existing stub's replay. Verdict: 22/22
   REPRODUCED, 13 strictly better.
3. **Materialisation** — `SchemeSweep --mode=materialize --field=Q
   --config=rectangular --baseFilter=2,4,4 --strategies=recomb --maxNodes=60000000`
   for the 14, `--config=simple --baseFilter=2,2,2 --strategies=recomb` for the
   8: 22 wins, 0 errors; each construction spot-checked before its stub was
   written. Then `VerifyOneScheme` on each committed stub (replay + strongest
   affordable check): **16 EXACT_SYMBOLIC proofs**, **6 RANDOM_SPOT_CHECK**
   (⟨13,23,31⟩, ⟨13,24,31⟩, ⟨13,30,31⟩, ⟨13,29,32⟩, ⟨13,31,32⟩, ⟨15,31,32⟩ —
   term count above the 30M exact cap; the spot-check is not an algebraic
   proof, per the optimality discipline).
4. Guards: `TestSweepSpotsSota.retains_issue8_lille_recipe_ranks` (rank ≤ ours
   for all 22), `TestSweepSpotsSota.compute_pipeline_reaches_6x14x25_1322_via_244_root`
   (the ⟨2,4,4⟩ root reps compose 1322 from disk leaves),
   `TestRootPoolContents.poolHasTwoContentDistinct244Roots`.

## Why the engine beats the recipes — and why it hadn't before

The recipe's rank is `Σ_k R(effective dims of outer product k)`; the
*support structure* of the outer ⟨2,4,4⟩:26 scheme decides which products
touch the small block on each axis. FMM's HK file and our `hk71` emitter share
one support structure; the catalog also holds **content-distinct** ⟨2,4,4⟩=26
schemes (Perminov-ZT, AlphaTensor-Z, Kauers-2026, flips) whose supports tile
differently — the Perminov-ZT rep lets more products land on the 4- and 5-blocks
(−14…−46 per shape). Our sweeps had used ⟨2,4,4⟩ only through the extended pool
(where `--baseFilter`-less thorough runs starve the per-base allocation budget)
or `--base=2x4x4`, which pins **one** file (`findWithSource`) — the wrong rep.
Fix shipped with this record: **both reps are now `rootPool()` bases** (user
request on issue #8: "register ⟨2,4,4⟩ as an additional base"); `poolContentKey`
keeps both because their supports differ. The 8 Strassen-outer shapes needed no
engine change — only the Perminov ⟨7,12,16⟩=876 / ⟨7,15,16⟩=1081 pieces
(synced 2026-09) and a re-materialisation, since the old stubs pin leaves by hash.

## Why not import the explicit files

Catalog policy (`CatalogPolicy.MATERIALISE_MAX_DIM`): above dim 16 a derived
scheme is stored as a lineage stub, replayed on demand. All 22 are
compositions of catalog pieces we hold; committing 321 MB of dense matrices
would duplicate what the stubs reproduce exactly, and the "prefer derivable
over imported" rule applies. The contributed files served as the independent
cross-check (their ranks are the bar the guards assert against) and stay on the
issue as the contributor's artifact. Attribution: the stubs are ours
(`derived`, lineage-pinned); the *recipes* are FMM-Lille's; the *first explicit
witnesses at the recipe ranks* are the contributor's — recorded here and in
REFERENCES.md [89], not in the stubs' `source`.

## Follow-ups

- The remaining 14 of the contributor's 36 (recipes whose listed pieces sum
  above the table rank; allocations outside the smallest-size rule; DIS09-style
  Kronecker-with-correction such as ⟨16,20,28⟩=4944) — the last group is our
  `references/fmm-report-16x20x28.md` topic; unchanged.
- ⟨6,14,25⟩: FMM index 1320 vs ours 1322 — index-only, no artifact; ask upstream.
- Regenerate `references/fmm-cross-check.md` / `docs/comparison/*` (stale for
  unrelated digest churn since 2026-07) — the ⟨10,·,·⟩ rows flip from WORSE to
  BETTER.
- Broader: re-run the closure over the 17–32 band with the two ⟨2,4,4⟩ root
  reps in the default pool (`SchemeSweep --mode=closure --field=Q
  --config=rectangular`) — every ⟨2,4,4⟩-outer recombination stub in the catalog
  was optimised against the hk71 rep only.
