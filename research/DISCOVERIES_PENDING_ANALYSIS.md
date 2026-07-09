# Discoveries pending analysis

Scratch file for experimental results that should NOT be lost if a
process is killed mid-investigation. Once a finding here is digested
into the catalog (cited-bounds.json / derived-bounds.json / a
materialised scheme JSON), move it from here to its permanent home.

---

## 2026-07-09 — INDEX-RECIPE decode #2: SEVEN more closed (two BEAT FMM) — most of the "absorbing class" was the POOL-CAP pattern in disguise; WORSE 16 → 9

**Key decode:** in FMM's sum displays, **leaf count = outer base rank**. The
"unreadable" 40/68/127/150-leaf displays are plain heterogeneous
recombinations over bases our `maxBaseDim=5` pool cap excludes:
⟨3,5,6⟩:68 (four rows!), ⟨5,5,7⟩:127, ⟨5,6,7⟩:150 — NOT shared-mul devices.
Closed via `--base` runs, all exact-verified where in term-cap:
- ⟨11,20,23⟩ **3035 — BEATS FMM 3039** (their own recipe re-priced with our
  leaves lands below their index!), exact proof;
- ⟨20,20,27⟩ **5986 — BEATS FMM 6006 by 20**, spot-check;
- ⟨12,20,23⟩ 3184 (tie, proof), ⟨11,20,24⟩ 3114 (tie, proof),
  ⟨15,25,29⟩ 6188 (tie, proof), ⟨20,23,28⟩ 7100 (tie, spot);
- ⟨3,29,29⟩ **1840** (tie, exact proof) — base ⟨2,5,5⟩:40 is in-pool but our
  hk71 rep has the wrong SUPPORT profile (per-product leaf dims); imported
  FMM's rep to bud-bases/section5 (the support-representative sibling of the
  bud-representative caveat). Guard rows for all four base families.
Cascade pass: 0 further (fixpoint). **WORSE 16 → 9, BETTER 1474 → 1476.**
Verification: 5 exact proofs + 2 spot-checks, 0 failures.

**The REAL remaining 9:** ⟨17,17,19⟩/⟨13,19,29⟩/⟨14,14,29⟩/⟨19,19,31⟩
(overlap-1 shared-mul family — the quantified ⟨17,17,19⟩ witness class),
⟨19,19,22⟩/⟨20,23,23⟩ (TA-family), ⟨7,11,30⟩ (overlap-2), ⟨8,27,30⟩
(dubious display, axis-mismatched recipe), ⟨9,11,22⟩ (display leaf-sorting
ambiguity, base undecoded — 38 leaves, sum exact; possibly ⟨3,4,4⟩:38 with
axis-permuted leaves, concentration run found nothing).

---

## 2026-07-09 — FmmCrossCheck now splits WORSE into artifact-backed (16) vs UPSTREAM-UNVERIFIED (20)

`FmmCrossCheck` consumes the new curated `references/fmm-artifact-audit.json`
(17 index_only + 3 placeholder shapes from the artifact audit) and moves those
rows out of `## WORSE` into a `## UPSTREAM-UNVERIFIED` section — named without
the word "WORSE" so the fmm-gap skill's section-scoped random picker never
selects them. Headline now honest: **WORSE=16 (real, artifact-backed),
UNVERIFIED=20, BETTER=1474**. Guard: `TestFmmArtifactAudit` (audit JSON
well-formed, keys sorted n≤m≤p — an unsorted key would silently fail to match
and fall back into WORSE). Skill §1 updated (re-audit rule on digest refresh).
The 16 real rows are the shared-correction-recombination frontier (see the
ANATOMY entry below) minus the already-blocked serendip stragglers.

---

## 2026-07-09 — CORRECTION (user catch): ⟨4,9,10⟩/⟨4,7,8⟩ "base-level mismatches" were OUR mirror artifacts; ⟨8,27,30⟩ recipe found unrealizable from FMM's published base

**User checked fmm.univ-lille.fr/4x9x10.html** (Perminov-credited, rank 250)
against our claim of an FMM artifact at 255. The 255 was PERMINOV'S GitHub
mirror copy — `tools/import_fmm_maple.py` prefers the mirror and only falls
back to FMM on 404. FMM's own artifacts, re-downloaded directly: ⟨4,9,10⟩ =
**250** and ⟨4,7,8⟩ = **161**, both exactly their index ranks — both
"base-level mismatch" lines RETRACTED from the report
(`fmm-artifact-audit.md` CORRECTION section). Provenance re-check: all 14
non-parametric PHANTOM rows 404 on the mirror → their audit counts came from
FMM's own files; the phantom classifications STAND. (Residual upstream-of-
upstream note: Perminov's own tensor/ mirror carries a 255 where FMM credits
him with 250.)

**⟨8,27,30⟩ follow-up with the genuine 250:** bud census (exact directions)
= **6 size-2 U-classes, zero V, zero W**. The index recipe
`(⟨4,9,10⟩:250 − 12) ⊗ ⟨2,3,3⟩:15 + 6·⟨4,3,3⟩:29 = 3744` needs six
shared-B (V-direction) classes — ⟨4,3,3⟩ is an n-axis enlargement; the
artifact's U-classes fuse to ⟨2,3,6⟩:34, a LOSS (vs 2·15=30). And no
orientation can swap the roles (all dims distinct; only (4,9,10) divides
(8,27,30) per-axis). So the display recipe is NOT realizable from FMM's own
published base — the ⟨16,20,28⟩ display-inconsistency class, sneakier
variant: the "−12" COUNT matches the artifact's bud census but the AXIS
doesn't. ⟨8,27,30⟩ moves from "flip-walk viable" to the dubious-display
report list. Bucket-3 remainder: ⟨20,24,25⟩, ⟨16,20,29⟩ only.

---

## 2026-07-09 — ABSORBING-CLASS ANATOMY: the device is CROSS-LEAF CORRECTION SHARING (first quantified witness, ⟨17,17,19⟩:3266)

**Method:** block-support fingerprint of downloaded FMM artifacts
(`tools/fingerprint_fmm_artifact.py`, NEW): partition a scheme's products by
which 2×2 blocks their U/V/W supports touch under the display's axis splits.

**⟨17,17,19⟩:3266 (display: 7 leaves summing 3267, so a 1-saving):** under
splits A=8/B=8/C=9 the artifact IS display-faithful — **7 main clusters**
sized 519/499/482/422×3/417, each at or BELOW its standalone leaf rank
(422 vs R(8,9,9)=430 ×3; 417 vs R(8,8,10)=427; 499/519 vs R(9,9,10)=534 ×2;
482 = R(8,9,10) exactly), plus **83 residue products** with cross/subset
signatures. The books close exactly: Σ(leafRank − clusterSize) = 24+10+35+15
= 84 products of leaf work replaced by 83 shared corrections — **saving
exactly the 1** separating 3267 from 3266. THE DEVICE: recombination whose
leaves are NOT independent sub-schemes; a small correction set is shared
across leaves (our engine prices leaves independently, hence can never
express this class). Some clusters combine 3 blocks per operand (not plain
Strassen — likely a GL-orbit/flip variant of a ⟨2,2,2⟩:7 base).

**⟨19,19,22⟩:4536:** NOT display-faithful under any split order (largest
class: 1460 products spanning 3 blocks on every operand) — basis-mixed
flip-graph output; the display is a post-hoc description. Artifact-anatomy
extraction only works on the faithful ones — fingerprint FIRST.

**Status / next steps (dedicated session):** (1) proper residue→leaf
attribution needs bilinear algebra (W-output attribution), not support
masks (my containment heuristic left 16 products ambiguous); (2) recover
the outer base (which ⟨2,2,2⟩:7 orbit member) and the per-cluster basis;
(3) formulate the sharing identity; (4) generic builder + search
integration ("shared-correction recombination"). Also drafted the user's
upstream report note for the inconsistent ⟨16,20,28⟩ display:
`references/fmm-report-16x20x28.md`.

---

## 2026-07-09 — fmm-gap INDEX-RECIPE sweep: SIX closed at exact ties (⟨6,25,32⟩ 2870, ⟨9,30,32⟩ 4851, ⟨8,21,21⟩ 2156, ⟨12,20,28⟩ 3782, ⟨18,20,21⟩ 4259, ⟨20,21,21⟩ 4967); WORSE 42 → 36; 20× sweep concurrency

**Method:** parsed FMM's INDEX-page `<td content=description>` cells (unquoted
attrs; fresher than the detail pages) for all 42 WORSE rows and priced every
recipe with OUR ingredient ranks. 9 serendipitous recipes priced exactly
FMM's rank; 7 plain sums priced consistent. Engine attempts: the 9 serendip
rows all failed on OUR base representatives (bud-representative gap) →
imported the 5 distinct FMM bases; TWO were useless (⟨4,9,10⟩ artifact 255
vs recipe's 250, ⟨4,7,8⟩ 164 vs 161 — MORE index-vs-artifact staleness, added
to the report material; our own catalog already holds 250/161). The THREE
matching bases (⟨3,5,8⟩:90, ⟨4,7,7⟩:144, ⟨3,5,7⟩:79 — bud-bases/, Kauers–Wood
era, attribution TBD from biblio) closed SIX rows at exact ties, one
(⟨20,21,21⟩) again requiring the ambiguous-orientation fix
(OrientAs(⟨5,3,7⟩)). ⟨3,5,7⟩:79 alone pays for three rows. The 7 sum rows: 0
wins (absorbing-class structures, as before). ⟨20,24,25⟩/⟨16,20,28⟩ still
blocked (known bud-representative / inconsistent-recipe cases). Verification:
1 exact proof (⟨8,21,21⟩) + 5 spot-checks (dense), 0 failures. Cascade pass:
0 (fixpoint). **WORSE 42 → 36.**

**Concurrency (user flagged low CPU):** three fixes, all guarded green —
(1) `bestViaAllocationOptimizer` now fans the per-base B&B over a parallel
stream (AtomicLong shared incumbent for cross-base pruning; deterministic
(rank, pool-index) pick); (2) `RecursiveClosureSota` made concurrent-safe
(ConcurrentHashMap cache — getRank is pure so duplicate computes are benign;
cycle-guard moved to ThreadLocal so sibling computations don't false-positive
to directRank); (3) `SchemeSweep.runMaterialize` submits heaviest-volume
shapes FIRST (the dim-27+ stragglers ran the tail near-single-threaded).
Measured: 72% → ~1700% CPU (0.7 → 17 of 18 cores) on the same 36-row sweep.

---

## 2026-07-09 — fmm-gap ⟨21,28,30⟩ + ⟨4,18,30⟩: CLOSED (ties 9473 / 1394) — NEW root cause: AMBIGUOUS-ORIENTATION masking in budBasesAt; ⟨16,20,28⟩ recipe proven inconsistent; WORSE 44 → 42

**Trigger:** user pointed at FMM's INDEX-page description for ⟨21,28,30⟩:
`(⟨7,7,5⟩:176 − 21) ⊗ ⟨3,4,6⟩:54 + ⟨15,4,6⟩:263 + 8·⟨6,4,6⟩:105 = 9473` —
arithmetic exact with OUR ingredient ranks (all match), base = the round-6
Kauers–Wood ⟨5,7,7⟩:176 import. Yet the engine stopped at 9477.

**Root cause (engine, silent):** ⟨7,7,5⟩ is reachable from a ⟨5,7,7⟩ file by
TWO distinct S₃ orientations (dims repeat) with DIFFERENT bud profiles —
probe: BCA orientation prices 9487, CBA prices **9473**. `budBasesAt` offered
only `orientAs`'s first match (BCA), so the σ-paying sibling was never fed to
the search. **Fix:** budBasesAt now offers EVERY S₃ orientation reaching the
target dims (`ORIENT_PERMS` loop); `BudBase` carries the exact perm; the pin
records an explicit `axisMap` OrientAs (new `durableLeafRef(…, exactPerm)`
overload) so replay is dim-repeat-deterministic — stub
`21x28x30-r9473-f794828` pins `OrientAs(axisMap=pmn, 5x7x7@a88ab92…)`.
Guards: `TestSerendipitousStubInner.ambiguous_orientation_offers_both_bud_profiles`
+ `…_stub_replays_exactly`; rows ⟨21,28,30⟩≤9473, ⟨4,18,30⟩≤1394.
Broad re-harvest (WORSE ∪ perm-multiples, 185 shapes): +1 collateral,
**⟨4,18,30⟩ = 1394** (tie with FMM index, BEATS their published artifact
1397). Verification: 1 exact proof (⟨4,18,30⟩) + 1 spot-check (⟨21,28,30⟩,
dense), 0 failures. **WORSE 44 → 42.**

**Audit revision:** "phantom" rows split — where an INDEX-page recipe
verifies, the right upstream report is "stale artifact", not "index
over-claim" (⟨21,28,30⟩ and ⟨4,18,30⟩ now reproduced by us; artifacts 9782 /
1397 are stale). See `references/fmm-artifact-audit.md` REVISION section.

**⟨16,20,28⟩ dig (user suggestion — recipe despite placeholder artifact):**
slice-rank analysis of FMM's published ⟨4,5,7⟩:104 base (exact rationals):
U-classes = 7×size-2 + 1×size-3 = **exactly 17 terms** (matches the "−17"),
but every class is slice-rank-TIGHT (the size-3 class's M_u has rank 3 → no
⟨1,1,2⟩ embedding exists), and 8×⟨4,4,8⟩ can hold at most 16 terms. So the
displayed recipe `(104−17)⊗48 + 8·96 = 4944` is TERM-COUNT INCONSISTENT and
unattainable from this base under bud semantics (best expressible: 4989 via
7×⟨4,4,8⟩ + 1×⟨4,4,12⟩). Either the display is corrupted (add to the FMM
report — its artifact is also a placeholder) or the true device fuses in the
PRODUCT space (base⊗inner term level — a finer mechanism than base-level
buds; unexplored).

---

## 2026-07-08 — fmm-gap BATCH: full WORSE-list gap-tree census; the ⟨2,2t,3t⟩ "HK family" rows are FMM INDEX-ONLY OVER-CLAIMS (their artifacts match or LOSE to ours)

**Scope:** user directive "more FMM gaps until all gaps are done" → recursive
gap-tree diagnosis of all 44 WORSE rows (leaves fetched from per-shape FMM
pages, compared to our ranks, deficient chains followed).
**Census tally (by root class): 23 OVERLAP** (sub-additive constructions —
absorbing-pad / staircase family, incl. the ⟨27,·,·⟩ run whose 2-leaf concats
save 18–26 below plain sums), **11 SERENDIP** (mostly span-compressed-bud
blocked, see previous entry), **4 TA-family, 4 unparsed (⟨27,28,·⟩ family
pages), 2 mis-classified "closable"** (⟨3,29,29⟩/⟨9,11,22⟩ — their 38/40-leaf
lists cannot come from any partition grid: an all-single-A-support rank-40
⟨2,5,5⟩ would direct-sum to ≥50; absorbing family again).

**The ⟨2,2t,3t⟩ t=6…10 + ⟨2,24,30⟩ rows (6 rows) — RESOLVED as
NOT-OUR-GAP:** all six FMM index ranks (333/452/588/743/915/1095) are exactly
the HK71 closed form `⌈(3pn+max)/2⌉`; their page "decompositions" are
fallback displays summing 1–5 HIGHER. Downloaded FMM's actual maple tensor
artifacts: `⟨2,12,18⟩` artifact = **334 = ours** (+1 over index);
`⟨2,16,24⟩` artifact = **592, WORSE than our 590** (+4 over index). Perminov
carries none of the six. So the index ranks are UNATTAINED formula citations
— consistent with our task-#9 impossibility results
(`research/hopcroft-kerr-2np/CONSTRUCTIVE_METHOD.md`: the g ≥ 6 triangle
family is formula-impossible within the local-atom framework; the arc-sum +
(3,3,bridge) theorems). FMM's own artifacts sitting at exactly our ranks is
empirical corroboration. **Action items:** (1) these six rows should be
EXCLUDED from actionable-WORSE (or FmmCrossCheck should annotate index-vs-
artifact discrepancies — extend `fmm-lille-discrepancies.md` tooling);
(2) worth reporting upstream to Sedoglavic (user's contact) — either the
index cites HK's paper bound optimistically, or an unpublished attaining
scheme exists and the artifact link is stale; (3) task #9's "resume points"
gain urgency only if upstream confirms an attaining scheme exists.

**Conclusion of the batch pass:** ZERO remaining WORSE rows are closable
with current engines — every row is blocked on one of: (a) absorbing-pad /
overlap constructions (~27+2 rows, biggest class), (b) span-compressed buds
(~11), (c) TA-family builder (4), (d) the 6 index-only HK rows (not real
gaps). "All gaps done" now means: build (a) and (b) as engine features, and
reclassify (d). Analysis imports were kept OUT of the catalog trees and
deleted after use per policy.

**Follow-up (same session) — FULL ARTIFACT AUDIT of all 44 rows**
(`references/fmm-artifact-audit.md`): downloaded every row's FMM tensor
artifact and compared its true rank to the index. **14 rows are PHANTOM**
(artifact ≥ ours; e.g. ⟨21,28,30⟩ artifact 9782 vs our 9477; the three
⟨27,28,·⟩ artifacts equal ours exactly) — the index over-claims and the
cross-check overstates our deficit. **≥3 rows are recipe-only placeholders**
(`Tensor:=Tensor:`, no data: ⟨16,20,28⟩/⟨16,20,29⟩/⟨20,24,25⟩). **≈25 rows
are REAL** (artifact attains index — explicit schemes exist that our engines
can't derive; the downloaded artifacts are the raw material for extracting
the absorbing-pad/overlap device). 5 downloads failed (⟨27,·,3x⟩ family +
⟨4,18,30⟩, likely placeholders too — re-check). Action items in the audit
doc: annotate phantoms in FmmCrossCheck, report upstream, mine the real
artifacts for the device.

---

## 2026-07-08 — fmm-gap ⟨16,20,29⟩→⟨16,20,28⟩: NOT closed — needs SPAN-COMPRESSED bud support in the serendipitous engine

**Target:** `Q⟨16,20,29⟩:m` (ours 5272, FMM 5264) = naive ⟨16,20,1⟩:320 +
`Q⟨16,20,28⟩:m` (ours 4986, FMM 4944, +42) — the whole gap is the ingredient.

**Diagnosis:** FMM's ⟨16,20,28⟩:4944 = `(⟨4,5,7⟩:104 − 17) ⊗ ⟨4,4,4⟩:48 +
8·⟨4,4,8⟩:96` (base = Kauers–Wood arXiv:2505.05896). 17 terms absorbed into
8 blocks of ⟨4,4,8⟩ — since size-2 buds save ZERO here (⟨4,4,8⟩=96=2·48), the
48-saving requires at least one 3-term class fused into a ⟨4,4,8⟩ block:
**span compression** (the class's enlarged width collapses below k·p₂ via
linear dependence among the factor vectors). Our `findBuds`/`costOf` are
span-BLIND (enlarged width = class SIZE × p₂), so: imported base's best
greedy σ → 4989; both catalog 104s → 4989; flip-graph walk (§4c.4 playbook,
3 seeds × 4 rngs × 30k steps, `FlipObjectives.serendipitous(Q,4,4,4)`)
plateaued at 4986 = incumbent. KEY NEGATIVE INSIGHT: the walk CANNOT succeed
while the objective itself is span-blind — even on FMM's exact representative
it would price 4989. Imported bud base purged per PURGE_REFCOUNT_POLICY
(produced nothing); `ProbeFlipWalkSerendip162028` kept as the negative record.

**Future work (well-defined engine feature):** span-compressed buds in
`SerendipitousBudProduct` — findBuds detects classes whose enlarged block
collapses (rank of the spanned factor subspace < class size), costOf prices
the collapsed width, productFromDecomposition builds the recombined block.
Would re-open ⟨16,20,28⟩ (−42), cascade to ⟨16,20,29⟩ (−8), and likely other
serendip-census rows.

---

## 2026-07-08 — fmm-gap ⟨7,14,32⟩ + ⟨14,17,30⟩: CLOSED (ties 2017 / 4180) — TWO more stub-blindness holes: diskBest leaf fetch + the overlay's oriented-hash pin; WORSE 46 → 44

**Targets:** `Q⟨7,14,32⟩:m` (ours 2020, FMM 2017) and collateral
`Q⟨14,17,30⟩:m` (4184 → 4180). All bounds. Re-scoped from ⟨17,17,19⟩ after
the overlap census (see the ⟨19,19,22⟩ entry below).

**Root cause #1 (diskBest stub-blind):** FMM's ⟨7,14,32⟩ = plain concat
`⟨7,14,2⟩:152 + ⟨7,14,30⟩:1865` — both ingredients ours (1865 is a round-7
stub). The scorer priced it (search incumbent 2018 visible in the round-7
logs) but the composed build threw `construct: missing sub-algorithm for
⟨7,14,30⟩`: `RecursiveLookup.find → diskBest → findWithSource` skips
lineage-only stubs — pricing includes stubs, the build couldn't fetch them
(the recombination-builder sibling of the serendipitous stub-blindness).
**Fix:** `diskBest` now detects stub-is-best (`findRank <` best dense file)
and resolves via the stub-replaying `resolveParentHit`, pinning the leaf with
`durableLeafRef`. Guards: `TestDiskBestStubLeaf` + row ⟨7,14,32⟩≤2017 in
`TestSweepSpotsSota`.

**Root cause #2 (dangling oriented-hash pin):** with #1 fixed, ⟨14,17,30⟩
built at 4180 but persist REFUSED: the overlay's corrupt-lineage fallback
pinned the ⟨8,9,12⟩ leaf by `contentHash(oriented alg)` — stamped on NO file
→ DANGLING (exactly the anti-pattern trySerendipitous's comment documents).
**Fix:** `diskBest` re-pins corrupt-lineage files by FILE hash via
`durableLeafRef` before the overlay ever sees them.

**Sweep to fixpoint over the 46-row list:** ⟨7,14,32⟩=2017 (tie, leaf
born-pinned to the 7x14x30 stub by full content hash), ⟨14,17,30⟩=4180
(tie); final pass 0 wins, 0 errors. Verification: **1 exact symbolic proof
(⟨7,14,32⟩) + 1 spot-check (⟨14,17,30⟩, dense), 0 failures**.

**Net cross-check: WORSE 46 → 44 (−2), BETTER 1474 (unchanged).**

**Addendum — overlap sub-census of the 31 block rows** (extends the census
below): 21 ZERO-overlap (plain leaf sums; after this round's fixes the
survivors are genuine INGREDIENT gaps — our best rank for some leaf exceeds
FMM's; each needs its own neighbour-diff chase). 10 POSITIVE-overlap
(sub-additive: the ⟨2,2t,3t⟩ t=6…10 run is the unported Hopcroft–Kerr
⟨2,m,n⟩ closed-form parities — porting them would close 5 rows at once;
⟨17,17,19⟩ / ⟨7,11,30⟩ are 1–2-saving overlaps of the same flavour as the
TA-family rows).

---

## 2026-07-08 — fmm-gap ⟨19,19,22⟩: NOT closed (gap 2 stands) — needs an ABSORBING-PAD Pan-pair builder; mechanism census of the WORSE list

**Target:** `Q⟨19,19,22⟩:m` — ours 4538, FMM 4536 (+2). All bounds.

**Diagnosis (decoded, construction not attempted):** FMM's 4536 =
`⟨9,9,12⟩:600 + ⟨9,9,10⟩:534 + ⟨10,10,12⟩:766 + ⟨10,9,12⟩:668 + ⟨9,10,12⟩:668
+ TA(⟨10,10,10⟩,⟨10,10,10⟩):1300` — 7 terms for a 2×2×2 grid's 8 blocks
(A=[9,10], B=[9,10], C=[12,10]). Our leaf ranks match ALL five explicit leaves
exactly. The TA pair (10³+3·10²=1300) covers the two cross-blocks
⟨9,10,10⟩/⟨10,9,10⟩ *padded to cubes with REAL data* such that the
over-computation absorbs the eighth block ⟨10,10,10⟩ — plain-sum accounting
(pair 1180 + solo 651) prices 5067, so the padding-absorption isload-bearing.
Note ⟨19,19,22⟩ = ⟨2n−1,2n−1,2n+2⟩ at n=10 — a parametric family.

**Why our engines can't reach it:** `PairFusedRecombination` is explicitly an
all-cubic MVP (uniform ⟨k,k,k⟩ leaves over cubic targets); the naive-grid
pairing sweep (`PairedSubProducts.applyPairing`) prices exact-rotation pairs
only — no padding, no absorption. `buildPeeledViaTa` covers only the
⟨N,N+s,N+s⟩ peel. Fired the cheap engines anyway: GL-orbit probe → tie 4538;
projection scatter → 0; recomb/serendip/proj thorough sweeps already visited
this shape twice in the 2026-07-08 rounds. Our 4538 = Strassen-style
`R[2x2x2; 9,10|10,9|12,10]`.

**Future work (logged, not started):** generalise
`RectangularTrilinearAggregation.build` to the deficient 2×2×2 grid
(absorbing-pad Pan pair) and register it as a `ConstructiveMethod`.
**Mechanism census of the 46-row WORSE list** (per-shape FMM pages):
**3 TA-family** (⟨19,19,22⟩, ⟨20,23,23⟩, ⟨22,23,23⟩ — all need this builder),
**32 block decompositions** (engine-reachable in principle; each needs
individual ingredient diagnosis), **10 serendipitous**, 1 other. Census
script + results: `$CLAUDE_JOB_DIR/tmp/worse8-mechanisms.txt` (job-local).

---

## 2026-07-08 — fmm-gap ⟨14,27,27⟩: CLOSED (tie 5862) — NEW root cause: allocation B&B STAGNATION-STARVED on wide axes; WORSE 85 → 46

**Target:** `Q⟨14,27,27⟩:m` — ours was 5889, FMM 5862 (+27). All bounds.
(Context: the CI merge `994d0bb3` refreshed the FMM digest — 80 shapes improved
UPSTREAM, which moved the cross-check from 57 to 85 WORSE / 1493→1469 BETTER
with no regression on our side. FmmCrossCheck now orders sections
shape-lexicographically, per user request — families read as runs.)

**Diagnosis:** FMM's 5862 = plain ⟨2,3,3⟩:15 grid `[8,6 | 9,9,9 | 9,9,9]` =
9·⟨8,9,9⟩:430 + 6·⟨6,9,9⟩:332 — and BOTH ingredients were already ours (the
332 arrived in the same CI merge). Pure closure gap, yet thorough AND
`--baseFilter=2x3x3` sweeps returned "no materialisable strategy".

**Root cause (engine, systematic):** `AllocationOptimizer`'s nested A→B→C
sweep with balance-first per-axis order + `ALLOC_STAGNATION=100k`. For a
⟨·,27,27⟩ target the B×C completion space under a SINGLE A-allocation is
C(26,2)² ≈ 105k nodes > the whole stagnation window — the B&B dies inside the
balanced A=[7,7] subtree (returns 5940) and A=[6,8] is NEVER tried. Every
target with two ~27-wide 3-part axes was starved the same way; the fresh
FMM ranks in the ⟨17–18,·,·⟩ band were exactly such grids. Diagnosis method
(fast, after the user pushed back on slow probing): call
`AllocationOptimizer.optimize` directly with a disk-only SotaResolver —
uncapped finds `[6,8|9³|9³]=5862 exhaustive=true` (1.38M nodes ≈ full space);
the sweep's stagnation budget stops at 100k with balanced 5940.

**Fix:** coordinate-descent SEEDING in `AllocationOptimizer.optimize` —
before the nested sweep, price every one-axis deviation from the greedily
updated incumbent (|uA|+|uB|+|uC| ≈ hundreds of cost() calls). Single-axis-
unbalanced optima (FMM's standard grid shape) become the incumbent regardless
of the stagnation cap, and the tightened threshold sharpens all downstream
pruning. The exact B&B is unchanged (pure seeding). NOT implemented (noted as
the principled alternative): global total-imbalance ordering via lazy
sorted-sums heap — it would break the nested-prefix lower-bound structure.
Guards: `TestAllocationOrdering` (hermetic, 0.07s — spells the per-axis visit
order literally + proves seeding finds a single-axis optimum under a 1-node
sweep budget), `TestAllocationOptimizer.stagnation_capped_run_finds_single_axis_unbalanced_optimum`,
rows ⟨14,27,27⟩≤5862 / ⟨17,22,29⟩≤6125 in `TestSweepSpotsSota`.

**Collateral (full-strategy sweep over the 85-row WORSE list, fixpoint at
pass 2):** target + **39 pass-1 wins** — the seeding-unlocked ⟨17–18,·,·⟩ band
en bloc (⟨17,22,29⟩ 6125, ⟨15,22,31⟩ 5844, ⟨22,27,31⟩ 9936, …), re-taking the
upstream-refreshed FMM ranks at parity. Pass 2: 0 (fixpoint). Zero errors.

**Net cross-check: WORSE 85 → 46 (−39), BETTER 1469 → 1474 (+5).**
Verification over the 40 new stubs: **3 exact symbolic proofs + 37
spot-checks (dense dim-27+ composites), 0 failures.** All born-pinned.

---

## 2026-07-08 — fmm-gap ⟨20,28,28⟩: CLOSED and beaten (8440 → 8434, FMM 8438) — NEW root cause: STUB-BLIND serendipitous search; WORSE 86 → 57

**Target:** `Q⟨20,28,28⟩:m` — ours was 8440, FMM 8438 (+2). All bounds.

**Diagnosis:** FMM's page states `⟨20,28,28⟩:8438 = (⟨5,7,7⟩:176 − 5) ⊗
⟨4,4,4⟩:48 + ⟨4,4,20⟩:230` — a serendipitous product whose base (Kauers–Wood
2025, arXiv:2505.05896) carries a size-5 ⟨1,1,5⟩ U-bud (5·48 − 230 = 10 saved).
Our catalog's rank-176 representative (perminov_c2535_ZT) tops out at
`U⟨1,1,4⟩` (saving 4 → 8444): the bud-representative caveat. Imported the FMM
base to `bud-bases/section7/fmm-lille_5x7x7_r176_a3315.json` (exact-verified;
discovery:false, attribution Kauers–Wood). Its profile prices **8434 < 8438**:
the U×5 bud (−10) PLUS a V×4 bud FMM's own recipe leaves unfused (−4 via
⟨4,4,16⟩:188, surfaced by the V-first orderings).

**Root cause (engine, silent):** even fed the base, `SerendipitousSearch`
returned EMPTY. It priced candidates with a `findWithSource`-based oracle and
built through `findWithSource` — both blind to lineage-only STUBS, so any
candidate whose enlarged fusion target exists only as a stub (⟨4,4,20⟩=230 is
a ConcatCols stub) was silently dropped (predict = −1) or threw at build.
Serendipitous search has been quietly weakened catalog-wide since stubs became
the norm for derived shapes. **Fix:** predict via `findRank` (stub-inclusive,
= `costOf`'s oracle); build via a new `SerendipitousBudProduct.InnerResolver`
hook — `RecursiveMaterialiser` passes its stub-replaying `resolveParentHit`,
`LineageReplayer` passes `resolveShape` (so serendip stubs with stub inners
replay), Phase-2 catches unbuildable candidates and falls through. Two
supporting fixes: `StampImportHashes` (NEW driver — stamps content `hash` on
8 hash-less old-convention imports so `durableLeafRef` can pin them) and
content-hash fallbacks in the persist audit's `findFileByHash` + the manifest's
`existingPinnedKeys` (both were filename-token-only → resolvable pins flagged
DANGLING). Guards: `TestSerendipitousStubInner` (3 tests: bestFor builds
through the stub inner ≤8434; default resolver doesn't crash; the persisted
stub replays to exactly 8434) + rows ⟨20,28,28⟩≤8434, ⟨21,25,28⟩≤8118 in
`TestSweepSpotsSota`.

**Collateral (serendip harvest over 141 candidates = perm(5,7,7)-multiples ∪
WORSE list, then full-strategy recomb/serendip/proj sweep over the 86-row
WORSE list, both to fixpoint at pass 2):** 37 collateral wins. Serendip pass 1:
10 (⟨10,21,28⟩ −36, ⟨21,25,28⟩ −24, ⟨21,28,30⟩ −13, ⟨16,27,32⟩ −10, ⟨20,21,28⟩
−9, ⟨14,20,28⟩ −7, … — several from PRE-EXISTING bases the stub-blindness had
been suppressing, e.g. 2x4x7 → ⟨10,21,28⟩). Full sweep pass 1: 27 more (the
⟨19,2x,2x⟩ recomb family en bloc, ⟨28,31,31⟩ 14028, ⟨25,29,31⟩ 12064,
⟨25,29,32⟩ 12410, ⟨3,29,29⟩ 1843, ⟨12,12,29⟩ 2404, ⟨9,15,28⟩ 2232, …).
Pass 2 both engines: 0 (fixpoint). Zero errors.

**Net cross-check: WORSE 86 → 57 (−29), BETTER 1482 → 1493 (+11).**
Verification over all 38 new stubs: **14 exact symbolic proofs + 24
spot-checks (dense, above the 30M-term exact cap), 0 failures.** All stubs
born-pinned (base = `5x7x7@a88ab92…` content hash).

**Follow-up:** `VerifyOneScheme` batch mode — verifying 38 stubs cost 38 JVM
starts × ~15s lookup-index builds (monothreaded, user flagged it); a
multi-file mode building the index once (+ optional thread pool) would cut
the batch ~10×.

---

## 2026-07-07 — fmm-gap ⟨7,14,24⟩: CLOSED (tie 1514) — NEW root cause: pool starvation of the per-base alloc budget; WORSE 95 → 86

**Target:** `Q⟨7,14,24⟩:m` — ours was 1519, FMM 1514 (+5). All bounds.

**Diagnosis:** FMM's 1514 = outer **⟨2,4,4⟩:26** with the two-axis-uneven alloc
`[4,3 | 4,4,3,3 | 6,6,6,6]` (leaves 8×⟨4,4,6⟩:73 + 8×⟨4,3,6⟩:54 + 7×⟨3,4,6⟩:54 +
3×⟨3,3,6⟩:40). The base is dim-4 — IN the default pool (both content-distinct 26s)
— so this is NOT the pool-cap pattern. Root cause: **pool starvation** — under the
full thorough pool (~3.3k entries) the per-base allocation B&B budget prunes deep
two-axis-uneven allocs; a `--baseFilter=2x4x4-r26`-concentrated run (6-entry pool)
finds 1514 in seconds (`R[2x4x4; 3,4 | 3,4,3,4 | 6,6,6,6]`, exact-verified,
born-pinned). Skill §4b updated with the caveat.

**Collateral (concentrated per-family sweeps over the 95-row list):**
`2x4x4-r26` → 8 wins (⟨5,22,23⟩ 1649, ⟨6,23,23⟩ 1940, ⟨7,23,25⟩ 2517, ⟨5,23,23⟩
1718, ⟨5,23,25⟩ 1865, ⟨7,14,32⟩ 2020, ⟨3,16,32⟩ 1136, ⟨7,23,23⟩ 2321);
`2x4x5-r32` → 1 (⟨14,20,29⟩ 4708); `2x3x3-r15`/`3x3x4-r29`/`3x4x4-r38`/`3x4x5-r47`
→ 0 (already exhausted); cascade pass → 0 (immediate fixpoint).

**Net cross-check: WORSE 95 → 86 (−9), BETTER 1476 → 1482 (+6).** Guard row
⟨7,14,24⟩≤1514. Verification: **9 exact symbolic proofs + 1 spot-check
(⟨14,20,29⟩, dense), 0 failures** over the 10 stubs.

---

## 2026-07-07 — fmm-gap ⟨3,20,30⟩: CLOSED (tie 1300); dim-6/7 base battery drops WORSE 128 → 95

**Target:** `Q⟨3,20,30⟩:m` — ours was 1320, FMM 1300 (+20). All bounds.

**Diagnosis:** FMM's 1300 = outer **⟨2,5,6⟩:47** with deficient A-split 3=2+1
(30×⟨2,4,5⟩:32 + 17×⟨1,4,5⟩:20 — 17 of 47 products isolate the width-1 A-block).
Fifth pool-cap instance; ⟨2,5,6⟩ was a base never yet swept across the list.

**Reaction:** target closed at **1300 (tie)** via `--base=2x5x6`
(`R[2x5x6; 2,1 | 4×5 | 5×6]`, our 47 carries the same 17-isolation profile).
Then the full 8-base battery of previously-unswept dim-6/7 bases
(`2x5x6,3x3x6,2x4x7,3x3x7,2x6x6,2x5x8,2x6x7,4x4x6`) + cascade pass with all
strategies over the 128-row WORSE list at 10g: **39 collateral wins** (pass 1:
20 — mostly the thin-A ⟨3,·,·⟩/⟨5,·,·⟩ band, incl. ⟨7,20,24⟩ 2088→2086 closing a
leftover; pass 2: 19 cascades; pass 3: fixpoint, 0). Zero errors.

**Net cross-check: WORSE 128 → 95 (−33), BETTER 1460 → 1476 (+16).** Guard row
⟨3,20,30⟩≤1300 added. Verification (`verifyAuto`-tiered batch over the 40 stubs):
**39 exact symbolic proofs + 1 spot-check, 0 failures** (thin-A stubs are sparse
enough for the exact term-map).

---

## 2026-07-07 — fmm-gap ⟨14,28,29⟩: CLOSED and beaten (6506 → 6494, FMM 6498); 66-stub collateral harvest drops WORSE 187 → 128

**Target:** `Q⟨14,28,29⟩:m` — ours was 6506, FMM 6498 (+8). All bounds.

**Diagnosis (from FMM's 14x28x29.html, leaf-multiset arithmetic):** their 6498 =
outer **⟨3,4,6⟩:54** block-recursive with TWO deficient axes — 14=[5,5,4],
28=[7,7,7,7], 29=[5,5,5,5,5,4] — leaves 40×⟨5,7,5⟩:127 + 10×⟨4,7,5⟩:104 +
2×⟨5,7,4⟩:104 + 2×⟨4,7,4⟩:85. A dim-6 outer base: FOURTH instance of the
`maxBaseDim=5` pool-cap pattern. We held every ingredient at identical ranks.

**Reaction:** `SchemeSweep --shape=14x28x29 --base=3x4x6,3x4x7,2x5x7` → **6494**
(same splits, our AlphaEvolve ⟨3,4,6⟩ representative prices richer isolation:
saving 364 vs FMM's 360). Exact-verified, born-pinned. **Beats FMM by 4.**

**Collateral (3 sweep passes over the 187-row WORSE list with
`--base=3x4x6,2x5x7[,3x4x7]`):** 66 further stubs. Pass 1 (recomb-only, 4g):
32 wins — 12 BEAT FMM (⟨15,28,31⟩ 7269→7134 vs 7266 = −132(!), ⟨9,31,32⟩ −63,
⟨9,32,32⟩ −61, ⟨13,25,28⟩ −49, ⟨13,19,28⟩ −34, ⟨13,18,24⟩ −26, ⟨14,18,24⟩ −18,
⟨11,16,26⟩ −15, ⟨21,26,28⟩ −12, ⟨13,18,23⟩ −6, ⟨11,28,29⟩ −3, ⟨11,16,23⟩ −1),
18 exact ties, 2 narrowed. Pass 2 (all strategies, 4g): 12 cascade wins but 92
worker OOMs (17 parallel proj/serendip workers at 4g — heap dump
java_pid62063.hprof). Pass 3 (same, 10g): 22 more wins, 0 errors (recovered the
OOM-aborted shapes, e.g. ⟨28,28,32⟩ 12503, ⟨22,23,23⟩ 6476, ⟨20,25,29⟩ 8024).

**Net cross-check: WORSE 187 → 128 (−59), BETTER 1437 → 1460 (+23).**
Guard row added (⟨14,28,29⟩≤6494). Verification (batch `verifyAuto`-tiered):
**40 exact symbolic proofs + 26 randomised spot-checks, 0 failures** over the 66
stubs (spot-check = dense dim-28+ composites whose exact term-map exceeds
`DEFAULT_MAX_EXACT_TERMS`; a real proof for those needs a streaming verifier).

**Verifier OOM lesson (same session):** unconditional `isExactNonCubic` on dense
dim-28+ composites builds a `HashMap<Long,BigInteger>` term-map of
`Σₖ|Uₖ||Vₖ||Wₖ|` entries → OOM even at 20g (22 GB dumps in target/oom-dumps/).
`Verifier.verifyAuto` (already existed) gates on `estimateExactTerms` and reports
the tier honestly; `VerifyOneScheme` now routes through it.

---

## 2026-07-07 — fmm-gap ⟨28,29,31⟩: CLOSED (tie 13091) via a 3-level upstream chain; 4 upstream WORSE rows closed too

**Target:** `Q⟨28,29,31⟩:m` — ours was 13097, FMM 13091 (+6). All ranks are bounds.
FMM's 13091 is an INDEX-ONLY bound (`rank == naive_rank`, no refs, no tensor; their
detail page documents a stale 13134) — nothing to import, pure closure gap.

**Diagnosis (math-first, from both catalogs' neighbour ranks):**
`FMM 13091 = ⟨28,28,29⟩:11319 + ⟨3,28,29⟩:1772` (p-concat 31=28+3). We lagged on
both ingredients, each with its own cause:
- `⟨28,28,29⟩` ours 11334 vs 11319 = `⟨28,28,28⟩:10535 + 784`: we HOLD 10535 (tie)
  — our concat closure was simply **stale** after 10535 landed. Same for
  `⟨28,28,31⟩` (12256 vs 12241 = 10535+1706).
- `⟨3,28,29⟩` ours 1782 vs 1772 = `252 + ⟨3,25,28⟩:1520`; the real root is
  `⟨3,25,28⟩` ours 1540 vs 1520 = **outer ⟨2,5,7⟩:55 with deficient A-split 3=2+1**
  (35×⟨2,5,4⟩:32 + 20×⟨1,5,4⟩:20) — ANOTHER dim-7 base the `maxBaseDim=5` pool
  cap excludes (the ⟨12,16,27⟩ pattern, third instance).

**Reaction:** 3-pass fixpoint sweep over the 5-shape chain with
`--base=2x5x7,3x4x7 --strategies=kron,concat,recomb,serendip,proj`:
pass 1 → ⟨3,25,28⟩ 1520, ⟨28,28,29⟩ 11319, ⟨28,28,31⟩ 12241; pass 2 → ⟨3,28,29⟩
1772; pass 3 → **⟨28,29,31⟩ 13091 (tie)**. All 5 stubs born-pinned; the two small
ones exact-verified; the three dim-28+ ones replay+spot-check verified (compose-time
write-gate) — exact `isExactNonCubic` OOMs at 16g on r≥11k dim-28 schemes (heap
attempt escalation recorded). Explicitable caveat: ⟨3,28,29⟩/⟨28,29,31⟩ read
`explicitable:false` because their DAG crosses the LEGACY stub `3x28x4-r252-a097e49`
(pre-pinning `ext[...] :: CANONICAL` base ref — replays as best-at-shape, so the
classifier is RIGHT to flag it). They flip to true when the planned legacy re-pin
pass rewrites that ancestor. Related classifier fix landed: parametric constructor
refs (`TA_lita(n=N)`, `DIS09Lemma4(n=N)`) are deterministic formula rebuilds and now
count as explicit (`LineageReplayer.isParametricRef`, shared with
`ComputeExplicitable`; guard `TestComputeExplicitableParametric`) — that alone
flipped ~130 TA-descendant stubs to true (⟨28,28,28⟩=10535 among them). Beating 13091 was checked and priced out: the m-split route
(10535+1706+868=13109) and all other 31-splits land at ≥13091; a strict win needs a
strictly better ⟨28,28,29⟩ or ⟨28,29,30⟩≤12278.

**Repro:** chain shapes in `tmp/chain-shapes.txt` pattern above; guards extended in
`TestSweepSpotsSota.retains_fmm_gap_dim7_base_wins` (⟨3,25,28⟩≤1520,
⟨28,29,31⟩≤13091). Lesson recorded: FMM index-only bounds come from their
combination generator over ingredients — when a gap has `rank==naive_rank`+no refs,
diff the NEIGHBOUR ranks first; the gap is usually upstream, and fixing the
deepest ingredient cascades for free.

---

## 2026-07-07 — fmm-gap ⟨12,16,27⟩: CLOSED and beaten (2988 → 2964, FMM 2984); 11 collateral wins via the dim-7 outer base ⟨3,4,7⟩

**Target:** `Q⟨12,16,27⟩:m` — ours was 2988, FMM 2984 (+4). All ranks below are
**bounds** (recomb closure output), not proven optima.

**Diagnosis (from FMM's 12x16x27.html, one read):** their 2984 = outer ⟨3,4,7⟩:63
block-recursive with uneven p-split 27=6·4+3 — 59×⟨4,4,4⟩:48 + 4×⟨4,4,3⟩:38
leaves (4 of 63 outer products isolate the deficient p-block; 3024−4·10). Our
2988 was the SAME mechanism mirrored — `Recombination(base=⟨4,4,4⟩:48,
allocC=[7,7,7,6], leaves ⟨3,4,7⟩:63/⟨3,4,6⟩:54)` = 3024−4·9 — because
`RecombinationPoolConfig` caps NON-cubic pool bases at `maxBaseDim=5`, so
⟨3,4,7⟩ (dim 7) could never be the OUTER base. Root cause: pool cap, not a
missing scheme or a representative gap.

**Reaction:** `SchemeSweep --mode=materialize --field=Q --shape=12x16x27
--config=thorough --strategies=recomb --base=3x4x7` → WIN 2964
(`Recombination(base=⟨3,4,7⟩:63 [DPS ac0e1ad], allocC=[4,4,4,4,4,4,3],
leaves ⟨4,4,3⟩:38/⟨4,4,4⟩:48)` = 3024−6·10). Exact-verified
(`VerifyOneScheme` → `isExactNonCubic=true`). **Beats FMM by 20.**

**Math-first support census of DPS ⟨3,4,7⟩:63** (one-sided isolation: a product
needs only a ⟨4,4,3⟩ leaf if EITHER its V p-support OR its W p-support sits in
the deficient block — V-side ⇒ zero padded column, W-side ⇒ truncated output):
per-axis one-sided counts n:15(!) m:3 p:7. The engine priced 6 of the 7
p-axis candidates → possible extra −10 (2954) if its isolation rule is the
one-sided OR; worth checking `AllocationOptimizer`'s leaf-pricing rule. The
n-axis count 15 is why the ⟨11,·,·⟩ shapes below fell so hard.

**Collateral (same base over the 203-shape WORSE list; 11 further wins, all
exact-verified):** ⟨16,28,28⟩ 6905→6860 (FMM 6902, −42), ⟨9,28,29⟩ 4221→4161
(FMM 4202, −41), ⟨12,18,20⟩ 2494→2460 (FMM 2493, −33), ⟨11,21,32⟩ 4383→4365
(−17), ⟨12,18,23⟩ 2854→2842 (−11), ⟨9,23,23⟩ 2841→2830 (−5), ⟨11,16,31⟩
3289→3287 (−1), ties with FMM at ⟨11,28,28⟩ 5069, ⟨11,16,28⟩ 2894, ⟨14,16,28⟩
3674; ⟨11,16,30⟩ 3174→3158 narrows (FMM 3126 still ahead, gap 48→32).
Net: 8 WORSE rows flipped to BETTER, 3 tied, 1 narrowed.

**Follow-ups:** (a) the same `--base` probe with other dim-6/7 non-cubic bases
the maxBaseDim=5 pool excludes (⟨3,3,6⟩, ⟨2,5,6⟩, ⟨3,4,8⟩, ⟨4,4,7⟩…) over the
WORSE list; (b) check whether the recomb leaf-pricing implements one-sided
isolation (potential ⟨12,16,27⟩ 2954); (c) consider a pool preset that admits
a curated dim-6/7 base list without the full cost of raising maxBaseDim.

**Repro:** commands above; guard `TestSweepSpotsSota.retains_fmm_gap_dim7_base_wins`
(disk-presence: lookup ⟨12,16,27⟩≤2964, ⟨11,16,28⟩≤2894 — NOT compute-pipeline rows,
since the default pool's maxBaseDim=5 cannot re-derive these; a CSV-row guard was
tried first and burned 70 min before failing at 3004).

**Lineage pinning (2026-07-07, user request):** the original stubs recorded the
outer base as the `--base` pool label (`base<3x4x7>=63`) — a best-at-shape ref
that re-resolves against future catalogs (`explicitable:false`). Fixed
`SchemeSweep.userBasePool` to stamp `originLineage` = `{storedShape}@{contentHash}`
(OrientAs-wrapped for non-native orientations), deleted + re-swept the 12 stubs
(identical ranks, base now `3x4x7@ac0e1ad…`), all re-verified exact and now
`explicitable:true`. Guard: `TestSchemeSweepUserBasePool`. **Going-forward pinning is now total**
(user 2026-07-07): rootPool's six synthetic AxisSplit/Naive entries pin as
`naive-NxMxP` (content-hash-identical to the naive schemes —
`TestRootPoolOrigins.axis_split_bases_are_content_identical_to_naive`), and
`RecursiveMaterialiser`'s null-origin fallback now pins by content hash when the
base resolves in the catalog, else WARNs loudly before writing a label ref.
Guards: `TestRootPoolOrigins` (root/default/thorough pools all carry origins).
The ~2.1k OLDER derived stubs with unpinned refs (bare `base=NxMxP`,
`ext[...] :: CANONICAL` labels) remain; agreed plan: a later **re-sweep driven by
the current partial lineage** — replay each stub's recorded base-shape/alloc/leaf
structure as search hints so the full pinned lineage regenerates quickly rather
than from scratch.

---

## 2026-07-06 — fmm-react ⟨20,24,25⟩: gap +8 NOT closed (FMM 6466 vs our 6474); 28 collateral serendipitous wins landed in the 20–25 band

**Target:** `Q⟨20,24,25⟩:m` — ours 6474 (stub `Serendipitous(⟨5,12,5⟩ ⊗ˢ ⟨4,2,5⟩)` = 204·32−54),
FMM 6466 (+8). Closure-gap diagnosis (max-dim 25, cross-check "small gaps = 0").

**What was tried (all Q, all failed to close the target):**
- `ProjectFmmGaps --shapes=20x24x25 --passes=2` → 0 wins (7 parent replays margin-pruned;
  all FMM parents ≥ 6924, projection hopeless — math confirms).
- `SchemeSweep --mode=materialize --field=Q --shape=20x24x25 --config=thorough
  --strategies=recomb,serendip,proj` → no improvement (18-base recomb pool; predictor
  correctly prunes ⟨2,2,2⟩: brute-force over ALL uneven Strassen block splits with
  today's catalog ranks bottoms out at 6750 = 5·⟨10,12,13⟩:990 + 2·⟨10,12,12⟩:900).
- `SerendipitousSweep --baseCap=12 --secondCap=5 --targetMin=20 --targetCap=25 --field=Q`
  → no ⟨20,24,25⟩ win, but 35 collateral target-wins (see below).
- `VerifyGLCandidate 20 24 25` (GL×GL×GL frontier of ⟨2,2,2⟩, 36 members) → tie:
  best eval 6528, materialises back to catalog 6474. No GL win here.

**FMM's 6466 construction (site came back mid-session — from their 20x24x25.html):**
`⟨20,24,25⟩:6466 = (⟨5,12,5⟩:204 − 42) ⊗ ⟨4,2,5⟩:32 + ⟨16,2,5⟩:126 + 2·⟨4,6,5⟩:90
+ 16·⟨4,4,5⟩:61` — the SAME serendipitous pairing as ours, cut deeper. In bud
language the removed 42 terms form `1×⟨4,1,1⟩ + 2×⟨1,3,1⟩ + 16×⟨1,2,1⟩` sub-matmuls
= bud profile **16×W2 + 2×W3 + 1×V4** (σ = 16·3 + 2·6 + 1·2 = 62; ours: 18×W2 → 54;
the σ-arithmetic reproduces both numbers exactly).

**Key finding — it's a bud-REPRESENTATIVE gap, not an engine gap:** we imported FMM's
own published ⟨5,5,12⟩:204 tensor (their page attributes it to Perminov FMM GitHub,
Feb 2026; verified `isExactNonCubic=true`; **0 shared triads** with our
`perminov_ZT-61a6cb7` 204) into
`bud-bases/section12/fmm-lille_5x5x12_r204_a2326.json`. Direction-class scans (Java
`SerendipitousBudProduct.findBuds`/`independentClassSizes` under ALL_ORDERINGS +
an independent Python up-to-scaling scan) agree: **both** 204 representatives carry
only 18×(size-2) classes — no size-3, no size-4, no second axis. Direction classes
are axis-permutation-invariant, so the representative FMM actually composed with is a
THIRD rank-204 scheme (flip/GL-moved) whose bud profile is richer than either
published one. Probe: `docs.explore.ProbeFmmSerendip202425`.

**Flip-walk follow-up (RAN, no win at this budget):**
`docs.explore.ProbeFlipWalkSerendip202425 30000 4` — flip-graph walk from both 204
seeds under the DIRECT objective `FlipObjectives.serendipitous(Q, 4,2,5)`. All
8 walks (2 seeds × 4 rngs × 30k steps, ~30 s each) plateaued at predicted 6474.
Bounded negative — NOT proof the richer representative is flip-unreachable. Next
escalations, in order: (a) longer walks / more seeds / larger `maxRankAbove`,
(b) seed from OTHER 204-adjacent schemes (5x5x12 at 207/208/213…, walking rank back
down), (c) parse FMM's composed `20x24x25_tensor.mpl` (site is back up) and factor
the base representative out of it directly — the 42 cut terms and their patch
blocks are identifiable in the composed scheme's term structure.
The imported base `bud-bases/section12/fmm-lille_5x5x12_r204_a2326.json` is retained
as the seed for these escalations (PURGE_REFCOUNT_POLICY: delete it if the follow-up
is abandoned).

**Collateral wins (LANDED, bounds not optima):** 28 serendipitous stubs persisted via
`MaterialiseSerendipitousWins --shapes=…` (driver extended this session to accept
`--shapes=NxMxP,…`), sections 20/21/24/25, e.g. `Q⟨20,20,25⟩:m` 5614→5611,
`Q⟨21,25,25⟩:m` 7304→7293, `Q⟨14,25,25⟩:m` 5039→5030, `Q⟨9,9,21⟩:m` 1060→1058,
`Q⟨4,20,14⟩:m` 755→736 (biggest single drop, −19). All spot-check-verified at compose
time; 7 further sweep hits resolved through orientation twins landed in the same run.
`TestSweepSpotsSota` green after landing.

**Repro:** commands above; the win list is the `< sota` lines of the
SerendipitousSweep run (baseCap=12 secondCap=5 targetMin=20 targetCap=25, field=Q).

**Stamper side-effect cleanup (2026-07-07):** the tree-wide stamper pass run for
the FMM import's backfill exposed two silent bugs, both fixed + regression-tested:
(1) `StampSource` on NEW-convention filenames (shape-first) fell back to the whole
stem and stamped 7,681 filename echoes (`"source": "10x3x3-R69-Derived-79abb2c"`) —
reverted by `docs.migrate.FixFilenameEchoSources`, stamper now skips shape-first
stems (`TestStampSource`); (2) `NarrowFields.expandFromBase` returned a hardcoded
chain for the strongest base tag, stripping VERIFIED F3 from 3,730 Q schemes with
power-of-2 denominators above the verify cap — files restored from HEAD, expansion
is now a union (`TestNarrowFields`). Kept: the legitimate widening of 88 Z-only
Perminov CI imports to `[F2,F3,Z,Q,R,C]`. Catalog count moved 9283→9233 (+28
session stubs, −78 `known/` duplicates that the widening exposed to the by-design
per-(format,field,comm,rank) dedup election — files remain on disk as witnesses).

**Status:** target OPEN (+8). Blocked on FMM site availability for their construction;
our in-scope closure mechanisms (projection / recomb / serendip / GL-frontier) are
exhausted at 6474.

---

## 2026-06-25 — VERIFIED R⟨17,19,20⟩=3780 (−17 vs catalog 3797) via GL-orbit frontier of ⟨2,2,2⟩

**Claim:** `R⟨17,19,20⟩ = 3780`, beating the catalog's 3797 by 17. The win comes
from a GL-orbit frontier member of the ⟨2,2,2⟩ Strassen base — a change-of-basis
support the live pool's AXIS_FLIP expansion cannot reach (it only does
row-reversals; this needs full `GL×GL×GL`).

**Status: REGISTERED to the catalog (2026-06-25).** Written as a lineage stub
`schemes/derived/section20/17x19x20-r3780-derived-1dd1c50.json` via
`docs.migrate.RegisterGLWin 17 19 20` (lineage = Recombination of base
`2x2x2@511df05…` = **Winograd_1971**, allocA=[9,8] allocB=[10,9] allocC=[12,8] +
Project leaves). `docs/catalog.json` regenerated (`findRank(17,19,20)=3780`).
Independently replay-verified from the COMMITTED catalog alone (no injected base):
`VerifyOneScheme` → replayed r=3780 (151 ms), `isExactNonCubic = true` (154 s).
`TestSweepSpotsSota` 16/16 still green. The superseded `…-r3797-…b4e28ce.json`
remains on disk (harmless — `findRank` picks the 3780 min). Manifest flags it
`explicitable=False` (conservative for >16 stubs) but the disk-only replay proves
it explicitable. Field-stamp note: tagged `[R,C]` (old 3797 was `[Q,R,C]`) — may
under-claim Q; not yet audited.

**Breakdown** (single-level ⟨2,2,2⟩ GL-member recombination, all sub-blocks
`find()`-buildable, maxDim ≤ 12 so not stub-blind):
`2×⟨9,10,12⟩=668 + 1×⟨9,9,8⟩=430 + 1×⟨8,9,12⟩=560 + 2×⟨8,10,8⟩=427 + 1×⟨9,9,12⟩=600 = 3780`.

**Repro:**
- `mvn -q -o -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.verify.VerifyGLCandidate -Dexec.args="17 19 20"`
- Discovered by `eu.solven.matmul.docs.explore.GLLargeSweep` (band 17–32, ⟨2,2,2⟩
  base, full GL frontier vs `FieldAwareLookup.findRank`): 2 strict wins + 480
  shapes where the single ⟨2,2,2⟩ GL frontier re-derives catalog SOTA exactly.

**Second candidate (UNVERIFIED): R⟨18,30,31⟩=8970 (−18 vs 8988).** Rests on
sub-blocks ⟨9,18,18⟩/⟨9,18,13⟩/⟨9,12,18⟩ that are `find()=none` — but maxDim=18>16,
so this is the `find` vs `findRank` stub-blindness case, NOT proven phantom.
Needs a materialiser build to confirm or reject. `VerifyGLCandidate 18 30 31`
materialise was too slow to finish in-session (huge >16 sub-blocks).

**To register R⟨17,19,20⟩=3780:** re-run the materialise with `writeNewSchemes=true`
into `schemes/derived/section20/`, lineage = recombination of the GL ⟨2,2,2⟩
member; then regenerate `docs/catalog.json` and run `TestSweepSpotsSota`.
Deferred pending user OK (autonomous session avoided catalog mutation).

**Methodology note — three bugs this run surfaced (see also memory):**
1. `SchemeIO.read` didn't dispatch the Perminov "reduced"/compact format → 168
   files (incl. NC integer bases ⟨4,4,4⟩=49, ⟨2,4,4⟩=26) unloadable. **FIXED**
   (`read(JsonNode)` now branches `isReduced`), +3 regression tests.
2. `Recombination.catalogResolver` returns CUBIC `a·b·c` for shapes outside its
   small hardcoded `KnownAlgorithmCatalog` (⟨9,9,9⟩→729 not 486) — invalidated the
   absolute numbers of every AllocationOptimizer probe (relative same-sota
   comparisons survive). Use `FieldAwareLookup.findRank`-backed sota instead.
3. `RecombinationMultisetOrbit.enumerateSampled` over-reports the frontier (226 vs
   exact 170 for ⟨2,3,3⟩) — canonicalization mismatch; harmless for win-finding
   (duplicates, build-verified) but wrong for any exhaustiveness claim. NOT fixed.

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
