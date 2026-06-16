# Sweep configuration ladder

The catalog-closure search has several **independent richness axes**. Each
"richer" setting is (almost always) a strict **superset** — more coverage at
more cost — so configs form a **priority ladder**: run the cheap / high-yield
tiers first to *match FMM fast*, and the expensive / thorough tier later to
*spot improvements beyond FMM*.

The **incumbent-bound prune** (improve mode) makes a richer tier cheap on the
shapes a cheaper tier already solved (it proves "no allocation beats what we
have" fast), so tiers **accumulate** rather than duplicate work.

## Axes (flag → meaning → cost)

| Axis | Flag | Cheap ⊂ … ⊂ Rich | Catches when richer |
| --- | --- | --- | --- |
| **Base pool** | `--config` | `simple` (cubicOnly) ⊂ `rectangular` ⊂ `includeDerived` ⊂ `thorough` | rectangular roots ⟨2,3,3⟩/⟨2,2,3⟩/⟨3,3,4⟩/⟨3,3,6⟩ (FMM's 2·3·3-style splits), then derived bases |
| **Allocation** | `--balancedOnly` (+`--maxImbalance`) | `true` ⊂ `false` | unbalanced block sizes (e.g. the (6,26) p-concat) |
| **Operator** | driver | compose (`SchemeSweep`) ⟂ project (`ProjectFmmGaps`) | downward projection catches FMM's projected schemes (complementary, not superset) |
| **Orbit** | `--config` orbitMode | CANONICAL ⊂ AXIS_FLIP ⊂ PERMUTATION_BOUNDED | flipped / permuted base orientations |
| **B&B budget** | `--stagnation` / `--maxNodes` | low (~10k) ⊂ high (100k / ∞) | exhaustive optimum vs anytime-near-optimal |
| **Mode** | `--improve=true` | (needed for re-sweep + incumbent-bound prune) | keeps only strict improvements |

## Priority ladder (tiers are DISJOINT — each does only NEW work)

Use `--baseShape=cubic|rectangular` to run **strictly** that base-shape class, so a
later tier never re-walks the bases an earlier tier already covered. The
incumbent-bound prune then means each tier only *writes* where its new bases beat
the running incumbent.

| Tier | Bases (`--baseShape`) | Alloc | Operator | Orbit | Budget | Catches | Cost |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **T1** | cubic only | unbalanced | compose | CANONICAL | 10k | cubic unbalanced splits | fast |
| **T2 — match-FMM** | **rectangular only** | **unbalanced** | compose | CANONICAL | 10k | **FMM 2·3·3-style rectangular splits** | medium |
| **T3** | (n/a) | — | **projection** | — | — | FMM's projection-derived schemes (complementary operator) | +medium |
| **T4 — full (later)** | all + derived | unbalanced | compose + project | AXIS_FLIP / PERM | high | flips, derived bases, exhaustive B&B → *new* improvements | slow |

## Commands (over the WORSE shapes, gap-ordered — `research/worse-by-gap.txt`)

```bash
# T1 — cubic bases only (fast baseline):
SchemeSweep --mode=materialize --improve=true --config=simple \
    --balancedOnly=false --maxImbalance=64 --stagnation=10000 --shape=<WORSE>

# T2 — match FMM: STRICTLY rectangular bases (disjoint from T1) → ⟨2,3,3⟩-style splits.
SchemeSweep --mode=materialize --improve=true --config=rectangular --baseShape=rectangular \
    --balancedOnly=false --maxImbalance=64 --stagnation=10000 --shape=<WORSE>

# T3 — projection complement (downward operator, parent-centric scatter):
ProjectFmmGaps --passes=2

# T4 — full, slow, may BEAT FMM (run once, overnight):
SchemeSweep --mode=materialize --improve=true --config=thorough \
    --balancedOnly=false --maxImbalance=128 --maxNodes=50000000 --shape=<all>
ProjectFmmGaps --passes=3
```

## Now vs later
- **Now (match FMM):** T2 then T3. Fast, closes the bulk of the gap with FMM's own recipes.
- **Later (beat FMM):** T4 — exhaustive budget + derived bases + orbit; slow, but the only tier that can find ranks *below* FMM by construction rather than by chance.

## Notes
- All sweeps are **anytime + incremental-persist** (crash-safe; interrupt keeps wins).
- The lineage write-guard + task-#91 fixes mean every win is bit-exactly replayable
  (no phantoms), so tiers compose safely.
- Verified: cubic-only ⟨22,27,29⟩=9447 vs `rectangular` =9392 (FMM 9369) via
  `base=2x3x3` — confirms T2's value over T1.
