# Purge & re-derivation policy — SOTA loss ≠ deletable; refcounts are the guard

*(2026-07-05; operational sibling of `PROJECTION_MARGIN_TRADEOFF.md` — that note
prices the (rank, μ) Pareto, this one says what catalog hygiene must do about it.
Status: policy + gap list; the motivating numbers are exact and catalog-dated.)*

## The problem

A scheme routinely loses rank-SOTA at its own shape while remaining
**load-bearing elsewhere** — as a projection parent, recombination base, Kron
factor, or concat block. Purge/rederivation passes that key on
"rank-dominated" alone therefore delete (or silently re-target) the very
schemes the downstream catalog stands on. The catalog must (a) count inbound
references per scheme, (b) surface the count honestly, and (c) recategorize
a dethroned scheme from *SOTA* to *base* rather than purging it.

## Motivating case (exact, 2026-07-05): LITA dethrones DIS09 at ⟨30,30,30⟩

| ⟨30,30,30⟩ cube | rank | proj. margin μ | proj → ⟨29,30,30⟩ |
|---|---|---|---|
| LITA (KGP 2026) | **12672** (cube SOTA) | 30 | 12642 |
| DIS09 Pan-TA | 12710 | **122** | **12588** = catalog & FMM SOTA |

LITA's local re-pairings save 38 cube products by merging exactly the
axis-localized correction products a coordinate drop would kill — so the new
cube SOTA is the *worse* projection parent by 54. "Regenerate rectangulars
from the new SOTA cube" would regress ⟨29,30,30⟩ (and its projection
descendants ⟨29,29,30⟩, …). Exchange-rate view
(`PROJECTION_MARGIN_TRADEOFF.md`): DIS09 pays δ=38 rank for Δμ=92 margin —
profitable, and with no a-priori ceiling on such trades, **this pattern will
recur every time a tighter cube family lands.**

And the current catalog.json (2026-07-05) mislabels it exactly as feared:

```
30x30x30-r12710-dis09  rank=12710  projection_margin=122  used_as_base: ABSENT
30x30x30-r12672 LITA   rank=12672  (rank-best)            used_as_base: ABSENT
29x30x30-r12588        lineage: Project(Atom "DIS09Lemma4(n=30)")  used_as_base: {by_projection: 3}
```

DIS09(30) reads as *rank-dominated + unreferenced* — a textbook purge
candidate — while actually parenting the 12588 SOTA and, transitively, its
three projection children. Two mechanisms hide the truth (see Footguns).

## What already exists (inventory)

- **`used_as_base`** — `GenerateCatalogManifest.stampBaseUsage` (+
  `catalog/BaseUsageStats`) stamps each catalog.json entry with inbound
  usage by op: `by_projection / by_recombination / by_kronecker / by_concat /
  by_sum / by_disjoint_sum / by_serendipitous / by_augment` + `total`.
  Attribution: pinned `shape@hash7` refs credit the exact scheme; bare
  `shape` refs credit the current rank-best entry of that shape.
- **`projection_margin`** (max of `axisMargins`) + `schemes/margin-bases/`
  as the home for margin-rich rank-worse finds.
- **`CatalogKeepClosure`** (docs.verify, advisory): keep = reachable from a
  keep-root (rank-best per shape) via lineage edges; everything else
  "DELETABLE".
- **`DedupDerivedSchemes`** (docs.migrate): deletes a derived scheme only if
  every field is strictly better-covered by a sibling AND its hash is pinned
  nowhere — never dangles a pin, never loses field coverage.
- **`RetireDerivableImports`** (docs.migrate): re-pins all inbound
  `shape@oldHash7` refs to the derived twin *before* deleting, fails loud on
  survivors. The template for any actual deletion.
- **`AuditLineageRefs`** (docs.verify) classifies refs (PINNED_OK /
  PINNED_DANGLING / BARE / PARAMETRIC / …); `RepinDanglingLineageRefs`
  repairs. Neither builds a per-scheme inbound count — `used_as_base` is the
  only reverse index.
- **`LineageReplayer.resolveLeaf`**: a dangling pin THROWS (since
  2026-07-05) — historically it WARNed and fell back to bare-shape
  (catalog-best) resolution ("phantom replay"; task #91 precedent:
  ⟨17,22,29⟩ 6129 → 6138). Guard:
  `TestReplayLeafFallback.dangling_pinned_ref_throws_instead_of_shape_best_substitution`.
  The legacy shape-extraction fallback for UNPINNED source-prefixed /
  canonical-key refs (no hash) remains — those refs never recorded exact
  content, so shape-best is their defined meaning.

## Footguns (why `used_as_base == 0` does NOT mean purgeable)

1. **Parametric parents are invisible.** `BaseUsageStats.baseKey` returns
   null for `DIS09Lemma4(n=30)` / `Strassen<2,2,2>=7` / `naive-*` refs — only
   `shape@hash7` and bare-shape refs are counted. The DIS09(30) *file* gets
   no credit for the 12588 child because the child pins the parametric
   constructor, not the file. (Replay is safe — parametric refs rebuild from
   code — but the catalog's refcount, and any purge pass reading it, is
   blind to the dependency.)
2. **Bare-ref credit migrates on SOTA change.** Bare-shape credit goes to
   the *current* rank-best (`stampBaseUsage.bestByShape`). The moment a new
   SOTA lands at a shape, the dethroned scheme's bare-ref credit silently
   transfers to the newcomer — its `used_as_base` drops exactly when its
   status question is being asked, and the newcomer inherits credit for
   derivations it never produced (whose replay through it may yield worse
   children — the LITA case, μ 30 vs 122).
3. **Phantom replay masks breakage** — CLOSED 2026-07-05: deleting a pinned
   parent now hard-fails children's rederivation (`LineageReplayer` throws on
   a dangling pin instead of silently re-resolving to shape-best). A purge
   "verified by successful replay" is now meaningful for pinned refs — but
   still proves nothing for bare/parametric refs (footguns 1–2).
4. **Keep-closure roots are rank-only.** `CatalogKeepClosure` roots at
   rank-best per shape; a (rank, μ)-Pareto-nondominated scheme that nothing
   pins (because its children reference it parametrically, or bare) is
   classified DELETABLE — DIS09(30) today.

## Policy

1. **Never purge on rank-dominance alone.** A deletion candidate must clear
   ALL of: (a) hash pinned nowhere (`AuditLineageRefs`); (b) `used_as_base`
   total 0 *including* bare-credit it would receive as shape-best and any
   parametric twin-ship; (c) not on the per-(shape, field) (rank, μ) Pareto
   frontier — some kept sibling has ≤ rank AND ≥ margin AND covers its
   fields; (d) outside the keep-closure.
2. **Recategorize, don't delete.** A scheme that loses shape-SOTA but keeps
   inbound refs or a nondominated margin changes *role* — SOTA → projection /
   recombination base. It stays in the catalog and in parent pools
   (`margin-bases/` is the home for new margin-rich finds; existing entries
   keep their place).
3. **Rederivation after a SOTA change re-prices, never swaps.** Do not
   mass-repin children (bare or pinned) onto the new SOTA parent: replay each
   child against both parents and keep the better result. A tighter cube is
   not a better base (exchange rate 1: δ rank buys only ≤ δ projected rank,
   and may cost far more margin).
4. **Any actual deletion follows the `RetireDerivableImports` discipline**:
   repin-or-verify every inbound ref first, delete, then fail loud if any
   old pin survives.

## Gaps to close (actionable)

- **G1 — count parametric refs.** Extend `BaseUsageStats.baseKey` (or a
  parallel key space) to attribute `DIS09Lemma4(n=…)`-style refs, and stamp
  them onto the matching catalog scheme when a file-twin exists (match by
  shape + rank + source). Today's catalog shows the failure live.
- **G2 — split pinned vs bare credit.** Emit
  `used_as_base: {…, pinned_total, bare_total}` so migrating credit is
  visible; a purge pass may trust `pinned_total`, never the merged total.
- **G3 — per-entry `role` field** in catalog.json, computed at manifest
  time: `sota` (rank-best for shape+field) / `base` (inbound refs > 0, or
  pinned, or parametric-twin) / `margin-pareto` ((rank, μ)-nondominated for
  its shape+field) / `historical` (everything else — the only purge-eligible
  class, and even then only via Policy §4). SPA badge to match.
- **G4 — keep-closure roots = Pareto, not rank-best.** `CatalogKeepClosure`
  should root at the per-(shape, field) (rank, μ) frontier plus everything
  G1/G2 counts, so "DELETABLE" stops flagging load-bearing parents.
- **G5 — strict-pin replay** — DONE 2026-07-05, and stricter than first
  scoped: hard-fail on dangling pins is the DEFAULT, not an opt-in mode
  (`LineageReplayer.resolveLeaf` throws; guard test in
  `TestReplayLeafFallback`). Corollary: `AuditLineageRefs` PINNED_DANGLING
  must stay 0 — a dangling pin is no longer a degradation but a replay
  breaker, so `RepinDanglingLineageRefs` must run before any catalog
  mutation lands.
- **G6 — regression guard** (per the CLAUDE.md silent-regression rule): a
  `TestSweepSpotsSota`-style row asserting ⟨29,30,30⟩ materialises ≤ 12588,
  so a future "regenerate from the LITA cube" sweep cannot silently regress
  the projection family.
- **G7 — legacy stem refs still resolve shape-best.** `AuditLineageRefs`
  snapshot 2026-07-05: PINNED_OK ×18532, PINNED_OK_STUB ×453,
  PINNED_DANGLING ×5, BARE ×17, PARAMETRIC ×96, NAIVE ×7757, and
  **UNRESOLVABLE ×16571 across 5866 files** — legacy source-prefixed /
  canonical-stem / `ext[…]` refs that replay via loose shape extraction to
  the *current* shape-best (deliberate: those refs never recorded exact
  content, and the replayer's javadoc'd fallback + `TestReplayLeafFallback`
  pin that contract). Extending throw-on-unclear to this class first
  requires a repin sweep (stems carrying a `-hash7` token can be upgraded to
  `shape@hash7` mechanically; the rest need shape-best-at-repin-time
  pinning) — until then, these 5866 files replay content that can drift
  with the catalog.
- **G8 — audit/repin pin handling is field-blind.** The 2026-07-05
  PINNED_DANGLING ×5 decomposed into: ×1 real (an unstamped stub —
  `31x31x31-r14878-dis09_Q` had no `"hash"` field, so `30x31x31-r14573`'s
  pin dangled) and ×4 false positives — F₂-only parents (`4x4x4@258e5b7`,
  `2x2x2@89207d8` = alphatensor-F2 files) invisible to the audit's
  `FieldAwareLookup("C")`, which does NOT index the F₂ universe.
  `AuditLineageRefs.classify` (and `RepinDanglingLineageRefs.isDangling`,
  which safely QUARANTINEs but shouldn't even flag them) should retry pins
  against an F₂/F₃ lookup before declaring PINNED_DANGLING.
  **Real case repaired 2026-07-05** via the chain that is now the template
  for hash-regime drift: `StampStubHashes` (new docs.migrate driver — replays
  each unstamped lineage-only stub, requires replay rank == declared rank,
  stamps `contentHash`; 31 stamped: the dis09_Q cubes n=4…32 plus two
  derived stubs, all named under the pre-`RehashRationalsComplex` legacy
  hash) → rename the 31 stale filename tokens to the new hash7 (filenames
  are labels, but `GenerateCatalogManifest.hash7Of` keys `used_as_base`
  credit off the token) → `RepinDanglingLineageRefs --apply --max-dim=31`
  (REPINNED_BITEXACT: the child replays to its stored hash with the
  re-stamped parent) → re-audit. Residual PINNED_DANGLING = the ×4 F₂ false
  positives. Un-replayable leftovers: 14 `rosowski_2019_thm2` stubs whose
  `RosowskiTheorem2(l,n,p)` parametric leaf is not wired into
  `LineageReplayer.resolveParametric` — wire it or accept them as
  non-replayable formula stubs.
