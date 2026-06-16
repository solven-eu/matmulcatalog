# Projection margin vs rank — the tradeoff, formalized

*(2026-06-12; the projection sibling of `SERENDIPITY_RANK_TRADEOFF.md`; see
also `BUD_STRUCTURE_THEORY.md` for the bud-count vs σ-value distinction and the
flip-graph move lifecycle. Status: the statements under "Exchange rate" and
"Free margin" are theorems; the scan numbers are catalog-dated.)*

## Catalog home: `schemes/margin-bases/`

Projection-rich schemes that are rank-worse/rank-tied for their own shape
live in `schemes/margin-bases/sectionN/` (sibling of `bud-bases/`, same
rationale: heuristic finds, recursed by every reader). Stamp:
`projection_margins: [μn, μm, μp]` — **target-free**, because margin is
INTRINSIC (supports only; no partner pricing, unlike serendipity's σ which
goes through R(inner-enlarged)). One scheme serves every target below it;
multi-drop margins recompute from the factors
(`ProjectionSearch.projectedRank`); the win question (proj vs R_cat(T)) is
catalog-relative and evaluated at query time, never stamped. Writer:
`FlipGraphSweep --objective=project|margin --write`.

## The metric

A base `B` (shape ⟨n,m,p⟩, rank r) used as a **projecting parent** for target
`T = ⟨n−a, m−b, p−c⟩` yields a T-scheme of rank

```
proj(B→T) = r − μ(B, drop)   (exact: min survivors over keep-combos — DCE)
```

A product dies under a drop iff its support on the dropped axis, in either
carrying matrix (n→{U,W}, m→{U,V}, p→{V,W}), collapses to dropped indices.
Code: `ProjectionSearch.projectedRank` (exact, target-priced — the real
currency), `ProjectionSearch.axisMargins` (per-axis triple `[μ_n,μ_m,μ_p]`,
best single drop each; the stamped scalar `projection_margin` is their max),
`FlipObjectives.projectedTo` (direct walk objective,
`FlipGraphSweep --objective=project --target=NxMxP`).

**Second-order margin (currently uncounted):** after dropping, surviving
products can become proportional (post-projection bud pairs) or zero-padded —
`survivorCount` counts deaths only. A merge-aware projected rank is a known
tightening, not yet implemented.

## Exchange rate is exactly 1 — no relaxation window

`proj = r − μ`, so one point of base rank trades against one point of margin:
rank+δ pays iff it buys `Δμ > δ`. The ceiling `μ ≤ r − R_lb(T)` also grows
with slope exactly 1 in r. Consequence: **unlike serendipity (rate R(inner),
ceiling slope σ-density ≪ R(inner) → finite window δ_max), projection has NO
a-priori limit on profitable rank relaxation.** Whether rank+5 with margin+6
exists is purely a realizability question — this is the regime where
"slightly higher rank, much better margin" walks make structural sense.

## Free margin — what a projection win actually requires

Margin is trivially manufacturable: glue the best T-scheme (rank R_T) with a
naive slice block (e.g. the dropped C-row computed as ⟨1,m,p⟩, mp products) —
a valid ⟨n,m,p⟩ base of rank `R_T + mp` and margin exactly `mp`, projecting
back to `R_T` on the nose. So "high μ" alone is worthless; the base must
encode T-knowledge **beyond the catalog**:

```
win ⟺ proj(B→T) < R_cat(T) ⟺ μ > r − R_cat(T)
```

and `proj(B→T) ≥ R(T)` always (the projection IS a T-scheme), so a target at
a **proven-optimal** catalog rank is permanently unwinnable — the gap may be
open but the door is closed. The valuable Pareto object is therefore
`(r ≈ SOTA(⟨n,m,p⟩), μ large)`: near-optimal big-shape schemes that also
project hard. This is the (rank, μ) Pareto axis the paper's §projmargin
gestures at, now priced per target.

## Scan findings (2026-06-12, dims ≤ 7, 1-drop, catalog-best parents)

`docs.explore.ScanProjectionFrontier`: **522 open gaps vs 18 matches**, gap
almost always exactly 1 (e.g. Laderman ⟨3,3,3⟩ → ⟨2,3,3⟩ projects to 16 vs
catalog 15; ⟨2,4,3⟩ → ⟨2,3,3⟩ same). Reading: at small dims the catalog
targets sit at/near proven optima (HK71 band), so those gaps are unwinnable
artifacts, NOT search failures. The actionable estate is where
`R_cat(T) − R_lb(T)` is large — the 17–32 band, where projection already
drives the catalog (FMM parents) and where a margin-walk on an imported
near-SOTA parent could realistically buy `Δμ > δ`. Caveat: the scan prices
the catalog-BEST representative per parent shape; margin-richer equal-rank
representatives (margin-bases, mirroring bud-bases/) would shift it.

## Margin potential μ_max(δ) — the ceiling as a function of rank slack

*(2026-06-12, paper-feed section. δ = r − R(⟨n,m,p⟩) is the base's rank slack;
ΔR = R(⟨n,m,p⟩) − R(T) the shape step's true rank increment. Where R is
unknown, catalog values give the honest bound-version of every statement.)*

**Ceiling theorem.** For ANY base at rank r = R(big) + δ:

```
μ ≤ r − R(T) = ΔR + δ          (the projected scheme IS a T-scheme)
```

So the potential at **tight rank is ΔR, not 0** — tightness does not kill
margin — and each unit of slack raises the ceiling by exactly 1.

**Nesting question (achievability at δ = 0).** μ = ΔR at tight rank ⟺ the
survivor set of some drop is an *optimal* T-scheme sitting inside an
*optimal* big-scheme — "optimal schemes nest". Worked invariants (exact,
catalog 2026-06-12):

| base | μ-triple | proj → target | ceiling | verdict |
|---|---|---|---|---|
| Strassen ⟨2,2,2⟩ r=7 (tight, de Groote-unique) | [3,3,3] | ⟨1,2,2⟩: 4 = R | 7−4=3 | **perfect nesting, every axis** — an invariant of THE rank-7 scheme |
| Laderman ⟨3,3,3⟩ r=23 | [7,6,6] | ⟨2,3,3⟩: 16 | 23−15=8 | nesting fails by 1 (this representative) |
| ⟨3,3,8⟩ r=55 (closure, rational) | [4,4,6] | ⟨3,3,7⟩: 49 = R_cat | 55−49=6 | **catalog-ceiling-tight** |

**Achievability at slack (block-diagonal).** Axis slices split cleanly
(C-columns partition with no cross terms), so `T-scheme ⊕ ⟨n,m,d⟩-scheme` is
a valid base of rank R(T)+R(slice) with proj = R(T) exactly. Hence the
ceiling is ALWAYS achieved from

```
δ_bd = R(T) + R(slice) − R(big)        (⟨3,3,8⟩: 49+9−55 = 3)
```

onwards. Between δ = 0 and δ_bd realizability is per-shape: the rank-55
⟨3,3,8⟩ scheme reaches proj = R_cat(T) at δ = 0, i.e. it "shares" 3 of the 9
products a concat would spend on the extra slice — non-trivial nested
structure beating block-diagonal economics.

**Win theorem (the punchline).** For any δ:

```
∃ base with proj(B→T) < R_cat(T)   ⟺   R(T) < R_cat(T)
```

(⇐ take any better T-scheme, block-diagonal it up; ⇒ the projection itself
is the better T-scheme.) **Projection-margin search is target-rank search in
a bigger parametrization** — never a way around the target's hardness, only
a different (sometimes structurally richer) place to look for it.
Trichotomy per (base, target): *catch-up* (proj > R_cat(T): real search room
with no new mathematics — close the gap to the tie; the 522 gap-1 cases),
*tie* (proj = R_cat(T): any further gain ⟺ improving R(T) itself), and
proj < R_cat(T) cannot persist (it would BE the new catalog value). This is
why the ⟨3,3,8⟩ walks froze at 49 in all three basins: the seed is
tie-tight, so the walk was implicitly asked to prove R⟨3,3,7⟩ ≤ 48.

**∃ vs ∀ duality with lower bounds.** Our μ is an ∃-scheme quantity (this
base projects well). Substitution-method LB proofs are the ∀-scheme dual:
"in EVERY rank-r scheme, ≥ k products must die on this slice" is exactly a
universal margin statement, and stepwise LB increments R(big) ≥ R(T) + L are
∀-margin theorems. The two meet at tight rank: L ≤ ΔR = the ∃-ceiling. Same
support structure, two quantifiers — a paper-ready framing of why projection
margins and substitution LBs are the same combinatorics.

**Second-order margin.** The computed μ counts DCE deaths only; a drop can
also make surviving products proportional (slicing CREATES bud pairs) or
zero. True projected rank ≤ survivor count, so every μ here is a lower bound
on the real margin — and merge-aware projection is the natural bridge from
the projection study to the serendipity study (drops manufacture buds).

## Walkability

Flips move supports (a flip u_i += s·u_j unions/cancels row support), so μ is
mobile at fixed rank in principle; splits buy rank for structure at the 1:1
rate. Empirically (2026-06-12): random ternary walks froze at the seed's
projected cost in ALL basins tried — block-diagonal r=58 (proj 49), Perminov
ZT r=56 (proj 51), and the scaled-lift rank-55 rational basin (proj 49; the
lift c·T trick makes rational bases walkable: clear per-product denominators,
equalize to a common c — flips preserve the sum and the projection objective
is support-based, hence scale-invariant). Structural reason: flips UNION
supports except at rare cancellations, so random flips almost never buy
margin. The win theorem above shows the freeze at tie-tight seeds is forced;
catch-up cases (proj > R_cat) are where directed support-concentrating moves
could genuinely pay. Harness: `--objective=project --target=…
--rank-above=δ --split-prob>0 --no-full-reduce`.
