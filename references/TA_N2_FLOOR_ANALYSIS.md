# The `15/4·N²` floor of trilinear aggregation — a four-lever analysis

*2026-08-12. Question posed by the user: given the recent SZ (2025) and LITA
(2026) results, can trilinear aggregation (TA) be improved further? This note
records the investigation and its conclusion: within the known TA toolkit, **no**
— and we now know precisely why, lever by lever. All work is over `Q`,
non-commutative, cubic `⟨N,N,N⟩`.*

## 0. What moved, and what didn't

The three modern TA closed forms (even `N`), lined up by coefficient:

| construction | `K⟨N,N,N⟩:m` | `N³` | `N²` | `N¹` | `N⁰` |
| --- | --- | --- | --- | --- | --- |
| Hadas–Schwartz 1982 | `(4N³+45N²+128N+108)/12` | **1/3** | **15/4** | 32/3 = 64/6 | 9 |
| Schwartz–Zwecher 2025 | `(4N³+45N²+122N+96)/12` | **1/3** | **15/4** | 61/6 | 8 |
| LITA 2026 (even) | `(4N³+45N²+116N+84)/12` | **1/3** | **15/4** | 29/3 = 58/6 | 7 |

- `N³/3` is **axis-count-bound** (fold 3 index-rotations into 1 product; only 3
  tensor factors) — no TA variant has ever moved it, and none can.
- `15/4·N²` has been **frozen since Hadas–Schwartz 1982**. SZ and LITA both leave
  it *exactly* untouched.
- All recent progress is in the linear + constant terms only, and it is startlingly
  regular: linear `64/6 → 61/6 → 58/6` (−`N/2` per paper), constant `9 → 8 → 7`.

So the whole question reduces to: **can the frozen `15/4·N²` move?**

## 1. Where the `2·N²` half lives (diagnostic)

Split `15/4 = 1.75 (correction) + 2.0 (aggregation tail)`. The aggregation families
(a) Ŝ and (b) Ṡ\Ṡ1 emit one product per index-triple `(i,j,k)∈[d]³`, `d=N/2+1`.
Bucketing those products by how many indices are **distinct** and fitting each
count exactly (a cubic in `d`) attributes the `N²` coefficient:

| family | index class | count in `d` | → `N²` coeff |
| --- | --- | --- | --- |
| (a) Ŝ | 3-distinct (generic) | `⅔d³ − 2d² + 4/3 d` | **0** |
| (a) Ŝ | 2-distinct (boundary) | `2d² − 2d` | **1/2** |
| (b) Ṡ | 3-distinct (generic) | `2d³ − 6d² + 4d` | **0** |
| (b) Ṡ | 2-distinct (boundary) | `6d² − 6d` | **3/2** |
| (b) Ṡ | 1-distinct (diagonal) | `d` | 0 (`O(N)`) |

**The generic triples carry the `N³/3` leading term with *zero* `N²` overhead. The
entire `2·N²` tail is the two-equal-index "boundary" triples** — the 2-D boundary
of the aggregation, one dimension up from the 1-D diagonal SZ optimized for the
*linear* term. `≈ 8` boundary products per ordered pair `(i,k)`, for only `≈ 3`
intrinsic trace terms each.

*Repro:* `docs.explore.ProbeAggTail`.

## 2. The four levers, and why each is closed

### Lever 1 — a better correction base (e.g. `⟨3,3,6⟩`). Provably WORSE.

The correction (family (c)) is a **2-D grid of blocks** (`Σ_{i≠j}`), not a 3-D
recursion. Tiling it with `⟨g,g,g⟩` blocks costs `(N/g)²·R_g = (R_g/g²)·N²`, so the
relevant efficiency is **rank-per-`g²`**, not the ω-style rank-per-`g³`:

| base | `R_g` (NC/Q) | `R_g/g²` |
| --- | --- | --- |
| `⟨2,2,2⟩` | 7 | **1.75** |
| `⟨3,3,3⟩` | 23 | 2.56 |
| `⟨4,4,4⟩` | 48 | 3.00 |
| `⟨6,6,6⟩` | 153 | 4.25 |

This is **provable**: the classical `rank⟨g,g,g⟩ ≥ 2g²−1` gives `R_g/g² ≥ 2 − 1/g²`,
strictly increasing in `g`, and `⟨2,2,2⟩ = 7` meets it exactly. **`g=2` (Strassen)
is optimal for the correction tiling** — bigger bases (`⟨3,3,6⟩` included) are
strictly worse. The `1.75·N²` correction half is essentially optimal.

*(Honesty tier: proven, modulo the block-tiling structure all TA uses.)*

### Lever 2 — harvest leftover unifications. Provably EMPTY.

A product unites with another (a −1) only if the two are **proportional on two of
the three axes** (`u⊗v⊗w + u⊗v⊗w' = u⊗v⊗(w+w')`). Measuring every product's
factor-direction in BOTH the delivered (pulled-back) basis and the transformed
(A*/B*/C*) basis, for SZ *and* LITA at `⟨20³⟩`/`⟨28³⟩`:

- **0** shared factor-directions on any axis — every product is distinct.
- **0** two-axis kin. Family (c)'s off-diagonal products are mutually independent
  and share nothing with the other families.

Both constructions are fully **unification-complete**: there is nothing to harvest,
in either space. (SZ's diagonal kin is already consumed into `R(i)`.)

*Repro:* `docs.explore.ProbeTaKinGraph` (modes `sz N`, `lita N`, `szsym N`).

### Lever 3 — re-decompose the boundary block. TIGHT for non-degenerate `N`.

Let `R_bd = Σ(boundary products)` — exactly what the boundary must produce with the
rest held fixed (since `Σ(all products) = T_matmul`). If `rank(R_bd) < #boundary`,
swapping in a minimal decomposition gives a strictly better *valid* scheme. Testing
via rank-continuation CP-ALS warm-started from the exact decomposition:

| `n₀` | γ | `B = #boundary` | `rank(R_bd)` | verdict |
| --- | --- | --- | --- | --- |
| 4 | −2 (dyadic) | 48 | **36** (bounded coefs → true rank) | reducible — but a **γ-degeneracy** |
| 6 | −5/4 | 96 | **96** (can't drop 1) | **tight** |
| 8 | −4/5 | 160 | **160** (can't drop 1) | **tight** |

The `n₀=4` reduction is real but degenerate (dyadic `γ=−2`, `d=3` creates
coincidental low rank); it vanishes at the first non-degenerate case. For `n₀ ≥ 6`
the boundary products are an essentially minimal decomposition of their own sum —
**the boundary is locally tight**; you cannot lower `N²` by re-decomposing it while
holding the rest fixed.

*(Honesty tier: strong empirical. CP-ALS is incomplete, but the method demonstrably
finds reductions when they exist — it found n₀=4's −12 — and finds none at n₀=6,8.)*

*Repro:* `docs.explore.ProbeBoundaryReducible n0`.

### Lever 4 — engineer NEW kin (the LITA-inspired lever). Structurally CLOSED.

SZ's diagonal trick (Cor 3.2) worked because a **tunable** sub-algorithm (the
diagonal `⟨2,2,2;7⟩`) was **co-located** with the target locus (same `{i,ī}` block),
so its coefficients could be tuned to make it *kin* with a diagonal aggregation
product → unite → −`O(N)` → the linear term. LITA systematized the search for such
re-pairings. The natural next step: do the same for the 2-D boundary → move `N²`.

It cannot be done, for two exactly-measured reasons (`n₀ = 6, 8, 10`):

- Every boundary product has **all `±1` coefficients** — rigid, nothing to tune on
  the boundary side (`96/96`, `160/160`, `240/240`).
- Every boundary product is **support-isolated**: it shares its 2-axis support with
  **zero** other products (`0/96`, `0/160`, `0/240`). Unification needs identical
  support on two axes; no coefficient search over φ or the correction `⟨2,2,2;7⟩`
  orbit can create it, because the variable **entries** are disjoint, not merely
  their coefficients.

**Structural reason.** SZ's mechanism needs a *tunable* block *co-located* with the
target. The 1-D diagonal has one (the diagonal correction); the 2-D boundary is
pure rigid aggregation with its own variables — no co-located tunable partner
exists. The trick that moved the linear term has **no analog** for the quadratic
term.

*(Honesty tier: proven per tested `n₀` — support-isolation is an exact computation;
a failed necessary condition, not a failed search.)*

*Repro:* `docs.explore.ProbeBoundaryKinPotential n0`.

## 3. Conclusion

Every lever in the SZ/LITA/flip-graph toolkit is closed:

1. better base → provably worse (2-D tiling favors `g=2`);
2. harvest leftovers → provably empty (0 shared directions);
3. compress the boundary → tight for `N ≥ 6`;
4. engineer kin → structurally closed (boundary rigid + support-isolated).

Combined with the empirical fact that **three independent constructions across 44
years (HS 1982, SZ 2025, LITA 2026) all land on exactly `15/4`**, this is strong
evidence — short of a formal lower-bound proof — that `15/4·N²` is a genuine **floor
for the trilinear-aggregation family**, not an oversight.

## 4. The one open path

The only thing left is a fundamentally different **aggregation skeleton** — one where
the 2-D boundary is handled by a *co-located tunable block* (giving it the structure
the 1-D diagonal has) instead of rigid aggregation, so the SZ kin-mechanism gains an
`O(N²)` analog. That is a from-scratch construction, not a search over any existing
parameters, and it is exactly what nobody — including LITA — has found. This is the
real, hard frontier for TA.

## 5. Instrumentation (all in `papers/schwartzzwecher2025/TaNew25Construction`)

Added for this analysis; `build(n0)` remains exact-verified throughout:

- `buildSymbolicForms(n0)` — transformed-space A*/B*/C* forms per product.
- `buildTagged(n0)` — scheme + per-product class (0 generic / 1 boundary /
  2 correction / 3 diagonal).
- `buildSymbolicTagged(n0)` — both of the above together.

Probes (in `docs.explore`): `ProbeAggTail`, `ProbeTaKinGraph`,
`ProbeBoundaryReducible`, `ProbeBoundaryKinPotential`.
