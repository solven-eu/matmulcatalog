# Trilinear aggregation / padded Kronecker

How fast matmul algorithms for one format are derived from another.
Companion to [SMALL_MATMUL_CATALOG.md](SMALL_MATMUL_CATALOG.md) §6
(recursive composition framework) and
[REFERENCES.md](../../REFERENCES.md) §[21] (Sedoglavic 2017, with local PDF).

---

## 1. Three composition strategies

For target `⟨n, m, p⟩` there are three principled ways to build an
algorithm from smaller pieces:

### A. Pure Kronecker product

`R(⟨n₁·n₂, m₁·m₂, p₁·p₂⟩) ≤ R(⟨n₁, m₁, p₁⟩) · R(⟨n₂, m₂, p₂⟩)`.

Multiplicatively combines two algorithms. Requires the target axes to
**factor exactly** — e.g. `n = n₁·n₂` for some integers. For prime `n`
(7, 11, 13, …) only `n = 1·n` works, which is trivial. So pure
Kronecker has nothing to say about prime-dimension cubes like
`⟨7,7,7⟩`.

### B. Naïve block decomposition

Split each axis into a sum, e.g. `7 = 3 + 4`. The matmul becomes an
outer 2×2 block tensor of mixed-shape sub-products. No sharing across
blocks. Always works but loses cross-block addition opportunities — see
§3 below for the `⟨7,7,7⟩` worked example.

### C. Sedoglavic's structured padding (this section)

Take `⟨n,n,n⟩` with `n` "almost" factorable. Split each axis into a
sum where one part has a great existing algorithm, and **pad the other
part with a zero row+col** so it matches a size with a great algorithm
too. Then apply Strassen-style composition over the resulting block
structure, and prove that many resulting multiplications hit purely
zero inputs and can be dropped.

This is **not** "enumerate all padding bases and search." It's a
structured construction tailored per target.

---

## 2. The `⟨7,7,7⟩ = 250` recipe (Sedoglavic 2017)

The construction documented in
[`references/papers/sedoglavic_2017_7x7x7_r250.pdf`](../../references/papers/sedoglavic_2017_7x7x7_r250.pdf):

1. **Split each `7` as `3 + 4`** along all three axes. The 7×7 matrices
   `A`, `B`, `C` become 2×2 block matrices with mixed-shape blocks:

   ```
       ⎡ 3×3  3×4 ⎤
   A = ⎢          ⎥
       ⎣ 4×3  4×4 ⎦
   ```

2. **Pad the `3` side to `4`** by inserting a zero row+col **in the
   middle** (at position 4 of the 7-axis). Now every "block" is the
   same shape 4×4, and the matrix lives in a `(4+4) × (4+4) = 8×8`
   structure with a known zero row+col at the inserted middle.

3. **Apply `Strassen × Strassen²`** over this 2×2-of-4×4 structure —
   that's Strassen³ structurally, with rank `7 · 49 = 343`
   multiplications **before** clipping.

4. **Drop the multiplications whose inputs are entirely supported on
   the inserted zero row+col.** Sedoglavic's analysis proves exactly
   `250` multiplications remain active.

The 21-multiplication saving over the naïve block decomposition (271,
see §3) comes from the same source as Strassen's saving over naïve
`⟨2,2,2⟩=8`: cross-block A-entry sums and C-entry sharing inside the
Strassen³ structure that the per-block decomposition cannot express.

---

## 2bis. Sedoglavic's closed-form identity

The whole calculation collapses to a one-line algebraic identity for
any cubic `⟨n,n,n⟩` with a split `n = u + v`, `u > v ≥ 1`:

> **`⟨u+v, u+v, u+v⟩ ≤ ⟨u,u,u⟩ + 3·⟨u,u,v⟩ + 3·⟨u,v,v⟩`**

(rank inequality; ranks are catalog values).

The naïve block decomposition of `⟨u+v,u+v,u+v⟩` would give 8
sub-products grouped by shape:

| shape | count in naïve block decomp | in Sedoglavic identity |
|---|---|---|
| `⟨u,u,u⟩` | 1 | 1 |
| `⟨u,u,v⟩` (× 3 axis-perms) | 3 | 3 |
| `⟨u,v,v⟩` (× 3 axis-perms) | 3 | 3 |
| `⟨v,v,v⟩` | 1 | **0** ← the saving |

The structured padding consolidates the `⟨v,v,v⟩` corner block into
sharing with the other terms — the missing 8th sub-product is the
exact savings.

**Equivalent constructive view (implemented in this repo):** the
saving comes from using **Strassen's 7-mult `⟨2,2,2⟩`** as the outer
2×2-block algorithm (instead of the naïve 8-mult block decomposition)
combined with **min-shape reduction per Strassen mult**. Strassen's 7
mults each see a different sub-shape:

| Strassen mult | sub-shape (with `u=4, v=3`) | sub-rank |
|---|---|---|
| `M0 = (A₁₁+A₂₂)(B₁₁+B₂₂)` | `⟨4,4,4⟩` | 49 |
| `M1 = (A₂₁+A₂₂)·B₁₁` | `⟨3,4,4⟩` | 38 |
| `M2 = A₁₁·(B₁₂−B₂₂)` | `⟨4,4,3⟩` | 38 |
| `M3 = A₂₂·(B₂₁−B₁₁)` | `⟨3,3,4⟩` | 29 |
| `M4 = (A₁₁+A₁₂)·B₂₂` | `⟨4,3,3⟩` | 29 |
| `M5 = (A₂₁−A₁₁)·(B₁₁+B₁₂)` | `⟨3,4,3⟩` | 29 |
| `M6 = (A₁₂−A₂₂)·(B₂₁+B₂₂)` | `⟨4,3,4⟩` | 38 |

Total: **49 + 3·38 + 3·29 = 250** — Sedoglavic's bound, reached by a
constructive `Recombination.constructWithAllocation(Strassen, lookup,
[4,3], [4,3], [4,3])` call. See
[`TestStrassenRecombination777`](../../src/test/java/io/cormoran/strassen/v3/catalog/TestStrassenRecombination777.java);
the resulting algorithm verifies under
{@code Verifier.isExactNonCubic}.

The trick: each Strassen mult's input (e.g. `A₁₁ + A₂₂`) involves
*blocks of mixed shape* (`4×4 + 3×3`); the smaller block is
implicitly zero-padded to `4×4` so the sum is well-defined. The
zero-padded positions contribute nothing to the output, so the
min-reduction in
[`processAdditions`](../../src/main/java/io/cormoran/strassen/v3/catalog/Recombination.java)
collapses each sub-mult to its truly-active shape.

Some numeric checks:

- `⟨7,7,7⟩` with `u=4, v=3` → `49 + 3·29 + 3·38 = 250` (naïve was 273;
  saving is `⟨3,3,3⟩=23`).
- `⟨5,5,5⟩` with `u=3, v=2` → `23 + 3·11 + 3·14 = 98` (naïve would be
  `23 + 33 + 42 + 8 = 106`; saving is `⟨2,2,2⟩=8`). Still worse than
  the direct AlphaEvolve `93`.
- `⟨7,7,7⟩` with `u=5, v=2` → `93 + 3·22 + 3·18 = 213` IF we had those
  ranks; in practice `⟨5,5,2⟩` best ≈ `25`, `⟨5,2,2⟩` ≈ `18`, giving
  `93 + 75 + 54 = 222`. Worse than the `u=4, v=3` split.

So **for each target `⟨n,n,n⟩`, try every `n = u + v` split with
`u > v`, pick the min RHS.** That's `⌊n/2⌋` evaluations, each a
catalog lookup. Cheap.

The Java class
[`BlockSplitSearch`](../../src/main/java/io/cormoran/strassen/v3/catalog/BlockSplitSearch.java)
automates this; `main` prints the per-target table. Once a split is
chosen, hand `(n, u, v)` to
[`Compose.blockSplitCubic`](../../src/main/java/io/cormoran/strassen/v3/catalog/Compose.java)
to construct the actual algorithm (Phase 2).

**Cross-field caveat.** The formula reads ranks from the catalog
without checking which field each comes from — so a Phase-1
prediction might mix `R_{F₂}` (AlphaTensor `⟨4,4,4⟩=47`) with `R_R`
(Laderman `⟨3,3,3⟩=23`). The result is a valid upper bound only over
whichever field is most restrictive. For a field-pure result, restrict
the lookup to schemes of the target field — see
[`feedback-state-field-explicitly`](../../../memory/feedback_state_field_explicitly.md).

## 2ter. Why zero placement doesn't matter

A natural question: if we pad `⟨7,7,7⟩` into `⟨8,8,8⟩`, does it matter
WHERE the zero row+col goes — at position 4 (Sedoglavic's middle
placement) vs position 8 (end) vs position 1?

**No** — the active-multiplication count is the same.

The matmul tensor and its decompositions are invariant under axis
permutations of A, B, C (independently). Permuting A's rows + C's rows
in the same way leaves both the algorithm and the count unchanged.
Different "zero placements" are just different choices of which
permutation to apply before clipping; they produce equivalent algorithms
with the same rank.

What placement DOES affect:
- The specific surviving multiplications' indices (a relabeling).
- Implementation details (memory layout, vectorization, cache).
- The visual structure of the resulting factor matrices.

Sedoglavic's middle placement is about getting a clean block-aligned
8×8 structure for the analysis, not about a better count.

## 3. Why `250 < 271`: cross-block sharing

Compare to the naïve block decomposition (Strategy B) of `⟨7,7,7⟩` with
the same 7 = 3+4 split:

| sub-product | format | best mults |
|---|---|---|
| `A₃₃·B₃₃` | ⟨3,3,3⟩ | 23 (Laderman) |
| `A₃₃·B₃₄`, `A₃₄·B₄₃` | ⟨3,3,4⟩, ⟨3,4,3⟩ | 29 each |
| `A₃₄·B₄₄`, `A₄₃·B₃₃` | ⟨3,4,4⟩, ⟨4,3,3⟩ | 38, 29 |
| `A₄₄·B₄₃`, `A₄₃·B₃₄`, `A₄₄·B₄₄` | ⟨4,4,3⟩, ⟨4,3,4⟩, ⟨4,4,4⟩ | 38, 38, 47 |

**Total: 271** — every sub-product is computed independently, no cross-
block sharing.

Sedoglavic's structured padding embeds the 7=3+4 split into an
Strassen³-on-8×8 framework where multiplications like
`(a₁₁ + a₂₂ + … + a₈₈)·(…)` combine A entries from BOTH blocks before
the multiplication. That single multiplication consolidates work that
the naïve decomposition would have to recompute per block.

---

## 3bis. Field discipline: once per field-class

The block-split / Sedoglavic identity holds **as a rank inequality
over any field**, but a *constructive* algorithm inherits the field
of whichever sub-algorithms it stitches together. You cannot mix
sub-algorithms from different fields in one construction — their
addition operations disagree (mod-2 XOR is not real-arithmetic
addition).

So the procedure is run **once per field-class**. For matmul, the
field-classes empirically collapse to three:

| field-class | typical UB at `⟨4,4,4⟩` | typical UB at `⟨5,5,5⟩` |
|---|---|---|
| `{Q, Z, R}` (single cluster) | 49 (Strassen²) | 93 (AlphaEvolve) |
| `C` | 48 (AlphaEvolve 2025) | 93 |
| `F₂` (and other `F_p`) | 47 (AlphaTensor 2022) | open |

### Why Q/Z/R cluster as one

All three are *characteristic-0 subfields of R*, and matmul-rank is the
same across them in practice for the following reasons:

1. **Q ⊆ R ⊆ C, so `R_R(T) ≥ R_Q(T)`** (a Q-algorithm is automatically
   an R-algorithm; a smaller-rank Q-algorithm would beat any R one).
   Conversely an R-algorithm uses real coefficients, which are at
   least as general as Q.
2. **Z ↔ Q by scaling.** Any rational-coefficient algorithm can be
   converted to integer coefficients by multiplying through by a
   common denominator — preserves rank (which counts multiplications,
   not coefficient magnitude). Conversely any integer algorithm
   trivially has rational coefficients.
3. **No known matmul algorithm beats Q with irrational coefficients.**
   All published-best algorithms over R (Strassen, Laderman, Smirnov,
   AlphaTensor's standard-arithmetic ones, AlphaEvolve's R-mode)
   use coefficients in `{-1, 0, 1, ±½}` — which are rational. There
   is no theorem saying R can't do strictly better than Q for matmul
   (such a strict gap is possible in principle), but for matmul
   tensors specifically no example is known.

So the rule of thumb: run the trilinear aggregation construction once
with `{Q, Z, R}` sub-algorithms (treat them as one pool), once with
`C` sub-algorithms (only useful when complex coefficients buy
something — known for `⟨4,4,4⟩` and `⟨5,5,5⟩` via AlphaEvolve), and
once with `F₂` sub-algorithms (mod-2 arithmetic; AlphaTensor's wins).

This matters for `BlockSplitSearch` Phase 1: if the catalog lookup
silently mixes fields (returns `R_{F₂}(⟨4,4,4⟩)=47` alongside
`R_Q(⟨3,3,3⟩)=23`), the formula sum is a *cross-field upper bound* —
mathematically valid but not realisable as a single-field algorithm.
See memory entry [`feedback-state-field-explicitly`](../../../memory/feedback_state_field_explicitly.md).

## 3ter. Cross-base comparison for `⟨7,7,7⟩`

The block-split construction works with **any** outer base, not just
Strassen `⟨2,2,2⟩`. For each base `⟨b,b,b⟩` and each non-degenerate
allocation summing to 7 over `b` blocks, compute the formula RHS and
pick the min. Result (from
[`Test777AcrossBases`](../../src/test/java/io/cormoran/strassen/v3/catalog/Test777AcrossBases.java)):

| outer base | base rank | best allocation | `⟨7,7,7⟩` rank |
|---|---|---|---|
| `⟨2,2,2⟩` Strassen | 7 | `[3,4]³` (or `[4,3]³`) | **250** ← Sedoglavic |
| `⟨3,3,3⟩` Laderman | 23 | `[2,2,3]³` | 281 |
| `⟨4,4,4⟩` (Strassen²) | 49 | `[1,2,2,2]³` | 280 |
| `⟨5,5,5⟩` AlphaEvolve | 93 | `[2,1,2,1,1]³` | 319 |
| `⟨6,6,6⟩` | 153 | `[1,1,1,1,2,1]³` | 279 |

So Strassen wins. Two non-obvious side observations:

- **The specific algorithm matters at the same rank.** Two equivalent
  rank-7 `⟨2,2,2⟩` algorithms can give different `⟨7,7,7⟩` totals
  under the same `[4,3]³` allocation, because `processAdditions`
  reduces each sub-mult's shape based on the **support pattern** of
  its `U/V/W` columns — which (block-row, block-col) positions are
  nonzero. Canonical Strassen has carefully balanced support
  (e.g. `M_2 = (A_{21}+A_{22})·B_{11}` touches row-block 1 only,
  giving small sub-row-dim under `[4,3]`); other rank-7 algorithms
  (e.g. `alphatensor-Z_2x2x2_r7_a22`) have different support and
  lose 5 mults to coarser sub-shapes:

  | algorithm | sub-shape sum | total |
  |---|---|---|
  | `strassen_2x2x2_r7_a18` | `49 + 3·38 + 3·29` | **250** |
  | `alphatensor-Z_2x2x2_r7_a22` | `29 + 49 + 29 + 38 + 38 + 49 + 23` | 255 |

  Both are valid `⟨7,7,7⟩` constructions; the second is just less
  efficient under min-reduction. **A complete search would enumerate
  not just `(base_format, allocation)` but `(specific_algorithm,
  allocation)`.** For now we pin the canonical file per format
  (Strassen for `⟨2,2,2⟩`, Laderman for `⟨3,3,3⟩`, etc.) when
  measuring.
- **Bigger outer bases don't necessarily win.** `⟨6,6,6⟩` (rank 153)
  gives 279 — better than `⟨3,3,3⟩` (281) — because the
  fine-grained 6-block allocation captures more block-level
  sharing. But Strassen `⟨2,2,2⟩` still wins overall because its
  7-mult outer structure is exceptionally well-aligned with the
  8-block decomposition tensor's symmetries.

## 3quart. Improvement propagation (systemic effect)

**Sedoglavic's formula is monotone in its sub-format inputs.** Concretely:

> If `R(⟨a,a,a⟩)`, `R(⟨a,a,b⟩)`, or `R(⟨a,b,b⟩)` improves for any
> small `a, b`, then **every** larger cubic target `⟨a+b, a+b, a+b⟩`
> immediately gets a (possibly) tighter upper bound — and
> recursively, every cubic target above it that uses these sizes in
> its split.

This is a fundamental reason to chase small-format improvements: a
better `⟨4,4,4⟩` algorithm doesn't just improve `⟨4,4,4⟩` — it
propagates to `⟨7,7,7⟩`, `⟨8,8,8⟩` (via `4+4`), then to
`⟨14,14,14⟩` (via `7+7` using improved `⟨7,7,7⟩`), then up to
`⟨16,16,16⟩` via `⟨4,4,4⟩²`, etc.

The
[`TestSedoglavicCoverageGap`](../../src/test/java/io/cormoran/strassen/v3/catalog/TestSedoglavicCoverageGap.java)
test makes this systemic. It scans every cubic target `⟨n,n,n⟩` for
`n ∈ [4, 32]` and reports:

- **Gaps**: targets where the formula bound is below the direct
  catalog rank — i.e. we can already do better by constructing.
- **Missing directs**: targets with no direct catalog scheme — the
  formula bound is then the SOTA we know.

The test runs in **~150 ms** (just catalog map lookups + `⌊n/2⌋`
formula evaluations per target). It serves as a **regression guard**:
introducing a new improved small-format scheme automatically surfaces
which big targets newly benefit. Wire it into the CI loop so any
catalog edit produces an updated coverage delta.

Current state (2026-05-27), **per field-class** — the test now scans
each of `{R, C, F₂}` independently via
`BlockSplitSearch.loadCatalogBestRanksForField(...)`:

| field | gaps | missing directs | notable |
|---|---|---|---|
| R (=Q=Z) | 0 | `⟨17..31⟩³` (15) | catalog covers `n ≤ 16` with direct schemes ≥ formula |
| C | 1 | `⟨17..31⟩³` (15) | **`⟨7,7,7⟩` direct=250, formula=249** via `48 + 3·38 + 3·29` using AlphaEvolve's `⟨4,4,4⟩=48` over C — a legitimate constructive improvement to materialise |
| F₂ | 0 | most `n ≥ 6` (catalog only has AlphaTensor `⟨2..5⟩³` directs) | formula gives the only F₂ bound we know for medium-sized targets |

The C-field gap is the most actionable: running
`Recombination.constructWithAllocation(Strassen, c_lookup, [4,3]³)`
would yield a verified `⟨7,7,7⟩` with **249 multiplications over C** —
beating both our direct entry (250) and the upper bound any R-field
scheme can offer (250 from Sedoglavic 2017).

Removing the field discipline (mixing across fields) gives an
unrealisable cross-field upper bound of 248 for `⟨7,7,7⟩` — but no
single-field algorithm can achieve this; see
[`docs/derived-from-cited-bounds.json`](../../docs/derived-from-cited-bounds.json) where the
formula uses `⟨4,4,4⟩=47` from F₂ + `⟨3,3,3⟩=23` from R, which
cannot be combined since their addition operations disagree.

## 3quint. Non-cubic targets — the real testing ground

The block-split / Strassen-recombine construction generalises to
**any** `⟨n, m, p⟩`. The Phase-1 search is
`BlockSplitSearch.findBestSplitNonCubic(n, m, p, strassen, sota)`:
for each axis independently, enumerate splits
`(u_n + v_n, u_m + v_m, u_p + v_p)` with all parts ≥ 1, evaluate via
`Recombination.recombineWithAllocation(Strassen, ...)`, pick the min.

This is where the procedure earns its keep: **most catalog gaps are
non-cubic**, because non-cubic targets are simultaneously
(a) the natural sub-products that bigger cubic constructions need and
(b) where AlphaTensor / AlphaEvolve direct schemes tend to be locally-optimal
without exploiting cross-block sharing.

The same `TestSedoglavicCoverageGap` test scans canonical non-cubic
targets `n ≤ m ≤ p` (excluding cubic). Current state, **R-field
only, dim ≤ 12** (runs in ~0.3 sec):

| signal | count | example |
|---|---|---|
| at-or-better | 253 | catalog matches/exceeds formula |
| **gaps** | **8** | `⟨8,8,10⟩` direct=441, formula=427 (Δ=14) via `[4,4]/[4,4]/[5,5]` |
| missing direct | 14 | `⟨2,9,9⟩`, `⟨4,9,9⟩`, … — formula gives the only SOTA |

Top actionable gaps:

| target | direct R | formula | Δ | recipe (Strassen allocs) |
|---|---|---|---|---|
| `⟨8,8,10⟩` | 441 | **427** | 14 | `[4,4] / [4,4] / [5,5]` |
| `⟨4,8,10⟩` | 230 | **224** | 6 | `[2,2] / [4,4] / [5,5]` |
| `⟨4,8,9⟩` | 209 | **206** | 3 | `[2,2] / [4,4] / [5,4]` |
| `⟨4,7,10⟩` | 206 | **203** | 3 | `[2,2] / [4,3] / [5,5]` |
| `⟨5,8,10⟩` | 286 | **284** | 2 | `[3,2] / [4,4] / [5,5]` |
| `⟨4,7,9⟩` | 187 | **186** | 1 | `[2,2] / [4,3] / [5,4]` |
| `⟨4,8,11⟩` | 253 | **252** | 1 | `[2,2] / [4,4] / [6,5]` |
| `⟨5,8,11⟩` | 313 | **312** | 1 | `[3,2] / [4,4] / [6,5]` |

Each row is an actionable
`Recombination.constructWithAllocation(canonicalStrassen,
RLookup, allocA, allocB, allocC)` call away from a verified
better algorithm to save to `section{N}/composed-strassen-mixed_NxMxP_rR.json`.

## 3sext. Catalog audit: AlphaTensor / AlphaEvolve non-cubic

We integrate both AT and AE non-cubic results — `tools/import_alphatensor.py`
covers **93 AT-Z** non-cubic formats up to `⟨12,12,12⟩` territory plus
**20 AT-F2** small non-cubic (`⟨2..5⟩³` + small mixed), and
`tools/import_alphaevolve.py` covers **15 AE** schemes (small non-cubic
+ `⟨3,3,3⟩=23`, `⟨4,4,4⟩=48`, `⟨5,5,5⟩=93`).

## 3sept. Extended trilinear scan (R / C / F₂, dim ≤ 16)

After integrating AT/AE, a per-field non-cubic block-split scan
(`TestSedoglavicCoverageGap`, ~1 sec) reveals where Strassen-recombine
beats or fills in the catalog. Summary:

| field | catalog ≥ formula | gaps (formula beats catalog) | missing direct |
|---|---|---|---|
| R | 625 | **14** | 26 |
| C | 602 | **37** | 26 |
| F₂ | 16 | 0 | **215** |

**Top R-field gaps** (actionable: `Recombination.constructWithAllocation` materialises them):

| target | direct R | formula | Δ | recipe (allocs) |
|---|---|---|---|---|
| `⟨13,16,16⟩` | 2038 | 2006 | **32** | `[8,5] / [8,8] / [8,8]` |
| `⟨8,8,10⟩` | 441 | 427 | 14 | `[4,4] / [4,4] / [5,5]` |
| `⟨7,16,16⟩` | 1164 | 1158 | 6 | `[4,3] / [8,8] / [8,8]` |
| `⟨4,8,10⟩` | 230 | 224 | 6 | `[2,2] / [4,4] / [5,5]` |

**C-field is even richer**: 37 gaps, all enabled by AlphaEvolve's
`⟨4,4,4⟩=48` propagating through composition. E.g. `⟨7,7,7⟩=249`
(vs direct R=250) — saving 1 mult by using AE's C-specific sub-block.

**F₂-field has 215 missing directs** for `n,m,p ≤ 16`. AT-F2 only
covers small cubes/non-cubics; the formula gives the first known F₂
bound for hundreds of medium targets, e.g. `⟨6,6,6⟩=167` via
`Strassen × AT-F2_4x4x4=47 + 3·AT-F2_4x4x2 + 3·AT-F2_4x2x2`.

All 169 above (cubic mixed-field + non-cubic per-field gaps + missing
directs) are now in [`docs/derived-from-cited-bounds.json`](../../docs/derived-from-cited-bounds.json)
and appear as italicised purple rows in the Pages browser.

See [`COVERAGE.md`](../../generated/COVERAGE.md) for the full per-format matrix
across fields, and `references/fmm-lille-discrepancies.md` for our
catalog vs fmm-lille comparison.

## 4. When padding can be WORSE

Padding is not magic. Direct algorithms can outperform any padded
construction.

### `⟨5,5,5⟩` — direct wins

- **Direct AlphaEvolve scheme** (`alphaevolve_5x5x5_r93`): **93 mults**.
- Padding into `⟨6,6,6⟩` (best base = 153): even after clipping, almost
  certainly > 93.
- Padding into `⟨8,8,8⟩` Strassen³ (343 unpadded): clipping helps but
  the base is far too oversized.

For `⟨5,5,5⟩` the direct scheme wins decisively. There's no Sedoglavic-
style structured padding that beats it.

### `⟨3,3,3⟩` — direct Laderman wins

`23` mults from Laderman 1976. No padded approach (e.g. from
`⟨4,4,4⟩=47` over F₂) beats it.

### General rule of thumb

| target | best strategy |
|---|---|
| factorable on every axis | pure Kronecker if both factor schemes are good |
| **prime / awkward, no direct scheme** | **Sedoglavic-style structured pad** (split + pad the small side, embed in next factorable size) |
| small + well-studied (`⟨2,2,2⟩`, `⟨3,3,3⟩`, `⟨4,4,4⟩`, `⟨5,5,5⟩`) | direct hand-crafted scheme |
| `⟨1, m, p⟩` or `⟨n, 1, p⟩` | trivial: `n·m·p / min(n,m,p)` mults |

---

## 5. Where the structured recipe scales

Sedoglavic's recipe generalises to any `⟨n,n,n⟩` where `n + small` is
factorable:

| target | split (per axis) | pad to | outer × inner |
|---|---|---|---|
| `⟨7,7,7⟩` | `7 = 3+4` | each block 4×4 | `Strassen × Strassen² = Strassen³` |
| `⟨5,5,5⟩` | `5 = 3+2` (pad 2→3) | each block 3×3 | `Strassen × Laderman` |
| `⟨5,5,5⟩` | `5 = 4+1` (pad 1→4) | each block 4×4 | `Strassen × Strassen²` |
| `⟨11,11,11⟩` | `11 = 6+5` (pad 5→6) | each block 6×6 | `Strassen × ⟨6,6,6⟩` |
| `⟨11,11,11⟩` | `11 = 8+3` (pad 3→8) | each block 8×8 | `Strassen × Strassen³` |
| `⟨13,13,13⟩` | `13 = 8+5` (pad 5→8) | each block 8×8 | `Strassen × Strassen³` |
| `⟨13,13,13⟩` | `13 = 6+7` (pad 7→8 too) | each block 8×8 | `Strassen × Strassen³` |

Each recipe gives **one** padded rank to compute (run
`Recombination.recombine` with the specific base + allocation). Then
compare across recipes, pick the min.

**This is targeted enumeration — not brute force over the whole
catalog.** A correct implementation should:
1. For each target `⟨n,n,n⟩`, propose a small set (≤ ~5) candidate
   splits like the table above.
2. For each candidate, construct the corresponding `Strassen × inner`
   base scheme by Kronecker product over our catalog.
3. Call `Recombination.recombine(target, target, target, base, sota)`
   to count active multiplications under the right allocation.
4. Pick the min.

(Non-cubic targets follow the same recipe per axis independently.)

---

## 6. Status & open work

- **`Recombination.recombine`** (Java) — does the structured rank count
  for any (target, base, allocation) tuple. The math is right.
- **Targeted recipe enumeration** for a single target — TODO. Should
  encode the table above as a small per-target candidate list, NOT
  brute-force over all catalog bases.
- **`FindBestPaddedCompositions.java`** — written but currently
  over-eagerly iterates every catalog scheme. Either rework into the
  targeted-recipe form, or restrict its scope to single-target probes.
- **Trivial `⟨1, m, p⟩` / `⟨n, 1, p⟩`** schemes — needed for the
  sub-problem lookup to be grounded (currently the SOTA resolver
  fallback covers `⟨1,1,1⟩=1` etc. by formula).

Tracked in [ROADMAP.md](../../ROADMAP.md).

---

## 7. References for further reading

- **Sedoglavic 2017** (this construction): [HAL](https://hal.science/hal-01572046v2/file/RFC1708.pdf)
  / [local PDF](../../references/papers/sedoglavic_2017_7x7x7_r250.pdf).
- **AlphaTensor 2022** generalised the recombination/aggregation idea
  with RL search over base allocations — see
  [REFERENCES.md §12](../../REFERENCES.md#12-alphatensor) and the Java port
  `io.cormoran.strassen.v3.catalog.Recombination` (which IS the
  `_process_additions` logic).
- **Schwartz–Zwecher 2025** ([§19](../../REFERENCES.md#19-schwartz-zwecher25))
  pushes feasible ω down to `O(n^{2.773203})` using a similar style of
  trilinear aggregation at large bases.
