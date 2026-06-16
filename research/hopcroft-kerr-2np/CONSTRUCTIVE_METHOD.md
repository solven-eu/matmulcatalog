# A constructive realisation of the Hopcroft–Kerr bound for ⟨2,p,n⟩

Status: **COMPLETE — integer schemes, exact-verified** (tasks #7/#8/#10/#11,
2026-06-11/12). 465 schemes over `3 ≤ p ≤ 32, p ≤ n ≤ 32` in
`schemes/constructed/`, 459 at the exact formula, all over ℤ, dual-certified.
This note records the complete constructive recipe and its proofs — the paper
section is `paper/sections/hk71.tex` (kept in sync). Code:
`src/main/java/eu/solven/matmul/papers/hopcroftkerr1971/` (`HopcroftKerr2bc`,
`HopcroftKerr2bcAsymmetric`, `LemmaOneAugmentation`); survey harness
`SurveyHk2npEmitterGaps`. Orientation/leftovers: [OVERVIEW.md](OVERVIEW.md);
the single open item is the g ≥ 6 family (task #9, see "Where we stopped"
at the bottom).

## Claim

For `⟨p,2,n⟩` (≅ `⟨2,p,n⟩` up to orientation), the construction below attains
the Hopcroft–Kerr 1971 bound

```
R(⟨2,p,n⟩) ≤ ⌈(3pn + max(p,n))/2⌉
```

**exactly**, with machine-verified INTEGER schemes (over ℤ, hence valid
over F₂/F₃/Q/R/C as well — all-unimodular Lemma-1, task #11), for **every**
shape in the swept range `3 ≤ p ≤ 32`, `p ≤ n ≤ 32` — band range
(`n ≤ 2p−1`: odd p via Case 1, even p via Case 2 + repaired Step 3) and
chained range (`n > 2p−1`, task #10) alike — **except** the six `g ≥ 6`
circulant shapes, which land at +1..+3 and are *provably*
formula-impossible within this framework (theorem below). Final tally:
**465 schemes, 459 at the exact formula**.

New-to-all-catalogs bounds produced (examples): `⟨2,10,15⟩=233` (every
catalog: 234), `⟨2,12,16⟩=296` (published best 298), `⟨2,11,12⟩=204`,
`⟨2,11,13⟩=221`, `⟨2,11,14⟩=238`, `⟨2,13,15⟩=300`, `⟨2,13,16⟩=320`,
`⟨2,10,16⟩=248`; beyond any published catalog range, e.g. `⟨2,13,25⟩=500`,
`⟨2,15,32⟩=736`. All claims are ℤ-field claims.

## The construction

Convention: `Y = Ā·X` with `Ā = M·A` (augmentation), `A ∈ K^{p×2}`,
`X ∈ K^{2×n}`, internal square-ish problem on `n` rows.

### Ingredients

1. **Diagonal methods (HK p9).** Each internal row `i ∈ [1,n]` gets a method
   `c(i) ∈ {1,2,3}`; each method computes the diagonal cell `y_ii` with 2
   products and leaves both products reusable:
   `1:{A=a₂(x₁+x₂), B=(a₁−a₂)x₁}`, `2:{C=(a₁−a₂)x₂, D=a₁(x₁+x₂)}`,
   `3:{E=a₂x₂, F=a₁x₁}`. Coloring: alternate `1,2,1,2,…` (even n) or
   `1,2,…,1,2,3` (odd n — the lone method-3 at position n).

2. **Different-method pairs (HK Lemma 2).** A pair of cells `(y_ij, y_ji)`
   with `c(i) ≠ c(j)` costs **3 new products** (reusing the four diagonal
   products of i and j). The `(1,2)` pair emits the cross products
   `E(a_i+a_j, x_i+x_j)`, `F(a_i+a_j, x_i+x_j)` and a mixed `G`.

3. **Same-method pairs with a bridge.** When `c(i) = c(j)`, the pair costs 3
   new products by routing through a third position `b` (the *bridge*),
   reusing (i) the diagonal products of `j`, (ii) the cross `E/F` products of
   the *pair* `(i,b)`, and (iii) subtracting the already-computed pair
   `(b,j)`. Sympy-verified identities exist for bridge methods opposite to the
   pair's: `(1,1,bridge-2)` and `(2,2,bridge-1)`.

4. **Bridge-selection lemma (ours — the key unblocking step).** In the cyclic
   distance-ordered processing (all pairs at distance d before distance d+1),
   a same-method pair `(i,j)` at distance d may use **any arc-interior
   position** `b = i+e (mod n)`, `1 ≤ e < d`, as its bridge: the pair `(i,b)`
   (distance e) and the pair `(b,j)` (distance d−e) are both already
   computed, and the `E/F` cross products of `(i,b)` exist for every
   different-method pair at any distance. With the alternating coloring, an
   interior position of the *opposite* method (1↔2) always exists:
   - `(1,1)` pairs: `i` odd ⇒ `i+1` even ⇒ method 2 (never the method-3
     position, which is odd-indexed). Take `e=1`.
   - `(2,2)` pairs: `i` even ⇒ `i+1` odd ⇒ method 1, except when `i+1 = n`
     (the method-3 position, odd n) — then `e=2` lands on position 1
     (method 1), and `d ≥ 3` holds for such pairs, so `e=2 < d` is legal.

   **Consequence: the `(2,2,bridge-3)` configuration never arises.** This
   sidesteps — rather than contradicts — the characteristic-0 impossibility
   theorem of `references/hopcroftkerr1971/README.md`: that theorem rules out
   3-product completions over a specific reusable set S; the schedule above
   simply never needs that case. (HK's "the other cases follow by symmetry"
   stays wrong as written; the *bound* is recovered by scheduling, not by the
   missing symmetry.)

5. **Lemma 1 augmentation, small-coefficient form (ours).** `M ∈ ℤ^{n×p}`,
   first p rows identity, every cyclic p-window nonsingular. Take augmented
   rows with entries in `{−2..2}` from a seeded deterministic stream, checked
   per window with **exact BigInteger Bareiss** determinants; retry on the
   rare singular draw. (The classical Vandermonde rows `(i+1)^k` satisfy the
   window property in theory but are computationally toxic: they overflow
   64-bit determinant checks and make the back-substitution numerically
   inexact from `n ≈ p+2`.)

6. **Exact back-substitution.** Per output column j, the internal band
   provides the cells `(band(j), j)`; with `M_j = M[band(j),:]` invertible,
   `A·x_j = M_j⁻¹ · (ĀX)_{band(j),j}`. Computed in **exact rational
   arithmetic** (BigFraction Gauss–Jordan), converted to double at the end —
   the output scheme is over ℚ — or over ℤ when the Lemma-1 matrix is
   **all-unimodular** (task #11, 2026-06-12): with every cyclic p-window
   determinant ±1, the window inverses are integer matrices and the
   back-substitution never leaves ℤ.

   **The Euclidean construction** (`LemmaOneAugmentation.buildUnimodular`;
   `sympy/unimodular/euclidean.py`). Window dets of `[I_p; B]` reduce to
   exactly three minor families of the augmented block B (m = n−p rows):
   leading k×k (cols 0..k−1), trailing w×w (rows m−w.., cols p−w..), and
   sliding m×m over contiguous column windows. Build B recursively:
   - **comb body**: `B[u][c] = 1 ⟺ c ≡ u (mod m)` on the first
     `p − (p mod m)` columns — sliding windows inside the body are
     permutation matrices and leading minors identity blocks (det ±1 by
     inspection);
   - **tail = transposed recursion**: for `r = p mod m > 0`, the last r
     columns are `B(m, r)ᵀ`. Eliminating the comb's permutation part shows
     the crossing-window and trailing minors of B equal (±) the
     leading/trailing/sliding minors of the tail block transposed — i.e.
     the SAME three-family problem at size `(m, r)`. The recursion descends
     like the gcd and terminates at `r = 0` (pure comb).

   Pure 0/1 entries, no search, O(p²); verified exhaustively over the whole
   band range `3 ≤ p ≤ 32` (`euclidean.py`: zero failures), with a
   belt-and-braces exact verification in the builder (fallback to the
   ternary draw — scheme stays ℚ — is loud, never silent). Discovered by
   examining seam solutions from an exact column-DFS
   (`comb_seam_coldfs.py`); the {−1,1}-dense rows of the previous draw can
   NEVER be unimodular for m ≥ 2 (any 2×2 ±1-minor is even) — which is
   exactly where the historical denominators 2/4/8/16 came from.

7. **Chained augmentation for n > 2p−1 (task #10).** One Lemma-1 band covers
   at most `2p−1` columns; beyond that, partition the n columns into segments
   of size in `[p, 2p−1]` and concatenate independent band constructions
   along the n axis (`buildChained` + `concatColumns`; ranks add, U is
   per-segment). Since `max(p,s) = s` on the segment range, per-segment
   formulas telescope to `⌈n(3p+1)/2⌉` whenever the ceiling slack vanishes;
   rather than hand-coding parity rules, the partition is a DP over the
   **achieved** ranks of the segment builds — which kills slack AND routes
   around degraded g ≥ 6 segment sizes automatically (`p=12, n=36` →
   16+20 = 666, not 18+18 = 668). All 196 chained shapes land at the formula.

### Case 1 — odd p = 2k+1

Internal: diagonals (2n products) + all cyclic pairs at distances `1..k`
(3 products each, n per distance) = `n(3k+2) = (3pn+n)/2` products. Band per
column = the symmetric window `[j−k, j+k]` (p rows). Back-substitute. ∎

### Case 2 — even p = 2k+2

Step 1: the odd band `1..k` as above (`n(3k+2)` products). Each column now has
`2k+1` cells and needs **one more**, at distance `k+1` (either side keeps the
window contiguous: `[j−k−1, j+k]` or `[j−k, j+k+1]`).

Step 2: a full pair `(i, i+k+1)` provides the *up* extra cell of column i AND
the *down* extra cell of column `i+k+1` — 3 products for two columns. Choosing
pairs = a perfect matching in the circulant `i ↔ i+(k+1) (mod n)`: alternate
edges along each orbit cycle of `+(k+1) mod n` (there are `g = gcd(n, k+1)`
cycles of length `n/g`).

Step 3 (boundary): each odd cycle leaves one column uncovered; its extra cell
is computed naively (2 products). Budget check: the formula allows
`⌈(3pn+n)/2⌉ − n(3k+2)` extra products = `3n/2` (n even) or `(3n+1)/2`
(n odd) — exactly `n/2` pairs, or `(n−1)/2` pairs + one naive cell. Hence the
formula is attained iff the circulant has ≤ 1 odd cycle (and n odd absorbs
that one); each additional odd cycle costs ~2 products over. The exceptional
`(n, k+1)` combos (e.g. `n=15, k+1=5`: five 3-cycles → `+2`) are exactly the
shapes where the published catalogs ALSO sit above the formula (`⟨2,10,15⟩`:
best known 234 vs formula 233) — suggesting HK's own Step 3 has structure we
(and the catalogs) haven't recovered yet.

## Verified results (final, 2026-06-12 emission)

| Family | Range swept | At formula | Exact-verified |
|---|---|---|---|
| Square `⟨2,n,n⟩` | n = 3..32 | all | all |
| Odd p band, p ≤ n ≤ 2p−1 | p ≤ 31 | all | all |
| Even p band, p ≤ n ≤ 2p−1 | p ≤ 32 | all but g ≥ 6 (6 shapes, +1..+3) | all |
| Chained, 2p ≤ n ≤ 32 | 3 ≤ p ≤ 16 | **196/196** | all |
| **Total** | 465 shapes | **459/465** | **465/465** |

**Headline: `⟨2,10,15⟩ = 233` over ℤ** — exact formula, machine-verified,
strictly below every published catalog (FMM-Lille and Perminov both hold
234). The only shapes above formula are the six `g ≥ 6` circulant cases
(`⟨2,12,18⟩+1, ⟨2,14,21⟩+1, ⟨2,16,24⟩+2, ⟨2,18,27⟩+2, ⟨2,20,30⟩+3,
⟨2,24,30⟩+1`) — provably formula-impossible within this framework (see the
arc-sum argument and the (3,3,·) theorem below), and beyond all published
catalogs regardless.

## The decisive step: bridge-3 over the TRUE reusable set

The published impossibility theorem (rank-1 atoms, char 0) rules out 3-product
`(2,2,bridge-3)` completions over `S = {Ci,Di,Cj,Dj,Ep,Fp}` — 6 products. But
the EMISSION has 12 reusables: the three diagonals AND every product of the
Lemma-2 pairs `(i,b)` and `(b,j)`, each individually weightable in W. Over
that true basis (`sympy/derive_bridge_true_reusables.py`, exact arithmetic):

- **`(2,2,bridge-3)`: SOLVABLE** — e.g. 3 new products `E(a_i+a_b−a_j)`,
  `F(a_i+a_b−a_j)`, `(a_i1+a_b2−a_j2)·(x_1b−x_1j−x_2i)`; explicit identity in
  `sympy/extract_bridge_identity.py` and implemented in
  `emitSameMethodPair_22_bridge3`.
- **`(1,1,bridge-3)`: SOLVABLE** — analogous identity, implemented.
- **`(3,3,bridge-1)` and `(3,3,bridge-2)`: NOT solvable** in a 120-atom
  candidate catalog even over the true basis — the robust gap. Consequently
  method-3 rows must stay pairwise more than k apart.
- **Upgrade (2026-06-11): the (3,3,·) gap is a theorem over ARBITRARY local
  atoms**, not just catalogs (`sympy/derive_33bridge_general.py` + `.out`).
  Block-decomposing the 6×6 monomial space (rows {a_i,a_b} vs a_j, cols
  {x_i,x_b} vs x_j) forces, for ANY 3 rank-1 atoms on rows/cols {i,b,j}: all
  three atoms full, α and β each supported on exactly 2 atoms (overlapping in
  one), and the (RB,CB) spill collapsing to two NONZERO dyads with factors in
  the virtual-row space U₀ = ⟨a_i+a_b⟩ and virtual-col space V₀ = ⟨x_i+x_b⟩,
  both required to lie in span(shared|_{RB×CB}). Exact computation:
  `(U₀⊗V₀) ∩ span(shared)` is 1-dimensional, generated by the **rank-2**
  virtual diagonal `(a_i1+a_b1)⊗(x_1i+x_1b) + (a_i2+a_b2)⊗(x_2i+x_2b)` —
  no rank-1 element, hence **no 3-product completion exists**, for either
  bridge method. Scope: the 9-product reusable set of the emitter pattern
  (bridge-sum targets); the remaining unexplored enrichment is the 12-product
  set (adding the (b,j)-pair products with free weights, which couples the
  blocks), non-local atoms, or ≥4-product completions with savings elsewhere.

This vindicates HK's "three additional multiplications" claim for the bridge-3
cases (the published impossibility holds for its narrow S, but the
construction was never confined to that S), while sharpening the genuine gap
to the `(3,3,·)` same-method case — which HK's own Lemma-3 sequence requires
and which no derivation attempt has closed.

**Triangle-family limit (arc-sum argument).** For `n = 3(k+1)` the circulant
splits into `k+1` triangles, forcing `g = k+1` leftover columns and `⌊g/2⌋`
Z-pairs, each needing one method-3 row. Three rows pairwise more than k apart
on `C_{3(k+1)}` would need their three arcs to sum to at least `3k+6 > n` —
impossible. Hence `g ≥ 6` (six such shapes in the full sweep range) degrades
by +1 per dropped Z-pair (2 naive cells instead of 3 shared products);
`g ≤ 5` closes exactly — including the `⟨2,10,15⟩` flagship (g=5: two
Z-pairs, feasible).

## HK's actual Case-2 Steps 2–3 (decoded from the paper, 2026-06-11)

With the HK71 PDF in hand (paper pages 7–13), their Case 2 differs from our
matching variant:

- **Step 2 (theirs)**: consecutive-interval distance-(k+1) pairs `(i, i+k+1)`
  for `i = 1..k+1` (and `i = 2k+3..n−k−1` when `3k+4 ≤ n`) — covering all but
  `m' = n−(2k+2)` (regime A) / `4k+4−n` (regime B) *deficient* columns
  `ℓ..ℓ+m'−1`.
- **Step 3 (theirs — the `Z`-trick)**: deficient columns are completed in
  PAIRS `(i₂, i₄)=(ℓ+2i, ℓ+1+2i)` via virtual-row aggregation: with
  `i₁ = i₂+(k+1)`, `i₃ = i₄−(k+1)`, treat `α = a_{i₂}+a_{i₃}` and
  `β = a_{i₁}+a_{i₄}` as virtual diagonal rows. Their method-j "diagonal
  products" ALREADY EXIST as the band pairs' cross-products (a Lemma-2 pair
  with methods {α,β} emits exactly the third method's products on the virtual
  row — the pattern (1,2)→E,F; (1,3)→C,D; (2,3)→A,B). The cross pair
  `Z(α,β), Z(β,α)` then costs only **3 new products** by Lemma 2, and the two
  missing cells fall out as `Z` minus already-computed cells — including
  cells like `y_{i₁,i₃}` at distance 2k+1 OUTSIDE the band, recovered as
  *rational combinations of fully-reconstructed columns* (a column with all p
  cells determines every `ĀX` entry in that column through the Lemma-1
  relation). Budget: 3 products per 2 columns ✓, one naive leftover absorbed
  by the ceiling when m' is odd ✓.

**But HK's method sequence (Lemma 3) breaks their own Lemma 2.** Their
assignment places method-3 rows at `ℓ+k+1+2i` — spacing 2 — so the band
contains **(3,3) same-method pairs**, a case the paper waves off ("several
cases exist…"). We searched for a 3-product `(3,3,bridge-1/2)` completion over
the RICH reusable set (diagonals of i, j, b + the (i,b) pair's virtual
cross-products + G — the same enrichment that solves `(2,2,bridge-3)`):
**no solution exists in an 80-atom candidate catalog** (exact Fraction
arithmetic; `references/hopcroftkerr1971/sympy/derive_33bridge.py` + `.out`).
This strengthens the known proof gap: HK's Theorem 1 as written relies on a
same-method case that has resisted every derivation attempt — and the public
catalogs sit above the formula at exactly the shapes (regime-A Case 2)
where the `Z`-trick forces clustered method-3 rows.

## The repaired Step 3 (implemented — this is what ships)

Combine the matching variant with the `Z`-trick, avoiding (3,3) pairs
(implemented in `HopcroftKerr2bc.planStep3` + `buildEvenBanded` and the
`buildEven` back-substitution):

1. Run the circulant matching (no method-3 rows at all); each odd orbit cycle
   leaves ONE leftover column, with **free choice of which vertex**.
2. Complete leftover columns in pairs via the `Z`-trick at any separation
   `δ ∈ [1, 2k+1]` (the band-distance constraints `|δ−(k+1)| ≤ k` allow it),
   placing each pair's single method-3 row so that all method-3 rows are
   pairwise more than k apart (the leftover-vertex freedom makes this a
   placement problem, not a forced clustering). One naive cell remains only
   when the leftover count g is odd — and `g odd ⟹ n odd`, absorbed by the
   ceiling. **This attains the formula for every (p, n), p ≤ n ≤ 2p−1, with
   g ≤ 5** — in particular `⟨2,10,15⟩ = 233`, below every published catalog
   (234). For `g ≥ 6` the placement is infeasible (arc-sum) and the plan
   degrades gracefully (+1 per dropped Z-pair).

## Reference emission (final: 2026-06-12, integer re-emission)

`schemes/constructed/`: **465 schemes** over the full range
`3 ≤ p ≤ 32`, `p ≤ n ≤ 32` — **459 at the exact HK formula**, 6 at +1..+3
(the g ≥ 6 family, beyond all catalogs regardless), **all over ℤ**
(fields `[F2,F3,Z,Q,R,C]`, checked from actual coefficients). Emission
history: 2026-06-11 band range (269 schemes, ℚ) with **197 strict
improvements over the union of FMM-Lille / Perminov / our prior catalog**
(margins up to ~29, e.g. `⟨2,24,27⟩ = 986` vs 1015); same-day chaining
extension (196 schemes, ALL at formula) adding **22 further strict
improvements** over everything previously held, including our own
recursive-closure constructions (e.g. `⟨2,15,32⟩ = 736` vs 741,
`⟨2,14,31⟩ = 667` vs 670); 2026-06-12 full re-emission with the
Euclidean-unimodular Lemma-1 upgrading every scheme ℚ → ℤ (identical
ranks, zero fallbacks). Each emission used the full quadruple gate; the
independent cross-language certificate reports **465/465 verified
exactly, 0 failures** (`verify_constructed_independent.py`).

## Verification policy (when do we check what)

These results beat published catalogs, so the certification story is explicit:

1. **Sweep (fast by default)** — `GenerateHk2npConstructed` sweeps the
   construction's entire applicability range
   (`3 ≤ p ≤ 32`, `p ≤ n ≤ 32`) and, by default, just BUILDS and
   WRITES — no in-sweep verification; schemes are stamped `verified: false`
   pending phase 2. The `--verify` flag opts into the full quadruple gate
   (20k-sample spot check; full residual over every tensor cell;
   exact-rational symbolic of the in-memory scheme; exact-rational symbolic
   of the **disk artifact** after a read round-trip — schemes whose exact
   coefficients exceed publishable range are then DEFERRED, never
   approximated). The reference emissions were produced WITH `--verify` —
   one-time certified registrations. Certified runs are **resumable**: a
   shape is skipped iff its on-disk file carries the IDENTICAL content hash
   (builders are deterministic) AND `verified: true` — those exact bits
   already passed the gate; `--force` redoes everything.
2. **Independent certificate (one-time, cross-language)** —
   `references/hopcroftkerr1971/verify_constructed_independent.py` re-verifies
   every published JSON with stdlib Python `Fraction` arithmetic, its own
   parser, and the documented disk conventions — sharing no code with the
   generator. Convention errors there can only produce false negatives, so a
   clean pass is an independent proof of the exact bilinear identity.
3. **After verification (steady state)** — once phase 2 passes, the schemes
   carry `verified: true` and are trusted like any other catalog entry: normal
   sweeps and consumers do NOT re-verify them. They participate in the
   standard catalog-wide verification flows (`VerifyAllSchemes`) like
   everything else. "Valid by construction" applies to *uses* of these
   schemes (Kron/concat/recombination of verified inputs), not to the
   emitter's own output — the emitter is new math, and this session caught
   three real bugs in it through exactly these gates.

## Where we stopped (2026-06-12) — how to resume

Everything below the formula bar is DONE and shipped: band construction
(tasks #7/#8), chaining (`buildChained`, task #10: DP over achieved segment
ranks — kills ceiling slack and routes around g ≥ 6 segments, e.g.
`p=12, n=36` → 16+20 = 666, not 18+18 = 668), and the ℤ upgrade
(Euclidean-unimodular Lemma-1, task #11). 465 schemes registered, dual
certified, paper section in sync. **One open item remains, plus two parked.**

1. **OPEN — task #9: the g ≥ 6 family, 6 shapes at +1..+3 over formula**:
   `⟨2,12,18⟩=334 (+1)`, `⟨2,14,21⟩=453 (+1)`, `⟨2,16,24⟩=590 (+2)`,
   `⟨2,18,27⟩=745 (+2)`, `⟨2,20,30⟩=918 (+3)`, `⟨2,24,30⟩=1096 (+1)`.
   These are formula-impossible within the current framework, by two
   stacked results: (a) the arc-sum argument — g ≥ 6 needs ≥ 3 method-3
   rows pairwise ≥ k+2 apart, geometrically impossible; (b) the
   `(3,3,bridge-1/2)` THEOREM (`sympy/derive_33bridge_general.py`): no
   3-product completion by ANY rank-1 atoms local to {i,b,j} exists over
   the emitter's 9-product reusable set — the would-be spill must be a
   rank-1 element of `(U₀⊗V₀) ∩ span(shared)`, which is 1-dimensional and
   generated by the RANK-2 virtual diagonal `(a_i+a_b)⊗(x_i+x_b)`.

   **Resume points, in order of promise:**
   - **12-product reusable set** (add the (b,j)-pair's three products with
     free weights; targets become raw `y_ij, y_ji`). Checked so far: the
     block forced-structure argument PARTIALLY survives — projecting out
     the a_b rows / x_b cols kills the shared set on the target blocks, so
     "all atoms full + α,β on exactly 2 atoms" still holds for the
     *projected* dyads — but the spill-collapse weakens to projected
     proportionality (free a_b/x_b components survive). What remains is a
     finite quadratic system after the structural substitution:
     Gröbner/elimination it. A solution closes all six shapes; another
     impossibility makes the +1 framework-optimal in the strongest local
     sense.
   - **Non-local atoms** (rows/cols beyond {i,b,j}) — no analysis yet.
   - **4-product (3,3) completions traded against a saving elsewhere** —
     budget bookkeeping says one saved product anywhere pays for one +1.
   - Note: all six shapes are already BELOW every published catalog — this
     gap is against our own formula target, not against the field.

2. **PARKED — lower bounds**: everything we emit is upper-bound tier; HK
   proved optimality only for `⟨2,2,n⟩` and `⟨2,3,3⟩`. Explicit user
   decision to skip for now.

3. **PARKED — housekeeping**: `references/hopcroftkerr1971/README.md` is
   the historical impossibility document (kept, with a status banner); the
   2026-06 working tree (schemes, manifest, paper, tests) was pending
   commit when we stopped.
