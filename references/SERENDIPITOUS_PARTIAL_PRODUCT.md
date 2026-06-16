# Serendipitous tensor product — Schönhage borrow-and-correct (#159)

> **CORRECTION (2026-06-05) — most of the model below is WRONG; read this first.**
> Perminov's draft (*Meta Flip Graph meets Serendipitous Product*, Def 2.9–2.12,
> §2.6) defines the serendipitous product **exactly and combinatorially** — NOT
> as subset-deletion (§3 here) and NOT as border rank (§4.6, `BORDER_RANK_LAYER.md`).
> The real mechanism:
> - A **bud** = a pair/group of rank-one terms sharing the SAME `u` vector (or
>   `v`, or `w`), where the shared vector is ANY linear combination (up to
>   scaling), not a basis vector.
> - Buds decompose a scheme into **elementary matmul tensors**:
>   `T = Σ Sᵢ⟨Nᵢ,Mᵢ,Pᵢ⟩`, `r = Σ Sᵢ·NᵢMᵢPᵢ`. `U`-bud→`⟨1,1,2⟩`, k-way→`⟨1,1,k⟩`;
>   `V`-bud→`⟨2,1,1⟩`; `W`-bud→`⟨1,2,1⟩`; lone term→`⟨1,1,1⟩`.
> - **Serendipitous product** (**Smith 2002, eq. (69)**; re-derived Perminov
>   Def 2.12 / Sedoglavic): `T₁ ⊗ ⟨n₂,m₂,p₂⟩ = Σ Sᵢ⟨Nᵢn₂,Mᵢm₂,Pᵢp₂⟩`, each
>   block at its BEST known rank ⇒ `r_s = Σ Sᵢ·R(Nᵢn₂,Mᵢm₂,Pᵢp₂)`.
>   This is the **same objective as the recombination/block-split search** —
>   `min Σ R(sub-shapes)` over a decomposition, given a SOTA oracle `R(·)`. The
>   recombination decomposition is a freely-chosen ALLOCATION (B&B-searched);
>   the serendipitous decomposition is the base's BUD STRUCTURE (non-unique —
>   we take the greedy one, so `r_s` is an upper bound). An optimal-bud-structure
>   search would be the dual of the allocation search.
>
> **Verified:** HK `⟨2,3,4⟩=20` has 4 proportional `U`-buds (`{1,19},{5,17},
> {8,16},{10,15}`) ⇒ `12·⟨1,1,1⟩ + 4·⟨1,1,2⟩`; `⊗⟨2,4,2⟩` →
> `12·⟨2,4,2⟩@14 + 4·⟨2,4,4⟩@26 = 272`. EXACT, combinatorial, no border rank.
>
> **Why my recognizers found nothing:** they searched for axis-aligned matmul
> *subtensors* (`u`=basis vector); buds need only `u^(i) ∝ u^(j)` as vectors.
>
> **Correct engine (replaces §1–§5 below):** (1) bud recognizer — group the r
> terms by proportional `u`/`v`/`w` columns → elementary `⟨N,M,P⟩` blocks;
> (2) §2.6 bud-block constructor (shared `u₁⊗u₂`, sum `v`/`w` across the k
> sub-blocks — Perminov Figs 1–3); (3) `Verifier.isExactNonCubic`. No ALS/SAT/ε.
> `SerendipitousProduct.build`'s contiguous-offset embedding is too restrictive
> (buds are non-axis-aligned) and needs the general bud-block construction.
> `BorrowBlockRecognizer`/`GeneralBorrowRecognizer` search the wrong object —
> replace with the bud recognizer. The text below is kept only as a record of
> the wrong turns.

---

Design note (2026-06-04, SUPERSEDED — see correction above). The constructive
engine for the "serendipitous tensor product" composition. Companion to
`KnownTauIdentities` (rank-prediction tier) and `references/MULTISET_FRONTIER.md`.

## 0. Provenance — say exactly what this is, and is NOT

- **Origin**: Schönhage 1981, *Partial and total matrix multiplication*
  (SIAM J. Comput. 10). Textbook: **BCS97 §14–16** (ref [65]). The
  asymptotic engine behind the τ-theorem / asymptotic-sum inequality, here
  used at **fixed small size**.
- **FMM naming**: "serendipitous tensor product" is Sedoglavic's (FMM-Lille)
  label. It appears verbatim on `fmm.univ-lille.fr/4x8x12.html`:
  `⟨4×8×12:272⟩ = (⟨2×3×4:20⟩ − 8) ⊗ ⟨2×4×2:14⟩ + 4⟨2×4×4:26⟩`.
- **NOT DIS09.** Drevet–Islam–Schost 2009 (ref [10]) is *sparseness-aware
  recursive padding* (prune zero/unneeded sub-products of a padded scheme),
  which the repo already implements (`Recombination` peel #87 + output-side
  zero masks #86). My earlier "DIS09 use partial products" was wrong —
  retracted. The two techniques are different: DIS09 *subtracts* useless work
  from an exact padded scheme; serendipitous *borrows* (deletes good work) and
  *adds back* a cheaper correction.

## 1. The identity

A bilinear algorithm `A` for `⟨a,b,c⟩` is the matmul tensor as a rank-`rA`
sum `T_A = Σ_{ℓ=1}^{rA} uℓ ⊗ vℓ ⊗ wℓ`. Kronecker of algorithms realises
`T_⟨a,b,c⟩ ⊗ T_⟨d,e,f⟩ = T_⟨ad,be,cf⟩` at rank `rA·rB`.

**Borrow-and-correct.** Pick a subset `S ⊆ {1..rA}`, `|S| = k`. Split
`T_A = T_A^{kept} + P_S` where `P_S = Σ_{ℓ∈S} uℓ⊗vℓ⊗wℓ` (rank ≤ k). Then for
any full `B = T_⟨d,e,f⟩`:

```
T_target = T_A ⊗ T_B
         = T_A^{kept} ⊗ T_B   +   P_S ⊗ T_B
           └── (rA−k)·rB ──┘     └─ the DEFECT D ─┘
```

This is an **exact identity, always true** — the content is purely the COST.
The construction wins iff the defect `D = P_S ⊗ T_B` admits a decomposition
cheaper than its naive `k·rB`:

```
R(target) ≤ (rA − k)·rB + cost(D),   serendipitous ⇔ cost(D) < k·rB.
```

The "serendipity": this only happens for special `(A, S, B)`. The cheapest
`cost(D)` arises when **`D` is a direct sum of matmul tensors**.

### 1.1 The key simplification

If `P_S` is itself a **direct sum of matmul tensors** in disjoint index blocks,
`P_S = ⊕ⱼ T_⟨aⱼ,bⱼ,cⱼ⟩`, then

```
D = P_S ⊗ T_B = ⊕ⱼ T_⟨aⱼd, bⱼe, cⱼf⟩      (direct sum of matmuls)
cost(D) = Σⱼ R(⟨aⱼd, bⱼe, cⱼf⟩).
```

So the whole construction reduces to: **find a subset of `A`'s products whose
removed part `P_S` is a direct sum of (smaller) matmul tensors.** That is the
"borrow block". The 4×8×12 instance: `A=⟨2,3,4⟩:20`, `S` with `k=8`,
`P_S = T_⟨2,?,?⟩`-block s.t. `P_S ⊗ ⟨2,4,2⟩ = 4·⟨2,4,4⟩` ⇒
`(20−8)·14 + 4·26 = 168 + 104 = 272 < 280`. (`⟨6,8,9⟩=296` is the same recipe
with `B=⟨3,3,2⟩`, aux `⟨3,3,4⟩` — already in `KnownTauIdentities.FMM_6_8_9_296`.)

## 2. The recognizer (the hard core)

**Problem.** Given the rank-`rA` decomposition `{uℓ⊗vℓ⊗wℓ}` of `A`, find a
subset `S` such that `P_S = Σ_{ℓ∈S} uℓ⊗vℓ⊗wℓ` is a direct sum of matmul
tensors, maximising the saving `k·rB − Σⱼ R(⟨aⱼd,bⱼe,cⱼf⟩)`.

**Why hard.** Recognising whether an arbitrary 3-tensor is (a direct sum of)
matmul tensors *up to basis change* is a structural decomposition problem
(matmul tensors are characterised, but recognition up to the GLₐ×GL_b×GL_c
action is nontrivial). Brute force over `2^{rA}` subsets is infeasible for
`rA ≥ 20`.

**Tractable restriction (do this first).** Require the blocks to be
**axis-aligned** in `A`'s natural index basis: a partition of the row index
`[a]`, mid index `[b]`, col index `[c]` into blocks, with `P_S` supported on a
union of `(rowblk × midblk × colblk)` cells that each carry a *complete* matmul
tensor. FMM's instances are axis-aligned in the natural Kronecker basis, so
this restriction is expected to cover them. Algorithm:

1. For each candidate sub-shape block `⟨a',b',c'⟩` aligned to index sub-ranges
   of `⟨a,b,c⟩`, compute the matmul subtensor `T_block` it would carry.
2. Test whether `T_block` is **spanned by a subset of `A`'s rank-1 terms**
   exactly (solve: is `T_block = Σ_{ℓ∈S'} uℓ⊗vℓ⊗wℓ` for some `S'`? — a
   support/linear-membership test over the term set).
3. Greedily/exactly select disjoint blocks maximising the saving; emit.

### 2.1 Worked structure of 4×8×12 (the recognizer's concrete target)

With `B = ⟨2,4,2⟩` ⇒ `(d,e,f) = (2,4,2)`. The correction `⟨2,4,4⟩ = ⟨aⱼd,bⱼe,cⱼf⟩`
forces `(aⱼ,bⱼ,cⱼ) = (1,1,2)` (since `1·2=2, 1·4=4, 2·2=4`). Multiplicity 4 ⇒
**the borrow is 4 disjoint `⟨1,1,2⟩` sub-matmuls of `⟨2,3,4⟩`**. Each `⟨1,1,2⟩`
has rank 2, so `k = 4·2 = 8` ✓, and each completes to `⟨1,1,2⟩⊗⟨2,4,2⟩ = ⟨2,4,4⟩`
(rank 26). Saving `= k·rB − Σ R = 8·14 − 4·26 = 112 − 104 = 8` ✓.

So the recognizer's job here is concrete and small: in a `⟨2,3,4⟩=20` base, find a
subset `S` of 8 terms that splits into 4 disjoint `⟨1,1,2⟩` matmul subtensors
(one `a`, one `b`, two `c`'s, each `C[a,c]=A[a,b]·B[b,c]`). **Caveat (the real
serendipity):** whether such clean sub-blocks exist depends on the SPECIFIC
`⟨2,3,4⟩` decomposition — a generic rank-20 scheme need not contain them. So the
search ranges over `⟨2,3,4⟩` base schemes (and their orbit) too; finding one with
the removable blocks IS the discovery. Self-certified by the symbolic verifier,
so a wrong guess fails verification rather than lying.

### 2.2 Self-certification (no FMM import — approach ii)

The constructor emits a standard flatten `(U,V,W)` `NonCubicBilinearAlgorithm`
(partial-A⊗B columns ++ offset-embedded correction columns). `Verifier.
isExactNonCubic` is an **independent oracle**: any recognizer mistake yields a
scheme that fails verification, never a silent wrong answer. The scheme is also
an SLP (R columns = R multiplications + U/V/W linear-combination adds), so the
#190 additive-complexity machinery applies directly. Ground-truth FMM data is
NOT needed.

**Honesty tier (CLAUDE.md optimality discipline).** The axis-aligned
recognizer is *optimal-within-scope* (axis-aligned blocks only); a full
up-to-basis recognizer would be needed to call any negative result "no
serendipitous form exists". Label outputs `optimality: "optimal-within:
axis-aligned-borrow-blocks"`.

## 3. Construct + verify layer (the unambiguous foundation — build first)

`SerendipitousProduct.build(A, S, B, List<Correction>)`:
- **Partial A**: drop columns `S` from `A.U/V/W` → `(rA−k)`-term factor set.
- **Kronecker** with `B` (reuse `Compose.kroneckerGeneral`) → a
  `(rA−k)·rB`-term partial scheme for `⟨ad,be,cf⟩` (computes `target − D`).
- **Corrections**: each `Correction` = (matmul scheme `Cⱼ` for
  `⟨aⱼd,bⱼe,cⱼf⟩`, row/mid/col block offset). Append its columns, mapped into
  the target `U/V/W` index space at the block offset (disjoint ⇒ no overlap).
- **Verify**: `Verifier.isExactNonCubic` on the assembled scheme. Self-check
  before any persistence.

This layer is testable independent of the recognizer: with `S = ∅` it must
reproduce `Compose.kroneckerGeneral(A,B)` exactly (regression anchor).

## 4. Lineage

New node (or reuse) for replay:
`BorrowAndCorrect(baseA, removedColumns[], baseB, corrections[(ref, blockOffset)])`.
Renders compact as `(A ⊖ k) ⊗ B ⊕ Σ cⱼ·Cⱼ`. Promotes the existing
`KnownTauIdentities.BorrowAndCorrect` **predictor** record into a **replayable
constructor** — the predicted rank and the constructed rank must agree (a
divergence is a catalog inconsistency, same contract as the
DisjointSum/Sedoglavic pair).

## 4.5 Findings (2026-06-04) — engine built, base is the bottleneck

The full engine is implemented and self-certifying:
- `SerendipitousProduct.build` (construct+verify) — `TestSerendipitousProduct` 4✓.
- `BorrowBlockRecognizer` (axis-aligned *confined* recognizer) — finds blocks via
  `Verifier.isExactNonCubic` on the extracted sub-scheme; on `naive ⟨2,3,4⟩` it
  rediscovers a verified `⟨4,12,8⟩=312 < 336` (zero FMM import).
- `GeneralBorrowRecognizer` (non-confined: solves `Σ_{ℓ∈S} tℓ = T_block` by
  least squares over the flattened tensor space; columns of a valid base are
  generically independent ⇒ unique 0/1 subset). Validated: finds the 4 `⟨1,1,2⟩`
  blocks of `naive ⟨2,2,2⟩`.

**Exact negative on our rank-20 bases.** Both `alphatensor_Z-2x3x4_m20` and
`alphatensor_F2-2x3x4_m20` contain **zero** borrowable matmul sub-blocks of ANY
shape (full scan a1×b1×c1). Their 20 products are fully entangled — no proper
subset equals an axis-aligned matmul subtensor. So **FMM's serendipitous
`⟨4,8,12⟩=272` is unreachable from our bases**: it requires FMM's specific
`⟨2,3,4⟩=20` basis whose products isolate `⟨1,1,2⟩` slices. The recognizer is
complete and correct; the bottleneck is now **base availability**, a separate
(harder) base-construction problem:
- (a) import FMM's `⟨2,3,4⟩=20` atom (one base, not the result) → recognizer
  rediscovers 272 from it;
- (b) search/synthesise a rank-20 `⟨2,3,4⟩` whose product set carries 4
  borrowable `⟨1,1,2⟩` blocks (ALS/SAT with block/sparsity constraints);
- (c) wire the engine into `BlockSplitSearch` so it opportunistically fires on
  any (target, base) where the catalog base DOES have borrowable blocks.

## 4.6 Boundary finding (2026-06-05) — the borrow is NOT column deletion

Fetched FMM's own `⟨2,3,4⟩=20` Maple tensor (`2x3x4_tensor.mpl`, attributed to
**Hopcroft-Kerr 1971**) and parsed it (`isExact=true`, so the parse is correct).
The general recognizer finds **zero** borrowable sub-blocks of ANY shape in it —
same as AT-Z/AT-F2.

This refutes the "delete a 0/1 subset" model of the serendipitous borrow. The
FMM identity `(⟨2,3,4⟩−8)⊗⟨2,4,2⟩ + 4⟨2,4,4⟩` forces (cancel `⊗T_B`)
`P_S = ⊕⁴ T_⟨1,1,2⟩`; but no subset of HK's 20 products sums to even one
`⟨1,1,2⟩`. So FMM's "`−8`" is a genuine **partial / border-rank** construction
(a degenerate ε-limit tensor in Schönhage/Bini's sense), NOT a sub-scheme of the
exact 20-product algorithm. Our subset-deletion engine (confined OR general)
cannot reach it.

**Consequence for #159.** The engine is correct and complete for the
*subset-deletion borrow* class (validated end-to-end on naive bases). Rediscovering
FMM's published serendipitous results requires modelling **partial matrix
multiplication / border rank** (Bini 1979, Schönhage 1981) — degenerate tensors
with an ε-parameter and a completion step — which is a substantially deeper
undertaking than column deletion. Tracked as the next layer, not closed here.
Files fetched to `target/fmm-maple/` (uncommitted).

## 5. Build order

1. **Construct+verify layer** (§3) + regression test (`S=∅` ≡ Kronecker).
2. **Reconstruct 4×8×12 = 272** end-to-end — needs the explicit `(A, S, B,
   corrections)` ground truth. Either (a) run the §2 recognizer on a
   `⟨2,3,4⟩=20` base to *discover* `S`, or (b) one targeted fetch of the
   `4x8x12` FMM Maple file for the published decomposition (a single page, not
   the full scan that was declined). Validating instance for the whole engine.
3. **Axis-aligned recognizer** (§2) + unit tests on synthetic `⊕ T_⟨a,b,c⟩`.
4. **Wire into `BlockSplitSearch.findBestStrategy`** as a candidate strategy
   (alongside Kron / Concat / recombination / `KnownTauIdentities`), gated by
   the saving criterion; emit via the §4 lineage.
5. **REFERENCES #195 tags** (already started): DIS09 = *sparseness-aware
   recursive padding*; this technique = *partial matrix multiplication
   (Schönhage borrow-and-correct)*; Schönhage 1981 + BCS97 §14–16 as sources.
