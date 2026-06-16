# Rosowski 2019/2020 — extracted algorithms

Source: `references/papers/rosowski_2019_commutative_matmul_arxiv1904.07683.pdf`
(arXiv:1904.07683, A. Rosowski, *Fast Commutative Matrix Algorithm*).

## Why these aren't in `src/main/resources/schemes/`

Rosowski's algorithms are **non-bilinear**: their products mix A and B
entries in both factors (e.g. `(a12 + b12)(a11 + b21)`). Our standard
`SchemeIO` scheme format requires strictly bilinear products
`(linear in A) × (linear in B)`, so these algorithms don't round-trip
through that format.

Two options to fix later:
- **Convert to bilinear**: expand each mixed product into 4 standard
  bilinear products. This inflates rank (e.g. ⟨3,3,3⟩=21 non-bilinear
  → ~80 bilinear). Loses the point.
- **Extend SchemeIO**: add a non-bilinear scheme JSON format with
  per-product `(u_a, u_b, v_a, v_b)` coefficient vectors. Recursive
  application would require special handling.

Until then this Markdown is the authoritative source.

---

## Algorithm 1 (Theorem 1) — `⟨n, 3, 3⟩` in `6n + 3` mults

Input: an `n × 3` vector/matrix `A` (here shown for a single row
`a = (a1, a2, a3)`) and a `3 × 3` matrix `B = (b_{ij})`.

```
p1 := (a2 + b12)(a1 + b21)
p2 := (a3 + b13)(a1 + b31)
p3 := (a3 + b23)(a2 + b32)
p4 := a1 · (b11 − b12 − b13 − a2 − a3)
p5 := a2 · (b22 − b21 − b23 − a1 − a3)
p6 := a3 · (b33 − b31 − b32 − a1 − a2)
p7 := b12 · b21
p8 := b13 · b31
p9 := b23 · b32

Output (the row a · B):
  (aB)_1 = p4 + p1 + p2 − p7 − p8
  (aB)_2 = p5 + p1 + p3 − p7 − p9
  (aB)_3 = p6 + p2 + p3 − p8 − p9
```

The 3 products `p7, p8, p9` involve only B entries and **can be reused
across all `n` rows of A**. So the full `⟨n,3,3⟩` cost is
`6 + 3·(n − 1) = 6n + 3`. For `n = 3` this gives the 21-mult
`⟨3,3,3⟩` algorithm below.

## Corollary 1 — `⟨3,3,3⟩` in 21 multiplications

For `A = (a_{ij})_{3×3}` and `B = (b_{ij})_{3×3}`:

```
p1  := (a12 + b12)(a11 + b21)
p2  := (a13 + b13)(a11 + b31)
p3  := (a13 + b23)(a12 + b32)
p4  := a11 · (b11 − b12 − b13 − a12 − a13)
p5  := a12 · (b22 − b21 − b23 − a11 − a13)
p6  := a13 · (b33 − b31 − b32 − a11 − a12)

p7  := (a22 + b12)(a21 + b21)
p8  := (a23 + b13)(a21 + b31)
p9  := (a23 + b23)(a22 + b32)
p10 := a21 · (b11 − b12 − b13 − a22 − a23)
p11 := a22 · (b22 − b21 − b23 − a21 − a23)
p12 := a23 · (b33 − b31 − b32 − a21 − a22)

p13 := (a32 + b12)(a31 + b21)
p14 := (a33 + b13)(a31 + b31)
p15 := (a33 + b23)(a32 + b32)
p16 := a31 · (b11 − b12 − b13 − a32 − a33)
p17 := a32 · (b22 − b21 − b23 − a31 − a33)
p18 := a33 · (b33 − b31 − b32 − a31 − a32)

p19 := b12 · b21
p20 := b13 · b31
p21 := b23 · b32

(AB)_{i,1} = p_{6i−2} + p_{6i−5} + p_{6i−4} − p19 − p20   (i = 1, 2, 3)
(AB)_{i,2} = p_{6i−1} + p_{6i−5} + p_{6i−3} − p19 − p21
(AB)_{i,3} = p_{6i}   + p_{6i−4} + p_{6i−3} − p20 − p21

Concretely:
AB = | p4 + p1 + p2  − p19 − p20    p5 + p1 + p3  − p19 − p21    p6 + p2 + p3  − p20 − p21 |
     | p10+ p7 + p8  − p19 − p20    p11+ p7 + p9  − p19 − p21    p12+ p8 + p9  − p20 − p21 |
     | p16+ p13+ p14 − p19 − p20    p17+ p13+ p15 − p19 − p21    p18+ p14+ p15 − p20 − p21 |
```

Total products: **21**. Notable: 6 of them (p4, p5, p6, p10, p11, p12,
p16, p17, p18) are of the form `a_{ij} · (linear in B AND A)` — the
non-bilinearity. The remaining 12 are products of two mixed
`(A entry + B entry)` linear forms (p1, p2, p3, p7, p8, p9, p13, p14,
p15) or pure B-only (p19, p20, p21).

## Theorem 2 — divisions-free `⟨l, n, m⟩` for even `n` (`n(lm+l+m−1)/2`)

For an `l × n` by `n × m` product over a commutative ring with **even**
contraction dimension `n`, Theorem 2 computes `AB` in
`n(lm + l + m − 1)/2` multiplications **without any divisions**.

The rank `n(lm+l+m−1)/2` is **not** Rosowski's — it was already known
(Waksman 1970 [19], Islam 2009 [10]), but those constructions divide by
2 (so they need 2 invertible: Q/R/C, not Z). Rosowski's contribution is
the *divisions-free* construction of the same rank, valid over **any**
commutative ring including `Z`, `F₂`, `F₃`.

Formulas (1-indexed, `k = 1..n/2`). For `i = 1..l`:

```
c_{i,1} = Σ_k a_{i,2k−1}(b_{2k−1,1} + a_{i,2k})
        + Σ_k a_{i,2k}  (b_{2k,1}   − a_{i,2k−1})
```

For `i = 1..l`, `j = 2..m`:

```
c_{i,j} = Σ_k (a_{i,2k−1} + b_{2k,j})(a_{i,2k} + b_{2k−1,1} + b_{2k−1,j})
        − Σ_k a_{i,2k−1}(b_{2k−1,1} + a_{i,2k})
        − Σ_k b_{2k,j}  (b_{2k−1,1} + b_{2k−1,j})
```

Four product families: `P1 = a(b+a)` (`l·n/2`), `P2 = a(b−a)` (`l·n/2`),
`S = b(b+b)` (B-only, `(m−1)·n/2`, shared across rows), and
`Q = (a+b)(a+a+b)` (`l(m−1)·n/2`); total `n(lm+l+m−1)/2`. The
`P1` products are reused in both `c_{i,1}` (added) and every `c_{i,j≥2}`
(subtracted).

**Specialised to `l = n = 2`, `m = p`** this is the `⟨2,2,p⟩ = 3p+1`
family. Implemented in
[`RosowskiTheorem2`](../src/main/java/eu/solven/matmul/papers/rosowski2019/RosowskiTheorem2.java)
(general even-contraction `build(l, n, m)` + `build22p(p)` shortcut),
materialised for `p = 3..16` by `MaterializeRosowskiTheorem2`
(field `Z`, `commutative: true`, `attribution_for_rank: Waksman 1970`,
`discovery: false`). These DO round-trip through `SchemeIO` via the
non-bilinear `scheme_type` (`SparseNonBilinearWriter`), unlike the
"TODO — bilinear conversion" note below which is about a *bilinear*
re-encoding (still pointless).

## Theorem 3 — odd-contraction case

Theorem 3 handles odd `n ≥ 3` (with `m ≥ 3`) by splitting
`A = [A₁ | A₂]`, `B = [B₁; B₂]` with `A₁ ∈ R^{l×3}`, computing `A₂B₂`
via Theorem 2 (the remaining even contraction `n−3`) and `A₁B₁` via an
Algorithm-1-style ⟨l,3,m⟩ block. Cost `n(lm+l+m−1)/2` for odd `m`, and
`(n(lm+l+m−1)+l−1)/2` for even `m`. Bound exposed via
`RosowskiBound.commutativeBoundBilinear`; the explicit ⟨2,2,p⟩ family is
even-contraction so it is fully covered by Theorem 2 above (no Theorem 3
constructor needed for that family).

## Theorem 4 / Theorem 5 — cubic `⟨n,n,n⟩`

Theorem 4 (n even) and Theorem 5 (n odd) give closed-form rank bounds
`n(n² + 3n + 1)/2` and `n(n² + 3n + 2)/2` respectively, via explicit
inductive construction over the n=3 base case (Algorithm 1). The
formulas are encoded in
[`RosowskiBound.nonBilinearRankBound`](../src/main/java/io/cormoran/strassen/v3/catalog/RosowskiBound.java);
the explicit per-n products are in Sections 3.1 and 3.2 of the paper.

**Why not extracted product-by-product here**: the Theorem 4/5
constructions are inductive (each n derives from n−2 via a structured
expansion of the n=3 base), so the explicit product list grows with n.
The recipe is mechanical; concrete products at large n are best
generated by code following Sections 3.1–3.2, not hand-copied.

## TODO — bilinear conversion

To express Rosowski's `⟨3,3,3⟩=21` algorithm as a standard bilinear
scheme, each non-bilinear product `(A_combo)(B_combo + A_combo')` must
be expanded into the 4 bilinear products it implies (after
distributing). Net cost: ~84 bilinear products — strictly worse than
Laderman's 23. So a "bilinearised Rosowski" wouldn't be useful as a
scheme entry. The non-bilinear formulation is the whole point.
