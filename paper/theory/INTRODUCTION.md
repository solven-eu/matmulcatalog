# Introduction — what is this project about, in plain English?

This repository studies **fast algorithms for multiplying two matrices**.
The headline finding of the field, from 1969, is that you can multiply
two `2×2` matrices using **7 multiplications** instead of the obvious
**8** — that's Strassen's algorithm. The repo catalogs every such
"trick" discovered since, lets you compose them recursively to handle
bigger matrices, and runs searches for new ones.

If you've never thought about this, the rest of this doc walks through
the prerequisites from scratch.

---

## What is a matrix?

A **matrix** is a rectangular grid of numbers. We write its size as
`m × n` ("m by n"): `m` rows, `n` columns.

```
A = [ 1  2  3 ]       ← a 2×3 matrix
    [ 4  5  6 ]

B = [ 7  ]            ← a 3×1 matrix
    [ 8  ]
    [ 9  ]
```

Matrices show up everywhere — image rotations, neural-network layers,
solving systems of equations, computer graphics, scientific simulation.

---

## What is matrix multiplication?

Given an `m × n` matrix `A` and an `n × p` matrix `B` (note: A's
columns must match B's rows), their product `C = A · B` is an `m × p`
matrix where:

```
C[i, j] = Σ A[i, k] · B[k, j]    over k from 1 to n
```

In English: **the (i, j) entry of C is the dot product of A's i-th row
with B's j-th column**.

Example with the matrices above:
```
C = A · B = [ 1·7 + 2·8 + 3·9 ]   =   [ 50 ]
            [ 4·7 + 5·8 + 6·9 ]       [ 122 ]
```

The "size" of a matmul problem is written `⟨m, n, p⟩` throughout this
repo — it's the **shape triple**.

### Why is matmul interesting?

It's the **fundamental bottleneck** in linear algebra. Neural network
training, image rotation, solving linear systems, eigenvalue
computation — all reduce to or are dominated by matmuls.

If you make matmul 1% faster, the entire deep-learning industry runs
1% faster.

---

## The naive algorithm

The definition above IS an algorithm: for each entry `C[i,j]`, compute
its dot product directly. For `n × n × n` square matmul this costs
**`n³` scalar multiplications** (and roughly `n³` additions).

```
for i in 1..m:
  for j in 1..p:
    C[i,j] = 0
    for k in 1..n:
      C[i,j] += A[i,k] · B[k,j]      ← n multiplications per output entry
```

For `2 × 2 × 2` (multiplying two `2×2` matrices) this is `2³ = 8`
multiplications.

---

## Strassen's discovery (1969)

Volker Strassen, a German mathematician, published a paper titled
*"Gaussian elimination is not optimal"* in 1969. He showed that
multiplying two `2×2` matrices needs only **7 multiplications** — not 8.

The trick is a clever rewriting:
```
M₁ = (A₁₁ + A₂₂)(B₁₁ + B₂₂)        ← 1 mult of two sums
M₂ = (A₂₁ + A₂₂)·B₁₁
M₃ = A₁₁·(B₁₂ − B₂₂)
M₄ = A₂₂·(B₂₁ − B₁₁)
M₅ = (A₁₁ + A₁₂)·B₂₂
M₆ = (A₂₁ − A₁₁)(B₁₁ + B₁₂)
M₇ = (A₁₂ − A₂₂)(B₂₁ + B₂₂)

C₁₁ = M₁ + M₄ − M₅ + M₇
C₁₂ = M₃ + M₅
C₂₁ = M₂ + M₄
C₂₂ = M₁ − M₂ + M₃ + M₆
```

You can verify by expanding: every term matches the naive formula.
**7 multiplications, 18 additions/subtractions.** The naive needs 8
multiplications and 4 additions. We traded 4 extra additions for 1
saved multiplication.

### Why does that matter?

For just two `2×2` matrices, saving one multiplication is silly. The
power is **recursion**: Strassen's trick works on `2×2` matrices whose
**entries are themselves matrices** (matrix blocks). So you can use it
to multiply huge matrices, applying the trick at every level.

A `2^k × 2^k` matrix multiplied via Strassen-recursion costs **`7^k`
multiplications** instead of the naive `8^k = (2^k)³`. For `k=10`
(so `1024 × 1024` matrices), naive is ~10⁹ mults; Strassen is ~3·10⁸ —
a 3× speedup. For larger matrices the asymptotic gain is even better.

This is the discovery that founded an entire research area:
**fast matrix multiplication**.

---

## Divide and Conquer

Strassen's algorithm is a **divide-and-conquer** algorithm — a
problem-solving pattern where you:

1. **Divide** the problem into smaller sub-problems of the same kind.
2. **Conquer** each sub-problem recursively (or directly if small enough).
3. **Combine** the results.

For matmul: split each `2n × 2n` matrix into four `n × n` blocks; apply
Strassen at the top level (7 sub-products of `n × n` matrices instead
of 8); recursively solve each sub-product the same way.

This pattern is why a "trick" saving 1 mult at the top level cascades
into massive savings at scale: each level of recursion multiplies the
savings.

### Why does commutativity matter for divide-and-conquer here?

When we recurse, the "scalar multiplications" become **block matrix
multiplications**. The matrix blocks generally **don't commute**
(`A · B ≠ B · A` in general — see the next section).

If Strassen's algorithm relied on `a · b = b · a`, applying it
recursively to blocks would give wrong answers. Strassen's `M₁` uses
`(A₁₁ + A₂₂)(B₁₁ + B₂₂)` — the order is fixed; it doesn't swap. The
algorithm works **regardless of whether the entries commute**, so it
recurses safely.

Algorithms that DO use commutativity (we call these "commutative
schemes") only work when the matrix ENTRIES are scalars (like real
numbers). They CAN'T be recursed into block matmuls. That distinction
is load-bearing across this repo — see the "Algebra" section below.

---

## Numbers live in fields, rings, and algebras

So far we've assumed entries are "numbers". But "number" comes in
flavors:

### Sets of numbers

- **ℤ** — the **integers**: `..., −2, −1, 0, 1, 2, ...`
- **ℚ** — the **rationals**: fractions of integers (`1/2`, `−7/3`, `8/1`)
- **ℝ** — the **reals**: all decimal numbers (rationals + irrationals
  like `π`, `√2`)
- **ℂ** — the **complex numbers**: pairs `a + b·i` where `i² = −1`
- **𝔽₂** — the **field of two elements**: just `{0, 1}` with arithmetic
  done modulo 2 (`1 + 1 = 0`, `1 · 1 = 1`). Used in cryptography,
  coding theory, hardware.
- **𝔽ₚ** — more generally, integers modulo a prime `p`
- **𝔽₃** = `{0, 1, 2}` modulo 3, etc.

### What are these algebraically?

A **ring** is a set where you can add, subtract, and multiply
sensibly (with the usual distributive law). Examples: ℤ, ℚ, ℝ, ℂ,
the set of `n × n` matrices over any ring.

A **field** is a ring where you can ALSO divide by non-zero elements.
Examples: ℚ, ℝ, ℂ, 𝔽₂, 𝔽ₚ. **Counter-example**: ℤ is NOT a field
(you can't divide 1 by 2 and stay in ℤ). The set of `n×n` matrices is
NOT a field (not every matrix has an inverse).

An **algebra** is a ring whose elements you can also multiply by
"scalars" from some underlying field — vectors of polynomials, for
instance.

The TLDR for this project: **the rank of a matmul algorithm depends
on which field the matrix entries live in**. AlphaTensor 2022 found
that `4×4` matmul over 𝔽₂ needs only **47 multiplications**, but over
ℝ the best known is still Strassen-recursion's **49**. AlphaEvolve
2025 then found an algorithm using only **48** multiplications over
ℂ. Same problem, different fields, different answers.

---

## Commutativity

Two numbers `a` and `b` **commute** if `a · b = b · a`. For ordinary
real numbers this is always true: `3 · 5 = 5 · 3`. So ℝ, ℚ, ℤ, ℂ,
𝔽₂, 𝔽ₚ are all **commutative rings/fields**.

### When does commutativity fail?

- **Matrices don't commute** in general:
  `[[1,2],[0,1]] · [[1,0],[1,1]] = [[3,2],[1,1]]` but
  `[[1,0],[1,1]] · [[1,2],[0,1]] = [[1,2],[1,3]]` — different.
- **3D rotations don't commute**: rotate a book 90° around the
  x-axis, then 90° around the y-axis — different orientation than
  doing y first, then x. Try it with a physical book.
- **Quaternions** (used in 3D graphics) don't commute.
- **Function composition** doesn't commute: applying `f` then `g`
  isn't generally the same as `g` then `f`.

### Why this matters for matmul

When you recurse a matmul algorithm to bigger sizes, the
"multiplications" between entries become multiplications between
sub-matrix blocks — and **matrix blocks don't commute in general**.

An algorithm that exploits `a · b = b · a` would compute wrong answers
when applied to non-commutative blocks. So such algorithms are
restricted to scalar matmul (matrices of plain numbers) and can't be
recursed.

This repo tags every scheme with both a field AND a commutative flag.
The combination is called an **Algebra** in the code — see
[`eu.solven.matmul.algebra.Algebra`](../../src/main/java/eu/solven/matmul/algebra/Algebra.java).

---

## What's "tensor rank" and why is it the same as multiplication count?

Skip this section unless you're curious about the mathematical
underpinning.

Every bilinear operation (like matmul) corresponds to a 3D tensor:
roughly, a 3-dimensional array of numbers `T[i, j, k]`. For matmul,
the entries are 1 when the indices match a valid product term, else 0.

A **rank-r decomposition** of a tensor is writing it as a sum of `r`
"simple" pieces (outer products of three vectors). The **tensor rank**
is the smallest such `r`.

The KEY fact: a rank-r decomposition of the matmul tensor IS exactly
an algorithm computing matmul with `r` multiplications. So:

**Finding a fast matmul algorithm ≡ finding a low-rank tensor
decomposition of the matmul tensor.**

This is why papers like AlphaTensor and AlphaEvolve use RL / search
to find tensor decompositions — that's literally the same thing as
finding fast matmul algorithms.

---

## What does this repo actually do?

1. **Catalog**: ~1900 known matrix-multiplication algorithms,
   organized by `(shape, field, commutativity)`, with their factor
   matrices on disk and a verifier that confirms each one actually
   computes matmul. See [`docs/`](../../docs) for the web-browsable
   version at [solven.eu/matmulcatalog](https://www.solven.eu/matmulcatalog/).

2. **Recursive composition**: given a small Strassen-like base and a
   target shape, compose recursively (with min-reduction tricks) to
   produce a verified bigger algorithm. See
   [`Recombination.java`](../../src/main/java/eu/solven/matmul/catalog/Recombination.java).

3. **Search**: reproduce and improve the Drevet–Islam–Schost 2009
   survey using modern catalog entries + multi-base + symmetry +
   Pan's trilinear aggregation. See [SEARCH_STRATEGY.md](SEARCH_STRATEGY.md)
   and [NEW_BOUNDS.md](../../research/NEW_BOUNDS.md) for what we've found.

4. **SAT-solver**: a separate pipeline that searches for new
   algorithms by encoding the tensor decomposition as a SAT problem.

5. **References + provenance**: every claim is tagged with its
   source paper. See [REFERENCES.md](../../REFERENCES.md).

---

## Where to go next

- **[README.md](../../README.md)** — top-level documentation map
- **[SMALL_MATMUL_CATALOG.md](SMALL_MATMUL_CATALOG.md)** — narrative
  walkthrough of the catalog
- **[RANK_KNOWLEDGE.md](RANK_KNOWLEDGE.md)** — theoretical
  rank bounds + field rules in depth
- **[OMEGA_HISTORY.md](OMEGA_HISTORY.md)** — the asymptotic ω race
  (Coppersmith–Winograd → Williams 2024 at 2.371552)
- **[NEW_BOUNDS.md](../../research/NEW_BOUNDS.md)** — bounds this repo establishes
- **[SEARCH_STRATEGY.md](SEARCH_STRATEGY.md)** — how to extend the
  search yourself

---

## Glossary cheat-sheet

- `⟨n, m, p⟩` — multiplying an `n×m` matrix by an `m×p` matrix
- **rank** of a matmul — minimum number of scalar multiplications
  needed; lower is better
- **scheme** — a concrete algorithm; in this repo, the U/V/W factor
  matrices encoding the bilinear products
- **field / algebra** — the arithmetic setting (ℝ, ℂ, 𝔽₂, ...) plus
  commutative-or-not
- **bilinear** — algorithm whose products are
  `(linear combo of A entries) · (linear combo of B entries)`
- **non-bilinear** — products can mix A and B entries on both sides;
  requires commutativity (see [Rosowski 2019](../../REFERENCES.md#28-rosowski19))
- **ω** — the matrix-multiplication exponent: how fast does matmul
  scale with `n`? Naive: `O(n³)`. Strassen: `O(n^2.807)`. Best known
  asymptotically: `O(n^2.371552)` (Williams 2024). True value
  conjecturally `2`.
