# Dense Matrix Multiplication: Algorithms Review

A focused survey on fast algorithms for dense matrix multiplication, centered on
Strassen-style divide-and-conquer schemes, their optimality, and the directions
in which the optimality result can be relaxed.

---

## 1. Landscape of algorithms

### 1.1 Naive algorithm
- Computes `C = A·B` with `N^3` scalar multiplications and `N^2(N-1)` additions.
- Exponent: `ω = 3`.
- Numerically the most stable; cache-friendly when properly blocked.

### 1.2 Strassen (1969)
Strassen, *"Gaussian elimination is not optimal"*, Numerische Mathematik 13.
- For `2×2` blocks, computes the product with **7 block multiplications** and 18
  block additions (instead of `8` and `4`).
- Recursing on `N = 2^k`: `T(N) = 7·T(N/2) + Θ(N^2)`, giving
  `ω ≤ log_2 7 ≈ 2.8074`.
- Works over **any ring** (entries need not commute), which is why it lifts to
  block matrices unchanged.

### 1.3 Winograd's variant (1971)
Same 7 multiplications, but only **15 additions/subtractions** (Strassen's
original uses 18). This is the variant actually implemented in most libraries
(GotoBLAS variants, BLIS, Intel MKL's Strassen path).

### 1.4 Subsequent asymptotic improvements
A long line of work pushed `ω` further down. None of these are Strassen-like
divide-and-conquer in the recursive `2×2` sense; they instead exploit
*aggregating, canceling, and uniting* (Pan), border rank (Bini), or the laser
method on large structured tensors (Coppersmith–Winograd and successors).

| Year | Author(s) | Bound on `ω` |
|---|---|---|
| 1969 | Strassen | 2.8074 |
| 1978 | Pan | 2.795 (then 2.78) |
| 1979 | Bini, Capovani, Romani, Lotti (border rank) | 2.7799 |
| 1981 | Schönhage (τ-theorem) | 2.522 |
| 1986 | Strassen (laser method) | 2.479 |
| 1990 | Coppersmith–Winograd | 2.376 |
| 2010 | Stothers | 2.3737 |
| 2012 | Vassilevska Williams | 2.3729 |
| 2014 | Le Gall | 2.37287 |
| 2021 | Alman–Williams | 2.37286 |
| 2022–24 | Duan–Wu–Zhou, Williams–Xu–Xu–Zhou | 2.3715… |

All of these post-Strassen results are **galactic**: the crossover with the
naive algorithm is at sizes vastly larger than anything physically realizable.
The only fast algorithm that beats naive on real hardware (around `N ≈ 500`–
`1000` for double precision) is Strassen/Winograd, possibly with 1 or 2 levels
of recursion.

---

## 2. Why 7 multiplications, and why not 6?

The relevant complexity measure is the **bilinear rank** (or just *rank*) of
the matrix-multiplication tensor `⟨n,m,p⟩`. Any algorithm of the form

> compute `r` products `m_k = (Σ α_{ijk} a_{ij}) · (Σ β_{ijk} b_{ij})`,
> then `c_{ij} = Σ γ_{ijk} m_k`

has *bilinear complexity* `r`. The minimum such `r` is the rank
`R(⟨n,m,p⟩)`. For `2×2` multiplication the tensor is `⟨2,2,2⟩` and Strassen
proves `R(⟨2,2,2⟩) ≤ 7`.

### 2.1 The matching lower bound: `R(⟨2,2,2⟩) = 7`
Two independent proofs in 1971:

- **Hopcroft & Kerr 1971**, *"On minimizing the number of multiplications
  necessary for matrix multiplication"*, SIAM J. Applied Math 20(1).
- **Winograd 1971**, *"On multiplication of 2×2 matrices"*, Linear Algebra and
  its Applications 4.

Both prove that 7 non-scalar multiplications are necessary, over any (infinite)
field. Winograd also showed that **15 additions** are necessary if one uses
exactly 7 multiplications — so Strassen–Winograd is additively optimal too.

A later, simpler proof: **de Groote 1978**, *"On varieties of optimal
algorithms for the computation of bilinear mappings"*, Theoretical CS 7.
De Groote classifies *all* optimal algorithms for `⟨2,2,2⟩`: they form a
single orbit under a known group action — every 7-multiplication algorithm is
obtained from Strassen's by a change of basis on the input/output spaces. There
is, in a precise sense, *one* algorithm.

### 2.2 What exactly does the proof assume?
Read carefully, the Hopcroft–Kerr / Winograd / de Groote bound says:

1. The model is **bilinear** (or, equivalently, *non-commutative* in the sense
   that entries of `A` and `B` are treated as independent indeterminates).
2. The algorithm is **exact** over a field — no limits, no approximations.
3. Only **multiplications of input-derived quantities** are counted; additions,
   subtractions, and scalar multiplications by field constants are free.
4. The product is `2×2` by `2×2`, square.

Each of these assumptions is a door to a relaxation.

---

## 3. Relaxing the constraints

The user's question — *what if we divide into 9 blocks instead of 4, or use
larger sums of smaller blocks?* — corresponds exactly to relaxing (4), and the
border-rank story corresponds to relaxing (2). Both have been deeply studied.

### 3.1 Relaxing the block shape: larger base cases
Replace the `2×2` base case with `n×n` (or rectangular). Recursing gives
`ω ≤ log_n R(⟨n,n,n⟩)`. To **beat Strassen** at base `n` you need

> `R(⟨n,n,n⟩) < n^{log_2 7}`.

| `n` | naive `n^3` | needed to beat Strassen | best known upper | best known lower |
|---|---|---|---|---|
| 2 | 8 | 7 | **7** (Strassen) | **7** (Hopcroft–Kerr) |
| 3 | 27 | < 21.85, i.e. ≤ 21 | **23** (Laderman 1976; many other 23-algs by Smirnov, Oh et al.) | **19** (Bläser 2003) |
| 4 | 64 | < 49 | **47** via 2×Strassen; **49** direct (Pan); AlphaTensor: 47 over `Z/2Z` | 34 (Landsberg) |
| 5 | 125 | < 110 | **100** (Smirnov 2017, approximate / specific) — exact ≤ ≈ 100 (Sedoglavic–Smirnov) | 47 (Massarenti–Raviolo) |

Key takeaway: **no exact algorithm at any `n ≥ 3` is known to beat Strassen
when applied recursively**. The 3×3 case has a 19/23 gap — open since the
1970s. Closing it to 21 or below would supersede Strassen.

References:
- **Laderman 1976**, *"A noncommutative algorithm for multiplying (3×3)
  matrices using 23 multiplications"*, Bull. AMS 82(1).
- **Smirnov 2013**, *"Bilinear complexity and practical algorithms for matrix
  multiplication"*, Comput. Math. Math. Phys. 53. Gives many `<n,m,p>`
  algorithms, including new 3×3 schemes.
- **Heun 1994** and **Bläser 2003**, *"On the complexity of the multiplication
  of matrices of small formats"*, J. Complexity 19 — lower bound
  `R(⟨3,3,3⟩) ≥ 19`.
- **Landsberg 2014**, *"New lower bounds for the rank of matrix
  multiplication"*, SIAM J. Comput. 43.
- **Fawzi et al. (AlphaTensor) 2022**, Nature 610. RL-discovered algorithms;
  improves several rectangular cases and modular-arithmetic cases, but does
  **not** beat Strassen for square exact `<n,n,n>` over `R` or `C`.

### 3.2 Rectangular base cases `⟨n,m,p⟩`
A more flexible relaxation: most progress on `ω` since Pan has actually come
from finding low-rank algorithms for *rectangular* multiplications, then
combining them via the Schönhage τ-theorem, which states

> `R(⟨n,m,p⟩) ≤ r`  implies  `ω ≤ 3 · log_{nmp} r`.

So a great `⟨2,3,4⟩` or `⟨3,3,6⟩` algorithm can yield ω improvements without
ever directly improving `⟨2,2,2⟩`. This is the route Pan, Bini, and Strassen
took in the late 70s / early 80s. The user's intuition — *"sums of many
smaller blocks"* — is exactly this idea: don't insist on a square base case.

### 3.3 Relaxing exactness: border rank
**Bini, Capovani, Romani, Lotti 1979** introduced *border rank* `R̲(T)`: the
smallest `r` such that `T` is the *limit* of rank-`r` tensors. Equivalently,
algorithms parameterized by an indeterminate `ε`, valid as `ε → 0`.

Crucially:

> `R̲(⟨2,2,2⟩) = 7` (Landsberg 2006) — border rank is also 7 here.
> But `R̲(⟨2,2,3⟩) = 14`, and the analogous *exact* rank is 15.

So for the very smallest square case, border rank gives nothing extra; the
gains appear in rectangular cases and at larger formats. This is the
mathematical reason "approximate algorithms for `2×2`" do not exist below 7
multiplications. The often-cited "Bini 6 multiplications for `2×2`" actually
refers to a partial `2×2` product (one output entry suppressed), not the full
`<2,2,2>` tensor.

### 3.4 Relaxing the ring: commutativity, characteristic, finite fields
- **Commutative entries** (e.g. `A`, `B` over a commutative ring viewed as
  scalars, not block matrices): allows divisions of input variables.
  **Pan 1980** and **Makarov 1970** gave commutative algorithms for `3×3` with
  fewer multiplications, but these *cannot be applied recursively* because
  block matrices do not commute. So commutativity does not help asymptotically.
- **Finite-field tricks (`Z/2Z`)**: AlphaTensor found 47-multiplication
  algorithms for `<4,4,4>` over `Z/2Z`. Over `Z/2Z`, addition = subtraction =
  XOR, and the proof of `R(⟨2,2,2⟩) ≥ 7` still holds, but rectangular cases can
  diverge. Use case: cryptographic linear algebra, not numerical computing.
- **Characteristic 0** vs positive characteristic: the lower bound proof for
  `<2,2,2>` does not care, so 7 is tight over `R`, `C`, `Q`, `Z/pZ` alike.

### 3.5 Relaxing what we count: divisions, scaling, preprocessing
- Allowing **divisions** does not help asymptotically (Strassen 1973 — every
  division-using arithmetic circuit can be simulated without divisions with
  only polynomial blowup).
- Allowing **arbitrary preprocessing on `A`** (e.g. precomputing once,
  multiplying by many `B`'s): the relevant measure becomes
  `R(⟨n,1,n⟩)`-style and is well understood; modest gains.
- Allowing **probabilistic / randomized** algorithms: no improvement over the
  deterministic bilinear lower bound.

### 3.6 Practical relaxations Strassen users actually care about
- **Numerical stability**: Strassen's componentwise error bound is worse than
  naive's by a factor depending on the recursion depth. *Higham 1990*,
  *"Exploiting fast matrix multiplication within the Level 3 BLAS"*. In
  practice this caps recursion at 1–2 levels.
- **Memory**: naive needs no scratch; Strassen needs `O(N^2)` extra. There are
  *in-place* and *low-memory* variants (Boyer–Pernet–Trumbore; Huang–Smith).
- **Crossover size**: roughly `N ≈ 500`–`2000` depending on architecture for
  one level of recursion.

---

## 4. Synthesis: can Strassen's 7 be improved?

**For the literal question — `2×2` exact bilinear multiplication, any ring:
no.** The Hopcroft–Kerr / Winograd / de Groote lower bound is matching, and
de Groote even shows there is essentially one optimal algorithm up to
symmetry. Border rank does not help here. This is one of the most-settled
lower bounds in algebraic complexity.

**For "Strassen-style recursive divide-and-conquer" more broadly — open.** The
asymptotic exponent `ω` is conjectured to be 2, but the gap between
`ω ≤ 2.3715…` and `ω ≥ 2` is wide. The natural next step is a better algorithm
at base 3 (anywhere in `[19, 22]` would beat Strassen) — open for 50 years.

**For practical computation — Strassen–Winograd is the only fast algorithm
that wins on real hardware**, and only with shallow recursion. Everything past
Coppersmith–Winograd is theoretical.

The most promising relaxations are:
1. **Larger / rectangular base cases** — actively explored, including ML-aided
   search (AlphaTensor, FBHHRBNRSSSHK, OpenAI's recent work on tensor
   decomposition).
2. **Border rank / approximate algorithms** — the workhorse of all
   post-Strassen asymptotic results, but inapplicable below 7 for `<2,2,2>`
   itself.
3. **Specific rings (modular, GF(2))** — practical for crypto/coding theory,
   irrelevant for numerical linear algebra.

---

## 5. Annotated references

Foundational:
- V. Strassen, *Gaussian elimination is not optimal*, Numer. Math. 13 (1969),
  354–356. **The original.**
- S. Winograd, *On multiplication of 2×2 matrices*, Lin. Alg. Appl. 4 (1971),
  381–388. **Lower bound + 15-addition variant.**
- J. E. Hopcroft, L. R. Kerr, *On minimizing the number of multiplications
  necessary for matrix multiplication*, SIAM J. Appl. Math. 20 (1971),
  30–36. **Independent lower bound.**
- H. F. de Groote, *On varieties of optimal algorithms for the computation of
  bilinear mappings II*, Theor. Comput. Sci. 7 (1978), 127–148.
  **Classifies all 7-mult algorithms — essentially unique.**

Larger base cases:
- J. D. Laderman, *A noncommutative algorithm for multiplying (3×3) matrices
  using 23 multiplications*, Bull. AMS 82 (1976), 126–128.
- A. V. Smirnov, *Bilinear complexity and practical algorithms for matrix
  multiplication*, Comput. Math. Math. Phys. 53 (2013), 1781–1795.
- M. Bläser, *On the complexity of the multiplication of matrices of small
  formats*, J. Complexity 19 (2003), 43–60. **Lower bound 19 for `<3,3,3>`.**
- J. M. Landsberg, *New lower bounds for the rank of matrix multiplication*,
  SIAM J. Comput. 43 (2014), 144–149.

Border rank, asymptotic:
- D. Bini, M. Capovani, F. Romani, G. Lotti, *O(n^{2.7799}) complexity for
  n×n approximate matrix multiplication*, Inf. Process. Lett. 8 (1979),
  234–235.
- A. Schönhage, *Partial and total matrix multiplication*, SIAM J. Comput. 10
  (1981), 434–455.
- D. Coppersmith, S. Winograd, *Matrix multiplication via arithmetic
  progressions*, J. Symb. Comput. 9 (1990), 251–280.
- J. Alman, V. Vassilevska Williams, *A refined laser method and faster matrix
  multiplication*, SODA 2021.
- R. Duan, H. Wu, R. Zhou, *Faster matrix multiplication via asymmetric
  hashing*, FOCS 2023.

ML-discovered:
- A. Fawzi et al., *Discovering faster matrix multiplication algorithms with
  reinforcement learning*, Nature 610 (2022). **AlphaTensor.**

Practical:
- N. J. Higham, *Exploiting fast matrix multiplication within the Level 3
  BLAS*, ACM TOMS 16 (1990), 352–368. **Stability analysis.**
- P. D'Alberto, M. Bodrato, A. Nicolau, *Exploiting parallelism in matrix
  computation kernels for symmetric multiprocessor systems*, ACM TOMS 38
  (2011). **Real implementations.**

Surveys / textbooks:
- P. Bürgisser, M. Clausen, M. A. Shokrollahi, *Algebraic Complexity Theory*,
  Springer, 1997. **The standard reference; chapters 14–15 cover everything
  above with proofs.**
- M. Bläser, *Fast matrix multiplication*, Theory of Computing Graduate
  Surveys 5 (2013). **Modern, free, well-written survey.**
- V. Vassilevska Williams, *On the complexity of matrix multiplication*, PhD
  thesis, Berkeley 2008 (and her later survey talks).

---

## 6. Open problems worth investigating in this repo

1. **`R(⟨3,3,3⟩) ≤ 22`?** Would beat Strassen. The current upper bound has
   stood at 23 since 1976. Active ML search target.
2. **Tight border rank `R̲(⟨n,n,n⟩)` for small `n`** — even `n=3` is open
   (current: `15 ≤ R̲(⟨3,3,3⟩) ≤ 21`).
3. **Better numerically-stable variants of Strassen** — alternative bases of
   the 7-mult algorithm with smaller error constants (Ballard et al. have
   explored this; many in de Groote's orbit are untested).
4. **Implementation engineering**: most BLAS libraries do not ship Strassen by
   default; experimenting with the crossover point and memory layout is a
   tractable, valuable engineering project.
