# Strassen as a System of Equations, and How to Rediscover It

This note shows Strassen's 7-multiplication algorithm as a constrained system,
then walks through the *structural* ideas one can use to (re)discover it
without brute search.

---

## 1. The algorithm, for reference

Let `A, B, C ∈ R^{2×2}` with `C = A·B`. Strassen (1969):

```
M1 = (A11 + A22) · (B11 + B22)
M2 = (A21 + A22) · B11
M3 = A11       · (B12 - B22)
M4 = A22       · (B21 - B11)
M5 = (A11 + A12) · B22
M6 = (A21 - A11) · (B11 + B12)
M7 = (A12 - A22) · (B21 + B22)

C11 = M1 + M4 - M5 + M7
C12 = M3 + M5
C21 = M2 + M4
C22 = M1 - M2 + M3 + M6
```

Seven products of linear combinations, recombined by signed sums.

---

## 2. The bilinear-algorithm template

Any bilinear algorithm with `r` multiplications has the shape:

```
For k = 1..r:
    M_k = (Σ_{a,b} U[a,b,k] · A[a,b]) · (Σ_{c,d} V[c,d,k] · B[c,d])

For each output (i,j):
    C[i,j] = Σ_k W[i,j,k] · M_k
```

The unknowns are three coefficient tensors:

- `U ∈ R^{2×2×r}`  (how to combine entries of `A`)
- `V ∈ R^{2×2×r}`  (how to combine entries of `B`)
- `W ∈ R^{2×2×r}`  (how to recombine the products)

For Strassen, `r = 7`. So the search space has `3 × 4 × 7 = 84` real
(or integer / `{-1,0,1}`) unknowns.

---

## 3. The constraints as a system

Substitute the template into `C = A·B`. The product
`(Σ U · A)·(Σ V · B)` expands as `Σ U V · A·B` (treating `A[a,b]` and `B[c,d]`
as independent indeterminates):

```
C[i,j] = Σ_k W[i,j,k] · Σ_{a,b,c,d} U[a,b,k] · V[c,d,k] · A[a,b] · B[c,d]
       = Σ_{a,b,c,d} ( Σ_k U[a,b,k] · V[c,d,k] · W[i,j,k] ) · A[a,b] · B[c,d]
```

Matching this against the desired `C[i,j] = Σ_l A[i,l] · B[l,j]`, the
coefficient of `A[a,b]·B[c,d]` in `C[i,j]` must be `1` if `(a=i, b=c, d=j)`
and `0` otherwise. That gives **`2^6 = 64` cubic equations** in the 84
unknowns:

> For all `i, j, a, b, c, d ∈ {1,2}`:
>
>     `Σ_k=1..7  U[a,b,k] · V[c,d,k] · W[i,j,k]  =  δ(a=i) · δ(b=c) · δ(d=j)`

Equivalently, defining the matrix-multiplication tensor

```
T[(a,b), (c,d), (i,j)] = δ(a=i) · δ(b=c) · δ(d=j)            ∈ R^{4×4×4}
```

the constraint is a **tensor rank-7 decomposition**:

```
T = Σ_{k=1..7}  u_k ⊗ v_k ⊗ w_k
```

where `u_k = U[:,:,k] ∈ R^4`, similarly `v_k, w_k`. **Strassen's solution is
one particular rank-7 decomposition of `T = ⟨2,2,2⟩`.** Hopcroft–Kerr and
Winograd 1971 proved no rank-6 decomposition exists; de Groote 1978 proved
all rank-7 decompositions are equivalent under a known group action.

### 3.1 The 64 equations, written out

Group them by `(i,j)` — for each output entry there are 16 equations (one per
`(a,b,c,d)`):

`C[1,1] = A[1,1]B[1,1] + A[1,2]B[2,1]`:
- Coefficient of `A[1,1]B[1,1]`: `Σ_k U[1,1,k]·V[1,1,k]·W[1,1,k] = 1`
- Coefficient of `A[1,2]B[2,1]`: `Σ_k U[1,2,k]·V[2,1,k]·W[1,1,k] = 1`
- Coefficient of the other 14 `A[a,b]B[c,d]` pairs: `= 0`

…and similarly for `C[1,2], C[2,1], C[2,2]`. Total 64 equations.

### 3.2 As a matrix-flattening equation

A common rewriting that makes the structure visible: define
the `16 × 7` matrix `UV` with rows indexed by `(a,b,c,d)` and column `k` equal
to `U[a,b,k]·V[c,d,k]`. Define `W` as a `7 × 4` matrix indexed by
`k, (i,j)`. The constraint is

```
UV · W  =  T_flat
```

where `T_flat ∈ R^{16×4}` is the flattened matmul tensor, with rows the 16
`(a,b,c,d)` and columns the 4 outputs `(i,j)`. The left side is *quadratic*
in the unknowns `U, V` and linear in `W`. This is the form solved by
*alternating least squares* (ALS, see §5).

---

## 4. Why 7? A back-of-the-envelope sanity check

`T = ⟨2,2,2⟩` is a `4×4×4` tensor (64 entries, of which 8 are 1's and 56 are
0's). A rank-`r` decomposition has `r·(4+4+4) = 12r` parameters, but the
decomposition is invariant under independent scaling
`(u_k, v_k, w_k) → (α u_k, β v_k, γ w_k)` with `αβγ = 1` (`2` params per `k`),
and the rank-1 terms can be permuted (`r!` symmetry). So the *essential*
parameter count is `r·(12 - 2) = 10r`, minus discrete symmetry.

Dimension count vs equations:
- `r = 6`: 60 essential parameters vs 64 equations → generically *no
  solution* (and indeed there is none).
- `r = 7`: 70 essential parameters vs 64 equations → 6-dimensional family
  of solutions. De Groote: this family is a single `GL_2 × GL_2 × GL_2`
  orbit modulo the cyclic `Z/3` action.

The naive heuristic "70 > 64, so it might work" is consistent with the
existence proof but does *not* explain why 7 is tight (the lower bound is
much subtler).

---

## 5. Ways to rediscover Strassen without brute force

### 5.1 Counting and the "shared sub-expression" insight
The naive algorithm computes 8 products
`A[i,l]·B[l,j]` for `(i,l,j) ∈ {1,2}^3`. To save one product, we need a
single bilinear quantity `M = (sum of A entries)·(sum of B entries)` whose
expansion contains *two* of those 8 monomials with matching coefficients, so
that `M` substitutes for both.

The simplest such product: `M1 = (A11 + A22)(B11 + B22)` expands to
```
A11 B11  +  A11 B22  +  A22 B11  +  A22 B22
```
of which `A11 B11` is part of `C11` and `A22 B22` is part of `C22` — two
"useful" monomials at the cost of two "junk" monomials (`A11 B22, A22 B11`).
The junk has to be canceled by other products. This is the *single key idea*
behind Strassen. From here, a methodical search over which pair of diagonal
or off-diagonal terms to merge — and which compensating products cancel the
junk — converges on Strassen with a few hours of pencil work.

This was, by Strassen's own account (and Pan's later commentary), the route
to the original discovery.

### 5.2 Exploit symmetry: the cyclic `Z/3` action
The matmul tensor is *cyclically symmetric*:
`trace(A·B·C) = trace(B·C·A) = trace(C·A·B)`. Algebraically, `⟨n,n,n⟩` is
invariant under cyclic permutation of its three factors. So one may **search
for `Z/3`-symmetric decompositions**: pick orbits of size 1 or 3 under the
cyclic action, and look for products that are themselves invariant or come
in orbits of three.

Strassen's algorithm decomposes as **1 + 3 + 3 = 7** under this `Z/3` action:
- One product is "diagonal" / cyclically fixed: `M1 = (A11+A22)(B11+B22)`.
- Three products form one orbit (the "edge" products `M5, M6, M7` after
  relabeling).
- Three products form another orbit (`M2, M3, M4`).

Restricting the unknowns to be `Z/3`-equivariant cuts the parameter count
by a factor of `~3` and turns the search into something a SAT/SMT solver or
hand calculation can settle. **Heun 1994** and **Oh, Hopcroft 1979** used
this style to enumerate algorithms; **Smirnov** built families of `3×3`
algorithms via the same approach.

### 5.3 Alternating Least Squares (ALS)
Treat the constraint `UV·W = T_flat` (§3.2) as a fixed-point iteration:
1. Initialize `U, V, W` randomly.
2. Fix `U, V` → solve for `W` linearly (it's an over-determined least-squares
   problem).
3. Fix `U, W` → solve for `V` linearly.
4. Fix `V, W` → solve for `U` linearly.
5. Repeat until the residual vanishes.

For `r = 7` and `⟨2,2,2⟩`, ALS converges to a rank-7 decomposition from a
random start with non-trivial probability. **Smirnov** uses this routinely
for larger formats. The catch: ALS gets stuck in local minima of rank
`> 7`; restart many times. AlphaTensor essentially does a much smarter
version of this with RL.

### 5.4 Algebraic identity / polynomial factoring
Strassen's algorithm can be re-derived as a clever rewriting of
`trace((A − λI)·(B − μI))` or of
`det(A + tB)` expansions. Pan's 1980 book and several modern surveys present
Strassen via the identity
```
(a+d)(e+h) + d(f−e) + a(g−h) = ae + dh + … (rearrange)
```
This is more presentational than discovery-oriented but useful for teaching.

### 5.5 Restrict the alphabet
Strassen's coefficients all lie in `{-1, 0, +1}`. Empirically, this is true
for most known small-format algorithms. Restricting `U, V, W` entries to
`{-1, 0, 1}` (or `{-1, 0, 1, ±1/2}`) shrinks an infinite parameter space
to a finite one — `3^84 ≈ 10^40` candidates for `r=7`, still huge, but
amenable to branch-and-bound with symmetry pruning. This is essentially the
search space AlphaTensor explored, with neural-network-guided priors.

### 5.6 Border-rank / approximation as a guide
Even when an exact rank-7 decomposition is hard to find, *border-rank*
decompositions (algorithms in `ε` valid as `ε → 0`) are often easier; the
limit can sometimes be massaged into an exact algorithm. **Bini's 1979
algorithms** for `⟨2,2,3⟩` were found this way. For `⟨2,2,2⟩` itself, the
border rank already equals 7, so this trick gives nothing extra — but for
the *rediscovery exercise*, it is a useful warm-up.

### 5.7 Reinforcement learning (AlphaTensor 2022)
Frame the problem as a single-player game: start with the tensor `T`,
each move subtracts a rank-1 tensor `u ⊗ v ⊗ w`, lose when stuck above zero
in more moves than your budget. Train a value/policy network (AlphaZero
style) over actions = entries in a finite alphabet. This rediscovers
Strassen and finds new algorithms in larger formats, *without* hand-derived
symmetry — though using symmetry as a data-augmentation prior helps.
Reference: Fawzi et al., *Nature* 610 (2022).

---

## 6. Minimal verification (sketch)

A direct way to verify any candidate `(U, V, W)` is to expand `Σ_k (U_k·A)
(V_k·B)·W_k` symbolically and check it equals `A·B`. With 16 monomials
`A[a,b]B[c,d]` per output and 4 outputs, this is a 64-coefficient check that
a CAS (SymPy, Mathematica) does in milliseconds.

Pseudocode:

```python
import sympy as sp

A = sp.MatrixSymbol('A', 2, 2)
B = sp.MatrixSymbol('B', 2, 2)
# define M1..M7 as sympy expressions per §1
C_strassen = sp.Matrix(2,2, lambda i,j: ...)   # recombination
C_truth   = sp.Matrix(A) * sp.Matrix(B)
assert sp.simplify(C_strassen - C_truth) == sp.zeros(2,2)
```

This same script, given any candidate `(U, V, W)` tensor, validates it. A
search loop wrapping this with ALS or branch-and-bound is the practical
starting point for repo-level experiments.

---

## 7. Suggested next experiments in this repo

1. **Implement the verifier** in §6 — it is the foundation of every search.
2. **Implement ALS** for `⟨2,2,r⟩` with `r ∈ {5,6,7}`; observe that `r=6`
   never converges to residual 0 (numerical hint of the lower bound).
3. **Symbolic search with `{-1,0,1}` alphabet and `Z/3` symmetry** —
   small enough to enumerate, large enough to non-trivially rediscover
   Strassen.
4. **Try `⟨3,3,3⟩` with `r = 21`** and the analogous symmetry restriction —
   you'd be reproducing the frontier of an open problem (current upper
   bound is 23).
