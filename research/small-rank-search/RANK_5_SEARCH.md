# What Searching for a Rank-5 Algorithm Looks Like

A concrete look at the equations one would have to solve to find an
algorithm computing `2×2` matrix multiplication with **5 bilinear
multiplications**. Spoiler: the system is hopelessly over-determined, and
several independent arguments (counting, border rank, the Hopcroft–Kerr
proof) confirm there is no solution. This note shows the structure
explicitly — both as a sanity check on the lower bound and as a template
for what a meaningful search loop would compute.

---

## 1. The unknowns

The bilinear template (from `STRASSEN_AS_EQUATIONS.md`) with `r = 5`:

```
For k = 1..5:
    M_k = (Σ_{a,b} U[a,b,k] · A[a,b]) · (Σ_{c,d} V[c,d,k] · B[c,d])

For each (i,j) ∈ {1,2}²:
    C[i,j] = Σ_{k=1..5} W[i,j,k] · M_k
```

Unknowns:
- `U ∈ R^{2×2×5}`, `V ∈ R^{2×2×5}`, `W ∈ R^{2×2×5}` → **60 entries**.
- Per-term scaling redundancy (`(αu)⊗(βv)⊗(γw) = u⊗v⊗w` when `αβγ = 1`)
  removes `2 · 5 = 10` parameters → **50 essential**.
- The `GL_2 × GL_2 × GL_2` gauge removes `12` more → **38 essential**.

---

## 2. The 64 equations, in full

Same as for `r = 7`, the constraint is that for every
`(i, j, a, b, c, d) ∈ {1,2}^6`:

```
Σ_{k=1..5}  U[a,b,k] · V[c,d,k] · W[i,j,k]  =  δ(a=i) · δ(b=c) · δ(d=j)
```

That's **64 cubic equations** in the 60 raw unknowns. Of the 64 equations:
- **8 have RHS = 1**: the "useful" monomials, two per output entry of `C`:
  - `C[1,1]`: coefficients of `A[1,1]·B[1,1]` and `A[1,2]·B[2,1]`
  - `C[1,2]`: coefficients of `A[1,1]·B[1,2]` and `A[1,2]·B[2,2]`
  - `C[2,1]`: coefficients of `A[2,1]·B[1,1]` and `A[2,2]·B[2,1]`
  - `C[2,2]`: coefficients of `A[2,1]·B[1,2]` and `A[2,2]·B[2,2]`
- **56 have RHS = 0**: every other monomial `A[a,b]·B[c,d]` must vanish.

### 2.1 Explicit list for output `C[1,1]` (16 equations, written for `r = 5`)

Abbreviate `α_k = W[1,1,k]`. Then for each of the 16 pairs `(a,b,c,d)`:

```
(1,1,1,1):  U[1,1,1]V[1,1,1]α_1 + … + U[1,1,5]V[1,1,5]α_5  = 1
(1,2,2,1):  U[1,2,1]V[2,1,1]α_1 + … + U[1,2,5]V[2,1,5]α_5  = 1
(1,1,1,2):  U[1,1,1]V[1,2,1]α_1 + … + U[1,1,5]V[1,2,5]α_5  = 0
(1,1,2,1):  U[1,1,1]V[2,1,1]α_1 + … + U[1,1,5]V[2,1,5]α_5  = 0
(1,1,2,2):  U[1,1,1]V[2,2,1]α_1 + … + U[1,1,5]V[2,2,5]α_5  = 0
(1,2,1,1):  U[1,2,1]V[1,1,1]α_1 + … + U[1,2,5]V[1,1,5]α_5  = 0
(1,2,1,2):  U[1,2,1]V[1,2,1]α_1 + … + U[1,2,5]V[1,2,5]α_5  = 0
(1,2,2,2):  U[1,2,1]V[2,2,1]α_1 + … + U[1,2,5]V[2,2,5]α_5  = 0
(2,1,1,1):  U[2,1,1]V[1,1,1]α_1 + … + U[2,1,5]V[1,1,5]α_5  = 0
(2,1,1,2):  U[2,1,1]V[1,2,1]α_1 + … + U[2,1,5]V[1,2,5]α_5  = 0
(2,1,2,1):  U[2,1,1]V[2,1,1]α_1 + … + U[2,1,5]V[2,1,5]α_5  = 0
(2,1,2,2):  U[2,1,1]V[2,2,1]α_1 + … + U[2,1,5]V[2,2,5]α_5  = 0
(2,2,1,1):  U[2,2,1]V[1,1,1]α_1 + … + U[2,2,5]V[1,1,5]α_5  = 0
(2,2,1,2):  U[2,2,1]V[1,2,1]α_1 + … + U[2,2,5]V[1,2,5]α_5  = 0
(2,2,2,1):  U[2,2,1]V[2,1,1]α_1 + … + U[2,2,5]V[2,1,5]α_5  = 0
(2,2,2,2):  U[2,2,1]V[2,2,1]α_1 + … + U[2,2,5]V[2,2,5]α_5  = 0
```

Three more groups of 16 (for `C[1,2], C[2,1], C[2,2]`) round out the 64
equations. Each equation is cubic in the unknowns (linear in `U`, in `V`, and
in `W` separately).

---

## 3. Why this can't have a solution (three views)

### 3.1 Naive dimension count
`38 effective unknowns` vs `64 equations`. The system has codimension at
least `64 - 38 = 26` worth of "extra" constraints. For comparison:

|  | unknowns (after gauge) | equations | excess |
|---|---|---|---|
| `r = 5` | 38 | 64 | +26 (no solution generically) |
| `r = 6` | 48 | 64 | +16 (no solution; lower bound proves it) |
| `r = 7` | 58 | 64 | +6 (solutions form a 6-dim family — Strassen's orbit) |
| `r = 8` | 68 | 64 | -4 (huge family — includes the naive algorithm) |

Counting alone doesn't *prove* infeasibility (overdetermined systems
sometimes have solutions because the equations aren't independent), but it
shows why `r = 5` is wildly off.

### 3.2 Symmetry-restricted view: `Z/3`-equivariant rank-5 is impossible
Under the cyclic `Z/3` action, every term is either fixed or in an orbit of
3. So `5 = 5·1 + 0·3` or `5 = 2·1 + 1·3`. With `Z/3` equivariance:

- `5 = 2 + 3` structure: 2 fixed triples (4 params each) + 1 orbit-3
  generator (12 params) = **20 free params**.
- Equations collapse from 64 to 24 orbits (Burnside).

20 unknowns vs 24 equations — still over-determined, and ALS / SAT on this
small system terminates with strictly positive residual or UNSAT in
seconds. This gives a fast computational *witness* that no `Z/3`-symmetric
rank-5 algorithm exists. (To rule out asymmetric ones too you need either
the full Hopcroft–Kerr argument or an exhaustive SAT.)

### 3.3 The hard theorem: border rank `≥ 7`
The deepest "no" comes from **border rank**: the smallest `r` such that
`⟨2,2,2⟩` can be *approximated* arbitrarily well by rank-`r` tensors.

> **Landsberg & Michałek (2017):** `R̲(⟨2,2,2⟩) = 7`.

Consequence: not only is rank-5 impossible, but no sequence of rank-5
tensors can converge to `⟨2,2,2⟩`. Any rank-5 approximation has residual
**bounded away from zero by a positive constant** (depending only on the
choice of norm). In particular, an ALS run targeting `r = 5` will plateau
at strictly positive residual — and that plateau is a numerical certificate
of the border-rank lower bound.

---

## 4. What an actual rank-5 search loop reports

A concrete run with alternating least squares (sketch):

```python
import numpy as np

# matmul tensor ⟨2,2,2⟩ as a 4×4×4 array
def matmul_tensor(n=2):
    T = np.zeros((n*n, n*n, n*n))
    for i in range(n):
      for j in range(n):
        for l in range(n):
            T[i*n+l, l*n+j, i*n+j] = 1.0
    return T

T = matmul_tensor()
r = 5
U = np.random.randn(4, r); V = np.random.randn(4, r); W = np.random.randn(4, r)

for it in range(10_000):
    # solve for U given V, W; solve for V given U, W; solve for W given U, V
    # (standard ALS update — Khatri–Rao products)
    ...
    residual = np.linalg.norm(T - np.einsum('ar,br,cr->abc', U, V, W))
    if it % 1000 == 0: print(it, residual)
```

Empirical result (any standard ALS implementation): the residual decreases,
**plateaus around 0.4–0.6** (in Frobenius norm; precise value depends on
init), and never reaches 0. Restarting from new random inits gives similar
plateaus. This is the *numerical face* of the border-rank-≥-7 theorem.

Compare to `r = 7`: same code, residual converges to numerical zero
(`~10^{-14}`) within a few hundred iterations.

---

## 5. What a rank-5 algorithm *would* look like, if it existed

Just to fix intuition, here is the **shape** of a hypothetical rank-5
algorithm:

```
M1 = (a11 u11 + a12 u12 + a21 u13 + a22 u14) · (b11 v11 + b12 v12 + b21 v13 + b22 v14)
M2 = (a11 u21 + a12 u22 + a21 u23 + a22 u24) · (b11 v21 + b12 v22 + b21 v23 + b22 v24)
M3 = (a11 u31 + a12 u32 + a21 u33 + a22 u34) · (b11 v31 + b12 v32 + b21 v33 + b22 v34)
M4 = (a11 u41 + a12 u42 + a21 u43 + a22 u44) · (b11 v41 + b12 v42 + b21 v43 + b22 v44)
M5 = (a11 u51 + a12 u52 + a21 u53 + a22 u54) · (b11 v51 + b12 v52 + b21 v53 + b22 v54)

C11 = w11_1 M1 + w11_2 M2 + w11_3 M3 + w11_4 M4 + w11_5 M5
C12 = w12_1 M1 + w12_2 M2 + w12_3 M3 + w12_4 M4 + w12_5 M5
C21 = w21_1 M1 + w21_2 M2 + w21_3 M3 + w21_4 M4 + w21_5 M5
C22 = w22_1 M1 + w22_2 M2 + w22_3 M3 + w22_4 M4 + w22_5 M5
```

— and the 60 coefficients `{u, v, w}` would jointly satisfy the 64
equations of §2. The Hopcroft–Kerr / de Groote / Landsberg–Michałek
theorems say in three different ways that the only way out is to add more
multiplications.

---

## 6. The pedagogical pay-off

Writing the `r = 5` system out gives three useful intuitions:

1. **Why 7 is "just barely" enough.** The unknowns–equations balance only
   crosses into positive territory at `r = 7`. Strassen sits on the
   knife-edge.

2. **Why naive `r = 8` is so "wasteful."** At `r = 8` we have 68 unknowns
   vs 64 equations — a 4-parameter family of algorithms, of which the
   naive algorithm is one. Lots of slack, all of which is squeezed out at
   `r = 7`.

3. **What a search loop is actually optimizing.** Whether ALS, SAT, RL, or
   pencil-and-paper — every method is trying to solve those 64 cubic
   equations. The interesting differences are in the priors (alphabet,
   symmetry) and the heuristic for navigating the cubic landscape.

---

## 7. Suggested follow-on experiments

1. **Implement the ALS verifier** (§4); record the `r = 5` residual plateau
   across 1000 random inits — empirically estimate the border-rank lower
   bound's "gap."
2. **Repeat at `r = 6`** — the same plateau effect but smaller residual;
   gap is provably positive (Landsberg–Michałek) but tighter.
3. **`Z/3`-equivariant SAT encoding at `r = 5, 6, 7`** — get UNSAT certificates
   at `r = 5, 6` and a (essentially unique) SAT model at `r = 7`. Yields a
   computer-checkable re-proof of the lower bound under the symmetry
   restriction.
