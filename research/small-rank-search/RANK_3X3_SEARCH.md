# The `⟨3,3,3⟩` Problem: Searching at Rank 21 and 22

A focused note on the **3×3 matrix multiplication tensor** `T = ⟨3,3,3⟩`.
This is the smallest open frontier of the bilinear-complexity world: the
best known upper bound is **23 multiplications** (Laderman 1976; Smirnov
2013/2017; matched by AlphaTensor 2022), and the best known lower bound
is `19`. Closing the gap matters because each unit of rank at this format
moves the matrix-multiplication exponent `ω` directly.

The two interesting thresholds:

| rank `r` | `log_3(r)` | meaning |
|---|---|---|
| `23` | `2.854` | current state of the art (Laderman et al.) |
| **`22`** | **`2.814`** | **would beat the world record** |
| **`21`** | **`2.771`** | **would also beat recursive Strassen** (`log_2 7 ≈ 2.807`) |
| `19` | `2.680` | the best known lower bound (still open) |

So `r = 22` is the next milestone for the field; `r = 21` is the
"Strassen-crossing" milestone — recursive Strassen on 6×6 = (2×2)⊗(3×3)
would become strictly worse than recursing directly on 3×3 once an
`r=21` algorithm is known.

---

## 1. The algorithm template

Same bilinear template as in [STRASSEN_AS_EQUATIONS.md](../../paper/theory/STRASSEN_AS_EQUATIONS.md), but with all
indices ranging over `{1,2,3}`:

```
For k = 1..r:
    M_k = (Σ_{a,b ∈ {1,2,3}} U[a,b,k] · A[a,b]) · (Σ_{c,d} V[c,d,k] · B[c,d])

For each (i,j) ∈ {1,2,3}²:
    C[i,j] = Σ_{k=1..r} W[i,j,k] · M_k
```

Unknowns:
- `U, V, W ∈ R^{3×3×r}` — `3 · 9 · r = 27r` entries.
- For `r = 21`: **567 raw entries**.
- For `r = 22`: **594 raw entries**.
- For `r = 23`: **621 raw entries**.

---

## 2. The equations

Substituting the template into `C = A·B` and matching coefficients of
each `A[a,b] · B[c,d]` monomial in each `C[i,j]`:

> For all `i, j, a, b, c, d ∈ {1,2,3}`:
>
>     Σ_{k=1..r}  U[a,b,k] · V[c,d,k] · W[i,j,k]  =  δ(a=i) · δ(b=c) · δ(d=j)

That is `3^6 = **729** cubic equations` in `27r` raw unknowns.

Of the 729 equations:
- **27 have RHS = 1**: the "useful" monomials. Each of the `9` output
  cells `C[i,j]` contains exactly `3` products `A[i,l]·B[l,j]` for
  `l ∈ {1,2,3}`.
- **702 have RHS = 0**: every other `A[a,b]·B[c,d]` monomial must cancel
  in every output cell.

Equivalently, this is a rank-`r` tensor decomposition

```
T = ⟨3,3,3⟩ ∈ R^{9×9×9}      with         T = Σ_{k=1..r}  u_k ⊗ v_k ⊗ w_k
```

where each `u_k, v_k, w_k ∈ R^9` (flattened 3×3 slices).

---

## 3. Counting: unknowns vs equations

After the standard quotient by scaling and the `GL_3 × GL_3 × GL_3` gauge:

|  | raw params | after `-2r` scaling | after `-27` gauge | equations | excess |
|---|---|---|---|---|---|
| `r = 19` | 513 | 475 | 448 | 729 | +281 (no solution — open whether any) |
| `r = 21` | 567 | 525 | 498 | 729 | +231 |
| `r = 22` | 594 | 550 | 523 | 729 | +206 |
| `r = 23` | 621 | 575 | 548 | 729 | +181 (solutions exist — Laderman, Smirnov…) |
| `r = 27` (naive) | 729 | 675 | 648 | 729 | +81 |

`GL_3` has dimension `9`, so the three-factor gauge eats `27` continuous
degrees of freedom. Even at `r = 23` the system is heavily overdetermined
by counting — and yet it has solutions. **Counting alone is not predictive
here** the way it is for `⟨2,2,2⟩`: the equations are far from generic.
This is precisely why the 3×3 case is hard.

A useful comparison:

| format | rank | params | equations | excess |
|---|---|---|---|---|
| `⟨2,2,2⟩` | 7 | 58 | 64 | +6 (Strassen exists) |
| `⟨3,3,3⟩` | 23 | 548 | 729 | +181 (Laderman exists) |
| `⟨3,3,3⟩` | 22 | 523 | 729 | +206 (?) |
| `⟨3,3,3⟩` | 21 | 498 | 729 | +231 (?) |

The "+181 already has solutions" point is the empirical evidence that the
3×3 equations are *very* non-generic; the search problem at `r = 22, 21`
is far from numerically hopeless.

---

## 4. Symmetries of `⟨3,3,3⟩`

The symmetry group (from [SYMMETRIES.md](../../paper/theory/SYMMETRIES.md)) is

```
G  =  (GL_3 × GL_3 × GL_3) ⋊ S_3
```

— continuous dimension `27`, discrete factor of order `6`.

### 4.1 The cyclic `Z/3` action and orbit structures

A `Z/3`-equivariant decomposition partitions its `r` rank-1 terms into
**fixed points** (size-1 orbits where `u = v = w` up to the scaling
quotient) and **orbits of 3**. The allowed structures for the two
interesting ranks:

| `r` | fixed `+` orbits-of-3 |
|---|---|
| `21` | `0 + 7`, `3 + 6`, `6 + 5`, `9 + 4`, `12 + 3`, `15 + 2`, `18 + 1`, `21 + 0` |
| `22` | `1 + 7`, `4 + 6`, `7 + 5`, `10 + 4`, `13 + 3`, `16 + 2`, `19 + 1` |

(For `r = 23` Laderman's algorithm has structure `2 + 7` or `5 + 6`
depending on the variant; Smirnov enumerated several.)

`r = 22` forces **at least one fixed triple**, since `22 mod 3 = 1`. That
fixed triple is `(u, u, u)` with `u ∈ R^9` — a single 9-vector instead of
the 27 entries of an unconstrained triple, **plus** it must be cyclically
"self-correcting" against itself in the equation system. This is a strong
constraint and a useful search prior: enumerate plausible fixed-vector
choices `u` first, then search the orbit generators.

`r = 21` has the cleanest option: `0 + 7` — purely orbits-of-3, no fixed
point. Equivalent to searching for 7 generator triples, each contributing
3 multiplications.

### 4.2 Burnside count on the 729 equations

The 729 equations themselves carry a `Z/3` action via cyclic permutation
of the index pairs `(a,b) → (c,d) → (i,j) → (a,b)`. Fixed equations are
those with `(a,b) = (c,d) = (i,j)`, of which there are `3² = 9`. By
Burnside:

```
# orbits of equations under Z/3  =  (729 + 9 + 9) / 3  =  249
```

So a `Z/3`-equivariant search has **249 distinct equations** instead of
729.

### 4.3 Equivariant parameter counts

For a `f + 3·g` structure (where `f + 3g = r`):

- `f` fixed triples: each `u ∈ R^9`, so `9·f` params.
- `g` orbit generators: each a full triple `(u,v,w) ∈ R^{27}`, so `27·g`
  params.
- Total: `9f + 27g`.

| `r` | structure | raw equivariant params |
|---|---|---|
| 21 | `0 + 7` | `0 + 189 = 189` |
| 21 | `3 + 6` | `27 + 162 = 189` |
| 21 | `21 + 0` (all fixed) | `189` |
| 22 | `1 + 7` | `9 + 189 = 198` |
| 22 | `4 + 6` | `36 + 162 = 198` |
| 22 | `19 + 1` | `171 + 27 = 198` |

Pattern: `9f + 27g = 9(f + 3g) = 9r`. So **every equivariant structure at
fixed `r` has the same raw parameter count `9r`**. The structures differ
in their *discrete shape*, not their dimension.

After gauge fixing (`-27`) and scaling quotient (`-2r` — but careful:
fixed points only have 1 scaling parameter, not 2):

- `r = 21`: `189 - 2·21 - 27 = 120` essential parameters vs `249`
  equations → excess `+129`.
- `r = 22`: `198 - (1·1 + 2·21) - 27 = 128` (for `1 + 7` structure)
  vs `249` → excess `+121`.

Still heavily overdetermined under equivariance — but the constants are
nothing like the raw counts (e.g. ALS over 128 unknowns and 249 equations
is genuinely numerically tractable).

### 4.4 The transpose `Z/2` and `S_3`

`S_3` equivariance (cyclic plus transpose) further halves things but
risks excluding solutions: Strassen itself is *not* `S_3`-equivariant on
the nose. Smirnov's `⟨3,3,3⟩` constructions exist with various levels of
discrete symmetry; the most successful searches use partial symmetry
(`Z/3` always, `S_3` opportunistically).

### 4.5 Canonicalization vs equivariance — two different things

This is the conceptual point that's easy to miss. There are two ways to
use a symmetry group `H` acting on the search space:

**Canonicalization (orbit-breaking).** Pick one representative from each
`H`-orbit (e.g. lex-smallest). The candidate `(U, V, W)` is still a
*generic* rank-`r` decomposition — you visit the same shapes as before,
just once per orbit instead of `|H|` times.

- Acts on the **traversal**, not the **shape** of candidates.
- Search-space reduction: factor `|H|` (or `|H|/|stabilizer|`).
- Dimension of the space of candidates: **unchanged**.
- Risk of missing solutions: **none** — the canonical representative
  always exists in each orbit.

**Equivariance (shape restriction).** Demand that the candidate itself
is `H`-invariant — applying any `h ∈ H` to `(U, V, W)` gives back the
same multiset of rank-1 terms. The candidate is *no longer generic*; it
is self-symmetric.

- Acts on the **shape** of candidates.
- Search-space reduction: the parameter dimension drops. For `Z/3`-
  equivariant decompositions of `⟨3,3,3⟩`: `27r` entries → `9r` entries
  (factor of 3 *in dimension*, exponentially huge in size).
- Risk of missing solutions: **yes** — if no `H`-equivariant
  decomposition of rank `r` exists but a non-equivariant one does, the
  equivariant search returns nothing.

**Concrete example for `Z/3` on a rank-7 decomposition of `⟨2,2,2⟩`:**

- *Canonicalization* would let you visit one representative per
  `Z/3`-orbit of decompositions — you still consider all 7-tuples
  `{(u_k, v_k, w_k)}`, just modulo cyclic relabeling of each term.
- *Equivariance* forces the 7 terms themselves to partition into
  `Z/3`-orbits: 1 fixed triple + 2 orbits of 3 = 7 terms. The candidate
  is described by 1 vector (`u_1 = v_1 = w_1`) + 2 generator triples,
  total `4 + 24 = 28` entries instead of `84`.

Strassen happens to be equivariant — that's why the bet pays off at
`r = 7`. Whether some `r = 22` decomposition of `⟨3,3,3⟩` is
equivariant is unknown; searching only equivariant ones is a useful
heuristic, not a complete method.

The other symmetries (`S_r` on terms, scaling, `GL_3³` gauge) are
always used as **canonicalizations** — they never restrict the shape.
Only the `S_3`-on-factors family (in particular `Z/3`) is naturally
applied as **equivariance** when one wants to compress the dimension.

### 4.6 Full symmetry hierarchy for a search loop

All four families of symmetry stack. For a `{-1, 0, +1}`-alphabet
search at `r = 21`:

| layer | symmetry | size | how used | dimension cut? |
|---|---|---|---|---|
| 1 | `S_r` on rank-1 terms | `21! ≈ 5·10^{19}` | canonicalize (lex-order columns of `U`) | no |
| 2 | scaling per term `(α,β,γ)`, `αβγ=1` | continuous, but `(±1)^{21}` preserves alphabet | canonicalize (force first nonzero of `u_k` to `+1`) | no (gauge quotient at parameter-count level only) |
| 3a | monomial subgroup of `(GL_3)³`: permutations | `(3!)³ = 216` | canonicalize (lex-min `u_1`) | no |
| 3b | monomial subgroup of `(GL_3)³`: sign flips | `(2³)³ = 512` | canonicalize | no |
| 3c | full `(GL_3)³` (continuous) | dim `27` | parameter quotient | no (gauge quotient only) |
| 4 | `Z/3` on tensor factors | `3` | **equivariance** (shape restriction) | **yes: `27r → 9r`** |
| 5 | transpose `Z/2` (giving full `S_3` on factors) | `2` | optional equivariance | yes if used |

Layer 1 (`S_r`) is your "order columns of `U`" trick: it permutes
`{1, …, r}` independently of the tensor structure.

Layer 3a is the **monomial / permutation subgroup of the gauge**.
There are exactly three independent permutation freedoms `σ_1, σ_2, σ_3 ∈
S_3`, each pairing two of the six possible row/column axes of `(A, B,
C)` and leaving the third matrix untouched:

| | acts on | acts on | leaves alone |
|---|---|---|---|
| `σ_1` | rows of `A` | rows of `C` | `B` |
| `σ_2` | cols of `A` | rows of `B` | `C` |
| `σ_3` | cols of `B` | cols of `C` | `A` |

These three pairings come straight from the Kronecker-delta structure
of the matmul tensor:

```
T[(a,b), (c,d), (i,j)]  =  δ(a = i) · δ(b = c) · δ(d = j)
                              ─────     ─────     ─────
                              σ_1       σ_2       σ_3
```

Each `δ(·)` ties two index types together — and that's exactly the
three permutation freedoms. There is no fourth one: any other paired
permutation (e.g. "rows of `A` and cols of `B`") would break `C = A·B`.
Total count: `(3!)³ = 216` discrete elements. **Note this acts on all
`r` terms simultaneously** (it permutes the row/column structure inside
each `u_k`'s 3×3 matrix), so it is *not* the same as layer 1.

Layer 4 is the unique source of *dimensional* compression: the candidate
itself is constrained to be self-symmetric, and the parameter count
drops from `27r` to `9r`.

**Composition order matters.** Canonicalizations don't commute: applying
"lex-min over `S_r`" then "lex-min over the monomial gauge" generally
doesn't give the same canonical form as the reverse order, because each
operation can disturb the previous one's canonical form. The robust
recipe is to **canonicalize over the full combined group** at the end —
or, for prefix pruning during enumeration, accept that you may revisit
some orbit members and dedupe with a hash set on the final canonical
form.

---

## 5. The `{-1, 0, +1}` alphabet

All known small-format fast algorithms have entries in a tiny finite set.
For `⟨3,3,3⟩` the typical ambient alphabets are `{-1, 0, +1}` and
sometimes `{-1, 0, +1, ±1/2, ±2}`. Restricting to `{-1, 0, +1}`:

### 5.1 Canonicalization budget

The combined size of the canonicalization group `|S_r × scaling × monomial
gauge × transpose|` determines how much enumeration the symmetries save.

**Raw (non-equivariant) at `r = 21`:**

| group | size | `log₁₀` |
|---|---|---|
| `S_r` on terms | `21! ≈ 5·10^{19}` | 19.7 |
| Per-term `±1` scaling (`αβγ=1` → 4 choices/term) | `4^{21} ≈ 4·10^{12}` | 12.6 |
| Monomial gauge `(S_3 ⋉ {±1}^3)^3` | `48³ ≈ 1.1·10^5` | 5.0 |
| Transpose `Z/2` | `2` | 0.3 |
| **total** | | **37.6** |

**`Z/3`-equivariant at `r = 21`, structure `0 + 7`:**

| group | size | `log₁₀` |
|---|---|---|
| `S_7` on orbit generators | `5040` | 3.7 |
| Cyclic choice of generator within each orbit | `3^7 = 2187` | 3.3 |
| Per-generator `±1` scaling | `4^7 ≈ 1.6·10^4` | 4.2 |
| Monomial gauge | `48³` | 5.0 |
| Transpose `Z/2` | `2` | 0.3 |
| **total** | | **16.5** |

(The reduction is smaller in absolute terms in the equivariant setting
because the raw space is itself much smaller — the *ratio* relative to
the raw space is similar.)

### 5.2 Candidate counts after canonicalization

| setting | raw count | after canonicalization |
|---|---|---|
| `r = 21`, no equivariance | `3^{567} ≈ 10^{270}` | `≈ 10^{232}` |
| `r = 21`, `Z/3`-equivariant `0 + 7` | `3^{189} ≈ 10^{90}` | **`≈ 10^{73}`** |
| `r = 22`, no equivariance | `3^{594} ≈ 10^{283}` | `≈ 10^{243}` |
| `r = 22`, `Z/3`-equivariant `1 + 7` | `3^{198} ≈ 10^{94}` | **`≈ 10^{78}`** |

So the most-compressed honest count for `r = 21` is **~10⁷³ canonical
candidates** — still vastly beyond brute-force enumeration (the
observable universe contains ~10⁸⁰ atoms). The reduction from raw
ternary enumeration is dominated by the `Z/3`-equivariance dimension
cut (`~10^{270} → ~10^{90}` — eighteen orders of magnitude per layer of
the three factors), with canonicalization contributing a further
~17 orders of magnitude.

### 5.3 Practical implications

For a brute-force feasibility comparison with the much smaller `⟨2,2,2⟩` at
`r=7` problem (where the canonical-candidate count comes out to ~`10^{7.9}` —
seconds on a laptop, which is why this repo's `v1`/`v2` enumeration finishes),
see [SOLVING_STRATEGIES.md](../../paper/theory/SOLVING_STRATEGIES.md) §2.1.

`10^{73}` is "solar-system-of-atoms" territory — no enumeration strategy
can ever touch it directly. What makes the problem tractable in
practice is *not* further symmetry pruning but **propagation**: each of
the 249 cubic equations, once partially evaluated against a fixed
prefix of variables, eliminates large branches of the search tree. A
modern SAT/SMT solver, ALS run, or AlphaTensor-style policy network
navigates this `10^{73}`-candidate space without ever materializing
even a tiny fraction of it.

A useful mental anchor: AlphaTensor's RL agent visited ~`10^{12}`
candidate moves over its entire training run and rediscovered Laderman
at `r = 23`. The `10^{73}` count is the *size of the formal search
space*, not the work that needs to be done — but it does mean that any
"just enumerate, with smart pruning" approach has to discard at least
~60 orders of magnitude of the space *per evaluated candidate* to
finish in reasonable time. That's the bar a search heuristic has to
clear.

The takeaway: pure enumeration is infeasible at `r = 21, 22`, but
**guided search** in a `~10^{73}`-candidate space is the right mental
picture. Any successful approach needs:

1. equivariance to chop the symmetry orbit;
2. gauge fixing to canonicalize residual freedoms;
3. a learned or hand-coded heuristic to traverse the tree;
4. SAT/SMT or ALS to certify candidate completions.

---

## 6. Search pseudocode

### 6.1 Verifier (foundation for any search)

```python
import numpy as np

def matmul_tensor(n=3):
    T = np.zeros((n*n, n*n, n*n))
    for i in range(n):
        for j in range(n):
            for l in range(n):
                T[i*n + l, l*n + j, i*n + j] = 1.0
    return T

def verify(U, V, W, T):
    """U, V, W are shape (9, r). Returns Frobenius residual."""
    approx = np.einsum('ak,bk,ck->abc', U, V, W)
    return np.linalg.norm(T - approx)
```

### 6.2 `Z/3`-equivariant branch-and-bound at `r = 21`, structure `0 + 7`

```python
def equivariant_search_21(alphabet=(-1, 0, 1)):
    T = matmul_tensor(3)
    n_orbits = 7

    # Each orbit generator is a triple (u, v, w) ∈ R^9 × R^9 × R^9.
    # Use gauge fixing on the first generator to break GL_3³ orbit
    # redundancy (e.g. fix u_1 to a canonical form).

    for gen_seq in canonical_orbit_generators(n_orbits, alphabet):
        U, V, W = expand_orbits(gen_seq)             # each shape (9, 21)
        if verify(U, V, W, T) < 1e-10:
            yield (U, V, W)

def expand_orbits(generators):
    """Each generator (u, v, w) yields three terms under Z/3:
       (u, v, w), (v, w, u), (w, u, v)."""
    U_cols, V_cols, W_cols = [], [], []
    for (u, v, w) in generators:
        U_cols += [u, v, w]
        V_cols += [v, w, u]
        W_cols += [w, u, v]
    return (np.stack(U_cols, axis=1),
            np.stack(V_cols, axis=1),
            np.stack(W_cols, axis=1))
```

`canonical_orbit_generators` is where every practical optimization
lives — gauge fixing, lex-min canonicalization, prefix pruning against
partial 249-equation residuals, isomorph-rejection, etc.

### 6.3 `Z/3`-equivariant ALS at `r = 22`, structure `1 + 7`

```python
def equivariant_als_22(n_restarts=10_000, n_iter=5_000):
    T = matmul_tensor(3)
    best = (np.inf, None)
    for restart in range(n_restarts):
        # 1 fixed triple: vector u_fix ∈ R^9
        # 7 orbit generators: each (u, v, w) ∈ R^9 × R^9 × R^9
        u_fix = np.random.randn(9)
        gens  = [tuple(np.random.randn(9) for _ in range(3)) for _ in range(7)]

        for it in range(n_iter):
            # Standard ALS but projected onto Z/3-equivariant subspace.
            # i.e. solve for u_fix and gens such that
            # expand([u_fix-fixed] + [orbits]) reconstructs T.
            u_fix, gens = als_step_z3_equivariant(T, u_fix, gens)

        U, V, W = expand_full(u_fix, gens)
        r = verify(U, V, W, T)
        if r < best[0]:
            best = (r, (U, V, W))
            if r < 1e-10:
                print(f"FOUND rank-22 algorithm at restart {restart}")
                return best
    return best
```

Expected behavior, based on the literature:
- At `r = 23`: converges with non-trivial probability per restart.
- At `r = 22`: residual plateaus tantalizingly low (no proof it's
  impossible, but no known successful run either as of this writing).
- At `r = 21`: even harder, but the `0 + 7` orbit structure is the most
  studied target. If found, would simultaneously beat Laderman *and*
  recursive Strassen.

### 6.4 SAT/SMT formulation sketch

Over `{-1, 0, +1}`, each scalar unknown is encoded as 2 booleans
(`sign` + `nonzero`). The 249 orbit-equation polynomials are cubic; expand
them into CNF via Tseitin and feed to a modern SMT solver
(e.g. Z3, cvc5, OR-tools) or a finite-field-aware solver
(e.g. CryptoMiniSat over GF(2) for the parity-encoded version).

Variable counts (equivariant, `{-1,0,1}`-encoded, before propagation):
- `r = 21, 0+7`: `189 · 2 ≈ 380` booleans.
- `r = 22, 1+7`: `198 · 2 ≈ 400` booleans.

This is well within solver range *per problem instance*; the bottleneck
is the 249 cubic constraints, each expanding into ~thousands of CNF
clauses after Tseitin. Real-world experience (Heule et al.; Sedoglavic's
work) suggests these encodings can be made to terminate in hours-to-days
on a single machine for the right structural restrictions.

---

## 7. Why this is the open question

`⟨3,3,3⟩` is special:

1. It's the smallest format where the true rank is **not known**.
2. The 23 upper bound has stood since **1976**.
3. The Coppersmith–Winograd and laser-method machinery gives the best
   *asymptotic* `ω` bounds, but those don't yield small-format algorithms.
4. AlphaTensor (2022) rediscovered Laderman's 23 but did *not* improve it
   in standard arithmetic — strongly suggesting that the `r = 22`
   threshold is genuinely hard (or perhaps the rank actually is 23 and we
   need a lower-bound advance).

For this repo specifically, the natural progression of experiments is:

1. Implement the verifier (§6.1).
2. Reproduce Laderman's 23-mul algorithm from a hard-coded `(U, V, W)`
   table and verify it.
3. Implement `Z/3`-equivariant ALS targeting `r = 23` first
   (to confirm the framework works), then `r = 22`.
4. Implement the SAT/SMT encoding with `Z/3` equivariance and
   `{-1, 0, +1}` alphabet at `r = 22`. Even a **negative** result under
   restricted alphabet + symmetry is publishable; a **positive** result
   would be a 50-year breakthrough.

---

## 8. Comparison summary

| question | `⟨2,2,2⟩` (`STRASSEN_AS_EQUATIONS`, `RANK_5_SEARCH`) | `⟨3,3,3⟩` (this doc) |
|---|---|---|
| true rank | known = 7 (de Groote, Hopcroft–Kerr) | open: `19 ≤ R ≤ 23` |
| best known algorithm | Strassen 1969 | Laderman 1976, Smirnov 2013/17 |
| `# equations` | 64 | 729 |
| `# equations under Z/3` | 24 | 249 |
| raw params at known rank | 84 (`r=7`) | 621 (`r=23`) |
| equivariant params at known rank | 28 | `9·23 = 207` |
| beats-Strassen threshold | n/a | `r = 21` (`log_3 21 < log_2 7`) |
| beats-current-art threshold | n/a | `r = 22` |
| `{-1,0,1}` brute-force space | `3^{84} ≈ 10^{40}` | `3^{621} ≈ 10^{296}` |
| equivariant `{-1,0,1}` space | `3^{28} ≈ 10^{13}` | `3^{207} ≈ 10^{99}` |

The combinatorial blow-up between the two formats is roughly **`10^{86}`**
in equivariant `{-1,0,1}` enumeration, which is exactly why `⟨3,3,3⟩`
remained open while `⟨2,2,2⟩` fell in a single 1969 paper.
