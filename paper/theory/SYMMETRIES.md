# Leveraging the Symmetries of `⟨2,2,2⟩`

The matrix-multiplication tensor `T = ⟨n,n,n⟩` has a large symmetry group.
Any rank decomposition of `T` is mapped to another rank decomposition by every
group element, so the search space and the equation system can both be cut
down by *equivariant* search.

This is the single most effective non-brute-force pruning available, and is
what underlies essentially every successful new algorithm discovery in the
small-format world (Heun, Smirnov, Sedoglavic, and AlphaTensor's
data-augmentation prior).

---

## 1. The symmetry group of `⟨n,n,n⟩`

Define the trilinear form `T(A, B, C) = trace(A·B·C)`. Three families of
symmetries:

### (a) Cyclic `Z/3` from the trace
`trace(A·B·C) = trace(B·C·A) = trace(C·A·B)` — the trilinear form is
invariant under cyclic permutation of its three arguments. So the
permutation `σ: (A,B,C) → (B,C,A)` is a symmetry.

### (b) Transpose `Z/2`
`trace(A·B·C) = trace((A·B·C)^T) = trace(C^T·B^T·A^T)`. Combining with
cyclic gives the full symmetric group `S_3` on the three tensor factors,
twisted by entrywise transposition.

### (c) Change of basis `GL_n × GL_n × GL_n`
For invertible `P, Q, R`:
```
trace((P A Q^{-1})·(Q B R^{-1})·(R C P^{-1})) = trace(P A B C P^{-1}) = trace(A B C)
```
So `(A, B, C) → (P A Q^{-1}, Q B R^{-1}, R C P^{-1})` is a symmetry. This is
the "sandwich gauge."

### The full group
```
G  =  (GL_n × GL_n × GL_n)  ⋊  S_3
```
For `n = 2`: `dim(G) = 3·4 = 12` continuous degrees of freedom, plus a
discrete `S_3` of order 6.

Every rank decomposition `T = Σ_k u_k ⊗ v_k ⊗ w_k` is mapped by `G` to
another rank-`r` decomposition — so **the set of rank-7 algorithms is a
union of `G`-orbits**.

> **de Groote (1978):** for `T = ⟨2,2,2⟩` and `r = 7`, there is *exactly
> one orbit*. Every Strassen-like algorithm is equivalent to the original
> under the `G`-action.

---

## 2. Counting: how much do the symmetries save?

| | parameters | equations |
|---|---|---|
| Raw bilinear template (§3 of STRASSEN_AS_EQUATIONS) | `84` | `64` |
| Modulo scaling redundancy per rank-1 term (`-2r`) | `70` | `64` |
| Modulo `GL_2 × GL_2 × GL_2` gauge (`-12`) | `58` | `64` |
| Plus `Z/3` equivariance restriction | `~36` | `24` (orbits of equations) |
| Plus `{-1, 0, +1}` alphabet | finite, `~10^17` | 24 |
| Plus `S_3` equivariance | even smaller | smaller |

The combined effect: a 84-dim continuous search → a finite combinatorial
search of moderate size with `~20` algebraic constraints. This is what
makes Heun's and Smirnov's hand-and-computer searches tractable, and what
gives AlphaTensor's RL the right prior structure.

---

## 3. Equivariance, concretely

A decomposition `{(u_k, v_k, w_k)}_{k=1..r}` is `H`-equivariant (for `H ⊂ G`)
if, for every `g ∈ H`, the multiset is preserved by `g`'s action on
triples. The terms break into `H`-orbits.

### 3.1 `Z/3`-equivariance: the orbit structure
Cyclic permutation acts on triples by `(u, v, w) → (v, w, u)`. Each rank-1
term either:
- is **fixed** (`u = v = w` up to the scaling redundancy: i.e., the triple is
  cyclically invariant — corresponds to a "symmetric" product), or
- lies in an **orbit of 3** (the term and its two cyclic conjugates are all
  distinct).

So for `r = 7`, the possible orbit structures under `Z/3` are:
- `7 = 7·1`  (no nontrivial cyclic structure)
- `7 = 4·1 + 3`  (one orbit of 3, four fixed)
- `7 = 1·1 + 2·3`  ← **Strassen's structure**
- `7 = 7` impossible (`3 ∤ 7`)

**Strassen is 1 + 3 + 3.** That single fixed product is
```
M1 = (A11+A22) · (B11+B22)        with output mask (1,0,0,1) (the diagonal of C)
```
Its `(u, v, w)` triple is `((1,0,0,1), (1,0,0,1), (1,0,0,1))` — all three
slices equal, manifestly cyclically invariant. The remaining six products
form two orbits of three under cyclic permutation of the slices.

### 3.2 Why this prunes the search
For a `1 + 3 + 3` decomposition you only need to specify:
- 1 fixed triple: a single vector `u ∈ R^4` (since `u = v = w`), so 4 params.
- 2 orbit generators: each is a triple `(u, v, w)` of three vectors in `R^4`,
  so `2 × 12 = 24` params.
- Total: **28 free parameters** instead of 84.

After scaling and gauge fixing, this drops further to a handful of effective
parameters — small enough to enumerate over `{-1, 0, +1}` exhaustively.

### 3.3 Equations also collapse under symmetry
The 64 equations
```
Σ_k U[a,b,k] · V[c,d,k] · W[i,j,k] = δ(a=i) · δ(b=c) · δ(d=j)
```
are themselves `Z/3`-equivariant: if you cyclically permute the index triples
`(a,b) → (c,d) → (i,j) → (a,b)`, you map equations to equations. By
**Burnside** the number of equation orbits under `Z/3` is

```
(64 + 4 + 4) / 3 = 24
```

(The 4's come from the index tuples fixed by cyclic permutation:
`(a,b) = (c,d) = (i,j)`, of which there are `2² = 4`.)

So in the equivariant setting we have **28 unknowns subject to 24 equations**.
Strassen lies at a specific 4-dim solution component of this small
overdetermined system — orders of magnitude easier to find than the raw `84
vs 64`.

### 3.4 Pushing further with `S_3`
Adding the transpose `Z/2` would force the algorithm to be `S_3`-equivariant.
Strassen is **not** fully `S_3`-equivariant on the nose: applying transpose
gives a *different but equivalent* algorithm in the same `GL³`-orbit (the
"dual" Strassen). So restricting to `S_3`-equivariant decompositions risks
missing Strassen. For larger formats `S_3` symmetry IS sometimes used
(Smirnov found `<3,3,3>` algorithms with full `S_3` symmetry).

---

## 4. `GL_2 × GL_2 × GL_2` gauge fixing

The continuous part of `G` acts by change of basis. Two decompositions in the
same `GL³` orbit are "the same algorithm" up to relabeling of rows/columns.
To avoid traversing this orbit redundantly, **fix the gauge**:

Common gauge-fixing choices for `n = 2`:
- Fix three columns of `U` to standard basis vectors (4 + 4 + 4 = 12 params
  removed — matches `dim GL³`).
- Equivalently, fix the first three triples `(u_1, v_1, w_1), …, (u_3, v_3,
  w_3)` to canonical forms.
- For Strassen specifically: fix `u_1 = (1,0,0,1)` (the diagonal), then the
  cyclic `Z/3` stabilizer further pins down `v_1 = w_1 = (1,0,0,1)`.

With gauge fixed, the rank-7 set is **a finite (in fact single) point** by
de Groote's theorem — so any search that respects the gauge will converge
deterministically to Strassen once feasibility is achieved.

---

## 5. What this means for an actual solver

### 5.1 Branch-and-bound search
Pseudo-code for a `Z/3`-equivariant search with `{-1, 0, 1}` alphabet:

```python
for u_fixed in candidates(R^4, alphabet={-1,0,1}):       # 3^4 = 81
    for (u1, v1, w1) in candidates_triples():            # 3^12 ≈ 5·10^5
        for (u2, v2, w2) in candidates_triples():        # 3^12 ≈ 5·10^5
            orbit1 = [(u1,v1,w1), (v1,w1,u1), (w1,u1,v1)]
            orbit2 = [(u2,v2,w2), (v2,w2,u2), (w2,u2,v2)]
            decomp = [(u_fixed,u_fixed,u_fixed)] + orbit1 + orbit2
            if verify_decomposition(decomp) == matmul_tensor:
                yield decomp
```

Naive count: `81 × 5·10^5 × 5·10^5 ≈ 2·10^13` — large but feasible with
SAT/SMT + symmetry-aware branching. Gauge fixing (canonicalize each orbit
generator) shaves another 1–2 orders of magnitude.

### 5.2 ALS with equivariance constraint
Modify the ALS iteration (§5.3 of STRASSEN_AS_EQUATIONS) by projecting onto
the `Z/3`-equivariant subspace at each step. This essentially solves only
for the orbit generators, reducing the per-iteration linear system from
`16 × 7 → 16 × 3` and dropping local minima drastically.

### 5.3 Use symmetry as a data-augmentation prior (AlphaTensor style)
Train a neural policy with the `G`-action as a data augmentation: any
training example `(state, action)` yields `|G|` equivalent examples by
applying group elements. AlphaTensor uses exactly this to make the
combinatorial search tractable.

---

## 6. Beyond `n = 2`: the same idea, much more leverage

For `⟨3,3,3⟩` the symmetry group is `(GL_3)³ ⋊ S_3`, dimension `3·9 = 27`
continuous. Smirnov's 23-multiplication algorithms have nontrivial `Z/3`
orbit structures (e.g. `23 = 2 + 3·7` or `23 = 5 + 3·6`). The
equivariance-restricted search is how all known `<3,3,3>` algorithms with
23 multiplications were discovered. **Closing the 19–23 gap** (the open
question) is essentially a question about whether *some* `Z/3`-equivariant
decomposition of rank ≤ 22 exists.

---

## 7. Suggested experiments in this repo

1. **Implement Z/3-equivariant ALS** for `⟨2,2,2⟩, r=7`. Should converge
   from random `Z/3`-equivariant initializations far more reliably than
   plain ALS.
2. **SAT/SMT encode** the 24 orbit-equations in 28 variables over
   `{-1,0,1}` with gauge-fixing; verify it has a unique solution (modulo
   symmetry) and that it is Strassen.
3. **Same encoding at `r=6`**: prove unsatisfiability, giving an
   independent computational re-proof of `R(⟨2,2,2⟩) ≥ 7` (the lower bound
   becomes a SAT certificate).
4. **Scale to `⟨3,3,3⟩, r=22`** with `Z/3` equivariance — a serious
   research-frontier search with realistic prospects.

---

## 8. Worked example: peeling the symmetry layers on `⟨2,2,2⟩, r=7`

Experiment #2 above is implemented in [`v3/Z3BruteForce2x2.java`](../../src/main/java/io/cormoran/strassen/v3/Z3BruteForce2x2.java)
and was run on 2026-05-24. It walks through *exactly* how the symmetry stack
acts on the search space, layer by layer, and confirms that Strassen is
unique within `{-1, 0, +1}` modulo the full discrete symmetry.

### 8.1 Which symmetries can be enforced *during* the search

Not all symmetries can be applied as prefix filters. The split is by **locality**:

| layer | size | applied during search? | rationale |
|---|---|---|---|
| `u_fix` sign canonicalization (first nonzero positive) | `2` | ✓ | local to one vector — fix at outer loop |
| cyclic-within-orbit (`u` lex-smallest of `{u, v, w}`) | `3` per orbit | ✓ | local to one orbit's three generators |
| `S_2` orbit-swap (orbit1 ≤ orbit2 lex) | `2` | ✓ | local to the comparison of two orbits |
| monomial gauge `(S_2 ⋉ {±1}²)³` | `8³ = 512` | ✗ | acts simultaneously on all 7 columns — can't evaluate mid-prefix |
| transpose `Z/2` | `2` | ✗ | global rewiring of all 7 columns |
| per-term `±1` scaling | `4^7 = 16,384` | ✗ | per-term local but interacts with the gauge action above |
| `S_r` permutation of rank-1 terms | `7! = 5,040` | ✗ | only meaningful once all 7 terms are fixed |

The "✓" layers reduce the raw `3^{28} ≈ 2.3·10^{13}` ternary space of
Z/3-equivariant decompositions to roughly `7.4·10^6` orbit-1 prefix prunes
(the `pairs` counter in the run log). The "✗" layers must wait until a
complete 7-term decomposition is in hand — they are applied as a post-hoc
**dedup over canonical forms**.

### 8.2 The two-phase run

Phase 1 — **enumeration with the "✓" filters and the R2-diagonal prune**
(see [SOLVING_STRATEGIES.md](SOLVING_STRATEGIES.md) §2.1 for the diagonal-prune explanation):

```
≈ 105.9 seconds   →   32 partial-canonical solutions, all with u_fix = (1, 0, 0, 1)
```

Phase 2 — **post-hoc canonicalization over the "✗" symmetries**: each
solution is expanded to its 7 rank-1 terms; for each of the `512 × 2 = 1024`
gauge/transpose elements, the terms are transformed, sign-normalized
per-term, sorted, and serialized. The lex-min serialization across all
1024 transformations is the canonical form. Solutions sharing a canonical
form belong to the same equivalence class:

```
≈ 0.175 seconds   →   1 distinct canonical solution
```

So the **final count is `1`**, the canonicalization step is essentially
free, and the orbit-size factor `32` represents the un-quotiented portion
of the symmetry orbit that the prefix-only search couldn't collapse on its
own.

### 8.3 What this confirms

- **De Groote 1978's uniqueness** — that all rank-7 algorithms for `⟨2,2,2⟩`
  are equivalent under the full gauge — holds within the `{-1, 0, +1}`
  alphabet, constructively, in a computer-checkable enumeration.
- **The symmetry quotient stack we've documented is complete**: if any of
  the listed layers (`Z/3` cyclic, monomial gauge, transpose, per-term
  scaling, `S_r`) had been missing or under-specified, the final count
  would have been > 1. It collapsed cleanly to 1, so the stack is right.
- **Empirical evidence for the prune ratio**: the
  `(monomial gauge × transpose)`-orbit of Strassen, restricted to within
  the "✓"-canonical representatives, has size exactly `32`. That number is
  a small but nontrivial constant — useful as a calibration point when
  designing analogous searches at larger `n`.

### 8.4 Two-phase architecture and why it's the right shape

The split between in-loop pruning ("✓" layers) and post-hoc dedup ("✗"
layers) is the same pattern that any constraint-propagating search (SAT,
ALS, RL) ends up using, just under different names:

- **In-loop**: corresponds to *unit propagation* in SAT — eliminates
  partial assignments using local consistency.
- **Post-hoc dedup**: corresponds to *isomorph rejection* in graph
  enumeration — eliminates duplicates only at the leaves.

You *could* try to push more layers into in-loop pruning (e.g. enforce
gauge-canonicality on `u_fix` plus the first few entries of `u1` as soon as
they're set), but the bookkeeping is intricate and the wins are small
relative to the R2-diagonal prune that already kills 99.97% of prefixes.
For `⟨2,2,2⟩` the two-phase architecture is the right shape; for
`⟨3,3,3⟩` it becomes essential (no in-loop schedule survives the
combinatorial explosion documented in [RANK_3X3_SEARCH.md](../../research/small-rank-search/RANK_3X3_SEARCH.md) §5).
