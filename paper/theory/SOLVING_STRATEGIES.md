# Solving Strategies for Bilinear-Algorithm Search

How do we actually *find* a rank-`r` decomposition of a matmul tensor
`T = ⟨n,n,n⟩`? This note surveys the methods, what each one can and cannot
prove, and how they map onto the code in this repo.

The short answer to "is brute-force out of the picture?" — **yes, for any rank
≥ 7 over discrete alphabets, and trivially yes over `R`** (which is
uncountable). But the *reason* matters: brute force is killed by combinatorics,
not by the existence of "convergent" methods. ALS converges to a local
minimum and gives no global guarantee, so it does not replace exhaustive
search — it just happens to be the only thing that scales.

---

## 1. The search problem

For each candidate rank `r`, we want triples of matrices `(U, V, W) ∈
(R^{n²×r})³` that satisfy the `n^6` cubic equations

```
∀ a,b,c,i,j,d ∈ {1..n}:    Σ_{k=1..r} U[ab,k] · V[cd,k] · W[ij,k]
                        =  δ(a=i) · δ(b=c) · δ(d=j)
```

(see [STRASSEN_AS_EQUATIONS.md](STRASSEN_AS_EQUATIONS.md) and [RANK_3X3_SEARCH.md](../../research/small-rank-search/RANK_3X3_SEARCH.md)).

Two qualitatively different problem regimes:

| regime | unknowns | what we want |
|---|---|---|
| **continuous over `R`** | `27r` real numbers (for `n=3`) | any real `(U,V,W)` satisfying the equations |
| **discrete over a finite alphabet** | `27r` ternary digits (`{-1,0,1}`) | any alphabet-bounded `(U,V,W)` satisfying the equations |

The continuous problem is the one ALS and gradient methods attack. The discrete
problem is the one SAT/SMT and pure enumeration attack. Solutions in one
regime usually don't transfer to the other directly, but a continuous solution
with rational entries is the easy bridge from one to the other.

---

## 2. Brute-force enumeration

For `n=3, r=21` over `{-1, 0, +1}` the search space size, after every available
symmetry reduction, is roughly **`10^73` canonical candidates**
(see [RANK_3X3_SEARCH.md](../../research/small-rank-search/RANK_3X3_SEARCH.md) §5.2). That is "solar-system-of-atoms" territory —
no enumeration strategy, however clever the pruning, can ever materialize even
a tiny fraction of it.

Brute force is *only* tractable when:

1. **rank is tiny**: at `n=2, r=4` over `{-1, 0, +1}` you have `3^{48} ≈ 10^{22}`
   candidates raw, ~`10^{16}` after symmetry — borderline feasible with massive
   parallelism plus aggressive prefix pruning.
2. **the search space is heavily structured by hand**: the AlphaTensor team
   often constrains entries to come from a tiny pre-selected set (e.g. just
   `{0, ±1}`), with most entries pinned to `0` by a sparsity prior. This is
   still "search" but over a space much smaller than the full alphabet count.

For `n=3` at the interesting ranks (21, 22), brute force is **definitively
out of the picture**, full stop. Even SAT/SMT (§4) doesn't enumerate — it
*propagates*, which is a different beast.

The repo's `v1` and `v2` packages are essentially exhaustive enumerations
over the small `⟨2,2,2⟩` case — they finish only because `n=2` is small
enough.

### 2.1 Worked example: `⟨2,2,2⟩` at `r = 7` (Strassen-rank)

This is the regime where this repo's `v1`/`v2` packages live. Numbers
parallel the `r=21` table in [RANK_3X3_SEARCH.md](../../research/small-rank-search/RANK_3X3_SEARCH.md) §5.

**Raw parameter count.** Three factors `U, V, W ∈ R^{4×7}` → `3 · 4 · 7 = 84`
entries. Over `{−1, 0, +1}`:

```
raw candidates  =  3^84  ≈  1.6·10^{40}
```

That alone is universe-of-atoms territory and brute-forceable by nothing.

**`Z/3`-equivariant parametrization.** Strassen has structure `1 + 2`: one
fixed triple (single vector in `R^4`) plus two orbit generators (each a triple
of vectors in `R^4`). Entry count:

```
1 · 4   +   2 · (3 · 4)   =   28
equivariant candidates    =   3^{28}  ≈  2.3·10^{13}
```

A factor `~10^{27}` cut purely from the dimension reduction. Still 23
trillion candidates — too many to verify one by one (would be ~year-of-CPU
single-threaded at 1 μs each), but tractable with parallelism.

**Canonicalization stack.** Stacking the discrete symmetries (in the order from
[RANK_3X3_SEARCH.md](../../research/small-rank-search/RANK_3X3_SEARCH.md) §4.6) on the equivariant `1+2` structure for `n=2`:

| layer | group | size | `log₁₀` |
|---|---|---|---|
| `S_g` on orbit generators (`g=2`) | `2!` | `2` | 0.3 |
| Cyclic choice of generator within each orbit | `3^g` | `9` | 1.0 |
| Per-orbit `±1` scaling (`αβγ=1` → 4 choices) | `4^g` | `16` | 1.2 |
| Per-fixed scaling | `1` (no freedom: `α³=1 ⇒ α=1`) | `1` | 0 |
| Monomial gauge `(S_2 ⋉ {±1}²)³` | `8³` | `512` | 2.7 |
| Transpose `Z/2` | `2` | `2` | 0.3 |
| **total** | | **`294,912`** | **5.47** |

So:

```
canonical candidates  =  10^{13.36} / 10^{5.47}  ≈  10^{7.89}  ≈  7.8·10^7
```

**Wall-clock.** Each canonical candidate needs an equation check (64 cubic
equations for `n=2`). At ~10 ns per equation eval that's ~640 ns per
candidate, so:

```
7.8·10^7 candidates  ×  640 ns  ≈  50 seconds single-threaded
```

Embarrassingly parallel, so realistically **single-digit seconds on a laptop
with a dozen cores**. This is why the `v1`/`v2` enumeration ever finishes.

**Empirical confirmation** (run via [`v3/Z3BruteForce2x2.java`](../../src/main/java/io/cormoran/strassen/v3/Z3BruteForce2x2.java)):

| stage | count | wall-clock |
|---|---|---|
| Enumeration with partial canonicalization (`u_fix` sign, cyclic-within-orbit, S₂ orbit-swap) and the `R2`-diagonal prune | **32** | ~105 s |
| After also quotienting by the monomial gauge `(S_2 ⋉ {±1}²)³`, transpose `Z/2`, and per-term `±1` scaling | **1** | <1 s |

So in the `⟨2,2,2⟩`-at-`r=7` regime, Strassen is **literally unique** over
`{-1, 0, +1}` modulo the full discrete symmetry — there is no alternative
ternary-alphabet rank-7 algorithm. This is a constructive lattice-level proof
of de Groote 1978 within the alphabet bound, and a sanity check that the
symmetry quotient stack in §4 of [RANK_3X3_SEARCH.md](../../research/small-rank-search/RANK_3X3_SEARCH.md) is complete. (Had the
final count been > 1, we'd be missing a symmetry.)

The R2-diagonal prune is what makes this finish in minutes. The orbit-2
contribution at the diagonal point `(i, i, i)` equals `3 · u₂[i] · v₂[i] · w₂[i]`,
because the three cyclic shifts of a rank-1 orbit term collapse to the same
scalar product when all three indices coincide. Over `{-1, 0, +1}` that forces
`R2[i,i,i] ∈ {-3, 0, +3}` for all `i` — and that filter kills roughly 99.97% of
`(u_fix, orbit1)` prefixes before any orbit-2 candidate is tried.

### 2.2 Performance experiments — what worked, what didn't

While implementing the brute force in [`v3/Z3BruteForce2x2.java`](../../src/main/java/io/cormoran/strassen/v3/Z3BruteForce2x2.java)
we tried four optimizations beyond the algorithmic R2-diagonal prune. Three
helped; two of those helped a lot; one was equal; one was a slight regression.
Recorded so future code reviewers don't re-derive the negative results.

Baseline: post-hoc dedup over the 32-solution enumeration result.

| optimization | wall-clock | speedup vs baseline | verdict |
|---|---|---|---|
| Baseline: `double[][][]` tensors, post-hoc dedup | **96.8 s** | 1.0× | — |
| `double` → `int` tensors (exact compare, no epsilon) | **76.7 s** | **1.26×** | ✓ kept |
| `int` + in-loop per-term-sign canonicalization on orbit1 | **22.5 s** | **4.3×** | ✓ kept |
| `int` → `byte` tensors | 25.2 s (vs 22.5 s) | **0.89×** | ✗ reverted |
| Bit-packed 17 MB orbit2-contribution LUT (4 longs per tensor, XOR compare) | 24.9 s (vs 22.5 s) | **0.95×** | ✗ reverted |
| **Multithreaded via `ExecutorService.newWorkStealingPool`** (on top of `int` + in-loop sign canon) | **5.1 s** | **18.6×** cumulative | ✓ kept |

**Why `int` over `double` (+26%):** every value in the trilinear matmul tensor
and residuals is an exact integer in `[-3, 3]`. Floating-point comparisons
required epsilon checks like `Math.abs(d - 3) < 0.5`; integer compares are
exact, branch-predictable, and use the 32-bit ALU path which JIT optimizes
better. The original `double` was an unthinking inheritance from the ALS code
in `v3/Als.java` (which genuinely *needs* reals).

**Why in-loop sign canonicalization wins (+3.4×):** filtering each orbit1 by
its per-term-sign canonical representative *during* enumeration cuts the
inner-loop count by ~4× (1.99 B verifications instead of 8.03 B). The filter
itself is local and cheap. The matching gauge-canonicalization can't be
applied in-loop — the gauge acts on the whole `(u_fix, orbit1, orbit2)`
decomposition — so it stays in post-hoc dedup. See [SYMMETRIES.md](SYMMETRIES.md) §8 for
which symmetries are prefix-local and which aren't.

**Why `byte` lost (–11%):** Java has no native `byte` arithmetic. Every byte
operand is sign-extended to `int` on the operand stack, computed as `int`,
then cast back when written. For our L1-cache-resident hot data (256-byte
R2 tensor + tiny vectors), the cache-savings from `byte` are zero and the
implicit-widening overhead is pure cost. `byte[][][]` arrays also retain
the ~16-byte object header per leaf row, so the in-memory footprint is much
closer to `int[][][]` than the naive 4× ratio suggests.

**Why the bit-packed LUT lost (–5%):** the obvious win — precompute all
`81³ = 531k` orbit-2 contribution tensors, pack each as 4 longs (32 bytes),
and replace the inner loop's 64 multiply-compares with a 4-long XOR — fails
because the LUT is **17 MB**. That doesn't fit in L1 (32 KB) or L2 (~1 MB),
so every lookup is an L3-hit (~10 cycles), and the access pattern across
adjacent inner iterations is strided by ±6561 LUT entries (≈ ±200 KB),
defeating both spatial locality and the hardware prefetcher.

In the int version, by contrast, R2 (256 bytes), `u2/v2/w2` (48 bytes each),
and the per-index candidate arrays are all permanently resident in L1.
Per-iteration ops are slower than 4 long XORs (~25 cycles vs ~10), but every
load is ~4 cycles instead of ~10 — and the int version benefits from
JIT-friendly early exit in `orbitMatches` (most candidates fail on the
first few entries).

The lesson: **before optimizing arithmetic, check whether the data structure
escapes the cache.** A 33×-arithmetic-op-reduction is worth nothing if it
forces an extra L3 hit per iteration.

**Why multithreading wins (+5.2× on 16 cores):** the `(u_fix, u1)` grid
naturally factors into ~3000 independent work units, each scanning its own
`(v1, w1)` sub-grid and orbit-2 enumeration with no shared mutable state.
We use `Executors.newWorkStealingPool` (Java's `ForkJoinPool` under the
hood) rather than parallel streams to avoid lambda/iterator overhead in
the dispatch loop. The 5.2× / 16 cores ≈ 32% scaling efficiency is
limited by **load imbalance**: `u_fix = (1, 0, 0, 1)` carries ~99% of the
work, so only its ~81 `u1` sub-tasks contribute meaningfully to parallel
throughput. Work-stealing partially hides this, but the achievable
parallelism is bounded by the count of "heavy" tasks, not by the core
count.

To get closer to linear scaling, the per-task granularity would need to
shrink — e.g., parallelize across `(u_fix, u1, v1)` triples, giving
~250k tasks. That would expose enough parallelism but each task would be
microseconds long, and dispatch overhead would start to dominate. The
current granularity is the right trade-off for ~10-32 cores.

**What might actually win further** (not pursued — current 5 s is fine
for the experiment we wanted):

- Finer-grained task partitioning combined with task batching to keep
  dispatch overhead amortized.
- Shrink the bit-packed LUT below the L2 boundary by storing only
  cyclic-canonical orbit-2 contributions (`81³ / 3 ≈ 5.7 MB` — still too
  big for L2 on most CPUs). Combined with sign-canonicalization,
  `~1.4 MB` would fit L2. Requires extra per-iteration canonical-form
  mapping; net gain unclear.
- Java Vector API (incubator) for SIMD-style 4-long compares. JIT support
  is uneven; engineering effort high.
- Move to a compiled language (Rust/C++) where data layouts and SIMD are
  under direct control. Likely 3–5× from the same algorithm.

None of these would change the algorithmic story or affect what the
brute force *proves* — they'd just lower the wall-clock further.

**Compare `⟨3,3,3⟩` at `r=21`** (from [RANK_3X3_SEARCH.md](../../research/small-rank-search/RANK_3X3_SEARCH.md) §5.2):

| problem | raw `{−1,0,1}` | `Z/3`-equiv | canonical | brute-force? |
|---|---|---|---|---|
| `⟨2,2,2⟩, r=7` (Strassen) | `3^{84} ≈ 10^{40}` | `3^{28} ≈ 10^{13}` | `~10^{7.9}` | **yes — seconds** |
| `⟨3,3,3⟩, r=21` | `3^{567} ≈ 10^{270}` | `3^{189} ≈ 10^{90}` | `~10^{73}` | **no — never** |
| `⟨3,3,3⟩, r=22` | `3^{594} ≈ 10^{283}` | `3^{198} ≈ 10^{94}` | `~10^{78}` | **no — never** |
| `⟨3,3,3⟩, r=23` (Laderman) | `3^{621} ≈ 10^{296}` | `3^{207} ≈ 10^{99}` | `~10^{82}` | **no — never** |

The jump from feasible (`10^8`) to "more candidates than atoms in the
observable universe" (`10^{73}+`) happens over a single dimension increment:
`n=2` → `n=3`. Both Z/3-equivariance and canonicalization help by *constant*
log-factors; the dimensional explosion is in the exponent and cannot be cut
by any symmetry argument.

That's why the field's energy moved from enumeration (1969–1990s) to ALS,
SAT, and learned search (2000s–present) the moment the `n=3` frontier opened up.

---

## 3. ALS (Alternating Least Squares)

**What it is.** A converging iterative optimization algorithm for problems
that are non-convex in all variables jointly but **linear in each variable
separately**. Our objective `‖T − Σ_k u_k ⊗ v_k ⊗ w_k‖²` is exactly that
shape: cubic in `(U, V, W)` jointly, but quadratic-and-convex if you fix
any two of `(U, V, W)` and solve for the third.

Relationship to gradient descent: same family (both find local minima of
non-convex objectives) but different mechanic. Gradient descent takes many
small steps along `−∇L`; ALS takes few big steps each of which is the
**exact** minimum of a one-variable sub-problem. When the problem has
block-linear structure, ALS is dramatically faster per iteration. Like
gradient descent, ALS converges to *some* fixed point (provably, by
monotone descent), but that fixed point is only guaranteed to be local —
no global-min guarantee.

**Method.** Hold two of `(U, V, W)` fixed, solve for the third by linear least
squares; alternate. Each sub-problem is convex and has a closed-form solution
via a normal-equation solve. Iterate until the Frobenius residual `‖T − Σ_k
u_k ⊗ v_k ⊗ w_k‖` stops decreasing.

**Convergence.** ALS is a *contraction* near a residual-0 solution: if you
start near one, you snap to it in dozens of iterations (verified in this
repo's `TestAls.alsConvergesNearStrassen` and `alsConvergesNearLaderman`).
Globally, however, it converges only to a **local minimum** of the loss
landscape. The number of distinct local minima is enormous, and most are
"swamps" — plateaus at residual `~10^{-3}` where progress slows to nothing.

**What ALS proves.** A single ALS run that converges to residual 0 is a
*constructive proof* that a rank-`r` decomposition exists. A million ALS runs
that all fail to converge **prove nothing**: there may still be a rank-`r`
solution in a basin that random init never sampled.

**Cost per run.** For `n=3, r=23`: ~few ms per outer iteration on a laptop.
Restarts are embarrassingly parallel.

**Implementation in this repo.** See [`v3/Als.java`](../../src/main/java/io/cormoran/strassen/v3/Als.java).
Z/3-equivariant variant in [`v3/Z3Als.java`](../../src/main/java/io/cormoran/strassen/v3/Z3Als.java),
which restricts the search to decompositions that are themselves Z/3-symmetric
— cuts the parameter dimension from `27r` to `9r` and halves equation count
under cyclic averaging (see [RANK_3X3_SEARCH.md](../../research/small-rank-search/RANK_3X3_SEARCH.md) §4).

**What we know from the literature.** ALS easily finds `r = 23` for ⟨3,3,3⟩
(rediscovers Laderman and Smirnov-style variants). At `r = 22` and `r = 21`
ALS has been run for tens of CPU-years across multiple research groups with
**zero successful runs** (as of 2026). That's strong empirical evidence that
either (a) no Z/3-equivariant solution exists at those ranks, or (b) the
basins are too narrow for random init to find. We can't distinguish (a) from
(b) without a complete method.

---

## 4. SAT/SMT (constraint propagation)

### 4.0 What SAT and SMT are

**SAT = Boolean satisfiability.** Given a logical formula over Boolean
variables, decide whether some assignment makes the formula true. The
canonical NP-complete problem (Cook 1971), worst-case exponential but in
practice extraordinarily fast on structured inputs thanks to **Conflict-
Driven Clause Learning (CDCL)**: when the solver hits a contradiction it
*learns* a new clause that prevents revisiting any branch that would hit
the same contradiction. Industrial solvers (**Kissat**, **CryptoMiniSat**,
**Glucose**, **CaDiCaL**) routinely handle formulas with millions of
variables in seconds. A SAT solver returns either **SAT** (with a
satisfying assignment) or **UNSAT** (with a proof certificate — typically
a resolution refutation that's machine-checkable).

**SMT = Satisfiability Modulo Theories.** SAT extended with richer
non-Boolean variables: integers, reals, bit-vectors, arrays, strings,
etc., parameterized by a "theory" that says what constraints over those
variables mean. The solver runs SAT over the formula's Boolean skeleton,
delegates arithmetic consistency checks to a theory-specific solver, and
ping-pongs between them until both agree or one finds a contradiction.
The reference SMT solver is **Z3** (Microsoft Research, open-source);
others include **cvc5**, **Yices**, **MathSAT**. Most have Java/Python
bindings.

**Why SAT/SMT for matmul.** The 729 cubic equations for `⟨3,3,3⟩` have
no continuous structure to descend along, so ALS and gradient methods
have nothing to grab onto. But the equations *do* have rich discrete
structure (sparsity, value bounds, symmetry), and that's exactly what
modern SAT/SMT solvers exploit. Crucially, **SAT/SMT is the only method
on this list that can return UNSAT** — a certified proof that no
algorithm exists at a given rank within a given alphabet. ALS and RL
can show existence (by constructing an algorithm) but their *failure*
to find one proves nothing.

The catch: this complete-search power is **alphabet-bounded**. SAT
proper is decidable but only over finite alphabets (Booleans); SMT over
reals is in principle decidable but practically intractable at our
problem sizes. So SAT/SMT in this field is almost always SAT over
`{-1, 0, +1}` or similar small finite sets — see §4.1 for why.

### 4.1 The alphabet question

**In practice, yes — SAT/SMT for matmul search always works over a fixed
finite alphabet.** Pretty much every published attack on `⟨3,3,3⟩` of this
flavor — Heule, Sedoglavic, Courtois et al. — fixes the alphabet to
`{-1, 0, +1}` (sometimes extended to `{-1, 0, +1, ±1/2, ±2}`).

The reason is theoretical and practical:

- **SAT proper** is decision over Booleans. Each scalar unknown gets encoded
  as ~`⌈log₂(|alphabet|)⌉` Boolean variables. An alphabet must therefore be
  finite — over `R` there is literally no SAT encoding.
- **SMT over nonlinear real arithmetic (NRA)** *can* in principle handle
  continuous unknowns: Z3 supports `QF_NRA` via cylindrical algebraic
  decomposition (CAD), which is decidable but doubly-exponential in the
  number of variables. For our problem at `n=3, r=22` (198 real unknowns,
  249 cubic constraints under `Z/3` equivariance) this is astronomically
  out of reach — minutes or hours just to load the formula, no published
  termination on the full problem.
- **Gröbner-basis methods** (CoCoA, Macaulay2, msolve, etc.) can compute
  over `Q` exactly and would in principle certify (non-)existence over the
  algebraic closure. They scale even worse than NRA-SMT for this problem.
  The famous Hopcroft–Kerr 1971 lower-bound proof that `R(⟨2,2,2⟩) = 7` over
  `R` used Gröbner-style polynomial-ideal arguments — but that proof is for
  `n=2`, and nothing analogous has been pushed through for `n=3`.

So the honest summary: **continuous-alphabet SAT/SMT exists in theory but is
not a practical tool for this problem size.** Every reported computational
non-existence result in the literature is alphabet-bounded.

### 4.2 The actual method

Restrict to a discrete alphabet (typically `{-1, 0, +1}`), encode each scalar
unknown as ~2 booleans, expand the `n^6` cubic equations into CNF clauses via
Tseitin transformation, and feed to a modern SAT/SMT solver (Z3, cvc5,
CryptoMiniSat over GF(2)). The solver runs unit-propagation + DPLL +
learned-clause heuristics over a `~10^{73}` formal space without ever
materializing it.

**What SAT proves.**
- **SAT** (a model found): constructive existence proof. Same standing as
  ALS — it produces an actual `(U, V, W)`.
- **UNSAT**: a *certified* non-existence proof, **conditional on the
  alphabet**. This is the only method on this list that can prove non-existence,
  full stop.

**Why "conditional on the alphabet" matters.** UNSAT over `{-1, 0, +1}` does
**not** rule out an `r`-decomp over `{-1, 0, +1, ±1/2}`, over `Q`, or over
`R`. The Smirnov 2017 catalog of `r=23` algorithms for `⟨3,3,3⟩` includes
some with entries like `1/3` and `2/3`, which a `{-1, 0, +1}` SAT solver
would (correctly) call UNSAT — alphabet-bounded UNSAT just doesn't say what
people often hope it says.

This is the fundamental epistemic asymmetry of the field:

| claim | what proves it |
|---|---|
| "rank `r` solution exists over `R`" | one ALS run that converges, or one SAT model |
| "no rank-`r` solution exists over `{−1,0,+1}`" | certified UNSAT |
| "no rank-`r` solution exists over `R`" | **only Gröbner / NRA-SMT — infeasible for `n=3`** |

The last row is precisely why nobody has been able to prove `R(⟨3,3,3⟩) > 22`
in 50 years. We have lower bounds of `19` from structural arguments (border
rank, asymptotic-rank techniques) that don't depend on alphabet enumeration —
but they're loose.

### 4.3 Cost

Scales badly with alphabet size and `r·n²`. For `n=3, r=23` over `{-1, 0, +1}`,
terminations in hours-to-days have been reported (Heule, Sedoglavic). At
`r = 22` it has been pushed to weeks of compute without a result, but no
published certified UNSAT either.

**Not yet implemented in this repo.** Would need a JNI-style bridge to Z3 or
similar — substantial engineering.

---

### 4.4 Fractional coefficients and the integer/rational alphabet question

**Do they exist?** Yes. The matmul literature contains rank-`r` decompositions
with entries in `{-1, 0, +1, ±½, ±⅓, ±⅔, ...}`. Per the survey discussion of
Heun-Sedoglavic and Smirnov: "many schemes involve fractional coefficients
and therefore do not apply to arbitrary coefficient rings." Smirnov's 2017
catalog of `⟨3,3,3⟩` `r=23` algorithms in particular includes variants with
denominators 2 and 3.

**Are they equivalent to integer-coefficient algorithms?** Per-term scaling
`(α, β, γ)` with `αβγ = 1` lets you rotate denominators among the three slots
`(U, V, W)` of any rank-1 term — but you cannot make all denominators
*disappear* unless they were absent to start with. Concretely:

| original alphabet | LCM of denoms | cleared-integer form (per-term scaled) |
|---|---|---|
| `{0, ±½, ±1}` | `2` | `{0, ±1, ±2}` lands in *one* slot, other slots integer |
| `{0, ±⅓, ±1}` | `3` | `{0, ±1, ±3}` (one slot) |
| `{0, ±½, ±⅓, ±⅔}` | `6` | `{0, ±2, ±3, ±6}` — `±6` needed in *one* slot |
| `{0, ±⅕}` | `5` | `{0, ±1, ±5}` (one slot) |

So a rational rank-`r` decomposition with combined-denominator `d` is
*equivalent* to an integer rank-`r` decomposition (same rank, computing
`d·C` instead of `C` — the final scaling by `1/d` is free per the bilinear-
complexity model), but the entries are integers up to `d` in magnitude, not
necessarily small.

**Why this matters for search.** A SAT search over `{-1, 0, +1}` only
covers ternary algorithms with no per-term rescaling cleanup. A SAT search
over `{-3, ..., +3}` covers all the denominator-≤3 rational algorithms
*after* per-term clearing, but only when the clearing leaves entries in
`{-3..3}`. The two searches:

- explore differently-shaped prefix trees (different SAT propagations),
- can succeed or fail on different sub-problems,
- and an UNSAT result over one alphabet says *nothing* about the other.

**TODO: pin down concrete references and add them as runnable Java.**
We have not yet hard-coded a fractional-coefficient algorithm in
`src/main/java/io/cormoran/strassen/v3/`. The natural candidates to
transcribe (with verbatim formulas from the original papers, since paraphrase
is error-prone) are:

- Smirnov 2017 — `⟨3,3,3⟩` `r=23` variants with denominator-2 entries.
  Paper: A.V. Smirnov, *"Optimizing Matrix Multiplication Is NP-Hard"* and
  the supplemental catalogs in the Russian-language sequel; explicit
  `(U, V, W)` tables typically appear as paper appendices or supplementary
  GitHub repos.
- Pan 1980 — `⟨2,2,2⟩` trilinear-aggregation derivations with rational
  intermediate coefficients.
- AlphaTensor 2022 — small-format algorithms over `Z/2` (boolean
  multiplication only) and over the standard ring; the standard-ring
  catalogs were published as supplementary data (Nature 2022) but
  individual algorithms have integer entries with magnitudes up to `~5`.

Until these are in the codebase, the v3 verifier (which doesn't care about
alphabet) can verify *any* `(U, V, W)` triple including ones with rational
entries — so the addition is purely a matter of transcription. PRs welcome.

---

## 5. Gradient / second-order methods

**Method.** Treat `‖T − Σ_k u_k ⊗ v_k ⊗ w_k‖²` as an objective and run gradient
descent (Adam, L-BFGS, Levenberg–Marquardt, etc.) instead of alternating
linear solves.

**Trade-off vs ALS.** Often faster per iteration in low-precision regimes
(GPU-friendly). Worse near solutions (loses the closed-form convergence ALS
gets for free). In practice mostly used as a **drop-in alternative to ALS**
with restart-based search, not a fundamentally different strategy. Both find
local minima.

---

## 6. Reinforcement learning (AlphaTensor)

**What it is.** A 2022 result from DeepMind (the AlphaGo / AlphaZero team)
that recast bilinear-algorithm discovery as a single-player game: the agent
places rank-1 terms one by one onto a "board" representing the residual
tensor, and the game ends when the residual is zero. Then standard
deep-RL machinery (AlphaZero-style MCTS + neural policy/value networks)
searches the game tree, learning from self-play episodes which rank-1
terms tend to make progress. The big idea: replace hand-coded heuristics
(symmetry priors, structural priors, restart strategies) with a learned
policy network.

**Result:** found new matmul algorithms for `⟨4,5,5⟩`, `⟨5,5,5⟩` and several
other formats over `Z/2`, and rediscovered Laderman at `r = 23` for `⟨3,3,3⟩`,
but did **not** improve the standard-arithmetic state of the art at `n = 3`.
Reference: Fawzi et al., *Discovering faster matrix multiplication algorithms
with reinforcement learning*, Nature 610, 47–53 (2022).

**Method.** Train a policy network to construct a decomposition one rank-1
term at a time, with reward = `-1` per term and a terminal bonus for hitting
residual 0. Episodes are generated by self-play; symmetries are exploited as
a data-augmentation prior.

**What it proved.** Rediscovered Laderman at `r = 23`. Found new algorithms
over `Z/2` and other small finite fields, but **did not improve standard
arithmetic at `n = 3`**. That null result is itself meaningful evidence that
`r = 22` is genuinely hard.

**Cost.** ~`10^{12}` simulated moves over the full training run — orders of
magnitude more than ALS for the same problem, but searches a structurally
richer space (the network learns problem-specific heuristics ALS doesn't have).

---

## 7. Symmetry exploitation (orthogonal)

This is not a *separate* solver — it composes with every method above.

**Canonicalization (orbit-breaking).** Pick a representative per orbit of the
symmetry group acting on the search space. Search visits each shape once
instead of `|G|` times. Doesn't change the parameter dimension.

**Equivariance (shape restriction).** Demand that the candidate itself is
invariant under a chosen subgroup. For `Z/3` on `⟨3,3,3⟩` decompositions:
parameter dimension drops from `27r` to `9r` (factor of 3 *in dimension*,
exponentially huge in enumeration size). Risk: may exclude all solutions if
none happen to be `Z/3`-equivariant.

See [SYMMETRIES.md](SYMMETRIES.md) for the group structure and [RANK_3X3_SEARCH.md](../../research/small-rank-search/RANK_3X3_SEARCH.md) §4 for
quantitative reduction tables.

The repo's `Z3Als` is the equivariance instance applied to ALS. Layering
canonicalization on top of ALS is also possible (e.g. lex-min the columns of
`U` after each update) but expensive enough to rarely be worth it.

---

## 8. Putting it together

A realistic search pipeline for ⟨3,3,3⟩ at `r = 22`:

1. **Symmetry layer**: search over `Z/3`-equivariant decompositions only.
2. **Continuous search**: run ALS (or gradient) from `~10^5–10^7` random
   initializations. Catalogue all distinct local minima below some residual
   threshold.
3. **Rationalization**: for the lowest-residual candidates, attempt to snap
   numerical entries to small rationals. If snapping preserves the
   equations, you have a candidate exact algorithm.
4. **Discrete search**: in parallel, run SAT/SMT over `{-1, 0, +1}` with the
   equivariance constraint baked in. A certified UNSAT at this alphabet
   rules out a wide subclass of solutions.
5. **Negative-result publication**: even step 4 by itself, completed, would
   be a publishable bound on the structure of any `r = 22` decomposition.

No method on this list is "the answer" alone. The opportunity is in the
combination: ALS provides fast existence checks, SAT provides certified
non-existence (within an alphabet), symmetries cut everything down, and RL
or hand-tuned heuristics improve sample efficiency. Brute force underlies
none of it past `n = 2`.

---

## 9. Summary table

| method | finds solutions? | proves non-existence? | scales to `n=3, r=22`? |
|---|---|---|---|
| Brute force (raw alphabet) | yes (tiny problems only) | yes over alphabet (tiny problems only) | **no** |
| ALS | yes (locally) | no | yes per run, no global guarantee |
| Gradient / Adam | yes (locally) | no | same as ALS |
| SAT over `{−1,0,+1}` | yes | yes — **alphabet-bounded only** | borderline (weeks of compute) |
| SMT-NRA over `R` | yes (in principle) | yes over `R` (in principle) | **no** — astronomical |
| Gröbner basis over `Q` | yes | yes over `R` | **no** — astronomical |
| AlphaTensor-style RL | yes | no | yes (requires GPU farm) |
| Symmetry reduction | n/a (layer) | n/a (layer) | composes with all above |

The bottom line: search at `r = 22` and `r = 21` is currently bottlenecked
not by lack of tools but by **landscape complexity**. A method that
narrows the basin geometry — through better symmetry priors, better
heuristics, or both — is the next breakthrough.

---

## 10. Open follow-up: alternative approaches considered

The methods documented above (ALS, SAT, RL, brute force, symmetry,
gradient) are the established techniques. None has cracked
`R(⟨3,3,3⟩) ≤ 22` in 50 years of attempts. Below is a short menu of
four alternative or hybrid approaches that are **not yet implemented in
this repo** but are concrete enough to test on `⟨2,2,2⟩, r=7` (ground
truth: Strassen) first, then port to `⟨3,3,3⟩` if they pass validation.

Honest framing: none of these is realistically expected to crack the
50-year open problem on a single-developer codebase. They are
**infrastructure investments** — well-engineered tools that future
attempts can build on, and that could surface new `r=23` algorithm
variants or strengthen known negative results.

### 10.1 A1 — ALS + rationalization

**Method.** Run plain `Als` to convergence; for any solution with residual
`< 1e-3`, attempt to snap each entry to the nearest small-denominator
rational (`d ∈ {1, 2, 3, 4, 6}`); re-verify exactness against
`Verifier.residual`. Catches both integer (Strassen-style) and rational
(Smirnov-style) algorithms.

- **2x2 validation**: run ALS at r=7, snap to `{−1, 0, +1}`, verify Strassen.
- **3x3 prospects**: catches Laderman + plausibly new `r=23` variants over
  small rationals — Smirnov's 2017 catalog has entries with denominator 2
  that our current ALS would miss without rationalization.
- **Effort**: ~100 lines of Java. New utility class `Rationalize.java` in
  `v3/`.

### 10.2 B1 — MCTS without learning (poor-man's AlphaTensor)

**Method.** Apply Monte Carlo Tree Search to the AlphaTensor formulation:
tree nodes are partial decompositions, edges are "add rank-1 term", leaves
are complete decompositions. Use UCB1 for selection; rollout policy is
random rank-1 sampling from `{−1, 0, +1}`. Evaluation function is the
Frobenius residual. **No neural network** — uses MCTS as a generic
search heuristic.

- **2x2 validation**: should find Strassen reliably (small search tree;
  even random rollouts from a decent prior should hit it).
- **3x3 prospects**: probably plateaus at `r=23` like AlphaTensor's full
  pipeline did. Could surface `r=23` variants not in published catalogs.
- **Effort**: ~600 lines. Tree node structure, UCB1 selection, rank-1 term
  enumeration, parallel rollout, transposition table.

### 10.3 A2 — SAT solver with ALS warm-start, embeddings, symmetries, cube-and-conquer ⭐

**Method.** A staged SAT/SMT pipeline that combines structural priors with
parallel **cube-and-conquer** orchestration. Phased from simplest to hardest
target, so each phase validates the next.

#### 10.3.0 Phased target strategy

Each phase produces real algorithmic results AND validates the next:

| phase | target tensor | what we expect to find | empirical status |
|---|---|---|---|
| 1 | `⟨2,2,2⟩` over `Z/2`, `r=7` | Strassen-equivalent (Z/2 version) | **✓ done**: SAT in ~3.8s with SAT4J + hand-coded column lex-ordering; UNSAT at r=6 in ~2.6s |
| 1b | `⟨2,2,2⟩` triangular over `Z/2`, `r=4` | 4-mult triangular algorithm | **✓ done**: SAT in ~2ms; UNSAT at r=3 in ~1ms |
| 1.5 | **`⟨3,3,3⟩` diagonal-plus-one over `Z/2`** | restricted-position sanity check at n=3 scale (see §10.3.0.3) | **✓ done**: r=5 SAT and r=4 UNSAT both in milliseconds |
| 1.6 | **`⟨2,3,3⟩` over `Z/2`, `r=14`** | SAT-certified `R_{Z/2}(⟨2,3,3⟩) ≥ 15` (= confirm believed-tight LB), or *new* `r ≤ 14` algorithm | encoder ready (non-cubic refactor landed, `TestSatNonCubic#denseZ2_223_r11_isSat` validates); ready to run via `Phase16Runner` |
| 2 | `⟨3,3,3⟩` over `Z/2`, `r=23` | Laderman-equivalent (Z/2 version) | **✗ blocked**: needs BreakID and/or embedding clauses, see §10.3.0.2 |
| 3 | `⟨2,2,2⟩` over `R`/`{-1,0,+1}`, `r=7` | Strassen | not started |
| 4 | `⟨3,3,3⟩` triangular at various `r` | optimal triangular matmul algorithms | not started |
| 5 | `⟨3,3,3⟩` over `R`/`{-1,0,+1}`, `r=23` | Laderman / Smirnov variants | not started |
| 6 | **`⟨3,3,3⟩` over `R`/`{-1,0,+1}`, `r=22`** | **either constructive SAT or certified UNSAT** | the actual goal |

##### 10.3.0.1 Phase 2 empirical scaling

Calibration run 2026-05-25 on dense `⟨3,3,3⟩` over `Z/2`, encoder default
(hand-coded column lex-ordering, no BreakID, no embeddings):

| solver | r | result | wall-clock |
|---|---|---|---|
| SAT4J 2.3.6 | 30 (heavy over-rank) | (no result) | **>60 min, killed** |
| kissat 4.0.4 | 30 | SAT | **22.7 min** |
| kissat 4.0.4 | 27 (= n³ naive) | SAT | **3.9 min** |
| kissat 4.0.4 | 25 | (no result) | >32 min, killed |
| kissat 4.0.4 | 23 (Laderman) | (not attempted, projected hours) | — |
| kissat 4.0.4 | 22 (research goal) | (not attempted, projected days+) | — |

**Bottom line**: SAT4J cannot reach n=3 scale at all. kissat works at heavy
over-rank but slows dramatically approaching Laderman. The next-step
additions described in the next subsection are **required, not optional**,
to make Phase 2 (Z/2 Laderman reproduction) tractable.

##### 10.3.0.2 What unblocks Phase 2

In rough priority order:

1. **BreakID preprocessing** — `v3/sat/BreakIdBridge.java` is wired and
   ready. Currently blocked by the user's macOS Command Line Tools issue
   (`clang++` can't find `<string>`). Resolution: `xcode-select --install`,
   then build BreakID from `https://bitbucket.org/krr/breakid.git`.
2. **Sub-block embedding constraints** — see §10.3.1 point 2 and §10.7.
   Adds derived clauses for the 27 `⟨2,2,2⟩` sub-block embeddings (and
   optionally for `⟨2,3,3⟩`); should give strong unit-propagation
   acceleration.
3. **Cryptominisat alternative** — `brew install cryptominisat` provides
   a solver natively optimized for XOR clauses (which our Z/2 sum
   constraints are). Worth A/B-testing against kissat on the same CNF.
4. **Cardinality bounds** (sparsity-based) — see §10.3.1 point 4.
5. **ALS warm-start** — freeze stable ALS entries as unit clauses.

Until at least one of (1)–(3) lands, the v3/sat pipeline is effectively
capped at Z/2 `⟨2,2,2⟩` (Phase 1) for *unrestricted* search. **Restricted
problems at n=3 (Phase 1.5) are tractable even without preprocessing.**

##### 10.3.0.3 The "diagonal-plus-one" sanity-check strategy

Between Phase 1 (Z/2 `⟨2,2,2⟩`, easy) and Phase 2 (Z/2 `⟨3,3,3⟩`, hard),
there's a useful family of intermediate validation targets: restrict the
matmul inputs/outputs to a small subset of positions. The pipeline runs
at the full n=3 index space (catches any bugs that only appear at that
scale) but with a tiny search space.

**Construction**: pick a position mask `S ⊆ {0..n²−1}` (interpreted as
flatten indices `i·n + j`). Force `U[i, k] = V[i, k] = W[i, k] = 0` for
`i ∉ S` via unit clauses in the encoder. The matmul tensor is restricted
to those positions; everything else is 0.

**The smallest non-trivial restriction at n=3** is the diagonal plus one
off-diagonal entry: `S = {(0,0), (1,1), (2,2), (i, j)}` for some `(i, j)`
off the diagonal. Concrete worked example for `(0, 2)`:

| restriction | structure | provable rank | encoder vars (Z/2, r=5) | search space |
|---|---|---|---|---|
| pure diagonal | independent products | `n` (trivial) | `3·n·r` | `2^{3·n·r}` |
| diagonal + (0,2) | triangular `⟨2,2,2⟩` sub-block + isolated diag cell | **5** (= 4 + 1 by direct sum) | `3·4·5 = 60` | `2^{60} ≈ 10^{18}` |
| diagonal + (1,0) | same shape, different positions | 5 | 60 | `10^{18}` |
| upper-triangular | `n(n+1)/2` per matrix | open (≤ 11 for n=3) | `3·6·r` | larger |
| dense `⟨3,3,3⟩` | — | 23 / open ≥ 19 | `3·9·r` | `10^{83+}` |

**Why this works as a sanity check**:

- **Catches n=3 index-space bugs**: the encoder allocates `n² × r` U/V/W
  rows and `n⁶ × r` product variables — same structure as full `⟨3,3,3⟩`,
  just with most variables pinned to 0. Any indexing bug appears here.
- **Cheap**: SAT and UNSAT both finish in milliseconds (validated
  empirically: ~370 ms for SAT-at-r=5 + UNSAT-at-r=4 together).
- **Provable lower bound**: the direct-sum argument
  (`triangular ⟨2,2,2⟩ + isolated diagonal cell`) gives a tight,
  derivable-by-hand `R = 5`, so UNSAT at r=4 is a checkable correctness
  certificate.
- **Composes with toggleable encoder options**: every symmetry layer,
  every preprocessing toggle can be validated independently against this
  small-but-real n=3 problem before being trusted at full ⟨3,3,3⟩.

**The general "diagonal-plus-one" recipe** as a sanity-check strategy:
whenever a new encoder feature or solver back-end is added, run it
against this Phase 1.5 target first. If `R = 5` doesn't appear within
seconds, the new component has a bug. This is the analogue of the
`Z3BruteForce2x2` "32 → 1" canonicalization-check from
_memory: `project-strassen-uniqueness-2x2-r7`_, lifted to the SAT pipeline at
n=3 scale.

Implementation: `SatMatmulPipeline.findZ2RestrictedDecomposition(n, r,
target, allowedPositions)` plus the `diagonalPlusOne(n, row, col)`
helper. Tests in `TestSatRestricted3x3`.

##### 10.3.0.4 Embedding constraints — empirical structure of Laderman

Counting "active terms" per sub-block embedding for Laderman r=23
(computed by `v3/SubblockAnalyzer.java`, reported in
`LadermanSubblockReport`):

| family | bound | # embeddings | min active | tight ones | slack range |
|---|---|---|---|---|---|
| `⟨2,2,2⟩` | ≥ 7 | 27 | 7 | 4 | [0, 9] |
| **`⟨2,2,3⟩`** | **≥ 11** | **27** | **11** | **12** | **[0, 5]** |
| `⟨2,3,3⟩` | ≥ 15 | 9 | 16 | 0 | [1, 3] |

**Key finding**: nearly half (12 of 27) of the `⟨2,2,3⟩` embeddings on
Laderman are at the *tight* lower bound. **No `r=22` algorithm with
Laderman-like structure can lose a term without violating at least one
`⟨2,2,3⟩` constraint.** The `⟨2,2,3⟩`-family constraints are the
sharpest structural prior available; the SAT encoder should include all
three families (`⟨2,2,2⟩`, `⟨2,2,3⟩`, `⟨2,3,3⟩`) as independent
necessary conditions — they constrain different aspects of the term
distribution.

##### 10.3.0.5 Phase 1.6 — non-cubic ⟨2,3,3⟩ at r=14

**Status (2026-05-26)**: encoder refactor landed — `Z2CnfEncoder` now
takes separate `dimU, dimV, dimW` (with the cubic `(n, r, target)`
constructor delegating to it). `Verifier.intMatmulTensor(n, m, p)` and
`SatMatmulPipeline.findZ2NonCubicDecomposition(n, m, p, r, target)` are
the new entry points. `TestSatNonCubic#denseZ2_223_r11_isSat` validates
the refactor end-to-end on the smallest non-cubic case (`⟨2,2,3⟩` at
the tight upper bound). The standalone `Phase16Runner` main runs the
research target: `java … Phase16Runner [r]` (default `r=15` for SAT
sanity check; pass `r=14` for the UNSAT research target).

`⟨n_A, n_B, n_C⟩` dimension convention:
- U: `n · m` rows (A's flattened dimensions)
- V: `m · p` rows (B's flattened dimensions)
- W: `n · p` rows (C's flattened dimensions)

Expected runtime: ⟨2,3,3⟩ at r=14 over Z/2 should run in **minutes to
hours with kissat** — much faster than ⟨3,3,3⟩ r=22 because the CNF is
~half the size:

| target | scalar vars | product vars | total vars |
|---|---|---|---|
| `⟨2,3,3⟩` r=14 | `(6+9+6)·14 = 294` | `6·9·6·14 = 4,536` | `~9,300` |
| `⟨3,3,3⟩` r=22 | `27·22 = 594` | `729·22 = 16,038` | `~32,700` |

**Either outcome is publishable** (per [RANK_KNOWLEDGE.md](RANK_KNOWLEDGE.md) §6):
- SAT at r=14 → new upper bound `R_{Z/2}(⟨2,3,3⟩) ≤ 14`.
- UNSAT at r=14 → SAT-certified lower bound `R_{Z/2}(⟨2,3,3⟩) ≥ 15`
  (the believed-tight value, but without a published rigorous proof).

**Why start with Z/2**: encoding is dramatically simpler — single boolean per
scalar (no sign + zero indicators), cubic products become AND-XOR (SAT's
native structure), no cardinality awkwardness over mixed signs. Probably 3-5×
faster pipeline per cube. AlphaTensor (2022) found new `Z/2` algorithms over
several small formats, so `Z/2` is also of independent research interest. The
shared structural code (symmetries, embeddings, cubes, warm-start) transfers
unchanged from Z/2 to ternary; only the variable encoding and equation form
change.

#### 10.3.1 Five sources of search-narrowing power

1. **Equation encoding** — the `n^6` cubic equations of `T_⟨n,n,n⟩`. Over
   Z/2: XOR constraints, natively handled by modern solvers. Over ternary
   `{-1, 0, +1}`: 2 boolean vars per scalar, Tseitin clauses for the cubic
   products, pseudo-Boolean constraints for the signed-sum equations.
2. **Sub-format embedding propagation** — every `⟨2,2,2⟩` (and `⟨2,3,3⟩`,
   `⟨3,2,3⟩`) embedded sub-tensor gives a constraint that any valid
   r-decomposition automatically satisfies. Encoded as additional clauses
   for unit propagation; doesn't reduce solution count but accelerates
   search (see §10.7 conceptual derivation).
3. **Symmetry-breaking constraints** that quotient the `~10^{37}` raw
   symmetry group. Each family is a separate, individually-toggleable
   layer:
   - Per-term sign breaking (`αβγ=1`, ±1 only): factor `4^r`.
   - `S_r` lex-ordering on rank-1 terms: factor `r!`.
   - `Z/3` equivariance baked into the variable layout: factor `3^{18r}`-ish.
   - Monomial gauge breaking (`(S_n ⋉ {±1}^n)^3`): factor `48³` at `n=3`.
   - Transpose `Z/2`: factor 2.
   - **BreakID preprocessor** (Devriendt et al. 2016) for automatic
     symmetry detection of anything hand-coded misses.
4. **Cardinality + structural constraints** per cube (see §10.3.2):
   - Exact non-zero count.
   - `Z/3` structure `(f, g)` choice.
   - "No fully-zero rank-1 term" clause.
   - Lower bound on per-corner non-vanishing-term count from embeddings.
5. **ALS warm-start**: run `Als` until residual plateaus; identify entries
   that are stable across multiple restarts (consistently near small
   rationals); freeze them as unit clauses; let SAT complete or prove UNSAT.

#### 10.3.2 Cube-and-conquer orchestration

Partition the search space into **cubes** (each cube is a conjunction of
literal assignments fixing some structural parameters), then solve each cube
independently in parallel. This is the **cube-and-conquer** technique (Heule,
Kullmann, Wieringa 2011), best known for the Pythagorean triples UNSAT proof
(Heule, Kullmann, Marek 2016).

For matmul, the natural cube variables are:

| stratification axis | values | typical # cubes |
|---|---|---|
| total non-zero count (single value per cube) | `N ∈ [N_min, N_max]` | 50–250 |
| `Z/3` structure `(f, g)` at r=22 | `(1, 7), (4, 6), (7, 5), …, (19, 1)` | 7 |
| max entry magnitude (alphabet width) | `{−1, 0, +1}` vs `{−2..+2}` vs `{0, ±½, ±1}` | 3 |
| corner-embedding non-vanishing count per corner | per-corner ≥7 forced; can pin exact count | 4–15 |

**Cube design**: prefer **single-value strips** (one cube per integer `N`,
not a range) — simpler API, no error-prone range arithmetic, equally
parallel. Per-cube startup cost is microseconds; no amortization needed.

**Cartesian product** of axes gives a few thousand cubes; most are
mutually-exclusive trivially-UNSAT, so the actual working set is a few
hundred non-trivial cubes per phase.

**Why this matters**: each UNSAT cube is **publishable structural progress**
("no rank-22 algorithm for `⟨3,3,3⟩` with `(Z/3 structure = (4,6), card =
243, alphabet = {−1,0,+1})` exists"). A complete cube survey gives a map of
which structural configurations of `r=22` are possible.

#### 10.3.3 Toggleable component API

```java
Sat4jMatmulSolver solver = new Sat4jMatmulSolver(n, r)
    .targetTensor(TargetTensor.DENSE)    // or .TRIANGULAR, .CYCLIC...
    .arithmetic(Arithmetic.Z2)           // or .TERNARY, .RATIONAL_HALF
    .withEquationConstraints()           // always on
    .withSignBreaking(true)
    .withTermOrderingBreaking(true)
    .withZ3EquivarianceLayout(true)
    .withMonomialGaugeBreaking(true)
    .withTransposeBreaking(true)
    .withEmbeddingConstraints(true)      // sub-block + rectangular embeddings
    .withNoFullyZeroTermsClause(true)
    .withCubeCardinalityExact(238)       // pin this cube to N=238
    .withCubeZ3Structure(4, 6)           // pin this cube to (f=4, g=6)
    .withCubeAlphabet(Alphabet.TERNARY)
    .withWarmStartFrom(alsResult);       // optional

Optional<BilinearAlgorithm> solution = solver.solve();
```

A separate orchestrator runs the parallel cube survey:

```java
CubeOrchestrator orchestrator = new CubeOrchestrator(baseConfig);
List<Cube> cubes = orchestrator.generateCubes(/* axis specs */);
Map<Cube, SolveResult> results = orchestrator.solveAllParallel(cubes, parallelism);
StructuralSurveyReport report = orchestrator.summarize(results);
// report contains: which cubes returned SAT (and their solutions),
//                  which returned UNSAT (and how this constrains the r=22 space),
//                  which timed out (and what bounded them).
```

#### 10.3.4 Validation discipline

Every component layer + every cube combination is validated by the rule:

> **For `⟨2,2,2⟩` at `r=7`, the SAT survey must produce exactly 1 canonical
> solution after dedup.** If any toggle or cube spec breaks this, the bug
> is in that toggle.

This is the same discipline that `Z3BruteForce2x2` proved out for the
discrete enumeration; we re-apply it to every SAT layer.

#### 10.3.5 External dependencies

- **SAT4J** (`org.ow2.sat4j:org.ow2.sat4j.pb:2.3.6`) — pure-Java SAT/PB
  solver, no native binaries. Good enough for n=2 and Z/2-flavored n=3
  validation; slow but workable for early ternary at n=3.
- **Kissat** or **CryptoMiniSat5** binary on `$PATH` — for production
  ternary runs at n=3, r=22. We export DIMACS, call out via
  `Runtime.exec`. Optional; SAT4J fallback always available.
- **BreakID** binary on `$PATH` — optional preprocessing.

#### 10.3.6 Triangular matmul as side excursion

Triangular `⟨n,n,n⟩` is a sub-tensor of dense (entries below the diagonal
of A and B are forced to zero), with smaller rank: `R(⟨2,2,2⟩_△) = 4` vs
`R(⟨2,2,2⟩) = 7`.

**Reduction**: dense `AB = A_U B_U + A_U B_L + A_L B_U + A_L B_L`, so
`R(⟨n,n,n⟩) ≤ 4 · R(⟨n,n,n⟩_△)`. This 4× factor is too loose to be useful
as a lower-bound transfer for dense (`R(⟨3,3,3⟩) ≤ 4·11 = 44` is way
weaker than Laderman's 23).

**Why include it anyway**:
- Smaller search space (fewer non-zero target entries) → faster validation
  of the SAT pipeline at n=3.
- Triangular matmul ranks are less studied than dense; we'd likely produce
  *new* algorithmic facts (independent value).
- Provides an intermediate between Z/2 ⟨3,3,3⟩ and dense ⟨3,3,3⟩ in the
  phased plan, with a smaller jump in problem size.

#### 10.3.7 Effort

~700-800 lines for the base SAT pipeline (Z/2 + ternary, all layers).
Cube-and-conquer orchestration: +200 lines. Total ~1000-1500 lines
distributed across many small modules. Each module is individually
testable; the phased target strategy ensures we always have a recent
working pipeline at the previous phase before adding complexity for the
next.

### 10.4 C1 — Border rank with parametric families (Bini 1979 style)

#### What "border rank" actually means

A tensor `T` has **border rank** `r` if there's a parametric family of
**exact-rank-r** tensors `T_ε` that converges to `T` as `ε → 0`:

```
T  =  lim_{ε → 0}  T_ε      where each T_ε = Σ_{k=1..r} u_k(ε) ⊗ v_k(ε) ⊗ w_k(ε)
```

The catch: individual factor entries `u_k(ε), v_k(ε), w_k(ε)` typically
**blow up** as `ε → 0` (entries proportional to `1/ε`). Each `T_ε`
*exactly* equals a rank-`r` decomposition, but the decomposition needs
*unbounded arithmetic precision* to represent in the limit.

**Bini 1979** showed `R̃(⟨2,2,2⟩) ≤ 5`, two less than the actual rank 7.
**This does NOT give us a 5-multiplication algorithm for `2×2` matmul in
fixed-precision arithmetic.** In any practical sense (BLAS, floating-point,
machine-integer arithmetic), Strassen at r=7 remains the best — Bini's
algorithm has entries scaling as `1/ε` and is unusable.

**Why border rank matters anyway**: the matmul **exponent** `ω` (smallest
constant with `O(n^ω)` cost for n×n matmul) depends on border rank, not
exact rank:

```
ω  ≤  3 · log_n( R̃(⟨n,n,n⟩) )
```

When you recursively apply a border-rank-r algorithm for `⟨n,n,n⟩` to
larger matrices, the precision issues "wash out" in the asymptotic limit.
So Bini-style border-rank work contributes to **theoretical** `ω` bounds
but produces no algorithm anyone would use in practice.

The two ranks are different complexity measures:
- **Exact rank** matters for *practical* algorithms (BLAS, GEMM).
- **Border rank** matters for the *theoretical* matmul exponent `ω`.

#### Method

Search for **rank-(r−1) approximate** decompositions parameterized by ε,
where `lim_{ε → 0} = T`. Bini's 1979 paper found `⟨2,2,2⟩` at border rank
5 by this trick.

- **Why this matters for `ω`**: improving border rank of `⟨3,3,3⟩` would
  contribute directly to the asymptotic matmul-complexity frontier.
  Current best: `R̃(⟨3,3,3⟩) ≤ 21` (Schönhage 1981, via aggregating
  ε-parametrized smaller-format border-rank algorithms).
- **2x2 validation**: rediscover Bini's border-rank-5 algorithm. Concrete
  and testable.
- **3x3 prospects**: improving the bound `R̃(⟨3,3,3⟩) ≤ 21` would
  contribute to `ω` work even without an exact-rank improvement.
- **Effort**: ~500 lines. Requires **symbolic ε-arithmetic** (not
  floating-point) so we can take the limit cleanly. Either hand-rolled
  rational-polynomial-in-ε class or `commons-math3` `PolynomialFunction`
  composed with `Fraction`.

### 10.5 B2 — Flip-graph navigation (Kauers et al., ternary meta)

The current state of the art for *finding new* small-format algorithms,
distinct from ALS / SAT / RL in design philosophy. Several flavors of
increasing power:

- **Flip graphs** (Kauers & Moosbauer 2023): each vertex is a valid
  rank-r decomposition `T = Σₖ uₖ ⊗ vₖ ⊗ wₖ`. Edges are *flips* — local
  transformations that modify a few rank-1 terms while preserving the
  total sum and the rank. Walk the graph from a known algorithm to
  enumerate variants. "Merge"-flips occasionally find rank-(r−1)
  decompositions, surfacing genuinely new algorithms.
- **Meta flip graphs** (Kauers & Wood 2025, arXiv:2510.19787): extends
  the flip space with *cross-format* moves — embed `⟨n,m,p⟩` inside
  `⟨n',m',p'⟩`, aggregate sub-schemes, recombine. Massively expands the
  searchable space.
- **Ternary meta flip graphs** (2025, arXiv:2511.20317): restricts the
  meta-flips to those preserving the `{-1, 0, +1}` alphabet. Solves the
  practicality problem of `Z/2`-only flip-graph results, which often
  don't lift to usable integer algorithms.

**Method.** Implement vertex = `BilinearAlgorithm`, edge generators =
the ~5–10 standard flip patterns (per the Kauers-Moosbauer formal
catalog). BFS/A* the neighborhood of a seed algorithm (e.g., Laderman
or AlphaTensor's `factorizations_f2.npz` ⟨3,3,3⟩ scheme), watching for
merge-flips that drop rank.

- **Why it's distinct from A2 (SAT)**: SAT does exhaustive search with
  certified UNSAT. Flip graphs do *local* search from a known starting
  point — no UNSAT proof, but much faster at finding *new* algorithms
  near an existing one.
- **Why it's distinct from A1 (ALS)**: ALS does continuous optimization
  with random restarts. Flip graphs operate over discrete-coefficient
  schemes throughout — no rationalization step needed, ternary
  guaranteed by construction.
- **No symmetry breaking required**: the flips themselves traverse the
  symmetry orbit. The "neighborhood" of a seed naturally covers gauge-
  equivalent and gauge-distinct schemes.

**2x2 validation**: load Strassen, BFS the flip graph at r=7. Should
recover the entire `~10^7`-element gauge orbit of Strassen in seconds.
Confirms flip generators are correct.

**3x3 prospects**: this is precisely the framework that produced the
recent 58-addition rank-23 scheme for `⟨3,3,3⟩` (a known improvement
over Laderman's 62 additions, same rank). Finding `r=22` via flip
navigation from Laderman is the research target; nobody has succeeded
yet but the framework is the only one that's *plausibly close*. If a
flip-graph descent from Laderman reaches r=22, that's the breakthrough.

**Recent concrete results from this lineage**:
- Moosbauer-Poole 2024 — new `⟨5,5,5⟩` algorithms via flip-graph
  exploration of AlphaTensor-discovered schemes.
- Kauers-Wood 2025 — meta-flip cross-format moves; structural
  consequences of Moosbauer-Poole algorithms.
- The "58-addition rank-23 `⟨3,3,3⟩`" scheme (from Perminov's repo
  `dronperminov/FastMatrixMultiplication`, see §11.2) — found via
  flip-graph methodology.

**Effort**: ~800–1200 lines.
- Flip generators: ~300 lines for the standard catalog (per
  Kauers-Moosbauer 2023 §3-4).
- Vertex deduplication via canonical form: ~150 lines (can reuse the
  `Z3BruteForce2x2` canonicalization).
- BFS/A* with merge-flip detection: ~200 lines.
- Meta-flip cross-format support: +~300 lines (optional, second phase).
- Seed loaders (Laderman, AlphaTensor `.npz`): +~100 lines.

**External deps**: none (pure Java). Optionally, the open-source
[github.com/jakobmoosbauer/flips](https://github.com/jakobmoosbauer/flips)
repo provides a C++ reference implementation we could shell out to.

### 10.6 Priority order

Per the discussion that produced this list (2026-05-25 / -26), the five
approaches in order of payoff/effort ratio:

1. **A1 (ALS + rationalization)** — quickest infrastructure win, useful
   regardless of which other direction we pursue.
2. **B1 (MCTS without learning)** — most "genuinely new method" for the
   codebase; testable on 2x2.
3. **A2 (SAT + ALS warm-start)** — **highest research-frontier potential**;
   user-flagged as the priority direction. Phase 1 done; Phase 2 blocked
   on BreakID / embedding clauses / kissat scaling.
4. **B2 (flip-graph navigation)** — current state of the art for *finding
   new* small-format algorithms. Best entry point once we have a seed
   algorithm loader (AlphaTensor `.npz` or Smirnov catalog import).
5. **C1 (Bini-style border rank)** — most academically interesting since
   it touches the `ω` frontier, but the highest entry cost.

Not on the list (deferred / skipped):

- **Simulated annealing** on discrete alphabet — slow, well-studied,
  insufficient novelty.
- **Group algebra / Cohn–Umans constructions** — requires substantial
  group-theory infrastructure before any code; weeks of math reading.
- **SDP/Lasserre relaxations for lower bounds** — high impact long-term
  (a proof of `R(⟨3,3,3⟩) ≥ 22` would be publishable) but worst
  effort/payoff ratio for a first pass.

---

## 11. External algorithm catalogs and data sources

Curated list of online catalogs of fast-matmul algorithms (rank decompositions
for various formats). Useful as **validation ground truth** (compare our
pipeline's discoveries against known algorithms), as **warm-start seeds** for
SAT/ALS searches, and for **structural-pattern study** (analyze known r=23
schemes for sparsity / Z/3 structure / coefficient distributions).

### 11.1 The main aggregator

**[FMM Catalogue](https://fmm.univ-lille.fr)** — the central repository of
fast matrix-multiplication algorithms. Curated reference catalog covering
many formats including ⟨2,2,2⟩, ⟨3,3,3⟩, ⟨2,3,3⟩, ⟨4,4,4⟩, ⟨5,5,5⟩, etc.,
with rank, coefficient field, and provenance for each. **Start here when
looking for a known algorithm to validate against.**

### 11.2 Algorithm-source repositories

| source | URL | content |
|---|---|---|
| **AlphaTensor** (DeepMind 2022) | [github.com/google-deepmind/alphatensor/tree/main/algorithms](https://github.com/google-deepmind/alphatensor/tree/main/algorithms) | `factorizations_r.npz` (standard arithmetic), `factorizations_f2.npz` (Z/2), `nonequivalence/alphatensor_14236_factorizations.npz` (14,236 distinct ⟨4,4,4⟩ algorithms). NumPy format. |
| **AlphaEvolve** (DeepMind) | mathematical_results.ipynb on Google Colab | Evolutionary-search results for matmul and other math problems. |
| **Original Flip Graph** (Moosbauer) | [github.com/jakobmoosbauer/flips/tree/main/solutions](https://github.com/jakobmoosbauer/flips/tree/main/solutions) | Foundational flip-graph methodology and the algorithms it produces. |
| **Adaptive Flip Graph** (Yamato-Arai) | [github.com/Yamato-Arai/adap](https://github.com/Yamato-Arai/adap) | Enhanced flip-graph techniques. |
| **Symmetric Flip Graph** (Moosbauer) | [github.com/jakobmoosbauer/symmetric-flips](https://github.com/jakobmoosbauer/symmetric-flips) | Symmetry-preserving flip-graph variant — directly relevant to Z/3-equivariant search. |
| **Meta Flip Graph** (Kauers et al.) | [github.com/mkauers/matrix-multiplication](https://github.com/mkauers/matrix-multiplication) | Advanced flip-graph methodology by M. Kauers' group. |
| **FMM Add Reduction** (Werekorren) | [github.com/werekorren/fmm_add_reduction/tree/main/algorithms](https://github.com/werekorren/fmm_add_reduction/tree/main/algorithms) | Reducing additive complexity of known schemes (separate optimization axis from rank). |
| **Perminov / FastMatrixMultiplication** | [github.com/dronperminov/FastMatrixMultiplication](https://github.com/dronperminov/FastMatrixMultiplication) (GH username `dronperminov`, author Andrew I. Perminov; papers: arXiv:2603.02398, 2511.20317, 2512.13365, 2512.21980, 2606.02480) | Broad catalog ⟨2,2,2⟩ to ⟨16,16,16⟩, ternary `{-1, 0, +1}` focus, JSON + Maple formats. Includes the recent "58-addition rank-23 scheme" for ⟨3,3,3⟩ (arXiv:2512.21980). Its own `analyzed-schemes--data-sources` section is the index that produced this list. |

### 11.3 Individual algorithms in the academic literature

These don't live in a public repository; the explicit `(U, V, W)` tables are
appendices/tables in specific papers (sometimes in supplementary material).

- **Laderman (1976)** — `⟨3,3,3⟩` at r=23 over `{-1, 0, +1}`. First known
  beat of the naive n³=27. Reproduced in this repo as
  [`v3/Laderman23.java`](../../src/main/java/io/cormoran/strassen/v3/Laderman23.java);
  source: Heun et al. arXiv:1108.2830 §2.4.
- **Smirnov (2013)** — *Bilinear complexity and practical algorithms for matrix
  multiplication*, Computational Mathematics and Mathematical Physics 53(12).
  Contains explicit `⟨3,3,3⟩` r=23 over `{-1, 0, +1}` and several other
  small-format schemes.
- **Smirnov (2017)** — Russian-language follow-up extending the 2013 catalog
  with denominator-2 (and higher) rational-coefficient variants.
- **Pan (1978-1980)** — trilinear aggregation; produces algorithms via
  combining smaller-format schemes.
- **AlphaTensor** (Fawzi et al., Nature 610, 47-53, 2022) — new algorithms
  for `⟨4,4,4⟩`, `⟨4,5,5⟩`, `⟨5,5,5⟩` over `Z/2`; the supplementary data is
  in their GitHub repo (above).

### 11.4 How to use these catalogs in this codebase

For our `v3/sat` pipeline specifically:

1. **Validation**: pull a known algorithm from the AlphaTensor `.npz`,
   convert to a `BilinearAlgorithm`, run `SatMatmulPipeline.verifyZ2(alg,
   target)`. If our encoder/decoder is correct, residual = 0.
2. **Warm-starts**: convert a known algorithm's entries into DIMACS unit
   clauses for the corresponding SAT variables, feed to the solver as a
   partial assignment seed. Particularly useful when looking for *new*
   algorithms near a known one (e.g., enumerate small perturbations).
3. **Structural study**: compute sparsity, Z/3 orbit structure, monomial
   gauge orbit for each known algorithm. Calibrates our cardinality
   constraints and informs which `(f, g)` Z/3 structures actually appear.

A `KnownAlgorithms` utility class to load and verify these would be a useful
addition once we have a JSON or `.npz` importer.
