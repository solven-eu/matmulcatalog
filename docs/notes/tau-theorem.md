# The τ-theorem — what it is and why it matters

A working note for the matmul-catalog. Aimed at someone who knows tensor rank
+ Strassen-style algorithms but hasn't sat down with Pan 1980 / Schönhage 1981
/ DIS09 yet.

## TL;DR

The **τ-theorem** (Schönhage 1981) is the bridge that lets you compute *several*
independent matrix products *jointly* using **fewer total multiplications than
the sum of their individual ranks**. It generalises Pan's 1980 "two-product
pair-fuse" trick to arbitrary collections of matmul tensors. It's the
foundation of the laser method and every post-Strassen ω improvement
(Coppersmith–Winograd, Stothers, Le Gall, …) and — more relevant to us — it's
what FMM-Lille's `⟨17,17,17⟩=2930` recipe uses at the outer level.

## The setup

We have multiple matmul tensors at possibly different shapes:
`T_⟨n₁,m₁,p₁⟩, T_⟨n₂,m₂,p₂⟩, …, T_⟨n_q,m_q,p_q⟩`. Their **direct sum** is the
tensor we get by considering them as independent computations (disjoint
variable sets):

```
T_disj = T_⟨n₁,m₁,p₁⟩ ⊕ T_⟨n₂,m₂,p₂⟩ ⊕ … ⊕ T_⟨n_q,m_q,p_q⟩
```

Concretely: compute `q` independent matmuls `C_i = A_i · B_i` where each
`(A_i, B_i)` is its own input pair.

The **trivial bound** is `R(T_disj) ≤ Σ R(T_⟨n_i,m_i,p_i⟩)` — just run each
algorithm independently. The interesting question is: can we do better?

## Pan 1980 — the entry point

For just `q = 2` and cyclically-related shapes, Pan showed:

```
R(T_⟨a,b,c⟩ ⊕ T_⟨b,c,a⟩)  ≤  abc + ab + bc + ca
```

vs. the trivial `2 · R(⟨a,b,c⟩)`. At `a = b = c = k` this is `k³ + 3k²` vs.
`2k³` — strictly better when `k ≥ 5` (`k³` grows faster than `3k²`).

This is the "TA pair-fuse" we keep seeing. It's `PairedSubProducts.pairCost`
in our codebase.

## Schönhage 1981 — the generalisation

The τ-theorem (often stated as a partial-matmul-to-full-matmul lifting) says
roughly:

> **If** `R̃(T_⟨a,b,c⟩) ≤ r` for `q` shapes whose dimensions satisfy a degree
> identity (basically Σ n_i m_i p_i terms), **then** the matmul exponent ω
> satisfies a corresponding inequality.

In its rank form, the τ-theorem says you can convert a "good" disjoint-sum
algorithm into a "good" single-large-matmul algorithm via aggregation. The
mechanism is what Pan called **trilinear aggregation**: build a single tensor
that *contains* the disjoint sum as a sub-tensor and exploit the embedding.

For practical small-shape construction (what we care about for ⟨17,17,17⟩),
the τ-theorem isn't itself invoked directly. What we use are its
**constructive corollaries**:

1. **Pair-fuse identities** — Pan 1980 (the `k³ + 3k²` formula).
2. **Triple-fuse identities** — Pan 1982 generalisation. Less well-known
   exactly because (per our own scan) the savings are small or zero for the
   `R(⟨k,k,k⟩)` values our catalog supports.
3. **Disjoint-sum constructions for ⟨n,n,n⟩** at odd n — Drevet–Islam–Schost
   2009 §3 puts these into an explicit closed form. The construction
   decomposes `T_⟨n,n,n⟩` (for n odd) into a sum of cyclically-related
   sub-tensors whose joint rank is what `(n³ + 15n² + 14n − 6) / 3` gives.
4. **Kin-row reductions** — Schwartz–Zwecher 2025 §2.20. A post-processing
   pass that finds linearly-dependent rows in (U, V, W) and drops them.
   This sharpens the τ-theorem-derived bounds further.

## Why this matters for `BlockSplitSearch`

Our outer-decomposition search is **Strassen-recombination only**: pick an
outer base ⟨n₀, m₀, p₀⟩ from a pool, enumerate per-axis block allocations,
recurse on the sub-products.

The τ-theorem says: there's a *different* way to compute a matmul tensor.
Pick a multiset of smaller matmul sub-tensors, find linear-form embeddings
(U_α, V_α, W_α) for each, and the sum of their tensor products equals the
target. **This is NOT a Strassen recombination** — it's not parametrised by
an outer scheme + allocation. It's parametrised by *which sub-tensors to
include* + *how to embed them*.

FMM's `⟨17,17,17⟩=2930`:
```
⟨17,17,17⟩ = 4·⟨8,9,9⟩=430 + ⟨9,9,9⟩=486 + ⟨8,8,8⟩=336 + ⟨8,8,9⟩=388  = 2930
```
is exactly this. **Seven sub-tensors at various shapes**, summing to
`4·430 + 486 + 336 + 388 = 2930`. The four `⟨8,9,9⟩` instances are NOT
related by Strassen's 7-product identity — they're four *independent*
sub-products, each with its own `(U_α, V_α, W_α)` embedding chosen so the
sum equals `T_⟨17,17,17⟩` exactly.

Our search misses this because it only knows the Strassen-recombination
shape. We get 2940 (Strassen `(9,8)³` decomposition); the τ-theorem family
gets ≤ 2930.

## The constructive challenge

Implementing the τ-theorem family in our search is harder than
Strassen-recombination because:

1. **The shape multiset is itself part of the search.** Not just allocations.
   For ⟨17,17,17⟩ we'd need to consider `{⟨a,b,c⟩ × multiplicity}` collections
   subject to cover constraints. The space is much bigger than Strassen's
   allocation space.

2. **Finding the linear embeddings is the hard part.** For each candidate
   shape multiset, you need to *prove* a valid `(U_α, V_α, W_α)` exists and
   *find* it. This is essentially solving a tensor-decomposition feasibility
   problem.

3. **Cover constraints are necessary but not sufficient.** The
   `Σ n_α·m_α ≥ n·m` etc. constraints are easy to check. They DO rule out
   most invalid candidates. But many candidates satisfying them are still
   infeasible.

## Paths the literature offers us

- **DIS09 §3 odd-n closed form**: gives a constructive recipe for any odd n,
  with rank `(n³ + 15n² + 14n − 6)/3`. We have this as
  `PanTrilinearAggregation.build(int n)` in the codebase (Islam 2009 form).
  At n=17 it gives 3160 — much worse than 2930. So Islam's odd-n recipe
  is the foundation but isn't the FMM recipe.

- **DIS09 + Schwartz-Zwecher 2025 kin-row reduction**: sharpens
  Islam's bound to something tighter. For n=44 (even) SZ get 36110 vs Pan
  1982's 36133. For n=17 (odd) there's no published improvement formula —
  one would need to manually apply kin-row reduction to Islam's construction.

- **FMM-Lille's specific 2930/2934 recipes**: these are explicit
  hand-engineered constructions that don't appear in any closed-form formula
  we know. They might be the result of a search procedure FMM-Lille runs
  internally, not a formula evaluation.

## What we'd implement first

The most productive concrete step for task #88 (derivation side) is to:

1. **Verify our existing `PanTrilinearAggregation.build(17)`** actually
   produces a correct algorithm and matches the rank-3160 Islam bound.
2. **Add kin-row reduction as a post-pass** (task #103 generalised). This
   alone might reduce 3160 to something competitive — though almost certainly
   not all the way to 2930.
3. **Implement DIS09 §3 odd-n directly** if it's not what
   `PanTrilinearAggregation.build` currently does — that's the canonical
   constructive recipe.
4. **Long-term**: implement a search over (shape-multiset, linear-embedding)
   pairs with a feasibility check (task #102).

## References

- Pan 1980: "Strassen's algorithm is not optimal" (BIT 18)
- Schönhage 1981: "Partial and total matrix multiplication" (SIAM J. Comput. 10)
- Drevet–Islam–Schost 2009 (DIS09): `references/papers/drevet_islam_schost_2009_DrIsSc09.pdf`
- Pan 2014: review article, `references/papers/pan_2014_trilinear_apa_arxiv1412.1145.pdf`
- Schwartz–Zwecher 2025: `references/papers/schwartz_zwecher_2025_feasible_matmul_arxiv2508.01748.pdf`
