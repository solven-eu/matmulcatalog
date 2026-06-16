# ω-cartography per composition strategy

A working note for predicting which **asymptotic exponent** `ω(strategy)` a
given composition operator can deliver, given the seed atoms' own ω's.
Per the canonical convention:

```
ω(⟨n,m,p⟩ : r)  =  3 · log(r) / log(n·m·p)
```

(`= 3 · ln(r) / ln(n·m·p)`; base of log irrelevant by ratio).

## Why this matters

Search picks the min-rank candidate at a target shape. But for **planning**
— e.g. asking "if my atoms have ω ≈ 2.79, what's the best ω I can reach
at large cubic shapes?" — we don't need to run the search. The formula
below tells you the FLOOR of the strategy's achievable ω as a function
of the seed atoms.

## Strategy → ω closed forms

### Kronecker

`⟨n₁,m₁,p₁⟩=r₁ ⊗ ⟨n₂,m₂,p₂⟩=r₂  →  ⟨n₁n₂, m₁m₂, p₁p₂⟩ = r₁·r₂`

Implied ω at the product shape:

```
ω(Kron) = 3 · log(r₁·r₂) / log((n₁n₂)·(m₁m₂)·(p₁p₂))
       = 3 · (log r₁ + log r₂) / (log(n₁m₁p₁) + log(n₂m₂p₂))
       = SIZE-WEIGHTED AVERAGE of ω₁, ω₂
```

where the weight on each ω is `log(n_i · m_i · p_i)`. Concretely: if you
Kron a low-ω atom with a high-ω atom, the result's ω is **between** them,
biased toward whichever has the larger product volume.

**Consequence**: to push ω down via Kron alone, you compose your lowest-ω
atom with itself recursively. Pure self-Kronecker gives `ω(Kron^k(A)) = ω(A)`
exactly — Kronecker can't IMPROVE ω, only propagate it.

### Concat (axis-split)

`⟨n,m,p₁⟩=r₁ + ⟨n,m,p₂⟩=r₂  →  ⟨n, m, p₁+p₂⟩ = r₁+r₂`

Implied ω:

```
ω(Concat) = 3 · log(r₁+r₂) / log(n·m·(p₁+p₂))
```

This is NOT a simple weighted average — for r₁ = r₂ = r:

```
ω(Concat) = 3 · log(2r) / log(n·m·2·p) = ω(⟨n,m,p⟩:r) + correction
```

The correction depends on r vs n·m·p. Concat ALWAYS HURTS ω asymptotically
relative to either component (rank doubles for output volume that only
doubles on one axis). It's used for shape-fitting, not ω-improvement.

### Strassen-recombination (block split with outer base)

Outer base `⟨bn, bm, bp⟩ = r`. Target `⟨n, m, p⟩` allocated as
`⟨n, m, p⟩ = bn × ⟨n/bn, m/bm, p/bp⟩`. The r outer products become
sub-products on the inner shape.

Implied ω: depends on the **multiset of sub-product shapes** the outer
base induces. For pure self-recursion (Strassen ⟨2,2,2⟩=7 on
balanced `[k,k]³` blocks), the inner shape is `⟨k,k,k⟩` and:

```
ω(Strassen-recur) = ω(outer base)  exactly
```

Strassen³ at ⟨2k,2k,2k⟩ recursed: ω = log_2(7) = 2.807. Doesn't depend
on k.

For unbalanced allocations, the sub-product multiset is mixed:
```
total rank = Σ R(sub_shape_i)  weighted by the outer base's structure
ω(Strassen-recur, alloc) ≈ ω-of-(arg max over sub-shapes weighted by their volumes)
```

So **Strassen-recombination "follows the multiset"** as you put it — the
output ω is dominated by the sub-shape whose `(R · log V)` term is largest.

### Sedoglavic Prop 1

`⟨u+v,u+v,u+v⟩ ≤ ⟨u,u,u⟩ + 3⟨u,u,v⟩ + 3⟨v,v,u⟩` (u > v)

The implied ω at the cubic target:

```
ω(SedoglavicProp1) ≈ max(ω(⟨u,u,u⟩), ω(⟨u,u,v⟩), ω(⟨v,v,u⟩))
                      with size-weighted adjustment for the (1, 3, 3) multiplicities
```

Concretely at (u,v)=(4,3): the three sub-shapes give ω(⟨4,4,4⟩=48) ≈ 2.793,
ω(⟨4,4,3⟩=38) ≈ 2.659, ω(⟨3,3,4⟩=29) ≈ 2.694. The output:
ω(⟨7,7,7⟩=249) = 3·log(249)/log(343) ≈ 2.835. Slightly worse than each
individual sub-ω because the formula's 1+3+3=7 multiplicity adds rank
faster than the volume grows.

### Pan TA pair-fusion (doubling)

`⟨2k,2k,2k⟩ ≤ ⟨k,k,k⟩ + 3·pair_cost(k,k,k)` where pair_cost = k³+3k².

Implied ω: dominated by `3·k³` term (the cubic dominates the quadratic
correction asymptotically). So:

```
ω(Sedoglavic-doubling) → 3 as k → ∞    (asymptotically NAIVE!)
```

Pair-fusion at the doubling shape is useful for SMALL k only. For large
k, plain Kron `⟨2,2,2⟩=7 ⊗ ⟨k,k,k⟩` is better.

## Cartography summary

| Strategy | ω(output) as a function of ω(atoms) |
| --- | --- |
| Self-Kron `A ⊗ A` | `= ω(A)` (preserves, can't improve) |
| Cross-Kron `A ⊗ B` | weighted average of ω(A), ω(B) |
| Concat | strictly WORSE than max(ω(A), ω(B)) |
| Strassen-recur on balanced cubic | `= ω(outer base)` |
| Strassen-recombine general | dominated by `(R · log V)`-max sub-shape |
| Sedoglavic Prop 1 (u > v) | `≈ max sub-ω` plus multiplicity penalty |
| Sedoglavic doubling (u=v=k) | → 3 as k → ∞ (asymptotically naïve) |
| Pan TA pair-fusion in cubic | small wins; asymptotically naïve |
| Schönhage τ-aggregation (#159) | can improve below ω(atoms) — the key |

## Implication for the seed atom set

To push ω DOWN below the best atom's ω, you need an operator that
DOESN'T trivially preserve or worsen. Concretely:

- Kron preserves at best.
- Concat / Recombination preserve at best.
- Sedoglavic Prop 1 preserves at best (multiplicity overhead).
- **Schönhage τ-aggregation can break the floor** — the `(rA−k)·rB + c·rCaux`
  identity gives an ω that can be strictly less than each sub-shape's ω
  when the structural conditions are met. This is what FMM ⟨6,8,9⟩=296
  achieves: ω(296) < min(ω(20), ω(15), ω(29)).

This is why #159 (general τ-search) is the most consequential of the
pending algorithmic gaps: it's the only operator family that can
improve ω, vs the others which just propagate.

## Open question: cartography UI

Could be a SPA panel: "given atom set `{S, A, AE, …}`, what ω can each
strategy reach at target shape `⟨n,n,n⟩`?" with bars showing per-strategy
floors. Would visualise the "what's possible from the seed set" answer
without running the search.
