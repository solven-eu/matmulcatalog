# Lineage syntax

Every composed (non-atom) scheme in the catalog carries a **lineage**: a
structured record of *how it was built from leaves*. The lineage is what
makes a derived bound reproducible — a `LineageReplayer` can re-materialise
the factor matrices `(U, V, W)` from the lineage plus the on-disk atoms,
without storing the (often large) matrices themselves. This is the basis of
the two-tier catalog (atoms + lineage stubs).

This document is the canonical reference for the lineage grammar. It is
written to be liftable verbatim into the paper as an appendix; the code is
the source of truth (`eu.solven.matmul.catalog.Lineage`,
`eu.solven.matmul.catalog.Compose`).

## Conventions

A matmul format is written `⟨n,m,p⟩`: it computes `C[n×p] = A[n×m] · B[m×p]`.
So `n` indexes output rows, `p` indexes output columns, and `m` is the shared
*inner* (contracted) dimension. Rank (multiplication count) is the `=r` /
`:m` suffix. Always tag the field in prose (`R⟨3,3,6⟩=40`); the field is not
part of the lineage grammar (it is inferred from the leaves).

## Three serialisations

A lineage is one DAG, emitted in three complementary forms. Each scheme JSON
carries all three:

| JSON field | form | audience | dedup |
| --- | --- | --- | --- |
| `lineage_str` | function-call pretty-print | humans (tooltips, review) | plain-repeat |
| `lineage_compact` | one-line infix shorthand | catalog manifest, at-a-glance | plain-repeat |
| `lineage` | JSON DAG of `{"op": …}` nodes | machines (replay) | `@ref` back-refs |

"Plain-repeat" means a shared subtree is written out in full at each
occurrence. The `lineage` JSON instead tags the first occurrence with
`"id":"L<i>"` and emits `{"op":"@ref","id":"L<i>"}` at later occurrences, so
the DAG stays compact and shared work is visible.

## Operation vocabulary

Each op mirrors one primitive in the composition code, so lineage and code
stay in lockstep. Adding an op means extending the `Lineage.Node` sealed
hierarchy *and* the renderers — never repurpose an existing op.

| op | pretty form | compact | source primitive | meaning |
| --- | --- | --- | --- | --- |
| `Atom` | `ref` | `ref` | on-disk leaf | reference to a primitive scheme by filename stem |
| `KronProduct` | `KronProduct(outer, inner)` | `outer ⊗ inner` | `Compose.kroneckerGeneral` | tensor product; rank `r₁·r₂`, shape `⟨n₁n₂, m₁m₂, p₁p₂⟩` |
| `KronChain` | `KronChain(L₁, …, Lₖ)` | `L₁ ⊗ … ⊗ Lₖ` | folded k-way Kronecker | left-fold of `KronProduct` |
| `ConcatCols` | `ConcatCols(left, right)` | `left +p right` | `Compose.concatRight` | **p-axis tile**: `⟨n,m,p₁⟩+⟨n,m,p₂⟩→⟨n,m,p₁+p₂⟩`, shares A, `C=[C₁\|C₂]` |
| `ConcatRows` | `ConcatRows(top, bottom)` | `top +n bottom` | `Compose.concatBelow` | **n-axis tile**: `⟨n₁,m,p⟩+⟨n₂,m,p⟩→⟨n₁+n₂,m,p⟩`, shares B, `C=[C₁;C₂]` |
| `SumInner` | `SumInner(left, right)` | `left +m right` | `Compose.concatInner` | **m-axis sum**: `⟨n,m₁,p⟩+⟨n,m₂,p⟩→⟨n,m₁+m₂,p⟩`, shares neither, `C=C₁+C₂` |
| `Recombination` | `Recombination(base=…, allocA=…, allocB=…, allocC=…, leaves=[…])` | `R[base; aA \| aB \| aC]` | `Recombination.constructWithAllocation` | Strassen-style block recombination of a `base` scheme over per-axis block allocations |
| `RecombinationWithPair` | `RecombinationWithPair(base=…, pairs=…, solo=…, leaves=[…])` | `R*[base; …; pairs/solo]` | `RecombinationWithPair` | recombination with pair-fused sub-products |
| `AugmentSquareDiscard` | `AugmentSquareDiscard(p, n, square)` | `AS(p, n, square)` | `HopcroftKerr2bcAsymmetric.buildNaive` | embed a square scheme into an asymmetric shape, discard padding |
| `Transpose` | `Transpose(child, perm)` | `child^perm` | tensor-symmetry rewrite | axis relabel, e.g. `perm = "NMP->NPM"` |
| `AxisFlip` | `AxisFlip[approximate](../../child, mask)` | `child^J<mask>` | `SymmetryTransforms.axisFlipOrbit` | reverse index order on axes per 3-bit `mask` (bit0=A, bit1=B, bit2=C) |
| `AxisPermute` | `AxisPermute[approximate](../../child, permA, permB, permC)` | `child^π(…)` | `SymmetryTransforms.permutationOrbit` | per-axis row permutations of `child` |
| `DCE` | `DCE(child)` | `child` (elided) | dead-code-elimination pass | drop zero products; semantics-preserving |
| `DisjointSum` | `DisjointSum(c₁ + … + cₖ; TA-legs=…)` | `c₁ ⊕ … ⊕ cₖ` | Pan/Schönhage τ-theorem | sum of smaller matmul sub-tensors, optional trilinear-aggregation legs (rank-prediction marker; materialisation is a follow-up) |

### The three additive siblings (`+p` / `+n` / `+m`)

The matmul tensor `⟨n,m,p⟩` is a direct sum along **any** of its three modes,
each giving a rank-`r₁+r₂` composition. They differ in *which* operand is
reused and whether the output is tiled or summed:

| op | split axis | shares | output `C` | tile or sum |
| --- | --- | --- | --- | --- |
| `ConcatCols` (`+p`) | output columns | A | `[C₁ \| C₂]` | tile (horizontal) |
| `ConcatRows` (`+n`) | output rows | B | `[C₁ ; C₂]` | tile (vertical) |
| `SumInner` (`+m`) | inner / contraction | neither | `C₁ + C₂` | **sum** (accumulate) |

`SumInner` is the odd one out: it splits the shared inner dimension
(`A=[A₁\|A₂]`, `B=[B₁;B₂]`) so `C = A₁·B₁ + A₂·B₂` — the two sub-products
*accumulate into the same output* rather than tiling disjoint regions. It is
a genuine rank-`r₁+r₂` construction nonetheless (a direct sum along the
middle mode). Hence the name `SumInner` rather than `Concat*`.

## Compact-form grammar

```
Atom(ref)                          → ref
KronProduct(A, B)                  → A ⊗ B
KronChain(L₁,…,Lₖ)                 → L₁ ⊗ L₂ ⊗ … ⊗ Lₖ
ConcatCols(L, R)                   → L +p R
ConcatRows(T, B)                   → T +n B
SumInner(L, R)                     → L +m R
Recombination(B, aA, aB, aC, …)    → R[B; aA | aB | aC]      (leaves dropped)
RecombinationWithPair(…)           → R*[B; aA | aB | aC; pairs/solo]
Transpose(child, perm)             → child^perm
AxisFlip(child, mask)              → child^J<mask>
AxisPermute(…)                     → child^π(permA,permB,permC)
DisjointSum([c₁,…,cₖ], …)          → c₁ ⊕ c₂ ⊕ … ⊕ cₖ
DCE(child)                         → child                    (elided)
```

**Recombination drops its leaves in the compact form** — they are recoverable
by re-running the search with the recorded `(base, allocations)` tuple against
the catalog at replay time (≈3× disk savings vs. inlining them). The full
`lineage` JSON keeps the leaves for direct replay.

Base references are kept as plain strings (no single-letter codes): "Makarov"
alone is ambiguous between the 3×3×3 and 5×5×5 results, so the full scheme
reference stays in the lineage.

## `lineage` JSON DAG

Each node is `{"op": "<Op>", …fields…}`. Child fields are named per op
(`outer`/`inner`, `left`/`right`, `top`/`bottom`, `base`/`leaves`, `child`,
`children`). Shared subtrees: the first occurrence gets `"id":"L<i>"`; later
occurrences become `{"op":"@ref","id":"L<i>"}`.

```json
{"op":"ConcatCols",
 "left":{"op":"Atom","ref":"fmm_lille-3x3x6_m40_a862"},
 "right":{"op":"Atom","ref":"fmm_lille-3x3x12_m80_a1724"}}
```

### Legacy op aliases (read-compat)

Before 2026-06 the tile ops were named `ConcatRight` (p-axis) and
`ConcatBelow` (n-axis). The parser still **reads** those op strings (and the
`left`/`right`, `top`/`bottom` field names are unchanged), mapping them to
`ConcatCols` / `ConcatRows`. New materialisations **write** the current
names. On-disk files emitted before the rename therefore still load; a
catalog-wide backfill of the stored strings is an optional follow-up, not a
correctness requirement.

## Worked examples

**`R⟨3,3,18⟩=120` — p-axis tile of two FMM atoms.** `18 = 6 + 12`, so the
same A multiplies both column-blocks of B; outputs concatenated horizontally.
```
lineage_compact: fmm_lille-3x3x6_m40 +p fmm_lille-3x3x12_m80
lineage_str:     ConcatCols(fmm_lille-3x3x6_m40_a862, fmm_lille-3x3x12_m80_a1724)
```

**`R⟨3,3,18⟩=120` — Kronecker, alternative recipe.** `⟨3,3,18⟩ = ⟨1,1,3⟩ ⊗
⟨3,3,6⟩`, rank `3 · 40 = 120` (note the unit factor `⟨1,1,3⟩`, naive rank 3).
```
lineage_compact: naive-1x1x3 ⊗ fmm_lille-3x3x6_m40
```

**Strassen² `⟨4,4,4⟩=49` — Kronecker square.**
```
lineage_compact: strassen-2x2x2_m7 ⊗ strassen-2x2x2_m7
```

**Recombination at an unbalanced cubic.** A `⟨2,2,2⟩=7` base reused over a
`(9,8)`-style block allocation:
```
lineage_compact: R[Strassen<2,2,2>=7 :: CANONICAL; 1,1 | 8,9 | 1,1]
```

**`SumInner` (m-axis sum).** Building `⟨n,m₁+m₂,p⟩` by accumulating two
inner-split sub-products:
```
lineage_compact: scheme-Axm1xP +m scheme-Axm2xP
lineage_str:     SumInner(scheme-Axm1xP, scheme-Axm2xP)
```

## See also

- `eu.solven.matmul.catalog.Lineage` — node hierarchy + the three renderers.
- `eu.solven.matmul.catalog.Compose` — `kroneckerGeneral`, `concatRight`,
  `concatBelow`, `concatInner`.
- `eu.solven.matmul.search.LineageReplayer` — re-materialise `(U,V,W)` from a
  lineage + atoms.
- `eu.solven.matmul.search.ConcatSplitSearch` — discover the best single-axis
  additive split (all three modes).
