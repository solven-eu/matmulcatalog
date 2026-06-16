# Lower bounds on bilinear rank `R(⟨n,m,p⟩)`

Companion to the machine-readable `docs/lower-bounds.json`. This file
explains technique families, attribution chains, and surveys the state
of known lower-bound results for small matrix multiplication formats.

A bilinear rank lower bound `R(⟨n,m,p⟩) ≥ ℓ` says **no algorithm
using bilinear products can compute the `⟨n,m,p⟩` matrix product in
fewer than `ℓ` scalar multiplications, regardless of how additions are
arranged.** Field discipline matters: a bound over `R` is also a bound
over `Q ⊂ R`, but the converse needn't hold — `F₂`-specific bounds
(Wang 2026) can be strictly tighter than the field-agnostic ones.

## Companion files

| File | Purpose |
|---|---|
| `docs/lower-bounds.json` | Machine-readable LB table consumed by the SPA |
| `docs/border-ranks.json` | Border-rank (R̃) bounds, separate axis |
| `docs/BORDER_RANKS.md` | Human-readable companion for border ranks |

## Technique families

### 1. Substitution / pinning
The oldest technique: fix some input variables to specific values,
reducing the matmul to a sub-problem with known rank. Strassen 1973
and Hopcroft-Kerr 1971 use substitution to prove `R(⟨2,2,2⟩) ≥ 7`.

### 2. Strassen's number
Strassen 1983 introduced an algebraic invariant (the "tensor's rank
over the closure of the field extension") that gives a lower bound on
border rank, which in turn lower-bounds rank.

### 3. Degeneration arguments
A tensor `T₁` *degenerates* to `T₂` if `T₂` is in the orbit closure of
`T₁`. Border rank is monotone under degeneration, so if `T₂` has
border rank `r`, then `T₁` does too. Used by Landsberg, Michalek to
prove tight border-rank bounds.

### 4. Bläser's substitution-method-with-symmetry
Bläser 2003 extended substitution to exploit the S₃ symmetry of the
matmul tensor (cyclic and transpose permutations), achieving
`R(⟨3,3,3⟩) ≥ 19` — the current best general-field bound.

### 5. Certifying dynamic programming (Wang 2026)
Wang 2026 used field-specific dynamic-programming arguments over `F₂`
to push the LB to `R_{F₂}(⟨3,3,3⟩) ≥ 20`. Highlights the field
sensitivity of rank.

## Key results by format

### `⟨2,2,2⟩ = 7` (tight)

| field | LB | UB | gap |
|---|--:|--:|--:|
| all | 7 | 7 (Strassen 1969) | 0 |

**Bound**: Hopcroft-Kerr 1971 (Cornell TR 69-44, Sep 1969 — same month
as Strassen's UB) + Winograd 1971 (independent proof, same paper that
introduces the 15-additions Winograd-Strassen variant).

### `⟨3,3,3⟩` — open

| field | LB | UB | gap |
|---|--:|--:|--:|
| `R, Q, C` | **19** (Bläser 2003) | 23 (Laderman 1976) | 4 |
| `F₂` | **20** (Wang 2026) | 23 (Laderman, mod 2) | 3 |

**Field discipline note**: AlphaTensor 2022 did *not* improve
`R(⟨3,3,3⟩)` over `F₂`. The F2-LB-tightening from 19 to 20 by Wang
2026 closes the gap from 4 to 3 in that field specifically. The
20-vs-19 distinction is exactly the one you'd have to look up in this
document.

### `⟨4,4,4⟩` — open, field-dependent UB

| field | LB | UB | gap |
|---|--:|--:|--:|
| `F₂` | 33 (Smirnov 2017) | **47** (AlphaTensor 2022) | 14 |
| `C` | 33 | **48** (AlphaEvolve 2025, DPS 2025) | 15 |
| `R, Q` | 33 | 48 (DPS 2025) | 15 |

The famous "different best UB per field" canonical case. F₂ vs R
differ by 2 because mod-2 reduction has structural slack.

### Other notable small formats

- `R(⟨2,2,3⟩) = 11` (HK 1971, tight, all fields)
- `R(⟨2,3,3⟩) = 15` (HK 1971, tight, all fields)
- `R(⟨2,4,4⟩) ≥ 26` (Hopcroft-Kerr), UB 26 → tight
- `R(⟨2,2,n⟩) = ⌈(3·2n + max(2,n))/2⌉` = `⌈3.5n + 0.5⌉` for n ≥ 2

See JSON for the full table.

## Attribution chains for compositional results

When a bound is the result of a non-trivial argument combining several
earlier results, the `source` field in `lower-bounds.json` should
record the **paper that did the analysis**, with the leaf-result chain
recoverable from the `notes` field. Example: Bläser's `R(⟨3,3,3⟩) ≥ 19`
combines Strassen-1973-style substitution with S₃-symmetry observations
that go back to de Groote. The bound is attributed to Bläser 2003
(the analyzer); the chain is mentioned in the notes for archival.

## Conventions

- `tight: true` means we have a matching upper-bound algorithm in our
  catalog (a `lookup.find(...)` returns a scheme at exactly `lb`
  rank). When new UBs appear that match the LB, set `tight: true`.
- Permutation: `format` is canonical-sorted (rank is `S₃`-invariant
  via tensor symmetry). `⟨3,2,3⟩` and `⟨3,3,2⟩` collapse to `[2,3,3]`.
- Field-agnostic bounds use `field: "all"`; field-specific use the
  specific tag (`F2`, `R`, `C`).
