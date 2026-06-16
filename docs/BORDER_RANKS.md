# Border ranks `R̃(⟨n,m,p⟩)`

Companion to `docs/border-ranks.json`. The **border rank** `R̃(T)` of
a tensor `T` is the smallest `r` such that `T` is the **limit** of a
sequence of tensors of rank `r`, where the limit is taken in an
algebraic-geometric (ε → 0) sense:

```
R̃(T) = min { r : ∃ T_ε of rank ≤ r with lim_{ε → 0} T_ε = T }
```

Equivalently, `R̃` is the rank of `T` when one is allowed to use
*approximate* bilinear products that recover `T` exactly only in the
limit. Bini-Capalbo-Lotti 1979 introduced this notion precisely
because it can be **strictly smaller than the exact rank**.

## Why care about R̃ vs R

- **Asymptotic complexity**: `ω ≤ log_n R̃(⟨n,n,n⟩)`. Strassen's
  laser method (1987) lifts border-rank bounds into actual fast-matmul
  algorithms via *aggregation* — paying a polynomial overhead in
  exchange for the ε-limit cost.
- **Lower bounds**: `R̃ ≤ R` always, so a border-rank LB is a rank LB.
  The substitution / degeneration arguments (Strassen 1983, Landsberg)
  primarily bound `R̃`.

## Companion files

| File | Purpose |
|---|---|
| `docs/border-ranks.json` | Machine-readable R̃ table |
| `docs/lower-bounds.json` | Exact-rank lower bounds (related but distinct) |
| `docs/LOWER_BOUNDS.md` | Human-readable companion for exact rank |

## Why the SPA toggles between R and R̃

Border rank and rank are **different things**, and showing them mixed
on the same row would be misleading (e.g. `R(⟨2,2,3⟩) = 11` but
`R̃(⟨2,2,3⟩) ≤ 10` — both are correct, neither subsumes the other for
practical purposes). The catalog SPA renders R̃ entries in a separate
view toggle so the comparison is explicit.

## Key results

### `R̃(⟨2,2,2⟩) = 7` (tight, no gap to R)

The classical case where border rank and rank coincide. Landsberg
2006 proved the matching lower bound `R̃ ≥ 7`. There is no Bini-style
ε-decomposition for the 2×2 case.

### `R̃(⟨2,2,3⟩) ≤ 10` — first border-rank improvement

| | rank | border rank |
|---|--:|--:|
| `⟨2,2,3⟩` | 11 (Hopcroft-Kerr 1971) | **10** (Bini-Capalbo-Lotti 1979) |

This is the seminal example. The Bini ε-decomposition uses 10 bilinear
products `(linA + ε · linA')(linB + ε · linB')` that recover the matmul
output up to `O(ε)` — exactly recovered in the limit. This trick is the
basis for the entire laser-method line of work that drives ω.

### `R̃(⟨3,3,3⟩) ≤ 19`

| | rank LB | rank UB | border rank UB |
|---|--:|--:|--:|
| `⟨3,3,3⟩` | 19 (Bläser 03) | 23 (Laderman 76) | **19** (Smirnov 17) |

Strikingly, the border-rank UB **matches the rank LB**. So either
`R(⟨3,3,3⟩) = 19` (and there's a still-unknown 19-mult algorithm), or
`R(⟨3,3,3⟩) > R̃(⟨3,3,3⟩)` (and the gap between exact and border rank
is at least 4).

## On R̃ lower bounds (sketch)

The same techniques used for `R` LB (substitution, degeneration) often
yield `R̃` LB directly because they're *invariant under limits*. The
two lower bounds typically match for small formats and diverge only
at large `n` (where laser-method aggregation pushes ω down via border
rank that doesn't translate to exact rank).

## TODO (catalog stretch)

- Add Landsberg-Michalek 2016 border-rank LBs for `⟨n,n,n⟩` family.
- Pan TA aggregations also give R̃ improvements at moderate sizes —
  cross-reference with the R̃ UB column from `references/papers/pan_2014_trilinear_apa_arxiv1412.1145.pdf` §6.
- Wire up `border-ranks.json` to the SPA via a toggle (right now
  it's emitted but not rendered).
