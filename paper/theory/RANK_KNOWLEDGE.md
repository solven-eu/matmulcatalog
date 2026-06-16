# Matrix Multiplication Rank Knowledge

A reference summary of known ranks, border ranks, addition counts, and
asymptotic complexity bounds for matrix multiplication tensors. Companion
to `RANK_3X3_SEARCH.md` (which focuses on the open `⟨3,3,3⟩` problem) and
`SOLVING_STRATEGIES.md` §11 (external catalogs and data sources).

Complements
[Perminov's FastMatrixMultiplication §research-findings--status][perminov],
which provides a similar table focused on ternary `{-1, 0, +1}` schemes and
their additive complexities — see that repo for ~200 format-specific
entries with provenance.

## Notation

- **`⟨n, m, p⟩`** — matrix multiplication tensor: compute `C = A·B` where
  `A` is `n×m`, `B` is `m×p`, `C` is `n×p`. By symmetry up to relabeling,
  `R(⟨n,m,p⟩) = R(⟨m,n,p⟩) = R(⟨n,p,m⟩)` etc.
- **`R(T)`** — exact tensor rank: minimum `r` such that
  `T = Σ_{k=1..r} u_k ⊗ v_k ⊗ w_k`. **Default convention in matmul
  literature: non-commutative rank (see below).**
- **`R̃(T)`** — border rank: minimum `r` such that `T` lies in the Zariski
  closure of rank-`r` tensors (informally: limit of rank-`r` decompositions
  with entries that may diverge).
- **`R_{Z/2}(T)`** — rank over the field of two elements (boolean
  arithmetic, addition = XOR).
- **Field notation synonyms**: `Z/2` ≡ `Z/2Z` ≡ `F_2` ≡ `F₂` ≡ `GF(2)` — all
  the same field of two elements. Likewise `Z/p` ≡ `F_p` ≡ `GF(p)` for
  any prime `p`. This document uses `Z/2` (matching our codebase
  notation); papers cited often use `F_2` (matching the algebra/SAT-LB
  literature).
- **`ω`** — matrix multiplication exponent: smallest constant with
  `n × n` matmul having an `O(n^ω)` algorithm. Tied to rank by
  `ω ≤ 3·log_n(R(⟨n,n,n⟩))` via recursion, and to border rank by
  `ω ≤ 3·log_n(R̃(⟨n,n,n⟩))` asymptotically.

### Non-commutative vs commutative rank

A bilinear matmul algorithm computes `r` products of the form `M_k =
(linear in A entries) · (linear in B entries)`, then linearly combines
the `M_k` to produce `C` entries.

- **Non-commutative rank** `R(T)`: the minimum `r` **without assuming
  `a·b = b·a`** for the input scalars. This is the version used in
  practical matmul algorithms — it's what's required for the algorithm
  to work when applied **recursively to matrix entries** (since matrices
  don't commute). Strassen's `r=7` for `⟨2,2,2⟩`, Laderman's `r=23` for
  `⟨3,3,3⟩`, AlphaTensor results — all are non-commutative ranks.
- **Commutative rank** `R_c(T)`: the minimum `r` **assuming `a·b = b·a`**.
  Lower than the non-commutative rank for SOME formats (e.g.,
  `R_c(⟨3,3,3⟩) ≤ 21` per Rosowski 2019 vs `R(⟨3,3,3⟩) = 23` Laderman 1976),
  but **doesn't translate to fast matmul** because the recursion on matrix
  entries doesn't preserve commutativity. Note: for `⟨2,2,2⟩`
  commutativity gives NO advantage — Winograd 1971 proved
  `R_c(⟨2,2,2⟩) ≥ 7`, matching Strassen's UB, so `R_c(⟨2,2,2⟩) = 7`
  exactly.

**Default convention in this document, the matmul literature, and this
codebase: non-commutative rank.** Where context matters, commutative
ranks are flagged explicitly (see §1.3 below).

---

## 1. Small format ranks

Bounds are **over `R` (or `Q`)** unless otherwise noted.

### 1.1 Cubic formats

| `⟨n,n,n⟩` | naive `n³` | rank `R` | border rank `R̃` | status | references |
|---|---|---|---|---|---|
| `⟨1,1,1⟩` | 1 | 1 | 1 | trivially tight | — |
| **`⟨2,2,2⟩`** | 8 | **7** | **≤ 5** | rank tight; border-rank UB 5 by Bini | [Strassen 1969][s69], [Hopcroft-Kerr 1971][hk71], [Bini 1979][bini79] |
| **`⟨3,3,3⟩`** | 27 | **`19 ≤ R ≤ 23`** (over `R`/`Q`); **`20 ≤ R ≤ 23`** (over `F₂`) | **`R̃ ≤ 21`** | gap open since 1976; `F₂` LB tightened 2026 | [Laderman 1976][lad76], [Smirnov 2013][smir13], [Schönhage 1981][sch81], [Bläser 2003][bla03], [Wang 2026][wang26] |
| `⟨4,4,4⟩` | 64 | `R ≤ 47` (`Z/2`) `≤ 49` (`R`) | `R̃ ≤ 38` (approx) | open | [AlphaTensor 2022][alpha22], Strassen²=49 |
| `⟨5,5,5⟩` | 125 | `R ≤ 96` (`Z/2`) `≤ 100` (`R`) | open | open | AlphaTensor 2022, Moosbauer-Poole 2024 |

### 1.2 Commutative ranks (for reference — NOT applicable to practical matmul)

Listed for completeness only; these don't help build faster `n × n` matrix
multiplication algorithms.

| `⟨n,n,n⟩` | non-commutative `R` (used in practice) | commutative `R_c` | reference |
|---|---|---|---|
| `⟨2,2,2⟩` | **7** | **7** (= NC; Winograd 1971 proved LB ≥ 7 even with commutativity) | [Strassen 1969][s69] UB; [Winograd 1971][win71] LB |
| `⟨3,3,3⟩` | `[19, 23]` | **`R_c ≤ 21`** (Rosowski closed form: `n(lm+l+m−1)/2 = 3·14/2`); Makarov333 in DIS09 has 22 | [Rosowski 2019/2020][rosowski19], [Drevet–Islam–Schost 2009][drisc09] |
| `⟨4,4,4⟩` | `49` (R) | 46 (Waksman, DIS09); Rosowski's even-n formula also gives 46 | [drisc09], [rosowski19] |
| `⟨5,5,5⟩` | `93` (R) | **`R_c ≤ 85`** (Rosowski formula); 93 from Waksman per DIS09; Rosowski has approximate non-bilinear at 89 | [drisc09], [rosowski19] |

For non-square formats with similar UB/LB tables in the commutative
case, see Drevet–Islam–Schost 2009 §3-4.

**The commutative rank is a theoretical baseline.** In the bilinear
complexity model (which is what matters for practical algorithms),
**non-commutative rank is the relevant complexity measure** throughout.

### 1.2bis How field choice affects rank

A rank bound is always **field-specific**: `R_K(T)` is defined per field `K`.
The relationships between bounds across fields are NOT symmetric — pay
close attention to direction.

#### Upper bounds (algorithms) transfer DOWNWARD

If `K ⊆ K'` (K is a subfield of K'), then **every K-bilinear decomposition
is also a K'-bilinear decomposition** (the coefficients still make sense in
the larger field). Therefore:

> `R_{K'}(T) ≤ R_K(T)` whenever `K ⊆ K'`.

Practical consequences:

- `R_C(T) ≤ R_R(T) ≤ R_Q(T) ≤ R_Z(T)` — a richer field permits more
  algorithm-design flexibility, so the rank can only stay the same or
  drop. **AlphaEvolve's `R_C(⟨4,4,4⟩) ≤ 48` does not (yet) imply
  `R_R(⟨4,4,4⟩) ≤ 48`** — the 48-mult algorithm uses complex
  coefficients that don't (necessarily) realise as real numbers.
- Strassen's algorithm uses `{-1, 0, +1}` coefficients only, so it works
  over `Z`, `Q`, `R`, `C` *and* (after reducing mod 2) over `F_2`.

#### Practical fallback rule (used by Pages and lookup code)

The downward-transfer property dictates a clean **fallback chain** for any
caller that needs a scheme for field `K`:

1. Use the catalog's best entry tagged exactly `K`.
2. If no `K` entry exists, fall through to schemes in any subfield
   `K' ⊂ K` — they're automatically valid in `K`.
3. Never fall through to schemes in `K'' ⊃ K` (they may use coefficients
   `K` can't represent), and never cross to `F_p` (different arithmetic).

So in practice:

| requested field | valid catalog sources (in priority order) |
|---|---|
| `Z` (= `Z`) | `Z` only |
| `Q` | `Q`, then `Z` |
| `R` | `R`, then `Q`, then `Z` |
| `C` | `C`, then `R`, then `Q`, then `Z` |
| `F_2` | `F_2` only — no fallback to / from characteristic-0 |

This rule is encoded in:

- **`BlockSplitSearch.loadCatalogBestRanksForField`** — the `"C"` filter
  permits R/Q/Z files; the `"R"` filter does not permit C files.
- **`docs/catalog.js`** Pages filter — when the user picks `C`, R-class
  rows show up too, marked `R/Q/Z ↗ C` so it's clear they're
  fallback-valid rather than C-native.

The reverse (e.g. "I asked for `R` and got a `C` scheme") is **a bug**:
catalog state and the UI should never return a strict super-field
scheme when a sub-field was requested, because the algorithm may use
unrepresentable coefficients.

#### Same rule for the commutative ↔ non-commutative axis

The commutative-vs-non-commutative distinction (§ "Non-commutative vs
commutative rank" above) follows the **same fallback structure** as
the field hierarchy:

> A **non-commutative** algorithm makes no commutativity assumption
> about its scalar inputs, so it remains valid when applied to
> commuting scalars. Hence:
> **`R_c(T) ≤ R(T)`** and any non-commutative scheme is automatically
> valid for the commutative case.

The reverse does NOT hold — a commutative-only scheme exploits
`a·b = b·a` and breaks when applied recursively to matrix entries
(which don't commute). This is why the matmul literature defaults to
non-commutative rank and why **every scheme in this catalog is
non-commutative** (no separate commutative tag needed).

Practical fallback chain combining both axes (in priority order from
strictest):

| requested mode | valid catalog sources |
|---|---|
| **Z, non-commutative** | only Z, non-commutative |
| **Q, non-commutative** | Q, then Z (always non-commutative in this catalog) |
| **R, non-commutative** | R, then Q, then Z |
| **C, non-commutative** | C, then R, then Q, then Z |
| **F₂, non-commutative** | only F₂ |
| **commutative** (any field) | any non-commutative scheme of the same field, with the same field-fallback chain |

Implementation: [`FieldAwareLookup`](../../src/main/java/io/cormoran/strassen/v3/catalog/FieldAwareLookup.java)
covers the field axis; since the catalog is entirely non-commutative,
no separate commutativity filter is needed today. If commutative-only
schemes (e.g. from
[Drevet–Islam–Schost 2009](../../REFERENCES.md#10-drisc09)) are imported in
future, tag them with a `commutative: true` JSON field and add the
corresponding fallback skip-rule.

#### Axis-orientation fallback (S₃ tensor symmetry)

The matmul tensor `T_{n,m,p}` is invariant under the full
**symmetric group S₃** on its three slots — generated by:

- **Cyclic shift** (`T(A,B,C) = T(B,C,A)`): `⟨n,m,p⟩ ≅ ⟨m,p,n⟩ ≅ ⟨p,n,m⟩`
- **Transpose** (`(AB)^T = B^T A^T`): `⟨n,m,p⟩ ≅ ⟨p,m,n⟩`

The 6 axis orderings of `⟨a,b,c⟩` (sorted `a ≤ b ≤ c`) all have the
same tensor rank. So when looking up a scheme for `⟨9,12,11⟩` and only
`⟨9,11,12⟩` exists in the catalog, the scheme can be re-oriented via
the appropriate S₃ element (here: transpose ∘ cyclic²) to produce the
requested orientation with the same rank.

Implementation: [`NonCubicBilinearAlgorithm.orientAs`](../../src/main/java/io/cormoran/strassen/v3/NonCubicBilinearAlgorithm.java)
tries all 6 S₃ orbit elements (3 cyclic shifts of this, plus 3 of its
`transpose()`); [`FieldAwareLookup.find`](../../src/main/java/io/cormoran/strassen/v3/catalog/FieldAwareLookup.java)
uses `orientAs` to fold all axis orderings back to the canonical
sorted form when searching the catalog.

#### `F_p` is a quotient, not a subfield

Going from `Z` to `F_p` is **not** a field embedding — it's a quotient by
the ideal `(p)`. Coefficients reduce modulo `p`, and some terms may cancel
that did not cancel in `Z`. So:

> An algorithm over `Z` reduces to an algorithm over `F_p` with possibly
> fewer terms (or possibly more cancellation in the XOR sum), giving
> **`R_{F_p}(T) ≤ R_Z(T)`**.

But the reverse does NOT hold: a low-rank algorithm over `F_2` does not
lift to a low-rank algorithm over `Z` (or `R`). The lifted operations
involve coefficients that can be `0` or `1` modulo `2` but `2`, `3`, or
anything else in `Z` — a free lift is a Hobson's choice that typically
doesn't satisfy the integer equations.

This is why **AlphaTensor 2022's `R_{F_2}(⟨4,4,4⟩) ≤ 47` does not imply
`R(⟨4,4,4⟩) ≤ 47` over `R` or `Z`** — the AlphaTensor scheme is
F_2-specific.

#### Lower bounds (impossibility) are direction-aware

`R_K(T) ≥ L` means "no algorithm over `K` can use fewer than `L`
multiplications." This bound does **not** automatically transfer to other
fields.

| direction | does it transfer? |
|---|---|
| `R_K(T) ≥ L`, K ⊆ K' | **NO** — `K'` has more flexibility; `R_{K'}` could be lower. Example: `R_R(⟨3,3,3⟩) ≥ 19` (Bläser 2003) does not imply `R_C(⟨3,3,3⟩) ≥ 19`. |
| `R_K(T) ≥ L`, K ⊇ K' | YES, often — fewer coefficient options ↔ harder. Example: `R_C ≥ L` ⇒ `R_R ≥ L`. |
| `R_{F_p}(T) ≥ L` | **NO direct implication** for `R_R`, `R_Z`, etc. — different field, different problem. |

#### What Wang 2026's `R_{F_2}(⟨3,3,3⟩) ≥ 20` actually says

[Wang 2026][wang26] proves the bound **over `F_2` specifically**.

| claim | true after Wang? |
|---|---|
| `R_{F_2}(⟨3,3,3⟩) ≥ 20` | **YES** — what the paper proves |
| `R(⟨3,3,3⟩) ≥ 20` over `R`/`Q`/`Z` | **NO** — still `≥ 19` (Bläser 2003); Wang's technique may or may not extend |
| `R_C(⟨3,3,3⟩) ≥ 20` | **NO** — `C` is a richer field; LB over `R` doesn't imply LB over `C` |
| `R_{F_2}(⟨3,3,3⟩) ≤ 23` | YES — Laderman's `±1` algorithm reduces mod 2 |
| Therefore `R_{F_2}(⟨3,3,3⟩) ∈ [20, 23]` | YES — gap of 3 |
| And `R(⟨3,3,3⟩) ∈ [19, 23]` over `R` | YES — gap of 4, unchanged since 2003 |

#### Worked example: why Wang 2026's `⟨2,3,3⟩` result is news

`R(⟨2,3,3⟩) = 15` over `R`, `Q`, `Z`, `C` has been known tight since
[Hopcroft–Kerr 1971][hk71] — their argument uses characteristic-0 algebra.

But characteristic 0 doesn't cover `F_2`, and the bound there was an
**independent open question** for 50+ years, because:

1. `F_2` has the richest cancellation patterns of any small field — every
   nonzero element equals `1`, so `a + a = 0` for every `a`. AlphaTensor
   exploited exactly this to find `R_{F_2}(⟨4,4,4⟩) ≤ 47`, which does NOT
   lift to `R`.
2. So `F_2` was the most plausible field where some sub-15 algorithm could
   hide for `⟨2,3,3⟩`.

Wang shut that door: `R_{F_2}(⟨2,3,3⟩) ≥ 15`, matching the upper bound.

What's still open for **other prime fields** (`F_3`, `F_5`, `F_7`, ...):
each one's `⟨2,3,3⟩` lower bound is technically a separate problem,
not implied by either Hopcroft–Kerr or Wang. No improvement is known,
but no proof has been published either.

#### Practical takeaway for this codebase

When citing or comparing ranks: **always state the field**. The catalog
([SMALL_MATMUL_CATALOG.md](SMALL_MATMUL_CATALOG.md) §3) makes `field` an explicit column for
exactly this reason. A SAT proof of `UNSAT at r=22` over `F_2` would
tighten `R_{F_2}(⟨3,3,3⟩) ≥ 23` (closing the F_2 gap), but would not by
itself say anything about `R(⟨3,3,3⟩)` over `R`.

### 1.3 Non-cubic formats (most useful for `⟨3,3,3⟩` embedding analysis)

| `⟨n,m,p⟩` | naive `nmp` | rank `R` | status | references |
|---|---|---|---|---|
| `⟨2,2,3⟩` | 12 | **11** | proven tight | [Hopcroft-Kerr 1971][hk71], [Smirnov][smir13] |
| `⟨2,2,4⟩` | 16 | **14** | proven tight | [Pan 1978][pan78], [Hopcroft-Kerr 1971][hk71] |
| `⟨2,2,5⟩` | 20 | **18** | tight | [Hopcroft-Kerr 1971][hk71] |
| `⟨2,3,3⟩` | 18 | **15** | tight (UB Smirnov, LB folklore-tight) | [Smirnov 2013][smir13] |
| `⟨2,3,4⟩` | 24 | **20** | tight | [Smirnov 2013][smir13] |
| `⟨2,4,4⟩` | 32 | **26** (latest), `≤ 27` historically | open lower bound | [AlphaTensor 2022][alpha22] |
| `⟨3,3,4⟩` | 36 | **`29 ≤ R ≤ 29`** (= tight per recent work) | mostly tight | [Smirnov 2013][smir13] |
| `⟨3,3,5⟩` | 45 | `≤ 36` | open | [Smirnov 2013][smir13] |
| `⟨3,4,5⟩` | 60 | `≤ 47` | open | [Smirnov 2013][smir13] |
| `⟨4,5,5⟩` | 100 | `≤ 76` (`Z/2`) | open | [AlphaTensor 2022][alpha22] |

**Notes:**
- Lower bounds for non-square formats are often less explicit in the
  literature than for square formats. The "tight" entries above with
  named UB authors typically have implicit folklore tight lower bounds
  via direct counting / Hopcroft-Kerr substitution method.
- `⟨2,m,p⟩` formats have the most reliable bounds. `⟨3,m,p⟩` becomes
  open quickly.

---

## 2. Addition counts (separate optimization axis)

Beyond multiplications, a fast matmul algorithm has some number of
**additions** (linear combinations of inputs/outputs). For `⟨n,n,n⟩` at
fixed rank `r`, the addition count varies between algorithms in the same
gauge orbit and is a separate research target.

| algorithm | mults | additions | year | reference |
|---|---|---|---|---|
| Strassen `⟨2,2,2⟩` r=7 | 7 | 18 (original) → 15 (Winograd) | 1969/1971 | [Strassen 1969][s69], [Winograd 1971][win71] |
| Laderman `⟨3,3,3⟩` r=23 | 23 | 62 (Laderman) | 1976 | [Laderman 1976][lad76] |
| Smirnov `⟨3,3,3⟩` r=23 | 23 | 56 | 2013-2017 | [Smirnov 2013][smir13] |
| Heun et al. `⟨3,3,3⟩` r=23 | 23 | 60-something | 2019 | [Heule et al.][heun18] |
| Stapleton (no basis change) `⟨3,3,3⟩` r=23 | 23 | **60** | 2025-08 | [arXiv:2508.03857][add60] |
| Mårtensson–Stankovski Wagner–Stapleton (no basis change) `⟨3,3,3⟩` r=23 | 23 | **59** | 2025-12 | [arXiv:2601.05272][mss25] |
| Perminov (no basis change) `⟨3,3,3⟩` r=23 | 23 | **58** | 2025-12 | [arXiv:2512.21980][add58] |

The "best addition count" is an active research target. AlphaTensor and
flip-graph methods (Moosbauer, Kauers, etc.) have produced reduced-
addition variants of known rank-23 schemes.

---

## 3. Z/2-specific results (when different from `R`/`Q`)

Over the field of two elements (boolean arithmetic), rank can differ
from standard arithmetic. **AlphaTensor 2022** found algorithms with
lower `Z/2` rank than the best known `R`-arithmetic rank for several
small formats:

| format | `R(⟨n,m,p⟩)` over `R` | `R_{Z/2}(⟨n,m,p⟩)` | reduction |
|---|---|---|---|
| `⟨4,4,4⟩` | ≤ 49 (Strassen²) | **≤ 47** | −2 |
| `⟨4,5,5⟩` | ≤ 80 (recursive) | **≤ 76** | −4 |
| `⟨5,5,5⟩` | ≤ 100 (recursive) | **≤ 96** | −4 |

For `⟨2,2,2⟩` and `⟨3,3,3⟩`, no `Z/2`-specific reduction is known —
`R_{Z/2}` likely equals `R_Q` at these formats, but neither has been
proven tight separately.

Reference: [Fawzi et al., Nature 2022][alpha22].

---

## 4. The matmul exponent `ω`

Historical progression:

| year | author(s) | `ω` upper bound | method |
|---|---|---|---|
| pre-1969 | — | 3 | naive |
| 1969 | Strassen | **2.807** | rank `⟨2,2,2⟩` ≤ 7 |
| 1978 | Pan | 2.78 | trilinear aggregation |
| 1979 | Bini et al. | 2.78 | border rank, `⟨2,2,2⟩` ≤ 5 |
| 1981 | Schönhage | 2.522 | `τ`-theorem |
| 1987 | Coppersmith-Winograd | 2.376 | laser method |
| 2010 | Stothers | 2.3727 | refined CW |
| 2014 | Le Gall | 2.3729 | refined CW (a different way) |
| 2020 | Alman-Vassilevska Williams | 2.37286 | new analysis |
| **2024** | Duan-Wu-Zhou-Niwa-Williams | **2.371552** | current best |

Lower bound: `ω ≥ 2` (trivial — output size `n²`). Currently believed
likely `ω = 2` ("conjecture of high tensor rank"), but **no proof**.

**Practical algorithms in BLAS/numerical software use `ω = 3` (naive)
or Strassen's `ω = 2.807`** — the asymptotic-frontier algorithms are
not used in practice because their hidden constants and condition
numbers are astronomical.

---

## 5. Lower-bound techniques

For each format, the lower bound on `R` typically comes from one of:

| technique | best result for `⟨3,3,3⟩` | reference |
|---|---|---|
| Substitution method | `R(⟨3,3,3⟩) ≥ 19` | [Bläser 2003][bla03] |
| Orbit-DP + verifiable certificates | `R_{F₂}(⟨3,3,3⟩) ≥ 20` | [Wang 2026][wang26] — methodology directly relevant to the SAT pipeline in this repo |
| Border-rank lower bounds | `R̃(⟨3,3,3⟩) ≥ 15` (or similar) | various |
| Hopcroft-Kerr (small formats) | tight for `⟨2,2,k⟩` family | [HK 1971][hk71] |
| Asymptotic / laser method | `ω ≥ 2` only | — |

**Current `⟨3,3,3⟩` rank lower bound gap**: `[19, 23]`, open since
1976 (UB) / 2003 (LB).

---

## 6. The frontier

- **`R(⟨3,3,3⟩) = 22?`** — would be first improvement on Laderman in
  50 years. Doesn't directly improve `ω` (`log₃ 22 ≈ 2.814 > log₂ 7 = 2.807`)
  but is a major algebraic-complexity result.
- **`R(⟨3,3,3⟩) = 21?`** — would beat Strassen recursively
  (`log₃ 21 ≈ 2.771 < log₂ 7`), a 50-year `ω` improvement at the small-
  format level.
- **`R(⟨2,3,3⟩) = 14?`** — would improve a tight 1971 result. Believed
  tight at 15 but the lower-bound proof is less explicit than for
  `⟨2,2,2⟩`. SAT-certified `R_{Z/2}(⟨2,3,3⟩) ≥ 15` would be a
  computational re-validation.
- **`R̃(⟨3,3,3⟩) < 21`?** — would directly improve `ω` via Schönhage's
  `τ`-theorem.
- **Cohn-Umans group-theoretic conjecture**: `ω = 2` via specific group
  algebras. Not constructively verified.

---

## References

- <a name="s69"></a>[s69] V. Strassen, "Gaussian Elimination is not Optimal,"
  *Numerische Mathematik* 13:354–356, 1969.
- <a name="win71"></a>[win71] S. Winograd, "On Multiplication of 2×2 Matrices,"
  *Linear Algebra and Its Applications* 4(4):381–388, 1971.
- <a name="hk71"></a>[hk71] J. Hopcroft, L. Kerr, "On Minimizing the Number of
  Multiplications Necessary for Matrix Multiplication," *SIAM J. Applied Math.*
  20(1):30–36, 1971.
- <a name="lad76"></a>[lad76] J. Laderman, "A noncommutative algorithm for multiplying
  (3×3) matrices using 23 multiplications," *Bull. AMS* 82(1):126–128, 1976.
  ([reproduced in this repo as `v3/Laderman23.java`](../../src/main/java/io/cormoran/strassen/v3/Laderman23.java))
- <a name="pan78"></a>[pan78] V. Pan, "Strassen's Algorithm is Not Optimal: Trilinear
  Technique of Aggregating, Uniting and Canceling for Constructing Fast
  Algorithms for Matrix Operations," *FOCS 1978*.
- <a name="bini79"></a>[bini79] D. Bini, M. Capannini, F. Lotti, F. Romani, "O(n^{2.7799})
  Complexity for n×n Approximate Matrix Multiplication,"
  *Information Processing Letters* 8(5):234–235, 1979.
- <a name="sch81"></a>[sch81] A. Schönhage, "Partial and Total Matrix Multiplication,"
  *SIAM J. Computing* 10(3):434–455, 1981.
- <a name="cw87"></a>[cw87] D. Coppersmith, S. Winograd, "Matrix Multiplication via
  Arithmetic Progressions," *J. Symbolic Computation* 9(3):251–280, 1990.
- <a name="wang26"></a>[wang26] C. Wang, *Automated Lower Bounds for Small Matrix Multiplication
  Complexity over Finite Fields*, arXiv:2603.07280 (2026).
  [arXiv](https://arxiv.org/abs/2603.07280). Proves `R_{F₂}(⟨3,3,3⟩) ≥ 20`
  via orbit classification + dynamic programming + verifiable proof
  certificates. The methodology — generating compact certificates a SAT
  solver could check — is directly relevant to this repo's pipeline.
- <a name="bla03"></a>[bla03] M. Bläser, "On the Complexity of the Multiplication of
  Matrices of Small Formats," *J. Complexity* 19(1):43–60, 2003.
- <a name="smir13"></a>[smir13] A. V. Smirnov, "Bilinear complexity and practical
  algorithms for matrix multiplication," *Computational Mathematics and Mathematical
  Physics* 53(12):1781–1795, 2013.
- <a name="heun18"></a>[heun18] M. J. H. Heule, M. Kauers, M. Seidl, "Local Search for Fast
  Matrix Multiplication," *International Conference on Theory and Applications of
  Satisfiability Testing*, 2019. [arXiv:1903.11391](https://arxiv.org/abs/1903.11391)
- <a name="alpha22"></a>[alpha22] A. Fawzi et al., "Discovering faster matrix
  multiplication algorithms with reinforcement learning," *Nature* 610:47–53,
  October 2022.
  [Algorithms on GitHub](https://github.com/google-deepmind/alphatensor)
- <a name="add58"></a>[add58] A. I. Perminov, *A 58-Addition, Rank-23 Scheme for General 3×3 Matrix
  Multiplication*, arXiv:2512.21980, Dec 2025. Local PDF: `references/papers/perminov_2025_3x3x3_r23_a58_arxiv2512.21980.pdf`.
- <a name="add60"></a>[add60] J. Stapleton, *A 60-Addition, Rank-23 Scheme for Exact 3×3 Matrix
  Multiplication*, arXiv:2508.03857, Aug 2025. Local PDF: `references/papers/stapleton_2025_3x3x3_r23_a60_arxiv2508.03857.pdf`.
- <a name="mss25"></a>[mss25] E. Mårtensson, P. Stankovski Wagner, J. Stapleton,
  *A Rank 23 Algorithm for Multiplying 3 × 3 Matrices with an Arithmetic Complexity of 59*,
  arXiv:2601.05272, 2025. **59 additions** (no basis change), local PDF
  in `references/papers/`.
- <a name="rosowski19"></a>[rosowski19] A. Rosowski, *Fast Commutative Matrix Algorithm*,
  arXiv:1904.07683 (v2: 2020). **Commutative** `R_c(⟨3,3,3⟩) ≤ 21`
  (beats Makarov 22). Closed-form: `n(lm+l+m−1)/2` for n,m odd;
  `(n(lm+l+m−1)+l−1)/2` for n odd, m even. Local PDF in `references/papers/`.
- <a name="perminov"></a>[perminov] Perminov, Andrew I., *Fast Matrix Multiplication catalog*,
  [github.com/dronperminov/FastMatrixMultiplication](https://github.com/dronperminov/FastMatrixMultiplication#research-findings--status).
  GitHub username `dronperminov`; cite as "Perminov". The methods/results
  behind the repo are published in arXiv:2603.02398 (small-formats flip-graph
  framework), arXiv:2511.20317 (ternary meta flip graphs), arXiv:2512.13365
  (additive-complexity reduction) and arXiv:2512.21980 (`⟨3,3,3⟩` r23/a58) —
  see REFERENCES.md [75], [77], [78], [79].
- <a name="drisc09"></a>[drisc09] C. Drevet, M. Islam, É. Schost,
  *Optimization techniques for small matrix multiplication*, 2009.
  [Preprint](https://cs.uwaterloo.ca/~eschost/publications/DrIsSc09.pdf) —
  contains separate tables for commutative and non-commutative ranks across
  many small formats, with explicit bounds and references.

[s69]: #s69
[win71]: #win71
[hk71]: #hk71
[lad76]: #lad76
[pan78]: #pan78
[bini79]: #bini79
[sch81]: #sch81
[cw87]: #cw87
[bla03]: #bla03
[wang26]: #wang26
[smir13]: #smir13
[heun18]: #heun18
[alpha22]: #alpha22
[add58]: #add58
[add60]: #add60
[mss25]: #mss25
[rosowski19]: #rosowski19
[perminov]: #perminov
[drisc09]: #drisc09

---

## Caveats

- **Bounds evolve.** This document reflects the state of the literature
  as of 2026-05; AlphaTensor/AlphaEvolve, Moosbauer-Kauers flip-graph
  work, and follow-up SAT-search efforts produce new results regularly.
  Check the [FMM Catalogue](https://fmm.univ-lille.fr) and
  [Perminov's catalog][perminov] for the latest.
- **Some lower bounds are "folklore tight"** — established by multiple
  independent verifications but not always rigorously published as
  papers. The `⟨2,3,3⟩ ≥ 15` bound is one such case; SAT-certified
  re-verification would be useful.
- **`R(T)` vs `R_{Z/2}(T)`** can differ. AlphaTensor 2022 found `Z/2`
  algorithms with strictly lower rank than the best `R`-arithmetic
  algorithm for several formats.
- **Border rank `R̃` is not always known.** For most formats above,
  border rank is an active research area with looser bounds than rank.
- **Non-commutative vs commutative.** All bounds in §1.1 and §1.3 are
  non-commutative (the practical version). §1.2 lists the commutative
  ranks for the same formats — these are theoretical lower-bound
  baselines but don't apply to matmul speedups. See
  [Drevet–Islam–Schost 2009][drisc09] for the most comprehensive side-
  by-side comparison.
