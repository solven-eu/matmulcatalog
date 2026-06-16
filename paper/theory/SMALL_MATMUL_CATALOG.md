# Small Matrix-Multiplication Algorithm Catalog

An attempt at an up-to-date catalog of **fast bilinear matrix-multiplication
algorithms** for small formats up to `⟨32, 32, 32⟩`, with explicit per-entry
provenance (who, when, where, in which arithmetic setting). Designed as the
basis for recursive-composition experiments in the style of
[AlphaEvolve][alphaevolve].

> Companion to [RANK_KNOWLEDGE.md](RANK_KNOWLEDGE.md) (theoretical bounds and ω history) and
> [SOLVING_STRATEGIES.md](SOLVING_STRATEGIES.md) (search pipeline). Where the rank-knowledge file
> states *what's true*, this catalog records *who proved it, when, and via
> which algorithm*.

---

## 1. Scope

- **Formats**: cubic `⟨n,n,n⟩` for `n ∈ {2..32}`; key non-cubic `⟨n,m,p⟩` that
  embed productively (e.g. `⟨2,2,3⟩`, `⟨2,3,4⟩`, `⟨4,5,5⟩`).
- **Arithmetic settings tracked separately** (this matters — best known
  rank depends on the field):
  - `Z/2` — the field of two elements (synonyms: `F_2`, `F₂`, `GF(2)`,
    `Z/2Z`). Boolean arithmetic, addition = XOR.
  - `Z` / `Q` / `R` — "standard integer/rational/real" (typically same rank)
  - `C` — complex (can be strictly better; see `⟨4,4,4⟩`)
  - "commutative" — bilinear-complexity model that assumes `a·b = b·a`;
     theoretical baseline only ([RANK_KNOWLEDGE.md](RANK_KNOWLEDGE.md) §1.2).
- **Cross-field rank relationships** (which lower/upper bounds transfer
  between fields, and which do not) are explained in [RANK_KNOWLEDGE.md](RANK_KNOWLEDGE.md)
  §1.2bis — read this before drawing conclusions like "AlphaTensor's
  47-mult `⟨4,4,4⟩` algorithm proves `R(⟨4,4,4⟩) ≤ 47` over reals" (it
  does not — the result is `F_2`-specific).
- **What counts as a "first"**: the first published algorithm (proof,
  certificate, or RL-discovered scheme) achieving a particular rank
  upper bound. Lower bounds get a separate column.
- **What's *not* in scope**: asymptotic results (those live in
  [RANK_KNOWLEDGE.md](RANK_KNOWLEDGE.md) §4); border-rank-only results without an
  exact-rank counterpart; partial results (e.g. "≤ R + ε" for some
  unspecified ε).

## 2. Sources we are synthesizing

| key | catalog / source | scope | note |
|---|---|---|---|
| [perminov][perminov] | Perminov, A., *FastMatrixMultiplication* (GitHub: dronperminov) | up to `⟨16,16,16⟩` | GitHub-hosted, structured per-scheme data + provenance |
| [drisc09][drisc09] | Drevet–Islam–Schost 2009 | up to `⟨30,30,30⟩` | Separate commutative vs non-commutative tables. **Predates AlphaTensor / AlphaEvolve.** |
| [smirnov][smirnov] | A. V. Smirnov, *Bilinear complexity and practical algorithms for matrix multiplication* | many small non-cubic | The integer/{-1,0,+1} catalog |
| [alphatensor][alphatensor] | DeepMind AlphaTensor (Fawzi et al. 2022) | small Z/2 + a few standard-arithmetic | RL discovery, single-format optima |
| [alphaevolve][alphaevolve] | DeepMind AlphaEvolve (Novikov et al. 2025) | new `⟨4,4,4⟩` complex result; recursive composition framework | "first improvement after 56 years over Strassen in this setting" |

This catalog supersedes the per-format upper-bound tables in
[drisc09][drisc09] where AlphaTensor / AlphaEvolve have since improved them.

## 3. Per-entry schema

Each catalog row records:

| column | meaning |
|---|---|
| format | `⟨n,m,p⟩` |
| field | `Z/2`, `Z`/`Q`/`R` (collapsed when same), `C`, commutative |
| **rank** | best known **upper bound** (number of multiplications) |
| **LB** | best known **lower bound** (when known to be sharp, rank = LB) |
| year | publication year of the first scheme attaining `rank` |
| source | algorithm name + author(s) |
| ref | citation key (see §8) |
| notes | improvement story, alternates, dependencies on lower-level algorithms |

## 4. Cubic formats up to `⟨32, 32, 32⟩`

### 4.1 `⟨2,2,2⟩`

| field | rank | LB | year | source | ref | notes |
|---|---|---|---|---|---|---|
| any non-commutative | **7** | **7** | 1969 | Strassen | [strassen69][strassen69] | Tight forever — lower bound proved [hopcroft-kerr 71][hk71] |
| commutative | **6** | **6** | 1971 | Hopcroft, Winograd (independently) | [hk71][hk71], [winograd71][win71] | Tight but doesn't recurse to matrix entries |

### 4.2 `⟨3,3,3⟩`

| field | rank | LB | year | source | ref | notes |
|---|---|---|---|---|---|---|
| any non-commutative | **23** | 19 | 1976 | Laderman | [laderman76][lad76] | LB `19` from [bläser 2003][blaser03]; open gap `[19, 23]` |
| Z/2 | 23 | **20** | 1976 | Laderman | [lad76][lad76] | LB **20** from [wang 2026][wang26] (tightened from 19 via SAT-style orbit-DP certificates over `F₂`) |

### 4.3 `⟨4,4,4⟩`

| field | rank | LB | year | source | ref | notes |
|---|---|---|---|---|---|---|
| Z/2 | **47** | (open) | 2022 | AlphaTensor | [alphatensor][alphatensor] | Strassen² gave 49 (1969); AlphaTensor's main headline |
| C | **48** | (open) | 2025 | AlphaEvolve | [alphaevolve][alphaevolve] | First improvement on Strassen² (49) in this setting *"after 56 years"* — for complex matrices specifically |
| R, Q | 49 | (open) | 1969 | Strassen² (recursive) | [strassen69][strassen69] | AlphaEvolve's `48` is over `C`; whether `48` reaches `R` is open |

### 4.4 `⟨5,5,5⟩`

| field | rank | LB | year | source | ref | notes |
|---|---|---|---|---|---|---|
| Z/2 | **96** | (open) | 2022 | AlphaTensor | [alphatensor][alphatensor] | Improves Smirnov's prior best |
| Q / R | ~98 | (open) | various | Smirnov + improvements | [smirnov][smirnov] | Different sources cite 96-100 depending on field assumptions; verify against [perminov][perminov] |

### 4.5 Powers of 2 up to `⟨32, 32, 32⟩` (via recursive composition)

Pure recursion bounds, no algorithm-mixing yet:

| format | Strassen recursion `R = 7^k` | best mixed (lower) | notes |
|---|---|---|---|
| `⟨2,2,2⟩` | `7^1 = 7` | 7 | base |
| `⟨4,4,4⟩` | `7^2 = 49` | **47** (Z/2), **48** (C) | see §4.3 |
| `⟨8,8,8⟩` | `7^3 = 343` | **329** (Z/2: 47·7) | use AlphaTensor for the `⟨4,4,4⟩` level |
| `⟨16,16,16⟩` | `7^4 = 2,401` | **2,209** (Z/2: 47²) | use AlphaTensor at each `⟨4,4,4⟩` level |
| `⟨32,32,32⟩` | `7^5 = 16,807` | **15,463** (Z/2: 47²·7) | best mixed: `R(⟨16,16,16⟩) · R(⟨2,2,2⟩)` |

For complex matrices using AlphaEvolve's `⟨4,4,4⟩ = 48`:

| format | composition | rank | notes |
|---|---|---|---|
| `⟨8,8,8⟩` C | `48·7` | **336** | |
| `⟨16,16,16⟩` C | `48²` | **2,304** | |
| `⟨32,32,32⟩` C | `48²·7` | **16,128** | improves Strassen `16,807` by ~4% |

**Caveats**:
- These are *composition* bounds, not independently-verified ranks. A
  direct algorithm could in principle be lower.
- Mixed compositions (e.g. Strassen at outer level + AlphaTensor inner)
  rely on the inner algorithm being valid over the same field — for
  recursing through matrix entries, "valid over the entry field" is the
  requirement.
- The Z/2 catalog assumes the outer matrices are over Z/2. Recursing
  AlphaTensor `⟨4,4,4⟩ = 47` for matrix entries over **`Z` or `R`** is
  not valid because the 47-scheme is field-specific to Z/2.

### 4.6 Non-power-of-2 cubic (selected)

| format | field | rank | year | source | ref |
|---|---|---|---|---|---|
| `⟨6,6,6⟩` | any | 153 | 1980s | Smirnov / Pan | [smirnov][smirnov] |
| `⟨7,7,7⟩` | any | ~250 | various | catalogs | [perminov][perminov] |
| `⟨9,9,9⟩` | any | 521 (Laderman² gives 529) | 2017 | Smirnov | [smirnov][smirnov] |
| `⟨12,12,12⟩` | any | ~~? — verify against | [perminov][perminov] |  |
| `⟨27,27,27⟩` | any | Laderman³ gives `23³ = 12,167` | 1976 | composition | base bound |

> Entries marked `~`/`?` are placeholders — verify against
> [perminov][perminov]'s structured DB before citing. Pull requests welcome.

## 5. Key non-cubic formats

These are useful for embedding-based lower bounds (see [RANK_KNOWLEDGE.md](RANK_KNOWLEDGE.md) §1.3)
and as ingredients for non-square recursive composition.

| format | field | rank | LB | year | source | ref |
|---|---|---|---|---|---|---|
| `⟨2,2,3⟩` | non-commutative | **11** | 11 | 1971 | Hopcroft–Kerr (tight) | [hk71][hk71] |
| `⟨2,3,3⟩` | non-commutative | **15** | 15 | ~1973 | Hopcroft–Kerr / Pan | [hk71][hk71], [pan78][pan78] |
| `⟨2,3,4⟩` | non-commutative | **20** | 20 | various | (verify) | [perminov][perminov] |
| `⟨3,3,4⟩` | non-commutative | **29** | (open) | various | (verify) | [perminov][perminov] |
| `⟨4,5,5⟩` | Z/2 | **76** | (open) | 2022 | AlphaTensor (improved from 80) | [alphatensor][alphatensor] |
| `⟨4,4,5⟩` | Z/2 | **63** | (open) | 2022 | AlphaTensor | [alphatensor][alphatensor] |

## 6. Recursive composition framework

The rank of a composed matmul is multiplicative across nested factorizations.

**Theorem (folklore, easy to verify)**: if `⟨n,m,p⟩ = ⟨n_1·n_2, m_1·m_2, p_1·p_2⟩`,
then

```
R(⟨n,m,p⟩) ≤ R(⟨n_1,m_1,p_1⟩) · R(⟨n_2,m_2,p_2⟩)
```

— provided both subalgorithms are valid over the entry-level field
(this is why non-commutativity matters: the entries become matrices,
which don't commute, so the inner algorithm must work non-commutatively).

**Three composition styles**:

1. **Pure (single base × power)** — Strassen iterated: `7^k`.
   Simple, never optimal beyond `k=1` once better small algorithms exist.
2. **Two-level mixed** — outer one format, inner another. The `Z/2`
   bound `47²·7 = 15,463` for `⟨32,32,32⟩` is this style.
3. **Multi-level mixed / discovered** — AlphaEvolve's approach: search
   over compositions to find the optimal tree, including non-cubic
   intermediate factorizations.

AlphaEvolve's `⟨4,4,4⟩ = 48` over `C` was specifically targeted because
that level appears in many compositions toward `⟨32,32,32⟩`; improving
even one level cascades.

### 6.1 Implementation

The Kronecker construction is implemented in
`io.cormoran.strassen.v3.catalog.Compose`:

- `Compose.kronecker(outer, inner)` — composes two cubic
  {@link BilinearAlgorithm}s into one for the product format.
- `Compose.chain(List<BilinearAlgorithm>)` — left-fold a chain.
- `Compose.strassenPower(k)` — convenience: `k` levels of Strassen.

Tests in `TestCompose` exercise: `Strassen²→⟨4,4,4⟩=49`,
`Strassen³→⟨8,8,8⟩=343`, `Strassen⁴→⟨16,16,16⟩=2401`,
`Strassen⁵→⟨32,32,32⟩=16,807` (construction only — the verifier's
`O(n⁶·r)` cost makes re-verification impractical beyond `n=8`), and
`Strassen × Laderman → ⟨6,6,6⟩=161` (verified exact).

**Scope limitation**: cubic-to-cubic only. To compose with non-cubic
algorithms (e.g. mix `⟨2,2,3⟩` factors into the chain), a non-cubic
`BilinearAlgorithm` class is needed first. The encoder is already
non-cubic-aware ([SOLVING_STRATEGIES.md](SOLVING_STRATEGIES.md) §10.3.0.5), but the
algorithm container isn't yet.

## 7. Where this codebase fits

The strassen repo's role in this catalog:

- **Verification layer**: `Verifier.java` exactly checks any candidate
  decomposition. Every catalog entry should round-trip through a stored
  factor file (planned: `catalog/algorithms/` subdir).
- **Discovery layer**:
  - [SAT pipeline][solving_strategies] — small Z/2 over-rank certifying
    UNSAT (e.g. `⟨2,3,3⟩ r=14` is the current research target —
    _memory: `project-next-steps-a2-sat-als-warm-start`_).
  - [Z3Als][z3als] — equivariant ALS for raw `Z/3`-symmetric tensors.
- **Composition layer** (planned): an `io.cormoran.strassen.v3.catalog`
  Java package: `KnownAlgorithm`, `KnownAlgorithmCatalog`,
  `RecursiveComposition`. The package data model captures everything
  in §3 above; tooling can iterate over the catalog and compute mixed
  compositions automatically.

### 7.1 Open problems where this codebase could contribute

| problem | format | field | current | target | approach |
|---|---|---|---|---|---|
| Confirm `R_{Z/2}(⟨2,3,3⟩) = 15` | `⟨2,3,3⟩` | Z/2 | `[14, 15]` | UNSAT at r=14 | SAT (Phase 1.6, encoder ready) |
| Confirm `R_{F₂}(⟨3,3,3⟩) = 23` | `⟨3,3,3⟩` | Z/2 | `[20, 23]` after [wang26][wang26] | UNSAT at r=22 | SAT + embedding cubes; build on Wang's orbit-DP technique |
| Tighten `R(⟨3,3,3⟩)` over `R` | `⟨3,3,3⟩` | R / Q | `[19, 23]` | improve LB toward 23 | reproduce Wang's method over larger fields |
| Mixed composition for `⟨32,32,32⟩` over `R` | `⟨32,32,32⟩` | R | 16,807 (Strassen⁵) | < 16,807 | catalog-driven search |

## 8. References

- <a name="strassen69"></a>[strassen69] V. Strassen, *Gaussian elimination is not optimal*,
  Numerische Mathematik 13 (1969), 354–356.
- <a name="hk71"></a>[hk71] J. E. Hopcroft, L. R. Kerr, *On minimizing the number of
  multiplications necessary for matrix multiplication*, SIAM J. Appl. Math. 20
  (1971), 30–36.
- <a name="win71"></a>[winograd71] S. Winograd, *On multiplication of 2×2 matrices*,
  Linear Algebra Appl. 4 (1971), 381–388.
- <a name="lad76"></a>[laderman76] J. D. Laderman, *A noncommutative algorithm for
  multiplying (3 × 3) matrices using 23 multiplications*, Bull. AMS 82 (1976), 126–128.
- <a name="pan78"></a>[pan78] V. Y. Pan, *Strassen's algorithm is not optimal*, FOCS 1978.
- <a name="blaser03"></a>[blaser03] M. Bläser, *On the complexity of the multiplication
  of matrices of small formats*, J. Complexity 19 (2003), 43–60.
- <a name="smirnov"></a>[smirnov] A. V. Smirnov, *Bilinear complexity and practical
  algorithms for matrix multiplication*, Comput. Math. Math. Phys. 53 (2013).
- <a name="drisc09"></a>[drisc09] C. Drevet, M. Islam, É. Schost, *Optimization
  techniques for small matrix multiplication*, 2009. [PDF](https://cs.uwaterloo.ca/~eschost/publications/DrIsSc09.pdf).
- <a name="wang26"></a>[wang26] C. Wang, *Automated Lower Bounds for Small Matrix Multiplication
  Complexity over Finite Fields*, arXiv:2603.07280 (2026).
  [arXiv](https://arxiv.org/abs/2603.07280). Proves `R_{F₂}(⟨3,3,3⟩) ≥ 20`
  via orbit classification + dynamic programming + verifiable proof
  certificates — methodology directly relevant to this codebase's SAT pipeline.
- <a name="alphatensor"></a>[alphatensor] A. Fawzi et al. (DeepMind),
  *Discovering faster matrix multiplication algorithms with reinforcement
  learning*, Nature 610 (2022). [Article](https://www.nature.com/articles/s41586-022-05172-4),
  [code+data](https://github.com/deepmind/alphatensor).
- <a name="alphaevolve"></a>[alphaevolve] A. Novikov et al. (DeepMind),
  *AlphaEvolve: A coding agent for scientific and algorithmic discovery*,
  arXiv:2506.13131 (2025). [arXiv](https://arxiv.org/abs/2506.13131).
- <a name="perminov"></a>[perminov] Perminov, Andrew I., *FastMatrixMultiplication*,
  [GitHub](https://github.com/dronperminov/FastMatrixMultiplication#research-findings--status).
  GitHub username `dronperminov`; cite as "Perminov". Papers behind the repo:
  arXiv:2603.02398, arXiv:2511.20317, arXiv:2512.13365, arXiv:2512.21980,
  arXiv:2606.02480 (REFERENCES.md [75], [77], [78], [79], [86]).

[strassen69]: #strassen69
[hk71]: #hk71
[win71]: #win71
[lad76]: #lad76
[pan78]: #pan78
[blaser03]: #blaser03
[smirnov]: #smirnov
[drisc09]: #drisc09
[alphatensor]: #alphatensor
[alphaevolve]: #alphaevolve
[wang26]: #wang26
[perminov]: #perminov
[solving_strategies]: SOLVING_STRATEGIES.md
[z3als]: src/main/java/io/cormoran/strassen/v3/Z3Als.java

---

## 9. Working notes

- Verify and replace `~`-marked entries in §4.6 against the
  [perminov][perminov] structured DB before publishing.
- `⟨3,3,3⟩` over `R`/`Q`: the `≥ 19` lower bound is [bläser03][blaser03].
  Verify the current state of the `Z/2`-specific bound.
- AlphaTensor improved several non-cubic formats over Z/2 that
  aren't captured in §5 yet; pull from [alphatensor][alphatensor]'s
  supplementary tables.
- Pan 1978's "border-rank" results and Bini 1979 are documented in
  [RANK_KNOWLEDGE.md](RANK_KNOWLEDGE.md) §1.1 (border rank column) but the *practical*
  border-rank-derived algorithms (e.g. Pan's 1980 trilinear aggregation)
  are not catalogued here yet.
