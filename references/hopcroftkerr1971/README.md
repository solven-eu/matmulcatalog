# Hopcroft–Kerr 1971: an impossibility theorem for the `(2,2,bridge-3)` same-method pair

> ## ⚑ STATUS (2026-06-12): the constructive closure HAPPENED — this is now a historical document
>
> The goal stated below ("a constructive way to attain the HK bound") was
> achieved in full: **465 integer schemes** over `3 ≤ p ≤ 32, p ≤ n ≤ 32`,
> 459 at the exact formula, dual-certified, registered in
> `schemes/constructed/`. The gap table below CLOSED — `⟨2,10,15⟩=233`,
> `⟨2,10,16⟩=248`, `⟨2,12,16⟩=296`, all over ℤ, all strictly below the
> catalogs. See `research/hopcroft-kerr-2np/CONSTRUCTIVE_METHOD.md` (the
> living recipe + proofs) and `paper/sections/hk71.tex`.
>
> **Standing of this document's theorem**: correct as stated — and twice
> superseded in practice. (1) *Sidestepped*: arc-interior bridge SELECTION
> makes the `(2,2,bridge-3)` configuration never arise in the schedule.
> (2) *Reversed operationally*: over the TRUE 12-product reusable set (the
> theorem's S has only 6), `(2,2,b3)` and `(1,1,b3)` ARE 3-product-solvable
> (`sympy/derive_bridge_true_reusables.py`, identities implemented as
> emitters). The genuinely impossible case moved to **`(3,3,bridge-1/2)`**,
> now a theorem over ARBITRARY local rank-1 atoms for the 9-product set
> (`sympy/derive_33bridge_general.py`) — it gates the six remaining
> `g ≥ 6` shapes at +1..+3 (task #9, the only open item).
>
> The text below is kept verbatim as the historical record of the
> impossibility analysis (its mathematics is reused by the newer proofs).

## What we're looking for (synthetic) — *historical framing, resolved*

A **constructive way** to attain the Hopcroft–Kerr 1971 upper bound
`R(⟨2, n, p⟩) ≤ ⌈(3np + max(n,p))/2⌉` at arbitrary `(n, p)`.

HK's page-10 same-method-pair derivation works through ONE case explicitly
and dismisses five others "by symmetry"; four close by direct algebra,
the fifth — `(2, 2, bridge-3)` — has resisted, and we have a
Gröbner-basis impossibility theorem (reproducible in `sympy/` +
`mathematica/`) ruling out the natural rank-1 atom family over
characteristic 0.

**The published catalogs are not a full substitute.** Sedoglavic et al.
(FMM-Lille) and Perminov publish explicit factor matrices for many
specific `⟨2, n, p⟩`, and several attain the HK formula, but coverage is
incomplete — concrete shapes where even Perminov sits **strictly above**
the formula (all three closed by our construction, 2026-06):

| Shape  | HK formula | Best published | Gap | Ours (2026-06) |
|--------|-----------:|---------------:|----:|---------------:|
| ⟨2,10,15⟩ | 233 | 234 | +1 | **233 (ℤ)** |
| ⟨2,10,16⟩ | 248 | 249 | +1 | **248 (ℤ)** |
| ⟨2,12,16⟩ | 296 | 298 | +2 | **296 (ℤ)** |

A constructive closure **generates** an attaining scheme at any
`(n, p)` rather than relying on case-by-case imports — it composes
inside larger Kronecker constructions, generalises to `n, p` beyond what
any catalog publishes, and provides algorithmic provenance for the
existing opaque factor matrices. (All delivered; see the status banner.)

The remainder of this document develops the impossibility theorem in
the rank-1 / characteristic-0 case.

---

A self-contained writeup of an open derivation in Hopcroft & Kerr 1971,
"On Minimizing the Number of Multiplications Necessary for Matrix
Multiplication" (*SIAM J. Appl. Math.* 20(1), pp. 30–36), and a
**characteristic-0 impossibility theorem** ruling out one natural
candidate solution. Self-contained for review by a peer or external
solver agent; all scripts + outputs are in `sympy/`.

## Abstract

HK 1971 derives `R(⟨2,n,n⟩) ≤ (3n²+2n)/2` for odd `n ≥ 3` by tiling the
symmetric `n × n` output matrix with "pair products" (Lemma 2) plus a
**same-method-pair** correction (page 10). The page-10 derivation
covers the `(1,1,bridge-2)` case; the paper writes that "the other 5
cases follow by symmetry" and never derives them. Of the remaining
five, three follow by sympy-verified algebra
(`(2,2,bridge-1)`, `(1,1,bridge-3)`, plus two unused cases); a fourth,
**`(2,2,bridge-3)`**, has resisted enumeration over standard
shift-form atom catalogs.

We prove the following:

> **Theorem (informal).** No 3-product `(2,2,bridge-3)` derivation exists
> within the rank-1 atom family `{ L(a)·M(x) : L, M ∈ R^6 }` in
> characteristic 0.

Concretely: the would-be three new multiplications cannot be chosen as
rank-1 tensors `L_k(a)·M_k(x)` (with `L_k, M_k` arbitrary linear forms
in the available `a`-row and `x`-row scalars). Whether the
`(2,2,bridge-3)` derivation can be closed with multilinear atoms of
higher rank, with characteristic-2 atoms, or via an entirely different
framing is **open**.

A targeted literature search (see "Prior art and originality claim"
below) suggests the result is new, pending verification of four
hard-to-access sources (Pan 1984 LNCS 179, Heun 1994, De Groote 1983,
Smirnov 2019). Strong empirical corroboration: FMM-Lille publishes
every `⟨2, n, n⟩` for n = 9..16 strictly above the HK formula by
+1 to +5, all credited to Hopcroft–Kerr 1971 — i.e. the community
has never published a construction that *attains* the formula at the
shapes our impossibility predicts cannot be attained.

Practical consequence *(historical — every claim in this paragraph was
subsequently resolved by the constructive closure; see the status
banner)*: HK Case 1 (odd `p`) reached its claimed `(3pn+n)/2` bound only
for `p ∈ {3, 5}` via the then-existing emitters; for `p ≥ 7` and for
Case 2 (even `p`) Steps 1–3, the fallback `buildNaiveDCE` cost `+δ`
extra mults per bridge-3 pair. The gap shapes `⟨2, 10, 15⟩`,
`⟨2, 10, 16⟩`, `⟨2, 12, 16⟩` remained `1–4` mults above the HK formula
in our catalog at the time, as in FMM-Lille and Perminov — consistent
with the narrow-S result above, until bridge selection + the
true-reusable identities dissolved it.

## Prior art and originality claim

A targeted literature search (sources listed below) found **no paper
that revisits HK 1971's page-10 same-method-pair derivation at the
level of detail required to close or to fail to close the
`(2,2,bridge-3)` case**. The terminology used here ("same-method
pair", "bridge product", "rank-1 atom") is repo-internal; no external
paper, to our knowledge, uses this framing.

**The result claimed here therefore appears to be new** — pending
verification of three sources we could not fully access (see "Sources
not fully checked" below). Even granting a hidden prior closure in one
of those, the result is at minimum a self-contained, sympy- +
Mathematica-verifiable, constructively-provable impossibility in
explicit Gröbner-basis form, which is independently useful.

### Strong corroborating evidence

The most compelling indirect witness for the impossibility is the
state of the existing catalogs:

| Shape | HK formula | Catalog SOTA | Excess | Source |
|-------|-----------:|-------------:|-------:|--------|
| `⟨2, 9, 9⟩` | 126 | 127 | +1 | FMM-Lille (hopcroft:1971-credited) |
| `⟨2, 10, 10⟩` | 155 | 156 | +1 | FMM-Lille |
| `⟨2, 10, 15⟩` | 232 | 234 | +2 | FMM-Lille |
| `⟨2, 10, 16⟩` | 247 | 249 | +2 | FMM-Lille |
| `⟨2, 12, 16⟩` | 294 | 298 | +4 | FMM-Lille |
| `⟨2, 14, 14⟩` | 301 | 303 | +2 | FMM-Lille |
| `⟨2, 14, 15⟩` | 324 | 327 | +3 | FMM-Lille |
| `⟨2, 14, 16⟩` | 346 | 350 | +4 | FMM-Lille |
| `⟨2, 16, 16⟩` | 392 | 397 | +5 | FMM-Lille |

Every entry credits Hopcroft–Kerr 1971 for the algorithm idea, yet
sits **strictly above** the HK formula. If a constructive closure of
the same-method-pair gap were known, these entries would attain the
formula. The +1 to +5 excess matches exactly the `+δ` slack a
`buildNaiveDCE`-style fallback would cost at each unclosed bridge-3
pair. See `references/fmm-lille-discrepancies.md` for the full
discrepancy listing.

### Sources consulted (literature-search verdict table)

| Source | Year | Verdict | Notes |
|--------|-----:|---------|-------|
| **Hopcroft & Kerr**, *SIAM J. Appl. Math.* 20(1) | 1971 | Source of the gap | Page 10: "the other cases follow by symmetry" — never derived. |
| **Hopcroft & Musinski**, draft / *SIAM J. Comput.* 2(3) | 1973 | Adjacent but distinct | Reframes HK via bilinear duality; does not re-derive same-method-pair cases. |
| **Probert**, *SIAM J. Comput.* 5(2) | 1976 | Doesn't address | Additive complexity of ⟨2,2,2⟩, different setting. |
| **Pan**, FOCS proceedings | 1978 | Adjacent but distinct | Trilinear aggregation; asymptotic-complexity focus. |
| **Pan**, *How to Multiply Matrices Faster* (LNCS 179) | 1984 | **Unresolved** (paywalled) | Most-likely candidate for hidden coverage of HK's finite-n details; we could not retrieve the full text. |
| **De Groote**, *SIAM J. Comput.* 12(1):101–117 | 1983 | **Closest structural** | Treats `⟨n,2,n⟩` via division-algebra orbit theorems, but only for *minimal-rank* algorithms (rank `m+n-1`); HK's `(3n²+2n)/2` construction is strictly *over-rank* so De Groote's theorems do not directly imply our impossibility. **Worth re-reading once retrievable**. |
| **Heun**, *Fast komplexe Matrixmultiplikation* (German) | 1994 | **Unresolved** | Could not retrieve online. Scope (asymptotic/aggregation) makes hidden coverage unlikely. |
| **Bürgisser, Clausen, Shokrollahi**, *Algebraic Complexity Theory* (Grundlehren 315) | 1997 | Doesn't address | Cites HK 1971 for the `R(⟨2,2,2⟩) ≥ 7` lower bound and the upper-bound formula; no derivation audit. |
| **Smirnov**, *Comput. Math. Math. Phys.* 53 | 2013 | Doesn't address | Numerical / ALS catalog; treats HK 1971 as black-box bound. |
| **Alekseev & Smirnov family** | 2013–2023 | Adjacent but distinct | Exact bilinear complexity of `⟨2,2,m⟩`, `⟨2,m,2⟩`, `⟨2×7,7×2⟩`, etc. Same general program (small-format bilinear complexity) but never `⟨2,n,n⟩` for n ≥ 3 or HK's bridge construction. |
| **Smirnov**, *Moscow Univ. CMC* | 2019 | **Unresolved** (paywalled) | Title suggests `⟨2,2,m⟩` over finite fields. Worth retrieving to check whether his analysis uses a HK-style bridge. |
| **Conner–Harper–Landsberg** (border-rank LB on `⟨n,m,m⟩`) | 2020s | Adjacent but distinct | Lower-bound side; doesn't address upper-bound constructions. |
| **Kauers & Moosbauer**, flip-graph search | 2023 | Doesn't address | Heuristic / SAT rank improvements; no structural analysis of HK's derivation. |
| **Kauers & Wood**, meta-flip-graph | 2025 | Doesn't address | Same. |
| **Heule et al.**, *Fast Matrix Multiplication in Small Formats* | 2026 | Doesn't address | Same. |
| **AlphaTensor** (DeepMind, Fawzi et al.), *Nature* | 2022 | Doesn't address | RL search for small-format ranks; does not address derivation structure. |
| **AlphaEvolve** (DeepMind, follow-up) | 2024 | Doesn't address | Same. |
| **FMM-Lille catalog** (Sedoglavic et al., ongoing) | — | **Empirical witness** | Lists every `⟨2,n,n⟩` for n ≥ 9 strictly above HK formula (table above), all credited "hopcroft:1971". |
| **Perminov catalog** (ongoing) | — | **Empirical witness** | Same pattern. |
| **DIS09 archive** (DeGroote, Ibarra, Sahni 2009 republication) | 2009 | Doesn't address | Catalog of small-format algorithms; doesn't audit HK's derivation. |

### Sources not fully checked

Three sources where a hidden prior closure remains possible. Verifying
each would strengthen the originality claim:

1. **Pan 1984 monograph** "How to Multiply Matrices Faster", LNCS 179.
   Likely chapter to check: the section on Hopcroft–Kerr's `⟨2,n,n⟩`
   construction. Not retrievable via free web search; needs library
   or Springer access.
2. **Heun 1994** "Fast komplexe Matrixmultiplikation" (German,
   bilinear-complexity textbook). Not findable online; likely needs
   German library access.
3. **Smirnov 2019** *Moscow Univ. CMC* 4: bilinear complexity of
   `⟨2, 2, m⟩` over finite fields. Springer-paywalled. Possibly uses
   HK-style bridge in the analysis.
4. **De Groote 1983** *SIAM J. Comput.* 12(1):101–117 — paywalled.
   Closest structural result; would need a full re-read to confirm
   the orbit theorems don't imply our impossibility as a special case
   (we believe they don't, because HK's construction is over-rank,
   but verification by direct reading is warranted).

**Reader request**: if you have access to any of the above and can
either confirm a prior closure or confirm "doesn't address", please
open a PR against this section.

### Originality claim (precise)

**Pending the 4 unresolved sources, the following is new to our
literature search**:

> (i) An explicit Gröbner-basis-verifiable, characteristic-0
> impossibility result for the same-method-pair `(2, 2, bridge-3)`
> case of Hopcroft & Kerr 1971's Lemma 2, ruling out 3-product
> derivations within the rank-1 atom family
> `{ L(a)·M(x) : L, M ∈ R^6 }`.
>
> (ii) An empirical decomposition of the would-be solution variety
> `X = { (A, A', c) : g(A') = c·g(A) }` into exactly two irreducible
> components — the proportional locus `P̄` and the zero-garbage
> slice `R¹² × V × {0}` — proved via Gröbner saturation
> `(I : c^∞) = I_P`.
>
> (iii) Three independent algebraic routes (tangent-space, Gröbner
> saturation, linearisation + minor-cuts) reaching the same
> conclusion, plus a Mathematica cross-CAS verification scaffold
> (`mathematica/01_independent_saturation.wls`).

What the result does **not** claim:
- A lower bound on `R(⟨2, n, n⟩)`. The HK upper bound stands; our
  result only says HK's *specific construction* has a gap.
- Originality of the framework. The "atom + quotient" / residual-rank
  analysis is standard bilinear-complexity machinery.
- Coverage outside rank-1 atoms over characteristic 0 (see "What's
  NOT ruled out" below).

## How to reproduce

Dependencies: Python 3.8+, sympy 1.10+. No project imports; each
script under `sympy/` runs as `python3 sympy/<script>.py`.

```bash
# Smoke-check: HK's published (1,1,bridge-2) derivation
python3 sympy/page10_22b1.py             # ~1 s

# Establish quotient: dim Q = 2, explicit basis q1, q2
python3 sympy/quotient_analysis_v2.py    # ~5 s

# Show no single rank-1 atom has zero garbage AND non-zero (α, β)
python3 sympy/zero_garbage_variety.py    # ~10 s

# Impossibility proof: dim P = dim X = 14, so atoms with proportional
# garbage are proportional rank-1 tensors
python3 sympy/paired_garbage_analysis_v2.py   # ~11 s
```

Expected outputs are checked in alongside (`.out` siblings). Re-running
should reproduce them up to non-deterministic random-point picks (the
Jacobian-rank checks pass at every numeric point we have tried; the
script picks 5 random proportional points to guard against degenerate
samples).

## Setup

### Notation

Write the `⟨2, n, n⟩` matmul as `Y = X · A` with `X ∈ K^{2 × n}` and
`A ∈ K^{n × n}`. HK collapses each output-pair `(y_{r,t}, y_{r,t'})`
into a 6-scalar bilinear problem with row-index variables
`a_*_1, a_*_2, x_*_1, x_*_2` for `* ∈ {i, p, j}` (three rows of `A`
selected for the pair construction). The 36 monomials `a_α_β · x_α'_β'`
span the rank-1 tensor space inside which all atoms live:

```
a-vars: a_i_1, a_i_2, a_p_1, a_p_2, a_j_1, a_j_2
x-vars: x_i_1, x_i_2, x_p_1, x_p_2, x_j_1, x_j_2
```

### The six base methods

HK's per-position constructors:

| Name | Formula |
|------|---------|
| `A(a_1,a_2,x_1,x_2)` | `a_2 · (x_1 + x_2)` |
| `B(a_1,a_2,x_1,x_2)` | `(a_1 − a_2) · x_1` |
| `C(a_1,a_2,x_1,x_2)` | `(a_1 − a_2) · x_2` |
| `D(a_1,a_2,x_1,x_2)` | `a_1 · (x_1 + x_2)` |
| `E(a_1,a_2,x_1,x_2)` | `a_2 · x_2` |
| `F(a_1,a_2,x_1,x_2)` | `a_1 · x_1` |

"Method 1" uses `{A, B}`; "method 2" uses `{C, D}`; "method 3" uses
`{E, F}`. Lemma 2 shows that when two adjacent diagonal positions use
different methods, the bridge position can be computed with 3 (not 4)
new products — saving one mult per pair.

### The (2, 2, bridge-3) setup

Both `y_{ii}` and `y_{jj}` are assigned method 2 (so `C_{ii}, D_{ii},
C_{jj}, D_{jj}` are computed and reusable). The bridge position
`y_{i+1, i+1}` uses method 3 (`E_{i+1,i+1}, F_{i+1,i+1}` reusable). The
"shared" set is

```
S = { Ci, Di, Cj, Dj, Ep, Fp }
where  Ci=(a_i_1−a_i_2)·x_i_2     Di=a_i_1·(x_i_1+x_i_2)
       Cj=(a_j_1−a_j_2)·x_j_2     Dj=a_j_1·(x_j_1+x_j_2)
       Ep=a_p_2·x_p_2              Fp=a_p_1·x_p_1     (using p ≡ i+1)
```

The targets are

```
T1 = y_{ij} = a_i_2·x_i_1 + a_p_2·x_p_1 + a_j_2·x_j_1
T2 = y_{ji} = a_i_1·x_i_2 + a_p_1·x_p_2 + a_j_1·x_j_2
```

and HK's framework allots **3 new products** to compute `T1` and `T2`
given the reusables `S`. Each new product is a rank-1 tensor
`L(a) · M(x)` with arbitrary linear forms `L, M` over the 6+6 scalars.

## The proof (outline)

The argument has three steps, each with a sympy harness.

### Step 1: The quotient `Q` has dimension 2

`Q := span(T1, T2, S) / S` is the algebraic space the new products must
cover. By direct calculation:

```
rank(S) = 6
rank(S, T1, T2) = 8
⇒ dim Q = 8 − 6 = 2
```

Explicit basis (Euclidean orthogonal projection onto S⊥):

```
q1 = P_⊥ T1 = a_i_2·x_i_1 + a_p_2·x_p_1 + a_j_2·x_j_1
q2 = P_⊥ T2 = −⅓ a_i_1·x_i_1 + ⅓ a_i_1·x_i_2 + ⅓ a_i_2·x_i_2
              + a_p_1·x_p_2
              − ⅓ a_j_1·x_j_1 + ⅓ a_j_1·x_j_2 + ⅓ a_j_2·x_j_2
```

The fractional coefficients are an artifact of the projection (a
non-canonical basis choice). What is canonical is the *coset*: `T1 ↦
q1` and `T2 ↦ q2` in `Q`. Computed in `sympy/quotient_analysis_v2.py`.

### Step 2: Decomposition of a rank-1 atom

Write `L = Σ λ_k a_k`, `M = Σ μ_k x_k`. The atom `L·M` projects to S⊥
with two structural components:

- **Q-component**: `(α(λ,μ), β(λ,μ))` in `Q ≅ R²`. Closed-form
  - `α = (l_2 m_1 + l_4 m_3 + l_6 m_5) / 3`
  - `β = (−l_1 m_1 + l_1 m_2 + l_2 m_2 + 3 l_3 m_4 − l_5 m_5 + l_5 m_6 + l_6 m_6) / 5`
- **garbage**: a 28-D component in `S⊥ ∩ Q⊥`, generically non-zero.

A 3-product `(2, 2, bridge-3)` derivation requires choosing 3 atoms
`P_k = L_k·M_k` and output weights `W_1, W_2 ∈ R^3` such that

```
W_1 · (Q-components)  = (1, 0)   AND   W_1 · (garbage_k) = 0       (yields T1 mod S)
W_2 · (Q-components)  = (0, 1)   AND   W_2 · (garbage_k) = 0       (yields T2 mod S)
```

So `W_1, W_2` must lie in the left null space of the 28×3 garbage
matrix `U = [g(A_1) | g(A_2) | g(A_3)]`. This null space has dimension
`3 − rank(U)`; for the system to admit two linearly independent
solutions, we need `rank(U) ≤ 1`.

### Step 3a: No single atom has zero garbage with non-zero `(α, β)`

The variety `V = { (λ, μ) : g(λ, μ) = 0 }` decomposes into 8
axis-aligned affine components (one per minimal vertex cover of the
bipartite incidence graph `{l_1,l_2}×{m_3..m_6}, {l_3,l_4}×{m_1,m_2,m_5,m_6},
{l_5,l_6}×{m_1..m_4}` formed by the 24 monomial garbage equations like
`l_1·m_3 = 0`). Five components are "axis" (some of `l_*, m_*` forced
to 0 such that none of the α, β monomial-terms survives); three are
"block-pair" components where the residual non-monomial equations
reduce `α, β` to 0 in the Gröbner basis. **In all 8 components, the
`(α, β)`-image is the single point `(0, 0)`.**

So if we try to build a solution from atoms with zero garbage, all
three atoms contribute `(0, 0)` to Q ⇒ cannot span Q. Computed in
`sympy/zero_garbage_variety.py`.

### Step 3b: Atoms with proportional garbage are proportional tensors

Define

```
X = { (A, A', c) ∈ R^12 × R^12 × R : g(A') = c · g(A) }
P = { (A, A', c) ∈ X : (L', M') = (α·L, β·M) for some α, β ∈ R, c = α·β }
```

`P ⊆ X` is the "proportional locus" — atoms `A'` that are scalar
multiples of `A` (in the L⊗M sense).

- **`dim P = 14`**: the polynomial map
  `ψ : R^14 → R^25, (λ, μ, α, β) ↦ (λ, μ, α·λ, β·μ, α·β)` is
  generically injective with rank-14 Jacobian.
- **`dim X = 14`**: at five random proportional points the symbolic
  Jacobian of the constraint map
  `Φ : R^25 → R^36, (A, A', c) ↦ g(A') − c·g(A)` has uniform rank 11,
  so `dim X = 25 − 11 = 14` on the smooth locus through `P`.

Computed in `sympy/paired_garbage_analysis_v2.py`.

### The shape of `X`: `P` plus a `c = 0` family

A subsequent Gröbner-basis analysis (`sympy/x_equals_p_proof.py`)
shows that `X` is **not** exhausted by `P`. The fibre ideal over a
generic anchor `A` strictly contains the proportional ideal — in
particular, equations like `13·L'_1 − 2·L'_6 = 0` appear in the fibre
basis but not in the proportional ideal. The extra solutions concentrate
on the locus `c = 0`:

```
X ⊇ X_{c=0} := { (A, A', 0) ∈ X } = { (A, A', 0) : g(A') = 0 } = R^12 × V × {0}
```

where `V` is the zero-garbage variety from Step 3a.

**Saturation closes the gap.** The same script (`x_equals_p_proof.py`,
Step D) computes `I : c^∞` via Rabinowitsch (adjoin `1 − t·c`,
eliminate `t`). The saturated ideal has 11 generators and reduces to
the proportional ideal `I_P` **both ways**: `(I : c^∞) = I_P`. So as
algebraic sets,

```
V(I)  =  V(I_P)  ∪  V(I + ⟨c⟩)  =  P-fibre  ∪  ({λ} × V × {0}).
```

There are **no other components**. Substituting `c = 0` into the 34
fibre equations reproduces exactly the 34 garbage-vanishing equations
of `zero_garbage_variety.py` (24 pure monomials `L_i M_j = 0` plus 10
non-monomial residuals), which is the literal defining ideal of `V`.
Hence the full set-theoretic decomposition

```
X  =  P̄  ∪  X_{c=0}      (over any A with all λ_k, μ_j ≠ 0)
```

is **proved**, not merely conjectured.

**The impossibility theorem still holds**, but the argument is slightly
more subtle than `P = X`. Three atoms `(P_1, P_2, P_3)` with pairwise
rank-1 garbage matrix split into two cases:

1. **All three are in `P`** (i.e. pairwise proportional rank-1
   tensors). Their `(α, β)` projections lie on a single line in
   `Q ≅ R²`. Cannot span both `q1` and `q2`.
2. **At least one atom is in `X_{c=0}`** (i.e. has `g(A_k) = 0`,
   so `A_k ∈ V`). By Step 3a, such an atom has `(α, β) = (0, 0)`
   — it contributes nothing to `Q`. The remaining ≤ 2 atoms must
   span all of `Q`, but they are still in the rank-1 family and
   are pairwise-proportional (or one of them is also in `V`),
   so their projections also lie on a single line through the
   origin — cannot span 2-D `Q`.

In both cases the 3-atom team fails to cover `Q`, so the construction
cannot be completed. **No 3-product `(2, 2, bridge-3)` derivation
exists in the rank-1 atom family over a field of characteristic 0.**

No further components: the saturation `(I : c^∞) = I_P` plus the
identification `V(I + ⟨c⟩) = {λ} × V × {0}` together exhaust `V(I)`.
The impossibility theorem for `(2, 2, bridge-3)` in the rank-1 atom
family over characteristic-0 fields is therefore **unconditional** —
no smooth-locus, generic-A, or "possibly more components" caveats
remain.

## Implications for HK 1971's main theorem

HK's main result `R(⟨2,n,n⟩) ≤ (3n²+2n)/2` for odd `n ≥ 3` stands, but
**its proof has a gap** at the `(2,2,bridge-3)` same-method-pair
configuration. The bound itself is correct: the gap costs at most `+1`
mult per affected pair, and the FMM-Lille / Perminov publications for
the affected shapes (`⟨2,10,15⟩`, `⟨2,10,16⟩`, `⟨2,12,16⟩`) sit
`1–4` above the HK formula, consistent with this slack.

The impossibility theorem above says the gap **cannot be closed
within the rank-1 atom family in characteristic 0**. To close it (and
recover HK's claimed exact rank), one of the following must hold:

- A non-rank-1 multilinear atom (sum of two or more rank-1 tensors,
  treated as one "product") suffices. Unusual but not ruled out.
- A characteristic-2 (or other small-prime) construction exists that
  doesn't lift to characteristic 0. HK was working over fields of
  characteristic 0; an F_2-specific closure would be a refinement, not
  a correction.
- An entirely different decomposition of `T1, T2` modulo `S` — e.g.
  with different sharing of the bridge — was implicitly intended.

The narrative in HK page 10 ("the other cases follow by symmetry")
appears to be wrong, at least at face value. The natural involutions
on the variables (`(a_1↔a_2)`, `(x_1↔x_2)`) swap method 1 ↔ method 2
but fix method 3, so the `(1, 1, bridge-2)` derivation does not map
cleanly onto `(2, 2, bridge-3)` under any single involution. The
sympy harnesses confirm that no "naive symmetry" reduction works.

## What's NOT ruled out

The impossibility theorem applies only to rank-1 atoms over
characteristic 0. Open:

1. **Higher-rank atoms.** Each "product" in HK's framework is one
   scalar multiplication, so each atom should be rank-1; but the
   bilinear-complexity literature has explored "approximate" rank-1
   atoms (border rank). Not investigated here.
2. **Characteristic 2 (or other small primes).** Some bilinear
   constructions exist over F_2 that don't lift to characteristic 0
   (e.g. AlphaTensor's `⟨4,4,4⟩ = 47` over F_2 vs `48` over C). A
   careful check whether the `(2,2,bridge-3)` system has an F_2
   solution is a clean follow-up.
3. **Reorganising HK's pair decomposition.** The whole "same-method
   pair" framework is HK's specific tiling of the symmetric output
   matrix. A different tiling might avoid the bridge-3 case entirely.
4. **A different residual subspace.** Our `S` reuses the specific
   methods HK assigned. With different method assignments for
   `y_{ii}, y_{jj}, y_{i+1,i+1}`, the residual quotient `Q` changes
   and the impossibility doesn't transfer.

A peer with a different framing or solver may close some of these. The
sympy `verify.py` harness accepts any 3-atom candidate and verifies
it against `T1, T2 mod S` within seconds.

## Sympy harnesses (this directory)

All scripts live under [`sympy/`](sympy/). No project imports.

| Script | Purpose |
|--------|---------|
| [`sympy/quotient_analysis_v2.py`](sympy/quotient_analysis_v2.py) | Quotient `Q` setup, explicit `q1, q2`, generic `(α, β)` formula. (v1 had a conflation bug — see `quotient_analysis.py` for the audit trail.) |
| [`sympy/zero_garbage_variety.py`](sympy/zero_garbage_variety.py) | Step 3a: variety `V` of atoms with zero garbage, prime decomposition into 8 components, Gröbner reduction showing `(α, β)` is identically 0 on each. |
| [`sympy/paired_garbage_analysis_v2.py`](sympy/paired_garbage_analysis_v2.py) | Step 3b: `dim P = dim X = 14`. Verified via Jacobian rank at 5 random proportional points. (v1 had a linearisation bug giving `dim X = 13` — kept as audit trail.) |
| [`sympy/x_equals_p_proof.py`](sympy/x_equals_p_proof.py) | Gröbner-basis tightening: shows `X = P̄ ∪ X_{c=0}` literally, by proving `(I : c^∞) = I_P` (saturation kills the `c = 0` component) and `V(I + ⟨c⟩) = V`. Closes the impossibility theorem unconditionally. |
| [`sympy/analyze_tensor_rank.py`](sympy/analyze_tensor_rank.py) | Tensor-rank decomposition of the residual subspace — the framing that exposed the early "joint matrix rank" mistake. |
| [`sympy/page10_22b1.py`](sympy/page10_22b1.py) | Sanity check of HK's published `(1, 1, bridge-2)` derivation. |
| [`sympy/sanity_11bridge2.py`](sympy/sanity_11bridge2.py) | Independent sympy verification of the `(1, 1, bridge-2)` emitter in our Java code. |
| [`sympy/derive_22bridge1.py`](sympy/derive_22bridge1.py) | Successful sympy derivation of `(2, 2, bridge-1)` — template for attacking other cases. |
| [`sympy/derive_all.py`](sympy/derive_all.py) | Sweep over all 6 same-method-pair cases. |
| [`sympy/same_method_pair.py`](sympy/same_method_pair.py) | Skeleton enumerator factored from the case-specific scripts. |
| [`sympy/verify.py`](sympy/verify.py) | Generic verifier: take a candidate 3-atom solution as input, check it computes `T1, T2` exactly. |
| [`sympy/search_22bridge3_v4.py`](sympy/search_22bridge3_v4.py) | The brute-force enumeration that, in retrospect, was searching for something that doesn't exist (per the impossibility theorem). Kept as audit trail. |
| [`sympy/search_22bridge3_v5.py`](sympy/search_22bridge3_v5.py) | Same. |
| [`sympy/quotient_analysis.py`](sympy/quotient_analysis.py) | v1 (buggy) of quotient analysis. Kept as audit trail. |
| [`sympy/paired_garbage_analysis.py`](sympy/paired_garbage_analysis.py) | Independent re-derivation of the v2 result via a different route: linearise the 28 garbage-difference equations in `t_kj := L'_k·M'_j`, nullspace-solve at a generic A (kernel dim 9), then impose `rank(T) ≤ 1` via 2×2 minors and a Gröbner basis. Confirms `c≠0` ⇒ proportional (tangent dim 1 at proportional point, in 5/5 random trials); the `c=0` branches recover the `V` variety. Cross-check for `paired_garbage_analysis_v2.py`. |
| `*.out` | Outputs of the corresponding scripts; commit-friendly. |

Audit-trail policy: failed attempts are kept under their version
numbers so a later reader can see what was tried and why it didn't
work. The current "use this" scripts have no `vN` suffix or have
`_v2` to mark them as the canonical version. Future revisions should
keep this convention.

## Surrounding code (Java)

The HK constructor lives at:

* `src/main/java/eu/solven/matmul/papers/hopcroftkerr1971/HopcroftKerr2bc.java`
* `src/main/java/eu/solven/matmul/papers/hopcroftkerr1971/HopcroftKerr2bcAsymmetric.java`
* `src/main/java/eu/solven/matmul/papers/hopcroftkerr1971/LemmaOneAugmentation.java`

The same-method-pair dispatcher is `emitSameMethodPair` at
`HopcroftKerr2bc.java:405`. The three implemented cases are
`emitSameMethodPair_11_bridge2`, `_22_bridge1`, `_11_bridge3`. The
open case `emitSameMethodPair_22_bridge3` currently dispatches to the
naïve fallback (`buildNaiveDCE`) — given the impossibility result, no
3-product implementation can replace it within the rank-1 family.

Java-side verification of any candidate construction uses
`Verifier.passesRandomMatmulSpotCheckNB` (fast spot check) or
`Verifier.residualNonBilinear` (full tensor identity).

## Literature search — status

A targeted literature search was conducted (see "Prior art and
originality claim" above for the full verdict table). The result
appears to be **new** pending verification of four hard-to-access
sources, listed there. Help wanted: if you have institutional access
to any of Pan 1984 (LNCS 179), Heun 1994, De Groote 1983, or Smirnov
2019, please confirm whether they address the `(2, 2, bridge-3)`
same-method-pair case, and open a PR against the prior-art section.

## License + contribution notes

This note and the sympy harnesses are released under the repository's
top-level licence (see `../../LICENSE`). To propose:

- A constructive solution (in some atom family the impossibility
  theorem doesn't cover): open a PR adding the construction to
  `sympy/` and the corresponding Java emitter at
  `HopcroftKerr2bc.emitSameMethodPair_22_bridge3`.
- An independent verification of the impossibility proof (Gröbner
  bases re-run in a different CAS, or an algebraic-geometry-style
  primary-decomposition check): PR against this README's "The shape
  of `X`" section.
- An extension to families the current theorem doesn't cover (higher
  rank, characteristic 2, …): PR against "What's NOT ruled out".
- A reference we missed: PR against the "Literature-search task"
  section.

## Citation suggestion

If this note is useful in your work, please cite as:

```
matmul-catalog contributors. "Hopcroft–Kerr 1971: an impossibility
theorem for the (2,2,bridge-3) same-method pair." matmul-catalog
references, 2026. https://github.com/<repo>/references/hopcroftkerr1971/
```

(Replace `<repo>` once the repository is public.)
