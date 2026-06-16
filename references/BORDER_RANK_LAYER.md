> **SUPERSEDED (2026-06-05).** This whole scope was based on a misdiagnosis.
> Perminov's draft (*Meta Flip Graph meets Serendipitous Product*, Def 2.9–2.12)
> shows the serendipitous product is **exact and combinatorial**: decompose the
> base scheme into elementary matmul tensors by its **buds** (rank-one terms
> sharing a `u`/`v`/`w` vector, up to scaling), then tensor each block at its
> best known rank. **No border rank, no ALS, no SAT.** Verified: HK `⟨2,3,4⟩=20`
> has 4 `U`-buds → `12·⟨1,1,1⟩+4·⟨1,1,2⟩` → `⊗⟨2,4,2⟩` = `12·14+4·26 = 272`.
> Border rank remains a real but SEPARATE capability (Bini/Schönhage bounds, ω,
> τ-theorem) — not needed for FMM's exact serendipitous schemes. The Milestone-0
> ALS/SAT experiments answered a non-question. See the correction banner in
> `SERENDIPITOUS_PARTIAL_PRODUCT.md`. Kept only as a record.

# Scope: the "border-rank" layer (#159 next layer)

Scoping note (2026-06-05). Requested as "scope border-rank layer" after the
finding (§4.6 of `SERENDIPITOUS_PARTIAL_PRODUCT.md`) that FMM's serendipitous
`⟨4,8,12⟩=272` is NOT a subset-deletion of any `⟨2,3,4⟩=20` scheme.

## 0. Headline: it splits into TWO distinct capabilities

Grounding the math shows the request actually contains two *different* layers,
and the one needed to rediscover FMM's exact serendipitous results is **not**
border rank:

| | **Layer A — partial matrix multiplication (EXACT rank)** | **Layer B — border rank (ε-degenerate)** |
| --- | --- | --- |
| object | exact decomposition of a *partial* tensor `D' = T_⟨n,m,p⟩ − ⊕ blocks` | ε-polynomial decomposition `Σ uℓ(ε)vℓ(ε)wℓ(ε) = ε^q·T + O(ε^{q+1})` |
| produces | **exact** catalog schemes (rational factor matrices) | a **border-rank BOUND** (limit object), not an exact scheme |
| rediscovers FMM 272? | **YES** (see §1) | **No** — a border component yields only a border bound on ⟨4,8,12⟩ |
| feeds | the scheme catalog + additions | ω-history / `cited-bounds` / τ-theorem disjoint sums |
| reuses | `Als(n,r,target)`, SAT pipeline, `SerendipitousProduct.build`, `Verifier` | needs NEW `R[ε]`-ALS + completion infra |

**Why FMM 272 is Layer A, not B.** FMM publishes 272 as an *exact* scheme
(rational factor matrices, no ε). A border (ε) component tensored with B and
corrected still carries ε; killing it needs `lim_{ε→0}`, which gives the *border
rank* of ⟨4,8,12⟩, not an exact rank-272 scheme. So for an **exact** 272 the
borrowed part `D'` must have **exact rank ≤ 12** — a Layer-A fact. Border rank is
a real and valuable capability, but it is the wrong tool for reproducing FMM's
exact serendipitous entries. (The literature interleaves them — Schönhage 1981
covers both — which is why the request reads as "border rank".)

## 1. The reframing that unlocks Layer A (corrects the earlier engine)

The earlier recognizer searched for a 0/1 **subset of an existing `⟨2,3,4⟩=20`
scheme** equal to `⊕⁴⟨1,1,2⟩`. That was the wrong object — none exists (proven
for AT-Z/F2 and FMM's own HK base). The right object is a **fresh decomposition
of the partial tensor**:

```
D' = T_⟨2,3,4⟩  −  ⊕⁴ T_⟨1,1,2⟩            (matmul minus 4 chosen output-blocks)
```
If `R(D') ≤ 12` (exact), then
```
⟨4,8,12⟩ = D'⊗⟨2,4,2⟩ (= 12·14 = 168)  +  4·⟨2,4,4⟩ (= 104)  = 272,  exact.
```
`D'` is computed FRESH (ALS/SAT), independent of any catalogued ⟨2,3,4⟩ scheme.
The borrowed "−8" of FMM's notation is the rank drop `20 → 12` achieved by *not*
computing 8 of the bilinear directions — those become the 4 `⟨2,4,4⟩` corrections.

## 2. Milestone 0 — the decisive, cheap experiment (do first)

**Question:** does `R(D' = ⟨2,3,4⟩ − ⊕⁴⟨1,1,2⟩) ≤ 12` hold exactly?

- Build the `double[][][]` target `D'` (matmul tensor minus the 4 chosen
  output-blocks; pick the output-tiling of C=2×4 into 4 `⟨1,1,2⟩`).
- Run `Als(·, r=12, D')` (general ALS, already in `search/als/Als.java`); also
  try SAT (`SatMatmulPipeline`) for an exact rational/F2 witness.
- If a rank-12 decomposition is found → **Layer B is unnecessary for 272**; the
  whole engine reduces to Layer A (ALS/SAT + the existing `SerendipitousProduct.
  build` + `Verifier`). Big de-risk.
- If exact rank-12 does NOT exist but border rank 12 does → 272 is genuinely a
  border construction and Layer B is required (less likely given FMM lists it
  exact).

Effort: ~½ day (target builder + driver). Reuses existing ALS/SAT/verify. This
single result decides whether the rest is small (A) or large (B).

## 3. Layer A build (if M0 succeeds — expected path)

1. `PartialTensor` builder: `T_⟨n,m,p⟩` minus a chosen `⊕⟨aⱼ,bⱼ,cⱼ⟩` block set
   (output-tiling enumerator over the C grid).
2. `PartialRankSearch`: drive `Als`/SAT to find a min-rank exact decomposition of
   `D'`; anytime ⇒ tier **bound** unless SAT proves minimality (then **proven**).
3. Assemble: `SerendipitousProduct.build(D'-decomp-as-base, removed=∅, B,
   corrections)` — already implemented and verified; the partial decomposition
   plays the role of `(A⊖S)`. `Verifier.isExactNonCubic` certifies.
4. Search wrapper: for target `⟨n,m,p⟩`, enumerate factor splits `⟨a,b,c⟩⊗B`,
   correction block-tilings, and `R(D')` budgets; keep verified wins below the
   naive-Kron rank. Cost guard + progress/ETA per CLAUDE.md.
5. Lineage: reuse the `BorrowAndCorrect` node (promote predictor→constructor),
   render `(⟨a,b,c⟩ partial r) ⊗ B ⊕ Σ cⱼ·Cⱼ`.

Effort: ~2–3 days. No new heavy math; the risk is ALS convergence at the target
rank (warm-start from naive, multiple restarts; SAT fallback for small cases).

## 4. Layer B build (only if M0 forces it — large)

1. `R[ε]`-arithmetic (truncated polynomials mod `ε^{q+1}`) for U/V/W coefficients.
2. `BorderAls`: ALS minimising the residual order in `ε` — find rank-`r`
   `Σ uℓ(ε)vℓ(ε)wℓ(ε) = ε^q·T + O(ε^{q+1})`. Seed from known Bini/Schönhage forms.
3. Completion / ε-cancellation: combine border component with exact corrections;
   where ε cannot be cancelled exactly, the result is a **border bound** only.
4. Honesty (CLAUDE.md optimality + field discipline): border-rank outputs are
   **bounds**, tagged `R̃` distinctly from exact `R`; they feed `cited-bounds.json`
   / ω-history, NOT the exact scheme catalog. Never conflate `R̃` with `R`.

Anchors: Bini 1979 `R̃(⟨2,2,2⟩) ≤ 5` (ref [6]); Schönhage 1981 partial/total &
asymptotic-sum inequality; BCS97 §15 (border rank), §14–16 (partial MM).

Effort: ~1–2 weeks; new infra (R[ε] ALS) is the bulk. Lower priority than A
unless we specifically want border-rank ω contributions (overlaps task #159's
sibling #159-area and the τ-theorem strand #159/#160).

## 5. Recommendation

Do **Milestone 0** first — it's cheap and tells us whether we even need Layer B.
Expected outcome: `R(D')=12` exact ⇒ build **Layer A** only, which rediscovers
FMM's exact serendipitous results (272 and family) with the engine already in
place. Treat **Layer B (border rank)** as a separate, later capability for
border *bounds* / ω, not for exact catalog schemes — and only if a concrete
target needs it.

## 6. Infra inventory

- `search/als/Als.java` — general ALS, `Als(n, r, double[][][] target)`. Check/extend
  for rectangular (non-square) targets; `Z3Als` is the symmetric/cubic variant.
- `f2/sat/SatMatmulPipeline` + Kissat/CryptoMiniSat/Sat4j — exact F2 witnesses.
- `SerendipitousProduct.build` + `Verifier.isExactNonCubic` — combine + certify (done).
- `KnownTauIdentities.BorrowAndCorrect` — rank predictor to promote to constructor.
- No border-rank / `R[ε]` infra exists yet (Layer B is greenfield).
