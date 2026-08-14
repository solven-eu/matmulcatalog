# Porting spec — Schwartz–Zwecher 2025 TA-New25 cubic constructor

*(2026-08-06. Target chosen by user: the paper-faithful cubic TA-New25, to
regenerate our dense SZ imports ⟨28³⟩/⟨30³⟩/⟨32³⟩ as stubs and to build the
un-held even n₀ = 34…50. Source: arXiv:2508.01748 §3 + Appendix B; PDF held at
`references/papers/schwartz_zwecher_2025_feasible_matmul_arxiv2508.01748.pdf`.)*

## What it is (one paragraph)

TA-New25 is Pan's ⟨n₀,n₀,n₀; tPan⟩ trilinear-aggregation cubic (even n₀≠16, in
the Hadas–Schwartz [26] presentation), improved by **kin-row unification**: for
the `n₀/2+1` *diagonal* traces (i=j) of Pan's cancellation step, a **specific**
⟨2,2,2;7⟩ algorithm is chosen (Cor 3.2 / Claim 3.3) so its first product is
*kin* (identical U- and V-rows) to a product from the aggregation step; the kin
pair is then **united** (Lemma 2.20). Net saving `n₀/2+1`:
`tNew = tPan − (n₀/2 + 1) = n₀³/3 + 15/4·n₀² + 61/6·n₀ + 8`.

## The complete explicit construction (Appendix B)

Index convention: `ī = i + n₀/2+1 (mod n₀+2)`, `d = n₀/2+1`, `γ = 1 − 9/d`.

**1. Transform** (φ; padding to zero-sum, so blocks are `d×d`). With `u` the
all-ones vector of length `d`:
```
L = [ I_{n₀/2} ; −uᵀ ]              (d × n₀/2)
R = [ I_{n₀/2} − (1/d)uuᵀ ; −(1/d)u ] (n₀/2 × d)   [transpose per the paper]
φ(X) = (I₂ ⊗ L) · X · (I₂ ⊗ R)
A* = φ(A), B* = φ(B), C* = φ(C)
```
(This is the SAME zero-sum padding our `PanTrilinearAggregationBuilder.buildL/buildR`
already implement — reuse/verify against it.)

**2. Trace(ABC) = sum of four families** (each `(·)(·)(·)` triple is one product):
```
(a) aggregation-symmetric,  Ŝ = {(i,j,k) : i≤j<k or k<j≤i}:
    Σ (A*_{ij}+A*_{jk}+A*_{ki})(B*_{jk}+B*_{ki}+B*_{ij})(C*_{ki}+C*_{ij}+C*_{jk})

(b) aggregation-barred,  Ṡ\Ṡ₁ (Ṡ pairs (i,j,k),(ī,j̄,k̄); Ṡ₁ = i=j=k):
    Σ (−A*_{ij}+A*_{j̄k}+A*_{k,ī})(B*_{j,k̄}+B*_{ki}+B*_{ī,j})(−C*_{k̄,i}+C*_{i,j̄}+C*_{jk})

(c) off-diagonal cancellation,  S̃\S̃₁ = {(i,j): i≠j}, ANY ⟨2,2,2;7⟩:
    − d · Σ Trace( [[A*_{ij},A*_{īj}],[A*_{i,j̄},A*_{ī,j̄}]]
                    [[B*_{ij},B*_{īj}],[B*_{i,j̄},B*_{ī,j̄}]]
                    [[C*_{ij},−C*_{īj}],[−C*_{i,j̄},C*_{ī,j̄}]] )

(d) united diagonal,  R(i) for i∈[d]  — the 7-term family that REPLACES the
    i=j diagonal traces, already kin-united.  ⚠ THIS IS THE HAZARD: p18 gives
    R(i) as 7 explicit products with coefficients in {1, γ, 1/γ, 1/γ², d,
    d/γ, ±(γ±1)/γ, …}. It MUST be transcribed from clean source, not PDF text.
```

**3. Map back** the products (which compute `Trace(A*B*C*)`) to factor rows in
the ORIGINAL A/B/C spaces by composing with φ (U ← (I₂⊗L)ᵀ side, V ← (I₂⊗R) side;
C* → C via the padded→output remap our builder's `remapPaddedToOutput` does).

## Implementation plan

1. **New class** `papers/khoruzhii2026`-sibling: `papers/schwartzzwecher2025/TaNew25Construction.build(int n0)`
   returning `NonCubicBilinearAlgorithm` (rational `BigFraction` coeffs; γ, 1/γ
   are rational since d is integer). Follow `LitaTaConstruction` (BigFraction) as
   the port template.
2. Reuse `PanTrilinearAggregationBuilder.buildL/buildR/remapPaddedToOutput` for φ
   and the C-remap (verify they match L/R above first).
3. Emit families (a),(b),(c),(d) as `mul(u,v,w)` products; (c) with a fixed
   ⟨2,2,2;7⟩ (Winograd), (d) with the R(i) formula.
4. Wire `TrilinearAggregations.SZ.canBuild(n)=even,≠16,≥? ; .build(n)`.
5. **Verify**: `Verifier.passesRandomMatmulSpotCheck` (exact symbolic is
   O(n⁶r) — intractable ≥ n₀=20). Cross-check rank == `tNew`. Cross-check the
   built ⟨28³⟩/⟨30³⟩/⟨32³⟩ against the held dense imports (same rank; spot-check
   both agree on random A,B) — this is the strongest available correctness gate.
6. On green: `MaterialiseTaNew25Cubes` (like `MaterialiseLitaCubes`) → stubs;
   retire the dense SZ imports (they become regenerable); add even n₀=34…50.

## R(i) — the 7 united diagonal products (clean, from ar5iv; cross-checked vs PDF p18)

For each `i∈[d]`, with `ī=i+d`, emit these 7 products `(A-form)(B-form)(C-form)`:
```
1. (A*īi + A*iī − A*ii)(B*īi + B*iī + B*ii)
     (C*īī·d(1−γ)/γ − C*īi(γ−d)/γ − C*iī(−γ+d)/γ + C*ii(1−d))
2. (A*iī)(B*īī(−γ−1)/γ − B*īi/γ + B*iī(1−1/γ²) + B*ii(γ−1)/γ)
     (C*īī·d + C*īi·d + C*iī·d/γ + C*ii·d)
3. (A*iī + A*ii·γ)(B*īī(γ+1)/γ + B*īi(γ+1)/γ + B*iī/γ² + B*ii/γ)
     (C*iī·d/γ + C*ii·d)
4. (A*īi + A*ii(−γ−1))(B*īī + B*īi + B*iī/γ + B*ii)
     (C*iī·d/γ² + C*ii(d + d/γ))
5. (A*īī + A*īi − A*iī/γ − A*ii)(B*īī(−γ−1) − B*iī/γ)
     (C*īī·d(γ−1)/γ − C*īi·d/γ)
6. (A*īi − A*ii)(B*īī(−γ−1) − B*īi + B*iī(−γ−1)/γ − B*ii)
     (C*īī·d(1−γ)/γ + C*īi·d/γ − C*iī·d(γ−1)/γ² + C*ii·d/γ)
7. (A*īī + A*iī(−γ−1)/γ)(−B*īī + B*iī(γ−1)/γ)
     (C*īī·d/γ − C*īi(−d − d/γ))
```
γ = 1−9/d is rational (d integer) so all coeffs are BigFraction. UNBLOCKED.

## Verification safety net

n₀=4 (d=3, γ=1−9/3=−2) builds a valid ⟨4,4,4⟩ scheme (rank tNew(4)=130, not
competitive but CORRECT) — small enough for EXACT `SymbolicVerifier`
(O(4⁶·130)≈5·10⁵ ops). Debug the port against n₀=4 exact BEFORE scaling to
n₀=20/28 (spot-check + cross-check vs held dense imports). Note γ can be
negative/zero-ish for small d — n₀=4 gives γ=−2 (fine); n₀=16 gives d=9, γ=0
(excluded, matches paper's n₀≠16).

## RESOLVED (2026-08-06) — constructor verifies exactly, wired

`TaNew25Construction.build(n0)` reconstructs the matmul tensor **exactly**: n0=4/6/8
pass `Verifier.isExactNonCubic`; n0=10/20/28/30/32 pass the random spot-check; rank
== `tNew` at every tested n0 (28/30/32 → 10550/12688/15096, matching the held dense
imports). Wired into `TrilinearAggregations.SZ.build`/`canBuild` and replayable via
the `TA_sz(n=N)` lineage atom (`LineageReplayer`). Guards: `TestTaNew25Construction`
(exact-verify n0=4/6/8) + SZ rows in `TestTrilinearAggregations`.

**The two bugs (both found by diffing the arXiv LaTeX e-print, `arxiv.org/e-print/2508.01748`,
NOT the garbled PDF):**
1. **Transform must be uniform φ for all three matrices** — `A*=B*=C*=(I2⊗L)X(I2⊗R)`.
   The earlier "asymmetric" `B←L·B·Lᵀ` was wrong; it only looked better under the
   (also-wrong) recovery. `Tr(φA·φB·φC)=Tr(ABC)` holds because `R·L=I`.
2. **Family (b) C-form transpose** — the third term is `C*_{q,s}` (`C*_{j,k}`), not
   `C*_{s,q}`. This hit all 51 family-(b) products.

Verification of the OTHER families against the source: family (a), family (c)
Strassen-trace (all 7 Mk W-forms re-derived from `Tr((MN)·Cblk)`), and all 7 R(i)
products — checked term-by-term, all correct. Recovery is the φ-pullback shared by
U/V/W: coeff of C[x][y] in `Σγ_pq C*_pq` is `Σ_pq γ_pq·PL[p][x]·PR[y][q]`, onto
output `(AB)_{yx}` (trace dual).

**Follow-up (optional, provenance-only):** retire the 3 dense `known/` SZ imports
(28/30/32) → `TA_sz(n=N)` derived stubs (prefer-derived-over-imported). Deferred:
they are load-bearing lineage parents (pinned-to by 28x28x29 / 30x30x31 / 30x32x32),
so it needs the `RetireDerivableImports` re-pin dance, and LITA already dominates
every SZ shape on rank, so there is no SOTA gain — pure cleanup.

## Implementation status (2026-08-06, superseded — kept for the debug chronology)

`papers/schwartzzwecher2025/TaNew25Construction.java` — WIP, compiles, NOT wired.
- ✓ **Rank correct**: build(4) = 130 = tNew(4). Family structure/counts verified
  (Ŝ 16 + Ṡ\Ṡ1 51 + off-diag 42 + R(i) 21 = 130). Ŝ/Ṡ each emit BOTH the triple
  and its all-barred partner; Ṡ1 excludes only the 3 unbarred diagonals.
- ✓ **Transform pinned**: A←L·A·R, B←L·B·Lᵀ (asymmetric; uniform/swap strictly
  worse per residual, matches the verified `PanTrilinearAggregationBuilder`).
- ✓ **Recovery now CONSISTENT**: W uses the same φ-pullback as U/V's `star` —
  coeff of C[x][y] in Σγ_pq C*_pq is Σ_pq γ_pq·PL[p][x]·PR[y][q], onto output
  (AB)_{yx} (trace dual). The earlier integer `(I2⊗Lᵀ)C*(I2⊗L)` variant was
  INCONSISTENT with U/V (only looked better because Lᵀ/L are ±1). Fractional
  W-entries are legal; "fractions" was never the bug.
- ✗ **Tensor still wrong under the correct recovery**: residual 360/4096 (256
  fractional). Progress over the debug run: 848 → 416 → 376 → 332 (integer,
  inconsistent recovery) → 360 (fractional, CONSISTENT recovery). The
  consistent-recovery number is the trustworthy one.
- **Defect localised to the FAMILY FORMS, not the recovery.** Per-output-cell
  residual (`ProbeTaNew25`) is BLOCK-patterned: intra-half cells ~36–38 errors,
  cross-half ~8, near-uniform across all four families. Signature = a systematic
  transcription error in the core φ ⊗ block algebra (the `1/d` from PR not
  cancelling across a block), NOT one mistranscribed product. Residual-poking is
  exhausted.

### NEXT STEPS — residual-poking exhausted; two clean routes
1. **Coefficient-level diff** of families (a)/(b)/(c)/(d) — especially the R(i)
   7-product γ-powers and the `R = [I−(1/d)J ; −(1/d)u]` "[transpose per the
   paper]" orientation — against the clean arXiv LaTeX **e-print**
   (`arxiv.org/e-print/2508.01748`), which the PDF-text extraction garbles. This
   is the spec's day-one flagged hazard; it is the blocker.
2. **Pan-builder-delta route (recommended)**: instead of re-deriving Appendix B,
   start from the VERIFIED `PanTrilinearAggregationBuilder` (which already gets φ
   + the padded-C remap right) and add ONLY the SZ kin-row unification of the
   n0/2+1 diagonal ⟨2,2,2;7⟩ traces (Cor 3.2 / Lemma 2.20). Caveat: the repo Pan
   builder is a different variant (⟨4,4,4⟩=100 ≠ tPan), so the delta needs the
   Hadas–Schwartz presentation specifically — not free either.

Both routes are gated on the clean SZ/HS source. Value is **provenance-only**
(LITA strictly beats SZ on rank: ⟨28³⟩ 10535<10550, ⟨30³⟩ 12672<12688, ⟨32³⟩
15079<15096), so this is coverage/chronology, not a SOTA gain. Once n0=4 residual
= 0: spot-check build(20)/build(28), cross-check vs held dense SZ imports, then
wire `SZ.build`/`canBuild` + add `MaterialiseTaNew25Cubes`.

## The R(i) blocker before coding — RESOLVED (see R(i) above)

The `R(i)` family (2d) and the barred-index arithmetic must be transcribed from
**clean LaTeX**, not the garbled PDF text (fractions/subscripts mangle on
extract, and a single wrong γ-power fails silently under spot-check). Get the
arXiv source (`arxiv.org/e-print/2508.01748`) or the ar5iv HTML render, OR have
the user drop the source. Everything else in this spec is unambiguous and ready.
