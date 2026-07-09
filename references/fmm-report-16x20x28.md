# Report note — FMM-Lille ⟨16,20,28⟩:4944 (drafted 2026-07-09)

**Page:** https://fmm.univ-lille.fr/16x20x28.html
**Index rank:** 4944. **Displayed construction:**

> ⟨16×20×28:4944⟩ is serendipitous tensor product
> (⟨4×5×7:104⟩ − 17) ⊗ ⟨4×4×4:48⟩ + 8⟨4×4×8:96⟩

**Issue 1 — the construction's term arithmetic is inconsistent.**
The rank arithmetic checks out ((104−17)·48 + 8·96 = 4944), but the term
count does not: 8 fused blocks of ⟨4,4,8⟩ absorb at most 8×2 = 16 of the
base's products (each ⟨4,4,8⟩ replaces exactly 2 products of the ⟨4,4,4⟩
inner under a size-2 bud), while the formula removes **17**.

**Issue 2 — the published base cannot support the recipe.**
We analysed the published ⟨4,5,7⟩:104 tensor (4x5x7_tensor.mpl) over exact
rationals. Its factor-direction classes are:
- U-direction: 7 classes of size 2 + 1 class of size 3 (total 17 — matching
  the "−17", so the recipe surely intends fusing all U-classes);
- every class is slice-rank tight: for each class, rank(Σᵢ vᵢwᵢᵀ) equals the
  class size. In particular the size-3 class has slice rank 3, so it cannot
  be realised inside a ⟨4,4,8⟩ block (that requires slice rank ≤ 2); it
  needs ⟨4,4,12⟩ (best known 141), giving
  (104−17)·48 + 7·96 + 141 = **4989**, not 4944.
Combining all U/V/W classes of the published base never prices below 4989.

**Issue 3 — the tensor artifact is an empty placeholder.**
https://fmm.univ-lille.fr/16x20x28_tensor.mpl.bz2 decompresses to ~10 KB
containing only the symbolic A/B/C matrix declarations and a trailing
`Tensor:=Tensor:` — no tensor data. The same placeholder format affects
16x20x29 and 20x24x25.

**Questions for the maintainers:**
1. Is 4944 backed by an unpublished explicit scheme (if so, could the
   artifact link be refreshed)?
2. If the serendipitous recipe is the sole support, which ⟨4,5,7⟩:104
   representative does it use? Ours and the published one provably cannot
   realise 8×⟨4,4,8⟩ from 17 products (slice-rank argument above); a
   17-into-8 fusion would need a base whose size-3 class has slice rank 2 —
   which would imply a rank-103 ⟨4,5,7⟩ scheme by direct rewriting, so the
   displayed recipe seems unrealisable as stated.

**Related sharpened case — ⟨2,12,18⟩ (and the ⟨2,2t,3t⟩ family):** the index
claims 333 (= the Hopcroft–Kerr 1971 closed form ⌈(3mn+max)/2⌉); the
published artifact is a 334-product scheme that we verified to be locally
rigid — zero proportional triads, and every factor-direction class span-tight
(52×size-2 + 2×size-4 in both U and W roles, all full rank), so no
multiplication is mergeable by linear analysis. Our own constructive work on
HK71's even case (independent port, 465 machine-verified schemes) proves the
formula unattainable for this g ≥ 6 circulant family within the natural
local-atom framework. Question: what supports the 333 index entry? Same
pattern for the whole family: ⟨2,14,21⟩ 452/453, ⟨2,16,24⟩ 588/592,
⟨2,18,27⟩ 743/747, ⟨2,20,30⟩ 915/920, ⟨2,24,30⟩ 1095/1104 (index/artifact).

**Related case — ⟨8,27,30⟩:** the display recipe
`(⟨4,9,10⟩:250 − 12) ⊗ ⟨2,3,3⟩:15 + 6·⟨4,3,3⟩:29 = 3744` is not realizable
from the published ⟨4,9,10⟩:250 artifact: its bud census is exactly six
size-2 classes, but they share the A-factor (fusing to ⟨2,3,6⟩, a loss),
while the recipe's ⟨4,3,3⟩ blocks require six shared-B classes, of which the
artifact has none — the "−12" count matches, the axis does not.

**Context:** part of a broader index-vs-artifact audit
(`references/fmm-artifact-audit.md`): 19 of the 44 shapes where the FMM
index beats our catalog have artifacts that do NOT attain the index
(examples: 2x16x24 index 588 / artifact 592; 21x28x30 index 9473 /
artifact 9782 — the latter's index we since reproduced via its index-page
recipe, so there the artifact is merely stale).
