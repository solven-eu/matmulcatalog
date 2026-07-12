# Draft report to FMM-Lille maintainers (prepared 2026-07-09)

*Consolidates: fmm-artifact-audit.{md,json}, fmm-report-16x20x28.md, and the
2026-07-08/09 research-log entries. Numbers verified against the digest of
2026-07-08 and per-shape artifact downloads (provenance-checked: every count
below is from fmm.univ-lille.fr's own files, not the Perminov mirror).*

---

Subject: fmm.univ-lille.fr — index/artifact discrepancies and two
inconsistent construction descriptions (44-shape audit)

Dear Alexandre,

while systematically cross-checking our catalog against fmm.univ-lille.fr we
audited every shape where your index rank beats our catalog (44 shapes at
audit time), downloading each `_tensor.mpl.bz2` artifact and counting its
actual rank. Three groups of findings you may want to look at:

## 1. Index rank not attained by the published artifact (19 shapes)

| shape | index | artifact | | shape | index | artifact |
|---|---:|---:|---|---|---:|---:|
| 2x12x18 | 333 | 334 | | 22x28x28 | 9383 | 9552 |
| 2x14x21 | 452 | 453 | | 22x30x30 | 10534 | 10705 |
| 2x16x24 | 588 | 592 | | 27x28x28 | 10413 | 10442 |
| 2x18x27 | 743 | 747 | | 27x28x29 | 11169 | 11198 |
| 2x20x30 | 915 | 920 | | 27x28x30 | 11561 | 11590 |
| 2x24x30 | 1095 | 1104 | | 27x28x31 | 12072 | 12117 |
| 21x28x30 | 9473 | 9782 | | 27x28x32 | 12339 | 12394 |
| 21x28x32 | 10134 | 10477 | | 27x29x31 | 12890 | 12942 |
| 22x23x23 | 6435 | 6501 | | 27x29x32 | 13187 | 13234 |
| 4x18x30 | 1394 | 1397 | | | | |

(The 27x·x· and 4x18x30 artifacts are parametric — symbolic coefficients —
counted structurally, calibrated on 2x12x18.)

Some are surely just stale artifact links: we reproduced 21x28x30 = 9473 and
4x18x30 = 1394 ourselves from your index-page construction descriptions, so
those index ranks are correct. For the ⟨2,m,n⟩ family the situation seems
different: the index equals the Hopcroft–Kerr closed form ⌈(3mn+max)/2⌉, we
verified your 2x12x18 artifact (334) is locally rigid (no proportional
triads; every factor-direction class span-tight), and our own constructive
work on HK71's even case — 465 machine-verified schemes, exact at the
formula everywhere except the g ≥ 6 circulant family — includes an
impossibility theorem suggesting the formula is not attainable there with
local atoms. Question: is 333 backed by an unpublished scheme, or is it the
paper bound?

## 2. Placeholder artifacts (3 shapes)

16x20x28, 16x20x29, 20x24x25: the `.mpl` files contain only the symbolic
A/B/C declarations ending `Tensor:=Tensor:` — no tensor data.

## 3. Construction descriptions that don't check out (2 shapes)

- **16x20x28**: "(⟨4,5,7⟩:104 − 17) ⊗ ⟨4,4,4⟩:48 + 8⟨4,4,8⟩:96". The rank
  arithmetic works, but 8 blocks of ⟨4,4,8⟩ absorb at most 16 of the base's
  products, not 17 — and your published ⟨4,5,7⟩:104 cannot support it: its
  U-classes are 7×size-2 + 1×size-3 (total 17), every class slice-rank-tight,
  so the size-3 class needs ⟨4,4,12⟩ and the best realisable total is 4989.
  (A base whose size-3 class had slice rank 2 would imply ⟨4,5,7⟩ = 103 by
  direct rewriting.)
- **8x27x30**: "(⟨4,9,10⟩:250 − 12) ⊗ ⟨2,3,3⟩:15 + 6⟨4,3,3⟩:29". Your
  published ⟨4,9,10⟩:250's bud census is exactly six size-2 classes — but
  they share the A-factor (fusing to ⟨2,3,6⟩, a loss), while ⟨4,3,3⟩ blocks
  require six shared-B classes, of which the artifact has none. The "−12"
  count matches, the axis does not.

Also worth noting: several description cells omit the fusion that makes the
stated rank (e.g. 17x17x19 lists seven leaves summing 3267 for a 3266 claim;
its artifact realises the −1 as 83 cross-leaf products replacing 84).

We'd be happy to share our verification scripts. And thanks for the catalog —
the per-shape construction descriptions in particular have been invaluable:
several of them led us to constructions our own search had missed.

Best regards,
Benoit

---
*Prep notes (not part of the email): our catalog now stands WORSE=6 /
BETTER=1476 vs the audited digest; the six remaining are the four −1
diffuse-sharing shapes + 7x11x30 (−2) + 8x27x30 (−6, item 3 above).*
