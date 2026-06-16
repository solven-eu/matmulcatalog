# Published typos in matrix-multiplication literature

This file tracks **typographical and transcription errors** in published
papers and theses that affect the matmul-rank literature — both formal
errata acknowledged by authors and "silent" errors we have caught while
implementing the algorithms in this catalog.

**Why this file exists.** One of the project's load-bearing arguments is
that the small-matmul literature is *fragile*: results get lost
([Makarov 1986](case_studies/makarov_1986_recovery.md)), citations drift
([alphatensor-* re-attribution](../DISCOVERIES_PENDING_ANALYSIS.md)),
and formulas get re-typeset incorrectly between thesis → journal →
secondary citation. A single transcription error in a formula can make
the algorithm unreproducible without independent re-derivation. Keeping
a public catalog of caught typos:

- **Helps reproducers**: anyone implementing the same construction can
  consult this file before debugging their own correct code against an
  incorrect reference.
- **Demonstrates the value of the catalog**: each entry is a concrete
  case where the open + cross-referenced setup of this repo caught a
  problem that a one-pass reader of any individual source would have
  missed.
- **Distinguishes "literature said X" from "the algorithm actually
  does X"**: we can quote the corrected formula and link to the typo
  so downstream users aren't propagating a misreading.

If you implement an algorithm from a paper cited here and find a
discrepancy, add an entry — even one-line tentative notes are useful.

---

## Confirmed typos

### Drevet-Islam-Schost 2009 — Lemma 4, 7th `u`-correction

- **Source**: Drevet, Islam, Schost. *Optimization techniques for small
  matrix multiplication*. TCS 412(22):2219–2236, 2011 (preprint 2009).
  See [`drisc09`](../REFERENCES.md#10-drisc09).
- **Location**: Appendix A.1 (page 29 of the preprint), the explicit
  expression for `t(Ã, B̃, C̃)` in Lemma 4. Among the eight `u`-correction
  terms listed, the **7th** reads:

  ```
  − u(−Ã^{2,2}, B̃^{1,2}, −C̃^{1,2}, C̃^{1,2}, C̃^{2,2})
  ```

  The 3rd argument is `−C̃^{1,2}`.

- **Correct version** (Islam's 2009 MSc thesis, page 50, same Lemma):

  ```
  − u(−Ã^{2,2}, B̃^{1,2}, −C̃^{2,1}, C̃^{1,2}, C̃^{2,2})
  ```

  The 3rd argument is **`−C̃^{2,1}`**.

- **Why the thesis is right**: derive the term directly. The 7th
  `u`-correction comes from applying Lemma 1 to the off-diag-2 outer
  T-term, whose arguments are
  `(A, B, C, U, V, W, X, Y, Z) = (Ã^{1,2}, B̃^{2,2}, C̃^{1,2}, Ã^{2,1}, B̃^{1,2}, C̃^{2,2}, −Ã^{2,2}, B̃^{2,1}, −C̃^{2,1})`.
  The third correction is `U(X, V, Z, C, W)` with
  `(X, V, Z, C, W) = (−Ã^{2,2}, B̃^{1,2}, −C̃^{2,1}, C̃^{1,2}, C̃^{2,2})`.
  Slot 3 must hold `Z = −C̃^{2,1}` — the thesis version is consistent
  with Lemma 1 + the off-diag-2 argument mapping, the published version
  is not.

- **Likely cause**: copy-paste during re-typesetting for the journal
  version. The neighbouring slot 4 holds `C̃^{1,2}`; the typo replaces
  slot 3 with the same.

- **Impact**: a naive transcription of DIS09 Lemma 4 into a sympy
  implementation will leave a 256-term residual against
  `T(Ã, B̃, C̃)` for the smallest non-trivial case (m = 3). The fix
  is one character (`1,2` → `2,1`).

- **Cross-checks attempted**: Islam's Magma validation code was hosted
  at `http://www.csd.uwo.ca/~mislam63/` per DIS09's appendix preamble;
  the URL is no longer reachable as of 2026-05. **Recovered via Wayback
  Machine**:
  [`https://web.archive.org/web/20120223044300/http://www.csd.uwo.ca:80/~mislam63/TA.mgm`](https://web.archive.org/web/20120223044300/http://www.csd.uwo.ca:80/~mislam63/TA.mgm),
  archived locally at [`references/islam2009/magma/TA.mgm`](islam2009/magma/TA.mgm).
  **Line 263** of the Magma source reads `u1(-A22, B12, -C21, C12, C22)`
  — the 3rd argument is `-C21`, confirming Islam's thesis. The Magma
  is the actual algorithm Islam used to validate the DIS09 appendix
  claims (per the preamble note); it pre-dates the journal-version
  typo and is the authoritative source for the formula.

---

## Style for new entries

```
### <Author Year> — <where in the paper>

- **Source**: full citation + link to REFERENCES.md entry.
- **Location**: page, equation/lemma number, exact symbol or formula
  fragment that contains the typo.
- **As published**: quoted verbatim (use a code block).
- **Correct version**: quoted with the fix; cite the source that has
  it right (older draft, thesis, peer correction, our own derivation).
- **Why the correct version is right**: one or two sentences of math
  derivation OR a pointer to the script in this repo that verifies it
  (e.g. a sympy file in `references/<author><year>/sympy/`).
- **Likely cause**: optional — speculation about HOW the typo happened
  (re-typesetting, index-relabeling, OCR, …). Useful for spotting
  similar patterns elsewhere.
- **Impact**: what breaks if you trust the typo'd version.
- **Cross-checks attempted**: other sources you consulted to confirm
  (other editions, errata pages, reference implementations).
```

## Related drift / fragility files

- [`case_studies/makarov_1986_recovery.md`](case_studies/makarov_1986_recovery.md) — Makarov ⟨3,3,3⟩=22, recovered from Russian original after Islam thesis transcription typo + paywall + access-loss.
- [`../DISCOVERIES_PENDING_ANALYSIS.md`](../DISCOVERIES_PENDING_ANALYSIS.md) — attribution audit (which schemes are genuine "discoveries" by their bulk-import source vs re-derivations of older results).

## On using the Wayback Machine as a recovery channel

The Drevet–Islam–Schost case demonstrates a pattern worth noting: a
reference implementation that the original authors pointed to as the
canonical validator (DIS09 Appendix opens with *"an implementation in
Magma that validates the claims in this appendix is available at
http://www.csd.uwo.ca/~mislam63/"*) became unreachable some time after
publication. The university homepage hosting it was eventually
decommissioned along with the user account.

For research-grade reproducibility, **the Wayback Machine** is often
the only path back to such material. The Magma script we used to
confirm the DIS09 typo was recovered from a 2012-02-23 crawl. When
investigating a paper with broken supplementary-material URLs:

1. Construct the candidate Wayback URL: `https://web.archive.org/web/*/<original-url>`.
2. If a snapshot exists, use `https://web.archive.org/web/<timestamp>if_/<original-url>`
   (the `if_` flag returns the raw archived bytes without the Wayback
   navigation frame, suitable for `curl`).
3. Archive the recovered file locally under `references/<author><year>/`
   so subsequent verification doesn't depend on the archive's continued
   availability.

This is a generally useful method for the catalog: the matmul
literature is small and old enough that many "definitive" supplementary
artefacts (Maple files, Magma scripts, university-hosted tensors) have
moved or vanished, but the academic crawl coverage is generally good.
