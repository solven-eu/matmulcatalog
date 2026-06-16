# How a 1986 matrix-multiplication algorithm got effectively lost — and why a public catalog matters

A small case study, written during the recovery of Makarov 1986's
explicit `⟨3,3,3⟩ = 22` non-bilinear commutative scheme. The goal of
this note is to document concretely how a *cited, useful, never
superseded for >30 years* algorithm can be hard for a modern
researcher to actually use, and to argue from that observation that
a publicly maintained matmul-scheme catalog (this repo's mission)
has a real role to play.

## The result we wanted to add

Makarov 1986 — "An algorithm for multiplication of 3×3 matrices",
*USSR Comput. Math. Math. Phys.* 26(2):293–294 + 320 — gives an
explicit 22-multiplication algorithm for the `⟨3,3,3⟩` matmul
product. The 22-mult result is:

- **Commutative**: each rank-1 product mixes entries of `A` and `B`
  in a single factor, so the algorithm requires the scalar entries
  to commute. It does not lift to recursive matmul over
  non-commutative rings.
- **Better than Laderman 1976's `⟨3,3,3⟩ = 23`** (which is
  non-commutative and DOES lift). The −1 is real, just bought at
  the cost of commutativity.
- **The state of the art for commutative `⟨3,3,3⟩` for 33 years** —
  improved only in 2019/2020 by Rosowski's `⟨3,3,3⟩ = 21`.

A working 3-matrix product algorithm with 22 multiplications,
published in a serious journal, recited in every textbook
discussion of small-format matmul, has been *the* commutative
`⟨3,3,3⟩` baseline for three decades.

The factor matrices are 9 × 22 = ~200 small integer coefficients.
Trivial to encode once you have them.

## What it took to actually obtain those 200 numbers in 2026

1. **Original Russian text** (Makarov 1986, *Zh. Vychisl. Mat. i
   Mat. Fiz.* 26(2)). Cyrillic. Requires either reading Russian or
   trusting a translation. PDF available at
   [mathnet.ru](https://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=zvmmf&paperid=4056&option_lang=eng).
2. **English translation** (Sci-Direct,
   [DOI 10.1016/0041-5553(86)90203-X](https://doi.org/10.1016/0041-5553(86)90203-X)).
   **Paywalled.** ~$24 per article, $36 with rental, full
   institutional access for those affiliated.
3. **Modern references** mostly cite **Islam 2009** ("Algorithms
   for Strassen-Like Matrix Multiplication", MSc thesis, University
   of Western Ontario) as the convenient secondary source. Rosowski
   2019/2020's own paper (the one that improves on Makarov to 21)
   cites Islam, not Makarov directly.
4. **Islam 2009 is itself not casually available online.** It was
   distributed as a master's thesis, not formally archived in arXiv
   / DBLP / Google Scholar with a stable open PDF. Tracking it down
   is non-trivial; the thesis is referenced in many papers and yet
   does not pop up in a regular Google search of the form
   "Makarov 22 product algorithm". It needs Western Ontario library
   access or a personal contact.
5. **Islam's transcription of Makarov has a single-index typo.**
   Out of the 22 products, exactly one — `γ_18` (= Makarov's
   `M_18`) — is transcribed with the second factor `b_{3,2}`
   instead of the original `b_{2,3}`. The typo is a single subscript
   swap in one product of one paper, but it propagates: 4 of 9
   output combinations become incorrect, the verifier residual is
   non-zero, and a naive consumer who blindly types the formula in
   gets a broken algorithm.

So the situation: **the textbook fact "Makarov 1986 multiplies 3×3
matrices in 22 commutative steps" is widely repeated, but actually
producing a working set of factor matrices in 2026 requires (a)
reading Russian, OR (b) paying a journal paywall, OR (c) finding a
private copy of an unpublished master's thesis, AND (d) catching a
transcription typo in that thesis.**

A researcher without the Russian or the journal subscription or the
thesis is, in practice, stuck.

## The Rosowski 2019/2020 paradox

Rosowski's `⟨3,3,3⟩ = 21` paper is the modern improvement over
Makarov 22. Rosowski cites Makarov (via Islam) as the prior state of
the art. Rosowski's paper, like Makarov's, gives explicit factor
matrices — *for the new 21-product algorithm*. Rosowski does NOT
include Makarov 22 in his paper. Why would he? Makarov is old, and
Rosowski's own result is strictly better. So a reader who finds
Rosowski's improvement gets the 21-product algorithm but not the
22-product one — and Makarov 22 still goes uncollected.

This is the natural state of the literature: **every paper that
*could* publish Makarov 22's matrices has either better results to
present, or doesn't think of itself as a catalog**. The original
publication itself is in a paywalled journal in Russian; everything
downstream skips over re-publishing the explicit construction.

## The typo, demonstrated

Concrete demonstration that the typo matters. Islam's transcription
gives, for the 18th product:

```
γ_18 = (a_{3,2} − a_{2,3} − a_{3,3}) · b_{3,2}
```

The correct formula, recovered from the Russian original where
`M_18 = (b_3 − c_2 − c_3) · k_6` with `k_6 = b_{2,3}` in our
indexing, is:

```
γ_18 = (a_{3,2} − a_{2,3} − a_{3,3}) · b_{2,3}
```

A single subscript swap on the second factor (`b_{3,2}` →
`b_{2,3}`, equivalent to `k_8` → `k_6` in Makarov's own naming).

The downstream impact:

- With Islam's transcription, the Verifier residual is
  **≈ 4.9** (sum-of-squares = 24, so on the order of 24 unit-sized
  coefficient errors).
- Specifically: outputs `c_{1,2}`, `c_{1,3}`, `c_{3,2}`, `c_{3,3}`
  are all wrong; the other 5 are correct.
- The 24 wrong coefficients are exactly the 6 mixed products
  `{a_{2,3}, a_{3,2}, a_{3,3}} × {b_{2,3}, b_{3,2}}` in 4 outputs,
  with `c_{1,2}` and `c_{1,3}` sharing one error pattern, and
  `c_{3,2}`, `c_{3,3}` sharing the opposite pattern.
- After the typo fix (`b_{3,2}` → `b_{2,3}` in `γ_18`), the
  residual drops to **exactly 0**. All 9 outputs verify.

The transcription script + the single-flip / single-removal /
single-sign-change search harness used during the recovery are at
[`src/test/java/eu/solven/matmul/research/MakarovSearch.java`](../../src/test/java/eu/solven/matmul/research/MakarovSearch.java)
and
[`MakarovPerOutput.java`](../../src/test/java/eu/solven/matmul/research/MakarovPerOutput.java).
Crucially, **no single sign-flip or single removal/addition would
have caught the typo** — the search would have come up empty. The
typo is an *index* change, not a sign change. Locating it required
cross-referencing the original Russian.

## What this argues

For a result like Makarov 22:

- **The bibliographic ledger is intact** (everyone cites it).
- **The publication record is intact** (the journal exists, the DOI
  resolves).
- **The mathematical knowledge is intact** (the algorithm has been
  rediscovered, improved on, and incorporated into the field's
  oral / textbook tradition).
- **The actual concrete numerical artifact — the 200 integers that
  let you run the algorithm — is essentially inaccessible** to a
  practitioner who lacks (a) Russian + (b) journal access + (c) the
  Islam thesis + (d) the patience to verify a thesis transcription
  against a Russian original.

The practitioner who *just wants the algorithm* — to benchmark, to
recurse on, to compose into a larger scheme, to teach with, to
double-check a result against — is locked out of a 40-year-old
public result by an entirely accidental combination of language,
paywall, archival, and transcription gaps.

A small, public, maintained, **executable** catalog of explicit
factor matrices — with each scheme verified against a reproducible
test — directly closes this gap. Reading the Russian, locating the
Islam thesis, fixing the typo: these are one-time costs. After
they're paid, the resulting factor matrices live forever in an open
file, citable, runnable, verifiable, queryable across all the
researchers and tools that come after.

This is, fundamentally, why we maintain this catalog. The Makarov
recovery makes the case concretely: even single, isolated,
already-published results can effectively leave circulation when
the access ladder gets long enough. A public catalog patches that.

## Generalisation

Makarov is not unique. We have already encountered:

- **Heun 1994** (German bilinear-complexity textbook) — referenced
  by HK 1971 follow-up work, no online English version found.
- **Probert-Fischer 1976** — referenced as the canonical
  decomposition framework for small-format upper bounds, original
  paper hard to retrieve.
- **Smith 2002 tables** — referenced as the next iteration of
  Probert-Fischer's enumeration, original publication channel
  unclear in 2026.
- **Pan 1984 monograph** "How to Multiply Matrices Faster" (LNCS
  179) — Springer-paywalled, no open PDF found.

Each of these has the same shape as Makarov 1986: technically
public, practically scarce. Each is a candidate for the same
"recover the explicit construction, verify, archive in the
catalog, cite the source faithfully" treatment.

### A second pattern: thesis → published paper attribution drift

Makarov-like access barriers are one failure mode. A related but
distinct one is the **thesis → shared-paper attribution chain**.
Concrete example:

- **Islam 2009** (Md. Nazrul Islam MSc thesis, Schost-supervised,
  University of Western Ontario) introduces several results,
  including the closed-form Pan-style TA bound `(n³+15n²+14n−6)/3`
  for *odd* `n` (Proposition 1 — Islam writes verbatim *"we do not
  know of any previous mention of the case of n odd"*).
- **DIS09** (Drevet, Islam, Schost 2009) republishes Islam's thesis
  Chapters 4–5 as the formal shared paper. Islam is second author of
  his own underlying work, alongside his supervisor and a third
  co-author.
- **Modern citations** — including ours, until corrected — credit
  *"DIS09 §3"* for the formula. The thesis (single-author, less
  indexed) effectively falls out of the conversation.

This is normal academic flow: the published paper is the citable
artefact, the thesis is the underlying work. The cost is that the
proper first-author attribution of a *2009 result* (Islam's
introduction of the odd case) gets diluted into "DIS09" — exactly the
same dilution that obscures Makarov 1986's algorithm behind Islam's
re-transcription of it. In both cases, the **practical citation
attaches to the easier-to-find document**, not the original
contribution.

**And nobody is going to fix this from the inside.** Islam's MSc
thesis is not online. He does not appear to have pursued an
academic career past the master's. There is no faculty page, no
preprint server presence, no later papers that re-iterate "by
the way, the odd-case formula in DIS09 §3 is mine, see my 2009
thesis Prop 1". The original author is not around to push back on
citation drift. The drift becomes permanent by default — DIS09
becomes "the source" not because Islam-2009-the-thesis was
superseded, but because Islam-as-an-author dropped out of the
conversation.

This is exactly the case where a maintained external catalog has
the most leverage. We are not competing with active corrections from
the thesis's author — we are filling in for a quiet historical
record that nobody else has incentive to maintain. The marginal cost
of recording "this formula is Islam 2009 Prop 1, published as
DIS09 §3" in a one-line code comment is essentially zero; the
marginal benefit (future readers, students, citation accuracy)
accrues forever. A catalog that names sources precisely (Islam
thesis Prop 1 for the odd-case formula; Makarov 1986 for the
22-product algorithm; Drevet-Islam-Schost 2009 for the formal
Pan-TA paper) is one of the few places where this kind of
quiet-attribution drift can be reversed at all.

We have updated the attribution in
`src/main/java/eu/solven/matmul/papers/pan1978/PanTrilinearAggregation.java`
to credit Islam 2009 + DIS09 as the source for the formulas instead
of "DIS09 §3" alone. Other classes citing "DIS09" should be
re-audited along the same lines.

## Action items raised by this case study

1. **Errata reporting**: Islam 2009 §3.3.1's typo in γ₁₈ has been
   isolated. A short erratum could usefully be sent to the author
   or appended as a note in any future redistribution of the
   thesis. Tracked in `FUTURE_WORK.md`.
2. **Catalog tagging**: add a "rescued" or "recovered" tag (or a
   `provenance_path` field in the JSON) to schemes whose
   constructive form was effectively unavailable until our import.
   Makarov 22 is the first such entry; the field signals to a
   reader "this took some work to obtain, here are the citations".
3. **Russian / German / French language pipeline**: a non-trivial
   share of the historical matmul literature is non-English.
   Investing in OCR + machine translation for the next round of
   pre-2000 papers (Russian / Polish / German) likely yields more
   recoveries on the Makarov pattern.

## Citation

If this case study is useful to a reader, the catalog entry for
Makarov 1986 lives at
[`src/main/resources/schemes/section3/makarov-1986_3x3x3_r22_a105_commutative.json`](../../src/main/resources/schemes/section3/makarov-1986_3x3x3_r22_a105_commutative.json).
The Russian primary source is archived at
[`references/papers/makarov_1986_RU_zvmmf4056.pdf`](../papers/makarov_1986_RU_zvmmf4056.pdf).
The bibliographic ledger entry is
[`REFERENCES.md#29-makarov86`](../../REFERENCES.md#29-makarov86).
