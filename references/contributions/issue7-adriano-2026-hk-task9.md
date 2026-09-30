# Issue #7 — Marcos Adriano (2026): the six `g ≥ 6` ⟨2,p,n⟩ shapes at the exact Hopcroft–Kerr formula, over ℤ

**Source.** GitHub issue
[solven-eu/matmulcatalog#7](https://github.com/solven-eu/matmulcatalog/issues/7)
(opened 2026-09-24 by [@marcosadriano27](https://github.com/marcosadriano27),
"Task #9: the six g ≥ 6 ⟨2,p,n⟩ shapes at the exact Hopcroft-Kerr formula, with
integer coefficients"), attachment `hk-task9-schemes.zip`
(<https://github.com/user-attachments/files/32606745/hk-task9-schemes.zip>,
143 768 bytes, six JSON files `adriano-{shape}_m{rank}.json` in this repo's
sparse convention). The contributor states the work was done *with the help
of AI tools*. Imported 2026-09-30 (this note is the verification record).

**What it closes.** Task #9 of
[`research/hopcroft-kerr-2np/`](../../research/hopcroft-kerr-2np/OVERVIEW.md) —
the only shapes in the swept range `3 ≤ p ≤ 32, p ≤ n ≤ 32` where our own
emitter (`HopcroftKerr2bcAsymmetric`, `schemes/constructed/`) could not reach
`⌈(3pn + max(p,n))/2⌉`, provably so *within its framework* (arc-sum placement +
the `(3,3,bridge-1/2)` theorem over the 9-product reusable set). FMM-Lille lists
the formula rank for all six but backs none with an artifact (our cross-check
class: `index_only`, "UPSTREAM-UNVERIFIED").

| shape (ℤ) | HK formula | ours before (`constructed/`) | FMM index / artifact | **contributed** | catalog file |
|---|---:|---:|---:|---:|---|
| ⟨2,12,18⟩ | 333 | 334 (+1) | 333 / 334 | **333** | `known/section18/2x12x18-r333-adriano_2026-d046d2b.json` |
| ⟨2,14,21⟩ | 452 | 453 (+1) | 452 / 453 | **452** | `known/section21/2x14x21-r452-adriano_2026-365353f.json` |
| ⟨2,16,24⟩ | 588 | 590 (+2) | 588 / 592 | **588** | `known/section24/2x16x24-r588-adriano_2026-b4f0c4f.json` |
| ⟨2,18,27⟩ | 743 | 745 (+2) | 743 / 747 | **743** | `known/section27/2x18x27-r743-adriano_2026-47aa850.json` |
| ⟨2,20,30⟩ | 915 | 918 (+3) | 915 / 920 | **915** | `known/section30/2x20x30-r915-adriano_2026-aaf5f01.json` |
| ⟨2,24,30⟩ | 1095 | 1096 (+1) | 1095 / 1104 | **1095** | `known/section30/2x24x30-r1095-adriano_2026-67ab18d.json` |

(FMM artifact ranks from
[`references/fmm-report-16x20x28.md`](../fmm-report-16x20x28.md); `g = gcd(n, p/2)`
= 6, 7, 8, 9, 10, 6 respectively.)

## What we verified (and what we did not)

Every file went through `ImportContributedSchemes`
(`eu.solven.matmul.docs.migrate`) — dry-run first, then `--execute`:

1. **Exact symbolic proof** — `Verifier.isExactNonCubic`: every coefficient is
   scaled to a `BigInteger` by the lcm denominator (here 1) and the full
   bilinear identity `Σ_k U[a][k]·V[b][k]·W[c][k] = T[a][b][c]` is checked
   term-wise. No floating point in the decision. **6/6 exact.**
2. **Random matmul spot-check** — `Verifier.passesRandomMatmulSpotCheck` (an
   independent numeric path: run the scheme on random `A,B`). **6/6.**
3. **Fields from the coefficients, not from the contributor's `fields[]`** —
   all-integer ⇒ `Z`; exact mod-2 / mod-3 reductions
   (`Verifier.isExactNonCubicFp`) ⇒ `F2`, `F3`; `Q`, `R`, `C` by inclusion.
   **All six: `[F2, F3, Z, Q, R, C]`.** Not ternary (`zt:false`): max |c| =
   9 / 10 / 19 / 63 / 35 / 17 (the issue says "between 9 and 63" — confirmed).
4. **Canonical-disk round-trip** — the written file (`SchemeIO.write` +
   `MatrixJsonFormatter`) is re-read, must re-verify exactly and hash
   identically. **6/6.** `VerifyOneScheme` on the committed files reports
   `EXACT_SYMBOLIC = true` (34 550 terms for ⟨2,12,18⟩, 144 541 for
   ⟨2,24,30⟩; a few ms each).
5. **Gate sanity** — a copy of ⟨2,12,18⟩ with ONE coefficient sign flipped is
   rejected (`exact=false spot=false`; it still reduces mod 2, correctly).
   Regression tests: `TestImportContributedSchemes` (gate accepts / rejects /
   never trusts claimed fields; the six files are on disk at the formula) and
   `TestSweepSpotsSota.retains_issue7_hk_task9_formula_schemes` (rank ≤ formula
   via `FieldAwareLookup("Q")`).

Whole run: ~4.6 s for the six files (index build included).

**Not verified here:** the construction narrative (tori / eigenlines, the
tripartite graph, the virtual-vertex bridge, the rational pencil for
`2 ≤ p ≤ n ≤ 32` "496/496"). We hold explicit schemes for the six shapes only;
the structural claims below are reported as the contributor's, and labelled so
in the scheme JSON (`discovery_note`).

## Import command (reproducible)

```bash
# provenance (UTF-8 — Windows CLI args are Cp1252-decoded, so ⌈…⌉ / ⟨…⟩ go through --meta):
#   references/contributions/issue7-adriano-2026.meta.json
mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.ImportContributedSchemes \
    -Dexec.args="<unzipped hk-task9-schemes dir> --meta=references/contributions/issue7-adriano-2026.meta.json --field=Q --execute"
```

Metadata decisions (per CLAUDE.md "distinguish discoveries from re-discoveries"
and the field discipline):

- `source: "Adriano 2026"`, `year: 2026`, `source_url: …/issues/7`;
- `discovery: false` — the *rank value* is Hopcroft & Kerr 1971's Theorem-1
  formula, listed by FMM-Lille; `attribution_for_rank` names HK71 **and**
  records that this is the **first explicit scheme** at the shape. (Our own
  reading — [`LOGIC_AND_LIMITS.md`](../../research/hopcroft-kerr-2np/LOGIC_AND_LIMITS.md),
  Layer 3 — is that HK71's proof has a gap at exactly these shapes, so the
  explicit witness is what turns the 1971 *claim* into a certified bound
  here; the attribution string says so.)
- `commutative: false`; `verified: true`; standard Phase-2 metrics stamped on
  import (`additions`, `has_buds`, `projection_margin`, `zt`).
- The superseded `constructed/` files (334 … 1096) are **kept**: they are the
  emitter's honest output, still `+k`-labelled, and the catalog keeps every
  scheme per shape (chronology, not selectivity).

## The contributor's construction (as described in issue #7 — not re-derived)

1. *Methods as tori.* Each HK diagonal method is a pair of eigenlines (a
   torus in some basis). Lemma 2 costs 3 new products **iff the two tori share
   exactly one line**, 4 otherwise — checked exhaustively over F₂, F₃, F₅ (all
   torus pairs) and F₇, F₁₁, F₁₃ (standard torus vs all), "only if" proved
   over ℚ by Gröbner bases (also modulo every prime ≤ 10 000). HK's three
   methods are the three tori of a triangle of lines — hence equal-method
   pairs need a bridge.
2. *Five shapes without bridges.* The three HK methods form a unimodular
   trio; on a tripartite graph (one method per part, complete tripartite minus
   a matching) every edge is an integer Lemma 2. Lemma-1 windows need not be
   cyclic or contiguous — only the minors actually used must be invertible:
   `M = [I; B]`, `B ∈ {−1,0,1}`, found by local search so every used window
   has det ±1. This sidesteps our circulant arc-sum obstruction.
3. *⟨2,24,30⟩: a bridge through a virtual vertex, (3,3) included.* Degree 23
   on 30 vertices forces equal-method edges (Turán: four pairwise-compatible
   integer methods do not exist). For an equal-method pair `(i,j)` with bridge
   `b`: two products of the pair `(i,b)` already span the torus of the virtual
   row `ā_i + ā_b` (column `x_i − x_b`), whose eigenlines are those of `i` and
   `b`; that torus shares exactly one line with `j`'s, so Lemma 2 on
   `(virtual, j)` costs 3 new products and yields `y_ij + y_bj` and
   `±(y_ji − y_jb)`; the pair `(b,j)`, computed anyway, then gives `y_ij`,
   `y_ji`. Works for every equal-method pair, **(3,3) included**. Exhaustively
   over F₂, F₃, F₅: never with 6 or 9 reusable products (consistent with our
   theorem), always with **12**.
4. *Rationally, no bridges at all.* A pencil of tori through one common line
   gives arbitrarily many pairwise-compatible methods, so any `(p−1)`-regular
   graph reaches the formula — every ⟨2,p,n⟩ with `2 ≤ p ≤ n ≤ 32` exactly at
   the formula (496/496, exact checks per the contributor), but with
   denominators `det(P_i, P_j)`; hence the integer schemes use 2–3 instead.
5. Side remark: the Lemma-2 identity is not unique — over F_p the 3-product
   solutions form `(p−1)²` spaces, a two-parameter family (may help lower the
   coefficients).

**Relation to our theorems.** Nothing we proved is contradicted: the
`(3,3,bridge-1/2)` impossibility is scoped to the **9-product** reusable set
and the contributor's exhaustive F₂/F₃/F₅ search agrees (never 6 or 9). Our
first-listed resume point — *does a (3,3,bridge) completion exist over the
**12-product** set?* — is answered **yes** by construction 3 (the same
enlargement that had already flipped `(2,2,bridge-3)` from impossible to
solvable in this repo). The "Layer 3" fuzziness (the bound never proved at
these shapes) is settled *for the six shapes* by explicit witnesses; the
general statement (all `g ≥ 6`, or all `(p,n)`) rests on the contributor's
construction 4 and is not certified here.

## Follow-ups

- Ask the contributor for the generator (code or a write-up) so the
  construction can be ported into `HopcroftKerr2bcAsymmetric` and the whole
  `g ≥ 6` family beyond the six shapes emitted + certified by us; the rational
  pencil (item 4) would make the emitter formula-exact for every `(p,n)`.
- Re-run the closure sweep over shapes that use ⟨2,12,18⟩ … ⟨2,24,30⟩ as
  Kron / concat / recombination ingredients (`SchemeSweep --field=Q`): the
  derived stubs `2x16x24-r591`, `2x18x27-r746`, `2x20x30-r919`,
  `2x24x30-r1100`, and anything built on them, are now improvable.
- Regenerate `references/fmm-cross-check.md` / `docs/comparison/fmm-gap-report.md`
  (`FmmCrossCheck`, `GenerateFmmGapReport`): the six "UPSTREAM-UNVERIFIED"
  rows become ties with FMM's index.
- The paper (`paper/sections/hk71.tex`) now carries a 2026-09-30 checkpoint
  paragraph; the arXiv text should cite the issue (`adriano2026hk` in
  `paper/refs.bib`).
