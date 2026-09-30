# Hopcroft–Kerr ⟨2,p,n⟩ — constructive closure (SOLVED 2026-06-11/12)

**One-paragraph status.** The HK 1971 bound `R(⟨2,p,n⟩) ≤ ⌈(3pn+max(p,n))/2⌉`
is now attained **constructively, in integer (ℤ) schemes, machine-verified**,
for every shape `3 ≤ p ≤ 32`, `p ≤ n ≤ 32`: **465 schemes in
`schemes/constructed/`, 459 at the exact formula**, the only six exceptions
(+1..+3) being the `g ≥ 6` circulant family, which is *provably*
formula-impossible within the framework (theorem, see below). 197 + 22 of
these strictly improve on every published catalog (FMM-Lille, Perminov, our
own prior derivations). Full recipe, proofs and verification policy:
[CONSTRUCTIVE_METHOD.md](CONSTRUCTIVE_METHOD.md). Synthetic explainer — why
the construction works and why the `g ≥ 6` shapes' reachability was *unclear*
rather than merely unsolved: [LOGIC_AND_LIMITS.md](LOGIC_AND_LIMITS.md).

**Update 2026-09-30 — the six `g ≥ 6` shapes are closed, externally.**
[Issue #7](https://github.com/solven-eu/matmulcatalog/issues/7) (Marcos
Adriano) contributed explicit **integer** schemes at the exact formula for all
six (333 / 452 / 588 / 743 / 915 / 1095); we verified them exactly and imported
them as `schemes/known/*-adriano_2026-*` (record:
[`references/contributions/issue7-adriano-2026-hk-task9.md`](../../references/contributions/issue7-adriano-2026-hk-task9.md)).
His route is our first-listed resume point — a `(3,3,bridge)` completion over
the **12-product** reusable set (via a virtual row `ā_i+ā_b`), plus
non-cyclic unimodular Lemma-1 windows — so the framework theorem below stands
(9-product scope) and its boundary was, once again, not the truth boundary.
The catalog now holds the HK formula at **every** shape in the swept range;
the emitter itself still lands +1..+3 there (unchanged).

**What was built** (tasks #7, #8, #10, #11 — all complete):

- **Odd p (Case 1)** — cyclic band + exact back-substitution; at formula
  everywhere. Key unlock: arc-interior *bridge selection* (the
  published-impossible `(2,2,bridge-3)` case never arises by scheduling).
- **Even p (Case 2)** — band + circulant matching + repaired Step 3
  (operational Z-trick with free leftover-vertex placement); at formula
  except `g ≥ 6`. Key unlock: the `(2,2,b3)`/`(1,1,b3)` identities ARE
  derivable over the TRUE 12-product reusable set — the published
  impossibility holds only for its narrow 6-product S.
- **n > 2p−1 (chaining)** — DP-optimal concatenation of band segments over
  *achieved* segment ranks; at formula at every swept shape, all parities.
- **Integer schemes** — all-unimodular Lemma-1 matrices by a *Euclidean
  recursion* (comb body + transposed recursive tail, gcd descent); every
  window inverse is integer, so all 465 schemes are over ℤ
  (fields `[F2,F3,Z,Q,R,C]`).

**Theorems / conclusions worth remembering:**

1. `(3,3,bridge-1/2)` admits **no 3-product completion by ANY local rank-1
   atoms** over the emitter's 9-product reusable set (block-reduction proof,
   `references/hopcroftkerr1971/sympy/derive_33bridge_general.py`). This is
   why HK's own Lemma-3 sequence (which forces (3,3) pairs) cannot work as
   written — their Theorem-1 proof has a real gap; the bound survives by
   *avoiding* (3,3) pairs, which is possible iff `g ≤ 5` (arc-sum argument).
2. The six `g ≥ 6` shapes (`⟨2,12,18⟩+1, ⟨2,14,21⟩+1, ⟨2,16,24⟩+2,
   ⟨2,18,27⟩+2, ⟨2,20,30⟩+3, ⟨2,24,30⟩+1`) are the ONLY shapes in range
   where the *emitter* sits above the formula — and were still below every
   published catalog. (Closed at the catalog level by issue #7, 2026-09-30:
   explicit ℤ schemes at the exact formula, via the 12-product reusable set.)
3. Dense ±1 Lemma-1 rows can **never** be unimodular for `n−p ≥ 2` (every
   2×2 ±1-minor is even) — the origin of the denominators 2/4/8/16 in the
   first emission; zeros are necessary, and the Euclidean comb construction
   is the clean fix.
4. Verification: quadruple gate at emission + independent cross-language
   certificate (`verify_constructed_independent.py`): **465/465 exact**, the
   integer claim checked from actual coefficients.

**Leftovers (resume here):**

- **Task #9 (CLOSED externally, 2026-09-30)** — the six `g ≥ 6` shapes. The
  frontier question — does a `(3,3,bridge)` completion exist over the
  **12-product** reusable set? — is answered **yes** by issue #7's
  construction (virtual-row bridge; exhaustive F₂/F₃/F₅ per the contributor:
  never with 6 or 9 reusables, always with 12). Our quadratic-system
  elimination was never run; the theorem (9-product scope) is untouched.
  What remains *ours*: port the construction (and the reported rational
  pencil-of-tori variant, formula-exact for all `2 ≤ p ≤ n ≤ 32`) into the
  emitter so the family is emitted + certified by us beyond the six
  contributed shapes; re-run the closure sweep on shapes consuming them. See
  CONSTRUCTIVE_METHOD.md §"Where we stopped" and
  `references/contributions/issue7-adriano-2026-hk-task9.md`.
- **Lower bounds** — parked by explicit decision (HK proved optimality only
  for ⟨2,2,n⟩ and ⟨2,3,3⟩; everything we emit is upper-bound tier).
- `references/hopcroftkerr1971/README.md` is the *historical* impossibility
  document — kept with a status banner; its narrow-S theorem is correct but
  was sidestepped (scheduling) and then reversed (true reusable set).
- The 2026-06 working tree (schemes, manifest, paper §hk71, tests) was
  pending commit when we stopped.

**Historical notes** (superseded but kept for provenance): the first
"joint matrix rank ⇒ impossible" claim was wrong (the right invariant is
residual *tensor* rank); old task numbers #112/#113 correspond to what
became tasks #9/#10.

**See also.** `references/hopcroftkerr1971/` (impossibility + sympy
artifacts), `paper/sections/hk71.tex` (the paper section),
`src/main/java/eu/solven/matmul/papers/hopcroftkerr1971/` (the emitters).
