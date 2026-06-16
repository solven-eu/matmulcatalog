# SOTA conventions for the paper's comparison tables

The frontier-closure catalog has more entries than fit into a clean
comparison table; the table needs deterministic rules for which entry
to report per `(shape, field, commutativity)`. This file documents
the rules so the regeneration tooling can be deterministic.

## The working definition

> **SOTA** = the lowest known rank for `(⟨n,m,p⟩, K, commutative)`
> witnessed by either an explicit catalog entry, a cited bound, or a
> derived bound.

The cell shows the rank, the **earliest** known attribution (from
`attribution_for_rank`), and a dagger (`†`) if the cell is a cited
bound without an explicit on-disk scheme.

## Per-axis decisions

### Field promotion

A scheme valid over a wider field is valid over every narrower one.
We populate cells transitively:

- A Z-scheme is also a Q-scheme, an R-scheme, and a C-scheme.
- A Q-scheme is also an R-scheme and a C-scheme.
- An R-scheme is also a C-scheme.
- **F₂ does not auto-promote from Z** (cancellations land on
  reduction; a Z-scheme's rank may not be its mod-2 rank).

The Z-to-Q-to-R-to-C promotion is automatic: a scheme tagged Z
populates the Q/R/C cells with the same rank value if no strictly
better one exists for that cell.

The F₂ row is populated from F₂-native schemes only (AlphaTensor,
F₂-SAT results, etc.).

### Commutative-vs-non-commutative

These are separate columns; we never let a commutative result fill
a non-commutative cell. The table for `c-K` reports only schemes
tagged `"commutative": true` in the JSON; the `nc-K` table reports
only schemes tagged `"commutative": false` (or absent).

A commutative scheme is **not** a valid non-commutative scheme for
the same rank: it relies on `ab = ba` which fails for matrix-block
recursion. Hence no cross-population.

### Border rank

Border-rank bounds (Bini-type) live in a separate
`docs/border-bounds.json` and a separate appendix table. They do
**not** populate the main comparison tables.

### Cited bounds (the dagger)

A cited bound is included in the cell if no on-disk scheme is
strictly better. The cell carries a dagger and a footnote referring
to the source.

When a cited bound is later materialised, the dagger disappears
automatically (next regeneration).

### Discovery status

When the entry has `"discovery": true`, the attribution column shows
the importing source. When `"discovery": false`, the attribution
column shows `attribution_for_rank` (the earliest known source); the
importing source goes into a footnote.

When `"discovery": "TBD"`, the cell uses `attribution_for_rank` if
present, otherwise the importing source with a daggered footnote
flagging the unverified discovery status.

### Verification status

Schemes that have failed `Verifier.passesRandomMatmulSpotCheck` are
excluded from the tables -- they are still on disk under a
`broken-` prefix but do not populate.

A scheme that has not been spot-checked since the last lineage
change is treated as verified for the table; the catalog-wide
verification sweep is responsible for catching divergence.

## Open questions for the regeneration tooling

1. When two schemes from different fields tie at the same rank, do
   we prefer the wider-field one (lifts to more downstream uses) or
   the narrower-field one (more specific)? Working answer: prefer
   the wider one, so a `Z` entry beats an `R` entry at the same
   rank when populating the R column.

2. Commutative tables for C and F₂ are sparse. Do we include empty
   sub-tables for completeness or omit? Working answer: omit if
   zero entries; include in an appendix.

3. Per-format border-rank entries (Bini ⟨2,2,2⟩ border rank 5) --
   do they go into a per-field border table or a single combined
   border table? Working answer: per-field, mirroring the main
   tables.

4. When a hand-crafted scheme and a composed scheme tie at the same
   rank but have distinct lineages, the cell should still pick one
   for display. Working answer: prefer the **earlier** of the two
   (oldest `attribution_for_rank`); the composed entry's existence
   shows up in a footnote.
