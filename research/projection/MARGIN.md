# Projection — margin (μ) and the linear-projection gap (open)

**Projection margin μ.** For axis `a`, product `k` is `i`-private iff it is
confined to input-block `i` (U) **or** output-block `i` (W); ρₐ(i) counts
`i`-private products and μ = max over axes/indices. A single-index projection
drops `R → R − μ`. A higher-rank but *more structured* parent can project better
than the flat rank-best one (the SZ-vs-PanTA ⟨30,30,30⟩ example). μ is stamped
per scheme (`projection_margin`) and surfaced in the SPA. (Paper §projmargin.)

**Open item — stronger projection.** `ProjectionSearch.bestFor` enumerates only
**coordinate-subset** drops (keep 27 of 28 raw indices). FMM's `[[1,0],[0]]`-style
projections are **linear**: they fold an index into a linear combination of others
before dropping, collapsing more product supports (more dead code). Example gap:
our best single-index drop kills μ=114 where FMM's linear projection kills 143
(10442 vs 10413). Matching it needs a **linear-projection operator**, strictly
stronger than index selection.

**Next.** Add a linear-projection search (rank-1 row/col fold + DCE), score by the
realized product-support collapse, and wire it into the downward closure wave.
