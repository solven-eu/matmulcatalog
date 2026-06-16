# Serendipity — savings-based bud scoring (open)

**Mechanism.** A serendipitous product reuses a base's *bud* (a rank-one term
shared across products) over an *enlarged* inner: a U-bud of multiplicity k turns
inner ⟨n₂,m₂,p₂⟩ into ⟨n₂,m₂,k·p₂⟩, which can be sub-additive vs k separate copies.
The search (`SerendipitousSearch.bestFor`) enumerates every divisor-shape as base
(both Kronecker orders), bud-decomposes the base, and predicts rank via
`predictRank`. (See `references/SERENDIPITOUS_PARTIAL_PRODUCT.md`.)

**Open item — a *formula* bud score.** Today `budScore` is a structural proxy and
base generation is gated on `budScore(b) > 0`. We want a **BudEffectTable**: the
realizable saving `saving(type,k,inner) = k·R(inner) − R(enlarge(type,k,inner))`
tabulated from the catalog, so bases are ranked by *actual* serendipity value, not
a structural heuristic. This also surfaces when a bud-rich base is missed because
its buds sit on the wrong axis (e.g. the prime-perfect-match case: only a bud on
the matched axis wins).

**Next.** Build the saving table; replace the structural `budScore`; fold
serendipity into the closure cost engine (#159, #196).
