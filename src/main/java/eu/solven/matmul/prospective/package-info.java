/**
 * <b>Prospective (non-core) mechanics.</b> Code in this package is exploratory —
 * it is NOT part of the catalog's core construction/search pipeline and is not
 * wired into materialisation. It is kept for investigation and future work, and
 * its outputs must not be treated as constructive bounds without independent
 * verification.
 *
 * <p>Current residents:</p>
 * <ul>
 *   <li>{@link eu.solven.matmul.prospective.DisjointSumSearch} — a τ-theorem /
 *       disjoint-MM-sum cover search (Pan 1980 / Schönhage 1981) with an optional
 *       same-shape trilinear-aggregation post-pass. It enforces only loose
 *       area-cover constraints, so its predicted ranks are an OVER-OPTIMISTIC
 *       heuristic (e.g. ⟨4,4,4⟩→28, ⟨17,17,17⟩→2395 are unrealisable) — a
 *       shape-exploration tool, not a rank predictor. The project deliberately
 *       does not leverage the τ-theorem as a strategy it runs (no catalogued
 *       scheme uses disjoint-sum lineage); this package location makes that
 *       non-core, speculative status explicit.</li>
 * </ul>
 */
package eu.solven.matmul.prospective;
