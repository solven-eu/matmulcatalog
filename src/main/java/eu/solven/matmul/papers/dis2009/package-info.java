/**
 * Drevet–Islam–Schost 2009 ({@code DIS09}) constructions.
 *
 * <p>Contents currently scoped to the trilinear-aggregation (TA) cubic
 * matmul algorithm: {@link eu.solven.matmul.papers.dis2009.PanTrilinearAggregation}
 * exposes both the closed-form rank bound
 * {@code (n³+12n²+11n)/3} (even) /
 * {@code (n³+15n²+14n−6)/3} (odd) and the explicit factor-matrix
 * constructor that yields a {@code NonCubicBilinearAlgorithm} matching
 * the bound. Implementation is a port of Islam's Magma reference
 * ({@code references/islam2009/magma/TA.mgm}, recovered from the
 * Wayback Machine — see {@code references/typos.md} for provenance).</p>
 *
 * <p><strong>Naming</strong>: the technique is "Pan trilinear
 * aggregation", introduced by Pan in 1978/1980 and refined through
 * Pan/Sha 1992. The specific {@link eu.solven.matmul.papers.dis2009.PanTrilinearAggregation#build}
 * construction in this package is the variant from Islam's 2009 MSc
 * thesis (subsequently published as DIS09 §3 / Appendix Lemma 4): even
 * case improved over Pan 1992, odd case introduced for the first time.
 * The "Pan" prefix on the class name preserves the technique name; a
 * future {@code eu.solven.matmul.papers.pan1978} variant may host the
 * original Pan TA construction if/when we port that one too.</p>
 */
package eu.solven.matmul.papers.dis2009;
