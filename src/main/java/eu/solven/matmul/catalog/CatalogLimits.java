package eu.solven.matmul.catalog;

/**
 * Single source of truth for upper-bound enumeration limits across the
 * catalog. Any algorithm that iterates {@code 2..MAX_DIM} (formula
 * families, derived bounds, gap audits) must reference these constants
 * rather than hard-coding a number — so a single edit retunes the
 * entire pipeline.
 *
 * <p>{@link #MAX_DIM} is the inclusive upper bound on individual matrix
 * dimensions {@code n, m, p} for all "go wide" enumerations. Currently
 * 32 — chosen to align with the project's "small-matmul catalog up to
 * {@code ⟨32, 32, 32⟩}" target.</p>
 */
public final class CatalogLimits {

	private CatalogLimits() {}

	/**
	 * Inclusive upper bound on matrix-dimension enumeration for cubic
	 * shapes ⟨n, n, n⟩. Currently {@code 32}.
	 */
	public static final int MAX_DIM = 32;

	/** Alias of {@link #MAX_DIM} for callers that want a name disambiguating from non-cubic. */
	public static final int MAX_CUBIC_DIM = MAX_DIM;

	/**
	 * Inclusive upper bound on each axis for <strong>non-cubic</strong>
	 * shape sweeps. Currently {@code 8} — the small slice where direct
	 * factor matrices on disk dominate the ranks and the cost is bounded.
	 * Beyond {@code 8} the non-cubic shape space gets dense and is mostly
	 * dominated by Kronecker decomposition of larger cubic shapes, so
	 * routine sweeps cap there to stay tractable.
	 */
	public static final int MAX_NONCUBIC_DIM = 8;

	// ───────────────────────────────────────────────────────────────────────
	// On-disk JSON-format customisation (single home; #162 / task #187)
	// ───────────────────────────────────────────────────────────────────────

	/**
	 * Maximum maxDim (= max axis) for which <em>derived</em> schemes get their
	 * factor matrices materialised on disk. Above this, derived schemes live as
	 * lineage-only stubs (reproduced on demand by
	 * {@link eu.solven.matmul.search.LineageReplayer}); atoms are always
	 * materialised regardless. Currently {@code 16}. Canonical home for what was
	 * historically {@code CatalogPolicy.MATERIALISE_MAX_DIM} (kept as a
	 * delegating alias). See {@link CatalogPolicy} for the two-tier rationale.
	 */
	public static final int MATERIALISE_MAX_DIM = 16;

	/**
	 * Switching threshold between the dense {@code u}/{@code v}/{@code w}
	 * encoding and the compact sparse {@code *_sparse} encoding: dense for
	 * {@code max(n,m,p) < this}, sparse at or above.
	 *
	 * <p>Set to {@code 9} (2026-06): <strong>dense up to {@code 2³ = 8}, sparse
	 * from {@code 9}</strong>. The empirical crossover where sparse first beats
	 * dense is dim&nbsp;5 ({@code CompareDenseVsSparseSizes}: sp/dn 1.24 at
	 * dim&nbsp;4, 0.85 at dim&nbsp;5, 0.49 at dim&nbsp;8, 0.44 at dim&nbsp;9),
	 * but dims 5–8 stay human-readable as dense matrices on still-modest files
	 * (&le;223&nbsp;KB); from dim&nbsp;9 sparse saves &ge;55% on files that grow
	 * past 450&nbsp;KB. {@code 2³} is the natural cut: it is exactly one
	 * Strassen recursion level above the {@code ⟨2,2,2⟩} base, so every shape
	 * that is "one split from a 2×2 leaf" stays dense and readable. (Was
	 * {@code 16}.)</p>
	 */
	public static final int SPARSE_DIM_THRESHOLD = 9;
}
