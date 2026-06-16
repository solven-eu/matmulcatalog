package eu.solven.matmul.catalog;

/**
 * Project-wide catalog policy constants.
 *
 * <h2>The two-tier catalog model (#162)</h2>
 *
 * <p>The catalog is conceptually split into three layers:</p>
 *
 * <ol>
 *   <li><strong>Atoms</strong> ({@code src/main/resources/schemes/} files
 *       whose lineage is {@code Atom(ref)} — Strassen ⟨2,2,2⟩=7, Laderman
 *       ⟨3,3,3⟩=23, AT-F2 ⟨4,4,4⟩=47, etc.). Hand-curated; the source of
 *       truth. Always carry explicit factor matrices.</li>
 *   <li><strong>Lineage stubs</strong> (same directory, JSON files
 *       carrying only {@code {format, rank, additions, lineage_compact,
 *       source, attribution_for_rank}} — NO factor matrices on disk).
 *       Reproducibility comes from {@link
 *       eu.solven.matmul.search.LineageReplayer} replaying the lineage
 *       against the atom set.</li>
 *   <li><strong>Materialised</strong> — explicit factor matrices for a
 *       lineage stub. Either inlined in the JSON (small shapes) or
 *       deferred to runtime replay (large shapes). Whether a stub gets
 *       materialised on disk is governed by {@link #MATERIALISE_MAX_DIM}
 *       below.</li>
 * </ol>
 *
 * <h2>Why bother</h2>
 *
 * <p>Above {@link #MATERIALISE_MAX_DIM} the bilinear factor matrices grow
 * to multi-megabyte JSON files that bloat the catalog and slow every
 * checkout / clone / load. Lineage stubs are tiny (~100 bytes each) and
 * fully reproducible. The trade-off is one extra replay step at load
 * time for above-cap shapes, which is much cheaper than carrying the
 * data.</p>
 */
public final class CatalogPolicy {

	private CatalogPolicy() {}

	/**
	 * Maximum maxDim (= max axis) for which derived schemes get their
	 * factor matrices materialised on disk. Above this dim, derived
	 * schemes live as lineage-only stubs and the materialisation step
	 * is deferred to runtime replay via
	 * {@link eu.solven.matmul.search.LineageReplayer}.
	 *
	 * <p>Default {@code 16}: keeps the on-disk catalog roughly the
	 * current size while letting larger shapes (⟨17,17,17⟩, ⟨32,32,32⟩,
	 * etc.) stay as stubs.</p>
	 *
	 * <p>Atoms are unconditionally materialised (they ARE the
	 * factor-matrix data the stubs reference) — this cap applies only
	 * to derived schemes, not atoms.</p>
	 *
	 * <p><strong>Canonical value lives in
	 * {@link CatalogLimits#MATERIALISE_MAX_DIM}</strong> (single home for all
	 * on-disk JSON-format constants); this is a delegating alias kept for the
	 * existing call sites and the rationale documented above.</p>
	 */
	public static final int MATERIALISE_MAX_DIM = CatalogLimits.MATERIALISE_MAX_DIM;
}
