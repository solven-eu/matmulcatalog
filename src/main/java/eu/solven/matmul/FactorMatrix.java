package eu.solven.matmul;

/**
 * Read-only view of one factor matrix ({@code U}, {@code V} or {@code W}) of a
 * {@link NonCubicBilinearAlgorithm}, abstracting over dense vs sparse storage.
 *
 * <p>Layout convention (matching the historical {@code double[][]}): a
 * <strong>row</strong> is a flattened matrix entry (e.g. {@code i·m+j} for
 * {@code U}), a <strong>column</strong> is a product {@code k ∈ [0,r)}. So
 * {@code get(row, col)} is the coefficient of entry {@code row} in product
 * {@code col}.</p>
 *
 * <p>Factor matrices are immutable after construction (verified: no in-place
 * writes anywhere in the codebase), so this interface is read-only. The visitor
 * methods let hot loops traverse only the non-zeros of a sparse backing instead
 * of scanning the full dense grid; {@link #get} stays available for the
 * scattered random-access sites.</p>
 */
public interface FactorMatrix {

	/** Number of rows = flattened matrix entries (e.g. {@code n·m} for U). */
	int rows();

	/** Number of columns = products {@code r}. */
	int cols();

	/** Coefficient of entry {@code row} in product {@code col}. */
	double get(int row, int col);

	/** Visit each non-zero {@code (row, value)} of a single product {@code col}
	 *  (its support). This is the per-product hot path. */
	void forEachInColumn(int col, ColumnConsumer consumer);

	/** Visit every non-zero {@code (row, col, value)}. Iteration order is
	 *  unspecified (column-major for the sparse backing). */
	void forEachNonZero(EntryConsumer consumer);

	/** Number of stored non-zeros (e.g. for addition-count style queries). */
	int nonZeros();

	/** Sparse dot-product of column {@code col} with a dense vector:
	 *  {@code Σ over non-zeros (row,val) of col: val · dense[row]}. Allocation-free
	 *  — the primitive for hot per-product accumulation (e.g. {@code α_k}). */
	double dotColumn(int col, double[] dense);

	/** Sparse AXPY into a dense target: for each non-zero {@code (row,val)} of
	 *  column {@code col}, {@code target[row] += scalar · val}. Allocation-free —
	 *  the primitive for hot output accumulation (e.g. {@code C += W[:,k]·γ_k}). */
	void axpyColumn(int col, double scalar, double[] target);

	@FunctionalInterface
	interface ColumnConsumer {
		void accept(int row, double value);
	}

	@FunctionalInterface
	interface EntryConsumer {
		void accept(int row, int col, double value);
	}
}
