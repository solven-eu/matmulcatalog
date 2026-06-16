package eu.solven.matmul;

/**
 * Column-major (CSC-style) sparse {@link FactorMatrix}. For each column (product)
 * {@code k} it stores the sorted row indices and values of that product's
 * non-zeros. Chosen because the hot loops are per-product ({@code dotColumn} /
 * {@code axpyColumn} / {@code forEachInColumn}); a factor matrix is typically
 * ~3–5 % dense, so this is ~15–25× smaller than the dense {@code double[][]} and
 * the per-product traversals touch only non-zeros.
 *
 * <p>Immutable; built once from a dense {@code double[][]} via {@link #fromDense}.</p>
 */
public final class SparseFactorMatrix implements FactorMatrix {

	private final int rows;
	private final int cols;
	private final int[][] rowIdx;   // rowIdx[col] = sorted rows of non-zeros in column col
	private final double[][] vals;  // vals[col]   = matching coefficients
	private final int nnz;

	private SparseFactorMatrix(int rows, int cols, int[][] rowIdx, double[][] vals, int nnz) {
		this.rows = rows;
		this.cols = cols;
		this.rowIdx = rowIdx;
		this.vals = vals;
		this.nnz = nnz;
	}

	/** Compress a dense {@code [rows][cols]} matrix to column-major sparse. */
	public static SparseFactorMatrix fromDense(double[][] dense) {
		int rows = dense.length;
		int cols = rows == 0 ? 0 : dense[0].length;
		int[][] rowIdx = new int[cols][];
		double[][] vals = new double[cols][];
		int nnz = 0;
		int[] tmpIdx = new int[rows];
		double[] tmpVal = new double[rows];
		for (int col = 0; col < cols; col++) {
			int c = 0;
			for (int row = 0; row < rows; row++) {
				double v = dense[row][col];
				if (v != 0.0) {
					tmpIdx[c] = row;
					tmpVal[c] = v;
					c++;
				}
			}
			rowIdx[col] = java.util.Arrays.copyOf(tmpIdx, c);
			vals[col] = java.util.Arrays.copyOf(tmpVal, c);
			nnz += c;
		}
		return new SparseFactorMatrix(rows, cols, rowIdx, vals, nnz);
	}

	/**
	 * Build column-major from per-column COO arrays {@code colRows[c]} /
	 * {@code colVals[c]} (the non-zeros of product {@code c}). Each column is
	 * sorted ascending by row in place — the caller may pass them unsorted, as the
	 * Kronecker builder does (the composite row index is not monotone in the
	 * factors' row order). No dense materialisation.
	 */
	public static SparseFactorMatrix fromColumns(int rows, int[][] colRows, double[][] colVals) {
		int cols = colRows.length;
		int nnz = 0;
		for (int c = 0; c < cols; c++) {
			int[] ri = colRows[c];
			double[] rv = colVals[c];
			int len = ri.length;
			for (int i = 1; i < len; i++) {  // insertion sort (columns are short)
				int ki = ri[i];
				double kv = rv[i];
				int j = i - 1;
				while (j >= 0 && ri[j] > ki) {
					ri[j + 1] = ri[j];
					rv[j + 1] = rv[j];
					j--;
				}
				ri[j + 1] = ki;
				rv[j + 1] = kv;
			}
			nnz += len;
		}
		return new SparseFactorMatrix(rows, cols, colRows, colVals, nnz);
	}

	@Override
	public int rows() {
		return rows;
	}

	@Override
	public int cols() {
		return cols;
	}

	@Override
	public double get(int row, int col) {
		int[] idx = rowIdx[col];
		int lo = 0, hi = idx.length - 1;
		while (lo <= hi) {
			int mid = (lo + hi) >>> 1;
			int v = idx[mid];
			if (v < row) lo = mid + 1;
			else if (v > row) hi = mid - 1;
			else return vals[col][mid];
		}
		return 0.0;
	}

	@Override
	public void forEachInColumn(int col, ColumnConsumer consumer) {
		int[] idx = rowIdx[col];
		double[] val = vals[col];
		for (int i = 0; i < idx.length; i++) {
			consumer.accept(idx[i], val[i]);
		}
	}

	@Override
	public void forEachNonZero(EntryConsumer consumer) {
		for (int col = 0; col < cols; col++) {
			int[] idx = rowIdx[col];
			double[] val = vals[col];
			for (int i = 0; i < idx.length; i++) {
				consumer.accept(idx[i], col, val[i]);
			}
		}
	}

	@Override
	public double dotColumn(int col, double[] dense) {
		int[] idx = rowIdx[col];
		double[] val = vals[col];
		double s = 0.0;
		for (int i = 0; i < idx.length; i++) {
			s += val[i] * dense[idx[i]];
		}
		return s;
	}

	@Override
	public void axpyColumn(int col, double scalar, double[] target) {
		int[] idx = rowIdx[col];
		double[] val = vals[col];
		for (int i = 0; i < idx.length; i++) {
			target[idx[i]] += scalar * val[i];
		}
	}

	@Override
	public int nonZeros() {
		return nnz;
	}

	/**
	 * A copy with every row index {@code r} relabelled to {@code rowMap(r)} and the
	 * row count set to {@code newRows}; columns and values are unchanged. {@code rowMap}
	 * MUST be injective over the stored rows (it encodes an index permutation/transpose,
	 * as used by {@link NonCubicBilinearAlgorithm#cyclicShift()} /
	 * {@link NonCubicBilinearAlgorithm#transpose()}). Fully sparse: {@code O(nnz)} plus a
	 * small per-column re-sort to restore ascending row order. No dense materialisation.
	 */
	public SparseFactorMatrix mapRows(int newRows, java.util.function.IntUnaryOperator rowMap) {
		int[][] nidx = new int[cols][];
		double[][] nval = new double[cols][];
		for (int col = 0; col < cols; col++) {
			int[] idx = rowIdx[col];
			double[] val = vals[col];
			int len = idx.length;
			int[] ni = new int[len];
			double[] nv = new double[len];
			for (int i = 0; i < len; i++) {
				ni[i] = rowMap.applyAsInt(idx[i]);
				nv[i] = val[i];
			}
			// Insertion sort (columns are short); keep values aligned with row indices.
			for (int i = 1; i < len; i++) {
				int ki = ni[i];
				double kv = nv[i];
				int j = i - 1;
				while (j >= 0 && ni[j] > ki) {
					ni[j + 1] = ni[j];
					nv[j + 1] = nv[j];
					j--;
				}
				ni[j + 1] = ki;
				nv[j + 1] = kv;
			}
			nidx[col] = ni;
			nval[col] = nv;
		}
		return new SparseFactorMatrix(newRows, cols, nidx, nval, nnz);
	}

	/** Materialise a fresh dense {@code [rows][cols]} copy (for {@code denseU()}). */
	public double[][] toDense() {
		double[][] dense = new double[rows][cols];
		for (int col = 0; col < cols; col++) {
			int[] idx = rowIdx[col];
			double[] val = vals[col];
			for (int i = 0; i < idx.length; i++) {
				dense[idx[i]][col] = val[i];
			}
		}
		return dense;
	}
}
