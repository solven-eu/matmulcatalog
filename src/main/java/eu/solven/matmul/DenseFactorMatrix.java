package eu.solven.matmul;

/**
 * {@link FactorMatrix} backed by the historical dense {@code double[][]} layout
 * ({@code data[row][col]}). Used as the drop-in backing during the migration
 * away from direct array access; the sparse backing
 * ({@link SparseFactorMatrix}) replaces it once all consumers go through the
 * {@link FactorMatrix} API.
 */
public final class DenseFactorMatrix implements FactorMatrix {

	private final double[][] data; // [row][col]
	private final int cols;

	public DenseFactorMatrix(double[][] data) {
		this.data = data;
		this.cols = data.length == 0 ? 0 : data[0].length;
	}

	@Override
	public int rows() {
		return data.length;
	}

	@Override
	public int cols() {
		return cols;
	}

	@Override
	public double get(int row, int col) {
		return data[row][col];
	}

	@Override
	public void forEachInColumn(int col, ColumnConsumer consumer) {
		for (int row = 0; row < data.length; row++) {
			double v = data[row][col];
			if (v != 0.0) {
				consumer.accept(row, v);
			}
		}
	}

	@Override
	public void forEachNonZero(EntryConsumer consumer) {
		for (int row = 0; row < data.length; row++) {
			double[] r = data[row];
			for (int col = 0; col < r.length; col++) {
				if (r[col] != 0.0) {
					consumer.accept(row, col, r[col]);
				}
			}
		}
	}

	@Override
	public double dotColumn(int col, double[] dense) {
		double s = 0.0;
		for (int row = 0; row < data.length; row++) {
			double v = data[row][col];
			if (v != 0.0) {
				s += v * dense[row];
			}
		}
		return s;
	}

	@Override
	public void axpyColumn(int col, double scalar, double[] target) {
		for (int row = 0; row < data.length; row++) {
			double v = data[row][col];
			if (v != 0.0) {
				target[row] += scalar * v;
			}
		}
	}

	@Override
	public int nonZeros() {
		int nz = 0;
		for (double[] r : data) {
			for (double v : r) {
				if (v != 0.0) {
					nz++;
				}
			}
		}
		return nz;
	}
}
