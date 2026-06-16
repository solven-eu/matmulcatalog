package eu.solven.matmul.catalog;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Memory-compact backing for a {@link NonCubicBilinearAlgorithm}, used
 * by {@link FieldAwareLookup}'s parse cache. Most catalog schemes hold
 * tiny coefficient sets <em>and</em> substantial zero density — storing
 * the {@code double[][]} factors uses 8 bytes per entry where a dense
 * ternary or a sparse list would use a fraction.
 *
 * <p>Per-matrix encoding (auto-selected by {@link #of} based on which
 * is smallest):</p>
 * <ul>
 *   <li>{@code SPARSE_INT} — only small-integer non-zeros, low density;
 *       4-byte flat index + 1-byte byte-coef per nonzero. Best for
 *       outputs of Kronecker / concat composition.</li>
 *   <li>{@code SPARSE_DOUBLE} — same shape but coefs need full
 *       precision (rationals / irrationals); 4-byte index + 8-byte value.</li>
 *   <li>{@code SPARSE_DICT} — entries that aren't byte-integers but draw
 *       from a small set of unique non-zero values (rationals like
 *       ±1/2, ±1/4 …). Build a {@code double[]} dictionary, then store
 *       each non-zero as a {@code (flat-index, dict-index)} pair where
 *       the dict-indices are packed into a {@code long[]} at 1/2/4/8 bits
 *       depending on {@code |dict|}. For a typical Q-rational scheme with
 *       4-8 unique values, this is roughly half of {@code SPARSE_DOUBLE}.</li>
 *   <li>{@code TERNARY} — all entries ∈ {-1, 0, 1}, packed 2 bits each
 *       into a {@code long[]} (32× over raw double; beats sparse when
 *       density &gt; ~25%).</li>
 *   <li>{@code DOUBLE} — raw {@code double[][]}, last-resort fallback.</li>
 * </ul>
 *
 * <p>{@link #expand} materialises a fresh {@link NonCubicBilinearAlgorithm}
 * on demand. Callers expand → use → drop; the compact backing in the
 * cache stays small.</p>
 */
public final class CompactScheme {

	public final int n, m, p, r;
	private final Matrix uMat, vMat, wMat;

	private static final byte TERNARY = 0;
	private static final byte SPARSE_INT = 1;
	private static final byte SPARSE_DICT = 2;
	private static final byte SPARSE_DOUBLE = 3;
	private static final byte DOUBLE = 4;

	private record Matrix(byte kind, int rows, int cols,
			int[] idx, byte[] valByte, double[] valDouble,
			long[] packed, double[][] dense,
			double[] dict, int dictBits) {

		long byteSize() {
			return switch (kind) {
				case TERNARY -> packed.length * 8L;
				case SPARSE_INT -> idx.length * 5L;
				case SPARSE_DICT -> idx.length * 4L + packed.length * 8L + dict.length * 8L;
				case SPARSE_DOUBLE -> idx.length * 12L;
				default -> (long) rows * cols * 8L;
			};
		}

		double[][] expand() {
			double[][] out = new double[rows][cols];
			switch (kind) {
				case TERNARY -> {
					int bp = 0;
					for (int i = 0; i < rows; i++) {
						double[] row = out[i];
						for (int j = 0; j < cols; j++) {
							int code = (int) ((packed[bp >>> 5] >>> ((bp & 31) * 2)) & 0x3L);
							if (code == 1) row[j] = 1.0;
							else if (code == 2) row[j] = -1.0;
							bp++;
						}
					}
				}
				case SPARSE_DICT -> {
					long mask = (1L << dictBits) - 1L;
					int slotsPerLong = 64 / dictBits;
					for (int k = 0; k < idx.length; k++) {
						int code = (int) ((packed[k / slotsPerLong] >>> ((k % slotsPerLong) * dictBits)) & mask);
						out[idx[k] / cols][idx[k] % cols] = dict[code];
					}
				}
				case SPARSE_INT -> {
					for (int k = 0; k < idx.length; k++) {
						out[idx[k] / cols][idx[k] % cols] = valByte[k];
					}
				}
				case SPARSE_DOUBLE -> {
					for (int k = 0; k < idx.length; k++) {
						out[idx[k] / cols][idx[k] % cols] = valDouble[k];
					}
				}
				default -> {
					for (int i = 0; i < rows; i++) {
						System.arraycopy(dense[i], 0, out[i], 0, cols);
					}
				}
			}
			return out;
		}
	}

	private CompactScheme(int n, int m, int p, int r, Matrix u, Matrix v, Matrix w) {
		this.n = n; this.m = m; this.p = p; this.r = r;
		this.uMat = u; this.vMat = v; this.wMat = w;
	}

	public static CompactScheme of(NonCubicBilinearAlgorithm alg) {
		return new CompactScheme(alg.n, alg.m, alg.p, alg.r,
				bestPack(alg.denseU()), bestPack(alg.denseV()), bestPack(alg.denseW()));
	}

	public long byteSize() { return uMat.byteSize() + vMat.byteSize() + wMat.byteSize(); }

	public NonCubicBilinearAlgorithm expand() {
		return new NonCubicBilinearAlgorithm(n, m, p, uMat.expand(), vMat.expand(), wMat.expand());
	}

	/** Cheapest encoding wins. Single pass tallies nnz + ternary + byte-int + dict. */
	private static Matrix bestPack(double[][] mat) {
		int rows = mat.length, cols = mat[0].length;
		long total = (long) rows * cols;
		int nnz = 0;
		boolean ternaryOk = true, byteIntOk = true;
		// Dictionary scan: track unique non-zero values (limit 64 for cheap detection).
		double[] dictScratch = new double[64];
		int dictSize = 0;
		boolean dictOverflowed = false;
		for (double[] row : mat) {
			for (double v : row) {
				if (v == 0.0) continue;
				nnz++;
				if (v != 1.0 && v != -1.0) ternaryOk = false;
				if (v != Math.rint(v) || v < -128 || v > 127) byteIntOk = false;
				if (!dictOverflowed) {
					boolean seen = false;
					for (int k = 0; k < dictSize; k++) {
						if (dictScratch[k] == v) { seen = true; break; }
					}
					if (!seen) {
						if (dictSize == dictScratch.length) dictOverflowed = true;
						else dictScratch[dictSize++] = v;
					}
				}
			}
		}

		long denseBytes = total * 8L;
		long ternaryBytes = (total + 31) / 32 * 8L;
		long sparseIntBytes = (long) nnz * 5L;
		long sparseDblBytes = (long) nnz * 12L;
		// Dict bytes = idx[] + packed[] + dict[]
		int dictBits = !dictOverflowed && dictSize > 0
				? (dictSize <= 2 ? 1 : dictSize <= 4 ? 2 : dictSize <= 16 ? 4 : 8)
				: 0;
		long sparseDictBytes = dictBits > 0
				? (long) nnz * 4L + ((long) nnz * dictBits + 63) / 64 * 8L + (long) dictSize * 8L
				: Long.MAX_VALUE;

		long best = denseBytes;
		byte kind = DOUBLE;
		if (ternaryOk && ternaryBytes < best) { best = ternaryBytes; kind = TERNARY; }
		if (byteIntOk && sparseIntBytes < best) { best = sparseIntBytes; kind = SPARSE_INT; }
		if (sparseDictBytes < best) { best = sparseDictBytes; kind = SPARSE_DICT; }
		if (sparseDblBytes < best) { kind = SPARSE_DOUBLE; }

		return switch (kind) {
			case TERNARY -> packTernary(mat, rows, cols);
			case SPARSE_INT -> packSparseInt(mat, rows, cols, nnz);
			case SPARSE_DICT -> packSparseDict(mat, rows, cols, nnz,
					java.util.Arrays.copyOf(dictScratch, dictSize), dictBits);
			case SPARSE_DOUBLE -> packSparseDouble(mat, rows, cols, nnz);
			default -> new Matrix(DOUBLE, rows, cols, null, null, null, null, mat, null, 0);
		};
	}

	private static Matrix packSparseDict(double[][] mat, int rows, int cols, int nnz,
			double[] dict, int dictBits) {
		int[] idx = new int[nnz];
		int slotsPerLong = 64 / dictBits;
		long mask = (1L << dictBits) - 1L;
		long[] packed = new long[(nnz + slotsPerLong - 1) / slotsPerLong];
		int k = 0;
		for (int i = 0; i < rows; i++) {
			double[] row = mat[i];
			for (int j = 0; j < cols; j++) {
				double v = row[j];
				if (v == 0.0) continue;
				idx[k] = i * cols + j;
				int dictIdx = -1;
				for (int d = 0; d < dict.length; d++) {
					if (dict[d] == v) { dictIdx = d; break; }
				}
				packed[k / slotsPerLong] |= ((long) dictIdx & mask) << ((k % slotsPerLong) * dictBits);
				k++;
			}
		}
		return new Matrix(SPARSE_DICT, rows, cols, idx, null, null, packed, null, dict, dictBits);
	}

	private static Matrix packTernary(double[][] mat, int rows, int cols) {
		long total = (long) rows * cols;
		long[] out = new long[(int) ((total + 31) / 32)];
		int bp = 0;
		for (double[] row : mat) {
			for (double v : row) {
				int code = (v == 0.0) ? 0 : (v == 1.0) ? 1 : 2;
				out[bp >>> 5] |= (long) code << ((bp & 31) * 2);
				bp++;
			}
		}
		return new Matrix(TERNARY, rows, cols, null, null, null, out, null, null, 0);
	}

	private static Matrix packSparseInt(double[][] mat, int rows, int cols, int nnz) {
		int[] idx = new int[nnz];
		byte[] val = new byte[nnz];
		int k = 0;
		for (int i = 0; i < rows; i++) {
			double[] row = mat[i];
			for (int j = 0; j < cols; j++) {
				if (row[j] == 0.0) continue;
				idx[k] = i * cols + j;
				val[k] = (byte) row[j];
				k++;
			}
		}
		return new Matrix(SPARSE_INT, rows, cols, idx, val, null, null, null, null, 0);
	}

	private static Matrix packSparseDouble(double[][] mat, int rows, int cols, int nnz) {
		int[] idx = new int[nnz];
		double[] val = new double[nnz];
		int k = 0;
		for (int i = 0; i < rows; i++) {
			double[] row = mat[i];
			for (int j = 0; j < cols; j++) {
				if (row[j] == 0.0) continue;
				idx[k] = i * cols + j;
				val[k] = row[j];
				k++;
			}
		}
		return new Matrix(SPARSE_DOUBLE, rows, cols, idx, null, val, null, null, null, 0);
	}
}
