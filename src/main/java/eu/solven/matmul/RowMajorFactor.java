package eu.solven.matmul;

/**
 * Compact row-major (CSR) snapshot of a {@link FactorMatrix}, for hot
 * <em>by-row</em> access without densifying.
 *
 * <p>The sparse factor backing is column-major (fast per-product traversal —
 * {@code forEachInColumn} / {@code dotColumn}). Some consumers instead need to
 * read a whole <em>row</em> (one flattened matrix entry across all {@code r}
 * products) — notably the sampled verifier, which evaluates
 * {@code Σ_k U[a][k]·V[b][k]·W[c][k]} at random {@code (a,b,c)}. Densifying to
 * {@code double[rows][r]} for that costs {@code rows·r·8} bytes — ~137&nbsp;MB
 * per factor (~413&nbsp;MB for U/V/W) at ⟨32,32,32⟩, almost all zeros. This CSR
 * stores only the non-zeros: ~{@code nnz·12} bytes (a factor is ~3–5&nbsp;%
 * dense, so single-digit MB), an order-of-magnitude memory win.</p>
 *
 * <p>Within each row the column indices are <strong>ascending</strong> (both
 * {@link SparseFactorMatrix} and {@link DenseFactorMatrix} emit
 * {@link FactorMatrix#forEachNonZero} with columns ascending for a fixed row),
 * so two or three rows can be combined by a linear merge — see
 * {@link #triProduct}.</p>
 */
public final class RowMajorFactor {

	private final int[] rowStart;  // length rows+1; row r spans [rowStart[r], rowStart[r+1])
	private final int[] col;       // length nnz; ascending within each row
	private final double[] val;    // length nnz

	private RowMajorFactor(int[] rowStart, int[] col, double[] val) {
		this.rowStart = rowStart;
		this.col = col;
		this.val = val;
	}

	/** Build a CSR snapshot in {@code O(nnz)} time and memory, reading only the
	 *  non-zeros of {@code fm} (never materialises a dense grid). */
	public static RowMajorFactor of(FactorMatrix fm) {
		int rows = fm.rows();
		int nnz = fm.nonZeros();
		int[] rowStart = new int[rows + 1];
		// Pass 1: count non-zeros per row (offset by 1 for the prefix-sum below).
		fm.forEachNonZero((r, c, v) -> rowStart[r + 1]++);
		for (int i = 0; i < rows; i++) {
			rowStart[i + 1] += rowStart[i];
		}
		int[] col = new int[nnz];
		double[] val = new double[nnz];
		int[] cursor = rowStart.clone();
		// Pass 2: scatter into CSR. Column-ascending-within-row is preserved
		// because both backings visit a fixed row's columns in ascending order.
		fm.forEachNonZero((r, c, v) -> {
			int idx = cursor[r]++;
			col[idx] = c;
			val[idx] = v;
		});
		return new RowMajorFactor(rowStart, col, val);
	}

	public int rows() {
		return rowStart.length - 1;
	}

	/** Sum of {@code uRow·vRow·wRow} over the shared product-columns (k) of the
	 *  three given rows — i.e. {@code Σ_k U[a][k]·V[b][k]·W[c][k]} where the
	 *  three CSRs hold rows {@code a}, {@code b}, {@code c} respectively. Linear
	 *  three-pointer merge over ascending column indices; only columns where all
	 *  three rows are non-zero contribute (zeros elsewhere add nothing). */
	public static double triProduct(RowMajorFactor u, int a, RowMajorFactor v, int b, RowMajorFactor w, int c) {
		int pa = u.rowStart[a], ea = u.rowStart[a + 1];
		int pb = v.rowStart[b], eb = v.rowStart[b + 1];
		int pc = w.rowStart[c], ec = w.rowStart[c + 1];
		double sum = 0.0;
		while (pa < ea && pb < eb && pc < ec) {
			int ka = u.col[pa], kb = v.col[pb], kc = w.col[pc];
			if (ka == kb && kb == kc) {
				sum += u.val[pa] * v.val[pb] * w.val[pc];
				pa++;
				pb++;
				pc++;
			} else {
				int kmax = Math.max(ka, Math.max(kb, kc));
				if (ka < kmax) pa++;
				if (kb < kmax) pb++;
				if (kc < kmax) pc++;
			}
		}
		return sum;
	}

	/** F₂ / bitwise variant of {@link #triProduct}: XOR-accumulates
	 *  {@code round(U[a][k]) & round(V[b][k]) & round(W[c][k])} over shared
	 *  columns. Matches the dense F₂ sampled-residual semantics exactly (a
	 *  missing entry is 0, and {@code 0 & x & y == 0} contributes nothing). */
	public static int triAndXor(RowMajorFactor u, int a, RowMajorFactor v, int b, RowMajorFactor w, int c) {
		int pa = u.rowStart[a], ea = u.rowStart[a + 1];
		int pb = v.rowStart[b], eb = v.rowStart[b + 1];
		int pc = w.rowStart[c], ec = w.rowStart[c + 1];
		int acc = 0;
		while (pa < ea && pb < eb && pc < ec) {
			int ka = u.col[pa], kb = v.col[pb], kc = w.col[pc];
			if (ka == kb && kb == kc) {
				int ui = (int) Math.round(u.val[pa]);
				int vi = (int) Math.round(v.val[pb]);
				int wi = (int) Math.round(w.val[pc]);
				// Mask to bit 0: for coefficients outside {0,1} (e.g. -1) the upper
				// bits of ui&vi&wi are garbage and would pollute the XOR accumulator,
				// falsely rejecting mod-2-exact Z schemes. Bit 0 alone is the parity.
				acc ^= (ui & vi & wi) & 1;
				pa++;
				pb++;
				pc++;
			} else {
				int kmax = Math.max(ka, Math.max(kb, kc));
				if (ka < kmax) pa++;
				if (kb < kmax) pb++;
				if (kc < kmax) pc++;
			}
		}
		return acc;
	}
}
