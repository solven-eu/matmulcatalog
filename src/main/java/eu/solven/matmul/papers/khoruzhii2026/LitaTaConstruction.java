package eu.solven.matmul.papers.khoruzhii2026;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.math3.fraction.BigFraction;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Exact-rational (BigFraction) constructor for the Khoruzhii&ndash;Gel&szlig;&ndash;Pokutta
 * 2026 LITA cubic matmul schemes &mdash; a faithful port of the Maple generators
 * {@code KGP2026_odd.mpl} / {@code KGP2026_even.mpl}. Produces a
 * {@link NonCubicBilinearAlgorithm} for {@code <N,N,N>}, {@code N >= 19}.
 *
 * <p>Rank equals {@link LitaTrilinearAggregation#cubicRank(int)}. Factors are
 * computed exactly with {@link BigFraction} and converted to {@code double} only
 * at the very end. The odd branch reproduces the reference scheme up to a row
 * permutation; the even branch is a valid alternate decomposition at the same
 * rank (it computes the matmul exactly &mdash; see {@code eu.solven.matmul.verifiers.Verifier}).</p>
 *
 * <h3>Representation</h3>
 * A linear form is a {@code Map<Integer,BigFraction>} keyed by {@code row*N + col}
 * (0-based) into an {@code N x N} matrix. Each rank-1 product is a triple of forms
 * (one in A, one in B, one in X). Flatten: a U/A coefficient at {@code (i,j)} goes to
 * factor row {@code i*N+j}; a V/B coefficient at {@code (j,k)} to {@code j*N+k}; a W/X
 * coefficient at {@code (a,b)} to {@code b*N+a} (transpose-then-row-major, the
 * {@code sum_l U_l(A) V_l(B) W_l(X^T) = trace(A B X)} convention). Per-row scaling is
 * rebalanced across the three factors (even split of the rational content's prime
 * factorisation) so every stored denominator stays small.
 *
 * <p>Source: Kirill Khoruzhii, Patrick Gel&szlig;, Sebastian Pokutta,
 * <em>Local Improvements to Trilinear Aggregation</em>, 2026.</p>
 */
public final class LitaTaConstruction {

	private LitaTaConstruction() {}

	private static BigFraction f(int n) { return new BigFraction(n); }

	private static final BigFraction ZERO = BigFraction.ZERO;
	private static final BigFraction ONE = BigFraction.ONE;

	/** Build the LITA scheme for {@code <n,n,n>}, exact then double. */
	public static NonCubicBilinearAlgorithm build(int n) {
		if (n < LitaTrilinearAggregation.MIN_N) {
			throw new IllegalArgumentException(
					"LITA is defined only for N >= " + LitaTrilinearAggregation.MIN_N + ", got " + n);
		}
		return ((n & 1) == 1) ? new Odd(n).build() : new Even(n).build();
	}

	// ===================================================================
	// linear-form helpers (Map<Integer,BigFraction>)
	// ===================================================================
	private static void addInto(Map<Integer, BigFraction> dst, Map<Integer, BigFraction> src, BigFraction scale) {
		if (scale.getNumerator().signum() == 0) {
			return;
		}
		for (Map.Entry<Integer, BigFraction> e : src.entrySet()) {
			int k = e.getKey();
			BigFraction nv = dst.getOrDefault(k, ZERO).add(e.getValue().multiply(scale));
			if (nv.getNumerator().signum() == 0) {
				dst.remove(k);
			} else {
				dst.put(k, nv);
			}
		}
	}

	private static void addInto(Map<Integer, BigFraction> dst, Map<Integer, BigFraction> src) {
		addInto(dst, src, ONE);
	}

	/** Canonical string of a form (sorted), for exact equality keys. */
	private static String canon(Map<Integer, BigFraction> form) {
		TreeMap<Integer, BigFraction> t = new TreeMap<>(form);
		StringBuilder sb = new StringBuilder();
		for (Map.Entry<Integer, BigFraction> e : t.entrySet()) {
			BigFraction v = e.getValue();
			sb.append(e.getKey()).append(':').append(v.getNumerator()).append('/')
					.append(v.getDenominator()).append(';');
		}
		return sb.toString();
	}

	/** Unit form (divided by value at smallest key) + that scale. */
	private static Object[] unit(Map<Integer, BigFraction> form) {
		int minK = Integer.MAX_VALUE;
		for (int k : form.keySet()) {
			if (k < minK) {
				minK = k;
			}
		}
		BigFraction sc = form.get(minK);
		TreeMap<Integer, BigFraction> u = new TreeMap<>();
		for (Map.Entry<Integer, BigFraction> e : form.entrySet()) {
			u.put(e.getKey(), e.getValue().divide(sc));
		}
		return new Object[] { u, sc };
	}

	// ===================================================================
	// rank-1 row -> rebalanced double factors
	// ===================================================================
	private static final class Row {
		final Map<Integer, BigFraction> a, b, x; // x keyed by (row*N+col) in X coords

		Row(Map<Integer, BigFraction> a, Map<Integer, BigFraction> b, Map<Integer, BigFraction> x) {
			this.a = a;
			this.b = b;
			this.x = x;
		}
	}

	private static BigInteger gcd(BigInteger a, BigInteger b) {
		return a.gcd(b);
	}

	private static BigInteger lcm(BigInteger a, BigInteger b) {
		if (a.signum() == 0 || b.signum() == 0) {
			return BigInteger.ZERO;
		}
		return a.divide(a.gcd(b)).multiply(b);
	}

	/** Factor out the rational content; returns {primitive-integer map, content}. */
	private static Object[] content(Map<Integer, BigFraction> form) {
		BigInteger commonDen = BigInteger.ONE;
		for (BigFraction v : form.values()) {
			commonDen = lcm(commonDen, v.getDenominator());
		}
		Map<Integer, BigInteger> ints = new HashMap<>();
		BigInteger g = BigInteger.ZERO;
		for (Map.Entry<Integer, BigFraction> e : form.entrySet()) {
			BigFraction scaled = e.getValue().multiply(commonDen);
			BigInteger iv = scaled.getNumerator(); // denominator is 1 after *commonDen
			ints.put(e.getKey(), iv);
			g = gcd(g, iv.abs());
		}
		if (g.signum() == 0) {
			g = BigInteger.ONE;
		}
		Map<Integer, BigInteger> prim = new HashMap<>();
		for (Map.Entry<Integer, BigInteger> e : ints.entrySet()) {
			prim.put(e.getKey(), e.getValue().divide(g));
		}
		BigFraction c = new BigFraction(g, commonDen);
		return new Object[] { prim, c };
	}

	/** Split den (smooth, small) into three near-even factors. */
	private static BigInteger[] splitDen(BigInteger den) {
		long d = den.longValueExact();
		long[] parts = { 1, 1, 1 };
		for (long pr = 2; pr * pr <= d; pr++) {
			while (d % pr == 0) {
				int mi = (parts[0] <= parts[1] && parts[0] <= parts[2]) ? 0 : (parts[1] <= parts[2] ? 1 : 2);
				parts[mi] *= pr;
				d /= pr;
			}
		}
		if (d > 1) {
			int mi = (parts[0] <= parts[1] && parts[0] <= parts[2]) ? 0 : (parts[1] <= parts[2] ? 1 : 2);
			parts[mi] *= d;
		}
		return new BigInteger[] { BigInteger.valueOf(parts[0]), BigInteger.valueOf(parts[1]),
				BigInteger.valueOf(parts[2]) };
	}

	private static NonCubicBilinearAlgorithm emit(int n, List<Row> rows) {
		int dim = n * n;
		int r = rows.size();
		double[][] U = new double[dim][r];
		double[][] V = new double[dim][r];
		double[][] W = new double[dim][r];
		for (int k = 0; k < r; k++) {
			Row row = rows.get(k);
			Object[] ca = content(row.a);
			Object[] cb = content(row.b);
			Object[] cx = content(row.x);
			@SuppressWarnings("unchecked")
			Map<Integer, BigInteger> pa = (Map<Integer, BigInteger>) ca[0];
			@SuppressWarnings("unchecked")
			Map<Integer, BigInteger> pb = (Map<Integer, BigInteger>) cb[0];
			@SuppressWarnings("unchecked")
			Map<Integer, BigInteger> px = (Map<Integer, BigInteger>) cx[0];
			BigFraction s = ((BigFraction) ca[1]).multiply((BigFraction) cb[1]).multiply((BigFraction) cx[1]);
			BigInteger numS = s.getNumerator();
			BigInteger denS = s.getDenominator();
			BigInteger[] dd = splitDen(denS);
			for (Map.Entry<Integer, BigInteger> e : pa.entrySet()) {
				U[e.getKey()][k] = new BigFraction(e.getValue(), dd[0]).doubleValue();
			}
			for (Map.Entry<Integer, BigInteger> e : pb.entrySet()) {
				V[e.getKey()][k] = new BigFraction(e.getValue(), dd[1]).doubleValue();
			}
			for (Map.Entry<Integer, BigInteger> e : px.entrySet()) {
				int kx = e.getKey();
				int a = kx / n, b = kx % n;
				int wIdx = b * n + a; // transpose-then-row-major
				W[wIdx][k] = new BigFraction(e.getValue().multiply(numS), dd[2]).doubleValue();
			}
		}
		return new NonCubicBilinearAlgorithm(n, n, n, U, V, W);
	}

	// ===================================================================
	// ODD construction
	// ===================================================================
	private static final class Odd {
		final int N; // matmul dimension (= 2m-1)
		final int m;
		final int d;

		Odd(int n) {
			this.N = n;
			this.m = (n + 1) / 2;
			this.d = m + 1;
		}

		private Map<Integer, BigFraction> single(int row, int col) {
			Map<Integer, BigFraction> r = new HashMap<>();
			r.put(row * N + col, ONE);
			return r;
		}

		private int[] blockDims(int block) {
			switch (block) {
			case 1: return new int[] { 0, 0, m, m };
			case 2: return new int[] { 0, m, m, m - 1 };
			case 3: return new int[] { m, 0, m - 1, m };
			case 4: return new int[] { m, m, m - 1, m - 1 };
			default: throw new IllegalArgumentException("block " + block);
			}
		}

		private Map<Integer, BigFraction> raw(int block, int r, int c) {
			int[] dm = blockDims(block);
			int ro = dm[0], co = dm[1], rows = dm[2], cols = dm[3];
			if (r < 0 || r >= rows || c < 0 || c >= cols) {
				return new HashMap<>();
			}
			return single(ro + r, co + c);
		}

		private Map<Integer, BigFraction> aEntry(int block, int r, int c) {
			int rr = r, cc = c;
			if (block == 3 || block == 4) {
				if (rr == m - 1) return new HashMap<>();
				if (rr == m) rr = m - 1;
			}
			if (block == 2 || block == 4) {
				if (cc == m - 1) return new HashMap<>();
				if (cc == m) cc = m - 1;
			}
			int[] dm = blockDims(block);
			int rows = dm[2], cols = dm[3];
			if (rr < 0 || rr > rows || cc < 0 || cc > cols) return new HashMap<>();
			Map<Integer, BigFraction> rowSum = new HashMap<>();
			for (int t = 0; t < cols; t++) addInto(rowSum, raw(block, rr, t));
			Map<Integer, BigFraction> colSum = new HashMap<>();
			for (int t = 0; t < rows; t++) addInto(colSum, raw(block, t, cc));
			Map<Integer, BigFraction> total = new HashMap<>();
			for (int t = 0; t < rows; t++) for (int s = 0; s < cols; s++) addInto(total, raw(block, t, s));
			BigFraction inv = f(1).divide(f(cols + 1));
			Map<Integer, BigFraction> res = new HashMap<>();
			if (rr < rows && cc < cols) {
				addInto(res, raw(block, rr, cc));
				addInto(res, rowSum, inv.negate());
			} else if (rr < rows && cc == cols) {
				addInto(res, rowSum, inv.negate());
			} else if (rr == rows && cc < cols) {
				addInto(res, colSum, f(-1));
				addInto(res, total, inv);
			} else if (rr == rows && cc == cols) {
				addInto(res, total, inv);
			}
			return res;
		}

		private Map<Integer, BigFraction> bEntry(int block, int r, int c) {
			int rr = r, cc = c;
			if (block == 3 || block == 4) {
				if (rr == m - 1) return new HashMap<>();
				if (rr == m) rr = m - 1;
			}
			if (block == 2 || block == 4) {
				if (cc == m - 1) return new HashMap<>();
				if (cc == m) cc = m - 1;
			}
			int[] dm = blockDims(block);
			int rows = dm[2], cols = dm[3];
			if (rr < 0 || rr > rows || cc < 0 || cc > cols) return new HashMap<>();
			Map<Integer, BigFraction> rowSum = new HashMap<>();
			for (int t = 0; t < cols; t++) addInto(rowSum, raw(block, rr, t));
			Map<Integer, BigFraction> colSum = new HashMap<>();
			for (int t = 0; t < rows; t++) addInto(colSum, raw(block, t, cc));
			Map<Integer, BigFraction> total = new HashMap<>();
			for (int t = 0; t < rows; t++) for (int s = 0; s < cols; s++) addInto(total, raw(block, t, s));
			Map<Integer, BigFraction> res = new HashMap<>();
			if (rr < rows && cc < cols) {
				addInto(res, raw(block, rr, cc));
			} else if (rr < rows && cc == cols) {
				addInto(res, rowSum, f(-1));
			} else if (rr == rows && cc < cols) {
				addInto(res, colSum, f(-1));
			} else if (rr == rows && cc == cols) {
				addInto(res, total);
			}
			return res;
		}

		private Map<Integer, BigFraction> cEntry(int block, int r, int c) {
			int[] dm = blockDims(block);
			int ro = dm[0], co = dm[1], rows = dm[2], cols = dm[3];
			if (c < 0 || c >= rows || r < 0 || r >= cols) return new HashMap<>();
			int row = ro + c;
			int col = co + r;
			return single(col, row); // keyed (col,row) like the Maple transposed pullback
		}

		private Map<Integer, BigFraction> af(int[][] codes, int r, int c) {
			Map<Integer, BigFraction> res = new HashMap<>();
			for (int[] e : codes) addInto(res, aEntry(e[1], r, c), f(e[0]));
			return res;
		}

		private Map<Integer, BigFraction> bf(int[][] codes, int r, int c) {
			Map<Integer, BigFraction> res = new HashMap<>();
			for (int[] e : codes) addInto(res, bEntry(e[1], r, c), f(e[0]));
			return res;
		}

		private Map<Integer, BigFraction> cf(int[][] codes, int r, int c) {
			Map<Integer, BigFraction> res = new HashMap<>();
			for (int[] e : codes) addInto(res, cEntry(e[1], r, c), f(e[0]));
			return res;
		}

		private Map<Integer, BigFraction> negInner(Map<Integer, BigFraction> inner) {
			Map<Integer, BigFraction> r = new HashMap<>();
			addInto(r, inner, f(-1));
			return r;
		}

		// label = {kind, call, i, j, k}; kinds: 0 s0,1 s1,2 s2,3 u2p,4 u2,5 u1,6 u3,7 u4
		private Map<Integer, BigFraction> labelFactor(int[] L, int axis) {
			int kind = L[0];
			if (kind == 0) { // s0
				int call = L[1], i = L[2];
				int[][] ac, bc, cc;
				if (call == 0) {
					ac = new int[][] { { 1, 2 }, { -1, 1 }, { 1, 3 } };
					bc = new int[][] { { 1, 3 }, { 1, 2 }, { 1, 1 } };
					cc = new int[][] { { 1, 1 }, { -1, 2 }, { 1, 3 } };
				} else {
					ac = new int[][] { { 1, 2 }, { 1, 3 }, { -1, 4 } };
					bc = new int[][] { { 1, 4 }, { 1, 2 }, { 1, 3 } };
					cc = new int[][] { { 1, 2 }, { 1, 4 }, { -1, 3 } };
				}
				if (axis == 1) return af(ac, i, i);
				if (axis == 2) return bf(bc, i, i);
				return cf(cc, i, i);
			}
			if (kind == 1) { // s1
				int call = L[1], i = L[2], j = L[3], k = L[4];
				int[][] ac, bc, cc;
				if (call == 0) { ac = new int[][] { { 1, 1 } }; bc = new int[][] { { 1, 1 } }; cc = new int[][] { { 1, 1 } }; }
				else { ac = new int[][] { { 1, 4 } }; bc = new int[][] { { 1, 4 } }; cc = new int[][] { { 1, 4 } }; }
				if (axis == 1) {
					Map<Integer, BigFraction> r = af(ac, i, j); addInto(r, af(ac, j, k)); addInto(r, af(ac, k, i)); return r;
				}
				if (axis == 2) {
					Map<Integer, BigFraction> r = bf(bc, j, k); addInto(r, bf(bc, k, i)); addInto(r, bf(bc, i, j)); return r;
				}
				Map<Integer, BigFraction> r = cf(cc, k, i); addInto(r, cf(cc, i, j)); addInto(r, cf(cc, j, k)); return r;
			}
			if (kind == 2) { // s2
				int call = L[1], i = L[2], j = L[3], k = L[4];
				int[][] ac, bc, cc, uc, vc, wc, xc, yc, zc;
				if (call == 0) {
					ac = new int[][] { { -1, 1 } }; bc = new int[][] { { 1, 2 } }; cc = new int[][] { { -1, 2 } };
					uc = new int[][] { { 1, 3 } }; vc = new int[][] { { 1, 1 } }; wc = new int[][] { { 1, 3 } };
					xc = new int[][] { { 1, 2 } }; yc = new int[][] { { 1, 3 } }; zc = new int[][] { { 1, 1 } };
				} else {
					ac = new int[][] { { -1, 4 } }; bc = new int[][] { { 1, 3 } }; cc = new int[][] { { -1, 3 } };
					uc = new int[][] { { 1, 2 } }; vc = new int[][] { { 1, 4 } }; wc = new int[][] { { 1, 2 } };
					xc = new int[][] { { 1, 3 } }; yc = new int[][] { { 1, 2 } }; zc = new int[][] { { 1, 4 } };
				}
				if (axis == 1) {
					Map<Integer, BigFraction> r = af(ac, i, j); addInto(r, af(uc, j, k)); addInto(r, af(xc, k, i)); return r;
				}
				if (axis == 2) {
					Map<Integer, BigFraction> r = bf(bc, j, k); addInto(r, bf(vc, k, i)); addInto(r, bf(yc, i, j)); return r;
				}
				Map<Integer, BigFraction> r = cf(cc, k, i); addInto(r, cf(wc, i, j)); addInto(r, cf(zc, j, k)); return r;
			}
			if (kind == 3) { // u2p
				int call = L[1], i = L[2];
				int[][] ac = (call == 0) ? new int[][] { { 1, 1 } } : new int[][] { { 1, 4 } };
				if (axis == 1) return af(ac, i, i);
				if (axis == 2) return bf(ac, i, i);
				Map<Integer, BigFraction> inner = new HashMap<>();
				addInto(inner, cf(ac, i, i), f(d - 9));
				for (int k = 0; k < d; k++) { addInto(inner, cf(ac, k, i)); addInto(inner, cf(ac, i, k)); }
				return negInner(inner);
			}
			if (kind == 4) { // u2
				int call = L[1], i = L[2];
				int[][] ac, bc, dc, wc, zc;
				int[][][] T = {
						{ { -1, 1 }, { 1, 3 }, { -1, 2 }, { 1, 3 }, { 1, 1 } },
						{ { 1, 2 }, { 1, 1 }, { 1, 1 }, { -1, 2 }, { 1, 3 } },
						{ { 1, 3 }, { 1, 2 }, { 1, 3 }, { 1, 1 }, { -1, 2 } },
						{ { -1, 4 }, { 1, 2 }, { -1, 3 }, { 1, 2 }, { 1, 4 } },
						{ { 1, 3 }, { 1, 4 }, { 1, 4 }, { -1, 3 }, { 1, 2 } },
						{ { 1, 2 }, { 1, 3 }, { 1, 2 }, { 1, 4 }, { -1, 3 } } };
				int[][] row = T[call];
				ac = new int[][] { row[0] }; bc = new int[][] { row[1] }; dc = new int[][] { row[2] };
				wc = new int[][] { row[3] }; zc = new int[][] { row[4] };
				if (axis == 1) return af(ac, i, i);
				if (axis == 2) return bf(bc, i, i);
				Map<Integer, BigFraction> inner = new HashMap<>();
				for (int k = 0; k < d; k++) addInto(inner, cf(dc, k, i));
				addInto(inner, cf(wc, i, i), f(d));
				for (int k = 0; k < d; k++) addInto(inner, cf(zc, i, k));
				return negInner(inner);
			}
			if (kind == 5) { // u1
				int call = L[1], i = L[2], j = L[3];
				int[][][] T = {
						{ { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 } },
						{ { -1, 1 }, { 1, 3 }, { -1, 2 }, { 1, 3 }, { 1, 1 } },
						{ { 1, 2 }, { 1, 3 }, { 1, 2 }, { 1, 4 }, { -1, 3 } },
						{ { -1, 4 }, { 1, 2 }, { -1, 3 }, { 1, 2 }, { 1, 4 } },
						{ { 1, 2 }, { 1, 1 }, { 1, 1 }, { -1, 2 }, { 1, 3 } },
						{ { 1, 3 }, { 1, 2 }, { 1, 3 }, { 1, 1 }, { -1, 2 } },
						{ { 1, 3 }, { 1, 4 }, { 1, 4 }, { -1, 3 }, { 1, 2 } },
						{ { 1, 4 }, { 1, 4 }, { 1, 4 }, { 1, 4 }, { 1, 4 } } };
				int[][] row = T[call];
				int[][] ac = { row[0] }, bc = { row[1] }, cc = { row[2] }, wc = { row[3] }, zc = { row[4] };
				if (axis == 1) return af(ac, i, j);
				if (axis == 2) return bf(bc, i, j);
				Map<Integer, BigFraction> inner = new HashMap<>();
				addInto(inner, cf(wc, i, j), f(d));
				for (int k = 0; k < d; k++) { addInto(inner, cf(cc, k, i)); addInto(inner, cf(zc, j, k)); }
				return negInner(inner);
			}
			if (kind == 6) { // u3
				int call = L[1], i = L[2];
				int[][] ac, yc, cc;
				switch (call) {
				case 0: ac = new int[][] { { 1, 1 }, { 1, 2 } }; yc = new int[][] { { 1, 1 } }; cc = new int[][] { { 1, 1 } }; break;
				case 1: ac = new int[][] { { 1, 1 }, { 1, 2 } }; yc = new int[][] { { 1, 3 } }; cc = new int[][] { { 1, 2 } }; break;
				case 2: ac = new int[][] { { 1, 3 }, { 1, 4 } }; yc = new int[][] { { 1, 2 } }; cc = new int[][] { { 1, 3 } }; break;
				default: ac = new int[][] { { 1, 3 }, { 1, 4 } }; yc = new int[][] { { 1, 4 } }; cc = new int[][] { { 1, 4 } }; break;
				}
				if (axis == 1) return af(ac, i, d - 1);
				if (axis == 2) return bf(yc, i, d - 1);
				Map<Integer, BigFraction> inner = new HashMap<>();
				for (int k = 0; k < d; k++) addInto(inner, cf(cc, k, i));
				return negInner(inner);
			}
			if (kind == 7) { // u4
				int call = L[1], j = L[2];
				int[][] ac, yc, zc;
				switch (call) {
				case 0: ac = new int[][] { { 1, 1 } }; yc = new int[][] { { 1, 1 }, { -1, 3 } }; zc = new int[][] { { 1, 1 } }; break;
				case 1: ac = new int[][] { { 1, 2 } }; yc = new int[][] { { 1, 1 }, { -1, 3 } }; zc = new int[][] { { 1, 3 } }; break;
				case 2: ac = new int[][] { { 1, 3 } }; yc = new int[][] { { -1, 2 }, { 1, 4 } }; zc = new int[][] { { 1, 2 } }; break;
				default: ac = new int[][] { { 1, 4 } }; yc = new int[][] { { -1, 2 }, { 1, 4 } }; zc = new int[][] { { 1, 4 } }; break;
				}
				if (axis == 1) return af(ac, d - 1, j);
				if (axis == 2) return bf(yc, d - 1, j);
				Map<Integer, BigFraction> inner = new HashMap<>();
				for (int k = 0; k < d; k++) addInto(inner, cf(zc, j, k));
				return negInner(inner);
			}
			throw new IllegalStateException("kind " + kind);
		}

		private boolean propsize(int i, int j, int k) {
			int c = 0;
			if (i == d - 1) c++;
			if (j == d - 1) c++;
			if (k == d - 1) c++;
			return c < 2;
		}

		private List<int[]> groupFromTarget(int[] L, int z) {
			List<int[]> g = new ArrayList<>();
			g.add(L);
			if (L[0] == 1 && L[1] == 1) {
				if (L[2] == L[3] && L[4] == z) {
					g.add(new int[] { 3, 1, L[2], -1, -1 });
				} else if (L[4] == z) {
					g.add(new int[] { 5, 7, L[2], L[3], -1 });
				} else if (L[2] == z) {
					g.add(new int[] { 5, 7, L[3], L[4], -1 });
				}
			}
			return g;
		}

		private List<int[]> cycleTargets(int first, int second, int z) {
			List<int[]> t = new ArrayList<>();
			t.add(new int[] { 1, 1, first, first, z });
			t.add(new int[] { 4, 4, first, -1, -1 });
			t.add(new int[] { 4, 2, first, -1, -1 });
			t.add(new int[] { 2, 1, first, first, z });
			t.add(new int[] { 4, 3, first, -1, -1 });
			t.add(new int[] { 2, 1, first, z, first });
			t.add(new int[] { 1, 1, first, second, z });
			t.add(new int[] { 2, 1, second, z, first });
			t.add(new int[] { 5, 6, first, second, -1 });
			t.add(new int[] { 5, 5, first, second, -1 });
			t.add(new int[] { 2, 1, first, second, z });
			t.add(new int[] { 5, 3, first, second, -1 });
			t.add(new int[] { 1, 1, second, second, z });
			t.add(new int[] { 2, 1, second, z, second });
			t.add(new int[] { 4, 4, second, -1, -1 });
			t.add(new int[] { 4, 2, second, -1, -1 });
			t.add(new int[] { 4, 3, second, -1, -1 });
			t.add(new int[] { 2, 1, second, second, z });
			t.add(new int[] { 1, 1, z, second, first });
			t.add(new int[] { 5, 3, second, first, -1 });
			t.add(new int[] { 5, 5, second, first, -1 });
			t.add(new int[] { 5, 6, second, first, -1 });
			t.add(new int[] { 2, 1, second, first, z });
			t.add(new int[] { 2, 1, first, z, second });
			return t;
		}

		private List<List<int[]>> orderedCycleGroups(int first, int second, int z) {
			List<List<int[]>> g = new ArrayList<>();
			for (int[] tgt : cycleTargets(first, second, z)) {
				g.add(groupFromTarget(tgt, z));
			}
			return g;
		}

		private boolean labelEq(int[] a, int[] b) {
			for (int i = 0; i < 5; i++) if (a[i] != b[i]) return false;
			return true;
		}

		private boolean groupEq(List<int[]> a, List<int[]> b) {
			if (a.size() != b.size()) return false;
			for (int i = 0; i < a.size(); i++) if (!labelEq(a.get(i), b.get(i))) return false;
			return true;
		}

		private void appendGroup(List<List<int[]>> out, List<int[]> g) {
			for (List<int[]> e : out) if (groupEq(e, g)) return;
			out.add(g);
		}

		private void appendRange(List<List<int[]>> out, List<List<int[]>> groups, int a, int b) {
			for (int t = a; t <= b; t++) appendGroup(out, groups.get(t - 1));
		}

		private List<List<int[]>> orderedPath3Groups(int[] idx, int z) {
			List<List<int[]>> left = orderedCycleGroups(idx[0], idx[1], z);
			List<List<int[]>> right = orderedCycleGroups(idx[1], idx[2], z);
			List<List<int[]>> out = new ArrayList<>();
			appendRange(out, left, 1, 6);
			appendRange(out, left, 7, 12);
			appendRange(out, left, 13, 18);
			appendRange(out, left, 19, 24);
			appendRange(out, right, 7, 12);
			appendRange(out, right, 13, 18);
			appendRange(out, right, 19, 24);
			if (out.size() != 42) throw new IllegalStateException("path3 " + out.size());
			return out;
		}

		private List<List<int[]>> orderedPath4Groups(int[] idx, int z) {
			List<List<int[]>> first = orderedCycleGroups(idx[0], idx[1], z);
			List<List<int[]>> out = new ArrayList<>();
			appendRange(out, first, 1, 6);
			for (int t = 0; t < idx.length - 1; t++) {
				List<List<int[]>> cyc = orderedCycleGroups(idx[t], idx[t + 1], z);
				appendRange(out, cyc, 7, 12);
				appendRange(out, cyc, 13, 18);
				appendRange(out, cyc, 19, 24);
			}
			if (out.size() != 60) throw new IllegalStateException("path4 " + out.size());
			return out;
		}

		private List<int[]> indexBlocks() {
			int regular = m - 1;
			int q = regular / 4, rem = regular % 4;
			List<Integer> lengths = new ArrayList<>();
			if (rem == 0) {
				for (int t = 0; t < q; t++) lengths.add(4);
			} else if (rem == 1) {
				for (int t = 0; t < q - 1; t++) lengths.add(4);
				lengths.add(3); lengths.add(2);
			} else if (rem == 2) {
				for (int t = 0; t < q; t++) lengths.add(4);
				lengths.add(2);
			} else {
				for (int t = 0; t < q; t++) lengths.add(4);
				lengths.add(3);
			}
			List<int[]> blocks = new ArrayList<>();
			int start = 0;
			for (int ell : lengths) {
				int[] blk = new int[ell];
				for (int t = 0; t < ell; t++) blk[t] = start + t;
				blocks.add(blk);
				start += ell;
			}
			return blocks;
		}

		private Map<Integer, BigFraction> groupFactor(List<int[]> group, int axis) {
			if (axis == 1 || axis == 2) return labelFactor(group.get(0), axis);
			Map<Integer, BigFraction> res = new HashMap<>();
			for (int[] L : group) addInto(res, labelFactor(L, 3));
			return res;
		}

		private Map<Integer, BigFraction> localForm(List<List<int[]>> groups, List<int[]> tabG, List<BigFraction> tabC, int term, int axis) {
			Map<Integer, BigFraction> res = new HashMap<>();
			for (int i = 0; i < tabG.size(); i++) {
				int[] gt = tabG.get(i);
				if (gt[1] == term) {
					addInto(res, groupFactor(groups.get(gt[0] - 1), axis), tabC.get(i));
				}
			}
			return res;
		}

		@SuppressWarnings("unchecked")
		NonCubicBilinearAlgorithm build() {
			List<BigFraction> coeffs = new ArrayList<>();
			List<Map<Integer, BigFraction>> As = new ArrayList<>();
			List<Map<Integer, BigFraction>> Bs = new ArrayList<>();
			List<Map<Integer, BigFraction>> Xs = new ArrayList<>();

			BigFraction D = new BigFraction(d);

			// ---- BASE ----
			for (int i = 0; i < m; i++) {
				addTerm(coeffs, As, Bs, Xs, ONE, new int[] { 0, 0, i, -1, -1 });
				addTerm(coeffs, As, Bs, Xs, ONE, new int[] { 0, 1, i, -1, -1 });
			}
			for (int call = 0; call <= 1; call++)
				for (int i = 0; i < d; i++)
					for (int j = 0; j < d; j++)
						for (int k = 0; k < d; k++)
							if (((i <= j && j < k) || (k < j && j <= i)) && propsize(i, j, k))
								addTerm(coeffs, As, Bs, Xs, ONE, new int[] { 1, call, i, j, k });
			for (int call = 0; call <= 1; call++)
				for (int i = 0; i < d; i++)
					for (int j = 0; j < d; j++)
						for (int k = 0; k < d; k++)
							if (!(i == j && j == k) && propsize(i, j, k))
								addTerm(coeffs, As, Bs, Xs, ONE, new int[] { 2, call, i, j, k });
			for (int i = 0; i < m; i++) {
				addTerm(coeffs, As, Bs, Xs, ONE, new int[] { 3, 0, i, -1, -1 });
				addTerm(coeffs, As, Bs, Xs, ONE, new int[] { 3, 1, i, -1, -1 });
			}
			for (int q = 0; q <= 5; q++)
				for (int i = 0; i < m; i++)
					addTerm(coeffs, As, Bs, Xs, ONE, new int[] { 4, q, i, -1, -1 });
			for (int q = 0; q <= 7; q++)
				for (int i = 0; i < m; i++)
					for (int j = 0; j < m; j++)
						if (i != j)
							addTerm(coeffs, As, Bs, Xs, ONE, new int[] { 5, q, i, j, -1 });
			for (int q = 0; q <= 3; q++) {
				for (int i = 0; i < m; i++)
					addTerm(coeffs, As, Bs, Xs, ONE, new int[] { 6, q, i, -1, -1 });
				for (int j = 0; j < m; j++)
					addTerm(coeffs, As, Bs, Xs, ONE, new int[] { 7, q, j, -1, -1 });
			}

			// ---- LOCAL REPLACEMENTS ----
			int z = m - 1;
			for (int[] block : indexBlocks()) {
				List<List<int[]>> groups;
				List<int[]> ug, vg, wg;
				List<BigFraction> uc, vc, wc;
				int rk;
				if (block.length == 4) {
					groups = orderedPath4Groups(block, z);
					ug = path4UG(); uc = path4UC(D);
					vg = path4VG(); vc = path4VC(D);
					wg = path4WG(); wc = path4WC(D);
					rk = 57;
				} else if (block.length == 3) {
					groups = orderedPath3Groups(block, z);
					ug = path3UG(); uc = path3UC(D);
					vg = path3VG(); vc = path3VC(D);
					wg = path3WG(); wc = path3WC(D);
					rk = 40;
				} else if (block.length == 2) {
					groups = orderedCycleGroups(block[0], block[1], z);
					ug = cycleUG(); uc = cycleUC(D);
					vg = cycleVG(); vc = cycleVC(D);
					wg = cycleWG(); wc = cycleWC(D);
					rk = 23;
				} else {
					continue;
				}
				for (List<int[]> group : groups)
					for (int[] L : group)
						addTerm(coeffs, As, Bs, Xs, f(-1), L);
				for (int r = 1; r <= rk; r++) {
					Map<Integer, BigFraction> fu = localForm(groups, ug, uc, r, 1);
					Map<Integer, BigFraction> fv = localForm(groups, vg, vc, r, 2);
					Map<Integer, BigFraction> fw = localForm(groups, wg, wc, r, 3);
					coeffs.add(ONE); As.add(fu); Bs.add(fv); Xs.add(fw);
				}
			}

			// ---- exact UV kin unification (up to scalar) + zero deletion ----
			Map<String, Object[]> acc = new HashMap<>();
			for (int t = 0; t < coeffs.size(); t++) {
				Map<Integer, BigFraction> A = As.get(t), B = Bs.get(t), X = Xs.get(t);
				if (A.isEmpty() || B.isEmpty()) continue;
				Object[] ua = unit(A);
				Object[] ub = unit(B);
				@SuppressWarnings("unchecked")
				TreeMap<Integer, BigFraction> au = (TreeMap<Integer, BigFraction>) ua[0];
				@SuppressWarnings("unchecked")
				TreeMap<Integer, BigFraction> bu = (TreeMap<Integer, BigFraction>) ub[0];
				BigFraction sa = (BigFraction) ua[1], sb = (BigFraction) ub[1];
				String key = canon(au) + "|" + canon(bu);
				Object[] slot = acc.get(key);
				if (slot == null) {
					slot = new Object[] { au, bu, new HashMap<Integer, BigFraction>() };
					acc.put(key, slot);
				}
				BigFraction scale = coeffs.get(t).multiply(sa).multiply(sb);
				addInto((Map<Integer, BigFraction>) slot[2], X, scale);
			}
			List<Row> rows = new ArrayList<>();
			for (Object[] slot : acc.values()) {
				@SuppressWarnings("unchecked")
				Map<Integer, BigFraction> W = (Map<Integer, BigFraction>) slot[2];
				if (W.isEmpty()) continue;
				@SuppressWarnings("unchecked")
				TreeMap<Integer, BigFraction> au = (TreeMap<Integer, BigFraction>) slot[0];
				@SuppressWarnings("unchecked")
				TreeMap<Integer, BigFraction> bu = (TreeMap<Integer, BigFraction>) slot[1];
				rows.add(new Row(new HashMap<>(au), new HashMap<>(bu), W));
			}
			return emit(N, rows);
		}

		private void addTerm(List<BigFraction> coeffs, List<Map<Integer, BigFraction>> As,
				List<Map<Integer, BigFraction>> Bs, List<Map<Integer, BigFraction>> Xs,
				BigFraction coeff, int[] L) {
			coeffs.add(coeff);
			As.add(labelFactor(L, 1));
			Bs.add(labelFactor(L, 2));
			Xs.add(labelFactor(L, 3));
		}

		private static List<int[]> path4UG() { return java.util.Arrays.asList(
			new int[]{1,5},new int[]{1,9},new int[]{1,10},new int[]{1,51},new int[]{2,1},new int[]{2,2},new int[]{2,3},new int[]{2,4},new int[]{2,5},new int[]{2,6},new int[]{2,8},new int[]{2,49},new int[]{2,51},new int[]{6,1},new int[]{6,2},new int[]{6,3},new int[]{6,4},new int[]{6,5},new int[]{6,6},new int[]{6,49},new int[]{7,3},new int[]{7,33},new int[]{7,34},new int[]{7,49},new int[]{8,1},new int[]{8,2},new int[]{8,3},new int[]{8,49},new int[]{13,18},new int[]{13,34},new int[]{13,42},new int[]{14,7},new int[]{14,11},new int[]{14,12},new int[]{14,15},new int[]{14,16},new int[]{14,17},new int[]{14,18},new int[]{14,19},new int[]{14,52},new int[]{15,7},new int[]{15,11},new int[]{15,12},new int[]{15,15},new int[]{15,16},new int[]{15,17},new int[]{15,18},new int[]{15,19},new int[]{15,43},new int[]{15,52},new int[]{19,9},new int[]{19,12},new int[]{19,14},new int[]{19,50},new int[]{21,11},new int[]{21,12},new int[]{21,13},new int[]{21,50},new int[]{25,16},new int[]{25,20},new int[]{25,21},new int[]{25,52},new int[]{26,7},new int[]{26,15},new int[]{26,16},new int[]{26,52},new int[]{31,20},new int[]{31,23},new int[]{31,25},new int[]{31,54},new int[]{32,22},new int[]{32,23},new int[]{32,35},new int[]{32,36},new int[]{32,39},new int[]{32,40},new int[]{32,41},new int[]{32,46},new int[]{32,57},new int[]{33,22},new int[]{33,23},new int[]{33,24},new int[]{33,35},new int[]{33,36},new int[]{33,39},new int[]{33,40},new int[]{33,41},new int[]{33,46},new int[]{33,54},new int[]{33,57},new int[]{37,34},new int[]{37,36},new int[]{37,38},new int[]{37,53},new int[]{39,35},new int[]{39,36},new int[]{39,37},new int[]{39,53},new int[]{43,40},new int[]{43,44},new int[]{43,45},new int[]{43,57},new int[]{44,39},new int[]{44,40},new int[]{44,46},new int[]{44,57},new int[]{49,31},new int[]{49,44},new int[]{49,47},new int[]{49,56},new int[]{50,26},new int[]{50,27},new int[]{50,30},new int[]{50,31},new int[]{50,32},new int[]{51,26},new int[]{51,27},new int[]{51,30},new int[]{51,31},new int[]{51,32},new int[]{51,48},new int[]{51,56},new int[]{55,20},new int[]{55,27},new int[]{55,29},new int[]{55,55},new int[]{57,26},new int[]{57,27},new int[]{57,28},new int[]{57,55}); }
		private static List<BigFraction> path4UC(BigFraction D) { return java.util.Arrays.asList(
			f(1).negate(), f(1), f(1).negate(), f(1).negate(), f(1).negate(), f(1).negate(), f(1).negate(), D.negate(), D.negate(), f(1).negate(), f(1), f(1).negate(), f(1), f(1), f(1), f(1), f(1), f(1), f(1), f(1), f(1).negate(), f(1).negate(), f(1).negate(), f(1).divide(f(2)), f(1).negate(), D.subtract(f(1)).negate().divide(D), D.subtract(f(1)).negate().divide(D), f(1).negate(), f(1), f(1).negate(), f(1), f(1), f(1).divide(D), f(1).divide(D), f(1), f(1), f(1), f(1), f(1), f(1), f(1).negate(), f(1).negate().divide(D), f(1).negate().divide(D), f(1).negate(), f(1).negate(), D.negate(), D.negate(), f(1).negate(), f(1), f(1).negate(), f(1), f(1), f(1), f(2).negate(), D.subtract(f(1)).negate().divide(D), D.subtract(f(1)).negate().divide(D), f(1).negate(), f(1).negate(), f(1), f(1), f(1), f(2).negate(), f(1).negate(), D.subtract(f(1)).negate().divide(D), D.subtract(f(1)).negate().divide(D), f(1).negate(), f(1), f(1), f(1), f(2).negate(), f(1).negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), D, D, f(1), f(1).negate(), f(1).negate(), D.divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), f(1), f(1).divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), D.negate(), D.negate(), f(1).negate(), f(1), f(1), f(1), f(1).negate(), f(1), f(1), f(2).negate(), f(1), f(1), f(1), f(1), f(1), f(1), f(1), f(1), f(1).subtract(D), f(1).subtract(D), f(1), f(1), f(1), f(1), f(1), f(1), f(1).divide(D), f(1).divide(D), f(1), f(1), f(1), f(1).negate().divide(D), f(1).negate().divide(D), D.negate(), D.negate(), f(1).negate(), f(1), f(1).negate().divide(f(2)), f(1), f(1), f(1), f(4), D.subtract(f(1)).negate().divide(D), D.subtract(f(1)).negate().divide(D), f(1).negate(), f(1).negate()); }
		private static List<int[]> path4VG() { return java.util.Arrays.asList(
			new int[]{1,4},new int[]{1,5},new int[]{1,6},new int[]{1,51},new int[]{3,8},new int[]{3,9},new int[]{3,10},new int[]{3,14},new int[]{3,51},new int[]{4,9},new int[]{4,10},new int[]{4,14},new int[]{7,2},new int[]{7,3},new int[]{7,6},new int[]{7,49},new int[]{10,1},new int[]{10,33},new int[]{10,34},new int[]{10,38},new int[]{10,42},new int[]{10,49},new int[]{11,33},new int[]{11,34},new int[]{11,38},new int[]{11,42},new int[]{13,17},new int[]{13,18},new int[]{13,19},new int[]{16,42},new int[]{16,43},new int[]{19,11},new int[]{19,12},new int[]{19,19},new int[]{19,50},new int[]{20,13},new int[]{20,14},new int[]{20,50},new int[]{25,15},new int[]{25,16},new int[]{25,19},new int[]{25,52},new int[]{28,7},new int[]{28,20},new int[]{28,21},new int[]{28,25},new int[]{28,29},new int[]{28,52},new int[]{29,20},new int[]{29,21},new int[]{29,25},new int[]{29,29},new int[]{31,22},new int[]{31,23},new int[]{31,41},new int[]{31,54},new int[]{34,24},new int[]{34,25},new int[]{34,54},new int[]{37,35},new int[]{37,36},new int[]{37,41},new int[]{37,53},new int[]{38,37},new int[]{38,38},new int[]{38,53},new int[]{43,39},new int[]{43,40},new int[]{43,41},new int[]{43,57},new int[]{46,44},new int[]{46,45},new int[]{46,46},new int[]{46,47},new int[]{46,57},new int[]{47,44},new int[]{47,45},new int[]{47,47},new int[]{49,30},new int[]{49,31},new int[]{49,32},new int[]{49,56},new int[]{52,47},new int[]{52,48},new int[]{52,56},new int[]{55,26},new int[]{55,27},new int[]{55,32},new int[]{55,55},new int[]{56,28},new int[]{56,29},new int[]{56,55}); }
		private static List<BigFraction> path4VC(BigFraction D) { return java.util.Arrays.asList(
			f(1), f(1), f(1), f(1), f(1), f(1).negate(), D.negate(), f(1).negate().divide(D), f(1), f(1), f(1), f(1).divide(D), f(1), f(1), f(1), f(1), f(1), D.negate(), f(1).negate(), f(1).negate().divide(D), f(1).negate().divide(D), f(2), f(1), f(1), f(1).divide(D), f(1).divide(D), f(1), f(1), f(1), D.subtract(f(1)).negate().divide(D), f(1).negate(), f(1), f(1), f(1), f(2), f(1).negate(), D.subtract(f(1)).negate().divide(D), f(1).negate(), f(1), f(1), f(1), f(1), f(1), f(1).negate(), D.negate(), f(1).negate().divide(D), f(1).negate().divide(D), f(1), f(1), f(1), f(1).divide(D), f(1).divide(D), f(1), f(1), f(1), f(1), f(1).negate(), D.subtract(f(1)).negate().divide(D), f(1).negate(), f(1), f(1), f(1), f(1), f(1).negate(), D.subtract(f(1)).negate().divide(D), f(1).negate(), f(1), f(1), f(1), f(1), f(1).negate(), D.negate(), f(1), f(1).negate().divide(D), f(2), f(1), f(1), f(1).divide(D), f(1), f(1), f(1), f(1), D.subtract(f(1)).negate().divide(D), f(1).negate(), f(1).negate(), f(1), f(1), f(1), f(2), f(1).negate(), D.subtract(f(1)).negate().divide(D), f(1).negate()); }
		private static List<int[]> path4WG() { return java.util.Arrays.asList(
			new int[]{1,4},new int[]{1,5},new int[]{2,4},new int[]{2,5},new int[]{2,8},new int[]{2,9},new int[]{2,10},new int[]{2,11},new int[]{2,12},new int[]{2,13},new int[]{2,14},new int[]{2,50},new int[]{2,51},new int[]{3,8},new int[]{4,9},new int[]{4,10},new int[]{4,14},new int[]{5,9},new int[]{5,10},new int[]{5,14},new int[]{6,4},new int[]{6,5},new int[]{6,8},new int[]{6,9},new int[]{6,10},new int[]{6,11},new int[]{6,12},new int[]{6,13},new int[]{6,14},new int[]{6,50},new int[]{6,51},new int[]{7,2},new int[]{7,3},new int[]{8,2},new int[]{8,4},new int[]{8,5},new int[]{8,6},new int[]{8,8},new int[]{8,9},new int[]{8,10},new int[]{8,11},new int[]{8,12},new int[]{8,13},new int[]{8,14},new int[]{8,50},new int[]{8,51},new int[]{9,2},new int[]{9,4},new int[]{9,5},new int[]{9,6},new int[]{9,8},new int[]{9,9},new int[]{9,10},new int[]{9,11},new int[]{9,12},new int[]{9,13},new int[]{9,14},new int[]{9,50},new int[]{9,51},new int[]{10,1},new int[]{11,1},new int[]{11,2},new int[]{11,3},new int[]{11,4},new int[]{11,5},new int[]{11,6},new int[]{11,8},new int[]{11,9},new int[]{11,10},new int[]{11,14},new int[]{11,33},new int[]{11,49},new int[]{11,51},new int[]{12,1},new int[]{12,2},new int[]{12,3},new int[]{12,4},new int[]{12,5},new int[]{12,6},new int[]{12,8},new int[]{12,9},new int[]{12,10},new int[]{12,14},new int[]{12,49},new int[]{12,51},new int[]{13,17},new int[]{13,18},new int[]{14,7},new int[]{14,11},new int[]{14,15},new int[]{14,16},new int[]{14,17},new int[]{14,19},new int[]{14,20},new int[]{14,21},new int[]{14,22},new int[]{14,23},new int[]{14,24},new int[]{14,25},new int[]{14,26},new int[]{14,27},new int[]{14,28},new int[]{14,29},new int[]{14,30},new int[]{14,31},new int[]{14,32},new int[]{14,35},new int[]{14,39},new int[]{14,40},new int[]{14,41},new int[]{14,44},new int[]{14,45},new int[]{14,46},new int[]{14,47},new int[]{14,48},new int[]{14,52},new int[]{14,54},new int[]{14,55},new int[]{14,56},new int[]{14,57},new int[]{15,7},new int[]{15,11},new int[]{15,15},new int[]{15,16},new int[]{15,17},new int[]{15,19},new int[]{15,20},new int[]{15,21},new int[]{15,22},new int[]{15,23},new int[]{15,24},new int[]{15,25},new int[]{15,26},new int[]{15,27},new int[]{15,28},new int[]{15,29},new int[]{15,30},new int[]{15,31},new int[]{15,32},new int[]{15,35},new int[]{15,39},new int[]{15,40},new int[]{15,41},new int[]{15,44},new int[]{15,45},new int[]{15,46},new int[]{15,47},new int[]{15,48},new int[]{15,52},new int[]{15,54},new int[]{15,55},new int[]{15,56},new int[]{15,57},new int[]{16,43},new int[]{17,1},new int[]{17,2},new int[]{17,3},new int[]{17,4},new int[]{17,5},new int[]{17,6},new int[]{17,8},new int[]{17,9},new int[]{17,10},new int[]{17,14},new int[]{17,33},new int[]{17,34},new int[]{17,38},new int[]{17,42},new int[]{17,49},new int[]{17,51},new int[]{18,1},new int[]{18,2},new int[]{18,3},new int[]{18,4},new int[]{18,5},new int[]{18,6},new int[]{18,8},new int[]{18,9},new int[]{18,10},new int[]{18,14},new int[]{18,33},new int[]{18,34},new int[]{18,38},new int[]{18,42},new int[]{18,49},new int[]{18,51},new int[]{19,11},new int[]{19,12},new int[]{20,14},new int[]{21,13},new int[]{22,11},new int[]{25,15},new int[]{25,16},new int[]{26,7},new int[]{26,15},new int[]{26,16},new int[]{26,20},new int[]{26,21},new int[]{26,22},new int[]{26,23},new int[]{26,24},new int[]{26,25},new int[]{26,26},new int[]{26,27},new int[]{26,28},new int[]{26,29},new int[]{26,30},new int[]{26,31},new int[]{26,32},new int[]{26,35},new int[]{26,39},new int[]{26,40},new int[]{26,41},new int[]{26,44},new int[]{26,45},new int[]{26,46},new int[]{26,47},new int[]{26,48},new int[]{26,52},new int[]{26,54},new int[]{26,55},new int[]{26,56},new int[]{26,57},new int[]{27,7},new int[]{27,15},new int[]{27,16},new int[]{27,20},new int[]{27,21},new int[]{27,22},new int[]{27,23},new int[]{27,24},new int[]{27,25},new int[]{27,26},new int[]{27,27},new int[]{27,28},new int[]{27,29},new int[]{27,30},new int[]{27,31},new int[]{27,32},new int[]{27,35},new int[]{27,39},new int[]{27,40},new int[]{27,41},new int[]{27,44},new int[]{27,45},new int[]{27,46},new int[]{27,47},new int[]{27,48},new int[]{27,52},new int[]{27,54},new int[]{27,55},new int[]{27,56},new int[]{27,57},new int[]{28,7},new int[]{29,20},new int[]{29,21},new int[]{29,22},new int[]{29,23},new int[]{29,24},new int[]{29,25},new int[]{29,26},new int[]{29,27},new int[]{29,28},new int[]{29,29},new int[]{29,30},new int[]{29,31},new int[]{29,32},new int[]{29,35},new int[]{29,36},new int[]{29,37},new int[]{29,38},new int[]{29,39},new int[]{29,40},new int[]{29,41},new int[]{29,44},new int[]{29,45},new int[]{29,46},new int[]{29,47},new int[]{29,48},new int[]{29,53},new int[]{29,54},new int[]{29,55},new int[]{29,56},new int[]{29,57},new int[]{30,20},new int[]{30,21},new int[]{30,22},new int[]{30,23},new int[]{30,24},new int[]{30,25},new int[]{30,26},new int[]{30,27},new int[]{30,28},new int[]{30,29},new int[]{30,30},new int[]{30,31},new int[]{30,32},new int[]{30,35},new int[]{30,36},new int[]{30,37},new int[]{30,38},new int[]{30,39},new int[]{30,40},new int[]{30,41},new int[]{30,44},new int[]{30,45},new int[]{30,46},new int[]{30,47},new int[]{30,48},new int[]{30,53},new int[]{30,54},new int[]{30,55},new int[]{30,56},new int[]{30,57},new int[]{31,22},new int[]{31,23},new int[]{32,22},new int[]{32,26},new int[]{32,30},new int[]{32,31},new int[]{32,32},new int[]{32,35},new int[]{32,39},new int[]{32,40},new int[]{32,41},new int[]{32,44},new int[]{32,45},new int[]{32,46},new int[]{32,47},new int[]{32,48},new int[]{32,56},new int[]{32,57},new int[]{33,22},new int[]{33,26},new int[]{33,30},new int[]{33,31},new int[]{33,32},new int[]{33,35},new int[]{33,39},new int[]{33,40},new int[]{33,41},new int[]{33,44},new int[]{33,45},new int[]{33,46},new int[]{33,47},new int[]{33,48},new int[]{33,56},new int[]{33,57},new int[]{34,24},new int[]{35,22},new int[]{35,23},new int[]{35,24},new int[]{35,26},new int[]{35,27},new int[]{35,28},new int[]{35,29},new int[]{35,30},new int[]{35,31},new int[]{35,32},new int[]{35,35},new int[]{35,36},new int[]{35,37},new int[]{35,38},new int[]{35,39},new int[]{35,40},new int[]{35,41},new int[]{35,44},new int[]{35,45},new int[]{35,46},new int[]{35,47},new int[]{35,48},new int[]{35,53},new int[]{35,54},new int[]{35,55},new int[]{35,56},new int[]{35,57},new int[]{36,22},new int[]{36,23},new int[]{36,24},new int[]{36,25},new int[]{36,26},new int[]{36,27},new int[]{36,28},new int[]{36,29},new int[]{36,30},new int[]{36,31},new int[]{36,32},new int[]{36,35},new int[]{36,36},new int[]{36,37},new int[]{36,38},new int[]{36,39},new int[]{36,40},new int[]{36,41},new int[]{36,44},new int[]{36,45},new int[]{36,46},new int[]{36,47},new int[]{36,48},new int[]{36,53},new int[]{36,54},new int[]{36,55},new int[]{36,56},new int[]{36,57},new int[]{37,35},new int[]{37,36},new int[]{38,38},new int[]{39,37},new int[]{40,35},new int[]{43,39},new int[]{43,40},new int[]{44,26},new int[]{44,30},new int[]{44,31},new int[]{44,32},new int[]{44,39},new int[]{44,40},new int[]{44,44},new int[]{44,45},new int[]{44,46},new int[]{44,47},new int[]{44,48},new int[]{44,56},new int[]{44,57},new int[]{45,26},new int[]{45,30},new int[]{45,31},new int[]{45,32},new int[]{45,39},new int[]{45,40},new int[]{45,44},new int[]{45,45},new int[]{45,46},new int[]{45,47},new int[]{45,48},new int[]{45,56},new int[]{45,57},new int[]{46,46},new int[]{47,26},new int[]{47,27},new int[]{47,28},new int[]{47,29},new int[]{47,30},new int[]{47,31},new int[]{47,32},new int[]{47,44},new int[]{47,45},new int[]{47,47},new int[]{47,48},new int[]{47,55},new int[]{47,56},new int[]{48,26},new int[]{48,27},new int[]{48,28},new int[]{48,29},new int[]{48,30},new int[]{48,31},new int[]{48,32},new int[]{48,44},new int[]{48,45},new int[]{48,47},new int[]{48,48},new int[]{48,55},new int[]{48,56},new int[]{49,30},new int[]{49,31},new int[]{50,26},new int[]{50,30},new int[]{50,32},new int[]{51,26},new int[]{51,30},new int[]{51,32},new int[]{52,48},new int[]{53,26},new int[]{53,27},new int[]{53,28},new int[]{53,29},new int[]{53,30},new int[]{53,31},new int[]{53,32},new int[]{53,48},new int[]{53,55},new int[]{53,56},new int[]{54,26},new int[]{54,27},new int[]{54,28},new int[]{54,29},new int[]{54,30},new int[]{54,31},new int[]{54,32},new int[]{54,47},new int[]{54,48},new int[]{54,55},new int[]{54,56},new int[]{55,26},new int[]{55,27},new int[]{56,29},new int[]{57,28},new int[]{58,26}); }
		private static List<BigFraction> path4WC(BigFraction D) { return java.util.Arrays.asList(
			f(1), f(1).negate(), f(1), f(1).negate(), f(1).negate(), f(1).divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(3), f(2).negate(), f(1).divide(f(2)), D.negate().divide(D.subtract(f(1))), f(1).negate().divide(f(2)), f(1), f(1), D.negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), D.add(f(1)), D.negate(), D.negate(), D.divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), f(3).multiply(D), f(2).negate().multiply(D), D.divide(f(2)), D.pow(2).negate().divide(D.subtract(f(1))), D.negate().divide(f(2)), D, f(1), f(1).negate(), D.negate().divide(D.subtract(f(1))), D.pow(2).negate().divide(D.subtract(f(1))), D, D.divide(D.subtract(f(1))), D, D.negate().divide(D.subtract(f(1))), D.negate().divide(D.subtract(f(1))), f(3).negate().multiply(D), f(2).multiply(D), D.negate().divide(f(2)), D.pow(2).divide(D.subtract(f(1))), D.divide(f(2)), D.negate(), D.negate().divide(D.subtract(f(1))), D.negate().divide(D.subtract(f(1))), f(1), f(1).divide(D.subtract(f(1))), f(1), f(1).negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(3).negate(), f(2), f(1).negate().divide(f(2)), D.divide(D.subtract(f(1))), f(1).divide(f(2)), f(1).negate(), f(1).negate(), f(2).multiply(D), D.multiply(f(3).multiply(D).subtract(f(1))).divide(f(2).multiply(D.subtract(f(1)))), D.negate().divide(f(2)), D.pow(2).divide(D.subtract(f(1))), D.negate(), D.negate().divide(D.subtract(f(1))), D.negate(), D.divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), D.negate().divide(D.subtract(f(1))), f(1), D.negate(), D, f(2), f(3).multiply(D).subtract(f(1)).divide(f(2).multiply(D.subtract(f(1)))), f(1).negate().divide(f(2)), D.divide(D.subtract(f(1))), f(1).negate(), f(1).negate().divide(D.subtract(f(1))), f(1).negate(), f(1).divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).negate(), f(1), f(1).negate(), f(1), D.negate(), D.negate().divide(D.subtract(f(1))), D.negate().multiply(f(3).multiply(D).subtract(f(2))).divide(D.subtract(f(1))), f(2).multiply(D), f(1).negate().divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), f(2).multiply(D).divide(D.subtract(f(1))), f(2).negate().multiply(D).divide(D.subtract(f(1))), f(3).multiply(D), f(2).negate().multiply(D), D, f(2).negate().multiply(D.pow(2)).divide(D.subtract(f(1))), D.multiply(f(3).multiply(D).subtract(f(4))).divide(D.subtract(f(1))), f(4).negate().multiply(D), D.negate().divide(f(2)), f(2).negate().multiply(D.pow(2)).divide(D.subtract(f(1))), D.multiply(f(2).multiply(D).subtract(f(3))).divide(D.subtract(f(1))), f(2).negate().multiply(D), D.divide(D.subtract(f(1))), f(1), D.pow(2).negate().divide(D.subtract(f(1))), D, D.divide(D.subtract(f(1))), f(2).multiply(D).divide(D.subtract(f(1))), f(2).negate().multiply(D).divide(D.subtract(f(1))), f(2).multiply(D), f(2).negate().multiply(D.pow(2)).divide(D.subtract(f(1))), D, D, D.negate(), D.divide(f(2)), f(2).multiply(D), D.negate(), f(1).negate(), f(1).negate().divide(D.subtract(f(1))), f(3).multiply(D).subtract(f(2)).negate().divide(D.subtract(f(1))), f(2), f(1).negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(2).divide(D.subtract(f(1))), f(2).negate().divide(D.subtract(f(1))), f(3), f(2).negate(), f(1), f(2).negate().multiply(D).divide(D.subtract(f(1))), f(3).multiply(D).subtract(f(4)).divide(D.subtract(f(1))), f(4).negate(), f(1).negate().divide(f(2)), f(2).negate().multiply(D).divide(D.subtract(f(1))), f(2).multiply(D).subtract(f(3)).divide(D.subtract(f(1))), f(2).negate(), f(1).divide(D.subtract(f(1))), f(1).divide(D), D.negate().divide(D.subtract(f(1))), f(1), f(1).divide(D.subtract(f(1))), f(2).divide(D.subtract(f(1))), f(2).negate().divide(D.subtract(f(1))), f(2), f(2).negate().multiply(D).divide(D.subtract(f(1))), f(1), f(1), f(1).negate(), f(1).divide(f(2)), f(2), f(1).negate(), f(1).negate(), f(2).negate(), f(3).multiply(D).subtract(f(1)).negate().divide(f(2).multiply(D.subtract(f(1)))), f(1).divide(f(2)), D.negate().divide(D.subtract(f(1))), f(1), f(1).divide(D.subtract(f(1))), f(1), f(1).negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), f(1), f(1).negate(), f(2).negate().multiply(D), D.negate().multiply(f(3).multiply(D).subtract(f(1))).divide(f(2).multiply(D.subtract(f(1)))), D.divide(f(2)), D.pow(2).negate().divide(D.subtract(f(1))), D, D.divide(D.subtract(f(1))), D, D.negate().divide(D.subtract(f(1))), D.negate().divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), D.negate().divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), D, D.negate(), f(1).negate(), f(1), f(1), f(1), f(1).negate(), f(1).negate(), f(1), D, f(3).multiply(D), f(2).negate().multiply(D), f(2).negate().multiply(D).divide(D.subtract(f(1))), f(2).multiply(D).divide(D.subtract(f(1))), f(3).negate().multiply(D), f(2).multiply(D), D.negate(), f(2).multiply(D.pow(2)).divide(D.subtract(f(1))), D.negate().multiply(f(3).multiply(D).subtract(f(4))).divide(D.subtract(f(1))), f(4).multiply(D), D.divide(f(2)), f(2).multiply(D.pow(2)).divide(D.subtract(f(1))), D.negate().multiply(f(2).multiply(D).subtract(f(3))).divide(D.subtract(f(1))), f(2).multiply(D), D.negate().divide(D.subtract(f(1))), f(1).negate(), D.pow(2).divide(D.subtract(f(1))), D.negate(), D.negate().divide(D.subtract(f(1))), f(2).negate().multiply(D).divide(D.subtract(f(1))), f(2).multiply(D).divide(D.subtract(f(1))), f(2).negate().multiply(D), f(2).multiply(D.pow(2)).divide(D.subtract(f(1))), D.negate(), D.negate(), D, D.negate().divide(f(2)), f(2).negate().multiply(D), D, f(1), f(2), f(2).negate(), f(2).negate().divide(D.subtract(f(1))), f(2).divide(D.subtract(f(1))), f(3).negate(), f(2), f(1).negate(), f(2).multiply(D).divide(D.subtract(f(1))), f(3).multiply(D).subtract(f(4)).negate().divide(D.subtract(f(1))), f(4), f(1).divide(f(2)), f(2).multiply(D).divide(D.subtract(f(1))), f(2).multiply(D).subtract(f(3)).negate().divide(D.subtract(f(1))), f(2), f(1).negate().divide(D.subtract(f(1))), f(1).negate().divide(D), D.divide(D.subtract(f(1))), f(1).negate(), f(1).negate().divide(D.subtract(f(1))), f(2).negate().divide(D.subtract(f(1))), f(2).divide(D.subtract(f(1))), f(2).negate(), f(2).multiply(D).divide(D.subtract(f(1))), f(1).negate(), f(1).negate(), f(1), f(1).negate().divide(f(2)), f(2).negate(), f(1), f(1).negate(), D.negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(3).negate().multiply(D).divide(f(2)), D, D.negate().divide(f(2)), D.pow(2).divide(D.subtract(f(1))), D.negate().multiply(f(3).multiply(D).subtract(f(4))).divide(f(2).multiply(D.subtract(f(1)))), f(2).multiply(D), D.divide(f(4)), D.pow(2).divide(D.subtract(f(1))), D.negate().multiply(f(2).multiply(D).subtract(f(3))).divide(f(2).multiply(D.subtract(f(1)))), D, D.negate().divide(f(2).multiply(D.subtract(f(1)))), f(3).negate().multiply(D).divide(f(2)), D, D.negate().divide(f(2)), D, D.pow(2).divide(f(2).multiply(D.subtract(f(1)))), D.negate().divide(f(2)), D.negate().divide(f(2).multiply(D.subtract(f(1)))), D.negate().divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), D.negate(), D.pow(2).divide(D.subtract(f(1))), D.negate().divide(f(2)), D.divide(f(2)), D.divide(f(2)), D.negate().divide(f(4)), D.negate(), D.divide(f(2)), f(1).negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(3).negate().divide(f(2)), f(1), f(1).negate().divide(f(2)), D.divide(D.subtract(f(1))), f(3).multiply(D).subtract(f(4)).negate().divide(f(2).multiply(D.subtract(f(1)))), f(2), f(1).divide(f(4)), D.divide(D.subtract(f(1))), f(2).multiply(D).subtract(f(3)).negate().divide(f(2).multiply(D.subtract(f(1)))), f(1), f(1).negate().divide(f(2).multiply(D.subtract(f(1)))), f(3).negate().divide(f(2)), f(1), f(1).negate().divide(f(2)), f(1), D.divide(f(2).multiply(D.subtract(f(1)))), f(1).negate().divide(f(2)), f(1).negate().divide(f(2).multiply(D.subtract(f(1)))), f(1).negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1).negate(), D.divide(D.subtract(f(1))), f(1).negate().divide(f(2)), f(1).divide(f(2)), f(1).divide(f(2)), f(1).negate().divide(f(4)), f(1).negate(), f(1).divide(f(2)), f(1).negate(), f(1), f(1), D.negate().divide(D.subtract(f(1))), D.multiply(f(2).multiply(D).subtract(f(3))).divide(D.subtract(f(1))), f(2).negate().multiply(D), D.divide(D.subtract(f(1))), f(1), D.pow(2).negate().divide(D.subtract(f(1))), D, D.divide(D.subtract(f(1))), f(2).multiply(D).divide(D.subtract(f(1))), f(2).negate().multiply(D).divide(D.subtract(f(1))), f(2).multiply(D), f(2).negate().multiply(D.pow(2)).divide(D.subtract(f(1))), D, f(2).multiply(D), D.negate(), f(1), f(1).negate().divide(D.subtract(f(1))), f(2).multiply(D).subtract(f(3)).divide(D.subtract(f(1))), f(2).negate(), f(1).divide(D.subtract(f(1))), f(1).divide(D), D.negate().divide(D.subtract(f(1))), f(1), f(1).divide(D.subtract(f(1))), f(2).divide(D.subtract(f(1))), f(2).negate().divide(D.subtract(f(1))), f(2), f(2).negate().multiply(D).divide(D.subtract(f(1))), f(1), f(2), f(1).negate(), f(1).negate(), f(3).divide(f(2)), f(1).negate(), f(1).divide(f(2)), f(3).multiply(D).subtract(f(4)).divide(f(2).multiply(D.subtract(f(1)))), f(2).negate(), f(1).negate().divide(f(4)), f(1).negate(), f(2).multiply(D).subtract(f(3)).divide(f(2).multiply(D.subtract(f(1)))), f(1).negate(), f(1).divide(f(2).multiply(D.subtract(f(1)))), f(3).divide(f(2)), f(1).negate(), f(1).divide(f(2)), f(1).negate(), D.negate().divide(f(2).multiply(D.subtract(f(1)))), f(1).divide(f(2)), f(1).divide(f(2).multiply(D.subtract(f(1)))), f(1).divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1), D.negate().divide(D.subtract(f(1))), f(1).divide(f(2)), f(1).negate().divide(f(2)), f(1).negate().divide(f(2)), f(1).divide(f(4)), f(1), f(1).negate().divide(f(2)), f(3).multiply(D).divide(f(2)), D.negate(), D.divide(f(2)), D.negate(), D.multiply(f(3).multiply(D).subtract(f(4))).divide(f(2).multiply(D.subtract(f(1)))), f(2).negate().multiply(D), D.negate().divide(f(4)), D.negate(), D.multiply(f(2).multiply(D).subtract(f(3))).divide(f(2).multiply(D.subtract(f(1)))), D.negate(), D.divide(f(2).multiply(D.subtract(f(1)))), f(3).multiply(D).divide(f(2)), D.negate(), D.divide(f(2)), D.negate(), D.pow(2).negate().divide(f(2).multiply(D.subtract(f(1)))), D.divide(f(2)), D.divide(f(2).multiply(D.subtract(f(1)))), D.divide(D.subtract(f(1))), D.negate().divide(D.subtract(f(1))), D, D.pow(2).negate().divide(D.subtract(f(1))), D.divide(f(2)), D.negate().divide(f(2)), D.negate().divide(f(2)), D.divide(f(4)), D, D.negate().divide(f(2)), f(1).negate(), f(1), f(1), f(1).negate(), D.subtract(f(1)).divide(D), f(1).negate(), f(1), D.divide(D.subtract(f(1))), D.negate().multiply(f(2).multiply(D).subtract(f(3))).divide(D.subtract(f(1))), f(2).multiply(D), D.negate().divide(D.subtract(f(1))), D.add(f(1)), D.negate(), f(2).negate().multiply(D).divide(D.subtract(f(1))), f(2).multiply(D).divide(D.subtract(f(1))), f(2).negate().multiply(D), f(2).multiply(D.pow(2)).divide(D.subtract(f(1))), D.negate(), f(2).negate().multiply(D), D, f(1).divide(D.subtract(f(1))), f(2).multiply(D).subtract(f(3)).negate().divide(D.subtract(f(1))), f(2), f(1).negate().divide(D.subtract(f(1))), f(1), f(1).negate(), f(2).negate().divide(D.subtract(f(1))), f(2).divide(D.subtract(f(1))), f(2).negate(), f(2).multiply(D).divide(D.subtract(f(1))), f(1).negate(), f(2).negate(), f(1), f(1), D.negate().multiply(f(3).multiply(D).subtract(f(4))).divide(f(2).multiply(D.subtract(f(1)))), f(2).multiply(D), D.divide(f(4)), D, D.negate().multiply(f(2).multiply(D).subtract(f(3))).divide(f(2).multiply(D.subtract(f(1)))), D, D.negate().divide(f(2).multiply(D.subtract(f(1)))), D.negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), D.pow(2).divide(D.subtract(f(1))), D.negate().divide(f(2)), D.negate().divide(f(4)), D.negate(), f(3).multiply(D).subtract(f(4)).negate().divide(f(2).multiply(D.subtract(f(1)))), f(2), f(1).divide(f(4)), f(1), f(2).multiply(D).subtract(f(3)).negate().divide(f(2).multiply(D.subtract(f(1)))), f(1), f(1).negate().divide(f(2).multiply(D.subtract(f(1)))), f(1).negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), f(1).negate().divide(f(2)), f(1).negate().divide(f(4)), f(1).negate(), f(1).negate(), f(1), D.negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1).negate(), f(3).multiply(D).subtract(f(4)).divide(f(2).multiply(D.subtract(f(1)))), f(2).negate(), f(1).negate().divide(f(4)), f(1).negate(), f(2).multiply(D).subtract(f(3)).divide(f(2).multiply(D.subtract(f(1)))), f(1).negate(), f(1).divide(f(2).multiply(D.subtract(f(1)))), f(1).divide(f(2)), f(1).divide(f(4)), f(1), D.multiply(f(3).multiply(D).subtract(f(4))).divide(f(2).multiply(D.subtract(f(1)))), f(2).negate().multiply(D), D.negate().divide(f(4)), D.negate(), D.multiply(f(2).multiply(D).subtract(f(3))).divide(f(2).multiply(D.subtract(f(1)))), D.negate(), D.divide(f(2).multiply(D.subtract(f(1)))), D.negate(), D.divide(f(2)), D.divide(f(4)), D, f(1).negate(), f(1), f(1), f(1), f(1).negate()); }
		private static List<int[]> path3UG() { return java.util.Arrays.asList(
			new int[]{1,5},new int[]{1,9},new int[]{1,10},new int[]{1,37},new int[]{2,1},new int[]{2,2},new int[]{2,3},new int[]{2,4},new int[]{2,5},new int[]{2,6},new int[]{2,8},new int[]{2,35},new int[]{2,37},new int[]{6,1},new int[]{6,2},new int[]{6,3},new int[]{6,4},new int[]{6,5},new int[]{6,6},new int[]{6,35},new int[]{7,3},new int[]{7,26},new int[]{7,27},new int[]{7,35},new int[]{8,1},new int[]{8,2},new int[]{8,3},new int[]{8,35},new int[]{13,18},new int[]{13,27},new int[]{13,33},new int[]{14,7},new int[]{14,11},new int[]{14,12},new int[]{14,15},new int[]{14,16},new int[]{14,17},new int[]{14,18},new int[]{14,19},new int[]{14,38},new int[]{15,7},new int[]{15,11},new int[]{15,12},new int[]{15,15},new int[]{15,16},new int[]{15,17},new int[]{15,18},new int[]{15,19},new int[]{15,34},new int[]{15,38},new int[]{19,9},new int[]{19,12},new int[]{19,14},new int[]{19,36},new int[]{21,11},new int[]{21,12},new int[]{21,13},new int[]{21,36},new int[]{25,16},new int[]{25,20},new int[]{25,21},new int[]{25,38},new int[]{26,7},new int[]{26,15},new int[]{26,16},new int[]{26,38},new int[]{31,20},new int[]{31,23},new int[]{31,25},new int[]{31,40},new int[]{32,22},new int[]{32,23},new int[]{32,28},new int[]{32,29},new int[]{32,32},new int[]{33,22},new int[]{33,23},new int[]{33,24},new int[]{33,28},new int[]{33,29},new int[]{33,32},new int[]{33,40},new int[]{37,27},new int[]{37,29},new int[]{37,31},new int[]{37,39},new int[]{39,28},new int[]{39,29},new int[]{39,30},new int[]{39,39}); }
		private static List<BigFraction> path3UC(BigFraction D) { return java.util.Arrays.asList(
			f(1).negate(), f(1), f(1).negate(), f(1).negate(), f(1).negate(), f(1).negate(), f(1).negate(), D.negate(), D.negate(), f(1).negate(), f(1), f(1).negate(), f(1), f(1), f(1), f(1), f(1), f(1), f(1), f(1), f(1).negate(), f(1).negate(), f(1).negate(), f(1).divide(f(2)), f(1).negate(), D.subtract(f(1)).negate().divide(D), D.subtract(f(1)).negate().divide(D), f(1).negate(), f(1), f(1).negate(), f(1), f(1), f(1).divide(D), f(1).divide(D), f(1), f(1), f(1), f(1), f(1), f(1), f(1).negate(), f(1).negate().divide(D), f(1).negate().divide(D), f(1).negate(), f(1).negate(), D.negate(), D.negate(), f(1).negate(), f(1), f(1).negate(), f(1), f(1), f(1), f(2).negate(), D.subtract(f(1)).negate().divide(D), D.subtract(f(1)).negate().divide(D), f(1).negate(), f(1).negate(), f(1), f(1), f(1), f(2).negate(), f(1).negate(), D.subtract(f(1)).negate().divide(D), D.subtract(f(1)).negate().divide(D), f(1).negate(), f(1), f(1), f(1), f(2).negate(), f(1).negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1), D.divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), f(1), f(1).divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1).negate(), f(1), f(1).negate(), f(1), f(1), f(2).negate(), f(1), f(1), f(1), f(1)); }
		private static List<int[]> path3VG() { return java.util.Arrays.asList(
			new int[]{1,4},new int[]{1,5},new int[]{1,6},new int[]{1,37},new int[]{3,8},new int[]{3,9},new int[]{3,10},new int[]{3,14},new int[]{3,37},new int[]{4,9},new int[]{4,10},new int[]{4,14},new int[]{7,2},new int[]{7,3},new int[]{7,6},new int[]{7,35},new int[]{10,1},new int[]{10,26},new int[]{10,27},new int[]{10,31},new int[]{10,33},new int[]{10,35},new int[]{11,26},new int[]{11,27},new int[]{11,31},new int[]{11,33},new int[]{13,17},new int[]{13,18},new int[]{13,19},new int[]{16,33},new int[]{16,34},new int[]{19,11},new int[]{19,12},new int[]{19,19},new int[]{19,36},new int[]{20,13},new int[]{20,14},new int[]{20,36},new int[]{25,15},new int[]{25,16},new int[]{25,19},new int[]{25,38},new int[]{28,7},new int[]{28,20},new int[]{28,21},new int[]{28,25},new int[]{28,38},new int[]{29,20},new int[]{29,21},new int[]{29,25},new int[]{31,22},new int[]{31,23},new int[]{31,32},new int[]{31,40},new int[]{34,24},new int[]{34,25},new int[]{34,40},new int[]{37,28},new int[]{37,29},new int[]{37,32},new int[]{37,39},new int[]{38,30},new int[]{38,31},new int[]{38,39}); }
		private static List<BigFraction> path3VC(BigFraction D) { return java.util.Arrays.asList(
			f(1), f(1), f(1), f(1), f(1), f(1).negate(), D.negate(), f(1).negate().divide(D), f(1), f(1), f(1), f(1).divide(D), f(1), f(1), f(1), f(1), f(1), D.negate(), f(1).negate(), f(1).negate().divide(D), f(1).negate().divide(D), f(2), f(1), f(1), f(1).divide(D), f(1).divide(D), f(1), f(1), f(1), D.subtract(f(1)).negate().divide(D), f(1).negate(), f(1), f(1), f(1), f(2), f(1).negate(), D.subtract(f(1)).negate().divide(D), f(1).negate(), f(1), f(1), f(1), f(1), f(1), f(1).negate(), D.negate(), f(1).negate().divide(D), f(1), f(1), f(1), f(1).divide(D), f(1), f(1), f(1), f(1), f(1).negate(), D.subtract(f(1)).negate().divide(D), f(1).negate(), f(1), f(1), f(1), f(1), f(1).negate(), D.subtract(f(1)).negate().divide(D), f(1).negate()); }
		private static List<int[]> path3WG() { return java.util.Arrays.asList(
			new int[]{1,4},new int[]{1,5},new int[]{2,4},new int[]{2,5},new int[]{2,8},new int[]{2,9},new int[]{2,10},new int[]{2,11},new int[]{2,12},new int[]{2,13},new int[]{2,14},new int[]{2,36},new int[]{2,37},new int[]{3,8},new int[]{4,9},new int[]{4,10},new int[]{4,14},new int[]{5,9},new int[]{5,10},new int[]{5,14},new int[]{6,4},new int[]{6,5},new int[]{6,8},new int[]{6,9},new int[]{6,10},new int[]{6,11},new int[]{6,12},new int[]{6,13},new int[]{6,14},new int[]{6,36},new int[]{6,37},new int[]{7,2},new int[]{7,3},new int[]{8,2},new int[]{8,4},new int[]{8,5},new int[]{8,6},new int[]{8,8},new int[]{8,9},new int[]{8,10},new int[]{8,11},new int[]{8,12},new int[]{8,13},new int[]{8,14},new int[]{8,36},new int[]{8,37},new int[]{9,2},new int[]{9,4},new int[]{9,5},new int[]{9,6},new int[]{9,8},new int[]{9,9},new int[]{9,10},new int[]{9,11},new int[]{9,12},new int[]{9,13},new int[]{9,14},new int[]{9,36},new int[]{9,37},new int[]{10,1},new int[]{11,1},new int[]{11,2},new int[]{11,3},new int[]{11,4},new int[]{11,5},new int[]{11,6},new int[]{11,8},new int[]{11,9},new int[]{11,10},new int[]{11,14},new int[]{11,26},new int[]{11,35},new int[]{11,37},new int[]{12,1},new int[]{12,2},new int[]{12,3},new int[]{12,4},new int[]{12,5},new int[]{12,6},new int[]{12,8},new int[]{12,9},new int[]{12,10},new int[]{12,14},new int[]{12,35},new int[]{12,37},new int[]{13,17},new int[]{13,18},new int[]{14,7},new int[]{14,11},new int[]{14,15},new int[]{14,16},new int[]{14,17},new int[]{14,19},new int[]{14,20},new int[]{14,21},new int[]{14,22},new int[]{14,23},new int[]{14,24},new int[]{14,25},new int[]{14,28},new int[]{14,32},new int[]{14,38},new int[]{14,40},new int[]{15,7},new int[]{15,11},new int[]{15,15},new int[]{15,16},new int[]{15,17},new int[]{15,19},new int[]{15,20},new int[]{15,21},new int[]{15,22},new int[]{15,23},new int[]{15,24},new int[]{15,25},new int[]{15,28},new int[]{15,32},new int[]{15,38},new int[]{15,40},new int[]{16,34},new int[]{17,1},new int[]{17,2},new int[]{17,3},new int[]{17,4},new int[]{17,5},new int[]{17,6},new int[]{17,8},new int[]{17,9},new int[]{17,10},new int[]{17,14},new int[]{17,26},new int[]{17,27},new int[]{17,31},new int[]{17,33},new int[]{17,35},new int[]{17,37},new int[]{18,1},new int[]{18,2},new int[]{18,3},new int[]{18,4},new int[]{18,5},new int[]{18,6},new int[]{18,8},new int[]{18,9},new int[]{18,10},new int[]{18,14},new int[]{18,26},new int[]{18,27},new int[]{18,31},new int[]{18,33},new int[]{18,35},new int[]{18,37},new int[]{19,11},new int[]{19,12},new int[]{20,14},new int[]{21,13},new int[]{22,11},new int[]{25,15},new int[]{25,16},new int[]{26,7},new int[]{26,15},new int[]{26,16},new int[]{26,20},new int[]{26,21},new int[]{26,22},new int[]{26,23},new int[]{26,24},new int[]{26,25},new int[]{26,28},new int[]{26,32},new int[]{26,38},new int[]{26,40},new int[]{27,7},new int[]{27,15},new int[]{27,16},new int[]{27,20},new int[]{27,21},new int[]{27,22},new int[]{27,23},new int[]{27,24},new int[]{27,25},new int[]{27,28},new int[]{27,32},new int[]{27,38},new int[]{27,40},new int[]{28,7},new int[]{29,20},new int[]{29,21},new int[]{29,22},new int[]{29,23},new int[]{29,24},new int[]{29,25},new int[]{29,28},new int[]{29,29},new int[]{29,30},new int[]{29,31},new int[]{29,32},new int[]{29,39},new int[]{29,40},new int[]{30,20},new int[]{30,21},new int[]{30,22},new int[]{30,23},new int[]{30,24},new int[]{30,25},new int[]{30,28},new int[]{30,29},new int[]{30,30},new int[]{30,31},new int[]{30,32},new int[]{30,39},new int[]{30,40},new int[]{31,22},new int[]{31,23},new int[]{32,22},new int[]{32,28},new int[]{32,32},new int[]{33,22},new int[]{33,28},new int[]{33,32},new int[]{34,24},new int[]{35,22},new int[]{35,23},new int[]{35,24},new int[]{35,28},new int[]{35,29},new int[]{35,30},new int[]{35,31},new int[]{35,32},new int[]{35,39},new int[]{35,40},new int[]{36,22},new int[]{36,23},new int[]{36,24},new int[]{36,25},new int[]{36,28},new int[]{36,29},new int[]{36,30},new int[]{36,31},new int[]{36,32},new int[]{36,39},new int[]{36,40},new int[]{37,28},new int[]{37,29},new int[]{38,31},new int[]{39,30},new int[]{40,28}); }
		private static List<BigFraction> path3WC(BigFraction D) { return java.util.Arrays.asList(
			f(1), f(1).negate(), f(1), f(1).negate(), f(1).negate(), f(1).divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(3), f(2).negate(), f(1).divide(f(2)), D.negate().divide(D.subtract(f(1))), f(1).negate().divide(f(2)), f(1), f(1), D.negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), D.add(f(1)), D.negate(), D.negate(), D.divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), f(3).multiply(D), f(2).negate().multiply(D), D.divide(f(2)), D.pow(2).negate().divide(D.subtract(f(1))), D.negate().divide(f(2)), D, f(1), f(1).negate(), D.negate().divide(D.subtract(f(1))), D.pow(2).negate().divide(D.subtract(f(1))), D, D.divide(D.subtract(f(1))), D, D.negate().divide(D.subtract(f(1))), D.negate().divide(D.subtract(f(1))), f(3).negate().multiply(D), f(2).multiply(D), D.negate().divide(f(2)), D.pow(2).divide(D.subtract(f(1))), D.divide(f(2)), D.negate(), D.negate().divide(D.subtract(f(1))), D.negate().divide(D.subtract(f(1))), f(1), f(1).divide(D.subtract(f(1))), f(1), f(1).negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(3).negate(), f(2), f(1).negate().divide(f(2)), D.divide(D.subtract(f(1))), f(1).divide(f(2)), f(1).negate(), f(1).negate(), f(2).multiply(D), D.multiply(f(3).multiply(D).subtract(f(1))).divide(f(2).multiply(D.subtract(f(1)))), D.negate().divide(f(2)), D.pow(2).divide(D.subtract(f(1))), D.negate(), D.negate().divide(D.subtract(f(1))), D.negate(), D.divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), D.negate().divide(D.subtract(f(1))), f(1), D.negate(), D, f(2), f(3).multiply(D).subtract(f(1)).divide(f(2).multiply(D.subtract(f(1)))), f(1).negate().divide(f(2)), D.divide(D.subtract(f(1))), f(1).negate(), f(1).negate().divide(D.subtract(f(1))), f(1).negate(), f(1).divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).negate(), f(1), f(1).negate(), f(1), D.negate(), D.negate().divide(D.subtract(f(1))), D.negate().multiply(f(3).multiply(D).subtract(f(2))).divide(D.subtract(f(1))), f(2).multiply(D), f(1).negate().divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), f(2).multiply(D).divide(D.subtract(f(1))), f(2).negate().multiply(D).divide(D.subtract(f(1))), f(3).multiply(D), f(2).negate().multiply(D), D, f(2).negate().multiply(D.pow(2)).divide(D.subtract(f(1))), f(1), D.divide(D.subtract(f(1))), D, D.negate(), f(1).negate(), f(1).negate().divide(D.subtract(f(1))), f(3).multiply(D).subtract(f(2)).negate().divide(D.subtract(f(1))), f(2), f(1).negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(2).divide(D.subtract(f(1))), f(2).negate().divide(D.subtract(f(1))), f(3), f(2).negate(), f(1), f(2).negate().multiply(D).divide(D.subtract(f(1))), f(1).divide(D), f(1).divide(D.subtract(f(1))), f(1), f(1).negate(), f(1).negate(), f(2).negate(), f(3).multiply(D).subtract(f(1)).negate().divide(f(2).multiply(D.subtract(f(1)))), f(1).divide(f(2)), D.negate().divide(D.subtract(f(1))), f(1), f(1).divide(D.subtract(f(1))), f(1), f(1).negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), f(1), f(1).negate(), f(2).negate().multiply(D), D.negate().multiply(f(3).multiply(D).subtract(f(1))).divide(f(2).multiply(D.subtract(f(1)))), D.divide(f(2)), D.pow(2).negate().divide(D.subtract(f(1))), D, D.divide(D.subtract(f(1))), D, D.negate().divide(D.subtract(f(1))), D.negate().divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), D.negate().divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), D, D.negate(), f(1).negate(), f(1), f(1), f(1), f(1).negate(), f(1).negate(), f(1), D, f(3).multiply(D), f(2).negate().multiply(D), f(2).negate().multiply(D).divide(D.subtract(f(1))), f(2).multiply(D).divide(D.subtract(f(1))), f(3).negate().multiply(D), f(2).multiply(D), D.negate(), f(2).multiply(D.pow(2)).divide(D.subtract(f(1))), f(1).negate(), D.negate().divide(D.subtract(f(1))), D.negate(), D, f(1), f(2), f(2).negate(), f(2).negate().divide(D.subtract(f(1))), f(2).divide(D.subtract(f(1))), f(3).negate(), f(2), f(1).negate(), f(2).multiply(D).divide(D.subtract(f(1))), f(1).negate().divide(D), f(1).negate().divide(D.subtract(f(1))), f(1).negate(), f(1), f(1).negate(), D.negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(3).negate().multiply(D).divide(f(2)), D, D.negate().divide(f(2)), D.pow(2).divide(D.subtract(f(1))), f(3).negate().multiply(D).divide(f(2)), D, D.negate().divide(f(2)), D, D.negate().divide(f(2).multiply(D.subtract(f(1)))), D.divide(f(2)), D.divide(f(2)), f(1).negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(3).negate().divide(f(2)), f(1), f(1).negate().divide(f(2)), D.divide(D.subtract(f(1))), f(3).negate().divide(f(2)), f(1), f(1).negate().divide(f(2)), f(1), f(1).negate().divide(f(2).multiply(D.subtract(f(1)))), f(1).divide(f(2)), f(1).divide(f(2)), f(1).negate(), f(1), f(1), f(1), D.divide(D.subtract(f(1))), f(1), f(1).divide(D), f(1).divide(D.subtract(f(1))), f(1).negate(), f(3).divide(f(2)), f(1).negate(), f(1).divide(f(2)), f(3).divide(f(2)), f(1).negate(), f(1).divide(f(2)), f(1).negate(), f(1).divide(f(2).multiply(D.subtract(f(1)))), f(1).negate().divide(f(2)), f(1).negate().divide(f(2)), f(3).multiply(D).divide(f(2)), D.negate(), D.divide(f(2)), D.negate(), f(3).multiply(D).divide(f(2)), D.negate(), D.divide(f(2)), D.negate(), D.divide(f(2).multiply(D.subtract(f(1)))), D.negate().divide(f(2)), D.negate().divide(f(2)), f(1).negate(), f(1), f(1), f(1).negate(), D.subtract(f(1)).divide(D)); }
		private static List<int[]> cycleUG() { return java.util.Arrays.asList(
			new int[]{6,1},new int[]{6,2},new int[]{6,4},new int[]{6,5},new int[]{6,6},new int[]{14,10},new int[]{14,12},new int[]{14,14},new int[]{14,16},new int[]{14,22},new int[]{24,10},new int[]{24,12},new int[]{24,22},new int[]{8,1},new int[]{8,2},new int[]{5,4},new int[]{5,7},new int[]{5,8},new int[]{5,23},new int[]{13,14},new int[]{13,15},new int[]{13,17},new int[]{13,20},new int[]{13,21},new int[]{7,2},new int[]{7,3},new int[]{7,17},new int[]{7,18},new int[]{19,7},new int[]{19,11},new int[]{19,12},new int[]{19,13},new int[]{19,22},new int[]{15,10},new int[]{15,12},new int[]{15,14},new int[]{15,16},new int[]{15,19},new int[]{15,21},new int[]{15,22},new int[]{2,1},new int[]{2,2},new int[]{2,4},new int[]{2,5},new int[]{2,6},new int[]{2,9},new int[]{2,23}); }
		private static List<BigFraction> cycleUC(BigFraction D) { return java.util.Arrays.asList(
			f(1), D.negate(), f(1).divide(D), f(1).divide(D), f(1), f(1).negate(), D.negate(), f(1).negate(), f(1).negate(), f(1), f(1), D.subtract(f(1)), f(1).negate(), f(1).negate(), D.subtract(f(1)), f(1), f(1).negate(), f(1).negate(), f(1), f(1).negate(), f(1).negate(), f(1).negate(), f(1).negate(), f(1), f(1).negate(), f(1).negate(), f(1).negate(), f(1), f(1), D, f(1).negate(), f(1).negate(), f(1), f(1), D, D, f(1), f(1), f(1), f(1).negate(), f(1).negate(), D, f(1).negate(), f(1).negate(), f(1).negate(), f(1).negate(), f(1)); }
		private static List<int[]> cycleVG() { return java.util.Arrays.asList(
			new int[]{4,7},new int[]{4,8},new int[]{4,10},new int[]{4,11},new int[]{4,22},new int[]{11,17},new int[]{11,18},new int[]{11,19},new int[]{11,20},new int[]{11,21},new int[]{18,19},new int[]{18,20},new int[]{18,21},new int[]{23,10},new int[]{23,11},new int[]{23,22},new int[]{6,4},new int[]{6,5},new int[]{6,6},new int[]{6,23},new int[]{7,2},new int[]{7,3},new int[]{7,6},new int[]{13,14},new int[]{13,15},new int[]{13,16},new int[]{13,20},new int[]{13,21},new int[]{19,11},new int[]{19,12},new int[]{19,13},new int[]{19,16},new int[]{19,22},new int[]{10,1},new int[]{10,17},new int[]{10,18},new int[]{10,19},new int[]{10,20},new int[]{10,21},new int[]{5,7},new int[]{5,8},new int[]{5,9},new int[]{5,10},new int[]{5,11},new int[]{5,22},new int[]{5,23}); }
		private static List<BigFraction> cycleVC(BigFraction D) { return java.util.Arrays.asList(
			f(1), f(1), f(1), f(1), f(1), f(1), f(1), f(1), f(1), f(1), f(1).negate(), D.subtract(f(1)).negate().divide(D), f(1).negate(), f(1).negate(), D.subtract(f(1)).negate().divide(D), f(1).negate(), f(1), f(1), f(1), f(1), f(1), f(1), f(1), f(1), f(1), f(1), f(1), f(1), f(1).negate().divide(D), f(1), f(1), f(1), f(1).negate(), f(1), f(1).negate(), D.negate(), f(1).negate(), f(1).negate(), f(1).negate(), f(1).negate(), D.negate(), f(1), f(1).negate(), f(1).negate(), f(1).negate(), f(1)); }
		private static List<int[]> cycleWG() { return java.util.Arrays.asList(
			new int[]{3,9},new int[]{16,19},new int[]{10,1},new int[]{21,10},new int[]{5,7},new int[]{5,8},new int[]{5,10},new int[]{5,11},new int[]{5,12},new int[]{5,13},new int[]{5,14},new int[]{5,15},new int[]{5,16},new int[]{5,19},new int[]{5,20},new int[]{5,21},new int[]{5,22},new int[]{9,2},new int[]{9,3},new int[]{9,5},new int[]{9,6},new int[]{15,4},new int[]{15,5},new int[]{15,7},new int[]{15,8},new int[]{15,9},new int[]{15,10},new int[]{15,11},new int[]{15,12},new int[]{15,13},new int[]{15,14},new int[]{15,15},new int[]{15,16},new int[]{15,22},new int[]{15,23},new int[]{12,15},new int[]{12,17},new int[]{12,18},new int[]{12,20},new int[]{17,15},new int[]{17,20},new int[]{2,5},new int[]{20,10},new int[]{20,12},new int[]{20,13},new int[]{20,14},new int[]{20,15},new int[]{20,16},new int[]{20,19},new int[]{20,20},new int[]{20,21},new int[]{20,22},new int[]{22,4},new int[]{22,5},new int[]{22,7},new int[]{22,8},new int[]{22,9},new int[]{22,10},new int[]{22,11},new int[]{22,13},new int[]{22,22},new int[]{22,23},new int[]{7,3},new int[]{13,15},new int[]{1,4},new int[]{1,5},new int[]{19,13},new int[]{23,10},new int[]{23,11},new int[]{23,12},new int[]{23,13},new int[]{23,14},new int[]{23,15},new int[]{23,16},new int[]{23,19},new int[]{23,20},new int[]{23,21},new int[]{23,22},new int[]{11,15},new int[]{11,17},new int[]{11,18},new int[]{11,20},new int[]{8,2},new int[]{8,3},new int[]{8,5},new int[]{8,6},new int[]{24,4},new int[]{24,5},new int[]{24,7},new int[]{24,8},new int[]{24,9},new int[]{24,10},new int[]{24,11},new int[]{24,12},new int[]{24,13},new int[]{24,22},new int[]{24,23},new int[]{4,7},new int[]{4,8},new int[]{4,10},new int[]{4,11},new int[]{4,12},new int[]{4,13},new int[]{4,14},new int[]{4,15},new int[]{4,16},new int[]{4,19},new int[]{4,20},new int[]{4,21},new int[]{4,22},new int[]{14,4},new int[]{14,5},new int[]{14,7},new int[]{14,8},new int[]{14,9},new int[]{14,10},new int[]{14,11},new int[]{14,12},new int[]{14,13},new int[]{14,14},new int[]{14,15},new int[]{14,16},new int[]{14,22},new int[]{14,23}); }
		private static List<BigFraction> cycleWC(BigFraction D) { return java.util.Arrays.asList(
			f(1).negate(), f(1).negate(), f(1).negate(), f(1).negate(), f(1).negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1).negate(), f(1).divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), D.subtract(f(3)).divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1), f(1).negate(), f(1).negate(), f(1).negate(), f(1).divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1), f(2).negate(), f(1).negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1).negate(), f(1).negate(), f(1).divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), D.subtract(f(3)).divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).negate(), f(1).negate(), f(1).divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1), f(1).negate(), f(1).negate(), f(1), f(1).negate().divide(D.subtract(f(1))), D.subtract(f(2)).negate().divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1).negate(), f(1), f(1), f(1), f(1).negate(), f(2), f(1).divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), f(1), f(1), f(1).negate().divide(D.subtract(f(1))), D.subtract(f(2)).negate().divide(D.subtract(f(1))), f(1), f(1), f(1).negate(), f(1).negate(), f(1).negate(), f(1), f(1).negate(), D, f(1).negate(), D.negate().divide(D.subtract(f(1))), D.pow(2).subtract(f(3).multiply(D)).add(f(1)).negate().divide(D.subtract(f(1))), D.negate().divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), D.negate(), D, D, D, D.divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), D.negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), D.negate().divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), D.negate(), f(2).multiply(D), D.divide(D.subtract(f(1))), D.negate().divide(D.subtract(f(1))), D, D, D.negate().divide(D.subtract(f(1))), f(1).negate(), D.pow(2).subtract(f(3).multiply(D)).add(f(1)).negate().divide(D.subtract(f(1))), D, D, D.negate().divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), D.negate(), D.divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), D.multiply(D.subtract(f(3))).divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), D.negate().divide(D.subtract(f(1))), D.negate().divide(D.subtract(f(1))), D, D.negate(), D.negate(), D.negate(), D, f(2).negate().multiply(D), D.negate().divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), D.negate(), D.negate(), D.divide(D.subtract(f(1))), D.divide(D.subtract(f(1))), D.multiply(D.subtract(f(3))).divide(D.subtract(f(1))), f(1).divide(D.subtract(f(1))), f(1).negate().divide(D.subtract(f(1))), D.negate().divide(D.subtract(f(1))), D.negate(), D.negate()); }
	}

	// ===================================================================
	// EVEN construction
	// ===================================================================
	private static final class Even {
		final int N; // matmul dimension (= 2n)
		final int n;
		final int d;
		final int twoD;
		final BigFraction[][] PL; // (2d) x N
		final BigFraction[][] PR; // N x (2d)
		final Map<Integer, Map<Integer, BigFraction>> phiCache = new HashMap<>();

		Even(int nn) {
			this.N = nn;
			this.n = nn / 2;
			this.d = n + 1;
			this.twoD = 2 * d;
			BigFraction[][] L = new BigFraction[d][n];
			for (int r = 0; r < d; r++) for (int c = 0; c < n; c++) L[r][c] = ZERO;
			for (int r = 0; r < n; r++) L[r][r] = ONE;
			for (int c = 0; c < n; c++) L[n][c] = f(-1);
			BigFraction[][] R = new BigFraction[n][d];
			BigFraction invd = f(1).divide(f(d));
			for (int r = 0; r < n; r++) {
				for (int c = 0; c < n; c++) R[r][c] = (r == c ? ONE : ZERO).subtract(invd);
				R[r][n] = invd.negate();
			}
			PL = new BigFraction[twoD][N];
			for (int p = 0; p < twoD; p++) for (int i = 0; i < N; i++) PL[p][i] = ZERO;
			for (int p = 0; p < twoD; p++) {
				if (p < d) for (int i = 0; i < n; i++) PL[p][i] = L[p][i];
				else for (int i = 0; i < n; i++) PL[p][n + i] = L[p - d][i];
			}
			PR = new BigFraction[N][twoD];
			for (int j = 0; j < N; j++) for (int q = 0; q < twoD; q++) PR[j][q] = ZERO;
			for (int j = 0; j < N; j++) {
				if (j < n) for (int q = 0; q < d; q++) PR[j][q] = R[j][q];
				else for (int q = 0; q < d; q++) PR[j][d + q] = R[j - n][q];
			}
		}

		private Map<Integer, BigFraction> phi(int p, int q) {
			int key = p * twoD + q;
			Map<Integer, BigFraction> f0 = phiCache.get(key);
			if (f0 != null) return f0;
			Map<Integer, BigFraction> form = new HashMap<>();
			for (int i = 0; i < N; i++) {
				BigFraction lv = PL[p][i];
				if (lv.getNumerator().signum() == 0) continue;
				for (int j = 0; j < N; j++) {
					BigFraction rv = PR[j][q];
					if (rv.getNumerator().signum() == 0) continue;
					form.put(i * N + j, lv.multiply(rv));
				}
			}
			phiCache.put(key, form);
			return form;
		}

		private int foo(int i) { return i; }

		private int bar(int i) { return i + d; }

		// sum of sc * phi(p,q) over (p,q,sc) triples
		private Map<Integer, BigFraction> sum(int[][] terms) {
			Map<Integer, BigFraction> res = new HashMap<>();
			for (int[] tr : terms) addInto(res, phi(tr[0], tr[1]), f(tr[2]));
			return res;
		}

		private Map<Integer, BigFraction> form4(BigFraction[] v, int i) {
			Map<Integer, BigFraction> res = new HashMap<>();
			addInto(res, phi(foo(i), foo(i)), v[0]);
			addInto(res, phi(bar(i), foo(i)), v[1]);
			addInto(res, phi(foo(i), bar(i)), v[2]);
			addInto(res, phi(bar(i), bar(i)), v[3]);
			return res;
		}

		private BigFraction[][][] litaLocal7() {
			BigFraction D = f(d);
			BigFraction[] e = { f(-1), ZERO, ZERO, ONE };
			BigFraction[] u = {
					D.multiply(f(8).subtract(D)).divide(D.subtract(f(6))),
					f(-2).multiply(D).divide(D.subtract(f(6))),
					D.multiply(D.subtract(f(2))).divide(D.subtract(f(6))),
					ZERO };
			BigFraction[] v = { ZERO, f(-1), ONE, ZERO };
			BigFraction[] w = { f(1).divide(f(2)), ZERO, f(1).divide(f(2)), ZERO };
			BigFraction A = D.divide(D.subtract(f(6)));
			BigFraction B = D.subtract(f(6)).negate().divide(D);
			BigFraction C = D.subtract(f(3)).divide(D.multiply(f(2)));
			BigFraction Ap = D.multiply(D.subtract(f(7))).divide(D.subtract(f(6)));
			BigFraction Bp = f(3).divide(D.multiply(f(2)));
			BigFraction Cp = D.subtract(f(6)).divide(D);
			BigFraction lam = D.pow(2).subtract(D.multiply(f(11))).add(f(27)).divide(D);
			BigFraction[] a = new BigFraction[4], b = new BigFraction[4], c = new BigFraction[4];
			BigFraction[] ap = new BigFraction[4], bp = new BigFraction[4], cp = new BigFraction[4];
			BigFraction[] lame = new BigFraction[4];
			for (int t = 0; t < 4; t++) {
				a[t] = u[t].add(A.multiply(e[t]));
				b[t] = v[t].add(B.multiply(e[t]));
				c[t] = w[t].add(C.multiply(e[t]));
				ap[t] = u[t].negate().add(Ap.multiply(e[t]));
				bp[t] = w[t].add(Bp.multiply(e[t]));
				cp[t] = v[t].add(Cp.multiply(e[t]));
				lame[t] = lam.multiply(e[t]);
			}
			BigFraction[][][] raw = {
					{ lame, e, e }, { a, b, c }, { b, c, a }, { c, a, b },
					{ ap, bp, cp }, { bp, cp, ap }, { cp, ap, bp } };
			BigFraction[] D0 = { ONE, f(-1), f(-1), ONE };
			BigFraction[] D2 = { ONE, ONE, f(-1), f(-1) };
			BigFraction[][][] out = new BigFraction[7][3][4];
			for (int r = 0; r < 7; r++) {
				for (int t = 0; t < 4; t++) {
					out[r][0][t] = raw[r][0][t];
					out[r][1][t] = D0[t].multiply(raw[r][1][t]);
					out[r][2][t] = D2[t].multiply(raw[r][2][t]);
				}
			}
			return out;
		}

		NonCubicBilinearAlgorithm build() {
			List<BigFraction> coeffs = new ArrayList<>();
			List<Map<Integer, BigFraction>> As = new ArrayList<>();
			List<Map<Integer, BigFraction>> Bs = new ArrayList<>();
			List<Map<Integer, BigFraction>> Xs = new ArrayList<>();

			// first family (unbarred / barred) and its symmetric half
			for (int mp = 0; mp < 2; mp++) {
				for (int k = 0; k <= n; k++)
					for (int j = 0; j < k; j++)
						for (int i = 0; i <= j; i++)
							addFirst(coeffs, As, Bs, Xs, mp, i, j, k);
				for (int i = 0; i <= n; i++)
					for (int j = 0; j <= i; j++)
						for (int k = 0; k < j; k++)
							addFirst(coeffs, As, Bs, Xs, mp, i, j, k);
			}

			// second family
			for (int i = 0; i <= n; i++)
				for (int j = 0; j <= n; j++)
					for (int k = 0; k <= n; k++) {
						if (i == j && j == k) continue;
						coeffs.add(ONE);
						As.add(sum(new int[][] { { bar(j), foo(k), 1 }, { foo(k), bar(i), 1 }, { foo(i), foo(j), -1 } }));
						Bs.add(sum(new int[][] { { foo(j), bar(k), 1 }, { foo(k), foo(i), 1 }, { bar(i), foo(j), 1 } }));
						Xs.add(sum(new int[][] { { foo(i), bar(j), 1 }, { foo(j), foo(k), 1 }, { bar(k), foo(i), -1 } }));
					}
			for (int i = 0; i <= n; i++)
				for (int j = 0; j <= n; j++)
					for (int k = 0; k <= n; k++) {
						if (i == j && j == k) continue;
						coeffs.add(ONE);
						As.add(sum(new int[][] { { foo(j), bar(k), 1 }, { bar(k), foo(i), 1 }, { bar(i), bar(j), -1 } }));
						Bs.add(sum(new int[][] { { bar(j), foo(k), 1 }, { bar(k), bar(i), 1 }, { foo(i), bar(j), 1 } }));
						Xs.add(sum(new int[][] { { bar(i), foo(j), 1 }, { bar(j), bar(k), 1 }, { foo(k), bar(i), -1 } }));
					}

			// diagonal local rank-7
			BigFraction[][][] cf7 = litaLocal7();
			for (int i = 0; i <= n; i++)
				for (int r = 0; r < 7; r++) {
					coeffs.add(ONE);
					As.add(form4(cf7[r][0], i));
					Bs.add(form4(cf7[r][1], i));
					Xs.add(form4(cf7[r][2], i));
				}

			// off-diagonal Winograd correction, scaled by -d
			BigFraction md = f(-d);
			for (int i = 0; i <= n; i++)
				for (int j = 0; j <= n; j++) {
					if (i == j) continue;
					addProd(coeffs, As, Bs, Xs, md,
							new int[][] { { foo(i), foo(j), 1 }, { bar(i), foo(j), -1 } },
							new int[][] { { foo(i), foo(j), 1 }, { bar(i), foo(j), -1 } },
							new int[][] { { foo(i), foo(j), 1 }, { bar(i), foo(j), 1 } });
					addProd(coeffs, As, Bs, Xs, md,
							new int[][] { { foo(i), foo(j), 1 }, { foo(i), bar(j), 1 } },
							new int[][] { { foo(i), foo(j), 1 }, { foo(i), bar(j), 1 } },
							new int[][] { { foo(i), foo(j), 1 }, { foo(i), bar(j), -1 } });
					addProd(coeffs, As, Bs, Xs, md,
							new int[][] { { bar(i), bar(j), 1 } },
							new int[][] { { bar(i), bar(j), 1 } },
							new int[][] { { bar(i), bar(j), 1 } });
					addProd(coeffs, As, Bs, Xs, md,
							new int[][] { { bar(i), foo(j), -1 } },
							new int[][] { { bar(i), foo(j), 1 }, { bar(i), bar(j), 1 }, { foo(i), foo(j), -1 }, { foo(i), bar(j), -1 } },
							new int[][] { { foo(i), bar(j), 1 } });
					addProd(coeffs, As, Bs, Xs, md,
							new int[][] { { bar(i), foo(j), 1 }, { bar(i), bar(j), 1 }, { foo(i), foo(j), -1 }, { foo(i), bar(j), -1 } },
							new int[][] { { foo(i), bar(j), -1 } },
							new int[][] { { bar(i), foo(j), 1 } });
					addProd(coeffs, As, Bs, Xs, md,
							new int[][] { { foo(i), bar(j), 1 } },
							new int[][] { { bar(i), foo(j), 1 } },
							new int[][] { { foo(i), bar(j), 1 }, { bar(i), bar(j), 1 }, { foo(i), foo(j), -1 }, { bar(i), foo(j), -1 } });
					addProd(coeffs, As, Bs, Xs, md,
							new int[][] { { bar(i), foo(j), 1 }, { foo(i), foo(j), -1 }, { foo(i), bar(j), -1 } },
							new int[][] { { foo(i), foo(j), 1 }, { foo(i), bar(j), 1 }, { bar(i), foo(j), -1 } },
							new int[][] { { foo(i), foo(j), 1 }, { foo(i), bar(j), -1 }, { bar(i), foo(j), 1 } });
				}

			// plain merge of bit-identical triples
			Map<String, Object[]> acc = new HashMap<>();
			for (int t = 0; t < coeffs.size(); t++) {
				Map<Integer, BigFraction> A = As.get(t), B = Bs.get(t), X = Xs.get(t);
				if (A.isEmpty() || B.isEmpty() || X.isEmpty()) continue;
				String key = canon(A) + "|" + canon(B) + "|" + canon(X);
				Object[] slot = acc.get(key);
				if (slot == null) {
					slot = new Object[] { A, B, X, coeffs.get(t) };
					acc.put(key, slot);
				} else {
					slot[3] = ((BigFraction) slot[3]).add(coeffs.get(t));
				}
			}
			List<Row> rows = new ArrayList<>();
			for (Object[] slot : acc.values()) {
				BigFraction c = (BigFraction) slot[3];
				if (c.getNumerator().signum() == 0) continue;
				@SuppressWarnings("unchecked")
				Map<Integer, BigFraction> A = (Map<Integer, BigFraction>) slot[0];
				@SuppressWarnings("unchecked")
				Map<Integer, BigFraction> B = (Map<Integer, BigFraction>) slot[1];
				@SuppressWarnings("unchecked")
				Map<Integer, BigFraction> X = (Map<Integer, BigFraction>) slot[2];
				Map<Integer, BigFraction> Ac = new HashMap<>();
				addInto(Ac, A, c);
				rows.add(new Row(Ac, new HashMap<>(B), new HashMap<>(X)));
			}
			return emit(N, rows);
		}

		private void addFirst(List<BigFraction> coeffs, List<Map<Integer, BigFraction>> As,
				List<Map<Integer, BigFraction>> Bs, List<Map<Integer, BigFraction>> Xs,
				int mp, int i, int j, int k) {
			int fi = (mp == 0) ? foo(i) : bar(i);
			int fj = (mp == 0) ? foo(j) : bar(j);
			int fk = (mp == 0) ? foo(k) : bar(k);
			coeffs.add(ONE);
			As.add(sum(new int[][] { { fi, fj, 1 }, { fj, fk, 1 }, { fk, fi, 1 } }));
			Bs.add(sum(new int[][] { { fj, fk, 1 }, { fk, fi, 1 }, { fi, fj, 1 } }));
			Xs.add(sum(new int[][] { { fk, fi, 1 }, { fi, fj, 1 }, { fj, fk, 1 } }));
		}

		private void addProd(List<BigFraction> coeffs, List<Map<Integer, BigFraction>> As,
				List<Map<Integer, BigFraction>> Bs, List<Map<Integer, BigFraction>> Xs,
				BigFraction coeff, int[][] aT, int[][] bT, int[][] xT) {
			coeffs.add(coeff);
			As.add(sum(aT));
			Bs.add(sum(bT));
			Xs.add(sum(xT));
		}
	}
}
