package eu.solven.matmul.recombination;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Exact rational subspace algebra over {@code ℚ^dim} for the recombination-multiset
 * frontier ({@link RecombinationMultisetOrbit}). Subspaces are carried as an integer
 * basis (a {@code List<int[]>}); all arithmetic is exact ({@link Frac}, gcd-normalised),
 * so containment / intersection / sum are decided without floating error.
 *
 * <p>Why it exists: the structural frontier sweep needs candidate flag directions that
 * lie inside every relevant subspace of the footprint <b>arrangement</b> — not just the
 * footprints themselves but their pairwise (and iterated) <b>intersections</b> (for the
 * ascending-flag / {@code uMin} side) and <b>sums</b> (for the tail-flag / {@code wMin}
 * side). A dense integer cube accidentally supplies these but explodes at {@code dim=4};
 * this class supplies them exactly and cheaply (one generic vector per lattice member).
 */
final class SubspaceArrangement {
	private SubspaceArrangement() {}

	/** Defensive cap on a lattice closure (intersections/sums); closures here are tiny in practice. */
	static final int LATTICE_CAP = 4096;

	// ---- exact rationals -----------------------------------------------------

	private static final class Frac {
		final long num, den; // den > 0, gcd(|num|,den)=1
		Frac(long n, long d) {
			if (d == 0) throw new ArithmeticException("zero denominator");
			if (d < 0) { n = -n; d = -d; }
			long g = gcd(Math.abs(n), d);
			if (g == 0) g = 1;
			num = n / g; den = d / g;
		}
		static final Frac ZERO = new Frac(0, 1);
		boolean isZero() { return num == 0; }
		Frac sub(Frac o) { return new Frac(num * o.den - o.num * den, den * o.den); }
		Frac mul(Frac o) { return new Frac(num * o.num, den * o.den); }
		Frac div(Frac o) { return new Frac(num * o.den, den * o.num); }
	}

	private static long gcd(long a, long b) { a = Math.abs(a); b = Math.abs(b); while (b != 0) { long t = a % b; a = b; b = t; } return a; }

	// ---- core: reduced row echelon over ℚ ------------------------------------

	/** RREF of the integer rows (as ℚ), in place on a fresh Frac matrix; returns pivot columns in order. */
	private static List<Integer> rref(Frac[][] M, int cols) {
		List<Integer> pivots = new ArrayList<>();
		int row = 0;
		for (int c = 0; c < cols && row < M.length; c++) {
			int piv = -1;
			for (int i = row; i < M.length; i++) if (!M[i][c].isZero()) { piv = i; break; }
			if (piv < 0) continue;
			Frac[] tmp = M[piv]; M[piv] = M[row]; M[row] = tmp;
			Frac inv = new Frac(1, 1).div(M[row][c]);
			for (int j = 0; j < cols; j++) M[row][j] = M[row][j].mul(inv);
			for (int i = 0; i < M.length; i++) {
				if (i == row || M[i][c].isZero()) continue;
				Frac f = M[i][c];
				for (int j = 0; j < cols; j++) M[i][j] = M[i][j].sub(f.mul(M[row][j]));
			}
			pivots.add(c);
			row++;
		}
		return pivots;
	}

	private static Frac[][] toFrac(List<int[]> rows, int dim) {
		Frac[][] M = new Frac[Math.max(1, rows.size())][dim];
		if (rows.isEmpty()) { for (int j = 0; j < dim; j++) M[0][j] = Frac.ZERO; return new Frac[0][dim]; }
		for (int i = 0; i < rows.size(); i++) for (int j = 0; j < dim; j++) M[i][j] = new Frac(rows.get(i)[j], 1);
		return M;
	}

	/** Clear denominators of a rational vector and gcd-reduce to a primitive integer vector. */
	private static int[] primitive(Frac[] v) {
		long lcm = 1;
		for (Frac f : v) if (!f.isZero()) lcm = lcm / gcd(lcm, f.den) * f.den;
		int[] out = new int[v.length];
		long g = 0;
		for (int i = 0; i < v.length; i++) { out[i] = (int) (v[i].num * (lcm / v[i].den)); g = gcd(g, Math.abs(out[i])); }
		if (g > 1) for (int i = 0; i < out.length; i++) out[i] /= g;
		return out;
	}

	// ---- public subspace ops -------------------------------------------------

	/** Integer basis of the null space {@code {x ∈ ℚ^dim : row·x = 0 ∀ row}}. */
	static List<int[]> nullspace(List<int[]> rows, int dim) {
		Frac[][] M = toFrac(rows, dim);
		List<Integer> pivots = M.length == 0 ? new ArrayList<>() : rref(M, dim);
		boolean[] isPivot = new boolean[dim];
		for (int c : pivots) isPivot[c] = true;
		List<int[]> basis = new ArrayList<>();
		for (int free = 0; free < dim; free++) {
			if (isPivot[free]) continue;
			Frac[] x = new Frac[dim];
			for (int j = 0; j < dim; j++) x[j] = Frac.ZERO;
			x[free] = new Frac(1, 1);
			for (int pi = 0; pi < pivots.size(); pi++) {
				int pc = pivots.get(pi);
				x[pc] = Frac.ZERO.sub(M[pi][free]); // pivot row pi has 1 at pc; x[pc] = -coeff(free)
			}
			basis.add(primitive(x));
		}
		return basis;
	}

	/** Reduce a generating set to an independent integer basis of its span. */
	static List<int[]> basisOf(List<int[]> gens, int dim) {
		Frac[][] M = toFrac(gens, dim);
		if (M.length == 0) return new ArrayList<>();
		List<Integer> pivots = rref(M, dim);
		List<int[]> basis = new ArrayList<>();
		for (int i = 0; i < pivots.size(); i++) basis.add(primitive(M[i]));
		return basis;
	}

	/** Orthogonal complement (annihilator) {@code Aᗮ}: null space of A's basis-as-rows. */
	static List<int[]> annihilator(List<int[]> basis, int dim) { return nullspace(basis, dim); }

	/** Intersection {@code A ∩ B = (Aᗮ + Bᗮ)ᗮ}, returned as an integer basis. */
	static List<int[]> intersect(List<int[]> A, List<int[]> B, int dim) {
		List<int[]> ann = new ArrayList<>(annihilator(A, dim));
		ann.addAll(annihilator(B, dim));
		return nullspace(ann, dim);
	}

	/** Sum {@code A + B} (span of the union), returned as an independent integer basis. */
	static List<int[]> sum(List<int[]> A, List<int[]> B, int dim) {
		List<int[]> all = new ArrayList<>(A); all.addAll(B);
		return basisOf(all, dim);
	}

	/** Canonical key of a subspace (RREF of its basis), for dedup in a lattice closure. */
	static String key(List<int[]> basis, int dim) {
		List<int[]> red = basisOf(basis, dim);
		StringBuilder sb = new StringBuilder();
		for (int[] v : red) sb.append(java.util.Arrays.toString(v)).append(';');
		return sb.toString();
	}

	/**
	 * Lattice closure of {@code generators} under {@code ∩} (op="meet") or {@code +} (op="join"),
	 * each member returned as a basis. Deduped by {@link #key}; capped at {@link #LATTICE_CAP}.
	 */
	static List<List<int[]>> closure(List<List<int[]>> generators, int dim, boolean meet) {
		Map<String, List<int[]>> seen = new LinkedHashMap<>();
		List<List<int[]>> work = new ArrayList<>();
		for (List<int[]> g : generators) {
			List<int[]> b = basisOf(g, dim);
			if (b.isEmpty()) continue;
			String k = key(b, dim);
			if (seen.putIfAbsent(k, b) == null) work.add(b);
		}
		for (int i = 0; i < work.size() && seen.size() < LATTICE_CAP; i++) {
			for (int j = 0; j < i; j++) {
				List<int[]> c = meet ? intersect(work.get(i), work.get(j), dim) : sum(work.get(i), work.get(j), dim);
				if (c.isEmpty()) continue;
				String k = key(c, dim);
				if (seen.putIfAbsent(k, c) == null) { work.add(c); if (seen.size() >= LATTICE_CAP) break; }
			}
		}
		return new ArrayList<>(seen.values());
	}

	/**
	 * Lattice closure of {@code generators} under BOTH {@code ∩} and {@code +} — the full joint
	 * arrangement. Needed because the index is {@code max(uMin,wMin)}: a cell can mix an
	 * ascending-flag (A-side, intersection) constraint on some products with a tail-flag (S-side,
	 * sum) constraint on others, governed by subspaces like {@code (∩A_K) ∩ (ΣS_J)} that neither
	 * single-op closure produces. Deduped by {@link #key}; capped at {@link #LATTICE_CAP}.
	 */
	static List<List<int[]>> closureBoth(List<List<int[]>> generators, int dim) {
		Map<String, List<int[]>> seen = new LinkedHashMap<>();
		List<List<int[]>> work = new ArrayList<>();
		for (List<int[]> g : generators) {
			List<int[]> b = basisOf(g, dim);
			if (b.isEmpty()) continue;
			if (seen.putIfAbsent(key(b, dim), b) == null) work.add(b);
		}
		for (int i = 0; i < work.size() && seen.size() < LATTICE_CAP; i++) {
			for (int j = 0; j < i && seen.size() < LATTICE_CAP; j++) {
				for (List<int[]> c : List.of(intersect(work.get(i), work.get(j), dim), sum(work.get(i), work.get(j), dim))) {
					if (c.isEmpty()) continue;
					if (seen.putIfAbsent(key(c, dim), c) == null) work.add(c);
				}
			}
		}
		return new ArrayList<>(seen.values());
	}

	/**
	 * A generic integer vector inside the span of {@code basis}: combination with strictly
	 * increasing coprime-ish coefficients {@code (1, t, t², …)}, {@code t} chosen above the
	 * basis entries so the combination avoids the basis's own proper sub-arrangement.
	 */
	static int[] genericIn(List<int[]> basis, int dim, int t) {
		Frac[] acc = new Frac[dim];
		for (int j = 0; j < dim; j++) acc[j] = Frac.ZERO;
		long c = 1;
		for (int[] b : basis) {
			for (int j = 0; j < dim; j++) acc[j] = acc[j].sub(new Frac(-c * b[j], 1)); // acc += c*b
			c *= t;
		}
		return primitive(acc);
	}
}
