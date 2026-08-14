package eu.solven.matmul.docs.explore;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.papers.khoruzhii2026.LitaTaConstruction;
import eu.solven.matmul.papers.schwartzzwecher2025.TaNew25Construction;

/**
 * Kin/bud graph measurement for TA schemes — quantifies how much "local
 * unification" headroom the SZ/LITA constructions leave, and (for SZ, where we
 * know the family emission order) whether the frozen {@code 15/4·N²} off-diagonal
 * cancellation family (c) carries any shared-factor structure.
 *
 * <p><b>Two-axis kin</b>: products proportional on two of {U,V,W}. Such a pair is
 * TRIVIALLY mergeable ({@code u⊗v⊗w + u⊗v⊗w' = u⊗v⊗(w+w')}, a free −1), so a
 * locally-reduced scheme should have ~0. Any found = a merge left on the table.</p>
 *
 * <p><b>One-axis buds</b>: products sharing (proportionally) a single factor. NOT
 * trivially mergeable — the substrate for serendipitous products and higher-order
 * re-aggregation. Reported as a size histogram, per axis and (for SZ) per family.</p>
 *
 * <p>Args: {@code sz N} | {@code lita N} | {@code <scheme.json path>}.</p>
 */
public final class ProbeTaKinGraph {
	private ProbeTaKinGraph() {}

	public static void main(String[] args) throws Exception {
		String kind = args[0];

		// SYMBOLIC (transformed-space) mode: the space where SZ/LITA unification lives.
		if (kind.equals("szsym")) {
			int n = Integer.parseInt(args[1]);
			double[][][] sym = TaNew25Construction.buildSymbolicForms(n);   // [product][axis][dim²]
			int r = sym.length;
			int[] fb = szFamilyBoundaries(n, r);
			int rows = sym[0][0].length;
			double[][] A = axisMatrix(sym, 0, rows), B = axisMatrix(sym, 1, rows), C = axisMatrix(sym, 2, rows);
			String[] cu = canonAll(A, r), cv = canonAll(B, r), cw = canonAll(C, r);
			System.out.printf("%n=== SZ ⟨%d³⟩ TRANSFORMED-space (A*/B*/C* symbolic) : r=%d ===%n", n, r);
			System.out.println("(kin here = genuine unification opportunity, unlike the pulled-back scheme)");
			report(cu, cv, cw, r, fb);
			return;
		}

		NonCubicBilinearAlgorithm alg;
		int[] famBound = null;    // SZ only: {aEnd, bEnd, cEnd} (dEnd = r)
		String label;
		if (kind.equals("sz")) {
			int n = Integer.parseInt(args[1]);
			alg = TaNew25Construction.build(n);
			famBound = szFamilyBoundaries(n, alg.r);
			label = "SZ ⟨" + n + "³⟩";
		} else if (kind.equals("lita")) {
			int n = Integer.parseInt(args[1]);
			alg = LitaTaConstruction.build(n);
			label = "LITA ⟨" + n + "³⟩";
		} else {
			alg = SchemeIO.read(new File(kind));
			label = new File(kind).getName();
		}
		int r = alg.r;
		double[][] U = alg.denseU(), V = alg.denseV(), W = alg.denseW();
		String[] cu = canonAll(U, r), cv = canonAll(V, r), cw = canonAll(W, r);
		System.out.printf("%n=== %s : rank r=%d ===%n", label, r);
		report(cu, cv, cw, r, famBound);
	}

	/** Shared reporting: distinct dirs, two-axis kin, bud histograms, and (SZ) per-family. */
	private static void report(String[] cu, String[] cv, String[] cw, int r, int[] famBound) {
		// per-axis distinct proportionality classes (r − distinct = shared-factor redundancy)
		System.out.printf("distinct factor-directions: U=%d  V=%d  W=%d  (of r=%d)%n",
				distinct(cu), distinct(cv), distinct(cw), r);

		// two-axis kin (proportional on two axes) = directly unitable pairs (−1 each).
		System.out.printf("TWO-AXIS KIN (unitable, −1 each): U∧V=%d  U∧W=%d  V∧W=%d  products%n",
				kinExcess(cu, cv), kinExcess(cu, cw), kinExcess(cv, cw));

		// one-axis bud histograms
		System.out.println("ONE-AXIS BUD size histogram (groups of ≥2 sharing a factor):");
		System.out.println("  U: " + budHisto(cu));
		System.out.println("  V: " + budHisto(cv));
		System.out.println("  W: " + budHisto(cw));

		if (famBound != null) {
			int[] fam = new int[r];
			for (int l = 0; l < r; l++) {
				fam[l] = l < famBound[0] ? 0 : l < famBound[1] ? 1 : l < famBound[2] ? 2 : 3;
			}
			String[] fname = { "(a)Ŝ-symm", "(b)Ṡ-barred", "(c)off-diag[15/4·N²]", "(d)R(i)-diag" };
			int[] cnt = new int[4];
			for (int f : fam) cnt[f]++;
			System.out.println("SZ family sizes: " + fname[0] + "=" + cnt[0] + "  " + fname[1] + "=" + cnt[1]
					+ "  " + fname[2] + "=" + cnt[2] + "  " + fname[3] + "=" + cnt[3]);

			// the crux: restrict to family (c) and measure its INTERNAL bud/kin structure.
			System.out.println("--- family (c) off-diagonal, in isolation (the frozen 15/4·N² term) ---");
			String[] cuC = restrict(cu, fam, 2), cvC = restrict(cv, fam, 2), cwC = restrict(cw, fam, 2);
			int rc = cuC.length;
			System.out.printf("  products=%d  distinct dirs U=%d V=%d W=%d%n",
					rc, distinct(cuC), distinct(cvC), distinct(cwC));
			System.out.printf("  two-axis kin WITHIN (c): U∧V=%d U∧W=%d V∧W=%d%n",
					kinExcess(cuC, cvC), kinExcess(cuC, cwC), kinExcess(cvC, cwC));
			System.out.println("  one-axis buds WITHIN (c): U=" + budHisto(cuC)
					+ "  V=" + budHisto(cvC) + "  W=" + budHisto(cwC));

			// cross-family: do (c) products share a factor-direction with ANY other family?
			System.out.println("--- do family-(c) directions also appear in families a/b/d? ---");
			crossFamilyShare(cu, fam, "U");
			crossFamilyShare(cv, fam, "V");
			crossFamilyShare(cw, fam, "W");
		}
	}

	/** Reshape sym[product][axis][row] → [row][product] for the chosen axis (canonAll input). */
	private static double[][] axisMatrix(double[][][] sym, int axis, int rows) {
		int r = sym.length;
		double[][] M = new double[rows][r];
		for (int l = 0; l < r; l++) for (int t = 0; t < rows; t++) M[t][l] = sym[l][axis][t];
		return M;
	}

	/** SZ family boundaries {aEnd,bEnd,cEnd} from d=n/2+1; asserts they sum to r. */
	private static int[] szFamilyBoundaries(int n, int r) {
		int d = n / 2 + 1;
		int countA = 0;
		for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) for (int k = 0; k < d; k++) {
			if ((i <= j && j < k) || (k < j && j <= i)) countA++;
		}
		int a = 2 * countA;
		int b = 2 * d * d * d - d;
		int c = 7 * (d * d - d);
		int dd = 7 * d;
		if (a + b + c + dd != r) {
			throw new IllegalStateException("family counts " + a + "+" + b + "+" + c + "+" + dd
					+ "=" + (a + b + c + dd) + " != r=" + r);
		}
		return new int[] { a, a + b, a + b + c };
	}

	/** Canonicalise each product's factor column into a proportionality-class key. */
	private static String[] canonAll(double[][] F, int r) {
		int rows = F.length;
		String[] out = new String[r];
		for (int l = 0; l < r; l++) {
			int i0 = -1;
			for (int i = 0; i < rows; i++) if (Math.abs(F[i][l]) > 1e-9) { i0 = i; break; }
			if (i0 < 0) { out[l] = "0"; continue; }
			double s = F[i0][l];
			StringBuilder sb = new StringBuilder(i0 + ":");
			for (int i = i0; i < rows; i++) {
				double q = F[i][l] / s;
				long rq = Math.round(q * 1e7);
				if (rq == 0) rq = 0;   // kill -0
				sb.append(rq).append(',');
			}
			out[l] = sb.toString();
		}
		return out;
	}

	private static int distinct(String[] keys) {
		Map<String, Integer> m = new HashMap<>();
		for (String k : keys) m.merge(k, 1, Integer::sum);
		return m.size();
	}

	/** Σ(size−1) over classes of the (a,b) pair key with size≥2 = trivially-mergeable count. */
	private static int kinExcess(String[] a, String[] b) {
		Map<String, Integer> m = new HashMap<>();
		for (int l = 0; l < a.length; l++) m.merge(a[l] + "|" + b[l], 1, Integer::sum);
		int excess = 0;
		for (int c : m.values()) if (c >= 2) excess += c - 1;
		return excess;
	}

	private static String budHisto(String[] keys) {
		Map<String, Integer> m = new HashMap<>();
		for (String k : keys) m.merge(k, 1, Integer::sum);
		TreeMap<Integer, Integer> h = new TreeMap<>();
		for (int c : m.values()) if (c >= 2) h.merge(c, 1, Integer::sum);
		if (h.isEmpty()) return "{none}";
		StringBuilder sb = new StringBuilder("{");
		int inBuds = 0;
		for (var e : h.entrySet()) {
			sb.append("size").append(e.getKey()).append("×").append(e.getValue()).append(' ');
			inBuds += e.getKey() * e.getValue();
		}
		sb.append("| ").append(inBuds).append(" products in buds}");
		return sb.toString();
	}

	private static String[] restrict(String[] keys, int[] fam, int f) {
		List<String> out = new ArrayList<>();
		for (int l = 0; l < keys.length; l++) if (fam[l] == f) out.add(keys[l]);
		return out.toArray(new String[0]);
	}

	/** For each family-(c) direction, report how many are shared with a/b/d. */
	private static void crossFamilyShare(String[] keys, int[] fam, String axis) {
		// map direction -> set of families that use it
		Map<String, boolean[]> use = new HashMap<>();
		for (int l = 0; l < keys.length; l++) {
			use.computeIfAbsent(keys[l], x -> new boolean[4])[fam[l]] = true;
		}
		int cDirs = 0, cSharedWithOther = 0;
		for (var e : use.entrySet()) {
			if (!e.getValue()[2]) continue;   // not a (c) direction
			cDirs++;
			if (e.getValue()[0] || e.getValue()[1] || e.getValue()[3]) cSharedWithOther++;
		}
		System.out.printf("  %s: %d distinct (c)-directions, %d also used by a/b/d%n",
				axis, cDirs, cSharedWithOther);
	}
}
