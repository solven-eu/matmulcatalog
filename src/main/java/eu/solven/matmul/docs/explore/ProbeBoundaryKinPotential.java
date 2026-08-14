package eu.solven.matmul.docs.explore;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import eu.solven.matmul.papers.schwartzzwecher2025.TaNew25Construction;

/**
 * Tests the ceiling of the LITA/SZ-style kin-unification lever for the {@code 2·N²}
 * boundary. A product unites with another (a −1 saving) only if they are PROPORTIONAL
 * on two of the three axes — which REQUIRES identical support (same nonzero
 * transformed-variable indices) on those two axes. If a boundary product shares no
 * 2-axis support with ANY other product, it can never be united — no coefficient
 * tuning (e.g. of the correction ⟨2,2,2;7⟩ orbit) can help, because the variable
 * ENTRIES differ, not just their coefficients.
 *
 * <p>Also reports whether boundary products are coefficient-RIGID (all nonzero coeffs
 * ±1 → nothing to tune on the boundary side). Together these decide whether SZ's
 * diagonal mechanism (tunable + co-located correction) has ANY analog at the 2-D
 * boundary. Args: {@code n0}.</p>
 */
public final class ProbeBoundaryKinPotential {
	private ProbeBoundaryKinPotential() {}

	private static final String[] CLASS = { "generic", "boundary", "correction", "diagonal" };

	public static void main(String[] args) {
		int n0 = Integer.parseInt(args[0]);
		TaNew25Construction.SymbolicTagged st = TaNew25Construction.buildSymbolicTagged(n0);
		double[][][] F = st.forms();
		int[] cls = st.productClass();
		int r = F.length;

		int[] hist = new int[4];
		for (int c : cls) hist[c]++;
		System.out.printf("%n=== SZ ⟨%d³⟩ boundary kin-potential : r=%d ===%n", n0, r);
		System.out.printf("classes: generic=%d boundary=%d correction=%d diagonal=%d%n",
				hist[0], hist[1], hist[2], hist[3]);

		// (1) boundary coefficient rigidity: are all nonzero coeffs ±1 (nothing to tune)?
		int rigid = 0, bcount = 0;
		for (int l = 0; l < r; l++) {
			if (cls[l] != 1) continue;
			bcount++;
			if (allPmOne(F[l][0]) && allPmOne(F[l][1]) && allPmOne(F[l][2])) rigid++;
		}
		System.out.printf("boundary products with ALL ±1 coefficients (untunable): %d / %d%n", rigid, bcount);

		// (2) 2-axis support match: build maps over each axis-pair.
		// pair 0 = (A*,B*), 1 = (A*,C*), 2 = (B*,C*)
		int[][] PAIR = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
		List<Map<String, List<Integer>>> maps = new ArrayList<>();
		for (int p = 0; p < 3; p++) {
			Map<String, List<Integer>> m = new HashMap<>();
			for (int l = 0; l < r; l++) {
				String key = supp(F[l][PAIR[p][0]]) + "#" + supp(F[l][PAIR[p][1]]);
				m.computeIfAbsent(key, x -> new ArrayList<>()).add(l);
			}
			maps.add(m);
		}

		// for each boundary product, does any OTHER product share a 2-axis support key?
		int isolated = 0;
		long[] partnerByClass = new long[4];
		int withPartner = 0;
		for (int l = 0; l < r; l++) {
			if (cls[l] != 1) continue;
			boolean anyPartner = false;
			int[] pc = new int[4];
			for (int p = 0; p < 3; p++) {
				String key = supp(F[l][PAIR[p][0]]) + "#" + supp(F[l][PAIR[p][1]]);
				for (int other : maps.get(p).get(key)) {
					if (other == l) continue;
					anyPartner = true;
					pc[cls[other]]++;
				}
			}
			if (!anyPartner) isolated++;
			else {
				withPartner++;
				for (int c = 0; c < 4; c++) partnerByClass[c] += pc[c];
			}
		}
		System.out.printf("boundary products that share a 2-axis support with SOME other product: %d / %d%n",
				withPartner, bcount);
		System.out.printf("boundary products fully SUPPORT-ISOLATED (no unification possible, ever): %d / %d%n",
				isolated, bcount);
		if (withPartner > 0) {
			System.out.print("  partners by class (across axis-pairs): ");
			for (int c = 0; c < 4; c++) System.out.printf("%s=%d ", CLASS[c], partnerByClass[c]);
			System.out.println();
		}
		System.out.println("--------");
		if (isolated == bcount) {
			System.out.println("VERDICT: every boundary product is support-isolated ⟹ the kin-unification "
					+ "lever (SZ/LITA style) is STRUCTURALLY CLOSED for the 2·N² term — no φ / correction-orbit "
					+ "choice can create a boundary unification (the variable entries are disjoint).");
		} else {
			System.out.printf("VERDICT: %d boundary products have a support-compatible partner ⟹ a "
					+ "coefficient search over those pairs COULD engineer unifications — worth a targeted try.%n",
					withPartner);
		}
	}

	private static boolean allPmOne(double[] v) {
		for (double x : v) if (Math.abs(x) > 1e-9 && Math.abs(Math.abs(x) - 1.0) > 1e-9) return false;
		return true;
	}

	private static String supp(double[] v) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < v.length; i++) if (Math.abs(v[i]) > 1e-9) sb.append(i).append(',');
		return sb.toString();
	}
}
