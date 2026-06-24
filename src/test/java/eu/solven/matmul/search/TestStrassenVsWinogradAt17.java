package eu.solven.matmul.search;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Side-by-side comparison: for target ⟨17,17,17⟩ at allocation (9,8)³,
 * what 7 sub-products does each of (Strassen, Winograd-cousin) produce,
 * and at what catalog rank? This is the "why Winograd beats Strassen
 * here" view at the shape level.
 */
class TestStrassenVsWinogradAt17 {

	@Test
	void print_side_by_side() throws Exception {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		CitedBound sota = new CitedBound(lookup);

		NonCubicBilinearAlgorithm strassen = SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint(
				"src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		NonCubicBilinearAlgorithm winogradCousin = SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint(
				"src/main/resources/schemes/known/section2/solven_winograd_cousin_axflip1-2x2x2_m7_a24.json"));

		int[] alloc = { 9, 8 };
		Recombination.Result rsStr = Recombination.recombineWithAllocation(strassen, sota, alloc, alloc, alloc);
		Recombination.Result rsWin = Recombination.recombineWithAllocation(winogradCousin, sota, alloc, alloc, alloc);

		System.out.println();
		System.out.println("==============================================================================");
		System.out.println("⟨17,17,17⟩ via 2x2 Strassen-style recombination, allocation (9,8)³");
		System.out.println("==============================================================================");
		System.out.println();
		System.out.printf("%-8s | %-16s | %-7s || %-16s | %-7s%n",
				"Product", "Strassen shape", "rank", "Winograd shape", "rank");
		System.out.println("---------|------------------|---------||------------------|--------");
		long totalStr = 0, totalWin = 0;
		for (int k = 0; k < 7; k++) {
			int[] s = rsStr.smallMatrixSizes[k];
			int[] w = rsWin.smallMatrixSizes[k];
			long rStr = sota.getRank(s[0], s[1], s[2]);
			long rWin = sota.getRank(w[0], w[1], w[2]);
			totalStr += rStr;
			totalWin += rWin;
			System.out.printf("M%d       | ⟨%d,%d,%d⟩            | %5d   || ⟨%d,%d,%d⟩            | %5d%n",
					k + 1, s[0], s[1], s[2], rStr,
					w[0], w[1], w[2], rWin);
		}
		System.out.println("---------|------------------|---------||------------------|--------");
		System.out.printf("%-8s | %-16s | %5d   || %-16s | %5d%n",
				"TOTAL", "", totalStr, "", totalWin);
		System.out.println();

		System.out.println("Shape distribution (count by shape multiset):");
		System.out.println();
		System.out.println("Strassen (9,8)³:");
		printShapeCounts(rsStr);
		System.out.println();
		System.out.println("Winograd-cousin (axflip mask=1) (9,8)³:");
		printShapeCounts(rsWin);
		System.out.println();

		System.out.println("Catalog ranks used (best-known over Q for our sota):");
		int[][] shapes = { {9,9,9}, {8,9,9}, {9,8,9}, {9,9,8}, {8,8,9}, {9,8,8}, {8,9,8}, {8,8,8} };
		for (int[] s : shapes) {
			System.out.printf("  R(⟨%d,%d,%d⟩) = %d%n",
					s[0], s[1], s[2], sota.getRank(s[0], s[1], s[2]));
		}
		System.out.println();
		System.out.println("==============================================================================");
	}

	private static void printShapeCounts(Recombination.Result r) {
		java.util.TreeMap<String, Integer> counts = new java.util.TreeMap<>();
		for (int k = 0; k < r.smallMatrixSizes.length; k++) {
			int[] s = r.smallMatrixSizes[k];
			String key = "⟨" + s[0] + "," + s[1] + "," + s[2] + "⟩";
			counts.merge(key, 1, Integer::sum);
		}
		for (var e : counts.entrySet()) {
			System.out.printf("  %d × %s%n", e.getValue(), e.getKey());
		}
	}
}
