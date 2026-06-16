package eu.solven.matmul.io;

import java.io.File;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.condition.EnabledIf;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Diagnostic / investigation harness — does NOT assert anything, just
 * reports structural statistics about the ⟨17,17,17⟩ schemes (LRP-parsed
 * rank 2931 and raw-parsed rank 2945) to the test log.
 *
 * <p>Specifically:</p>
 * <ol>
 *   <li>How many product columns share the same {@code (U[:,k], V[:,k])}
 *       pair? (zero would mean kin-row reduction has already merged
 *       all shareable products; a small positive count would corroborate
 *       the "3 free multiplications = 3 duplicates" hypothesis)</li>
 *   <li>What is the joint distribution of {@code (rowExtent, colExtent)}
 *       for the A-side support of each product? — to cross-check FMM-Lille's
 *       advertised polynomial breakdown
 *       {@code 95×⟨4,4,4⟩ + 336×⟨3,3,3⟩ + ...} summing to 2931.</li>
 * </ol>
 */
public class TestMapleLRPInvestigation {

	private static final File LRP = new File("references/fmm-lille/17x17x17/17x17x17_LRP.mpl");
	// STRICT: this investigation targets the specific m2945 raw import (purged
	// 2026-06). Lax byHint fell back to an arbitrary same-shape file, making
	// rawJsonAvailable() true with the WRONG scheme. Strict returns null when the
	// import is absent, and the @EnabledIf gate disables the test honestly.
	private static final File RAW_JSON = eu.solven.matmul.catalog.SchemeResolver.byHintStrict(
			"src/main/resources/schemes/known/section17/fmm_lille_2025-17x17x17_m2945_a68812.json");

	static boolean lrpAvailable() { return LRP.isFile(); }
	static boolean rawJsonAvailable() { return RAW_JSON != null && RAW_JSON.isFile(); }

	@Test
	@EnabledIf("lrpAvailable")
	public void report_duplicate_UV_pairs_in_LRP() throws Exception {
		NonCubicBilinearAlgorithm alg = MapleLRPParser.parse(LRP, 17, 17, 17);
		int dups = MapleLRPParser.countDuplicateUVPairs(alg);
		System.out.println("\n=== LRP scheme (r=" + alg.r + ") duplicate (U,V) pair count: "
				+ dups + " ===");
	}

	@Test
	@EnabledIf("rawJsonAvailable")
	public void report_duplicate_UV_pairs_in_raw_2945() throws Exception {
		NonCubicBilinearAlgorithm alg = SchemeIO.read(RAW_JSON);
		int dups = MapleLRPParser.countDuplicateUVPairs(alg);
		System.out.println("\n=== Raw r2945 scheme duplicate (U,V) pair count: " + dups + " ===");
	}

	@Test
	@EnabledIf("lrpAvailable")
	public void report_row_extents_polynomial_LRP() throws Exception {
		NonCubicBilinearAlgorithm alg = MapleLRPParser.parse(LRP, 17, 17, 17);
		int[][] extsA = MapleLRPParser.rowExtentsA(alg);
		int[][] extsB = MapleLRPParser.rowExtentsB(alg);

		// Histogram of (rowsA, colsA) on the A side.
		java.util.Map<String, Integer> histA = new java.util.TreeMap<>();
		for (int[] e : extsA) {
			String k = "(" + e[0] + "," + e[1] + ")";
			histA.merge(k, 1, Integer::sum);
		}
		System.out.println("\n=== LRP A-side (rowExtent, colExtent) histogram ===");
		histA.forEach((k, v) -> System.out.println("  " + k + " : " + v));

		java.util.Map<String, Integer> histB = new java.util.TreeMap<>();
		for (int[] e : extsB) {
			String k = "(" + e[0] + "," + e[1] + ")";
			histB.merge(k, 1, Integer::sum);
		}
		System.out.println("\n=== LRP B-side (rowExtent, colExtent) histogram ===");
		histB.forEach((k, v) -> System.out.println("  " + k + " : " + v));

		// Joint product-shape histogram: for each k, the triple
		// (rowsA × colsA × colsB) approximates the sub-MM shape ⟨a, b, c⟩.
		// (We use colsA = rowsB necessarily; here we use colsA as the inner
		// dim and rowsA × colsB for the outer dims.)
		java.util.Map<String, Integer> histShape = new java.util.TreeMap<>();
		for (int k = 0; k < alg.r; k++) {
			int a = extsA[k][0], b = extsA[k][1], c = extsB[k][1];
			String key = "⟨" + a + "," + b + "," + c + "⟩";
			histShape.merge(key, 1, Integer::sum);
		}
		System.out.println("\n=== LRP per-product sub-shape ⟨rowsA, colsA, colsB⟩ histogram (top 20) ===");
		histShape.entrySet().stream()
				.sorted(java.util.Map.Entry.<String, Integer>comparingByValue().reversed())
				.limit(20)
				.forEach(e -> System.out.println("  " + e.getKey() + " : " + e.getValue()));

		// Quick sanity: total count by inner-dim bucket
		int[] byInner = new int[18];
		for (int k = 0; k < alg.r; k++) byInner[extsA[k][1]]++;
		System.out.println("\n=== LRP per-product A-colExtent (= inner-dim of sub-MM) sum ===");
		for (int i = 0; i <= 17; i++) if (byInner[i] > 0)
			System.out.println("  colsA=" + i + " : " + byInner[i] + " products");
	}
}
