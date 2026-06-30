package eu.solven.matmul.search;

import org.junit.jupiter.api.Tag;


import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.isotropy.PairedSubProducts;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Prototype: scan ALL 2×2-block allocations (n₀, m₀, p₀) of the target
 * via Strassen ⟨2,2,2⟩=7, compute the resulting 7 sub-product shapes
 * for each allocation, apply same-shape TA pair-fuse where profitable,
 * and report the top allocations by total rank.
 *
 * <p>Focused on the user's question: "stay with 2×2 Strassen split,
 * find the allocation that gives multiple same-cubic-shape products
 * enabling TA." For ⟨17,17,17⟩ FMM's recipe has 2× ⟨8,8,8⟩ + 2× ⟨9,9,9⟩;
 * this prototype asks whether ANY Strassen 2×2 allocation produces such
 * clusters AND beats the baseline 2940.</p>
 */
@Tag("slow")
class TestStrassen2x2AllocationScan {

	@Test
	void scan_17_17_17() throws Exception {
		scanTargetCubic(17);
	}

	@Test
	void scan_8_8_8() throws Exception {
		scanTargetCubic(8);
	}

	@Test
	void scan_9_9_9() throws Exception {
		scanTargetCubic(9);
	}

	@Test
	void scan_3_3_3() throws Exception {
		scanTargetCubic(3);
	}

	@Test
	void scan_6_6_6() throws Exception {
		scanTargetCubic(6);
	}

	@Test
	void scan_2_2_3_noncubic() throws Exception {
		scanTargetNonCubic(2, 2, 3);
	}

	@Test
	void scan_3_4_5_noncubic() throws Exception {
		scanTargetNonCubic(3, 4, 5);
	}

	private static void scanTargetNonCubic(int n, int m, int p) throws Exception {
		File strassenFile = eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json");
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(strassenFile);
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		CitedBound sota = new CitedBound(lookup);

		List<Candidate> candidates = new ArrayList<>();

		for (int n0 = 1; n0 < n; n0++) {
			for (int m0 = 1; m0 < m; m0++) {
				for (int p0 = 1; p0 < p; p0++) {
					int n1 = n - n0, m1 = m - m0, p1 = p - p0;
					if (n1 <= 0 || m1 <= 0 || p1 <= 0) continue;
					int[] allocA = { n0, n1 };
					int[] allocB = { m0, m1 };
					int[] allocC = { p0, p1 };
					Recombination.Result r;
					try {
						r = Recombination.recombineWithAllocation(strassen, sota, allocA, allocB, allocC);
					} catch (Exception e) {
						continue;
					}
					long noTA = r.totalRank;
					long withTA = PairedSubProducts.applyPairing(r.smallMatrixSizes, sota);
					boolean taHelped = withTA < noTA;
					candidates.add(new Candidate(n0, m0, p0, r.smallMatrixSizes, noTA, withTA, taHelped));
				}
			}
		}
		candidates.sort(Comparator.comparingLong(c -> c.bestRank()));

		System.out.printf("%n=== Strassen 2×2 allocation scan for ⟨%d,%d,%d⟩ ===%n", n, m, p);
		System.out.printf("Total allocations scanned: %d%n", candidates.size());
		System.out.printf("%nTop 5 by total rank:%n");
		for (int i = 0; i < Math.min(5, candidates.size()); i++) {
			Candidate c = candidates.get(i);
			System.out.printf("  alloc=(%d,%d)x(%d,%d)x(%d,%d)  rank: noTA=%d, withTA=%d%s  shapes=%s%n",
					c.n0, n - c.n0, c.m0, m - c.m0, c.p0, p - c.p0,
					c.noTARank, c.withTARank, c.taHelped ? " (TA SAVED)" : "", countShapes(c.shapes));
		}
		long anyTA = candidates.stream().filter(c -> c.taHelped).count();
		System.out.printf("Allocations where TA helps: %d / %d%n", anyTA, candidates.size());
	}

	private static void scanTargetCubic(int n) throws Exception {
		File strassenFile = eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json");
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(strassenFile);
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		CitedBound sota = new CitedBound(lookup);

		List<Candidate> candidates = new ArrayList<>();
		long baselineRank = -1;  // (⌈n/2⌉, ⌊n/2⌋) on each axis = "default" Strassen

		for (int n0 = 1; n0 < n; n0++) {
			for (int m0 = 1; m0 < n; m0++) {
				for (int p0 = 1; p0 < n; p0++) {
					int n1 = n - n0, m1 = n - m0, p1 = n - p0;
					int[] allocA = { n0, n1 };
					int[] allocB = { m0, m1 };
					int[] allocC = { p0, p1 };
					Recombination.Result r;
					try {
						r = Recombination.recombineWithAllocation(strassen, sota, allocA, allocB, allocC);
					} catch (Exception e) {
						continue;
					}
					long noTA = r.totalRank;
					long withTA = PairedSubProducts.applyPairing(r.smallMatrixSizes, sota);
					long best = Math.min(noTA, withTA);
					boolean taHelped = withTA < noTA;
					candidates.add(new Candidate(n0, m0, p0, r.smallMatrixSizes, noTA, withTA, taHelped));
					if (n0 == (n + 1) / 2 && m0 == (n + 1) / 2 && p0 == (n + 1) / 2) {
						baselineRank = best;
					}
				}
			}
		}

		candidates.sort(Comparator.comparingLong(c -> c.bestRank()));

		System.out.printf("%n=== Strassen 2×2 allocation scan for ⟨%d,%d,%d⟩ ===%n", n, n, n);
		System.out.printf("Baseline (⌈%d/2⌉, ⌊%d/2⌋)^3 rank: %d%n", n, n, baselineRank);
		System.out.printf("Total allocations scanned: %d%n", candidates.size());
		System.out.printf("%n--- Top 10 by total rank ---%n");
		int show = Math.min(10, candidates.size());
		for (int i = 0; i < show; i++) {
			Candidate c = candidates.get(i);
			Map<String, Integer> shapeCounts = countShapes(c.shapes);
			System.out.printf("  alloc=(%d,%d)x(%d,%d)x(%d,%d)  rank: noTA=%d, withTA=%d%s  shapes=%s%n",
					c.n0, n - c.n0, c.m0, n - c.m0, c.p0, n - c.p0,
					c.noTARank, c.withTARank, c.taHelped ? " (TA SAVED!)" : "",
					shapeCounts);
		}

		System.out.printf("%n--- Allocations where TA pair-fuse strictly helps ---%n");
		long taHelpedCount = candidates.stream().filter(c -> c.taHelped).count();
		System.out.printf("Count: %d / %d%n", taHelpedCount, candidates.size());
		int shown = 0;
		for (Candidate c : candidates) {
			if (!c.taHelped) continue;
			if (shown >= 8) break;
			Map<String, Integer> shapeCounts = countShapes(c.shapes);
			System.out.printf("  alloc=(%d,%d)x(%d,%d)x(%d,%d)  noTA=%d → withTA=%d (saved %d)  shapes=%s%n",
					c.n0, n - c.n0, c.m0, n - c.m0, c.p0, n - c.p0,
					c.noTARank, c.withTARank, c.noTARank - c.withTARank, shapeCounts);
			shown++;
		}

		System.out.printf("%n--- Allocations with ≥2 cubic same-shape products (TA candidates) ---%n");
		int sameCubic = 0;
		for (Candidate c : candidates) {
			Map<String, Integer> shapeCounts = countCubicSameShape(c.shapes);
			boolean hasPair = shapeCounts.values().stream().anyMatch(v -> v >= 2);
			if (!hasPair) continue;
			sameCubic++;
			if (sameCubic > 8) continue;
			System.out.printf("  alloc=(%d,%d)x(%d,%d)x(%d,%d)  rank=%d  cubic-clusters=%s%n",
					c.n0, n - c.n0, c.m0, n - c.m0, c.p0, n - c.p0, c.bestRank(), shapeCounts);
		}
		System.out.printf("Total allocations with ≥2 same-cubic products: %d%n", sameCubic);
		System.out.println("(If best of these still ≥ baseline, Strassen-2×2 cannot reach FMM via simple TA on same-cubic clusters.)");
	}

	private static Map<String, Integer> countShapes(int[][] shapes) {
		Map<String, Integer> out = new HashMap<>();
		for (int[] s : shapes) {
			out.merge("⟨" + s[0] + "," + s[1] + "," + s[2] + "⟩", 1, Integer::sum);
		}
		return out;
	}

	private static Map<String, Integer> countCubicSameShape(int[][] shapes) {
		Map<String, Integer> out = new HashMap<>();
		for (int[] s : shapes) {
			if (s[0] == s[1] && s[1] == s[2]) {
				out.merge("⟨" + s[0] + "," + s[1] + "," + s[2] + "⟩", 1, Integer::sum);
			}
		}
		return out;
	}

	private record Candidate(int n0, int m0, int p0, int[][] shapes,
			long noTARank, long withTARank, boolean taHelped) {
		long bestRank() { return Math.min(noTARank, withTARank); }
	}
}
