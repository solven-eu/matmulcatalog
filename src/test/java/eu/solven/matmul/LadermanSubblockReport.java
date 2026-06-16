package eu.solven.matmul;

import java.util.List;

import eu.solven.matmul.BilinearAlgorithm;
import eu.solven.matmul.papers.laderman1976.Laderman23;
import eu.solven.matmul.SubblockAnalyzer;

/**
 * Prints the sub-block embedding analysis for Laderman's r=23 ⟨3,3,3⟩
 * algorithm. Run via:
 *   mvn -q test-compile && java -cp target/classes:target/test-classes \
 *     eu.solven.matmul.LadermanSubblockReport
 */
public class LadermanSubblockReport {

	public static void main(String[] args) {
		BilinearAlgorithm laderman = Laderman23.get();
		System.out.println("Laderman ⟨3,3,3⟩ rank-23 algorithm.");
		System.out.printf("Total non-zero entries across U, V, W: %d/%d%n",
				countNonZero(laderman), 27 * 23);
		System.out.println();

		List<SubblockAnalyzer.Subblock> all = SubblockAnalyzer.enumerateAll222Subblocks(3);
		System.out.printf("All %d ⟨2,2,2⟩ sub-block embeddings of ⟨3,3,3⟩:%n%n", all.size());

		SubblockAnalyzer.Analysis analysis = SubblockAnalyzer.analyze(laderman, all);

		// 1. Per-subblock counts: must each be ≥ 7 (Hopcroft–Kerr)
		System.out.println("=== Per-sub-block: how many Laderman terms have non-vanishing restriction? ===");
		System.out.println("(R(⟨2,2,2⟩) = 7 forces every count ≥ 7)");
		System.out.println();
		System.out.printf("%-50s %10s%n", "sub-block", "# active");
		int minCount = Integer.MAX_VALUE, maxCount = 0;
		for (int s = 0; s < all.size(); s++) {
			int c = analysis.termsPerSubblock[s];
			minCount = Math.min(minCount, c);
			maxCount = Math.max(maxCount, c);
			System.out.printf("%-50s %10d%n", all.get(s).label(), c);
		}
		System.out.printf("%nRange across 27 sub-blocks: [min=%d, max=%d]%n", minCount, maxCount);
		System.out.printf("Mean: %.1f terms per sub-block%n%n",
				java.util.Arrays.stream(analysis.termsPerSubblock).average().orElse(0));

		// 2. Per-term counts: how many sub-blocks does each term contribute to?
		System.out.println("=== Per-term: how many of the 27 sub-blocks does each Laderman term contribute to? ===");
		System.out.printf("%5s  %10s%n", "term", "# subbl");
		int minSpt = Integer.MAX_VALUE, maxSpt = 0;
		for (int k = 0; k < laderman.r; k++) {
			int c = analysis.subblocksPerTerm[k];
			minSpt = Math.min(minSpt, c);
			maxSpt = Math.max(maxSpt, c);
			System.out.printf("M%-4d  %10d%n", k + 1, c);
		}
		System.out.printf("%nRange across 23 terms: [min=%d, max=%d]%n", minSpt, maxSpt);
		System.out.printf("Total non-vanishing (term, subblock) pairs: %d%n",
				java.util.Arrays.stream(analysis.subblocksPerTerm).sum());
		System.out.printf("Sanity: total = sum of per-subblock = sum of per-term: %d ?= %d%n",
				java.util.Arrays.stream(analysis.termsPerSubblock).sum(),
				java.util.Arrays.stream(analysis.subblocksPerTerm).sum());

		// 3. Histogram of "entanglement degree"
		System.out.println("\n=== Histogram: terms vs # sub-blocks they contribute to ===");
		int[] degHist = new int[28];
		for (int k = 0; k < laderman.r; k++) {
			degHist[analysis.subblocksPerTerm[k]]++;
		}
		for (int d = 0; d < degHist.length; d++) {
			if (degHist[d] > 0) {
				System.out.printf("%2d sub-blocks: %d terms %s%n", d, degHist[d],
						"*".repeat(degHist[d]));
			}
		}

		// 4. Corner sub-blocks specifically
		System.out.println("\n=== Just the 4 'corner' 2x2 sub-blocks ===");
		List<SubblockAnalyzer.Subblock> corners = SubblockAnalyzer.enumerateCornerSubblocks();
		SubblockAnalyzer.Analysis cornerAnalysis = SubblockAnalyzer.analyze(laderman, corners);
		String[] names = { "TL (rows 0,1; cols 0,1)", "(rows 0,1; cols 1,2)",
				"(rows 1,2; cols 0,1)", "BR (rows 1,2; cols 1,2)" };
		for (int s = 0; s < corners.size(); s++) {
			System.out.printf("  %s: %d active terms%n", names[s], cornerAnalysis.termsPerSubblock[s]);
		}

		// 5. ⟨2,2,3⟩-family embeddings (rank lower bound 11 each)
		System.out.println("\n=== ⟨2,2,3⟩-family (27 embeddings, each ≥ 11 active terms required) ===");
		List<SubblockAnalyzer.Subblock> sb223 = SubblockAnalyzer.enumerateAll223FamilySubblocks(3);
		SubblockAnalyzer.Analysis a223 = SubblockAnalyzer.analyze(laderman, sb223);
		printActivityHistogram(a223, 11);

		// 6. ⟨2,3,3⟩-family embeddings (rank lower bound 15 each)
		System.out.println("\n=== ⟨2,3,3⟩-family (9 embeddings, each ≥ 15 active terms required) ===");
		List<SubblockAnalyzer.Subblock> sb233 = SubblockAnalyzer.enumerateAll233FamilySubblocks(3);
		SubblockAnalyzer.Analysis a233 = SubblockAnalyzer.analyze(laderman, sb233);
		printActivityHistogram(a233, 15);
	}

	private static void printActivityHistogram(SubblockAnalyzer.Analysis a, int lowerBound) {
		int[] hist = new int[24];
		for (int s = 0; s < a.subblocks.size(); s++) hist[a.termsPerSubblock[s]]++;
		int minC = Integer.MAX_VALUE, maxC = 0;
		for (int s = 0; s < a.subblocks.size(); s++) {
			minC = Math.min(minC, a.termsPerSubblock[s]);
			maxC = Math.max(maxC, a.termsPerSubblock[s]);
		}
		System.out.printf("  range=[%d, %d], mean=%.1f, lower bound=%d, slack=[%d, %d]%n",
				minC, maxC,
				java.util.Arrays.stream(a.termsPerSubblock).average().orElse(0),
				lowerBound, minC - lowerBound, maxC - lowerBound);
		System.out.println("  histogram (active count → # embeddings):");
		for (int c = 0; c < hist.length; c++) {
			if (hist[c] > 0) {
				System.out.printf("    %2d: %d %s%n", c, hist[c], "*".repeat(hist[c]));
			}
		}
	}

	private static int countNonZero(BilinearAlgorithm alg) {
		int count = 0;
		for (int k = 0; k < alg.r; k++) {
			for (int i = 0; i < alg.n * alg.n; i++) {
				if (alg.U[i][k] != 0) count++;
				if (alg.V[i][k] != 0) count++;
				if (alg.W[i][k] != 0) count++;
			}
		}
		return count;
	}
}
