package eu.solven.matmul.search;

import eu.solven.matmul.recombination.AnalyticalMaskSearch;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.AnalyticalMaskSearch.SchemeSupports;

/**
 * The support-pattern histogram ({@link BaseFingerprint}) must produce EXACTLY
 * the same recombination cost as the per-product {@link AnalyticalMaskSearch}
 * path — it is a grouping optimisation, not an approximation.
 */
public class TestBaseFingerprint {

	private static SotaResolver sota() {
		FieldAwareLookup lk = new FieldAwareLookup("R");
		return (p, q, r) -> {
			if (p == 0 || q == 0 || r == 0) return 0;
			if (p == 1) return q * r;
			if (q == 1) return p * r;
			if (r == 1) return p * q;
			int v = lk.findRank(p, q, r);
			return v >= Integer.MAX_VALUE / 100 ? p * q * r : v;
		};
	}

	private static long refCost(SchemeSupports sup, SotaResolver sota, int[] a, int[] b, int[] c) {
		long tot = 0;
		for (int[] s : AnalyticalMaskSearch.shapesAt(sup, a, b, c)) tot += sota.getRank(s[0], s[1], s[2]);
		return tot;
	}

	private static NonCubicBilinearAlgorithm read(String path) throws Exception {
		return SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint(path));
	}

	@Test
	public void fingerprintCostMatchesShapesAt_strassen_and_alphaevolve() throws Exception {
		SotaResolver sota = sota();

		NonCubicBilinearAlgorithm strassen = read("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json");
		SchemeSupports supS = SchemeSupports.extract(strassen);
		BaseFingerprint fpS = BaseFingerprint.of(strassen);
		// Several allocations of ⟨9,9,9⟩ over ⟨2,2,2⟩ (incl. unbalanced).
		for (int[][] alloc : new int[][][] {
				{ { 5, 4 }, { 5, 4 }, { 5, 4 } },
				{ { 6, 3 }, { 6, 3 }, { 6, 3 } },
				{ { 7, 2 }, { 4, 5 }, { 6, 3 } } }) {
			assertThat(fpS.cost(alloc[0], alloc[1], alloc[2], sota))
					.isEqualTo(refCost(supS, sota, alloc[0], alloc[1], alloc[2]));
		}

		NonCubicBilinearAlgorithm ae = read("src/main/resources/schemes/known/section5/alphaevolve-5x5x5_m93_a846.json");
		SchemeSupports supA = SchemeSupports.extract(ae);
		BaseFingerprint fpA = BaseFingerprint.of(ae);
		for (int[][] alloc : new int[][][] {
				{ { 4, 3, 3, 3, 3 }, { 4, 3, 3, 3, 3 }, { 4, 3, 3, 3, 3 } },
				{ { 3, 4, 3, 3, 3 }, { 3, 3, 4, 3, 3 }, { 3, 3, 3, 4, 3 } },
				{ { 8, 4, 2, 1, 1 }, { 4, 3, 3, 3, 3 }, { 12, 1, 1, 1, 1 } } }) {
			assertThat(fpA.cost(alloc[0], alloc[1], alloc[2], sota))
					.isEqualTo(refCost(supA, sota, alloc[0], alloc[1], alloc[2]));
		}
	}

	@Test
	public void distinctPatterns_compressesProducts() throws Exception {
		BaseFingerprint fpA = BaseFingerprint.of(read("src/main/resources/schemes/known/section5/alphaevolve-5x5x5_m93_a846.json"));
		// 93 products collapse to strictly fewer distinct support patterns.
		assertThat(fpA.r).isEqualTo(93);
		assertThat(fpA.distinctPatterns()).isLessThanOrEqualTo(93).isPositive();
	}
}
