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
 * {@link PerAxisFingerprint} must (1) reproduce the {@link AnalyticalMaskSearch}
 * cost exactly, and (2) actually compress on at least one axis for a dense base
 * (where the full-tuple {@link BaseFingerprint} did not).
 */
public class TestPerAxisFingerprint {

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
	public void costMatchesShapesAt() throws Exception {
		SotaResolver sota = sota();
		NonCubicBilinearAlgorithm ae = read("src/main/resources/schemes/known/section5/alphaevolve-5x5x5_m93_a846.json");
		SchemeSupports sup = SchemeSupports.extract(ae);
		PerAxisFingerprint fp = PerAxisFingerprint.of(ae);
		for (int[][] alloc : new int[][][] {
				{ { 4, 3, 3, 3, 3 }, { 4, 3, 3, 3, 3 }, { 4, 3, 3, 3, 3 } },
				{ { 3, 4, 3, 3, 3 }, { 3, 3, 4, 3, 3 }, { 3, 3, 3, 4, 3 } },
				{ { 8, 4, 2, 1, 1 }, { 4, 3, 3, 3, 3 }, { 12, 1, 1, 1, 1 } } }) {
			assertThat(fp.cost(alloc[0], alloc[1], alloc[2], sota))
					.isEqualTo(refCost(sup, sota, alloc[0], alloc[1], alloc[2]));
		}
	}

	@Test
	public void perAxisCompressesWhereFullTupleDidNot() throws Exception {
		NonCubicBilinearAlgorithm ae = read("src/main/resources/schemes/known/section5/alphaevolve-5x5x5_m93_a846.json");
		PerAxisFingerprint fp = PerAxisFingerprint.of(ae);
		assertThat(fp.r).isEqualTo(93);
		// Each axis has at most r groups, and at least one axis strictly compresses.
		assertThat(fp.distinctA()).isBetween(1, 93);
		assertThat(fp.distinctB()).isBetween(1, 93);
		assertThat(fp.distinctC()).isBetween(1, 93);
		int minD = Math.min(fp.distinctA(), Math.min(fp.distinctB(), fp.distinctC()));
		assertThat(minD).as("at least one axis compresses below r=93").isLessThan(93);
	}
}
