package eu.solven.matmul.search;

import eu.solven.matmul.recombination.AnalyticalMaskSearch;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.SymmetryTransforms;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.AnalyticalMaskSearch.MaskCandidate;
import eu.solven.matmul.recombination.AnalyticalMaskSearch.SchemeSupports;

/**
 * Tests that {@link AnalyticalMaskSearch} matches brute-force
 * {@code Recombination.recombineWithAllocation(applyAxisFlip(scheme, mask), alloc)}
 * for all 8 masks, then asserts the headline findings (⟨17,17,17⟩ winners).
 */
class TestAnalyticalMaskSearch {

	private final FieldAwareLookup lookup = newLookup();
	private final CitedBound sota = new CitedBound(lookup);

	private static FieldAwareLookup newLookup() {
		try {
			return new FieldAwareLookup("Q");
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	private static NonCubicBilinearAlgorithm load(String name) throws Exception {
		return SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/" + name));
	}

	/**
	 * Cross-check #1: for every mask × every allocation we try, the analytical
	 * shape multiset must equal (after sort) the brute-force shape multiset
	 * produced by applying the mask to the scheme then calling
	 * recombineWithAllocation with the original allocation.
	 *
	 * <p>This is the correctness contract — without it, no other test result is
	 * meaningful.
	 */
	@Test
	void crossCheck_strassen_matches_bruteForce() throws Exception {
		NonCubicBilinearAlgorithm strassen = load("strassen-2x2x2_m7_a18.json");
		crossCheckScheme(strassen, new int[][] {
				{9, 8}, {8, 9}, {10, 7}, {16, 16}, {17, 0}, {1, 16}
		});
	}

	@Test
	void crossCheck_winograd_matches_bruteForce() throws Exception {
		NonCubicBilinearAlgorithm winograd = load("winograd_1971-2x2x2_m7_a24.json");
		crossCheckScheme(winograd, new int[][] {
				{9, 8}, {8, 9}, {10, 7}, {16, 16}, {12, 5}
		});
	}

	@Test
	void crossCheck_alphatensorZ_matches_bruteForce() throws Exception {
		NonCubicBilinearAlgorithm at = load("alphatensor_Z-2x2x2_m7_a22.json");
		crossCheckScheme(at, new int[][] {
				{9, 8}, {8, 9}, {16, 16}
		});
	}

	private void crossCheckScheme(NonCubicBilinearAlgorithm canonical, int[][] allocs) {
		SotaResolver resolver = sota::getRank;
		SchemeSupports supports = SchemeSupports.extract(canonical);
		for (int[] allocA : allocs) {
			for (int[] allocB : allocs) {
				for (int[] allocC : allocs) {
					for (int mask = 0; mask < 8; mask++) {
						int[] aA = ((mask & 1) != 0) ? reverse(allocA) : allocA;
						int[] aB = ((mask & 2) != 0) ? reverse(allocB) : allocB;
						int[] aC = ((mask & 4) != 0) ? reverse(allocC) : allocC;
						int[][] analytical = AnalyticalMaskSearch.shapesAt(supports, aA, aB, aC);
						long analyticalCost = AnalyticalMaskSearch.costOf(analytical, resolver);

						NonCubicBilinearAlgorithm masked = applyMaskByOrbitIndex(canonical, mask);
						Recombination.Result bf = Recombination.recombineWithAllocation(
								masked, resolver, allocA, allocB, allocC);
						assertThat(sortedMultiset(analytical))
								.as("scheme=%s allocA=%s allocB=%s allocC=%s mask=%d",
										labelOf(canonical), Arrays.toString(allocA),
										Arrays.toString(allocB), Arrays.toString(allocC), mask)
								.isEqualTo(sortedMultiset(bf.smallMatrixSizes));
						assertThat(analyticalCost)
								.as("cost mismatch mask=%d", mask)
								.isEqualTo(bf.totalRank);
					}
				}
			}
		}
	}

	/**
	 * Mask index → (swapA, swapB, swapC) bits using
	 * {@link SymmetryTransforms#axisFlipOrbit}'s convention: bit 0 = A, 1 = B,
	 * 2 = C.
	 */
	private static NonCubicBilinearAlgorithm applyMaskByOrbitIndex(NonCubicBilinearAlgorithm alg, int mask) {
		// axisFlipOrbit iterates mask 0..7 in this exact bit order, so the
		// list index equals the mask.
		List<NonCubicBilinearAlgorithm> orbit = SymmetryTransforms.axisFlipOrbit(alg);
		// orbit may have de-duplicated some masks (e.g. on fully-symmetric
		// schemes), so we can't just index by mask. Rebuild via the bits.
		boolean sA = (mask & 1) != 0;
		boolean sB = (mask & 2) != 0;
		boolean sC = (mask & 4) != 0;
		return applyAxisFlipReflectively(alg, sA, sB, sC);
	}

	/** Mirror of SymmetryTransforms.applyAxisFlip (which is package-private). */
	private static NonCubicBilinearAlgorithm applyAxisFlipReflectively(NonCubicBilinearAlgorithm alg,
			boolean swapA, boolean swapB, boolean swapC) {
		int a = alg.n, b = alg.m, c = alg.p, r = alg.r;
		double[][] srcU = alg.denseU();
		double[][] srcV = alg.denseV();
		double[][] srcW = alg.denseW();
		double[][] U2 = new double[a * b][r];
		double[][] V2 = new double[b * c][r];
		double[][] W2 = new double[a * c][r];
		for (int i = 0; i < a; i++) for (int j = 0; j < b; j++) {
			int iP = swapA ? (a - 1 - i) : i;
			int jP = swapB ? (b - 1 - j) : j;
			for (int k = 0; k < r; k++) U2[iP * b + jP][k] = srcU[i * b + j][k];
		}
		for (int j = 0; j < b; j++) for (int l = 0; l < c; l++) {
			int jP = swapB ? (b - 1 - j) : j;
			int lP = swapC ? (c - 1 - l) : l;
			for (int k = 0; k < r; k++) V2[jP * c + lP][k] = srcV[j * c + l][k];
		}
		for (int i = 0; i < a; i++) for (int l = 0; l < c; l++) {
			int iP = swapA ? (a - 1 - i) : i;
			int lP = swapC ? (c - 1 - l) : l;
			for (int k = 0; k < r; k++) W2[iP * c + lP][k] = srcW[i * c + l][k];
		}
		return new NonCubicBilinearAlgorithm(a, b, c, U2, V2, W2);
	}

	private static String labelOf(NonCubicBilinearAlgorithm a) {
		return String.format("⟨%d,%d,%d⟩ r=%d", a.n, a.m, a.p, a.r);
	}

	private static List<String> sortedMultiset(int[][] shapes) {
		String[] s = new String[shapes.length];
		for (int i = 0; i < shapes.length; i++) {
			int[] x = shapes[i];
			s[i] = "⟨" + x[0] + "," + x[1] + "," + x[2] + "⟩";
		}
		Arrays.sort(s);
		return Arrays.asList(s);
	}

	private static int[] reverse(int[] a) {
		int n = a.length;
		int[] r = new int[n];
		for (int i = 0; i < n; i++) r[i] = a[n - 1 - i];
		return r;
	}

	// ─────────────────────────────────────────────────────────────────────
	// Headline findings
	// ─────────────────────────────────────────────────────────────────────

	@Test
	void winograd_at_17x17x17_98_finds_mask_1_with_2930() throws Exception {
		NonCubicBilinearAlgorithm winograd = load("winograd_1971-2x2x2_m7_a24.json");
		int[] alloc = {9, 8};
		List<MaskCandidate> top = AnalyticalMaskSearch.topKMasks(
				winograd, alloc, alloc, alloc, sota::getRank, 3);
		assertThat(top).isNotEmpty();
		assertThat(top.get(0).cost).isEqualTo(2930L);
	}

	@Test
	void strassen_at_17x17x17_98_canonical_mask_is_best_at_2940() throws Exception {
		NonCubicBilinearAlgorithm strassen = load("strassen-2x2x2_m7_a18.json");
		int[] alloc = {9, 8};
		List<MaskCandidate> top = AnalyticalMaskSearch.topKMasks(
				strassen, alloc, alloc, alloc, sota::getRank, 8);
		// Strassen has multiple distinct multisets under axis-flip even at
		// (9,8)³: the min(U_view, W_view) reduction lets some masks produce
		// a ⟨8,8,8⟩ sub-product, but they come paired with duplicated ⟨9,9,9⟩
		// → net higher cost (2944). Mask=0 (canonical) remains the best at
		// 2940.
		assertThat(top).isNotEmpty();
		assertThat(top.get(0).mask).isEqualTo(0);
		assertThat(top.get(0).cost).isEqualTo(2940L);
		// Verify the secondary cluster exists and is worse
		assertThat(top.size()).isGreaterThanOrEqualTo(2);
		assertThat(top.get(1).cost).isGreaterThan(2940L);
	}

	@Test
	void winograd_at_16x16x16_88_balanced_all_masks_equivalent() throws Exception {
		NonCubicBilinearAlgorithm winograd = load("winograd_1971-2x2x2_m7_a24.json");
		int[] alloc = {8, 8};
		List<MaskCandidate> top = AnalyticalMaskSearch.topKMasks(
				winograd, alloc, alloc, alloc, sota::getRank, 8);
		// Balanced (8,8) — all 8 masks produce identical multiset.
		assertThat(top).hasSize(1);
		// 7 × R(⟨8,8,8⟩) = 7 × 336 = 2352
		assertThat(top.get(0).cost).isEqualTo(7L * 336);
	}
}
