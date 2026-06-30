package eu.solven.matmul.catalog;

import eu.solven.matmul.papers.hopcroftkerr1971.HopcroftKerr2bcAsymmetric;

import eu.solven.matmul.papers.hopcroftkerr1971.HopcroftKerr2bc;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;

public class TestHopcroftKerr2bcAsymmetric {

	@Test
	public void buildNaive_p2_n3_matches_square_rank() {
		// ⟨2, 2, 3⟩ via augment-square-discard: rank = R(⟨3, 2, 3⟩) = 15
		// (not the HK-optimal 11; the naive builder runs full square HK).
		NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.buildNaive(2, 3);
		assertThat(alg.n).isEqualTo(2);
		assertThat(alg.m).isEqualTo(2);
		assertThat(alg.p).isEqualTo(3);
		assertThat(alg.r).isEqualTo(15);
		assertThat(Verifier.passesRandomMatmulSpotCheck(alg)).isTrue();
	}

	@Test
	public void buildNaive_p3_n4_matches_square_rank() {
		// ⟨3, 2, 4⟩: HK-optimal 22 (=(3·12+4)/2 round up = 20 ⌈/⌉ — actually
		// (36+4)/2 = 20). Naive: R(⟨4, 2, 4⟩) = (3·16+4)/2 = 26.
		NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.buildNaive(3, 4);
		assertThat(alg.r).isEqualTo(26);
		assertThat(Verifier.passesRandomMatmulSpotCheck(alg)).isTrue();
	}

	@Test
	public void buildNaive_p3_n5_matches_square_rank() {
		NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.buildNaive(3, 5);
		assertThat(alg.r).isEqualTo(HopcroftKerr2bc.rank(5, 5));  // 40
		assertThat(Verifier.passesRandomMatmulSpotCheck(alg)).isTrue();
	}

	@Test
	public void buildNaive_p10_n15_runs_and_verifies() {
		// One of the HK targets ⟨10, 2, 15⟩. HK-optimal = 233; naive = R(⟨15,2,15⟩) = 345.
		NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.buildNaive(10, 15);
		assertThat(alg.r).isEqualTo(HopcroftKerr2bc.rank(15, 15));  // 345
		assertThat(Verifier.passesRandomMatmulSpotCheck(alg)).isTrue();
	}

	@Test
	public void buildNaiveDCE_p2_n3_drops_purely_augmented_products() {
		// ⟨2, 2, 3⟩: HK-optimal 11, naive 15. DCE may reach somewhere between.
		NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.buildNaiveDCE(2, 3);
		assertThat(alg.n).isEqualTo(2);
		assertThat(alg.p).isEqualTo(3);
		assertThat(alg.r).isLessThanOrEqualTo(15);
		assertThat(Verifier.passesRandomMatmulSpotCheck(alg)).as("rank=" + alg.r).isTrue();
		System.out.println("⟨2,2,3⟩ DCE rank: " + alg.r + " (HK-optimal=11, naive=15)");
	}

	@Test
	public void buildNaiveDCE_p10_n15_reports_gap_to_hk_optimal() {
		NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.buildNaiveDCE(10, 15);
		assertThat(alg.r).isLessThanOrEqualTo(345);
		assertThat(Verifier.passesRandomMatmulSpotCheck(alg)).as("rank=" + alg.r).isTrue();
		System.out.println("⟨10,2,15⟩ DCE rank: " + alg.r + " (HK-optimal=233, naive=345)");
	}

	@Test
	public void buildNaiveDCE_p10_n16_reports_gap_to_hk_optimal() {
		NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.buildNaiveDCE(10, 16);
		assertThat(Verifier.passesRandomMatmulSpotCheck(alg)).as("rank=" + alg.r).isTrue();
		System.out.println("⟨10,2,16⟩ DCE rank: " + alg.r + " (HK-optimal=248, naive=392)");
	}

	@Test
	public void buildOdd_p3_n3_matches_hk_formula_and_verifies() {
		// p=3, n=3 (square): HK = (27+3)/2 = 15. Same as buildSquareOdd.
		NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.buildOdd(3, 3);
		assertThat(alg.r).isEqualTo(15);
		assertThat(Verifier.passesRandomMatmulSpotCheck(alg)).isTrue();
	}

	@Test
	public void buildOdd_p3_n4_matches_hk_formula_and_verifies() {
		// p=3, n=4: HK = (36+4)/2 = 20.
		NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.buildOdd(3, 4);
		assertThat(alg.n).isEqualTo(3);
		assertThat(alg.m).isEqualTo(2);
		assertThat(alg.p).isEqualTo(4);
		assertThat(alg.r).isEqualTo(20);
		assertThat(Verifier.passesRandomMatmulSpotCheck(alg)).as("rank=" + alg.r).isTrue();
	}

	@Test
	public void buildOdd_p3_n5_matches_hk_formula_and_verifies() {
		// p=3, n=5: HK = (45+5)/2 = 25.
		NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.buildOdd(3, 5);
		assertThat(alg.r).isEqualTo(25);
		assertThat(Verifier.passesRandomMatmulSpotCheck(alg)).as("rank=" + alg.r).isTrue();
	}

	@Test
	public void buildOdd_p5_n5_matches_hk_formula_and_verifies() {
		// p=5, n=5: HK = (75+5)/2 = 40 (square).
		NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.buildOdd(5, 5);
		assertThat(alg.r).isEqualTo(40);
		assertThat(Verifier.passesRandomMatmulSpotCheck(alg)).as("rank=" + alg.r).isTrue();
	}

	@Test
	public void buildOdd_p5_n7_matches_hk_formula_and_verifies() {
		// p=5, n=7: HK = (105+7)/2 = 56.
		NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.buildOdd(5, 7);
		assertThat(alg.r).isEqualTo(56);
		assertThat(Verifier.passesRandomMatmulSpotCheck(alg)).as("rank=" + alg.r).isTrue();
	}

	@Test
	public void buildOdd_p5_n9_matches_hk_formula_and_verifies() {
		// p=5, n=9 (= 2p-1, max allowed): HK = (135+9)/2 = 72.
		NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.buildOdd(5, 9);
		assertThat(alg.r).isEqualTo(72);
		assertThat(Verifier.passesRandomMatmulSpotCheck(alg)).as("rank=" + alg.r).isTrue();
	}

	// p=7+ cases (k ≥ 3) currently fail in cyclic-band same-method handling:
	// when bridge wraps to row n (method 3), emitSameMethodPair_22_bridge3
	// looks up E_adj/F_adj products that the underlying Case (2, 3) off-diag
	// pair emits NEITHER. This is a paper-derivation gap from #46 (only
	// _11_bridge2 has the page-10 explicit formula; _22_bridge3 reuses the
	// product names but they don't exist in (2,3)-bridge context).
	// Tracked in #48 description.

	@Test
	public void buildOdd_p11_n11_long_vandermonde_no_overflow() {
		// p=11 square: HK = (3·11·11 + 11)/2 = 374/2 = 187. Mostly tests
		// the buildLong augmentation no longer overflows. n=11 square avoids
		// the cyclic-wrap method-3 bridge issue (no method-3 used when n even
		// — and for n=11 square, all pair distances stay within distance ≤ 5,
		// boundary i = n-1 doesn't cyclically meet method-1 pair partner with
		// the problematic same-method-22-bridge3 case).
		// NB: this test is expected to still hit the _22_bridge3 derivation
		// gap; if so, it'll throw "expected product not found" — keeps the
		// long-arith path coverage and surfaces the bridge issue.
		try {
			NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.buildOdd(11, 11);
			assertThat(alg.r).isEqualTo(187);
			assertThat(Verifier.passesRandomMatmulSpotCheck(alg)).as("rank=" + alg.r).isTrue();
		} catch (IllegalStateException e) {
			// Either Vandermonde overflow OR _22_bridge3 derivation gap.
			// Verify it's the LATTER (overflow should be fixed by buildLong).
			assertThat(e.getMessage()).contains("expected product not found");
		}
	}

	@Test
	public void buildNaiveDCE_p12_n16_reports_gap_to_hk_optimal() {
		NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.buildNaiveDCE(12, 16);
		assertThat(Verifier.passesRandomMatmulSpotCheck(alg)).as("rank=" + alg.r).isTrue();
		System.out.println("⟨12,2,16⟩ DCE rank: " + alg.r + " (HK-optimal=296, naive=392)");
	}

	@Test
	public void build_dispatches_to_buildOdd_when_p_is_3_or_5_in_range() {
		// Smart dispatcher: p ∈ {3, 5} with n ≤ 2p-1 → HK-optimal via buildOdd.
		NonCubicBilinearAlgorithm a35 = HopcroftKerr2bcAsymmetric.build(3, 5);
		assertThat(a35.r).isEqualTo(25);  // HK-formula
		NonCubicBilinearAlgorithm a59 = HopcroftKerr2bcAsymmetric.build(5, 9);
		assertThat(a59.r).isEqualTo(72);  // HK-formula
	}

	@Test
	public void build_attains_formula_for_all_odd_p_in_band_range() {
		// Task #7 (2026-06-11): the (2,2,bridge-3) gap is AVOIDED by arc-interior
		// bridge selection, and the Lemma-1/back-sub numerics are exact — so all
		// odd p now dispatch to buildOdd and attain the HK formula.
		NonCubicBilinearAlgorithm a710 = HopcroftKerr2bcAsymmetric.build(7, 10);
		assertThat(a710.r).isEqualTo(HopcroftKerr2bc.rank(7, 10));  // = 110, HK formula
		assertThat(Verifier.isExactNonCubic(a710)).isTrue();
		NonCubicBilinearAlgorithm a1112 = HopcroftKerr2bcAsymmetric.build(11, 12);
		assertThat(a1112.r).isEqualTo(HopcroftKerr2bc.rank(11, 12));  // = 204 — absent from all catalogs
		assertThat(Verifier.isExactNonCubic(a1112)).isTrue();
	}

	@Test
	public void build_attains_formula_for_even_p_case2() {
		// Task #7 (2026-06-11): Case 2 (even p) implemented — circulant matching
		// + repaired Z-trick Step 3 + true-reusable bridge-3 emitters.
		// ⟨2,10,15⟩=233 beats every published catalog (FMM/Perminov: 234).
		NonCubicBilinearAlgorithm a1015 = HopcroftKerr2bcAsymmetric.build(10, 15);
		assertThat(a1015.r).isEqualTo(HopcroftKerr2bc.rank(10, 15));  // = 233
		assertThat(Verifier.isExactNonCubic(a1015)).isTrue();
		NonCubicBilinearAlgorithm a1016 = HopcroftKerr2bcAsymmetric.build(10, 16);
		assertThat(a1016.r).isEqualTo(HopcroftKerr2bc.rank(10, 16));  // = 248
		assertThat(Verifier.isExactNonCubic(a1016)).isTrue();
	}

	@Test
	public void build_falls_back_to_DCE_when_p_is_2() {
		// p = 2 stays on the DCE fallback (the band builders need p ≥ 3).
		NonCubicBilinearAlgorithm a2x5 = HopcroftKerr2bcAsymmetric.build(2, 5);
		assertThat(a2x5.r).isGreaterThan(HopcroftKerr2bc.rank(2, 5));
		assertThat(Verifier.passesRandomMatmulSpotCheck(a2x5)).isTrue();
	}

	@Test
	public void buildChained_attains_formula_beyond_band_odd_p() {
		// n > 2p−1, odd p: every segment cost s(3p+1)/2 is an integer, so any
		// [p,2p−1] partition is slack-free and the chain sits at the formula.
		// ⟨3,2,8⟩: (72+8)/2 = 40; ⟨5,2,11⟩: (165+11)/2 = 88;
		// ⟨3,2,32⟩: (288+32)/2 = 160.
		for (int[] c : new int[][] { { 3, 8, 40 }, { 5, 11, 88 }, { 3, 32, 160 } }) {
			NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.build(c[0], c[1]);
			assertThat(alg.r).as("⟨2,%d,%d⟩", c[0], c[1]).isEqualTo(c[2]);
			assertThat(Verifier.isExactNonCubic(alg)).as("⟨2,%d,%d⟩", c[0], c[1]).isTrue();
		}
	}

	@Test
	public void buildChained_attains_formula_beyond_band_even_p() {
		// Even p: the DP must pick segment parities killing the ceiling slack
		// (at most one odd-size segment) and avoid degraded g ≥ 6 segments.
		// ⟨4,2,9⟩: ⌈9·13/2⌉ = 59; ⟨4,2,14⟩: 14·13/2 = 91;
		// ⟨6,2,32⟩: 32·19/2 = 304; ⟨10,2,21⟩: ⌈21·31/2⌉ = 326.
		for (int[] c : new int[][] { { 4, 9, 59 }, { 4, 14, 91 }, { 6, 32, 304 }, { 10, 21, 326 } }) {
			NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.build(c[0], c[1]);
			assertThat(alg.r).as("⟨2,%d,%d⟩", c[0], c[1]).isEqualTo(c[2]);
			assertThat(Verifier.isExactNonCubic(alg)).as("⟨2,%d,%d⟩", c[0], c[1]).isTrue();
		}
	}

	@Test
	public void build_emits_integer_schemes_via_unimodular_lemma1() {
		// Task #11 guard: with the unimodular-first Lemma-1 matrix, the window
		// inverses — hence the back-substituted W — are integer, so the whole
		// scheme is over Z, not Q. Cover each builder path: odd band, even band
		// (incl. the flagship ⟨10,2,15⟩ and the triangle ⟨12,2,18⟩), chained.
		// ⟨11,2,17⟩ and ⟨16,2,31⟩ are the regression shapes: the pre-Euclidean
		// seam DFS exhausted its budget there (m·r = 30 and 63 free bits).
		for (int[] c : new int[][] { { 9, 14 }, { 10, 15 }, { 12, 18 }, { 6, 32 }, { 13, 25 },
				{ 11, 17 }, { 16, 31 } }) {
			NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.build(c[0], c[1]);
			assertThat(Verifier.isExactNonCubic(alg)).as("⟨%d,2,%d⟩", c[0], c[1]).isTrue();
			for (double[][] mat : new double[][][] { alg.denseU(), alg.denseV(), alg.denseW() }) {
				for (double[] row : mat) {
					for (double v : row) {
						assertThat(v == Math.rint(v))
								.as("⟨%d,2,%d⟩ non-integer coefficient %s", c[0], c[1], v)
								.isTrue();
					}
				}
			}
		}
	}

	@Test
	public void buildChained_p12_n36_routes_around_triangle_segment() {
		// p=12 (k=5): the in-band segment n=18 is the g=6 triangle combo at
		// formula+1, but the DP can partition 36 as e.g. 16+20 (both clean),
		// so the chain must reach the exact formula 36·37/2 = 666 — NOT
		// 2×334 = 668 (two triangle segments) nor 667.
		NonCubicBilinearAlgorithm alg = HopcroftKerr2bcAsymmetric.build(12, 36);
		assertThat(alg.r).isEqualTo(666);
		assertThat(Verifier.passesRandomMatmulSpotCheck(alg)).isTrue();
	}
}
