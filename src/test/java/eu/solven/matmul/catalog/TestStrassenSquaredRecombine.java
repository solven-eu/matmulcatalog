package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Recombination;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Smallest possible recombination test: ⟨4,4,4⟩ via Strassen[2,2]³ —
 * the textbook Strassen² recipe should give a 49-product algorithm
 * that passes the spot-check.
 *
 * <p>If this passes, recombination is sound for the simple case and
 * the n=21 bug lies in non-balanced / mixed-shape allocations. If
 * this FAILS, recombination itself has a fundamental issue.</p>
 */
public class TestStrassenSquaredRecombine {

	@Test
	public void strassen_squared_444_via_recombine_verifies() throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("R");

		NonCubicBilinearAlgorithm result = Recombination.constructWithAllocation(
				strassen, lookup, new int[]{2, 2}, new int[]{2, 2}, new int[]{2, 2});

		System.out.printf("Constructed ⟨%d,%d,%d⟩ with rank %d%n", result.n, result.m, result.p, result.r);
		assertThat(result.n).isEqualTo(4);
		assertThat(result.m).isEqualTo(4);
		assertThat(result.p).isEqualTo(4);
		assertThat(result.r).isEqualTo(49); // 7 outer × 7 inner

		assertThat(Verifier.passesRandomMatmulSpotCheck(result))
				.as("Strassen[2,2]³ recombination must produce a valid ⟨4,4,4⟩=49 algorithm")
				.isTrue();
	}

	/**
	 * Non-balanced uniform allocation: ⟨5,5,5⟩ via Strassen[3,2]³. Each
	 * axis split as (3,2). Inner sub-products are various shapes from
	 * {⟨3,3,3⟩, ⟨3,3,2⟩, ⟨3,2,3⟩, ⟨2,3,3⟩, ⟨3,2,2⟩, ⟨2,3,2⟩, ⟨2,2,3⟩,
	 * ⟨2,2,2⟩}; if Recombination handles them correctly, verifies.
	 */
	@Test
	public void strassen_555_via_recombine_3_2_uniform_verifies() throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("R");

		NonCubicBilinearAlgorithm result = Recombination.constructWithAllocation(
				strassen, lookup, new int[]{3, 2}, new int[]{3, 2}, new int[]{3, 2});

		System.out.printf("⟨5,5,5⟩ Strassen[3,2]³ → rank %d%n", result.r);
		assertThat(result.n).isEqualTo(5);

		boolean ok = Verifier.passesRandomMatmulSpotCheck(result);
		System.out.println("Spot-check: " + (ok ? "PASS" : "FAIL"));
		assertThat(ok).isTrue();
	}

	@Test
	public void strassen_777_via_recombine_4_3_uniform_verifies() throws Exception {
		// Sedoglavic's ⟨7,7,7⟩=250 recipe.
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("R");
		NonCubicBilinearAlgorithm result = Recombination.constructWithAllocation(
				strassen, lookup, new int[]{4, 3}, new int[]{4, 3}, new int[]{4, 3});
		System.out.printf("⟨7,7,7⟩ Strassen[4,3]³ → rank %d%n", result.r);
		assertThat(Verifier.passesRandomMatmulSpotCheck(result))
				.as("Sedoglavic ⟨7,7,7⟩=250 recipe must verify")
				.isTrue();
	}

	@Test
	public void strassen_777_via_recombine_3_4_uniform_verifies() throws Exception {
		// Same as [4,3]³ but axis-flipped — exercises transpose-orientation lookup.
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("R");
		NonCubicBilinearAlgorithm result = Recombination.constructWithAllocation(
				strassen, lookup, new int[]{3, 4}, new int[]{3, 4}, new int[]{3, 4});
		System.out.printf("⟨7,7,7⟩ Strassen[3,4]³ → rank %d%n", result.r);
		assertThat(Verifier.passesRandomMatmulSpotCheck(result))
				.as("⟨7,7,7⟩ axis-flipped allocation [3,4]³ must verify")
				.isTrue();
	}

	@Test
	public void scan_uniform_alloc_a_b_triples() throws Exception {
		// For each (a, b) with a<b, build ⟨a+b, a+b, a+b⟩ via Strassen[a,b]³.
		// Find the smallest (a, b) where the construction fails.
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("R");
		for (int a = 2; a <= 10; a++) {
			for (int b = a + 1; b <= a + 5; b++) {
				int n = a + b;
				try {
					NonCubicBilinearAlgorithm r = Recombination.constructWithAllocation(
							strassen, lookup, new int[]{a, b}, new int[]{a, b}, new int[]{a, b});
					boolean ok = Verifier.passesRandomMatmulSpotCheck(r);
					System.out.printf("[%d,%d]³ n=%d rank=%d %s%n", a, b, n, r.r, ok ? "PASS" : "FAIL");
				} catch (RuntimeException ex) {
					System.out.printf("[%d,%d]³ n=%d MISSING-INNER: %s%n", a, b, n, ex.getMessage());
				}
			}
		}
	}

	@Test
	public void strassen_151515_via_recombine_uniform_7_8_verifies() throws Exception {
		// Same n=15 but with [7,8] on ALL three axes — should work since it's "symmetric".
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("R");
		NonCubicBilinearAlgorithm result = Recombination.constructWithAllocation(
				strassen, lookup, new int[]{7, 8}, new int[]{7, 8}, new int[]{7, 8});
		System.out.printf("⟨15,15,15⟩ Strassen[7,8]³ → rank %d%n", result.r);
		assertThat(Verifier.passesRandomMatmulSpotCheck(result))
				.as("⟨15,15,15⟩ Strassen[7,8]³ uniform — spot-check")
				.isTrue();
	}

	@Test
	public void debug_n15_failing_case_symbolic_diff() throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("R");
		NonCubicBilinearAlgorithm result = Recombination.constructWithAllocation(
				strassen, lookup, new int[]{7, 8}, new int[]{7, 8}, new int[]{8, 7});
		System.out.printf("⟨15,15,15⟩ Strassen[7,8]/[7,8]/[8,7] → rank %d%n", result.r);

		// Originally this test caught a recombination bug — discrepancies were
		// EXPECTED until the bug was fixed. The underlying construction now
		// produces a correct scheme; the test inverts to "should NOT regress".
		java.util.List<Verifier.SymbolicDiff> diffs = Verifier.symbolicDiscrepancies(result, 5, 1e-10);
		if (!diffs.isEmpty()) {
			System.out.println("Regression detected — first " + diffs.size() + " discrepancies:");
			for (Verifier.SymbolicDiff d : diffs) System.out.println("  " + d);
		}
		assertThat(diffs).as("recombination should be regression-free").isEmpty();
	}

	@Test
	public void scan_growing_sizes_C_axis_flipped() throws Exception {
		// Bisect: for each n=4..21, find an allocation where A=B=[a,b] but C=[b,a].
		// Identify the smallest n where this construction fails.
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("R");
		int firstFailing = -1;
		for (int n = 4; n <= 21; n++) {
			int a = n / 2;
			int b = n - a;
			if (a == b) {
				// symmetric, axis-flip is the same — skip
				System.out.printf("n=%2d: skip (a=b=%d)%n", n, a);
				continue;
			}
			NonCubicBilinearAlgorithm result;
			try {
				result = Recombination.constructWithAllocation(
						strassen, lookup, new int[]{a, b}, new int[]{a, b}, new int[]{b, a});
			} catch (RuntimeException e) {
				System.out.printf("n=%2d: missing sub-shape: %s%n", n, e.getMessage());
				continue;
			}
			boolean ok = Verifier.passesRandomMatmulSpotCheck(result);
			System.out.printf("n=%2d [%d,%d]/[%d,%d]/[%d,%d] rank=%d  %s%n",
					n, a, b, a, b, b, a, result.r, ok ? "PASS" : "FAIL");
			if (!ok && firstFailing == -1) firstFailing = n;
		}
		System.out.println("First failing n with axis-flipped C: " + firstFailing);
	}

	@Test
	public void strassen_777_via_recombine_3_4_with_C_flipped_verifies() throws Exception {
		// Same pattern as n=21 failing case but at tiny n: A,B use one allocation
		// while C uses a DIFFERENT one. Strassen[3,4]/[3,4]/[4,3] on ⟨7,7,7⟩.
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("R");
		NonCubicBilinearAlgorithm result = Recombination.constructWithAllocation(
				strassen, lookup, new int[]{3, 4}, new int[]{3, 4}, new int[]{4, 3});
		System.out.printf("⟨7,7,7⟩ Strassen[3,4]/[3,4]/[4,3] → rank %d%n", result.r);
		boolean ok = Verifier.passesRandomMatmulSpotCheck(result);
		System.out.println("  spot-check: " + (ok ? "PASS" : "FAIL"));
		if (!ok) {
			java.util.List<Verifier.SymbolicDiff> diffs = Verifier.symbolicDiscrepancies(result);
			System.out.println("  first " + Math.min(10, diffs.size()) + " discrepancies:");
			for (int i = 0; i < Math.min(10, diffs.size()); i++) {
				System.out.println("    " + diffs.get(i));
			}
		}
		assertThat(ok).as("⟨7,7,7⟩ Strassen[3,4]/[3,4]/[4,3] must verify").isTrue();
	}

	@Test
	public void strassen_212121_via_recombine_balanced_11_10_verifies() throws Exception {
		// Balanced [11,10]³ — should give Sedoglavic-like result.
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("R");
		NonCubicBilinearAlgorithm result = Recombination.constructWithAllocation(
				strassen, lookup, new int[]{11, 10}, new int[]{11, 10}, new int[]{11, 10});
		System.out.printf("⟨21,21,21⟩ Strassen[11,10]³ → rank %d%n", result.r);
		assertThat(Verifier.passesRandomMatmulSpotCheck(result))
				.as("⟨21,21,21⟩ Strassen[11,10]³ balanced — spot-check")
				.isTrue();
	}

	@Test
	public void strassen_212121_via_recombine_nonbal_9_12_verifies() throws Exception {
		// Non-balanced [9,12]³ — symmetric across axes but uses [9,12] instead of [11,10].
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("R");
		NonCubicBilinearAlgorithm result = Recombination.constructWithAllocation(
				strassen, lookup, new int[]{9, 12}, new int[]{9, 12}, new int[]{9, 12});
		System.out.printf("⟨21,21,21⟩ Strassen[9,12]³ → rank %d%n", result.r);
		assertThat(Verifier.passesRandomMatmulSpotCheck(result))
				.as("⟨21,21,21⟩ Strassen[9,12]³ non-balanced symmetric — spot-check")
				.isTrue();
	}

	@Test
	public void strassen_141414_via_recombine_7_7_verifies() throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("R");
		NonCubicBilinearAlgorithm result = Recombination.constructWithAllocation(
				strassen, lookup, new int[]{7, 7}, new int[]{7, 7}, new int[]{7, 7});
		System.out.printf("⟨14,14,14⟩ Strassen[7,7]³ → rank %d%n", result.r);
		assertThat(Verifier.passesRandomMatmulSpotCheck(result))
				.as("⟨14,14,14⟩ Strassen[7,7]³ recombination must verify")
				.isTrue();
	}

	/**
	 * Non-balanced + non-symmetric allocation: ⟨5,5,5⟩ via Strassen[3,2]/[3,2]/[2,3].
	 * One axis allocation is reversed — checks that the per-axis allocation
	 * ordering matters and is handled correctly.
	 */
	@Test
	public void strassen_555_via_recombine_mixed_axis_orders_verifies() throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		FieldAwareLookup lookup = new FieldAwareLookup("R");

		NonCubicBilinearAlgorithm result = Recombination.constructWithAllocation(
				strassen, lookup, new int[]{3, 2}, new int[]{3, 2}, new int[]{2, 3});

		System.out.printf("⟨5,5,5⟩ Strassen[3,2]/[3,2]/[2,3] → rank %d%n", result.r);
		boolean ok = Verifier.passesRandomMatmulSpotCheck(result);
		System.out.println("Spot-check: " + (ok ? "PASS" : "FAIL"));
		assertThat(ok).isTrue();
	}
}
