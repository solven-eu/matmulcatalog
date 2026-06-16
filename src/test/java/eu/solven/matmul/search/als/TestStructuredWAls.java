package eu.solven.matmul.search.als;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;

/**
 * Regression guard for a SILENT false-negative (2026-06-12): a fixed ridge
 * (1e-8) in the ALS floored the residual near 1e-7, so even a scheme
 * warm-started from EXACT factors reported {@code solved()=false} — and any
 * de-novo "NO convergence" verdict was untrustworthy. The ridge must decay
 * with the residual so exact solutions pass the 1e-9 gate.
 */
public class TestStructuredWAls {

	/** Naive ⟨2,2,2⟩ (r=8): products sharing the output cell (i,k) share their
	 *  W direction → 4 tied classes of size 2. Warm-exact MUST solve. */
	@Test
	public void warm_exact_tied_classes_must_solve() {
		NonCubicBilinearAlgorithm naive = NonCubicBilinearAlgorithm.naive(2, 2, 2);
		// Product order is (i,j,k) lexicographic; class = i*2+k.
		int[] classOf = new int[8];
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				for (int k = 0; k < 2; k++) {
					classOf[(i * 2 + j) * 2 + k] = i * 2 + k;
				}
			}
		}
		double[][] warmW = new double[4][4];
		for (int i = 0; i < 2; i++) {
			for (int k = 0; k < 2; k++) {
				warmW[i * 2 + k][i * 2 + k] = 1;
			}
		}
		StructuredWAls.Result r = StructuredWAls.solve(2, 2, 2, classOf, 1, 50,
				naive.denseU(), naive.denseV(), warmW);
		assertThat(r.residual()).as("warm-exact residual must beat the solved() gate")
				.isLessThan(1e-9);
		assertThat(r.solved()).isTrue();
		assertThat(Verifier.passesRandomMatmulSpotCheck(
				StructuredWAls.expand(2, 2, 2, classOf, r))).isTrue();
	}

	/** Tied-U / tied-V go through cyclic tensor rotation — a wrong index
	 *  transpose in the back-map would yield a "solved" result that fails the
	 *  exactness spot-check. Fully asymmetric shape so every transpose bites;
	 *  warm-exact start so the test is deterministic (cold-start POWER is a
	 *  separate, calibrated concern — see the unconstrained-control protocol). */
	@Test
	public void tied_slot_rotation_round_trips() {
		record Rotated(int n, int m, int p) {}
		var rotatedOf = java.util.Map.of(
				eu.solven.matmul.catalog.SerendipitousBudProduct.BudType.W, new Rotated(2, 3, 4),
				eu.solven.matmul.catalog.SerendipitousBudProduct.BudType.U, new Rotated(3, 4, 2),
				eu.solven.matmul.catalog.SerendipitousBudProduct.BudType.V, new Rotated(4, 2, 3));
		for (var e : rotatedOf.entrySet()) {
			Rotated rot = e.getValue();
			NonCubicBilinearAlgorithm naive = NonCubicBilinearAlgorithm.naive(rot.n(), rot.m(), rot.p());
			int[] classOf = new int[24];
			for (int l = 0; l < 24; l++) {
				classOf[l] = l;
			}
			StructuredWAls.Result r = StructuredWAls.solve(rot.n(), rot.m(), rot.p(), classOf,
					1, 50, naive.denseU(), naive.denseV(), naive.denseW());
			assertThat(r.solved()).as("warm-exact on rotated shape %s must solve", rot).isTrue();
			assertThat(Verifier.passesRandomMatmulSpotCheck(
					StructuredWAls.expandTied(e.getKey(), 2, 3, 4, classOf, r)))
					.as("expandTied(%s) must map factors back to ⟨2,3,4⟩ exactly", e.getKey())
					.isTrue();
		}
	}

	/** De-novo all-singleton at the trivial rank must also converge — guards
	 *  against the ridge decay breaking cold-start behaviour. */
	@Test
	public void cold_start_naive_rank_converges() {
		int[] classOf = { 0, 1, 2, 3, 4, 5, 6, 7 };
		boolean any = false;
		for (long seed = 0; seed < 10 && !any; seed++) {
			any = StructuredWAls.solve(2, 2, 2, classOf, seed, 2_000, null, null, null).solved();
		}
		assertThat(any).as("rank-8 all-singleton ⟨2,2,2⟩ must solve within 10 restarts").isTrue();
	}
}
