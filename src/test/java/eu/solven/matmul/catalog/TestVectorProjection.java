package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.Random;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;

/**
 * {@link VectorProjection}: projection along a direction (FMM-Lille's {@code ProjOfAlongMode})
 * as a change of basis + coordinate drop. Guards the 2026-11 finding behind issue #23's
 * ⟨9,17,k⟩ family: every ⟨3,3,6⟩=40 representative we hold has coordinate projection margin
 * {@code [0,0,0]}, so the Kronecker stub ⟨9,18,18⟩=1600 was (correctly) pruned and never
 * projected, while FMM priced the same product down to ⟨9,17,18⟩:1591 and ⟨9,17,17⟩:1574 —
 * reachable once the factor is rebased.
 */
public class TestVectorProjection {

	private static final File STRASSEN = new File("src/main/resources/schemes/known/section2/2x2x2-r7-strassen-db11bcc.json");
	private static final File SMIRNOV_336 = new File("src/main/resources/schemes/known/section6/3x3x6-r40-smirnov13-d246191.json");

	/** The group action is exact for any invertible basis, on every axis (and composed). */
	@Test
	public void change_of_basis_on_any_axis_is_an_exact_scheme() throws Exception {
		NonCubicBilinearAlgorithm s = SchemeIO.readBilinear(STRASSEN);
		Random rnd = new Random(42);
		NonCubicBilinearAlgorithm cur = s;
		for (int axis = 0; axis < 3; axis++) {
			double[][] G = { { 1 + rnd.nextInt(3), rnd.nextInt(3) - 1 }, { rnd.nextInt(3) - 1, 1 + rnd.nextInt(2) } };
			if (Math.abs(G[0][0] * G[1][1] - G[0][1] * G[1][0]) < 1e-9) {
				G[0][1] += 1; // keep it invertible
			}
			cur = VectorProjection.changeBasis(cur, axis, G);
			assertThat(Verifier.isExactNonCubic(cur)).as("axis " + axis).isTrue();
			assertThat(cur.r).isEqualTo(7);
		}
	}

	/**
	 * Smirnov's ⟨3,3,6⟩=40: coordinate margin 0 everywhere, yet rank-1 factors share directions —
	 * the rebase exposes them as a coordinate margin, verified exact, same rank.
	 */
	@Test
	public void smirnov_3x3x6_rebases_from_margin_zero_to_a_positive_coordinate_margin() throws Exception {
		NonCubicBilinearAlgorithm a = SchemeIO.readBilinear(SMIRNOV_336);
		assertThat(ProjectionSearch.axisMargins(a)).containsExactly(0, 0, 0);
		int[] direction = VectorProjection.directionMargins(a);
		assertThat(direction[0] + direction[1] + direction[2]).as("some axis carries a common direction")
				.isGreaterThan(0);
		NonCubicBilinearAlgorithm r = VectorProjection.rebaseAll(a);
		assertThat(r).isNotSameAs(a);
		assertThat(r.r).isEqualTo(40);
		assertThat(Verifier.isExactNonCubic(r)).isTrue();
		int[] margins = ProjectionSearch.axisMargins(r);
		for (int axis = 0; axis < 3; axis++) {
			assertThat(margins[axis]).as("axis " + axis).isGreaterThanOrEqualTo(direction[axis]);
		}
	}

	/**
	 * The rebase is a change of representative, not of scheme: the Kronecker square of the
	 * rebased factor is still an exact ⟨9,18,18⟩=1600. (It does NOT reach FMM's ⟨9,17,18⟩:1591 —
	 * FMM's own June-2026 ⟨9,18,18⟩ artifact has direction margin 0 on every axis too, so
	 * {@code ProjOfAlongMode} saves its 9 by something beyond killing common-direction rank-1
	 * products; see research/DISCOVERIES_PENDING_ANALYSIS.md 2026-11-09.)
	 */
	@Test
	public void rebased_kronecker_square_is_still_the_exact_product() throws Exception {
		NonCubicBilinearAlgorithm r = VectorProjection.rebaseAll(SchemeIO.readBilinear(SMIRNOV_336));
		NonCubicBilinearAlgorithm k = Compose.kroneckerGeneral(r, r.orientAs(3, 6, 3).orElseThrow());
		assertThat(k.r).isEqualTo(1600);
		assertThat(k.n).isEqualTo(9);
		assertThat(k.m).isEqualTo(18);
		assertThat(k.p).isEqualTo(18);
		assertThat(Verifier.isExactNonCubic(k)).isTrue();
	}
}
