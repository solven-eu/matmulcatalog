package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;

/**
 * {@link SubstitutionProjection} — FMM-Lille's {@code ProjOfAlongMode} lemma: dropping one
 * index of the p / m / n axis removes at least m / n / p products, from ANY representative.
 * Guards the 2026-11 finding behind issue #23's 44 FMM rows (all re-pricings by this lemma of
 * products we held) and the ranks it reaches, SOTA-or-better.
 */
public class TestSubstitutionProjection {

	private static final File STRASSEN = new File("src/main/resources/schemes/known/section2/2x2x2-r7-strassen-db11bcc.json");
	private static final File SMIRNOV_336 = new File("src/main/resources/schemes/known/section6/3x3x6-r40-smirnov13-d246191.json");

	@Test
	public void each_axis_drop_is_exact_and_removes_the_guaranteed_count() throws Exception {
		NonCubicBilinearAlgorithm s = SchemeIO.readBilinear(SMIRNOV_336); // ⟨3,3,6⟩ = 40
		int[][] expectShape = { { 2, 3, 6 }, { 3, 2, 6 }, { 3, 3, 5 } };
		int[] kills = { 6, 3, 3 }; // n-drop: p ; m-drop: n ; p-drop: m
		for (int axis = 0; axis < 3; axis++) {
			for (int idx = 0; idx < (axis == 2 ? 6 : 3); idx++) {
				NonCubicBilinearAlgorithm d = SubstitutionProjection.drop(s, axis, idx);
				assertThat(new int[] { d.n, d.m, d.p }).as("axis " + axis + " idx " + idx).containsExactly(expectShape[axis]);
				assertThat(d.r).as("axis " + axis + " idx " + idx).isLessThanOrEqualTo(40 - kills[axis]);
				assertThat(Verifier.isExactNonCubic(d)).as("axis " + axis + " idx " + idx).isTrue();
			}
			assertThat(SubstitutionProjection.guaranteedKills(s, axis)).isEqualTo(kills[axis]);
		}
	}

	@Test
	public void strassen_drops_to_an_exact_2x2x1_at_most_5() throws Exception {
		NonCubicBilinearAlgorithm s = SchemeIO.readBilinear(STRASSEN);
		NonCubicBilinearAlgorithm d = SubstitutionProjection.drop(s, 2, 1);
		assertThat(d.p).isEqualTo(1);
		assertThat(d.r).isLessThanOrEqualTo(5);
		assertThat(Verifier.isExactNonCubic(d)).isTrue();
	}

	/**
	 * The issue-#23 family, from the product we already held (coordinate projection: 0 saved on
	 * every axis). FMM's order m-then-p gives 1591 then 1574; p-then-m gives 1573. SOTA-or-better.
	 */
	@Test
	public void kronecker_9x18x18_projects_to_fmm_ranks_or_better() throws Exception {
		NonCubicBilinearAlgorithm a = SchemeIO.readBilinear(SMIRNOV_336);
		NonCubicBilinearAlgorithm k = Compose.kroneckerGeneral(a, a.orientAs(3, 6, 3).orElseThrow());
		assertThat(k.r).isEqualTo(1600);
		assertThat(ProjectionSearch.projectedRank(k, 9, 17, 18, 1)).as("coordinate projection sees nothing").isEqualTo(1600);

		NonCubicBilinearAlgorithm m17 = SubstitutionProjection.drop(k, 1, 17); // ⟨9,17,18⟩
		assertThat(m17.r).isLessThanOrEqualTo(1591);
		assertThat(Verifier.isExactNonCubic(m17)).isTrue();
		NonCubicBilinearAlgorithm fmmOrder = SubstitutionProjection.drop(m17, 2, 17); // ⟨9,17,17⟩
		assertThat(fmmOrder.r).isLessThanOrEqualTo(1574);
		assertThat(Verifier.isExactNonCubic(fmmOrder)).isTrue();

		NonCubicBilinearAlgorithm pFirst = SubstitutionProjection.drop(SubstitutionProjection.drop(k, 2, 17), 1, 17);
		assertThat(pFirst.r).as("p-then-m beats FMM's order").isLessThanOrEqualTo(1573);
		assertThat(Verifier.isExactNonCubic(pFirst)).isTrue();
	}
}
