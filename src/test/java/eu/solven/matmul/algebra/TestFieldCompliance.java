package eu.solven.matmul.algebra;

import eu.solven.matmul.papers.laderman1976.Laderman23;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.List;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Validates the {@link FieldCompliance} coefficient-field checker.
 */
public class TestFieldCompliance {

	@Test
	public void strassen_222_passes_Z() throws Exception {
		NonCubicBilinearAlgorithm s = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		assertThat(FieldCompliance.isCompliant(s, Field.Z)).isTrue();
		assertThat(FieldCompliance.isCompliant(s, Field.Q)).isTrue();
		assertThat(FieldCompliance.isCompliant(s, Field.R)).isTrue();
		assertThat(FieldCompliance.isCompliant(s, Field.F2)).isTrue();
	}

	@Test
	public void at_q_3x4x11_m103_a708_passes_Q_fails_Z() throws Exception {
		// The re-fetched Q scheme uses 0.5 coefficients — passes Q, fails Z.
		NonCubicBilinearAlgorithm s = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section11/alphatensor_Q-3x4x11_m103_a708.json"));
		assertThat(FieldCompliance.isCompliant(s, Field.Q)).isTrue();
		assertThat(FieldCompliance.isCompliant(s, Field.R)).isTrue();
		// 0.5 is not in Z — should fail with a precise location.
		List<FieldCompliance.Discrepancy> diffs =
				FieldCompliance.checkAllInField(s, Field.Z, 5);
		assertThat(diffs).as("Q scheme with 0.5 must fail Z compliance").isNotEmpty();
		System.out.println("First Q-coefficient violations under Z claim:");
		for (FieldCompliance.Discrepancy d : diffs) System.out.println("  " + d);
		// And not in F2 either.
		assertThat(FieldCompliance.isCompliant(s, Field.F2)).isFalse();
	}

	@Test
	public void laderman_passes_Z() {
		NonCubicBilinearAlgorithm l = NonCubicBilinearAlgorithm.fromCubic(
				eu.solven.matmul.papers.laderman1976.Laderman23.get());
		assertThat(FieldCompliance.isCompliant(l, Field.Z)).isTrue();
	}
}
