package eu.solven.matmul;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.List;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.AxisSplitBases;
import eu.solven.matmul.papers.laderman1976.Laderman23;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.SymmetryTransforms;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Verifies the S₃ symmetry transforms preserve correctness of bilinear
 * algorithms (each orbit element computes the relevant matmul).
 */
public class TestSymmetryTransforms {

	private static NonCubicBilinearAlgorithm loadStrassen() throws Exception {
		return SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
	}

	@Test
	public void transpose_twice_is_identity() throws Exception {
		NonCubicBilinearAlgorithm s = loadStrassen();
		NonCubicBilinearAlgorithm s_tt = s.transpose().transpose();
		assertThat(s_tt.n).isEqualTo(s.n);
		assertThat(s_tt.m).isEqualTo(s.m);
		assertThat(s_tt.p).isEqualTo(s.p);
		assertThat(s_tt.r).isEqualTo(s.r);
		double[][] sU = s.denseU();
		double[][] sV = s.denseV();
		double[][] sW = s.denseW();
		double[][] ttU = s_tt.denseU();
		double[][] ttV = s_tt.denseV();
		double[][] ttW = s_tt.denseW();
		assertThat(ttU).isDeepEqualTo(sU);
		assertThat(ttV).isDeepEqualTo(sV);
		assertThat(ttW).isDeepEqualTo(sW);
	}

	@Test
	public void transpose_of_strassen_is_exact_matmul() throws Exception {
		NonCubicBilinearAlgorithm s = loadStrassen();
		NonCubicBilinearAlgorithm st = s.transpose();
		assertThat(st.n).isEqualTo(s.p);
		assertThat(st.m).isEqualTo(s.m);
		assertThat(st.p).isEqualTo(s.n);
		assertThat(Verifier.isExactNonCubic(st)).isTrue();
	}

	@Test
	public void s3_orbit_strassen_contains_six_or_fewer() throws Exception {
		NonCubicBilinearAlgorithm s = loadStrassen();
		List<NonCubicBilinearAlgorithm> orbit = SymmetryTransforms.s3Orbit(s);
		assertThat(orbit).hasSizeBetween(1, 6);
		for (NonCubicBilinearAlgorithm a : orbit) {
			assertThat(Verifier.isExactNonCubic(a)).isTrue();
			assertThat(a.n).isEqualTo(s.n);
			assertThat(a.m).isEqualTo(s.m);
			assertThat(a.p).isEqualTo(s.p);
			assertThat(a.r).isEqualTo(s.r);
		}
	}

	@Test
	public void s3_orbit_laderman_all_verify() throws Exception {
		NonCubicBilinearAlgorithm laderman = NonCubicBilinearAlgorithm.fromCubic(Laderman23.get());
		List<NonCubicBilinearAlgorithm> orbit = SymmetryTransforms.s3Orbit(laderman);
		assertThat(orbit).hasSizeBetween(1, 6);
		System.out.println("Laderman S₃ orbit size: " + orbit.size());
		for (NonCubicBilinearAlgorithm a : orbit) {
			assertThat(a.n).isEqualTo(3);
			assertThat(a.m).isEqualTo(3);
			assertThat(a.p).isEqualTo(3);
			assertThat(a.r).isEqualTo(23);
			assertThat(Verifier.isExactNonCubic(a)).isTrue();
		}
	}

	@Test
	public void mul211_orbit_includes_mul121_and_mul112() {
		NonCubicBilinearAlgorithm mul211 = AxisSplitBases.mul211();
		List<NonCubicBilinearAlgorithm> orbit = SymmetryTransforms.s3Orbit(mul211);
		// Among the orbit, formats (2,1,1), (1,2,1), (1,1,2) must all appear.
		boolean has211 = orbit.stream().anyMatch(a -> a.n == 2 && a.m == 1 && a.p == 1);
		boolean has121 = orbit.stream().anyMatch(a -> a.n == 1 && a.m == 2 && a.p == 1);
		boolean has112 = orbit.stream().anyMatch(a -> a.n == 1 && a.m == 1 && a.p == 2);
		assertThat(has211).as("orbit contains ⟨2,1,1⟩").isTrue();
		assertThat(has121).as("orbit contains ⟨1,2,1⟩").isTrue();
		assertThat(has112).as("orbit contains ⟨1,1,2⟩").isTrue();
		for (NonCubicBilinearAlgorithm a : orbit) {
			assertThat(Verifier.isExactNonCubic(a)).isTrue();
		}
	}
}
