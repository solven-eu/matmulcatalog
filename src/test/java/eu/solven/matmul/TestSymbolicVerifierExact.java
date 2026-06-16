package eu.solven.matmul;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.SymbolicVerifier.Algebra;
import eu.solven.matmul.SymbolicVerifier.Result;
import eu.solven.matmul.catalog.SchemeIO;
import tools.jackson.databind.JsonNode;

/**
 * Exercises the per-algebra exact symbolic verifier
 * ({@link SymbolicVerifier}), which checks the trilinear identity in
 * exact BigInteger / Gaussian-rational arithmetic rather than via the
 * legacy floating-point spot check.
 *
 * <p>Canonical cases per CLAUDE.md "Field discipline":</p>
 * <ul>
 *   <li>Strassen ⟨2,2,2⟩=7 verifies over Z;</li>
 *   <li>tampering one Strassen coefficient produces a Z-FAIL;</li>
 *   <li>AlphaTensor-F2 ⟨4,4,4⟩=47 verifies over F2 and FAILS over Z
 *       (because in F2 1+1=0, but in Z 1+1=2);</li>
 *   <li>AlphaEvolve ⟨4,4,4⟩=48 (complex) verifies over C.</li>
 * </ul>
 */
public class TestSymbolicVerifierExact {

	private static final File STRASSEN_222 = eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json");

	private static final File AT_F2_444 = eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section4/alphatensor_F2-4x4x4_m47_a340.json");

	private static final File AE_C_444 = eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section4/alphaevolve-4x4x4_m48_a1264.json");

	@Test
	public void strassen_222_verifies_over_Z() throws Exception {
		JsonNode root = SchemeIO.parseJson(STRASSEN_222);
		Result r = SymbolicVerifier.verify(root);
		assertThat(r.verified()).as("Strassen should verify; got: " + r).isTrue();
		assertThat(r.algebra()).isEqualTo(Algebra.Z);
	}

	@Test
	public void corrupted_strassen_fails_symbolically() throws Exception {
		NonCubicBilinearAlgorithm s = SchemeIO.read(STRASSEN_222);
		double[][] srcU = s.denseU();
		double[][] srcV = s.denseV();
		double[][] srcW = s.denseW();
		double[][] badW = new double[srcW.length][];
		for (int i = 0; i < srcW.length; i++) badW[i] = srcW[i].clone();
		badW[0][0] += 1.0;
		NonCubicBilinearAlgorithm bad = new NonCubicBilinearAlgorithm(s.n, s.m, s.p, srcU, srcV, badW);

		Result r = SymbolicVerifier.verifyBilinear(bad, Algebra.Z);
		assertThat(r.verified()).as("corrupted Strassen must fail; got: " + r).isFalse();
		assertThat(r.algebra()).isEqualTo(Algebra.Z);
	}

	@Test
	public void alphatensor_F2_444_verifies_over_F2() throws Exception {
		JsonNode root = SchemeIO.parseJson(AT_F2_444);
		Result r = SymbolicVerifier.verify(root);
		assertThat(r.verified())
				.as("AlphaTensor-F2 ⟨4,4,4⟩=47 should verify over F2; got: " + r)
				.isTrue();
		assertThat(r.algebra()).isEqualTo(Algebra.F2);
	}

	@Test
	public void alphatensor_F2_444_fails_over_Z() throws Exception {
		// Same algorithm, treated as a Z scheme, fails: 1+1+1 = 3 ≠ 1 in Z
		// where it equals 1 in F2 (XOR). Load-bearing case per CLAUDE.md.
		NonCubicBilinearAlgorithm alg = SchemeIO.read(AT_F2_444);
		Result r = SymbolicVerifier.verifyBilinear(alg, Algebra.Z);
		assertThat(r.verified())
				.as("F2 scheme must NOT verify over Z; got: " + r)
				.isFalse();
	}

	@Test
	public void alphaevolve_C_444_verifies_over_C() throws Exception {
		JsonNode root = SchemeIO.parseJson(AE_C_444);
		Result r = SymbolicVerifier.verify(root);
		assertThat(r.verified())
				.as("AlphaEvolve ⟨4,4,4⟩=48 should verify over C; got: " + r)
				.isTrue();
		assertThat(r.algebra()).isEqualTo(Algebra.C);
	}
}
