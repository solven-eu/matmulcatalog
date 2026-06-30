package eu.solven.matmul;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.List;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.papers.laderman1976.Laderman23;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Validates {@link Verifier#symbolicDiscrepancies}: returns empty for
 * correct schemes, returns precise pinpoint diagnostics for broken
 * ones.
 */
public class TestSymbolicVerifier {

	@Test
	public void strassen_has_no_symbolic_discrepancies() throws Exception {
		NonCubicBilinearAlgorithm s = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		List<Verifier.SymbolicDiff> diffs = Verifier.symbolicDiscrepancies(s);
		assertThat(diffs).as("Strassen 1969 must have zero symbolic discrepancies").isEmpty();
	}

	@Test
	public void laderman_has_no_symbolic_discrepancies() {
		NonCubicBilinearAlgorithm l = NonCubicBilinearAlgorithm.fromCubic(Laderman23.get());
		List<Verifier.SymbolicDiff> diffs = Verifier.symbolicDiscrepancies(l);
		assertThat(diffs).as("Laderman 1976 must have zero symbolic discrepancies").isEmpty();
	}

	@Test
	public void corrupted_strassen_emits_precise_diagnostic() throws Exception {
		NonCubicBilinearAlgorithm s = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		double[][] srcU = s.denseU();
		double[][] srcV = s.denseV();
		double[][] srcW = s.denseW();
		double[][] badW = new double[srcW.length][];
		for (int i = 0; i < srcW.length; i++) badW[i] = srcW[i].clone();
		// Flip a single coefficient: this will skew specific bilinear terms.
		badW[0][0] += 1.0;
		NonCubicBilinearAlgorithm bad = new NonCubicBilinearAlgorithm(s.n, s.m, s.p, srcU, srcV, badW);

		List<Verifier.SymbolicDiff> diffs = Verifier.symbolicDiscrepancies(bad);
		assertThat(diffs).as("corrupted Strassen must produce at least one discrepancy").isNotEmpty();
		System.out.println("First discrepancies in corrupted Strassen:");
		for (Verifier.SymbolicDiff d : diffs.subList(0, Math.min(5, diffs.size()))) {
			System.out.println("  " + d);
		}
	}
}
