package eu.solven.matmul.additions;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

public class TestSchemeAdditiveComplexity {

	private static NonCubicBilinearAlgorithm load(String name) throws Exception {
		return SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/" + name));
	}

	@Test
	public void naive_matches_VerifierAdditionCount() throws Exception {
		NonCubicBilinearAlgorithm strassen = load("strassen-2x2x2_m7_a18.json");
		SchemeAdditiveComplexity.Result r = SchemeAdditiveComplexity.analyse(strassen);
		// The naive total must equal the existing convention exactly.
		assertThat(r.naive()).isEqualTo(Verifier.additionCount(strassen));
		System.out.println("Strassen ⟨2,2,2⟩=7 : " + r);
	}

	@Test
	public void minimal_never_exceeds_naive() throws Exception {
		for (String f : new String[] {
				"strassen-2x2x2_m7_a18.json",
				"winograd_1971-2x2x2_m7_a24.json",
				"alphatensor_Z-2x2x2_m7_a22.json" }) {
			NonCubicBilinearAlgorithm alg = load(f);
			SchemeAdditiveComplexity.Result r = SchemeAdditiveComplexity.analyse(alg);
			assertThat(r.minimal()).as("minimal ≤ naive for %s", f).isLessThanOrEqualTo(r.naive());
			// The registered SLP must reconstruct the scheme's factor matrices exactly.
			assertThat(r.reconstructs(alg)).as("SLP reconstructs %s", f).isTrue();
			System.out.printf("%-40s %s%n", f, r);
		}
	}

	@Test
	public void cse_reduces_winograd_below_its_naive() throws Exception {
		// Winograd's ⟨2,2,2⟩ is the canonical case where shared sub-expressions
		// matter: its 24-addition naive form collapses toward the famous 15.
		NonCubicBilinearAlgorithm winograd = load("winograd_1971-2x2x2_m7_a24.json");
		SchemeAdditiveComplexity.Result r = SchemeAdditiveComplexity.analyse(winograd);
		assertThat(r.naive()).isEqualTo(24);
		assertThat(r.minimal()).isLessThan(24);
		System.out.println("Winograd ⟨2,2,2⟩=7 : " + r);
	}
}
