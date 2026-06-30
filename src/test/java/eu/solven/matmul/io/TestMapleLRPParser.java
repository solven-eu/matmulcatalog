package eu.solven.matmul.io;

import org.junit.jupiter.api.Tag;

import org.junit.jupiter.api.Disabled;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.condition.DisabledIfSystemProperty;
import org.junit.jupiter.api.condition.EnabledIf;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.SymbolicVerifier;
import eu.solven.matmul.verifiers.Verifier;

/**
 * Sanity tests for {@link MapleLRPParser}: parses the FMM-Lille
 * {@code 17x17x17_LRP.mpl} file, asserts the rank announced by the
 * file (2931 — three less than the publicly cited 2934 — see the
 * import report for the gap analysis), and runs the
 * {@link SymbolicVerifier} as a hard sanity check that the
 * {@code (U, V, W)} we produced really do compute {@code ⟨17,17,17⟩}
 * exactly over {@code Q}.
 *
 * <p>The LRP file lives outside the git tree (≈5 MB Maple expression);
 * the tests are skipped if the file is not present.</p>
 */
@Tag("slow")
public class TestMapleLRPParser {

	private static final File LRP = new File("references/fmm-lille/17x17x17/17x17x17_LRP.mpl");

	static boolean lrpAvailable() { return LRP.isFile(); }

	@Test
	@EnabledIf("lrpAvailable")
	public void parses_dimensions() throws Exception {
		MapleLRPParser.LRPMatrices mats = MapleLRPParser.parseMatrices(LRP);
		assertThat(mats.rank()).isEqualTo(2931);
		assertThat(mats.dimA()).isEqualTo(17 * 17);
		assertThat(mats.dimB()).isEqualTo(17 * 17);
		assertThat(mats.dimC()).isEqualTo(17 * 17);
	}

	@Test
	@EnabledIf("lrpAvailable")
	public void parses_to_algorithm_with_rank_2931() throws Exception {
		NonCubicBilinearAlgorithm alg = MapleLRPParser.parse(LRP, 17, 17, 17);
		assertThat(alg.n).isEqualTo(17);
		assertThat(alg.m).isEqualTo(17);
		assertThat(alg.p).isEqualTo(17);
		assertThat(alg.r).isEqualTo(2931);
	}

	@Test
	@EnabledIf("lrpAvailable")
	public void numeric_spot_check_passes() throws Exception {
		NonCubicBilinearAlgorithm alg = MapleLRPParser.parse(LRP, 17, 17, 17);
		// Numeric spot-check (fast). If this fails, the parser produced
		// wrong U/V/W matrices and the symbolic verifier will also fail.
		boolean ok = Verifier.passesRandomMatmulSpotCheck(alg);
		assertThat(ok).as("numeric spot check on parsed ⟨17,17,17⟩ scheme").isTrue();
	}

	@Test
	@EnabledIf("lrpAvailable")
	@DisabledIfSystemProperty(named = "skip.slow", matches = "true")
	public void symbolic_verifier_passes_on_emitted_json() throws Exception {
		// Parse → write → SymbolicVerifier.verify(File). This mirrors the
		// production import path (ImportFmmLille17LRP) end-to-end.
		NonCubicBilinearAlgorithm alg = MapleLRPParser.parse(LRP, 17, 17, 17);
		int adds = Verifier.additionCount(alg);
		File out = new File("target/test-output/fmm-lille-LRP-r" + alg.r + "-a" + adds + ".json");
		out.getParentFile().mkdirs();
		ImportFmmLille17LRP.writeWithMetadata(alg, out, adds);

		SymbolicVerifier.Result r = SymbolicVerifier.verify(out);
		System.out.println("LRP-parsed ⟨17,17,17⟩: verified=" + r.verified()
				+ ", algebra=" + r.algebra() + ", reason=" + r.reason());
		assertThat(r.verified())
				.as("LRP-parsed ⟨17,17,17⟩ should verify symbolically over Q")
				.isTrue();
		assertThat(r.algebra()).isEqualTo(SymbolicVerifier.Algebra.Q);
	}

	@Test
	@EnabledIf("lrpAvailable")
	public void counts_duplicate_UV_pairs() throws Exception {
		// Sanity: the LRP file (kin-row-reduced) should have NO duplicate
		// (U[:,k], V[:,k]) products — reduction would have merged them.
		// Report the count for the investigation log.
		NonCubicBilinearAlgorithm alg = MapleLRPParser.parse(LRP, 17, 17, 17);
		int dupPairs = MapleLRPParser.countDuplicateUVPairs(alg);
		System.out.println("LRP-parsed duplicate (U,V) pairs: " + dupPairs);
		assertThat(dupPairs).as("kin-row-reduced LRP should have no duplicate products")
				.isEqualTo(0);
	}
}
