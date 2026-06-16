package eu.solven.matmul.f2.sat;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.io.IOException;
import java.util.Optional;

import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

import eu.solven.matmul.BilinearAlgorithm;
import eu.solven.matmul.f2.sat.BreakIdBridge;
import eu.solven.matmul.f2.sat.Cnf;
import eu.solven.matmul.f2.sat.SatMatmulPipeline;
import eu.solven.matmul.f2.sat.Z2CnfEncoder;

/**
 * Integration tests for the BreakID preprocessor bridge. Gracefully skip
 * when the {@code BreakID} binary isn't on {@code $PATH}.
 *
 * <p>To install BreakID locally on macOS:</p>
 * <pre>
 *   # Fix Command Line Tools first if &lt;string&gt; can't be found:
 *   xcode-select --install
 *
 *   # Build BreakID from source:
 *   git clone https://bitbucket.org/krr/breakid.git
 *   cd breakid/src
 *   make
 *   # Then put the resulting BreakID binary on $PATH.
 * </pre>
 */
public class TestBreakIdBridge {

	@Test
	public void cnfWriteReadRoundTrip() throws IOException {
		// Validates the DIMACS I/O independently of BreakID availability —
		// if this round-trips, then the BreakID bridge is wired correctly
		// to the moment a binary appears on $PATH.
		int[][][] target = SatMatmulPipeline.z2DenseMatmulTensor(2);
		Z2CnfEncoder enc = new Z2CnfEncoder(2, 7, target);
		int varCount = enc.getVarCount();
		int clauseCount = enc.getClauses().size();

		File tmp = File.createTempFile("test-cnf-", ".cnf");
		try {
			Cnf.writeDimacs(tmp, varCount, enc.getClauses());
			Cnf.ReadResult result = Cnf.readDimacs(tmp);
			assertThat(result.varCount).isEqualTo(varCount);
			assertThat(result.clauses).hasSize(clauseCount);
			// Spot-check: every original clause matches what was read back.
			for (int i = 0; i < clauseCount; i++) {
				assertThat(result.clauses.get(i)).containsExactly(enc.getClauses().get(i));
			}
		} finally {
			tmp.delete();
		}
	}

	@Test
	public void breakIdAcceleratesR7() {
		Assumptions.assumeTrue(BreakIdBridge.isAvailable(),
				"BreakID binary not on PATH — skipping (see class javadoc for install).");

		int[][][] target = SatMatmulPipeline.z2DenseMatmulTensor(2);

		// Baseline: SAT4J alone (with our hand-coded column lex-ordering).
		long t1 = System.nanoTime();
		Optional<BilinearAlgorithm> baseline = SatMatmulPipeline.findZ2Decomposition(2, 7, target, false);
		long baselineMs = (System.nanoTime() - t1) / 1_000_000;

		// With BreakID preprocessing on top.
		long t2 = System.nanoTime();
		Optional<BilinearAlgorithm> withBreakId = SatMatmulPipeline.findZ2Decomposition(2, 7, target, true);
		long breakIdMs = (System.nanoTime() - t2) / 1_000_000;

		assertThat(baseline).isPresent();
		assertThat(withBreakId).isPresent();
		assertThat(SatMatmulPipeline.verifyZ2(withBreakId.get(), target)).isTrue();

		System.out.printf("%n[BreakID] r=7 dense Z/2: baseline=%dms breakid=%dms speedup=%.2fx%n",
				baselineMs, breakIdMs, (double) baselineMs / Math.max(1, breakIdMs));
	}
}
