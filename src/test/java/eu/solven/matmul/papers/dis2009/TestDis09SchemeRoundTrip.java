package eu.solven.matmul.papers.dis2009;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Round-trip smoke test: build a small Pan TA scheme, write it to a temp
 * file via {@link SchemeIO}, read it back, and verify the recovered
 * algorithm still passes the matmul-tensor exactness check. Guards against
 * encoder/decoder bugs that the
 * {@code TestPanTrilinearAggregationBuild} unit tests (in-memory only)
 * would miss.
 */
class TestDis09SchemeRoundTrip {

	@Test
	void roundTripN4() throws Exception {
		NonCubicBilinearAlgorithm built = PanTrilinearAggregation.build(4);
		assertThat(Verifier.isExactNonCubic(built)).isTrue();

		Path tmp = Files.createTempFile("dis09-roundtrip-", ".json");
		try {
			Lineage.Node lineage = new Lineage.Atom("DIS09Lemma4(n=4)");
			SchemeIO.write(built, tmp.toFile(), lineage);

			File f = tmp.toFile();
			assertThat(f).exists();
			assertThat(f.length()).isGreaterThan(1000);

			NonCubicBilinearAlgorithm reloaded = SchemeIO.read(f);
			assertThat(reloaded.n).isEqualTo(4);
			assertThat(reloaded.r).isEqualTo(built.r);
			assertThat(Verifier.isExactNonCubic(reloaded)).isTrue();
		} finally {
			Files.deleteIfExists(tmp);
		}
	}
}
