package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;
import java.util.stream.Stream;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Tag;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Sampled verification for fmm-lille schemes imported via
 * {@code tools/import_fmm_maple.py}. The full residual at max-dim ≥ 12 is
 * minutes-to-hours; the 50k-sample check takes &lt;1s and is conclusive in
 * practice for Kronecker-derived schemes.
 */
@Tag("catalog-iterating")
public class TestFmmMapleImport {

	@Test
	public void fmm_lille_16x16x16_m2304_passes_sampled_verification() throws IOException {
		File f = findFmmScheme(16, 16, 16, 2304);
		if (f == null) return; // Skip when not yet imported (test runs in CI without --details).
		NonCubicBilinearAlgorithm alg = SchemeIO.read(f);
		assertThat(alg.n).isEqualTo(16);
		assertThat(alg.r).isEqualTo(2304);
		int wrong = Verifier.residualSampled(alg, 50_000, 0xCAFEBABE);
		assertThat(wrong).as("sampled residual on fmm-lille ⟨16,16,16⟩ r=2304").isEqualTo(0);
	}

	@Test
	public void fmm_lille_2x2x2_m7_matches_strassen_verifies_exact() throws IOException {
		File f = findFmmScheme(2, 2, 2, 7);
		if (f == null) return;
		NonCubicBilinearAlgorithm alg = SchemeIO.read(f);
		assertThat(alg.n).isEqualTo(2);
		assertThat(alg.r).isEqualTo(7);
		// Full residual check (fast at ⟨2,2,2⟩).
		assertThat(Verifier.isExactNonCubic(alg)).isTrue();
	}

	/** Locate `fmm-lille_{n}x{m}x{p}_r{rank}_a*.json` under the schemes tree. */
	private static File findFmmScheme(int n, int m, int p, int rank) throws IOException {
		String prefix = String.format("fmm-lille_%dx%dx%d_r%d_", n, m, p, rank);
		Path root = Path.of("src/main/resources/schemes");
		try (Stream<Path> s = Files.walk(root)) {
			Optional<Path> hit = s.filter(pp -> pp.getFileName().toString().startsWith(prefix))
					.findFirst();
			return hit.map(Path::toFile).orElse(null);
		}
	}
}
