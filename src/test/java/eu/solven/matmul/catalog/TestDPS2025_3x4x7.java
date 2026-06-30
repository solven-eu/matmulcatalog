package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.io.IOException;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;

/**
 * Verifies the Dumas-Pernet-Sedoglavic 2025 scheme for {@code Q⟨3,4,7⟩:m=63},
 * the rational-coefficient counterpart of AlphaEvolve's complex
 * {@code C⟨3,4,7⟩:m=63} from May 2025.
 *
 * <p>Source: HAL preprint hal-05121550 (June 2025).
 * Maple factor matrices mirrored at
 * {@code https://fmm.univ-lille.fr/3x4x7_tensor.mpl.bz2}.</p>
 */
public class TestDPS2025_3x4x7 {

	@Test
	public void dps_2025_3x4x7_m63_verifies_exact_non_cubic() throws IOException {
		// Resolve by content identity (shape + source note), not a pinned path —
		// the 2026-06 migration made filenames content-driven and they re-hash.
		File f = SchemeResolver.byHint("dumas_pernet_sedoglavic_2025-3x4x7_m63.json");
		assertThat(f).as("DPS-2025 ⟨3,4,7⟩=63 scheme JSON present").exists();
		NonCubicBilinearAlgorithm alg = SchemeIO.read(f);
		assertThat(alg.n).isEqualTo(3);
		assertThat(alg.m).isEqualTo(4);
		assertThat(alg.p).isEqualTo(7);
		assertThat(alg.r).isEqualTo(63);
		// Full residual check across all input entries of A ∈ Q^{3x4} and B ∈ Q^{4x7}.
		assertThat(Verifier.isExactNonCubic(alg))
				.as("DPS-2025 ⟨3,4,7⟩=63 over Q must exactly compute matmul")
				.isTrue();
	}
}
