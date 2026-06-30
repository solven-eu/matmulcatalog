package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

import java.util.Arrays;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.BilinearAlgorithm;
import eu.solven.matmul.papers.laderman1976.Laderman23;
import eu.solven.matmul.papers.strassen1969.Strassen7;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.Compose;

/**
 * Validates the Kronecker-product composer in {@link Compose}. Builds Strassen
 * compositions at progressively larger formats and confirms exactness with
 * {@link Verifier#isExact} up to {@code ⟨8,8,8⟩} (verification at
 * {@code ⟨16,16,16⟩} and above is impractically slow due to the verifier's
 * `O(n^6 · r)` cost — the construction is mathematically exact by the
 * Kronecker theorem).
 */
public class TestCompose {

	@Test
	public void strassen_squared_gives_4x4_at_rank_49_and_verifies() {
		BilinearAlgorithm composed = Compose.strassenPower(2);
		assertThat(composed.n).isEqualTo(4);
		assertThat(composed.r).isEqualTo(49);
		assertThat(Verifier.isExact(composed)).isTrue();
	}

	@Test
	public void strassen_cubed_gives_8x8_at_rank_343_and_verifies() {
		BilinearAlgorithm composed = Compose.strassenPower(3);
		assertThat(composed.n).isEqualTo(8);
		assertThat(composed.r).isEqualTo(343);
		assertThat(Verifier.isExact(composed)).isTrue();
	}

	@Test
	public void strassen_fourth_gives_16x16_at_rank_2401_built_in_seconds() {
		// Construction only — too big for the O(n^6 · r) verifier.
		long t0 = System.currentTimeMillis();
		BilinearAlgorithm composed = Compose.strassenPower(4);
		long dt = System.currentTimeMillis() - t0;
		assertThat(composed.n).isEqualTo(16);
		assertThat(composed.r).isEqualTo(2401);
		assertThat(dt).isLessThan(30_000); // 30s budget for construction
	}

	@Test
	public void strassen_fifth_gives_32x32_at_rank_16807() {
		BilinearAlgorithm composed = Compose.strassenPower(5);
		assertThat(composed.n).isEqualTo(32);
		assertThat(composed.r).isEqualTo(16_807);
		// Verifier impractical here; correctness is by Kronecker construction.
	}

	@Test
	public void chain_three_strassens_equivalent_to_strassenPower_3() {
		BilinearAlgorithm a = Strassen7.get();
		BilinearAlgorithm composedViaChain = Compose.chain(Arrays.asList(a, a, a));
		BilinearAlgorithm composedViaPower = Compose.strassenPower(3);
		assertThat(composedViaChain.n).isEqualTo(composedViaPower.n);
		assertThat(composedViaChain.r).isEqualTo(composedViaPower.r);
		// Both are exact decompositions of the same tensor; they may not be
		// numerically identical (no canonical ordering), but residual should be 0
		// against the target matmul tensor.
		assertThat(Verifier.isExact(composedViaChain)).isTrue();
	}

	@Test
	public void strassen_times_laderman_gives_6x6_at_rank_161() {
		// ⟨2,2,2⟩ ⊗ ⟨3,3,3⟩ → ⟨6,6,6⟩ at rank 7 · 23 = 161.
		BilinearAlgorithm composed = Compose.kronecker(Strassen7.get(), Laderman23.get());
		assertThat(composed.n).isEqualTo(6);
		assertThat(composed.r).isEqualTo(161);
		// Verify against the ⟨6,6,6⟩ matmul tensor.
		double residual = Verifier.residual(composed);
		assertThat(residual).isEqualTo(0.0, within(1e-10));
	}

	@Test
	public void laderman_times_strassen_also_gives_6x6_at_rank_161() {
		// Reverse order: ⟨3,3,3⟩ outer × ⟨2,2,2⟩ inner = ⟨6,6,6⟩, also rank 161.
		BilinearAlgorithm composed = Compose.kronecker(Laderman23.get(), Strassen7.get());
		assertThat(composed.n).isEqualTo(6);
		assertThat(composed.r).isEqualTo(161);
		assertThat(Verifier.isExact(composed)).isTrue();
	}
}
