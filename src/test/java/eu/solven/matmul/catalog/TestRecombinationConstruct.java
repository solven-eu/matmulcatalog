package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.papers.strassen1969.Strassen7;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.Compositions;
import eu.solven.matmul.catalog.Recombination;
import eu.solven.matmul.catalog.Recombination.AlgorithmLookup;

/**
 * Validates the constructive recombination: given a base algorithm + a lookup
 * for sub-algorithms, build an actual {@link NonCubicBilinearAlgorithm} for
 * the target format and verify it.
 */
public class TestRecombinationConstruct {

	private static final NonCubicBilinearAlgorithm STRASSEN =
			NonCubicBilinearAlgorithm.fromCubic(Strassen7.get());

	/**
	 * Lookup that returns a fixed mapping from {@code (n,m,p)} to verified
	 * sub-algorithms. Used to drive the construction over real-valued schemes.
	 */
	private static AlgorithmLookup mapLookup(Map<String, NonCubicBilinearAlgorithm> map) {
		return (n, m, p) -> Optional.ofNullable(map.get(n + "," + m + "," + p));
	}

	@Test
	public void construct_222_via_111_lookup_recovers_strassen_rank_7() throws IOException {
		// Target ⟨2,2,2⟩ recombined via Strassen base; sub-problems all become ⟨1,1,1⟩.
		Map<String, NonCubicBilinearAlgorithm> map = new HashMap<>();
		map.put("1,1,1", trivial111());
		NonCubicBilinearAlgorithm composed =
				Recombination.construct(2, 2, 2, STRASSEN, mapLookup(map));
		assertThat(composed.n).isEqualTo(2);
		assertThat(composed.r).isEqualTo(7);
		assertThat(Verifier.isExactNonCubic(composed)).isTrue();
	}

	@Test
	public void construct_444_via_strassen_with_222_lookup_yields_49() throws IOException {
		// Target ⟨4,4,4⟩ via Strassen base: each base mult sees sub-format ⟨2,2,2⟩.
		// With Strassen as the sub-lookup, total rank = 7·7 = 49.
		Map<String, NonCubicBilinearAlgorithm> map = new HashMap<>();
		map.put("2,2,2", STRASSEN);
		NonCubicBilinearAlgorithm composed =
				Recombination.construct(4, 4, 4, STRASSEN, mapLookup(map));
		assertThat(composed.n).isEqualTo(4);
		assertThat(composed.r).isEqualTo(49);
		assertThat(Verifier.isExactNonCubic(composed)).isTrue();
	}

	@Test
	public void construct_666_via_strassen_with_333_laderman_yields_161() throws IOException {
		// Target ⟨6,6,6⟩ via Strassen base: each base mult sees sub-format ⟨3,3,3⟩.
		// With Laderman lookup: total rank = 7·23 = 161.
		NonCubicBilinearAlgorithm laderman = Compositions.loadScheme("laderman_1976-3x3x3_m23.json");
		Map<String, NonCubicBilinearAlgorithm> map = new HashMap<>();
		map.put("3,3,3", laderman);
		NonCubicBilinearAlgorithm composed =
				Recombination.construct(6, 6, 6, STRASSEN, mapLookup(map));
		assertThat(composed.n).isEqualTo(6);
		assertThat(composed.r).isEqualTo(161);
		assertThat(Verifier.isExactNonCubic(composed)).isTrue();
	}

	@Test
	public void construct_888_via_strassen_with_444_lookup_yields_recursive_strassen() throws IOException {
		// Target ⟨8,8,8⟩ via Strassen base; each base mult sees ⟨4,4,4⟩.
		// With Strassen² (=49) lookup: total = 7·49 = 343 (matches Strassen³).
		NonCubicBilinearAlgorithm strassenSquared = Compose.kroneckerGeneral(STRASSEN, STRASSEN);
		Map<String, NonCubicBilinearAlgorithm> map = new HashMap<>();
		map.put("4,4,4", strassenSquared);
		NonCubicBilinearAlgorithm composed =
				Recombination.construct(8, 8, 8, STRASSEN, mapLookup(map));
		assertThat(composed.n).isEqualTo(8);
		assertThat(composed.r).isEqualTo(343);
		assertThat(Verifier.isExactNonCubic(composed)).isTrue();
	}

	@Test
	public void missing_sub_algorithm_throws() {
		// No lookup entries → construct should fail.
		Map<String, NonCubicBilinearAlgorithm> map = new HashMap<>();
		try {
			Recombination.construct(2, 2, 2, STRASSEN, mapLookup(map));
		} catch (IllegalStateException expected) {
			assertThat(expected.getMessage()).contains("missing sub-algorithm");
			return;
		}
		assertThat(true).as("should have thrown").isFalse();
	}

	/** Trivial scalar matmul {@code ⟨1,1,1⟩}: rank 1, one multiplication. */
	private static NonCubicBilinearAlgorithm trivial111() {
		double[][] U = { { 1.0 } };
		double[][] V = { { 1.0 } };
		double[][] W = { { 1.0 } };
		return new NonCubicBilinearAlgorithm(1, 1, 1, U, V, W);
	}
}
