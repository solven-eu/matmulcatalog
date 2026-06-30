package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.search.ConstructiveMethod.Prediction;

/**
 * The LITA {@link ConstructiveMethod}: cubic-only, n≥19 domain, a replayable
 * {@code TA_lita(n=N)} lineage, and registry membership so the search treats it
 * as an alternative-rank candidate alongside TA_dis.
 */
class TestLitaTrilinearAggregationMethod {

	private final LitaTrilinearAggregationMethod method = new LitaTrilinearAggregationMethod();

	@Test
	void predicts_cubic_at_or_above_19_with_replayable_lineage() {
		Optional<Prediction> p21 = method.predict(21, 21, 21, null);
		assertThat(p21).isPresent();
		assertThat(p21.get().predictedRank()).isEqualTo(5198L);
		assertThat(p21.get().lineageCompact()).isEqualTo("TA_lita(n=21)");
		assertThat(p21.get().verified()).isTrue();
	}

	@Test
	void does_not_apply_off_domain() {
		assertThat(method.predict(18, 18, 18, null)).isEmpty(); // below n=19
		assertThat(method.predict(21, 21, 22, null)).isEmpty(); // non-cubic
	}

	@Test
	void registered_in_method_catalog() {
		assertThat(MethodCatalog.all())
				.anyMatch(m -> m instanceof LitaTrilinearAggregationMethod);
	}

	@Test
	void construct_builds_an_exact_scheme_at_the_predicted_rank() {
		// Odd N=19 (sparse) — exact-verifiable and fast.
		Optional<NonCubicBilinearAlgorithm> alg = method.construct(19, 19, 19, null);
		assertThat(alg).isPresent();
		assertThat((long) alg.get().r).isEqualTo(4016L);
		assertThat(eu.solven.matmul.verifiers.Verifier.isExactNonCubic(alg.get())).isTrue();
	}
}
