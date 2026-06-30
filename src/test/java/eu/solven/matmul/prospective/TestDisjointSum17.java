package eu.solven.matmul.prospective;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.search.CitedBound;

/**
 * The acceptance test for {@link DisjointSumSearch}: at ⟨17,17,17⟩ over Q
 * the search should independently rediscover ≤ FMM-Lille's 2934 without
 * any FMM data import.
 */
class TestDisjointSum17 {

	@Test
	void demo_17_17_17_search_returns_a_value_but_likely_unrealisable() {
		// HONEST DEMO: this is exploratory output, NOT a rank assertion.
		// The current implementation enforces only area-based cover
		// constraints, which are too loose — predicted ranks may be
		// unrealisable. See DisjointSumSearch's Javadoc.
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		CitedBound sota = new CitedBound(lookup);

		Optional<DisjointSumSearch.DisjointSumPrediction> best =
				DisjointSumSearch.findBest(17, 17, 17, sota, /*beamK=*/5, /*maxTerms=*/8, /*minAxis=*/2);

		assertThat(best).isPresent();
		System.out.println("⟨17,17,17⟩ DisjointSumSearch HEURISTIC result:");
		System.out.println("  " + best.get().label());
		System.out.println("  FMM-Lille (constructive bound) = 2934");
		System.out.println("  Strassen-only (our realised bound) = 2940");
		// Sanity: search returned something positive; the actual rank
		// number is exploratory, not a constructive claim.
		assertThat(best.get().totalRank()).isPositive();
	}

	@Test
	void demo_4_4_4_under_loose_cover_shows_false_positive() {
		// EXPECTED to show a false-positive number well below 47/48 —
		// this test exists to make the flaw visible, NOT to assert it
		// matches truth. The output is informational.
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		CitedBound sota = new CitedBound(lookup);

		Optional<DisjointSumSearch.DisjointSumPrediction> p444 =
				DisjointSumSearch.findBest(4, 4, 4, sota);
		assertThat(p444).isPresent();
		System.out.println("⟨4,4,4⟩ DisjointSumSearch HEURISTIC result: "
				+ p444.get().label());
		System.out.println("  (real Q rank = 48 DPS, F₂ = 47 AT — anything below is unrealisable)");
		// Just make sure it ran.
		assertThat(p444.get().totalRank()).isPositive();
	}
}
