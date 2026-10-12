package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.Optional;

import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.LowerBoundRegistry.Bound;
import eu.solven.matmul.catalog.LowerBoundRegistry.Model;
import org.junit.jupiter.api.Test;

/**
 * Unit-tests the two applicability axes of {@link LowerBoundRegistry}: which
 * FIELDS a published floor binds, and which MODEL (commutative vs not).
 */
public class TestLowerBoundRegistry {

	private final LowerBoundRegistry registry = new LowerBoundRegistry();

	@Test
	public void winograd_222_binds_the_commutative_model_too() {
		// The whole point of the 'model' key: Winograd 1971 proves 7 minimal even
		// when commutativity is allowed, so a commutative ⟨2,2,2⟩ claim of 6 is
		// refuted by this entry — which is what the catalog asserted for years.
		Optional<Bound> cmt = registry.binding(2, 2, 2, Field.R, true);
		assertThat(cmt).isPresent();
		assertThat(cmt.get().lb()).isEqualTo(7);
		assertThat(cmt.get().model()).isEqualTo(Model.QUADRATIC);

		assertThat(registry.lowerBound(2, 2, 2, Field.R, false)).isEqualTo(7);
		assertThat(registry.lowerBound(2, 2, 2, Field.F2, true)).isEqualTo(7);
	}

	@Test
	public void bilinear_floors_do_not_bind_commutative_claims() {
		// Bläser's R(⟨3,3,3⟩) ≥ 19 is a non-commutative bilinear floor. Since
		// R_c ≤ R, it says nothing about a commutative-only algorithm.
		assertThat(registry.lowerBound(3, 3, 3, Field.R, false)).isEqualTo(19);
		assertThat(registry.binding(3, 3, 3, Field.R, true)).isEmpty();
	}

	@Test
	public void floors_do_not_transfer_across_characteristics() {
		// Wang 2026's F₂ floors are characteristic-2 only …
		assertThat(registry.lowerBound(3, 4, 4, Field.F2, false)).isEqualTo(29);
		assertThat(registry.binding(3, 4, 4, Field.R, false)).isEmpty();
		// … and a characteristic-0 floor doesn't reach F₂/F₃ either.
		assertThat(registry.lowerBound(4, 4, 4, Field.R, false)).isEqualTo(34);
		assertThat(registry.binding(4, 4, 4, Field.F2, false)).isEmpty();
	}

	@Test
	public void a_floor_binds_subfields_but_not_extensions() {
		// Landsberg–Michałek proved ⟨4,4,4⟩ ≥ 34 over R. A Z or Q scheme is an R
		// scheme, so R_Z ≥ R_Q ≥ R_R ≥ 34 — the floor binds downward …
		assertThat(registry.lowerBound(4, 4, 4, Field.Z, false)).isEqualTo(34);
		assertThat(registry.lowerBound(4, 4, 4, Field.Q, false)).isEqualTo(34);
		// … but NOT upward: rank can drop over C (AlphaEvolve's 48 < 49 is exactly
		// that), so an R floor is no floor for C.
		assertThat(registry.binding(4, 4, 4, Field.C, false)).isEmpty();
	}

	@Test
	public void absent_model_key_defaults_to_bilinear() {
		assertThat(registry.bounds())
				.filteredOn(b -> b.model() == Model.QUADRATIC)
				.hasSize(1);
		assertThat(registry.bounds()).allSatisfy(b -> assertThat(b.lb()).isPositive());
	}

	@Test
	public void unpublished_shapes_report_no_floor() {
		assertThat(registry.binding(19, 23, 29, Field.Q, false)).isEmpty();
		assertThat(registry.lowerBound(19, 23, 29, Field.Q, false)).isZero();
	}
}
