package eu.solven.matmul.docs.migrate;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;

import org.junit.jupiter.api.Test;

/**
 * Guards the 2026-07-07 silent bug: {@link NarrowFields#expandFromBase} returned a
 * hardcoded containment chain for the strongest base tag, so a Q scheme carrying a
 * VERIFIED F3 membership (e.g. [F3,Q,R,C] — all denominators powers of 2, coprime
 * to 3) mechanically "expanded" to [Q,R,C], stripping F3 from 3.7k schemes above
 * the verify cap. No crash — just quietly narrower field sets in catalog.json.
 * Mechanical expansion must be a union: it may only ADD tags, never remove one.
 */
public class TestNarrowFields {

	@Test
	public void expansion_never_drops_a_verified_prime_field_tag() {
		// The exact shape of the regression: ½-denominator Q scheme, F3-verified.
		assertThat(NarrowFields.expandFromBase(List.of("F3", "Q", "R", "C")))
				.containsExactly("F3", "Q", "R", "C");
		// Symmetric case: ⅓-denominator Q scheme, F2-verified.
		assertThat(NarrowFields.expandFromBase(List.of("F2", "Q", "R", "C")))
				.containsExactly("F2", "Q", "R", "C");
	}

	@Test
	public void expansion_still_widens_along_the_containment_chain() {
		assertThat(NarrowFields.expandFromBase(List.of("Z")))
				.containsExactly("F2", "F3", "Z", "Q", "R", "C");
		assertThat(NarrowFields.expandFromBase(List.of("Q")))
				.containsExactly("Q", "R", "C");
		assertThat(NarrowFields.expandFromBase(List.of("R")))
				.containsExactly("R", "C");
	}

	@Test
	public void prime_fields_have_no_char0_inclusion() {
		assertThat(NarrowFields.expandFromBase(List.of("F2"))).containsExactly("F2");
		assertThat(NarrowFields.expandFromBase(List.of("F3"))).containsExactly("F3");
	}

	@Test
	public void expansion_is_idempotent_on_already_full_sets() {
		List<String> full = List.of("F2", "F3", "Z", "Q", "R", "C");
		assertThat(NarrowFields.expandFromBase(full)).containsExactlyElementsOf(full);
	}
}
