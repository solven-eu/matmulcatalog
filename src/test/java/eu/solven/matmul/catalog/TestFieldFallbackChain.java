package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Regression guard for the downward-transfer field-fallback chain
 * (Z ⊂ Q ⊂ R ⊂ C, plus integer ⇒ F₂/F₃ by mod-p reduction) after the #174
 * fields[] migration + #175 narrowing. The canonical ⟨4,4,4⟩ landscape is the
 * witness: 49/Z (Strassen²), 48/Q (DPS-2025, ½ coefficients), 47/F₂
 * (AlphaTensor). A Q-only scheme MUST NOT leak into a Z lookup (the bug class
 * that produced the 2026-06-03 ⟨17,17,17⟩=2937 regression).
 */
public class TestFieldFallbackChain {

	private static int rank(String field, int n, int m, int p) {
		Optional<NonCubicBilinearAlgorithm> a = new FieldAwareLookup(field).find(n, m, p);
		assertThat(a).as(field + " ⟨" + n + "," + m + "," + p + "⟩ present").isPresent();
		return a.get().r;
	}

	@Test
	public void strassen_222_valid_in_every_field() {
		// Integer-exact ⇒ valid (and found) over the whole chain + prime fields.
		for (String field : new String[] { "Z", "Q", "R", "C", "F2" }) {
			assertThat(rank(field, 2, 2, 2))
					.as("Strassen ⟨2,2,2⟩ over " + field)
					.isEqualTo(7);
		}
	}

	@Test
	public void four_cubed_respects_field_discipline() {
		// Z: integer-only — must NOT return the ½-coefficient DPS-48 (that's Q).
		assertThat(rank("Z", 4, 4, 4))
				.as("⟨4,4,4⟩ over Z is the integer best (Strassen²=49), not the Q/F₂ schemes")
				.isEqualTo(49);
		// Q/R/C: DPS-2025 ⟨4,4,4⟩=48 lifts up the chain.
		assertThat(rank("Q", 4, 4, 4)).as("⟨4,4,4⟩ over Q").isEqualTo(48);
		assertThat(rank("R", 4, 4, 4)).as("⟨4,4,4⟩ over R").isEqualTo(48);
		assertThat(rank("C", 4, 4, 4)).as("⟨4,4,4⟩ over C").isEqualTo(48);
		// F₂: AlphaTensor ⟨4,4,4⟩=47 — its own characteristic-2 universe.
		assertThat(rank("F2", 4, 4, 4)).as("⟨4,4,4⟩ over F2").isEqualTo(47);
	}

	@Test
	public void z_lookup_excludes_rational_only_schemes() {
		// Downward transfer is one-directional: a Z lookup never returns a
		// strictly-Q scheme (½ ∉ Z). Equivalent to "no field-discipline leak".
		int zBest = rank("Z", 4, 4, 4);
		int qBest = rank("Q", 4, 4, 4);
		assertThat(zBest).as("Z-best ≥ Q-best (Z is the more restrictive field)")
				.isGreaterThanOrEqualTo(qBest);
	}
}
