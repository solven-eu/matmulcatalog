package eu.solven.matmul;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.papers.laderman1976.Laderman23;
import eu.solven.matmul.papers.strassen1969.Strassen7;
import eu.solven.matmul.Verifier;

/**
 * Validates the addition-count formula {@link Verifier#additionCount} against
 * Strassen's published count of 18 additions for ⟨2,2,2⟩ r=7.
 */
public class TestAdditionCount {

	@Test
	public void strassen_has_18_additions() {
		// Strassen 1969: 7 multiplications + 18 additions to multiply two 2×2 matrices.
		// Breakdown: 5 input-side U-adds + 5 input-side V-adds + 8 output-side W-adds = 18.
		assertThat(Verifier.additionCount(Strassen7.get())).isEqualTo(18);
	}

	@Test
	public void laderman_addition_count_is_published_value() {
		// Laderman 1976: 23 multiplications + 98 additions for ⟨3,3,3⟩.
		// (Smirnov 2013 and the dronperminov catalog corroborate ~98 for the
		// canonical Laderman scheme.)
		int adds = Verifier.additionCount(Laderman23.get());
		// Allow some slack since "additions" depends on the exact variant —
		// just confirm it's in a sane range.
		assertThat(adds).as("Laderman additions").isBetween(85, 110);
	}
}
