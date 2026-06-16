package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

/**
 * Pins the field-aware sub-shape rank lookup that a Kronecker SOTA must use.
 *
 * <p>The catalog holds two {@code ⟨3,3,8⟩} schemes that differ <b>by field</b>:
 * a rank-<b>55</b> one that is F₃/Q-only (½-symmetrisation, {@code 1/8}
 * coefficients) and a rank-<b>56</b> one that is genuinely integer (Z). Computing
 * the rank of {@code ⟨6,6,16⟩ = ⟨2,2,2⟩ ⊗ ⟨3,3,8⟩} therefore depends on the
 * target field:</p>
 *
 * <ul>
 *   <li><b>over Z</b>: the {@code ⟨3,3,8⟩} factor must be the Z-valid <b>56</b>
 *       (the 55 is not valid over Z), so {@code Z⟨6,6,16⟩ = 7·56 = 392};</li>
 *   <li><b>over Q</b> (or F₃): the 55 is valid, so {@code Q⟨6,6,16⟩ = 7·55 = 385}.</li>
 * </ul>
 *
 * <p>The bug this guards: a field-blind lookup picks the lowest rank (55)
 * regardless of target field, producing a {@code 385} scheme that then gets
 * tagged Z — an over-claim, since no Z-valid ⟨3,3,8⟩ of rank 55 exists. The
 * genuine {@code Z⟨6,6,16⟩} SOTA is Perminov's 392. (Companion to
 * {@link TestKronStubFieldConsistency}, which catches the over-claim in the
 * committed catalog; this one pins the lookup that should prevent producing it.)</p>
 */
public class TestKronFieldAwareSota {

	@Test
	public void z_lookup_picks_the_integer_3x3x8_56_not_the_rational_55() {
		FieldAwareLookup z = new FieldAwareLookup("Z");
		// The 55-rank ⟨3,3,8⟩ is F₃/Q-only; over Z the lookup must skip it.
		assertThat(z.findRank(3, 3, 8))
				.as("Z⟨3,3,8⟩ must be the integer 56, not the F₃/Q-only 55")
				.isEqualTo(56);
		assertThat(z.findRank(2, 2, 2))
				.as("Z⟨2,2,2⟩ = Strassen 7")
				.isEqualTo(7);
	}

	@Test
	public void q_lookup_may_use_the_rational_3x3x8_55() {
		FieldAwareLookup q = new FieldAwareLookup("Q");
		assertThat(q.findRank(3, 3, 8))
				.as("Q⟨3,3,8⟩ may use the rational 55")
				.isEqualTo(55);
	}

	@Test
	public void kron_sota_for_6x6x16_differs_by_field() {
		FieldAwareLookup z = new FieldAwareLookup("Z");
		FieldAwareLookup q = new FieldAwareLookup("Q");
		// ⟨6,6,16⟩ = ⟨2,2,2⟩ ⊗ ⟨3,3,8⟩: rank is the product of the field-correct
		// sub-shape ranks. Over Z that is 7·56 = 392 (= Perminov's Z SOTA); over Q
		// it is 7·55 = 385. A field-blind lookup would wrongly yield 385 for Z.
		int zKron = z.findRank(2, 2, 2) * z.findRank(3, 3, 8);
		int qKron = q.findRank(2, 2, 2) * q.findRank(3, 3, 8);
		assertThat(zKron).as("Z⟨6,6,16⟩ via Kron = 7·56").isEqualTo(392);
		assertThat(qKron).as("Q⟨6,6,16⟩ via Kron = 7·55").isEqualTo(385);
		assertThat(zKron).as("the Z Kron must NOT collapse onto the F₃/Q 385").isGreaterThan(qKron);
	}
}
