package eu.solven.matmul.docs.migrate;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;

/**
 * Guards the flip-graph {@code .exp} parser against the two silent-drop hazards
 * that distinguish flip-graph files from the {@code ±1}-only Kauers-2026 form:
 * integer coefficients ({@code 2*a11}) and un-parenthesised single-variable
 * factors ({@code (-a21)*(-b11+b13)*c12}). Either bug would parse a non-matmul
 * tensor that fails verification — so "parses to ⟨2,2,3⟩ r=11 AND verifies" is a
 * single end-to-end assertion covering both.
 */
public class TestImportMoosbauerFlips {

	/** Verbatim {@code solutions/223-11-mod0.exp} from jakobmoosbauer/flips.
	 *  Has {@code 2*}/{@code -2*} coefficients AND bare {@code c12}/{@code a12}/
	 *  {@code b21}/{@code c11} factors. */
	private static final String FLIPS_223_11_MOD0 = String.join("\n",
			"(-a21)*(-b11+b13)*c12",
			"a12*b21*(c11+c12)",
			"(-a11+a12+a21-a22)*b12*(c21+2*c22)",
			"(a12-a22)*(-b21+b23)*c12",
			"(-2*a11+2*a12+a21-a22)*(b22+b23)*(c21+c22)",
			"a11*b13*(-c11+2*c12+c31+2*c32)",
			"a12*b23*(-c21-c22+c31+c32)",
			"(a12-a22)*(b12+b22)*(c12-c22+c32)",
			"(-2*a11+a21)*(b12+b13+b22+b23)*(c12+c32)",
			"(2*a11-a12-a21+a22)*(b12+b22+b23)*(c12+c21+c22+c32)",
			"a11*(b11+b13)*c11");

	@Test
	public void parses_integer_coefficients_and_bare_factors() {
		NonCubicBilinearAlgorithm alg = ImportMoosbauerFlips.parse(FLIPS_223_11_MOD0);

		assertThat(alg).isNotNull();
		// true (unsorted) orientation recovered from index ranges, not the sorted name
		assertThat(alg.n).isEqualTo(2);
		assertThat(alg.m).isEqualTo(2);
		assertThat(alg.p).isEqualTo(3);
		assertThat(alg.r).isEqualTo(11);

		// the payoff: it really computes ⟨2,2,3⟩ matmul. Fails if the 2*/-2*
		// coefficient is dropped OR a bare single-variable factor is dropped.
		assertThat(Verifier.isExactNonCubic(alg)).isTrue();
		// integer ⇒ also exact over F₂ (exercises the mod2 verify path).
		assertThat(Verifier.isExactNonCubicF2(alg)).isTrue();
	}

	@Test
	public void handles_outer_negated_whole_product() {
		// 555-97-mod0.exp line 43 wraps the entire rank-one term in a negation:
		// `-(a-form * b-form * c-form)`. A depth-0 split sees no top-level '*', so a
		// naive parser collapses it to a single factor and silently drops the term
		// (97→96 → fails to verify). The peel-and-fold-sign path keeps all 3 terms.
		String body = String.join("\n",
				"a11*b11*c11",
				"-(a12*b21*c11)",                  // outer-negated whole product
				"a12*b21*(c11+c12)");
		NonCubicBilinearAlgorithm alg = ImportMoosbauerFlips.parse(body);
		assertThat(alg).isNotNull();
		assertThat(alg.r).isEqualTo(3);  // all three terms kept, none dropped

		// the negated term's a-form coefficient must carry the −1.
		// product 1 (k=0): +a12·b21·c11 ; the a12 entry (i=1,j=2 → U row 0*m+1) is +1 at k=2,
		// and −1 at the negated k=1.
		assertThat(alg.denseU()[1][1]).isEqualTo(-1d);
		assertThat(alg.denseU()[1][2]).isEqualTo(1d);
	}

	@Test
	public void dropping_the_coefficient_would_break_verification() {
		// Sanity that the coefficient is load-bearing: the SAME scheme with every
		// `2*` weakened to `1*` must NOT verify — proving the assertion above is not
		// vacuously passing on a coefficient-blind parse.
		String coeffStripped = FLIPS_223_11_MOD0.replace("2*", "");
		NonCubicBilinearAlgorithm broken = ImportMoosbauerFlips.parse(coeffStripped);
		assertThat(Verifier.isExactNonCubic(broken)).isFalse();
	}
}
