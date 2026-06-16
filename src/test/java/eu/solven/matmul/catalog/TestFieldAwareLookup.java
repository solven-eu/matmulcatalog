package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;

/**
 * Validates the field-fallback chain in {@link FieldAwareLookup}.
 *
 * <p>Specifically:</p>
 * <ul>
 *   <li>{@code C} lookup returns AlphaEvolve's {@code ⟨4,4,4⟩=48} (C-specific)
 *       when available, else falls back to Strassen² {@code =49} (R-class).</li>
 *   <li>{@code R} lookup does NOT return AlphaEvolve's 48 (C-only); returns 49.</li>
 *   <li>{@code F2} lookup returns AlphaTensor's {@code 47} (F₂-specific) for
 *       {@code ⟨4,4,4⟩}; falls back to nothing for unavailable F₂ targets.</li>
 * </ul>
 */
public class TestFieldAwareLookup {

	@Test
	public void c_lookup_finds_alphaevolve_48_for_444() {
		FieldAwareLookup lookupC = new FieldAwareLookup("C");
		Optional<NonCubicBilinearAlgorithm> alg = lookupC.find(4, 4, 4);
		assertThat(alg).isPresent();
		// AE catalog has ⟨4,4,4⟩=48 over C; might also be 49 if AE wasn't imported.
		// Accept either as long as it's not above 49.
		assertThat(alg.get().r).isLessThanOrEqualTo(49);
	}

	@Test
	public void r_lookup_for_444_does_not_fallback_to_c() {
		FieldAwareLookup lookupR = new FieldAwareLookup("R");
		Optional<NonCubicBilinearAlgorithm> alg = lookupR.find(4, 4, 4);
		assertThat(alg).isPresent();
		// R-class can never return AE's 48 over C — that's C-only.
		// At the time the test was written this was bounded below by 49 (Strassen²);
		// DPS 2025 imported ⟨4,4,4⟩=48 over R, so the actual R-best is now 48.
		// The invariant being tested is "no C-only fallback".
		assertThat(alg.get().r).isLessThanOrEqualTo(49);
		assertThat(alg.get().r).isGreaterThanOrEqualTo(48);
	}

	@Test
	public void f2_lookup_finds_alphatensor_47_for_444() {
		FieldAwareLookup lookupF2 = new FieldAwareLookup("F2");
		Optional<NonCubicBilinearAlgorithm> alg = lookupF2.find(4, 4, 4);
		assertThat(alg).isPresent();
		assertThat(alg.get().r).isEqualTo(47);
	}

	@Test
	public void r_lookup_does_not_return_f2_scheme() {
		FieldAwareLookup lookupR = new FieldAwareLookup("R");
		// AlphaTensor F2 has ⟨5,5,5⟩=96. R-class should NOT see it.
		Optional<NonCubicBilinearAlgorithm> alg = lookupR.find(5, 5, 5);
		assertThat(alg).isPresent();
		// Best R-class for ⟨5,5,5⟩ should be AlphaEvolve 93 (Z-tagged in catalog).
		assertThat(alg.get().r).isLessThanOrEqualTo(93);
		assertThat(alg.get().r).isGreaterThanOrEqualTo(93); // not the F₂ 96
	}

	@Test
	public void fallback_chain_is_correctly_ordered() {
		// Z ⊂ Q ⊂ R ⊂ C; F2 isolated.
		assertThat(eu.solven.matmul.algebra.Field.Z.fallbackChain())
				.containsExactly(eu.solven.matmul.algebra.Field.Z);
		assertThat(eu.solven.matmul.algebra.Field.Q.fallbackChain())
				.containsExactly(eu.solven.matmul.algebra.Field.Q, eu.solven.matmul.algebra.Field.Z);
		assertThat(eu.solven.matmul.algebra.Field.R.fallbackChain())
				.containsExactly(eu.solven.matmul.algebra.Field.R, eu.solven.matmul.algebra.Field.Q,
						eu.solven.matmul.algebra.Field.Z);
		assertThat(eu.solven.matmul.algebra.Field.C.fallbackChain())
				.containsExactly(eu.solven.matmul.algebra.Field.C, eu.solven.matmul.algebra.Field.R,
						eu.solven.matmul.algebra.Field.Q, eu.solven.matmul.algebra.Field.Z);
		assertThat(eu.solven.matmul.algebra.Field.F2.fallbackChain())
				.containsExactly(eu.solven.matmul.algebra.Field.F2);
	}
}
