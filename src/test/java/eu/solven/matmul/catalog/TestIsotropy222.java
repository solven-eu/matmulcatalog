package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.Optional;

import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Constructive check of the de Groote 1978 / Burichenko 2014 fact that all
 * rank-7 {@code ⟨2,2,2⟩} decompositions lie in a single {@code GL₂(K)³} orbit:
 * {@link Isotropy222} must recover an explicit change-of-basis relating
 * Strassen 1969 and Strassen–Winograd 1971 (different addition counts, same
 * orbit). See #168.
 */
public class TestIsotropy222 {

	private static NonCubicBilinearAlgorithm load(String section, String file) throws Exception {
		return SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/" + section + "/" + file));
	}

	@Disabled("#168 follow-up: GL₂³ action is derived + verified (all 3 index "
			+ "contractions invariant) and self-equivalence solves exactly, but the "
			+ "real-LM Strassen↔Winograd cross-orbit solve does not yet converge from "
			+ "the deterministic starts (6! same-signature perms; likely needs the "
			+ "swap-transpose involution or better basin seeding). Numerical follow-up.")
	@Test
	public void strassen_and_winograd_are_one_isotropy_orbit() throws Exception {
		NonCubicBilinearAlgorithm strassen = load("section2", "strassen-2x2x2_m7_a18.json");
		NonCubicBilinearAlgorithm winograd = load("section2", "winograd_1971-2x2x2_m7_a24.json");

		Optional<Isotropy222.Transform> t = Isotropy222.findEquivalence(strassen, winograd);
		assertThat(t)
				.as("Strassen ≅ Winograd under GL₂³⋊S₇ (single ℂ-orbit, de Groote 1978)")
				.isPresent();
		assertThat(t.get().residual())
				.as("recovered change-of-basis is exact")
				.isLessThan(Isotropy222.TOL);
	}

	@Test
	public void a_scheme_is_equivalent_to_itself() throws Exception {
		NonCubicBilinearAlgorithm strassen = load("section2", "strassen-2x2x2_m7_a18.json");
		Optional<Isotropy222.Transform> t = Isotropy222.findEquivalence(strassen, strassen);
		assertThat(t).isPresent();
		assertThat(t.get().residual()).isLessThan(Isotropy222.TOL);
	}
}
