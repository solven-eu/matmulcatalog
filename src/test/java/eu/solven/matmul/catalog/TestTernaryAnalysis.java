package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.Map;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;

/**
 * Guards the {@link SchemeAnalysis.Ternary} enrichment analysis — the
 * dense-representation computation of the {@code zt} ("ternary integer",
 * Perminov's sub-class of {@code Z} whose every U/V/W coefficient is in
 * {@code {-1,0,1}}) flag.
 *
 * <p>Regression context: the entire {@code derived/} tree (7784 files) was
 * 0/7784 stamped for {@code zt} because the flag was only ever written by a
 * separate {@code MaterialiseZT} pass that ran over {@code known/} alone. So a
 * genuinely-ternary derived scheme (e.g. {@code 2x2x3-r11-derived}, a projection
 * of a ⟨2,2,4⟩ atom with all coefficients in {@code {-1,0,1}}) showed
 * {@code zt = null} in the manifest and was invisible under the SPA's ZT
 * selector. Folding the computation into {@link SchemeAnalysis#defaults()} makes
 * every enrichment pass stamp {@code zt} from the dense {@code alg} it already
 * expanded, so no scheme can silently miss it again.</p>
 */
public class TestTernaryAnalysis {

	/** The exact derived scheme that surfaced the gap: ⟨2,2,3⟩=11, all coeffs ∈ {-1,0,1}. */
	private static final File DERIVED_223 =
			new File("src/main/resources/schemes/derived/section3/2x2x3-r11-derived-29e94cf.json");

	/** A rational (non-integer) scheme: ZT is meaningless, so nothing is stamped. */
	private static final File ROSOWSKI_223 =
			new File("src/main/resources/schemes/constructed/section3/2x2x3-r10-rosowski_2019_thm2-576d992.json");

	@Test
	public void ternary_is_in_the_default_enrichment_set() {
		assertThat(SchemeAnalysis.defaults())
				.as("Ternary must run in the single-expansion enrichment pass")
				.anyMatch(a -> a instanceof SchemeAnalysis.Ternary);
	}

	@Test
	public void derived_223_is_flagged_zt_true() throws Exception {
		NonCubicBilinearAlgorithm alg = SchemeIO.read(DERIVED_223);
		assertThat(SchemeIO.isTernary(alg))
				.as("the ⟨2,2,3⟩=11 derived scheme has all coeffs in {-1,0,1}")
				.isTrue();

		Map<String, Object> stamped = new SchemeAnalysis.Ternary().analyse(alg, Field.Z);
		assertThat(stamped).containsEntry("zt", Boolean.TRUE);
	}

	@Test
	public void non_Z_field_stamps_nothing() throws Exception {
		// ZT is a sub-class of Z only — for any non-integer field the analysis must
		// emit NOTHING (not even zt:false), matching MaterialiseZT's omit/clear path.
		NonCubicBilinearAlgorithm alg = SchemeIO.read(DERIVED_223);
		assertThat(new SchemeAnalysis.Ternary().analyse(alg, Field.Q))
				.as("zt is meaningless off Z")
				.isEmpty();
		assertThat(new SchemeAnalysis.Ternary().analyse(alg, Field.C)).isEmpty();
		assertThat(new SchemeAnalysis.Ternary().analyse(alg, Field.F2)).isEmpty();
	}

	@Test
	public void non_ternary_Z_scheme_is_flagged_zt_false() throws Exception {
		// A scheme with an integer coefficient outside {-1,0,1} (e.g. a 2) over Z
		// must stamp zt:false — the flag is present but negative, not omitted.
		NonCubicBilinearAlgorithm base = SchemeIO.read(DERIVED_223);
		double[][] u = base.denseU();
		double[][] v = base.denseV();
		double[][] w = base.denseW();
		double[][] u2 = new double[u.length][];
		for (int i = 0; i < u.length; i++) u2[i] = u[i].clone();
		u2[0][0] = 2.0; // push one coefficient outside {-1,0,1}
		NonCubicBilinearAlgorithm alg = new NonCubicBilinearAlgorithm(base.n, base.m, base.p, u2, v, w);

		assertThat(SchemeIO.isTernary(alg)).isFalse();
		assertThat(new SchemeAnalysis.Ternary().analyse(alg, Field.Z))
				.containsEntry("zt", Boolean.FALSE);
	}
}
