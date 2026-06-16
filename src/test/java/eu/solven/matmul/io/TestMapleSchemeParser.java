package eu.solven.matmul.io;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Parses the FMM-Lille {@code 17x17x17_raw.mpl} scheme and validates it
 * against the matmul tensor for ⟨17,17,17⟩ over Q.
 */
public class TestMapleSchemeParser {

	private static final File RAW = new File("references/fmm-lille/17x17x17/17x17x17_raw.mpl");

	@Test
	public void parses_17x17x17_raw_with_correct_dimensions() throws Exception {
		assertThat(RAW).as("FMM-Lille raw.mpl must exist for this test").exists();
		NonCubicBilinearAlgorithm alg = MapleSchemeParser.parseRawFmmLille(RAW, 17, 17, 17);
		assertThat(alg.n).isEqualTo(17);
		assertThat(alg.m).isEqualTo(17);
		assertThat(alg.p).isEqualTo(17);
		// The raw file has 2945 products; the published rank 2934 reflects
		// a post-processing kin-row reduction not applied here.
		assertThat(alg.r).isEqualTo(2945);
	}

	/**
	 * Fast randomised correctness check ({@code O(samples · r · (nm + mp + np))}).
	 * For ⟨17,17,17⟩ at rank 2945 with Q-rational coefficients this is decisive:
	 * any structural error in the parse would explode the Frobenius residual.
	 */
	@Test
	public void fmm_lille_scheme_computes_matmul() throws Exception {
		assertThat(RAW).exists();
		NonCubicBilinearAlgorithm alg = MapleSchemeParser.parseRawFmmLille(RAW, 17, 17, 17);
		// epsilon scaled relative to dimensions and rank — Q-rational coefs with
		// 1/8 denominators only contribute trivial rounding, so a tight bound works.
		double eps = 1e-6 * Math.sqrt(alg.n * alg.p) * Math.sqrt(alg.r);
		assertThat(Verifier.passesRandomMatmulSpotCheck(alg, 5, eps))
				.as("FMM-Lille ⟨17,17,17⟩ raw.mpl scheme must compute matrix multiplication")
				.isTrue();
	}

	@Test
	public void addition_count_is_reported() throws Exception {
		assertThat(RAW).exists();
		NonCubicBilinearAlgorithm alg = MapleSchemeParser.parseRawFmmLille(RAW, 17, 17, 17);
		int adds = Verifier.additionCount(alg);
		// Sanity: a non-trivial ⟨17,17,17⟩ scheme at rank 2945 must have a large
		// but bounded addition count. Just guard the plausibility envelope here —
		// the exact value goes into the JSON filename downstream.
		assertThat(adds).as("addition count must be positive").isPositive();
		System.out.println("[MapleSchemeParser] FMM-Lille ⟨17,17,17⟩ raw.mpl: r="
				+ alg.r + ", additions=" + adds);
	}

	/**
	 * Once the FMM-Lille JSON has been written to {@code section17/}, reading
	 * it back through {@link SchemeIO} must produce a scheme that still
	 * computes ⟨17,17,17⟩ — guards the import + sparse-format round-trip.
	 */
	@Test
	public void section17_catalog_json_round_trips() throws Exception {
		File catalog = eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section17/fmm_lille_2025-17x17x17_m2945_a68812.json");
		if (!catalog.isFile()) {
			// Catalog write is driven by ImportFmmLille17; skip cleanly if not yet emitted.
			System.out.println("[MapleSchemeParser] catalog JSON not present, skipping round-trip");
			return;
		}
		NonCubicBilinearAlgorithm alg = SchemeIO.read(catalog);
		assertThat(alg.n).isEqualTo(17);
		assertThat(alg.m).isEqualTo(17);
		assertThat(alg.p).isEqualTo(17);
		assertThat(alg.r).isEqualTo(2945);
		double eps = 1e-6 * Math.sqrt(alg.n * alg.p) * Math.sqrt(alg.r);
		assertThat(Verifier.passesRandomMatmulSpotCheck(alg, 5, eps))
				.as("FMM-Lille catalog JSON must round-trip to a correct ⟨17,17,17⟩ scheme")
				.isTrue();
	}
}
