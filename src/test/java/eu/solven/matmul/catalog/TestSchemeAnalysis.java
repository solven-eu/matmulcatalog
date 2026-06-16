package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;

/**
 * The expanded-scheme analyses (#159 pipeline): verify / additions / buds run
 * off a single expanded scheme and stamp JSON fields. Validated against
 * committed catalog schemes.
 */
public class TestSchemeAnalysis {

	private static NonCubicBilinearAlgorithm read(String path) throws Exception {
		return SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint(path));
	}

	/** Run the default analysis set and merge all stamped fields. */
	private static Map<String, Object> runAll(NonCubicBilinearAlgorithm alg, Field field) {
		Map<String, Object> out = new HashMap<>();
		for (SchemeAnalysis a : SchemeAnalysis.defaults()) {
			out.putAll(a.analyse(alg, field));
		}
		return out;
	}

	@Test
	public void laderman_3x3x3_verifies_with_98_additions() throws Exception {
		NonCubicBilinearAlgorithm a =
				read("src/main/resources/schemes/known/section3/laderman_1976-3x3x3_m23_a98.json");
		Map<String, Object> r = runAll(a, Field.Q);
		assertThat(r.get("verified")).isEqualTo(true);
		assertThat(r.get("additions")).isEqualTo(98);
		assertThat(r).containsKey("has_buds");
	}

	@Test
	public void strassen_2x2x2_has_no_buds() throws Exception {
		NonCubicBilinearAlgorithm a =
				read("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json");
		SchemeAnalysis.Buds buds = new SchemeAnalysis.Buds();
		Map<String, Object> r = buds.analyse(a, Field.Q);
		// Rank-minimal Strassen mixes all entries → zero buds (paper §serendipitous).
		assertThat(r.get("has_buds")).isEqualTo(false);
		assertThat(r).doesNotContainKey("buds");
	}

	@Test
	public void verify_field_routing_and_spotcheck() throws Exception {
		// Strassen ⟨2,2,2⟩=7 verifies over Q (spot-check) and over F2 (exact mod-2).
		NonCubicBilinearAlgorithm a =
				read("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json");
		assertThat(new SchemeAnalysis.Verify().analyse(a, Field.Q).get("verified")).isEqualTo(true);
		assertThat(new SchemeAnalysis.Verify().analyse(a, Field.F2).get("verified")).isEqualTo(true);
	}
}
