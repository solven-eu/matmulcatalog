package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Sanity for the human-friendly renderer (#188): Strassen ⟨2,2,2⟩=7 must yield
 * 7 multiplication lines and 4 element lines in the Perminov-style syntax.
 */
public class TestHumanScheme {

	@Test
	public void strassen_222_renders_7_products_and_4_elements() throws Exception {
		NonCubicBilinearAlgorithm strassen =
				SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		HumanScheme.Readable hr = HumanScheme.of(strassen);

		assertThat(hr.multiplications()).hasSize(7);
		assertThat(hr.elements()).hasSize(4);
		// Shape of a product line: "mK = (… a …)*(… b …)".
		assertThat(hr.multiplications().get(0))
				.startsWith("m1 = (")
				.contains(")*(")
				.contains("a")
				.contains("b");
		// Element lines reference products and the right output cells c11..c22.
		assertThat(hr.elements().get(0)).startsWith("c11 = ").contains("m");
		assertThat(hr.elements()).anyMatch(e -> e.startsWith("c22 = "));
	}
}
