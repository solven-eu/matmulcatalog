package eu.solven.matmul;

import org.junit.jupiter.api.Tag;


import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.catalog.SchemeIO;
import tools.jackson.databind.JsonNode;

/**
 * Live catalog spot-check for {@link SymbolicVerifier}.
 *
 * <p>Walks the small-shape slice of the catalog and asserts every
 * non-stub bilinear scheme verifies symbolically in its declared
 * algebra. This is the audit pass requested by task #97 — failures
 * pinpoint files whose declared field is too tight for their actual
 * coefficient ring (e.g. {@code _Q.json} but contains {@code i}).</p>
 *
 * <p>The slice covers section2..section10. Larger sections are
 * mostly stubs whose verification requires lineage materialisation,
 * which is out of scope for this test — see
 * {@code VerifyAllSchemes --strict} for the full audit.</p>
 */
@Tag("slow")
public class TestSymbolicVerifierCatalogSpotCheck {

	/**
	 * Files we ALREADY know are mistagged and are tracked separately. Names use
	 * the current {@code _m{rank}} convention (#120/#173); the half-integer
	 * recombine schemes have since gained {@code fields:["Q","R","C"]} via the
	 * #174/#175 migration so most now verify — the list is kept as a tolerant
	 * safety net rather than a known-failing set.
	 */
	private static final List<String> KNOWN_MISTAGS = List.of(
			"derived_strassen_recombine-4x7x9_m186.json",
			"derived_strassen_recombine-4x8x9_m206.json",
			"derived_strassen_recombine-4x7x10_m203.json",
			"derived_strassen_recombine-5x8x10_m284.json",
			"derived_strassen_recombine-4x8x10_m224.json"
	);

	@Test
	public void smallCatalogSlice_verifies_in_declared_algebra() throws IOException {
		List<Path> slice = new ArrayList<>();
		// Slice limited to sections 2..10 (small-shape, all non-stub).
		// Wider sweep is delegated to {@code VerifyAllSchemes --strict}
		// which has shape-weighted parallel iteration.
		for (String section : new String[]{
				"section2", "section3", "section4", "section5", "section6",
				"section7", "section8", "section9", "section10"}) {
			Path dir = Path.of("src/main/resources/schemes", section);
			if (!Files.isDirectory(dir)) continue;
			try (Stream<Path> walk = Files.list(dir)) {
				walk.filter(p -> p.toString().endsWith(".json")).forEach(slice::add);
			}
		}
		assertThat(slice).as("expected catalog files in section2..section10").isNotEmpty();

		List<String> unknownFailures = new ArrayList<>();
		List<String> knownFailures = new ArrayList<>();
		for (Path p : slice) {
			File f = p.toFile();
			JsonNode root = SchemeIO.parseJson(f);
			if (SchemeIO.isStub(root)) continue;
			if (SchemeIO.isNonBilinear(root)) continue;
			SymbolicVerifier.Result r = SymbolicVerifier.verify(root, f.getName());
			if (!r.verified()) {
				String line = p.getFileName() + " : " + r;
				if (KNOWN_MISTAGS.contains(f.getName())) {
					knownFailures.add(line);
				} else {
					unknownFailures.add(line);
				}
			}
		}

		// Known mistags (derived-strassen-recombine half-integer schemes
		// that should be tagged Q but lack any field marker) are recorded
		// here; rename / add "field": "Q" to clear them.
		if (!knownFailures.isEmpty()) {
			System.out.println("[catalog-audit] known mistags reproduced ("
					+ knownFailures.size() + "):");
			knownFailures.forEach(System.out::println);
		}
		// Anything else is a regression and must be investigated.
		assertThat(unknownFailures)
				.as("Unknown field-discipline failures (besides the tracked "
						+ KNOWN_MISTAGS.size() + " derived-strassen-recombine files)")
				.isEmpty();
	}
}
