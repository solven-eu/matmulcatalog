package eu.solven.matmul.docs.migrate;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Guards the importers' idempotence index. The silent regression it pins
 * (2026-09-30): the Perminov importer keyed "already held" on {@code (shape, rank)}
 * parsed from {@code -perminov_} filenames, so an upstream <b>ternary</b> scheme at a
 * rank we held only as a non-ternary integer scheme was skipped without being
 * downloaded — ⟨2,13,15⟩=304 ZT and six more never reached the catalog while the
 * sync job reported "0 imported".
 */
public class TestKnownSchemeKeys {

	private static final List<String> ALL = List.of("F2", "F3", "Z", "Q", "R", "C");

	/** A dense ⟨1,1,1⟩ "scheme" of the given rank — the scan reads, it does not verify. */
	private static void dense(Path dir, String name, int rank, String u, String extra) throws Exception {
		String ones = "[" + "[1], ".repeat(rank - 1) + "[1]]";
		Files.createDirectories(dir);
		Files.writeString(dir.resolve(name), "{\"n\": [1, 1, 1], \"m\": " + rank + ", \"u\": " + u + ", \"v\": " + ones
				+ ", \"w\": " + ones + ", " + extra + "}", StandardCharsets.UTF_8);
	}

	@Test
	public void a_non_ternary_integer_scheme_does_not_cover_an_upstream_ZT_file(@TempDir Path known) throws Exception {
		// Coefficient 2 → integer but not ternary: the ⟨2,13,15⟩=304 situation.
		dense(known.resolve("section1"), "whatever-name.json", 3, "[[1], [2], [2]]", "\"fields\": [\"Z\", \"Q\", \"R\", \"C\"]");

		Set<String> keys = KnownSchemeKeys.scan(known);

		assertThat(KnownSchemeKeys.covers(keys, 1, 1, 1, 3, "Q")).isTrue();
		assertThat(KnownSchemeKeys.covers(keys, 1, 1, 1, 3, "Z")).isTrue();
		assertThat(KnownSchemeKeys.covers(keys, 1, 1, 1, 3, "ZT"))
				.as("a ternary upstream scheme at this (shape, rank) is NEW — it must be imported")
				.isFalse();
	}

	@Test
	public void class_is_read_from_content_and_the_ladder_is_ZT_Z_Q(@TempDir Path known) throws Exception {
		Path dir = known.resolve("section1");
		// Ternary integer → covers all three classes. Written by the canonical writer
		// at a NON-sorted orientation: the key is canonical (sorted dims).
		File ternary = known.resolve("section3").resolve("label-says-nothing.json").toFile();
		Files.createDirectories(ternary.getParentFile().toPath());
		SchemeIO.write(NonCubicBilinearAlgorithm.naive(3, 1, 2), ternary);
		SchemeIO.addFields(ternary, Map.of("fields", ALL), true);
		// Rational → Q only.
		dense(dir, "rational.json", 2, "[[\"1/2\"], [\"1/2\"]]", "\"fields\": [\"Q\", \"R\", \"C\"]");
		// Commutative-only and F2-native schemes say nothing about the NC char-0 classes.
		dense(dir, "commutative.json", 5, "[[1], [1], [1], [1], [1]]",
				"\"fields\": [\"Z\", \"Q\", \"R\", \"C\"], \"commutative\": true");
		dense(dir, "f2.json", 4, "[[1], [1], [1], [1]]", "\"fields\": [\"F2\"]");
		// Not a scheme at all.
		Files.writeString(dir.resolve("index.json"), "{\"total\": 3}", StandardCharsets.UTF_8);

		Set<String> keys = KnownSchemeKeys.scan(known);

		assertThat(keys).containsExactlyInAnyOrder(
				"1x2x3-r6-Q", "1x2x3-r6-Z", "1x2x3-r6-ZT",
				"1x1x1-r2-Q");
		assertThat(KnownSchemeKeys.covers(keys, 2, 3, 1, 6, "ZT")).as("any orientation").isTrue();
	}

	@Test
	public void upstream_tags_map_to_classes() {
		assertThat(KnownSchemeKeys.classOfTag("ZT")).isEqualTo("ZT");
		assertThat(KnownSchemeKeys.classOfTag("ZT_reduced")).isEqualTo("ZT");
		assertThat(KnownSchemeKeys.classOfTag("Z")).isEqualTo("Z");
		assertThat(KnownSchemeKeys.classOfTag("Q")).isEqualTo("Q");
	}

	/** Fast real-catalog probe (two small sections): content-driven keys exist for the classics. */
	@Test
	public void real_catalog_sections_yield_the_classic_keys() throws Exception {
		Path known = Path.of("src/main/resources/schemes/known");
		Set<String> keys = new java.util.HashSet<>(KnownSchemeKeys.scan(known.resolve("section2")));
		keys.addAll(KnownSchemeKeys.scan(known.resolve("section3")));

		assertThat(keys).as("Strassen ⟨2,2,2⟩=7 is ternary; Laderman-rank ⟨3,3,3⟩=23 is held over Z")
				.contains("2x2x2-r7-ZT", "3x3x3-r23-Z");
	}
}
