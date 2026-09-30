package eu.solven.matmul.docs.migrate;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Stream;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.docs.migrate.ImportContributedSchemes.Gate;

/**
 * Guards for the contributed-scheme import gate (issue #7, the six HK task-#9
 * shapes). The gate is the ONLY thing standing between a third-party zip and
 * {@code schemes/known/}; it must accept a genuine scheme, reject a
 * single-coefficient sabotage (the failure mode is silent — a wrong scheme has
 * the right shape/rank and still parses), and stamp fields from the
 * coefficients rather than from the contributor's own claims.
 */
public class TestImportContributedSchemes {

	private static final Path SCHEMES = Path.of("src/main/resources/schemes");
	/** Small committed HK scheme: ⟨2,3,4⟩=20, integer, fast to verify exactly. */
	private static final File HK_2X3X4 =
			SCHEMES.resolve("constructed/section4/2x3x4-r20-hk71-89186e6.json").toFile();

	@Test
	public void gate_accepts_a_genuine_integer_scheme_and_stamps_fields_from_coefficients() throws Exception {
		Gate g = ImportContributedSchemes.gate(HK_2X3X4, null);
		assertThat(g.passed()).isTrue();
		assertThat(g.exactSymbolic()).isTrue();
		assertThat(g.spotCheck()).isTrue();
		assertThat(g.alg().r).isEqualTo(20);
		assertThat(g.integer()).isTrue();
		assertThat(g.fields()).containsExactly("F2", "F3", "Z", "Q", "R", "C");
		assertThat(g.catalogBest()).as("no lookup → no catalog comparison").isNull();
	}

	@Test
	public void gate_rejects_a_single_flipped_coefficient(@TempDir Path tmp) throws Exception {
		// Flip the sign of the first non-zero U coefficient: same shape, same rank,
		// still parses — only the bilinear identity breaks. Both the exact
		// symbolic proof AND the independent numeric spot-check must say no.
		String json = Files.readString(HK_2X3X4.toPath(), StandardCharsets.UTF_8);
		// Dense layout: the first U row is "[0, 0, 0, 1, 0, 0]".
		String sabotaged = json.replaceFirst("\\[0, 0, 0, 1, 0, 0\\]", "[0, 0, 0, -1, 0, 0]");
		assertThat(sabotaged).as("the sabotage must actually change the file").isNotEqualTo(json);
		File f = tmp.resolve("flipped.json").toFile();
		Files.writeString(f.toPath(), sabotaged, StandardCharsets.UTF_8);

		Gate g = ImportContributedSchemes.gate(f, null);
		assertThat(g.passed()).isFalse();
		assertThat(g.exactSymbolic()).isFalse();
		assertThat(g.spotCheck()).isFalse();
	}

	@Test
	public void gate_never_trusts_the_contributors_fields_claim(@TempDir Path tmp) throws Exception {
		// A contributor claiming ["Z"] on a scheme that is not integer must get the
		// coefficient-driven answer (no Z), not the claim.
		String json = Files.readString(HK_2X3X4.toPath(), StandardCharsets.UTF_8);
		// Rationals are stored as quoted "p/q" tokens (the canonical writer's form).
		String halved = json.replaceFirst("\\[0, 0, 0, 1, 0, 0\\]", "[0, 0, 0, \"1/2\", 0, 0]");
		assertThat(halved).as("the edit must actually change the file").isNotEqualTo(json);
		File f = tmp.resolve("halved.json").toFile();
		Files.writeString(f.toPath(), halved, StandardCharsets.UTF_8);

		Gate g = ImportContributedSchemes.gate(f, null);
		assertThat(g.integer()).isFalse();
		assertThat(g.fields()).doesNotContain("Z");
		// (it also no longer computes matmul — the gate is coefficient-blind on
		// purpose: it reports fields for what it read, and passed() separately)
		assertThat(g.passed()).isFalse();
	}

	@Test
	public void meta_file_is_read_as_utf8_and_flags_override_it(@TempDir Path tmp) throws Exception {
		// Windows CLI args are Cp1252-decoded; the ⌈…⌉ / ⟨…⟩ glyphs the notes use
		// only survive through the UTF-8 --meta file.
		Path meta = tmp.resolve("meta.json");
		Files.writeString(meta, "{\"note\":\"someone_2026\",\"source\":\"Someone 2026\",\"year\":2026,"
				+ "\"discovery\":\"false\",\"attribution_for_rank\":\"HK71 ⌈(3pn+max)/2⌉ for ⟨2,p,n⟩\"}",
				StandardCharsets.UTF_8);
		ImportContributedSchemes.Options o = ImportContributedSchemes.parseArgs(new String[] {
				"some/dir", "--meta=" + meta, "--year=2025" });
		assertThat(o.note()).isEqualTo("someone_2026");
		assertThat(o.source()).isEqualTo("Someone 2026");
		assertThat(o.year()).as("explicit flag overrides the meta file").isEqualTo(2025);
		assertThat(o.discovery()).isEqualTo("false");
		assertThat(o.attribution()).isEqualTo("HK71 ⌈(3pn+max)/2⌉ for ⟨2,p,n⟩");
		assertThat(o.execute()).isFalse();
	}

	@Test
	public void discovery_false_requires_an_attribution_when_executing() {
		assertThatThrownBy(() -> ImportContributedSchemes.parseArgs(new String[] {
				"x", "--note=a_2026", "--source=A 2026", "--discovery=false", "--execute" }))
				.isInstanceOf(IllegalArgumentException.class)
				.hasMessageContaining("--attribution");
	}

	/**
	 * The issue-#7 import itself: six {@code adriano_2026} files on disk, each at
	 * the exact HK formula {@code ⌈(3pn+max(p,n))/2⌉}, integer (Z ∈ fields[]),
	 * verified, and attributed (discovery:false + attribution_for_rank). Losing
	 * any of them would silently regress the six task-#9 shapes back to the
	 * emitter's +1..+3 — see TestSweepSpotsSota for the rank-level guard.
	 */
	@Test
	public void issue7_hk_task9_schemes_are_on_disk_at_the_formula() throws Exception {
		int[][] shapes = { { 2, 12, 18 }, { 2, 14, 21 }, { 2, 16, 24 }, { 2, 18, 27 }, { 2, 20, 30 }, { 2, 24, 30 } };
		for (int[] s : shapes) {
			int p = s[1], n = s[2];
			int formula = (int) Math.ceil((3.0 * p * n + Math.max(p, n)) / 2.0);
			Path dir = SCHEMES.resolve("known/section" + n);
			List<Path> hits;
			try (Stream<Path> ls = Files.list(dir)) {
				hits = ls.filter(f -> f.getFileName().toString()
						.startsWith("2x" + p + "x" + n + "-r" + formula + "-adriano_2026-")).toList();
			}
			assertThat(hits).as("⟨2,%d,%d⟩ adriano_2026 file at r=%d", p, n, formula).hasSize(1);
			var root = SchemeIO.parseJson(hits.get(0).toFile());
			assertThat(root.get("m").asInt()).isEqualTo(formula);
			assertThat(SchemeIO.fieldTags(root)).contains("Z", "Q");
			assertThat(SchemeIO.isVerified(root)).isTrue();
			assertThat(root.get("discovery").asBoolean()).isFalse();
			assertThat(root.get("attribution_for_rank").asString()).contains("Hopcroft").contains("Adriano");
			assertThat(root.get("source").asString()).isEqualTo("Adriano 2026");
		}
	}
}
