package eu.solven.matmul.docs.migrate;

import static org.assertj.core.api.Assertions.assertThat;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Optional;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.LineageReplayer;
import tools.jackson.databind.JsonNode;

/**
 * {@link RestampStubRanks} repairs a stub whose recorded rank drifted from its replay:
 * the rewritten stub records the replayed rank, keeps the lineage, fields and
 * provenance keys, is renamed to the new rank/hash, and the stale file is gone.
 */
public class TestRestampStubRanks {
	private static FieldAwareLookup lookup;
	private static LineageReplayer replayer;

	@BeforeAll
	static void setup() {
		lookup = new FieldAwareLookup("Q");
		replayer = LineageReplayer.withDefaultPool(lookup);
	}

	private static Lineage.Node projectionOfStrassen() throws Exception {
		// A hash-stamped explicit rank-7 ⟨2,2,2⟩ (the lookup's first file may be a legacy
		// record without a `hash`).
		Path strassen = lookup.findFiles(2, 2, 2).stream().filter(p -> {
			try {
				JsonNode r = SchemeIO.parseJson(p.toFile());
				return r.has("hash") && r.path("m").asInt() == 7 && !SchemeIO.isStub(r) && !SchemeIO.isNonBilinear(r)
						&& !SchemeIO.isComplex(r);
			} catch (Exception e) {
				return false;
			}
		}).findFirst().orElseThrow();
		String hash = SchemeIO.parseJson(strassen.toFile()).get("hash").asString();
		return new Lineage.Project(new Lineage.Atom("2x2x2@" + hash), new int[] { 0, 1 }, new int[] { 0, 1 },
				new int[] { 0 });
	}

	@Test
	public void an_over_claiming_stub_is_rewritten_at_its_replayed_rank(@TempDir Path dir) throws Exception {
		Lineage.Node lineage = projectionOfStrassen();
		NonCubicBilinearAlgorithm honest = replayer.replay(lineage);
		Path stub = dir.resolve("2x2x1-r" + (honest.r - 1) + "-derived-0000000.json");
		SchemeIO.writeStub(honest, stub.toFile(), lineage, List.of("Z", "Q", "R", "C"));
		String overClaimed = Files.readString(stub).replace("\"m\": " + honest.r + ",", "\"m\": " + (honest.r - 1) + ",");
		Files.writeString(stub, overClaimed);
		SchemeIO.addFields(stub.toFile(), java.util.Map.of("source", "test-provenance", "year", 2026), true);

		Optional<RestampStubRanks.Drift> drift = RestampStubRanks.driftOf(stub, replayer);
		assertThat(drift).isPresent();
		assertThat(drift.get().overClaim()).isTrue();
		assertThat(drift.get().replayed()).isEqualTo(honest.r);

		Path out = RestampStubRanks.restamp(stub, replayer);

		assertThat(Files.exists(stub)).as("the stale file is deleted").isFalse();
		assertThat(out.getFileName().toString()).startsWith("2x2x1-r" + honest.r + "-derived-");
		JsonNode fixed = SchemeIO.parseJson(out.toFile());
		assertThat(fixed.get("m").asInt()).isEqualTo(honest.r);
		assertThat(SchemeIO.isStub(fixed)).isTrue();
		assertThat(SchemeIO.fieldTags(fixed)).containsExactly("Z", "Q", "R", "C");
		assertThat(fixed.get("source").asString()).as("provenance keys carried over").isEqualTo("test-provenance");
		assertThat(fixed.get("year").asInt()).isEqualTo(2026);
		assertThat(RestampStubRanks.driftOf(out, replayer)).as("honest after the rewrite").isEmpty();
	}

	/**
	 * A loose stub whose honest version already exists (a sibling replaying to the same
	 * content at the same rank) is a pure duplicate: the stale copy goes, the existing
	 * file — lineage, provenance, computed metrics — is not touched, and the result points
	 * at the existing file so phase 2 re-points the duplicate's pins to its hash. (The
	 * 2026-10-18 bulk run overwrote 162 such files before this rule existed.)
	 */
	@Test
	public void a_loose_duplicate_of_an_existing_honest_stub_is_removed_not_overwritten(@TempDir Path dir)
			throws Exception {
		Lineage.Node lineage = projectionOfStrassen();
		NonCubicBilinearAlgorithm honest = replayer.replay(lineage);
		Path existing = dir.resolve(SchemeIO.canonicalName(honest, "derived"));
		SchemeIO.writeStub(honest, existing.toFile(), lineage, List.of("Z", "Q", "R", "C"));
		SchemeIO.addFields(existing.toFile(), java.util.Map.of("source", "the-honest-one", "additions", 12), true);
		String before = Files.readString(existing);

		Path loose = dir.resolve("2x2x1-r" + (honest.r + 3) + "-derived-1234567.json");
		SchemeIO.writeStub(honest, loose.toFile(), lineage, List.of("Z", "Q", "R", "C"));
		Files.writeString(loose,
				Files.readString(loose).replace("\"m\": " + honest.r + ",", "\"m\": " + (honest.r + 3) + ","));
		assertThat(RestampStubRanks.driftOf(loose, replayer)).isPresent();

		RestampStubRanks.Restamped r = RestampStubRanks.restamp(loose, replayer, null, true);

		assertThat(Files.exists(loose)).as("the stale duplicate is gone").isFalse();
		assertThat(r.newFile()).isEqualTo(existing);
		assertThat(Files.readString(existing)).as("the honest file is untouched").isEqualTo(before);
		assertThat(r.newHash()).isEqualTo(SchemeIO.readHash(SchemeIO.parseJson(existing.toFile())));
	}

	/**
	 * A restamped stub gets a new stamped hash, so every lineage that pinned the OLD hash
	 * ({@code shape@hash7} or the full hash) would dangle — the bulk run re-points those
	 * pins (exact by construction: the stub's lineage is unchanged, so the pinning file
	 * already replayed through it to the same matrices). Other pins are untouched and
	 * the rewritten file stays canonically formatted.
	 */
	@Test
	public void pins_to_a_restamped_stubs_old_hash_are_rewritten(@TempDir Path dir) throws Exception {
		String oldHash = "a".repeat(64);
		String newHash = "b".repeat(64);
		Path child = dir.resolve("section3").resolve("2x2x3-r11-derived-0000000.json");
		Files.createDirectories(child.getParent());
		String json = "{\n  \"n\": [2, 2, 3],\n  \"m\": 11,\n  \"scheme_type\": \"stub\",\n"
				+ "  \"lineage\": {\"op\": \"ConcatCols\", \"children\": ["
				+ "{\"op\": \"Atom\", \"ref\": \"2x2x1@" + oldHash.substring(0, 7) + "\"}, "
				+ "{\"op\": \"Atom\", \"ref\": \"2x2x2@" + oldHash + "\"}, "
				+ "{\"op\": \"Atom\", \"ref\": \"2x2x2@db11bcc\"}]}\n}";
		Files.writeString(child, MatrixJsonFormatter.format(json));
		Path untouched = dir.resolve("section3").resolve("3x3x3-r23-derived-1111111.json");
		Files.writeString(untouched, MatrixJsonFormatter.format("{\n  \"n\": [3, 3, 3],\n  \"m\": 23,\n"
				+ "  \"lineage\": {\"op\": \"Atom\", \"ref\": \"3x3x3@db11bcc\"}\n}"));
		String untouchedBefore = Files.readString(untouched);

		// Transitive: the stub was a duplicate of a twin (old → twin's hash) and the twin
		// itself was restamped in the same run (twin's hash → newHash); pins must land on the
		// final hash, not on an intermediate one that no file carries any more.
		String twinHash = "c".repeat(64);
		int rewritten = RestampStubRanks.rewritePins(dir, List.of(
				new RestampStubRanks.Restamped(dir.resolve("old.json"), dir.resolve("twin.json"), oldHash, twinHash),
				new RestampStubRanks.Restamped(dir.resolve("twin.json"), dir.resolve("new.json"), twinHash, newHash)));

		assertThat(rewritten).isEqualTo(1);
		String after = Files.readString(child);
		assertThat(after).contains("2x2x1@" + newHash.substring(0, 7)).contains("2x2x2@" + newHash)
				.contains("2x2x2@db11bcc").doesNotContain(oldHash.substring(0, 7)).doesNotContain(twinHash.substring(0, 7));
		assertThat(after).as("still canonically formatted")
				.isEqualTo(MatrixJsonFormatter.format(SchemeIO.parseJson(child.toFile())));
		assertThat(Files.readString(untouched)).isEqualTo(untouchedBefore);
	}

	@Test
	public void note_token_is_taken_from_the_old_name() {
		assertThat(RestampStubRanks.noteOf("13x13x17-r1876-derived-81d3a84.json")).isEqualTo("derived");
		assertThat(RestampStubRanks.noteOf("16x20x28-r4960-derived_serendip-818b551.json")).isEqualTo("derived_serendip");
		assertThat(RestampStubRanks.noteOf("legacy-name.json")).isEqualTo("derived");
	}
}
