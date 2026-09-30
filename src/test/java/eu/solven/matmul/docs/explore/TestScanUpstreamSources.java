package eu.solven.matmul.docs.explore;

import static org.assertj.core.api.Assertions.assertThat;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import eu.solven.matmul.CanonicalShape;
import eu.solven.matmul.Shape;

/**
 * Guards the 2026-06 → 2026-09 silent blindness of the weekly upstream scan:
 * {@code loadLocalRanks} parsed "ours" out of the FILENAME with a legacy
 * {@code {source}-{shape}_m{rank}} regex; after the catalog rename to
 * {@code {n}x{m}x{p}-r{rank}-{note}-{hash7}} it matched 7 of ~11k files, so every
 * upstream shape read as MISSING (never a GAP) and the scan reported
 * "Strict gaps: 0" for three months — through the 2026-09 KGP LITA refresh that
 * put ~530 shapes below the catalog. No crash, just an empty report.
 */
public class TestScanUpstreamSources {

	private static void write(Path dir, String name, String body) throws Exception {
		Files.createDirectories(dir);
		Files.writeString(dir.resolve(name), body, StandardCharsets.UTF_8);
	}

	@Test
	public void local_ranks_are_read_from_content_not_filenames(@TempDir Path root) throws Exception {
		// New-convention name, tiny dense scheme body (content is what must be read).
		write(root.resolve("known/section2"), "2x2x2-r7-strassen-db11bcc.json",
				"{\"n\": [2, 2, 2], \"m\": 7, \"u\": [], \"fields\": [\"F2\", \"F3\", \"Z\", \"Q\", \"R\", \"C\"]}");
		// A name that carries NO rank/shape token at all: still read.
		write(root.resolve("derived/section3"), "whatever.json",
				"{\"n\": [3, 2, 2], \"m\": 11, \"fields\": [\"Z\", \"Q\", \"R\", \"C\"]}");
		// Same shape, worse rank → min wins.
		write(root.resolve("derived/section3"), "2x2x3-r12-derived-0000000.json",
				"{\"n\": [2, 2, 3], \"m\": 12, \"fields\": [\"Q\", \"R\", \"C\"]}");
		Map<CanonicalShape, Integer> local = ScanUpstreamSources.loadLocalRanks(root);
		assertThat(local).containsEntry(Shape.of(2, 2, 2).canonical(), 7);
		assertThat(local).containsEntry(Shape.of(2, 2, 3).canonical(), 11);
	}

	@Test
	public void local_ranks_exclude_incomparable_schemes(@TempDir Path root) throws Exception {
		// Commutative-only (Waksman) must never "close" a non-commutative digest gap.
		write(root.resolve("constructed/section2"), "2x2x2-r7-waksman_1970-0d1209a.json",
				"{\"n\": [2, 2, 2], \"m\": 7, \"commutative\": true, \"fields\": [\"Z\", \"Q\", \"R\", \"C\"]}");
		// F₂-native (AlphaTensor ⟨4,4,4⟩=47) is not a Q rank.
		write(root.resolve("known/section4"), "4x4x4-r47-alphatensor_F2-0000000.json",
				"{\"n\": [4, 4, 4], \"m\": 47, \"fields\": [\"F2\"]}");
		write(root.resolve("known/section4"), "4x4x4-r49-alphatensor_Z-0000000.json",
				"{\"n\": [4, 4, 4], \"m\": 49, \"fields\": [\"F2\", \"F3\", \"Z\", \"Q\", \"R\", \"C\"]}");
		// A non-scheme JSON (e.g. a manifest) is skipped, not fatal.
		write(root, "catalog.json", "{\"schemes\": []}");
		Map<CanonicalShape, Integer> local = ScanUpstreamSources.loadLocalRanks(root);
		assertThat(local).doesNotContainKey(Shape.of(2, 2, 2).canonical());
		assertThat(local).containsEntry(Shape.of(4, 4, 4).canonical(), 49);
	}

	/** The real catalog: thousands of shapes, ⟨2,2,2⟩=7 — a regex-blind loader returns ~7 shapes. */
	@Test
	public void real_catalog_yields_thousands_of_shapes() throws Exception {
		Map<CanonicalShape, Integer> local = ScanUpstreamSources.loadLocalRanks(Path.of("src/main/resources/schemes"));
		assertThat(local.size()).as("content-driven load must see the whole catalog").isGreaterThan(4000);
		assertThat(local).containsEntry(Shape.of(2, 2, 2).canonical(), 7);
	}
}
