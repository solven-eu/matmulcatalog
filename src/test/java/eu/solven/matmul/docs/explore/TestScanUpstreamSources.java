package eu.solven.matmul.docs.explore;

import static org.assertj.core.api.Assertions.assertThat;

import java.net.InetSocketAddress;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Instant;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import com.sun.net.httpserver.HttpServer;

import eu.solven.matmul.CanonicalShape;
import eu.solven.matmul.Shape;
import eu.solven.matmul.docs.explore.ScanUpstreamSources.FmmMeta;
import eu.solven.matmul.docs.explore.ScanUpstreamSources.Gap;
import eu.solven.matmul.docs.explore.ScanUpstreamSources.GapRow;
import eu.solven.matmul.docs.explore.ScanUpstreamSources.SourceScan;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

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

	/**
	 * 2026-10-19: {@code fetch} wrote the download over the existing file WITHOUT
	 * truncating it ({@code BodyHandlers.ofFile(dst)} opens with CREATE+WRITE only), so a
	 * destination longer than the download kept its old tail — a digest that no longer
	 * parsed. On CI the file had never shrunk; locally the CRLF working copy is always
	 * longer than the LF download, and the refreshed digest ended in half of an old
	 * {@code "format"} line.
	 */
	@Test
	public void fetch_truncates_a_longer_existing_file(@TempDir Path root) throws Exception {
		String body = "{\"entries\": [1, 2, 3]}\n";
		HttpServer server = HttpServer.create(new InetSocketAddress("127.0.0.1", 0), 0);
		server.createContext("/digest.json", ex -> {
			byte[] bytes = body.getBytes(StandardCharsets.UTF_8);
			ex.sendResponseHeaders(200, bytes.length);
			try (var os = ex.getResponseBody()) {
				os.write(bytes);
			}
		});
		server.start();
		try {
			Path dst = root.resolve("digest.json");
			Files.writeString(dst, "{\"entries\": [9, 9, 9, 9, 9, 9, 9, 9, 9, 9]}\n".repeat(50), StandardCharsets.UTF_8);
			long before = Files.size(dst);
			assertThat(before).isGreaterThan(body.length());

			ScanUpstreamSources.fetch("http://127.0.0.1:" + server.getAddress().getPort() + "/digest.json", dst, true);

			assertThat(Files.readString(dst, StandardCharsets.UTF_8)).as("no stale tail after the download")
					.isEqualTo(body);
		} finally {
			server.stop(0);
		}
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

	/**
	 * Sedoglavic switched his digest from a {@code {"NxMxP": {"rank": R}}} map to
	 * {@code {"entries": [{"format": [n,m,p], "rank": R}]}} on 2026-10-09. The map-only
	 * parser then found no shape key and reported 0 gaps / 0 better / 0 missing — silently.
	 * Both layouts must read the same, and an unreadable digest must yield NO entries (the
	 * caller then warns loudly) rather than a clean zero.
	 */
	@Test
	public void sedoglavic_digest_is_read_in_both_layouts() throws Exception {
		var mapper = tools.jackson.databind.json.JsonMapper.builder().build();
		Map<CanonicalShape, Integer> legacy = ScanUpstreamSources.sedoglavicRanks(mapper.readTree(
				"{\"description\":\"d\",\"2x2x2\":{\"rank\":7,\"exponent\":2.8},\"3x2x2\":{\"rank\":11}}"));
		Map<CanonicalShape, Integer> current = ScanUpstreamSources.sedoglavicRanks(mapper.readTree(
				"{\"description\":\"d\",\"createdAt\":\"2026-10-09\",\"entries\":[{\"format\":[2,2,2],\"rank\":7},"
						+ "{\"format\":[3,2,2],\"rank\":11,\"exponent\":2.89}]}"));
		assertThat(legacy).containsEntry(Shape.of(2, 2, 2).canonical(), 7).containsEntry(Shape.of(2, 2, 3).canonical(), 11)
				.hasSize(2);
		assertThat(current).isEqualTo(legacy);
		assertThat(ScanUpstreamSources.sedoglavicRanks(mapper.readTree("{\"description\":\"d\",\"table\":[]}")))
				.as("an unknown layout yields no entry, never a clean zero-gap scan").isEmpty();
	}

	/** The committed digest (refreshed 2026-10-13) is the new layout and reads completely. */
	@Test
	public void committed_sedoglavic_digest_reads_completely() throws Exception {
		var mapper = tools.jackson.databind.json.JsonMapper.builder().build();
		Map<CanonicalShape, Integer> ranks = ScanUpstreamSources.sedoglavicRanks(
				mapper.readTree(Files.readString(Path.of("references/sedoglavic-fmm-sota.json"), StandardCharsets.UTF_8)));
		assertThat(ranks.size()).isGreaterThan(5000);
		assertThat(ranks).containsEntry(Shape.of(2, 2, 2).canonical(), 7);
	}

	/**
	 * The Z and ZT axes (issue #20): Z from fields[] (stubs included), ZT from the stored flag
	 * or the coefficients of an explicit scheme — never from a stub, which cannot be judged.
	 */
	@Test
	public void local_z_and_zt_axes_follow_fields_and_coefficients(@TempDir Path root) throws Exception {
		// Stamped ternary integer scheme: on all three axes.
		write(root.resolve("known/section2"), "a.json",
				"{\"n\": [2, 2, 2], \"m\": 7, \"zt\": true, \"u\": [[1]], \"fields\": [\"F2\", \"F3\", \"Z\", \"Q\", \"R\", \"C\"]}");
		// Unstamped explicit integer scheme with a coefficient 2: Z yes, ZT no (read from the coefficients).
		write(root.resolve("known/section3"), "b.json",
				"{\"n\": [1, 1, 1], \"m\": 1, \"u\": [[2]], \"v\": [[1]], \"w\": [[1]], \"fields\": [\"Z\", \"Q\", \"R\", \"C\"]}");
		// A stub with Z among its fields: Z yes, ZT unknowable → not counted.
		write(root.resolve("derived/section4"), "c.json",
				"{\"n\": [2, 2, 4], \"m\": 14, \"scheme_type\": \"stub\", \"fields\": [\"Z\", \"Q\", \"R\", \"C\"]}");
		// A rational scheme: Q only.
		write(root.resolve("known/section4"), "d.json",
				"{\"n\": [4, 4, 4], \"m\": 48, \"u\": [], \"fields\": [\"Q\", \"R\", \"C\"]}");

		ScanUpstreamSources.LocalRanks local = ScanUpstreamSources.loadLocalRanksByClass(root);

		assertThat(local.q()).containsEntry(Shape.of(2, 2, 2).canonical(), 7).containsEntry(Shape.of(4, 4, 4).canonical(), 48);
		assertThat(local.z()).containsEntry(Shape.of(1, 1, 1).canonical(), 1).containsEntry(Shape.of(2, 2, 4).canonical(), 14)
				.doesNotContainKey(Shape.of(4, 4, 4).canonical());
		assertThat(local.zt()).containsEntry(Shape.of(2, 2, 2).canonical(), 7)
				.doesNotContainKey(Shape.of(1, 1, 1).canonical()).doesNotContainKey(Shape.of(2, 2, 4).canonical());
	}

	/**
	 * The machine-readable hand-off (2026-11-09, issue #23): the weekly scan writes the gap
	 * list the local closure consumes, merging the sources per (axis, shape) and classifying
	 * each row from FMM's recipe cell — a closure candidate (a recipe over formats we hold:
	 * projection / Kronecker / block sum) vs a base candidate (a cited explicit scheme we
	 * lack). Nothing in CI runs the closure, so the file is what a human/agent starts from.
	 */
	@Test
	public void gap_file_merges_sources_and_classifies_from_the_fmm_recipe() throws Exception {
		CanonicalShape s91717 = Shape.of(9, 17, 17).canonical();
		CanonicalShape s444 = Shape.of(4, 4, 4).canonical();
		CanonicalShape s7717 = Shape.of(7, 7, 17).canonical();
		CanonicalShape s555 = Shape.of(5, 5, 5).canonical();
		CanonicalShape s236 = Shape.of(2, 3, 6).canonical();
		SourceScan sed = new SourceScan(List.of(new Gap(26, s91717, 1600, 1574)), List.of(), List.of(),
				List.of(new Gap(1, s7717, 573, 572)));
		SourceScan fmm = new SourceScan(List.of(new Gap(26, s91717, 1600, 1574), new Gap(1, s444, 49, 48)),
				List.of(), List.of(), List.of(new Gap(1, s7717, 573, 572)));
		SourceScan perm = new SourceScan(List.of(new Gap(1, s444, 49, 48)), List.of(), List.of(), List.of(),
				List.of(), List.of(new Gap(Integer.MAX_VALUE, s555, -1, 93)));
		SourceScan mizopl = new SourceScan(List.of(), List.of(), List.of());
		Map<CanonicalShape, FmmMeta> meta = Map.of(
				s91717, new FmmMeta("ProjOfAlongMode(⟨9×17×18:1591⟩,2)", List.of()),
				s444, new FmmMeta("Novikov et ali (AlphaEvolve 2025) [pdf]", List.of("alphaevolve")),
				s7717, new FmmMeta("⟨7×7×8:…⟩ + ⟨7×7×9:…⟩", List.of()),
				s236, new FmmMeta("⟨1×1×2:2⟩ ⊗ ⟨2×3×3:15⟩", List.of()));

		List<GapRow> rows = ScanUpstreamSources.mergeGaps(sed, perm, fmm, mizopl, meta);

		GapRow proj = rows.stream().filter(r -> r.shape().equals(s91717)).findFirst().orElseThrow();
		assertThat(proj.bySource()).containsEntry("fmm", 1574).containsEntry("sedoglavic", 1574).hasSize(2);
		assertThat(proj.kind()).isEqualTo("closure-candidate:projection");
		assertThat(proj.gap()).isEqualTo(26);
		GapRow base = rows.stream().filter(r -> r.shape().equals(s444)).findFirst().orElseThrow();
		assertThat(base.kind()).isEqualTo("base-candidate");
		assertThat(base.bySource()).containsKeys("fmm", "perminov");
		assertThat(rows.stream().filter(r -> r.shape().equals(s7717)).findFirst().orElseThrow().kind())
				.isEqualTo("upstream-unverified");
		GapRow zt = rows.stream().filter(r -> r.axis().equals("ZT")).findFirst().orElseThrow();
		assertThat(zt.shape()).isEqualTo(s555);
		assertThat(zt.ours()).isEqualTo(-1);
		assertThat(zt.kind()).isEqualTo("coefficient-class");
		// Q axis first, biggest gap first; one row per (axis, shape)
		assertThat(rows.get(0).shape()).isEqualTo(s91717);
		assertThat(rows).hasSize(4);
		assertThat(ScanUpstreamSources.classifyKind(meta.get(s236))).isEqualTo("closure-candidate:kron");
		assertThat(ScanUpstreamSources.classifyKind(null)).isEqualTo("unclassified");

		JsonNode json = JsonMapper.builder().build()
				.readTree(ScanUpstreamSources.gapFileJson(rows, Instant.parse("2026-11-09T05:00:00Z")));
		assertThat(json.get("generated_at").asText()).isEqualTo("2026-11-09T05:00:00Z");
		assertThat(json.get("counts").get("closure-candidate:projection").asInt()).isEqualTo(1);
		assertThat(json.get("counts").get("upstream-unverified").asInt()).isEqualTo(1);
		JsonNode first = json.get("gaps").get(0);
		assertThat(first.get("shape").asText()).isEqualTo("9x17x17");
		assertThat(first.get("ours").asInt()).isEqualTo(1600);
		assertThat(first.get("upstream").asInt()).isEqualTo(1574);
		assertThat(first.get("fmm_recipe").asText()).contains("ProjOfAlongMode");
		assertThat(first.get("by_source").get("sedoglavic").asInt()).isEqualTo(1574);

		String txt = ScanUpstreamSources.shapeList(rows);
		List<String> shapeLines = txt.lines().filter(l -> !l.startsWith("#")).toList();
		// unverified and ZT rows left out; leading token = the shape, as both drivers read it
		assertThat(shapeLines).containsExactly(
				"9x17x17 1600 1574 26 closure-candidate:projection fmm,sedoglavic",
				"4x4x4 49 48 1 base-candidate fmm,perminov");
		assertThat(txt).doesNotContain("7x7x17").doesNotContain("5x5x5");
	}

	/** The shape list is read unchanged by the projection driver (and SchemeSweep's reader takes
	 *  the same leading token): a row with trailing columns is one shape, comments are skipped. */
	@Test
	public void shape_list_is_read_by_the_projection_driver_unchanged(@TempDir Path root) throws Exception {
		Path f = root.resolve("upstream-gaps.txt");
		Files.writeString(f, "# Columns: NxMxP ours upstream gap kind sources\n"
				+ "9x17x17 1600 1574 26 closure-candidate:projection fmm,sedoglavic\n"
				+ "13,14,14\n"
				+ "6 21 23\n", StandardCharsets.UTF_8);
		List<int[]> shapes = ProjectFmmGaps.readShapeFile(f.toString());
		assertThat(shapes).hasSize(3);
		assertThat(shapes.get(0)).startsWith(9, 17, 17);
		assertThat(shapes.get(1)).startsWith(13, 14, 14);
		assertThat(shapes.get(2)).startsWith(6, 21, 23);
	}
}
