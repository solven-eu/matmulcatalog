package eu.solven.matmul.docs.explore;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import eu.solven.matmul.CanonicalShape;
import eu.solven.matmul.Shape;
import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * Scan upstream catalogs for ranks that beat what we have on disk — Java port of
 * the former {@code tools/scan_upstream_sources.py}, so the GitHub-Actions
 * process is the same one anyone can replay locally with a single Maven command.
 *
 * <p>Sources:</p>
 * <ul>
 *   <li><b>Sedoglavic fmm_digest</b> — manually-curated SOTA, refreshed from
 *       GitHub into {@code references/sedoglavic-fmm-sota.json} (a versioned cache).</li>
 *   <li><b>Perminov status.json</b> — fetched fresh to a temp file.</li>
 *   <li><b>FMM Université de Lille</b> — read from the local versioned snapshot
 *       {@code references/catalogs/fmm-lille-catalog.json} (no URL refresh; updated manually).</li>
 * </ul>
 *
 * <p>Output: a Markdown report listing per-source "we should investigate" rows.
 * Non-destructive — writes no scheme files. Intended for periodic CI execution.</p>
 *
 * <pre>
 *   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.explore.ScanUpstreamSources
 *   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.explore.ScanUpstreamSources \
 *       -Dexec.args="--refresh --out generated/UPSTREAM_SCAN.md"
 * </pre>
 */
@Slf4j
public final class ScanUpstreamSources {

	private static final String SEDOGLAVIC_URL =
			"https://raw.githubusercontent.com/sedoglavic/fmm_digest/master/fmm_sota.json";
	private static final Path SEDOGLAVIC_LOCAL = Path.of("references/sedoglavic-fmm-sota.json");
	private static final String PERMINOV_URL =
			"https://raw.githubusercontent.com/dronperminov/FastMatrixMultiplication/master/schemes/status.json";
	private static final Path PERMINOV_LOCAL =
			Path.of(System.getProperty("java.io.tmpdir"), "perminov_status.json");
	private static final Path FMM_LILLE_LOCAL = Path.of("references/catalogs/fmm-lille-catalog.json");

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");
	private static final Path DEFAULT_OUT = Path.of("generated/UPSTREAM_SCAN.md");

	private static final JsonMapper MAPPER = JsonMapper.builder().build();

	private static final Pattern SHAPE_KEY = Pattern.compile("^(\\d+)x(\\d+)x(\\d+)$");

	private ScanUpstreamSources() {}

	/** A shape where an upstream beats (or only-has) what we hold. */
	private record Gap(int gap, CanonicalShape shape, int ours, int theirs) {}

	private record Missing(CanonicalShape shape, int theirs) {}

	private record SourceScan(List<Gap> gaps, List<Missing> missing, List<Gap> weBetter) {}

	public static void main(String[] args) throws IOException, InterruptedException {
		boolean refresh = false;
		Path out = DEFAULT_OUT;
		for (int i = 0; i < args.length; i++) {
			String a = args[i];
			if (a.equals("--refresh")) {
				refresh = true;
			} else if (a.equals("--out")) {
				out = Path.of(args[++i]);
			} else if (a.startsWith("--out=")) {
				out = Path.of(a.substring("--out=".length()));
			} else {
				throw new IllegalArgumentException("unknown arg: " + a);
			}
		}

		fetch(SEDOGLAVIC_URL, SEDOGLAVIC_LOCAL, refresh);
		fetch(PERMINOV_URL, PERMINOV_LOCAL, refresh);

		Map<CanonicalShape, Integer> local = loadLocalRanks();
		log.info("Loaded {} local-best shapes", local.size());

		SourceScan sed = scanSedoglavic(local);
		SourceScan perm = scanPerminov(local);
		SourceScan fmm = scanFmmLille(local);

		String md = reportMd(sed, perm, fmm);
		if (out.getParent() != null) {
			Files.createDirectories(out.getParent());
		}
		Files.writeString(out, md);
		log.info("Wrote {}", out);
		System.out.println(md);
	}

	/** Download {@code url} to {@code dst} when forced, missing, or suspiciously
	 *  small (&lt;1&nbsp;KB) — matching the Python heuristic. */
	private static void fetch(String url, Path dst, boolean force) throws IOException, InterruptedException {
		boolean stale = force
				|| !Files.exists(dst)
				|| Files.size(dst) < 1000;
		if (!stale) {
			return;
		}
		log.info("fetching {}", url);
		if (dst.getParent() != null) {
			Files.createDirectories(dst.getParent());
		}
		HttpClient client = HttpClient.newBuilder()
				.followRedirects(HttpClient.Redirect.NORMAL)
				.build();
		HttpRequest req = HttpRequest.newBuilder(URI.create(url)).GET().build();
		HttpResponse<Path> resp = client.send(req, HttpResponse.BodyHandlers.ofFile(dst));
		if (resp.statusCode() / 100 != 2) {
			throw new IOException("HTTP " + resp.statusCode() + " fetching " + url);
		}
	}

	/** Best (lowest) NON-COMMUTATIVE, Q-valid rank we hold per canonical (sorted)
	 *  shape, read from scheme CONTENT ({@code n}, {@code m}, {@code fields[]},
	 *  {@code commutative}) under {@link #SCHEMES_ROOT} — the comparison scope of
	 *  every upstream digest (FMM / Perminov / Sedoglavic are NC, char-0). */
	private static Map<CanonicalShape, Integer> loadLocalRanks() throws IOException {
		return loadLocalRanks(SCHEMES_ROOT);
	}

	/**
	 * Content-driven local best per shape. Until 2026-09-30 this parsed the rank
	 * out of the FILENAME with {@link #FILE_SHAPE_RANK} (legacy
	 * {@code {source}-{shape}_m{rank}} names); the 2026-06 rename to
	 * {@code {n}x{m}x{p}-r{rank}-{note}-{hash7}} left that regex matching 7 of
	 * ~11k files, so every upstream shape read as MISSING (not a gap) and the
	 * weekly scan reported "Strict gaps: 0" for three months — including the
	 * 2026-09 KGP LITA refresh that put ~530 shapes below the catalog. Reading
	 * {@code n}/{@code m} from the JSON (per CLAUDE.md "content-driven, never
	 * the filename") cannot drift with naming. Commutative-only schemes and
	 * schemes not valid over Q (F₂-native, C-only) are excluded so the scan
	 * never "closes" a Q gap with an incomparable scheme.
	 */
	static Map<CanonicalShape, Integer> loadLocalRanks(Path schemesRoot) throws IOException {
		Map<CanonicalShape, Integer> out = new LinkedHashMap<>();
		if (!Files.isDirectory(schemesRoot)) {
			return out;
		}
		List<Path> files;
		try (Stream<Path> walk = Files.walk(schemesRoot)) {
			files = walk.filter(p -> p.getFileName().toString().endsWith(".json")).toList();
		}
		int skipped = 0;
		for (Path p : files) {
			try {
				JsonNode root = eu.solven.matmul.catalog.SchemeIO.parseJson(p.toFile());
				JsonNode dims = root.get("n");
				JsonNode rank = root.get("m");
				if (dims == null || !dims.isArray() || dims.size() != 3 || rank == null || !rank.isInt()) {
					skipped++;
					continue;
				}
				if (root.path("commutative").asBoolean(false)) {
					continue;
				}
				List<String> fields = eu.solven.matmul.catalog.SchemeIO.fieldTags(root);
				if (!fields.isEmpty() && !fields.contains("Q")) {
					continue; // F₂-native / C-only: not comparable to the NC char-0 digests
				}
				CanonicalShape sh = Shape.of(dims.get(0).asInt(), dims.get(1).asInt(), dims.get(2).asInt()).canonical();
				out.merge(sh, rank.asInt(), Math::min);
			} catch (Exception e) {
				skipped++;
				log.debug("loadLocalRanks: skipping {}: {}", p.getFileName(), e.toString());
			}
		}
		log.info("local ranks: {} shapes from {} files ({} skipped: non-scheme JSON / unreadable)",
				out.size(), files.size(), skipped);
		return out;
	}

	private static SourceScan scanSedoglavic(Map<CanonicalShape, Integer> local) {
		JsonNode sed = readJson(SEDOGLAVIC_LOCAL);
		List<Gap> gaps = new ArrayList<>();
		List<Missing> missing = new ArrayList<>();
		List<Gap> weBetter = new ArrayList<>();
		for (Map.Entry<String, JsonNode> e : sed.properties()) {
			CanonicalShape sh = parseShapeKey(e.getKey());
			if (sh == null) {
				continue;
			}
			int sr = e.getValue().get("rank").asInt();
			Integer our = local.get(sh);
			if (our == null) {
				missing.add(new Missing(sh, sr));
			} else if (sr < our) {
				gaps.add(new Gap(our - sr, sh, our, sr));
			} else if (our < sr) {
				weBetter.add(new Gap(sr - our, sh, our, sr));
			}
		}
		return new SourceScan(gaps, missing, weBetter);
	}

	private static SourceScan scanFmmLille(Map<CanonicalShape, Integer> local) {
		List<Gap> gaps = new ArrayList<>();
		List<Missing> missing = new ArrayList<>();
		if (!Files.exists(FMM_LILLE_LOCAL)) {
			return new SourceScan(gaps, missing, List.of());
		}
		JsonNode root = readJson(FMM_LILLE_LOCAL);
		for (JsonNode e : root.path("entries")) {
			JsonNode fmt = e.get("format");
			CanonicalShape sh = Shape.of(fmt.get(0).asInt(), fmt.get(1).asInt(), fmt.get(2).asInt()).canonical();
			int r = e.get("rank").asInt();
			Integer our = local.get(sh);
			if (our == null) {
				missing.add(new Missing(sh, r));
			} else if (r < our) {
				gaps.add(new Gap(our - r, sh, our, r));
			}
		}
		return new SourceScan(gaps, missing, List.of());
	}

	private static SourceScan scanPerminov(Map<CanonicalShape, Integer> local) {
		JsonNode status = readJson(PERMINOV_LOCAL);
		List<Gap> gaps = new ArrayList<>();
		List<Missing> missing = new ArrayList<>();
		for (Map.Entry<String, JsonNode> e : status.properties()) {
			CanonicalShape sh = parseShapeKey(e.getKey());
			if (sh == null) {
				continue;
			}
			JsonNode ranks = e.getValue().path("ranks");
			int best = Integer.MAX_VALUE;
			for (JsonNode r : ranks) {
				if (r.isInt() || r.isIntegralNumber()) {
					best = Math.min(best, r.asInt());
				}
			}
			if (best == Integer.MAX_VALUE) {
				continue;
			}
			Integer our = local.get(sh);
			if (our == null) {
				missing.add(new Missing(sh, best));
			} else if (best < our) {
				gaps.add(new Gap(our - best, sh, our, best));
			}
		}
		return new SourceScan(gaps, missing, List.of());
	}

	private static CanonicalShape parseShapeKey(String key) {
		Matcher m = SHAPE_KEY.matcher(key);
		if (!m.matches()) {
			return null;
		}
		return Shape.of(Integer.parseInt(m.group(1)), Integer.parseInt(m.group(2)), Integer.parseInt(m.group(3)))
				.canonical();
	}

	private static JsonNode readJson(Path p) {
		try {
			return MAPPER.readTree(Files.readString(p));
		} catch (IOException ex) {
			throw new UncheckedIOException("reading " + p, ex);
		}
	}

	/** Descending by gap, then by shape — the "biggest opportunities first" order. */
	private static final Comparator<Gap> BY_GAP_DESC = Comparator
			.comparingInt(Gap::gap).reversed()
			.thenComparing(g -> g.shape().toString());

	private static String reportMd(SourceScan sed, SourceScan perm, SourceScan fmm) {
		StringBuilder out = new StringBuilder();
		out.append("# Upstream source scan\n\n");
		out.append("Generated by `eu.solven.matmul.docs.explore.ScanUpstreamSources`.\n\n");

		out.append("## Sedoglavic fmm_digest\n\n");
		out.append("- Strict gaps (Sedoglavic better): **").append(sed.gaps().size()).append("**\n");
		out.append("- Missing from us (Sedoglavic only): **").append(sed.missing().size()).append("**\n");
		out.append("- We-better: **").append(sed.weBetter().size()).append("**\n");
		appendGapTable(out, "Top 20 strict gaps", "Sedoglavic", sed.gaps());

		out.append("\n## FMM Université de Lille (local snapshot)\n\n");
		out.append("- Strict gaps: **").append(fmm.gaps().size()).append("**\n");
		out.append("- Missing from us: **").append(fmm.missing().size()).append("**\n");
		appendGapTable(out, "Top 20 FMM strict gaps", "FMM", fmm.gaps());

		out.append("\n## Perminov status.json\n\n");
		out.append("- Strict gaps: **").append(perm.gaps().size()).append("**\n");
		out.append("- Missing from us: **").append(perm.missing().size()).append("**\n");
		appendGapTable(out, "Top 20 Perminov strict gaps", "Perminov", perm.gaps());

		return out.toString();
	}

	private static void appendGapTable(StringBuilder out, String title, String theirCol, List<Gap> gaps) {
		if (gaps.isEmpty()) {
			return;
		}
		out.append("\n### ").append(title).append("\n\n");
		out.append("| shape | ours | ").append(theirCol).append(" | gap |\n");
		out.append("|---|--:|--:|--:|\n");
		gaps.stream().sorted(BY_GAP_DESC).limit(20).forEach(g -> {
			CanonicalShape s = g.shape();
			out.append("| ⟨").append(s.n()).append(',').append(s.m()).append(',').append(s.p())
					.append("⟩ | ").append(g.ours()).append(" | ").append(g.theirs())
					.append(" | ").append(g.gap()).append(" |\n");
		});
	}
}
