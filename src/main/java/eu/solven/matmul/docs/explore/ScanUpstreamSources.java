package eu.solven.matmul.docs.explore;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import eu.solven.matmul.CanonicalShape;
import eu.solven.matmul.Shape;

import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ArrayNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Scan upstream catalogs for ranks that beat what we have on disk — Java port of
 * the former {@code tools/scan_upstream_sources.py}, so the GitHub-Actions
 * process is the same one anyone can replay locally with a single Maven command.
 *
 * <p>Sources:</p>
 * <ul>
 *   <li><b>Sedoglavic fmm_digest</b> — manually-curated SOTA, refreshed from
 *       GitHub into {@code references/sedoglavic-fmm-sota.json} (a versioned cache).</li>
 *   <li><b>Perminov status.json</b> — fetched fresh to a temp file. Compared on three
 *       axes: the char-0 rank (min over his Q/Z/ZT, against our Q-valid best), the
 *       integer rank ({@code Z}) and the ternary-integer rank ({@code ZT}, explicit
 *       schemes only — a stub's ternarity is not known without a replay, exactly as the
 *       catalog's own {@code zt} flag). Until 2026-11-02 only the first axis existed, so
 *       a ZT record we lacked at a shape whose Q rank we matched was invisible
 *       (the ten mizopl ZT records of issue #20).</li>
 *   <li><b>mizopl / we*star status.json</b> — the same Perminov layout, fetched fresh;
 *       the origin of the schemes Perminov mirrors under {@code known/mizopl/}.</li>
 *   <li><b>FMM Université de Lille</b> — read from the local versioned snapshot
 *       {@code references/catalogs/fmm-lille-catalog.json} (no URL refresh; updated manually).</li>
 * </ul>
 *
 * <p>Output: a Markdown report listing per-source "we should investigate" rows, and —
 * with {@code --gaps-out} / {@code --shapes-out} — the machine-readable hand-off to the
 * closure: one JSON row per (axis, shape) merging every source that beats us, classified
 * from FMM's recipe cell ({@link #classifyKind}: closure candidate vs base candidate vs
 * upstream-unverified), and a plain shape list whose leading token is what
 * {@code SchemeSweep --shape-file} and {@code ProjectFmmGaps --shapes-file} read. The scan
 * only <em>detects</em>: nothing in CI runs the closure (RESYNC_PLAYBOOK step 4 — CPU-heavy,
 * writes schemes), so issue #23's 44 FMM rows sat open until a human ran it; the files make
 * that hand-off a download plus two commands.
 * Non-destructive — writes no scheme files. Intended for periodic CI execution.</p>
 *
 * <pre>
 *   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.explore.ScanUpstreamSources
 *   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.explore.ScanUpstreamSources \
 *       -Dexec.args="--refresh --out generated/UPSTREAM_SCAN.md \
 *                    --gaps-out generated/upstream-gaps.json --shapes-out generated/upstream-gaps.txt"
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
	private static final String MIZOPL_URL =
			"https://raw.githubusercontent.com/mizopl/we-star/main/solutions/fast-matrix-multiplication/status.json";
	private static final Path MIZOPL_LOCAL =
			Path.of(System.getProperty("java.io.tmpdir"), "mizopl_status.json");
	private static final Path FMM_LILLE_LOCAL = Path.of("references/catalogs/fmm-lille-catalog.json");

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");
	private static final Path DEFAULT_OUT = Path.of("generated/UPSTREAM_SCAN.md");

	private static final JsonMapper MAPPER = JsonMapper.builder().build();

	private static final Pattern SHAPE_KEY = Pattern.compile("^(\\d+)x(\\d+)x(\\d+)$");

	private ScanUpstreamSources() {}

	/** A shape where an upstream beats (or only-has) what we hold. */
	record Gap(int gap, CanonicalShape shape, int ours, int theirs) {}

	record Missing(CanonicalShape shape, int theirs) {}

	/** {@code unverified}: FMM rows whose index rank is not backed by FMM's own artifact
	 *  (manual audit, {@link eu.solven.matmul.docs.verify.FmmArtifactAudit}) — reported,
	 *  never counted as a strict gap. Empty for the other sources. {@code zGaps} /
	 *  {@code ztGaps}: the integer and ternary-integer axes (Perminov-layout sources only). */
	record SourceScan(List<Gap> gaps, List<Missing> missing, List<Gap> weBetter, List<Gap> unverified,
			List<Gap> zGaps, List<Gap> ztGaps) {
		SourceScan(List<Gap> gaps, List<Missing> missing, List<Gap> weBetter) {
			this(gaps, missing, weBetter, List.of(), List.of(), List.of());
		}

		SourceScan(List<Gap> gaps, List<Missing> missing, List<Gap> weBetter, List<Gap> unverified) {
			this(gaps, missing, weBetter, unverified, List.of(), List.of());
		}
	}

	/** Our best rank per canonical shape on the three comparison axes (see {@link #loadLocalRanksByClass}). */
	public record LocalRanks(Map<CanonicalShape, Integer> q, Map<CanonicalShape, Integer> z,
			Map<CanonicalShape, Integer> zt) {}

	public static void main(String[] args) throws IOException, InterruptedException {
		boolean refresh = false;
		Path out = DEFAULT_OUT;
		Path gapsOut = null;
		Path shapesOut = null;
		for (int i = 0; i < args.length; i++) {
			String a = args[i];
			if (a.equals("--refresh")) {
				refresh = true;
			} else if (a.equals("--out")) {
				out = Path.of(args[++i]);
			} else if (a.startsWith("--out=")) {
				out = Path.of(a.substring("--out=".length()));
			} else if (a.equals("--gaps-out")) {
				gapsOut = Path.of(args[++i]);
			} else if (a.startsWith("--gaps-out=")) {
				gapsOut = Path.of(a.substring("--gaps-out=".length()));
			} else if (a.equals("--shapes-out")) {
				shapesOut = Path.of(args[++i]);
			} else if (a.startsWith("--shapes-out=")) {
				shapesOut = Path.of(a.substring("--shapes-out=".length()));
			} else {
				throw new IllegalArgumentException("unknown arg: " + a);
			}
		}

		fetch(SEDOGLAVIC_URL, SEDOGLAVIC_LOCAL, refresh);
		fetch(PERMINOV_URL, PERMINOV_LOCAL, refresh);
		fetch(MIZOPL_URL, MIZOPL_LOCAL, refresh);

		LocalRanks local = loadLocalRanksByClass(SCHEMES_ROOT);
		log.info("Loaded {} local-best shapes ({} with an integer scheme, {} with a ternary one)", local.q().size(),
				local.z().size(), local.zt().size());

		SourceScan sed = scanSedoglavic(local.q());
		SourceScan perm = scanPerminovLayout(PERMINOV_LOCAL, local);
		SourceScan mizopl = scanPerminovLayout(MIZOPL_LOCAL, local);
		SourceScan fmm = scanFmmLille(local.q());

		// The closure hand-off: every row, merged across sources and classified from FMM's recipe.
		Map<CanonicalShape, FmmMeta> meta = Files.exists(FMM_LILLE_LOCAL) ? fmmMeta(readJson(FMM_LILLE_LOCAL))
				: Map.of();
		List<GapRow> rows = mergeGaps(sed, perm, fmm, mizopl, meta);
		if (gapsOut != null) {
			writeText(gapsOut, gapFileJson(rows, Instant.now()));
		}
		if (shapesOut != null) {
			writeText(shapesOut, shapeList(rows));
		}

		String md = reportMd(sed, perm, fmm, mizopl, rows, gapsOut, shapesOut);
		writeText(out, md);
		System.out.println(md);
	}

	private static void writeText(Path out, String text) throws IOException {
		if (out.getParent() != null) {
			Files.createDirectories(out.getParent());
		}
		Files.writeString(out, text);
		log.info("Wrote {}", out);
	}

	/** Download {@code url} to {@code dst} when forced, missing, or suspiciously
	 *  small (&lt;1&nbsp;KB) — matching the Python heuristic.
	 *
	 *  <p>The destination is TRUNCATED before the body is written.
	 *  {@code BodyHandlers.ofFile(dst)} alone opens the file with {@code CREATE, WRITE} only,
	 *  so an existing file longer than the download kept its old tail after the new
	 *  content — a JSON that no longer parsed ("Unrecognized token 'mat'": the old
	 *  {@code "format"} line, split). Seen 2026-10-19 on a Windows checkout, where the CRLF
	 *  working copy is longer than the LF download; on CI it only waited for the day the
	 *  upstream digest shrinks, when the weekly scan would have crashed instead of
	 *  reporting (it did, 2026-10-26 and 2026-11-02, before the fix was merged).</p> */
	static void fetch(String url, Path dst, boolean force) throws IOException, InterruptedException {
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
		HttpResponse<Path> resp = client.send(req, HttpResponse.BodyHandlers.ofFile(dst,
				StandardOpenOption.CREATE, StandardOpenOption.WRITE, StandardOpenOption.TRUNCATE_EXISTING));
		if (resp.statusCode() / 100 != 2) {
			throw new IOException("HTTP " + resp.statusCode() + " fetching " + url);
		}
	}

	/**
	 * Content-driven local best per shape (the char-0 / Q axis). Until 2026-09-30 this
	 * parsed the rank out of the FILENAME with a legacy {@code {source}-{shape}_m{rank}}
	 * regex; the 2026-06 rename to {@code {n}x{m}x{p}-r{rank}-{note}-{hash7}} left that
	 * regex matching 7 of ~11k files, so every upstream shape read as MISSING (not a
	 * gap) and the weekly scan reported "Strict gaps: 0" for three months — including
	 * the 2026-09 KGP LITA refresh that put ~530 shapes below the catalog. Reading
	 * {@code n}/{@code m} from the JSON (per CLAUDE.md "content-driven, never the
	 * filename") cannot drift with naming. Commutative-only schemes and schemes not
	 * valid over Q (F₂-native, C-only) are excluded so the scan never "closes" a Q gap
	 * with an incomparable scheme.
	 */
	static Map<CanonicalShape, Integer> loadLocalRanks(Path schemesRoot) throws IOException {
		return loadLocalRanksByClass(schemesRoot).q();
	}

	/**
	 * {@link #loadLocalRanks} on all three axes: {@code q} every Q-valid NC scheme,
	 * {@code z} those with {@code Z} among their {@code fields[]} (stubs included — field
	 * membership is stamped on them), {@code zt} the explicit integer schemes whose every
	 * coefficient is in {-1,0,1} ({@code SchemeIO.ternaryOf}: the stored {@code zt} flag,
	 * else the coefficients; a stub has neither, so it never counts as ternary — the same
	 * reading as the catalog manifest and the SPA's ZT selector).
	 */
	static LocalRanks loadLocalRanksByClass(Path schemesRoot) throws IOException {
		Map<CanonicalShape, Integer> q = new LinkedHashMap<>();
		Map<CanonicalShape, Integer> z = new LinkedHashMap<>();
		Map<CanonicalShape, Integer> zt = new LinkedHashMap<>();
		if (!Files.isDirectory(schemesRoot)) {
			return new LocalRanks(q, z, zt);
		}
		List<Path> files;
		try (Stream<Path> walk = Files.walk(schemesRoot)) {
			files = walk.filter(p -> p.getFileName().toString().endsWith(".json")).toList();
		}
		int skipped = 0, ternaryReads = 0;
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
				q.merge(sh, rank.asInt(), Math::min);
				if (fields.contains("Z")) {
					z.merge(sh, rank.asInt(), Math::min);
					Integer ztBest = zt.get(sh);
					// Only a candidate for the ZT best needs the (possibly matrix-reading) verdict.
					if (ztBest == null || rank.asInt() < ztBest) {
						if (eu.solven.matmul.catalog.SchemeIO.readZT(root) == null
								&& !eu.solven.matmul.catalog.SchemeIO.isStub(root)) {
							ternaryReads++;
						}
						if (Boolean.TRUE.equals(eu.solven.matmul.catalog.SchemeIO.ternaryOf(root))) {
							zt.merge(sh, rank.asInt(), Math::min);
						}
					}
				}
			} catch (Exception e) {
				skipped++;
				log.debug("loadLocalRanks: skipping {}: {}", p.getFileName(), e.toString());
			}
		}
		log.info("local ranks: {} shapes from {} files ({} skipped: non-scheme JSON / unreadable; {} integer, "
				+ "{} ternary; {} coefficient reads for an unstamped zt)", q.size(), files.size(), skipped, z.size(),
				zt.size(), ternaryReads);
		return new LocalRanks(q, z, zt);
	}

	private static SourceScan scanSedoglavic(Map<CanonicalShape, Integer> local) {
		JsonNode sed = readJson(SEDOGLAVIC_LOCAL);
		List<Gap> gaps = new ArrayList<>();
		List<Missing> missing = new ArrayList<>();
		List<Gap> weBetter = new ArrayList<>();
		Map<CanonicalShape, Integer> ranks = sedoglavicRanks(sed);
		if (ranks.isEmpty()) {
			// A digest the parser does not understand must not read as "no gaps".
			log.warn("[layout] {} yields no (format, rank) entry — the digest layout changed; this scan is blind "
					+ "to Sedoglavic until the parser is updated", SEDOGLAVIC_LOCAL);
			System.out.println("::warning title=Sedoglavic digest layout changed::ScanUpstreamSources parsed 0 entries");
		}
		// Sedoglavic's digest IS the FMM-Lille index, so the artifact audit applies here too.
		eu.solven.matmul.docs.verify.FmmArtifactAudit audit = loadAudit();
		List<Gap> unverified = new ArrayList<>();
		for (Map.Entry<CanonicalShape, Integer> e : ranks.entrySet()) {
			CanonicalShape sh = e.getKey();
			int sr = e.getValue();
			Integer our = local.get(sh);
			if (our == null) {
				missing.add(new Missing(sh, sr));
			} else if (sr < our) {
				Gap gap = new Gap(our - sr, sh, our, sr);
				if (audit.coveringClass(sh.n(), sh.m(), sh.p(), sr).isPresent()) unverified.add(gap);
				else gaps.add(gap);
			} else if (our < sr) {
				weBetter.add(new Gap(sr - our, sh, our, sr));
			}
		}
		return new SourceScan(gaps, missing, weBetter, unverified);
	}

	private static eu.solven.matmul.docs.verify.FmmArtifactAudit loadAudit() {
		try {
			return eu.solven.matmul.docs.verify.FmmArtifactAudit.load();
		} catch (IOException ex) {
			throw new UncheckedIOException("reading the FMM artifact audit", ex);
		}
	}

	/**
	 * Sedoglavic's {@code fmm_digest} in either layout: the original map
	 * {@code {"NxMxP": {"rank": R, …}}} and, since 2026-10-09 (file {@code createdAt}),
	 * {@code {"entries": [{"format": [n,m,p], "rank": R, …}]}}. The switch went unnoticed
	 * for a refresh: the map-only parser found no shape key and reported 0 gaps / 0
	 * better / 0 missing — the same silent shape as the June filename-regex regression.
	 */
	static Map<CanonicalShape, Integer> sedoglavicRanks(JsonNode sed) {
		Map<CanonicalShape, Integer> out = new LinkedHashMap<>();
		JsonNode entries = sed.get("entries");
		if (entries != null && entries.isArray()) {
			for (JsonNode e : entries) {
				JsonNode fmt = e.get("format");
				JsonNode rank = e.get("rank");
				if (fmt == null || !fmt.isArray() || fmt.size() != 3 || rank == null || !rank.isInt()) continue;
				out.put(Shape.of(fmt.get(0).asInt(), fmt.get(1).asInt(), fmt.get(2).asInt()).canonical(), rank.asInt());
			}
			return out;
		}
		for (Map.Entry<String, JsonNode> e : sed.properties()) {
			CanonicalShape sh = parseShapeKey(e.getKey());
			if (sh == null) continue;
			JsonNode rank = e.getValue().get("rank");
			if (rank != null && rank.isInt()) out.put(sh, rank.asInt());
		}
		return out;
	}

	private static SourceScan scanFmmLille(Map<CanonicalShape, Integer> local) {
		List<Gap> gaps = new ArrayList<>();
		List<Missing> missing = new ArrayList<>();
		if (!Files.exists(FMM_LILLE_LOCAL)) {
			return new SourceScan(gaps, missing, List.of());
		}
		JsonNode root = readJson(FMM_LILLE_LOCAL);
		// The artifact audit (same reading as FmmCrossCheck): an index rank FMM's own
		// published artifact does not back is an upstream claim, not a gap to alert on —
		// before this (2026-10-13) the weekly scan would have opened an `upstream-gap`
		// issue every Monday for the same thirteen audited rows.
		eu.solven.matmul.docs.verify.FmmArtifactAudit audit = loadAudit();
		List<Gap> unverified = new ArrayList<>();
		for (JsonNode e : root.path("entries")) {
			JsonNode fmt = e.get("format");
			CanonicalShape sh = Shape.of(fmt.get(0).asInt(), fmt.get(1).asInt(), fmt.get(2).asInt()).canonical();
			int r = e.get("rank").asInt();
			Integer our = local.get(sh);
			if (our == null) {
				missing.add(new Missing(sh, r));
			} else if (r < our) {
				Gap gap = new Gap(our - r, sh, our, r);
				if (audit.coveringClass(sh.n(), sh.m(), sh.p(), r).isPresent()) {
					unverified.add(gap);
				} else {
					gaps.add(gap);
				}
			}
		}
		return new SourceScan(gaps, missing, List.of(), unverified);
	}

	/**
	 * A Perminov-layout {@code status.json} ({@code {"NxMxP": {"ranks": {"Q": r, "Z": r,
	 * "ZT": r}, …}}}) on the three axes. Missing status files (an optional source that could
	 * not be fetched) read as an empty scan, with a warning.
	 */
	private static SourceScan scanPerminovLayout(Path statusFile, LocalRanks local) {
		List<Gap> gaps = new ArrayList<>();
		List<Missing> missing = new ArrayList<>();
		List<Gap> zGaps = new ArrayList<>();
		List<Gap> ztGaps = new ArrayList<>();
		if (!Files.exists(statusFile)) {
			log.warn("{} is absent — this source reads as empty", statusFile);
			return new SourceScan(gaps, missing, List.of(), List.of(), zGaps, ztGaps);
		}
		JsonNode status = readJson(statusFile);
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
			Integer our = local.q().get(sh);
			if (our == null) {
				missing.add(new Missing(sh, best));
			} else if (best < our) {
				gaps.add(new Gap(our - best, sh, our, best));
			}
			axisGap(ranks.get("Z"), local.z().get(sh), sh, zGaps);
			axisGap(ranks.get("ZT"), local.zt().get(sh), sh, ztGaps);
		}
		return new SourceScan(gaps, missing, List.of(), List.of(), zGaps, ztGaps);
	}

	/** A strict gap on one coefficient-class axis; a shape where we hold nothing of that class is a gap too. */
	private static void axisGap(JsonNode theirs, Integer ours, CanonicalShape sh, List<Gap> into) {
		if (theirs == null || !theirs.isIntegralNumber()) {
			return;
		}
		int t = theirs.asInt();
		if (ours == null || t < ours) {
			into.add(new Gap(ours == null ? Integer.MAX_VALUE : ours - t, sh, ours == null ? -1 : ours, t));
		}
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

	// ---- the machine-readable gap file (closure input) -------------------------------

	/** FMM's own account of a format — the index row's recipe cell ({@code ProjOfAlongMode(⟨9×17×18:1591⟩,2)},
	 *  a block sum, a Kronecker product) or a citation when the row is an imported explicit scheme —
	 *  plus its biblio keys. Read from the committed digest ({@code SyncReferenceCatalogs --fmm}). */
	record FmmMeta(String description, List<String> references) {}

	static Map<CanonicalShape, FmmMeta> fmmMeta(JsonNode fmmDigest) {
		Map<CanonicalShape, FmmMeta> out = new LinkedHashMap<>();
		for (JsonNode e : fmmDigest.path("entries")) {
			JsonNode fmt = e.get("format");
			if (fmt == null || !fmt.isArray() || fmt.size() != 3) {
				continue;
			}
			CanonicalShape sh = Shape.of(fmt.get(0).asInt(), fmt.get(1).asInt(), fmt.get(2).asInt()).canonical();
			JsonNode d = e.get("description");
			String desc = d == null || d.isNull() ? null : d.asText();
			List<String> refs = new ArrayList<>();
			for (JsonNode r : e.path("references")) {
				refs.add(r.asText());
			}
			out.put(sh, new FmmMeta(desc, refs));
		}
		return out;
	}

	static final String KIND_CLOSURE = "closure-candidate";
	static final String KIND_BASE = "base-candidate";
	static final String KIND_UNVERIFIED = "upstream-unverified";
	static final String KIND_MISSING = "missing";
	static final String KIND_AXIS = "coefficient-class";
	static final String KIND_UNCLASSIFIED = "unclassified";

	/**
	 * What kind of work a Q-axis gap row calls for, read off FMM's recipe cell. A recipe over
	 * other formats ({@code ⟨…⟩} tokens: a projection, a Kronecker product, a block sum, a
	 * serendipitous product) is a <b>closure candidate</b>: FMM reached the rank with operators
	 * we also run (RESYNC_PLAYBOOK step 4) from ingredients we hold or derive — the 2026-11
	 * wave (61 formats re-priced by {@code ProjOfAlongMode}, issue #23) was entirely this kind,
	 * and ours sat at exactly the unprojected parents. A citation (no format token) is a
	 * <b>base candidate</b>: FMM holds an explicit scheme we lack, to import as a reaction base
	 * (fmm-gap skill §3), never swept for. No FMM row (or an empty cell: the digest predates
	 * the 2026-11-09 scraper fix) → unclassified.
	 */
	static String classifyKind(FmmMeta meta) {
		if (meta == null || meta.description() == null || meta.description().isBlank()) {
			return KIND_UNCLASSIFIED;
		}
		String d = meta.description();
		if (!d.contains("⟨")) {
			return KIND_BASE;
		}
		String lower = d.toLowerCase(Locale.ROOT);
		if (lower.contains("projof") || lower.contains("projection")) {
			return KIND_CLOSURE + ":projection";
		}
		if (lower.contains("serendip")) {
			return KIND_CLOSURE + ":serendipitous";
		}
		if (d.contains("⊗") || lower.contains("otimes") || lower.contains("kronecker")) {
			return KIND_CLOSURE + ":kron";
		}
		if (d.contains("+")) {
			return KIND_CLOSURE + ":block-sum";
		}
		return KIND_CLOSURE;
	}

	/** One row of the gap file: a (axis, shape) where some upstream is strictly below us (or holds
	 *  what we lack), the best upstream rank, the sources publishing one and the work it calls for. */
	record GapRow(CanonicalShape shape, String axis, int ours, int upstream, Map<String, Integer> bySource,
			String kind, FmmMeta fmm) {
		int gap() {
			return ours < 0 ? Integer.MAX_VALUE : ours - upstream;
		}

		String shapeKey() {
			return shape.n() + "x" + shape.m() + "x" + shape.p();
		}
	}

	/** Q axis first, then the biggest gap, then the shape — the closure's "best opportunity first". */
	private static final Comparator<GapRow> BY_ROW = Comparator
			.comparingInt((GapRow r) -> r.axis().equals("Q") ? 0 : r.axis().equals("Z") ? 1 : 2)
			.thenComparing(Comparator.comparingInt(GapRow::gap).reversed())
			.thenComparing(GapRow::shapeKey);

	/**
	 * Every row of every source, merged per (axis, shape): strict Q gaps (classified from the FMM
	 * recipe), the FMM-index rows the artifact audit marks unverified, shapes we lack, and the
	 * Z / ZT coefficient-class gaps of the Perminov-layout sources. A shape that is strict for one
	 * source and audited-unverified for another keeps the strict reading (a sync skew between the
	 * two FMM-index digests, not a verdict).
	 */
	static List<GapRow> mergeGaps(SourceScan sed, SourceScan perm, SourceScan fmm, SourceScan mizopl,
			Map<CanonicalShape, FmmMeta> meta) {
		Map<String, GapRow> rows = new LinkedHashMap<>();
		addGaps(rows, "Q", "sedoglavic", sed.gaps(), meta, null);
		addGaps(rows, "Q", "fmm", fmm.gaps(), meta, null);
		addGaps(rows, "Q", "perminov", perm.gaps(), meta, null);
		addGaps(rows, "Q", "mizopl", mizopl.gaps(), meta, null);
		addGaps(rows, "Q", "sedoglavic", sed.unverified(), meta, KIND_UNVERIFIED);
		addGaps(rows, "Q", "fmm", fmm.unverified(), meta, KIND_UNVERIFIED);
		addMissing(rows, "sedoglavic", sed.missing(), meta);
		addMissing(rows, "fmm", fmm.missing(), meta);
		addMissing(rows, "perminov", perm.missing(), meta);
		addMissing(rows, "mizopl", mizopl.missing(), meta);
		addGaps(rows, "Z", "perminov", perm.zGaps(), meta, KIND_AXIS);
		addGaps(rows, "Z", "mizopl", mizopl.zGaps(), meta, KIND_AXIS);
		addGaps(rows, "ZT", "perminov", perm.ztGaps(), meta, KIND_AXIS);
		addGaps(rows, "ZT", "mizopl", mizopl.ztGaps(), meta, KIND_AXIS);
		List<GapRow> out = new ArrayList<>(rows.values());
		out.sort(BY_ROW);
		return out;
	}

	private static void addGaps(Map<String, GapRow> rows, String axis, String source, List<Gap> gaps,
			Map<CanonicalShape, FmmMeta> meta, String fixedKind) {
		for (Gap g : gaps) {
			String key = axis + " " + g.shape();
			GapRow prev = rows.get(key);
			Map<String, Integer> by = prev == null ? new TreeMap<>() : new TreeMap<>(prev.bySource());
			by.put(source, g.theirs());
			int upstream = prev == null ? g.theirs() : Math.min(prev.upstream(), g.theirs());
			String kind;
			if (fixedKind == null) {
				kind = classifyKind(meta.get(g.shape())); // a strict reading wins over an unverified one
			} else {
				kind = prev == null ? fixedKind : prev.kind();
			}
			rows.put(key, new GapRow(g.shape(), axis, g.ours(), upstream, by, kind, meta.get(g.shape())));
		}
	}

	private static void addMissing(Map<String, GapRow> rows, String source, List<Missing> missing,
			Map<CanonicalShape, FmmMeta> meta) {
		for (Missing m : missing) {
			String key = "Q " + m.shape();
			GapRow prev = rows.get(key);
			Map<String, Integer> by = prev == null ? new TreeMap<>() : new TreeMap<>(prev.bySource());
			by.put(source, m.theirs());
			int upstream = prev == null ? m.theirs() : Math.min(prev.upstream(), m.theirs());
			rows.put(key, new GapRow(m.shape(), "Q", -1, upstream, by, KIND_MISSING, meta.get(m.shape())));
		}
	}

	/** The JSON gap file: provenance, a legend of the kinds, per-kind counts and one object per row. */
	static String gapFileJson(List<GapRow> rows, Instant now) {
		ObjectNode root = MAPPER.createObjectNode();
		root.put("generated_by", "eu.solven.matmul.docs.explore.ScanUpstreamSources");
		root.put("generated_at", now.toString());
		root.put("axes", "Q = non-commutative characteristic-0 rank (our Q-valid best vs the FMM/Sedoglavic index "
				+ "and the Perminov-layout min over Q/Z/ZT); Z / ZT = Perminov-layout coefficient-class axes "
				+ "(best integer / ternary-integer explicit scheme; ours null = we hold none)");
		ObjectNode kinds = root.putObject("kinds");
		kinds.put(KIND_CLOSURE + "[:mechanism]", "FMM reached the rank with a recipe over other formats (projection"
				+ " / kron / block-sum / serendipitous) — run the closure over the .txt (RESYNC_PLAYBOOK step 4)");
		kinds.put(KIND_BASE, "FMM cites an explicit scheme we lack — import it as a reaction base (fmm-gap skill §3)");
		kinds.put(KIND_UNVERIFIED, "FMM index rank not backed by FMM's own artifact (references/fmm-artifact-audit.json)"
				+ " — report upstream, not a target");
		kinds.put(KIND_MISSING, "an upstream holds a shape we have no Q-valid scheme for");
		kinds.put(KIND_AXIS, "Z / ZT axis: the Q rank may be matched, the integer / ternary-integer explicit scheme is not");
		kinds.put(KIND_UNCLASSIFIED, "no FMM recipe for the shape (Perminov-layout source only, or a digest without"
				+ " descriptions)");
		Map<String, Integer> counts = new TreeMap<>();
		for (GapRow r : rows) {
			counts.merge(r.kind(), 1, Integer::sum);
		}
		ObjectNode countsNode = root.putObject("counts");
		counts.forEach((k, v) -> countsNode.put(k, v.intValue()));
		ArrayNode arr = root.putArray("gaps");
		for (GapRow r : rows) {
			ObjectNode o = arr.addObject();
			o.put("shape", r.shapeKey());
			o.putArray("n").add(r.shape().n()).add(r.shape().m()).add(r.shape().p());
			o.put("axis", r.axis());
			if (r.ours() < 0) {
				o.putNull("ours");
			} else {
				o.put("ours", r.ours());
			}
			o.put("upstream", r.upstream());
			if (r.gap() == Integer.MAX_VALUE) {
				o.putNull("gap");
			} else {
				o.put("gap", r.gap());
			}
			ObjectNode by = o.putObject("by_source");
			r.bySource().forEach((s, v) -> by.put(s, v.intValue()));
			o.put("kind", r.kind());
			if (r.fmm() != null) {
				if (r.fmm().description() == null) {
					o.putNull("fmm_recipe");
				} else {
					o.put("fmm_recipe", r.fmm().description());
				}
				ArrayNode refs = o.putArray("fmm_references");
				r.fmm().references().forEach(refs::add);
			}
		}
		return MAPPER.writerWithDefaultPrettyPrinter().writeValueAsString(root) + "\n";
	}

	/**
	 * The closure drivers' input: one actionable Q-axis row per line, leading token {@code NxMxP}
	 * — the only token {@code SchemeSweep --shape-file} reads, and what
	 * {@code ProjectFmmGaps --shapes-file} reads too — then ours, best upstream, gap, kind and
	 * sources. Left out: upstream-unverified rows (report-to-FMM material, not targets) and the
	 * Z / ZT rows (a Q sweep cannot address a coefficient-class gap; they stay in the JSON).
	 */
	static String shapeList(List<GapRow> rows) {
		StringBuilder sb = new StringBuilder();
		sb.append("# Upstream gaps - closure input, written by ScanUpstreamSources --shapes-out.\n");
		sb.append("# Columns: NxMxP ours upstream gap kind sources. Only the leading token is read by\n");
		sb.append("# SchemeSweep --shape-file / ProjectFmmGaps --shapes-file; the rest is for the human.\n");
		for (GapRow r : rows) {
			if (!r.axis().equals("Q") || r.kind().equals(KIND_UNVERIFIED)) {
				continue;
			}
			sb.append(r.shapeKey()).append(' ')
					.append(r.ours() < 0 ? "-" : String.valueOf(r.ours())).append(' ')
					.append(r.upstream()).append(' ')
					.append(r.gap() == Integer.MAX_VALUE ? "inf" : String.valueOf(r.gap())).append(' ')
					.append(r.kind()).append(' ')
					.append(String.join(",", r.bySource().keySet())).append('\n');
		}
		return sb.toString();
	}

	/** Descending by gap, then by shape — the "biggest opportunities first" order. */
	private static final Comparator<Gap> BY_GAP_DESC = Comparator
			.comparingInt(Gap::gap).reversed()
			.thenComparing(g -> g.shape().toString());

	private static String reportMd(SourceScan sed, SourceScan perm, SourceScan fmm, SourceScan mizopl,
			List<GapRow> rows, Path gapsOut, Path shapesOut) {
		StringBuilder out = new StringBuilder();
		out.append("# Upstream source scan\n\n");
		out.append("Generated by `eu.solven.matmul.docs.explore.ScanUpstreamSources`.\n\n");

		// The hand-off summary. (Deliberately no "Strict gaps … **N**" line: the workflow's
		// issue step sums those across the per-source sections.)
		out.append("## Closure input\n\n");
		Map<String, Integer> counts = new TreeMap<>();
		rows.forEach(r -> counts.merge(r.kind(), 1, Integer::sum));
		long actionable = rows.stream().filter(r -> r.axis().equals("Q") && !r.kind().equals(KIND_UNVERIFIED))
				.count();
		out.append("- Rows: ").append(rows.size()).append(" (").append(actionable)
				.append(" actionable on the Q axis, listed in the shape file)\n");
		out.append("- By kind: ");
		if (counts.isEmpty()) {
			out.append("none");
		} else {
			out.append(String.join(", ", counts.entrySet().stream().map(e -> e.getKey() + " " + e.getValue()).toList()));
		}
		out.append('\n');
		if (gapsOut != null || shapesOut != null) {
			out.append("- Files: ");
			if (gapsOut != null) {
				out.append('`').append(gapsOut.toString().replace('\\', '/')).append('`');
			}
			if (gapsOut != null && shapesOut != null) {
				out.append(", ");
			}
			if (shapesOut != null) {
				out.append('`').append(shapesOut.toString().replace('\\', '/')).append('`');
			}
			out.append(" — the `.txt` feeds `SchemeSweep --shape-file` / `ProjectFmmGaps --shapes-file`"
					+ " (RESYNC_PLAYBOOK step 4, run by a human/agent, never in CI)\n");
		}

		out.append("\n## Sedoglavic fmm_digest\n\n");
		out.append("- Strict gaps (Sedoglavic better): **").append(sed.gaps().size()).append("**\n");
		out.append("- Missing from us (Sedoglavic only): **").append(sed.missing().size()).append("**\n");
		out.append("- We-better: **").append(sed.weBetter().size()).append("**\n");
		out.append("- Upstream-unverified (same audit as FMM-Lille): **").append(sed.unverified().size()).append("**\n");
		appendGapTable(out, "Top 20 strict gaps", "Sedoglavic", sed.gaps());

		out.append("\n## FMM Université de Lille (local snapshot)\n\n");
		out.append("- Strict gaps: **").append(fmm.gaps().size()).append("**\n");
		out.append("- Missing from us: **").append(fmm.missing().size()).append("**\n");
		out.append("- Upstream-unverified (index rank not backed by FMM's artifact, "
				+ "`references/fmm-artifact-audit.json`): **").append(fmm.unverified().size()).append("**\n");
		appendGapTable(out, "Top 20 FMM strict gaps", "FMM", fmm.gaps());
		appendGapTable(out, "Upstream-unverified FMM rows (audited, not gaps)", "FMM", fmm.unverified());

		appendPerminovLayout(out, "Perminov status.json", "Perminov", perm);
		appendPerminovLayout(out, "mizopl we*star status.json", "mizopl", mizopl);

		return out.toString();
	}

	/** The three-axis section of a Perminov-layout source. Every "Strict gaps" line is
	 *  counted by the workflow's issue step, so a Z / ZT record we lack alerts like a Q one. */
	private static void appendPerminovLayout(StringBuilder out, String title, String col, SourceScan s) {
		out.append("\n## ").append(title).append("\n\n");
		out.append("- Strict gaps: **").append(s.gaps().size()).append("**\n");
		out.append("- Missing from us: **").append(s.missing().size()).append("**\n");
		out.append("- Strict gaps (Z axis, integer schemes): **").append(s.zGaps().size()).append("**\n");
		out.append("- Strict gaps (ZT axis, ternary-integer explicit schemes): **").append(s.ztGaps().size())
				.append("**\n");
		appendGapTable(out, "Top 20 " + col + " strict gaps", col, s.gaps());
		appendGapTable(out, "Top 20 " + col + " Z-axis gaps (ours = best integer scheme, −1 = none)", col + " Z",
				s.zGaps());
		appendGapTable(out, "Top 20 " + col + " ZT-axis gaps (ours = best explicit ternary scheme, −1 = none)",
				col + " ZT", s.ztGaps());
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
					.append(" | ").append(g.gap() == Integer.MAX_VALUE ? "∞" : String.valueOf(g.gap())).append(" |\n");
		});
	}
}
