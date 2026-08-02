package eu.solven.matmul.docs.verify;

import java.io.IOException;
import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.SerializationFeature;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ArrayNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Refresh the external reference-rank digests we compare our catalog against —
 * Java port of {@code tools/scrape_perminov.py} + {@code tools/scrape_fmm_lille.py},
 * so the {@code sync-reference-catalogs} GitHub Action is one Maven command anyone
 * can replay locally.
 *
 * <p>Outputs (under {@code references/}):</p>
 * <ul>
 *   <li>{@code perminov-catalog.json} — per-format best NC char-0 rank = min over
 *       {Q, Z, ZT} from Perminov's {@code status.json}. Consumed by
 *       {@link GenerateCatalogManifest} (3-catalog comparison + {@code solven_discovery}).</li>
 *   <li>{@code fmm-lille-catalog.json} — every {@code (format, rank, refs)} row from
 *       the FMM Université de Lille main page. Consumed by {@link FmmCrossCheck},
 *       {@link FmmSearchComparison}, {@link GenerateFmmGapReport},
 *       {@link GenerateCatalogManifest} and {@link ScanUpstreamSources}.</li>
 *   <li>{@code fmm-lille-biblio.json} — bibliography entries (parity with the old
 *       scraper; not consumed by any code today).</li>
 * </ul>
 *
 * <p>Fault-tolerance: on a failed fetch the existing committed digest is left
 * untouched and the process exits non-zero (improving on the Python, which wrote
 * an empty digest on FMM downtime). The CI step guards the FMM refresh with
 * {@code || echo "keeping existing"} so a site outage doesn't fail the job.</p>
 *
 * <pre>
 *   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.verify.SyncReferenceCatalogs -Dexec.args="--perminov"
 *   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.verify.SyncReferenceCatalogs -Dexec.args="--fmm"
 *   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.verify.SyncReferenceCatalogs   # both
 * </pre>
 */
@Slf4j
public final class SyncReferenceCatalogs {

	private static final String PERMINOV_STATUS_URL =
			"https://raw.githubusercontent.com/dronperminov/FastMatrixMultiplication/master/schemes/status.json";
	private static final Path PERMINOV_OUT = Path.of("references/catalogs/perminov-catalog.json");
	private static final Path PERMINOV_SEREND_OUT = Path.of("references/catalogs/perminov-serendipitous-catalog.json");

	/** Human-browsable GitHub blob base — for per-entry {@code source_scheme_url}. */
	private static final String PERMINOV_BLOB_BASE =
			"https://github.com/dronperminov/FastMatrixMultiplication/blob/master/";
	/** The serendipitous 17–32 results are published in Perminov's June-2026 paper
	 *  "Meta Flip Graph meets Serendipitous Product" — NOT status.json. */
	private static final String PERMINOV_SEREND_SOURCE = "Perminov 2026 (serendipitous)";
	private static final String PERMINOV_SEREND_PAPER = "https://arxiv.org/abs/2606.02480";

	private static final String FMM_BASE = "https://fmm.univ-lille.fr";
	private static final Path FMM_CATALOG_OUT = Path.of("references/catalogs/fmm-lille-catalog.json");
	private static final Path FMM_BIBLIO_OUT = Path.of("references/fmm-lille-biblio.json");
	private static final String USER_AGENT =
			"strassen-catalog-scraper/0.1 (+https://github.com/solven-eu/strassen)";

	private static final JsonMapper MAPPER =
			JsonMapper.builder().enable(SerializationFeature.INDENT_OUTPUT).build();

	private static final Pattern PERMINOV_KEY = Pattern.compile("(\\d+)x(\\d+)x(\\d+)$");
	private static final Pattern CATALOG_ROW = Pattern.compile(
			"<a href=\"(\\d+x\\d+x\\d+\\.html)\">&lang;(\\d+)&times;(\\d+)&times;(\\d+):(\\d+)&rang;</a>"
					+ "(.*?)(?=<a href=\"\\d+x\\d+x\\d+\\.html\"|</table>|</section>)",
			Pattern.DOTALL);
	private static final Pattern NAIVE = Pattern.compile("#view\\d+x\\d+x\\d+\">(\\d+)</a>");
	// The index row's block-sum recipe, e.g. "4 ⟨5×5×8:144⟩ + 4 ⟨5×6×8:170⟩ + …".
	// This is the FRESH construction (the per-shape *.html pages lag the index), so it is
	// the authoritative source for diagnosing an FMM gap (which blocks, which allocation).
	private static final Pattern DESC_TD =
			Pattern.compile("<td[^>]*content=\"?description\"?[^>]*>(.*?)</td>", Pattern.DOTALL);
	private static final Pattern BIBLIO_KEY = Pattern.compile("biblio\\.html#([\\w:_-]+)");
	private static final Pattern BIBLIO_ENTRY = Pattern.compile(
			"class=\"bibtexnumber\">\\s*\\[<a name=\"([^\"]+)\">\\d+</a>\\].*?"
					+ "<td[^>]*class=\"bibtexitem\">(.*?)(?=<td[^>]*class=\"bibtexnumber\"|</table>|$)",
			Pattern.DOTALL);
	private static final Pattern HREF = Pattern.compile("href=\"([^\"]+)\"");
	private static final Pattern TAG = Pattern.compile("<[^>]+>");
	private static final Pattern WS = Pattern.compile("\\s+");
	private static final Pattern YEAR = Pattern.compile("(19[5-9]\\d|20[0-9]\\d|21\\d\\d)");

	private SyncReferenceCatalogs() {}

	public static void main(String[] args) throws Exception {
		boolean doPerminov = false;
		boolean doFmm = false;
		for (String a : args) {
			switch (a) {
				case "--perminov" -> doPerminov = true;
				case "--fmm" -> doFmm = true;
				case "--all" -> { doPerminov = true; doFmm = true; }
				default -> throw new IllegalArgumentException("unknown arg: " + a);
			}
		}
		if (!doPerminov && !doFmm) {
			doPerminov = true;
			doFmm = true;
		}
		if (doPerminov) {
			syncPerminov();
		}
		if (doFmm) {
			syncFmmLille();
		}
	}

	// ---- Perminov -----------------------------------------------------------

	private static void syncPerminov() throws IOException, InterruptedException {
		log.info("fetching {}", PERMINOV_STATUS_URL);
		JsonNode status = MAPPER.readTree(fetch(PERMINOV_STATUS_URL));

		List<ObjectNode> entries = new ArrayList<>();
		for (Map.Entry<String, JsonNode> e : status.properties()) {
			Matcher m = PERMINOV_KEY.matcher(e.getKey());
			if (!m.matches()) {
				continue;
			}
			JsonNode ranks = e.getValue().path("ranks");
			int n = Integer.parseInt(m.group(1)), mm = Integer.parseInt(m.group(2)), p = Integer.parseInt(m.group(3));
			// Emit PER FIELD-CLASS so the field-scoped comparison routes correctly
			// (ZT ⇒ Z ⇒ Q ⇒ R ⇒ C; an integer scheme reduces over F2/F3 too):
			//   - integer best = min(Z, ZT) — valid over EVERY field. ZT ⊂ Z so
			//     Z-best ≤ ZT-best and the Z entry dominates ZT everywhere; tag it
			//     ZT only when no plain-Z rank exists.
			//   - Q best — valid over Q/R/C ONLY; emit only when it STRICTLY beats
			//     the integer best (a rational scheme the integers can't match,
			//     e.g. ⟨4,4,4⟩ 48/Q vs 49/Z), else the integer entry already covers
			//     Q/R/C at that rank.
			// (The old code emitted just min over {Q,Z,ZT}, always tagged Q — which
			// starved the Z/F2/F3 tables of Perminov's integer schemes.)
			Integer zr = intRank(ranks.get("Z"));
			Integer ztr = intRank(ranks.get("ZT"));
			Integer qr = intRank(ranks.get("Q"));
			Integer intBest = null;
			String intField = null;
			if (zr != null) { intBest = zr; intField = "Z"; }
			if (ztr != null && (intBest == null || ztr < intBest)) { intBest = ztr; intField = "ZT"; }
			if (intBest != null) {
				entries.add(perminovEntry(n, mm, p, intBest, intField));
			}
			if (qr != null && (intBest == null || qr < intBest)) {
				entries.add(perminovEntry(n, mm, p, qr, "Q"));
			}
		}
		writePerminovDigest(PERMINOV_OUT, entries,
				"best NC rank per format, per field-class: an integer entry "
				+ "(field Z/ZT, valid over F2/F3/Z/Q/R/C) and — only when strictly better — "
				+ "a Q entry (valid over Q/R/C only). ZT ⇒ Z ⇒ Q ⇒ R ⇒ C. From status.json "
				+ "(maxdim ≤ 16); the 17–32 band is a SEPARATE catalog (perminov-serendipitous-catalog.json).");

		// The 17–32 band is kept as a SEPARATE additional catalog (not merged into
		// the status.json digest): distinct provenance (his serendipitous paper),
		// and the auto-synced digest stays a faithful mirror of status.json. The
		// consumers (loadExternalBest, GenerateSourceComparison) read both Perminov
		// files and fold them into one "perminov" external source.
		List<ObjectNode> serend = buildSerendipitousEntries();
		if (!serend.isEmpty()) {
			writePerminovDigest(PERMINOV_SEREND_OUT, serend,
					"Perminov's serendipitous 17–32 catalog (from "
					+ "references/perminov-serendipitous-17-32.json; his status.json omits this band). "
					+ "Same per-field-class convention; ranks are s1 ⊗ˢ s2 serendipitous products. "
					+ "Published in Perminov's June-2026 paper 'Meta Flip Graph meets "
					+ "Serendipitous Product' (arXiv:2606.02480).",
					PERMINOV_SEREND_PAPER);
		}
	}

	private static void writePerminovDigest(Path out, List<ObjectNode> entries, String note) throws IOException {
		writePerminovDigest(out, entries, note, null);
	}

	private static void writePerminovDigest(Path out, List<ObjectNode> entries, String note, String paperUrl)
			throws IOException {
		entries.sort(Comparator
				.<ObjectNode>comparingInt(SyncReferenceCatalogs::maxFormat)
				.thenComparing(SyncReferenceCatalogs::formatTuple));
		ObjectNode root = MAPPER.createObjectNode();
		root.put("source", "https://github.com/dronperminov/FastMatrixMultiplication");
		if (paperUrl != null) {
			root.put("source_paper_url", paperUrl);
		}
		root.put("note", note);
		ArrayNode arr = root.putArray("entries");
		entries.forEach(arr::add);
		writeJson(out, root);
		log.info("wrote {}: {} formats", out, entries.size());
	}

	/** Build digest entries for Perminov's serendipitous 17–32 catalog (same
	 *  per-field-class convention: ZT/Z ⇒ integer "Z" valid everywhere, Q ⇒
	 *  rational). Emitted to a SEPARATE file, not merged into the status.json digest. */
	private static List<ObjectNode> buildSerendipitousEntries() throws IOException {
		List<ObjectNode> out = new ArrayList<>();
		Path f = Path.of("references/perminov-serendipitous-17-32.json");
		if (!Files.isRegularFile(f)) {
			return out;
		}
		JsonNode root = MAPPER.readTree(Files.readString(f));
		for (Map.Entry<String, JsonNode> e : root.properties()) {
			String[] d = e.getKey().split("x");
			if (d.length != 3) {
				continue;
			}
			JsonNode v = e.getValue();
			JsonNode rk = v.get("serendipitous_rank");
			if (rk == null || !rk.isIntegralNumber()) {
				continue;
			}
			String path = v.has("path") ? v.get("path").asString() : "";
			String tag = path.isEmpty() ? "ZT"
					: path.substring(path.lastIndexOf('/') + 1).replace(".json", "").replaceAll(".*_", "");
			String field = tag.equals("Q") ? "Q" : "Z";
			out.add(serendipitousEntry(
					Integer.parseInt(d[0]), Integer.parseInt(d[1]), Integer.parseInt(d[2]), rk.asInt(), field, path));
		}
		return out;
	}

	/** A serendipitous-band entry: attributed to the June-2026 paper (NOT status.json),
	 *  with a {@code source_scheme_url} to the base scheme's file in Perminov's repo. */
	private static ObjectNode serendipitousEntry(int n, int m, int p, int rank, String field, String path) {
		ObjectNode node = MAPPER.createObjectNode();
		node.putArray("format").add(n).add(m).add(p);
		node.put("rank", rank);
		node.put("field", field);
		node.put("source", PERMINOV_SEREND_SOURCE);
		node.put("source_paper_url", PERMINOV_SEREND_PAPER);
		if (path != null && !path.isEmpty()) {
			node.put("source_scheme_url", PERMINOV_BLOB_BASE + path);
		}
		return node;
	}

	/** Parse a status.json rank cell into an int, or null when absent/non-integral. */
	private static Integer intRank(JsonNode r) {
		return (r != null && r.isIntegralNumber()) ? r.asInt() : null;
	}

	/** Build one {@code {format, rank, field, source}} digest entry. */
	private static ObjectNode perminovEntry(int n, int m, int p, int rank, String field) {
		ObjectNode node = MAPPER.createObjectNode();
		node.putArray("format").add(n).add(m).add(p);
		node.put("rank", rank);
		node.put("field", field);
		node.put("source", "perminov:status.json");
		return node;
	}

	private static int maxFormat(ObjectNode e) {
		ArrayNode f = (ArrayNode) e.get("format");
		return Math.max(f.get(0).asInt(), Math.max(f.get(1).asInt(), f.get(2).asInt()));
	}

	private static String formatTuple(ObjectNode e) {
		ArrayNode f = (ArrayNode) e.get("format");
		return String.format("%05d-%05d-%05d", f.get(0).asInt(), f.get(1).asInt(), f.get(2).asInt());
	}

	// ---- FMM Université de Lille -------------------------------------------

	private static void syncFmmLille() throws IOException, InterruptedException {
		// Catalog (consumed by 5 Java tools) — refresh first; bail before writing
		// if it fails, so the committed digest survives a site outage.
		log.info("fetching FMM catalog {}", FMM_BASE + "/");
		String catalogHtml = fetch(FMM_BASE + "/");
		List<ObjectNode> rows = new ArrayList<>();
		Matcher m = CATALOG_ROW.matcher(catalogHtml);
		while (m.find()) {
			int n = Integer.parseInt(m.group(2)), mm = Integer.parseInt(m.group(3)), p = Integer.parseInt(m.group(4));
			int rank = Integer.parseInt(m.group(5));
			String tail = m.group(6);
			ObjectNode node = MAPPER.createObjectNode();
			node.putArray("format").add(n).add(mm).add(p);
			node.put("rank", rank);
			Matcher nm = NAIVE.matcher(tail);
			if (nm.find()) {
				node.put("naive_rank", Integer.parseInt(nm.group(1)));
			} else {
				node.putNull("naive_rank");
			}
			ArrayNode refs = node.putArray("references");
			for (String key : dedup(BIBLIO_KEY, tail)) {
				refs.add(key);
			}
			// Block-sum recipe from the index row (authoritative; per-shape pages lag).
			Matcher dm = DESC_TD.matcher(tail);
			if (dm.find()) {
				String desc = TAG.matcher(dm.group(1)).replaceAll("");
				desc = desc.replace("&lang;", "⟨").replace("&rang;", "⟩")
						.replace("&times;", "×").replace("&InvisibleTimes;", "")
						.replace("&plus;", "+").replace("&minus;", "−").replace("&nbsp;", " ");
				desc = WS.matcher(desc).replaceAll(" ").trim();
				node.put("description", desc.isEmpty() ? null : desc);
			} else {
				node.putNull("description");
			}
			node.put("details_url", FMM_BASE + "/" + m.group(1));
			rows.add(node);
		}
		if (rows.isEmpty()) {
			throw new IOException("FMM catalog parse yielded 0 rows — refusing to overwrite "
					+ FMM_CATALOG_OUT + " (site changed or down)");
		}
		ObjectNode catalog = MAPPER.createObjectNode();
		catalog.put("source", FMM_BASE + "/");
		ArrayNode catArr = catalog.putArray("entries");
		rows.forEach(catArr::add);
		writeJson(FMM_CATALOG_OUT, catalog);
		log.info("wrote {}: {} catalog rows", FMM_CATALOG_OUT, rows.size());

		// Biblio (parity; consumed by nothing today) — best-effort, never fatal.
		try {
			log.info("fetching FMM biblio {}", FMM_BASE + "/biblio.html");
			String biblioHtml = fetch(FMM_BASE + "/biblio.html");
			List<ObjectNode> biblio = parseBiblio(biblioHtml);
			ObjectNode root = MAPPER.createObjectNode();
			root.put("source", FMM_BASE + "/biblio.html");
			ArrayNode arr = root.putArray("entries");
			biblio.forEach(arr::add);
			writeJson(FMM_BIBLIO_OUT, root);
			log.info("wrote {}: {} biblio entries", FMM_BIBLIO_OUT, biblio.size());
		} catch (IOException | RuntimeException ex) {
			log.warn("biblio refresh failed (keeping existing {}): {}", FMM_BIBLIO_OUT, ex.toString());
		}
	}

	private static List<ObjectNode> parseBiblio(String html) {
		List<ObjectNode> entries = new ArrayList<>();
		Matcher m = BIBLIO_ENTRY.matcher(html);
		while (m.find()) {
			String key = m.group(1);
			String content = m.group(2);
			ObjectNode node = MAPPER.createObjectNode();
			node.put("key", key);
			String raw = WS.matcher(TAG.matcher(content).replaceAll(" ")).replaceAll(" ").trim();
			node.put("raw", raw);
			Matcher ym = YEAR.matcher(raw);
			if (ym.find()) {
				node.put("year", Integer.parseInt(ym.group(1)));
			}
			int dot = raw.indexOf(". ");
			node.put("authors", (dot >= 0 ? raw.substring(0, dot) : raw).trim());
			ArrayNode links = node.putArray("links");
			Matcher hm = HREF.matcher(content);
			Set<String> seen = new LinkedHashSet<>();
			while (hm.find()) {
				if (seen.add(hm.group(1))) {
					links.addObject().put("url", hm.group(1));
				}
			}
			entries.add(node);
		}
		return entries;
	}

	private static List<String> dedup(Pattern p, String s) {
		Set<String> out = new LinkedHashSet<>();
		Matcher m = p.matcher(s);
		while (m.find()) {
			out.add(m.group(1));
		}
		return new ArrayList<>(out);
	}

	// ---- HTTP / IO ----------------------------------------------------------

	/** Fetch text with bounded retries + back-off (2s, 4s), matching the Python. */
	private static String fetch(String url) throws IOException, InterruptedException {
		HttpClient client = HttpClient.newBuilder()
				.followRedirects(HttpClient.Redirect.NORMAL)
				.connectTimeout(Duration.ofSeconds(30))
				.build();
		HttpRequest req = HttpRequest.newBuilder(URI.create(url))
				.header("User-Agent", USER_AGENT)
				.timeout(Duration.ofSeconds(60))
				.GET()
				.build();
		IOException last = null;
		for (int attempt = 0; attempt <= 2; attempt++) {
			try {
				HttpResponse<String> resp = client.send(req, HttpResponse.BodyHandlers.ofString());
				if (resp.statusCode() / 100 == 2) {
					return resp.body();
				}
				last = new IOException("HTTP " + resp.statusCode() + " for " + url);
			} catch (IOException e) {
				last = e;
			}
			if (attempt < 2) {
				Thread.sleep(2000L * (attempt + 1));
			}
		}
		throw new IOException("fetch failed after 3 attempts: " + url, last);
	}

	private static void writeJson(Path out, ObjectNode root) throws IOException {
		if (out.getParent() != null) {
			Files.createDirectories(out.getParent());
		}
		Files.writeString(out, MAPPER.writeValueAsString(root) + "\n");
	}
}
