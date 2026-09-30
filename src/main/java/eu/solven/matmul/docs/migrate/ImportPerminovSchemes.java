package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Duration;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * Import scheme files (factor matrices) from
 * <a href="https://github.com/dronperminov/FastMatrixMultiplication">dronperminov/FastMatrixMultiplication</a>
 * — the Java port of {@code tools/import_dronperminov.py}.
 *
 * <p>Perminov publishes the actual U/V/W per format under {@code schemes/results/<FIELD>/}
 * (2.5k files; {@code <n>x<m>x<p>_m<rank>_<tag>.json}). We rely on those <b>direct
 * scheme JSONs</b> (the {@code status.json} {@code path} pointers are stale / 404),
 * and {@link SchemeIO#read} parses Perminov's native shape ({@code n}, {@code m}=rank,
 * {@code u}/{@code v}/{@code w}) directly — so importing is: list the git tree →
 * download → {@link SchemeIO#read} → {@link Verifier#isExactNonCubic} → write the
 * canonical {@code known/section{maxdim}/<shape>-r<rank>-perminov_<tag>-<hash7>.json}.</p>
 *
 * <p>Idempotent: a {@code (shape, rank)} we already carry under {@code known/}
 * (any perminov file) is skipped <em>without downloading</em>. A NEW or BETTER
 * (lower-rank) upstream scheme is a fresh {@code (shape, rank)} → imported. Run the
 * digest sync ({@code SyncReferenceCatalogs --perminov}) first so the comparison
 * reflects the same upstream snapshot.</p>
 *
 * <p>Run (defaults: maxDim 32, minDim 2, no limit):</p>
 * <pre>MAVEN_OPTS="-Xmx2g" mvn -q -ntp exec:java \
 *   -Dexec.mainClass=eu.solven.matmul.docs.migrate.ImportPerminovSchemes \
 *   -Dexec.args="--max-dim=16"</pre>
 *
 * <p>After importing, stamp + regenerate: {@code StampFields}, {@code StampAdditions},
 * {@code GenerateCatalogManifest} (the latter computes {@code zt} from coefficients).</p>
 */
@Slf4j
public final class ImportPerminovSchemes {

	private ImportPerminovSchemes() {}

	private static final String TREE_URL =
			"https://api.github.com/repos/dronperminov/FastMatrixMultiplication/git/trees/master?recursive=1";
	private static final String RAW_BASE =
			"https://raw.githubusercontent.com/dronperminov/FastMatrixMultiplication/master/";
	/** Human-browsable GitHub blob base for {@code source_scheme_url} (the scheme's
	 *  file in Perminov's own repo — distinct from {@code source_paper_url}, the paper). */
	static final String BLOB_BASE =
			"https://github.com/dronperminov/FastMatrixMultiplication/blob/master/";
	private static final String USER_AGENT = "solven-matmul-catalog/perminov-import";
	private static final Path KNOWN = Path.of("src/main/resources/schemes/known");
	/** Reaction bases (bud-rich representatives; rank need not be the catalog best). */
	private static final Path BUD_BASES = Path.of("src/main/resources/schemes/bud-bases");
	/**
	 * Perminov's serendipitous-product BASES ({@code {shape}_m{rank}_{sha}_{FIELD}.json}):
	 * several CONTENT-DISTINCT schemes per {@code (shape, rank)} — e.g. six
	 * {@code 2x3x11_m55_*} — differing in BUD structure, each the base of published
	 * {@code s1 ⊗ˢ s2} ranks in the 17–32 band (perminov-serendipitous-catalog.json names
	 * the exact file per format). The generic {@code (shape, rank)} skip below kept at
	 * most ONE of them (whichever came first, not necessarily the bud-rich one a recipe
	 * needs), which is why 621 formats of that band sat below our catalog while the sync
	 * job reported "0 imported" (2026-09-30). These are keyed by upstream PATH instead
	 * and land in {@code bud-bases/} as reaction bases.
	 */
	private static final String SER_BASE_DIR = "schemes/results/serendipitous_base/";
	private static final String SER_PAPER_URL = "https://arxiv.org/abs/2606.02480";
	/** Upstream basename: {@code <n>x<m>x<p>_m<rank>_<tag>.json}. */
	private static final Pattern NAME =
			Pattern.compile("^(\\d+)x(\\d+)x(\\d+)_m(\\d+)_(.+)\\.json$");

	private static final JsonMapper MAPPER = JsonMapper.builder().build();

	public static void main(String[] args) throws Exception {
		int minDim = intArg(args, "--min-dim", 2);
		int maxDim = intArg(args, "--max-dim", 32);
		int limit = intArg(args, "--limit", 0);
		boolean overwrite = flag(args, "--overwrite");

		// Existing (shape, rank) we already carry as a Perminov import — skip those
		// without re-downloading. Keyed "NxMxP-rRANK".
		Set<String> have = existingPerminovKeys();
		log.info("already carry {} perminov (shape,rank) keys under known/", have.size());

		log.info("listing upstream scheme tree …");
		JsonNode tree = MAPPER.readTree(fetch(TREE_URL)).get("tree");
		List<String> paths = new ArrayList<>();
		for (JsonNode t : tree) {
			String p = t.path("path").asString();
			if (p.startsWith("schemes/results/") && p.endsWith(".json")) paths.add(p);
		}
		log.info("found {} upstream result schemes", paths.size());
		// Upstream paths of serendipitous bases we already hold (bud-bases/ + the known/
		// files at the same (shape, rank) keys) — the idempotence key for that directory.
		Set<String> havePaths = existingSerendipitousBasePaths(paths);
		log.info("already carry {} serendipitous-base upstream paths", havePaths.size());

		int wrote = 0, skipExisting = 0, skipRange = 0, fail = 0, processed = 0;
		long t0 = System.nanoTime();
		for (String path : paths) {
			String base = path.substring(path.lastIndexOf('/') + 1);
			Matcher m = NAME.matcher(base);
			if (!m.matches()) continue;
			int n = Integer.parseInt(m.group(1)), mm = Integer.parseInt(m.group(2)), p = Integer.parseInt(m.group(3));
			int rank = Integer.parseInt(m.group(4));
			String tag = m.group(5);
			int maxd = Math.max(n, Math.max(mm, p));
			if (maxd < minDim || maxd > maxDim) { skipRange++; continue; }
			String key = n + "x" + mm + "x" + p + "-r" + rank;
			boolean serBase = path.startsWith(SER_BASE_DIR);
			// Rank results: one Perminov file per (shape, rank) is enough. Serendipitous
			// bases: every content-distinct file matters (bud profile) → key by path.
			if (!overwrite && (serBase ? havePaths.contains(path) : have.contains(key))) {
				skipExisting++;
				continue;
			}

			try {
				String body = fetch(RAW_BASE + path);
				// "_reduced" schemes use Perminov's sparse-list JSON shape (index/value
				// lists), read by SchemeIO.readReduced; the rest are dense u/v/w.
				NonCubicBilinearAlgorithm alg = base.contains("reduced")
						? SchemeIO.readReduced(MAPPER.readTree(body))
						: SchemeIO.read(body);
				if (alg.n != n || alg.m != mm || alg.p != p || alg.r != rank) {
					log.warn("[SKIP] {} content {}x{}x{} r{} != filename", base, alg.n, alg.m, alg.p, alg.r);
					fail++;
					continue;
				}
				if (!Verifier.isExactNonCubic(alg)) {
					log.warn("[FAIL] {} did not verify as exact matmul", base);
					fail++;
					continue;
				}
				String hash = SchemeIO.contentHash(alg);
				String hash7 = hash.substring(0, 7);
				if (serBase && !overwrite && existsWithHash(n, mm, p, rank, hash7)) {
					// Identical content already in the catalog (typically the ONE variant the
					// (shape, rank) rule imported into known/ earlier) — nothing to add.
					havePaths.add(path);
					skipExisting++;
					continue;
				}
				// Clean cosmetic label: perminov_{ZT|Z|Q} (the raw tag carries cr/cn/
				// hash cruft); content + metadata are authoritative, the name is a label.
				String note = (serBase ? "perminov_serbase_" : "perminov_") + fieldOf(path, tag);
				Path dir = (serBase ? BUD_BASES : KNOWN).resolve("section" + maxd);
				Files.createDirectories(dir);
				File out = dir.resolve(n + "x" + mm + "x" + p + "-r" + rank + "-" + note + "-" + hash7 + ".json").toFile();
				if (out.exists() && !overwrite) { skipExisting++; continue; }
				SchemeIO.write(alg, out);
				// Metadata the importer is authoritative for; fields stamped from the
				// upstream field-class (integer ⇒ all fields; rational ⇒ Q/R/C). zt and
				// additions are (re)computed downstream by GenerateCatalogManifest.
				Map<String, Object> meta = new LinkedHashMap<>();
				// Attribute by Perminov's OWN directory layout, not blindly to "Perminov
				// 2023": his schemes/known/<sub> subtree re-hosts others' work (e.g.
				// known/meta_flip_graph = Kauers & Wood 2025). This loop only lists
				// schemes/results/* (line ~94 — Perminov's own), so forPath returns
				// PERMINOV_OWN here; the routing is defensive so a widened filter stays
				// correct, and SKIP_FRESH_IMPORT mirrors (tensor=FMM, matmulcatalog=ours)
				// are never freshly pulled.
				var attr = eu.solven.matmul.catalog.PerminovKnownAttribution.forPath(path)
						.orElse(new eu.solven.matmul.catalog.PerminovKnownAttribution.Attribution(
								"Perminov 2023",
								eu.solven.matmul.catalog.PerminovKnownAttribution.Disposition.PERMINOV_OWN));
				if (attr.disposition()
						== eu.solven.matmul.catalog.PerminovKnownAttribution.Disposition.SKIP_FRESH_IMPORT) {
					log.info("[skip-mirror] {} — obtained by other means ({})", base, attr.source());
					skipExisting++;
					continue;
				}
				if (serBase) {
					// A REACTION base, not a rank claim: the bud-rich representative behind
					// Perminov's published serendipitous products. Stamp the content hash —
					// durableLeafRef refuses to pin an un-hashed base when a reaction wins
					// (the ⟨20,28,28⟩=8434 persist failure, 2026-07-08).
					meta.put("source", "Perminov 2026 (serendipitous)");
					meta.put("year", 2026);
					meta.put("source_paper_url", SER_PAPER_URL);
					meta.put("hash", hash);
					meta.put("discovery", false);
					meta.put("attribution_for_rank", "reaction base only (bud-rich representative of a "
							+ "serendipitous product, Perminov arXiv:2606.02480) — the rank at this shape is "
							+ "attributed by the catalog's rank-best scheme, not by this file");
				} else {
					meta.put("source", attr.source());
					if (!attr.isPerminovOwn()) {
						meta.put("imported_via", "Perminov FastMatrixMultiplication");
					}
				}
				meta.put("original_source_path", path);
				// Clear pointer to the scheme's file in Perminov's own repo (the file,
				// not the author/paper) — carried into catalog.json and the SPA.
				meta.put("source_scheme_url", BLOB_BASE + path);
				meta.put("commutative", false);
				meta.put("verified", true);
				meta.put("fields", fieldsForField(fieldOf(path, tag)));
				SchemeIO.addFields(out, meta, /* apply */ true);
				have.add(key);
				if (serBase) havePaths.add(path);
				wrote++;
				if (wrote % 50 == 0) {
					long ms = (System.nanoTime() - t0) / 1_000_000L;
					log.info("[progress] {} imported ({} skipped-existing), {}ms elapsed", wrote, skipExisting, ms);
				}
			} catch (RuntimeException | IOException e) {
				log.warn("[ERR] {}: {}", base, e.toString());
				fail++;
			}
			if (limit > 0 && wrote >= limit) break;
			processed++;
		}
		long ms = (System.nanoTime() - t0) / 1_000_000L;
		log.info("Done: {} imported, {} skipped-existing, {} out-of-range, {} failed ({}ms). "
				+ "Next: StampFields / StampAdditions / GenerateCatalogManifest.",
				wrote, skipExisting, skipRange, fail, ms);
	}

	/** Upstream field-class for an import, from the {@code schemes/results/<FIELD>/} path
	 *  segment (falls back to the filename tag). {@code Z}/{@code ZT} ⇒ integer. */
	private static String fieldOf(String path, String tag) {
		String[] seg = path.split("/");
		if (seg.length >= 3) {
			String f = seg[2];
			if (f.equals("Z") || f.equals("ZT") || f.equals("Q")) return f;
		}
		String t = tag.toUpperCase(java.util.Locale.ROOT);
		if (t.contains("ZT")) return "ZT";
		if (t.contains("Q")) return "Q";
		return "Z"; // integer default (serendipitous_base etc. are ZT/Z)
	}

	/** Integer schemes (Z/ZT) reduce over every field; rational (Q) only Q/R/C.
	 *  Conservative — never claims F2/F3 from a rational (that needs verification). */
	private static List<String> fieldsForField(String field) {
		return field.equals("Q")
				? List.of("Q", "R", "C")
				: List.of("F2", "F3", "Z", "Q", "R", "C");
	}

	/** Scan {@code known/} for existing perminov imports → set of "NxMxP-rRANK". */
	private static Set<String> existingPerminovKeys() throws IOException {
		Set<String> keys = new HashSet<>();
		if (!Files.isDirectory(KNOWN)) return keys;
		Pattern fn = Pattern.compile("^(\\d+x\\d+x\\d+-r\\d+)-perminov_.*\\.json$");
		try (Stream<Path> w = Files.walk(KNOWN)) {
			w.filter(Files::isRegularFile).forEach(p -> {
				Matcher m = fn.matcher(p.getFileName().toString());
				if (m.matches()) keys.add(m.group(1));
			});
		}
		return keys;
	}

	/**
	 * Upstream {@code serendipitous_base/} paths already imported: the
	 * {@code original_source_path} of every file under {@code bud-bases/}, plus those of
	 * the {@code known/} Perminov files at a {@code (shape, rank)} that occurs among the
	 * upstream serendipitous bases (the single variant the pre-2026-09-30 rule let in).
	 * Bounded read: only files at those keys are opened.
	 */
	private static Set<String> existingSerendipitousBasePaths(List<String> upstreamPaths) throws IOException {
		Set<String> serKeys = new HashSet<>();
		for (String path : upstreamPaths) {
			if (!path.startsWith(SER_BASE_DIR)) continue;
			Matcher m = NAME.matcher(path.substring(path.lastIndexOf('/') + 1));
			if (m.matches()) serKeys.add(m.group(1) + "x" + m.group(2) + "x" + m.group(3) + "-r" + m.group(4));
		}
		Set<String> out = new HashSet<>();
		Pattern fn = Pattern.compile("^(\\d+x\\d+x\\d+-r\\d+)-.*\\.json$");
		Pattern src = Pattern.compile("\"original_source_path\"\\s*:\\s*\"([^\"]+)\"");
		for (Path root : new Path[] { BUD_BASES, KNOWN }) {
			if (!Files.isDirectory(root)) continue;
			List<Path> files;
			try (Stream<Path> w = Files.walk(root)) {
				files = w.filter(Files::isRegularFile).toList();
			}
			for (Path f : files) {
				Matcher m = fn.matcher(f.getFileName().toString());
				if (!m.matches() || !serKeys.contains(m.group(1))) continue;
				Matcher s = src.matcher(Files.readString(f));
				if (s.find() && s.group(1).startsWith(SER_BASE_DIR)) out.add(s.group(1));
			}
		}
		return out;
	}

	/** Is a scheme with this exact content (hash7 in the canonical filename) already on
	 *  disk at {@code (shape, rank)}, under {@code known/} or {@code bud-bases/}? */
	private static boolean existsWithHash(int n, int m, int p, int rank, String hash7) throws IOException {
		int maxd = Math.max(n, Math.max(m, p));
		String prefix = n + "x" + m + "x" + p + "-r" + rank + "-";
		for (Path root : new Path[] { BUD_BASES, KNOWN }) {
			Path dir = root.resolve("section" + maxd);
			if (!Files.isDirectory(dir)) continue;
			try (Stream<Path> ls = Files.list(dir)) {
				if (ls.anyMatch(f -> {
					String name = f.getFileName().toString();
					return name.startsWith(prefix) && name.endsWith("-" + hash7 + ".json");
				})) {
					return true;
				}
			}
		}
		return false;
	}

	private static int intArg(String[] args, String key, int dflt) {
		for (String a : args) {
			if (a.startsWith(key + "=")) return Integer.parseInt(a.substring(key.length() + 1));
		}
		return dflt;
	}

	private static boolean flag(String[] args, String key) {
		for (String a : args) if (a.equals(key) || a.equals(key + "=true")) return true;
		return false;
	}

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
				if (resp.statusCode() / 100 == 2) return resp.body();
				last = new IOException("HTTP " + resp.statusCode() + " for " + url);
			} catch (IOException e) {
				last = e;
			}
			if (attempt < 2) Thread.sleep(2000L * (attempt + 1));
		}
		throw last;
	}
}
