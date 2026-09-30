package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.LocalDate;
import java.time.ZoneOffset;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;

import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * Pull new schemes from Merlijn S. Witteveen's dataset
 * <a href="https://github.com/MerlijnW70/fmm-schemes">MerlijnW70/fmm-schemes</a>
 * (doi:10.5281/zenodo.22943995, <b>CC BY 4.0</b>) — REFERENCES.md [90].
 *
 * <p>Why a dedicated channel: the dataset is cited by FMM-Lille (⟨11,13,15⟩,
 * ⟨11,14,14⟩) and mirrored by Perminov under {@code schemes/known/MerlijnW70_fmm_schemes/},
 * but neither route reached the catalog — we import no FMM artifact automatically,
 * and {@link ImportPerminovSchemes} only lists Perminov's own {@code schemes/results/}.
 * Release 1.0 (2026-09-24, with a new best over every ring at ⟨11,13,15⟩ and
 * ⟨11,14,14⟩) therefore sat un-imported for six days behind green sync jobs. The
 * origin repo is also ahead of its mirrors (release 1.2's ⟨7,11,15⟩=772 was nowhere
 * else on 2026-09-30).</p>
 *
 * <p>Per upstream {@code schemes/{n}x{m}x{p}_m{rank}_{tag}.json}: skipped without
 * downloading when the catalog already holds that {@code (shape, rank, class)}
 * ({@link KnownSchemeKeys}, content-driven, any source); otherwise downloaded and put
 * through the contributed-scheme gate ({@link ImportContributedSchemes#gate}: exact
 * symbolic proof, random spot-check, fields from coefficients incl. the mod-2 / mod-3
 * reductions) and writer (canonical re-encoding, Phase-2 metrics, disk round-trip).
 * A file that fails the gate is logged and skipped — it never lands.</p>
 *
 * <p>{@code discovery} is stamped {@code "TBD"}: an automated import can see that a
 * rank beats <em>this catalog</em>, not that it is a record — the dataset's own table
 * (or FMM / Perminov) has to be read for that. The verdict against the catalog at
 * import time is recorded in {@code discovery_note}.</p>
 *
 * <pre>mvn -q -ntp exec:java \
 *   -Dexec.mainClass=eu.solven.matmul.docs.migrate.ImportWitteveenSchemes \
 *   [-Dexec.args="--dry-run --max-dim=32 --limit=0"]</pre>
 *
 * <p>Runs in {@code sync-reference-catalogs.yml} (every 3 h) after the Perminov
 * import; the job's manifest regeneration, field-claim gate and commit cover what it
 * writes.</p>
 */
@Slf4j
public final class ImportWitteveenSchemes {
	private ImportWitteveenSchemes() {}

	private static final String REPO = "MerlijnW70/fmm-schemes";
	private static final String TREE_URL = "https://api.github.com/repos/" + REPO + "/git/trees/main?recursive=1";
	private static final String RAW_BASE = "https://raw.githubusercontent.com/" + REPO + "/main/";
	private static final String BLOB_BASE = "https://github.com/" + REPO + "/blob/main/";
	private static final String CITATION = "Merlijn S. Witteveen, New fast matrix multiplication schemes, "
			+ "doi:10.5281/zenodo.22943995, github.com/" + REPO + " — licensed CC BY 4.0";

	static final String NOTE = "witteveen_2026";
	static final String SOURCE = "Witteveen 2026";

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");
	private static final Path KNOWN = SCHEMES_ROOT.resolve("known");
	/** Upstream layout: {@code schemes/<n>x<m>x<p>_m<rank>_<tag>.json}. */
	private static final Pattern NAME = Pattern.compile("^schemes/(\\d+)x(\\d+)x(\\d+)_m(\\d+)_([A-Za-z0-9]+)\\.json$");

	private static final JsonMapper MAPPER = JsonMapper.builder().build();

	/** An upstream scheme file the catalog does not hold yet. */
	record Candidate(String path, int n, int m, int p, int rank, String cls) {
	}

	public static void main(String[] args) throws Exception {
		int minDim = intArg(args, "--min-dim", 2);
		int maxDim = intArg(args, "--max-dim", 32);
		int limit = intArg(args, "--limit", 0);
		boolean dryRun = List.of(args).contains("--dry-run");

		Set<String> have = KnownSchemeKeys.scan(KNOWN);
		log.info("already carry {} (shape, rank, class) keys under known/", have.size());

		log.info("listing {} …", REPO);
		List<String> paths = new ArrayList<>();
		for (JsonNode t : MAPPER.readTree(ImportPerminovSchemes.fetch(TREE_URL)).get("tree")) {
			paths.add(t.path("path").asString());
		}
		List<Candidate> todo = plan(paths, have, minDim, maxDim);
		long upstream = paths.stream().filter(p -> NAME.matcher(p).matches()).count();
		log.info("{} upstream schemes, {} not held yet", upstream, todo.size());
		if (upstream == 0) {
			// The layout moved (renamed directory / branch / naming): say so, do not read it as "nothing new".
			log.warn("[layout] no file under {} matches schemes/<n>x<m>x<p>_m<rank>_<tag>.json — "
					+ "the upstream layout changed; this importer is blind until NAME is updated", REPO);
			System.out.println("::warning title=fmm-schemes layout changed::ImportWitteveenSchemes matched 0 upstream files");
		}
		if (todo.isEmpty()) {
			return;
		}
		if (dryRun) {
			todo.forEach(c -> log.info("[dry-run] would import {} ({})", c.path(), c.cls()));
			return;
		}

		FieldAwareLookup lookup = new FieldAwareLookup(Field.fromTag("Q"), SCHEMES_ROOT);
		int wrote = 0, fail = 0, already = 0;
		long t0 = System.nanoTime();
		for (Candidate c : todo) {
			Path tmp = Files.createTempFile("witteveen-", ".json");
			try {
				Files.writeString(tmp, ImportPerminovSchemes.fetch(RAW_BASE + c.path()), StandardCharsets.UTF_8);
				ImportContributedSchemes.Gate g = ImportContributedSchemes.gate(tmp.toFile(), lookup);
				NonCubicBilinearAlgorithm alg = g.alg();
				if (!KnownSchemeKeys.key(alg.n, alg.m, alg.p, alg.r, "").equals(KnownSchemeKeys.key(c.n(), c.m(), c.p(), c.rank(), ""))) {
					log.warn("[SKIP] {} content ⟨{},{},{}⟩ r{} != filename", c.path(), alg.n, alg.m, alg.p, alg.r);
					fail++;
					continue;
				}
				if (!g.passed()) {
					log.warn("[FAIL] {} did not pass the gate (exact={}, spot-check={})", c.path(), g.exactSymbolic(),
							g.spotCheck());
					fail++;
					continue;
				}
				String actual = g.ternary() ? KnownSchemeKeys.ZT : g.integer() ? KnownSchemeKeys.Z : KnownSchemeKeys.Q;
				if (!actual.equals(c.cls())) {
					log.warn("{}: filename says {} but the coefficients are {} — stamping from the coefficients", c.path(),
							c.cls(), actual);
				}
				int maxd = Math.max(alg.n, Math.max(alg.m, alg.p));
				File out = KNOWN.resolve("section" + maxd).resolve(SchemeIO.canonicalName(alg, NOTE)).toFile();
				if (out.isFile()) {
					// Same content already imported (its class differs from the filename tag).
					already++;
					continue;
				}
				out = ImportContributedSchemes.write(g, options(g));
				Map<String, Object> provenance = new LinkedHashMap<>();
				provenance.put("original_source_path", c.path());
				provenance.put("source_scheme_url", BLOB_BASE + c.path());
				SchemeIO.addFields(out, provenance, /* apply */ true);
				KnownSchemeKeys.add(have, alg.n, alg.m, alg.p, alg.r, actual);
				wrote++;
				log.info("[progress] {}/{} imported — {} r={} {} catalog(Q): {} ({}ms elapsed)", wrote, todo.size(),
						out.getName(), alg.r, actual, g.vsCatalog(), (System.nanoTime() - t0) / 1_000_000L);
			} catch (RuntimeException | IOException e) {
				log.warn("[ERR] {}: {}", c.path(), e.toString());
				fail++;
			} finally {
				Files.deleteIfExists(tmp);
			}
			if (limit > 0 && wrote >= limit) break;
		}
		log.info("Done: {} imported, {} already held by content, {} failed ({}ms). Next: GenerateCatalogManifest.", wrote,
				already, fail, (System.nanoTime() - t0) / 1_000_000L);
	}

	/**
	 * Which upstream files to pull: those named like a scheme, within the dimension
	 * window, whose {@code (shape, rank, class)} is not already held. Pure — the tests
	 * drive it with a literal listing.
	 */
	static List<Candidate> plan(List<String> upstreamPaths, Set<String> have, int minDim, int maxDim) {
		List<Candidate> out = new ArrayList<>();
		for (String path : upstreamPaths) {
			Matcher m = NAME.matcher(path);
			if (!m.matches()) continue;
			int n = Integer.parseInt(m.group(1)), mm = Integer.parseInt(m.group(2)), p = Integer.parseInt(m.group(3));
			int rank = Integer.parseInt(m.group(4));
			int maxd = Math.max(n, Math.max(mm, p));
			if (maxd < minDim || maxd > maxDim) continue;
			String cls = KnownSchemeKeys.classOfTag(m.group(5));
			if (KnownSchemeKeys.covers(have, n, mm, p, rank, cls)) continue;
			out.add(new Candidate(path, n, mm, p, rank, cls));
		}
		return out;
	}

	/** Provenance stamped on an automated import (see the class doc for {@code discovery: "TBD"}). */
	private static ImportContributedSchemes.Options options(ImportContributedSchemes.Gate g) {
		String today = LocalDate.now(ZoneOffset.UTC).toString();
		String discoveryNote = CITATION + ". Auto-imported " + today + " by ImportWitteveenSchemes (exact BigInteger"
				+ " symbolic proof, random spot-check, mod-2 / mod-3 reductions, canonical-disk round-trip). Against this"
				+ " catalog over Q at import time: " + g.vsCatalog() + ". Record status (any ring / Z / ZT class) not yet"
				+ " audited against the dataset's own table — `discovery` stays TBD until it is.";
		String notes = "Re-encoded from the upstream Perminov-style JSON (n, m, u, v, w) into the catalog's canonical"
				+ " format — no coefficient was changed. CC BY 4.0: credit Merlijn S. Witteveen, doi:10.5281/zenodo.22943995.";
		return new ImportContributedSchemes.Options(List.of(), NOTE, SOURCE, 2026, "https://github.com/" + REPO, "TBD",
				null, discoveryNote, notes, "known", "Q", true);
	}

	private static int intArg(String[] args, String key, int dflt) {
		for (String a : args) {
			if (a.startsWith(key + "=")) return Integer.parseInt(a.substring(key.length() + 1));
		}
		return dflt;
	}
}
