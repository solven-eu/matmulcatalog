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

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.docs.migrate.UpstreamSchemeRepos.UpstreamSchemeRepo;

import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * Pull new schemes from every repository in {@link UpstreamSchemeRepos} (or one of
 * them with {@code --repo=key}).
 *
 * <p>Per upstream scheme file: skipped without downloading when the catalog already
 * holds that {@code (shape, rank, class)} ({@link KnownSchemeKeys}, content-driven, any
 * source); otherwise downloaded and put through the contributed-scheme gate
 * ({@link ImportContributedSchemes#gate}: exact symbolic proof, random spot-check,
 * fields from coefficients incl. the mod-2 / mod-3 reductions) and writer (canonical
 * re-encoding, Phase-2 metrics, disk round-trip). A file that fails the gate is
 * logged and skipped — it never lands.</p>
 *
 * <p>{@code discovery} is stamped {@code "TBD"}: an automated import can see that a
 * rank beats <em>this catalog</em>, not that it is a record — the dataset's own table
 * (or FMM / Perminov) has to be read for that; {@link RestampContributedProvenance}
 * is the audit step. The verdict against the catalog at import time is recorded in
 * {@code discovery_note}.</p>
 *
 * <pre>mvn -q -ntp exec:java \
 *   -Dexec.mainClass=eu.solven.matmul.docs.migrate.ImportUpstreamSchemeRepos \
 *   [-Dexec.args="--repo=witteveen --dry-run --max-dim=32 --limit=0"]</pre>
 *
 * <p>Runs in {@code sync-reference-catalogs.yml} (every 3 h) after the Perminov
 * import; the job's manifest regeneration, field-claim gate and commit cover what it
 * writes.</p>
 */
@Slf4j
public final class ImportUpstreamSchemeRepos {
	private ImportUpstreamSchemeRepos() {}

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");
	private static final Path KNOWN = SCHEMES_ROOT.resolve("known");
	private static final JsonMapper MAPPER = JsonMapper.builder().build();

	/** An upstream scheme file the catalog does not hold yet. */
	public record Candidate(String path, int n, int m, int p, int rank, String cls) {
	}

	public static void main(String[] args) throws Exception {
		int minDim = intArg(args, "--min-dim", 2);
		int maxDim = intArg(args, "--max-dim", 32);
		int limit = intArg(args, "--limit", 0);
		boolean dryRun = List.of(args).contains("--dry-run");
		String only = null;
		for (String a : args) if (a.startsWith("--repo=")) only = a.substring("--repo=".length());

		Set<String> have = KnownSchemeKeys.scan(KNOWN);
		log.info("already carry {} (shape, rank, class) keys under known/", have.size());
		FieldAwareLookup lookup = null;
		int totalWrote = 0, totalFail = 0;
		for (UpstreamSchemeRepo repo : only == null ? UpstreamSchemeRepos.ALL : List.of(UpstreamSchemeRepos.byKey(only))) {
			log.info("=== {} ({}) ===", repo.key(), repo.repo());
			List<String> paths = new ArrayList<>();
			for (JsonNode t : MAPPER.readTree(ImportPerminovSchemes.fetch(repo.treeUrl())).get("tree")) {
				paths.add(t.path("path").asString());
			}
			List<Candidate> todo = plan(repo, paths, have, minDim, maxDim);
			long upstream = paths.stream().filter(p -> repo.parse(p).isPresent()).count();
			log.info("{}: {} upstream schemes, {} not held yet", repo.key(), upstream, todo.size());
			if (upstream == 0) {
				// The layout moved (renamed directory / branch / naming): say so, do not read it as "nothing new".
				log.warn("[layout] no file under {} matches {} — the upstream layout changed; this importer is blind "
						+ "until the registry pattern is updated", repo.repo(), repo.name());
				System.out.println("::warning title=" + repo.repo() + " layout changed::ImportUpstreamSchemeRepos matched 0 files");
			}
			if (todo.isEmpty()) continue;
			if (dryRun) {
				todo.forEach(c -> log.info("[dry-run] {} would import {} ({})", repo.key(), c.path(), c.cls()));
				continue;
			}
			if (lookup == null) lookup = new FieldAwareLookup(Field.fromTag("Q"), SCHEMES_ROOT);
			int[] counts = importAll(repo, todo, have, lookup, limit);
			totalWrote += counts[0];
			totalFail += counts[1];
		}
		log.info("ImportUpstreamSchemeRepos: {} imported, {} failed. Next: GenerateCatalogManifest.", totalWrote, totalFail);
	}

	/** Download, gate and write each candidate; returns {@code {written, failed}}. */
	static int[] importAll(UpstreamSchemeRepo repo, List<Candidate> todo, Set<String> have, FieldAwareLookup lookup,
			int limit) throws IOException, InterruptedException {
		int wrote = 0, fail = 0, already = 0;
		long t0 = System.nanoTime();
		for (Candidate c : todo) {
			Path tmp = Files.createTempFile(repo.key() + "-", ".json");
			try {
				Files.writeString(tmp, ImportPerminovSchemes.fetch(repo.rawBase() + c.path()), StandardCharsets.UTF_8);
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
				if (!actual.equals(c.cls()) && !KnownSchemeKeys.Q.equals(c.cls())) {
					log.warn("{}: filename says {} but the coefficients are {} — stamping from the coefficients", c.path(),
							c.cls(), actual);
				}
				int maxd = Math.max(alg.n, Math.max(alg.m, alg.p));
				File out = KNOWN.resolve("section" + maxd).resolve(SchemeIO.canonicalName(alg, repo.note())).toFile();
				if (out.isFile()) {
					// Same content already imported (its class is finer than the filename said).
					already++;
					continue;
				}
				out = ImportContributedSchemes.write(g, options(repo, g));
				Map<String, Object> provenance = new LinkedHashMap<>();
				provenance.put("original_source_path", c.path());
				provenance.put("source_scheme_url", repo.blobBase() + c.path());
				SchemeIO.addFields(out, provenance, /* apply */ true);
				KnownSchemeKeys.add(have, alg.n, alg.m, alg.p, alg.r, actual);
				wrote++;
				log.info("[progress] {} {}/{} imported — {} r={} {} catalog(Q): {} ({}ms elapsed)", repo.key(), wrote,
						todo.size(), out.getName(), alg.r, actual, g.vsCatalog(), (System.nanoTime() - t0) / 1_000_000L);
			} catch (RuntimeException | IOException e) {
				log.warn("[ERR] {}: {}", c.path(), e.toString());
				fail++;
			} finally {
				Files.deleteIfExists(tmp);
			}
			if (limit > 0 && wrote >= limit) break;
		}
		log.info("{}: {} imported, {} already held by content, {} failed ({}ms)", repo.key(), wrote, already, fail,
				(System.nanoTime() - t0) / 1_000_000L);
		return new int[] { wrote, fail };
	}

	/**
	 * Which upstream files to pull: those the repo's pattern names, within the
	 * dimension window, whose {@code (shape, rank, class)} is not already held. Pure —
	 * the tests drive it with a literal listing.
	 */
	public static List<Candidate> plan(UpstreamSchemeRepo repo, List<String> upstreamPaths, Set<String> have, int minDim,
			int maxDim) {
		List<Candidate> out = new ArrayList<>();
		for (String path : upstreamPaths) {
			Candidate c = repo.parse(path).orElse(null);
			if (c == null) continue;
			int maxd = Math.max(c.n(), Math.max(c.m(), c.p()));
			if (maxd < minDim || maxd > maxDim) continue;
			if (KnownSchemeKeys.covers(have, c.n(), c.m(), c.p(), c.rank(), c.cls())) continue;
			out.add(c);
		}
		return out;
	}

	/** Provenance stamped on an automated import (see the class doc for {@code discovery: "TBD"}). */
	private static ImportContributedSchemes.Options options(UpstreamSchemeRepo repo, ImportContributedSchemes.Gate g) {
		String today = LocalDate.now(ZoneOffset.UTC).toString();
		String discoveryNote = repo.citation() + ". Auto-imported " + today + " by ImportUpstreamSchemeRepos (exact"
				+ " BigInteger symbolic proof, random spot-check, mod-2 / mod-3 reductions, canonical-disk round-trip)."
				+ " Against this catalog over Q at import time: " + g.vsCatalog() + ". Record status (any ring / Z / ZT"
				+ " class) not yet audited against the dataset's own table — `discovery` stays TBD until it is"
				+ " (RestampContributedProvenance).";
		String notes = "Re-encoded from the upstream Perminov-style JSON (n, m, u, v, w) into the catalog's canonical"
				+ " format — no coefficient was changed. " + repo.citation() + ".";
		return new ImportContributedSchemes.Options(List.of(), repo.note(), repo.source(), repo.year(),
				"https://github.com/" + repo.repo(), "TBD", null, discoveryNote, notes, "known", "Q", true);
	}

	private static int intArg(String[] args, String key, int dflt) {
		for (String a : args) {
			if (a.startsWith(key + "=")) return Integer.parseInt(a.substring(key.length() + 1));
		}
		return dflt;
	}
}
