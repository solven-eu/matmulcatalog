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
import eu.solven.matmul.additions.SlpJson;
import eu.solven.matmul.additions.SlpJson.SchemeSlp;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.docs.migrate.UpstreamSchemeRepos.UpstreamSchemeRepo;
import eu.solven.matmul.verifiers.Verifier;

import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ObjectNode;

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
 * <p><b>Addition-reduced programs</b> ({@link UpstreamSchemeRepo#slpGroup}): an upstream
 * file in Perminov's {@code _reduced} layout carries a straight-line program, so the
 * "already held" test keys on the <em>program's source</em> instead
 * ({@link KnownSchemeKeys#scanSlpSources}). The program is read with
 * {@link SlpJson#fromFreshVars}, replayed against the gated matrices (a program that
 * does not reconstruct them is refused), then attached to the scheme file we hold when
 * the matrices are identical (same content hash, any label) or to the freshly written
 * file otherwise — as {@code min_additions} + our object-layout {@code slp}, tier
 * {@code bound}, {@code slp.source} = the repository. An existing program at a count
 * at most the contributor's is kept ({@code [slp-kept]}).</p>
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

	/** An upstream scheme file the catalog does not hold yet ({@code reduced}: it carries an SLP). */
	public record Candidate(String path, int n, int m, int p, int rank, String cls, boolean reduced) {
		/** A plain rank-scheme candidate. */
		public Candidate(String path, int n, int m, int p, int rank, String cls) {
			this(path, n, m, p, rank, cls, false);
		}
	}

	public static void main(String[] args) throws Exception {
		int minDim = intArg(args, "--min-dim", 2);
		int maxDim = intArg(args, "--max-dim", 32);
		int limit = intArg(args, "--limit", 0);
		boolean dryRun = List.of(args).contains("--dry-run");
		String only = null;
		for (String a : args) if (a.startsWith("--repo=")) only = a.substring("--repo=".length());

		Set<String> have = KnownSchemeKeys.scan(KNOWN);
		Map<String, Set<String>> slpHave = KnownSchemeKeys.scanSlpSources(KNOWN);
		log.info("already carry {} (shape, rank, class) keys under known/, {} with a sourced SLP", have.size(),
				slpHave.size());
		FieldAwareLookup lookup = null;
		int totalWrote = 0, totalFail = 0;
		for (UpstreamSchemeRepo repo : only == null ? UpstreamSchemeRepos.ALL : List.of(UpstreamSchemeRepos.byKey(only))) {
			log.info("=== {} ({}) ===", repo.key(), repo.repo());
			List<String> paths = new ArrayList<>();
			for (JsonNode t : MAPPER.readTree(ImportPerminovSchemes.fetch(repo.treeUrl())).get("tree")) {
				paths.add(t.path("path").asString());
			}
			List<Candidate> todo = plan(repo, paths, have, slpHave, minDim, maxDim);
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
				todo.forEach(c -> log.info("[dry-run] {} would import {} ({}{})", repo.key(), c.path(), c.cls(),
						c.reduced() ? ", SLP" : ""));
				continue;
			}
			if (lookup == null) lookup = new FieldAwareLookup(Field.fromTag("Q"), SCHEMES_ROOT);
			int[] counts = importAll(repo, todo, have, slpHave, lookup, limit);
			totalWrote += counts[0];
			totalFail += counts[1];
		}
		log.info("ImportUpstreamSchemeRepos: {} imported, {} failed. Next: GenerateCatalogManifest.", totalWrote, totalFail);
	}

	/** Download, gate and write each candidate; returns {@code {written-or-attached, failed}}. */
	static int[] importAll(UpstreamSchemeRepo repo, List<Candidate> todo, Set<String> have,
			Map<String, Set<String>> slpHave, FieldAwareLookup lookup, int limit) throws IOException, InterruptedException {
		int wrote = 0, attached = 0, kept = 0, fail = 0, already = 0;
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
				SchemeSlp program = null;
				if (c.reduced()) {
					JsonNode root = SchemeIO.parseJson(tmp.toFile());
					if (!SlpJson.isFreshVarLayout(root)) {
						log.warn("[SKIP] {} is named as an addition-reduced file but carries no fresh variables", c.path());
						fail++;
						continue;
					}
					program = SlpJson.fromFreshVars(root);
					if (!program.reconstructs(alg)) {
						log.warn("[FAIL] {}: the program does not replay to its own matrices — refused", c.path());
						fail++;
						continue;
					}
					int declared = root.path("complexity").path("reduced").asInt(-1);
					if (declared >= 0 && declared != program.additions()) {
						log.warn("{}: declared complexity.reduced={} but the program costs {} additions — using the replay",
								c.path(), declared, program.additions());
					}
				}
				int maxd = Math.max(alg.n, Math.max(alg.m, alg.p));
				File out = KNOWN.resolve("section" + maxd).resolve(SchemeIO.canonicalName(alg, repo.note())).toFile();
				File holder = out.isFile() ? out : program == null ? null : heldWithSameContent(alg, maxd);
				if (holder != null) {
					if (program == null) {
						// Same content already imported (its class is finer than the filename said).
						already++;
						continue;
					}
					if (attachProgram(holder, alg, program, programMethod(repo), repo.source())) {
						attached++;
						slpHave.computeIfAbsent(KnownSchemeKeys.key(alg.n, alg.m, alg.p, alg.r, actual),
								k -> new java.util.HashSet<>()).add(repo.source());
						log.info("[progress] {} attached {} ({} additions, {}) to {} ({}ms elapsed)", repo.key(), c.path(),
								program.additions(), actual, holder.getName(), (System.nanoTime() - t0) / 1_000_000L);
					} else {
						kept++;
					}
					continue;
				}
				out = ImportContributedSchemes.write(g, options(repo, g));
				Map<String, Object> provenance = new LinkedHashMap<>();
				provenance.put("original_source_path", c.path());
				provenance.put("source_scheme_url", repo.blobBase() + c.path());
				SchemeIO.addFields(out, provenance, /* apply */ true);
				if (program != null) {
					attachProgram(out, alg, program, programMethod(repo), repo.source());
					slpHave.computeIfAbsent(KnownSchemeKeys.key(alg.n, alg.m, alg.p, alg.r, actual),
							k -> new java.util.HashSet<>()).add(repo.source());
				}
				KnownSchemeKeys.add(have, alg.n, alg.m, alg.p, alg.r, actual);
				wrote++;
				log.info("[progress] {} {}/{} imported — {} r={} {}{} catalog(Q): {} ({}ms elapsed)", repo.key(), wrote,
						todo.size(), out.getName(), alg.r, actual, program == null ? "" : " +" + program.additions() + " adds",
						g.vsCatalog(), (System.nanoTime() - t0) / 1_000_000L);
			} catch (RuntimeException | IOException e) {
				log.warn("[ERR] {}: {}", c.path(), e.toString());
				fail++;
			} finally {
				Files.deleteIfExists(tmp);
			}
			if (limit > 0 && wrote + attached >= limit) break;
		}
		log.info("{}: {} imported, {} programs attached, {} programs not better than what we hold, "
				+ "{} already held by content, {} failed ({}ms)", repo.key(), wrote, attached, kept, already, fail,
				(System.nanoTime() - t0) / 1_000_000L);
		return new int[] { wrote + attached, fail };
	}

	/** {@link #plan(UpstreamSchemeRepo, List, Set, Map, int, int)} without an SLP index (rank schemes only). */
	public static List<Candidate> plan(UpstreamSchemeRepo repo, List<String> upstreamPaths, Set<String> have, int minDim,
			int maxDim) {
		return plan(repo, upstreamPaths, have, Map.of(), minDim, maxDim);
	}

	/**
	 * Which upstream files to pull: those the repo's pattern names, within the
	 * dimension window, whose {@code (shape, rank, class)} is not already held — or, for
	 * an addition-reduced file, whose program from this repository is not already attached
	 * to a held scheme at that key. Pure — the tests drive it with a literal listing.
	 */
	public static List<Candidate> plan(UpstreamSchemeRepo repo, List<String> upstreamPaths, Set<String> have,
			Map<String, Set<String>> slpHave, int minDim, int maxDim) {
		List<Candidate> out = new ArrayList<>();
		for (String path : upstreamPaths) {
			Candidate c = repo.parse(path).orElse(null);
			if (c == null) continue;
			int maxd = Math.max(c.n(), Math.max(c.m(), c.p()));
			if (maxd < minDim || maxd > maxDim) continue;
			if (c.reduced()) {
				Set<String> sources = slpHave.get(KnownSchemeKeys.key(c.n(), c.m(), c.p(), c.rank(), c.cls()));
				if (sources != null && sources.contains(repo.source())) continue;
			} else if (KnownSchemeKeys.covers(have, c.n(), c.m(), c.p(), c.rank(), c.cls())) {
				continue;
			}
			out.add(c);
		}
		return out;
	}

	/**
	 * Attach a replayed program to a scheme file as {@code min_additions} + object-layout
	 * {@code slp} (tier {@code bound}), unless the file already holds a program or schedule
	 * at a count at most {@code program}'s. Returns whether the file was written.
	 */
	static boolean attachProgram(File f, NonCubicBilinearAlgorithm alg, SchemeSlp program, String method, String source)
			throws IOException {
		ObjectNode root = (ObjectNode) SchemeIO.parseJson(f);
		int current = root.path("min_additions").asInt(Integer.MAX_VALUE);
		JsonNode slp = root.get("slp");
		if (slp != null && slp.isObject() && slp.path("additions").isInt()) {
			current = Math.min(current, slp.get("additions").asInt());
		}
		if (root.path("scheduled_additions").isInt()) {
			current = Math.min(current, root.get("scheduled_additions").asInt());
		}
		if (program.additions() >= current) {
			log.info("[slp-kept] {} already holds a program at {} ≤ {} additions ({})", f.getName(), current,
					program.additions(), source);
			return false;
		}
		if (!root.path("additions").isInt()) {
			root.put("additions", Verifier.additionCount(alg));
		}
		root.put("min_additions", program.additions());
		// Honesty flag (CLAUDE.md optimality discipline): a contributor's reducer is a
		// heuristic too — an upper bound, never a proven minimum.
		root.put("additions_optimal", false);
		root.set("slp", SlpJson.toJson(MAPPER, program, SlpJson.OPTIMALITY_BOUND, method, source));
		MatrixJsonFormatter.write(f, root);
		return true;
	}

	/** The catalog file (any label, same orientation) whose matrices are {@code alg}'s, or null. */
	static File heldWithSameContent(NonCubicBilinearAlgorithm alg, int maxd) {
		String hash = SchemeIO.contentHash(alg);
		for (String category : new String[] { "known", "constructed", "curated" }) {
			File dir = SCHEMES_ROOT.resolve(category).resolve("section" + maxd).toFile();
			File[] files = dir.listFiles((d, name) -> name.endsWith(".json"));
			if (files == null) continue;
			for (File f : files) {
				try {
					JsonNode root = SchemeIO.parseJson(f);
					JsonNode n = root.get("n");
					if (n == null || !n.isArray() || n.size() != 3 || root.path("m").asInt(-1) != alg.r
							|| n.get(0).asInt() != alg.n || n.get(1).asInt() != alg.m || n.get(2).asInt() != alg.p
							|| SchemeIO.isStub(root) || root.path("complex").asBoolean(false)) {
						continue;
					}
					if (hash.equals(SchemeIO.contentHash(SchemeIO.readBilinear(f)))) {
						return f;
					}
				} catch (Exception e) {
					log.debug("heldWithSameContent: skipping {}: {}", f.getName(), e.toString());
				}
			}
		}
		return null;
	}

	private static String programMethod(UpstreamSchemeRepo repo) {
		return "contributor's addition-reduced program (" + repo.repo()
				+ "), replayed exactly against the catalog matrices on import by ImportUpstreamSchemeRepos";
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
