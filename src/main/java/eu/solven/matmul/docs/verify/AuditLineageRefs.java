package eu.solven.matmul.docs.verify;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.atomic.AtomicLong;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.LineageVerifier;

import lombok.extern.slf4j.Slf4j;

import tools.jackson.databind.JsonNode;

/**
 * READ-ONLY lineage-integrity audit over the whole catalog (known/derived/curated).
 *
 * <p>Motivated by the 2026-06 purge/renames: ~2.7k derived files carry
 * {@code lineage_str} display strings referencing purged {@code fmm_lille-*}
 * filenames, and structured refs may dangle (pinned {@code shape@hash} whose
 * hash no longer exists) or drift (bare {@code NxMxP} refs resolve to
 * catalog-best, which changed under the purge — the bit-exact-replay risk).
 * This auditor quantifies the damage WITHOUT writing anything, so the repair
 * (re-pin / regenerate / quarantine) can be sized before it is run.</p>
 *
 * <p>Per-ref classification (refs are read straight from the lineage JSON —
 * every {@code {"op":"Atom","ref":…}} node):</p>
 * <ul>
 *   <li><b>PINNED_OK</b> — {@code shape@hash} resolves via the catalog index.</li>
 *   <li><b>PINNED_DANGLING</b> — pinned hash matches nothing (purged/re-hashed):
 *       replay is broken or silently re-targeted. The repair set.</li>
 *   <li><b>BARE</b> — shape-only ref; resolves to catalog-best <i>today</i>, which
 *       may differ from build time.</li>
 *   <li><b>LEGACY_SOURCE_PREFIXED</b> — old {@code source-NxMxP_m…} stem; the
 *       replayer falls back to shape resolution (same risk class as BARE) and the
 *       display string is stale.</li>
 *   <li><b>NAMED / PARAMETRIC / NAIVE</b> — pool names, formulas, naive blocks;
 *       resolved from code, not files.</li>
 *   <li><b>UNRESOLVABLE</b> — nothing at the shape at all.</li>
 * </ul>
 *
 * <p>Also: counts stale {@code lineage_str} mentions of purged stems, and runs a
 * sampled {@link LineageVerifier#verifyFile} pass to estimate replay health.</p>
 */
@Slf4j
public class AuditLineageRefs {

	private static final Path ROOT = Path.of("src/main/resources/schemes");
	/** Sampled-verify only replays shapes up to this max dim — replaying ⟨30,30,30⟩
	 *  stubs is exactly the off-heap blowup the long-running discipline warns about. */
	private static final int VERIFY_MAX_DIM = 12;
	private static final Pattern PINNED = Pattern.compile("^(\\d+)x(\\d+)x(\\d+)@([0-9a-f]{6,})$");
	private static final Pattern BARE = Pattern.compile("^(naive-)?(\\d+)x(\\d+)x(\\d+)(-direct)?$");
	private static final Pattern NAMED = Pattern.compile("^[A-Za-z].*<\\d+,\\d+,\\d+>=\\d+$");
	private static final Pattern PARAMETRIC = Pattern.compile("^[A-Za-z0-9]+\\(.*\\)$");
	/** Legacy source-prefixed stem, e.g. {@code fmm_lille-2x3x4_m11_a40} or {@code smirnov13-3x3x6_m40}. */
	private static final Pattern LEGACY = Pattern.compile("^[a-zA-Z][\\w.]*[-_](\\d+)x(\\d+)x(\\d+).*$");

	public static void main(String[] args) throws Exception {
		long start = System.nanoTime();
		FieldAwareLookup lookup = new FieldAwareLookup("C"); // broadest: indexes everything

		List<Path> all;
		try (Stream<Path> walk = Files.walk(ROOT)) {
			all = walk.filter(p -> p.toString().endsWith(".json")).sorted().toList();
		}
		log.info("Auditing lineage refs across {} scheme files", all.size());

		Map<String, AtomicLong> refClass = new TreeMap<>();
		Map<String, List<String>> samples = new TreeMap<>();
		long withLineage = 0, staleLineageStr = 0, parseErrors = 0;
		List<String> danglingFiles = new ArrayList<>();
		List<Path> lineageFiles = new ArrayList<>();

		int done = 0;
		for (Path p : all) {
			done++;
			if (done % 2000 == 0) {
				log.info("[progress] {}/{} files scanned, {} with lineage, {} dangling-pinned so far, {} ms",
						done, all.size(), withLineage, count(refClass, "PINNED_DANGLING"),
						(System.nanoTime() - start) / 1_000_000);
			}
			// RSS discipline: do NOT Jackson-parse every file — the dense imports
			// (sections 16–32) are multi-MB and tree-parsing them got this audit
			// SIGKILL'ed by macOS memory pressure (exit 144, the CLAUDE.md off-heap
			// pattern). A substring probe on the raw text decides whether a parse
			// is needed at all; the lineage-bearing files (stubs/derived) are small.
			String raw;
			try {
				raw = Files.readString(p);
			} catch (Exception e) {
				parseErrors++;
				continue;
			}
			// Stale display string? (purged fmm_lille stems are the known case)
			if (raw.contains("fmm_lille")) {
				staleLineageStr++;
			}
			if (!raw.contains("\"lineage\"")) {
				continue;
			}
			JsonNode root;
			try {
				root = SchemeIO.parseJson(p.toFile());
			} catch (Exception e) {
				parseErrors++;
				continue;
			}
			JsonNode lineage = root.get("lineage");
			if (lineage == null || lineage.isNull()) continue;
			withLineage++;
			// Track shape so the sampled-verify phase can skip heavyweight replays.
			int maxDim = 0;
			JsonNode shapeArr = root.get("n");
			if (shapeArr != null && shapeArr.isArray()) {
				for (JsonNode d : shapeArr) maxDim = Math.max(maxDim, d.asInt());
			}
			if (maxDim <= VERIFY_MAX_DIM) lineageFiles.add(p);

			List<String> refs = new ArrayList<>();
			collectAtomRefs(lineage, refs);
			boolean dangling = false;
			for (String ref : refs) {
				String cls = classify(ref, lookup);
				refClass.computeIfAbsent(cls, k -> new AtomicLong()).incrementAndGet();
				List<String> s = samples.computeIfAbsent(cls, k -> new ArrayList<>());
				if (s.size() < 8) s.add(ref + "   [in " + p.getFileName() + "]");
				if (cls.equals("PINNED_DANGLING") || cls.equals("UNRESOLVABLE")) dangling = true;
			}
			if (dangling) danglingFiles.add(p.toString());
		}

		log.info("==================================================================");
		log.info(" Lineage-ref audit ({} files, {} with structured lineage, {} parse errors)",
				all.size(), withLineage, parseErrors);
		log.info("==================================================================");
		for (var e : refClass.entrySet()) {
			log.info("  {}  ×{}", String.format("%-24s", e.getKey()), e.getValue().get());
			for (String s : samples.getOrDefault(e.getKey(), List.of())) {
				log.info("      e.g. {}", s);
			}
		}
		log.info("  stale lineage_str (mentions purged fmm_lille stems): {}", staleLineageStr);
		log.info("  files with ≥1 dangling/unresolvable ref: {}", danglingFiles.size());
		for (String f : danglingFiles.stream().limit(20).toList()) {
			log.info("      {}", f);
		}
		if (danglingFiles.size() > 20) log.info("      … and {} more", danglingFiles.size() - 20);

		// Sampled compositional replay-health check (maxDim ≤ VERIFY_MAX_DIM only).
		int sampleEvery = Math.max(1, lineageFiles.size() / 150);
		log.info("Sampled compositional verification (every {}th of {} lineage-bearing files with maxDim ≤ {}):",
				sampleEvery, lineageFiles.size(), VERIFY_MAX_DIM);
		LineageVerifier verifier = new LineageVerifier(lookup);
		long vOk = 0, vFail = 0;
		List<String> vFails = new ArrayList<>();
		for (int i = 0; i < lineageFiles.size(); i += sampleEvery) {
			Path p = lineageFiles.get(i);
			try {
				LineageVerifier.Result r = verifier.verifyFile(p.toFile());
				if (r.certified()) vOk++;
				else { vFail++; if (vFails.size() < 12) vFails.add(p.getFileName() + " — " + r.detail()); }
			} catch (Exception | StackOverflowError e) {
				vFail++;
				String msg = String.valueOf(e.getMessage());
				if (vFails.size() < 12) vFails.add(p.getFileName() + " — " + e.getClass().getSimpleName()
						+ ": " + msg.substring(0, Math.min(140, msg.length())));
			}
			if ((vOk + vFail) % 30 == 0) {
				log.info("[progress] sampled-verify {}: {} OK, {} FAILED, {} ms",
						vOk + vFail, vOk, vFail, (System.nanoTime() - start) / 1_000_000);
			}
		}
		log.info("  sampled verify: {} OK, {} FAILED", vOk, vFail);
		for (String f : vFails) log.info("      FAIL {}", f);
		log.info("Done in {} ms (read-only — nothing written).", (System.nanoTime() - start) / 1_000_000);
	}

	private static long count(Map<String, AtomicLong> m, String k) {
		AtomicLong v = m.get(k);
		return v == null ? 0 : v.get();
	}

	/** Collect every {@code {"op":"Atom","ref":…}} ref in a lineage JSON tree. */
	public static void collectAtomRefs(JsonNode node, List<String> out) {
		if (node == null) return;
		if (node.isObject()) {
			JsonNode op = node.get("op");
			JsonNode ref = node.get("ref");
			if (op != null && "Atom".equals(op.asString()) && ref != null && ref.isTextual()) {
				out.add(ref.asString());
			}
			node.properties().forEach(e -> collectAtomRefs(e.getValue(), out));
		} else if (node.isArray()) {
			node.forEach(child -> collectAtomRefs(child, out));
		}
	}

	/** stub file → stamped hash (classification cache). */
	private static final Map<java.nio.file.Path, String> STAMPED = new java.util.concurrent.ConcurrentHashMap<>();

	static String classify(String ref, FieldAwareLookup lookup) {
		Matcher m = PINNED.matcher(ref);
		if (m.matches()) {
			int n = Integer.parseInt(m.group(1));
			int mm = Integer.parseInt(m.group(2));
			int p = Integer.parseInt(m.group(3));
			if (lookup.findByHash(n, mm, p, m.group(4)).isPresent()) return "PINNED_OK";
			// findByHash only hashes materialised candidates; a pin to a STUB's hash
			// resolves via its stamped "hash" field (LineageReplayer does the same) —
			// without this check the auditor misclassifies valid stub pins as dangling.
			String hash = m.group(4);
			for (java.nio.file.Path c : lookup.findFiles(n, mm, p)) {
				String stamped = STAMPED.computeIfAbsent(c, f -> {
					try {
						String h = SchemeIO.readHash(SchemeIO.parseJson(f.toFile()));
						return h == null ? "" : h;
					} catch (Exception e) {
						return "";
					}
				});
				if (!stamped.isEmpty() && stamped.startsWith(hash)) return "PINNED_OK_STUB";
			}
			return "PINNED_DANGLING";
		}
		Matcher b = BARE.matcher(ref);
		if (b.matches()) {
			if (b.group(1) != null) return "NAIVE"; // naive-NxMxP — constructed, always resolvable
			int n = Integer.parseInt(b.group(2));
			int mm = Integer.parseInt(b.group(3));
			int p = Integer.parseInt(b.group(4));
			return lookup.findFiles(n, mm, p).isEmpty() ? "UNRESOLVABLE" : "BARE";
		}
		if (NAMED.matcher(ref).matches()) return "NAMED";
		if (PARAMETRIC.matcher(ref).matches()) return "PARAMETRIC";
		Matcher l = LEGACY.matcher(ref);
		if (l.matches()) {
			int n = Integer.parseInt(l.group(1));
			int mm = Integer.parseInt(l.group(2));
			int p = Integer.parseInt(l.group(3));
			return lookup.findFiles(n, mm, p).isEmpty() ? "UNRESOLVABLE" : "LEGACY_SOURCE_PREFIXED";
		}
		return "UNRESOLVABLE";
	}
}
