package eu.solven.matmul.docs.verify;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Stream;

import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.LineageReplayer;
import eu.solven.matmul.search.LineageReplayer.LeafClass;
import eu.solven.matmul.verifiers.LineageVerifier;

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
 * every {@code {"op":"Atom","ref":…}} node). <b>The classifier is the replayer's own
 * leaf grammar</b> ({@link LineageReplayer#classifyLeaf}), evaluated under the
 * file's own {@code fields[]}: a lookup for one field only sees the files valid over
 * it, so an F₂ scheme's pin to the AlphaTensor F₂ ⟨4,4,4⟩=47 is dangling for a
 * char-0 lookup and perfectly fine for an F₂ one. (Until 2026-10-18 the audit kept
 * its own regexes under one C lookup and reported 17 859 UNRESOLVABLE refs that all
 * replayed — canonical filename stems, {@code TA_lita(n=N)},
 * {@code Strassen<2,2,2>=7 :: CANONICAL} — plus four "dangling" F₂ pins.)</p>
 * <ul>
 *   <li><b>PINNED_OK</b> — {@code shape@hash} resolves to a materialised file.</li>
 *   <li><b>PINNED_OK_STUB</b> — {@code shape@hash} resolves to a stub's stamped hash.</li>
 *   <li><b>PINNED_DANGLING</b> — pinned hash matches nothing under any of the file's
 *       fields: the replay throws. The repair set ({@code RepinDanglingLineageRefs}).</li>
 *   <li><b>BARE</b> — shape-only ref; resolves to catalog-best <i>today</i>, which
 *       may differ from build time (the stub-rank-drift mechanism).</li>
 *   <li><b>LEGACY_STEM</b> — any other ref a shape can be pulled out of (old
 *       {@code source-NxMxP_m…} stems, pre-2026-10-18 canonical stems
 *       {@code NxMxP-r…-note-hash7}, {@code Name<N,M,P>=R :: …} labels): resolved by
 *       shape, same risk class as BARE.</li>
 *   <li><b>NAMED / PARAMETRIC / NAIVE</b> — pool names, formulas, naive blocks;
 *       resolved from code, not files.</li>
 *   <li><b>UNRESOLVABLE</b> — nothing at the shape at all, or no shape in the ref.</li>
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

	/** One replayer per field, built on demand (each lookup indexes the catalog once). */
	private static final Map<Field, LineageReplayer> REPLAYERS = new EnumMap<>(Field.class);

	static synchronized LineageReplayer replayerFor(Field f) {
		return REPLAYERS.computeIfAbsent(f, k -> LineageReplayer.withDefaultPool(new FieldAwareLookup(k)));
	}

	/**
	 * The lookups a file's leaves must be resolved under: C covers every characteristic-0
	 * tag (its fallback chain is C ⊃ R ⊃ Q ⊃ Z), F₂ and F₃ are their own universes. A file
	 * without {@code fields[]} (a data bug) is tried under all three.
	 */
	public static List<Field> fieldsToTry(List<String> tags) {
		List<Field> out = new ArrayList<>();
		boolean char0 = tags.stream().anyMatch(t -> t.equals("Z") || t.equals("Q") || t.equals("R")
				|| t.equals("C") || t.equals("ZT"));
		if (char0 || tags.isEmpty()) out.add(Field.C);
		if (tags.contains("F2") || tags.contains("Z2") || tags.isEmpty()) out.add(Field.F2);
		if (tags.contains("F3") || tags.contains("Z3") || tags.isEmpty()) out.add(Field.F3);
		return out;
	}

	public static void main(String[] args) throws Exception {
		long start = System.nanoTime();

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

			List<LineageReplayer> replayers = new ArrayList<>();
			for (Field f : fieldsToTry(SchemeIO.fieldTags(root))) replayers.add(replayerFor(f));

			List<String> refs = new ArrayList<>();
			collectAtomRefs(lineage, refs);
			boolean dangling = false;
			for (String ref : refs) {
				String cls = classify(ref, replayers);
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
		LineageVerifier verifier = new LineageVerifier(replayerFor(Field.C));
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

	/**
	 * Classify a ref under one or more replayers (one per field the file is valid over):
	 * the first replayer that resolves it decides; otherwise the "least broken" verdict
	 * (a dangling pin over a shapeless unresolvable).
	 */
	static String classify(String ref, List<LineageReplayer> replayers) {
		LeafClass best = null;
		for (LineageReplayer r : replayers) {
			LeafClass c = r.classifyLeaf(ref).cls();
			if (c.resolves()) return label(c);
			if (best == null || (c == LeafClass.DANGLING_PIN && best == LeafClass.UNRESOLVABLE)) best = c;
		}
		return label(best == null ? LeafClass.UNRESOLVABLE : best);
	}

	/** Single-replayer convenience (tests, ad-hoc probes). */
	static String classify(String ref, LineageReplayer replayer) {
		return classify(ref, List.of(replayer));
	}

	/** The audit's label for a replayer leaf class (kept stable for the report tables). */
	static String label(LeafClass c) {
		return switch (c) {
			case PINNED -> "PINNED_OK";
			case PINNED_STUB -> "PINNED_OK_STUB";
			case DANGLING_PIN -> "PINNED_DANGLING";
			case SHAPE_BEST -> "BARE";
			case LEGACY_SHAPE_BEST -> "LEGACY_STEM";
			case NAMED, PARAMETRIC, NAIVE, UNRESOLVABLE -> c.name();
		};
	}
}
