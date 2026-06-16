package eu.solven.matmul.docs.verify;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;

import lombok.extern.slf4j.Slf4j;

import tools.jackson.databind.JsonNode;

/**
 * READ-ONLY catalog keep-closure analysis (the catalog-build "what to keep / what to
 * delete" pass; deletion itself is a separate opt-in step).
 *
 * <p>The same dependency DAG {@code VerifyScheme} walks: every scheme's lineage refs
 * are edges to the SPECIFIC schemes it depends on. A scheme is KEPT iff it is reachable
 * from a <b>keep-root</b> — the best (lowest-rank) scheme for some shape — by following
 * those reference edges. Everything else (a worse-rank alternate that nothing pins) is
 * <b>DELETABLE</b>: removing it changes no best rank and breaks no pinned lineage.</p>
 *
 * <p>Edge resolution per lineage ref:</p>
 * <ul>
 *   <li>{@code shape@hash7} — a bit-exact pin → the exact file (must keep it, even if it
 *       is NOT the rank-best for its shape, e.g. a bud-rich projection/recombination base).</li>
 *   <li>named canonical ref ({@code perminov_ZT-2x6x7_m66_a797}, …) → that exact file.</li>
 *   <li>bare / {@code -direct} / {@code @sota} / {@code naive} — resolves to best-at-shape
 *       (already a keep-root) or the elementary scheme (no file); contributes no extra keep.</li>
 * </ul>
 *
 * <p>Writes NOTHING to the catalog. Emits the deletable list to
 * {@code target/catalog-deletable.txt} so a later step can act on it.</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.verify.CatalogKeepClosure</pre>
 */
@Slf4j
public final class CatalogKeepClosure {

	private CatalogKeepClosure() {}

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");
	/** {@code <shape>-r<rank>-<note>-<hash7>} — shape, rank and hash straight from the label. */
	private static final Pattern FN =
			Pattern.compile("^(\\d+x\\d+x\\d+)-r(\\d+)-.+-([0-9a-f]{4,})$");

	public static void main(String[] args) throws IOException {
		List<Path> files;
		try (Stream<Path> w = Files.walk(SCHEMES_ROOT)) {
			files = w.filter(p -> p.toString().endsWith(".json")).sorted().toList();
		}
		log.info("Keep-closure over {} scheme files", files.size());

		// Indices: pinned-key (shape@hash7) → file; canonical named key → file; and the
		// per-shape best (lowest rank) file = the keep-roots.
		Map<String, Path> byPinnedKey = new HashMap<>();
		Map<String, Path> byCanonKey = new HashMap<>();
		Map<String, Path> bestByShape = new HashMap<>();
		Map<String, Integer> bestRankByShape = new HashMap<>();
		Map<Path, String> shapeOf = new HashMap<>();

		for (Path p : files) {
			String stem = p.getFileName().toString().replaceFirst("\\.json$", "");
			byCanonKey.putIfAbsent(Lineage.canonicalKey(stem), p);
			Matcher m = FN.matcher(stem);
			if (!m.matches()) continue;            // non-canonical name: indexed for canon only
			String shape = m.group(1);
			int rank = Integer.parseInt(m.group(2));
			String hash7 = m.group(3).substring(0, Math.min(7, m.group(3).length()));
			byPinnedKey.put(shape + "@" + hash7, p);
			String sortedShape = sortShape(shape);
			shapeOf.put(p, sortedShape);
			Integer prev = bestRankByShape.get(sortedShape);
			if (prev == null || rank < prev) {
				bestRankByShape.put(sortedShape, rank);
				bestByShape.put(sortedShape, p);
			}
		}
		log.info("Indexed {} shapes; {} pinned-key entries", bestByShape.size(), byPinnedKey.size());

		// Edges: file → the SPECIFIC files its lineage references (pinned or named).
		Map<Path, List<Path>> edges = new HashMap<>();
		int danglingPins = 0;
		for (Path p : files) {
			List<Path> targets = new ArrayList<>();
			for (String ref : lineageRefs(p)) {
				Path t = resolveRefToFile(ref, byPinnedKey, byCanonKey);
				if (t != null) {
					targets.add(t);
				} else if (ref.contains("@") && !ref.endsWith("@naive") && !ref.endsWith("@sota")) {
					danglingPins++;   // a shape@hash pin whose file is gone (a real integrity bug)
				}
			}
			if (!targets.isEmpty()) edges.put(p, targets);
		}

		// BFS the keep-closure from the per-shape best schemes.
		Set<Path> keep = new HashSet<>();
		Deque<Path> queue = new ArrayDeque<>(bestByShape.values());
		keep.addAll(bestByShape.values());
		while (!queue.isEmpty()) {
			for (Path t : edges.getOrDefault(queue.poll(), List.of())) {
				if (keep.add(t)) queue.add(t);
			}
		}

		int keptBest = bestByShape.size();
		int keptReferenced = keep.size() - keptBest;
		List<Path> deletable = new ArrayList<>();
		for (Path p : files) if (!keep.contains(p)) deletable.add(p);

		// Write the deletable list (NO deletion here).
		Path out = Path.of("target", "catalog-deletable.txt");
		Files.createDirectories(out.getParent());
		StringBuilder sb = new StringBuilder();
		sb.append("# ").append(deletable.size()).append(" deletable scheme files (not best-per-shape,")
				.append(" not referenced by any kept scheme). Generated by CatalogKeepClosure — review")
				.append(" before any deletion.\n");
		for (Path p : deletable) sb.append(SCHEMES_ROOT.relativize(p)).append('\n');
		Files.writeString(out, sb.toString());

		// Deletable breakdown by section dir (where the cruft lives).
		Map<String, Integer> byDir = new TreeMap<>();
		for (Path p : deletable) {
			Path rel = SCHEMES_ROOT.relativize(p);
			String dir = rel.getNameCount() >= 2 ? rel.subpath(0, 2).toString() : rel.toString();
			byDir.merge(dir, 1, Integer::sum);
		}

		log.info("================ KEEP-CLOSURE ================");
		log.info("  total files        : {}", files.size());
		log.info("  KEEP               : {}  ({} best-per-shape + {} referenced)",
				keep.size(), keptBest, keptReferenced);
		log.info("  DELETABLE          : {}  → target/catalog-deletable.txt", deletable.size());
		log.info("  dangling @hash pins: {}  (referenced file missing — integrity bug, NOT deletable)",
				danglingPins);
		log.info("  -- deletable by dir (top 12) --");
		byDir.entrySet().stream()
				.sorted((a, b) -> Integer.compare(b.getValue(), a.getValue()))
				.limit(12)
				.forEach(e -> log.info("     {}  {}", e.getValue(), e.getKey()));
	}

	/** Resolve a lineage ref to the SPECIFIC file it pins (pinned {@code shape@hash} or a
	 *  named canonical key), or {@code null} for best-at-shape / naive / unresolved refs. */
	private static Path resolveRefToFile(String ref, Map<String, Path> byPinnedKey,
			Map<String, Path> byCanonKey) {
		int at = ref.indexOf('@');
		if (at > 0) {
			String shape = ref.substring(0, at);
			String tail = ref.substring(at + 1);
			if (tail.equals("naive") || tail.equals("sota")) return null; // not a specific file
			String hash7 = tail.substring(0, Math.min(7, tail.length()));
			Path direct = byPinnedKey.get(shape + "@" + hash7);
			if (direct != null) return direct;
			return null; // dangling pin (counted by the caller)
		}
		// A named canonical ref (perminov_…, alphatensor_…) → exact file via canonical key.
		Path named = byCanonKey.get(Lineage.canonicalKey(ref));
		return named; // null for bare-shape / -direct / naive refs (best-at-shape, no extra keep)
	}

	private static List<String> lineageRefs(Path f) {
		List<String> out = new ArrayList<>();
		try {
			collectRefs(SchemeIO.parseJson(f.toFile()).get("lineage"), out);
		} catch (Exception e) {
			// unreadable lineage → no edges (a separate verify concern)
		}
		return out;
	}

	private static void collectRefs(JsonNode node, List<String> out) {
		if (node == null) return;
		if (node.isObject()) {
			JsonNode ref = node.get("ref");
			if (ref != null && ref.isString()) out.add(ref.asString());
			for (var e : node.properties()) collectRefs(e.getValue(), out);
		} else if (node.isArray()) {
			for (JsonNode c : node) collectRefs(c, out);
		}
	}

	private static String sortShape(String shape) {
		String[] parts = shape.split("x");
		int[] d = { Integer.parseInt(parts[0]), Integer.parseInt(parts[1]), Integer.parseInt(parts[2]) };
		java.util.Arrays.sort(d);
		return d[0] + "x" + d[1] + "x" + d[2];
	}
}
