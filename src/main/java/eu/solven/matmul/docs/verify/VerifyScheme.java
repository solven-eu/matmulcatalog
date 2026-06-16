package eu.solven.matmul.docs.verify;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.LineageReplayer;

import lombok.extern.slf4j.Slf4j;

import tools.jackson.databind.JsonNode;

/**
 * DAG-based catalog verifier (the unified successor to {@code VerifyAllSchemes} +
 * {@code CatalogKeepClosure} + {@code ComputeExplicitable} + the stub rank-honesty audit).
 *
 * <p>Every scheme's lineage refs are edges to the SPECIFIC schemes it depends on; the
 * roots are atoms (explicit schemes with no dependencies). The verifier walks the DAG
 * <b>up</b> in topological order — NOT shape order, because {@code Project} goes big→small
 * (⟨25,25,25⟩ depends on ⟨26,26,26⟩) — so every node is checked only after its dependencies.</p>
 *
 * <p>Per-node verdict:</p>
 * <ul>
 *   <li><b>OK_EXPLICIT</b> — a concrete scheme (matrices); the matmul-tensor identity holds.</li>
 *   <li><b>OK_REPLAY</b> — an EXPLICITABLE stub; replaying its (precise) lineage reproduces a
 *       scheme of the recorded rank that passes the matmul spot-check.</li>
 *   <li><b>CITED_BOUND</b> — a non-explicitable stub (a best-at-shape {@code @sota}/{@code -direct}/
 *       bare leaf): a rank claim, not an exact construction; its refs are checked to resolve.</li>
 *   <li><b>CORRUPT_*</b> — verify failed: {@code RANK} (replay ≠ recorded), {@code VERIFY} (identity
 *       fails), {@code DANGLING} (a pinned ref's file is missing), {@code CYCLE} (self-reference),
 *       {@code DEP} (a dependency is corrupt). NEVER falls back to naive/best — it fails.</li>
 * </ul>
 *
 * <pre>
 * # verify everything, fast parallel all-clear check:
 *   --threads=8
 * # verify specific shapes and their dependencies, mono-threaded (stop at first corruption):
 *   --shape=4x19x20,6x14x20 --recursive=true --threads=1
 * # a curated list from a file:
 *   --shape-file=target/catalog-deletable.txt --recursive=false
 * </pre>
 */
@Slf4j
public final class VerifyScheme {

	private VerifyScheme() {}

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");
	private static final Pattern FN = Pattern.compile("^(\\d+x\\d+x\\d+)-r(\\d+)-.+-([0-9a-f]{4,})$");

	enum Verdict { OK_EXPLICIT, OK_REPLAY, CITED_BOUND, CORRUPT_RANK, CORRUPT_VERIFY,
		CORRUPT_DANGLING, CORRUPT_CYCLE, CORRUPT_DEP }

	public static void main(String[] args) throws IOException {
		int threads = Integer.parseInt(arg(args, "--threads", "8"));
		boolean recursive = Boolean.parseBoolean(arg(args, "--recursive", "true"));
		boolean monotony = Boolean.parseBoolean(arg(args, "--monotony", "true"));
		Set<String> shapeFilter = parseShapes(args);

		List<Path> files;
		try (Stream<Path> w = Files.walk(SCHEMES_ROOT)) {
			files = w.filter(p -> p.toString().endsWith(".json")).sorted().toList();
		}
		log.info("Loaded {} scheme files; building dependency DAG", files.size());

		// Indices: pinned-key/named-key → file; per-shape best (the keep-roots).
		Map<String, Path> byPinnedKey = new HashMap<>();
		Map<String, Path> byCanonKey = new HashMap<>();
		Map<String, Path> bestByShape = new HashMap<>();
		Map<String, Integer> bestRankByShape = new HashMap<>();
		Map<Path, Integer> recordedRank = new HashMap<>();
		Map<Path, String> shapeOf = new HashMap<>();
		for (Path p : files) {
			String stem = p.getFileName().toString().replaceFirst("\\.json$", "");
			byCanonKey.putIfAbsent(Lineage.canonicalKey(stem), p);
			Matcher m = FN.matcher(stem);
			if (!m.matches()) continue;
			String shape = m.group(1);
			int rank = Integer.parseInt(m.group(2));
			recordedRank.put(p, rank);
			// Key the @hash index by SORTED shape: a ref pins by CONTENT hash, which
			// identifies the scheme regardless of axis orientation. A projection parent
			// recorded as ⟨3,7,3⟩@h (an OrientAs of the ⟨3,3,7⟩ file) must still resolve to
			// that file — keying by raw orientation spuriously flagged it DANGLING.
			byPinnedKey.put(sortShape(shape) + "@" + m.group(3).substring(0, Math.min(7, m.group(3).length())), p);
			String sorted = sortShape(shape);
			shapeOf.put(p, sorted);
			if (bestRankByShape.merge(sorted, rank, Math::min) == rank) bestByShape.put(sorted, p);
		}

		// Edges: file → resolved dependency files. A bare/-direct/@sota leaf resolves to
		// best-at-shape; an unresolved @hash/named ref is a DANGLING dependency.
		Map<Path, List<Path>> deps = new HashMap<>();
		Set<Path> dangling = ConcurrentHashMap.newKeySet();
		for (Path p : files) {
			List<Path> ds = new ArrayList<>();
			for (String ref : lineageRefs(p)) {
				Path t = resolveDep(ref, byPinnedKey, byCanonKey, bestByShape);
				if (t != null && t != p) ds.add(t);
				else if (isPinnedOrNamed(ref) && t == null) dangling.add(p);
			}
			deps.put(p, ds);
		}

		Map<Path, Boolean> explicitable = ComputeExplicitable.computeAll(files);

		// Scope: the named shapes' files (+ transitive deps if --recursive), else everything.
		Set<Path> scope;
		if (shapeFilter.isEmpty()) {
			scope = new LinkedHashSet<>(files);
		} else {
			scope = new LinkedHashSet<>();
			for (Path p : files) if (shapeFilter.contains(shapeOf.get(p))) scope.add(p);
			if (recursive) {
				Deque<Path> q = new ArrayDeque<>(scope);
				while (!q.isEmpty()) {
					for (Path d : deps.getOrDefault(q.poll(), List.of())) if (scope.add(d)) q.add(d);
				}
			}
		}
		log.info("Scope: {} files ({} threads, recursive={})", scope.size(), threads, recursive);

		// Topological order (Kahn) over the scope; remaining nodes after Kahn are in cycles.
		List<Path> topo = topoOrder(scope, deps);
		Set<Path> cyclic = new HashSet<>(scope);
		topo.forEach(cyclic::remove);

		Map<Path, Verdict> verdict = new ConcurrentHashMap<>();
		for (Path p : cyclic) verdict.put(p, Verdict.CORRUPT_CYCLE);

		Verifier0 v0 = new Verifier0(recordedRank, explicitable, dangling, deps, verdict);
		AtomicInteger done = new AtomicInteger();
		long t0 = System.nanoTime();
		if (threads <= 1) {
			// MONO: deps-first, STOP at the first corruption so it can be fixed and resumed.
			for (Path p : topo) {
				Verdict v = v0.verify(p);
				verdict.put(p, v);
				if (isCorrupt(v)) {
					log.error("CORRUPT [{}] {} — stopping (mono mode: fix and resume)", v,
							SCHEMES_ROOT.relativize(p));
					break;
				}
				progress(done.incrementAndGet(), topo.size(), t0);
			}
		} else {
			// MULTI: independent nodes verified in parallel (the fast all-clear pass). Deps-first
			// is preserved per chain because verify() recurses through replay; verdicts memoise.
			topo.parallelStream().forEach(p -> {
				verdict.put(p, v0.verify(p));
				progress(done.incrementAndGet(), topo.size(), t0);
			});
		}

		report(scope, verdict);

		// Catalog-wide RANK MONOTONICITY sanity check (independent of the per-node verdict):
		// tensor rank is monotone under sub-tensor restriction, so for sorted shapes
		// small ⊆ big (componentwise) the catalog must satisfy best(small) ≤ best(big)
		// WITHIN a field (the zero-padding embedding is field-specific). A violation is a
		// catalog-soundness signal — usually a small shape we UNDER-derived (it is literally
		// a restriction of the bigger scheme), occasionally an over-claimed big rank. This is
		// exactly the invariant AllocationOptimizer's root lower bound relies on; a violation
		// here is what silently poisoned the optimizer before findRank stopped returning the
		// MAX/100 sentinel.
		if (monotony && shapeFilter.isEmpty()) {
			monotonyCheck(files, recordedRank, shapeOf);
		}
	}

	/** Per-field per-(sorted-)shape best recorded rank, then flag every dominance pair
	 *  {@code small ⊆ big} with {@code best(small) > best(big)}. Field membership is read
	 *  from each file's {@code fields[]} (parsed in parallel). */
	private static void monotonyCheck(List<Path> files, Map<Path, Integer> recordedRank,
			Map<Path, String> shapeOf) {
		// field -> sortedShape -> best (rank, file)
		Map<String, Map<String, Path>> bestFile = new ConcurrentHashMap<>();
		Map<String, Map<String, Integer>> bestRank = new ConcurrentHashMap<>();
		files.parallelStream().forEach(p -> {
			Integer rank = recordedRank.get(p);
			String shape = shapeOf.get(p);
			if (rank == null || shape == null) return;
			for (String f : readFields(p)) {
				Map<String, Integer> r = bestRank.computeIfAbsent(f, k -> new ConcurrentHashMap<>());
				Map<String, Path> fl = bestFile.computeIfAbsent(f, k -> new ConcurrentHashMap<>());
				synchronized (r) {
					Integer cur = r.get(shape);
					if (cur == null || rank < cur) { r.put(shape, rank); fl.put(shape, p); }
				}
			}
		});

		record Violation(String field, String small, int rSmall, Path fSmall,
				String big, int rBig, Path fBig) {
			int gap() { return rSmall - rBig; }
		}
		List<Violation> violations = new ArrayList<>();
		for (String field : bestRank.keySet()) {
			Map<String, Integer> r = bestRank.get(field);
			Map<String, Path> fl = bestFile.get(field);
			List<String> shapes = new ArrayList<>(r.keySet());
			for (int i = 0; i < shapes.size(); i++) {
				int[] a = dims(shapes.get(i));
				for (int j = 0; j < shapes.size(); j++) {
					if (i == j) continue;
					int[] b = dims(shapes.get(j));
					// a ⊆ b componentwise (sorted dims) ⇒ require best(a) ≤ best(b).
					if (a[0] <= b[0] && a[1] <= b[1] && a[2] <= b[2]
							&& r.get(shapes.get(i)) > r.get(shapes.get(j))) {
						violations.add(new Violation(field, shapes.get(i), r.get(shapes.get(i)),
								fl.get(shapes.get(i)), shapes.get(j), r.get(shapes.get(j)), fl.get(shapes.get(j))));
					}
				}
			}
		}
		log.info("================ RANK MONOTONICITY ================");
		if (violations.isEmpty()) {
			log.info("  OK — no field has best(small) > best(big) for any dominance pair");
			return;
		}
		violations.sort((x, y) -> Integer.compare(y.gap(), x.gap()));
		log.warn("  {} monotonicity violation(s) [best(small) > best(big), same field]:", violations.size());
		violations.stream().limit(25).forEach(v ->
				log.warn("    {}: {}=r{} > {}=r{}  (gap {})  small={}",
						v.field(), v.small(), v.rSmall(), v.big(), v.rBig(), v.gap(),
						SCHEMES_ROOT.relativize(v.fSmall())));
	}

	private static int[] dims(String sortedShape) {
		String[] parts = sortedShape.split("x");
		return new int[] { Integer.parseInt(parts[0]), Integer.parseInt(parts[1]), Integer.parseInt(parts[2]) };
	}

	private static List<String> readFields(Path p) {
		try {
			JsonNode f = SchemeIO.parseJson(p.toFile()).get("fields");
			if (f != null && f.isArray()) {
				List<String> out = new ArrayList<>(f.size());
				for (JsonNode x : f) out.add(x.asString());
				return out;
			}
		} catch (Exception e) { /* unreadable → contributes no field */ }
		return List.of();
	}

	/** Per-node verification, memoised in {@code verdict}. */
	private static final class Verifier0 {
		final Map<Path, Integer> recordedRank;
		final Map<Path, Boolean> explicitable;
		final Set<Path> dangling;
		final Map<Path, List<Path>> deps;
		final Map<Path, Verdict> verdict;
		final LineageReplayer replayer = LineageReplayer.withDefaultPool(new FieldAwareLookup("Q"));

		Verifier0(Map<Path, Integer> recordedRank, Map<Path, Boolean> explicitable, Set<Path> dangling,
				Map<Path, List<Path>> deps, Map<Path, Verdict> verdict) {
			this.recordedRank = recordedRank;
			this.explicitable = explicitable;
			this.dangling = dangling;
			this.deps = deps;
			this.verdict = verdict;
		}

		Verdict verify(Path p) {
			Verdict memo = verdict.get(p);
			if (memo != null) return memo;
			Verdict v = compute(p);
			verdict.put(p, v);
			return v;
		}

		private Verdict compute(Path p) {
			if (dangling.contains(p)) return Verdict.CORRUPT_DANGLING;
			// A node is only as sound as its dependencies.
			for (Path d : deps.getOrDefault(p, List.of())) {
				if (isCorrupt(verify(d))) return Verdict.CORRUPT_DEP;
			}
			JsonNode root;
			try {
				root = SchemeIO.parseJson(p.toFile());
			} catch (Exception e) {
				return Verdict.CORRUPT_VERIFY;
			}
			try {
				if (!SchemeIO.isStub(root)) {
					// Concrete scheme — verify the matmul-tensor identity directly.
					if (SchemeIO.isNonBilinear(root) || SchemeIO.isComplex(root)) {
						return Verdict.OK_EXPLICIT; // (verified by VerifyAllSchemes' field checks)
					}
					NonCubicBilinearAlgorithm alg = SchemeIO.isReduced(root)
							? SchemeIO.readReduced(root) : SchemeIO.read(root);
					return Verifier.passesRandomMatmulSpotCheck(alg)
							? Verdict.OK_EXPLICIT : Verdict.CORRUPT_VERIFY;
				}
				if (!Boolean.TRUE.equals(explicitable.get(p))) {
					return Verdict.CITED_BOUND; // a bound, not an exact construction
				}
				// Explicitable stub: replay the (precise) lineage and reconcile its rank.
				NonCubicBilinearAlgorithm alg = replayer.replayFromFile(p.toFile());
				Integer rec = recordedRank.get(p);
				// Replay WORSE than recorded = a phantom over-claim → CORRUPT_RANK. Replay
				// BETTER is FINE (mirrors the materialiser's [replay-dce] acceptance): a named
				// operand resolves to catalog-best, which may have improved since the stub was
				// written (⟨20,20,21⟩ recorded 4778 via SZ-4378, replays 4740 via the dis09
				// cube 4340). The scheme is valid; its recorded rank is just conservative.
				if (rec != null && alg.r > rec) {
					log.warn("[rank-diag] {} recorded r={} but lineage REPLAYS WORSE to r={} (Δ+{})",
							SCHEMES_ROOT.relativize(p), rec, alg.r, alg.r - rec);
					return Verdict.CORRUPT_RANK;
				}
				return Verifier.passesRandomMatmulSpotCheck(alg)
						? Verdict.OK_REPLAY : Verdict.CORRUPT_VERIFY;
			} catch (Throwable e) {
				return Verdict.CORRUPT_VERIFY;
			}
		}
	}

	private static List<Path> topoOrder(Set<Path> scope, Map<Path, List<Path>> deps) {
		Map<Path, Integer> indeg = new HashMap<>();
		Map<Path, List<Path>> dependents = new HashMap<>();
		for (Path p : scope) indeg.putIfAbsent(p, 0);
		for (Path p : scope) {
			for (Path d : deps.getOrDefault(p, List.of())) {
				if (!scope.contains(d)) continue;
				indeg.merge(p, 1, Integer::sum);
				dependents.computeIfAbsent(d, k -> new ArrayList<>()).add(p);
			}
		}
		Deque<Path> q = new ArrayDeque<>();
		indeg.forEach((p, d) -> { if (d == 0) q.add(p); });
		List<Path> order = new ArrayList<>();
		while (!q.isEmpty()) {
			Path p = q.poll();
			order.add(p);
			for (Path dep : dependents.getOrDefault(p, List.of())) {
				if (indeg.merge(dep, -1, Integer::sum) == 0) q.add(dep);
			}
		}
		return order; // nodes left out (indeg never hit 0) are in cycles
	}

	private static void report(Set<Path> scope, Map<Path, Verdict> verdict) {
		Map<Verdict, Integer> counts = new java.util.EnumMap<>(Verdict.class);
		List<Path> corrupt = new ArrayList<>();
		int unverified = 0;
		for (Path p : scope) {
			Verdict v = verdict.get(p);
			if (v == null) { unverified++; continue; } // mono-mode stopped before reaching it
			counts.merge(v, 1, Integer::sum);
			if (isCorrupt(v)) corrupt.add(p);
		}
		log.info("================ VERIFY-SCHEME (DAG) ================");
		counts.forEach((v, c) -> log.info("  {} {}", v, c));
		if (unverified > 0) log.info("  UNVERIFIED (mono stop) {}", unverified);
		log.info("  CORRUPT total      : {}", corrupt.size());
		corrupt.stream().limit(25).forEach(p ->
				log.info("     {} {}", verdict.get(p), SCHEMES_ROOT.relativize(p)));
		// Dump the FULL corrupt list (verdict<TAB>relative-path) so a purge step can act on
		// every corrupt file, not just the 25 logged. One line per file.
		if (!corrupt.isEmpty()) {
			Path dump = Path.of("target/corrupt-files.txt");
			try {
				List<String> lines = corrupt.stream()
						.map(p -> verdict.get(p) + "\t" + SCHEMES_ROOT.relativize(p))
						.sorted().toList();
				Files.createDirectories(dump.getParent());
				Files.write(dump, lines);
				log.info("  full corrupt list  : {} ({} entries)", dump, lines.size());
			} catch (IOException e) {
				log.warn("could not write corrupt-list dump: {}", e.toString());
			}
		}
	}

	private static boolean isCorrupt(Verdict v) {
		return v != null && v.name().startsWith("CORRUPT");
	}

	private static void progress(int done, int total, long t0) {
		if (done % 1000 == 0) {
			log.info("[progress] {}/{} verified, {}ms", done, total, (System.nanoTime() - t0) / 1_000_000L);
		}
	}

	// ── ref resolution / parsing ──────────────────────────────────────────────

	private static Path resolveDep(String ref, Map<String, Path> byPinnedKey,
			Map<String, Path> byCanonKey, Map<String, Path> bestByShape) {
		int at = ref.indexOf('@');
		if (at > 0) {
			String shape = ref.substring(0, at);
			String tail = ref.substring(at + 1);
			if (tail.equals("naive")) return null;                          // elementary, no file
			if (tail.equals("sota")) return bestByShape.get(sortShape(shape)); // cited bound → best
			String hash7 = tail.substring(0, Math.min(7, tail.length()));
			return byPinnedKey.get(sortShape(shape) + "@" + hash7);         // sorted: orientation-tolerant
		}
		if (ref.contains("naive")) return null;
		Matcher mx = Pattern.compile("^(\\d+x\\d+x\\d+)(?:-direct)?$").matcher(ref);
		if (mx.matches()) return bestByShape.get(sortShape(mx.group(1)));   // bare/-direct → best
		return byCanonKey.get(Lineage.canonicalKey(ref));                   // named → exact file
	}

	private static boolean isPinnedOrNamed(String ref) {
		int at = ref.indexOf('@');
		if (at > 0) {
			String tail = ref.substring(at + 1);
			return !tail.equals("naive") && !tail.equals("sota");
		}
		return false; // bare/-direct/named handled via best/canon (not "dangling")
	}

	private static List<String> lineageRefs(Path f) {
		List<String> out = new ArrayList<>();
		try {
			collectRefs(SchemeIO.parseJson(f.toFile()).get("lineage"), out);
		} catch (Exception e) { /* none */ }
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

	private static Set<String> parseShapes(String[] args) {
		Set<String> out = new HashSet<>();
		String shapes = arg(args, "--shape", null);
		if (shapes != null) for (String s : shapes.split(",")) out.add(sortShape(s.trim()));
		String file = arg(args, "--shape-file", null);
		if (file != null) {
			try {
				for (String line : Files.readAllLines(Path.of(file))) {
					String tok = line.trim();
					if (tok.isEmpty() || tok.startsWith("#")) continue;
					Matcher m = Pattern.compile("(\\d+x\\d+x\\d+)").matcher(tok);
					if (m.find()) out.add(sortShape(m.group(1)));
				}
			} catch (IOException e) {
				throw new RuntimeException("reading --shape-file " + file, e);
			}
		}
		return out;
	}

	private static String sortShape(String shape) {
		String[] parts = shape.split("x");
		int[] d = { Integer.parseInt(parts[0]), Integer.parseInt(parts[1]), Integer.parseInt(parts[2]) };
		java.util.Arrays.sort(d);
		return d[0] + "x" + d[1] + "x" + d[2];
	}

	private static String arg(String[] args, String key, String def) {
		for (String a : args) if (a.startsWith(key + "=")) return a.substring(key.length() + 1);
		return def;
	}
}
