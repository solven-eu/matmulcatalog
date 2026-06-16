package eu.solven.matmul.docs.generate;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.Semaphore;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Stream;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ObjectNode;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeAnalysis;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.LineageReplayer;
import eu.solven.matmul.util.ProgressMonitor;
import lombok.extern.slf4j.Slf4j;

/**
 * Phase-2 enrich (#68/#190): expand every bilinear scheme ONCE and stamp the
 * derived metrics computed by {@link SchemeAnalysis#defaults()} —
 * {@code verified, additions, has_buds/buds, projection_margin} — into the
 * scheme JSON body, so the catalog manifest can FORWARD them (no re-expansion)
 * and the filename no longer has to carry them.
 *
 * <p>This is the unified single-replay procedure: μ and buds (and verify, adds)
 * share one expansion. Stubs (maxDim&gt;16) are replayed; explicit-matrix
 * schemes are read directly. Heavy schemes are throttled to bound peak heap.</p>
 *
 * <p>Dry-run by default (reports coverage + would-change counts, writes
 * nothing); pass {@code --apply=true} to stamp. {@code --max-dim=N} bounds the
 * run to small shapes for a quick pass.</p>
 *
 * <p>Structure: {@link #main} only parses the CLI and reports; the enrichment
 * engine is {@link #enrich(Options)} (callable from tests/tools), which returns
 * a {@link Result}. Per-file work is {@link #enrichOne}.</p>
 *
 * <pre>{@code
 * mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.generate.EnrichSchemeMetrics \
 *     -Dexec.args="--apply=false --max-dim=12"
 * }</pre>
 */
@Slf4j
public final class EnrichSchemeMetrics {
	private EnrichSchemeMetrics() {}

	private static final Path ROOT = Path.of("src/main/resources/schemes");
	/** Self-heal keys: cleared when a previously-corrupted stub replays successfully. */
	private static final List<String> CORRUPT_KEYS = List.of("corrupted", "corrupted_reason");
	private static final JsonMapper MAPPER = JsonMapper.builder().build();
	/** Weight (r·(nm+mp+np)) above which a scheme is throttled (heavy expansion). */
	private static final long HEAVY_WEIGHT = 10_000_000L;

	/** Parsed CLI options. */
	public record Options(boolean apply, int maxDim, boolean revalidateStubs) {
		public static Options defaults() {
			return new Options(false, Integer.MAX_VALUE, false);
		}
	}

	/** Aggregated outcome of an enrichment run. */
	public record Result(int stamped, int normalized, int unchanged, int skipped,
			int errors, int corruptedStamped, int healed) {}

	// ------------------------------------------------------------------ CLI

	public static void main(String[] args) throws Exception {
		Options opts = parseArgs(args);
		Result r = enrich(opts);
		log.info("=== Enrich {}: {} {}, {} {} (formatting only), {} already-current, "
				+ "{} skipped (non-bilinear/maxDim), {} errors ({} flagged corrupted, {} healed) ===",
				opts.apply() ? "APPLIED" : "DRY-RUN",
				r.stamped(), opts.apply() ? "stamped" : "would-stamp",
				r.normalized(), opts.apply() ? "normalised" : "would-normalise",
				r.unchanged(), r.skipped(), r.errors(),
				r.corruptedStamped(), r.healed());
	}

	/** Parse the CLI flags into {@link Options}. The only CLI-aware code path. */
	static Options parseArgs(String[] args) {
		boolean apply = false;
		boolean revalidateStubs = false;
		int maxDim = Integer.MAX_VALUE;
		for (String a : args) {
			if (a.startsWith("--apply=")) {
				apply = Boolean.parseBoolean(a.substring("--apply=".length()));
			} else if (a.startsWith("--max-dim=")) {
				maxDim = Integer.parseInt(a.substring("--max-dim=".length()));
			} else if (a.equals("--revalidate-stubs")) {
				// Re-replay every STUB even if it already carries projection_margin, so a
				// stub that has BECOME unreplayable (its lineage now reconstructs a
				// different shape because a dependency changed) is re-flagged corrupted —
				// the metric sentinel can't see that a once-good stub rotted. Explicit-matrix
				// files keep the fast path (no replay risk). Pair with --apply=true to stamp.
				revalidateStubs = true;
			} else {
				throw new IllegalArgumentException("unknown arg " + a);
			}
		}
		return new Options(apply, maxDim, revalidateStubs);
	}

	// ---------------------------------------------------------------- engine

	/** Mutable per-run counters (thread-safe), folded into a {@link Result} at the end. */
	private static final class Counters {
		final AtomicInteger stamped = new AtomicInteger();
		final AtomicInteger normalized = new AtomicInteger();
		final AtomicInteger unchanged = new AtomicInteger();
		final AtomicInteger skipped = new AtomicInteger();
		final AtomicInteger errors = new AtomicInteger();
		final AtomicInteger corruptedStamped = new AtomicInteger();
		final AtomicInteger healed = new AtomicInteger();

		Result toResult() {
			return new Result(stamped.get(), normalized.get(), unchanged.get(), skipped.get(),
					errors.get(), corruptedStamped.get(), healed.get());
		}
	}

	/**
	 * Run the enrichment over the whole catalog under {@code opts} and return the
	 * aggregated counts. No CLI parsing and no summary line — pure engine, so it
	 * can be driven from a test or another tool.
	 */
	public static Result enrich(Options opts) throws Exception {
		List<Path> files;
		try (Stream<Path> w = Files.walk(ROOT)) {
			files = w.filter(p -> p.getFileName().toString().endsWith(".json")).sorted().toList();
		}
		int nThreads = Math.max(1, Runtime.getRuntime().availableProcessors() - 1);
		Semaphore heavyGate = new Semaphore(Math.max(1, nThreads / 4));
		ExecutorService pool = Executors.newFixedThreadPool(nThreads);
		FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);
		// One shared replayer across all worker threads (it is documented thread-safe
		// and its named-base pool is immutable). Building one per file raced on shared
		// pool construction → sporadic ConcurrentModificationException.
		LineageReplayer replayer = LineageReplayer.withDefaultPool(lookup);
		List<SchemeAnalysis> analyses = SchemeAnalysis.defaults();
		Counters c = new Counters();

		ProgressMonitor monitor = ProgressMonitor.builder().total(files.size())
				.logEveryMillis(5_000).logger(log).label("enrich-metrics").build();
		monitor.start();
		log.info("Enrich {} scheme files (apply={}, maxDim≤{}) on {} threads",
				files.size(), opts.apply(), opts.maxDim(), nThreads);

		List<Future<?>> fs = new ArrayList<>();
		for (int i = 0; i < files.size(); i++) {
			Path p = files.get(i);
			int idx = i;
			fs.add(pool.submit(() -> enrichOne(p, idx, opts, replayer, analyses, c, heavyGate, monitor)));
		}
		for (Future<?> f : fs) {
			try {
				f.get();
			} catch (Exception ignored) {
				// per-file failures are already counted/logged inside enrichOne
			}
		}
		pool.shutdown();
		monitor.done();
		return c.toResult();
	}

	/** Enrich a single scheme file (one worker task). Mutates {@code c}; never throws. */
	private static void enrichOne(Path p, int idx, Options opts, LineageReplayer replayer,
			List<SchemeAnalysis> analyses, Counters c, Semaphore heavyGate, ProgressMonitor monitor) {
		boolean heavy = false;
		boolean wasCorrupted = false;
		boolean isStubFile = false;
		try {
			String existing = Files.readString(p);
			ObjectNode root = (ObjectNode) MAPPER.readTree(existing);
			wasCorrupted = SchemeIO.isCorrupted(root);
			isStubFile = SchemeIO.isStub(root);
			int[] dims = shapeOf(p.getFileName().toString());
			if (dims == null || Math.max(dims[0], Math.max(dims[1], dims[2])) > opts.maxDim()) {
				c.skipped.incrementAndGet();
				return;
			}
			if (SchemeIO.isNonBilinear(root) || SchemeIO.isComplex(root)) {
				c.skipped.incrementAndGet();
				return;
			}
			// Fast path: a file already carrying the full default-metric set
			// (sentinel: projection_margin — all metrics are stamped together)
			// needs NO re-expansion. Run an empty-field addFields: it re-emits
			// through the canonical MatrixJsonFormatter, so a file whose
			// formatting had drifted (e.g. an older 1-space hand-appended field)
			// is normalised, while an already-canonical file is a pure no-op.
			if (root.has("projection_margin") && !(opts.revalidateStubs() && isStubFile)) {
				// Reuse the node + text we already read — no second read/parse.
				if (SchemeIO.addFields(p.toFile(), root, existing, Map.of(), opts.apply())) {
					c.normalized.incrementAndGet();
				} else {
					c.unchanged.incrementAndGet();
				}
				return;
			}
			heavy = weight(dims, root) > HEAVY_WEIGHT;
			if (heavy) heavyGate.acquireUninterruptibly();
			NonCubicBilinearAlgorithm alg = expand(p, root, replayer);
			Field field = fieldOf(root);
			// Collect the derived metrics ABSENT from the JSON; SchemeIO.addFields
			// adds them and re-emits the whole file via the canonical
			// MatrixJsonFormatter (THE single writer — never a hand-built string),
			// so on-disk style stays consistent with every other writer.
			java.util.LinkedHashMap<String, Object> toAdd = new java.util.LinkedHashMap<>();
			for (SchemeAnalysis an : analyses) {
				for (Map.Entry<String, Object> e : an.analyse(alg, field).entrySet()) {
					if (!root.has(e.getKey())) toAdd.put(e.getKey(), e.getValue());
				}
			}
			// Success: also clear any stale corrupted flag (self-heal — a stub
			// that now replays is no longer corrupted). Reuse the parsed node + text.
			boolean changed = SchemeIO.updateFields(p.toFile(), root, existing, toAdd, CORRUPT_KEYS, opts.apply());
			if (wasCorrupted) c.healed.incrementAndGet();
			if (!toAdd.isEmpty()) {
				c.stamped.incrementAndGet();
			} else if (changed) {
				c.normalized.incrementAndGet();
			} else {
				c.unchanged.incrementAndGet();
			}
		} catch (Throwable t) {
			c.errors.incrementAndGet();
			log.warn("enrich failed {}: {}", p.getFileName(), t.toString());
			// A stub whose recipe can't be replayed into matrices is flagged
			// corrupted so the search-gating lookup treats it as ABSENT (it must
			// not shadow a fresh, verifiable discovery). Non-stub read failures
			// are not flagged (their matrices are on disk; the failure is elsewhere).
			if (isStubFile && !wasCorrupted) {
				String reason = t.toString();
				if (reason.length() > 180) reason = reason.substring(0, 180);
				try {
					if (SchemeIO.updateFields(p.toFile(),
							Map.of("corrupted", true, "corrupted_reason", reason),
							List.of(), opts.apply())) {
						c.corruptedStamped.incrementAndGet();
					}
				} catch (Exception ex) {
					log.warn("could not flag corrupted {}: {}", p.getFileName(), ex.toString());
				}
			}
		} finally {
			if (heavy) heavyGate.release();
			monitor.tick(idx);
		}
	}

	// --------------------------------------------------------------- helpers

	/** Expand to explicit matrices: replay stubs, read explicit otherwise. */
	private static NonCubicBilinearAlgorithm expand(Path p, JsonNode root, LineageReplayer replayer) {
		if (SchemeIO.isStub(root)) {
			return replayer.replayFromFile(p.toFile());
		}
		try {
			return SchemeIO.isReduced(root) ? SchemeIO.readReduced(root) : SchemeIO.read(root);
		} catch (Exception e) {
			throw new RuntimeException("read " + p.getFileName(), e);
		}
	}

	private static int[] shapeOf(String name) {
		var m = java.util.regex.Pattern.compile("(\\d+)x(\\d+)x(\\d+)").matcher(name);
		if (!m.find()) return null;
		return new int[] { Integer.parseInt(m.group(1)), Integer.parseInt(m.group(2)),
				Integer.parseInt(m.group(3)) };
	}

	private static long weight(int[] d, JsonNode root) {
		long n = d[0], mm = d[1], p = d[2];
		long r = root.has("m") && root.get("m").isInt() ? root.get("m").asInt() : n * mm * p;
		return r * (n * mm + mm * p + n * p);
	}

	private static Field fieldOf(JsonNode root) {
		List<String> tags = SchemeIO.fieldTags(root);
		for (String t : new String[] { "Z", "Q", "R", "C", "F3", "F2" }) {
			if (tags.contains(t)) return Field.fromTag(t);
		}
		return Field.Q;
	}
}
