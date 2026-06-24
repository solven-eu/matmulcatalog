package eu.solven.matmul.docs.explore;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.search.CitedBound;
import eu.solven.matmul.search.PoolConfig;
import lombok.extern.slf4j.Slf4j;

/**
 * (a) — Mechanism map of the "behind" band. For every FMM-sourced format where
 * our committed best &gt; FMM, re-run {@link BlockSplitSearch#findBestStrategy}
 * with the UNBALANCED {@code includeDerived} pool and classify:
 *
 * <ul>
 *   <li><b>CLOSEABLE</b>: predicted ≤ FMM — the unbalanced recombination now
 *       matches/beats FMM (re-materialise to capture; this is bucket (0)).</li>
 *   <li><b>PARTIAL</b>: predicted improves but still &gt; FMM.</li>
 *   <li><b>HOLDOUT</b>: no improvement &amp; still &gt; FMM — the rotation-fuse
 *       (TA) candidates that need mechanism (b).</li>
 * </ul>
 *
 * Reads {@code docs/catalog.json} for the external baseline, runs in parallel,
 * and writes the CLOSEABLE shape list to {@code target/behind-closeable.txt} so
 * (0) can re-materialise exactly those.
 */
@Slf4j
public class MapBehindBand {

	record Behind(int n, int m, int p, long ourBest, long ext) {}

	public static void main(String[] args) throws Exception {
		JsonMapper mapper = JsonMapper.builder().build();
		JsonNode root = mapper.readTree(new File("docs/catalog.json"));
		// Per canonical format: our best COMPARABLE rank + external_best (FMM-sourced only).
		java.util.Map<String, long[]> best = new java.util.HashMap<>(); // key -> [our,ext]
		for (JsonNode s : root.get("schemes")) {
			if (!s.has("external_best_rank")) continue;
			if (!"fmm-lille".equals(s.path("external_best_source").asString(""))) continue;
			JsonNode fmt = s.get("format");
			int[] d = { fmt.get(0).asInt(), fmt.get(1).asInt(), fmt.get(2).asInt() };
			java.util.Arrays.sort(d);
			String key = d[0] + "x" + d[1] + "x" + d[2];
			long r = s.get("rank").asLong();
			long ext = s.get("external_best_rank").asLong();
			long[] cur = best.get(key);
			if (cur == null) best.put(key, new long[] { r, ext });
			else cur[0] = Math.min(cur[0], r);
		}
		List<Behind> behind = new ArrayList<>();
		for (var e : best.entrySet()) {
			long[] v = e.getValue();
			if (v[0] > v[1]) {
				String[] t = e.getKey().split("x");
				behind.add(new Behind(Integer.parseInt(t[0]), Integer.parseInt(t[1]),
						Integer.parseInt(t[2]), v[0], v[1]));
			}
		}
		log.info("FMM-sourced behind formats: {}", behind.size());

		FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);
		CitedBound sota = new CitedBound(lookup);
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.buildPool(PoolConfig.includeDerived());
		log.info("pool={} (includeDerived, unbalanced)", pool.size());

		AtomicInteger done = new AtomicInteger();
		AtomicInteger closeable = new AtomicInteger(), partial = new AtomicInteger(), holdout = new AtomicInteger();
		ConcurrentLinkedQueue<String> closeableShapes = new ConcurrentLinkedQueue<>();
		ConcurrentLinkedQueue<String> holdoutRows = new ConcurrentLinkedQueue<>();
		long t0 = System.currentTimeMillis();
		int nThreads = Math.max(1, Runtime.getRuntime().availableProcessors() - 1);
		var exec = Executors.newFixedThreadPool(nThreads);
		List<java.util.concurrent.Future<?>> futs = new ArrayList<>();
		for (Behind b : behind) {
			futs.add(exec.submit(() -> {
				Optional<BlockSplitSearch.NonCubicStrategy> best2 = BlockSplitSearch.findBestStrategy(
						b.n, b.m, b.p, pool, sota, false, PoolConfig.UNBOUNDED_IMBALANCE,
						PoolConfig.UNBOUNDED_COMBINATIONS, 0, Long.MAX_VALUE);
				long pred = best2.map(BlockSplitSearch.NonCubicStrategy::rank).orElse(-1L);
				if (pred > 0 && pred <= b.ext) {
					closeable.incrementAndGet();
					closeableShapes.add(b.n + "x" + b.m + "x" + b.p);
				} else if (pred > 0 && pred < b.ourBest) {
					partial.incrementAndGet();
				} else {
					holdout.incrementAndGet();
					holdoutRows.add(String.format("%dx%dx%d our=%d fmm=%d pred=%d [%s]",
							b.n, b.m, b.p, b.ourBest, b.ext, pred,
							best2.map(BlockSplitSearch.NonCubicStrategy::label).orElse("none")));
				}
				int n = done.incrementAndGet();
				if (n % 50 == 0 || n == behind.size()) {
					log.info("[progress] {}/{} (closeable={} partial={} holdout={}) {}ms",
							n, behind.size(), closeable.get(), partial.get(), holdout.get(),
							System.currentTimeMillis() - t0);
				}
			}));
		}
		for (var f : futs) f.get();
		exec.shutdown();

		log.info("=== MECHANISM MAP (FMM-derived behind, unbalanced includeDerived) ===");
		log.info("CLOSEABLE (pred ≤ FMM, re-materialise): {}", closeable.get());
		log.info("PARTIAL   (improved, still > FMM):      {}", partial.get());
		log.info("HOLDOUT   (no improvement, TA candidates): {}", holdout.get());

		Path out = Path.of("target/behind-closeable.txt");
		Files.writeString(out, closeableShapes.stream().sorted().collect(Collectors.joining("\n")));
		log.info("wrote {} closeable shapes → {}", closeable.get(), out);
		Path hout = Path.of("target/behind-holdouts.txt");
		Files.writeString(hout, holdoutRows.stream().sorted().collect(Collectors.joining("\n")));
		log.info("wrote {} holdout rows → {}", holdout.get(), hout);
	}
}
