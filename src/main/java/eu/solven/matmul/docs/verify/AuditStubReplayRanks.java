package eu.solven.matmul.docs.verify;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.search.LineageReplayer;

import lombok.extern.slf4j.Slf4j;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * READ-ONLY audit of stub <b>rank honesty</b>: for every lineage-only stub
 * (maxDim &gt; 16, no on-disk matrices), replay its lineage via
 * {@link LineageReplayer#replayFromFile} and compare the <i>replayed</i> rank
 * against the <i>recorded</i> rank in {@code docs/catalog.json}.
 *
 * <p>Motivation: the {@code SchemeSweep} ⟨4,19,20⟩ investigation (2026-06-15)
 * found stubs whose recorded rank is LOWER than what their lineage actually
 * rebuilds (e.g. a {@code Serendipitous} ⟨4,16,20⟩=822 that replays to 832),
 * which inflates every parent that consumes them and the FMM comparison. The
 * materialiser now FAILS LOUD on such a worse rebuild
 * ({@code RecursiveMaterialiser.assertRebuildNotWorse}); this audit SIZES the
 * pre-existing population so the purge can be scoped before it is run.</p>
 *
 * <p>Per-stub verdict:</p>
 * <ul>
 *   <li><b>OVER_CLAIM</b> — replayed rank &gt; recorded rank. The recorded rank is
 *       a phantom; the stub (and the manifest) over-state the result. The purge set.</li>
 *   <li><b>EXACT</b> — replayed == recorded. Honest.</li>
 *   <li><b>BETTER</b> — replayed &lt; recorded (projection DCE). Honest-but-pessimistic;
 *       the catalog could carry the improved rank.</li>
 *   <li><b>REPLAY_FAIL</b> — replay threw / cycle / dangling ref. A different bucket
 *       (see {@link AuditLineageRefs} / {@link DetectCyclicStubs}).</li>
 * </ul>
 *
 * <p>Writes NOTHING. Args: {@code --field=Q} (default), {@code --sample=N}
 * (stride-sample 1 of every N stubs; default 1 = all), {@code --limit=N} (stop
 * after N audited; default unbounded).</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.verify.AuditStubReplayRanks \
 *     -Dexec.args="--field=Q --sample=20"</pre>
 */
@Slf4j
public final class AuditStubReplayRanks {

	private AuditStubReplayRanks() {}

	private static final String SCHEMES_ROOT = "src/main/resources/schemes";

	enum Verdict { OVER_CLAIM, EXACT, BETTER, REPLAY_FAIL }

	public static void main(String[] args) throws Exception {
		String field = arg(args, "--field", "Q");
		int sample = Integer.parseInt(arg(args, "--sample", "1"));
		int limit = Integer.parseInt(arg(args, "--limit", String.valueOf(Integer.MAX_VALUE)));

		JsonMapper mapper = JsonMapper.builder().build();
		JsonNode catalog = mapper.readTree(Files.readString(Path.of("docs/catalog.json")));
		JsonNode schemes = catalog.isArray() ? catalog : catalog.get("schemes");

		// Best (lowest) recorded rank per sorted shape — to flag over-claimers that are
		// the catalog-BEST for their shape (those directly corrupt the FMM comparison).
		Map<String, Integer> bestByShape = new java.util.HashMap<>();
		for (JsonNode s : schemes) {
			JsonNode f = s.get("format");
			if (f == null || f.size() != 3 || s.get("rank") == null) continue;
			int[] d = { f.get(0).asInt(), f.get(1).asInt(), f.get(2).asInt() };
			java.util.Arrays.sort(d);
			String key = d[0] + "x" + d[1] + "x" + d[2];
			bestByShape.merge(key, s.get("rank").asInt(), Math::min);
		}

		// Collect stub entries (maxDim>16, has a file).
		List<JsonNode> stubs = new ArrayList<>();
		for (JsonNode s : schemes) {
			JsonNode f = s.get("format");
			if (f == null || f.size() != 3 || s.get("rank") == null || s.get("file") == null) continue;
			int maxDim = Math.max(f.get(0).asInt(), Math.max(f.get(1).asInt(), f.get(2).asInt()));
			if (maxDim > 16) stubs.add(s);
		}
		log.info("Stub population (maxDim>16, recorded rank + file): {}", stubs.size());

		FieldAwareLookup lookup = new FieldAwareLookup(field);
		LineageReplayer replayer = LineageReplayer.withDefaultPool(lookup);

		Map<Verdict, Integer> counts = new LinkedHashMap<>();
		for (Verdict v : Verdict.values()) counts.put(v, 0);
		Map<String, int[]> byOp = new TreeMap<>();        // op -> {overClaim, total}
		Map<Integer, Integer> magnitudeHist = new TreeMap<>(); // delta -> #stubs
		List<String> worst = new ArrayList<>();           // "shape recorded->replayed (+d) [op] BEST?"
		int overClaimBest = 0;

		long t0 = System.nanoTime();
		int audited = 0;
		for (int i = 0; i < stubs.size() && audited < limit; i++) {
			if (sample > 1 && (i % sample) != 0) continue;
			JsonNode s = stubs.get(i);
			JsonNode f = s.get("format");
			int recorded = s.get("rank").asInt();
			int[] d = { f.get(0).asInt(), f.get(1).asInt(), f.get(2).asInt() };
			java.util.Arrays.sort(d);
			String shapeKey = d[0] + "x" + d[1] + "x" + d[2];
			File file = new File(SCHEMES_ROOT, s.get("file").asString());
			String op = topOp(mapper, file);

			Verdict v;
			int replayed = -1;
			try {
				NonCubicBilinearAlgorithm rep = replayer.replayFromFile(file);
				replayed = rep.r;
				if (replayed > recorded) v = Verdict.OVER_CLAIM;
				else if (replayed == recorded) v = Verdict.EXACT;
				else v = Verdict.BETTER;
			} catch (Throwable e) {
				v = Verdict.REPLAY_FAIL;
			}
			counts.merge(v, 1, Integer::sum);
			byOp.computeIfAbsent(op, k -> new int[2])[1]++;
			if (v == Verdict.OVER_CLAIM) {
				byOp.get(op)[0]++;
				int delta = replayed - recorded;
				magnitudeHist.merge(delta, 1, Integer::sum);
				boolean isBest = Integer.valueOf(recorded).equals(bestByShape.get(shapeKey));
				if (isBest) overClaimBest++;
				worst.add(String.format("⟨%d,%d,%d⟩ %d→%d (+%d) [%s]%s",
						d[0], d[1], d[2], recorded, replayed, delta, op, isBest ? " BEST" : ""));
			}
			audited++;
			if (audited % 50 == 0) {
				long ms = (System.nanoTime() - t0) / 1_000_000L;
				log.info("[progress] {} audited ({} over-claim / {} exact / {} better / {} fail) {}ms",
						audited, counts.get(Verdict.OVER_CLAIM), counts.get(Verdict.EXACT),
						counts.get(Verdict.BETTER), counts.get(Verdict.REPLAY_FAIL), ms);
			}
		}

		log.info("================ STUB RANK-HONESTY AUDIT ({}; sample=1/{}; audited={}) ================",
				field, sample, audited);
		counts.forEach((v, c) -> log.info("  {} {}", v, c));
		log.info("  OVER_CLAIM that are catalog-BEST for their shape (corrupt the comparison): {}", overClaimBest);
		log.info("  -- over-claim by top-level op (overClaim/total audited of that op) --");
		byOp.forEach((op, a) -> { if (a[0] > 0) log.info("     {} {}/{}", op, a[0], a[1]); });
		log.info("  -- over-claim magnitude histogram (delta -> #stubs) --");
		magnitudeHist.forEach((delta, c) -> log.info("     +{} {}", delta, c));
		worst.sort(Comparator.comparingInt(AuditStubReplayRanks::deltaOf).reversed());
		log.info("  -- worst 25 over-claimers --");
		worst.stream().limit(25).forEach(w -> log.info("     {}", w));
		if (sample > 1) {
			log.info("  NOTE: sampled 1/{} — multiply counts by ~{} for a full-catalog estimate.", sample, sample);
		}
	}

	/** Top-level lineage op of a scheme file (e.g. {@code Serendipitous}, {@code SumInner},
	 *  {@code Recombination}, {@code Project}), or {@code ?} if unreadable. */
	private static String topOp(JsonMapper mapper, File f) {
		try {
			JsonNode root = mapper.readTree(Files.readString(f.toPath()));
			JsonNode lin = root.get("lineage");
			if (lin != null && lin.get("op") != null) return lin.get("op").asString();
			String str = root.path("lineage_str").asString("");
			int paren = str.indexOf('(');
			if (paren > 0) return str.substring(0, paren);
		} catch (Exception e) {
			// fall through
		}
		return "?";
	}

	private static int deltaOf(String worstLine) {
		int a = worstLine.indexOf("(+");
		int b = worstLine.indexOf(')', a);
		try {
			return Integer.parseInt(worstLine.substring(a + 2, b));
		} catch (RuntimeException e) {
			return 0;
		}
	}

	private static String arg(String[] args, String key, String def) {
		for (String a : args) {
			if (a.startsWith(key + "=")) return a.substring(key.length() + 1);
		}
		return def;
	}
}
