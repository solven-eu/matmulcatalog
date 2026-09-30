package eu.solven.matmul.docs.explore;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.SymmetryTransforms;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.AllocationOptimizer;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.search.LineageReplayer;
import eu.solven.matmul.search.SearchBudget;

import lombok.extern.slf4j.Slf4j;

/**
 * Issue #8 probe: can OUR engine re-derive the 22 FMM-Lille "recipe" ranks that
 * Marcos Adriano executed as explicit schemes (recombination over unequal blocks,
 * outer ⟨2,4,4⟩=26 HK or Strassen ⟨2,2,2⟩=7, every product at the smallest size
 * its two sides need)? Per target, prints:
 * <ol>
 *   <li>our catalog cost for every recipe piece vs the recipe's claimed piece rank
 *       (a piece gap = missing import, not an engine gap);</li>
 *   <li>the recipe sum re-costed with OUR leaves;</li>
 *   <li>the exact allocation optimizer's best rank over EVERY ⟨2,4,4⟩=26 (resp.
 *       ⟨2,2,2⟩=7) scheme in the catalog, in every orientation of that shape —
 *       different schemes at the same (shape, rank) have different product
 *       supports and tile differently (the hk71-vs-alphatensor ⟨2,4,4⟩ precedent);</li>
 *   <li>what the existing on-disk stub replays to (lineage leaves are hash-pinned,
 *       so a newer, better piece is NOT picked up by replay).</li>
 * </ol>
 * Read-only. {@code --no-replay} skips (4); {@code --shape=NxMxP} restricts.
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.explore.ProbeIssue8Recipes -Dexec.args="--field=Q"</pre>
 */
@Slf4j
public final class ProbeIssue8Recipes {
	private ProbeIssue8Recipes() {}

	record Target(int n, int m, int p, int recipeRank, String outer, int[] a, int[] b, int[] c, String recipe) {}

	private static Target t(int n, int m, int p, int r, String outer, int[] a, int[] b, int[] c, String recipe) {
		return new Target(n, m, p, r, outer, a, b, c, recipe);
	}

	/** The 22 recipes exactly as stated in the contributed files' lineage_str (issue #8). */
	static final List<Target> TARGETS = List.of(
			t(6, 14, 25, 1322, "2x4x4", new int[] { 3, 3 }, new int[] { 3, 3, 4, 4 }, new int[] { 6, 6, 6, 7 },
					"2 <3x3x7:49> + 6 <3x4x7:63> + 9 <3x3x6:40> + 9 <3x4x6:54>"),
			t(10, 19, 31, 3532, "2x4x4", new int[] { 5, 5 }, new int[] { 5, 5, 4, 5 }, new int[] { 8, 8, 8, 7 },
					"<4x5x7:104> + 4 <4x5x8:118> + 4 <5x5x7:127> + 17 <5x5x8:144>"),
			t(10, 22, 25, 3308, "2x4x4", new int[] { 5, 5 }, new int[] { 5, 5, 6, 6 }, new int[] { 6, 7, 6, 6 },
					"4 <5x5x7:127> + 4 <5x6x7:150> + 7 <5x5x6:110> + 11 <5x6x6:130>"),
			t(10, 22, 29, 3795, "2x4x4", new int[] { 5, 5 }, new int[] { 5, 5, 6, 6 }, new int[] { 7, 8, 7, 7 },
					"4 <5x5x8:144> + 4 <5x6x8:170> + 7 <5x5x7:127> + 11 <5x6x7:150>"),
			t(10, 23, 23, 3183, "2x4x4", new int[] { 5, 5 }, new int[] { 6, 6, 5, 6 }, new int[] { 6, 6, 6, 5 },
					"<5x5x5:93> + 8 <5x5x6:110> + 17 <5x6x6:130>"),
			t(10, 23, 26, 3571, "2x4x4", new int[] { 5, 5 }, new int[] { 6, 6, 5, 6 }, new int[] { 6, 6, 7, 7 },
					"2 <5x5x6:110> + 3 <5x5x7:127> + 9 <5x6x6:130> + 12 <5x6x7:150>"),
			t(10, 23, 27, 3688, "2x4x4", new int[] { 5, 5 }, new int[] { 6, 6, 5, 6 }, new int[] { 7, 7, 7, 6 },
					"<5x5x6:110> + 4 <5x5x7:127> + 4 <5x6x6:130> + 17 <5x6x7:150>"),
			t(10, 23, 30, 4076, "2x4x4", new int[] { 5, 5 }, new int[] { 6, 6, 5, 6 }, new int[] { 7, 7, 8, 8 },
					"2 <5x5x7:127> + 3 <5x5x8:144> + 9 <5x6x7:150> + 12 <5x6x8:170>"),
			t(10, 23, 31, 4193, "2x4x4", new int[] { 5, 5 }, new int[] { 6, 6, 5, 6 }, new int[] { 8, 8, 8, 7 },
					"<5x5x7:127> + 4 <5x5x8:144> + 4 <5x6x7:150> + 17 <5x6x8:170>"),
			t(10, 25, 26, 3864, "2x4x4", new int[] { 5, 5 }, new int[] { 6, 6, 6, 7 }, new int[] { 6, 7, 7, 6 },
					"4 <5x7x7:176> + 7 <5x6x6:130> + 15 <5x6x7:150>"),
			t(10, 26, 26, 4028, "2x4x4", new int[] { 5, 5 }, new int[] { 6, 6, 7, 7 }, new int[] { 7, 7, 6, 6 },
					"4 <5x6x6:130> + 8 <5x7x7:176> + 14 <5x6x7:150>"),
			t(10, 26, 27, 4172, "2x4x4", new int[] { 5, 5 }, new int[] { 6, 6, 7, 7 }, new int[] { 7, 7, 6, 7 },
					"2 <5x6x6:130> + 12 <5x6x7:150> + 12 <5x7x7:176>"),
			t(10, 26, 29, 4482, "2x4x4", new int[] { 5, 5 }, new int[] { 6, 6, 7, 7 }, new int[] { 7, 8, 7, 7 },
					"4 <5x6x8:170> + 4 <5x7x8:204> + 7 <5x6x7:150> + 11 <5x7x7:176>"),
			t(10, 27, 27, 4322, "2x4x4", new int[] { 5, 5 }, new int[] { 7, 7, 6, 7 }, new int[] { 7, 7, 7, 6 },
					"<5x6x6:130> + 8 <5x6x7:150> + 17 <5x7x7:176>"),
			t(13, 23, 31, 5396, "2x2x2", new int[] { 6, 7 }, new int[] { 11, 12 }, new int[] { 15, 16 },
					"<6x11x16:684> + <6x12x15:686> + <6x12x16:736> + <7x11x15:777> + <7x11x16:822> + <7x12x15:815> + <7x12x16:876>"),
			t(13, 23, 32, 5552, "2x2x2", new int[] { 6, 7 }, new int[] { 11, 12 }, new int[] { 16, 16 },
					"<6x11x16:684> + 2 <6x12x16:736> + 2 <7x11x16:822> + 2 <7x12x16:876>"),
			t(13, 24, 31, 5540, "2x2x2", new int[] { 6, 7 }, new int[] { 12, 12 }, new int[] { 15, 16 },
					"<6x12x15:686> + 2 <6x12x16:736> + 2 <7x12x15:815> + 2 <7x12x16:876>"),
			t(13, 25, 32, 6008, "2x2x2", new int[] { 6, 7 }, new int[] { 12, 13 }, new int[] { 16, 16 },
					"<6x12x16:736> + 2 <6x13x16:798> + 2 <7x12x16:876> + 2 <7x13x16:962>"),
			t(13, 29, 32, 6910, "2x2x2", new int[] { 6, 7 }, new int[] { 14, 15 }, new int[] { 16, 16 },
					"<6x14x16:864> + 2 <6x15x16:920> + 2 <7x14x16:1022> + 2 <7x15x16:1081>"),
			t(13, 30, 31, 6913, "2x2x2", new int[] { 6, 7 }, new int[] { 15, 15 }, new int[] { 16, 15 },
					"<6x15x16:920> + <7x15x15:1032> + 2 <6x15x15:859> + 3 <7x15x16:1081>"),
			t(13, 31, 32, 7322, "2x2x2", new int[] { 6, 7 }, new int[] { 15, 16 }, new int[] { 16, 16 },
					"<6x15x16:920> + 2 <6x16x16:972> + 2 <7x15x16:1081> + 2 <7x16x16:1148>"),
			t(15, 31, 32, 8185, "2x2x2", new int[] { 7, 8 }, new int[] { 16, 15 }, new int[] { 16, 16 },
					"<7x16x16:1148> + <8x15x16:1185> + 2 <7x15x16:1081> + 3 <8x16x16:1230>"));

	private static final Pattern PIECE = Pattern.compile("(?:(\\d+)\\s*)?<(\\d+)x(\\d+)x(\\d+):(\\d+)>");

	public static void main(String[] args) throws Exception {
		String fieldTag = "Q";
		boolean replay = true;
		String only = null;
		long maxNodes = 20_000_000L;
		for (String a : args) {
			if (a.startsWith("--field=")) fieldTag = a.substring("--field=".length());
			else if (a.equals("--no-replay")) replay = false;
			else if (a.startsWith("--shape=")) only = a.substring("--shape=".length());
			else if (a.startsWith("--max-nodes=")) maxNodes = Long.parseLong(a.substring("--max-nodes=".length()));
		}
		FieldAwareLookup lk = new FieldAwareLookup(fieldTag);
		SotaResolver sota = (p, q, r) -> {
			if (p == 0 || q == 0 || r == 0) return 0;
			if (p == 1) return q * r;
			if (q == 1) return p * r;
			if (r == 1) return p * q;
			return lk.findRank(p, q, r);
		};
		Map<String, List<NamedVariant>> basesByOuter = new LinkedHashMap<>();
		basesByOuter.put("2x4x4", orientedVariants(lk, 2, 4, 4, 26));
		basesByOuter.put("2x2x2", orientedVariants(lk, 2, 2, 2, 7));
		for (var e : basesByOuter.entrySet()) {
			log.info("outer {}: {} oriented base variant(s): {}", e.getKey(), e.getValue().size(),
					e.getValue().stream().map(NamedVariant::label).toList());
		}
		LineageReplayer replayer = replay ? LineageReplayer.withDefaultPool(lk) : null;

		List<String> summary = new ArrayList<>();
		for (Target t : TARGETS) {
			String shape = t.n() + "x" + t.m() + "x" + t.p();
			if (only != null && !only.equals(shape)) continue;
			long t0 = System.nanoTime();
			int catalog = lk.findRank(t.n(), t.m(), t.p());
			log.info("================ ⟨{},{},{}⟩ recipe={} catalog({})={} outer={} blocks A{} B{} C{}",
					t.n(), t.m(), t.p(), t.recipeRank(), fieldTag, catalog, t.outer(),
					java.util.Arrays.toString(t.a()), java.util.Arrays.toString(t.b()), java.util.Arrays.toString(t.c()));

			// (1)+(2) pieces: our leaf cost vs claimed; recipe sum with OUR leaves.
			long sumOurs = 0, sumClaimed = 0;
			int pieceGaps = 0;
			Matcher mt = PIECE.matcher(t.recipe());
			StringBuilder pieces = new StringBuilder();
			while (mt.find()) {
				int cnt = mt.group(1) == null ? 1 : Integer.parseInt(mt.group(1));
				int a = Integer.parseInt(mt.group(2)), b = Integer.parseInt(mt.group(3)), c = Integer.parseInt(mt.group(4));
				int claimed = Integer.parseInt(mt.group(5));
				int ours = sota.getRank(a, b, c);
				sumOurs += (long) cnt * ours;
				sumClaimed += (long) cnt * claimed;
				if (ours > claimed) pieceGaps++;
				pieces.append(String.format("%d×⟨%d,%d,%d⟩ ours=%d claimed=%d%s  ", cnt, a, b, c, ours, claimed,
						ours > claimed ? " ◄GAP" : (ours < claimed ? " (better)" : "")));
			}
			log.info("  pieces: {}", pieces);
			log.info("  recipe sum: claimed={} withOurLeaves={}{}", sumClaimed, sumOurs,
					pieceGaps > 0 ? "  (" + pieceGaps + " piece gap(s) — import first)" : "");

			// (3) exact optimizer over every base variant of the outer shape.
			long best = Long.MAX_VALUE;
			String bestLabel = null, bestAlloc = null;
			boolean bestExhaustive = false;
			for (NamedVariant nv : basesByOuter.get(t.outer())) {
				long s = System.nanoTime();
				SearchBudget budget = SearchBudget.upTo(catalog).withMaxNodes(maxNodes);
				AllocationOptimizer.Result r = AllocationOptimizer.optimize(nv.alg(), sota, t.n(), t.m(), t.p(), budget, null);
				String alloc = java.util.Arrays.toString(r.allocA()) + "|" + java.util.Arrays.toString(r.allocB()) + "|"
						+ java.util.Arrays.toString(r.allocC());
				log.info("    base {}: rank={} improved={} exhaustive={} nodes={}/{} alloc={} ({} ms)", nv.label(),
						r.rank(), r.improvedOnBound(), r.exhaustive(), r.nodes(), r.fullSpace(), alloc,
						(System.nanoTime() - s) / 1_000_000L);
				if (r.improvedOnBound() && r.rank() < best) {
					best = r.rank();
					bestLabel = nv.label();
					bestAlloc = alloc;
					bestExhaustive = r.exhaustive();
				}
			}

			// (4) what the on-disk stub replays to today.
			String replayed = "n/a";
			if (replayer != null) {
				try {
					for (Path f : lk.findFiles(t.n(), t.m(), t.p())) {
						var root = SchemeIO.parseJson(f.toFile());
						if (!SchemeIO.isStub(root)) continue;
						int stored = root.get("m").asInt();
						if (stored != catalog) continue;
						NonCubicBilinearAlgorithm re = replayer.replayFromFile(f.toFile());
						replayed = f.getFileName() + " → " + re.r;
						break;
					}
				} catch (Throwable ex) {
					replayed = "REPLAY FAILED: " + ex;
				}
			}
			String verdict = best == Long.MAX_VALUE ? "NO IMPROVEMENT over catalog " + catalog
					: (best <= t.recipeRank() ? "REPRODUCED " + best : "PARTIAL " + best + " (+" + (best - t.recipeRank()) + " vs recipe)");
			String line = String.format("⟨%d,%d,%d⟩ recipe=%d catalog=%d ourLeavesSum=%d | optimizer: %s via %s alloc=%s%s | stub replay: %s (%d ms)",
					t.n(), t.m(), t.p(), t.recipeRank(), catalog, sumOurs, verdict, bestLabel, bestAlloc,
					bestExhaustive ? "" : " [node cap]", replayed, (System.nanoTime() - t0) / 1_000_000L);
			log.info("  ==> {}", line);
			summary.add(line);
		}
		log.info("==================== SUMMARY ====================");
		for (String s : summary) log.info("{}", s);
	}

	record NamedVariant(String label, NonCubicBilinearAlgorithm alg) {}

	/** Every catalog scheme of shape ⟨n,m,p⟩ (any axis order) at exactly {@code rank}, in every
	 *  S₃ orientation that lands on ⟨n,m,p⟩, deduplicated by content hash. */
	private static List<NamedVariant> orientedVariants(FieldAwareLookup lk, int n, int m, int p, int rank) {
		List<NamedVariant> out = new ArrayList<>();
		java.util.Set<String> seen = new java.util.HashSet<>();
		for (Path f : lk.findFiles(n, m, p)) {
			try {
				var root = SchemeIO.parseJson(f.toFile());
				if (SchemeIO.isStub(root) || root.get("m").asInt() != rank) continue;
				NonCubicBilinearAlgorithm alg = SchemeIO.read(root);
				int idx = 0;
				for (NonCubicBilinearAlgorithm v : SymmetryTransforms.s3Orbit(alg)) {
					if (v.n != n || v.m != m || v.p != p) continue;
					String h = SchemeIO.contentHash(v);
					if (!seen.add(h)) continue;
					out.add(new NamedVariant(f.getFileName().toString().replace(".json", "") + "#" + (idx++), v));
				}
			} catch (Exception e) {
				log.warn("skip {}: {}", f.getFileName(), e.toString());
			}
		}
		return out;
	}
}
